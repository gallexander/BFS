#include "project.h"

double mytime(void){
    struct timeval now;
    gettimeofday(&now,NULL);
    return (double) ((long long)now.tv_usec+(long long)now.tv_sec*1000000);
}

int main(){
    uint64_t nodes = pow(2,SCALE);
    uint64_t edges = nodes*EDGEFACTOR;
    uint64_t root = ROOT;
    uint64_t *startVertex = (uint64_t *) calloc(edges, I64_BYTES);
    uint64_t *endVertex = (uint64_t *) calloc(edges, I64_BYTES);
    uint64_t *parents = (uint64_t *) calloc(nodes, I64_BYTES);
    unsigned char *level = (unsigned char *) calloc(1, nodes / BITS); //LEVEL BUFFER FOR MASTER
    float initiator[] = {0.25,0.25,0.25,0.25};
    
    struct edge **node_edge_list = (struct edge **) calloc(nodes, 8);
	uint64_t *count_edges_per_node = (uint64_t *) calloc(nodes, I64_BYTES);

    

    generate_graph(SCALE, EDGEFACTOR, initiator, startVertex, endVertex);
    create_node_edge_lists(nodes, edges, startVertex, endVertex, node_edge_list, count_edges_per_node);
	//MAYBE CREATING BUFFERS WITH REALLOC, SO THERE ARE NOT LISTS, create_node_edge_lists, NECESSARY

	uint64_t buffer_size = 0;
	uint64_t *buffer = lists_to_buffer(&buffer_size, node_edge_list, count_edges_per_node, 0, nodes-1);

	//SCATTER WHOLE BUFFER TO CHILD PROCESSES

	//GET_EDGE_BUFFER FOR EACH PROC
	//GET_COUNT_BUFFER FOR EACH PROC

	//TRANSFER COUNT_BUFFER TO INDEX_BUFFER FOR EACH PROC
	uint64_t i;
	uint64_t prev = count_edges_per_node[0];
	uint64_t tmp;
	count_edges_per_node[0] = 0;
	for (i = 1; i < nodes; i++){
		tmp = count_edges_per_node[i];
		count_edges_per_node[i] = prev + count_edges_per_node[i-1];
		prev = tmp;
	}
	
    //SET ROOT LEVEL
    level[(ROOT/BITS)] = level[(ROOT/BITS)] | (unsigned char) pow(2,(ROOT % BITS));
    //printf("level: %i\n", level[(ROOT/BITS)]);

    //SCATTER LEVEL BUFFER

    //BFS
    double time = mytime();
    bfs(level, buffer, buffer_size, count_edges_per_node, nodes, 0);

	//OUTPUT
    /*uint64_t j;
    struct edge *p;
    for (j = 0; j < nodes; j++){
        printf("%llu:%llu Elements:", (unsigned long long) j, (unsigned long long) count_edges_per_node[j]);
        p = node_edge_list[j];
        if (p != NULL){
            printf(" %llu", (unsigned long long) p->end);
            while (p->next != NULL){
                p = p->next;
                printf(" %llu", (unsigned long long) p->end);
            }
        }
        printf("\n");
    }
	printf("\nBuffer: ");
	for (j = 0; j < buffer_size; j++){
		printf("%llu,", (unsigned long long) buffer[j]);
	}
	printf("\n\n");

	for (j = 0; j < nodes; j++){
		printf("%llu,", (unsigned long long) count_edges_per_node[j]);
	}
	printf("\n");*/

    time = mytime() - time;
    printf("Time: %f\n", time/1000000);

    free(startVertex);
    free(endVertex);
    free(parents);
	free(buffer);
    free(level);
	free(count_edges_per_node);
    freelists(nodes, node_edge_list);
    return 0;
}

void bfs(unsigned char *level, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, uint64_t nodes_owned, int proc){
    unsigned char *next_level = (unsigned char *) calloc(pow(2,SCALE) / BITS, 1);
    unsigned char *visited = (unsigned char *) calloc(nodes_owned / BITS, 1);
    char oneChildisVisited = 1;
    uint64_t i;
    unsigned char position;
    while (oneChildisVisited){
        oneChildisVisited = 0;
        //printf("Visit this round:");
        for (i = 0; i < nodes_owned; i++){
            position = (unsigned char) pow(2,(i % BITS));
            if (position & level[(i / BITS)] & ~visited[(i / BITS)]){
                visited[(i / BITS)] = visited[(i / BITS)] | position;
                //printf("%llu,", (unsigned long long) i);
                uint64_t j = index_of_node[i];
                if (i == nodes_owned -1){
                    for (; j < buffer_size; j++){
                        next_level[(buffer[j]/BITS)] = next_level[(buffer[j]/BITS)] | (unsigned char) pow(2,(buffer[j] % BITS));
                        oneChildisVisited = 1;
                    }
                }else{
                    for (; j < buffer_size && j < index_of_node[i+1]; j++){
                        next_level[(buffer[j]/BITS)] = next_level[(buffer[j]/BITS)] | (unsigned char) pow(2,(buffer[j] % BITS));
                        oneChildisVisited = 1;
                    }
                }
            }
        }
        //printf("\n");
        /*printf("Next level:\n");
        for (i = 0; i < nodes_owned; i++){
            if (((unsigned char) pow(2,(i % BITS))) & next_level[(i / BITS)]){
                printf("%llu,",(unsigned long long) i);
            }
        }
        printf("\n");*/

        // SEND MESSAGE THAT THERE ARE CHILDS TO EVALUATE, CAN BE ONE BYTE FROM ALL PROCS
        // AFTER SEND LEVEL BUFFER, ALLTOALL

        // NEXT STEPS ARE JUST FOR SEQ ALGO
        memcpy((void *)level,(void *)next_level,(nodes_owned/BITS));
        memset(next_level, 0, (nodes_owned/BITS));
    }

    free(visited);
    free(next_level);
}


void create_node_edge_lists(uint64_t nodes, uint64_t edges, uint64_t *startVertex, uint64_t *endVertex, struct edge **node_edge_list, uint64_t *count_edges_per_node){
    uint64_t i;
    char notAlreadyExists = 1;
    struct edge *p;
    for (i = 0; i < edges; i++){
        if (startVertex[i] != endVertex[i]){
            if (node_edge_list[startVertex[i]] == NULL){
                node_edge_list[startVertex[i]] = (struct edge *) calloc(1, sizeof(struct edge));
                node_edge_list[startVertex[i]]->end = endVertex[i];
				count_edges_per_node[startVertex[i]] = 1;
            }else{
                notAlreadyExists = 1;
                p = node_edge_list[startVertex[i]];
                if (p->end == endVertex[i]){
                    notAlreadyExists = 0;
                }
                while (p->next != NULL && notAlreadyExists){
                    p = p->next;
                    if (p->end == endVertex[i]){
                        notAlreadyExists = 0;
                    }
                }
                if (notAlreadyExists){
                    p->next = (struct edge *) calloc(1, sizeof(struct edge));
                    (p->next)->end = endVertex[i];
					count_edges_per_node[startVertex[i]] +=1;
                }
            }
        }
    }
}

uint64_t calculate_size(uint64_t *count_edges_per_node, uint64_t first, uint64_t last){
	uint64_t size;

	uint64_t i;
	for (i = first; i <= last; i++){
		size += count_edges_per_node[i];
	}
	return size;
}

uint64_t *lists_to_buffer(uint64_t *size, struct edge **node_edge_list, uint64_t *count_edges_per_node, uint64_t first, uint64_t last){
	struct edge *p = NULL;
	uint64_t *buffer = NULL;
	
	uint64_t i,j;
	for (i = first; i <= last; i++){
		(*size) += count_edges_per_node[i];
	}
	if ((*size)){
		buffer = (uint64_t *) calloc((*size), I64_BYTES);
		j = 0;
		for (i = first; i <= last; i++){
			if (node_edge_list[i] != NULL){
				p = node_edge_list[i];
				buffer[j++] = p->end;				
				while (p->next != NULL){
					p = p->next;
					buffer[j++] = p->end;
				}	
			}
		}
	}
	return buffer;
}

void freelists(uint64_t nodes, struct edge **node_edge_list){
    uint64_t i;
    for (i = 0; i < nodes; i++){
        freelist(node_edge_list[i]);
    }
    free(node_edge_list);
}

void freelist(struct edge *pp){
    if (pp != NULL){
        freelist(pp->next);
        free(pp);
    }
}
