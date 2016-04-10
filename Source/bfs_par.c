#include "project.h"
#include "mpi.h"

double mytime(void){
    struct timeval now;
    gettimeofday(&now,NULL);
    return (double) ((long long)now.tv_usec+(long long)now.tv_sec*1000000);
}

int main(int argc, char *argv[]){
    int my_rank, procs, tag=0;
    uint64_t nodes = pow(2,SCALE);
    uint64_t edges = nodes*EDGEFACTOR;
    uint64_t root = ROOT;
        
    MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &procs);

    uint64_t *buffer = NULL;
    // WHY IS INT * AT SCATTERV NECESSARY
    int *buffer_send_counts = NULL; 
    int *buffer_displs = NULL;
    uint64_t buffer_recv_size = 0;
    unsigned char *level = NULL;    
    uint64_t *count_edges_per_node = NULL;
    uint64_t *buffer_recvbuf = NULL;
    unsigned char *level_recvbuf = (unsigned char *) calloc(nodes / BITS / procs, 1);
    uint64_t *count_edges_per_node_recvbuf = (uint64_t *) calloc(nodes / procs, I64_BYTES);
    if (my_rank == 0){
        double time = mytime();
        uint64_t *startVertex = (uint64_t *) calloc(edges, I64_BYTES);
        uint64_t *endVertex = (uint64_t *) calloc(edges, I64_BYTES);
        buffer_send_counts = (int *) calloc(procs, sizeof(int));
        buffer_displs = (int *) calloc(procs, sizeof(int));

        level = (unsigned char *) calloc(nodes / BITS, 1); //LEVEL BUFFER FOR MASTER
        
        float initiator[] = {0.25,0.25,0.25,0.25};
        
        struct edge **node_edge_list = (struct edge **) calloc(nodes, 8);
	    count_edges_per_node = (uint64_t *) calloc(nodes, I64_BYTES);

        generate_graph(SCALE, EDGEFACTOR, initiator, startVertex, endVertex);
        create_node_edge_lists(nodes, edges, startVertex, endVertex, node_edge_list, count_edges_per_node);
	    //MAYBE CREATING BUFFERS WITH REALLOC, SO THERE ARE NOT LISTS, create_node_edge_lists, NECESSARY

	    uint64_t buffer_size = 0;
	    buffer = lists_to_buffer(&buffer_size, node_edge_list, count_edges_per_node, 0, nodes-1);

        //SCATTER HOW MANY EDGES EACH NODES HAS
        int j;
        for (j = 0; j < procs; j++){
            buffer_send_counts[j] = getSumOfEdges(count_edges_per_node+j*(nodes / procs), nodes / procs);
            if (j){
                buffer_displs[j] = buffer_displs[j-1] + buffer_send_counts[j-1];            
            }else{
                buffer_displs[j] = 0;
            }
        }

        MPI_Scatter((void *)count_edges_per_node,nodes / procs, MPI_UINT64_T, (void *)count_edges_per_node_recvbuf,nodes / procs, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        
        buffer_recv_size = getSumOfEdges(count_edges_per_node_recvbuf, nodes / procs);
        buffer_recvbuf = (uint64_t *) calloc(buffer_recv_size, I64_BYTES);
        
        //SCATTER WHOLE BUFFER TO CHILD PROCESSES
        MPI_Scatterv((void *) buffer, buffer_send_counts, buffer_displs, MPI_UINT64_T, (void *) buffer_recvbuf, buffer_recv_size, MPI_UINT64_T, 0, MPI_COMM_WORLD);
	       	    
	    //TRANSFER COUNT_BUFFER TO INDEX_BUFFER FOR ROOT PROC
        uint64_t i;
	    uint64_t prev = count_edges_per_node_recvbuf[0];
	    uint64_t tmp;
	    count_edges_per_node_recvbuf[0] = 0;
	    for (i = 1; i < nodes / procs; i++){
		    tmp = count_edges_per_node_recvbuf[i];
		    count_edges_per_node_recvbuf[i] = prev + count_edges_per_node_recvbuf[i-1];
		    prev = tmp;
	    }
	
        //SET ROOT LEVEL
        level[(ROOT/BITS)] = level[(ROOT/BITS)] | (unsigned char) pow(2,(ROOT % BITS));

        //SCATTER LEVEL BUFFER
        
        MPI_Scatter((void *)level,nodes / BITS / procs,MPI_UNSIGNED_CHAR, (void *)level_recvbuf, nodes / BITS / procs, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        
        //BFS
        time = mytime() - time;
        printf("Time for generating and scattering: %f\n", time/1000000);
        time = mytime();
        bfs(level_recvbuf, buffer_recvbuf, buffer_recv_size, count_edges_per_node_recvbuf, (nodes / procs), procs);

	    //OUTPUT
        /*uint64_t j;
        struct edge *p;
        for (j = 0; j < (nodes / procs); j++){
            printf("%llu:%llu Elements:", (unsigned long long) j, (unsigned long long) count_edges_per_node_recvbuf[j]);
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
	    for (j = 0; j < buffer_recv_size; j++){
		    printf("%llu,", (unsigned long long) buffer_recvbuf[j]);
	    }
	    printf("\n\n");
        
	    for (j = 0; j < (nodes / procs); j++){
		    printf("%llu,", (unsigned long long) count_edges_per_node_recvbuf[j]);
	    }
	    printf("\n");
        */
        time = mytime() - time;
        printf("Time for bfs searching: %f\n", time/1000000);

        free(buffer_send_counts);
        free(buffer_displs);
        free(startVertex);
        free(endVertex);
	    free(buffer);
        free(level);
	    free(count_edges_per_node);
        freelists(nodes, node_edge_list);
    }else{
        // RECEIVE COUNT OF EDGES PER NODE
        MPI_Scatter((void *)count_edges_per_node,nodes / procs * I64_BYTES, MPI_BYTE, (void *)count_edges_per_node_recvbuf, nodes / procs * I64_BYTES, MPI_BYTE, 0, MPI_COMM_WORLD);

        // CALCULATE SUM OF EDGES THE PROC RECEIVES
        buffer_recv_size = getSumOfEdges(count_edges_per_node_recvbuf, nodes / procs);
        // ALLOCATE THE BUFFER TO RECEIVE THE DATA
        buffer_recvbuf = (uint64_t *) calloc(buffer_recv_size, I64_BYTES);
        // RECEIVE THE GRAPH DATA
        MPI_Scatterv((void *) buffer, buffer_send_counts, buffer_displs, MPI_UINT64_T, buffer_recvbuf, buffer_recv_size, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        //TRANSFER COUNT_BUFFER TO INDEX_BUFFER FOR SPEC PROC
        uint64_t i;
	    uint64_t prev = count_edges_per_node_recvbuf[0];
	    uint64_t tmp;
	    count_edges_per_node_recvbuf[0] = 0;
	    for (i = 1; i < nodes / procs; i++){
		    tmp = count_edges_per_node_recvbuf[i];
		    count_edges_per_node_recvbuf[i] = prev + count_edges_per_node_recvbuf[i-1];
		    prev = tmp;
	    }
        // GET THE FIRST LEVEL
        MPI_Scatter((void *)level,nodes / BITS / procs,MPI_UNSIGNED_CHAR, (void *)level_recvbuf, nodes / BITS / procs, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        bfs(level_recvbuf, buffer_recvbuf, buffer_recv_size, count_edges_per_node_recvbuf, (nodes / procs), procs);
    }
    free(level_recvbuf);
    free(count_edges_per_node_recvbuf);
    free(buffer_recvbuf);
   
    MPI_Finalize ();
    return 0;
}

void bfs(unsigned char *level, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, uint64_t nodes_owned, int procs){
    unsigned char *next_level = (unsigned char *) calloc(pow(2,SCALE) / BITS, 1);
    unsigned char *level_alltoall = (unsigned char *) calloc(pow(2,SCALE) / BITS, 1);
    unsigned char *visited = (unsigned char *) calloc(nodes_owned / BITS, 1);
    char oneChildisVisited = 1;
    uint64_t i;
    int rounds = 0;
    unsigned char position;
    while (oneChildisVisited){
        oneChildisVisited = 0;
        for (i = 0; i < nodes_owned; i++){
            position = (unsigned char) pow(2,(i % BITS));
            if (position & level[(i / BITS)] & ~visited[(i / BITS)]){
                visited[(i / BITS)] = visited[(i / BITS)] | position;
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
        /*printf("Next level:\n");
        for (i = 0; i < pow(2,SCALE); i++){
            if (((unsigned char) pow(2,(i % BITS))) & next_level[(i / BITS)]){
                printf("%llu,",(unsigned long long) i);
            }
        }
        printf("\n");*/

        // SEND MESSAGE THAT THERE ARE CHILDS TO EVALUATE, CAN BE ONE BYTE FROM ALL procs
        char isVisited_reduced = 0;
        MPI_Allreduce((void *) &oneChildisVisited, (void *) &isVisited_reduced, 1, MPI_CHAR, MPI_BOR, MPI_COMM_WORLD);
        oneChildisVisited = isVisited_reduced;
        // AFTER SEND LEVEL BUFFER, ALLTOALL
        if (oneChildisVisited){
            MPI_Alltoall((void *) next_level, pow(2,SCALE) / BITS / procs, MPI_UNSIGNED_CHAR, (void *) level_alltoall, pow(2,SCALE) / BITS / procs, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);

            memset(level, 0, nodes_owned / BITS);
            for (i = 0; i < (pow(2,SCALE) / BITS); i++){
                level[(i % (nodes_owned / BITS))] = level[(i % (nodes_owned / BITS))] | level_alltoall[i];
            }
       
            memset(next_level, 0, (pow(2,SCALE) / BITS));
        }
        rounds++;
    }
    printf("round: %i\n", rounds);
    free(level_alltoall);
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

uint64_t getSumOfEdges(uint64_t *count_edges_per_node, uint64_t nodes){
    uint64_t sum = 0;
    uint64_t i;
    for (i = 0; i < nodes; i++){
        sum += count_edges_per_node[i];
    }
    return sum;
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
