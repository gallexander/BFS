#include "project.h"
#include "mpi.h"

/*
 author: Alexander Gallauner
*/

int main(int argc, char *argv[]){
    int my_rank, procs;
    uint64_t nodes = pow(2,SCALE);
    uint64_t edges = nodes*EDGEFACTOR;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &procs);

    uint64_t *startVertex = NULL;
    uint64_t *endVertex = NULL;
    double timer;
    struct result1 result;
    if (my_rank == 0){
        startVertex = (uint64_t *) calloc(edges, I64_BYTES);
        endVertex = (uint64_t *) calloc(edges, I64_BYTES);
        timer = mytime();
        read_graph(SCALE, EDGEFACTOR, startVertex, endVertex);
        timer = mytime() - timer;
        printf("Time for reading the graph: %f\n", timer/1000000);
        timer = mytime();
    }
    kernel_1(startVertex, endVertex, edges, procs, my_rank, &result);
    if (my_rank == 0){
        timer = mytime() - timer;
        printf("Time for generating edge buffer, sorting and scattering: %f\n", timer/1000000);
    }
    int i = 0;
    uint64_t *roots = NULL;
    if (my_rank == 0){
        roots = (uint64_t *) calloc(64, sizeof(uint64_t));
    }
    uint64_t root = 0;
    uint64_t root_local = 0;
    uint64_t answer = 0;
    srand(time(NULL));
    while (i < SEARCHKEY_CNT){
        if (my_rank == 0){
            root = rand() % ((uint64_t)pow(2,result.scale));
        }
        MPI_Bcast(&root, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        if (pow(2,result.scale)/procs*my_rank <= root && pow(2,result.scale)/procs*(my_rank+1) > root){
            root_local = root % (uint64_t)(pow(2,result.scale)/procs);
            if (result.index_of_node[root_local] < result.index_of_node[root_local+1]){
                //printf("i = %i, root = %llu, my_rank = %i\n", i, (unsigned long long) root, my_rank);
                answer = 1;
            }else{
                answer = 0;
            }
        }else{
            answer = 0;
        }
        MPI_Allreduce(&answer, &root_local, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
        if (root_local == 1){
            if (my_rank == 0){
                roots[i] = root;
            }
            i++;
        }else{
        }
    }
    if (my_rank == 0){
        timer = mytime();
    }
    uint64_t traversed_edges = kernel_2(result.buffer, result.index_of_node, my_rank, procs, result.scale, startVertex, endVertex, roots);
    if (my_rank == 0){
        timer = mytime() - timer;
        printf("Time for bfs searching: %f\n", timer/1000000);
        printf("Traversed edges: %llu\n", (unsigned long long) traversed_edges);
        printf("GTEPS: %f\n", traversed_edges/(timer/1000000.0)/1000000000);
        free(roots);
    }
    free(startVertex);
    free(endVertex);
    MPI_Finalize ();
    return 0;
}

void kernel_1(uint64_t *startVertex, uint64_t *endVertex, uint64_t edges, int procs, int my_rank, struct result1 *result){
    int scale = 0;
    int *edgelist_send_counts = NULL;
    int *edgelist_send_displs = NULL;
    int edgelist_counts_recvbuf = 0;
    uint64_t *startVertex_recvbuf = (uint64_t *) calloc(edges / procs, I64_BYTES);
    uint64_t *endVertex_recvbuf = (uint64_t *) calloc(edges / procs, I64_BYTES);
    
    MPI_Scatter((void *) startVertex, edges / procs, MPI_UINT64_T, startVertex_recvbuf, edges / procs, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Scatter((void *) endVertex, edges / procs, MPI_UINT64_T, endVertex_recvbuf, edges / procs, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    
    //SORTING THE EDGE LIST
    sort(startVertex_recvbuf, endVertex_recvbuf, 0, (edges / procs) -1);
    tree_merge(startVertex, endVertex, startVertex_recvbuf, endVertex_recvbuf, edges, my_rank, procs);
    
    if (my_rank == 0){
        scale = (int) (log2(startVertex[edges-1])+1);
    }
    
    MPI_Bcast(&scale, 1, MPI_INT, 0, MPI_COMM_WORLD);
    uint64_t nodes = pow(2,scale);
    
    //FINDING OUT THE BOUNDS OF THE EDGE LIST FOR EACH PROC
    if (my_rank == 0){
        edgelist_send_counts = (int *) calloc(procs, sizeof(int));
        edgelist_send_displs = (int *) calloc(procs, sizeof(int));
        int j;
        int last_node_number = 0;
        int core_count = 0;
        for (j = 0; j < procs; j++){
            last_node_number = nodes / procs * (j+1) - 1;
            core_count = (edges / procs * (j+1)) - 1;
            if (j < procs -1){
                while (startVertex[core_count] <= last_node_number) {
                    core_count++;
                }
                while (startVertex[core_count] > last_node_number){
                    core_count--;
                }
                if (j){
                    edgelist_send_counts[j] = core_count - edgelist_send_displs[j] + 1;
                }else{
                    edgelist_send_counts[j] = core_count + 1;
                }
                edgelist_send_displs[j+1] = core_count + 1;
            }else{
                edgelist_send_displs[0] = 0;
                edgelist_send_counts[j] = edges - edgelist_send_displs[j];
            }
        }
    }
    
    MPI_Scatter((void *) edgelist_send_counts, 1, MPI_INT, &edgelist_counts_recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    startVertex_recvbuf = (uint64_t *) calloc(edgelist_counts_recvbuf, I64_BYTES);
    endVertex_recvbuf = (uint64_t *) calloc(edgelist_counts_recvbuf, I64_BYTES);
    
    MPI_Scatterv((void *) startVertex, edgelist_send_counts, edgelist_send_displs, MPI_UINT64_T, (void *) startVertex_recvbuf, edgelist_counts_recvbuf, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Scatterv((void *) endVertex, edgelist_send_counts, edgelist_send_displs, MPI_UINT64_T, (void *) endVertex_recvbuf, edgelist_counts_recvbuf, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    
    result->index_of_node = create_buffer_from_edgelist(startVertex_recvbuf, endVertex_recvbuf, nodes / procs, edgelist_counts_recvbuf, my_rank);
    
    result->buffer = startVertex_recvbuf;
    result->scale = scale;
    
    free(edgelist_send_counts);
    free(edgelist_send_displs);
    free(endVertex_recvbuf);
}

uint64_t kernel_2(uint64_t *buffer, uint64_t *index_of_node, int my_rank, int procs, int scale, uint64_t *startVertex, uint64_t *endVertex, uint64_t *roots){
    uint64_t root;
    uint64_t nodes = pow(2, (scale));
    uint64_t *level;
    uint64_t *parent_array;
    uint64_t count = 0;
    uint64_t j;
    for (j = 0; j < SEARCHKEY_CNT; j++){
        level = (uint64_t *) calloc(nodes / BITS, sizeof(uint64_t));
        parent_array = (uint64_t *) calloc(pow(2,scale), sizeof(uint64_t));
        if (my_rank == 0){ 
            root = j;           
            //root = roots[j];
            level[(root/BITS)] = level[(root/BITS)] | (uint64_t) pow(2,(root % BITS));
        }
        MPI_Bcast((void *)level, nodes / BITS, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        bfs(level, buffer, index_of_node[nodes/procs], index_of_node, my_rank, procs, scale, parent_array);
        if (my_rank == 0){
            parent_array[root] = root + 1;
            /*uint64_t i;
            for (i = 0; i < pow(2,scale)*EDGEFACTOR; i++){
                if (parent_array[startVertex[i]] != 0){
                    count++;
                }
            }*/
        }
        free(parent_array);
        free(level);
    }
    free(buffer);
    free(index_of_node);
    return count;
}

void bfs(uint64_t *level, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, int my_rank, int procs, int scale, uint64_t *parent_array){
    uint64_t nodes_owned = pow(2,scale) / procs;
    uint64_t *visited = (uint64_t *) calloc(nodes_owned / BITS, sizeof(uint64_t));
    char oneChildisVisited = 1;
    int level_count = 0;
    uint64_t i;
    uint64_t position;
    while (oneChildisVisited){
        oneChildisVisited = 0;
        for (i = 0; i < nodes_owned; i++){
            position = (uint64_t) pow(2,(i % BITS));
            //checks if the current node in the iteration is in the current level and not visited
            if (position & level[((nodes_owned*my_rank+i) / BITS)] & ~visited[(i / BITS)]){
                visited[(i / BITS)] = visited[(i / BITS)] | position;
                uint64_t j = index_of_node[i];
                //differentiate between the last node of the owned nodes and one node in the middle
                if (i == nodes_owned -1){
                    for (; j < buffer_size; j++){
                        position = (uint64_t) pow(2, (buffer[j] % BITS));
                        if (position & ~level[(buffer[j]/BITS)]){
                            parent_array[buffer[j]] = (uint64_t) nodes_owned*my_rank+i+1;
                            level[(buffer[j]/BITS)] = level[(buffer[j]/BITS)] | position;
                            oneChildisVisited = 1;
                        }
                    }
                }else{
                    for (; j < buffer_size && j < index_of_node[i+1]; j++){
                        position = (uint64_t) pow(2, (buffer[j] % BITS));
                        if (position & ~level[(buffer[j]/BITS)]){
                            parent_array[buffer[j]] = (uint64_t) nodes_owned*my_rank+i+1;
                            level[(buffer[j]/BITS)] = level[(buffer[j]/BITS)] | position;
                            oneChildisVisited = 1;
                        }
                    }
                }
            }
        }

        // SEND MESSAGE THAT THERE ARE CHILDS TO EVALUATE, CAN BE ONE BYTE FROM ALL procs
        MPI_Allreduce(MPI_IN_PLACE, (void *) &oneChildisVisited, 1, MPI_CHAR, MPI_BOR, MPI_COMM_WORLD);
        // AFTER SEND LEVEL BUFFER, ALLTOALL
        if (oneChildisVisited){
            MPI_Allreduce(MPI_IN_PLACE, (void *)level, (pow(2,scale) / BITS), MPI_UINT64_T, MPI_BOR, MPI_COMM_WORLD);
            level_count++;
        }
    }
    if (my_rank){
        MPI_Reduce((void *) parent_array, NULL, pow(2, scale), MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
    }else{
        MPI_Reduce(MPI_IN_PLACE, (void *) parent_array, pow(2, scale), MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
        //printf("rounds: %i\n", level_count);
    }
    
    free(visited);
}

uint64_t *create_buffer_from_edgelist(uint64_t *startVertex, uint64_t *endVertex, uint64_t nodes, uint64_t edges, uint64_t proc_number){
    uint64_t i,j,k, duplicates;
    uint64_t *index_of_node = (uint64_t *) calloc(nodes+1, I64_BYTES);
    j = 0;
    k = 0;
    duplicates = 0;
    index_of_node[k++] = j;
    for (i = nodes*proc_number; i < nodes*(proc_number+1); i++){
        while (startVertex[j] == i){
            if ((j > index_of_node[k-1] && endVertex[j] == startVertex[j-duplicates-1]) || (startVertex[j] == endVertex[j])){
                duplicates++;
            }else{
                startVertex[j-duplicates] = endVertex[j];
            }
            j++;
        }
        index_of_node[k++] = j-duplicates;
    }
    //REALLOC startVertex
    return index_of_node;
}

void read_graph(int scale, int edgefactor, uint64_t *startVertex, uint64_t *endVertex){
    uint64_t i;
    FILE *fp;
    fp = fopen(GRAPHFILE, "r");
    for (i = 0; i < pow(2,scale)*edgefactor; i++){
        fscanf(fp, "%llu %llu\n", (unsigned long long *)(startVertex+i), (unsigned long long *)(endVertex+i));
    }
    fclose(fp);
}

void sort(uint64_t *startVertex, uint64_t *endVertex, int64_t l, int64_t r){
    int64_t j;
    
    if(l < r){
        j = partition(startVertex, endVertex, l, r);
        sort(startVertex, endVertex, l, j-1);
        sort(startVertex, endVertex, j+1, r);
    }
}

void tree_merge(uint64_t *endresult_start, uint64_t *endresult_end, uint64_t *startVertex, uint64_t *endVertex, uint64_t edges, int my_rank, int procs){
    int step = 1;
    MPI_Status status;
    uint64_t s = edges / procs;
    uint64_t *recv_buf_start = NULL;
    uint64_t *recv_buf_end = NULL;
    uint64_t *result_start = startVertex;
    uint64_t *result_end = endVertex;
    uint64_t *chunk_start = NULL;
    uint64_t *chunk_end = NULL;
    while (step < procs){
        if (my_rank % (2 * step) == 0){
            if (my_rank + step < procs){
                recv_buf_start = (uint64_t *) calloc(s, sizeof(uint64_t));
                recv_buf_end = (uint64_t *) calloc(s, sizeof(uint64_t));
                MPI_Recv((void *)recv_buf_start, s, MPI_UINT64_T, my_rank + step, 0, MPI_COMM_WORLD, &status);
                MPI_Recv((void *)recv_buf_end, s, MPI_UINT64_T, my_rank + step, 0, MPI_COMM_WORLD, &status);
                chunk_start = result_start;
                chunk_end = result_end;
                result_start = (uint64_t *) calloc(2*s, sizeof(uint64_t));
                result_end = (uint64_t *) calloc(2*s, sizeof(uint64_t));
                merge(result_start, result_end, chunk_start, chunk_end, s, recv_buf_start, recv_buf_end, s);
            }
        }else{
            int near = my_rank - step;
            MPI_Send((void*) result_start, s, MPI_UINT64_T, near, 0, MPI_COMM_WORLD);
            MPI_Send((void*) result_end, s, MPI_UINT64_T, near, 0, MPI_COMM_WORLD);
            break;
        }
        step = step * 2;
        s = s * 2;
    }
    if (my_rank == 0){
        memcpy(endresult_start, result_start, edges*sizeof(uint64_t));
        memcpy(endresult_end, result_end, edges*sizeof(uint64_t));
    }
    free(result_start);
    free(result_end);
}

void merge(uint64_t *result_start, uint64_t *result_end, uint64_t *v1_start, uint64_t *v1_end, uint64_t n1, uint64_t *v2_start, uint64_t *v2_end, uint64_t n2){
    uint64_t i = 0, j = 0, k = 0;
    while (i < n1 && j < n2){
        if (v1_start[i] < v2_start[j] || (v1_start[i] == v2_start[j] && v1_end[i] < v2_end[j])){
            result_start[k] = v1_start[i];
            result_end[k] = v1_end[i];
            i++;
            k++;
        }else{
            result_start[k] = v2_start[j];
            result_end[k] = v2_end[j];
            j++;
            k++;
        }
    }
    if (i == n1){
        while(j < n2){
            result_start[k] = v2_start[j];
            result_end[k] = v2_end[j];
            j++;
            k++;
        }
    }else{
        while (i < n1){
            result_start[k] = v1_start[i];
            result_end[k] = v1_end[i];
            i++;
            k++;
        }
    }
    free(v1_start);
    free(v1_end);
    free(v2_start);
    free(v2_end);
}

int64_t partition(uint64_t *startVertex, uint64_t *endVertex, int64_t l, int64_t r){
    int64_t pivot_start, pivot_end, i, j, t;
    pivot_start = startVertex[l];
    pivot_end = endVertex[l];
    i = l; j = r+1;
    
    while(1){
        do ++i; while( (startVertex[i] < pivot_start || (startVertex[i] == pivot_start && endVertex[i] <= pivot_end)) && i <= r );
        do --j; while( (startVertex[j] > pivot_start || (startVertex[j] == pivot_start && endVertex[j] > pivot_end)) );
        if( i >= j ) break;
        t = startVertex[i]; startVertex[i] = startVertex[j]; startVertex[j] = t;
        t = endVertex[i]; endVertex[i] = endVertex[j]; endVertex[j] = t;
    }
    t = startVertex[l]; startVertex[l] = startVertex[j]; startVertex[j] = t;
    t = endVertex[l]; endVertex[l] = endVertex[j]; endVertex[j] = t;
    return j;
}
