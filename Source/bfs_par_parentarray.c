#include "project.h"
#include "mpi.h"

/*
 author: Alexander Gallauner
*/

int main(int argc, char *argv[]){
    int my_rank, procs, tag=0;
    uint64_t nodes = pow(2,SCALE);
    uint64_t edges = nodes*EDGEFACTOR;
    uint64_t root = ROOT;
        
    MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &procs); //SHOULD BE POWER OF TWO

    uint64_t *startVertex = NULL;
    uint64_t *endVertex = NULL;
    /* MUST BE INT BECAUSE OF MPI RESTRICTION */
    int *edgelist_send_counts = NULL;
    int *edgelist_send_displs = NULL;
    uint64_t *level = NULL;
    uint64_t *startVertex_recvbuf = NULL;
    uint64_t *endVertex_recvbuf = NULL;
    uint64_t *index_of_node = NULL;
    uint64_t *level_recvbuf = (uint64_t *) calloc(nodes / procs, sizeof(uint64_t));
    int edgelist_counts_recvbuf = 0;
    if (my_rank == 0){
        level = (uint64_t *) calloc(nodes, sizeof(uint64_t));
        startVertex = (uint64_t *) calloc(edges, I64_BYTES);
        endVertex = (uint64_t *) calloc(edges, I64_BYTES);
        edgelist_send_counts = (int *) calloc(procs, sizeof(int));
        edgelist_send_displs = (int *) calloc(procs, sizeof(int));

        read_graph(SCALE, EDGEFACTOR, startVertex, endVertex);
        
        double time = mytime();
        //SORTING THE EDGE LIST
        sort(startVertex, endVertex, 0, edges-1);
        
        //FINDING OUT THE BOUNDS OF THE EDGE LIST FOR EACH PROC
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
        
        MPI_Scatter((void *) edgelist_send_counts, 1, MPI_INT, &edgelist_counts_recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        startVertex_recvbuf = (uint64_t *) calloc(edgelist_counts_recvbuf, I64_BYTES);
        endVertex_recvbuf = (uint64_t *) calloc(edgelist_counts_recvbuf, I64_BYTES);
        
        MPI_Scatterv((void *) startVertex, edgelist_send_counts, edgelist_send_displs, MPI_UINT64_T, (void *) startVertex_recvbuf, edgelist_counts_recvbuf, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        
        MPI_Scatterv((void *) endVertex, edgelist_send_counts, edgelist_send_displs, MPI_UINT64_T, (void *) endVertex_recvbuf, edgelist_counts_recvbuf, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        
        index_of_node = create_buffer_from_edgelist(startVertex_recvbuf, endVertex_recvbuf, nodes / procs, edgelist_counts_recvbuf, my_rank);
	
        //SET ROOT LEVEL
        level[ROOT] = ROOT+1;

        //SCATTER LEVEL BUFFER
        MPI_Scatter((void *)level,nodes / procs,MPI_UINT64_T, (void *) level_recvbuf, nodes / procs, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        
        /*for (i = 0; i < index_of_node[(nodes/procs)]; i++){
            printf("%llu = %llu\n", (unsigned long long) buffer_recvbuf[i], (unsigned long long) startVertex_recvbuf[i]);
        }
        
        for (i = 0; i < nodes / procs; i++){
            printf("%llu = %llu\n", (unsigned long long) count_edges_per_node_recvbuf[i], (unsigned long long) index_of_node[i]);
        }*/
        
        //BFS
        time = mytime() - time;
        printf("Time for reading, generating edge buffer and scattering: %f\n", time/1000000);
        time = mytime();
        bfs(level_recvbuf, startVertex_recvbuf, index_of_node[nodes/procs], index_of_node, my_rank, procs, scale);

        time = mytime() - time;
        printf("Time for bfs searching: %f\n", time/1000000);

        free(edgelist_send_counts);
        free(edgelist_send_displs);
        free(startVertex);
        free(endVertex);
        free(level);
    }else{
        MPI_Scatter((void *) edgelist_send_counts, 1, MPI_INT, &edgelist_counts_recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        startVertex_recvbuf = (uint64_t *) calloc(edgelist_counts_recvbuf, I64_BYTES);
        endVertex_recvbuf = (uint64_t *) calloc(edgelist_counts_recvbuf, I64_BYTES);
        
        MPI_Scatterv((void *) startVertex, edgelist_send_counts, edgelist_send_displs, MPI_UINT64_T, (void *) startVertex_recvbuf, edgelist_counts_recvbuf, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        
        MPI_Scatterv((void *) endVertex, edgelist_send_counts, edgelist_send_displs, MPI_UINT64_T, (void *) endVertex_recvbuf, edgelist_counts_recvbuf, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        
        index_of_node = create_buffer_from_edgelist(startVertex_recvbuf, endVertex_recvbuf, nodes / procs, edgelist_counts_recvbuf, my_rank);
        
        // GET THE FIRST LEVEL
        MPI_Scatter(NULL, nodes / procs, MPI_UINT64_T, (void *)level_recvbuf, nodes / procs, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        
        bfs(level_recvbuf, startVertex_recvbuf, index_of_node[nodes/procs], index_of_node, my_rank, procs, scale);
    }
    free(level_recvbuf);
    free(startVertex_recvbuf);
    free(endVertex_recvbuf);
    free(index_of_node);
   
    MPI_Finalize ();
    return 0;
}

void bfs(uint64_t *level, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, int my_rank, int procs, int scale){
    uint64_t nodes_owned = pow(2,SCALE) / procs;
    uint64_t *parent_array = (unsigned long *) calloc(nodes_owned, sizeof(uint64_t));
    uint64_t *next_level = (unsigned long *) calloc(pow(2, SCALE), sizeof(uint64_t));
    double time_sum = 0;
    double time;
    int *recvcounts = (int *) calloc(procs, sizeof(int));
    char oneChildisVisited = 1;
    uint64_t position_level = nodes_owned*my_rank;
    uint64_t i;
    for (i = 0; i < procs; i++){
        recvcounts[i] = nodes_owned;
    }
    while (oneChildisVisited){
        oneChildisVisited = 0;
        for (i = 0; i < nodes_owned; i++){
            if (level[i] != 0 && parent_array[i] == 0){ //checks if the current node in the iteration is in the current level and not visited
                parent_array[i] = level[i];
                uint64_t j = index_of_node[i];
                if (i == nodes_owned -1){ //differentiate between the last node of the owned nodes and one node in the middle
                    for (; j < buffer_size; j++){
                        next_level[buffer[j]] = position_level+i+1;
                        oneChildisVisited = 1;
                    }
                }else{
                    for (; j < buffer_size && j < index_of_node[i+1]; j++){
                        next_level[buffer[j]] = position_level+i+1;
                        oneChildisVisited = 1;
                    }
                }
            }
        }

        // SEND MESSAGE THAT THERE ARE CHILDS TO EVALUATE, CAN BE ONE BYTE FROM ALL procs
        MPI_Allreduce(MPI_IN_PLACE, (void *) &oneChildisVisited, 1, MPI_CHAR, MPI_BOR, MPI_COMM_WORLD);
        // AFTER SEND LEVEL BUFFER, ALLTOALL
        if (oneChildisVisited){
            memset((void *) level, 0, nodes_owned*sizeof(uint64_t));
            time = mytime();
            MPI_Reduce_scatter((void *) next_level, (void *)level, recvcounts, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
            time_sum += ((mytime() - time)/1000000);
            memset((void *) next_level, 0, pow(2,SCALE)*sizeof(uint64_t));
        }
    }
    /*printf("parent array: ");
    for (i = 0; i < nodes_owned; i++){
        printf("%llu ", (unsigned long long) parent_array[i]);
    }
    printf("\n");*/
    printf("Time for reduce_scattering: %f\n", time_sum);
    free(parent_array);
    free(next_level);
    free(recvcounts);
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
