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
    uint64_t *roots = NULL;
    if (my_rank == 0){
        roots = (uint64_t *) calloc(SEARCHKEY_CNT, sizeof(uint64_t));
        FILE *fp = fopen(SEARCHKEYFILE, "r");
        uint64_t i;
        for (i = 0; i < SEARCHKEY_CNT; i++){
            fscanf(fp, "%llu\n", (unsigned long long *)(roots+i));
            //printf("%llu, ", (unsigned long long)(roots[i]));
        }
        printf("\n");
        fclose(fp);
    }
    /*uint64_t root = 0;
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
    }*/
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
    uint64_t *parent_array;
    uint64_t count = 0;
    uint64_t j;
    
    double time_allwork = 0;
    double time_alltoall = 0;
    double time_parentgather = 0;
    for (j = 0; j < SEARCHKEY_CNT; j++){
        if (my_rank == 0){
            parent_array = (uint64_t *) calloc(pow(2,scale), sizeof(uint64_t));
        }else{
            parent_array = (uint64_t *) calloc(pow(2,scale)/procs, sizeof(uint64_t));
        }
        if (my_rank == 0){                     
            root = roots[j];
        }
        MPI_Bcast((void *)&root, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        bfs_alltoall(root, buffer, index_of_node[nodes/procs], index_of_node, my_rank, procs, scale, parent_array, &time_allwork, &time_alltoall, &time_parentgather);
        if (my_rank == 0){
            if (SEARCHKEY_CNT == 1){
                FILE *fp;
                fp = fopen(PARENTFILE, "w");
                uint64_t i;
                for (i = 0; i < pow(2,scale); i++){
                    fprintf(fp, "%llu\n", (unsigned long long) parent_array[i]);            
                }
                fclose(fp);
            }
            /*uint64_t i;
            for (i = 0; i < pow(2,scale)*EDGEFACTOR; i++){
                if (parent_array[startVertex[i]] != 0){
                    count++;
                }
            }*/
        }
        free(parent_array);
    }
    if (my_rank == 0)
        printf("Time for all-to-all: %f\nTime for work without all-to-all: %f\nTime for gather parent array: %f\n", time_alltoall/1000000, (time_allwork-time_alltoall)/1000000, time_parentgather/1000000);
    free(buffer);
    free(index_of_node);
    return count;
}

void bfs_alltoall(uint64_t root, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, int my_rank, int procs, int scale, uint64_t *parent_array, double *time_allwork, double *time_alltoall, double *time_parentgather){
    uint64_t nodes_owned = pow(2,scale) / procs;
    uint64_t *visited = (uint64_t *) calloc(pow(2,scale) / BITS, sizeof(uint64_t));
    int *sdispls = (int *) calloc(procs, sizeof(int));
    int *rcounts = (int *) calloc(procs, sizeof(int));
    int *rdispls = (int *) calloc(procs, sizeof(int));
    uint64_t *prepare_edges[procs];
    int prepare_edges_count[procs];
    uint64_t *edges_received;
    uint64_t edges_count = 0;
    uint64_t position = 0;
    int proc_nr = -1;
    
    char oneChildisVisited = 1;
    int level_count = 1;
    uint64_t intern_node_number = 0;
    
    double time = 0;
    double time_work = 0;
    double time_parent = 0;
    char this_round_active = 0;
    if ((root / nodes_owned) == my_rank){
        edges_count = 2;
        edges_received = (uint64_t *) calloc(2, sizeof(uint64_t));
        edges_received[0] = root;
        edges_received[1] = root;
        this_round_active = 1;     
    }
    uint64_t i;
    for (i = 0; i < procs; i++){
        prepare_edges[i] = (uint64_t *) calloc(nodes_owned*2, sizeof(uint64_t));
    }
    if (my_rank == 0){
        time_work = mytime();
    }
    while (oneChildisVisited){
        oneChildisVisited = 0;
        //initialize the buffers that are sent to the procs at the end
        for (i = 0; i < procs; i++){
            memset((void *) prepare_edges[i], 0, sizeof(uint64_t)*(nodes_owned*2));
            prepare_edges_count[i] = 0;
        }
        memset((void *) sdispls, 0, sizeof(int)*procs);
        memset((void *) rcounts, 0, sizeof(int)*procs);
        memset((void *) rdispls, 0, sizeof(int)*procs);
        for (i = 0; i < edges_count; i=i+2){
            intern_node_number = edges_received[i] % nodes_owned;
            if (parent_array[intern_node_number] == 0){
                parent_array[intern_node_number] = edges_received[i+1]+1;
                uint64_t j = index_of_node[intern_node_number];
                for (; j < index_of_node[intern_node_number+1]; j++){
                    position = (uint64_t) pow(2, (buffer[j] % BITS));
                    if (position & ~visited[(buffer[j]/BITS)]){
                        visited[(buffer[j]/BITS)] = visited[(buffer[j]/BITS)] | position;
                        proc_nr = buffer[j] / nodes_owned;
                        prepare_edges[proc_nr][prepare_edges_count[proc_nr]] = buffer[j];
                        prepare_edges[proc_nr][prepare_edges_count[proc_nr]+1] = edges_received[i];
                        prepare_edges_count[proc_nr] += 2;
                        oneChildisVisited = 1;
                    }
                }
            }
        }
        // SEND MESSAGE THAT THERE ARE CHILDS TO EVALUATE, CAN BE ONE BYTE FROM ALL procs
        MPI_Allreduce(MPI_IN_PLACE, (void *) &oneChildisVisited, 1, MPI_CHAR, MPI_BOR, MPI_COMM_WORLD);
        // AFTER SEND LEVEL BUFFER, ALLTOALL
        if (my_rank == 0){
            time = mytime();
        }
        if (this_round_active){
            free(edges_received);
        }
        this_round_active = 1;
        if (oneChildisVisited){
            uint64_t *sendbuf;
            uint64_t sdispls_sum = 0;            
            edges_count = 0;
            MPI_Alltoall((void *)prepare_edges_count, 1, MPI_INT, (void *)rcounts, 1, MPI_INT, MPI_COMM_WORLD);
            for (i = 0; i < procs; i++){
                rdispls[i] = edges_count;
                sdispls[i] = sdispls_sum;
                edges_count += rcounts[i];
                sdispls_sum += prepare_edges_count[i];
            }
            sendbuf = (uint64_t *) calloc(sdispls_sum, sizeof(uint64_t));
            for (i = 0; i < procs; i++){
                memcpy((void *)(sendbuf+sdispls[i]), (void *) prepare_edges[i], prepare_edges_count[i]*sizeof(uint64_t)); 
            }
            edges_received = (uint64_t *) calloc(edges_count, sizeof(uint64_t));

            MPI_Alltoallv((void *)sendbuf, prepare_edges_count, sdispls, MPI_UINT64_T, edges_received, rcounts, rdispls, MPI_UINT64_T, MPI_COMM_WORLD);
            
            level_count++;
            free(sendbuf);
        }
        if (my_rank == 0){
            time = mytime() - time;
            *time_alltoall += time;
        }
    }
    if (my_rank == 0){
        time_work = mytime() - time_work;
        *time_allwork += time_work;
        time_parent = mytime();
        MPI_Gather(MPI_IN_PLACE, pow(2,scale)/procs, MPI_UINT64_T, (void *)parent_array, pow(2,scale)/procs, MPI_UINT64_T, 0, MPI_COMM_WORLD);
        time_parent = mytime() - time_parent;
        *time_parentgather += time_parent;
    }else{
        MPI_Gather((void *)parent_array, pow(2,scale)/procs, MPI_UINT64_T, NULL, pow(2,scale)/procs, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    }
    for (i = 0; i < procs; i++){
        free(prepare_edges[i]);
    }
    free(visited);
    free(sdispls);
    free(rcounts);
    free(rdispls);
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
