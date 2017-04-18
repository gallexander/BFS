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
    uint64_t *roots;
    double timer;
    struct result1 result;
    if (my_rank == 0){
        init_KISS(my_rank+1);
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
        roots = (uint64_t *) calloc(SEARCHKEY_CNT, sizeof(uint64_t));
        FILE *fp;
        if (!VALID_CHECKING){
            int i = 0;
            uint64_t root;
            fp = fopen(SEARCHKEYFILE, "w");
            while (i < SEARCHKEY_CNT){
                root = JKISS32() % ((uint64_t)pow(2,result.scale));
                if (result.index_of_node[root] < result.index_of_node[root+1]){
                    roots[i++] = root;
                    fprintf(fp, "%llu\n", (unsigned long long) root);
                }
            }
        }else{
            uint64_t i;
            fp = fopen(SEARCHKEYFILE, "r");
            for (i = 0; i < SEARCHKEY_CNT; i++){
                fscanf(fp, "%llu\n", (unsigned long long *)(roots+i));
            }
        }
        fclose(fp);
        timer = mytime();
        uint64_t traversed_edges = kernel_2(result.buffer, result.index_of_node, my_rank, procs, result.scale, startVertex, endVertex, roots);
        timer = mytime() - timer;
        printf("Time for bfs searching: %f\n", timer/1000000);
        printf("Traversed edges: %llu\n", (unsigned long long) traversed_edges);
        printf("GTEPS: %f\n", traversed_edges/(timer/1000000.0)/1000000000);
        free(roots);
        free(startVertex);
        free(endVertex);
    }
    MPI_Finalize ();
    return 0;
}

void kernel_1(uint64_t *startVertex, uint64_t *endVertex, uint64_t edges, int procs, int my_rank, struct result1 *result){
    int scale = 0;
    uint64_t *startVertex_recvbuf = (uint64_t *) calloc(edges / procs, I64_BYTES);
    uint64_t *endVertex_recvbuf = (uint64_t *) calloc(edges / procs, I64_BYTES);
    
    MPI_Scatter((void *) startVertex, edges / procs, MPI_UINT64_T, startVertex_recvbuf, edges / procs, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    MPI_Scatter((void *) endVertex, edges / procs, MPI_UINT64_T, endVertex_recvbuf, edges / procs, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    
    //SORTING THE EDGE LIST
    sort(startVertex_recvbuf, endVertex_recvbuf, 0, (edges / procs) -1);
    tree_merge(startVertex, endVertex, startVertex_recvbuf, endVertex_recvbuf, edges, my_rank, procs);
    
    if (my_rank == 0){
        scale = (int) (log2(startVertex[edges-1])+1);
        uint64_t nodes = (uint64_t) pow(2,scale);
        startVertex_recvbuf = (uint64_t *) calloc(edges, I64_BYTES);
        endVertex_recvbuf = (uint64_t *) calloc(edges, I64_BYTES);
        memcpy((void*)startVertex_recvbuf, (void*)startVertex, sizeof(uint64_t)*edges);
        memcpy((void*)endVertex_recvbuf, (void*)endVertex, sizeof(uint64_t)*edges);

        result->index_of_node = create_buffer_from_edgelist(startVertex_recvbuf, endVertex_recvbuf, nodes, edges, my_rank);
    
        result->buffer = startVertex_recvbuf;
        result->scale = scale;
        free(endVertex_recvbuf);
    }
}

uint64_t kernel_2(uint64_t *buffer, uint64_t *index_of_node, int my_rank, int procs, int scale, uint64_t *startVertex, uint64_t *endVertex, uint64_t *roots){
    uint64_t root;
    uint64_t nodes = pow(2, (scale));
    uint64_t *parent_array;
    int *distance_array = NULL;
    uint64_t count = 0;
    uint64_t j;
    for (j = 0; j < SEARCHKEY_CNT; j++){
        parent_array = (uint64_t *) calloc(pow(2,scale), sizeof(uint64_t));
        if (VALID_CHECKING){
            distance_array = (int *) calloc(pow(2,scale), sizeof(int));        
        }       
        root = roots[j];
        bfs_seq(root, buffer, index_of_node[nodes], index_of_node, scale, parent_array, distance_array);
        uint64_t i;
        uint64_t parent = 0;
        char valid = 1;
        uint64_t positions = 0;
        if (SEARCHKEY_CNT == 1){
            FILE *fp = NULL;
            if (!VALID_CHECKING){
                fp = fopen(PARENTFILE, "w");
                for (i = 0; i < pow(2,scale); i++){
                    fprintf(fp, "%llu\n", (unsigned long long) parent_array[i]);            
                }
            }else{
                fp = fopen(PARENTFILE, "r");
                for (i = 0; i < pow(2,scale); i++){
                    fscanf(fp, "%llu\n", (unsigned long long *)(&parent));
                    if (parent != parent_array[i]){
                        if (distance_array[parent-1] != distance_array[parent_array[i]-1]){
                            valid = 0;
                            positions++;
                            printf("%llu--> %llu:%i  != %llu:%i \n", (unsigned long long) i, (unsigned long long) parent-1, distance_array[parent-1], (unsigned long long) parent_array[i]-1, distance_array[parent_array[i]-1]);                   
                        }
                    }            
                }
                if (valid){
                    printf("The parent array from the file is VALID.\n");            
                }else{
                    printf("THE parent array from the file is NOT VALID.\n");
                    printf("NOT VALID ON %llu positions\n", (unsigned long long) positions);
                }
            }
            fclose(fp);
        }   
        /*uint64_t i;
        for (i = 0; i < pow(2,scale)*EDGEFACTOR; i++){
            if (parent_array[startVertex[i]] != 0){
                count++;
            }
        }*/
        free(parent_array);
        if (VALID_CHECKING){
            free(distance_array);    
        }
    }
    free(buffer);
    free(index_of_node);
    return count;
}

void bfs_seq(uint64_t root, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, int scale, uint64_t *parent_array, int *distance_array){
    uint64_t nodes = pow(2,scale);
    uint64_t FS_count = 0;
    uint64_t NS_count = 0;
    int level = 1;
    uint64_t *FS = (uint64_t *) calloc(pow(2,scale), sizeof(uint64_t));
    uint64_t *NS = (uint64_t *) calloc(pow(2,scale), sizeof(uint64_t));
    FS[FS_count++] = root;
    parent_array[root] = root + 1;
    if (VALID_CHECKING){
        distance_array[root] = 0;        
    }
    while (FS_count){
        uint64_t i;
        for (i = 0; i < FS_count; i++){
            uint64_t j = index_of_node[FS[i]];
            if (FS[i] == nodes -1){
                for (; j < buffer_size; j++){
                    if (parent_array[buffer[j]] == 0){
                        NS[NS_count++] = buffer[j];
                        parent_array[buffer[j]] = FS[i] + 1;
                        if (VALID_CHECKING){
                            distance_array[buffer[j]] = level;
                        }
                    }
                }
            }else{
                for (; j < buffer_size && j < index_of_node[FS[i]+1]; j++){
                    if (parent_array[buffer[j]] == 0){
                        NS[NS_count++] = buffer[j];
                        parent_array[buffer[j]] = FS[i] + 1;
                        if (VALID_CHECKING){
                            distance_array[buffer[j]] = level;
                        }
                    }
                }
            }
        }
        memset((void *)FS, 0, sizeof(uint64_t)*nodes);
        memcpy((void *)FS, (void *)NS, sizeof(uint64_t)*NS_count);
        FS_count = NS_count;
        NS_count = 0;
        level++;
    }
    free(FS);
    free(NS);
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
