#include "project.h"
#include "mpi.h"

int main(int argc, char *argv[]){
    int my_rank, procs, tag=0;
    //float initiator[] = {0.7,0.1,0.1,0.1};
    float initiator[] = {0.57,0.19,0.19,0.05};
    MPI_Status status;
    
    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &procs);
    uint64_t nodes = pow(2,SCALE);
    uint64_t edges = nodes*EDGEFACTOR;
    uint64_t *startVertex = NULL;
    uint64_t *endVertex = NULL;
    uint64_t *index_buffer = (uint64_t *) calloc(nodes, sizeof(uint64_t));
    if (my_rank == 0){
        //INDEX OF NODES ARE SHUFFLED, SO NO LOCALITY CAN BE EXPLOITED
        shuffle(index_buffer, nodes);

        startVertex = (uint64_t *) calloc(edges, sizeof(uint64_t));
        endVertex = (uint64_t *) calloc(edges, sizeof(uint64_t));
    }
    else{
        startVertex = (uint64_t *) calloc((edges / procs), sizeof(uint64_t));
        endVertex = (uint64_t *) calloc((edges / procs), sizeof(uint64_t));
    }
    MPI_Bcast((void *) index_buffer, nodes, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    double time = mytime();
    generate_graph(SCALE, EDGEFACTOR, initiator, startVertex, endVertex, procs, my_rank, index_buffer);
    
    if (my_rank == 0){
        MPI_Gather((void*) startVertex, (edges / procs), MPI_UINT64_T, (void *) startVertex, (edges / procs), MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Gather((void*) endVertex, (edges / procs), MPI_UINT64_T, (void *) endVertex, (edges / procs), MPI_UINT64_T, 0, MPI_COMM_WORLD);
        write_graph(SCALE, EDGEFACTOR, startVertex, endVertex);
        time = mytime() - time;
        printf("Time for generating and writing the graph to a file: %f\n", time/1000000);
        /*time = mytime();
        read_graph(SCALE, EDGEFACTOR, startVertex, endVertex);
        time = mytime() - time;
        printf("Time for reading graph: %f\n", time/1000000);*/
    }else{
        MPI_Gather((void *) startVertex, (edges / procs), MPI_UINT64_T, NULL, (edges / procs), MPI_UINT64_T, 0, MPI_COMM_WORLD);
        MPI_Gather((void *) endVertex, (edges / procs), MPI_UINT64_T, NULL, (edges / procs), MPI_UINT64_T, 0, MPI_COMM_WORLD);
    }
    free(startVertex);
    free(endVertex);
    free(index_buffer);
    
    MPI_Finalize ();
    return 0;
}

void shuffle(uint64_t *index_buffer, uint64_t nodes){
    srand(time(NULL));
    uint64_t i;
    for (i = 0; i < nodes; i++){
        index_buffer[i] = i;
    }
    
    uint64_t w;
    uint64_t t;
    for (i = nodes - 1; i > 0; --i) {
        // generate random index
        w = rand()%i;
        // swap items
        t = index_buffer[i];
        index_buffer[i] = index_buffer[w];
        index_buffer[w] = t;
    }
}

void write_graph(int scale, int edgefactor, uint64_t *startVertex, uint64_t *endVertex){
    uint64_t i;
    FILE *fp;
    fp = fopen(GRAPHFILE,"w");
    for (i = 0; i < pow(2,scale)*edgefactor; i++){
        fprintf(fp, "%llu %llu\n", (unsigned long long)(startVertex[i]), (unsigned long long)(endVertex[i]));
    }
    fclose(fp);
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

void generate_graph(int scale, int edgefactor, float *initiator, uint64_t *startVertex, uint64_t *endVertex, int procs, int my_rank, uint64_t *index_buffer){
    srand(time(NULL)*(my_rank+1));
    int i;
    int help_initiator[BLOCKS];
    int next;
    uint64_t row = 0;
    uint64_t col = 0;
    for (i = 0; i < BLOCKS; i++){
        help_initiator[i] = (float)(initiator[i]*100.0);
        if (i > 0){
            help_initiator[i] += help_initiator[i-1];
        }
    }
    uint64_t h;
    for (h = 0; h < (((pow(2,scale)*edgefactor)/procs)); h=h+2){
        col = 0;
        row = 0;
        for (i = scale-1; i >= 0; i--){
            next = rand() % 100;
            int j;
            for (j = 0; j < BLOCKS; j++){
                if (next < help_initiator[j]){
                    next = j;
                    j = BLOCKS;
                }
            }
            if (next == 0){
            }else if (next == 1){
                col += pow(2,i);
            }else if (next == 2){
                row += pow(2,i);
            }else if (next == 3){
                col +=pow(2,i);
                row +=pow(2,i);
            }
        }
        startVertex[h] = index_buffer[row];
        startVertex[h+1] = index_buffer[col];
        endVertex[h] = index_buffer[col];
        endVertex[h+1] = index_buffer[row];
    }
}
