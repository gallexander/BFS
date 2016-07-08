#include "project.h"

double mytime(void){
    struct timeval now;
    gettimeofday(&now,NULL);
    return (double) ((long long)now.tv_usec+(long long)now.tv_sec*1000000);
}

int main(){
    uint64_t nodes = pow(2,SCALE);
    uint64_t edges = nodes*EDGEFACTOR;
    uint64_t *startVertex = malloc(edges*I64_BYTES);
    uint64_t *endVertex = malloc(edges*I64_BYTES);
    float initiator[] = {0.25,0.25,0.25,0.25};
    //float initiator[] = {0.57,0.19,0.19,0.05};
    double time = mytime();
    generate_graph(SCALE, EDGEFACTOR, initiator, startVertex, endVertex);
    write_graph(SCALE, EDGEFACTOR, startVertex, endVertex);
    time = mytime() - time;
    printf("Time for generating and writing the graph to a file: %f\n", time/1000000);
    free(startVertex);
    free(endVertex);
    return 0;
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

void generate_graph(int scale, int edgefactor, float *initiator, uint64_t *startVertex, uint64_t *endVertex){
    srand(time(NULL));
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
    for (h = 0; h < (pow(2,scale)*edgefactor); h++){
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
        startVertex[h] = row;
        endVertex[h] = col;
    }
}
