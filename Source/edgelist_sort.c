#include "project.h"

int main(){
    uint64_t *startVertex = (uint64_t *) calloc(pow(2,SCALE)*EDGEFACTOR, I64_BYTES);
    uint64_t *endVertex = (uint64_t *) calloc(pow(2,SCALE)*EDGEFACTOR, I64_BYTES);
    read_graph(SCALE, EDGEFACTOR, startVertex, endVertex);
    
    double time = mytime();
    
    sort(startVertex, endVertex, 0, pow(2,SCALE)*EDGEFACTOR-1);
    
    time = mytime() - time;
    printf("Time for sorting the edge list: %f\n", time/1000000);
    
    write_graph(SCALE, EDGEFACTOR, startVertex, endVertex);
    return 0;
}

void sort(uint64_t *startVertex, uint64_t *endVertex, uint64_t l, uint64_t r){
    uint64_t j;
    
    if(l < r){
        j = partition(startVertex, endVertex, l, r);
        sort(startVertex, endVertex, l, j-1);
        sort(startVertex, endVertex, j+1, r);
    }
    
}

uint64_t partition(uint64_t *startVertex, uint64_t *endVertex, uint64_t l, uint64_t r){
    uint64_t pivot, i, j, t;
    pivot = startVertex[l];
    i = l; j = r+1;
    
    while(1){
        do ++i; while( startVertex[i] <= pivot && i <= r );
        do --j; while( startVertex[j] > pivot );
        if( i >= j ) break;
        t = startVertex[i]; startVertex[i] = startVertex[j]; startVertex[j] = t;
        t = endVertex[i]; endVertex[i] = endVertex[j]; endVertex[j] = t;
    }
    t = startVertex[l]; startVertex[l] = startVertex[j]; startVertex[j] = t;
    t = endVertex[l]; endVertex[l] = endVertex[j]; endVertex[j] = t;
    return j;
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