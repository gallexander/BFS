#include "project.h"

int main(){
    int procs = 4;
    uint64_t *startVertex = (uint64_t *) calloc(pow(2,SCALE)*EDGEFACTOR, I64_BYTES);
    uint64_t *endVertex = (uint64_t *) calloc(pow(2,SCALE)*EDGEFACTOR, I64_BYTES);
    int *edgelist_send_displs = (int *) calloc(procs, sizeof(int));
    int *edgelist_send_counts = (int *) calloc(procs, sizeof(int));
    uint64_t nodes = pow(2,SCALE);
    uint64_t edges = nodes*EDGEFACTOR;
    read_graph(SCALE, EDGEFACTOR, startVertex, endVertex);
    
    double time = mytime();
    
    sort(startVertex, endVertex, 0, edges-1);
    
    //FINDING OUT THE BOUNDS OF THE EDGE LIST FOR EACH PROC
    int j;
    uint64_t last_node_number = 0;
    uint64_t core_count = 0;
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
    
    time = mytime() - time;
    printf("Time for sorting the edge list: %f\n", time/1000000);
    
    write_graph(SCALE, EDGEFACTOR, startVertex, endVertex);
    free(startVertex);
    free(endVertex);
    free(edgelist_send_displs);
    free(edgelist_send_counts);
    return 0;
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