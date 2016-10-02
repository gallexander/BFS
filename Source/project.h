#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>

#define SCALE 25    // 2^30 MALLOC MAX = SCALE 26 + EDGEFACTOR 16
#define EDGEFACTOR 8

#define I64_BYTES 8
#define BLOCKS 4
#define BITS 64
#define GRAPHFILE "graph"

#define ROOT 0

struct edge {
    uint64_t end;
    struct edge *next;
};

struct result1 {
    uint64_t *buffer;
    uint64_t *index_of_node;
    int scale;
};

double mytime(void){
    struct timeval now;
    gettimeofday(&now,NULL);
    return (double) ((long long)now.tv_usec+(long long)now.tv_sec*1000000);
}

//KERNELS
void kernel_1(uint64_t *startVertex, uint64_t *endVertex, uint64_t edges, int procs, int my_rank, struct result1 *result);
void kernel_2(uint64_t *buffer, uint64_t *index_of_node, int my_rank, int procs, int scale);

//BFS
//void bfs(unsigned long *level, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, uint64_t nodes_owned, int procs);
void bfs(uint64_t *level, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, int my_rank, int procs, int scale);

//GENERATING
void generate_graph(int scale, int edgefactor, float *initiator, uint64_t *startVertex, uint64_t *endVertex, int procs, int my_rank, uint64_t *index_buffer);
void shuffle(uint64_t *index_buffer, uint64_t nodes);
uint64_t *create_buffer_from_edgelist(uint64_t *startVertex, uint64_t *endVertex, uint64_t nodes, uint64_t edges, uint64_t proc_number);
uint64_t calculate_size(uint64_t *count_edges_per_node, uint64_t first, uint64_t last);
uint64_t getSumOfEdges(uint64_t *count_edges_per_node, uint64_t nodes);

void write_graph(int scale, int edgefactor, uint64_t *startVertex, uint64_t *endVertex);
void read_graph(int scale, int edgefactor, uint64_t *startVertex, uint64_t *endVertex);

// SORTING
void tree_merge(uint64_t *endresult_start, uint64_t *endresult_end, uint64_t *startVertex, uint64_t *endVertex, uint64_t edges, int my_rank, int procs);
void merge(uint64_t *result_start, uint64_t *result_end, uint64_t *v1_start, uint64_t *v1_end, uint64_t n1, uint64_t *v2_start, uint64_t *v2_end, uint64_t n2);
void sort(uint64_t *startVertex, uint64_t *endVertex, int64_t l, int64_t r);
int64_t partition(uint64_t *startVertex, uint64_t *endVertex, int64_t l, int64_t r);

//WORKING WITH LISTS
void create_node_edge_lists(uint64_t nodes, uint64_t edges, uint64_t *startVertex, uint64_t *endVertex, struct edge **node_edge_list, uint64_t *count_edges_per_node);
uint64_t *lists_to_buffer(uint64_t *size, struct edge **node_edge_list, uint64_t *count_edges_per_node, uint64_t first, uint64_t last);
void freelists(uint64_t nodes, struct edge **node_edge_list);
void freelist(struct edge *pp);
