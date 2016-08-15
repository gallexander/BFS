#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>

#define SCALE 25    // 2^30 MALLOC MAX = SCALE 26 + EDGEFACTOR 16
#define EDGEFACTOR 4

#define I64_BYTES 8
#define BLOCKS 4
#define BITS 64
#define GRAPHFILE "graph"

#define ROOT 0

struct edge {
    uint64_t end;
    struct edge *next;
};

double mytime(void){
    struct timeval now;
    gettimeofday(&now,NULL);
    return (double) ((long long)now.tv_usec+(long long)now.tv_sec*1000000);
}

//void bfs(unsigned long *level, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, uint64_t nodes_owned, int procs);
void bfs(uint64_t *level, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, int my_rank, int procs);
void generate_graph(int scale, int edgefactor, float *initiator, uint64_t *startVertex, uint64_t *endVertex);
uint64_t *create_buffer_from_edgelist(uint64_t *startVertex, uint64_t *endVertex, uint64_t nodes, uint64_t edges, uint64_t proc_number);
uint64_t calculate_size(uint64_t *count_edges_per_node, uint64_t first, uint64_t last);
uint64_t getSumOfEdges(uint64_t *count_edges_per_node, uint64_t nodes);

void write_graph(int scale, int edgefactor, uint64_t *startVertex, uint64_t *endVertex);
void read_graph(int scale, int edgefactor, uint64_t *startVertex, uint64_t *endVertex);

void sort(uint64_t *startVertex, uint64_t *endVertex, int64_t l, int64_t r);
int64_t partition(uint64_t *startVertex, uint64_t *endVertex, int64_t l, int64_t r);

//WORKING WITH LISTS
void create_node_edge_lists(uint64_t nodes, uint64_t edges, uint64_t *startVertex, uint64_t *endVertex, struct edge **node_edge_list, uint64_t *count_edges_per_node);
uint64_t *lists_to_buffer(uint64_t *size, struct edge **node_edge_list, uint64_t *count_edges_per_node, uint64_t first, uint64_t last);
void freelists(uint64_t nodes, struct edge **node_edge_list);
void freelist(struct edge *pp);
