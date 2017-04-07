#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>

#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))
#define SCALE 18
#define EDGEFACTOR 32
#define SEARCHKEY_CNT 1
#define VALID_CHECKING 0
#define PARENTFILE "parent"
#define GRAPHFILE "graph"
#define SEARCHKEYFILE "keys"
#define BITS 64

#define I64_BYTES 8
#define BLOCKS 4


#define ROOT 0

const int tab64[64] = {
    63,  0, 58,  1, 59, 47, 53,  2,
    60, 39, 48, 27, 54, 33, 42,  3,
    61, 51, 37, 40, 49, 18, 28, 20,
    55, 30, 34, 11, 43, 14, 22,  4,
    62, 57, 46, 52, 38, 26, 32, 41,
    50, 36, 17, 19, 29, 10, 13, 21,
    56, 45, 25, 31, 35, 16,  9, 12,
    44, 24, 15,  8, 23,  7,  6,  5};

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

int log2_64(uint64_t);

//KERNELS
void kernel_1(uint64_t *startVertex, uint64_t *endVertex, uint64_t edges, int procs, int my_rank, struct result1 *result);
uint64_t kernel_2(uint64_t *buffer, uint64_t *index_of_node, int my_rank, int procs, int scale, uint64_t *startVertex, uint64_t *endVertex, uint64_t *roots);

//BFS
//void bfs(unsigned long *level, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, uint64_t nodes_owned, int procs);
void bfs(uint64_t *level, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, int my_rank, int procs, int scale, uint64_t *parent_array, double *time_allwork, double *time_allreduce, double *time_parentreduce);
void bfs_seq(uint64_t root, uint64_t *buffer, uint64_t buffer_size, uint64_t *index_of_node, int scale, uint64_t *parent_array, int *distance_array);

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
