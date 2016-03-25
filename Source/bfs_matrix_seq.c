#include "project.h"

int isNotZero(char *array);
void transpose_matrix(char matrix[][NODES]);
void bfs_2d(char matrix[][NODES], char *level, char *parents);
void matrix_vector_multiplication(char matrix[][NODES], char *vector, char *result);
void negate(char *array);
void vector_multiplication_elementwise(char *vector1, char *vector2);
void vector_addition_elementwise(char *vector1, char *vector2);

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
    float initiator[] = {0.57,0.19,0.19,0.05};

    generate_graph(SCALE, EDGEFACTOR, startVertex, endVertex,initiator);

    char matrix[NODES][NODES] = {{0,1,0,0,0,0,0},
                                 {0,0,1,0,0,0,0},
                                 {0,0,0,1,0,0,0},
                                 {0,0,0,0,1,0,0},
                                 {0,0,0,0,0,1,0},
                                 {0,0,0,0,0,0,1},
                                 {0,0,0,0,0,0,0}};
    char parents[] = {0,0,0,0,0,0,0};
    char level[] = {0,0,0,0,0,0,0};

    double time = mytime();

    transpose_matrix(matrix);

    bfs_2d(matrix, level, parents);

    time = mytime() - time;
    printf("Time: %f\n", time/1000000);

    return 0;
}

void bfs_2d(char matrix[][NODES], char *level, char *parents){
    level[0] = 1;
    parents[0] = 1;
    char result[] = {0,0,0,0,0,0,0};
    while (isNotZero(level)){
        matrix_vector_multiplication(matrix, level, result);
        negate(parents);
        vector_multiplication_elementwise(result, parents);
        negate(parents);
        vector_addition_elementwise(parents, result);
        memcpy((void*)level, (void*)result, NODES);
        int i;
        for (i = 0; i < NODES; i++){
            printf("%i,", level[i]);
        }
        printf("\n");
        for (i = 0; i < NODES; i++){
            printf("%i,", parents[i]);
        }
        printf("\n\n");
    }
}

void matrix_vector_multiplication(char matrix[][NODES], char *vector, char *result){
    int i,j;
    int sum = 0;
    for (i = 0; i < NODES; i++){
        for (j = 0; j < NODES; j++){
            sum += matrix[i][j] * vector[j];
        }
        result[i] = sum;
        sum = 0;
    }
}

void vector_multiplication_elementwise(char *vector1, char *vector2){
    int i;
    for (i = 0; i < NODES; i++){
        vector1[i] = vector1[i] * vector2[i];
    }
}

void vector_addition_elementwise(char *vector1, char *vector2){
    int i;
    for (i = 0; i < NODES; i++){
        vector1[i] = vector1[i] + vector2[i];
    }
}

int isNotZero(char *array){
    int r = 0;
    int i;
    for (i = 0; i < NODES; i++){
        if (array[i]==1){
            r = 1;
        }
    }
    return r;
}

void negate(char *array){
    int i;
    for (i = 0; i < NODES; i++){
        if (array[i]==0){
            array[i] = 1;
        }else{
            array[i] = 0;
        }
    }
}

void transpose_matrix(char matrix[][NODES]){
    char matrix_transposed[NODES][NODES];
    int i,j;
    for (i = 0; i < NODES; i++){
        for (j = 0; j < NODES; j++){
            matrix_transposed[i][j] = matrix[j][i];
        }
    }
    memcpy((void*)matrix,(void*)matrix_transposed,NODES*NODES);
}
