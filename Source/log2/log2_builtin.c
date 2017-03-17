#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

double mytime(void){
    struct timeval now;
    gettimeofday(&now,NULL);
    return (double) ((long long)now.tv_usec+(long long)now.tv_sec*1000000);
}

const int tab64[64] = {
    63,  0, 58,  1, 59, 47, 53,  2,
    60, 39, 48, 27, 54, 33, 42,  3,
    61, 51, 37, 40, 49, 18, 28, 20,
    55, 30, 34, 11, 43, 14, 22,  4,
    62, 57, 46, 52, 38, 26, 32, 41,
    50, 36, 17, 19, 29, 10, 13, 21,
    56, 45, 25, 31, 35, 16,  9, 12,
    44, 24, 15,  8, 23,  7,  6,  5};

int log2_64(uint64_t value);

int main(){
    FILE *fp;
    fp = fopen("numbers","r");
    int i;
    uint64_t x;
    uint64_t position;
    uint64_t counter = 0;
    double time = mytime();
    for (i = 0; i < pow(2,22); i++){
        fscanf(fp, "%llu\n", (unsigned long long *) (&x));
        while (x){
            position = LOG2(x);
            x = x & ~((unsigned long long) pow(2,position));
            counter++;
        }
    }
    time = mytime() -time;
    printf("Time for log2: %f counter: %llu\n", time/1000000, (unsigned long long) counter);
    fclose(fp);
    return 0;
}

int log2_64 (uint64_t value)
{
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    value |= value >> 32;
    return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}
