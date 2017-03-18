#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

double mytime(void){
    struct timeval now;
    gettimeofday(&now,NULL);
    return (double) ((long long)now.tv_usec+(long long)now.tv_sec*1000000);
}

int main(){
    FILE *fp;
    fp = fopen("numbers","r");
    int i;
    uint64_t x;
    uint64_t position;
    uint64_t counter = 0;
    double time = mytime();
    for (i = 0; i < pow(2,23); i++){
        fscanf(fp, "%llu\n", (unsigned long long *) (&x));
        int j;
        for (j = 0; j < 64; j++){
            position = pow(2,j);
            if (position & x){
                counter++;
            }
        }
    }
    time = mytime() -time;
    printf("Time for log2: %f counter: %llu\n", time/1000000, (unsigned long long) counter);
    fclose(fp);
    return 0;
}
