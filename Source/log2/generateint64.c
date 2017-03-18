#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

int main(){
    srand(time(NULL));
    FILE *fp;
    fp = fopen("numbers","w");
    uint64_t i;
    uint64_t r;
    for (i = 0; i < pow(2,23); i++){
        r = rand();
        r <<= 32;
        r = r | rand();
        fprintf(fp, "%llu\n", (unsigned long long) r);  
    }
    fclose(fp);
    return 0;
}
