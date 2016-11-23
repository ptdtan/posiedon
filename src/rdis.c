#include "bam.h"

#define LOOP(i,s,e) for(int i = s, i < e, i++)

typedef struct {
    int32_t beg;
    int32_t end;
    int32_t tid;
} reg;

reg *init_reg(int32_t beg, int32_t end, int32_t tid){
    reg *ireg = (reg*)malloc(sizeof(reg));
    ireg->beg = beg;
    ireg->end = end;
    ireg->tid = tid;
    return ireg;
}

void destroy_reg(reg *preg){
    free(preg);
}

int main(int argc, char *argv[]){
    int32_t beg = 14983070;
    int32_t end = 14983071;
    int32_t tid = 20;
    
    bamFile fbam = bam_open(argv[1], "r");
    bam_index_t *idx = bam_index_load(argv[1]);
    
    Read_collect *rcollection = init_reads();
    extract_read(fbam, idx, tid, beg, end, rcollection);
    printf("number of reads %d", rcollection->n);
    return 0;
}