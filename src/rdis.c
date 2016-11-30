//
// Created by ptdtan on 30/11/2016.
//  Some functions and libraries were imported from work_with_bam project of Phuong Thao
//

#include <stdint.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include "../minilib/bam.h"

#ifndef LOOP
#define LOOP(i,s,e) for(int i = s; i < e; i++)
#endif

#ifndef SHIFT_LEFT
#define SHIFT_LEFT(cigar,l) \
  memmove(&cigar[1], &cigar[0], (l-1)*sizeof(unsigned int); \
  --l;
#endif

#ifndef SHIFT_NLEFT
#define SHIFT_NLEFT(cigar, l , k) \
  memmove(&cigar[k-1], &cigar[k], (l-k)*sizeof(unsigned int)); \
  --l;
#endif

#ifndef MERGE_SOFTCLIP_LEFT
#define MERGE_SOFTCLIP_LEFT(cigar, l, pos, e_dis) \
  unsigned int l_clip = cigar[0]>>BAM_CIGAR_SHIFT; \
  unsigned int l_match = cigar[1]>>BAM_CIGAR_SHIFT; \
  unsigned int m_match = l_clip + l_match; \
  cigar[1] = bam_cigar_gen(m_match, bam_cigar_op('M')); \
  SHIFT_LEFT(cigar, l); \
  pos -= l_clip; \
  e_dis += l_clip;
#endif

#ifndef MERGE_SOFTCLIP_RIGHT
#define MERGE_SOFTCLIP_RIGHT(cigar, l, e_dis) \
  unsigned int l_clip = cigar[l-1]>>BAM_CIGAR_SHIFT; \
  unsigned int l_match = cigar[l-2]>>BAM_CIGAR_SHIFT; \
  unsigned int m_match = l_clip + l_match; \
  cigar[l-2] = bam_cigar_gen(m_match, BAM_CMATCH); \
  e_dis += l_clip; \
  cigar[l-1] = 0; \
  --l;
#endif

#ifndef INDEL_AT_SPLICE
#define INDEL_AT_SPLICE(cigar, flag, l) \
  LOOP(int i, 0, l-1){ \
    if(((bam_cigar_opchr(cigar[i])=='I' | bam_cigar_opchr(cigar[i])=='D' ) && bam_cigar_opchr(cigar[i+1])=='N')\
        | (bam_cigar_opchr(cigar[i])=='N' && (bam_cigar_opchr(cigar[i+1])=='I' | bam_cigar_opchr(cigar[i+1])=='D'))) \
      { flag++; \
        break; }\
  }
#endif

#ifndef INDEL_RIGHT
#define INDEL_RIGHT(cigar, l, i, t, l_id) \
  l_id = 0; \
  int j; \
  LOOP(j, i, l-1){ \
    if(bam_cigar_opchr(cigar[j])=='M' && (bam_cigar_opchr(cigar[j+1])=='I' | bam_cigar_opchr(cigar[j+1])=='D')) \
      { t = bam_cigar_opchr(cigar[j+1]); \
        l_id = cigar[j]>>BAM_CIGAR_SHIFT;\
        break; }\
        i = j;\
  }
#endif

#ifndef MERGE_INDEL
#define MERGE_INDEL(cigar, l, e_dis) \
  unsigned int l_indel = 0; \
  char tmp_c = '\0'; \
  int i = 0;\
  while(true) { \
    INDEL_RIGHT(cigar, l, i, tmp_c, l_indel) \
    if(!l_indel) \
      break; \
    if(tmp_c=='D'){ \
      l_indel+=cigar[i]>>BAM_CIGAR_SHIFT; \
      cigar[i+1] = bam_cigar_gen(l_indel, BAM_CMATCH); \
      SHIFT_NLEFT(cigar, l, i+1); \
    } else if(tmp_c=='I'){ \
        cigar[i+1] = bam_cigar_gen(l_indel, BAM_CMATCH); \
        SHIFT_NLEFT(cigar, l, i+1); \
    } else {\
      fprintf(stderr, "Didn't reconigze BAM CIGAR character %s", &tmp_c); \
      exit(EXIT_FAILURE); \
    }\
}
#endif

#ifndef MERGE_MATCH
#define MERGE_MATCH(cigar, l, i) \
  unsigned l_m1 = cigar[i]>>BAM_CIGAR_SHIFT; \
  unsigned l_m2 = cigar[i+1]>>BAM_CIGAR_SHIFT; \
  cigar[i+1] = bam_cigar_gen(l_m1+l_m2, BAM_CMATCH)\
  SHIFT_NLEFT(cigar, l, i+1);
#endif

#ifndef CONSECUTIVE_MATCH
#define CONSECUTIVE_MATCH(cigar, l, i, flag) \
  flag = 0; \
  LOOP(i, 0, l-1){ \
    if(bam_cigar_opchr(cigar[i])=='M' && bam_cigar_opchr(cigar[i+1])=='M'){ \
      flag++; \
      break; \
    } \
}
#endif

#ifndef MERGE_MATCH_ITER
#define MERGE_MATCH_ITER(cigar, l) \
  int i=0, flag; \
  while(true){ \
    CONSECUTIVE_MATCH(cigar, l, i, flag); \
    if(flag){ \
      MERGE_MATCH(cigar, l , i)\
    } else break; \
}
#endif

#ifndef DEBUG_CIGAR
#define DEBUG_CIGAR(cigar, l) \
  int i; \
  unsigned int cigar_c; \
  LOOP(i, 0, l){ \
    cigar_c = bam_cigar_opchr(cigar[i]); \
    printf("%d%s", cigar[i]>>BAM_CIGAR_SHIFT, &cigar_c); \
    } \
  printf("\n");
#endif

#ifndef MAX_SPLIT_AROUND
#define MAX_SPLIT_AROUND 200
#endif

typedef struct {
  int32_t beg;
  int32_t end;
  int32_t tid;
} Reg;

typedef struct {
  int32_t *cover;
  int32_t *pos;
} Bucket;

Reg *init_reg(int32_t beg, int32_t end, int32_t tid){
  Reg *ireg = (Reg*)malloc(sizeof(Reg));
  ireg->beg = beg;
  ireg->end = end;
  ireg->tid = tid;
  return ireg;
}

Bucket *init_bucket(void){
  Bucket *bucket = (Bucket*)malloc(sizeof(Bucket));
  bucket->cover = calloc(MAX_SPLIT_AROUND, sizeof(int32_t));
  memset(bucket->cover, 0, MAX_SPLIT_AROUND*sizeof(int32_t));
  memset(bucket->pos, 0, MAX_SPLIT_AROUND*sizeof(int32_t));
  bucket->pos = calloc(MAX_SPLIT_AROUND, sizeof(int32_t));
  return bucket;
}

void destroy_bucket(Bucket *bucket){
  free(bucket->cover);
  free(bucket->pos);
  free(bucket);
}

void destroy_reg(Reg *preg){
  free(preg);
}

void fill_bucket(bam1_t *read, Reg *region, Bucket *bucket){
  unsigned int *cigar;
  unsigned int r_i = 0, i;
  unsigned int c, l, p;

  cigar = bam1_cigar(read);
  LOOP(i, 0, read->core.n_cigar){
    c = bam_cigar_opchr(cigar[i]);
    l = cigar[i]>>BAM_CIGAR_SHIFT;
  }
}

int main(int argc, char *argv[]){
  int beg = atoi(argv[1]);
  int end = atoi(argv[2]);
  int tid = 20, it, c_i;
  unsigned int *cigar;

  bamFile fbam = bam_open(argv[3], "r");
  bam_index_t *idx = bam_index_load(argv[3]);

  Read_collect *rcollection = init_reads();
  extract_read(fbam, idx, tid, beg, end, rcollection);
  printf("number of reads in %d and %d %d\n", beg, end, rcollection->n);
  LOOP(it, 0, rcollection->n) {
    cigar = bam1_cigar(&rcollection->r[it]);
    DEBUG_CIGAR(cigar, rcollection->r[it].core.n_cigar)
  }
  return 0;
}