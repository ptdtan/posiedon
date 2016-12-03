//
// Created by ptdtan on 30/11/2016.
//  Some functions and libraries were imported from work_with_bam project of Phuong Thao
//

#include <stdint.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include "bam.h"

#ifndef MAX_EDIT_DISTANCE
#define MAX_EDIT_DISTANCE 2
#endif

#ifndef DEBUG
#define DEBUG(l,i) printf("l: %d, i: %d\n", l, i);
#endif

#ifndef LOOP
#define LOOP(i,s,e) for(i = s; i < e; i++)
#endif

#ifndef SHIFT_LEFT
#define SHIFT_LEFT(cigar,l) \
  memmove(&cigar[0], &cigar[1], (l-1)*sizeof(unsigned int)); \
  --l;
#endif

#ifndef SHIFT_NLEFT
#define SHIFT_NLEFT(cigar, l , k) \
  memmove(&cigar[k-1], &cigar[k], (l-k)*sizeof(unsigned int)); \
  --l;
#endif

#ifndef MERGE_SOFTCLIP_LEFT
#define MERGE_SOFTCLIP_LEFT(cigar, l, pos, l_clip, l_match, e_dis) \
  l_clip = cigar[0]>>BAM_CIGAR_SHIFT; \
  l_match = cigar[1]>>BAM_CIGAR_SHIFT; \
  cigar[1] = bam_cigar_gen(l_clip + l_match, BAM_CMATCH); \
  SHIFT_LEFT(cigar, l); \
  pos -= l_clip; \
  e_dis+= l_clip;
#endif

#ifndef MERGE_SOFTCLIP_RIGHT
#define MERGE_SOFTCLIP_RIGHT(cigar, l, l_clip, l_match, e_dis) \
  l_clip = cigar[l-1]>>BAM_CIGAR_SHIFT; \
  l_match = cigar[l-2]>>BAM_CIGAR_SHIFT; \
  cigar[l-2] = bam_cigar_gen(l_clip + l_match, BAM_CMATCH); \
  cigar[l-1] = 0; \
  --l; \
  e_dis+= l_clip;
#endif

#ifndef INDEL_AT_SPLICE
#define INDEL_AT_SPLICE(cigar, l, i, flag) \
  flag = 0; \
  LOOP(i, 0, l-1){ \
    if(((bam_cigar_opchr(cigar[i])=='I' | bam_cigar_opchr(cigar[i])=='D' ) && bam_cigar_opchr(cigar[i+1])=='N')\
        | (bam_cigar_opchr(cigar[i])=='N' && (bam_cigar_opchr(cigar[i+1])=='I' | bam_cigar_opchr(cigar[i+1])=='D'))) \
      { flag = EXIT_FAILURE; \
        break; }\
  }
#endif

#ifndef INDEL_RIGHT
#define INDEL_RIGHT(cigar, l, i, t, l_id) \
  l_id = 0; \
  LOOP(i, 0, l-1){ \
    if(bam_cigar_opchr(cigar[i])=='M' && (bam_cigar_opchr(cigar[i+1])=='I' | bam_cigar_opchr(cigar[i+1])=='D')) \
      { t = bam_cigar_opchr(cigar[i+1]); \
        l_id = cigar[i+1]>>BAM_CIGAR_SHIFT;\
        break; }\
  }
#endif

#ifndef MERGE_INDEL
#define MERGE_INDEL(cigar, l, i, l_indel, tmp_c) \
  l_indel = 0; \
  while(1) { \
    INDEL_RIGHT(cigar, l, i, tmp_c, l_indel) \
    if(!l_indel) \
      break; \
    if(tmp_c=='D'){ \
      l_indel+=cigar[i]>>BAM_CIGAR_SHIFT; \
      cigar[i+1] = bam_cigar_gen(l_indel, BAM_CMATCH); \
      SHIFT_NLEFT(cigar, l, i+1); \
    } else if(tmp_c=='I'){ \
        cigar[i+1] = cigar[i]; \
        SHIFT_NLEFT(cigar, l, i+1) \
    } else {\
      fprintf(stderr, "Didn't reconigze BAM CIGAR character %s", &tmp_c); \
      exit(EXIT_FAILURE); \
    }\
}
#endif

#ifndef MERGE_MATCH
#define MERGE_MATCH(cigar, l, i, l_m1, l_m2) \
  l_m1 = cigar[i]>>BAM_CIGAR_SHIFT; \
  l_m2 = cigar[i+1]>>BAM_CIGAR_SHIFT; \
  cigar[i+1] = bam_cigar_gen(l_m1+l_m2, BAM_CMATCH);\
  SHIFT_NLEFT(cigar, l, i+1)
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
#define MERGE_MATCH_ITER(cigar, l, i, l_m1, l_m2, flag) \
  while(1){ \
    CONSECUTIVE_MATCH(cigar, l, i, flag) \
    if(flag){\
      MERGE_MATCH(cigar, l , i, l_m1, l_m2) \
    } else break; \
}
#endif

#ifndef EDIT_DISTANCE
#define EDIT_DISTANCE(cigar, l, i, e_dis) \
  e_dis = 0; \
  LOOP(i, 0, l){ \
    if(bam_cigar_opchr(cigar[i])=='S' | bam_cigar_opchr(cigar[i])=='D' | bam_cigar_opchr(cigar[i])=='I') \
      e_dis+=cigar[i]>>BAM_CIGAR_SHIFT;\
  }
#endif

#ifndef DEBUG_CIGAR
#define DEBUG_CIGAR(cigar, l, i, cigar_c) \
  LOOP(i, 0, l){ \
    cigar_c = bam_cigar_opchr(cigar[i]); \
    printf("%d%s", cigar[i]>>BAM_CIGAR_SHIFT, &cigar_c); \
    } \
  printf("\n");
#endif

#ifndef FIND_JUNCTIONS
#define FIND_JUNCTIONS(cigar, l, i, js, jidx, pos) \
  LOOP(i,0,l){ \
    if(bam_cigar_opchr(cih))\
}
#endif

#ifndef CHECK_SKIPPING
#define CHECK_SKIPPING(cigar, l, i, flag) \
  flag = 0; \
  LOOP(i, 0, l){ \
    if(bam_cigar_opchr(cigar[i])!='N') {\
      flag=EXIT_FAILURE;\
      break;\
      }\
  }
#endif

#ifndef MAX_SPLIT_AROUND
#define MAX_SPLIT_AROUND 200
#endif

typedef struct
{
  int32_t beg;
  int32_t end;
  int32_t tid;
} Junc;

typedef struct
{
  int32_t *cover;
} Bucket;

Junc *init_reg(int32_t beg, int32_t end, int32_t tid)
{
  Junc *pjunc = (Junc*)malloc(sizeof(Junc));
  pjunc->beg = beg;
  pjunc->end = end;
  pjunc->tid = tid;
  return pjunc;
}

Bucket *init_bucket(void)
{
  Bucket *bucket = (Bucket*)malloc(sizeof(Bucket));
  bucket->cover = calloc(MAX_SPLIT_AROUND, sizeof(int32_t));
  memset(bucket->cover, 0, MAX_SPLIT_AROUND*sizeof(int32_t));
  return bucket;
}

void destroy_bucket(Bucket *bucket)
{
  free(bucket->cover);
  free(bucket);
}

void destroy_reg(Junc *pjunc)
{
  pjunc->beg = 0;
  pjunc->end = 0;
  pjunc->tid = 0;
  free(pjunc);
}

int fill_bucket(bam1_t *read, Junc **junctions, Bucket **bucket)
{
  unsigned int *cigar;
  unsigned int i, cigar_c, flag;
  unsigned int l, e_dis, lm1, lm2, l_indel;
  int32_t p;

  cigar = bam1_cigar(read);
  l = read->core.n_cigar;
  DEBUG_CIGAR(cigar, l, i, cigar_c)
  CHECK_SKIPPING(cigar, l ,i, flag)
  if(flag)
    return EXIT_FAILURE;
  INDEL_AT_SPLICE(cigar, l, i, flag)
  if(flag)
    return EXIT_FAILURE;
  EDIT_DISTANCE(cigar, l, i, e_dis)
  if(e_dis>MAX_EDIT_DISTANCE)
    return EXIT_FAILURE;

  p = read->core.pos;
  DEBUG_CIGAR(cigar, l, i, cigar_c)
  if(bam_cigar_opchr(cigar[l-1])=='S')
  {
    MERGE_SOFTCLIP_RIGHT(cigar, l, lm1, lm2, e_dis)
  }
  if(bam_cigar_opchr(cigar[0])=='S')
  {
    MERGE_SOFTCLIP_LEFT(cigar, l, p, lm1, lm2, e_dis)
  }
  printf("after soft-clip remove: ");
  DEBUG_CIGAR(cigar, l, i, cigar_c)
  MERGE_INDEL(cigar, l, i , l_indel, cigar_c)
  printf("after indel remove: ");
  DEBUG_CIGAR(cigar, l, i, cigar_c)
  return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
  int beg = atoi(argv[1]);
  int end = atoi(argv[2]);
  int tid = 20, it, j, i, flag, l_id;
	unsigned int len, e_dis = 0, cigar_c, pos, l_indel, l_m1, l_m2;
  uint32_t *cigar_tmp;
  Bucket **buckets;
  Junc **juncs;
  bamFile fbam = bam_open(argv[3], "r");
  bam_index_t *idx = bam_index_load(argv[3]);
  bam1_t *reads;

  Read_collect *rcollection = init_reads();
  extract_read(fbam, idx, tid, beg, end, rcollection);
  printf("number of reads in %d and %d %d\n", beg, end, rcollection->n);
  reads = rcollection->r;
  LOOP(i, 0, rcollection->n){
    //int a = fill_bucket(&rcollection->r[i], juncs, buckets);
    //if(a == EXIT_FAILURE){
      //printf("accepted read: \n");
      cigar_tmp = bam1_cigar(&reads[i]);
      printf("%d \n", reads[i].core.n_cigar);
      //DEBUG_CIGAR(cigar_tmp, rcollection->r[i].core.n_cigar, j, cigar_c)
    }

  /*unsigned int match = bam_cigar_gen(51, BAM_CMATCH); //800
  unsigned int softc =  bam_cigar_gen(52, BAM_CSOFT_CLIP); //804
  unsigned int del =  bam_cigar_gen(53, BAM_CDEL); //802
  unsigned int ins = bam_cigar_gen(54, BAM_CINS); // 801
  unsigned int skip = bam_cigar_gen(100, BAM_CREF_SKIP); // 803
  unsigned int *cigar;

  len = 10;
  cigar = (unsigned int *)calloc(len, sizeof(unsigned int));
  memcpy(cigar, (unsigned int[]){match, match, skip, match, ins, match, del, ins, match, softc}, len*sizeof(unsigned int));
  EDIT_DISTANCE(cigar, len, i , e_dis)
  DEBUG_CIGAR(cigar, len, i, cigar_c)
  MERGE_SOFTCLIP_RIGHT(cigar, len, l_m1, l_m2)
  INDEL_AT_SPLICE(cigar, flag, len, i)
  printf("indel at splice flag: %d\n", flag);
  INDEL_RIGHT(cigar, len, i, cigar_c, l_id)
  printf("indel %s at %d with length %d\n", &cigar_c, i, l_id);
  MERGE_INDEL(cigar, len, i , l_indel, cigar_c)
  DEBUG_CIGAR(cigar, len, i, cigar_c)
  MERGE_MATCH_ITER(cigar, len, i, flag, l_m1, l_m2)
  DEBUG_CIGAR(cigar, len, i, cigar_c)
  printf("total edit-distance: %d\n", e_dis);
  DEBUG(len, i)
  free(cigar);*/
  return 0;
}
