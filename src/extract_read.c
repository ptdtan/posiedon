//
// Created by ptdtan on 30/11/2016.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../minilib/bam.h"
#define LOOP(i, s, e) for(i = s; i < e; i++)

typedef struct {
    unsigned int node[2];		// Exon order of two terminal
    unsigned short count;		// Number of read split here
} Connect;

typedef struct {
    unsigned int len;			// gene length
    unsigned short *cover;		// Read cover for each exon
    unsigned short n_connect;	// Number of exon connection
    unsigned short n_extra;	// Number of exon connection
    Connect *connections;		// Store connections between exons
    Connect *extra;
} Read_count;

Read_count *init_readCount(unsigned int gene_len)
{
  Read_count *read_count = malloc(sizeof(Read_count));
  read_count->len = gene_len;
  read_count->cover = calloc(gene_len, sizeof(short));
  memset(read_count->cover, 0, gene_len*sizeof(short));
  read_count->n_connect = 0;
  read_count->connections = calloc(1, sizeof(Connect));
  read_count->n_extra = 0;
  read_count->extra = calloc(1, sizeof(Connect));
  return read_count;
}

void destroy_readCount(Read_count *r)
{
  free(r->cover);
  free(r->connections);
  free(r);
}

// Check which exon a pos belong to
short get_exon(unsigned int pos, unsigned int *e_start,
               unsigned int *e_end, unsigned short n_exon)
{
  unsigned short i;
  LOOP(i, 0, n_exon){
    if (pos >= e_start[i] && pos <= e_end[i]){
      return i;
    }
  }
  return -1;
}

// Add extra connects that not connect two exons
void add_extra(unsigned int p1, unsigned int p2, Read_count *data)
{

  /* Local var */
  unsigned short i;

  LOOP(i, 0, data->n_extra){
    if (p1 > data->extra[i].node[0] - 10 && p1 < data->extra[i].node[0] + 10 &&
        p2 > data->extra[i].node[1] - 10 && p2 < data->extra[i].node[1] + 10){
      ++data->extra[i].count;
      return;
    }
  }

  ++data->n_extra;
  data->extra = realloc(data->extra, data->n_extra*sizeof(Connect));

  data->extra[i].node[0] = p1;
  data->extra[i].node[1] = p2;
  data->extra[i].count = 1;
}

// Store connect between exons containing position p1/p2
void add_connect(unsigned int p1, unsigned int p2, Read_count *data,
                 unsigned int *e_start, unsigned int *e_end, unsigned short n_exon)
{

  /* Local var */
  short ex1, ex2; 	//Exon at two end of connection
  unsigned short i;

  ex1 = get_exon(p1, e_start, e_end, n_exon);
  ex2 = get_exon(p2, e_start, e_end, n_exon);
  if (ex1 < 0 || ex2 < 0){
    add_extra(p1, p2, data);
    return;
  }

  LOOP(i, 0, data->n_connect){
    if (ex1 == data->connections[i].node[0] &&
        ex2 == data->connections[i].node[1]){
      ++data->connections[i].count;
      return;
    }
  }

  ++data->n_connect;
  data->connections = realloc(data->connections,
                              data->n_connect*sizeof(Connect));

  data->connections[i].node[0] = ex1;
  data->connections[i].node[1] = ex2;
  data->connections[i].count = 1;
}

//Count read cover per base
void cover_count(bam1_t *Read, Read_count *data, unsigned int g_start,
                 unsigned int g_end, unsigned int *e_start,
                 unsigned int *e_end, unsigned short n_exon)
{
  /* Local var */
  unsigned int *cigar;		//Cigar string
  unsigned short i, r_i, k;	//iterator
  unsigned int c, l, p;		// type and len of each cigar block


  // Count cover for each base in read match/missmatch with gene's base
  // If read is split, make a connection between 2 site
  cigar = bam1_cigar(Read);
  r_i = 0;
  LOOP(i, 0, Read->core.n_cigar){
    c = bam_cigar_opchr(cigar[i]);
    l = cigar[i]>>BAM_CIGAR_SHIFT;
    if (c == 'M' || c == 'X'){
      LOOP(k, r_i, r_i + l){
        p = Read->core.pos + k;
        if (p < g_start || p > g_end) continue;
        ++data->cover[p - g_start];
      }
    } else if (c == 'N'){
      add_connect(Read->core.pos + r_i, Read->core.pos + r_i + l,
                  data, e_start, e_end, n_exon);
    }

    if (c != 'I' && c != 'S') r_i += l;
  }
}

Read_count *get_split_count(bamFile in_bam, bam_index_t *idx,
                            unsigned short n_exon, unsigned int *e_start, unsigned int *e_end,
                            unsigned int start, unsigned int end, unsigned short chr)
{
  /* Local var */
  Read_collect *read_inf;		//Collection of read in query region
  Read_count *ret_data;		//Return information
  unsigned int flag;			//Read flag
  unsigned int i;

  //Extract read from file-chromosome order-start-end
  //Return reads collection {n: number of read, r: store read_information}
  read_inf = init_reads();
  extract_read(in_bam, idx, chr, start, end, read_inf);

  // Count cover
  ret_data = init_readCount(end-start+1);
  LOOP(i, 0, read_inf->n){
    flag = read_inf->r[i].core.flag;
    if ((flag&BAM_FUNMAP) != 0 || (flag&BAM_FSECONDARY) != 0 ||
        read_inf->r[i].core.n_cigar == 0) continue;
    cover_count(&read_inf->r[i], ret_data, start,
                end, e_start, e_end, n_exon);
  }

  destroy_reads(read_inf);
  return ret_data;
}

int main(){
  char file[50] = "/home/ginny/SMO.bam";
  unsigned short n_exon = 12;
  unsigned int e_start[12] =
          {128828712,128843224,128845043,128845450,128845990,128846304,
           128848599,128849129,128850203,128850805,128851476,128851864};
  unsigned int e_end[12] =
          {128829323,128843430,128845253,128845623,128846210,128846428,
           128848692,128849238,128850389,128850954,128851611,128853385};
  unsigned int start = 128826713;
  unsigned int end = 128855383;
  unsigned short chr = 6; // ORDER of chromatid in chromosome (6: chr7)

  // Open file bam + bai
  bamFile in_bam = bam_open(file, "r");
  bam_index_t *idx = bam_index_load(file);

  Read_count *ret_data = get_split_count(in_bam, idx, n_exon,
                                         e_start, e_end, start, end, chr);

  // Destroy index & close file
  bam_index_destroy(idx);
  bam_close(in_bam);

  /* Temporary output file for drawing */
  FILE *out = fopen("sashimi", "w");
  unsigned int i = 0;
  fprintf(out, "var sashimi = {\nCover: [");
  LOOP(i, 0, ret_data->len){
    fprintf(out, "%u", ret_data->cover[i]);
    if (i < ret_data->len-1) fprintf(out, ", ");
  }
  fprintf(out, "],\nConnect:[\n");
  LOOP(i, 0, ret_data->n_connect){
    fprintf(out, "{n1: %u, n2: %u, cover: %u}", ret_data->connections[i].node[0], ret_data->connections[i].node[1], ret_data->connections[i].count);
    if (i < ret_data->n_connect - 1)fprintf(out, ", ");
  }

  fprintf(out, "],\nExtra:[\n");
  LOOP(i, 0, ret_data->n_extra){
    if (ret_data->extra[i].count > 1){
      fprintf(out, "{n1: %u, n2: %u, cover: %u}", ret_data->extra[i].node[0], ret_data->extra[i].node[1], ret_data->extra[i].count);
      if (i < ret_data->n_extra - 1)fprintf(out, ", ");
    }
  }

  fprintf(out, "],\ngene_start: %u,\ngene_end: %u,\nexon_start: [", start, end);
  LOOP(i, 0, n_exon){
    fprintf(out, "%u", e_start[i]);
    if (i < n_exon -1) fprintf(out, ", ");
  }
  fprintf(out, "],\nexon_end: [");
  LOOP(i, 0, n_exon){
    fprintf(out, "%u", e_end[i]);
    if (i < n_exon -1) fprintf(out, ", ");
  }
  fprintf(out, "]\n}");

  destroy_readCount(ret_data);
  return 0;
}

