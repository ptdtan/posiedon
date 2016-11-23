#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>

#include "bam.h"
#include "khash.h"
#include "ksort.h"
#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif

KSORT_INIT(off, pair64_t, pair64_lt)
KHASH_MAP_INIT_INT(i, bam_binlist_t)

struct _bam_index_t{
	int32_t n;
	uint64_t n_no_coor; // unmapped reads without coordinate
	khash_t(i) **index;
	bam_lidx_t *index2;
};

int bam_is_be = 0, bam_verbose = 2, bam_no_B = 0;
/**********************
 * Work with bam file *
 **********************/

void bam_iter_destroy(bam_iter_t iter)
{
	if (iter) { free(iter->off); free(iter); }
}

uint32_t bam_calend(const bam1_core_t *c, const uint32_t *cigar)
{
	int k, end = c->pos;
	for (k = 0; k < c->n_cigar; ++k) {
		int op  = bam_cigar_op(cigar[k]);
		int len = bam_cigar_oplen(cigar[k]);
		if (op == BAM_CBACK) { // move backward
			int l, u, v;
			if (k == c->n_cigar - 1) break; // skip trailing 'B'
			for (l = k - 1, u = v = 0; l >= 0; --l) {
				int op1  = bam_cigar_op(cigar[l]);
				int len1 = bam_cigar_oplen(cigar[l]);
				if (bam_cigar_type(op1)&1) { // consume query
					if (u + len1 >= len) { // stop
						if (bam_cigar_type(op1)&2) v += len - u;
						break;
					} else u += len1;
				}
				if (bam_cigar_type(op1)&2) v += len1;
			}
			end = l < 0? c->pos : end - v;
		} else if (bam_cigar_type(op)&2) end += bam_cigar_oplen(cigar[k]);
	}
	return end;
}

int bam_remove_B(bam1_t *b)
{
	int i, j, end_j, k, l, no_qual;
	uint32_t *cigar, *new_cigar;
	uint8_t *seq, *qual, *p;
	// test if removal is necessary
	if (b->core.flag & BAM_FUNMAP) return 0; // unmapped; do nothing
	cigar = bam1_cigar(b);
	for (k = 0; k < b->core.n_cigar; ++k)
		if (bam_cigar_op(cigar[k]) == BAM_CBACK) break;
	if (k == b->core.n_cigar) return 0; // no 'B'
	if (bam_cigar_op(cigar[0]) == BAM_CBACK) goto rmB_err; // cannot be removed
	// allocate memory for the new CIGAR
	if (b->data_len + (b->core.n_cigar + 1) * 4 > b->m_data) { // not enough memory
		b->m_data = b->data_len + b->core.n_cigar * 4;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
		cigar = bam1_cigar(b); // after realloc, cigar may be changed
	}
	new_cigar = (uint32_t*)(b->data + (b->m_data - b->core.n_cigar * 4)); // from the end of b->data
	// the core loop
	seq = bam1_seq(b); qual = bam1_qual(b);
	no_qual = (qual[0] == 0xff); // test whether base quality is available
	i = j = 0; end_j = -1;
	for (k = l = 0; k < b->core.n_cigar; ++k) {
		int op  = bam_cigar_op(cigar[k]);
		int len = bam_cigar_oplen(cigar[k]);
		if (op == BAM_CBACK) { // the backward operation
			int t, u;
			if (k == b->core.n_cigar - 1) break; // ignore 'B' at the end of CIGAR
			if (len > j) goto rmB_err; // an excessively long backward
			for (t = l - 1, u = 0; t >= 0; --t) { // look back
				int op1  = bam_cigar_op(new_cigar[t]);
				int len1 = bam_cigar_oplen(new_cigar[t]);
				if (bam_cigar_type(op1)&1) { // consume the query
					if (u + len1 >= len) { // stop
						new_cigar[t] -= (len - u) << BAM_CIGAR_SHIFT;
						break;
					} else u += len1;
				}
			}
			if (bam_cigar_oplen(new_cigar[t]) == 0) --t; // squeeze out the zero-length operation
			l = t + 1;
			end_j = j; j -= len;
		} else { // other CIGAR operations
			new_cigar[l++] = cigar[k];
			if (bam_cigar_type(op)&1) { // consume the query
				if (i != j) { // no need to copy if i == j
					int u, c, c0;
					for (u = 0; u < len; ++u) { // construct the consensus
						c = bam1_seqi(seq, i+u);
						if (j + u < end_j) { // in an overlap
							c0 = bam1_seqi(seq, j+u);
							if (c != c0) { // a mismatch; choose the better base
								if (qual[j+u] < qual[i+u]) { // the base in the 2nd segment is better
									bam1_seq_seti(seq, j+u, c);
									qual[j+u] = qual[i+u] - qual[j+u];
								} else qual[j+u] -= qual[i+u]; // the 1st is better; reduce base quality
							} else qual[j+u] = qual[j+u] > qual[i+u]? qual[j+u] : qual[i+u];
						} else { // not in an overlap; copy over
							bam1_seq_seti(seq, j+u, c);
							qual[j+u] = qual[i+u];
						}
					}
				}
				i += len, j += len;
			}
		}
	}
	if (no_qual) qual[0] = 0xff; // in very rare cases, this may be modified
	// merge adjacent operations if possible
	for (k = 1; k < l; ++k)
		if (bam_cigar_op(new_cigar[k]) == bam_cigar_op(new_cigar[k-1]))
			new_cigar[k] += new_cigar[k-1] >> BAM_CIGAR_SHIFT << BAM_CIGAR_SHIFT, new_cigar[k-1] &= 0xf;
	// kill zero length operations
	for (k = i = 0; k < l; ++k)
		if (new_cigar[k] >> BAM_CIGAR_SHIFT)
			new_cigar[i++] = new_cigar[k];
	l = i;
	// update b
	memcpy(cigar, new_cigar, l * 4); // set CIGAR
	p = b->data + b->core.l_qname + l * 4;
	memmove(p, seq, (j+1)>>1); p += (j+1)>>1; // set SEQ
	memmove(p, qual, j); p += j; // set QUAL
	memmove(p, bam1_aux(b), b->l_aux); p += b->l_aux; // set optional fields
	b->core.n_cigar = l, b->core.l_qseq = j; // update CIGAR length and query length
	b->data_len = p - b->data; // update record length
	return 0;

rmB_err:
	b->core.flag |= BAM_FUNMAP;
	return -1;
}

static void swap_endian_data(const bam1_core_t *c, int data_len, uint8_t *data)
{
	uint8_t *s;
	uint32_t i, *cigar = (uint32_t*)(data + c->l_qname);
	s = data + c->n_cigar*4 + c->l_qname + c->l_qseq + (c->l_qseq + 1)/2;
	for (i = 0; i < c->n_cigar; ++i) bam_swap_endian_4p(&cigar[i]);
	while (s < data + data_len) {
		uint8_t type;
		s += 2; // skip key
		type = toupper(*s); ++s; // skip type
		if (type == 'C' || type == 'A') ++s;
		else if (type == 'S') { bam_swap_endian_2p(s); s += 2; }
		else if (type == 'I' || type == 'F') { bam_swap_endian_4p(s); s += 4; }
		else if (type == 'D') { bam_swap_endian_8p(s); s += 8; }
		else if (type == 'Z' || type == 'H') { while (*s) ++s; ++s; }
		else if (type == 'B') {
			int32_t n, Bsize = bam_aux_type2size(*s);
			memcpy(&n, s + 1, 4);
			if (1 == Bsize) {
			} else if (2 == Bsize) {
				for (i = 0; i < n; i += 2)
					bam_swap_endian_2p(s + 5 + i);
			} else if (4 == Bsize) {
				for (i = 0; i < n; i += 4)
					bam_swap_endian_4p(s + 5 + i);
			}
			bam_swap_endian_4p(s+1); 
		}
	}
}

int bam_read1(bamFile fp, bam1_t *b)
{
	bam1_core_t *c = &b->core;
	int32_t block_len, ret, i;
	uint32_t x[8];

	assert(BAM_CORE_SIZE == 32);
	if ((ret = bam_read(fp, &block_len, 4)) != 4) {
		if (ret == 0) return -1; // normal end-of-file
		else return -2; // truncated
	}
	if (bam_read(fp, x, BAM_CORE_SIZE) != BAM_CORE_SIZE) return -3;
	if (bam_is_be) {
		bam_swap_endian_4p(&block_len);
		for (i = 0; i < 8; ++i) bam_swap_endian_4p(x + i);
	}
	c->tid = x[0]; c->pos = x[1];
	c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
	c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
	c->l_qseq = x[4];
	c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];
	b->data_len = block_len - BAM_CORE_SIZE;
	if (b->m_data < b->data_len) {
		b->m_data = b->data_len;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	if (bam_read(fp, b->data, b->data_len) != b->data_len) return -4;
	b->l_aux = b->data_len - c->n_cigar * 4 - c->l_qname - c->l_qseq - (c->l_qseq+1)/2;
	if (bam_is_be) swap_endian_data(c, b->data_len, b->data);
	if (bam_no_B) bam_remove_B(b);
	return 4 + block_len;
}

/************************
 * Work with index file *
 ************************/
static bam_index_t *bam_index_load_core(FILE *fp)
{
	int i;
	char magic[4];
	bam_index_t *idx;
	if (fp == 0) {
		fprintf(stderr, "[bam_index_load_core] fail to load index.\n");
		return 0;
	}
	fread(magic, 1, 4, fp);
	if (strncmp(magic, "BAI\1", 4)) {
		fprintf(stderr, "[bam_index_load] wrong magic number.\n");
		fclose(fp);
		return 0;
	}
	idx = (bam_index_t*)calloc(1, sizeof(bam_index_t));	
	fread(&idx->n, 4, 1, fp);
	if (bam_is_be) bam_swap_endian_4p(&idx->n);
	idx->index = (khash_t(i)**)calloc(idx->n, sizeof(void*));
	idx->index2 = (bam_lidx_t*)calloc(idx->n, sizeof(bam_lidx_t));
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index;
		bam_lidx_t *index2 = idx->index2 + i;
		uint32_t key, size;
		khint_t k;
		int j, ret;
		bam_binlist_t *p;
		index = idx->index[i] = kh_init(i);
		// load binning index
		fread(&size, 4, 1, fp);
		if (bam_is_be) bam_swap_endian_4p(&size);
		for (j = 0; j < (int)size; ++j) {
			fread(&key, 4, 1, fp);
			if (bam_is_be) bam_swap_endian_4p(&key);
			k = kh_put(i, index, key, &ret);
			p = &kh_value(index, k);
			fread(&p->n, 4, 1, fp);
			if (bam_is_be) bam_swap_endian_4p(&p->n);
			p->m = p->n;
			p->list = (pair64_t*)malloc(p->m * 16);
			fread(p->list, 16, p->n, fp);
			if (bam_is_be) {
				int x;
				for (x = 0; x < p->n; ++x) {
					bam_swap_endian_8p(&p->list[x].u);
					bam_swap_endian_8p(&p->list[x].v);
				}
			}
		}
		// load linear index
		fread(&index2->n, 4, 1, fp);
		if (bam_is_be) bam_swap_endian_4p(&index2->n);
		index2->m = index2->n;
		index2->offset = (uint64_t*)calloc(index2->m, 8);
		fread(index2->offset, index2->n, 8, fp);
		if (bam_is_be)
			for (j = 0; j < index2->n; ++j) bam_swap_endian_8p(&index2->offset[j]);
	}
	if (fread(&idx->n_no_coor, 8, 1, fp) == 0) idx->n_no_coor = 0;
	if (bam_is_be) bam_swap_endian_8p(&idx->n_no_coor);
	return idx;
}

bam_index_t *bam_index_load_local(const char *_fn)
{
	FILE *fp;
	char *fnidx, *fn;

	fn = strdup(_fn);
	fnidx = (char*)calloc(strlen(fn) + 5, 1);
	strcpy(fnidx, fn); strcat(fnidx, ".bai");
	fp = fopen(fnidx, "rb");

	free(fnidx); free(fn);
	if (fp) {
		bam_index_t *idx = bam_index_load_core(fp);
		fclose(fp);
		return idx;
	} else return 0;
}

bam_index_t *bam_index_load(const char *fn)
{
	bam_index_t *idx;
	idx = bam_index_load_local(fn);
	if (idx == 0) fprintf(stderr, "[bam_index_load] fail to load BAM index.\n");
	return idx;
}

static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[BAM_MAX_BIN])
{
	int i = 0, k;
	if (beg >= end) return 0;
	if (end >= 1u<<29) end = 1u<<29;
	--end;
	list[i++] = 0;
	for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
	for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
	for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
	for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
	for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
	return i;
}

// bam_fetch helper function retrieves 
bam_iter_t bam_iter_query(const bam_index_t *idx, int tid, int beg, int end)
{
	uint16_t *bins;
	int i, n_bins, n_off;
	pair64_t *off;
	khint_t k;
	khash_t(i) *index;
	uint64_t min_off;
	bam_iter_t iter = 0;

	if (beg < 0) beg = 0;
	if (end < beg) return 0;
	// initialize iter
	iter = calloc(1, sizeof(struct __bam_iter_t));
	iter->tid = tid, iter->beg = beg, iter->end = end; iter->i = -1;
	//
	bins = (uint16_t*)calloc(BAM_MAX_BIN, 2);
	n_bins = reg2bins(beg, end, bins);
	index = idx->index[tid];
	if (idx->index2[tid].n > 0) {
		min_off = (beg>>BAM_LIDX_SHIFT >= idx->index2[tid].n)? idx->index2[tid].offset[idx->index2[tid].n-1]
			: idx->index2[tid].offset[beg>>BAM_LIDX_SHIFT];
		if (min_off == 0) { // improvement for index files built by tabix prior to 0.1.4
			int n = beg>>BAM_LIDX_SHIFT;
			if (n > idx->index2[tid].n) n = idx->index2[tid].n;
			for (i = n - 1; i >= 0; --i)
				if (idx->index2[tid].offset[i] != 0) break;
			if (i >= 0) min_off = idx->index2[tid].offset[i];
		}
	} else min_off = 0; // tabix 0.1.2 may produce such index files
	for (i = n_off = 0; i < n_bins; ++i) {
		if ((k = kh_get(i, index, bins[i])) != kh_end(index))
			n_off += kh_value(index, k).n;
	}
	if (n_off == 0) {
		free(bins); return iter;
	}
	off = (pair64_t*)calloc(n_off, 16);
	for (i = n_off = 0; i < n_bins; ++i) {
		if ((k = kh_get(i, index, bins[i])) != kh_end(index)) {
			int j;
			bam_binlist_t *p = &kh_value(index, k);
			for (j = 0; j < p->n; ++j)
				if (p->list[j].v > min_off) off[n_off++] = p->list[j];
		}
	}
	free(bins);
	if (n_off == 0) {
		free(off); return iter;
	}
	{
		bam1_t *b = (bam1_t*)calloc(1, sizeof(bam1_t));
		int l;
		ks_introsort(off, n_off, off);
		// resolve completely contained adjacent blocks
		for (i = 1, l = 0; i < n_off; ++i)
			if (off[l].v < off[i].v)
				off[++l] = off[i];
		n_off = l + 1;
		// resolve overlaps between adjacent blocks; this may happen due to the merge in indexing
		for (i = 1; i < n_off; ++i)
			if (off[i-1].v >= off[i].u) off[i-1].v = off[i].u;
		{ // merge adjacent blocks
#if defined(BAM_TRUE_OFFSET) || defined(BAM_VIRTUAL_OFFSET16)
			for (i = 1, l = 0; i < n_off; ++i) {
#ifdef BAM_TRUE_OFFSET
				if (off[l].v + BAM_MIN_CHUNK_GAP > off[i].u) off[l].v = off[i].v;
#else
				if (off[l].v>>16 == off[i].u>>16) off[l].v = off[i].v;
#endif
				else off[++l] = off[i];
			}
			n_off = l + 1;
#endif
		}
		bam_destroy1(b);
	}
	iter->n_off = n_off; iter->off = off;
	return iter;
}

int bam_validate1(const bam_header_t *header, const bam1_t *b)
{
	char *s;

	if (b->core.tid < -1 || b->core.mtid < -1) return 0;
	if (header && (b->core.tid >= header->n_targets || b->core.mtid >= header->n_targets)) return 0;

	if (b->data_len < b->core.l_qname) return 0;
	s = memchr(bam1_qname(b), '\0', b->core.l_qname);
	if (s != &bam1_qname(b)[b->core.l_qname-1]) return 0;

	// FIXME: Other fields could also be checked, especially the auxiliary data

	return 1;
}

static inline int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b)
{
	uint32_t rbeg = b->core.pos;
	uint32_t rend = b->core.n_cigar? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + 1;
	return (rend > beg && rbeg < end);
}

int bam_iter_read(bamFile fp, bam_iter_t iter, bam1_t *b)
{
	int ret;
	if (iter && iter->finished) return -1;
	if (iter == 0 || iter->from_first) {
		ret = bam_read1(fp, b);
		if (ret < 0 && iter) iter->finished = 1;
		return ret;
	}
	if (iter->off == 0) return -1;
	for (;;) {
		if (iter->curr_off == 0 || iter->curr_off >= iter->off[iter->i].v) { // then jump to the next chunk
			if (iter->i == iter->n_off - 1) { ret = -1; break; } // no more chunks
			if (iter->i >= 0) assert(iter->curr_off == iter->off[iter->i].v); // otherwise bug
			if (iter->i < 0 || iter->off[iter->i].v != iter->off[iter->i+1].u) { // not adjacent chunks; then seek
				bam_seek(fp, iter->off[iter->i+1].u, SEEK_SET);
				iter->curr_off = bam_tell(fp);
			}
			++iter->i;
		}
		if ((ret = bam_read1(fp, b)) >= 0) {
			iter->curr_off = bam_tell(fp);
			if (b->core.tid != iter->tid || b->core.pos >= iter->end) { // no need to proceed
				ret = bam_validate1(NULL, b)? -1 : -5; // determine whether end of region or error
				break;
			}
			else if (is_overlap(iter->beg, iter->end, b)) return ret;
		} else break; // end of file or error
	}
	iter->finished = 1;
	return ret;
}

void bam_index_destroy(bam_index_t *idx)
{
	khint_t k;
	int i;
	if (idx == 0) return;
	for (i = 0; i < idx->n; ++i) {
		khash_t(i) *index = idx->index[i];
		bam_lidx_t *index2 = idx->index2 + i;
		for (k = kh_begin(index); k != kh_end(index); ++k) {
			if (kh_exist(index, k))
				free(kh_value(index, k).list);
		}
		kh_destroy(i, index);
		free(index2->offset);
	}
	free(idx->index); free(idx->index2);
	free(idx);
}

int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)
{
	int ret;
	bam_iter_t iter;
	bam1_t *b;
	b = bam_init1();
	iter = bam_iter_query(idx, tid, beg, end);
	while ((ret = bam_iter_read(fp, iter, b)) >= 0) func(b, data);
	bam_iter_destroy(iter);
	bam_destroy1(b);
	return (ret == -1)? 0 : ret;
}

/******************
 * Extra function *
 ******************/
Read_collect *init_reads(){
	Read_collect *r = malloc(sizeof(Read_collect));
	r->n = 0;
	r->r = calloc(1, sizeof(bam1_t));
	return r;
}

void destroy_reads(Read_collect *r)
{
	unsigned int i;
	for (i = 0; i < r->n; i++)
		free(r->r[i].data);
	free(r->r);
	free(r);
}

static int get_read(const bam1_t *b, Read_collect *data)
{
	unsigned int i;
	++data->n;
	i = data->n-1;
	data->r = realloc(data->r, data->n*sizeof(bam1_t));
	memcpy(data->r + i, b, sizeof(bam1_t));
	data->r[i].data = calloc(b->data_len, sizeof(uint8_t));
	memcpy(data->r[i].data, b->data, b->data_len*sizeof(uint8_t));
}

unsigned int extract_read(bamFile in_bam, const bam_index_t *idx, int chr,
										 int beg, int end, Read_collect *ret_data)
{
	if (idx == 0) return;
	bam_fetch(in_bam, idx, chr, beg, end, ret_data, get_read);
	return 0;
}
