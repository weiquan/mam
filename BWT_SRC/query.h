#ifndef __QUERY_H
#define __QUERY_H
/*-----------------------------------------------------------------------------
 *seed in C part:                       
 *  flag = 0
 *  offset = 0
 *  pos = seed pos in Ref 
 *
 *seed in R part:
 *  flag = 1
 *  offset = first R pos in seed
 *  pos = Ri
 *-----------------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <stdio.h>
//#include "kstring.h"
#include "kseq.h"
#include "kvec.h"

KSEQ_INIT(gzFile, gzread)

#define cigar_t uint32_t
#define MAX_UNIQUE_HITS 1024
#define MAX_STR_LEN 1024
#define STRAND_FORWARD 0
#define STRAND_BACKWARD 1
typedef struct{
    uint32_t bwtidx;
    uint32_t vs_idx;//variant_seq idx
    const char* var_preseq;//pre seq
    uint32_t l_aln;//non-variant align length
    const char* extend_seed;
    kstring_t template_seq;
    uint32_t edit_distance;
    uint32_t k, l;//template seq bwt range
    uint32_t pos;
} template_seed_t;
static inline void template_seed_clean(template_seed_t *p)
{
    p->bwtidx = (uint32_t)-1;
    p->vs_idx = (uint32_t) -1;
    p->var_preseq = NULL;
    p->l_aln = (uint32_t) -1;
    p->extend_seed = NULL;
    p->template_seq.l = p->template_seq.m = 0;
    p->template_seq.s = NULL;
    p->pos = (uint32_t )-1;
    p->k = (uint32_t) -1;
    p->l = (uint32_t) -1;
    p->pos = (uint32_t) -1;
    p->edit_distance = (uint32_t) -1;
}
static  void log_template_seed(template_seed_t *ts)
{
    fprintf(stderr, "\nbwtindex = %u\n", ts->bwtidx);
    fprintf(stderr, "vs_idx = %u\n", ts->vs_idx);
    fprintf(stderr, "l_aln = %u\n", ts->l_aln);
    fprintf(stderr, "template_seq = %s\n", ts->template_seq.s);
    fprintf(stderr, "template seed bwt index= %u-%u\n", ts->k, ts->l);
    fprintf(stderr, "pos=%u\n", ts->pos);
}
typedef struct{
    int offset;//seed start pos on read
    int l;//seed length
    uint32_t pos;//reference position
    int n;
    template_seed_t *ts;
} seed_t;
typedef struct{
    uint32_t pos;
    char *cigar;
} aln_t;
typedef struct{
    //Seq info
    char *name, *comment;
    int l_seq;  
    uint8_t *seq, *rseq, *qual;    
    int n_N;
    //seed info
    int n_seed;
    seed_t *seed_list;
    //alignment info
    uint32_t seq_start, seq_end;
    uint32_t pos;//primary alignment
    int l_cigar;
    char *cigar;
    int n_alt;
    aln_t *alt;//alternative alignments
    //sam
    //kstring_t *sam;
} query_t;


typedef struct{
    gzFile fp;
    kseq_t *kseq;
} queryio_t;


//query_t *query_init();
void query_destroy(query_t *query);
queryio_t *query_open(const char *fn_fa);
void query_seq_reverse(int len, uint8_t *seq, int is_comp);
void query_close(queryio_t *qs);

void query_gen_cigar(uint32_t l_ref, const uint32_t *mixRef, query_t *query);

int query_read_seq(queryio_t *qs, query_t *query);
int query_read_multiSeqs(queryio_t *qs, int n_seq, query_t *multiSeqs);
int query_read_multiPairedSeqs(queryio_t *qs[], int n_seq, query_t *multiSeqs);

uint32_t gen_mapq(uint32_t b0, uint32_t b1);
void print_query(query_t *q);
#endif
