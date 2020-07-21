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
//#include "kstring.h"
//#include "aln.h"
#include "kseq.h"
#include "kvec.h"

KSEQ_INIT(gzFile, gzread)

#define cigar_t uint32_t
#define MAX_UNIQUE_HITS 1024
#define MAX_STR_LEN 1024
#define STRAND_FORWARD 0
#define STRAND_BACKWARD 1

typedef struct{
  int sc;
  uint32_t ref_bg, ref_ed;
  uint32_t tb, te;
  int qb, qe;
  int strand;
} aln_t;


typedef kvec_t(uint32_t) vec_uint32_t;
typedef struct{
    char *name, *comment;
    int l_seq;  
    uint8_t *seq, *rseq, *qual;    
    uint32_t seq_start, seq_end;
    uint32_t ref_start, ref_end;
    //candidate_t candidate;
    //candidate_t rcandidate;

    int b0, b1;//best and secondary  hit score
    //best hit
    int flag;
    uint32_t pos;
    kstring_t *cigar;
    int strand;
    uint8_t mapq;
    int n_ambiguous;
    kvec_t(aln_t) xa;
    //multi hits
    //hits_t hits[2];
    //vec_uint_t ordered_index;
    //vec_uint_t unique_index;
    kstring_t *sam;
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

//void query_set_hits(query_t *query, int max_hits, hits_t *hits0, hits_t* hits1);
uint32_t gen_mapq(uint32_t b0, uint32_t b1);

#endif
