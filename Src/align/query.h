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
#include <stdio.h>
#include <zlib.h>
//#include "kstring.h"
//#include "bwt.h"
#include "kseq.h"
//#include "kvec.h"

KSEQ_INIT(gzFile, gzread)

#define cigar_t uint32_t
#define MAX_UNIQUE_HITS 1024
#define MAX_STR_LEN 1024
#define STRAND_FORWARD 0
#define STRAND_BACKWARD 1
/*  
typedef struct {
    uint32_t pos;
    uint8_t n_diff;
    uint8_t is_gap;//1:have gap 2:no gap
    uint16_t strand;//1 forwatd/0 backward
} hit_t;
*/
//typedef kvec_t(hit_t) hits_t;

//typedef kvec_t(uint32_t) vec_uint32_t;
typedef struct{
    //seq info
    char *name, *comment;
    int l_seq, n_ambiguous, l_trim3, l_trim5;  
    uint8_t *seq, *rseq, *qual, *apos;    

} qseq_t;   
typedef struct{
    //seq info
    /*  
    char *name, *comment;
    int l_seq;  
    uint8_t *seq, *rseq, *qual, *apos;    
    */
    qseq_t *qseq;
    char *name, *comment;
    int l_seq, n_ambiguous, l_trim3, l_trim5;  
    uint8_t *seq, *rseq, *qual, *apos, *target;    
    int8_t is_rev;
    int fr_flg[2];  
    uint8_t *read_seq; 
    int query_err;
    uint32_t (*target_idx_f)[2];
    uint32_t (*target_idx_r)[2];
    uint32_t (*target_idx)[2];
    int m_off;
    int l_off;
    int r_off;
    uint8_t q_flg[LEN_READ+10];    
    int8_t max_ext;
/*  
    int8_t max_ext_L[2];
    int8_t max_ext_R[2];
*/
    //best hit
    int flag;
    uint32_t pos;
    
    //alignment info
    int b0, b1;//best and secondary  hit score
    uint32_t ref_start, ref_end;//tlen
    uint32_t seq_start, seq_end;
    int n_cigar;
    uint32_t *cigar;
 
    int strand;
    uint8_t n_diff;
    uint8_t is_gap;
    int n_indel;
    int n_mismatch;
    uint8_t mapq;
    
    //multi hits
    //hits_t hits[2];
    int candi_thres;   
    int best_thres;   
    int error_thres;   
    int n_olp_for, n_olp_rev; 
    int rep_num; 
    //kstring_t *sam;
    int seed_num;
    int seed_add;

    int16_t *hash_idx[2];
    int16_t hash_pos[2][2*LEN_READ]; 
    int32_t hash_bak[2][2*LEN_READ]; 
 

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

int query_read_seq(queryio_t *qs, qseq_t *query);
int query_read_multiSeqs(queryio_t *qs, int n_seq, qseq_t *multiSeqs);


int query_read_multiPairedSeqs(queryio_t *qs[], int n_seq, qseq_t *multiSeqs);

//void query_set_hits(query_t *query, int max_hits, hits_t *hits0, hits_t* hits1);
uint32_t gen_mapq(uint32_t b0, uint32_t b1);

#endif
