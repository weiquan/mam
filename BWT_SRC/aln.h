#ifndef __ALN_H
#define __ALN_H
/*
 * =====================================================================================
 *
 *       Filename:  alnseed.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/27/2012 11:20:04 AM
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#include <pthread.h>
#include "kstring.h"
#include "utils.h"
#include "query.h"
#include "index.h"


#define N_SEQS 0x40000
//#define N_SEQS 1
#define EXTEND_LV 0
#define EXTEND_SW 1


typedef struct{
    int se;//1:SE;0:PE
    unsigned int max_tlen, min_tlen;
    int n_mismatch;
    int n_diff;

    unsigned long seed;

    int n_threads;
    char *fn_index;
    char *fn_read1;
    char *fn_read2; 
    char *rg_id;
    int l_read;
    int l_seed;
    int l_overlap;
    int ref;
    int use_sw_extend;
    int print_xa_cigar;
    int print_nm_md;
    kstring_t *cmd;
    uint32_t max_walk; 
    uint32_t max_seed;
    uint32_t max_locate;
    uint32_t max_hits;
    int extend_algo;
    int mismatch_penalty;
    int gapop_penalty;
    int gapext_penalty;

} opt_t;
typedef struct{
    int n_threads;
    int max_diff;
    int max_hits;
    //locationg arguments
    uint32_t max_locate; //max locations per bwt range
    //seeding arguments 
    int seed_only_ref;
    int extend_algo;
    uint32_t max_seed; //max locations per bwt range
    int l_seed;
    int l_overlap;
    //chaining 
    int max_chain_gap;
    int max_extend;

    //for PE
    uint32_t max_tlen, min_tlen;
    // smith-watemam arguments 
    int gap_op;
    int gap_ex;
    uint8_t flag;
    int filters;
    int filterd;
    int thres_score;
    //sam format arguments;
    int print_xa_cigar;
    int print_nm_md;
    char *rg_id;

} aln_opt_t;


typedef struct{
    uint32_t sp;
    uint32_t ep;
    uint32_t offset;
} sai_t;


typedef kvec_t(sai_t) vec_sai_t;
#ifdef HAVE_THREAD
typedef struct {
    int tid;
    idx_t *index;
    int n_query;
    query_t *query;
    aln_opt_t *aln_opt;


} thread_aux_t;
#endif
static inline aln_opt_t* aln_opt_init(const opt_t *opt){
    aln_opt_t *aln_opt = calloc(1, sizeof(aln_opt_t));
    aln_opt->n_threads = opt->n_threads;
    aln_opt->l_seed = opt->l_seed;
    aln_opt->max_diff = opt->n_diff; 
    
    aln_opt->min_tlen = opt->min_tlen;
    aln_opt->max_tlen = opt->max_tlen; 
    aln_opt->seed_only_ref = opt->ref;
    aln_opt->l_overlap = opt->l_overlap;
    aln_opt->max_locate = opt->max_locate; 
    aln_opt->max_chain_gap = 16; 
    aln_opt->max_extend = 50; 
    aln_opt->max_seed = opt->max_seed; 
    aln_opt->max_hits = 5;
    aln_opt->extend_algo = opt->extend_algo;
    
    
    aln_opt->gap_op = 3;
    aln_opt->gap_ex = 1;
    aln_opt->flag = 2;
//    opt->filters = 150;
    aln_opt->filters = 0;
    aln_opt->filterd = 20;
//    opt->thres_score = 150;
    aln_opt->thres_score = 50;
    aln_opt->print_nm_md = opt->print_nm_md;
    aln_opt->print_xa_cigar = opt->print_xa_cigar;
    aln_opt->rg_id = opt->rg_id;

    return aln_opt;

}
static inline void aln_opt_destroy(aln_opt_t *aln_opt){
    free(aln_opt); 
}

int usage();





int alnse_core(const opt_t *opt);

int aln_main(int argc, char *argv[]);


#endif
