/*
 * =====================================================================================
 *
 *       Filename:  aln.c
 *
 *    Description:  O
 *
 *        Version:  1.0
 *        Created:  2017年10月12日 10时54分32秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (wq), wquanhit@gmail.com
 *   Organization:  
 *
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include "utils.h"
#include "index.h"
#include "query.h"
#include "bsearch.h"
#include "kvec.h"
#include "kstring.h"
#define DEBUG
typedef kvec_t(aln_t) aln_vec_t;
typedef struct{
    uint8_t *seq;
    int n;
       
} seed_candidate_t;
typedef struct{
    int trim_left, trim_right;//trim <int> bases from 5'(left)or 3'(right) of each read before alignment
    int seed_intv;//seed extract interval 
} aux_t;


typedef struct{
    int n_lext0, n_lext1, n_rext;
    ext_t *local_lext0, *local_rext;
    uint32_t *local_lext1;
    uint32_t k, l;
    bp_t *map;//auxiliary bit map
    uint32_t *seq_buf[2];
} working_data_t;

static int aln_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "vmap    [opt]  <index.prefix>  <Reads>\n\n"); 

    fprintf(stderr, "       -h\n"); 
    fprintf(stderr, "\n");
    return 1;
}

extern unsigned char nst_nt4_table[256];

const int L0 = 12;
const int L1 = 20;
const int L_SEED = 36;

int is_variant(var_t *var, int l_seq, const uint8_t *seq)
{
    char alt[128] = {};
    const uint8_t *tail;
    int i, l;
    for(i =0; i < l_seq; ++i) LOG("%c", "ACGT"[seq[i]]);
    LOG("\t%s\n", var->var);
    strcpy(alt, var->var);
    const char*vt = strtok(alt, "\\");
    l = strlen(vt);
    tail = seq+l_seq-l;
    for(i=0; i <l; ++i ){
        if(tail[i] != nst_nt4_table[(uint8_t)vt[i]]) {
            LOG("alt %s miss!\n", vt);
            break;
        }
    }
    if(i == l){ 
        
        LOG("alt %s hit!\n", vt);
        return l;
    }
    while((vt=strtok(NULL, "\\"))!= NULL){
        l = strlen(vt); 
        tail = seq+l_seq-l;
        for(i=0; i <l; ++i ){
            if(tail[i] != nst_nt4_table[(uint8_t)vt[i]]){ 
                LOG("alt %s miss!\n", vt);
                break;
            }
        }
        if(i == l) {
            LOG("alt %s hit!\n", vt);
            return l;
        }
    }

    fprintf(stderr, "no variant hit %s\n", vt);
    return -1;
}
int is_suffix(uint32_t seq16, int skip_ref, int skip_seed, int l_seq, uint8_t *seq)
{
    if(l_seq > 16) return 0;
    int i;
    seq16 >>=skip_ref*2; 
    for(i = l_seq-skip_seed-1; i >=0; --i){
        uint8_t c0, c1;
        c0 = seq16&3;
        seq16 >>=2;
        c1 = seq[i];

        if(c0 != c1) {return 0;}
             
    } 
    return 1;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  seeding0
 *  Description:  20bp variant-aware seeding(seed should be longer than 28bp)
 *  |                       |0-11 12 mer fastmap 
 *  19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0
 *                          
 *  sa_range[0] = bwt[11..0]   n_ext = 0
 *  sa_range[1] = bwt[12..0]   n_ext = 1
 *  ...
 *---------------------------------------------------------------------------------
 *                      alignment table
 *  variant preseq (6bp)| variant (2bp)|align length | template seed bwt index| 
 * =====================================================================================
 */
typedef struct{
    uint8_t preseq[8];/* variant preseq include 2bp variants */
    uint8_t var;
    int l_aln;
    uint32_t k, l;/* template seq bwt index */
    uint32_t k0, l0;/* 36bp bwt index */
    uint32_t pos;
} seed_entry_t;
typedef kvec_t(seed_entry_t) table_t;
static void log_seed_entry(seed_entry_t *e)
{
    fprintf(stderr, "[seed entry]: %u\n", e->pos);
    fprintf(stderr, "[seed entry]: k, l %u, %u\n", e->k, e->l);

}


static inline uint32_t __POPCOUNT(uint32_t v)
{
    v = v-((v>>1)&0x55555555);
    v = (v&0x33333333)+((v>>2)&0x33333333);
    return ((v+(v>>4)&0xF0F0F0F)*0x1010101)>>24;
}


#define N_TABLE 1024
int hamming_distance(const uint32_t seq16, const int l_seq, const uint8_t *seq)
{
    if(l_seq < 16) return -1;
    int i, d=0;

    for(i = 0; i < 16; ++i){
        if(__get_seq16(seq16, i) != seq[i])  ++d;      
    }
    return d;
}
int hamming_distance_alt(uint32_t seq16, const int l_seq, const uint8_t *seq)
{
    if(l_seq < 16) return -1;
    int i;
    uint32_t x = 0;
    for(i = 0; i < 16; ++i){
        x <<=2; 
        x |= (uint32_t)seq[i];
    }
    x &= seq16;
    x = x>>1&x&0x55555555;
    return 16-__POPCOUNT(x);
}
#define MAX_DIFF 2
void log_bwt_seq(bwt_t *bwt, bwtint_t k, bwtint_t l, uint8_t *pac, int n)
{
    int j;
    bwtint_t i;
    for(i = k; i <= l; ++i){
        bwtint_t pos = bwt_sa(bwt, i);
        fprintf(stderr, "%u\t", pos);       
        for(j =0; j <n;++j) fprintf(stderr, "%c", "ACGT"[__get_pac(pac, pos+j)]);
        fprintf(stderr, "\n");
    } 
}
void log_seq16(uint32_t seq16)
{
    int i;
    for(i = 0; i < 16; ++i) fprintf(stderr, "%c", "ACGT"[__get_seq16(seq16, i)]);
    fprintf(stderr, "\n");
}
void seq16tobt(uint32_t seq16, uint8_t *seq)
{
    int i;
    for(i = 0; i < 16; ++i) seq[i] = __get_seq16(seq16, i);
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_local_ext
 *  Description:  find start and end of extend index  
 * =====================================================================================
 */
int find_local_ext(idx_t *idx, bwtint_t k, uint32_t ext[][2] )
{

    int i = bp_rank(idx->is_multiseeds, k)-1; 
    assert(idx->n_tot > i);
    ext[0][0] = idx->lext_idx[i];
    ext[0][1] = idx->lext_idx[i+1];
    ext[1][0] = idx->lext1_idx[i];
    ext[1][1] = idx->lext1_idx[i+1];
    ext[2][0] = idx->rext_idx[i];
    ext[2][1] = idx->rext_idx[i+1];
    
    return 0;
    
}		/* -----  end of function find_local_ext  ----- */
int find_local_ext_alt(idx_t *idx, bwtint_t k, bwtint_t l, working_data_t *w )
{
    int i = bp_rank(idx->is_multiseeds, k)-1; 
    assert(idx->n_tot > i);
   
    w->n_lext0 = idx->lext_idx[i+1]-idx->lext_idx[i];
    w->local_lext0 = idx->lext0+idx->lext_idx[i]; 
    w->n_lext1 = idx->lext1_idx[i+1]-idx->lext1_idx[i];
    w->local_lext1 = idx->lext1+idx->lext1_idx[i]; 
    w->n_rext = idx->rext_idx[i+1]-idx->rext_idx[i];
    w->local_rext = idx->rext0+idx->rext_idx[i]; 
    w->k = k; w->l = l;
     
    return 0;
    
}		/* -----  end of function find_local_ext  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  find_pos
 *  Description:  locate pos in [k, l]
 * =====================================================================================
 */
int find_pos(bwt_t *bwt, bwtint_t k, bwtint_t l, int l_seq, uint8_t *seq, int offset, aln_vec_t *aln)
{
    bwtint_t i, x;
    for(i = k; i <l+1; ++i){
        if((x=bwt_sa(bwt, i))>(bwtint_t)offset){ 
            aln_t *p = kv_pushp(aln_t, *aln);
            p->pos = x; 
        }
    } 
    return 0;
}

uint32_t bttoseq16(const int l_seq, const uint8_t *seq)
{
    int i;
    uint32_t r = 0;
    for(i = 0; i< l_seq; ++i){
        __set_seq16(r, i, seq[i]);
    }
    return r;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  gen_ext_seed
 *  Description:  
 * =====================================================================================
 */
int gen_ext_seed ( const int l_seq[2], const uint8_t *seq[2], uint32_t ext[][2], ext_t *pext[2], const uint32_t *lext1, uint32_t *ext_buf[2], bp_t *map, int n_num, uint32_t pout_data[][8] )
{
    int n_seed = 0;
    int i, l, n_buf0=0, n_buf1=0, n_lext, n_rext, n_l2r;
    const uint8_t *q;
    /* check left 16bp seq */ 
    n_lext = ext[0][1]-ext[0][0];
    ext_t *lext = pext[0]+ext[0][0]; 
    l = l_seq[0]; q = seq[0];
    for(i = 0; i < n_lext; ++i){
        if(hamming_distance(lext[i].seq, l, q) < 2) ext_buf[0][n_buf0++] = i; 
    }
    /* check right 16bp seq*/ 
    l = l_seq[1]; q= seq[1];
    n_rext = ext[2][1]-ext[2][0];
    ext_t *rext = pext[1]+ext[2][0];
    for(i = 0; i < n_rext; ++i){
        if(hamming_distance(rext[i].seq, l, q) < 2) {
            //uint32_t seq16 = bttoseq16(l, q);
            ext_buf[1][n_buf1++] = i;
            bp_set1(map, i); 
        }
    }
    /* pairing left and right seq */
    const uint32_t *l2r = lext1+ext[1][0];
    n_l2r = ext[1][1]-ext[1][0];
    uint32_t x ,j;
    for(i = 0; i < n_buf0; ++i){
        x = ext_buf[0][i]; 
        for(j = lext[x].idx; j < (x+1 ==(uint32_t)n_lext?n_l2r:lext[x+1].idx); ++j){//[fix me]:x+1 out of bounds?
            
            if(bp_get(map, l2r[j]) == 1) {
                if(n_seed == n_num) return n_seed;//out of range
                uint32_t ri = l2r[j];
                pout_data[n_seed][0] = lext[i].seq;//left extend seq 16
                pout_data[n_seed][1] = rext[ri].seq; //right extend seq 16
                pout_data[n_seed][2] = rext[ri].idx; //sa interval of right 16 bp extend 
                pout_data[n_seed][3] = rext[ri+1].idx-1;//out of bound?
                n_seed++;
            }//add to extend seed  
        }
    }
    
    for(i = 0; i < n_buf1; ++i){ bp_set0(map, ext_buf[1][i]); }
    return n_seed;
}		/* -----  end of function gen_ext_seed  ----- */
#define MAX_DISTANCE 0
int gen_ext_seed_alt( const int l_seq[2], const uint8_t *seq[2], working_data_t *working_pool, int n_num, uint32_t (*pout_data)[8] )
{
    int n_seed = 0;
    int i, l, n_buf0=0, n_buf1=0, n_lext, n_rext, n_l2r;
    ext_t *lext, *rext;uint32_t *l2r; 
    uint32_t *ext_buf[2];
    ext_buf[0] = working_pool->seq_buf[0];
    ext_buf[1] = working_pool->seq_buf[1];
    bp_t *map = working_pool->map;
    const uint8_t *q;
    /* check left 16bp seq */ 
    n_lext = working_pool->n_lext0;
    lext = working_pool->local_lext0; 
    l = l_seq[0]; q = seq[0];
    for(i = 0; i < n_lext; ++i){
        if(hamming_distance(lext[i].seq, l, q) <= MAX_DISTANCE) ext_buf[0][n_buf0++] = i; 
    }
    /* check right 16bp seq*/ 
    l = l_seq[1]; q= seq[1];
    n_rext = working_pool->n_rext;
    rext = working_pool->local_rext;
    for(i = 0; i < n_rext; ++i){
        if(hamming_distance(rext[i].seq, l, q) <=MAX_DISTANCE) {
            //uint32_t seq16 = bttoseq16(l, q);
            ext_buf[1][n_buf1++] = i;
            bp_set1(map, i); 
        }
    }
    /* pairing left and right seq */
    n_l2r = working_pool->n_lext1;
    l2r = working_pool->local_lext1;

    uint32_t x, j;
    for(i = 0; i < n_buf0; ++i){
        x = ext_buf[0][i]; 
        for(j = lext[x].idx; j < (x+1 ==(uint32_t)n_lext?n_l2r:lext[x+1].idx); ++j){//[fix me]:x+1 out of bounds?
            if(bp_get(map, l2r[j]) == 1) {
                if(n_seed == n_num) return n_seed;//out of range
                uint32_t ri = l2r[j];
                pout_data[n_seed][0] = lext[i].seq;//left extend seq 16
                pout_data[n_seed][1] = rext[ri].seq; //right extend seq 16
                pout_data[n_seed][2] = rext[ri].idx; //sa interval of right 16 bp extend 
                pout_data[n_seed][3] = ri==n_rext-1?working_pool->l:rext[ri+1].idx-1;//out of bound?
                n_seed++;
            }//add to extend seed  
        }
    }
    for(i = 0; i < n_buf1; ++i){
        bp_set0(map, ext_buf[1][i]); 
    }
    return n_seed;
}		/* -----  end of function gen_ext_seed  ----- */



int test_seeding_approx(idx_t *idx, int l_seed, uint8_t *seed, table_t *table, int start, int end)
{
    bp_t *m = idx->is_multiseeds;
    uint32_t *iter_seq16 = idx->isa2seq16; 
    uint32_t *seq16s = idx->seq16s;
    int64_t i, j, d;
    bwtint_t k ,l, iter;
    for(i = start; i <end; ++i){
        seed_entry_t *e = table->a+i;
        log_seed_entry(e);
        k = e->k, l = e->l; 
        
        if(k == l) {//uniq seed
            fprintf(stderr, "%u-16\t", bwt_sa(idx->bwt, k)); 
            continue;
        }
        if(l < k) continue;//not hit
        //log_bwt_seq(idx->bwt, k-1, l+1, idx->pac, 20);
         
        fprintf(stderr, "multi seeds [%u ,%u]\t", k, l); 
        if(bp_get(m, k)!= 1) fprintf(stderr, "%s multi_seeds flag is wrong!!!\n", __func__ );
        /* extend left 16mer */ 
        iter = bp_rank(m, k)-1;
        for(j = iter_seq16[iter]; j < (int)iter_seq16[iter+1]; ++j){
            uint32_t seq16 = seq16s[j];
            log_seq16(seq16);
            d = hamming_distance(seq16, 16, seed); 
            if(d >2) continue;
            /* print seed36 */ 
            uint8_t seq16bt[16] = {};
            for(i = 0; i < 16; ++i) seq16bt[i] = __get_seq16(seq16, i);      

            bwtint_t __k = k, __l= l;
            bwt_match_exact_alt(idx->bwt, 16, seq16bt, &__k, &__l);
            for(i = __k; i <= __l; ++i){
                fprintf(stderr, "%u\t", bwt_sa(idx->bwt, i)); 
            } 
        }
    
       
    }
    fprintf(stderr, "\n"); 
    return 0;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  seeding_approx
 *  Description:  l_seed = 36 extend 20bp seed to 36bp(left 16bp) seed
 * =====================================================================================
 */
int seeding_approx(idx_t *idx, int l_seed, uint8_t *seed, table_t *table, int start, int end)
{
    bp_t *m = idx->is_multiseeds;
    uint32_t *iter_seq16 = idx->isa2seq16; 
    uint32_t *seq16s = idx->seq16s;
    int64_t i, j, d;
    bwtint_t k ,l, iter;
    for(i = start; i <end; ++i){
        seed_entry_t *e = table->a+i;
        k = e->k, l = e->l; 
        if(k >= l) {continue;}/* hit number <=1 */
        if(bp_get(m, k)!= 1) fprintf(stderr, "%s multi_seeds flag is wrong!!!\n", __func__ );
        /* extend left 16mer */ 
        iter = bp_rank(m, k)-1;
        for(j = iter_seq16[iter]; j < (int)iter_seq16[iter+1]; ++j){
            uint32_t seq16 = seq16s[j];
            //log_seq16(seq16);
            d = hamming_distance(seq16, 16, seed); 
            if(d >2) continue;
            /* print seed36 */ 
            uint8_t seq16bt[16] = {};
            for(i = 0; i < 16; ++i) seq16bt[i] = __get_seq16(seq16, i);      
            bwtint_t __k = k, __l= l;
            bwt_match_exact_alt(idx->bwt, 16, seq16bt, &__k, &__l);
        }
    }
    return 0;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  seeding_approx_alt
 *  Description:  l_seed = 52
 * =====================================================================================
 */
#define MAX_SEED_NUM 500000
#define MAX_ALN 10000
extern void log_bt2nt(int l, uint8_t *bt);
static inline void working_data_init(working_data_t *w)
{
    w->n_lext0 = 0;
    w->n_lext1 = 0;
    w->n_rext = 0;
    w->local_rext = NULL;
    w->local_lext0 = NULL;
    w->local_lext1 = NULL;

    w->seq_buf[0] = (uint32_t *)calloc(MAX_SEED_NUM, sizeof(uint32_t));
    w->seq_buf[1] = (uint32_t *)calloc(MAX_SEED_NUM, sizeof(uint32_t));
    w->map = bp_init(MAX_SEED_NUM);//auxiliary bit map to pair left extend seq and right extend seq
    return;
}
static inline void working_data_destroy(working_data_t *w)
{
    free(w->seq_buf[0]);
    free(w->seq_buf[1]);
    bp_destroy(w->map);
}
int seeding_approx_alt(idx_t *idx, int l_seed, uint8_t *seed, table_t *table, int start, int end)
{
    int64_t i, j;
    bwtint_t k ,l, iter;
    bp_t *m = idx->is_multiseeds;
    
    aln_vec_t aln;//alignment  
    working_data_t working_pool; 
    uint32_t (*extend_seeds)[8]; 

    kv_init(aln); kv_resize(aln_t, aln, MAX_ALN); 
    working_data_init(&working_pool);
    extend_seeds =  (uint32_t(*)[8])malloc(MAX_SEED_NUM*sizeof(uint32_t)*8);//extend_seed
    
    for(i = start; i <end; ++i){
        int l_seq[2] = {16, 16};
        uint8_t *seq[2];
      
        seed_entry_t *e = table->a+i;
        k = e->k, l = e->l; 
        if(k >= l) {
            //fprintf(stderr, "[sa interval]:\t%u\t%u\n", k, l);
            find_pos(idx->bwt, k, l, 0, NULL, 0, &aln);
            continue;
        }/* hits number <=1 */
        if(bp_get(m, k)!= 1) fprintf(stderr, "%s multi_seeds flag is wrong!!!\n", __func__ );
        /* extend left 16mer */ 
        iter = bp_rank(m, k)-1;
        find_local_ext_alt(idx, k, l, &working_pool);
        seq[0] = seed; seq[1] = seed+36;
        int n = gen_ext_seed_alt(l_seq, seq, &working_pool, MAX_SEED_NUM, extend_seeds);
        for(j =0; j < n; ++j) {
            uint8_t seq[52] = {};
            find_pos(idx->bwt, extend_seeds[j][2], extend_seeds[j][3], 0, NULL, 0, &aln);
            fprintf(stderr, "36-range(%u, %u)\n", extend_seeds[j][2], extend_seeds[j][3]);
            uint32_t __k= extend_seeds[j][2], __l=extend_seeds[j][3];
            bwt_match_exact_alt(idx->bwt, 16, seed, &__k, &__l);
            fprintf(stderr, "52-range(%u, %u)\n", __k, __l);
            memcpy(seq, seed+16, 20);
        }
    }
    if(aln.n != 0) fprintf(stderr, "\n[locate pos]:\t");
    for(i = 0; i < aln.n; ++i){
        //fprintf(stderr, "%u\t", aln.a[i].pos);
    }
    fprintf(stderr, "\n");
    
    kv_destroy(aln);
    free(extend_seeds);
    working_data_destroy(&working_pool);
    return 0;
}

void snpaware_extend(idx_t *idx, int l_seed, uint8_t *seed, const int n_ext, const bwtint_t k, const bwtint_t l, table_t *table)
{
    int i, j, n;
    uint8_t c;
    uint32_t preseq = 0;
    int preseq_start = l_seed-n_ext-L0-8;
    int preseq_end = preseq_start+6;
    if(preseq_start < 0) {
        fprintf(stderr, "preseq_start <0\n");
        exit(1);
    }
    for(j = preseq_start; j <preseq_end; ++j ){
        preseq <<=2;
        preseq |= seed[j];
    }
    /* find pattern like [6mer][NN][n_ext][12mer] */
    bs_iter_t x, y,  x1, y1, __i; 
    x = idx->pmap[preseq];
    y = idx->pmap[preseq+1];
    if(x >=y)  return;
    x1 = lower_bound(x, y, idx->sai, k);
    y1 = uper_bound(x, y, idx->sai, l);
    LOG("bwt index [%u, %u]\n", k, l);
    for(__i = x; __i < y; ++__i){
        LOG("%u\t",idx->sai[__i]); 
    }
    LOG("\n%u %u\n", x1, y1); 
    if(x1 >= y1) return;
    /* fill NN*/
    uint8_t var_map[16] = {};
    memset(var_map, 0, sizeof(uint8_t)*16);
    for(n= 0, __i = x1;__i< y1; ++__i ){
        c = idx->refseq[__i]; 
        if(var_map[c] == 0) {
            var_map[c] = 1;
            ++n;
            if(n >= 16) break;
        }
    }
    var_map[seed[preseq_end]<<2|seed[preseq_end+1]] = 0;
    for(i=0; i<16; ++i){
        if(var_map[i] == 0) continue;
        seed_entry_t *iter_table = kv_pushp(seed_entry_t, *table);
        iter_table->l_aln = n_ext+L0;
        iter_table->var = i;
        /* generate template seq before  */  
        uint8_t *pseq = iter_table->preseq;
        for(n =0, j = preseq_start; j < preseq_end; ++j, ++n){//move to the loop may be more efficient
            pseq[n] = seed[j];
        }
        pseq[n++] = i>>2;
        pseq[n] = i&3;
        //bwt_match_exact_alt(bwt, 8, pseq, &k, &l);
        bwtint_t ok = k, ol = l;
        bwt_match_exact_alt(idx->bwt, 8-n_ext, pseq+n_ext, &ok, &ol);
        iter_table->k = ok;
        iter_table->l = ol;
        iter_table->pos = bwt_sa(idx->bwt, ok); 
    }
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  seeding_snpaware
 *  Description:  find 20 bp snp-aware seed
 *      @l_seed:  length of seed, >=28, (20 bp snpaware seed and 8 bp preseq)
 *      @seed  :  tail 20 bp is snpaware seed
 *      @table :  output table
 * =====================================================================================
 */
int seeding_snpaware(idx_t *idx, int l_seed, uint8_t *seed, aux_t *aux, table_t *table)
{
    lkt_t *fastmap = idx->fastmap;/* fast map 12mer suffix */
    bwt_t *bwt  = idx->bwt;

    seed_entry_t *iter_table;  
    bwtint_t k, l, ok, ol;
    bwtint_t sa_range[16] = {};/* sa_range[0] store bwt range of last 12mer seq*/
    uint32_t n_ext;
    uint8_t c;
    int i;   
    
    k = 1; l = bwt->seq_len;
    if(lkt_lookup_sa(fastmap, seed, l_seed-L0, l_seed-1, &k, &l) > 0){/*tail 12-mer exact match*/
        //fprintf(stderr, "[lookuptable]: [%u, %u]->bwt[%u, %u]\n", l_seed-L0, l_seed-1, k, l); 
        /*if(bwt_match_exact_alt(bwt, L1-L0, seed+l-L1, &k, &l)>0){ }*/
        for(i = l_seed-L0-1, n_ext=0; i >= l_seed-L1; --i) {//extend 12-mer to 20mer seed
            sa_range[n_ext*2] = k;
            sa_range[n_ext*2+1] = l; 
            c = seed[i];
            bwt_2occ(bwt, k-1, l, c, &ok, &ol);
            k = bwt->L2[c] + ok +1;
            l = bwt->L2[c] + ol;
            if(k > l) break;
            ++n_ext;

        }
        //fprintf(stderr, "[backward search]: [%u, %u]->bwt[%u, %u]\n", l_seed-L0-1, l_seed-L0-n_ext, k, l); 
        if(n_ext == L1-L0){//20-mer exact match
            //add seed to template seeds vector 
            iter_table = kv_pushp(seed_entry_t, *table);
            iter_table->l_aln = L1;
            iter_table->k = k; iter_table->l = l;
            iter_table->pos = bwt_sa(bwt, k);
            //printf("%u\n", iter_table->pos);
        } else{//20-mer variant-aware match
            for(i = n_ext; i >= 0; --i){
                k = sa_range[i*2]; l = sa_range[i*2+1]; 
                //fprintf(stderr, "[backtrack %u]: bwt[%u, %u]\n", L0+i,  k, l); 
                /* look up variants table and find i(which is the first sai higher than k), j(which is the first sai lower than l)*/
                snpaware_extend(idx, l_seed, seed,i, k, l, table); 
            }   
        }
    } else{//12-mer exact match fail
    
    }
    return table->n;
}		/* -----  end of function seed  ----- */


static inline void accumulate_cnt(int n, uint8_t *cnt_2nt)
{
    int i;
    for(i = 1; i < n-1; ++i) cnt_2nt[i] += cnt_2nt[i-1]; 
    for(i = n-1; i >0; --i) cnt_2nt[i] = cnt_2nt[i-1];
    cnt_2nt[0] = 0;
}
void gen_cnt(uint8_t *bwt, int n_data, uint8_t cnt_2nt[][17])
{
    uint8_t rot, ch, bg_row, ed_row, i;
    for(rot=0; rot<8; ++rot){
        for(ch=0; ch<17; ++ch) {cnt_2nt[rot][ch] = 0;} 
    }
    for(rot =0; rot<8; ++rot){
        bg_row = rot*((n_data+1)/2);
        ed_row = bg_row +n_data/2; 
        for(i = bg_row; i < ed_row; ++i){
            ch = bwt[i]>>4;
            ++cnt_2nt[rot][ch];
            ch = bwt[i]&0xF;
            ++cnt_2nt[rot][ch]; 
        }
        if(n_data%2 != 0){
            ch = bwt[ed_row]>>4;
            ++cnt_2nt[rot][ch];
        }
        accumulate_cnt(17, cnt_2nt[rot]);
    }
}
int count_bwt(uint8_t *bwt, int k, int l, uint8_t nt2)
{
    int i, n= 0;
    for(i = (k+1)/2; i < (l-1)/2; ++i){
        n += ((bwt[i]&0xF)==nt2);
        n += ((bwt[i]>>4)==nt2);
    } 
    if(k%2 != 0) n += ((bwt[k/2]&0xF)==nt2);  
    if(l%2 == 0) n += ((bwt[l/2]>>4)==nt2);  
    return n;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  align_128
 *  Description:  aln n_data < 128
 *           
 * =====================================================================================
 */
void align_128(uint8_t *bwt, uint8_t seq_2nt[8], int seq_st, uint8_t cnt_2nt[][17], int n_data)
{
    int stat_flag, mid_num, top_num, bot_num, top_row, bot_row, row;
    uint32_t bg_idx, ed_idx, i, j;
    uint8_t ch_seq, ch;
    bg_idx = 0; ed_idx = n_data;
    i = seq_st;
    int cur_len = 1;
    int rot = (8-seq_st)%8;
    while(1){
        top_row = (bg_idx+1)/2;
        bot_row = (ed_idx -1)/2;
        i %= 8;
        ch_seq = seq_2nt[i*2]<<4| seq_2nt[i*2+1];
        mid_num =0; top_num = 0; bot_num = 0; 
        if(bg_idx <= n_data-ed_idx){
            stat_flag = 1;
            if(ed_idx%2==0){
                ch = bwt[ed_idx/2]>>4;
                if(ch_seq == ch) ++mid_num;
            }
            for(j = (bg_idx+1)/2; j < (ed_idx-1)/2; ++j){
                ch = bwt[j] >>4;
                if(ch_seq == ch) ++mid_num;
                ch = bwt[j]&0xF;
                if(ch_seq == ch) ++mid_num;
                if(n_data == mid_num) break;
            }
            if(bg_idx%2>0){
                ch = bwt[bg_idx/2];
                if(ch_seq==ch) ++mid_num;
                ch = bwt[bg_idx/2]>>4;
                if(ch_seq == ch) ++top_num;
                --top_row;//???
            }
            if(n_data > mid_num && mid_num > 0) {
                for(row = top_row-1; row>=0; --row){
                    ch = bwt[row] >>4;
                    if(ch_seq==ch) ++top_num;
                    ch = bwt[row] &0xF;
                    if(ch_seq==ch) ++top_num; 
                
                } 
            }
        } else{ //bg_idx > n_data-ed_idx
            stat_flag = 2;
            if(bg_idx%2 != 0){
                ch = bwt[bg_idx/2]&0xF;
                if(ch_seq == ch) ++mid_num;
            }
            for(j = (bg_idx+1)/2; j < (ed_idx-1)/2; ++j){
                ch = bwt[j] >>4;
                if(ch_seq == ch) ++mid_num;
                ch = bwt[j]&0xF;
                if(ch_seq == ch) ++mid_num;
                if(n_data == mid_num) break;
            }
            if(ed_idx%2>0){
                ch = bwt[bg_idx/2]>>4;
                if(ch_seq==ch) ++mid_num;
                ch = bwt[bg_idx/2]&0xF;
                if(ch_seq == ch) ++bot_num;
                ++bot_row;//???
            }
            if(n_data > mid_num && mid_num > 0) {
                for(row = top_row-1; row>=0; --row){
                    ch = bwt[row] >>4;
                    if(ch_seq==ch) ++bot_num;
                    ch = bwt[row] &0xF;
                    if(ch_seq==ch) ++bot_num; 
                
                }
                if(n_data%2>0){
                    ch = bwt[row]>>4;
                    if(ch_seq==ch) ++bot_num;    
                } 
            }
       
        }
        if(mid_num >0 && cur_len <8){
            if(stat_flag == 1){
                bg_idx = cnt_2nt[rot][ch]+top_num;
                ed_idx = cnt_2nt[rot][ch+1]+top_num+mid_num;
            } else if(stat_flag == 2){
                bg_idx = cnt_2nt[rot][ch]-bot_num-mid_num+1;
                ed_idx = cnt_2nt[rot][ch+1] -bot_num;
            
            }
            bg_idx += rot*(n_data+1)/2;
            ed_idx += rot*(n_data+1)/2; 
            rot = (rot+1)%8;
            ++cur_len; 
        }
        if(mid_num ==0){
            if(cur_len == 7){
            
            } else if(cur_len <7){
            
            } else if(cur_len ==8){
            
            }
        
        
        } 
    }//end while(1) 
    
}

int alnse(idx_t *idx, query_t *q, aux_t *aux)
{
    int i, j;   
    // print_query(q);
    
    /* seeding */
    table_t *seed_table = (table_t *)calloc(1, sizeof(table_t));
    kv_resize(seed_entry_t, *seed_table, N_TABLE);
    fprintf(stderr, "\n***************************************************\n");    
    fprintf(stderr, "seq name:\t%s\n", q->name);    
    //fprintf(stderr, "seq\t%s\n", q->seq);    
    for(i=0; i < 52; ++i) fprintf(stderr, "%c", "ACGT"[q->seq[i]]);

    fprintf(stderr, "\n");
    for(i = aux->trim_left; i < q->l_seq-aux->trim_right-(L_SEED+16)+1; i+= aux->seed_intv){
        fprintf(stderr, "seeding offset\t%u\n",i);    
        /* seq */
        int start = seed_table->n;
        seeding_snpaware(idx, L_SEED, q->seq+i, aux, seed_table);//right 20bp snpaware seeding
        //seeding_approx(idx, L_SEED, q->seq+i, seed_table, start, seed_table->n);//left 16bp approximate seeding
        for(j =start; j < seed_table->n; ++j) printf("20 range(%u, %u)\n", seed_table->a[j].k, seed_table->a[j].l);
        seeding_approx_alt(idx, 52, q->seq+i, seed_table, start, seed_table->n);//left 16bp approximate seeding
        
        /* reverse complement */
        start = seed_table->n;    
        seeding_snpaware(idx, L_SEED, q->rseq+i, aux, seed_table);//right 20bp snpaware seeding
        //seeding_approx(idx, L_SEED, q->rseq+i, seed_table, start, seed_table->n);//left 16bp approximate seeding
        seeding_approx_alt(idx, 52, q->rseq+i, seed_table, start, seed_table->n);//left 16bp approximate seeding
    }
    kv_destroy(*seed_table); 
    free(seed_table);
    /* extending */
    /* aln to sam */
    return 0;
}
int aln_main(int argc, char *argv[])
{
    int c;
    while((c = getopt(argc, argv, "k:h"))>=0){
        switch(c)
        {
            case 'h':
                return aln_usage();
            default:
                return 1;
        }
    }
    if (argc - optind  != 2 ){
        return aln_usage();
    }
    const char* fn_idx = argv[optind];
    const char* fn_seq = argv[optind+1];
    idx_t *idx = idx_restore(fn_idx);
    queryio_t *fp_seq = query_open(fn_seq); 
    query_t *multi_seqs = (query_t *)calloc(N_SEQS, sizeof(query_t)); 
    int n, i;
    aux_t aux;
    aux.trim_left = 0;
    aux.trim_right = 0;
    aux.seed_intv = 10;
    while((n = query_read_multiSeqs(fp_seq, N_SEQS, multi_seqs))>0){
        for(i =0; i < n; ++i){
            alnse(idx, multi_seqs+i, &aux); 
        }
        for(i =0; i < n; ++i) query_destroy(multi_seqs+i); 
    } 
    query_close(fp_seq);
    free(multi_seqs); 
    idx_destroy(idx);
    return 0;
}
#ifdef ALN_MAIN
int main(int argc, char *argv[])
{

    return aln_main(argc, argv);
}

#endif
