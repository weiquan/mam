/*
 * =====================================================================================
 *
 *       Filename:  index1.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/11/2013 02:29:22 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <time.h>
#include <assert.h>

#include "utils.h"
#include "khash.h"
#include "index.h"
#include "kvec.h"
#include "ksort.h"
#include "setFileName.h"
#include "malloc_wrap.h"
#define USE_MALLOC_WRAPPERS
#define MAX_NAME 128
const int SA_INTV = 1;

int *global_stat_20seed0;
int *global_stat_20seed1;
int *global_stat_36seed;

KHASH_MAP_INIT_STR(str, int)
typedef kvec_t(uint32_t) vec_uint_t;
typedef kvec_t(ext_t) vec_ext_t;
//#define __sort_lt_ext(a, b)((a).seq == (a).seq?(a).idx<(b).idx:(a).seq < (b).seq ) 
#define __sort_lt_ext(a, b)((a).seq < (b).seq||((a).seq==(b).seq&&(a).idx < (b).idx )) 

KSORT_INIT_GENERIC(uint32_t)
KSORT_INIT(ext, ext_t, __sort_lt_ext)


#define MAX_COUNT 800000
#define SWAP(type_t, a, b) do{type_t x=(a); (a)=(b); (b)=x;} while(0)
#define __set_bwt(bwt, l, c) ((bwt)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define __get_bwt(bwt, l) ((bwt)[(l)>>2]>>((~(l)&3)<<1)&3)
#define __set_bwt2(bwt, l, c) ((bwt)[(l)>>1] |= (c)<<((~(l)&1)<<2))
#define __get_bwt2(bwt, l) ((bwt)[(l)>>1] >> ((~(l)&1)<<2)&15)


#define __get_col(seq, i) (((seq)>>((i)*4))&0xF)

static int index_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "index    [OPT]    <ref.fa>    <variants.txt>    <prefix>\n\n"); 
    fprintf(stderr, "         -h\n"); 
    fprintf(stderr, "\n");
    return 1;
}
/*****************
 * for debug
 * **************/
int log_sa(bwtint_t n_sa, bwtint_t *sa)
{
    bwtint_t i;
    for(i = 0; i < n_sa; ++i){
        printf("%u\t%u\n", i, sa[i]);
    
    }
    return 0;
}
void log_bt2nt(int l, uint8_t *bt)
{
    int i;
    for(i = 0; i < l; ++i){
        printf("%c", "ACGT"[bt[i]]); 
    }
    printf("\n");
}
void log_seq162nt(uint32_t x)
{
    uint32_t i;
    char seq[17] = {};
    for(i = 0; i < 16; ++i){
        seq[i] = "ACGT"[x>>(30-i*2)&3] ;
    } 
    fprintf(stderr, "%s\n", seq);
}

void restore_extend_idx(idx_t *idx, const char *prefix)
{
    int i;
    char fn_rext0[MAX_NAME], fn_totext[MAX_NAME];
    char fn_lext0[MAX_NAME], fn_lext1[MAX_NAME];
    strcat(strncpy(fn_lext0, prefix, MAX_NAME), ".lext0");
    strcat(strncpy(fn_lext1, prefix, MAX_NAME), ".lext1");
    strcat(strncpy(fn_rext0, prefix, MAX_NAME), ".rext0");
    strcat(strncpy(fn_totext, prefix, MAX_NAME), ".totext");
    FILE *fp_lext0 = xopen(fn_lext0, "r"); 
    FILE *fp_lext1 = xopen(fn_lext1, "r"); 
    FILE *fp_rext0 = xopen(fn_rext0, "r"); 
    FILE *fp_totext = xopen(fn_totext, "r"); 
    //restore totext
    fseek(fp_totext, 0, SEEK_END);
    int n_tot = idx->n_tot = ftell(fp_totext)/sizeof(int)/3; fseek(fp_totext, 0, SEEK_SET); 
    idx->lext_idx = (uint32_t *)calloc(n_tot, sizeof(uint32_t)); 
    idx->rext_idx = (uint32_t *)calloc(n_tot, sizeof(uint32_t)); 
    idx->lext1_idx = (uint32_t *)calloc(n_tot, sizeof(uint32_t)); 
    for(i = 0; i < n_tot; ++i){
        fread(idx->lext_idx+i, sizeof(uint32_t), 1, fp_totext); 
        fread(idx->lext1_idx+i, sizeof(uint32_t), 1, fp_totext); 
        fread(idx->rext_idx+i, sizeof(uint32_t), 1, fp_totext);  
    }
    fseek(fp_rext0, 0, SEEK_END);
    int n_rext = idx->n_rext= ftell(fp_rext0)/sizeof(ext_t); fseek(fp_rext0, 0, SEEK_SET);
    idx->rext0 = (ext_t *)calloc(n_rext, sizeof(ext_t));
    fread(idx->rext0, sizeof(ext_t), n_rext, fp_rext0); 

    fseek(fp_lext0, 0, SEEK_END);
    int n_lext0 = idx->n_lext0 = ftell(fp_lext0)/sizeof(ext_t); fseek(fp_lext0, 0, SEEK_SET);
    idx->lext0 = (ext_t *)calloc(n_lext0, sizeof(ext_t));
    fread(idx->lext0, sizeof(ext_t), n_lext0, fp_lext0); 
    
    fseek(fp_lext1, 0, SEEK_END);
    int n_lext1 = idx->n_lext1 = ftell(fp_lext1)/sizeof(uint32_t); fseek(fp_lext1, 0, SEEK_SET);
    idx->lext1 = (uint32_t *)calloc(n_lext1, sizeof(uint32_t));
    fread(idx->lext1, sizeof(uint32_t), n_lext1, fp_lext1); 
    
    
    
    err_fclose(fp_lext0);
    err_fclose(fp_lext1);
    err_fclose(fp_rext0);
    err_fclose(fp_totext);


}
idx_t *idx_restore(const char *prefix)
{
    idx_t *idx = (idx_t *)calloc(1, sizeof(idx_t)); 
    
    char fn0[MAX_NAME];
  

    /* restore pac and ann */ 

    idx->bns = bns_restore(prefix);  
    idx->pac = (uint8_t *)calloc(idx->bns->l_pac/4+1, sizeof(uint8_t));

    fread(idx->pac, sizeof(uint8_t), idx->bns->l_pac/4+1, idx->bns->fp_pac);
    //err_fclose(fp_pac);
    /* restore FM-index */
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".bwt");
    idx->bwt = bwt_restore_bwt(fn0); 
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".sa");
    bwt_restore_sa(fn0, idx->bwt); 
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".isa");
    //bwt_restore_isa(fn0, idx->bwt); 
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".lkt");
    idx->fastmap = lkt_restore(fn0); 
    /* restore snp-aware index */ 
/*
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".vseq");
    FILE *fp_varseq = xopen(fn0, "r");
    fread(&idx->n_pmap, sizeof(uint32_t), 1, fp_varseq);
    idx->pmap = (uint32_t *)calloc(idx->n_pmap, sizeof(uint32_t));
    fread(idx->pmap, sizeof(uint32_t), idx->n_pmap, fp_varseq);
    fread(&idx->n_var, sizeof(uint32_t), 1, fp_varseq);
    idx->sai = (uint32_t *) calloc(idx->n_var, sizeof(uint32_t));
    idx->refseq = (uint8_t *) calloc(idx->n_var, sizeof(uint8_t));
    fread(idx->sai, sizeof(uint32_t), idx->n_var, fp_varseq);
    fread(idx->refseq, sizeof(uint8_t), idx->n_var, fp_varseq);
    err_fclose(fp_varseq);
*/   
    //strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".repeat");
    //idx->is_multiseeds = bp_restore(fn0, 1);
    //restore_extend_idx(idx, prefix);
/*
    FILE *fp;   
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".idxseq16");
    fp = xopen(fn0, "r");
    fread(&idx->n_isa2seq16, sizeof(int), 1, fp);
    idx->isa2seq16 = (uint32_t*) calloc(idx->n_isa2seq16, sizeof(uint32_t));
    fread(idx->isa2seq16, sizeof(uint32_t), idx->n_isa2seq16, fp);
    fclose(fp);
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".seq16");
    fp = xopen(fn0, "r");
    idx->n_seq16s = idx->isa2seq16[idx->n_isa2seq16-1];
    idx->seq16s = (uint32_t*) calloc(idx->n_seq16s, sizeof(uint32_t));
    fread(idx->seq16s, sizeof(uint32_t), idx->n_seq16s, fp); 
    fclose(fp); 
*/
    return idx;
}
void idx_destroy(idx_t *idx)
{
    if(idx->fastmap != NULL ) lkt_destroy(idx->fastmap);
    bwt_destroy(idx->bwt);
    bns_destroy(idx->bns);
    if(idx->pac != NULL ) free(idx->pac);
/*
    free(idx->pmap);
    free(idx->refseq);
    free(idx->sai);

    free(idx->lext0);
    free(idx->lext1);
    free(idx->rext0);    
    free(idx->lext_idx);
    free(idx->rext_idx);
    free(idx->lext1_idx); 
    bp_destroy(idx->is_multiseeds);
*/
    //free(idx->isa2seq16);
    //free(idx->seq16s);
    //bp_destroy(idx->is_var);
    //free(idx->var_seqs);
    //variants_destroy(idx->var);
    free(idx);
}
/* 
 *  cut 16 bp seqs from reference , begin pos in [SA[i] +offset for i in range(k, l)]
 *
 *
 * */

int cut_ext_seq(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, vec_uint_t *seq_buf)
{
    bwtint_t i, pos, st, ed;
    uint32_t seq16;
    //fprintf(stderr, "\n"); 
    /* fetch 16 mer preseq in [k, l) and rm duplicate preseq */
    for(i = k; i < l; ++i){
        pos = bwt_sa(bwt, i);
        st = pos+offset; ed = st+16;
        if(offset < 0 && pos < (bwtint_t)-offset) st = 0;
        if(offset > 0 && ed > bwt->seq_len) ed = bwt->seq_len; 
        seq16 = bns_extract_seq16(pac, st, ed);//???
        //fprintf(stderr, "%u\t", seq16); 
        kv_push(uint32_t, *seq_buf, seq16);
    }
    //fprintf(stderr, "\n");
    //for(i = seq16_st; i <seq16s->n; ++i) fprintf(stderr, "%u\t", seq16s->a[i]);
    //fprintf(stderr, "\n");
    return seq_buf->n;
}

int rem_euqal(vec_uint_t *seq_buf, bwtint_t k, vec_ext_t *rext)
{
    int i, n_uniq;
    uint32_t *b = seq_buf->a;
    
    ext_t *p = kv_pushp(ext_t, *rext);
    p->seq = b[0];
    p->idx = k;
 
    n_uniq = 1; 
    for(i = 1; i < seq_buf->n; ++i){
        if(b[i] > b[i-1]) {
            assert(rext->m >n_uniq*2+1);
            ext_t *p = kv_pushp(ext_t, *rext);
            p->seq = b[i];
            p->idx = k+i;
            //a[n_uniq*2] = b[i];
            //a[n_uniq*2+1] = k+i;          
            ++n_uniq;
        } 
    }
    return n_uniq;
}


int filter_ext_seq(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, vec_uint_t *seq16s, bp_t *is_repeat)
{
    
    int seq16_st = seq16s->n;
    bwtint_t i, pos;
    uint32_t seq16;
    //fprintf(stderr, "\n"); 
    /* fetch 16 mer preseq in [k, l) and rm duplicate preseq */
    for(i = k; i < l; ++i){
        pos = bwt_sa(bwt, i);
        if(pos <16) continue; 
        seq16 = bns_extract_seq16(pac, pos-16, pos);//???
        //fprintf(stderr, "%u\t", seq16); 
        if(bp_get(is_repeat, seq16) == 1) continue;//skip if occured before
        bp_set1(is_repeat, seq16);
        kv_push(uint32_t, *seq16s, seq16);
    }
    ks_introsort(uint32_t, seq16s->n-seq16_st, seq16s->a+seq16_st);/* sort seq16 */
    //fprintf(stderr, "\n");
    //for(i = seq16_st; i <seq16s->n; ++i) fprintf(stderr, "%u\t", seq16s->a[i]);
    //fprintf(stderr, "\n");
    for(i=seq16_st; i < seq16s->n; ++i) bp_set0(is_repeat, seq16s->a[i]);/* reset is_repeat */
    ++global_stat_20seed0[l-k+1]; 
    ++global_stat_36seed[seq16s->n-seq16_st]; 
    global_stat_20seed1[l-k+1] += seq16s->n -seq16_st+1; 
    return seq16s->n;
}
#define MASK_2NT 0xF
#define MASK_6NT 0xFFF
#define L_MID_SEED 20
const int HASH_SIZE = 4097;


static void cat_filename(char *str, const char *prefix, const char *suffix)
{
    strcat(strncpy(str, prefix, MAX_NAME), suffix);

}
static void dump_cnt(FILE *fp, int tot_lext0, int tot_lext1, int tot_rext)
{
    fwrite(&tot_lext0, sizeof(uint32_t), 1, fp);
    fwrite(&tot_lext1, sizeof(uint32_t), 1, fp);
    fwrite(&tot_rext, sizeof(uint32_t), 1, fp);
 
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_fmidx
 *  Description:  build fm-index (include sa and isa)
 *       @fn_fa:  file name of reference
 *      @prefix:  file name prefix of index
 * =====================================================================================
 */
void build_fmidx(const char *fn_fa, const char *prefix)
{
    double t;
    char fn0[MAX_NAME], fn1[MAX_NAME];

    gzFile fp_fa;
    t = realtime();
    fp_fa = gzopen(fn_fa, "r");
    fprintf(stderr, "\n[variant_index]: convert genome %s  to packedfile!\n",fn_fa); 
    bns_fasta2bntseq(fp_fa, prefix);
    fprintf(stderr, "fa2pac %.2lf seconds elapse.\n", realtime() - t);
    gzclose(fp_fa);   
    
    t = realtime(); 

    strcat(strncpy(fn0, prefix, MAX_NAME), ".pac");
    strcat(strncpy(fn1, prefix, MAX_NAME), ".bwt");
    fprintf(stderr, "\n[va_index]: convert pac %s to bwt %s\n", fn0, fn1);

    bwt_bwtgen(fn0, fn1); 
    fprintf(stderr, "pac2bwt %.2lf seconds elapse.\n", realtime() - t);
    {    
        bwt_t *bwt;
        t = realtime();
        strcat(strncpy(fn0, prefix, MAX_NAME), ".bwt");
        fprintf(stderr, "[va_index]: Update BWT... \n");
        bwt = bwt_restore_bwt(fn0);
        bwt_bwtupdate_core(bwt);
        bwt_dump_bwt(fn0, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "update bwt %.2f sec\n", (float)(realtime() - t));
    }
    {   
        char fn_lkt[MAX_NAME];
        t = realtime();
        fprintf(stderr, "[va_index]: generate look up table... \n");
        strcat(strncpy(fn_lkt, prefix, MAX_NAME), ".lkt");
        lkt_t *lkt = lkt_init(12);
        strcat(strncpy(fn0, prefix, MAX_NAME), ".pac");
        lkt_build_lookuptable(fn0, lkt); 
        lkt_dump(fn_lkt, lkt);
        lkt_destroy(lkt); 
        fprintf(stderr, "generate lookup table %.2f sec\n", (float)(realtime() - t));
    }

    {   
        bwt_t *bwt;
        t = realtime();
        fprintf(stderr, "[va_index]: generate SA and ISA from BWT\n");
        strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".bwt");
        bwt = bwt_restore_bwt(fn0);
        fprintf(stderr, "[va_index]: call SA\n");
        bwt_cal_sa(bwt, SA_INTV);
        fprintf(stderr, "[va_index]: call ISA\n");
        bwt_cal_isa(bwt);
        strcat(strncpy(fn0, prefix, MAX_NAME), ".sa");
        fprintf(stderr, "[va_index]: dump SA\n");
        bwt_dump_sa(fn0, bwt);
        strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".isa");
        bwt_dump_isa(fn0, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "cal SA and ISA %.2f sec\n", (float)(realtime() - t));
    }

}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_snpaware_idx
 *  Description:  build snpaware index for 20bp snpaware seed
 *       @fn_av:  file name of variants(UCSC format)
 *      @prefix:  file name prefix of index
 * =====================================================================================
 */
void build_snpaware_idx(const char *fn_av, const char *prefix)
{
    double t;
    char fn0[MAX_NAME];
    
    strcat(strncpy(fn0, prefix, MAX_NAME), ".bwt");
    bwt_t *bwt = bwt_restore_bwt(fn0);
    strcat(strncpy(fn0, prefix, MAX_NAME), ".sa");
    bwt_restore_sa(fn0, bwt);
    strcat(strncpy(fn0, prefix, MAX_NAME), ".isa");
    bwt_restore_isa(fn0, bwt);
 
    fprintf(stderr, "[va_index]: build variant map\n");
    t = realtime();        
    bwtint_t *isa, i;        
    isa= bwt->isa;
    all_var_t *av = variants_build(fn_av);    
    var_t *v = av->var;
    bp_t *is_variant = bp_init(bwt->seq_len+1);  

    uint32_t vi;
    for(vi = 0; vi < av->n_var; ++vi){
        uint32_t pos = v[vi].end;// 
        bp_set1(is_variant, isa[pos]);
        LOG(stderr, "\nvariant pos %u-> bwt index %u\n", pos, isa[pos]);
    }
    bp_gen_occ(is_variant);        
    if(bp_rank(is_variant, bwt->seq_len) != av->n_var){
        fprintf(stderr, "[%s ]:variant num wrong %u != %u ", __func__, bp_rank(is_variant, bwt->seq_len), av->n_var); 
        //exit(1);
    }
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".vmap"); 
    bp_dump(is_variant, fn0);
    fprintf(stderr, "build variants map %.2f sec\n", (float)(realtime() - t));
    
    t = realtime(); 
    fprintf(stderr, "[va_index]: generate variant seqs...");
    
    bwtint_t *sa; bntseq_t *bntseq; uint8_t *pac;
    uint32_t *varseq;        
    sa = bwt->sa;  
    bntseq = bns_restore(prefix);         
    varseq = (uint32_t *)calloc(av->n_var*2, sizeof(uint32_t));        
    //fill variant seqs
    pac = (uint8_t *)calloc(bntseq->l_pac/4+1, sizeof(uint8_t));
    fread(pac, bntseq->l_pac/4+1, sizeof(uint8_t), bntseq->fp_pac);
   
    for(i = 0; i <= bwt->seq_len; ++i){
        if(bp_get(is_variant, i) == 0) continue;
        vi = bp_rank(is_variant, i)-1;
        uint32_t s, e;
        e = sa[i];
        //fprintf(stderr, "[debug]: variant pos %u\n", e);
        assert(e >= 16);
        s = e-16;
        varseq[vi*2] = bns_extract_seq16(pac, s, e);            
        varseq[vi*2+1] = i;        }
    
    uint32_t preseq, k;
    //uint32_t *preseq = (uint32_t *)calloc(av->n_var, sizeof(uint32_t)); 
    uint8_t *origin = (uint8_t *)calloc(av->n_var, sizeof(uint32_t)); 
    uint32_t *sai = (uint32_t *)calloc(av->n_var, sizeof(uint32_t)); 
    uint32_t *preseq_cumulate = (uint32_t *)calloc(HASH_SIZE, sizeof(uint32_t));
    for(i =0; i < av->n_var; ++i){
        preseq = (varseq[i*2]>>4)&MASK_6NT;
        ++preseq_cumulate[preseq]; 
    }
    uint32_t tot, last_tot;
    for(i =0, tot=0; i < HASH_SIZE; ++i){
        last_tot = tot;
        tot += preseq_cumulate[i];
        preseq_cumulate[i] = last_tot;     
    }        
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".vseq");
    FILE *fp_varseq = xopen(fn0, "w");
    fwrite(&HASH_SIZE, sizeof(int), 1, fp_varseq);         
    fwrite(preseq_cumulate, sizeof(uint32_t), HASH_SIZE, fp_varseq);
    for(i =0; i < av->n_var; ++i){
        preseq = (varseq[i*2]>>4)&MASK_6NT;//trim 2 trail nt
        k = preseq_cumulate[preseq]++;
        origin[k] =  varseq[i*2]&MASK_2NT;
        sai[k]  = varseq[i*2+1];
    }

    fwrite(&av->n_var, sizeof(uint32_t), 1, fp_varseq);         
    fwrite(sai, sizeof(uint32_t), av->n_var, fp_varseq); 
    fwrite(origin, sizeof(uint8_t), av->n_var, fp_varseq); 
    
    err_fclose(fp_varseq);
    free(preseq_cumulate);
    free(sai);free(origin);
    free(varseq);
    fprintf(stderr, "build variants seq %.2f sec\n", (float)(realtime() - t) );
    bp_destroy(is_variant); 
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".var");
    variants_dump(av, fn0); 
    variants_destroy(av);
    free(pac);
    bns_destroy(bntseq); 
    bwt_destroy(bwt);


}
static inline void accumulate_cnt(int n, int *cnt_2nt)
{
    int i;
    for(i = 1; i < n-1; ++i) cnt_2nt[i] += cnt_2nt[i-1]; 
    for(i = n-1; i >0; --i) cnt_2nt[i] = cnt_2nt[i-1];
    cnt_2nt[0] = 0;
 
}
void build_local_bwt(ext_t *lext, uint32_t l, uint32_t r, uint32_t *sa0, uint32_t *sa1, uint8_t *bwt[])
{
    int i, j;
    int cnt_2nt[17] = {};
    uint32_t *last_sa, *cur_sa;
    lext += l;
    for(j = 0; j < r-l; ++j){
        uint32_t x = lext[j].seq&0xF;
        __set_bwt(bwt[0], j*2+1, x&3);
        __set_bwt(bwt[0], j*2, (x>>2)&3);
        ++cnt_2nt[x];
    }
    for(j =0; j < r-l; ++j) sa0[j] = j;
    last_sa = sa0; cur_sa = sa1;
    for(i = 2; i < 16; i+=2){
        accumulate_cnt(17, cnt_2nt);
        //sort cur_sa 
        for(j = 0; j < r-l; ++j){
            uint32_t x = lext[last_sa[j]].seq&0xF;
            cur_sa[cnt_2nt[x]++]  = last_sa[j];
        }
        //generate bwt from cur_sa        
        memset(cnt_2nt, 4*17, 0);
        for(j = 0; j < r-l; ++j){
            uint32_t x = lext[cur_sa[j]].seq&0xF;
            __set_bwt(bwt[i], j*2+1, x&3);
            __set_bwt(bwt[i], j*2, (x>>2)&3);
            ++cnt_2nt[x];
        }

        SWAP(uint32_t *, last_sa, cur_sa); 
    }
}
void build_occ(uint8_t *bwt[], int n_lseqs)
{
    if(n_lseqs <=128){
        return; 
    } else if( n_lseqs <=256 ){
        
    } else if(n_lseqs <= 256*8){
    
    } else if(n_lseqs <= 256*256) {
    
    } else {
    
    }
  
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_seq16_bwt
 *  Description:  
 * =====================================================================================
 */
void build_seq16_bwt (const char*prefix)
{
    int64_t i, j;
    /* reload extend index from disk */
    idx_t idx;
    restore_extend_idx(&idx, prefix); 
    for(i = 0; i < idx.n_tot; ++i){
        int n_lseqs = idx.lext_idx[i+1] - idx.lext_idx[i];
        uint8_t *bwt[8];
        for(j = 0; j < 8; ++j) bwt[j] = (uint8_t *)calloc((n_lseqs+1)/2, sizeof(uint8_t));  
        uint32_t *sa0 = (uint32_t *)calloc(n_lseqs, sizeof(uint32_t));
        uint32_t *sa1 = (uint32_t *)calloc(n_lseqs, sizeof(uint32_t));
        build_local_bwt(idx.lext0, idx.lext_idx[i], idx.lext_idx[i+1], sa0, sa1, bwt); 
        build_occ(bwt, n_lseqs); 
        
        
       //dump bwt
        free(sa0);free(sa1);
    
        for(j = 0; j<8;++j) free(bwt[j]);  
        
    } 
    return;
}		/* -----  end of function build_seq16_bwt  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_bwt_is_repeat
 *  Description:  
 * =====================================================================================
 */
#define N_EXT 5
#define LEN_NEXT_IDX 1024
#define LEN_SEED 20
#define LEN_EXT 16
#define IDX_MID_SIZE (254*256)
#define IS_SMLSIZ 16
void create_jmpmod(const char*prefix, bwt_t *bwt, int ext_i, uint32_t *jump_idx, uint8_t *mod_idx, uint32_t *cap_pos) 
{

    char fn[1024] = {};
    uint32_t seq_len = bwt->seq_len;
    sprintf(fn, "%s_file_next_%u.txt", prefix, ext_i);
    FILE *fp_next = fopen64(fn, "r");
    fseek(fp_next, 0, SEEK_END);
    uint32_t n_next_idx = ftell(fp_next)/(sizeof(uint32_t)*3);  
    fseek(fp_next, 0, SEEK_SET); 
    uint32_t (*next_idx)[3] = calloc(n_next_idx, 3*sizeof(uint32_t));
    fread(next_idx, 3*sizeof(uint32_t), n_next_idx, fp_next); 
    jump_idx = calloc((seq_len+1+255)/256 ,sizeof(uint32_t));//+1 for $ 
    mod_idx = calloc(n_next_idx, sizeof(uint8_t)); 
    cap_pos = calloc(n_next_idx, sizeof(uint32_t)); 
    uint32_t row_jump = 0;
    uint32_t last_idx_q = 0;
    //uint32_t jump_idx[LEN_NEXT_IDX+1]; 
    uint32_t row, i;
    for(row=0; row<n_next_idx; ++row){
        cap_pos[row] = next_idx[row][2]; 
        uint32_t idx = next_idx[row][0];
        uint32_t idx_q = idx/256;
        mod_idx[row] = idx%256;
        if(idx_q == last_idx_q){
            continue;
        } else{
            for(i= last_idx_q+1; i <= idx_q; ++i){
                jump_idx[i] = row; 
            }
            last_idx_q = idx_q;
        }
    
    } 
    jump_idx[(n_next_idx+255)/256] = n_next_idx; 
    //free(jump_idx); free(mod_idx);free(next_idx); 
    fclose(fp_next);


}
#define ROOT_PATH "../Index"
//define ROOT_PATH "Index"
void build_struct_idx(uint8_t *bwt_is_repeat, bwt_t *bwt, int ext_i, const char*prefix)
{
    int i, j;
    char fn[1024]={};
    FILE *fp_next;
/*
    FILE *fp_pos[3], *fp_idx[3], *fp_next;
  
    for(i=0; i<3; ++i){
        sprintf(fn, "%s_file_pos_%u_%u.txt", prefix, ext_i, i);
        fp_pos[i] = fopen64(fn, "wb"); 
        sprintf(fn, "%s_file_idx_%u_%u.txt", prefix, ext_i, i);
        fp_idx[i] = fopen64(fn, "wb");
    }
*/
    //sprintf(fn, "%s_file_next_%u.txt", prefix, ext_i);
    //sprintf(fn, "%s/data/seedidx/seedidx_%u.txt", ROOT_PATH, ext_i);
    sprintf(fn, "%s/seedidx_%u.txt", ROOT_PATH,ext_i);

    fp_next = fopen64(fn, "wb");
    if(fp_next == NULL){
        printf("open file %s fail\n", fn);
        printf("%s\n", strerror(errno));
        exit(1);
    }

    const uint32_t *sa = bwt->sa;
    int stt_val;
    uint8_t pre_flag_idx = 0;
    uint32_t glb_idx = 0;
    uint32_t num_idx =0;
    uint32_t bg_idx, ed_idx, sum_idx, bg_pos; 

    uint32_t l_ext_seed = LEN_SEED+2*LEN_EXT*ext_i;
    uint32_t len_idx[3] = {0,0,0};
    uint32_t next_idx[3] = {0,0,0};
    
    while(glb_idx <= bwt->seq_len+1){
//fprintf(stderr, "glb_idx = %u\n", glb_idx);        
        stt_val = 0;
        uint8_t cur_flag_idx = bwt_is_repeat[glb_idx];  
        if(pre_flag_idx ==0){
            if(cur_flag_idx == 2){
                bg_idx = glb_idx;
                num_idx = 1;
                //pre_flag_idx = 2; 
            } else if(cur_flag_idx == 1){
fprintf(stderr, "bg_idx = %u, ed_idx = %u\n", bg_idx, ed_idx);
fprintf(stderr, "glb_idx = %u, flag = %u %u\n", glb_idx, bwt_is_repeat[glb_idx], bwt_is_repeat[glb_idx-1]);

                printf("699, %s, should't be here!\n", __func__);
                exit(1);
                        
            
            }


        } else if(pre_flag_idx ==1){
            if(cur_flag_idx == 0){
                ed_idx = glb_idx-1;
                sum_idx = num_idx;
                stt_val = 1;
            } else if(cur_flag_idx==1){
                ++num_idx;
            } else if(cur_flag_idx ==2){
                sum_idx = num_idx;
                stt_val = 1;
                ed_idx = glb_idx-1;
                //bg_idx = glb_idx;
                num_idx = 1;
                //pre_flag_idx = 2;
            }


        } else if(pre_flag_idx ==2){
            if(cur_flag_idx ==0){
fprintf(stderr, "bg_idx = %u, ed_idx = %u\n", bg_idx, ed_idx);
fprintf(stderr, "glb_idx = %u, flag = %u %u\n", glb_idx, bwt_is_repeat[glb_idx], bwt_is_repeat[glb_idx-1]);

                printf("728, %s, should't be here!\n", __func__);
                exit(1);
            } else if(cur_flag_idx ==1){
                ++num_idx;
            } else if(cur_flag_idx == 2){
fprintf(stderr, "bg_idx = %u, ed_idx = %u\n", bg_idx, ed_idx);
fprintf(stderr, "glb_idx = %u, flag = %u %u\n", glb_idx, bwt_is_repeat[glb_idx], bwt_is_repeat[glb_idx-1]);

                printf("736, %s, should't be here!\n", __func__);
                exit(1);
            }
        

        } 
        pre_flag_idx = cur_flag_idx;
        ++glb_idx;
        

        if(stt_val == 0) continue;


        bg_idx = ed_idx+1-sum_idx;
        uint32_t idx = bg_idx;
        next_idx[0] = bg_idx;
        next_idx[1] = sum_idx;
        //if(sum_idx >IS_SMLSIZ) 
        fwrite(next_idx, sizeof(uint32_t), 2, fp_next);
    }
//+++++++++++++++++++++++++++++       
//以下代码是最后一个三元组数据生成的代码
/*  
        if(ext_i==NUM_EXT){
            if(sum_idx < 255){
                next_idx[2] = len_idx[0];
                len_idx[0] += sum_idx;
            } else{
                next_idx[2] = len_idx[0];
                len_idx[0] += ++sum_idx;
            } 
            fwrite(next_idx, sizeof(uint32_t), 3, fp_next);
            continue;
        }
//++++++++++++++++++++++++++++++
        if(sum_idx <= IS_SMLSIZ){
            next_idx[2] = len_idx[0];
            len_idx[0] += sum_idx;
        } else{
            next_idx[2] = len_idx[1]++;
        }            
        fwrite(next_idx, sizeof(uint32_t), 3, fp_next);
    }
*/
    
    
/*  
    for(i=0; i<3; ++i){
        fclose(fp_pos[i]);
        fclose(fp_idx[i]);
    }
*/
    fclose(fp_next); 
    return;
}
//+++++++++++++++++++++++++++++++
//以下代码已被上一段代码取代

/*  
        if(sum_idx <= IS_SMLSIZ){
            //next_idx[2] = len_idx[0]++;
            next_idx[2] = len_idx[0];
            len_idx[0] += sum_idx;
            fwrite(sa+bg_idx, sizeof(uint32_t), sum_idx, fp_pos[0]); 
            fwrite(&idx, sizeof(uint32_t), 1, fp_idx[0]); 
        } else if(sum_idx <= IDX_MID_SIZE){
            next_idx[2] = len_idx[1]++;
            fwrite(sa+bg_idx, sizeof(uint32_t), sum_idx, fp_pos[1]); 
            fwrite(&idx, sizeof(uint32_t), 1, fp_idx[1]); 
        } else{//>IDX_MID_SIZE
            next_idx[2] = len_idx[2]++;
            fwrite(sa+bg_idx, sizeof(uint32_t), sum_idx, fp_pos[2]); 
            fwrite(&idx, sizeof(uint32_t), 1, fp_idx[2]); 
        }

        printf("[next_]:bg_idx,ed_idx, sum_idx=(%u,%u, %u)\n",bg_idx, ed_idx, sum_idx);
//+++++++++++++++++++++++++++++++++++++++++++++++++
         
        printf("for(i = bg_idx; i <= bg_idx+sum_idx-1; ++i), begin\n");
        for(i = bg_idx; i <= bg_idx+sum_idx-1; ++i){
            uint32_t pos = bwt_sa(bwt, i);
            if(pos+l_ext_seed > bwt->seq_len){
printf("if(pos+l_ext_seed > bwt->seq_len), pos = %u, l_ext_seed=%u\n", pos, l_ext_seed); 
            
            } 
        
        }
//+++++++++++++++++++++++++++++++++++++++++++++++++
        
        fwrite(next_idx, sizeof(uint32_t), 3, fp_next);
    }

 
    for(i=0; i<3; ++i){
        fclose(fp_pos[i]);
        fclose(fp_idx[i]);
    }
    fclose(fp_next); 

    
}
*/
void build_bwt_is_repeat_old( const char *prefix, const bwt_t *bwt, const uint8_t *pac, uint8_t* ref_is_repeat, uint8_t *bwt_is_repeat[])
{
    uint8_t *bwt_repeat0, *bwt_repeat1; 
    bwt_repeat0 = bwt_is_repeat[0];
    bwt_repeat1 = bwt_is_repeat[1];
    int64_t i, j, ext_i;

    for(i =0; i < bwt->seq_len; ++i) {
        if(ref_is_repeat[i] == 2) ref_is_repeat[i] = 0;
    }
    for(ext_i =1; ext_i <= N_EXT; ++ext_i){
        int stt_val;
        uint8_t pre_flag_idx = 0;
        uint32_t glb_idx = 0;
        uint32_t num_idx =0;
        uint32_t bg_idx, ed_idx, sum_idx, bg_pos; 
        //while(glb_idx < bwt->seq_len ){
        uint32_t l_ext_seed = LEN_SEED+2*LEN_EXT*ext_i;
        while(glb_idx <= bwt->seq_len ){
            stt_val = 0;
	        uint8_t cur_flag_idx = bwt_repeat0[glb_idx];  
            if(pre_flag_idx ==0){
                if(cur_flag_idx == 2){
                    bg_idx = glb_idx;
                    num_idx = 1;
                    //pre_flag_idx = 2; 
                }

            } else if(pre_flag_idx ==1){
                if(cur_flag_idx == 0){
                    ed_idx = glb_idx-1;
                    sum_idx = num_idx;
                    stt_val = 1;
                } else if(cur_flag_idx==1){
                    ++num_idx;
                } else if(cur_flag_idx ==2){
                    sum_idx = num_idx;
                    stt_val = 1;
                    ed_idx = glb_idx-1;
                    //bg_idx = glb_idx;
                    num_idx = 1;
                    //pre_flag_idx = 2;
                }
            } else if(pre_flag_idx ==2){
                if(cur_flag_idx ==0){
                    fprintf(stderr, "should't be here!\n");
                    exit(1);
                } else if(cur_flag_idx ==1){
                    ++num_idx;
                } else if(cur_flag_idx == 2){
                    fprintf(stderr, "should't be here!\n");
                    exit(1);
                }
            
            } 
            pre_flag_idx = cur_flag_idx;
            ++glb_idx;
            if(stt_val == 0) continue;
        
            //bg_idx = ed_idx-sum_idx;
            bg_idx = ed_idx+1-sum_idx;
 
 
            fprintf(stderr, "bg_idx, ed_idx=%u, %u\n", bg_idx, ed_idx);
            uint32_t cur_idx = bg_idx;        
            while(cur_idx < ed_idx){
                uint32_t cur_pos = bwt_sa(bwt, cur_idx)+l_ext_seed-1;            
                if(cur_pos >= bwt->seq_len){
                    ++cur_idx;
                    continue;
                }
                for(i=0; i< 2*LEN_EXT; ++i){
                    //if(ref_is_repeat[cur_pos-i]>0) break;
                    if(ref_is_repeat[cur_pos-i]==1) break;
                }
                if(i < 2*LEN_EXT){
                    ref_is_repeat[cur_pos] =1;
                    cur_idx++;
                    continue;
                }
                //bg_pos = cur_pos-LEN_SEED-l_ext_seed;//fix?
                bg_pos = cur_pos-l_ext_seed+1;//fix?
                //l_seq = LEN_SEED+l_ext_seed;//fix?
                int l_seq = l_ext_seed;//fix?
                uint8_t seed[256] = {};
                for(i = 0; i < l_seq; ++i){ seed[i] = __get_pac(pac, bg_pos+i);}
                uint32_t k =0, l=bwt->seq_len;
                bwt_match_exact_alt(bwt, l_seq, seed, &k, &l);
                if(k < l) {
           	/*
		         if(k != cur_idx){
                        fprintf(stderr, "k != cur_idx!, %u, %u\t, %u\n", k, l, cur_idx);
                    }
                    if(l != ed_idx){
                        fprintf(stderr, "l != ed_idx!, %u, %u\n", l, ed_idx);
                    }
		*/
			        bwt_repeat1[k] =2;
                    //j = cur_idx+l;
                    for(i=k+1; i <=l;++i){
                        bwt_repeat1[i] =1;
                    }
                    for(i = k; i <=l; ++i){
                        cur_pos = bwt_sa(bwt, i)+l_ext_seed-1;
                        //ref_is_repeat[i] = 2;
                        ref_is_repeat[cur_pos] = 2; 
                    }
                    //cur_idx = cur_idx+l-k+1;
                    cur_idx = l+1;
                } else{
                    ref_is_repeat[cur_pos] = 1;
                    ++cur_idx; 
                }
            }//end while(cur_idx < ed_idx)
        }// end while(idx < LEN_REF_SEQ)
        fprintf(stderr, "\nseed len = %u\n", l_ext_seed);
        //
        //for(i=0; i< bwt->seq_len; ++i){ fprintf(stderr, "%u\t%u\n", i, bwt_repeat1[i]);}
        for(i =0; i < bwt->seq_len; ++i) {
            if(ref_is_repeat[i] == 2) ref_is_repeat[i] = 0;
        }
        /* for debug */
/*  */        
        fprintf(stderr, "for debug\n");
        i = 0;
        while(i<= bwt->seq_len){
	        fprintf(stderr, "idx = %u\n", i);
            uint32_t pos = bwt_sa(bwt, i);  
            if(pos > bwt->seq_len-l_ext_seed) {
                ++i;
                continue;
            }
            int l_seq = l_ext_seed;//fix?
            uint8_t seed[256] = {};
            for(j = 0; j < l_seq; ++j){ seed[j] = __get_pac(pac, pos+j);}
            uint32_t k =0, l=bwt->seq_len;
            bwt_match_exact_alt(bwt, l_seq, seed, &k, &l);
            if(k > l) {
                fprintf(stderr, "should't be here!\n");
                return;
            } else if(k == l){
                if(bwt_repeat1[k] != 0) fprintf(stderr, "[Error]:flag =0, seed_len= %u\t repeat[%u]=%u\n", l_ext_seed, k, bwt_repeat1[k]);
               
            } else{

                if(bwt_repeat1[k] != 2) fprintf(stderr, "[Error]:flag =2, seed_len= %u\t repeat[%u]=%u\n", l_ext_seed, k, bwt_repeat1[k]);            
                for(j = k+1; j <=l; ++j) {
                    if(bwt_repeat1[j] != 1) fprintf(stderr, "[Error]:flag =1, seed_len= %u\t repeat[%u]=%u\n", l_ext_seed, j, bwt_repeat1[j]);
                }
              
            }  
            i += l-k+1;
        }
        //for(i=0; i<=bwt->seq_len; ++i) printf("[%u %u]:\t%u\n", l_ext_seed, i, bwt_repeat1[i]); 
        build_struct_idx(bwt_repeat1, bwt, ext_i-1, prefix); 
        SWAP(uint8_t *, bwt_repeat0, bwt_repeat1);
        memset(bwt_repeat1, 0, bwt->seq_len+1);
    }
    
    
    
    return;
}		/* -----  end of function build_bwt_is_repeat  ----- */

/*
//need more test
int cut_ext_seq2(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, uint8_t *bwt_repeat)
{
    bwtint_t i, bg_idx, ed_idx;
    uint64_t last_seq32, cur_seq32;
    bg_idx= k;
    while(bg_idx <=l){
        //caculate [bg_idx, ed_idx)
        last_seq32 = fetch_seq32(bwt, bg_idx);
        ed_idx = bg_idx+1;
       
        while(ed_idx <= l){
            uint32_t pos = bwt(bwt, ed_idx);
             
            if(pos+32>bwt->seq_len  || cur_seq32 != last_seq32) {
                break;
            }
            ++ed_idx;
        }
        //set bwt_repeat[bg_idx:ed_idx] 
        if(ed_idx - bg_idx == 1){
            bwt_repeat[bg_idx] = 0;
        
        } else{
            bwt_repeat[bg_idx] = 2;
            for(i = bg_idx+1; i < ed_idx; ++i) bwt_repeat[i] = 1;
        }
    }
    return 0;

}
*/

int cut_ext_seq0(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, uint32_t sort_flag[])
{
    
//fprintf(stderr, "k = %u, l = %u\n", k, l);
    bwtint_t i, pos, st, ed;
    uint32_t *ItoPos = bwt->sa;
    /* fetch 16 mer preseq in [k, l) and rm duplicate preseq */
    uint32_t max_pos, min_pos;

    if(offset<0){
        printf("%s\n","Err: offset can not be less than 0 !");
        exit(1);
    }
    max_pos = bwt->seq_len-offset-8;        
    //max_pos = bwt->seq_len-20;        
    uint32_t c;
    uint32_t ii;
    uint32_t r = 0, r0;
    uint8_t ch, jj, kk;

    uint32_t old_seq16, seq16;
    uint32_t sum = 0 ;
    int j = 0;
    int old_flag = 0;
    int cur_flag = 0;
    sort_flag[0] = sum; 
    pos = bwt_sa(bwt, k);       
    if(pos <= max_pos){
        st = pos + offset; 
        r =0;
        ch = pac[st/4];
        jj = (st%4)*2;
        kk = (8-jj); 
        r = (uint32_t )((0xFF>>jj)&ch);
        //for(ii = st/4+1; ii <st/4+4 ; ++ii){
        for(ii = st/4+1; ii <st/4+2 ; ++ii){
            r<<=8;
            r0 = (uint32_t)pac[ii];  
            r|= r0;
        }
        if(st%4>0){
            r <<=jj;
            r0 = (uint32_t)pac[ii];
            r0>>=kk;
            r|=r0;
        }
        seq16 = r;

        old_seq16 = seq16;
        old_flag = 0;
    }else{
        //fprintf(stderr, "pos = %u, K = %u, sum = %u\n", pos, k, sum);
        old_seq16 = 0; 
        old_flag = 1;
//fprintf(stderr, "1060, pos = %u\n", pos);

    }
//fprintf(stderr, "seq16 = %u, pos = %u\n", seq16, pos);
    for(i = k+1; i <= l; ++i){
        j = i-k;
        pos = bwt_sa(bwt, i);       
        if(pos <= max_pos){
            st = pos + offset; 
            r =0;
            ch = pac[st/4];
            jj = (st%4)*2;
            kk = (8-jj); 
            r = (uint32_t )((0xFF>>jj)&ch);
            for(ii = st/4+1; ii <st/4+2 ; ++ii){
                r<<=8;
                r0 = (uint32_t)pac[ii];  
                r|= r0;
            }
            if(st%4>0){
                r <<=jj;
                r0 = (uint32_t)pac[ii];
                r0>>=kk;
                r|=r0;
            }
            seq16 = r;
            cur_flag = 0;
        }else{
            ++sum;
            cur_flag = 1;
        }
//fprintf(stderr, "pos = %u, seq16 = %u, flag = %u\n", pos, seq16, cur_flag);
        if(cur_flag == 0){
            if(old_flag == 1 || old_seq16 < seq16){
                ++sum;
            } 
        }
        old_flag = cur_flag;
        old_seq16 = seq16; 
        sort_flag[j] = sum;
//fprintf(stderr, "pos = %u, seq16 = %u, flag = %u, sort_flag[%u]=%u\n", pos, seq16, cur_flag, j, sort_flag[j]);

//fprintf(stderr, "1104, sum = %u, pos = %u\n", sum, pos);
    }
    sort_flag[++j] = ++sum;
    sort_flag[++j] = ++sum;
    return j;
}

//uint32_t *sort_flag = calloc(MAX_COUNT, sizeof(uint32_t));

int cut_ext_seq1(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, uint32_t sort_flag[])
{
    bwtint_t i, pos, st, ed;
    uint32_t *ItoPos = bwt->sa;
    /* fetch 16 mer preseq in [k, l) and rm duplicate preseq */
    uint32_t max_pos, min_pos;

    if(offset<0){
    	printf("%s\n","Err: offset can not be less than 0 !");
    	exit(1);
    }
    max_pos = bwt->seq_len-offset-16;        
    uint32_t c;
    uint32_t ii;
    uint32_t r = 0, r0;
    uint8_t ch, jj, kk;

    uint32_t old_seq16, seq16;
    uint32_t sum = 0 ;
    int j = 0;
    int old_flag = 0;
    int cur_flag = 0;
    sort_flag[0] = sum; 
    pos = bwt_sa(bwt, k);       
    if(pos <= max_pos){
        st = pos + offset; 
        r =0;
        ch = pac[st/4];
        jj = (st%4)*2;
        kk = (8-jj); 
        r = (uint32_t )((0xFF>>jj)&ch);
        for(ii = st/4+1; ii <st/4+4 ; ++ii){
            r<<=8;
            r0 = (uint32_t)pac[ii];  
            r|= r0;
        }
        if(st%4>0){
            r <<=jj;
            r0 = (uint32_t)pac[ii];
            r0>>=kk;
            r|=r0;
        }
        seq16 = r;

        old_seq16 = seq16;
        old_flag = 0;

    }else{

        //fprintf(stderr, "pos = %u, K = %u, sum = %u\n", pos, k, sum);
        old_seq16 = 0; 
        old_flag = 1;

    }

 
    for(i = k+1; i <= l; ++i){
        j = i-k;
        pos = bwt_sa(bwt, i);       
		if(pos <= max_pos){
			st = pos + offset; 
            r =0;
            ch = pac[st/4];
            jj = (st%4)*2;
            kk = (8-jj); 
            r = (uint32_t )((0xFF>>jj)&ch);
            for(ii = st/4+1; ii <st/4+4 ; ++ii){
                r<<=8;
                r0 = (uint32_t)pac[ii];  
                r|= r0;
            }
            if(st%4>0){
                r <<=jj;
                r0 = (uint32_t)pac[ii];
                r0>>=kk;
                r|=r0;
            }
            seq16 = r;
            cur_flag = 0;
        }else{
            cur_flag = 1;
            ++sum;
        } 
        if(cur_flag == 0) {
            if(old_flag == 1 || old_seq16 < seq16) {
                sum++;
            }
        }
        old_flag = cur_flag;
        old_seq16 = seq16; 
        sort_flag[j] += sum;
    }
    sort_flag[++j] += ++sum;
    sort_flag[++j] += ++sum;
    return j;
}

//=============================================================================================================
int get_12mer_correct(uint32_t hash_boundry[], uint32_t end_kmer)
{
    int i;
    //uint32_t top = hash_boundry[0], bot = hash_boundry[LEN_SEED-1];
    uint32_t top = 0, bot = 11;
    uint32_t mid;
    int flag = 0;
    while(top < bot) {
        mid = (top+bot)/2;
        if(end_kmer == hash_boundry[mid]){
            flag = 1;
            break;
        } else if(end_kmer < hash_boundry[mid]){
            bot = mid;
        } else{
            top = mid+1;
        }
    }
    return flag;

}
void build_bwt_is_repeat ( const char *prefix, const idx_t *idx, const uint8_t *pac)
{
    clock_t t = clock();

    int64_t i, j, ext_i;
    bwt_t *bwt = idx->bwt;
    lkt_t *fastmap = idx->fastmap; 
 
printf("bwt->seq_len = %u\n", bwt->seq_len+2);
    uint8_t *bwt_is_repeat[2];
    bwt_is_repeat[0] = (uint8_t *)calloc(bwt->seq_len+2, sizeof(uint8_t));/* indicate repeat seed on bwt index */
    
    uint32_t hash_boundry[12] = {};
    uint8_t seed[12] = {};
    for(i =0; i < 12; ++i){
        seed[i] = __get_pac(pac, bwt->seq_len-12+i);

    }
    uint32_t k = 0, l = 0;
    for(i =0; i < 11; ++i){
        int n = bwt_match_exact_alt(bwt, 1, seed+12-1-i, &k, &l);
        hash_boundry[i] = k;
    }
    for(i = 0; i < 11; ++i){
        for(j = i+1; j < 11; ++j){
            if(hash_boundry[i] > hash_boundry[j]){
                uint32_t tmp;
                tmp = hash_boundry[i];
                hash_boundry[i] = hash_boundry[j];
                hash_boundry[j] = tmp;
            }
        }
    }

    uint32_t MAX_12_count = fastmap->item[0];
    for(i = 0; i < idx->fastmap->n_item+1; ++i){
        uint32_t count = fastmap->item[i+1] - fastmap->item[i];
        if(count > MAX_12_count) MAX_12_count = count;
    }
 
printf("sizeof(sort_flag_12) = %u\n", MAX_12_count);
printf("time elapsed time :%u sec\n", clock()/CLOCKS_PER_SEC);
    uint32_t stat_8[2] = {};
    uint32_t *sort_flag = calloc(MAX_12_count+3, sizeof(uint32_t));
   
    bwt_is_repeat[1] = (uint8_t *)calloc(bwt->seq_len+2, sizeof(uint8_t));/* indicate repeat seed on bwt index */
    uint8_t *bwt_repeat0, *bwt_repeat1; 
    bwt_repeat0 = bwt_is_repeat[0];
    bwt_repeat1 = bwt_is_repeat[1];


//+++++++++++++++++++++++++++++++++++++++++++++++++++
    uint32_t bg_idx, ed_idx, sum_idx, bg_pos; 
    uint32_t glb_idx = 0;
    bg_idx = ed_idx= glb_idx;
    //for(glb_idx = 0; glb_idx <= fastmap->n_item; ++glb_idx){
    //for(glb_idx =0; glb_idx < fastmap->n_item; ++glb_idx){
    uint32_t MAX_20_count = 0;
    for(glb_idx =0; glb_idx < fastmap->n_item; ++glb_idx){
        
        bg_idx = fastmap->item[glb_idx];
        ed_idx = fastmap->item[glb_idx+1];
//fprintf(stderr, "[glb_idx=%u]:bg_idx = %u, ed_idx =%u\n", glb_idx, bg_idx, ed_idx);
        ed_idx -= get_12mer_correct(hash_boundry, ed_idx-1);


        //ed_idx = (glb_idx == fastmap->n_item-1)?bwt->seq_len:fastmap->item[glb_idx+1];
//fprintf(stderr, "[glb_idx=%u]:bg_idx = %u, ed_idx =%u\n", glb_idx, bg_idx, ed_idx);
        if(bg_idx >= bwt->seq_len) break;
        uint32_t num_idx;
        if(bg_idx +1>= ed_idx) {

                ++stat_8[0]; 
		continue;}
        else{ 
            num_idx = ed_idx - bg_idx;
            ed_idx -= 1;
        }
/*        
        for(i = 0; i <= num_idx; ++i){
            sort_flag[i] = 0; 
        }
*/
        int offset = 12;
        cut_ext_seq0(pac, bwt, bg_idx, ed_idx, offset, sort_flag);
        uint32_t i_st = 0; 
//fprintf(stderr, "num_idx = %u, bg_idx= %u, ed_idx = %u\n", num_idx, bg_idx, ed_idx);
        for(i = 0; i < num_idx; ++i){
//fprintf(stderr, "sort_flag[%u] = %u\n", i, sort_flag[i]);
            if(sort_flag[i] < sort_flag[i+1]){
//=====================================================
                if(i - i_st >= 1){
//fprintf(stderr, "[1248, glb_idx=%u]\n", glb_idx);
//fprintf(stderr, "[1249, idx=%u]\n", bg_idx+i_st);

//50521244
                    bwt_repeat0[bg_idx+i_st] = 2; 
//if(bg_idx+i_st == 844453254) fprintf(stderr, "[%s, %u]bg_idx = %u, num = %u\n", __func__,__LINE__, bg_idx+i_st, i-i_st);
                    for(j = i_st+1; j <= i; ++j) bwt_repeat0[bg_idx+j] = 1;
                    if(i-i_st > MAX_20_count) MAX_20_count = i-i_st;
                }
                i_st = i+1;                 
                ++stat_8[1]; 
            }      
        }
    }// end while(glb_idx < bwt->seq_len)

printf("sizeof(sort_flag_20) = %u\n", MAX_20_count);
    free(sort_flag);
//exit(1);
    sort_flag = calloc(MAX_20_count+3, sizeof(uint32_t));
/* 
    int flag = bwt_repeat0[0];
    for(glb_idx = 1; glb_idx < bwt->seq_len;++glb_idx){
        if(flag==2 && bwt_repeat0[glb_idx]!= 1){
            fprintf(stderr, "glb_idx = %u, bwt_repeat = %u, %u\n", glb_idx, bwt_repeat0[glb_idx-1], bwt_repeat0[glb_idx]);
        }
        flag = bwt_repeat0[glb_idx];
    }

*/
    //处理右边界问题
/*
    uint8_t last_seq[LEN_SEED];
    for(j = 0; j < LEN_SEED; ++j) {last_seq[j] = __get_pac(pac, bwt->seq_len-19+j);}
    k = 0; l = 0;
    for(i =1; i <=LEN_SEED; ++i){
        bwt_match_exact_alt(bwt, 1, last_seq+LEN_SEED-i, &k, &l);
        if(bwt_repeat0[k] >0) {
            //if(bwt_repeat0[0] == 2)
            bwt_repeat[k] = 0;
        } 
    }
*/
//++++++++++++++++++++++++++++++++++++++++++++++++
fprintf(stderr, "for debug 12\n");
fprintf(stderr, "stat[0] = %u, stat[1] =%u\n\n", stat_8[0], stat_8[1]);
    

//exit(1);
    
    int n_error[3] ={};
/*
    i = 0;
    while(i<= bwt->seq_len){
        //fprintf(stderr, "idx = %u\n", i);
        uint32_t pos = bwt_sa(bwt, i);  
        if(pos > bwt->seq_len-12) {
            ++i;
            continue;
        }
        int l_seq = 12;//fix?
        uint8_t seed[256] = {}; uint32_t seq12 = 0;
        for(j = 0; j < 12; ++j){ 
            seed[j] = __get_pac(pac, pos+j);
            seq12 <<= 2;
            seq12 |= seed[j];
        }
        uint32_t k0 = fastmap->item[seq12], l0 = fastmap->item[seq12+1]; 
        l0 = l0-1 -get_12mer_correct(hash_boundry, l0-1);
        uint32_t k =0, l=bwt->seq_len;
        bwt_match_exact_alt(bwt, l_seq, seed, &k, &l);
        if(k0 != k || l0 != l){
            fprintf(stderr, "[%u] 12mer(%u)error: k0, k = %u,%u   l0, l = %u, %u\n", i, seq12+1, k0, k, l0, l);

        }
        i += l-k+1;
    }
*/    
    for(i =0; i <12; ++i) fprintf(stderr, "hash = %u\n", hash_boundry[i]);

//exit(1);
//++++++++++++++++++++++++++++++++++++++++++++++
    fprintf(stderr, "for debug 20\n");
    
/*  
//exit(1);
    //int n_error[3] ={};
    //i = 2961267;
  
    i = 0;
    while(i<= bwt->seq_len){
        //fprintf(stderr, "idx = %u\n", i);
        uint32_t pos = bwt_sa(bwt, i);  
        if(pos > bwt->seq_len-20) {
            ++i;
            continue;
        }
        int l_seq = 20;//fix?
        uint8_t seed[256] = {};
        for(j = 0; j < l_seq; ++j){ seed[j] = __get_pac(pac, pos+j);}
        uint32_t k =0, l=bwt->seq_len;
    
        if(bwt_match_exact_alt(bwt, l_seq, seed, &k, &l)==0) {
            fprintf(stderr, "should't be here!\n");
            exit(1);
        } else if(k == l){
            if(bwt_repeat0[k] != 0) {
                fprintf(stderr, "[Error]:flag =0, seed_len= %u\t repeat[%u]=%u\n", 20, k, bwt_repeat0[k]);
exit(1);
++n_error[0]; 
            }

        } else{

            if(bwt_repeat0[k] != 2) {
                fprintf(stderr, "[Error]:k= %u, l= %u\n", k,l);            
++n_error[2];
                
                fprintf(stderr, "[Error]:flag =2, seed_len= %u num = %u\t repeat[%u]=%u\n", 20,l-k+1, k, bwt_repeat0[k]);            
            
exit(1);

            }
            for(j = k+1; j <=l; ++j) {

                if(bwt_repeat0[j] != 1){
                    fprintf(stderr, "[Error]:k= %u, l= %u\n", k,l);            
                    
                    fprintf(stderr, "[Error]:flag =1, seed_len= %u\t repeat[%u]=%u\n", 20, j, bwt_repeat0[j]);
exit(1);
++n_error[1];

                }

            }
          
        }  
        i += l-k+1;
    }
fprintf(stderr, "error = %u, %u, %u\n", n_error[0], n_error[1], n_error[2]);

*/

    //memset(bwt_repeat1, 0, bwt->seq_len+2);
//+++++++++++++++++++++++++++++++++++++++++++++++++++
    build_struct_idx(bwt_repeat0, bwt, 0, prefix); 
fprintf(stderr, "dump 20 bp seed data!\n");
    for(ext_i =1; ext_i <= N_EXT; ++ext_i){
        t = clock();
        uint8_t pre_flag_idx = 0;
        uint32_t num_idx =0;
        uint32_t l_ext_seed = LEN_SEED+2*LEN_EXT*ext_i;
        glb_idx = 0;
        bg_idx = ed_idx= glb_idx;
int n_interval = 0;

        while(glb_idx <= bwt->seq_len ){
++n_interval;
            
            while(glb_idx <=bwt->seq_len && bwt_repeat0[glb_idx] ==0){ glb_idx++;}
            bg_idx = glb_idx++; 
            while(glb_idx <=bwt->seq_len && bwt_repeat0[glb_idx] ==1){ glb_idx++;}
            ed_idx = glb_idx-1;
            
            num_idx = ed_idx+1 - bg_idx;
            for(i = 0; i <= num_idx; ++i){
                sort_flag[i] = 0; 
            }
            int offset = LEN_SEED+2*LEN_EXT*(ext_i-1);
            cut_ext_seq1(pac, bwt, bg_idx, ed_idx, offset, sort_flag);
            cut_ext_seq1(pac, bwt, bg_idx, ed_idx, offset+16, sort_flag);
            
            uint32_t i_st = 0; 
            //if(num_idx == 2) fprintf(stderr, "sort_flag = %u %u %u\n", sort_flag[0], sort_flag[1], sort_flag[2]);
            for(i = 0; i < num_idx; ++i){
                if(sort_flag[i] < sort_flag[i+1]){
//===================================================================
                    if(i - i_st < 1){
                        //ref_is_repeat[cur_pos] = 3;
                    } else{
                        bwt_repeat1[bg_idx+i_st] = 2; 
                        for(j = i_st+1; j <= i; ++j) bwt_repeat1[bg_idx+j] = 1;
                    }
                    i_st = i+1;
 
                }                 
            }

if(n_interval%1000000 == 0) printf("[%u]: glb_idx = %u, bg_idx = %u, ed_idx = %u\n",n_interval,  glb_idx, bg_idx, ed_idx);
        }// end while(glb_idx < bwt->seq_len)
        fprintf(stderr, "for debug\n");
/*  
        if(ext_i > 1) continue; 
        i = 0;
        while(i<= bwt->seq_len){
            //fprintf(stderr, "idx = %u\n", i);
            uint32_t pos = bwt_sa(bwt, i);  
            if(pos > bwt->seq_len-l_ext_seed) {
                ++i;
                continue;
            }
            int l_seq = l_ext_seed;//fix?
            uint8_t seed[256] = {};
            for(j = 0; j < l_seq; ++j){ seed[j] = __get_pac(pac, pos+j);}
            uint32_t k =0, l=bwt->seq_len;
        
            if(bwt_match_exact_alt(bwt, l_seq, seed, &k, &l)==0) {
                fprintf(stderr, "should't be here!\n");
                exit(1);
                return;
            } else if(k == l){
                if(bwt_repeat1[k] != 0) {
                    fprintf(stderr, "[Error]:flag =0, seed_len= %u\t repeat[%u]=%u\n", l_ext_seed, k, bwt_repeat1[k]);
                    exit(1);
                }

            } else{

                if(bwt_repeat1[k] != 2) {
                    fprintf(stderr, "[Error]:flag =2, seed_len= %u num = %u\t repeat[%u]=%u\n", l_ext_seed,l-k+1, k, bwt_repeat1[k]);            
                    exit(1);

                }
                for(j = k+1; j <=l; ++j) {
                    if(bwt_repeat1[j] != 1){
                        fprintf(stderr, "[Error]:flag =1, seed_len= %u\t repeat[%u]=%u\n", l_ext_seed, j, bwt_repeat1[j]);
                        exit(1);
                    }

                }
              
            }  
            i += l-k+1;
        }
*/
/*   
        //for(i=0; i<=bwt->seq_len; ++i) printf("[%u %u]:\t%u\n", l_ext_seed, i, bwt_repeat1[i]); 
        uint32_t __glb_idx=0;
        while(__glb_idx < bwt->seq_len){
            while(__glb_idx <=bwt->seq_len && bwt_repeat1[__glb_idx] ==0){ __glb_idx++;}
            bg_idx = __glb_idx++; 
            while(__glb_idx <=bwt->seq_len && bwt_repeat1[__glb_idx] ==1){ __glb_idx++;}
            ed_idx = __glb_idx-1;
            for(i = bg_idx; i <= ed_idx; ++i){
                uint32_t pos = bwt_sa(bwt, i);
                if(pos+l_ext_seed > bwt->seq_len){
//printf("if(pos+l_ext_seed > bwt->seq_len), pos = %u, l_ext_seed=%u\n", pos, l_ext_seed); 
                
                }            
            }
        } 
*/
        fprintf(stderr, "ext = %u, len = %u ++++++++++++++\n", ext_i, 20+ext_i*32);
        build_struct_idx(bwt_repeat1, bwt, ext_i, prefix); 
	fprintf(stderr, "dump data finish\n");
        SWAP(uint8_t *, bwt_repeat0, bwt_repeat1);
        memset(bwt_repeat1, 0, bwt->seq_len+2);
    
printf("ext(%u)time elapsed time :%u sec\n", ext_i, clock()/CLOCKS_PER_SEC);
    }
    
    
    free(sort_flag);
    free(bwt_is_repeat[1]);
    free(bwt_is_repeat[0]);
fprintf(stderr, "1706\n");
fprintf(stderr, "%u %s return\n", __LINE__, __func__);
    return;
}		/* -----  end of function build_bwt_is_repeat  ----- */
void build_extend_idx_alt(const char *prefix)
{
/*
    global_stat_36seed = (int *)calloc(MAX_COUNT, sizeof(int));
    global_stat_20seed0 = (int *)calloc(MAX_COUNT, sizeof(int));
    global_stat_20seed1 = (int *)calloc(MAX_COUNT, sizeof(int));
*/  
    char fn0[MAX_NAME];
    double t = realtime();
    
    /* reload FM-index and packed reference */
    fprintf(stderr, "reload FM-index..\n");
    int64_t i;
    bwt_t *bwt; bntseq_t *bntseq; uint8_t *pac;
    bwtint_t ed, st, k, l, pos, isa;

/*
    cat_filename(fn0, prefix, ".bwt"); 
    bwt = bwt_restore_bwt(fn0);   
    cat_filename(fn0, prefix, ".sa"); 
    bwt_restore_sa(fn0, bwt);
    bntseq = bns_restore(prefix);
    pac = (uint8_t *)calloc(bntseq->l_pac/4+1, sizeof(uint8_t));
    fread(pac, bntseq->l_pac/4+1, sizeof(uint8_t), bntseq->fp_pac);
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".lkt");
    lkt_t *fastmap= lkt_restore(fn0); 
*/
    idx_t *idx = idx_restore(prefix);
    bwt = idx->bwt;
    pac = idx->pac;
    
    
    fprintf(stderr, "%.2f sec\n", realtime()-t); 
    /* generate repeat seeds and bwt range */
    fprintf(stderr, "gen repeat on bwt\n");
    //bp_t *ref_is_repeat; 
    //uint8_t *ref_is_repeat; 

    uint8_t seed[L_MID_SEED];    
    
    //ref_is_repeat = (bp_t *)bp_init(bwt->seq_len);/*indicate repeat seeds on reference*/
  
    //ref_is_repeat = (uint8_t *)calloc(bwt->seq_len+2, sizeof(uint8_t));/* indicate repeat seed on bwt index */
fprintf(stderr, "bwt->seq_len = %u\n", bwt->seq_len);
    int n_multi_seeds = 0;//number of multi_seeds(occurence of seed >= 2)
    int n_seedid = 0;
/* 
    for(ed =bwt->seq_len-1; ed>=L_MID_SEED-1; --ed){
        if(ref_is_repeat[ed] ==2 ) continue;//skip repeat seed if not 1st occurence
        st = ed+1-L_MID_SEED;
        for(i = 0; i < L_MID_SEED; ++i){ seed[i] = __get_pac(pac, st+i);}
        
        k = 0, l = bwt->seq_len;
        uint32_t __iter = lkt_seq2LktItem(seed, 12,20);
        k = fastmap->item[__iter]; 
        l = fastmap->item[__iter+1];
        if(bwt_match_exact_alt(bwt, 8, seed, &k, &l) <=1) {
		    ref_is_repeat[ed]=1;//set end pos of seed   
		    continue;//skip uniq seeds 
        }
	    //log_bt2nt(20, seed);
*/
/*  
        if(st >= 16&& ed+16<bwt->seq_len){
            uint8_t __seed[128] = {};
            for(i = 0; i<= 52; ++i) __seed[i] = __get_pac(pac, st-16+i);
            uint32_t __k1=0, __k2=0, __k3=0, __l1=bwt->seq_len, __l2=bwt->seq_len, __l3=bwt->seq_len;
            bwt_match_exact_alt(bwt, 20, __seed+16, &__k1, &__l1); 
            bwt_match_exact_alt(bwt, 36, __seed+16, &__k2, &__l2); 
            bwt_match_exact_alt(bwt, 52, __seed, &__k3, &__l3); 

            //printf(">seq%u_(%u,%u)_(%u,%u)_(%u,%u)\n", n_seedid, __k1, __l1, __k2, __l2, __k3, __l3);
            //log_bt2nt(52, __seed);
            ++n_seedid;
        } 
*/
/*
        while((pos = bwt_sa(bwt, k)) + L_MID_SEED-1 > bwt->seq_len){
            k++; 
        }
        if(k +1 <=l){
            for(isa = k; isa <=l; ++isa){
                pos = bwt_sa(bwt, k); 
                ref_is_repeat[pos+L_MID_SEED-1]=2;//set end pos of seed   
                bwt_is_repeat[0][isa] = isa==k?2:1;
            }
            ++n_multi_seeds;
        }
        
    }
*/
/*
    for(ed =bwt->seq_len; ed < bwt->seq_len+L_MID_SEED; ++ed){
        st = ed+1-L_MID_SEED;
        memset(seed, 0, L_MID_SEED);
        for(i = 0; i < bwt->seq_len-st; ++i){ seed[i] = __get_pac(pac, st+i);}
        k = 0, l = bwt->seq_len;
        int n = bwt_match_exact_alt(bwt, 20, seed, &k, &l);
        if(n == 1){
        
        }
	    
        } 
    }
*/
    build_bwt_is_repeat(prefix, idx, idx->pac);

//fprintf(stderr, "1702\n");
    //bp_destroy(ref_is_repeat);
    //free(ref_is_repeat);
//fprintf(stderr, "1709\n");
    idx_destroy(idx);

//fprintf(stderr, "1712\n");
/*
    lkt_destroy(fastmap);
    bwt_destroy(bwt);
    bns_destroy(bntseq);
    free(pac);
*/
/*  
    printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    for(i = 0; i < MAX_COUNT; ++i){
        printf("%u\t%u\n", i, global_stat_20seed0[i]); 
    }
    for(i = 0; i < MAX_COUNT; ++i){
        printf("%u\t%u\n", i, global_stat_20seed1[i]); 
    }

    for(i = 0; i < MAX_COUNT; ++i){
        printf("%u\t%u\n", i, global_stat_36seed[i]); 
    }

*/
/*
    free(global_stat_36seed);
    free(global_stat_20seed0);
    free(global_stat_20seed1);
*/

}


#define  LEN_SEED  20 
#define  OFF_SEED  10 
#define  LEN_EXT   16 
#define  LEN_READ  150
#define  NUM_EXT  (((LEN_READ- LEN_SEED )/2)+OFF_SEED+ LEN_EXT-1)/LEN_EXT  
#define  MAX_SEED_NUM 600000
#define  LEN_FILE_NAME 100 
#define  IS_SMLSIZ 16
#define  LEN_EXT_IDX  
//LEN_EXT_IDX 是所有ExtIdx[][3]数组长度中最大的，可以用文件大小来计算。

//--------------------------------------------------------------
//扩展种子Index区间生成过程：
//主BWT索引已经创建，包括:
//Bwt后缀前置序列数组SA[]，
//累加数组Rank[]，
//Index到Pos的映射数组sa[]数组，
//参考基因组序列RefSeq[]数组，
//-----------------------------------------------------------------

typedef struct{
	uint32_t relat;	
	uint8_t  smbwt;
	//uint32_t smpos;
    uint32_t nxtcap;	

	uint32_t num_seqL ;  // 当前块中seqL数据个数；
	uint32_t num_seqR ;  // 当前块中seqR数据个数；
	uint32_t num_relat;  // 当前块中relat数据个数；	
	uint32_t num_pos ; 	 // 当前块中pos数据个数；

} CapInfo;

typedef struct{  

	uint32_t (*seqL_buf)[2];  //截取的数据缓冲区
	uint32_t (*seqR_buf)[2];
//+++++++++++++++++++++++++++++++++++++++++
	uint32_t *sortL; //排序后的唯一序列
	uint32_t *sortR;
	uint32_t *relat ;	
	uint32_t *L2rel ;
	uint32_t (*idxR)[2]; // sortR[]的每一个序列的(bgnIdx,ednIdx)

	//注：idxR[][2]的数据存放在;seqL_buf[][2]的物理空间。
//++++++++++++++++++++++++++++++++++++++++++++++++++++
	uint32_t (*extIdx)[2]; // sortR[]的每一个序列的(bgnIdx,ednIdx)
	//注：extIdx[][2]的数据存放在;seqR_buf[][2]的物理空间。	
	uint32_t (*seqL_idxR)[3];  // relat[]的每个行上左序列和右序列对应的idx,即(seqL,bgnIdxR,endIdxR)	

	uint32_t *pos_buf;  // 保存当前快中的pos数据； 	
	uint32_t *nxt_idx;  // 保存当前快中的nxt_idx数据；
	uint8_t  *nxt_flg;  // 保存当前快中的nxt_flg数据；


	uint8_t  *bwt_seqL;   
	uint8_t  *bwt_sumL; 
	uint8_t  *bwt_seqR;   
	uint8_t  *bwt_sumR;
	uint8_t  *smbwt  ;

	uint32_t num_seqL ;  // 当前块中seqL数据个数；
	uint32_t num_seqR ;  // 当前块中seqR数据个数；
	uint32_t num_relat;  // 当前块中relat数据个数；	
	uint32_t num_pos ; 	 // 当前块中pos数据个数；

	uint8_t  num_bwtseqL;   
	uint8_t  num_bwtsumL; 
	uint8_t  num_bwtseqR;   
	uint8_t  num_bwtsumR;

	uint32_t num_bwtL;
	uint32_t num_bwtR;
	uint32_t len_capidx ;  // 当前块前为止capidx数组中的全局capidx数据长度；
	uint32_t len_relat  ;  // 当前块前为止relat数组中的全局relat数据长度；	
	uint32_t len_smbwt  ;  // 当前块前为止smbwt数组中的全局smbwt数据长度；	
	uint32_t len_nxtpnt ;  // 当前块前为止nxtpnt数组中的全局nxtpnt数据长度；
	uint32_t len_nxtflg ;  // 当前块前为止nxtflg数组中的全局flg数据长度；
	uint32_t len_smpos  ;  // 当前块前为止smpos数组中的全局pos数据长度；

	CapInfo cap[1];
	//uint32_t off_pos ; 
    //uint32_t blck_id;

} SubBuf;  


void setFileName(char *in_fname, FileName *out_fname,int flg );
uint32_t getlen(uint32_t LenExtIdx[],FileName *f);
void InitSubBuf(SubBuf *sub_buf);
void CrtJmpMod(int n, uint32_t cur_seedidx[][3], uint32_t* jmp_idx, uint8_t *mod_idx, uint32_t *cap_pos);

void OpenSeedIdxfiles(FileName *f, uint32_t (*sidx)[3], uint32_t LenExtIdx[],uint32_t repNum,FILE *fpw[6]);
void relatNxtCapIdx(SubBuf *subBuf,uint32_t *jmp_idx, uint8_t *jmp_mod,uint32_t *cap_pos,uint32_t *sa);
void comFile(FileName *rf, FileName *wf);
long getFileSize(char* file);
void getBlckData(uint32_t sidx[2],uint32_t num_ext, SubBuf *sub_buf, idx_t *idx);

void buldSmBwt(SubBuf *sub,int flg);
int build_smbwt_idx(idx_t *fm_idx){
    uint32_t *sa=fm_idx->bwt->sa;
	//+++++++++++++++++++++++++++++++++++++++++++
	//初始化
	char in_f[1024];
	FileName  f;
	strcpy(in_f,"seedidx");
	setFileName(in_f,&f,0);

	uint32_t LenExtIdx[NUM_EXT];
	uint32_t MaxLen; 
	uint32_t i = getlen(LenExtIdx,&f);  
	//LenExtIdx是所有ExtIdx[][3]数组长度中最大的，
	//可以用文件大小来计算,返回最大值序号。
	MaxLen = LenExtIdx[i];

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//内存初始化------------------------------------------------------
	uint32_t (*cur_seedidx)[3], (*ord_seedidx)[3], (*nxt_seedidx)[3];
	uint32_t *jmp_idx, *cap_pos;
	uint8_t  *mod_idx; 

	if(NULL == (cur_seedidx = (uint32_t(*)[3])malloc(MaxLen*3*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}

	if(NULL == (nxt_seedidx = (uint32_t(*)[3])malloc(MaxLen*3*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}

	if(NULL == (jmp_idx = (uint32_t*)malloc(((MaxLen+255)/256)*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}


	if(NULL == (mod_idx = (uint8_t*)malloc(MaxLen*sizeof(uint8_t)))){
	    perror("error...");
	    exit(1);
	}

	if(NULL == (cap_pos = (uint32_t*)malloc(MaxLen*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}

	//==============================================================
	//
	//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//数组初始化----------------------------------------------------
	fprintf(stderr, "数组初始化！\n");
    SubBuf sub;
	InitSubBuf(&sub);

	CapInfo  *cap;
	cap = sub.cap;

	FILE *fp;
	if((fp=fopen64(f.seedidx[0],"r"))==NULL){
	    printf("can't open file %s\n", f.seedidx[0]);
	    exit(0);
	}
	uint32_t numread=fread(cur_seedidx,3*sizeof(uint32_t),LenExtIdx[0],fp);
	//printf("Number of items read=%d\n",numread);
	fclose(fp);
	//================================================================
	//计算cur_seedidx的jmp_idx,mod_idx, cap_pos数组，并输出到FlgIdxFile_1
	CrtJmpMod(LenExtIdx[0], cur_seedidx, jmp_idx, mod_idx,cap_pos);
	//FileName  jf;
	strcpy(in_f,"jmpmod");
	setFileName(in_f,&f,1);

	if((fp=fopen64(f.jmpmod,"w"))==NULL){//判断是否打开文件
	    printf("can't open file\n");
	    exit(0);
	}
	fwrite(cap_pos,sizeof(uint32_t),LenExtIdx[0],fp);  
	fwrite(jmp_idx,sizeof(uint32_t),(LenExtIdx[0]+512-1)/256,fp);  
	fwrite(mod_idx,sizeof(uint32_t),(LenExtIdx[0]+3)/4,fp);  
	fclose(fp); 

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//开始循环处理----------------------------------------------------

	strcpy(in_f,"idxfile");
	setFileName(in_f,&f,1);
	//return 0;
    FILE *fpw[6], *fp_capidx, *fp_relat, *fp_smbwt, *fp_nxtpnt, *fp_nxtflg, *fp_smpos;
	fprintf(stderr, "loop!\n");
 	uint8_t repNum;
    for(repNum=0; repNum < NUM_EXT; ++repNum){
		OpenSeedIdxfiles( &f, nxt_seedidx,LenExtIdx,repNum, fpw);
		fp_relat  = fpw[0];
		fp_smbwt  = fpw[1];
		fp_nxtpnt = fpw[2];
		fp_nxtflg = fpw[3];
		fp_capidx = fpw[4];
		fp_smpos  = fpw[5];

		CrtJmpMod(LenExtIdx[repNum], nxt_seedidx, jmp_idx, mod_idx,cap_pos);

		uint32_t bgn,end,num,j,idx,len_seedidx ;


		len_seedidx = LenExtIdx[repNum];

		for(idx =0; idx<len_seedidx; ++idx){
			bgn = cur_seedidx[idx][0];
			num = cur_seedidx[idx][1];
			end = bgn + num-1; 

			if(num <= IS_SMLSIZ){
				sub.num_pos = 0 ;
				for(j=bgn; j<=end; j++) sub.pos_buf[sub.num_pos++]=sa[j]; 
                fwrite(sub.pos_buf,sizeof(uint32_t), sub.num_pos,fp_smpos);
                continue;
				//小规模数据不做任何处理
				//比对过程中可以通过自身的jmp_idx,mod_idx信息，
				//可以找到相应的SmPos[]数组的行号
			} 
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//以下是当num_idx>IS_SMLSIZ时进行------------------------- 
			//生成左扩展序列集合SeqL[]和右扩展序列集合SeqR[]，
			//计算关联关系矩阵数组Relat[ ]和映射数组LtoRel[ ], 

			uint32_t sidx[2];
			sidx[0] = bgn;
			sidx[1] = num;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//预编译期间临时屏蔽
	        getBlckData(sidx,(uint32_t)repNum,&sub,fm_idx);  
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	        

			//Relat[  ]和LtoRel[ ]分段输出到同一个文件RelatFile_i。

			fwrite(sub.relat,sizeof(uint32_t),sub.num_relat,fp_relat);
			// buf_relat是当前Relat[]数据的缓冲区
			// len_relat是当前buf_relat数据的个数
			fwrite(sub.L2rel,sizeof(uint32_t),sub.num_seqL,fp_relat);
			// buf_L_Rel是当前LtoRel[ ]数据的缓冲区
			// num_seq_L是当前SeqL数据的个数
			//num_seq_L, num_seq_R，len_relat, 分别保存在cap_idx的相应分量中；
		
			//---------------------------------------------------------	
			//计算SeqL[]的SmBwt,生成buf_bwt_seq_L和buf_bwt_sum_L;
			//len_bwtseq_L = (num_seq_L+1)/2;
			//len_bwtsum_L = getbwtsumsize(num_seq_L);


	        buldSmBwt(&sub,0); //0生成左Bwt，1生成右Bwt


			

			fwrite(sub.bwt_seqL,sizeof(uint8_t),sub.num_bwtseqL,fp_smbwt);
			fwrite(sub.bwt_sumL,sizeof(uint8_t),sub.num_bwtsumL,fp_smbwt);

			//计算SeqR[]的SmBwt,生成buf_bwt_seq_R和buf_bwt_sum_R;
			//len_bwtseq_R = (num_seq_R+1)/2;
			//len_bwtsum_R = getbwtsumsize(num_seq_R);


	        buldSmBwt(&sub,1); //0生成左Bwt，1生成右Bwtt
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


			fwrite(sub.bwt_seqR,sizeof(uint8_t),sub.num_bwtseqR,fp_smbwt);
			fwrite(sub.bwt_sumR,sizeof(uint8_t),sub.num_bwtsumR,fp_smbwt);
			
			//计算smbwt数组的指针增量  
			//cur_num_smbwt = len_bwtseq_L+ len_bwtsum_L + len_bwtseq_R+ len_bwtsum_R
				
			//-------------------------------------------------------
			//配对左序列seq_L[]与相关联的右序列seq_L[]的Indx区间
			//关联数据保存在buf_rel_idx[][]数组	

			relatNxtCapIdx(&sub,jmp_idx,mod_idx,cap_pos,sa);

			//输出nxt_idx[i],nxt_idx[i],和apIdx[]数据
			fwrite(sub.pos_buf,sizeof(uint32_t), sub.num_pos,fp_smpos);
            fwrite(sub.nxt_idx,sizeof(uint32_t),sub.len_relat,fp_nxtpnt);
			fwrite(sub.nxt_flg,sizeof(uint8_t), sub.len_relat,fp_nxtflg);

			cap[0].num_seqL  = sub.num_seqL;
			cap[0].num_seqR  = sub.num_seqR;
			cap[0].num_relat = sub.num_relat;

			fwrite(cap,sizeof(CapInfo),1,fp_capidx);	

	        cap[0].relat += sub.num_relat + sub.num_seqL ;	 
			cap[0].smbwt += sub.num_bwtL + sub.num_bwtR ; 
			cap[0].nxtcap+= sub.num_relat;
  
	        //计算relat数组的指针增量 
			sub.len_capidx++;

			//idx++;
		} // End: while(idx<len_seedidx) +++++++++++++++++++++++++++
			//---------------------------------------------------------
		fclose(fp_nxtpnt);
		fclose(fp_nxtflg);
		fclose(fp_capidx);
		fclose(fp_relat);
		fclose(fp_smbwt);
		fclose(fp_smpos);

		ord_seedidx = cur_seedidx;
		cur_seedidx = nxt_seedidx ;
		//nxt_seedidx = NULL;
		nxt_seedidx = ord_seedidx;
	}
	strcpy(in_f,"comfile");
	setFileName(in_f,&f,1);
    comFile(&f,&f);
	return 0;
}


//==================================================================
void comFile(FileName *rf, FileName *wf){
	//+++++++++++++++++++++++++++++++++++++++++++++
	FILE *fp_com;
	FILE *fp;
	int f_id;
	uint32_t num = 0;
	char  buf[2];

    int fsize[NUM_FILES+1][NUM_EXT];
    int head[NUM_FILES+1][NUM_EXT];

    char *fName[NUM_FILES][NUM_EXT] ;

	fName[0][0]	=  rf->capidx[0];
	fName[1][0]	=  rf->nxtpnt[0];
	fName[2][0]	=  rf->nxtflg[0];
	fName[3][0]	=  rf->relat[0];
	fName[4][0]	=  rf->smbwt[0];
	fName[5][0]	=  rf->smpos[0];

    int i,j;
	for(i=0;i<NUM_EXT;i++){
		fsize[1][i] =  getFileSize(rf->capidx[i]);
		fsize[2][i] =  getFileSize(rf->nxtpnt[i]);
		fsize[3][i] =  getFileSize(rf->nxtflg[i]);
		fsize[4][i] =  getFileSize(rf->relat[i]);
		fsize[5][i] =  getFileSize(rf->smbwt[i]);		
		fsize[6][i] =  getFileSize(rf->smpos[i]);
	}
	for(i=0;i<NUM_EXT;i++){
		head[0][i] = 5;
		head[1][i] = ((fsize[1][i]+3)/4)*4 ; 
		head[2][i] = ((fsize[2][i]+3)/4)*4 ; 
		head[3][i] = ((fsize[3][i]+3)/4)*4 ; 
		head[4][i] = ((fsize[4][i]+3)/4)*4 ; 
		head[5][i] = ((fsize[5][i]+3)/4)*4 ; 
		head[6][i] = ((fsize[6][i]+3)/4)*4 ; 
	}
	int num_ext = 0 ;
	while(num_ext<NUM_EXT){
		if((fp_com=fopen64(wf->comfile[num_ext],"w"))==NULL){
		    printf("can't open file\n");
		    exit(0);
		}
		f_id = 0;
		while(f_id<NUM_FILES){
			//++++++++++++++++++++++++++++++++++++++++++
			//file_capidx数据输出到file_extidx
			if((fp=fopen64(fName[f_id][num_ext],"r"))==NULL){
			    printf("can't open file\n");
			    exit(0);
			}
			num = 0;
			//whiel(num < head[i][1]){
			
			while((buf[0] = fgetc(fp)) != EOF){
				num++;
				fwrite(buf,sizeof(uint8_t),1,fp_com);
			}
			if(num != fsize[num_ext][1]){
				    printf("file read Err\n");
				    exit(0);
			}

			for(j=0; j<(head[num_ext][1]-fsize[num_ext][1]);j++ ){
				buf[0] = 0;
				fwrite(buf,sizeof(uint8_t),1,fp_com);
			}
			fclose(fp);
			f_id++;
			//++++++++++++++++++++++++++++++++++++++++++++
		}//End:  while(num_file<NumFiles) +++++++++++++++++++
		fclose(fp_com);
		num_ext++;
	}
	return;
} 
// End： void comDataFile(FileName *file) +++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void OpenSeedIdxfiles(FileName *f, uint32_t (*sidx)[3], uint32_t LenExtIdx[],uint32_t repNum,FILE *fpw[]){
	FILE *fp;
	if((fp=fopen64(f->seedidx[repNum+1],"r"))==NULL){//判断是否打开文件
	    fprintf(stderr, "can't open file %s\n", f->seedidx[repNum+1]);
	    exit(0);
	}
	uint32_t  numread=fread(sidx,3*sizeof(uint32_t),LenExtIdx[repNum+1],fp);
	//printf("Number of items read=%d\n",numread);
	fclose(fp);
	
	if((fpw[0]=fopen64(f->relat[repNum+1],"w"))==NULL){
	    fprintf(stderr, "can't open file %s\n", f->relat[repNum+1]);
	    exit(0);
	}
	
	if((fpw[1]=fopen64(f->smbwt[repNum+1],"w"))==NULL){
	    fprintf("can't open file %s \n", f->smbwt[repNum+1]);
	    exit(0);
	}

	if((fpw[2]=fopen64(f->nxtpnt[repNum+1],"w"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}

	if((fpw[3]=fopen64(f->nxtflg[repNum+1],"w"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}

	if((fpw[4]=fopen64(f->capidx[repNum+1],"w"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}

	if((fpw[5] = fopen64(f->smpos[repNum+1],"w"))== NULL){
	    printf("can't open file\n");
	    exit(0);
	}
	return;
} // End: void Openfiles( ) ++++++++++++++++++++++++++++++++

//*************************************************
//*************************************************


void CrtJmpMod(int n, uint32_t cur_seedidx[][3], uint32_t* jump_idx, 
	uint8_t *mod_idx, uint32_t *cap_pos){

    //uint32_t row_jump = 0;
    uint32_t last_idx_q = 0;
    //uint32_t jump_idx[LEN_NEXT_IDX+1]; 
    uint32_t row, i;
    for(row=0; row<n; ++row){
        cap_pos[row] = cur_seedidx[row][2]; 
        uint32_t idx = cur_seedidx[row][0];
        uint32_t idx_q = idx/256;
        mod_idx[row] = idx%256;
        if(idx_q == last_idx_q){
            continue;
        } else{
            for(i= last_idx_q+1; i <= idx_q; ++i){
                jump_idx[i] = row; 
            }
            last_idx_q = idx_q;
        }
    
    } 
    jump_idx[(n+255)/256] = n; 
    return;
}
/*
int cut_ext_seq(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, vec_uint_t *seq_buf)
{
    bwtint_t i, pos, st, ed;
    uint32_t seq16;
    //fprintf(stderr, "\n"); 
    //fetch 16 mer preseq in [k, l) and rm duplicate preseq
    for(i = k; i < l; ++i){
        pos = bwt_sa(bwt, i);
        st = pos+offset; ed = st+16;
        if(offset < 0 && pos < (bwtint_t)-offset) st = 0;
        if(offset > 0 && ed > bwt->seq_len) ed = bwt->seq_len; 
        seq16 = bns_extract_seq16(pac, st, ed);//???
        //fprintf(stderr, "%u\t", seq16); 
        kv_push(uint32_t, *seq_buf, seq16);
    }
    //fprintf(stderr, "\n");
    //for(i = seq16_st; i <seq16s->n; ++i) fprintf(stderr, "%u\t", seq16s->a[i]);
    //fprintf(stderr, "\n");
    return seq_buf->n;
}*/
int rm_euqal(vec_uint_t *seq_buf, bwtint_t k, vec_ext_t *rext)
{
    int i, n_uniq;
    uint32_t *b = seq_buf->a;
    
    ext_t *p = kv_pushp(ext_t, *rext);
    p->seq = b[0];
    p->idx = k;
 
    n_uniq = 1; 
    for(i = 1; i < seq_buf->n; ++i){
        if(b[i] > b[i-1]) {
            assert(rext->m >n_uniq*2+1);
            ext_t *p = kv_pushp(ext_t, *rext);
            p->seq = b[i];
            p->idx = k+i;
            //a[n_uniq*2] = b[i];
            //a[n_uniq*2+1] = k+i;          
            ++n_uniq;
        } 
    }
    return n_uniq;
}


void getBlckData(uint32_t sidx[2],uint32_t num_ext, 
        	SubBuf *sub_buf, idx_t *idx/*pubParm pParm*/ ){
//
#define MAX_COUNT 800000
    vec_uint_t seq_buf; 
    vec_ext_t rext, lext;
    kv_init(seq_buf); kv_init(rext);kv_init(lext);
    kv_resize(uint32_t, seq_buf, MAX_COUNT);
    kv_resize(ext_t, rext, MAX_COUNT);
    kv_resize(ext_t, lext, MAX_COUNT);
    ext_t tmp; 
    uint32_t *lext_uniq_idx = (uint32_t *)calloc(MAX_COUNT, sizeof(ext_t));       
    int tot_lext0 = 0, tot_lext1 = 0, tot_rext = 0;
 
    uint8_t *pac = idx->pac;
    bwt_t *bwt = idx->bwt;    
    //rext.n = 0; lext.n =0; seq_buf.n = 0;
    uint32_t k = sidx[0], l = k+sidx[1];
    cut_ext_seq(pac, bwt, k, l,  20, &seq_buf);//right seq buf
    rem_euqal(&seq_buf, k, &rext);//gen right extend seq 
    uint32_t i;
    for(i =0; i < rext.n; ++i) { 
        seq_buf.n = 0;
        int st, ed;
        st = rext.a[i].idx;
        ed = i+1 >= rext.n?l:rext.a[i+1].idx;
        cut_ext_seq(pac, bwt, st, ed, -16, &seq_buf);//left seq buf 
        ks_introsort(uint32_t, seq_buf.n, seq_buf.a); 
        rem_euqal(&seq_buf, st, &lext);
    }
    int row_r, row_l;
    /*
        fprintf(stderr, "----------------------\n");
 
        for(row_l = 0; row_l < lext.n; ++row_l){
            fprintf(stderr, "[lext]:\t%u\t", lext.a[row_l].idx);
            log_seq162nt(lext.a[row_l].seq); 
        }
        for(row_r = 0; row_r < rext.n; ++row_r){
            fprintf(stderr, "[rext]:\t%u\t", rext.a[row_r].idx);
            log_seq162nt(rext.a[row_r].seq); 
        }
    */  
    row_r = 1, row_l = 0; 
    while(row_l < lext.n && row_r < rext.n){
        if(lext.a[row_l].idx  < rext.a[row_r].idx) lext.a[row_l++].idx = row_r-1;
        else row_r++; 
    } 
    if(row_l < lext.n) for(i = row_l; i < lext.n; ++i) lext.a[i].idx = rext.n-1; 
    ks_introsort(ext, lext.n, lext.a);

    int n_uniq = 0;
    uint32_t last_seq16 = lext.a[0].seq;
    lext_uniq_idx[n_uniq++] = 0;
    for(i = 1; i < lext.n; ++i){
        if(lext.a[i].seq!= last_seq16) {
            //lext_uniq_idx[n_uniq++] = i+tot_lext0; 
            lext_uniq_idx[n_uniq++] = i; 
            last_seq16 = lext.a[i].seq;
        }
    }
    //dump_cnt(fp_totext, tot_lext0, tot_lext1, tot_rext);
    fprintf(stderr, "%u\t%u\t%u\n", tot_lext0, tot_lext1, tot_rext); 
    tot_lext1 += lext.n; 
    tot_lext0 += n_uniq;
    tot_rext += rext.n; 
    //fwrite(rext.a, sizeof(ext_t), rext.n, fp_rext0);         
    
    //for(i = 0; i < rext.n; ++i) log_seq162nt(rext.a[i].seq);

    for(i = 0; i < rext.n; ++i) {
        sub_buf->sortR[i] = rext.a[i].seq;
        sub_buf->idxR[i][0] = rext.a[i].idx;
        sub_buf->idxR[i][1] = i == rext.n?l:rext.a[i].idx-1;
    }
    for(i = 0; i < n_uniq; ++i) {

        sub_buf->sortL[i] = lext.a[lext_uniq_idx[i]].seq;
        sub_buf->L2rel[i] = lext_uniq_idx[i]; 
        //fwrite(&tmp, sizeof(ext_t), 1, fp_lext0);
        //fwrite(&lext.a[lext_uniq_idx[i]].seq, sizeof(uint32_t), 1, fp_lext0); 
        //fwrite(lext_uniq_idx+i, sizeof(uint32_t), 1, fp_lext0);
        fprintf(stderr, "%u\t", lext_uniq_idx[i]);
        log_seq162nt(lext.a[lext_uniq_idx[i]].seq);

    } 
    for(i = 0; i < lext.n; ++i) sub_buf->relat[i] = lext.a[i].idx;
    
    kv_destroy(seq_buf);
    kv_destroy(lext);
    kv_destroy(rext);
    return;
}

void build_local_bwt_alt(uint32_t *sort_seq, int l, int r, uint32_t *sa0, uint32_t *sa1, uint8_t *bwt)
{
    int i, j;
    uint32_t cnt_2nt[17];
    uint32_t *last_sa, *cur_sa;
    sort_seq += l;
    int n = r-l;
    for(j=0; j< 17; ++j) cnt_2nt[j] = 0;
    for(j = 0; j < n; ++j){
        uint32_t x = sort_seq[j]&0xF;
        __set_bwt(bwt, j*2+1, x&3);
        __set_bwt(bwt, j*2, (x>>2)&3);
        //__set_bwt2(bwt, j, x);
        ++cnt_2nt[x];
    }
    for(j =0; j < n; ++j) sa0[j] = j;
    last_sa = sa0; cur_sa = sa1;
    for(i = 1; i < 8; ++i){
        //fprintf(stderr, "\n==Iter = %u\n", i-1);

        //log_array(17, cnt_2nt);
        accumulate_cnt(17, cnt_2nt);
        //log_array(17, cnt_2nt);
        //log_array(r-l, last_sa);
        //log_array(r-l, cur_sa);
        //sort cur_sa 
        for(j = 0; j < n; ++j){
            //uint32_t x = (lext[last_sa[j]].seq>>((i-1)*4))&0xF;
            uint32_t x = __get_col(sort_seq[last_sa[j]], i-1);
            //fprintf(stderr, "%u->%u\n", x, cnt_2nt[x]);
            cur_sa[cnt_2nt[x]++]  = last_sa[j];
        }
        //generate bwt from cur_sa        
        for(j=0; j< 17; ++j) cnt_2nt[j] = 0;
        for(j = 0; j < n; ++j){
            //uint32_t x = (lext[cur_sa[j]].seq>>(i*4))&0xF;
            uint32_t x = __get_col(sort_seq[cur_sa[j]], i);
            //__set_bwt(bwt, i*n*2+j*2+1, x&3);
            //__set_bwt(bwt, i*n*2+j*2, (x>>2)&3);
            int n_2nt = (n+1)/2*2;
            __set_bwt2(bwt, i*n_2nt+j, x);
            ++cnt_2nt[x];
        }
        SWAP(uint32_t *, last_sa, cur_sa); 
    }
    fprintf(stderr, "\n==Iter = %u\n", i-1);
    //log_seq16array(n, last_sa, lext, i-1);
}

uint8_t *build_seq16_bwt_alt(int n, uint32_t *sort_seq)
{
   
    uint8_t *bwt; uint32_t *sa0, *sa1;
    bwt = (uint8_t *)calloc((n+1)/2*8, sizeof(uint8_t));  
    sa0 = (uint32_t *)calloc(n, sizeof(uint32_t));
    sa1 = (uint32_t *)calloc(n, sizeof(uint32_t));
    build_local_bwt_alt(sort_seq, 0, n, sa0, sa1, bwt); 
    free(sa0);free(sa1);
    build_occ(bwt, n); 
        
        
       //dump bwt
        
       
 
        

    return bwt;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define MAX_BLCK_SIZE 500000
#define MAX_SEED_SIZE 500000 
#define ALL_BUFS_SIZE 12  // All bufs size in SubBuf



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void InitSubBuf(SubBuf *sub){
/*  
	SubBuf *sub ;
    if(NULL == (sub = (SubBuf*)malloc(sizeof(SubBuf)))){
	    perror("error...");
	    exit(1);
	}
*/
	uint32_t *buf;
	if(NULL == (buf = (uint32_t*)malloc(MAX_BLCK_SIZE*ALL_BUFS_SIZE*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
	sub->seqL_buf = (uint32_t(*)[2])buf;
	sub->seqR_buf = (uint32_t(*)[2])(buf + MAX_BLCK_SIZE*2);

	sub->sortL    = buf + MAX_BLCK_SIZE*4;
	sub->sortR    = buf + MAX_BLCK_SIZE*5;

	sub->relat    = buf + MAX_BLCK_SIZE*6;
	sub->L2rel    = buf + MAX_BLCK_SIZE*7;
	sub->seqL_idxR= (uint32_t(*)[3])(buf + MAX_BLCK_SIZE*8);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	sub->idxR     = (uint32_t(*)[2])buf; //注：idxR[][2]的数据存放在;seqL_buf[][2]的物理空间。
	sub->extIdx   = (uint32_t(*)[2])(buf + MAX_BLCK_SIZE*2); //注：extIdx[][2]的数据存放在;seqR_buf[][2]的物理空间。

	

	sub->pos_buf = buf + MAX_BLCK_SIZE*9; 	
	sub->nxt_idx = buf + MAX_BLCK_SIZE*10;
	sub->nxt_flg = (uint8_t*)(buf + MAX_BLCK_SIZE*11);


//+++++++++++++++++++++++++++++++++++
	uint8_t *bwt;
	if(NULL == (bwt = (uint8_t*)malloc(12*MAX_BLCK_SIZE*sizeof(uint8_t)))){
	    perror("error...");
	    exit(1);
	}
	sub->smbwt    = bwt;
	sub->bwt_seqL = bwt;   
	sub->bwt_sumL = bwt + 4*MAX_BLCK_SIZE;
	sub->bwt_seqR = bwt + 6*MAX_BLCK_SIZE;    
	sub->bwt_sumR = bwt + 10*MAX_BLCK_SIZE;  

//+++++++++++++++++++++++++++++++++++++++++++++++++
	sub->num_seqL = 0;
	sub->num_seqR = 0;
	sub->num_relat= 0;	
	//sub->num_ext  = 0;

	//sub->off_pos  = 0; 
    //sub->blck_id  = 0; 

	sub->len_capidx= 0;
	sub->len_relat = 0;	
	sub->len_smbwt = 0;	
	sub->len_nxtpnt= 0;
	sub->len_nxtflg= 0;
	sub->len_smpos = 0;
	sub->num_pos   = 0;
	sub->num_bwtL  = 0;
	sub->num_bwtR  = 0;


	sub->cap[0].relat    = 0;	
	sub->cap[0].smbwt    = 0;
    sub->cap[0].nxtcap   = 0;	
	sub->cap[0].num_seqL = 0;
	sub->cap[0].num_seqR = 0;
	sub->cap[0].num_relat= 0;	
	sub->cap[0].num_pos  = 0;

//+++++++++++++++++++++++++++++++++++++++++++++++++
	//sub_buf = sub;

	return ;
}


void buldSmBwt(SubBuf *sub,int flg){
    uint8_t *bwt; uint32_t *sa0, *sa1;
    //bwt = (uint8_t *)calloc((n+1)/2*8, sizeof(uint8_t));  
    int n;
    uint32_t *sort_seq;
    bwt = flg==0?sub->bwt_seqL:sub->bwt_seqR;
    if(flg == 0 ){
        bwt = sub->bwt_seqL; 
        n = sub->num_seqL; 
        sort_seq = sub->sortL;
    } else{
        bwt = sub->bwt_seqR; 
        n = sub->num_seqR;
        sort_seq = sub->sortR;
    }
    sa0 = (uint32_t *)calloc(n, sizeof(uint32_t));
    sa1 = (uint32_t *)calloc(n, sizeof(uint32_t));
    build_local_bwt_alt(sort_seq, 0, n, sa0, sa1, bwt); 
    free(sa0);free(sa1);
    build_occ(bwt, n); 

        
       //dump bwt
        
       
 
        

   

  
	return;
}

//**************************************************
//*************************************************************
void relatNxtCapIdx(SubBuf *subBuf,uint32_t *jmp_idx, uint8_t *mod_idx,uint32_t *cap_pos,uint32_t *sa){
	uint32_t  i, j,	bgn,end, num,  bgn_q, end_q, i_r, 
			  seq, pos, row, j_idx;

	uint32_t buf_idx[5];
    SubBuf *sub;
	sub = subBuf;
    sub->num_pos = 0;
    for(i=0;i<sub->num_seqL;i++){
		bgn = sub->L2rel[i];
		end = sub->L2rel[i+1];
		for(j=bgn; j<end; j++){
			row = sub->relat[j];
			sub->seqL_idxR[j][0] = sub->sortL[i];
			sub->seqL_idxR[j][1] = sub->idxR[row][0]; 
			sub->seqL_idxR[j][2] = sub->idxR[row][1];										
		}				
	}
	//-------------------------------------------------------
	//buf_rel_idx[][]中的每一条数据，进行Bwt比对获得(bgnIdx,endIdx)的数据，
	//保存在sub.extIdx[][]
	for(i=0;i<sub->num_relat;i++){
		seq = sub->seqL_idxR[i][0];
		buf_idx[0] = sub->seqL_idxR[i][1];
		buf_idx[1] = sub->seqL_idxR[i][2];
//++++++++++++++++++++++++++++++++++++++++++++++
//预编译期间屏蔽
		//AlgnSeq(pParm,seq,buf_idx); //seq 是整数
//++++++++++++++++++++++++++++++++++++++++++++++++
		sub->extIdx[i][0]= buf_idx[3]; //返回的bgnIdx
		sub->extIdx[i][1]= buf_idx[4]; //返回的endIdx
	}
	//---------------------------------------------------------
	//buf_rel_idx[][]中的每一条数据，进行Bwt比对获得(bgnIdx,endIdx)的数据，
	//保存在sub.extIdx[][]
	sub->num_pos = 0;
	for(i=0;i<sub->num_relat;i++){
		bgn = sub->extIdx[i][0];
		end = sub->extIdx[i][1];
		num = end - bgn + 1;
		if(num==1){
			pos = sa[bgn];
			sub->nxt_idx[i] = pos;
			sub->nxt_flg[i] = num;
			continue;
		}
		//if(num_idx>1)+++++++++++++++++++++++++++++++++++++++++++++++++++
		bgn_q = jmp_idx[i/256];
		end_q = jmp_idx[(i/256)+1];
		i_r = i%256;
		for(j_idx=bgn_q; j_idx<end_q; j_idx++){
			if(i_r == mod_idx[j_idx]){
				break;
			} 
		}		
		sub->nxt_idx[i] = cap_pos[j_idx]; //指向下一级CapIdx[]的行号
		if(num < 255){
			sub->nxt_flg[i] = (uint8_t)num ; //指向CapIdx[]的类型
		}else{
			sub->nxt_flg[i] = 0xFF ;
		}
		if(num<=IS_SMLSIZ){
			for(j=bgn; j<end; j++){
				pos = sa[j]; 
				sub->pos_buf[sub->num_pos++] = pos; 
			}
		}
	} // End : for(i=0;i<len_relat;i++)
	return;
}

uint32_t getlen(uint32_t l_extidx[], FileName *f)
{
	int i;
	int max_len = 0, max_idx = -1;
	for(i =0; i < NUM_EXT; ++i){
        long len = getFileSize(f->seedidx[i]);
        l_extidx[i] = (uint32_t)len/(3*sizeof(uint32_t));
		if(len > max_len) {
			max_len = len;
			max_idx = i; 
		}

	}
	return max_idx;

}
long getFileSize(char *fn)
{
	FILE *fp = fopen64(fn, "r");
	fseek(fp, 0, SEEK_END);
	long len = ftell(fp);
	fclose(fp);
    return len; 


}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_extend_idx
 *  Description:  build extend index(extend 20bp snpaware seed to 52bp seed)
 *      @prefix:  prefix of index file name 
 * =====================================================================================
 */

void build_extend_idx(const char *prefix)
{

    global_stat_36seed = (int *)calloc(MAX_COUNT, sizeof(int));
    global_stat_20seed0 = (int *)calloc(MAX_COUNT, sizeof(int));
    global_stat_20seed1 = (int *)calloc(MAX_COUNT, sizeof(int));
    /* 
     *
     *
     * */
   
    char fn0[MAX_NAME];
    double t = realtime();
    
    /* reload FM-index and packed reference */
    fprintf(stderr, "reload FM-index..\n");
    int64_t i;
    bwt_t *bwt; bntseq_t *bntseq; uint8_t *pac;
    bwtint_t ed, st, k, l, pos, isa;

    cat_filename(fn0, prefix, ".bwt"); 
    bwt = bwt_restore_bwt(fn0);   
    cat_filename(fn0, prefix, ".sa"); 
    bwt_restore_sa(fn0, bwt);
    bntseq = bns_restore(prefix);
    pac = (uint8_t *)calloc(bntseq->l_pac/4+1, sizeof(uint8_t));
    fread(pac, bntseq->l_pac/4+1, sizeof(uint8_t), bntseq->fp_pac);
    fprintf(stderr, "%.2f sec\n", realtime()-t); 
    /* generate repeat seeds and bwt range */
    fprintf(stderr, "gen repeat on bwt\n");
    bp_t *ref_is_repeat; 
    uint8_t *bwt_is_repeat;
    uint8_t seed[L_MID_SEED];    
    
    ref_is_repeat = (bp_t *)bp_init(bwt->seq_len);/*indicate repeat seeds on reference*/
    bwt_is_repeat = (uint8_t *)calloc(bwt->seq_len+1, sizeof(uint8_t));/* indicate repeat seed on bwt index */
    /*  backward search mid seed on refernece and generate bwt_is_repeat array;
     *  bwt_is_repeat[i] == 0 if uniq
     *  bwt_is_repeat[i] == 1 if sequence Ref[SA[i]:SA[i]+L_MID_SEED] is occur more than once on ref and i is not begin idx on bwt
     *  bwt_is_repeat[i] == 2 if sequence Ref[SA[i]:SA[i]+L_MID_SEED] is occur more than once on ref and i is begin idx on bwt
     * */ 
    int n_multi_seeds = 0;//number of multi_seeds(occurence of seed >= 2)
    int n_seedid = 0;
    for(ed =bwt->seq_len-1; ed>=L_MID_SEED-1; --ed){
        if(bp_get(ref_is_repeat, ed) ==1 ) continue;//skip repeat seed if not 1st occurence
        st = ed+1-L_MID_SEED;
        for(i = 0; i < L_MID_SEED; ++i){ seed[i] = __get_pac(pac, st+i);}
        k = 0, l = bwt->seq_len;
        if(bwt_match_exact_alt(bwt, L_MID_SEED, seed, &k, &l) <=1) continue;//skip uniq seeds 
        //log_bt2nt(20, seed);
        //printf("%u %u\n", k, l);
        if(st >= 16&& ed+16<bwt->seq_len){
            uint8_t __seed[128] = {};
            for(i = 0; i<= 52; ++i) __seed[i] = __get_pac(pac, st-16+i);
            uint32_t __k1=0, __k2=0, __k3=0, __l1=bwt->seq_len, __l2=bwt->seq_len, __l3=bwt->seq_len;
            bwt_match_exact_alt(bwt, 20, __seed+16, &__k1, &__l1); 
            bwt_match_exact_alt(bwt, 36, __seed+16, &__k2, &__l2); 
            bwt_match_exact_alt(bwt, 52, __seed, &__k3, &__l3); 

            printf(">seq%u_(%u,%u)_(%u,%u)_(%u,%u)\n", n_seedid, __k1, __l1, __k2, __l2, __k3, __l3);
            log_bt2nt(52, __seed);
            ++n_seedid;
        } 
   
        for(isa = k; isa <=l; ++isa){
            pos = bwt_sa(bwt, isa); 
            bp_set1(ref_is_repeat, pos+L_MID_SEED-1);//set end pos of seed   
            bwt_is_repeat[isa] = isa==k?2:1;
        }
        ++n_multi_seeds;
    }
    fprintf(stderr, "%.2f sec\n", realtime()-t); 
    fprintf(stderr, "gen uniq 16mer\n");
    
    /* fetch repeat seeds from bwt_is_repeat and generate uniq preseq  */
    char fn_rext0[MAX_NAME], fn_totext[MAX_NAME], fn_lext0[MAX_NAME], fn_lext1[MAX_NAME];
    FILE *fp_lext0, *fp_lext1, *fp_rext0, *fp_totext; 
    strcat(strncpy(fn_lext0, prefix, MAX_NAME), ".lext0");
    strcat(strncpy(fn_lext1, prefix, MAX_NAME), ".lext1");
    strcat(strncpy(fn_rext0, prefix, MAX_NAME), ".rext0");
    strcat(strncpy(fn_totext, prefix, MAX_NAME), ".totext");
    fp_lext0 = xopen(fn_lext0, "w"); 
    fp_lext1 = xopen(fn_lext1, "w"); 
    fp_rext0 = xopen(fn_rext0, "w"); 
    fp_totext = xopen(fn_totext, "w"); 
    
    vec_uint_t seq_buf; 
    vec_ext_t rext, lext;
    kv_init(seq_buf); kv_init(rext);kv_init(lext);
    kv_resize(uint32_t, seq_buf, MAX_COUNT);
    kv_resize(ext_t, rext, MAX_COUNT);
    kv_resize(ext_t, lext, MAX_COUNT);
    
    ext_t tmp;       
    uint32_t *lext_uniq_idx = (uint32_t *)calloc(MAX_COUNT, sizeof(ext_t));       
    int tot_lext0 = 0, tot_lext1 = 0, tot_rext = 0;
    k = 0; l = 0;//extend seq [k, l)
    while(k < bwt->seq_len){
        if(bwt_is_repeat[k] ==0 ){
            ++k;
            continue;
        } else if(bwt_is_repeat[k] == 2){
            l = k + 1;
            while(l < bwt->seq_len && bwt_is_repeat[l]==1 ){++l;}
            /* generate right extend seq  and left extend seq in [l, r)*/ 
            rext.n = 0; lext.n =0; seq_buf.n = 0;
            cut_ext_seq(pac, bwt, k, l,  L_MID_SEED, &seq_buf);//right seq buf
            rem_euqal(&seq_buf, k, &rext);//gen right extend seq 
            for(i =0; i < rext.n; ++i) { 
                seq_buf.n = 0;
                int st, ed;
                st = rext.a[i].idx;
                ed = i+1 >= rext.n?l:rext.a[i+1].idx;
                cut_ext_seq(pac, bwt, st, ed, -16, &seq_buf);//left seq buf 
                ks_introsort(uint32_t, seq_buf.n, seq_buf.a); 
                rem_euqal(&seq_buf, st, &lext);
            }
            fprintf(stderr, "----------------------\n");
            int row_r, row_l;
            for(row_l = 0; row_l < lext.n; ++row_l){
                fprintf(stderr, "[lext]:\t%u\t", lext.a[row_l].idx);
                log_seq162nt(lext.a[row_l].seq); 
            }
            for(row_r = 0; row_r < rext.n; ++row_r){
                fprintf(stderr, "[rext]:\t%u\t", rext.a[row_r].idx);
                log_seq162nt(rext.a[row_r].seq); 
            }
            
            row_r = 1, row_l = 0; 
            while(row_l < lext.n && row_r < rext.n){
                if(lext.a[row_l].idx  < rext.a[row_r].idx) lext.a[row_l++].idx = row_r-1;
                else row_r++; 
            } 
            if(row_l < lext.n) for(i = row_l; i < lext.n; ++i) lext.a[i].idx = rext.n-1; 
            ks_introsort(ext, lext.n, lext.a);

            int n_uniq = 0;
            uint32_t last_seq16 = lext.a[0].seq;
            lext_uniq_idx[n_uniq++] = 0;
            for(i = 1; i < lext.n; ++i){
                if(lext.a[i].seq!= last_seq16) {
                    //lext_uniq_idx[n_uniq++] = i+tot_lext0; 
                    lext_uniq_idx[n_uniq++] = i; 
                    last_seq16 = lext.a[i].seq;
                }
            }
            dump_cnt(fp_totext, tot_lext0, tot_lext1, tot_rext);
            fprintf(stderr, "%u\t%u\t%u\n", tot_lext0, tot_lext1, tot_rext); 
            tot_lext1 += lext.n; 
            tot_lext0 += n_uniq;
            tot_rext += rext.n; 
            fwrite(rext.a, sizeof(ext_t), rext.n, fp_rext0);         
            for(i = 0; i < rext.n; ++i) log_seq162nt(rext.a[i].seq);

            for(i = 0; i < n_uniq; ++i) {

                tmp.seq = lext.a[lext_uniq_idx[i]].seq;
                tmp.idx = lext_uniq_idx[i]; 
                fwrite(&tmp, sizeof(ext_t), 1, fp_lext0);
                //fwrite(&lext.a[lext_uniq_idx[i]].seq, sizeof(uint32_t), 1, fp_lext0); 
                //fwrite(lext_uniq_idx+i, sizeof(uint32_t), 1, fp_lext0);
                fprintf(stderr, "%u\t", lext_uniq_idx[i]);
                log_seq162nt(lext.a[lext_uniq_idx[i]].seq);

            }               
            for(i = 0; i < lext.n; ++i) {
                fwrite(&lext.a[i].idx, sizeof(uint32_t), 1, fp_lext1); 
                fprintf(stderr, "%u\t", lext.a[i].idx);
            }
            fprintf(stderr, "\n");
            //kv_push(uint32_t, isa2seq16, n_seqs);
            k = l;
        } else{/* should't be here! */}
    }
    
    tmp.seq = 0; tmp.idx = bwt->seq_len+1;
    fwrite(&tmp, sizeof(ext_t), 1, fp_rext0);
    dump_cnt(fp_totext, tot_lext0, tot_lext1, tot_rext);
    
    err_fclose(fp_lext0);
    err_fclose(fp_lext1);
    err_fclose(fp_rext0);
    err_fclose(fp_totext);
    
    free(lext_uniq_idx); 
    kv_destroy(seq_buf);
    kv_destroy(lext); kv_destroy(rext);   
    fprintf(stderr, "%.2f sec\n", realtime()-t); 
    //bp_destroy(bitmap);
    /* use bwt_is_repeat to generate bwt_is_repeat2*/
    bp_t *bwt_is_repeat2 = bp_init(bwt->seq_len+1); /* bwt_is_repeat2[i] if Ref[sa[i]:sa[i]+20] is a multi seed*/
    for(st= 0; st < bwt->seq_len; ++st){
        if(bwt_is_repeat[st] == 2) bp_set1(bwt_is_repeat2, st);
    }
    
    cat_filename(fn0, prefix, ".repeat"); 
    bp_gen_occ(bwt_is_repeat2);
    bp_dump(bwt_is_repeat2, fn0);
    
    bp_destroy(bwt_is_repeat2);
    free(bwt_is_repeat);
    bp_destroy(ref_is_repeat);
    bwt_destroy(bwt);
    bns_destroy(bntseq);
    free(pac);
/*  
    printf("\n++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    for(i = 0; i < MAX_COUNT; ++i){
        printf("%u\t%u\n", i, global_stat_20seed0[i]); 
    }
    for(i = 0; i < MAX_COUNT; ++i){
        printf("%u\t%u\n", i, global_stat_20seed1[i]); 
    }

    for(i = 0; i < MAX_COUNT; ++i){
        printf("%u\t%u\n", i, global_stat_36seed[i]); 
    }

*/
    free(global_stat_36seed);
    free(global_stat_20seed0);
    free(global_stat_20seed1);


}

int idx_build_core(const char *fn_fa, const char *fn_av, const char *prefix)
{


    fprintf(stderr, "build fm-idx\n");
    build_fmidx(fn_fa, prefix);
    //build_snpaware_idx(fn_av, prefix);
    fprintf(stderr, "build_extend_idx_alt\n");
    build_extend_idx_alt(prefix);
return 0;
    //build_extend_idx(prefix);/*

    fprintf(stderr, "begin restore idx\n");
    idx_t *idx = idx_restore(prefix);
/*
//++++++++++++++++++++++++++++++++++++++++++++++
uint32_t __k=0, __l=idx->bwt->seq_len;


    uint8_t ch[36] = {3,3,3,3,0,0,0,0,
                     3,3,3,3,3,3,3,3,
                     0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,3,3};
 

int is_aln = bwt_match_exact_alt(idx->bwt, 36, ch, &__k, &__l);
printf("aln = %d, %u, %u\n", is_aln, __k, __l);
*/
//+++++++++++++++++++++++++++++++++++++++++++++
    fprintf(stderr, "begin build_smbwt_idx!\n");
    build_smbwt_idx(idx);
    idx_destroy(idx);

    /*         
*/
           return 0;
}
void test_idx_reload(const char *prefix)
{
    bwtint_t i, j;
   
    
    idx_t *idx = idx_restore(prefix);
    bwt_t *bwt = idx->bwt;
    fprintf(stderr, "\nTest idx reload!\n"); 

    fprintf(stderr, "====FM-index info====\n"); 
    fprintf(stderr, "Bwt seq length:%u\n", idx->bwt->seq_len); 
    fprintf(stderr, "Bwt sa INTV :%u\n", idx->bwt->sa_intv); 
    
    fprintf(stderr, "====Reference info====\n"); 
    fprintf(stderr, "chromsome num :%u\n", idx->bns->n_seqs); 
    fprintf(stderr, "chromsome1 name :%s\n", idx->bns->anns[0].name); 
    
    fprintf(stderr, "\n********************************************\n"); 
    fprintf(stderr, "log information of snpaware index\n"); 
    /*
    for(i = 0; i < 0xFFF; ++i){
        log_seq162nt(i);
        int st = idx->pmap[i];
        int  ed = idx->pmap[i+1];
        for(j = st; j < ed; ++j){
            fprintf(stderr, "%u\t", idx->sai[j]); 
        }
        fprintf(stderr, "\n"); 
    }
    
    */
    
    int n = 36;
    fprintf(stderr, "\n********************************************\n"); 
    fprintf(stderr, "log information of approx index\n"); 
    bp_t *is_multi = idx->is_multiseeds;
    /* print BWM and multi seeds */
    for(i = 0; i <= idx->bwt->seq_len; ++i){
        bwtint_t pos = bwt_sa(bwt, i);
        if(pos == -1) pos = bwt->seq_len;
        fprintf(stderr, "%u\t%u\t%u\t", i, bp_get(is_multi, i), pos);       
        for(j =0; j <n;++j){ 
            if(pos + j >= bwt->seq_len) break; 
            fprintf(stderr, "%c", "ACGT"[__get_pac(idx->pac, pos+j)]);
        }
        fprintf(stderr, "\t");
        for(j = pos<16?0:pos-16; j <pos;++j){ 
            if(j >= bwt->seq_len) break; 
            fprintf(stderr, "%c", "ACGT"[__get_pac(idx->pac, j)]);
        }

        
        
        fprintf(stderr, "\n");
    }
    /* 
    //print preseq
    uint32_t *seq16_idx = idx->isa2seq16;
    uint32_t *seq16s = idx->seq16s;
    for(i = 0; i <= idx->bwt->seq_len; ++i){
        if(bp_get(is_multi, i) == 0 ) continue;
        int k = bp_rank(is_multi, i)-1;
        fprintf(stderr, "%u\t", i);
        for(j = seq16_idx[k]; j < seq16_idx[k+1]; ++j){
            log_seq162nt(seq16s[j]);
        } 
    }
    */
    
    for(i=0; i < idx->n_tot-1; ++i){
        fprintf(stderr, "[multi seed id %u]\n", i);
        for(j = idx->lext_idx[i]; j < idx->lext_idx[i+1]; ++j){
            fprintf(stderr, "[%u]\t", idx->lext0[j].idx);
            log_seq162nt(idx->lext0[j].seq); 
        }
        for(j = idx->lext1_idx[i]; j <idx->lext1_idx[i+1]; ++j){
            //fprintf(stderr, "[lext1]: %u\t", idx->lext1[j]); 
        }
        fprintf(stderr, "\n");
        for(j = idx->rext_idx[i]; j < idx->rext_idx[i+1]; ++j){
            fprintf(stderr, "[%u]\t", idx->rext0[j].idx);
            log_seq162nt(idx->rext0[j].seq); 
        } 
    
    }
    





    idx_destroy(idx);
}
int index_main ( int argc, char *argv[] )
{
    int c;
    while((c = getopt(argc, argv, "k:h"))>=0){
        switch(c){
            case 'h':
                return index_usage();
            default:
         return 1;
              
        }
    
    }
    if (argc - optind  != 3 ){
        return index_usage();
    }
 
    
    
    
    idx_build_core(argv[optind], argv[optind+1], argv[optind+2]);

    
    return 0;
}				/* ----------  end of function main  ---------- */


#ifdef MAIN_INDEX
int main(int argc, char *argv[]){

    return index_main(argc, argv);
}
#endif
