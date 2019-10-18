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

#define MAX_NAME 128
const int SA_INTV = 1;

int *global_stat_20seed0;
int *global_stat_20seed1;
int *global_stat_36seed;

KHASH_MAP_INIT_STR(str, int)
typedef kvec_t(uint32_t) vec_uint_t;
KSORT_INIT_GENERIC(uint32_t)
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
        putchar("ACGT"[bt[i]]); 
    }
    putchar('\n');
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
    bwt_restore_isa(fn0, idx->bwt); 
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".lkt");
    idx->fastmap = lkt_restore(fn0); 
    /* restore snp-aware index */ 
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".vseq");
    FILE *fp_varseq = err_xopen_core(__func__, fn0, "r");
    fread(&idx->n_pmap, sizeof(uint32_t), 1, fp_varseq);
    idx->pmap = (uint32_t *)calloc(idx->n_pmap, sizeof(uint32_t));
    fread(idx->pmap, sizeof(uint32_t), idx->n_pmap, fp_varseq);
    fread(&idx->n_var, sizeof(uint32_t), 1, fp_varseq);
    idx->sai = (uint32_t *) calloc(idx->n_var, sizeof(uint32_t));
    idx->refseq = (uint8_t *) calloc(idx->n_var, sizeof(uint8_t));
    fread(idx->sai, sizeof(uint32_t), idx->n_var, fp_varseq);
    fread(idx->refseq, sizeof(uint8_t), idx->n_var, fp_varseq);
    err_fclose(fp_varseq);
    
    /* restore approximate index 
    FILE *fp;   
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".repeat");

    idx->is_multiseeds = bp_restore(fn0, 1);
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

    /* restore exntend index left 16 and right 16 */
    return idx;
}
void idx_destroy(idx_t *idx)
{
    if(idx->fastmap != NULL ) lkt_destroy(idx->fastmap);
    bwt_destroy(idx->bwt);
    bns_destroy(idx->bns);
    if(idx->pac != NULL ) free(idx->pac);
    free(idx->pmap);
    free(idx->refseq);
    free(idx->sai);

    bp_destroy(idx->is_multiseeds);
    free(idx->isa2seq16);
    free(idx->seq16s);
    //bp_destroy(idx->is_var);
    //free(idx->var_seqs);
    //variants_destroy(idx->var);
    free(idx);
}

int cut_ext_seq(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, vec_uint_t *seq16s, bp_t *is_repeat)
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


static void gen_filename(char *str, const char *prefix, const char *suffix)
{
    strncpy(str, prefix, MAX_NAME);strcat(str, suffix);

}
#define MAX_COUNT 800000
int idx_build_core(const char *fn_fa, const char *fn_av, const char *prefix)
{
    double t; 

    char fn0[MAX_NAME], fn1[MAX_NAME];
    /* fa2pac */
    t = realtime();
    gzFile fp_fa;
    fp_fa = gzopen(fn_fa, "r");
    fprintf(stderr, "\n[variant_index]: convert genome %s  to packedfile!\n",fn_fa); 
    bns_fasta2bntseq(fp_fa, prefix);
    gzclose(fp_fa);
    fprintf(stderr, "fa2pac %.2lf seconds elapse.\n", realtime() - t);

    /* pac2bwt */
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".pac");
    strncpy(fn1, prefix, MAX_NAME);strcat(fn1, ".bwt");
    t = realtime(); 
    fprintf(stderr, "\n[va_index]: convert pac %s to bwt %s\n", fn0, fn1);
    bwt_bwtgen(fn0, fn1); 
    fprintf(stderr, "pac2bwt %.2lf seconds elapse.\n", realtime() - t);
 

    {   /* update bwt */       
        strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".bwt");
        bwt_t *bwt;
        t = realtime();
        fprintf(stderr, "[va_index]: Update BWT... \n");
        bwt = bwt_restore_bwt(fn0);
        bwt_bwtupdate_core(bwt);
        bwt_dump_bwt(fn0, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "update bwt %.2f sec\n", (float)(realtime() - t));
    }
    {   /* lookup table */
        t = realtime();
        fprintf(stderr, "[va_index]: generate look up table... \n");
        char fn_lkt[MAX_NAME];
        strncpy(fn_lkt, prefix, MAX_NAME);strcat(fn_lkt, ".lkt");
        lkt_t *lkt = lkt_init(12);
        strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".pac");
        lkt_build_lookuptable(fn0, lkt); 
        lkt_dump(fn_lkt, lkt);
        lkt_destroy(lkt); 
        fprintf(stderr, "generate lookup table %.2f sec\n", (float)(realtime() - t));
    }
    {   /* cal sa and generate isa */
        bwt_t *bwt;
        t = realtime();
        fprintf(stderr, "[va_index]: generate SA and ISA from BWT\n");
        strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".bwt");
        bwt = bwt_restore_bwt(fn0);
        bwt_cal_sa(bwt, SA_INTV);
        bwt_cal_isa(bwt);
        strcat(strncpy(fn0, prefix, MAX_NAME), ".sa");
        bwt_dump_sa(fn0, bwt);
        strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".isa");
        bwt_dump_isa(fn0, bwt);
        fprintf(stderr, "cal SA and ISA %.2f sec\n", (float)(realtime() - t));
       
        /* generate snp-aware index 
         *1. generate is_variant array
         *2. generate vairiant preseq seqs 
         *3. 
         *
         *         |  varseq |seq|
         *         |preseq|orign|seq|
         * */

        /*generate is_variant map */
        fprintf(stderr, "[va_index]: build variant map\n");
        t = realtime();        
        bwtint_t *isa, i;        
        isa= bwt->isa;
        all_var_t *av = variants_build(fn_av);    
        var_t *v = av->var;
        bp_t *is_variant = bp_init(bwt->seq_len+1);/* if is_variant[i] == 1, there is a variant before sa[i] */       

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
        uint32_t *varseq;/* variant pre sequence*/
        
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
            varseq[vi*2] = bns_extract_seq16(pac, s, e);/* variant preseq include variants */
            varseq[vi*2+1] = i; /* bwt index */
        }
        
        uint32_t preseq, k;
        //uint32_t *preseq = (uint32_t *)calloc(av->n_var, sizeof(uint32_t)); 
        uint8_t *origin = (uint8_t *)calloc(av->n_var, sizeof(uint32_t)); 
        uint32_t *sai = (uint32_t *)calloc(av->n_var, sizeof(uint32_t)); 
        uint32_t *preseq_cumulate = (uint32_t *)calloc(HASH_SIZE, sizeof(uint32_t));/*preseq_cumulate[preseq] is subscript of sai and origin */

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
    {
        global_stat_36seed = (int *)calloc(MAX_COUNT, sizeof(int));
        global_stat_20seed0 = (int *)calloc(MAX_COUNT, sizeof(int));
        global_stat_20seed1 = (int *)calloc(MAX_COUNT, sizeof(int));
        /* 
         *
         *
         * */
       
        t = realtime();
        
        /* reload FM-index and packed reference */
        fprintf(stderr, "reload FM-index..\n");
        int64_t i;
        bwt_t *bwt;
        bntseq_t *bntseq; uint8_t *pac;
        bwtint_t ed, st, k, l, pos, isa;

        
        gen_filename(fn0, prefix, ".bwt"); 
        bwt = bwt_restore_bwt(fn0);   
        gen_filename(fn0, prefix, ".sa"); 
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
        for(ed =bwt->seq_len-1; ed>=L_MID_SEED-1; --ed){
            if(bp_get(ref_is_repeat, ed) ==1 ) continue;//skip repeat seed if not 1st occurence
            st = ed+1-L_MID_SEED;
            for(i = 0; i < L_MID_SEED; ++i){ seed[i] = __get_pac(pac, st+i);}
            k = 0, l = bwt->seq_len;
            if(bwt_match_exact_alt(bwt, L_MID_SEED, seed, &k, &l) <=1) continue;//skip uniq seeds 
            //log_bt2nt(20, seed);
            //printf("%u %u\n", k, l);
            for(isa = k; isa <=l; ++isa){
                pos = bwt_sa(bwt, isa); 
                bp_set1(ref_is_repeat, pos+L_MID_SEED-1);//set end pos of seed   
                bwt_is_repeat[isa] = isa==k?2:1;
            }
            ++n_multi_seeds;
        }
        fprintf(stderr, "%.2f sec\n", realtime()-t); 
        fprintf(stderr, "gen uniq 16mer\n");
        
        char fn_is_repeat[MAX_NAME], fn_preseq16[MAX_NAME], fn_isa2seq16[MAX_NAME];
        strncpy(fn_is_repeat, prefix, MAX_NAME);strcat(fn_is_repeat, ".repeat");
        strncpy(fn_isa2seq16, prefix, MAX_NAME);strcat(fn_isa2seq16, ".idxseq16");
        strncpy(fn_preseq16, prefix, MAX_NAME);strcat(fn_preseq16, ".seq16");
        
        FILE *fp_preseq16 = xopen(fn_preseq16, "w");

        /* fetch repeat seeds from bwt_is_repeat and generate uniq preseq  */
        vec_uint_t isa2seq16;//bwt index to seq16 index
        vec_uint_t seq16s;
        kv_init(isa2seq16);
        isa2seq16.a = (uint32_t *)calloc(n_multi_seeds+2, sizeof(uint32_t));
        isa2seq16.m = n_multi_seeds +2;
        kv_init(seq16s);
        seq16s.a = (uint32_t *)calloc(MAX_COUNT, sizeof(uint32_t));        
        seq16s.m = MAX_COUNT;
        kv_push(uint32_t, isa2seq16, 0);
        
        int n_seqs = 0;
        k = 0; l = 0;//extend seq [k, l)
        bp_t *bitmap = (bp_t *) bp_init(0x100000000); //bitmap to rm repeat preseq       
        while(k < bwt->seq_len){
            if(bwt_is_repeat[k] ==0 ){
                ++k;
                continue;
            } else if(bwt_is_repeat[k] == 2){
                l = k + 1;
                while(l < bwt->seq_len && bwt_is_repeat[l]==1 ){++l;}
                seq16s.n = 0;
                n_seqs += cut_ext_seq(pac, bwt, k, l,  &seq16s, bitmap);
                fwrite(seq16s.a, sizeof(uint32_t), seq16s.n, fp_preseq16);
                kv_push(uint32_t, isa2seq16, n_seqs);
                k = l;
            } else{/* should't be here! */}
        } 
   
        fprintf(stderr, "%.2f sec\n", realtime()-t); 
        bp_destroy(bitmap);
        /* use bwt_is_repeat to generate bwt_is_repeat2*/
        bp_t *bwt_is_repeat2 = bp_init(bwt->seq_len+1); /* bwt_is_repeat2[i] if Ref[sa[i]:sa[i]+20] is a multi seed*/

        for(st= 0; st < bwt->seq_len; ++st){
            if(bwt_is_repeat[st] == 2) bp_set1(bwt_is_repeat2, st);
        }
        bp_gen_occ(bwt_is_repeat2);
         
        /* 
         * dump bwt_is_repeat2
         * dump seq16s
         * dump isa2seq16
         * */
        
        fclose(fp_preseq16);     
        
        gen_filename(fn0, prefix, ".repeat"); 
        bp_dump(bwt_is_repeat2, fn0);
        FILE *fp;
        gen_filename(fn0, prefix, ".idxseq16");
        fp = xopen(fn0, "w");
        fwrite(&isa2seq16.n, sizeof(int), 1, fp); 
        fwrite(isa2seq16.a, sizeof(uint32_t), isa2seq16.n, fp); 
        fclose(fp); 
        
        kv_destroy(seq16s);
        kv_destroy(isa2seq16); 
        bp_destroy(bwt_is_repeat2);
        free(bwt_is_repeat);
        bp_destroy(ref_is_repeat);
        bwt_destroy(bwt);
        bns_destroy(bntseq);
        free(pac);
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


        free(global_stat_36seed);
        free(global_stat_20seed0);
        free(global_stat_20seed1);
    }
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
    for(i = 0; i < 0xFFF; ++i){
        log_seq162nt(i);
        int st = idx->pmap[i];
        int  ed = idx->pmap[i+1];
        for(j = st; j < ed; ++j){
            fprintf(stderr, "%u\t", idx->sai[j]); 
        }
        fprintf(stderr, "\n"); 
    }
    
    
    
    int n = 20;
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
    /* print preseq */
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
    //test_idx_reload(argv[optind+2]); 
    
    return 0;
}				/* ----------  end of function main  ---------- */


#ifdef MAIN_INDEX
int main(int argc, char *argv[]){

    return index_main(argc, argv);
}
#endif
