/*
 * =====================================================================================
 *
 *       Filename:  index.c
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
#include <getopt.h>
#include "utils.h"
#include "khash.h"
#include "index.h"
#include "kvec.h"
#include "ksort.h"
#include "setFileName.h"


#define USE_MALLOC_WRAPPERS
#define MAX_NAME 128

#define FLAG_UNIQ 0
#define FLAG_REP_ST 2
#define FLAG_REP 1



const int SA_INTV = 4;
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
    fprintf(stderr, "index\t[OPT]\t<ref.fa>\t<prefix>\n\n"); 
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

idx_t *fbwt_fmidx_restore(const char *prefix)
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
  return idx;
}
idx_t *fbwt_hier_idx_restore(idx_t *idx, const char *prefix)
{
  int i, j, x;
  int n_ext = 5;
  //idx->hier_idx = (hidx_t **)calloc(n_ext, sizeof(hidx_t *)); 
  char fn[128];
  //restore sai to local index
  sprintf(fn, "%s_jmpmod.txt", prefix);
  FILE *fp_jmp = xopen(fn, "rb"); 
  fread(&idx->n_jmp, sizeof(uint32_t), 1, fp_jmp);
  fread(&idx->n_mod, sizeof(uint32_t), 1, fp_jmp);
  fread(&idx->max_buf, sizeof(uint32_t), 1, fp_jmp);
  idx->jmp_idx = calloc(idx->n_jmp, sizeof(uint32_t));
  idx->mod_idx = calloc(idx->n_mod, sizeof(uint8_t));
  fread(idx->jmp_idx, sizeof(uint32_t), idx->n_jmp, fp_jmp);
  fread(idx->mod_idx, sizeof(uint8_t) , idx->n_mod, fp_jmp);  
  err_fclose(fp_jmp);
  for(i = 0; i < n_ext; ++i) {
    fprintf(stderr, "[%s] gen hier%u idx!\n", __func__, i);
    sprintf(fn, "%s_relat_%u.txt", prefix, i);
    FILE *fp_relat = xopen(fn, "rb"); 
    sprintf(fn, "%s_smbwt_%u.txt", prefix, i);
    FILE *fp_smbwt = xopen(fn, "rb"); 

    int n_local_idx = 0, n;
    fread(&n_local_idx, sizeof(uint32_t), 1, fp_relat);
    idx->n_hier[i] = n_local_idx; 
    idx->hier_idx[i] = calloc(n_local_idx, sizeof(hidx_t)); 
    hidx_t *hidx = idx->hier_idx[i];
    for(j = 0; j < n_local_idx; ++j) {
      if(j % 100000 == 0)fprintf(stderr, "[%s] %u local idx has been reload!\n", __func__, j);
      fread(&hidx[j].bgn, sizeof(uint32_t), 1, fp_relat);
      fread(&hidx[j].num, sizeof(uint32_t), 1, fp_relat);
      //fprintf(stderr, "bgn = %u, num = %u\n", hidx[j].bgn, hidx[j].num);
      fread(&n, sizeof(uint32_t), 1, fp_relat);
      //fprintf(stderr, "n_L = %u\n", n);
      hidx[j].n_L = n;
      hidx[j].L2rel = calloc(n, sizeof(uint32_t));
      fread(hidx[j].L2rel,sizeof(uint32_t), n, fp_relat);
      
      fread(&n, sizeof(uint32_t), 1, fp_relat);
      //fprintf(stderr, "n_R = %u\n", n);
      hidx[j].n_R = n;
     /*  
      hidx[j].R2rel = calloc(n, sizeof(uint32_t));
      fread(hidx[j].R2rel,sizeof(uint32_t), n, fp_relat);
      free(hidx[j].R2rel);
    */
      fread(&n, sizeof(uint32_t), 1, fp_relat);
      hidx[j].relat_bg = n;
      fread(&n, sizeof(uint32_t), 1, fp_relat);
      hidx[j].n_relat = n;
      //fprintf(stderr, "n_relat = %u\n", n);
      hidx[j].relat = calloc(n, sizeof(uint32_t));
      fread(hidx[j].relat, sizeof(uint32_t), n, fp_relat);
      hidx[j].nxt_idx = calloc(n, sizeof(uint32_t));
      fread(hidx[j].nxt_idx, sizeof(uint32_t), n, fp_relat);
      hidx[j].nxt_flg = calloc(n, sizeof(uint8_t));
      fread(hidx[j].nxt_flg, sizeof(uint8_t), n, fp_relat);
    }

    for(j = 0; j < n_local_idx; ++j) {
      hidx[j].bwt0 = rbwt_bwt_restore(fp_smbwt);
      hidx[j].bwt1 = rbwt_bwt_restore(fp_smbwt);
    }
    err_fclose(fp_relat);
    err_fclose(fp_smbwt);
  } 
  idx->cnt_table = calloc(256, sizeof(uint32_t));
  rbwt_gen_cnt_table(idx->cnt_table); 
  // initialize scoring matrix
	int sa = 1, sb = 2;
  int k;
  uint8_t *mat = idx->mat;
  for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? sa : -sb;
		mat[k++] = 0; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = 0;
  gen_hash_boundry(idx->fastmap_correct, idx->bwt, idx->pac);
  return idx;
}
void fbwt_hier_idx_destroy(idx_t *idx)
{
  int i, j, x;
  int n_ext = 5;
  //idx->hier_idx = (hidx_t **)calloc(n_ext, sizeof(hidx_t *)); 
  free(idx->jmp_idx);
  free(idx->mod_idx);
  free(idx->cnt_table);
  for(i = 0; i < n_ext; ++i) {
    hidx_t *hidx = idx->hier_idx[i];
    for(j = 0; j < idx->n_hier[i]; ++j) {
      free(hidx[j].L2rel);
      //free(hidx[j].R2rel);
      free(hidx[j].relat);
      free(hidx[j].nxt_idx);
      free(hidx[j].nxt_flg);
    }
    free(hidx);
  } 
  
}




void fbwt_fmidx_destroy(idx_t *idx)
{
    if(idx->fastmap != NULL ) lkt_destroy(idx->fastmap);
    bwt_destroy(idx->bwt);
    bns_destroy(idx->bns);
    if(idx->pac != NULL ) free(idx->pac);
    free(idx);
}

#define MASK_2NT 0xF
#define MASK_6NT 0xFFF
#define L_MID_SEED 20
const int HASH_SIZE = 4097;


static inline void cat_filename(char *str, const char *prefix, const char *suffix)
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
    fprintf(stderr, "\n[fm-idx]:  convert genome %s  to packedfile!\n",fn_fa); 
    bns_fasta2bntseq(fp_fa, prefix);
    fprintf(stderr, "[fm-idx]:  fa2pac %.2lf seconds elapse.\n", realtime() - t);
    gzclose(fp_fa);   
    
    t = realtime(); 

    strcat(strncpy(fn0, prefix, MAX_NAME), ".pac");
    strcat(strncpy(fn1, prefix, MAX_NAME), ".bwt");
    fprintf(stderr, "\n[fm-idx]:  convert pac %s to bwt %s\n", fn0, fn1);

    bwt_bwtgen(fn0, fn1); 
    fprintf(stderr, "[fm-idx]:  pac2bwt %.2lf seconds elapse.\n", realtime() - t);
    {    
        bwt_t *bwt;
        t = realtime();
        strcat(strncpy(fn0, prefix, MAX_NAME), ".bwt");
        fprintf(stderr, "[fm-idx]:  Update BWT... \n");
        bwt = bwt_restore_bwt(fn0);
        bwt_bwtupdate_core(bwt);
        bwt_dump_bwt(fn0, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "[fm-idx]:  update bwt %.2f sec\n", (float)(realtime() - t));
    }
    {   
        char fn_lkt[MAX_NAME];
        t = realtime();
        fprintf(stderr, "[fm-idx]:  generate look up table... \n");
        strcat(strncpy(fn_lkt, prefix, MAX_NAME), ".lkt");
        lkt_t *lkt = lkt_init(12);
        strcat(strncpy(fn0, prefix, MAX_NAME), ".pac");
        lkt_build_lookuptable(fn0, lkt); 
        lkt_dump(fn_lkt, lkt);
        lkt_destroy(lkt); 
        fprintf(stderr, "[fm-idx]:  generate lookup table %.2f sec\n", (float)(realtime() - t));
    }

    {   
        bwt_t *bwt;
        t = realtime();
        fprintf(stderr, "[fm-idx]:  generate SA and ISA from BWT\n");
        strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".bwt");
        bwt = bwt_restore_bwt(fn0);
        fprintf(stderr, "[fm-idx]:  call SA\n");
        bwt_cal_sa(bwt, SA_INTV);
        fprintf(stderr, "[fm-idx]:  call ISA\n");
        bwt_cal_isa(bwt);
        strcat(strncpy(fn0, prefix, MAX_NAME), ".sa");
        fprintf(stderr, "[fm-idx]:  dump SA\n");
        bwt_dump_sa(fn0, bwt);
        strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".isa");
        bwt_dump_isa(fn0, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "[fm-idx]:  cal SA and ISA %.2f sec\n", (float)(realtime() - t));
    }

}
void build_sa(const char *fn_fa, const char *prefix)
{
    double t;
    char fn0[MAX_NAME], fn1[MAX_NAME];

    gzFile fp_fa;
    t = realtime();
    fp_fa = gzopen(fn_fa, "r");
    fprintf(stderr, "\n[fm-idx]:  convert genome %s  to packedfile!\n",fn_fa); 

    fprintf(stderr, "[fm-idx]:  fa2pac %.2lf seconds elapse.\n", realtime() - t);
    gzclose(fp_fa);   
    
    t = realtime(); 

    strcat(strncpy(fn0, prefix, MAX_NAME), ".pac");
    strcat(strncpy(fn1, prefix, MAX_NAME), ".bwt");
    fprintf(stderr, "\n[fm-idx]:  convert pac %s to bwt %s\n", fn0, fn1);

    
    

    {   
        bwt_t *bwt;
        t = realtime();
        fprintf(stderr, "[fm-idx]:  generate SA and ISA from BWT\n");
        strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".bwt");
        bwt = bwt_restore_bwt(fn0);
        fprintf(stderr, "[fm-idx]:  call SA\n");
        bwt_cal_sa(bwt, SA_INTV);
        strcat(strncpy(fn0, prefix, MAX_NAME), ".sa");
        fprintf(stderr, "[fm-idx]:  dump SA\n");
        bwt_dump_sa(fn0, bwt);
        bwt_destroy(bwt);
        fprintf(stderr, "[fm-idx]:  cal SA and ISA %.2f sec\n", (float)(realtime() - t));
    }

}


#define N_EXT 5
#define LEN_NEXT_IDX 1024
#define LEN_SEED 20
#define LEN_EXT 16
#define IDX_MID_SIZE (254*256)

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
//generate bwt idx for every repetive seeds of i-th extension 
void gen_seedidx(uint8_t *bwt_is_repeat, idx_t *idx, int ext_i, const char*prefix)
{
    //__check_bwt_is_repeat(bwt_is_repeat, ext_i, idx);
    bwt_t *bwt = idx->bwt;
    char fn[1024]={};
    sprintf(fn, "%s_seedidx_%u.txt", prefix, ext_i);
    FILE *fp_seedidx = xopen(fn, "wb");
    int is_rep = 0;
    int last_flag = 0, cur_flag = 0;
    uint32_t bg_idx, ed_idx, tot_idx; 
    uint32_t next_idx[3] = {0,0,0};
    uint32_t glb_idx = 0;   
    int n_seedidx = 0;
    while(glb_idx <= bwt->seq_len+1) {
      cur_flag = bwt_is_repeat[glb_idx];  
      /*  
      if(last_flag == FLAG_UNIQ){
        if(cur_flag == FLAG_REP_ST){
          bg_idx = glb_idx;
          num_idx = 1;
          //last_flag = 2; 
        } else if(cur_flag == FLAG_REP){
          fprintf(stderr, "bg_idx = %u, ed_idx = %u\n", bg_idx, ed_idx);
          fprintf(stderr, "glb_idx = %u, flag = %u %u\n", glb_idx, bwt_is_repeat[glb_idx], bwt_is_repeat[glb_idx-1]);
          fprintf(stderr, "%u, %s, should't be here!\n", __LINE__, __func__);
          exit(1);
        }
      } else if(last_flag ==FLAG_REP){
        if(cur_flag == FLAG_UNIQ){
          ed_idx = glb_idx-1;
          tot_idx = num_idx;
          is_rep = 1;
        } else if(cur_flag==FLAG_REP){
          ++num_idx;
        } else if(cur_flag ==FLAG_REP_ST){
          tot_idx = num_idx;
          is_rep = 1;
          ed_idx = glb_idx-1;
          //bg_idx = glb_idx;
          num_idx = 1;
          //last_flag = 2;
        }
      } else if(last_flag ==FLAG_REP_ST){
        if(cur_flag ==FLAG_UNIQ){
          fprintf(stderr, "bg_idx = %u, ed_idx = %u\n", bg_idx, ed_idx);
          fprintf(stderr, "glb_idx = %u, flag = %u %u\n", glb_idx, bwt_is_repeat[glb_idx], bwt_is_repeat[glb_idx-1]);
          fprintf(stderr, "%u, %s, should't be here!\n", __LINE__, __func__);
          exit(1);
        } else if(cur_flag ==FLAG_REP){
          ++num_idx;
        } else if(cur_flag == FLAG_REP_ST){
          fprintf(stderr, "bg_idx = %u, ed_idx = %u\n", bg_idx, ed_idx);
          fprintf(stderr, "glb_idx = %u, flag = %u %u\n", glb_idx, bwt_is_repeat[glb_idx], bwt_is_repeat[glb_idx-1]);
          fprintf(stderr, "%u, %s, should't be here!\n", __LINE__, __func__);
          exit(1);
        }
      } 
      last_flag = cur_flag;
      ++glb_idx;
      */
      if(cur_flag == FLAG_UNIQ) { 
        bg_idx = glb_idx;
        ed_idx = glb_idx+1;  
      } else if(cur_flag == FLAG_REP_ST) {
        bg_idx = glb_idx;
        ed_idx = glb_idx+1;  
        int next_flag = ed_idx <= bwt->seq_len+1?bwt_is_repeat[ed_idx]:FLAG_UNIQ; 
        while(ed_idx <= bwt->seq_len+1 &&  next_flag == FLAG_REP){
          ed_idx++; 
          next_flag = ed_idx <= bwt->seq_len+1?bwt_is_repeat[ed_idx]:FLAG_UNIQ; 
        } 
      } else{
          fprintf(stderr, "[%s]:  cur_flg = %u, glb_idx = %u\n", __func__,cur_flag, glb_idx);
      
      }
      tot_idx = ed_idx  - bg_idx;
      glb_idx = ed_idx;

      if(cur_flag == FLAG_UNIQ) continue; 
      //print begin index and total number of index
#ifdef DEBUG
      if(!__check_seedidx(bg_idx, tot_idx, 20+32*ext_i, idx)){
        fprintf(stderr, "[%s]:  bg_idx = %u, ed_idx = %u, tot_idx = %u, cur_flag = %u\n", __func__, bg_idx, ed_idx, tot_idx, cur_flag);
        bwtint_t i;
        for(i = bg_idx; i < ed_idx; ++i) {
            fprintf(stderr, "[%s]:  idx = %u, flag = %u\n", __func__, i, bwt_is_repeat[i]);
        }
        exit(1);
      }
#endif
      if(n_seedidx % 10000 == 0) fprintf(stderr, "[%s]: n_seedidx = %d\n", __func__, n_seedidx);
      n_seedidx++;
      next_idx[0] = bg_idx;
      next_idx[1] = tot_idx;
      fwrite(next_idx, sizeof(uint32_t), 2, fp_seedidx);
    }
    fclose(fp_seedidx); 
    return;
}

int gen_rank(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, uint32_t rank[])
{
  bwtint_t i,ii, j, pos, st, ed;

  uint32_t max_pos, min_pos;
  if(offset<0){
      printf("%s\n","Err: offset can not be less than 0 !");
      exit(1);
  }
  max_pos = bwt->seq_len-offset-8;        
  //max_pos = bwt->seq_len-20;        
  uint32_t last_seq8, seq8;
  uint32_t n = 0;
  int wrong_pos = 0;
  
  rank[0] = n; 
  pos = bwt_sa(bwt, k);       
  if(pos <= max_pos){
    st = pos + offset; 
    seq8 = 0;
    for(ii = st; ii < st+8; ++ii) {
      seq8 |= __get_pac(pac, ii);
      seq8 <<= 2;
    }
    last_seq8 = seq8;
    wrong_pos = 0;
  }else{ // pos > reference length
    last_seq8 = 0; 
    wrong_pos = 1;
  }
  for(i = k+1; i <= l; ++i){
    j = i-k;
    pos = bwt_sa(bwt, i);       
    if(pos <= max_pos){
      st = pos + offset; 
      seq8 = 0;
      for(ii = st; ii < st+8; ++ii) {
        seq8 |= __get_pac(pac, ii);
        seq8 <<= 2;
      } 
      if(wrong_pos == 1 || last_seq8 < seq8){ ++n; }
      wrong_pos = 0;
    }else{
      ++n;
      wrong_pos = 1;
    }

    last_seq8 = seq8; 
    rank[j] = n;
  }
  rank[++j] = ++n;
  rank[++j] = ++n;
  return j;
}

int cut_ext_seq0(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, uint32_t sort_flag[])
{
    bwtint_t i,j,ii, pos, st, ed;

    uint32_t max_pos, min_pos;
    if(offset<0){
        printf("%s\n","Err: offset can not be less than 0 !");
        exit(1);
    }
    max_pos = bwt->seq_len-offset-8;        
    uint32_t last_seq16, seq16;
    uint32_t n = 0;
    int last_flag = 0, cur_flag = 0;
    sort_flag[0] = n; 
    pos = bwt_sa(bwt, k);       
    if(pos <= max_pos){
        st = pos + offset; 
        seq16 = 0;
        for(ii = st; ii < st+16; ++ii) {
          seq16 |= __get_pac(pac, ii);
          seq16 <<= 2;
        }
        last_seq16 = seq16;
        last_flag = 0;
    }else{
        last_seq16 = 0; 
        last_flag = 1;
    }
    for(i = k+1; i <= l; ++i){
        j = i-k;
        pos = bwt_sa(bwt, i);       
        if(pos <= max_pos){
            st = pos + offset; 
            seq16 = 0;
            for(ii = st; ii < st+16; ++ii) {
              seq16 |= __get_pac(pac, ii);
              seq16 <<= 2;
            } 
            cur_flag = 0;
        }else{
            ++n;
            cur_flag = 1;
        }
        if(cur_flag == 0){
            if(last_flag == 1 || last_seq16 < seq16){
                ++n;
            } 
        }
        last_flag = cur_flag;
        last_seq16 = seq16; 
        sort_flag[j] = n;
    }
    sort_flag[++j] = ++n;
    sort_flag[++j] = ++n;
    return j;
}
int cut_ext_seq1(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, uint32_t sort_flag[])
{
  bwtint_t i, j, ii, pos, st, ed;
  uint32_t max_pos, min_pos;
  if(offset<0){
    printf("%s\n","Err: offset can not be less than 0 !");
    exit(1);
  }
  max_pos = bwt->seq_len-offset-16;
  uint32_t last_seq16, seq16;
  uint32_t n = 0;
  int last_flag = 0, cur_flag = 0;
  sort_flag[0] = n;
  pos = bwt_sa(bwt, k);
  if(pos <= max_pos){
    st = pos + offset; 
    seq16 = 0;
    for(ii = st; ii < st+16; ++ii) {
      seq16 <<= 2;
      seq16 |= __get_pac(pac, ii);
    }
    last_seq16 = seq16;
    last_flag = 0;
  }else{
    //fprintf(stderr, "pos = %u, K = %u, n = %u\n", pos, k, n);
    last_seq16 = 0; 
    last_flag = 1;
  }
  for(i = k+1; i <= l; ++i){
    j = i-k;
    pos = bwt_sa(bwt, i);       
    if(pos <= max_pos){
      st = pos + offset;
      seq16 = 0;
      for(ii = st; ii < st+16; ++ii) {
        seq16 <<= 2;
        seq16 |= __get_pac(pac, ii);
      }
      cur_flag = 0;
    }else{
      cur_flag = 1;
      ++n;
    } 
    if(cur_flag == 0) {
      if(last_flag == 1 || last_seq16 < seq16) {
        n++;
      }
    }
    last_flag = cur_flag;
    last_seq16 = seq16; 
    sort_flag[j] += n;
  }
  sort_flag[++j] += ++n;
  sort_flag[++j] += ++n;
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


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  build_bwt_is_repeat
 *  Description:  
 * =====================================================================================
 */
int gen_hash_boundry(uint32_t hash_boundry[12], bwt_t *bwt, uint8_t *pac)
{
  uint32_t k = 0, l = 0, i, j;
  uint8_t seed[12] = {};
  for(i =0; i < 12; ++i){
    seed[i] = __get_pac(pac, bwt->seq_len-12+i);
  }

  for(i =0; i < 11; ++i){
    bwt_match_exact_alt(bwt, 1, seed+12-1-i, &k, &l);
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
  return 0;
}
void build_bwt_is_repeat ( const char *prefix, const idx_t *idx, const uint8_t *pac)
{
  int64_t i, j, ext_i;
  bwt_t *bwt = idx->bwt;
  lkt_t *fastmap = idx->fastmap; 
  uint8_t *bwt_is_repeat[2];
  bwt_is_repeat[0] = (uint8_t *)calloc(bwt->seq_len+2, sizeof(uint8_t));/* indicate repetive seed on bwt index */
  bwt_is_repeat[1] = (uint8_t *)calloc(bwt->seq_len+2, sizeof(uint8_t));/* indicate repetive seed on bwt index */ 
  uint32_t *hash_boundry = idx->fastmap_correct;
  gen_hash_boundry(hash_boundry, idx->bwt, idx->pac);
  uint32_t MAX_12_count = fastmap->item[0];
  for(i = 0; i < idx->fastmap->n_item+1; ++i){
    uint32_t count = fastmap->item[i+1] - fastmap->item[i];
    if(count > MAX_12_count) MAX_12_count = count;
  }
  

  uint32_t *rank = calloc(MAX_12_count+3, sizeof(uint32_t));
  uint8_t *bwt_repeat0 = bwt_is_repeat[0];
  uint8_t *bwt_repeat1 = bwt_is_repeat[1];

  //+++++++++++++++++++++++++++++++++++++++++++++++++++
  uint32_t bg_idx = 0, ed_idx = 0, glb_idx = 0; 
  uint32_t MAX_20_count = 0;
  for(glb_idx =0; glb_idx < fastmap->n_item; ++glb_idx){
    bg_idx = fastmap->item[glb_idx];
    ed_idx = fastmap->item[glb_idx+1];
    ed_idx -= get_12mer_correct(hash_boundry, ed_idx-1);
    if(bg_idx >= bwt->seq_len) break;
    uint32_t num_idx;
    if(bg_idx +1>= ed_idx) { 
      continue;
    } else{ 
      num_idx = ed_idx - bg_idx;
      ed_idx -= 1;
    }
    int offset = 12;
    gen_rank(pac, bwt, bg_idx, ed_idx, offset, rank);
    uint32_t i_st = 0; 
    for(i = 0; i < num_idx; ++i){
      if(rank[i] < rank[i+1]){
        if(i - i_st >= 1){
          bwt_repeat0[bg_idx+i_st] = FLAG_REP_ST; 
          for(j = i_st+1; j <= i; ++j) bwt_repeat0[bg_idx+j] = FLAG_REP;
          if(i-i_st > MAX_20_count) MAX_20_count = i-i_st;
        }
        i_st = i+1;                 
      }      
    }
  }
 
  /*
    fprintf(stderr, "for debug 12\n");
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
 

  gen_seedidx(bwt_repeat0, idx, 0, prefix); 
  fprintf(stderr, "dump 20 bp seed data!\n");

  for(ext_i =1; ext_i <= N_EXT; ++ext_i){
    uint32_t num_idx =0;
    glb_idx = 0;
    bg_idx = ed_idx= glb_idx;
    int n_interval = 0;
    while(glb_idx <= bwt->seq_len ){
      ++n_interval;
      while(glb_idx <=bwt->seq_len && bwt_repeat0[glb_idx] == FLAG_UNIQ){ glb_idx++;}
      bg_idx = glb_idx++; 
      while(glb_idx <=bwt->seq_len && bwt_repeat0[glb_idx] == FLAG_REP){ glb_idx++;}
      ed_idx = glb_idx-1;
      num_idx = ed_idx+1 - bg_idx;
      for(i = 0; i <= num_idx; ++i){
          rank[i] = 0; 
      }
      int offset = LEN_SEED+2*LEN_EXT*(ext_i-1);
      cut_ext_seq1(pac, bwt, bg_idx, ed_idx, offset, rank);
      cut_ext_seq1(pac, bwt, bg_idx, ed_idx, offset+16, rank);
      uint32_t i_st = 0; 
      for(i = 0; i < num_idx; ++i){

        fprintf(stderr, "[%s]: i = %u, i_st = %u, bg_idx+i_st = %u, rank[i] = %u, rank[i+1] = %u\n", __func__, i, i_st, bg_idx+i_st, rank[i], rank[i+1]);
        if(rank[i] < rank[i+1]){
          //ref_is_repeat[cur_pos] = 3;
          if(i - i_st >= 1){
            //fprintf(stderr, "[%s]: i = %u, i_st = %u, bg_idx+i_st = %u\n", __func__, i, i_st, bg_idx+i_st);
            bwt_repeat1[bg_idx+i_st] = FLAG_REP_ST; 
            for(j = i_st+1; j <= i; ++j) bwt_repeat1[bg_idx+j] = FLAG_REP;
          }
          i_st = i+1;
        }                 
      }
    }// end while(glb_idx < bwt->seq_len)
    
   
    fprintf(stderr, "ext = %u, len = %u ++++++++++++++\n", ext_i, 20+ext_i*32);
    
    //__check_bwt_is_repeat(bwt_repeat1, ext_i, idx);
    gen_seedidx(bwt_repeat1, idx, ext_i, prefix);// generate bg-idx and ed-idx for repetive seeds 
    fprintf(stderr, "dump data finish\n");
    SWAP(uint8_t *, bwt_repeat0, bwt_repeat1);
    memset(bwt_repeat1, FLAG_UNIQ, bwt->seq_len+2);
  } // end for
  
  free(rank);
  free(bwt_is_repeat[1]);
  free(bwt_is_repeat[0]);
  return;
}		/* -----  end of function build_bwt_is_repeat  ----- */
int __check_bwt_is_repeat(uint8_t *bwt_repeat1, int ext_i, idx_t *idx)
{

  fprintf(stderr, "[%s]:  check bwt repeat for generation %d\n", __func__, ext_i);
  bwt_t *bwt = idx->bwt;
  uint8_t *pac = idx->pac;
  int l_ext_seed = 20 + 32*ext_i;

  bwtint_t i = 0, j, k, l;
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
  return 1;

}

/*  
{
      fprintf(stderr, "for debug\n");

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

}
*/


void build_extend_idx_alt(const char *prefix, int init_k, int ext_k)
{
  char fn0[MAX_NAME];
  double t = realtime();
  /* reload FM-index and packed reference */
  fprintf(stderr, "[%s]:  reload FM-index..\n", __func__);
  idx_t *idx = fbwt_fmidx_restore(prefix);
  fprintf(stderr, "[%s]:  %.2f sec\n", __func__, realtime()-t); 
  /* generate repeat seeds and bwt range */
  fprintf(stderr, "[%s]:  gen repeat on bwt\n", __func__);
  //bp_t *ref_is_repeat; 
  //uint8_t *ref_is_repeat; 

  //ref_is_repeat = (bp_t *)bp_init(bwt->seq_len);/*indicate repeat seeds on reference*/
  //ref_is_repeat = (uint8_t *)calloc(bwt->seq_len+2, sizeof(uint8_t));/* indicate repeat seed on bwt index */
  build_bwt_is_repeat(prefix, idx, idx->pac);
  fbwt_fmidx_destroy(idx);
  return;
}
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
            ++n_dup_seeds;
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

#define  LEN_SEED  20 
#define  OFF_SEED  10 
#define  LEN_EXT   16 
#define  LEN_READ  150
#define  NUM_EXT  (((LEN_READ- LEN_SEED )/2)+OFF_SEED+ LEN_EXT-1)/LEN_EXT  
#define  MAX_SEED_NUM 600000
#define  LEN_FILE_NAME 100 

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
	uint32_t (*relat_sai)[2]; // sortR[]的每一个序列的(bgnIdx,ednIdx)
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



uint32_t getlen(uint32_t LenExtIdx[],FileName *f);

void CrtJmpMod(int n, uint32_t cur_seedidx[][3], uint32_t* jmp_idx, uint8_t *mod_idx, uint32_t *cap_pos);

void OpenSeedIdxfiles(FileName *f, uint32_t (*sidx)[3], uint32_t LenExtIdx[],uint32_t repNum,FILE *fpw[6]);

void comFile(FileName *rf, FileName *wf);
long getFileSize(char* file);
void getBlckData(uint32_t sidx[2],uint32_t num_ext, SubBuf *sub_buf, idx_t *idx);


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define MAX_BLCK_SIZE 500000
#define MAX_SEED_SIZE 500000 
#define ALL_BUFS_SIZE 12  // All bufs size in SubBuf





long getFileSize(char *fn)
{
	FILE *fp = xopen(fn, "r");
	fseek(fp, 0, SEEK_END);
	long len = ftell(fp);
	err_fclose(fp);
  return len; 
}

int fbwt_rbwt_build(idx_t *idx)
{
  idx->n_rot = LEN_EXT;
  idx->cnt_table = calloc(256, sizeof(uint32_t));
  rbwt_gen_cnt_table(idx->cnt_table);
}



int idx_build_core(const char *fn_fa, const char *prefix, int init_k, int ext_k)
{
  fprintf(stderr, "[%s]:  build fm-idx.\n", __func__);
  build_fmidx(fn_fa, prefix);
  build_sa(fn_fa, prefix);
  fprintf(stderr, "[%s]:  build extend idx.\n", __func__);
  build_extend_idx_alt(prefix, init_k, ext_k);
  fprintf(stderr, "[%s]:  build hier index.\n", __func__);
  build_hier_idx(prefix); 

  return 0;
}
void test_idx_reload(const char *prefix)
{
    bwtint_t i, j;
   
    idx_t *idx = fbwt_fmidx_restore(prefix);
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
    fbwt_fmidx_destroy(idx);
}
int index_main ( int argc, char *argv[] )
{
    int c;
    int k = 20;
    while((c = getopt(argc, argv, "k:h"))>=0){
        switch(c){
            case 'h':
                return index_usage();
            case 'k':
                k = atoi(optarg);
                break; 
            default:
         return 1;
              
        }
    
    }
    if (argc - optind  != 2 ){
        return index_usage();
    }
 
    
    
    
    idx_build_core(argv[optind], argv[optind+1], k, 16);

    
    return 0;
}				/* ----------  end of function main  ---------- */


#ifdef MAIN_INDEX
int main(int argc, char *argv[]){

    return index_main(argc, argv);
}
#endif
