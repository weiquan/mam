/*
 * =====================================================================================
 *
 *       Filename:  index.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017年10月09日 17时23分22秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (wq), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef __INDEX_H
#define __INDEX_H


#include <stdint.h>
#include "utils.h"
//#include "variants.h"
#include "bwt.h"
#include "bntseq.h"
//#include "bitmap.h"
#include "lookup.h"
#include "rbwt_alt.h"
typedef struct{
    uint32_t seq, idx;
} ext_t;
typedef struct{
  int n;
  int n_L, n_R, n_relat;
  int relat_bg;//unused?
  uint32_t bgn, num;
  uint32_t *L2rel;
  uint32_t *R2rel;
  uint32_t *relat;
  uint32_t *nxt_idx;
  uint8_t *nxt_flg;
  rbwt_t *bwt0, *bwt1;
  //setbwt0, setbwt1
  //relation array
} hidx_t;
typedef struct{
  //global FM-index
  uint32_t fastmap_correct[12];
  lkt_t *fastmap;/*fastmap 12bp seq*/
  bwt_t *bwt;/*FM index*/
  bntseq_t *bns;
  uint8_t *pac;/* packed reference seq */
  //
  int n_rot;
  uint32_t *cnt_table; 
  //rbwt_t *rbwt; 
  
  //local FM-index
  int n_hier[5];
  hidx_t *hier_idx[5];  
  //sai to local FM-index, only used in generation 0
  int n_jmp, n_mod, max_buf;
  uint32_t *jmp_idx;  
  uint8_t *mod_idx; 
 
  //matrix for extending
  uint8_t mat[25]; 
  /* snp-aware seed index */
  int n_pmap;
  uint32_t *pmap;/* preseq map */ 
  int n_var;
  uint8_t *refseq;/*2bp refseq at variant site*/
  uint32_t *sai;/* sa index */
  
  /* approximate match */
  //bp_t *is_multiseeds;/*is_multiseeds[i] = 1 if Reference[sa[i]:sa[i]+20] is multi location seeds */
  uint32_t n_isa2seq16;
  uint32_t* isa2seq16; /* isa2seq16[i] = start of seq16s, i is occ of 1 in is_multiseeds before isa*/
  uint32_t n_seq16s;
  uint32_t *seq16s;

  int n_tot;
  uint32_t *lext_idx;
  uint32_t *lext1_idx;
  uint32_t *rext_idx; 
  
  int n_rext, n_lext0, n_lext1;
  ext_t *rext0, *lext0;
  uint32_t *lext1;
  /* 
  bp_t *is_var;// variants map
  uint32_t n_var;
  uint32_t *var_seqs;// var_seqs[2i] = variant pre seq; var_seqs[2i+1] = variant index 
  all_var_t *var;*/



} idx_t;
#define IS_SMLSIZ 50
#define LEN_INIT_SEED 20
#define LEN_EXT_SEQ 16
#define MAX_EXT 5

idx_t *fbwt_fmidx_restore(const char *prefix);
void fbwt_fmidx_destroy(idx_t *idx);
idx_t *fbwt_hier_idx_restore(idx_t *idx, const char *prefix);
void fbwt_hier_idx_destroy(idx_t *idx);
uint32_t aln_mam(int l_seq, uint8_t *seq, int s_bg, int s_ed, bwtint_t *bg, bwtint_t *ed, idx_t *idx, int max_loc);
uint32_t aln_mem(idx_t *idx, int l_seq, uint8_t *seq, int *s_k, int *s_l, bwtint_t *bg, bwtint_t *ed, int max_loc);
#endif
