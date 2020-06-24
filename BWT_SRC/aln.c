/*
 * =====================================================================================
 *
 *       Filename:  aln.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/20/2020 09:38:57 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <getopt.h>
#include <stdio.h>
#include "index.h"



int aln_main ( int argc, char *argv[] )
{
    int c;
    int k = 20;
    while((c = getopt(argc, argv, "k:h"))>=0){
        switch(c){
            case 'h':
                //return index_usage();
                break;
            case 'k':
                k = atoi(optarg);
                break; 
            default:
         return 1;
              
        }
    
    }
    if (argc - optind  != 1 ){
        //return index_usage();
      fprintf(stderr, "cmd is wrong!!!");
      exit(1);
    }
 
    
    
    
    aln_core(argv[optind]);

    
    return 0;
}				/* ----------  end of function main  ---------- */

int check_hier_idx(idx_t *idx)
{
  int i, j, k, x, n;
  int n_ext = 5;
  //idx->hier_idx = (hidx_t **)calloc(n_ext, sizeof(hidx_t *)); 
  n = 0;
  for(i = 0; i < idx->n_jmp; ++i) {
    for(j = idx->jmp_idx[i]; j < idx->jmp_idx[i+1]; ++j){
      fprintf(stderr, "[jmp %u]: bg = %u\n", n++, i*256+idx->mod_idx[j]); 
    } 
  
  }
  
  
  for(i = 0; i < n_ext; ++i) {
    fprintf(stderr, "gen = %u, n_hier = %u\n", i, idx->n_hier[i]);
    hidx_t *hidx = idx->hier_idx[i];
    for(j = 0; j < idx->n_hier[i]; ++j) {
      fprintf(stderr, "j = %u, bgn = %u, num = %u\n", j, hidx[j].bgn, hidx[j].num);
      fprintf(stderr, "n_L = %d, n_R = %d, n_rel = %d\n", hidx[j].n_L, hidx[j].n_R, hidx[j].n_relat);
      for(k =0 ; k < hidx[j].n_relat; ++k) {
        fprintf(stderr, "relat = %d, nxt_idx = %d, nxt_flg = %d\n", hidx[j].relat[k], hidx[j].nxt_idx[k], hidx[j].nxt_flg[k]);
      
      }

    }

  }
  fprintf(stderr, "check finish!\n");
  return EXIT_SUCCESS;
}



int check_exact_match(idx_t *idx)
{
  int i, j, k, x;
  int n_ext = 5;
  //idx->hier_idx = (hidx_t **)calloc(n_ext, sizeof(hidx_t *)); 
  for(i = 0; i < n_ext; ++i) {
    fprintf(stderr, "gen = %u, n_hier = %u\n", i, idx->n_hier[i]);
    hidx_t *hidx = idx->hier_idx[i];
    for(j = 0; j < idx->n_hier[i]; ++j) {
      fprintf(stderr, "j = %u, bgn = %u, num = %u\n", j, hidx[j].bgn, hidx[j].num);
      fprintf(stderr, "n_L = %d, n_R = %d, n_rel = %d\n", hidx[j].n_L, hidx[j].n_R, hidx[j].n_relat);
      for(k =0 ; k < hidx[j].n_relat; ++k) {
        fprintf(stderr, "relat = %d, nxt_idx = %d, nxt_flg = %d\n", hidx[j].relat[k], hidx[j].nxt_idx[k], hidx[j].nxt_flg[k]);
      }

    }

  }
  fprintf(stderr, "check finish!\n");
  return EXIT_SUCCESS;
}



/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
int aln_core(char *prefix)
{
  idx_t idx;
  fbwt_fmidx_restore(&idx, prefix); 
  fbwt_hier_idx_restore(&idx, prefix);
  //check_hier_idx(&idx);
  idx.cnt_table = calloc(256, sizeof(uint32_t));
  rbwt_gen_cnt_table(idx.cnt_table); 
  bwtint_t i, j;
  for(i = 200; i < idx.bwt->seq_len-512; ++i){
    uint8_t seq[256] = {}; 
    for(j = 0; j < 256; ++j) {
      seq[j] = __get_pac(idx.pac, i+j); 
    }

    uint32_t init_bg = 16*5; 
    uint32_t init_ed = init_bg + 20; 
    uint32_t k = 0, l = idx.bwt->seq_len;
    bwt_match_exact_alt(idx.bwt, init_ed-init_bg, seq+init_bg, &k, &l);
    fprintf(stderr, "[%s:%u]: pos = %u, init = [%u, %u) (k, l)= (%u, %u)\n", __func__, __LINE__, i, init_bg, init_ed, k, l); 
    aln_mem(256, seq, &init_bg, &init_ed, &k, &l, &idx, IS_SMLSIZ);       
    fprintf(stderr, "[%s:%u]: pos = %u, init = [%u, %u) (k, l)= (%u, %u)\n", __func__, __LINE__, i, init_bg, init_ed, k, l); 
    uint32_t k1 = 0, l1 = idx.bwt->seq_len;   
    bwt_match_exact_alt(idx.bwt, init_ed-init_bg, seq+init_bg, &k1, &l1);
    if(k1 != k || l1 != l) {
      fprintf(stderr, "[%s:%u]: pos = %u, init = [%u, %u) (k, l)= (%u, %u) != (%u, %u)!\n", __func__, __LINE__, i, init_bg, init_ed, k, l, k1, l1); 
      exit(1); 
    }




  }  
  
  return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
#ifdef MAIN_ALN
int main(int argc, char *argv[]){

    return aln_main(argc, argv);
}
#endif
