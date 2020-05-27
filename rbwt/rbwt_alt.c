/*
 * =====================================================================================
 *
 *       Filename:  rbwt.c
 *
 *    Description:  rotation bwt 
 *
 *        Version:  1.0
 *        Created:  05/19/2020 05:14:23 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "ksort.h"
//#include "bwt.h"
#include "rbwt_alt.h"
#define USE_MALLOC_WRAPPERS
KSORT_INIT_GENERIC(uint32_t)
#define SWAP(type_t, a, b) do{type_t x=(a); (a)=(b); (b)=x;} while(0)
void rbwt_gen_cnt_table(uint32_t cnt_table[])
{
  int i, j;
  for (i = 0; i != 256; ++i) {
    uint32_t x = 0;
    for (j = 0; j != 4; ++j)
      x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
    cnt_table[i] = x;
  }
}


static inline void rbwt_accumulate_C(int n, int C[])
{
  int i;
  for(i = 1; i < n-1; ++i) C[i] += C[i-1]; 
  for(i = n-1; i >0; --i) C[i] = C[i-1];
  C[0] = 0;
}
rbwt_t *rbwt_bwt_init(int n_seqs, int l_seq)
{
  int rot;
  rbwt_t *rbwt = calloc(1, sizeof(rbwt_t));
  rbwt->n_seqs = n_seqs;
  rbwt->n_rot = l_seq;
  rbwt->cnt_table = calloc(256, sizeof(uint32_t));
  rbwt->rbwt = calloc(l_seq, sizeof(__rbwt_t)); 
  for(rot = 0; rot < l_seq; ++rot){
    rbwt->rbwt[rot].bwt = calloc((n_seqs*BITS_PER_NT+NT_PER_BYTE-1)/NT_PER_BYTE, sizeof(uint8_t)); 
    rbwt->rbwt[rot].Occ = calloc((n_seqs+OCC_INTV-1)/OCC_INTV, 4*sizeof(uint32_t)); 
  }
  return rbwt;
}
rbwt_t *rbwt_bwt_destroy(rbwt_t *rbwt)
{
  int rot;
  for(rot = 0; rot < rbwt->n_rot; ++rot){
    free(rbwt->rbwt[rot].bwt); 
    free(rbwt->rbwt[rot].Occ); 
  }
  free(rbwt->rbwt);
  free(rbwt->cnt_table);
  free(rbwt);

}



rbwt_t *rbwt_bwt_build(int n_seqs, uint32_t *sorted_seqs, uint32_t *sa0, uint32_t *sa1)
{
  uint32_t i, j;
  uint32_t *last_sa, *cur_sa;
  rbwt_t *rbwt = rbwt_bwt_init(n_seqs, LEN_SEQ);
  rbwt_gen_cnt_table(rbwt->cnt_table);
  //rbwt_t *rbwt = calloc(1, sizeof(rbwt_t));
  //rbwt->bwt = calloc((n_seqs+3)/4*16, sizeof(uint8_t));
  //uint8_t *bwt = rbwt->bwt;
  for(j=0; j< NT_SIZE; ++j) rbwt->rbwt[0].C[j] = 0;
  for(j = 0; j < n_seqs; ++j){
    uint32_t c = sorted_seqs[j]&3;
    __set_bwt(rbwt->rbwt[0].bwt, j, c);
    rbwt->rbwt[0].C[c]++;
  }
  rbwt_accumulate_C(NT_SIZE, rbwt->rbwt[0].C);
  for(j =0; j < n_seqs; ++j) sa0[j] = j;
  last_sa = sa0; cur_sa = sa1;
  for(i = 1; i < rbwt->n_rot; ++i){
    uint32_t C[NT_SIZE];
    for(j = 0; j < NT_SIZE; ++j) C[j] = rbwt->rbwt[i-1].C[j];
    for(j = 0; j < n_seqs; ++j){
      uint32_t c = __get_nt(sorted_seqs[last_sa[j]], i-1);
      cur_sa[C[c]++] = last_sa[j];
    }
    //generate bwt from cur_sa
    for(j=0; j< NT_SIZE; ++j) rbwt->rbwt[i].C[j] = 0;
    for(j = 0; j < n_seqs; ++j){
      uint32_t c = __get_nt(sorted_seqs[cur_sa[j]], i);
      //int n_alned = (n_seqs+NT_PER_BYTE-1)/NT_PER_BYTE*NT_PER_BYTE;
      __set_bwt(rbwt->rbwt[i].bwt, j, c);
      rbwt->rbwt[i].C[c]++;
    }
    rbwt_accumulate_C(NT_SIZE, rbwt->rbwt[i].C);
    SWAP(uint32_t *, last_sa, cur_sa); 
  }
  //fprintf(stderr, "\n==Iter = %u\n", i-1);
  //log_seq16array(n, last_sa, lext, i-1);
  return rbwt;
}


void rbwt_gen_C_8(uint8_t *bwt, int n_seqs, uint32_t C[][NT_SIZE])
{ 
    uint8_t c;
    uint32_t bg, ed, i, rot;
    for(rot=0; rot<LEN_SEQ; ++rot){
        for(c=0; c<NT_SIZE; ++c) {C[rot][c] = 0;} 
    }
    int n_alned = (n_seqs+NT_PER_BYTE-1)/NT_PER_BYTE*NT_PER_BYTE;
    for(rot =0; rot<LEN_SEQ; ++rot){
      bg = rot*n_alned;
      ed = bg + n_seqs; 
      for(i = bg; i < ed; ++i){
        c = __get_bwt(bwt, i);  
        ++C[rot][c];
      }
      //log_array2(17, rot, C_2nt);
      rbwt_accumulate_C(NT_SIZE, C[rot]);
      //log_array2(17, rot, C_2nt);
    }
}

void rbwt_gen_C_32(uint8_t *bwt, int n_seqs, uint32_t C[][NT_SIZE], uint32_t (*Occ[])[NT_SIZE])
{
  int rot; 
  uint32_t n_nt[NT_SIZE];
  uint32_t bg, ed, i;
  uint8_t c;
  for(rot=0; rot<LEN_SEQ; ++rot){
    for(c=0; c<NT_SIZE; ++c) {C[rot][c] = 0;} 
  }
  int n_alned = (n_seqs+NT_PER_BYTE-1)/NT_PER_BYTE*NT_PER_BYTE;
  
  for(rot =0; rot<LEN_SEQ; ++rot){
    uint32_t (*occ)[NT_SIZE] = Occ[rot];
    bg = rot*n_alned;
    ed = bg + n_seqs; 
    for(c = 0; c < NT_SIZE; ++c) n_nt[c] = 0;
    for(i = bg; i < ed; ++i){
      uint32_t ii = i - bg;
      if(ii%OCC_INTV == 0) {
        for(c =0 ; c < NT_SIZE; ++c) occ[ii/OCC_INTV][c] = n_nt[c]; 
      }
      c = __get_bwt(bwt, i);  
      ++n_nt[c];
    }
    for(c = 0; c < NT_SIZE; ++c) C[rot][c] = n_nt[c];
    //log_array2(17, rot, C_2nt);
    rbwt_accumulate_C(NT_SIZE, C[rot]);
    //log_array2(17, rot, C_2nt);
  }
}



/*  
fprintf(stderr, "\n----------------------\n");
if( DataNum <= 254*256 ) siz_C = 2 ;
if( DataNum > 254*256 ) siz_C = 4 ;
int rng_num = ((DataNum+253) /254) ;
int len_sum = (8+ siz_C) *( rng_num /8) + rng_num %8 + siz_C ;
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        fprintf(stderr, "rbwt_gen_C_all = %u   \t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            fprintf(stderr, "%u\t", C[__j+__l]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "-------------------\n\n");
}
fprintf(stderr, "\n+++++++++++++++++\n");
}*/
//+++++++++++++++++++++++++++++++
void rbwt_gen_C_all(const uint8_t* bwt, int n_data, uint8_t *C)
{
    //fprintf(stderr, "[rbwt_gen_C_all]\n");
    /*  
    if(n_data<=NO_BWT_SUM_SIZE) {
        fprintf(stderr, "Error in rbwt_gen_C_all!!\n"); 
        exit(1);
    } else if(n_data <= 255) rbwt_gen_C_8(bwt, n_data, (uint8_t (*)[17])C);
    else rbwt_gen_C2(bwt, n_data, C);
    */
    if(n_data <= 255) rbwt_gen_C_8(bwt, n_data, (uint8_t (*)[17])C);
    return;
}

//++++++++++++++++++++++++++++++++++++++
/*
if(n_data == 419){
int n_data = n_data;
int siz_C; 
uint32_t __i, __j, __k, __l;
fprintf(stderr, "\n----------------------\n");
if( n_data <= 254*256 ) siz_C = 2 ;
if( n_data > 254*256 ) siz_C = 4 ;
int rng_num = ((n_data+253) /254) ;
int len_sum = (8+ siz_C) *( rng_num /8) + rng_num %8 + siz_C ;
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        fprintf(stderr, "rbwt_gen_C2 =  %4u   \t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            fprintf(stderr, "%u\t", C[__j+__l]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "-------------------\n\n");
}
fprintf(stderr, "\n+++++++++++++++++\n");
}
*/
//++++++++++++++++++++++++++++++++++++++
//生成序列
uint32_t __dna_count(uint8_t *bwt, uint32_t k , uint32_t l, uint8_t c)
{
  uint32_t i, n;
  for(n=0, i =k; i <l; ++i){
    n += (__get_bwt(bwt, i)==c)?1:0;
  }
  return n;

}
#define NT_PER_BYTE 4
uint32_t __dna_count2(uint8_t *bwt, uint32_t k , uint32_t l, uint8_t c, uint32_t cnt_table[])
{
  uint32_t i, n;
  uint32_t k0 = (k+NT_PER_BYTE-1)/NT_PER_BYTE*NT_PER_BYTE; 
  uint32_t l0 = l/NT_PER_BYTE*NT_PER_BYTE; 
  n = 0;
  for(i =k; i <k0; ++i){
    n += (__get_bwt(bwt, i)==c)?1:0;
  }
  for(i =k0; i <l0; i += 4){
    n += (cnt_table[bwt[i>>2]]>>(c<<3))&255;
  }
  for(i =l0; i <l; ++i){
    n += (__get_bwt(bwt, i)==c)?1:0;
  }
  int n0 = __dna_count(bwt, k, l, c);
  if(n != n0) {   
    printf("n = %u, n0 = %u, k = %d, k0 = %d, l0 = %d, l =%d\n", n, n0, k, k0, l0, l); 
    exit(1);
  }
  return n;
}



/*  
void rbwt_test_smbwt(struct SubBuf *sub, int flg)
{
    uint8_t *seq_buf, *seqs0, *bwt;
    uint32_t num, *sort_seq;
    if(flg == 0){
        bwt = sub->bwt_seqL;
        num = sub->num_seqL;
        sort_seq = sub->sortL;
    } else if(flg == 1){
        bwt = sub->bwt_seqR;
        num = sub->num_seqR;
        sort_seq = sub->sortR;
    }
    seq_buf = malloc(16*num*2);  
    rbwt_test_rbwt(bwt, num, seq_buf);
    int i, j;
    seqs0 = seq_buf+num*16;
    for(i = 0; i < num; ++i){
        uint8_t ch;
        for(j = 0; j < 16; ++j){
            seqs0[i*16+j] = (sort_seq[i]>>((15-j)*2))&3;      
        }
    }

    for(i = 0; i < num; ++i){
        for(j = 0; j < 16; ++j){
            //printf("[%u, %u]: %u %u\n", i, j, seq_buf[i*16+j], seqs0[i*16+j]);
            if(seq_buf[i*16+j] != seqs0[i*16+j]){
                fprintf(stderr, "test_smbwt error!!!!\n"); 
                exit(1);
            
            }  
        }
       
    
    
    }
    free(seq_buf);
    return;
hi IndentGuidesEven ctermbg=lightgrey
}
*/
void rbwt_test_build()
{
  int n_seqs = 500;
  uint32_t *pac_seqs = calloc(n_seqs, sizeof(uint32_t));
  
  int i, j;
  srand(time(NULL));
  for(i = 0; i < n_seqs; ++i) {
    pac_seqs[i] = rand()%(uint32_t)(-1); 
  }
  ks_introsort_uint32_t(n_seqs, pac_seqs);
  uint32_t *sa0 = calloc(n_seqs, sizeof(uint32_t));
  uint32_t *sa1 = calloc(n_seqs, sizeof(uint32_t));
  rbwt_t *rbwt = rbwt_bwt_build(n_seqs, pac_seqs, sa0, sa1);
  uint8_t *seqs1 = calloc(n_seqs*LEN_SEQ, sizeof(uint8_t));
  rbwt_test_rbwt(rbwt, seqs1);
  rbwt_bwt_destroy(rbwt); 
  
  uint8_t *seqs0 = calloc(n_seqs*LEN_SEQ, sizeof(uint8_t));
  for(i = 0; i < n_seqs; ++i){
    for(j = 0; j < LEN_SEQ; ++j){
        seqs0[i*LEN_SEQ+j] = (pac_seqs[i]>>((LEN_SEQ-1-j)*2))&3;      
    }
  }
  for(i = 0; i < n_seqs; ++i){
    int eq = 1;
    for(j = 0; j < LEN_SEQ; ++j){
      if(seqs1[i*LEN_SEQ+j] != seqs0[i*16+j]){
        eq = 0;
      }  
    }
    printf("%c\n", "!="[eq]);
    for(j = 0; j < LEN_SEQ; ++j) printf("%d", seqs0[i*LEN_SEQ+j]);
    putchar('\n');
    for(j = 0; j < LEN_SEQ; ++j) printf("%d", seqs1[i*LEN_SEQ+j]);
    putchar('\n');
    putchar('\n');
  }
  free(pac_seqs);  
  free(sa0);
  free(sa1);
  free(seqs0);
  free(seqs1);
}
uint32_t rbwt_occ(rbwt_t *rbwt, int rot, uint32_t i, uint8_t c)
{
  /*
  int n_seqs = rbwt->n_seqs;
  uint32_t n_alned = (n_seqs+NT_PER_BYTE-1)/NT_PER_BYTE*NT_PER_BYTE;
  uint32_t bg = n_alned*rot;
  uint32_t occ = rbwt->Occ[rot][i/256][c] + __dna_count2(bwt, bg + i/OCC_INTV*OCC_INTV, bg+i, c, rbwt->cnt_table);
  */
  uint32_t occ = rbwt->rbwt[rot].Occ[i/256][c] + __dna_count2(rbwt->rbwt[rot].bwt, i/OCC_INTV*OCC_INTV, i, c, rbwt->cnt_table);
  return occ;
}
void rbwt_test_rbwt(rbwt_t *rbwt, uint8_t *seq_buf)
{
  //恢复序列算法
  int i, rot;
  rot = 0; 
  uint8_t seq[LEN_SEQ];
  for(i = 0; i < rbwt->n_seqs; ++i){
    uint32_t k = i;
    for(rot=0; rot< rbwt->n_rot; ++rot){
      //uint32_t n_alned = (n_seqs+NT_PER_BYTE-1)/NT_PER_BYTE*NT_PER_BYTE;
      //uint32_t rot_bg = n_alned*rot;
      //uint32_t ii = rot_bg + k;
      uint8_t c = __get_bwt(rbwt->rbwt[rot].bwt, k); 
      seq[LEN_SEQ-1-rot] = c; 
      uint32_t occ = __dna_count2(rbwt->rbwt[rot].bwt, 0, k, c, rbwt->cnt_table);
      k = rbwt->rbwt[rot].C[c] + occ;
    }
    memcpy(seq_buf+i*16, seq, 16); 
  }    

  return;
}
#ifdef RBWT_MAIN
int main()
{
  rbwt_test_build();
}
#endif
