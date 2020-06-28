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
#include "rbwt_alt.h"

//KSORT_INIT_GENERIC(uint32_t)
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
int rbwt_bwt_dump(rbwt_t *rbwt, FILE *fp)
{
  int i;
  int n_seqs = rbwt->n_seqs;
  int n_rot = rbwt->n_rot;
  int bwt_size = (n_seqs*BITS_PER_NT+NT_PER_BYTE-1)/NT_PER_BYTE;
  //int occ_size =  (n_seqs+OCC_INTV-1)/OCC_INTV;
  int occ_size =  rbwt->n_seqs/OCC_INTV + 1;
  fwrite(&n_seqs, sizeof(int), 1, fp); 
  fwrite(&n_rot, sizeof(int), 1, fp); 
  for(i = 0; i < n_rot; ++i) {
    fwrite(rbwt->rbwt[i].bwt, sizeof(uint8_t), bwt_size, fp); 
    fwrite(rbwt->rbwt[i].Occ, 4*sizeof(uint32_t), occ_size, fp); 
    fwrite(rbwt->rbwt[i].C, sizeof(uint32_t), NT_SIZE, fp); 
  }
  return 0;
}
rbwt_t *rbwt_bwt_restore(FILE *fp)
{
  int i;
  rbwt_t *rbwt = calloc(1, sizeof(rbwt_t));
  fread(&rbwt->n_seqs, sizeof(int), 1, fp); 
  fread(&rbwt->n_rot, sizeof(int), 1, fp); 
  rbwt->rbwt = calloc(rbwt->n_rot, sizeof(__rbwt_t)); 
  int bwt_size = (rbwt->n_seqs*BITS_PER_NT+NT_PER_BYTE-1)/NT_PER_BYTE;
  //int occ_size =  (rbwt->n_seqs+OCC_INTV-1)/OCC_INTV;
  int occ_size =  rbwt->n_seqs/OCC_INTV + 1;
  for(i = 0; i < rbwt->n_rot; ++i) {
    rbwt->rbwt[i].bwt = calloc(bwt_size, sizeof(uint8_t)); 
    rbwt->rbwt[i].Occ = calloc(occ_size, 4*sizeof(uint32_t)); 
    fread(rbwt->rbwt[i].bwt, sizeof(uint8_t), bwt_size, fp); 
    fread(rbwt->rbwt[i].Occ, 4*sizeof(uint32_t), occ_size, fp); 
    fread(rbwt->rbwt[i].C, sizeof(uint32_t), NT_SIZE, fp); 
  }
  return rbwt;
}


rbwt_t *rbwt_bwt_init(int n_seqs, int l_seq)
{
  int rot;
  rbwt_t *rbwt = calloc(1, sizeof(rbwt_t));
  rbwt->rbwt = calloc(l_seq, sizeof(__rbwt_t)); 
  rbwt->n_seqs = n_seqs;
  rbwt->n_rot = l_seq;
  for(rot = 0; rot < l_seq; ++rot){
    rbwt->rbwt[rot].bwt = calloc((n_seqs*BITS_PER_NT+NT_PER_BYTE-1)/NT_PER_BYTE, sizeof(uint8_t)); 
    rbwt->rbwt[rot].Occ = calloc(n_seqs/OCC_INTV+1, 4*sizeof(uint32_t)); 
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
  free(rbwt);

}
#define LEN_EXT 16

rbwt_t *rbwt_bwt_build(int n_seqs, int n_rot, uint32_t *sorted_seqs, uint32_t *sa0, uint32_t *sa1)
{
  uint32_t i, j, k;
  uint32_t *last_sa, *cur_sa;
  rbwt_t *rbwt = rbwt_bwt_init(n_seqs, LEN_SEQ);
  for(j=0; j< NT_SIZE; ++j) rbwt->rbwt[0].C[j] = 0;
  for(j = 0; j < n_seqs; ++j){
    uint32_t c = sorted_seqs[j]&3;
    __set_bwt(rbwt->rbwt[0].bwt, j, c);
    if(j % OCC_INTV == 0) {
      for(k = 0; k < NT_SIZE; ++k) rbwt->rbwt[0].Occ[j/OCC_INTV][k] = rbwt->rbwt[0].C[k];
    }
    rbwt->rbwt[0].C[c]++;
  }
  if(j % OCC_INTV == 0){ 
    for(k = 0; k < NT_SIZE; ++k) 
      rbwt->rbwt[0].Occ[j/OCC_INTV][k] = rbwt->rbwt[0].C[k];
  }
  rbwt_accumulate_C(NT_SIZE, rbwt->rbwt[0].C);
  for(j =0; j < n_seqs; ++j) sa0[j] = j;
  last_sa = sa0; cur_sa = sa1;
  for(i = 1; i < n_rot; ++i){
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
      __set_bwt(rbwt->rbwt[i].bwt, j, c);
      if(j % OCC_INTV == 0){ 
        for(k = 0; k < NT_SIZE; ++k) 
          rbwt->rbwt[i].Occ[j/OCC_INTV][k] = rbwt->rbwt[i].C[k];
      }
      rbwt->rbwt[i].C[c]++;
    }
    if(j % OCC_INTV == 0){ 
      for(k = 0; k < NT_SIZE; ++k) 
        rbwt->rbwt[i].Occ[j/OCC_INTV][k] = rbwt->rbwt[i].C[k];
    }

    rbwt_accumulate_C(NT_SIZE, rbwt->rbwt[i].C);
    SWAP(uint32_t *, last_sa, cur_sa); 
  }
  return rbwt;
}

uint32_t __dna_count(uint8_t *bwt, uint32_t k , uint32_t l, uint8_t c)
{
  uint32_t i, n;
  for(n=0, i =k; i < l; ++i){
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
  for(i =l0; i < l; ++i){
    n += (__get_bwt(bwt, i)==c)?1:0;
  }
  int n0 = __dna_count(bwt, k, l, c);
  if(n != n0) {   
    printf("n = %u, n0 = %u, k = %d, k0 = %d, l0 = %d, l =%d\n", n, n0, k, k0, l0, l); 
    exit(1);
  }
  return n;
}



uint32_t rbwt_occ(rbwt_t *rbwt, int rot, uint32_t i, uint8_t c, uint32_t cnt_table[])
{
  if(i == 0) return 0;
  uint32_t occ = rbwt->rbwt[rot].Occ[i/OCC_INTV][c] + __dna_count2(rbwt->rbwt[rot].bwt, i/OCC_INTV*OCC_INTV, i, c, cnt_table);
  //fprintf(stderr, "occ0= %u, occ1 = %u, occ = %u, i = %u\n", rbwt->rbwt[rot].Occ[i/OCC_INTV][c], __dna_count2(rbwt->rbwt[rot].bwt, i/OCC_INTV*OCC_INTV, i, c, rbwt->cnt_table),__dna_count2(rbwt->rbwt[rot].bwt, 0, i, c, rbwt->cnt_table),i );
  return occ;
}
void rbwt_bwt_to_seqs(rbwt_t *rbwt, int n_rot, uint8_t *seq_buf, uint32_t *cnt_table)
{
  //恢复序列算法
  int i, rot;
  rot = 0; 
  uint8_t seq[LEN_SEQ];
  for(i = 0; i < rbwt->n_seqs; ++i){
    uint32_t k = i;
    for(rot=0; rot< n_rot; ++rot){
      uint8_t c = __get_bwt(rbwt->rbwt[rot].bwt, k); 
      seq[LEN_SEQ-1-rot] = c; 
      uint32_t occ = rbwt_occ(rbwt, rot, k, c, cnt_table);
#ifdef DEBUG  
      uint32_t occ1 = __dna_count2(rbwt->rbwt[rot].bwt, 0, k, c, cnt_table); 
      if(occ != occ1) {
        fprintf(stderr, "[occ != occ1], occ = %u, occ1 = %u, k = %u\n", occ, occ1, k);
        exit(1); 
      }
#endif
      k = rbwt->rbwt[rot].C[c] + occ;
    }
    memcpy(seq_buf+i*16, seq, 16); 
  }    

  return;
}
uint32_t rbwt_exact_match(rbwt_t *rbwt, int n_rot, uint32_t *cnt_table, int l_seq, uint8_t *seq, uint32_t *k, uint32_t *l)
{
  int i, rot;
  rot = 0; 
  uint32_t k0 = *k, l0 = *l;
  while(k0 < l0 && rot < n_rot) {
    uint8_t c = seq[l_seq-1-rot];
    if(c > 3) return 0;
    k0 = rbwt->rbwt[rot].C[c] + rbwt_occ(rbwt, rot, k0, c, cnt_table);
    l0 = rbwt->rbwt[rot].C[c] + rbwt_occ(rbwt, rot, l0, c, cnt_table);
    rot++; 
  }
  if(k0 < l0) {
    *k = k0;
    *l = l0;
    return l0 - k0; 
  }

  return 0;
}
uint32_t rbwt_appr_match(rbwt_t *rbwt, int n_rot, uint32_t *cnt_table, int l_seq, uint8_t *seq, uint32_t *k, uint32_t *l)
{
  int i, rot;
  rot = 0; 
  uint32_t k0 = *k, l0 = *l, k1, l1;
  while(k1 < l1 && rot < n_rot) {
    uint8_t c = seq[l_seq-1-rot];
    k0 = k1, l0 = l1;
    k1 = rbwt->rbwt[rot].C[c] + rbwt_occ(rbwt, rot, k0, c, cnt_table);
    l1 = rbwt->rbwt[rot].C[c] + rbwt_occ(rbwt, rot, l0, c, cnt_table);
    rot++; 
  }
  if(rot >= n_rot-1){ goto end; }
  if(rot < n_rot -1) {
    int ri;
    for(ri = 1; ri < 16; ++ri){
      int n_step = 0;
      rot = ri;
      k0 = *k, l0 = *l, k1, l1;
      while(k1 < l1 && n_step < n_rot) {
        uint8_t c = seq[l_seq-1-rot];
        k0 = k1, l0 = l1;
        k1 = rbwt->rbwt[rot].C[c] + rbwt_occ(rbwt, rot, k0, c, cnt_table);
        l1 = rbwt->rbwt[rot].C[c] + rbwt_occ(rbwt, rot, l0, c, cnt_table);
        rot = (rot+1)%16;
        n_step++; 
      }   
      if(n_step >= n_rot-1){ goto end; }
    } 
  
  }
end:
  if(k1 < l1) {
    *k = k1;
    *l = l1;
    return l0 - k0; 
  } else {
    *k = k0;
    *l = l0;
    return l0 - k0; 
  }
  return 0;
}


void rbwt_test_build()
{
  int n_seqs = 500;
  uint32_t *pac_seqs = calloc(n_seqs, sizeof(uint32_t));
  
  int i, j;
  srand(time(NULL));
  for(i = 0; i < n_seqs; ++i) {
    pac_seqs[i] = rand()%(uint32_t)(-1); 
  }
  ks_introsort(uint32_t, n_seqs, pac_seqs);
  uint32_t *sa0 = calloc(n_seqs, sizeof(uint32_t));
  uint32_t *sa1 = calloc(n_seqs, sizeof(uint32_t));
  
  uint32_t *cnt_table = calloc(256, sizeof(uint32_t)); 
  rbwt_gen_cnt_table(cnt_table);
  rbwt_t *rbwt = rbwt_bwt_build(n_seqs, LEN_SEQ, pac_seqs, sa0, sa1);
  for(i = 0; i < n_seqs; ++i) {
    uint8_t seq[16] = {};
    for(j = 0; j < LEN_SEQ; ++j){
        seq[j] = (pac_seqs[i]>>((LEN_SEQ-1-j)*2))&3;      
    }
    uint32_t k = 0, l = n_seqs;
    rbwt_exact_match(rbwt, 16, cnt_table, 16, seq, &k, &l);
    if(k != i || l != i+1) {
        fprintf(stderr, "[rbwt_exact:], i = %u, k= %u, l = %u\n", i, k, l);
        exit(1);  
    }
  }
  
  uint8_t *seqs1 = calloc(n_seqs*LEN_SEQ, sizeof(uint8_t));
  rbwt_bwt_to_seqs(rbwt, LEN_SEQ, seqs1, cnt_table);
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
  free(cnt_table);
  free(pac_seqs);  
  free(sa0);
  free(sa1);
  free(seqs0);
  free(seqs1);
}


#ifdef RBWT_MAIN
int main()
{
  rbwt_test_build();
}
#endif
