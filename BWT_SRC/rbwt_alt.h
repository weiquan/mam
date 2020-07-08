#ifndef __ROT_BWT_H
#define __ROT_BWT_H

/*
 * =====================================================================================
 *
 *       Filename:  rbwt.h
 *
 *    Description:  rotation bwt
 *
 *        Version:  1.0
 *        Created:  05/19/2020 05:14:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#define __set_bwt(bwt, l, c) ((bwt)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define __get_bwt(bwt, l) ((bwt)[(l)>>2]>>((~(l)&3)<<1)&3)
#define __set_bwt2(bwt, l, c) ((bwt)[(l)>>1] |= (c)<<((~(l)&1)<<2))
#define __get_bwt2(bwt, l) ((bwt)[(l)>>1] >> ((~(l)&1)<<2)&15)
#define __get_nt(seq, i) (((seq)>>((i)*2))&0x3)
#define MIN_BWT_SIZE 0
#define NO_BWT_SUM_SIZE 64
#define NT_SIZE 4
#define BITS_PER_NT 2
#define LEN_SEQ 16
#define NT_PER_BYTE 4
#define OCC_INTV 256
typedef struct{
  uint8_t *bwt;
  uint32_t (*Occ)[NT_SIZE];
  uint32_t C[NT_SIZE];
} __rbwt_t;
typedef struct{
  int n_seqs;
  int n_rot;
  __rbwt_t *rbwt;
} rbwt_t;

rbwt_t *rbwt_bwt_build(int n_seqs, int n_rot, uint32_t *sorted_seqs, uint32_t *sa0, uint32_t *sa1);
rbwt_t *rbwt_bwt_restore(FILE *fp);
int rbwt_bwt_dump(rbwt_t *rbwt, FILE *fp);
rbwt_t *rbwt_bwt_destroy(rbwt_t *rbwt);
uint32_t rbwt_exact_match(rbwt_t *rbwt, int n_rot, uint32_t *cnt_table, int l_seq, uint8_t *seq, uint32_t *k, uint32_t *l);
#endif
