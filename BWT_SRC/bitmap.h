/*
 * =====================================================================================
 *
 *       Filename:  bitmap.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017年10月03日 17时09分08秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (wq), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef __BITMAP_H
#define __BITMAP_H
#include <stdint.h>
#define size_t int64_t
typedef struct {
    size_t l, m, n_a;
    uint32_t *a;//n_a = m*(sizeof(uint32_t))
    int intv;
    int n_occ;
    size_t *occ;
} bp_t;



void bp_dump(bp_t *bp, const char *fn);
bp_t *bp_restore(const char *fn, int restore_occ);
bp_t *bp_init (size_t l);
void bp_destroy ( bp_t *bp );

int bp_get(bp_t *bp, size_t i);
int bp_set1(bp_t *bp, size_t i);
int bp_set0(bp_t *bp, size_t i);
void bp_reset(bp_t *bp);
uint32_t bp_rank(bp_t *bp, size_t i);
void bp_gen_occ ( bp_t* bp);
#endif
