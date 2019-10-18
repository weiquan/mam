#ifndef __LOOLUPTABLE_H
#define __LOOLUPTABLE_H
/*
 * =====================================================================================
 *
 *       Filename:  LookUpTable.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  01/07/2013 02:59:17 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */
#include <stdint.h>
typedef struct{
    uint32_t maxLookupLen;
    uint32_t n_item;
    uint32_t *item;
} lkt_t;

lkt_t *lkt_init(const int maxLookupLen);
void lkt_destroy(lkt_t *lkt);
lkt_t *lkt_restore(const char *fn_lkt);
uint32_t lkt_seq2LktItem(const uint8_t *seq, int from, int to);
void lkt_build_lookuptable(const char *fn_pac, lkt_t *lkt);
void lkt_dump(const char *fn_lkt, lkt_t *lkt);
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  lkt_lookup_sa
 *  Description:  return SA interval [l, k] of seq[from:to] 
 * =====================================================================================
 */
static inline int lkt_lookup_sa(lkt_t *lkt, const uint8_t *seq, int from, int to, uint32_t *l, uint32_t *k)
{
    uint32_t i;
    
    i = lkt_seq2LktItem(seq, from, to);
    if(i == 0xFFFFFFFF) {
        *l = 1;
        *k = 0;
        return 0;
    }

    *l = lkt->item[i];
    *k = lkt->item[i+1]-1;
    if(*l > *k) return 0;
    else return *k -*l +1;
}

#endif
