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
#define MAX_NAME 1024
#define USE_MALLOC_WRAPPERS
		/* ----------  end of function main  ---------- */
bwtint_t locate_step(const bwt_t *bwt, bwtint_t k)
{
	bwtint_t sa = 0;
	while (k % 16 != 0) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	//printf("pos :%d\n", sa + bwt->sa[k/bwt->sa_intv]);
    return sa;
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
    return idx;
}
void idx_destroy(idx_t *idx)
{

    free(idx->bwt->bwt);
    free(idx->bwt);
    //bwt_destroy(idx->bwt);
    bns_destroy(idx->bns);
    if(idx->pac != NULL ) free(idx->pac);
    free(idx);
}

#ifdef MAIN_LOCATE
int main(int argc, char *argv[]){
    const *prefix = argv[1];
    idx_t *idx = idx_restore(prefix);
    uint32_t i;
    uint64_t tot = 0;
    for(i =0; i < idx->bwt->seq_len; ++i) {
        int l = locate_step(idx->bwt, i);
        tot += l;
        printf("%u\t%u\n", i, l); 
    }
    printf("tot = %ul, average = %u\n", tot, tot/idx->bwt->seq_len);
    return;
}
#endif
