/*
 * =====================================================================================
 *
 *       Filename:  bitmap.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017年10月03日 17时08分24秒
 *       Revision:  none 
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (wq), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */


#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "bitmap.h"

#define INTV_BP 256
#define BITS_PER_WORD 32
#define BP_WORD(i) ((i)/BITS_PER_WORD) 
#define BP_MASK(i) (1U<<(i%BITS_PER_WORD))
static uint32_t __BITS_COUNT(uint32_t v)
{
    v = v-((v>>1)&0x55555555);
    v = (v&0x33333333)+((v>>2)&0x33333333);
    return ((v+(v>>4)&0xF0F0F0F)*0x1010101)>>24;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  bp_init
 *  Description: Init a bit set  capable of holding l integers, [0, l) 
 * =====================================================================================
 */

bp_t *bp_init (size_t l)
{
    bp_t *bp = (bp_t *)calloc(1, sizeof(bp_t));
    bp->l = l;
    bp->n_a = (l-1)/BITS_PER_WORD +1;
    bp->m = bp->n_a *BITS_PER_WORD; 
    bp->intv = INTV_BP;
    bp->n_occ = bp->l / bp->intv+2;
    bp->a = (uint32_t *)calloc(bp->n_a, sizeof(uint32_t));
    bp->occ = (size_t *)calloc(bp->n_occ, sizeof(size_t)); 
    return bp;
}		/* -----  end of function bp_init  ----- */



/*  
 * ===  FUNCTION  ======================================================================
 *         Name:  bp_destroy
 *  Description:  
 * ====================================================================================
 */
void bp_destroy ( bp_t *bp )
{
    free(bp->a);
    free(bp->occ);
    free(bp);
}		/* -----  end of function bp_destroy  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  bp_set1
 *  Description:  
 * =====================================================================================
 */

int bp_set1(bp_t *bp, size_t i)
{
    if(i > bp->l) {
        fprintf(stderr, "[bitmap]: index %ld out of range!\n", i);
        exit(1);
    }
    //uint32_t *word = bp->a+i/BITS_PER_WORD;
    bp->a[BP_WORD(i)] |= BP_MASK(i);//from left to right      
    return 0;
}		/* -----  end of function bp_set1  ----- */
int bp_set0(bp_t *bp, size_t i)
{
    if(i > bp->l) {
        fprintf(stderr, "[bitmap]: index %ld out of range!\n", i);
        exit(1);
    }
    //uint32_t *word = bp->a+i/BITS_PER_WORD;
    //*word &= ~(1<<(i%BITS_PER_WORD));//from left to right      
    bp->a[BP_WORD(i)] &= ~BP_MASK(i);
    return 0;
}		/* -----  end of function bp_set1  ----- */


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  bp_get
 *  Description:  
 * =====================================================================================
 */

bp_t *bp_restore(const char *fn, int restore_occ)
{
    FILE *fp = fopen(fn, "r");
    bp_t *bp = (bp_t *)calloc(1, sizeof(bp_t));
    fread(&bp->l, sizeof(size_t), 1,  fp); 
    fread(&bp->m, sizeof(size_t), 1, fp); 
    fread(&bp->intv, sizeof(int), 1, fp); 
   
    fread(&bp->n_a, sizeof(size_t), 1, fp); 
    bp->a = (uint32_t *)calloc(bp->n_a, sizeof(uint32_t));
    fread(bp->a,sizeof(uint32_t), bp->n_a, fp ); 
    if(restore_occ == 0) return bp;  
    fread(&bp->n_occ, sizeof(int), 1, fp); 
    bp->occ = (size_t *)calloc(bp->n_occ, sizeof(size_t));
    fread(bp->occ,sizeof(size_t), bp->n_occ, fp ); 
    fclose(fp);
    return bp;
}
void bp_dump(bp_t *bp, const char *fn)
{
    FILE *fp = fopen(fn, "w");
    fwrite(&bp->l, sizeof(size_t), 1, fp); 
    fwrite(&bp->m, sizeof(size_t), 1, fp); 
    fwrite(&bp->intv, sizeof(int), 1, fp); 
    fwrite(&bp->n_a, sizeof(size_t), 1, fp); 
    fwrite(bp->a, sizeof(uint32_t), bp->n_a, fp ); 
    fwrite(&bp->n_occ, sizeof(int), 1, fp); 
    fwrite(bp->occ,sizeof(size_t), bp->n_occ, fp ); 
    fclose(fp);
}

int bp_get(bp_t *bp, size_t i)
{
    if(i >= bp->l) {
        fprintf(stderr, "[bitmap]: index %ld out of range!\n", i);
        exit(1);
    }
    //uint32_t *word = bp->a+i/BITS_PER_WORD;
    //return 1&(*word >>(i%BITS_PER_WORD));      
    return (bp->a[BP_WORD(i)]&BP_MASK(i))!= 0;
}		/* -----  end of function bp_get  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  bp_reset
 *  Description:  
 * =====================================================================================
 */

void bp_reset(bp_t *bp)
{

    memset(bp, 0, bp->n_a);  
}		/* -----  end of function bp_reset  ----- */

/* 
 * ===  FUNCTION  =============================================,=========================
 *         Name:  bp_gen_occ
 *  Description:  
 * =====================================================================================
 */

void bp_gen_occ ( bp_t* bp)
{
    size_t i,c=0;
    uint32_t *word = bp->a;
    size_t *occ = bp->occ;
    for(i=0; i < bp->n_a; ++i){
        if((i*BITS_PER_WORD%INTV_BP)== 0){
            occ[i*BITS_PER_WORD/INTV_BP] = c;
            c = 0;
        }
        c += __BITS_COUNT(word[i]); 
    
    }
    for(i=1; i < bp->n_occ; ++i){
        occ[i] += occ[i-1]; 
        //printf("%u\n", occ[i]); 
    }
}		/* -----  end of function bp_gen_occ  ----- */

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  bp_gen_occ
 *  Description: count number of 1 before i(i included) 
 * =====================================================================================
 */

uint32_t bp_rank(bp_t *bp, size_t i)
{
    uint32_t j, c= 0;
    const uint32_t *word = bp->a;
    for(j = i/INTV_BP*INTV_BP/BITS_PER_WORD;j < i/BITS_PER_WORD;++j){
         c += __BITS_COUNT(word[j]);
    }

    c+= __BITS_COUNT(word[j]<<(BITS_PER_WORD-(i%BITS_PER_WORD)-1));//shift arithmetic left?

    return bp->occ[i/INTV_BP] + c;
}		/* -----  end of function bp_gen_occ  ----- */
/*
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================*/

#ifdef MAIN_BITMAP
#define MAX_NUM 100
static inline int __cmp(const void *x, const void* y)
{
    uint32_t a = *(uint32_t *)x;
    uint32_t b = *(uint32_t *)y;
    if(a<b) return -1;
    else if(a == b) return 0;
    else return 1;
}
int main ( int argc, char *argv[] )
{
    srand(time(NULL));
    uint32_t i;
    uint32_t a[MAX_NUM] = {};
    bp_t *bp = bp_init(1000);
    for(i=0; i< MAX_NUM; ++i){
        a[i] = rand()%1000;
        bp_set1(bp, a[i]); 
    }
    for(i=0;i < bp->l; ++i){
        if(bp_get(bp, i)==1) printf("%u\t", i); 
    }
    printf("\n");
    bp_gen_occ(bp);
    for(i=0;i < bp->l; ++i){
        printf("%u:%u\n", i, bp_rank(bp, i)); 
    }
    printf("\n");


    qsort(a, MAX_NUM, sizeof(uint32_t), __cmp);
    for(i=0;i < MAX_NUM; ++i){
        printf("%u\t", a[i]); 
    }
    printf("\n"); 
    bp_destroy(bp);
    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
#endif
