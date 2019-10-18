/*
 * =====================================================================================
 *
 *       Filename:  bsearch.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2017年10月19日 10时33分54秒
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
#include "bsearch.h"

bs_iter_t lower_bound(bs_iter_t first, bs_iter_t last, bs_val_t *a, bs_val_t val)
{
    bs_iter_t half, len, middle;
    len = last-first; 
    while(len >0){
        half = len>>1; 
        middle = first;
        middle += half;
        if(a[middle] < val){
            first = middle;
            ++first;
            len = len-half-1;
        
        } else{ len = half;}
    
    }  
    return first;
}
bs_iter_t uper_bound(bs_iter_t first, bs_iter_t last, bs_val_t *a, bs_val_t val)
{
    bs_iter_t half, len, middle;
    len = last-first;
    while(len >0){
        half = len>>1;
        middle = first;
        middle += half;
        if(val < a[middle]) len = half; 
        else{
            first = middle;
            ++first;
            len = len-half-1;
        
        }
    }
    return first;
}
int test_bsearch()
{
    int i, k;
    uint32_t A[] = {1,2,3,3,3,5,8};
    
    for(i = 0; i < 7; ++i){
        printf("A[%d]: %u\n", i, A[i]); 
    }
    printf("\nlower bound search\n");
    for(i = 0; i < 10; ++i){
        k = lower_bound(0, 7, A, i);
        printf("Searching for %d: Result: index = %d, A[%u]== %u\n", i, k, k, A[k]); 
    }

    printf("\nuper bound search\n");
    for(i = 0; i < 10; ++i){
        k = uper_bound(0, 7, A, i);
        printf("Searching for %d: Result: index = %d, A[%u]== %u\n", i, k, k, A[k]); 
    }


    return 1;
}

#ifdef BSEARCH_MAIN

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
int main ( int argc, char *argv[] )
{
    test_bsearch();
    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
#endif
