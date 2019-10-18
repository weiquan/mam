/*
 * =====================================================================================
 *
 *       Filename:  test.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年03月01日 12时17分20秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (wq), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>

#include	<stdlib.h>

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>

#define ROOT_PATH "../Index"
int main ( int argc, char *argv[] )
{
    char fn[1024];
    sprintf(fn, "%s/data/seedidx/seedidx_%u.txt", ROOT_PATH, 1);
    //FILE *fp = fopen64("../Index/data/seedidx/seed_idx0", "wb");
    FILE *fp = fopen(fn, "wb");
    if(fp != NULL) { 
        fprintf(stderr, "file %s ok\n", fn);
        fclose(fp);}
    else{fprintf(stderr, "path wrong\n");}
    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
