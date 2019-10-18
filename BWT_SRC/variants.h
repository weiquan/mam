/*
 * =====================================================================================
 *
 *       Filename:  hapmap.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/14/2012 05:23:38 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#ifndef __VARIANTS_H_
#define __VARIANTS_H_
#include <stdint.h>
#include <stdio.h>
#include "khash.h"
#include "kvec.h"
KHASH_MAP_INIT_STR(id, int)
typedef struct{
    uint32_t start;//pack reference pos
    uint32_t end;//variant len,eg: ins = 0, sub =1
    char varclass;
    char* orign; 
    char* var;

} var_t;
typedef struct{
    char *name;
    uint32_t offset; 
    uint32_t len;
} chromvar_t;

typedef struct{
    uint32_t n_chrom;    
    chromvar_t *chromvar;
    uint32_t n_var;   
    khash_t(id) *h;//rname to rid 
    var_t *var;
    FILE *fp;
} all_var_t;

all_var_t *variants_restore(char *fn);
void variants_dump(all_var_t *av, const char *fn);
all_var_t *variants_build(const char *fp);
void variants_destroy(all_var_t *av);
#endif
