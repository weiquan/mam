/*
 * =====================================================================================
 *
 *       Filename:  hapmap.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/14/2012 08:04:55 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "khash.h"
#include "variants.h"

KHASH_MAP_INIT_STR(str, int)



#define MAX_CHAR 128

all_var_t *variants_build(const char*fn)
{

    FILE *fp = fopen(fn, "r");
    char line[MAX_CHAR];
    uint32_t n_var=0, n_chrom=0;
    

    if(fp == NULL){
        fprintf(stderr, "[variants_build] : fp == NULL!\n");
        exit(1);
    }
    int k, r;
    khash_t(str) *h = kh_init(str);//hash[rname] = count
    khash_t(id) *h1 = kh_init(id); //hash[rname] = rid
    while(fgets(line, MAX_CHAR, fp)!= NULL){
        char *rname = strdup((strtok(line, "\t")));
        k = kh_put(str, h, rname, &r); 

        if(r==0){ 
            kh_val(h, k) += 1;
            free(rname);
        }
        else if(r==1) {
            kh_val(h, k) = 1;
  
            k = kh_put(id, h1, rname, &r); 
            kh_val(h1, k) = n_chrom++;
        }
        ++n_var;
    }    
    all_var_t *av = (all_var_t*)calloc(1, sizeof(all_var_t));       
    av->n_chrom = n_chrom; 
    av->n_var = n_var;
    av->var = (var_t *)calloc(n_var, sizeof(var_t)); 
    av->chromvar = (chromvar_t *)calloc(n_chrom, sizeof(chromvar_t)); 

    av->fp = fp;
    av->h = h1;    
    //fill chromvar

    chromvar_t *chromvar = av->chromvar;
    const char *rname;uint32_t rid;
    /* test hash 
    kh_foreach(h1, rname,rid, (printf("%s\t%u\n", rname, rid))); 
    uint32_t count;
    kh_foreach(h, rname,count, (printf("%s\t%u\n", rname, count))); 
    */
    kh_foreach(h1, rname,rid, (chromvar[rid].name = (char *)rname)); 
    for(rid=0; rid< n_chrom; ++rid){
        chromvar_t *vars = chromvar+rid;  
        k = kh_get(str, h, vars->name);
        vars->len = kh_val(h, k); 
        if(rid == 0) vars->offset=0;
        else vars->offset = (vars-1)->offset+(vars-1)->len;
        //printf("%s\t%u\t%u\n", vars->name, vars->offset, vars->len);
    } 

    kh_clear(str, h); 
    //fill var 
    fseek(fp, 0, SEEK_SET);    
    uint32_t i=0;
    while(fgets(line, MAX_CHAR, fp)!= NULL){
       
        char *rname = strtok(line, "\t"); 
        uint32_t start = strtol(strtok(NULL, "\t"), NULL, 10); 
        uint32_t end =  strtol(strtok(NULL, "\t"), NULL, 10);  
        char varclass = strtok(NULL,"\t")[0]; 
        char *orign = strtok(NULL,"\t"); 
        char *var = strtok(NULL,"\t\n"); 

  
        k = kh_put(str, h, rname, &r);          
        if(r!=0) {
            k =  kh_get(id, h1, rname); 
            rid = kh_val(h1, k);
            i = av->chromvar[rid].offset;
        } 
        av->var[i].start = start;
        av->var[i].end = end;
        av->var[i].varclass = varclass;
        av->var[i].orign = strdup(orign);
        av->var[i].var = strdup(var);
        ++i;
    } 
    kh_destroy(str, h); 
    return av;
}
void variants_dump(all_var_t *av, const char *fn)
{
    uint32_t i;
    FILE *fp = fopen(fn, "w");
    fwrite(&(av->n_chrom), sizeof(uint32_t), 1, fp); 
    for(i = 0; i < av->n_chrom; ++i){
        chromvar_t *x = av->chromvar+i; 
        uint32_t l = strlen(x->name);
        fwrite(&l, sizeof(uint32_t), 1, fp); 
        fwrite(x->name, sizeof(char), l, fp); 
        fwrite(&(x->offset), sizeof(uint32_t), 1, fp);
        fwrite(&(x->len), sizeof(uint32_t), 1, fp);
    
    }
    fwrite(&(av->n_var), sizeof(uint32_t), 1, fp); 
    for(i =0; i < av->n_var; ++i){
        var_t *x = av->var+i;
        fwrite(&(x->start), sizeof(uint32_t), 1, fp);
        fwrite(&(x->end), sizeof(uint32_t), 1, fp);
        fwrite(&(x->varclass), sizeof(char), 1, fp);
        uint32_t l = strlen(x->orign);
        fwrite(&l, sizeof(uint32_t), 1, fp); 
        fwrite(x->orign, sizeof(char), l, fp); 
        l = strlen(x->var);
        fwrite(&l, sizeof(uint32_t), 1, fp); 
        fwrite(x->var, sizeof(char), l, fp); 
    } 
    
    fclose(fp);
}
all_var_t *variants_restore(char *fn)
{
    FILE *fp = fopen(fn, "r");
    all_var_t *av = (all_var_t *)calloc(1, sizeof(all_var_t)); 


    uint32_t i; 
    fread(&(av->n_chrom), sizeof(uint32_t), 1, fp); 
    av->chromvar = (chromvar_t *)calloc(av->n_chrom, sizeof(chromvar_t));
    for(i = 0; i < av->n_chrom; ++i){
        chromvar_t *x = av->chromvar+i; 
        uint32_t l;
        fread(&l, sizeof(uint32_t), 1, fp); 
        x->name = (char *)calloc(l+1 , sizeof(char));
        fread(x->name, sizeof(char), l, fp); 
        
        fread(&(x->offset), sizeof(uint32_t), 1, fp);
        fread(&(x->len), sizeof(uint32_t), 1, fp);
    
    }
    fread(&(av->n_var), sizeof(uint32_t), 1, fp); 
    av->var = (var_t *)calloc(av->n_var, sizeof(var_t));
    for(i =0; i < av->n_var; ++i){
        uint32_t l;
        var_t *x = av->var+i;
        fread(&(x->start), sizeof(uint32_t), 1, fp);
        fread(&(x->end), sizeof(uint32_t), 1, fp);
        fread(&(x->varclass), sizeof(char), 1, fp);

        fread(&l, sizeof(uint32_t), 1, fp); 
        x->orign = (char *)calloc(l, sizeof(char));
        fread(x->orign, sizeof(char), l, fp); 
        fread(&l, sizeof(uint32_t), 1, fp); 
        x->var = (char *)calloc(l, sizeof(char));
        fread(x->var, sizeof(char), l, fp); 
    } 
    khash_t(id) *h = kh_init(id); //hash[rname] = rid
    for(i = 0; i < av->n_chrom; ++i){
        chromvar_t *x = av->chromvar+i; 
        int r; uint32_t k;
        k = kh_put(id, h, x->name, &r);     
        kh_val(h, k) = i; 
    }
    av->h = h;
    av->fp = fp; 
    return av;
}
void variants_destroy(all_var_t *av)
{
    
    fclose(av->fp);   
    uint32_t i; 
    for(i =0; i < av->n_chrom; ++i){
        if(av->chromvar[i].name != NULL ) free(av->chromvar[i].name);
    }
    free(av->chromvar);
    for(i =0; i< av->n_var; ++i){
        free(av->var[i].orign);   
        free(av->var[i].var);   
    }
    free(av->var);
    kh_destroy(id, av->h);
    free(av);

}

#ifdef MAIN_VR
int main()
{
    uint32_t i, j, k; 
    all_var_t *vr;
    vr = variants_build("variants.txt");
    printf("\n");
    printf("chrom informations n_chrom = %u\n", vr->n_chrom);
    for(i=0; i < vr->n_chrom; ++i){
        chromvar_t *x = vr->chromvar;
        k = kh_get(id, vr->h, x[i].name);
        if(i != kh_val(vr->h, k)) printf("id %u != %u\n", i, kh_val(vr->h, k));
        printf(">>>%u\t%s\t%u\t%u\n", i, x[i].name, x[i].offset, x[i].len);
        for(j = x[i].offset; j <x[i].offset+x[i].len; ++j){
            var_t *v = vr->var+j;
            printf("%s\t%u\t%u\t%c\t%s\t%s\n", x[i].name, v->start, v->end, v->varclass, v->orign, v->var); 
        
        }
    }

    variants_destroy(vr);

    // 
    return 0;
}
#endif
