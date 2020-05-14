/*
 * =====================================================================================
 *
 *       Filename:  aln.c
 *
 *    Description:  O
 *
 *        Version:  1.0
 *        Created:  2017年10月12日 10时54分32秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (wq), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include "index.h"
#include "query.h"

static int test_seed_usage()
{
    fprintf(stderr, "generate seeds seq include variant\n");
    fprintf(stderr, "seed    [OPT]   <index.prefix>\n\n"); 

    fprintf(stderr, "       -h\n"); 
    fprintf(stderr, "\n");
    return 1;
}

void print_reads(const uint8_t*refseq, int var_start, int var_end, const char*alt, uint32_t ref_pos, bwtint_t *isa)
{

    int i, j;
    for(i = 0; i < 8; ++i){
        printf(">%u_%u_%d_%d_%s_%d\n", isa[ref_pos], ref_pos, var_start, var_end, alt, i);
        for(j = var_start-8; j < var_start; ++j){
            putchar("ACGT"[refseq[j]]); 
        }
        int l_alt = strlen(alt);
        if(strcmp(alt, "-") == 0) l_alt = 0;
        else printf("%s", alt);
        //printf("%s", alt);
        for(j=var_end ; j < var_end+12+i; ++j){
            putchar("ACGT"[refseq[j]]); 
        }
        putchar('\n'); 
    }

}
#define MAX_NAME 128
int gen_seed(const char *prefix)
{
    char fn_isvariant[MAX_NAME], fn_var[MAX_NAME];

    strncpy(fn_isvariant, prefix, MAX_NAME);strcat(fn_isvariant, ".vmap"); 
    strncpy(fn_var, prefix, MAX_NAME);strcat(fn_var, ".var");
 
    idx_t *idx = idx_restore(prefix);
    uint8_t *pac = idx->pac;
    bwt_t *bwt = idx->bwt;
    bwtint_t *isa = bwt->isa;
    bwtint_t *sa = bwt->sa; 
    bp_t *is_variant = bp_restore(fn_isvariant, 1);
    all_var_t *av = variants_restore(fn_var); 
    var_t *var = av->var;
    uint32_t *varseq = (uint32_t *)calloc(av->n_var*2, sizeof(uint32_t));        
    //fill variant seqs
    bwtint_t i, j; uint32_t vi,vseqi; 
    for(vi = 0; vi < av->n_var; ++vi){
        uint32_t pos = var[vi].end;// 
        vseqi = bp_rank(is_variant, isa[pos])-1;
        varseq[vseqi*2+1] = vi; 
        //fprintf(stderr, "variant pos %u-> bwt index %u\n", pos, isa[pos]);
    }
 
    
    for(i = 0; i < bwt->seq_len; ++i){
        if(bp_get(is_variant, i) == 0) continue;
        vseqi = bp_rank(is_variant, i)-1;
        uint32_t s, e;
        e = sa[i];
        //fprintf(stderr, "[debug]: variant pos %u\n", e);
        s = e-16;
        varseq[vseqi*2] = bns_extract_seq16(pac, s, e);/* variant preseq include variants */
    }

    for(i = 0; i < bwt->seq_len; ++i){
        uint8_t refseq[32] = {}; 
        if(bp_get(is_variant, i) == 0) continue; 
        uint32_t vseqi = bp_rank(is_variant, i) -1;
        uint32_t x = bwt_sa(bwt, i);
        //fprintf(stderr, "pos = %u variants [%u, %u]\n", x, var[varseq[vseqi*2+1]].start, var[varseq[vseqi*2+1]].end);
        if( x < 8|| x+20>bwt->seq_len) continue;//ignore reference head and tail
        for(j = 0; j < 28; ++j) {
            refseq[j] = __get_pac(pac, x-8+j); 
        }
        //for alt in variants
        var_t *v = var+varseq[vseqi*2+1];
        const char*alt = strtok(v->var, "\\"); 
        print_reads(refseq, v->start-(x-8), v->end-(x-8), alt, v->start, isa);
        while((alt =strtok(NULL, "\\")) != NULL){
            print_reads(refseq, v->start-(x-8), v->end-(x-8), alt, v->start, isa);
        } 
    } 
    idx_destroy(idx);
    bp_destroy(is_variant);
    variants_destroy(av);
    free(varseq); 
    return 0;
}

int gen_seeds_main(int argc, char *argv[])
{
    int c;
    while((c = getopt(argc, argv, "k:h"))>=0){
        switch(c){
            case 'h':
                return test_seed_usage();
            default:
                return 1;
              
        }
    
    }
    if (argc - optind  != 1 ){
        return test_seed_usage();
    }
    const char* fn_idx = argv[optind];
    gen_seed(fn_idx);   
    
    return 0;
}
int main(int argc, char *argv[])
{
    gen_seeds_main(argc, argv);
    return 1;
}
