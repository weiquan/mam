#ifndef __BWT_H
#define __BWT_H
/*
#include <stdint.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <zlib.h>
#include <time.h>
#include <assert.h>

//#include "utils.h"
//#include "khash.h"
//#include "index.h"
//#include "kvec.h"
//#include "ksort.h"
//#include "setFileName.h"
//#include "malloc_wrap.h"
//#define USE_MALLOC_WRAPPERS
#define MAX_NAME 128
// requirement: (OCC_INTERVAL%16 == 0)
#define OCC_INTERVAL 0x80

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef unsigned char ubyte_t;
#endif
typedef uint32_t bwtint_t;
typedef struct {
  int64_t offset;
  int32_t len;
  int32_t n_ambs;
  uint32_t gi;
  char *name, *anno;
} bntann1_t;

typedef struct {
  int64_t offset;
  int32_t len;
  char amb;
} bntamb1_t;

typedef struct {
  int64_t l_pac;
  int32_t n_seqs;
  uint32_t seed;
  bntann1_t *anns; // n_seqs elements
  int32_t n_holes;
  bntamb1_t *ambs; // n_holes elements
  FILE *fp_pac;
} bntseq_t;

typedef struct {
    bwtint_t primary; // S^{-1}(0), or the primary index of BWT
    bwtint_t L2[5]; // C(), cumulative count
    bwtint_t seq_len; // sequence length
    bwtint_t bwt_size; // size of bwt, about seq_len/4
    uint32_t *bwt; // BWT
    // occurance array, separated to two parts
    uint32_t cnt_table[256];
    // suffix array
    int sa_intv;
    bwtint_t n_sa;
    bwtint_t *sa;
    //inverse suffix array
    bwtint_t *isa;
} bwt_t;

typedef struct {
    size_t l, m, n_a;
    uint32_t *a;//n_a = m*(sizeof(uint32_t))
    int intv;
    int n_occ;
    size_t *occ;
} bp_t;
typedef struct{
    uint32_t seq, idx;
} ext_t;
typedef struct{

    //lkt_t *fastmap;//fastmap 12bp seq/
    
    bwt_t *bwt;//FM index/
    bntseq_t *bns;
    uint8_t *pac;// packed reference seq /
    int n_pmap;
    uint32_t *pmap;// preseq map / 
    int n_var;
    uint8_t *refseq;/2bp refseq at variant site/
    uint32_t *sai;
    
    bp_t *is_multiseeds;/is_multiseeds[i] = 1 if Reference[sa[i]:sa[i]+20] is multi location seeds /
    uint32_t n_isa2seq16;
    uint32_t* isa2seq16; // isa2seq16[i] = start of seq16s, i is occ of 1 in is_multiseeds before isa/
    uint32_t n_seq16s;
    uint32_t *seq16s;

    int n_tot;
    uint32_t *lext_idx;
    uint32_t *lext1_idx;
    uint32_t *rext_idx; 
    
    int n_rext, n_lext0, n_lext1;
    ext_t *rext0, *lext0;
    uint32_t *lext1;

    bp_t *is_var;// variants map
    uint32_t n_var;
    uint32_t *var_seqs;// var_seqs[2i] = variant pre seq; var_seqs[2i+1] = variant index 
    all_var_t *var;
} idx_t;

idx_t *idx_restore(const char *prefix, idx_t *idx);

void idx_destroy(idx_t *idx);

const int SA_INTV = 1;
*/
//int *global_stat_20seed0;
//int *global_stat_20seed1;
//int *global_stat_36seed;

//KHASH_MAP_INIT_STR(str, int)
//typedef kvec_t(uint32_t) vec_uint_t;
//typedef kvec_t(ext_t) vec_ext_t;
//#define __sort_lt_ext(a, b)((a).seq == (a).seq?(a).idx<(b).idx:(a).seq < (b).seq ) 
//#define __sort_lt_ext(a, b)((a).seq < (b).seq||((a).seq==(b).seq&&(a).idx < (b).idx )) 

//KSORT_INIT_GENERIC(uint32_t)
//KSORT_INIT(ext, ext_t, __sort_lt_ext)

#include <stdint.h>
#include "lookup.h"
#define size_t int64_t
#define SA_INTV 12
#define ISA_INTV 12


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

#define MAX_COUNT 800000
#define SWAP(type_t, a, b) do{type_t x=(a); (a)=(b); (b)=x;} while(0)
#define __set_bwt(bwt, l, c) ((bwt)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define __get_bwt(bwt, l) ((bwt)[(l)>>2]>>((~(l)&3)<<1)&3)
#define __set_bwt2(bwt, l, c) ((bwt)[(l)>>1] |= (c)<<((~(l)&1)<<2))
#define __get_bwt2(bwt, l) ((bwt)[(l)>>1] >> ((~(l)&1)<<2)&15)


#define __get_col(seq, i) (((seq)>>((i)*4))&0xF)
#define __set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define __get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
#define __set_seq16(word, l, c) ((word)|=((c)<<((~(l)&15)<<1)))
#define __get_seq16(word, l) ((word)>>((~(l)&15)<<1)&3)
static int index_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "index    [OPT]    <ref.fa>    <variants.txt>    <prefix>\n\n"); 
    fprintf(stderr, "         -h\n"); 
    fprintf(stderr, "\n");
    return 1;
}


#define bwt_bwt(b, k) ((b)->bwt[(k)/OCC_INTERVAL*12 + 4 + (k)%OCC_INTERVAL/16])

#define bwt_B0(b, k) (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)

#define bwt_occ_intv(b, k) ((b)->bwt + (k)/OCC_INTERVAL*12)

#define bwt_invPsi(bwt, k)                                              \
    (((k) == (bwt)->primary)? 0 :                                       \
     ((k) < (bwt)->primary)?                                            \
     (bwt)->L2[bwt_B0(bwt, k)] + bwt_occ(bwt, k, bwt_B0(bwt, k))        \
     : (bwt)->L2[bwt_B0(bwt, (k)-1)] + bwt_occ(bwt, k, bwt_B0(bwt, (k)-1)))


    //void bwt_dump_bwt(const char *fn, const bwt_t *bwt);
    //void bwt_dump_sa(const char *fn, const bwt_t *bwt);
    //void bwt_dump_isa(const char *fn, const bwt_t *bwt);
    
    bwt_t *bwt_restore_bwt(const char *fn);
    //void bwt_restore_sa(const char *fn, bwt_t *bwt);
    void bwt_restore_isa(const char *fn, int intv, bwt_t *bwt);

    void bwt_destroy(bwt_t *bwt);

    void bwt_bwtgen(const char *fn_pac, const char *fn_bwt); // from BWT-SW
    void bwt_cal_sa(bwt_t *bwt, int intv);
    void bwt_cal_isa(bwt_t *bwt);

    void bwt_bwtupdate_core(bwt_t *bwt);

    static inline bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c);

    static inline void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4]);
    bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k);

    // more efficient version of bwt_occ/bwt_occ4 for retrieving two close Occ values
    void bwt_gen_cnt_table(bwt_t *bwt);
    void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol);
    static inline void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4]);

    int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end);
    int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0);




//#ifndef LH3_UTILS_H
//#define LH3_UTILS_H

#include <stdint.h>
#include <stdio.h>

//#include "malloc_wrap.h"
//#ifdef __GNUC__
// Tell GCC to validate printf format string and args
//#define ATTRIBUTE(list) __attribute__ (list)
//#else
//#define ATTRIBUTE(list)
//#endif

//#ifdef DEBUG
//#define LOG(format, args...) fprintf(stderr, (format), ##args)
//#else
//#define LOG(format, args...)
//#endif

//#define err_fatal_simple(msg) _err_fatal_simple(__func__, msg)
//#define err_fatal_simple_core(msg) _err_fatal_simple_core(__func__, msg)

//#define fopen(fn, mode) err_fopen_core(__func__, fn, mode)
//#define xreopen(fn, mode, fp) err_xreopen_core(__func__, fn, mode, fp)
//#define xzopen(fn, mode) err_xzopen_core(__func__, fn, mode)

//#define xassert(cond, msg) if ((cond) == 0) _err_fatal_simple_core(__func__, msg)

//#define FSYNC_ON_FLUSH

//#include <stdio.h>
//#include <stdarg.h>
//#include <stdlib.h>
//#include <string.h>
//#include <errno.h>
//#ifdef FSYNC_ON_FLUSH
//#include <sys/types.h>
//#include <sys/stat.h>
//#include <unistd.h>
//#endif
//#include <sys/time.h>

#include "ksort.h"
//#define pair64_lt(a, b) ((a).x < (b).x || ((a).x == (b).x && (a).y < (b).y))
//KSORT_INIT(128, pair64_t, pair64_lt)
//KSORT_INIT(64,  uint64_t, ks_lt_generic)


//void bns_destroy(bntseq_t *bns);
uint32_t bns_extract_seq16(uint8_t *pac, uint32_t start, uint32_t end);

uint32_t bns_extract_seq16(uint8_t *pac, uint32_t start, uint32_t end)
{
    
    uint32_t c;
    uint32_t i;
    uint32_t r = 0;
    assert(start+16 >= end);
    //fprintf(stderr, "\n[seq]%u-%u\n", start, end);
    for(i = start; i < end; ++i){

        //c = pac[i>>2] >> ((3-(i&3))<<1);
        //c &= 3;

        //fprintf(stderr, "%c\t", "acgt"[c]);
        c = __get_pac(pac, i);
        //__set_seq16(r, i-start, c);
        r |= (c<<((15+start-i)<<1));
    }
    //fprintf(stderr, "\n");
    return r;
}


//===============================================================
bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[1024];
	FILE *fp;
	//bntseq_t *bns;
	long long xx;
	int i;
	bntseq_t *bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = fopen(ann_filename, "r");
		fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			fscanf(fp, "%u%s", &p->gi, str);
			p->name = strdup(str);
			// read fasta comments 
			while ((c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			*q = 0;
			if (q - str > 1) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			p->offset = xx;
		}
		fclose(fp);
	}
	{ // read .amb
		int64_t l_pac;
		int32_t n_seqs;
		fp = fopen(amb_filename, "r");
		fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		l_pac = xx;
		
		bns->ambs = (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t));
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			p->offset = xx;
			p->amb = str[0];
		}
		fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = fopen(pac_filename, "rb");
	}
	return bns;
}

bntseq_t *bns_restore(const char *prefix)
{  
	char ann_filename[1024], amb_filename[1024], pac_filename[1024];
	strcat(strcpy(ann_filename, prefix), ".ann");
	strcat(strcpy(amb_filename, prefix), ".amb");
	strcat(strcpy(pac_filename, prefix), ".pac");
	return bns_restore_core(ann_filename, amb_filename, pac_filename);
}
int bns_coor_pac2real(const bntseq_t *bns, int64_t pos, int len, int32_t *real_seq)
{
	uint32_t pac_coor = pos;
    int left, mid, right, nn;
	if (pac_coor >= bns->l_pac) {
        printf("bns_coor_pac2real bug! Coordinate is longer than sequence (%lld>=%lld).", pac_coor, bns->l_pac);
        exit(1);
    }
    // binary search for the sequence ID. Note that this is a bit different from the following one...
	left = 0; mid = 0; right = bns->n_seqs;
	while (left < right) {
		mid = (left + right) >> 1;
		if (pac_coor >= bns->anns[mid].offset) {
			if (mid == bns->n_seqs - 1) break;
			if (pac_coor < bns->anns[mid+1].offset) break;
			left = mid + 1;
		} else right = mid;
	}
	*real_seq = mid;
	// binary search for holes
	left = 0; right = bns->n_holes; nn = 0;
	while (left < right) {
		int64_t mid = (left + right) >> 1;
		if (pac_coor >= bns->ambs[mid].offset + bns->ambs[mid].len) left = mid + 1;
		else if (pac_coor + len <= bns->ambs[mid].offset) right = mid;
		else { // overlap
			if (pac_coor >= bns->ambs[mid].offset) {
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pac_coor + len?
					bns->ambs[mid].offset + bns->ambs[mid].len - pac_coor : len;
			} else {
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pac_coor + len?
					bns->ambs[mid].len : len - (bns->ambs[mid].offset - pac_coor);
			}
			break;
		}
	}
	return nn;
}
void bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else {
		int i;
		if (bns->fp_pac) fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i) {
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}

idx_t *idx_restore(const char *prefix,idx_t *idx)
{
    //idx_t *idx = (idx_t *)calloc(1, sizeof(idx_t)); 
    
    char fn0[MAX_NAME];
  

    /* restore pac and ann */ 

//idx->bns = bns_restore(prefix);  
// 权威改写上面的函数，
// ++++++++++++++++++++++++++++++
    //bntseq_t *bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
    idx->bns = bns_restore(prefix);

//FILE *fp;
//fp = fopen("./XXXX","r");
//idx->bns->fp_pac = fp;
//+++++++++++++++++++++++++++++++++++++++++++
    idx->pac = (uint8_t *)calloc(idx->bns->l_pac/4+1, sizeof(uint8_t));

    fread(idx->pac, sizeof(uint8_t), idx->bns->l_pac/4+1, idx->bns->fp_pac);
    //err_fclose(fp_pac);
    /* restore FM-index */

    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".bwt");
    printf("[log]: restore bwt...\n");
    
    idx->bwt = bwt_restore_bwt(fn0); 
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".sa");
    //bwt_restore_sa(fn0, idx->bwt); 
    printf("[log]: restore sa...\n");
    bwt_restore_sa(fn0, SA_INTV, idx->bwt); 
    //if(idx->bwt->sa_intv != 1) test_bwt_sa(idx->bwt);
    printf("[log]: restore isa...\n");
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".isa");
    //bwt_restore_isa(fn0, idx->bwt); 
    bwt_restore_isa(fn0, ISA_INTV, idx->bwt); 

    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".lkt");
    idx->fastmap = lkt_restore(fn0); 
    /* restore snp-aware index */ 
    /*
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".vseq");
    FILE *fp_varseq = fopen(fn0, "r");
    fread(&idx->n_pmap, sizeof(uint32_t), 1, fp_varseq);
    idx->pmap = (uint32_t *)calloc(idx->n_pmap, sizeof(uint32_t));
    fread(idx->pmap, sizeof(uint32_t), idx->n_pmap, fp_varseq);
    fread(&idx->n_var, sizeof(uint32_t), 1, fp_varseq);
    idx->sai = (uint32_t *) calloc(idx->n_var, sizeof(uint32_t));
    idx->refseq = (uint8_t *) calloc(idx->n_var, sizeof(uint8_t));
    fread(idx->sai, sizeof(uint32_t), idx->n_var, fp_varseq);
    fread(idx->refseq, sizeof(uint8_t), idx->n_var, fp_varseq);
    err_fclose(fp_varseq);
    */
    //strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".repeat");
    //idx->is_multiseeds = bp_restore(fn0, 1);
    //restore_extend_idx(idx, prefix);
/*
    FILE *fp;   
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".idxseq16");
    fp = fopen(fn0, "r");
    fread(&idx->n_isa2seq16, sizeof(int), 1, fp);
    idx->isa2seq16 = (uint32_t*) calloc(idx->n_isa2seq16, sizeof(uint32_t));
    fread(idx->isa2seq16, sizeof(uint32_t), idx->n_isa2seq16, fp);
    fclose(fp);
    strncpy(fn0, prefix, MAX_NAME);strcat(fn0, ".seq16");
    fp = fopen(fn0, "r");
    idx->n_seq16s = idx->isa2seq16[idx->n_isa2seq16-1];
    idx->seq16s = (uint32_t*) calloc(idx->n_seq16s, sizeof(uint32_t));
    fread(idx->seq16s, sizeof(uint32_t), idx->n_seq16s, fp); 
    fclose(fp); 
*/
    return idx;
}
void idx_destroy(idx_t *idx)
{
    if(idx->fastmap != NULL ) lkt_destroy(idx->fastmap);
    bwt_destroy(idx->bwt);
    bns_destroy(idx->bns);
    if(idx->pac != NULL ) free(idx->pac);
    //free(idx->pmap);
    //free(idx->refseq);
    //free(idx->sai);

    //free(idx->lext0);
    //free(idx->lext1);
    //free(idx->rext0);    
    //free(idx->lext_idx);
    //free(idx->rext_idx);
    //free(idx->lext1_idx); 
    //bp_destroy(idx->is_multiseeds);
    //free(idx->isa2seq16);
    //free(idx->seq16s);
    //bp_destroy(idx->is_var);
    //free(idx->var_seqs);
    //variants_destroy(idx->var);
    //free(idx);
}
/*  
bwtint_t bwt_get_idx(bwt_t *bwt, uint32_t pos)
{
    if(pos < bwt->seq_len ){
        if(pos%2 == 0 ) return bwt->isa[pos/2]; 
        else {
            uint32_t idx = bwt->isa[(pos+1)/2];
            idx = bwt_invPsi(bwt, idx); 
            return idx;
        }
    } else{
        return 0;
    }

}
*/
/*  
bwtint_t bwt_get_idx(bwt_t *bwt, uint32_t pos)
{
    int isa_intv = bwt->isa_intv;
    if(pos < bwt->seq_len ){
        if(pos%isa_intv == 0 ) return bwt->isa[pos/isa_intv]; 
        else {
            uint32_t idx = bwt->isa[(pos+1)/isa_intv];
            idx = bwt_invPsi(bwt, idx); 
            return idx;
        }
    } else{
        return 0;
    }

}
*/
bwtint_t bwt_get_idx(bwt_t *bwt, uint32_t pos)
{
    int i, isa_intv = bwt->isa_intv;
    bwtint_t idx;
    if(pos%isa_intv == 0) {
        idx = bwt->isa[pos/isa_intv];
    } else {
        idx = bwt->isa[pos/isa_intv+1];
        for(i = 0; i < isa_intv - pos%isa_intv; ++i) {
            idx = bwt_invPsi(bwt, idx); 
        }
    }
    return idx;  
}



void bwt_restore_isa1(const char *fn, bwt_t *bwt)
{
    FILE *fp;
    int intv = bwt->isa_intv;
    fp = fopen(fn, "rb");
    
    fread(&bwt->n_isa, sizeof(bwtint_t), 1, fp);
    bwt->n_isa = ((bwt->n_isa)+intv-1)/intv;
    bwt->isa  = (bwtint_t *)calloc ( bwt->n_isa, sizeof(bwtint_t) );
    
    if ( bwt->isa==NULL ) {
        fprintf ( stderr, "\ndynamic memory allocation failed\n" );
        exit (EXIT_FAILURE);
    }
    int max_num = 0xFFFF, i, n;
    uint32_t tot = 0;
    bwtint_t *buf = calloc(max_num, sizeof(bwtint_t));
printf("\ndynamic memory allocation\n" );
    while( (n = fread(buf, sizeof(bwtint_t), max_num, fp)) >0 ){
        for(i = 0; i < n; ++i, ++tot){
            if(tot %intv != 0) continue;
            bwt->isa[tot/intv] = buf[i]; 
        } 
    
    }
printf("\ndynamic memory allocation\n" );
    free(buf);
    fclose(fp);
}
void bwt_restore_isa(const char *fn, int intv, bwt_t *bwt)
{

  if(intv == 1) {
    bwt->isa_intv = 1;
    bwt_restore_isa1(fn, bwt);
    return; 
  } else if(intv < 1){
    printf("error!\n");
    exit(1); 
  } 


  char prefix[256] = {}, fn_intv_isa[256];
  strncpy(prefix, fn, strlen(fn)-4);
  sprintf(fn_intv_isa, "%s_%d.isa", prefix, intv);

 
  FILE *fp_intv_isa = NULL;
  if(fp_intv_isa = fopen(fn_intv_isa, "r")) {
    fread(&bwt->isa_intv, sizeof(uint32_t), 1, fp_intv_isa); 
    fread(&bwt->n_isa, sizeof(uint32_t), 1, fp_intv_isa);
printf("isa_intv = %d, n_isa = %u\n", bwt->isa_intv, bwt->n_isa); 
    bwt->isa = (bwtint_t*)calloc(bwt->n_isa, sizeof(bwtint_t));
    fread(bwt->isa, sizeof(bwtint_t), bwt->n_isa, fp_intv_isa); 
    fclose(fp_intv_isa); 
  } else {
    bwt->isa_intv = intv;
    //bwt->n_isa = (bwt->seq_len+intv-1)/intv + 1; 
    //bwt->isa = (bwtint_t*)calloc(bwt->n_isa, sizeof(bwtint_t));
    
    fp_intv_isa = fopen(fn_intv_isa, "w");
    bwt_restore_isa1(fn, bwt);
    fwrite(&bwt->isa_intv, sizeof(uint32_t), 1, fp_intv_isa); 
    fwrite(&bwt->n_isa, sizeof(uint32_t), 1, fp_intv_isa);
    fwrite(bwt->isa, sizeof(bwtint_t), bwt->n_isa, fp_intv_isa); 
    fclose(fp_intv_isa);
   
  
  }
  return; 
}
void bwt_restore_sa_1(const char *fn, bwt_t *bwt)
{
  char skipped[256];
  FILE *fp;
  bwtint_t primary;

  fp = fopen(fn, "rb");
  fread(&primary, sizeof(bwtint_t), 1, fp);
  //xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
  fread(skipped, sizeof(bwtint_t), 4, fp); // skip
  fread(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
  fread(&primary, sizeof(bwtint_t), 1, fp);
  //xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

  bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
  bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
  bwt->sa[0] = (bwtint_t)-1;

  fread(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
  fclose(fp);
}
void bwt_restore_sa(const char *fn, int sa_intv, bwt_t *bwt)
{
  if(sa_intv == 1) {
    bwt_restore_sa_1(fn, bwt);
    return; 
  } else if(sa_intv < 1){
    printf("error!\n");
    exit(1); 
  } 

  //bwt_restore_sa_1(fn, bwt);
  bwt->n_pos = (bwt->seq_len+sa_intv-1)/sa_intv + 1;
  bwt->n_mod = (bwt->seq_len+63)/64 + 1;
  bwt->n_flg = (bwt->seq_len+63)/64 + 1;
 
  
  bwt->sa_pos = (bwtint_t*)calloc(bwt->n_pos, sizeof(bwtint_t));
  bwt->sa_mod = (bwtint_t*)calloc(bwt->n_mod, sizeof(bwtint_t));
  bwt->sa_flg = (uint64_t *)calloc(bwt->n_flg, sizeof(uint64_t));


  bwtint_t *sa_mod, *sa_pos;
  uint64_t *sa_flg;
  sa_mod = bwt->sa_mod; // #pos before idx
  sa_pos = bwt->sa_pos;
  sa_flg = bwt->sa_flg;
  char skipped[256], fn_mod[256], fn_flg[256], fn_pos[256];
  FILE *fp;
  bwtint_t primary;
  char prefix[256] = {};
  strncpy(prefix, fn, strlen(fn)-3);
  sprintf(fn_pos, "%s_%d.sa", prefix, sa_intv);
  sprintf(fn_mod, "%s_%d.mod", prefix, sa_intv);
  sprintf(fn_flg, "%s_%d.flg", prefix, sa_intv);
  FILE *fp_pos = NULL, *fp_mod = NULL,  *fp_flg = NULL;
  int open_flg = 0;
  if(fp_pos = fopen(fn_pos, "r")) {
    open_flg += 1; 
  }
  if(fp_mod = fopen(fn_mod, "r")) {
    open_flg += 1; 
  }
  if(fp_flg = fopen(fn_flg, "r")) {
    open_flg += 1; 
  }
  if(open_flg == 3) {
    //read pos, mod, flg
    bwt->sa_intv = sa_intv;     
    fread(&bwt->n_pos, sizeof(uint32_t), 1, fp_pos);
    fread(bwt->sa_pos, sizeof(bwtint_t), bwt->n_pos, fp_pos); 
    fread(&bwt->n_mod, sizeof(uint32_t), 1, fp_mod);
    fread(bwt->sa_mod, sizeof(bwtint_t), bwt->n_mod, fp_mod); 
    fread(&bwt->n_flg, sizeof(uint32_t), 1, fp_flg);
    fread(bwt->sa_flg, sizeof(uint64_t), bwt->n_flg, fp_flg); 
    fclose(fp_pos); 
    fclose(fp_mod); 
    fclose(fp_flg); 
    return;
  }


 
  
  fp = fopen(fn, "rb");
  fread(&primary, sizeof(bwtint_t), 1, fp);
  //xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
  fread(skipped, sizeof(bwtint_t), 4, fp); // skip
  fread(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
  fread(&primary, sizeof(bwtint_t), 1, fp);
  //xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");
  if(bwt->sa_intv != 1) {
    printf("%u, error!\n", __LINE__);
    exit(1); 
  
  }
  bwtint_t pos, pos_num = 0;
#define __set_flg(flg, p) ((flg)[p/64] |= (uint64_t)1 <<(63-p%64))
#define __get_flg(flg, p) (((flg)[p/64] >> (63-p%64)) & 1)
  uint32_t idx = 1;
  while(fread(&pos, sizeof(bwtint_t), 1, fp)){
    uint32_t idx_q, idx_r;
    if(idx % 64 == 0) { 
        //idx_q = idx/64;
        //idx_r = idx%64;
        sa_mod[idx/64] = pos_num;        
       
        
    }
    if(pos % sa_intv == 0) {
        sa_pos[pos_num] = pos;  
        idx_q = idx/64;
        idx_r = idx%64;
        sa_flg[idx_q] |= (uint64_t)1 << (63-idx_r);
        //set_flg(sa_flg, idx);
        pos_num++;
    } 
    idx++; 
  }
  fclose(fp);
  
  bwt->sa_intv = sa_intv;     
  open_flg = 0;
  if(fp_pos = fopen(fn_pos, "w")) {
    open_flg += 1; 
  }
  if(fp_mod = fopen(fn_mod, "w")) {
    open_flg += 1; 
  }
  if(fp_flg = fopen(fn_flg, "w")) {
    open_flg += 1; 
  }
  if(open_flg == 3) {
    //write sa, mod, flg
    fwrite(&bwt->n_pos, sizeof(uint32_t), 1, fp_pos);
    fwrite(bwt->sa_pos, sizeof(bwtint_t), bwt->n_pos, fp_pos); 
    fwrite(&bwt->n_mod, sizeof(uint32_t), 1, fp_mod);
    fwrite(bwt->sa_mod, sizeof(bwtint_t), bwt->n_mod, fp_mod); 
    fwrite(&bwt->n_flg, sizeof(uint32_t), 1, fp_flg);
    fwrite(bwt->sa_flg, sizeof(uint64_t), bwt->n_flg, fp_flg); 
     
      

    fclose(fp_pos); 
    fclose(fp_mod); 
    fclose(fp_flg); 

    return;
  } else {
    printf("%u, error !\n", __LINE__); 
    exit(1);
  
  }

  

  /*  
  bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
  bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
  bwt->sa[0] = (bwtint_t)-1;
  */
  //fread(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
  //fclose(fp);
}
int bwt_lcp(uint8_t *pac, bwt_t *bwt, uint32_t idx0, uint32_t idx1)
{
    int len = 0;
    uint8_t c0, c1;
    c0 = __get_pac(pac, bwt_sa(bwt, idx0));
    c1 = __get_pac(pac, bwt_sa(bwt, idx1));
    while(c0 == c1) {
        len++;  
        c0 = __get_pac(pac, bwt_sa(bwt, idx0) + len);
        c1 = __get_pac(pac, bwt_sa(bwt, idx1) + len);
    }
    return len;
}
bwt_t *bwt_restore_bwt(const char *fn)
{
  bwt_t *bwt;
  FILE *fp;

  bwt = (bwt_t*)calloc(1, sizeof(bwt_t));

  fp = fopen(fn, "rb");
  fseek(fp, 0, SEEK_END);
  bwt->bwt_size = (ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
//+++++++++++++++++++++++++++
//bwt->bwt_size =100;
//printf("bwt->bwt_size = %d\n", bwt->bwt_size);
//----------------------------
  bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
  

  fseek(fp, 0, SEEK_SET);
  fread(&bwt->primary, sizeof(bwtint_t), 1, fp);
  fread(bwt->L2+1, sizeof(bwtint_t), 4, fp);
  fread(bwt->bwt, 4, bwt->bwt_size, fp);
  bwt->seq_len = bwt->L2[4];
  fclose(fp);
  bwt_gen_cnt_table(bwt);

  return bwt;
}

void bwt_destroy(bwt_t *bwt)
{
  if(bwt == 0) return;
  if(bwt->isa) free(bwt->isa);
  if(bwt->sa) free(bwt->sa);
  else {
    free(bwt->sa_pos);
    free(bwt->sa_mod);
    free(bwt->sa_flg);
  } 
  free(bwt->bwt);
  free(bwt);
}

void bwt_gen_cnt_table(bwt_t *bwt)
{
    int i, j;
    for (i = 0; i != 256; ++i) {
        uint32_t x = 0;
        for (j = 0; j != 4; ++j)
            x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
        bwt->cnt_table[i] = x;
    }
}
/*
void bwt_gen_cnt_table(bwt_t *bwt)
{
  int i, j;
  for (i = 0; i != 256; ++i) {
    uint32_t x = 0;
    for (j = 0; j != 4; ++j)
      x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
    bwt->cnt_table[i] = x;
  }
}
*/
// bwt->bwt and bwt->occ must be precalculated
void bwt_cal_sa(bwt_t *bwt, int intv)
{
  bwtint_t isa, sa, i; // S(isa) = sa

  //xassert(bwt->bwt, "bwt_t::bwt is not initialized.");

  if (bwt->sa) free(bwt->sa);
  bwt->sa_intv = intv;
  bwt->n_sa = (bwt->seq_len + intv) / intv;
  bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
  // calculate SA value
  isa = 0; sa = bwt->seq_len;
  for (i = 0; i < bwt->seq_len; ++i) {
    if (isa % intv == 0) bwt->sa[isa/intv] = sa;
    --sa;
    isa = bwt_invPsi(bwt, isa);
  }
  if (isa % intv == 0) bwt->sa[isa/intv] = sa;
  bwt->sa[0] = (bwtint_t)-1; // before this line, bwt->sa[0] = bwt->seq_len
}
void bwt_cal_isa(bwt_t *bwt)
{
        bwtint_t i;
        bwtint_t *sa = bwt->sa;
        bwtint_t *isa = (bwtint_t *)calloc(bwt->n_sa, sizeof(bwtint_t));

        isa[bwt->seq_len] = 0;
        for(i=1; i < bwt->n_sa; ++i){
            isa[sa[i]] = i;
        } 
        bwt->isa = isa; 
}
bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t idx)
{
    if(idx > bwt->seq_len) {
        printf("%u, (idx = %u )> (seq_len = %u) !!!\n", __LINE__, idx, bwt->seq_len);
        exit(1);
    }
    if(bwt->sa_intv == 1) {
        return bwt->sa[idx]; 
    } else if(bwt->sa_intv < 1){
        printf("%u, error !!!\n", __LINE__);
        exit(1); 
    } 
    bwtint_t *sa_mod = bwt->sa_mod, *sa_pos = bwt->sa_pos;
    uint64_t *sa_flg = bwt->sa_flg;
    uint32_t n_mod = bwt->n_mod, n_pos = bwt->n_pos, n_flg = bwt->n_flg; 
    
    bwtint_t sa = 0;
    //printf("%u\n", __LINE__); 
    int flg = (sa_flg[idx/64] >> (63-idx%64))&1;       
    //while (__get_flg(sa_flg, idx) == 0) {
    while (flg == 0) {

//#define __get_flg(flg, p) (((flg)[p/64] >> (63-p%64)) & 1)
    //printf("idx = %u, idx_q = %u, idx_r = %u, n_flg = %u\n", idx, idx/64, idx%64, bwt->n_flg);
    //printf("idx_flg = %d\n", __get_flg(sa_flg, idx));

        ++sa;
        idx = bwt_invPsi(bwt, idx);
        flg = (sa_flg[idx/64] >> (63-idx%64))&1;       
    }
    if(sa > bwt->sa_intv) {
        printf("%u, sa = %u, bwt->sa_intv = %u, error!!!\n", __LINE__, sa, bwt->sa_intv);
        exit(1); 
    }
    //count 1 before idx in flg
    
    //printf("%u\n", __LINE__); 
    uint32_t idx_q = idx/64; 
    uint32_t idx_r = idx%64; 
    uint32_t n = sa_mod[idx_q];
    int i;
    for(i = 0; i <= idx_r; ++i) {
//printf("i = %u, idx_q = %u, idx_r = %u, n_flg = %u, n_mod = %u\n", i, idx_q, idx_r, bwt->n_flg, bwt->n_mod);

//printf("flg[idx_q] = %d\n", sa_flg[idx_q]);

//#define __get_flg(flg, p) (((flg)[p/64] >> (63-p%64)) & 1)
       int flg = (sa_flg[idx_q] >> (63-i))&1;       
//printf("flg = %d\n", flg);
//printf("get_flg = %d\n", __get_flg(sa_flg, idx_q*64+i));

        //n += __get_flg(sa_flg, idx_q*64+i); 
        n += flg; 
    }
    return sa + bwt->sa_pos[n-1];
}
int test_bwt_sa(bwt_t *bwt)
{
    int64_t k, pos;

    int64_t i, n = 0;
    for(i = 1; i < bwt->seq_len; ++i){
        
        if(i % 1000000 == 0) printf("%u pos checked!\n", i); 
        uint32_t pos0, pos1; 
        pos0 = bwt->sa[i];  

        if(1) {  
      
            pos1 = bwt_sa(bwt, i);
            if(pos0 != pos1) {
                printf("%u != %u, idx = %u\n", pos0, pos1, i); 
                exit(1);
            } 
        }
        /*   
        if(i %64 == 0) {
            if(n != bwt->sa_mod[i/64]) {
                printf("%u, %u != %u, idx = %u\n", __LINE__, n, bwt->sa_mod[i/64], i); 
                exit(1);
            }
        }
        

        uint32_t pos0 = bwt->sa[i];
        if(pos0 % SA_INTV == 0) {
            int flg = (bwt->sa_flg[i/64] >> (63-i%64))&1;       
            if(flg == 0) {
                printf("%u, flg = %u, idx = %u\n", __LINE__, flg, i); 
                //exit(1);
            }
            if(pos0 != bwt->sa_pos[n]) {
                printf("%u, %u != %u, idx = %u\n", __LINE__, pos0, bwt->sa_pos[n], i); 
                exit(1);
            }
            ++n;
        }
        */
    }
    /*  
    for(pos = bwt->seq_len, k = 0; pos >= 0; --pos){
        if(pos % 1 == 0 ){
            printf("%u, pos = %u, k = %u\n", __LINE__, pos, k); 
        }
        if(k != 0) {
            uint32_t pos1 = bwt_sa(bwt, k);
            if(pos1 != pos) {
                printf("%u != %u, k = %u\n", pos, pos1, k); 
                exit(1);
            } 
        }
        k = bwt_invPsi(bwt, k);

    }
    */
    return 1;
}


bwtint_t bwt_sa0(const bwt_t *bwt, bwtint_t k)
{
  bwtint_t sa = 0;
  while (k % bwt->sa_intv != 0) {
    ++sa;
    k = bwt_invPsi(bwt, k);
  }
  /* without setting bwt->sa[0] = -1, the following line should be
     changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
  //printf("pos :%d\n", sa + bwt->sa[k/bwt->sa_intv]);
    return sa + bwt->sa[k/bwt->sa_intv];
}

static inline int __occ_aux(uint64_t y, int c)
{
  // reduce nucleotide counting to bits counting
  y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
  // count the number of 1s in y
  y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
  return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

static inline bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
  bwtint_t n, l, j;
  uint32_t *p;

  if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
  if (k == (bwtint_t)(-1)) return 0;
  if (k >= bwt->primary) --k; // because $ is not in bwt

  // retrieve Occ at k/OCC_INTERVAL
  n = (p = bwt_occ_intv(bwt, k))[c];
  p += 4; // jump to the start of the first BWT cell

  // calculate Occ up to the last k/32
  j = k >> 5 << 5;
  for (l = k/OCC_INTERVAL*OCC_INTERVAL; l < j; l += 32, p += 2)
    n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);

  // calculate Occ
  n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
  if (c == 0) n -= ~k&31; // corrected for the masked bits

  return n;
}

// an analogy to bwt_occ() but more efficient, requiring k <= l
void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol)
{
  bwtint_t _k, _l;
  if (k == l) {
    *ok = *ol = bwt_occ(bwt, k, c);
    return;
  }
  _k = (k >= bwt->primary)? k-1 : k;
  _l = (l >= bwt->primary)? l-1 : l;
  if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
    *ok = bwt_occ(bwt, k, c);
    *ol = bwt_occ(bwt, l, c);
  } else {
    bwtint_t m, n, i, j;
    uint32_t *p;
    if (k >= bwt->primary) --k;
    if (l >= bwt->primary) --l;
    n = (p = bwt_occ_intv(bwt, k))[c];
    p += 4;
    // calculate *ok
    j = k >> 5 << 5;
    for (i = k/OCC_INTERVAL*OCC_INTERVAL; i < j; i += 32, p += 2)
      n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
    m = n;
    n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
    if (c == 0) n -= ~k&31; // corrected for the masked bits
    *ok = n;
    // calculate *ol
    j = l >> 5 << 5;
    for (; i < j; i += 32, p += 2)
      m += __occ_aux((uint64_t)p[0]<<32 | p[1], c);
    m += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~l&31)<<1)) - 1), c);
    if (c == 0) m -= ~l&31; // corrected for the masked bits
    *ol = m;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define __occ_aux4(bwt, b)                      \
  ((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]   \
   + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

static inline void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
  bwtint_t l, j, x;
  uint32_t *p;
  if (k == (bwtint_t)(-1)) {
    memset(cnt, 0, 4 * sizeof(bwtint_t));
    return;
  }
  if (k >= bwt->primary) --k; // because $ is not in bwt
  p = bwt_occ_intv(bwt, k);
  memcpy(cnt, p, 16);
  p += 4;
  j = k >> 4 << 4;
  for (l = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; l < j; l += 16, ++p)
    x += __occ_aux4(bwt, *p);
  x += __occ_aux4(bwt, *p & ~((1U<<((~k&15)<<1)) - 1)) - (~k&15);
  cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

// an analogy to bwt_occ4() but more efficient, requiring k <= l
static inline void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
  bwtint_t _k, _l;
  if (k == l) {
    bwt_occ4(bwt, k, cntk);
    memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
    return;
  }
  _k = (k >= bwt->primary)? k-1 : k;
  _l = (l >= bwt->primary)? l-1 : l;
  if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
    bwt_occ4(bwt, k, cntk);
    bwt_occ4(bwt, l, cntl);
  } else {
    bwtint_t i, j, x, y;
    uint32_t *p;
    int cl[4];
    if (k >= bwt->primary) --k; // because $ is not in bwt
    if (l >= bwt->primary) --l;
    cl[0] = cl[1] = cl[2] = cl[3] = 0;
    p = bwt_occ_intv(bwt, k);
    memcpy(cntk, p, 4 * sizeof(bwtint_t));
    p += 4;
    // prepare cntk[]
    j = k >> 4 << 4;
    for (i = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; i < j; i += 16, ++p)
      x += __occ_aux4(bwt, *p);
    y = x;
    x += __occ_aux4(bwt, *p & ~((1U<<((~k&15)<<1)) - 1)) - (~k&15);
    // calculate cntl[] and finalize cntk[]
    j = l >> 4 << 4;
    for (; i < j; i += 16, ++p) y += __occ_aux4(bwt, *p);
    y += __occ_aux4(bwt, *p & ~((1U<<((~l&15)<<1)) - 1)) - (~l&15);
    memcpy(cntl, cntk, 16);
    cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
    cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
  }
}

int bwt_match_exact(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end)
{
  bwtint_t k, l, ok, ol;
  int i;
  k = 0; l = bwt->seq_len;
  for (i = len - 1; i >= 0; --i) {
    ubyte_t c = str[i];
    if (c > 3) return 0; // no match
//printf("bwt_2occ(bwt, k - 1, l, c, &ok, &ol) = %s\n", " ok1" );
    bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
//printf("bwt_2occ(bwt, k - 1, l, c, &ok, &ol) = %s\n", " ok2" );    
    k = bwt->L2[c] + ok + 1;
    l = bwt->L2[c] + ol;
    if (k > l) break; // no match
  }
  if (k > l) return 0; // no match
  if (sa_begin) *sa_begin = k;
  if (sa_end)   *sa_end = l;
  return l - k + 1;
}
int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0)
{
   int i;
    int ret =1;
  bwtint_t k, l, ok, ol;
  k = *k0; l = *l0;
  for (i = len - 1; i >= 0; --i) {
    ubyte_t c = str[i];
    if (c > 3){             
        ret = 0; // there is an N here. no match
        break;
    }
    bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
    k = bwt->L2[c] + ok + 1;
    l = bwt->L2[c] + ol;
    if (k > l) {
            ret = 0; // no match
          break;
        }
    }
   if(ret ){
      *k0 = k; *l0 = l;
        return l - k + 1;
    }else{
        return ret;
    }


}




#endif

