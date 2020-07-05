/*
 * =====================================================================================
 *
 *       Filename:  alnse.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/01/2013 02:08:59 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT
 *
 * =====================================================================================
 */

//cbwt should use forward-search
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <pthread.h>
#include "index.h"
#include "ksort.h"
#include "khash.h"
#include "ksw.h"
#include "kvec.h"
#include "aln.h"
#include "sam.h" 
typedef struct{
  int n;
  uint32_t bg, ed;
} chain_t;
typedef struct{
  int sc;
  uint32_t tb, te;
  int qb, qe;
  int strand;
} aln_t;

#define __sai_lt(a, b) ((a).ep -(a).sp < (b).ep - (b).sp)
//KHASH_MAP_INIT_INT(32, unsigned char)
//KSORT_INIT_GENERIC(uint32_t)
KSORT_INIT(sai, sai_t, __sai_lt)	


#define __chain_lg(a, b) ((a).n == (b).n? (a).ed-(a).bg<(b).ed-(b).bg:(a).n>(b).n)
KSORT_INIT(chain, chain_t, __chain_lg)


#define __MAX(a, b) ((a)>(b)?(a):(b))
#define __MIN(a,b) ((a)<(b)?(a):(b) )



int g_n_seqs;
pthread_rwlock_t rwlock;



#define LITTLE_SEED_MAX 25
#define PRINT_TIME() \
do{\
    fprintf(stderr, "[aln_main]aln reads ...");\
    fprintf(stderr, "%.2f sec\n", (float)(clock()-tot_aln)/CLOCKS_PER_SEC);\
} while(0)





void alnse_core1(idx_t *index, int n, query_t *multi_query, aln_opt_t* aln_opt);

static inline int __lt(uint32_t *a, uint32_t *b)
{
    return *a <*b ; 

}


//#define OVERLAP_INTV 5
//void alnse_seed_overlap(index_t *index, const uint8_t *seq,uint32_t l_seq, int l_seed, candidate_t *p_can)
#define candi_t kvec_t(uint32_t)

void alnse_seed_overlap(idx_t *index, uint32_t l_seq, const uint8_t *seq, aln_opt_t *opt, candi_t *candi)
{
    //int m;
    int l_seed = opt->l_seed;
    int n;       

    int l_init = l_seed; 

    int OVERLAP_INTV = opt->l_overlap;
    bwt_t *bwt = index->bwt;
    //rbwt_t *rbwt0 = index->rbwt2->rbwt0; //forward


    int seed_start, seed_end; 
    for(seed_start= 0; (unsigned int)seed_start < l_seq-l_init+1; ++seed_start )
    {   
      if (seed_start % OVERLAP_INTV != 0) continue;
      seed_end = seed_start+l_init;
      uint32_t k, l, p, i,n, init_bg = seed_start, init_ed = seed_end; 
      uint8_t c;
      /* non-snp seeding*/
      /*  
      k = 1; //k=0???
      l = bwt->seq_len;
      LKT_lookup_sa(lkt, seq, seed_end -l_lkt+1, seed_end, &k, &l);
      */
      k = 1; //k=0???
      l = bwt->seq_len;
      uint32_t iter = lkt_seq2LktItem(seq, seed_start+8, seed_end-1);
#ifdef DEBUG
      uint32_t __k = 1, __l = bwt->seq_len, __n, __k0, __l0;
      __n = bwt_match_exact_alt(index->bwt, 12, seq+init_bg+8, &__k, &__l);

      if(__n >0 && iter == 0xFFFFFFFF){
        __k0 = index->fastmap->item[iter];
        __l0 = index->fastmap->item[iter+1];
        __l0 -= get_12mer_correct(index->fastmap_correct, __l0-1);
        fprintf(stderr, "12 mer search  fail!!\n");
        fprintf(stderr, "seed_st = %u, k = %u, l = %u, __k =%u, __l = %u, [%u, %u],\n", seed_start,__k0, __l0, __k, __l, init_bg, init_ed);
        exit(1);
      } 
    

#endif 


      if(iter ==0xFFFFFFFF) continue; 
      k = index->fastmap->item[iter];
      l = index->fastmap->item[iter+1];
      l -= get_12mer_correct(index->fastmap_correct, l-1)+1;

     
      if(k > l) continue;
      n = bwt_match_exact_alt(index->bwt, 8, seq+init_bg, &k, &l);
      if(n == 0) continue;
      if(n > IS_SMLSIZ) {
        aln_mem_alt(index, l_seq, seq, &init_bg, &init_ed, &k, &l, IS_SMLSIZ);       
      }
#ifdef DEBUG
      fprintf(stderr, "seed_st = %u, k = %u, l = %u, [%u, %u]\n", seed_start, k, l, init_bg, init_ed);
#endif 
      //if(k + opt->max_locate > l) {
      for(i = k; i <= __MIN(l, k+opt->max_locate); ++i) {
        p = bwt_sa(index->bwt, i);
#ifdef DEBUG
        
        print_real_coor(index->bns, p, l_seq);
#endif
        kv_push(uint32_t, *candi, p>init_bg?p-init_bg:0); // append bgn
        kv_push(uint32_t, *candi, p+l_seq-init_bg); // append end
      }
      //}         
    }    
    ks_introsort(uint32_t, candi->n, candi->a);
}
void alnse_chain(candi_t *candi, int max_gap, int max_len, kvec_t(chain_t) *chain)
{
  uint32_t i, j, p0, p1;
  uint32_t *a = candi->a;
  chain_t c = {0, 0, 0};
  if(candi->n == 0) return;
  kv_init(*chain);
  kv_push(chain_t, *chain, c);
  p0 = a[0];
  chain->a[0].bg = p0; 
  chain->a[0].ed = p0;
  chain->a[0].n = 1; 
  for(i = 1, j = 0; i < candi->n; ++i){
    p1 = a[i];
    if(p0 + max_gap >= p1 && p1 - chain->a[j].bg < max_len) {
      chain->a[j].ed = p1;
      ++chain->a[j].n; 
    } else {
      kv_push(chain_t, *chain, c);
      j = chain->n - 1;
      chain->a[j].bg = p1; 
      chain->a[j].ed = p1;
      chain->a[j].n = 1; 
    }
    p0 = p1;
  }
  ks_introsort(chain, chain->n, chain->a);
}
kswr_t alnse_extend(uint8_t *pac, uint32_t ref_bg, uint32_t ref_ed, int l_seq, const uint8_t *seq, kswq_t **q, int8_t mat[25], int minsc)
{
  int i;
  uint8_t ref[1024];
  int gapo = 5, gape = 2, xtra = KSW_XSTART;
  if(ref_ed - ref_bg >= 1024) {
    fprintf(stderr, "Ref seq len (%u) too long!", ref_ed - ref_bg); 
  }
  for(i = 0; i < ref_ed - ref_bg; ++i) {
    ref[i] = __get_pac(pac, ref_bg+i); 
  }
  //kswq_t *q[2] = {0, 0};
  xtra |= KSW_XBYTE|minsc;
  kswr_t r = ksw_align(l_seq, seq, ref_ed -ref_bg, ref,5, mat, gapo, gape, xtra, q);
  //free(q[0]);
  return r; 
}
int alnse_global(uint8_t *pac, uint32_t ref_bg, uint32_t ref_ed, int l_seq, const uint8_t *seq, int8_t mat[25], int minsc)
{
  int i;
  uint8_t ref[1024];
  int gapo = 5, gape = 2, xtra = KSW_XSTART;
  for(i = 0; i < ref_ed - ref_bg; ++i) {
    ref[i] = __get_pac(pac, ref_bg+i); 
  }
  int n_cigar;
  uint32_t *cigar = NULL;
  ksw_global(l_seq, seq, ref_ed-ref_bg, ref, 5, mat, gapo, gape, 100, &n_cigar, &cigar);
  fprintf(stderr, "sc = %d, cigar:", minsc);
  for(i =0; i < n_cigar; ++i) {
    fprintf(stderr, "%u%c", cigar[i]>>4, "MID"[cigar[i]&0xF]); 
  } 
  fprintf(stderr, "\n");
}
int gen_cigar(uint8_t *pac, uint32_t ref_bg, uint32_t ref_ed, int l_seq, const uint8_t *seq, int8_t mat[25], int minsc, int *n_cigar, uint32_t **cigar)
{
  int i;
  uint8_t ref[1024];
  int gapo = 5, gape = 2, xtra = KSW_XSTART;
  for(i = 0; i < ref_ed - ref_bg; ++i) {
    ref[i] = __get_pac(pac, ref_bg+i); 
  }

  int w = ((l_seq>>1)*mat[0]-gapo)/gape;
  ksw_global(l_seq, seq, ref_ed-ref_bg, ref, 5, mat, gapo, gape, w, n_cigar, cigar);
  /*  
  fprintf(stderr, "sc = %d, cigar:", minsc);
  for(i =0; i < *n_cigar; ++i) {
    fprintf(stderr, "%u%c", (*cigar)[i]>>4, "MID"[(*cigar)[i]&0xF]); 
  } 
  fprintf(stderr, "\n");
  */
}


void alnse_core_single(idx_t *index, query_t *query, aln_opt_t* aln_opt)
{
  int i, j;   
  uint32_t max_locate = aln_opt->max_locate;
  int max_diff = aln_opt->max_diff; 
 


  int l_seq = query->l_seq;
  uint8_t *seq= query->seq;
  uint8_t *rseq = query->rseq;
  //fprintf(stderr, "[%s]: seeding\n", query->name);
  //seeding 
  candi_t candi0, candi1;
  kv_init(candi0);  kv_init(candi1);   

  alnse_seed_overlap(index, l_seq, seq, aln_opt, &candi0);
  alnse_seed_overlap(index, l_seq, rseq, aln_opt, &candi1);


  //chaining
  kvec_t(chain_t) chain0, chain1;
  kv_init(chain0);  kv_init(chain1);   
  alnse_chain(&candi0, aln_opt->max_chain_gap, 1000, &chain0);
  alnse_chain(&candi1, aln_opt->max_chain_gap, 1000, &chain1);
  //extending
#ifdef DEBUG
  for(j = 0; j < candi0.n; ++j) {
    fprintf(stderr, "%u\t",candi0.a[j]);
  } 
  fprintf(stderr, "\n[%s]: backward\n", query->name);
  for(j = 0; j < candi1.n; ++j) {
    fprintf(stderr, "%u\t",candi1.a[j]);
  } 
  fprintf(stderr, "\n");

  fprintf(stderr, "[%s]: chaining\n", query->name);
  for(j = 0; j < chain0.n; ++j) {
    fprintf(stderr, "[bg, ed, n, l] = [%u, %u, %u, %u]\n", chain0.a[j].bg, chain0.a[j].ed, chain0.a[j].n, chain0.a[j].ed - chain0.a[j].bg);
    print_real_coor(index->bns, chain0.a[j].bg, query->l_seq);
    print_real_coor(index->bns, __MIN(chain0.a[j].ed, index->bwt->seq_len-1), query->l_seq);
  } 
  fprintf(stderr, "\n[%s]: backward\n", query->name);
  for(j = 0; j < chain1.n; ++j) {
    fprintf(stderr, "[bg, ed, n, l] = [%u, %u, %u, %u]\n", chain1.a[j].bg, chain1.a[j].ed, chain1.a[j].n, chain1.a[j].ed - chain1.a[j].bg);
    print_real_coor(index->bns, chain1.a[j].bg, query->l_seq);
    print_real_coor(index->bns, __MIN(chain1.a[j].ed, index->bwt->seq_len-1), query->l_seq);

  } 
  
#endif   
  // 
  uint32_t ref_bg, ref_ed;
  int n_extend = aln_opt->max_extend < chain0.n? aln_opt->max_extend:chain0.n;
  kvec_t(aln_t) aln;
  kv_init(aln);
  int thres = query->l_seq/3;
  int maxsc = thres, maxi = -1; 
  kswq_t *q[2] = {0, 0};
  for(j =0; j < n_extend; ++j) {
    ref_bg = chain0.a[j].bg; 
    ref_ed = __MIN(chain0.a[j].ed, index->bwt->seq_len-1); 
    kswr_t r = alnse_extend(index->pac, ref_bg, ref_ed, query->l_seq, seq, &q[0], index->mat, maxsc);
#ifdef DEBUG
    fprintf(stderr, "[bg, ed, n, l, sc] = [%u, %u, %u, %u, %u, %u, %u, %u, %u]\n", chain0.a[j].bg, chain0.a[j].ed, chain0.a[j].n, chain0.a[j].ed - chain0.a[j].bg, r.score, r.tb, r.te, r.qb, r.qe);
    print_real_coor(index->bns, ref_bg, query->l_seq);
    print_real_coor(index->bns, ref_ed, query->l_seq);
 
#endif
    if(r.score > thres) {
      aln_t x;
      x.strand = 0;
      x.sc = r.score;
      x.tb = r.tb+ref_bg;
      x.te = r.te+ref_bg+1;
      x.qb = r.qb;
      x.qe = r.qe+1;
      kv_push(aln_t, aln, x); // append bgn
      if(r.score > maxsc){ 
        maxsc = r.score;
        maxi = aln.n-1;    
      }
 
    }
  }
  n_extend = aln_opt->max_extend < chain1.n? aln_opt->max_extend:chain1.n;
  for(j =0; j < n_extend; ++j) {
    ref_bg = chain1.a[j].bg; 
    ref_ed = __MIN(chain1.a[j].ed, index->bwt->seq_len-1);  
    kswr_t r = alnse_extend(index->pac, ref_bg, ref_ed, query->l_seq, rseq, &q[1], index->mat, maxsc);
#ifdef DEBUG
    fprintf(stderr, "[bg, ed, n, l, sc] = [%u, %u, %u, %u, %u, %u, %u, %u, %u]\n", chain1.a[j].bg, chain1.a[j].ed, chain1.a[j].n, chain1.a[j].ed - chain1.a[j].bg, r.score, r.tb, r.te, r.qb, r.qe);
    print_real_coor(index->bns, ref_bg, query->l_seq);
    print_real_coor(index->bns, ref_ed, query->l_seq);
 
#endif
    if(r.score > thres) {
      aln_t x;
      x.strand = 1;
      x.sc = r.score;
      x.tb = r.tb+ref_bg;
      x.te = r.te+ref_bg+1;
      x.qb = r.qb;
      x.qe = r.qe+1;
      kv_push(aln_t, aln, x); // append bgn
      if(r.score > maxsc){ 
        maxsc = r.score;
        maxi = aln.n-1;    
      }
    }
  }
  free(q[0]);
  free(q[1]);
  if(maxi != -1) {
    query->pos = query->ref_start = aln.a[maxi].tb;
    query->b0 = maxsc;
    query->strand = aln.a[maxi].strand;
    query->ref_start = aln.a[maxi].tb;
    query->ref_end = aln.a[maxi].te;
    query->seq_start = aln.a[maxi].qb;
    query->seq_end = aln.a[maxi].qe;

    //alnse_global(index->pac, aln.a[maxi].tb, aln.a[maxi].te, l_seq, aln.a[maxi].strand==0?seq:rseq, index->mat, maxsc);
  } else {
    query->pos = 0xFFFFFFFF;
    
  }
 

  //destroy
  kv_destroy(candi0);
  kv_destroy(candi1);
  kv_destroy(chain0);
  kv_destroy(chain1);
  kv_destroy(aln);
  return;    
}
void alnse_core1(idx_t *index, int n, query_t *multi_query, aln_opt_t* aln_opt)
{
  int i, j;   
  for(i = 0; i < n; ++i) {
    alnse_core_single(index, multi_query+i, aln_opt); 
  }
}
/*  
void alnse_core1(idx_t *index, int n, query_t *multi_query, aln_opt_t* aln_opt)
{
  int i, j;   
  uint32_t max_locate = aln_opt->max_locate;
  int max_diff = aln_opt->max_diff; 
 
  for(i = 0; i < n; ++i) {
    query_t *query = multi_query + i; 
    int l_seq = query->l_seq;
    uint8_t *seq= query->seq;
    uint8_t *rseq = query->rseq;
    //fprintf(stderr, "[%s]: seeding\n", query->name);
    //seeding 
    candi_t candi0, candi1;
    kv_init(candi0);  kv_init(candi1);   

    alnse_seed_overlap(index, l_seq, seq, aln_opt, &candi0);
    alnse_seed_overlap(index, l_seq, rseq, aln_opt, &candi1);


    //chaining
    kvec_t(chain_t) chain0, chain1;
    kv_init(chain0);  kv_init(chain1);   
    alnse_chain(&candi0, aln_opt->max_chain_gap, 512, &chain0);
    alnse_chain(&candi1, aln_opt->max_chain_gap, 512, &chain1);
    //extending
#ifdef DEBUG
    for(j = 0; j < candi0.n; ++j) {
      fprintf(stderr, "%u\t",candi0.a[j]);
    } 
    fprintf(stderr, "\n[%s]: backward\n", query->name);
    for(j = 0; j < candi1.n; ++j) {
      fprintf(stderr, "%u\t",candi1.a[j]);
    } 
    fprintf(stderr, "\n");

    fprintf(stderr, "[%s]: chaining\n", query->name);
    for(j = 0; j < chain0.n; ++j) {
      fprintf(stderr, "[bg, ed, n, l] = [%u, %u, %u, %u]\n", chain0.a[j].bg, chain0.a[j].ed, chain0.a[j].n, chain0.a[j].ed - chain0.a[j].bg);
    } 
    fprintf(stderr, "\n[%s]: backward\n", query->name);
    for(j = 0; j < chain1.n; ++j) {
      fprintf(stderr, "[bg, ed, n, l] = [%u, %u, %u, %u]\n", chain1.a[j].bg, chain1.a[j].ed, chain1.a[j].n, chain1.a[j].ed - chain1.a[j].bg);
    } 
    
#endif   
    // 
    uint32_t ref_bg, ref_ed;
    int n_extend = aln_opt->max_extend < chain0.n? aln_opt->max_extend:chain0.n;
    kvec_t(aln_t) aln;
    kv_init(aln);
    int thres = query->l_seq/3;
    int maxsc = thres, maxi = -1; 
    for(j =0; j < n_extend; ++j) {
      ref_bg = chain0.a[j].bg; 
      ref_ed = chain0.a[j].ed; 
      kswr_t r = alnse_extend(index->pac, ref_bg, ref_ed, query->l_seq, seq, index->mat, maxsc);
#ifdef DEBUG
      fprintf(stderr, "[bg, ed, n, l, sc] = [%u, %u, %u, %u, %u, %u, %u, %u, %u]\n", chain0.a[j].bg, chain0.a[j].ed, chain0.a[j].n, chain0.a[j].ed - chain0.a[j].bg, r.score, r.tb, r.te, r.qb, r.qe);
#endif
      if(r.score > thres) {
        aln_t x;
        x.strand = 0;
        x.sc = r.score;
        x.tb = r.tb+ref_bg;
        x.te = r.te+ref_bg+1;
        x.qb = r.qb;
        x.qe = r.qe+1;
        kv_push(aln_t, aln, x); // append bgn
        if(r.score > maxsc){ 
          maxsc = r.score;
          maxi = aln.n-1;    
        }
   
      }
    }
    n_extend = aln_opt->max_extend < chain1.n? aln_opt->max_extend:chain1.n;
    for(j =0; j < n_extend; ++j) {
      ref_bg = chain1.a[j].bg; 
      ref_ed = chain1.a[j].ed;  
      kswr_t r = alnse_extend(index->pac, ref_bg, ref_ed, query->l_seq, rseq, index->mat, maxsc);
#ifdef DEBUG
      fprintf(stderr, "[bg, ed, n, l, sc] = [%u, %u, %u, %u, %u, %u, %u, %u, %u]\n", chain1.a[j].bg, chain1.a[j].ed, chain1.a[j].n, chain1.a[j].ed - chain1.a[j].bg, r.score, r.tb, r.te, r.qb, r.qe);
#endif
      if(r.score > thres) {
        aln_t x;
        x.strand = 1;
        x.sc = r.score;
        x.tb = r.tb+ref_bg;
        x.te = r.te+ref_bg+1;
        x.qb = r.qb;
        x.qe = r.qe+1;
        kv_push(aln_t, aln, x); // append bgn
        if(r.score > maxsc){ 
          maxsc = r.score;
          maxi = aln.n-1;    
        }
      }
    }
    if(maxi != -1) {
      query->pos = query->ref_start = aln.a[maxi].tb;
      query->b0 = maxsc;
      query->strand = aln.a[maxi].strand;
      query->ref_start = aln.a[maxi].tb;
      query->ref_end = aln.a[maxi].te;
      query->seq_start = aln.a[maxi].qb;
      query->seq_end = aln.a[maxi].qe;

      //alnse_global(index->pac, aln.a[maxi].tb, aln.a[maxi].te, l_seq, aln.a[maxi].strand==0?seq:rseq, index->mat, maxsc);
    } else {
      query->pos = 0xFFFFFFFF;
      
    }


    //destroy
    
    
    
    kv_destroy(candi0);
    kv_destroy(candi1);
    kv_destroy(chain0);
    kv_destroy(chain1);
    kv_destroy(aln);
  }
  return;    
}
*/

#ifdef HAVE_THREAD

static void *alnse_worker(void *data)
{
    thread_aux_t *d = (thread_aux_t*)data;
    //alnse_core1(d->tid, d->index, d->n_query, d->query, d->aln_opt, d->aux); 
    alnse_core_thread(d->tid, d->index, d->n_query, d->query, d->aln_opt); 
    return 0;
}
#endif
#define MAX_N_PERSEQ 200
void alnse_core_thread(int tid, idx_t *index, int n_query, query_t *multi_query, aln_opt_t *aln_opt)
{
  int i;
  while(1){
    pthread_rwlock_wrlock(&rwlock);
    i = g_n_seqs++; 
    pthread_rwlock_unlock(&rwlock);
    if(i >= n_query) break;
    query_t *query = multi_query+i;
    if(query->n_ambiguous > MAX_N_PERSEQ) continue;
    alnse_core_single(index, query, aln_opt); 
    aln_samse(index, query, aln_opt);
 }  
}






int alnse_core(const opt_t *opt)
{
    fprintf(stderr, "[alnse_core]:  Start single end alignment!\n");
    double t_real = realtime();     
    aln_opt_t *aln_opt = aln_opt_init(opt);
    fprintf(stderr, "[%s]:  Reload index...\n", __func__); 
    idx_t *idx = fbwt_fmidx_restore(opt->fn_index);
    fbwt_hier_idx_restore(idx, opt->fn_index);
    fprintf(stderr, "%lf sec escaped.\n", realtime()-t_real);
    t_real = realtime();
  

#ifdef HAVE_THREAD
    pthread_attr_t attr;
    pthread_t *tid; 
    thread_aux_t *thread_data;
    pthread_rwlock_init(&rwlock, NULL); 
    if(opt->n_threads<=1){
    }else{
      pthread_attr_init(&attr);
      pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
      tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
      thread_data = calloc(opt->n_threads, sizeof(thread_aux_t)); 
      int j;
      for(j=0; j < opt->n_threads;++j){
          thread_data[j].tid = j; thread_data[j].index = idx;
          thread_data[j].aln_opt= aln_opt;
      } 
    }       
#else
#endif

    queryio_t *qs = query_open(opt->fn_read1);
    query_t *multiSeqs = calloc(N_SEQS, sizeof(query_t));
    if(multiSeqs == NULL){
        fprintf(stderr, "[%s]:  alloc mem fail!\n", __func__);
        exit(1);
    }
    fprintf(stderr, "[%s]:  Begin to align short reads...\n", __func__); 
    aln_samhead(opt, idx->bns);
    int n, i, n_tot=0;
    while( (n = query_read_multiSeqs(qs, N_SEQS, multiSeqs)) >0){
        n_tot += n;
        //if(n_tot < 3000000) continue;
#ifdef HAVE_THREAD
        g_n_seqs = 0;
        if(opt->n_threads <=1){
            alnse_core1(idx, n, multiSeqs, aln_opt); 
            for(i = 0; i< n; ++i){
                query_t *query = multiSeqs+i;
                aln_samse(idx, query, opt);
            }
        }else{
            int j;
            for (j = 0; j < opt->n_threads; ++j) {
                thread_data[j].query = multiSeqs;
                thread_data[j].n_query = n; 
                pthread_create(&tid[j], &attr, alnse_worker, thread_data + j);
            }
            for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0); 
        }
#else            
        alnse_core1(idx, n, multiSeqs, aln_opt); 
        for(i = 0; i< n; ++i){
            query_t *query = multiSeqs+i;
            aln_samse(idx, query, opt);
        } 

#endif       
        for(i = 0; i< n; ++i){
            query_t *query = multiSeqs+i;
            //aln_samse(idx, query, opt);
            //printf("%s\n", query->sam->s);
            puts(query->sam->s);
            query_destroy(query);
        } 
        fflush(stdout); 
        memset(multiSeqs, 0, N_SEQS*sizeof(query_t));

        //if(n_tot%10000==0 && n_tot!=0) fprintf(stderr, "%d reads have been aligned!\n", n_tot);
        fprintf(stderr, "%d reads have been aligned!\n", n_tot);
    }
    //fprintf(stderr, "[%s]: total %.2f sec escaped\n", __func__, (float)(clock()-t)/CLOCKS_PER_SEC);
    fprintf(stderr, "[%s]: total %lf sec escaped\n", __func__, realtime()-t_real);

    
#ifdef HAVE_THREAD
    pthread_rwlock_destroy(&rwlock);
    if(opt->n_threads<=1){
    }else{
        free(tid);
        free(thread_data);
    }
#else
#endif    
    free(multiSeqs);
    query_close(qs);
    aln_opt_destroy(aln_opt);



    return EXIT_SUCCESS;



}
