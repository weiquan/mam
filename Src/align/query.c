/*
 * =====================================================================================
 *
 *       Filename:  query.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/28/2012 10:11:20 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Quan, Wei (mn), wquanhit@gmail.com
 *        Company:  BIC, HIT, China
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ksw.h"
//#include "kstring.h"

//#include "bwt.h"
//#include "editdistance.h"
//#include "aln.h"
#include "query.h"
#include "SeedExtStruct.h"
#include "debug.h"
//#define UINT_MAX 0xFFFFFFFF
#define SCORE_MATCH 1
#define SCORE_MISMATCH -2
#define SCORE_GAP -5
#define CIGAR_MATCH 0
#define CIGAR_INS 1
#define CIGAR_DEL 2


unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

//extern unsigned char nst_nt4_table[256];

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  query_seq_reverse
 *  Description:  from bwa seqio.c
 *                if comp == 1 return reverse complement
 *                if comp == 0 return mirror seq 
: * =====================================================================================
 */
void query_seq_reverse(int len, uint8_t *seq, int is_comp)
{
	
    int i;
	if (is_comp) {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			if (tmp < 4) tmp = 3 - tmp;
			seq[len-1-i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
			seq[i] = tmp;
		}
		if (len&1) seq[i] = (seq[i] >= 4)? seq[i] : 3 - seq[i];
	} else {
		for (i = 0; i < len>>1; ++i) {
			char tmp = seq[len-1-i];
			seq[len-1-i] = seq[i]; seq[i] = tmp;
		}
	}
}
#define FLAG_UNMAPPED 0x0004
#define FLAG_REV 0x0010

#define __get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
/*  
void query_add_best(query_t *query, idx_t *idx, int mat[25], uint32_t pos, struct SubBuf *sub)
{
    int i;
    //query->pos = pos;
    //uint8_t *target = (uint8_t *)malloc(query->l_seq);
    uint8_t target[LEN_READ+50];
    uint32_t l_pos = pos > MAX_CLIP?pos-MAX_CLIP:0;
    uint32_t r_pos = pos+query->l_seq+MAX_CLIP;
    for(i =0; i < r_pos - l_pos; ++i) target[i] = __get_pac(idx->pac, l_pos+i);
    
    int xtra = KSW_XSTART;
    kswq_t *q[2] = {0, 0};
	kswr_t r = ksw_align(query->l_seq, query->is_rev?query->rseq:query->seq, r_pos-l_pos, target, 5, mat, 6, 1, xtra, &q[0], sub->kswq_R);
fprintf(stderr, "pos = %u, r.score = %d, query->b0 = %d, query->pos = %u\n", pos, r.score, query->b0, query->pos);
    if(r.score > query->b0) {
        query->strand = query->is_rev;
        query->b0 = r.score;
        query->pos = pos;

        query->tlen = r.te+1-r.tb;
        query->ref_start = r.tb;
        query->seq_start = r.qb;
        query->seq_end = r.qe+1;
        query->n_cigar = 0;
        //if(query->cigar) free(query->cigar);

        //int score = ksw_global(r.qe-r.qb, (query->is_rev?query->rseq:query->seq)+r.qb, r.te-r.tb, target+r.tb, 5, mat, 6, 1, 20, &query->n_cigar, &query->cigar, &sub->kswgb);
        //fprintf(stderr, "score = %d, n_cigar = %d\n", score, query->n_cigar);
        //for(i = 0; i < query->n_cigar; ++i) fprintf(stderr, "%u%c", query->cigar[i]>>4, "MID"[query->cigar[i]&0xF]); 

    }


    //free(target);
    return;
}
*/
void query_aln2sam1(query_t *query,  idx_t *idx, struct SubBuf *sub,  FILE *fp)
{
    bntseq_t *bns = idx->bns;

    int i, j;
    unsigned short flag = 0;
    int nm = 255; 
//fflush(fp);
fprintf(fp, "%s\t", query->qseq->name); 
    //flag Rname pos mapq cigar Rnext Pnext tlen
    if(query->strand == 1) flag |= FLAG_REV;
    if(query->pos ==(uint32_t)-1) {
        flag |= FLAG_UNMAPPED;
fprintf(fp, "%u\t*\t0\t255\t*\t*\t0\t0\t",flag); 
    } else{
        uint8_t target[LEN_READ+50];

        uint32_t pos = query->pos;
        uint32_t l_pos = pos > MAX_CLIP?pos-MAX_CLIP:0;
        uint32_t r_pos = pos+query->l_seq+MAX_CLIP;
        for(i =0; i < r_pos - l_pos; ++i) target[i] = __get_pac(idx->pac, l_pos+i);
        uint8_t *seq = query->strand==1?query->rseq:query->seq;
        kswq_t **qry, **kswq; 
        if(query->strand) {
            qry = &sub->qry_r[0];
            kswq = sub->kswq_r; 
        } else{
            qry = &sub->qry_f[0];
            kswq = sub->kswq_f; 
        }
  
        kswr_t r = ksw_align2(query->seq_end-1, seq, query->ref_end-1, target, 5, sub->mat, 6, 1, 0, qry, kswq, query->b0);
        query->ref_start = r.tb;
        query->seq_start = r.qb;

        int score = ksw_global(query->seq_end-query->seq_start, seq+query->seq_start, query->ref_end-query->ref_start, target+query->ref_start, 5, sub->mat, 6, 1, 50, &query->n_cigar, &query->cigar, &sub->kswgb);
        
        uint8_t *x0 = target+query->ref_start, *x1 = seq+query->seq_start;

        nm = 0;
        for(i = 0; i < query->n_cigar; ++i){
            int op = query->cigar[i] & 0xF;
            int len = query->cigar[i] >>4;
            if(op == 0) {
                for(j =0; j < len; ++j, ++x0, ++x1) {
                    if(*x0 != *x1) ++nm; 
                }
            } else if(op == 1) {
                nm += len;
                x1 += len; 
            } else{
                nm += len;
                x0 += len;
            }
        }

        
        if(query->ref_start <= MAX_CLIP) {
            query->pos -= MAX_CLIP-query->ref_start; 
            query->ref_start = 0;
        } else{
            query->pos += query->ref_start-MAX_CLIP;
            query->ref_start -= MAX_CLIP; 
        }
        
        uint32_t rid;
bns_coor_pac2real(bns, query->pos, query->l_seq, &rid); 

fprintf(fp, "%u\t%s\t%u\t255\t",flag, bns->anns[rid].name, query->pos-bns->anns[rid].offset+1);
if(query->seq_start != 0) fprintf(fp, "%uS", query->seq_start);
        for(i = 0; i < query->n_cigar; ++i){ 
            if((i == 0|| i == query->n_cigar-1) && (query->cigar[i] &3) ==2 ) {
                query->b0 += 7;
                continue;
            }
fprintf(fp, "%u%c", query->cigar[i]>>4, "MID"[query->cigar[i]&3]);
        }
if(query->seq_end != query->l_seq) fprintf(fp, "%uS", query->l_seq-query->seq_end);
fprintf(fp, "\t*\t0\t0\t");
        //if(query->cigar != NULL )free(query->cigar);


    }
    if(query->strand == 1) {
        for(i = 0; i < query->l_seq; ++i){ 
fprintf(fp, "%c", "ACGTN"[query->rseq[i]]);
        }
fprintf(fp, "\t");
        for(i = query->l_seq-1; i >=0; --i) {
fprintf(fp, "%c", query->qual[i]);
        }
    } else{ 
        for(i = 0; i < query->l_seq; ++i) {
fprintf(fp, "%c", "ACGTN"[query->seq[i]]);
        }
fprintf(fp, "\t");
        for(i = 0; i < query->l_seq; ++i) {
fprintf(fp, "%c", query->qual[i]);
        }
    }
    
fprintf(fp, "\tAS:i:%d\tNM:i:%d", query->b0, nm);
    //+++++++++++++++++++++++++++++++++++++
//fprintf(stderr, "%u, %s\n", __LINE__, __func__ );

    uint32_t (*out_buf)[2] = sub->aln_out->out_buf; 
    int (*found)[4] = sub->aln_out->found; 
    int out_len = sub->aln_out->len;
    int out_num = sub->aln_out->num; 
 
    int fnd_i, fnd_num, f_score, f_bg, f_ed, f_num, p, cur_row, cur_num;
    uint32_t (*pos_out)[2] = sub->aln_out->out_buf+out_len;

    fnd_i = 0;
    //int delta = sub->delta;
    int delta = 1;
    for(i = 0; i < delta; ++i){
        if(found[i][1] > 0) {
            found[fnd_i][0] = found[i][0];  
            found[fnd_i][1] = found[i][1];  
            found[fnd_i][2] = found[i][2];  
            found[fnd_i][3] = found[i][3];  
            fnd_i++;
        }
    }
    fnd_num = fnd_i;
    int pos_i = 0;
    for(fnd_i = 0; fnd_i < fnd_num; ++fnd_i){
        //f_score = found[fnd_i][0];
        //f_num = found[fnd_i][1];
        //f_ed = found[fnd_i][3]; 
        f_bg = found[fnd_i][2];
        p = f_bg;
        while(1){
            cur_num = out_buf[p][0];
            for(j = 0; j < cur_num; ++j) {
                pos_out[pos_i][0] = out_buf[p+j+1][0];  
                pos_out[pos_i][1] = out_buf[p+j+1][1];  
                pos_i++;   
            }
            if(out_buf[p][1] > 0) {
                p = out_buf[p][1]; 
            } else {
                break;
            }
        }
    }
    //------------------------------------- 
    int n_aln = 0;    
    for(i = 0; i < pos_i; ++i) { 
        uint32_t pos = pos_out[i][0]; 
        int sc = pos_out[i][1]%1024; 
        if( pos != query->pos) {
            if(n_aln == 0) fprintf(fp, "\tXA:Z:");
            int rid;
            bns_coor_pac2real(bns, pos, query->l_seq, &rid); 
fprintf(fp, "%s,%u,%d;",bns->anns[rid].name, pos-bns->anns[rid].offset+1,sc);
            //fprintf(stderr,"pos = %u, ref = %s,pos = %u\n",aln_pos[i], bns->anns[rid].name, aln_pos[i]-bns->anns[rid].offset+1);
            //fprintf(fp, "%u;", aln_pos[i]+1);
            n_aln++;
        }
    }
fprintf(fp, "\n");
    return;
}
int print_pos(idx_t *fm_idx, int n_pos, uint32_t *pos_buf, int len)
{
return 0;
    int i, rid;
fprintf(stderr, "%u, pos_i = %d\n", __LINE__, n_pos);
    for(i = 0; i < n_pos; ++i) {
        uint32_t pos = pos_buf[i];
        bns_coor_pac2real(fm_idx->bns, pos, len, &rid); 
        fprintf(stderr, "pos_i = %d, pos = %u, chr = %s, pos = %u\n",i, pos, fm_idx->bns->anns[rid].name, pos-fm_idx->bns->anns[rid].offset+1);
    }
}
int query_add_md(uint8_t *target, int tb, int te, uint8_t *read_seq, int qb, int qe, int n_cigar, int *cigar) 
{
    int i, j;    
    int nm = 0, sc = 0, match = 1, mis = -4, gap_o = 6, gap_e = 1;
    uint8_t *x0 = target+tb, *x1 = read_seq+qb;
    for(i = 0; i < n_cigar; ++i){
        int op = cigar[i] & 0xF;
        int len = cigar[i] >>4;
        fprintf(stderr, "%d %c\n", len, "MID"[op]);
        if(op == 0) {
            for(j =0; j < len; ++j, ++x0, ++x1) {
                
                if(*x0 != *x1) {
fprintf(stderr, "%3u\t%3u\t%d\t%d\tS\t%d\n", (x0-target), (x1-read_seq), *x0, *x1,sc);              
                    ++nm;
                    sc += mis;
                } else {
fprintf(stderr, "%3u\t%3u\t%d\t%d\t=\t%d\n", (x0-target), (x1-read_seq), *x0, *x1, sc);              
                    sc += match;
                } 
            }
        } else if(op == 1) {
            for(j = 0; j < len; ++j, ++x1) {
              
fprintf(stderr, "%3u\t%3u\t%c\t%d\tI\t%d\n", (x0-target), (x1-read_seq), '-', *x1, sc);              
                if(j == 0) sc -= gap_o;
                sc -= gap_e;
            }
            nm += len;
            //x1 += len; 
        } else{
    
            for(j = 0; j < len; ++j, ++x0) {
              
fprintf(stderr, "%3u\t%3u\t%d\t%c\tD\t%d\n", (x0-target), (x1-read_seq), *x0, '-', sc);              
                if(j == 0) sc -= gap_o;
                sc -= gap_e;
            }

            nm += len;
            //x0 += len;
        }
    }


}
int query_aln2sam3(query_t *query,  idx_t *idx, struct SubBuf *sub,  FILE *fp_sam)
{
    int len_str = 0;
    
    bntseq_t *bns = idx->bns;

    int i, j;
    unsigned short flag = 0;
    int nm = 255; 
//fflush(fp);
    fprintf(fp_sam, "%s\t", query->qseq->name); 
    //flag Rname pos mapq cigar Rnext Pnext tlen
    if(query->strand == 1) flag |= FLAG_REV;
    if(query->pos ==(uint32_t)-1) {
        flag |= FLAG_UNMAPPED;
        fprintf(fp_sam, "%u\t*\t0\t255\t*\t*\t0\t0\t",flag); 
    } else{
        uint8_t target[LEN_READ+50];

        uint32_t pos = query->pos;
        uint32_t l_pos = pos > MAX_CLIP?pos-MAX_CLIP:0;
        uint32_t r_pos = pos+query->l_seq+MAX_CLIP;
        for(i =0; i < r_pos - l_pos; ++i) target[i] = __get_pac(idx->pac, l_pos+i);
        uint8_t *seq = query->strand==1?query->rseq:query->seq;
        kswq_t **qry, **kswq; 
        if(query->strand) {
            qry = &sub->qry_r[0];
            kswq = sub->kswq_r; 
        } else{
            qry = &sub->qry_f[0];
            kswq = sub->kswq_f; 
        }
        if(query->seq_start == -1 && query->ref_start == -1) {
            kswr_t r = ksw_align2(query->seq_end-1, seq, query->ref_end-1, target, 5, sub->mat, 6, 1, 0, qry, kswq, query->b0);
            query->ref_start = r.tb;
            query->seq_start = r.qb;
        }

        int score = ksw_global(query->seq_end-query->seq_start, seq+query->seq_start, query->ref_end-query->ref_start, target+query->ref_start, 5, sub->mat, 6, 1, 50, &query->n_cigar, &query->cigar, &sub->kswgb);
        
        uint8_t *x0 = target+query->ref_start, *x1 = seq+query->seq_start;

        nm = 0;
        for(i = 0; i < query->n_cigar; ++i){
            int op = query->cigar[i] & 0xF;
            int len = query->cigar[i] >>4;
            if(op == 0) {
                for(j =0; j < len; ++j, ++x0, ++x1) {
                    if(*x0 != *x1) ++nm; 
                }
            } else if(op == 1) {
                nm += len;
                x1 += len; 
            } else{
                nm += len;
                x0 += len;
            }
        }

        
        if(query->ref_start <= MAX_CLIP) {
            query->pos -= MAX_CLIP-query->ref_start; 
            query->ref_start = 0;
        } else{
            query->pos += query->ref_start-MAX_CLIP;
            query->ref_start -= MAX_CLIP; 
        }
        
        uint32_t rid;
        bns_coor_pac2real(bns, query->pos, query->l_seq, &rid); 
        fprintf(fp_sam, "%u\t%s\t%u\t255\t",flag, bns->anns[rid].name, query->pos-bns->anns[rid].offset+1);
        if(query->seq_start != 0) fprintf(fp_sam, "%uS", query->seq_start);
        for(i = 0; i < query->n_cigar; ++i){ 
            if((i == 0|| i == query->n_cigar-1) && (query->cigar[i] &3) ==2 ) {
                query->b0 += 7;
                continue;
            }
            fprintf(fp_sam, "%u%c", query->cigar[i]>>4, "MID"[query->cigar[i]&3]);
        }
        if(query->seq_end != query->l_seq) fprintf(fp_sam, "%uS", query->l_seq-query->seq_end);
        fprintf(fp_sam, "\t*\t0\t0\t");
        //if(query->cigar != NULL )free(query->cigar);


    }
    if(query->strand == 1) {
        for(i = 0; i < query->l_seq; ++i){ fprintf(fp_sam, "%c", "ACGTN"[query->rseq[i]]);}
        fprintf(fp_sam, "\t");
        for(i = query->l_seq-1; i >=0; --i) { fprintf(fp_sam, "%c", query->qual[i]);}
    } else{ 
        for(i = 0; i < query->l_seq; ++i) { fprintf(fp_sam, "%c", "ACGTN"[query->seq[i]]); }
        fprintf(fp_sam, "\t");
        for(i = 0; i < query->l_seq; ++i) { fprintf(fp_sam, "%c", query->qual[i]); }
    }
    fprintf(fp_sam, "\tAS:i:%d\tNM:i:%d", query->b0, nm);
    //+++++++++++++++++++++++++++++++++++++
//fprintf(stderr, "%u, %s\n", __LINE__, __func__ );

    uint32_t (*out_buf)[2] = sub->aln_out->out_buf; 
    int (*found)[4] = sub->aln_out->found; 
    int out_len = sub->aln_out->len;
    int out_num = sub->aln_out->num; 
 
    int fnd_i, fnd_num, f_score, f_bg, f_ed, f_num, p, cur_row, cur_num;
    uint32_t (*pos_out)[2] = sub->aln_out->out_buf+out_len;
    /*  
    fnd_i = 0;
    for(i = 0; i < sub->delta; ++i){
        if(found[i][1] > 0) {
            found[fnd_i][0] = found[i][0];  
            found[fnd_i][1] = found[i][1];  
            found[fnd_i][2] = found[i][2];  
            found[fnd_i][3] = found[i][3];  
            fnd_i++;
        }
    }
    fnd_num = fnd_i;
    int pos_i = 0;
    for(fnd_i = 0; fnd_i < fnd_num; ++fnd_i){
        //f_score = found[fnd_i][0];
        //f_num = found[fnd_i][1];
        //f_ed = found[fnd_i][3]; 
        f_bg = found[fnd_i][2];
        p = f_bg;
        while(1){
            cur_num = out_buf[p][0];
            for(j = 0; j < cur_num; ++j) {
                pos_out[pos_i][0] = out_buf[p+j+1][0];  
                pos_out[pos_i][1] = out_buf[p+j+1][1];  
                pos_i++;   
            }
            if(out_buf[p][1] > 0) {
                p = out_buf[p][1]; 
            } else {
                break;
            }
        }
    }
    //------------------------------------- 
    int n_aln = 0;    
    for(i = 0; i < pos_i; ++i) { 
        uint32_t pos = pos_out[i][0]; 
        if( pos != query->pos) {
            if(n_aln == 0) fprintf(fp_sam, "\tXA:Z:");
            int rid;
            bns_coor_pac2real(bns, pos, query->l_seq, &rid); 
            fprintf(fp_sam, "%s,%u;",bns->anns[rid].name, pos-bns->anns[rid].offset+1);
            //fprintf(stderr,"pos = %u, ref = %s,pos = %u\n",aln_pos[i], bns->anns[rid].name, aln_pos[i]-bns->anns[rid].offset+1);
            //fprintf(fp, "%u;", aln_pos[i]+1);
            n_aln++;
        }
    }
    */
    fprintf(fp_sam, "\n");
    //printf("%s", str);
    return len_str;
}

int query_aln2sam2(query_t *query,  idx_t *idx, struct SubBuf *sub,  char *str)
//int query_aln2sam2(query_t *query,  idx_t *idx, struct SubBuf *sub, FILE *fp_sam)
{
    //char str_buf[4000] = {};
    //char *str = str_buf; 
    int len_str = 0;
    bntseq_t *bns = idx->bns;

    int i, j;
    unsigned short flag = 0;
    int nm = 255; 
//fflush(fp);
    len_str += sprintf(str, "%s\t", query->qseq->name); 
    //flag Rname pos mapq cigar Rnext Pnext tlen
    uint32_t pos = query->pos;
    if(query->strand == 1) flag |= FLAG_REV;
    if(query->pos ==(uint32_t)-1) {
        flag |= FLAG_UNMAPPED;
        len_str += sprintf(str+len_str, "%u\t*\t0\t255\t*\t*\t0\t0\t",flag); 
    } else{
        uint8_t target[LEN_READ+100];


        uint32_t l_pos = pos > MAX_CLIP?pos-MAX_CLIP:0;
        uint32_t r_pos = pos+query->l_seq+MAX_CLIP;
fprintf(stderr, "%u, ref_len = %u, pos = %u, clip = %u, query->l_seq = %d, l_pos = %u, r_pos = %u\n", __LINE__, idx->bwt->seq_len, pos, MAX_CLIP, query->l_seq, l_pos, r_pos);
        for(i =0; i < r_pos - l_pos; ++i) target[i] = __get_pac(idx->pac, l_pos+i);
        uint8_t *seq = query->strand==1?query->rseq:query->seq;
        kswq_t **qry, **kswq; 
        if(query->strand) {
            qry = &sub->qry_r[0];
            kswq = sub->kswq_r; 
        } else{
            qry = &sub->qry_f[0];
            kswq = sub->kswq_f; 
        }
fprintf(stderr, "%u, ref_st = %d, ref_end = %d, seq_st = %d, seq_end = %d\n", __LINE__, query->ref_start, query->ref_end, query->seq_start, query->seq_end); 
//fprintf(stderr, "%u, kswq = %p\n", __LINE__, kswq); 

        if(query->seq_start == -1 && query->ref_start == -1) {
        //if(1) {
            kswr_t r = ksw_align2(query->seq_end-1, seq, query->ref_end-1, target, 5, sub->mat, 6, 1, 0, qry, kswq, query->b0);
            query->ref_start = r.tb;
            query->seq_start = r.qb;
        }
fprintf(stderr, "%u, ref_st = %d, ref_end = %d, seq_st = %d, seq_end = %d\n", __LINE__, query->ref_start, query->ref_end, query->seq_start, query->seq_end); 
        int score = ksw_global(query->seq_end-query->seq_start, seq+query->seq_start, query->ref_end-query->ref_start, target+query->ref_start, 5, sub->mat, 6, 1, 50, &query->n_cigar, &query->cigar, &sub->kswgb);
        
        if(score != query->b0) {
            fprintf(stderr, "%u, b0 = %d, sc = %d\n", __LINE__, query->b0, score); 
        }
        //query_add_md(target, query->ref_start, query->ref_end, seq, query->seq_start, query->seq_end, query->n_cigar, query->cigar);
        uint8_t *x0 = target+query->ref_start, *x1 = seq+query->seq_start;

        nm = 0;
        for(i = 0; i < query->n_cigar; ++i){
            int op = query->cigar[i] & 0xF;
            int len = query->cigar[i] >>4;
            if(op == 0) {
                for(j =0; j < len; ++j, ++x0, ++x1) {
                    if(*x0 != *x1) ++nm; 
                }
            } else if(op == 1) {
                nm += len;
                x1 += len; 
            } else{
                nm += len;
                x0 += len;
            }
        }

   
        if(query->ref_start <= MAX_CLIP) {
            pos -= MAX_CLIP-query->ref_start; 
            query->ref_start = 0;
        } else{
            pos += query->ref_start-MAX_CLIP;
            query->ref_start -= MAX_CLIP; 
        }
        
        uint32_t rid;
        //pos = bns->l_pac - 1 - pos - (query->seq_end - query->seq_start);
        bns_coor_pac2real(bns, pos, query->l_seq, &rid); 
        len_str += sprintf(str+len_str, "%u\t%s\t%u\t255\t",flag, bns->anns[rid].name, pos-(uint32_t)bns->anns[rid].offset+1);
        if(query->seq_start != 0) len_str += sprintf(str+len_str, "%uS", query->seq_start);
        for(i = 0; i < query->n_cigar; ++i){ 
            if((i == 0|| i == query->n_cigar-1) && (query->cigar[i] &3) ==2 ) {
                //query->b0 += 7;
                //continue;
            }
            len_str += sprintf(str+len_str, "%u%c", query->cigar[i]>>4, "MID"[query->cigar[i]&3]);
        }
        if(query->seq_end != query->l_seq) len_str += sprintf(str+len_str, "%uS", query->l_seq-query->seq_end);
        len_str += sprintf(str+len_str, "\t*\t0\t0\t");
        //if(query->cigar != NULL )free(query->cigar);


    }
    if(query->strand == 1) {
        for(i = 0; i < query->l_seq; ++i){ len_str += sprintf(str+len_str, "%c", "ACGTN"[query->rseq[i]]);}
        len_str += sprintf(str+len_str, "\t");
        for(i = query->l_seq-1; i >=0; --i) { len_str += sprintf(str+len_str, "%c", query->qual[i]);}
    } else{ 
        for(i = 0; i < query->l_seq; ++i) { len_str += sprintf(str+len_str, "%c", "ACGTN"[query->seq[i]]); }
        len_str += sprintf(str+len_str, "\t");
        for(i = 0; i < query->l_seq; ++i) { len_str += sprintf(str+len_str, "%c", query->qual[i]); }
    }
    len_str += sprintf(str+len_str, "\tAS:i:%d\tNM:i:%d", query->b0, nm);
    //+++++++++++++++++++++++++++++++++++++
fprintf(stderr, "%u, %s\n", __LINE__, __func__ );

    uint32_t (*out_buf)[2] = sub->aln_out->out_buf; 
    int (*found)[4] = sub->aln_out->found; 
    int out_len = sub->aln_out->len;
    int out_num = sub->aln_out->num; 
 
    int fnd_i, fnd_num, f_score, f_bg, f_ed, f_num, p, cur_row, cur_num;
    uint32_t (*pos_out)[2] = sub->aln_out->out_buf+out_len;
    
    int delta = sub->delta;
    //int delta = 1;
    fnd_i = 0;
    for(i = 0; i < delta; ++i){
        if(found[i][1] > 0) {
            found[fnd_i][0] = found[i][0];  
            found[fnd_i][1] = found[i][1];  
            found[fnd_i][2] = found[i][2];  
            found[fnd_i][3] = found[i][3];  
            fnd_i++;
        }
    }
    fnd_num = fnd_i;
    int pos_i = 0;
    for(fnd_i = 0; fnd_i < fnd_num; ++fnd_i){
        //f_score = found[fnd_i][0];
        //f_num = found[fnd_i][1];
        //f_ed = found[fnd_i][3]; 
        f_bg = found[fnd_i][2];
        p = f_bg;
        while(1){
            cur_num = out_buf[p][0];
            for(j = 0; j < cur_num; ++j) {
                pos_out[pos_i][0] = out_buf[p+j+1][0];  
                pos_out[pos_i][1] = out_buf[p+j+1][1];  
                pos_i++;   
            }
            if(out_buf[p][1] > 0) {
                p = out_buf[p][1]; 
            } else {
                break;
            }
        }
    }
    //------------------------------------- 
    int n_aln = 0;    
fprintf(stderr, "%u, pos_i = %d\n", __LINE__, pos_i);

    for(i = 0; i < pos_i; ++i) { 
if(n_aln >= 499) break;
        uint32_t pos = pos_out[i][0]; 
        int sc = pos_out[i][1]%1024; 
        if( pos != query->pos) {
            if(n_aln == 0) len_str += sprintf(str+len_str, "\tXA:Z:");
            int rid;
            //pos = bns->l_pac - 1 - pos - query->l_seq;
            bns_coor_pac2real(bns, pos, query->l_seq, &rid); 
            len_str += sprintf(str+len_str, "%s,%u,%d;",bns->anns[rid].name, (uint32_t)(pos-bns->anns[rid].offset+1), sc);
            //fprintf(stderr,"pos = %u, ref = %s,pos = %u\n",aln_pos[i], bns->anns[rid].name, aln_pos[i]-bns->anns[rid].offset+1);
            //fprintf(fp, "%u;", aln_pos[i]+1);
            n_aln++;
        }
    }
   
    len_str += sprintf(str+len_str, "\n");
    //fprintf(fp_sam, "%s", str_buf);
    return len_str;
}

void query_aln2sam(query_t *query, int n_pos, uint32_t *aln_pos,  idx_t *idx, struct SubBuf *sub,  FILE *fp)
{
    bntseq_t *bns = idx->bns;
    int i;
    unsigned short flag = 0;
    int nm = 255; 
    fflush(fp);
//fprintf(fp, "%s\t", query->qseq->name); 
    //flag Rname pos mapq cigar Rnext Pnext tlen
    if(query->strand == 1) flag |= FLAG_REV;
    if(query->pos ==(uint32_t)-1) {
        flag |= FLAG_UNMAPPED;
        fprintf(fp, "%u\t*\t0\t255\t*\t*\t0\t0\t",flag); 
    } else{
        uint8_t target[LEN_READ+50];

        uint32_t pos = query->pos;
        uint32_t l_pos = pos > MAX_CLIP?pos-MAX_CLIP:0;
        uint32_t r_pos = pos+query->l_seq+MAX_CLIP;
        for(i =0; i < r_pos - l_pos; ++i) target[i] = __get_pac(idx->pac, l_pos+i);

        
       
        uint8_t *seq = query->strand==1?query->rseq:query->seq;
//fprintf(stderr, "%u, pos = %u\n", __LINE__, query->pos);
        int score = ksw_global(query->seq_end-query->seq_start, seq+query->seq_start, query->ref_end - query->ref_start, target+query->ref_start, 5, sub->mat, 6, 1, 50, &query->n_cigar, &query->cigar, &sub->kswgb);
        
        uint8_t *x0 = target+query->ref_start, *x1 = seq+query->seq_start;
        int j;
        nm = 0;
        for(i = 0; i < query->n_cigar; ++i){
            int op = query->cigar[i] & 0xF;
            int len = query->cigar[i] >>4;
            if(op == 0) {
                for(j =0; j < len; ++j, ++x0, ++x1) {
                    if(*x0 != *x1) ++nm; 
                }
            } else if(op == 1) {
                nm += len;
                x1 += len; 
            } else{
                nm += len;
                x0 += len;
            }
        }

        
        if(query->ref_start <= MAX_CLIP) {
            query->pos -= MAX_CLIP-query->ref_start; 
            query->ref_start = 0;
        } else{
            query->pos += query->ref_start-MAX_CLIP;
            query->ref_start -= MAX_CLIP; 
        }
        
        uint32_t rid;
        bns_coor_pac2real(bns, query->pos, query->l_seq, &rid); 

//fprintf(fp, "%u\t%s\t%u\t255\t",flag, bns->anns[rid].name, query->pos-bns->anns[rid].offset+1);
//if(query->seq_start != 0) fprintf(fp, "%uS", query->seq_start);
        for(i = 0; i < query->n_cigar; ++i){ 
            if((i == 0|| i == query->n_cigar-1) && (query->cigar[i] &3) ==2 ) {
                query->b0 += 7;
                continue;
            }
//fprintf(fp, "%u%c", query->cigar[i]>>4, "MID"[query->cigar[i]&3]);
        }
//if(query->seq_end != query->l_seq) fprintf(fp, "%uS", query->l_seq-query->seq_end);
//fprintf(fp, "\t*\t0\t0\t");
        //if(query->cigar != NULL )free(query->cigar);


    }
    if(query->strand == 1) {
        for(i = 0; i < query->l_seq; ++i) {
//fprintf(fp, "%c", "ACGTN"[query->rseq[i]]);
        }
//fprintf(fp, "\t");
        for(i = query->l_seq-1; i >=0; --i){ 
//fprintf(fp, "%c", query->qual[i]);
        }
    } else{ 
        for(i = 0; i < query->l_seq; ++i) {
//fprintf(fp, "%c", "ACGTN"[query->seq[i]]);
        }
//fprintf(fp, "\t");
        for(i = 0; i < query->l_seq; ++i) {
//fprintf(fp, "%c", query->qual[i]);
        }
    }
//fprintf(fp, "\tAS:i:%d\tNM:i:%d", query->b0, nm);
    int n_aln = 0;    
    //fprintf(stderr, "%u, n_pos = %u\n", __LINE__, n_pos);
    for(i = 0; i < n_pos; ++i) { 
        if( aln_pos[i] != query->pos) {
//if(n_aln == 0) fprintf(fp, "\tXA:Z:");
            int rid;
            bns_coor_pac2real(bns, aln_pos[i], query->l_seq, &rid); 
//fprintf(fp, "%s,%u;",bns->anns[rid].name, aln_pos[i]-bns->anns[rid].offset+1);
            //fprintf(stderr,"pos = %u, ref = %s,pos = %u\n",aln_pos[i], bns->anns[rid].name, aln_pos[i]-bns->anns[rid].offset+1);
            //fprintf(fp, "%u;", aln_pos[i]+1);
            n_aln++;
        }
    }
//fprintf(fp, "\n");
    return;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  query_init
 *  Description:  query init func
 * =====================================================================================
 */


void query_destroy(query_t *query)
{
    //fprintf(stderr, "%s\n", __func__);

    if(query->qseq->name != NULL) free(query->qseq->name);
    if(query->qseq->comment != NULL) free(query->qseq->comment);
    if(query->seq != NULL) free(query->seq);
    if(query->rseq != NULL) free(query->rseq);
    if(query->qual != NULL) free(query->qual);
/*
    if(query->cigar != NULL) {
        free(query->cigar->s);
        free(query->cigar);
    }
    if(query->sam != NULL ) {
        free(query->sam->s);
        free(query->sam);
    }
    kv_destroy(query->hits[0]);
    kv_destroy(query->hits[1]);
*/

    //free(query->hits[0].a);
    //free(query->hits[1].a);
    //kv_destroy(query->ordered_index);
    //kv_destroy(query->unique_index);
}
queryio_t *query_open(const char *fn_fa)
{
    //fprintf(stderr, "%s\n", __func__);

    gzFile fp = 0;
    kseq_t *ks;
    queryio_t *qs = (queryio_t *)calloc(1, sizeof(queryio_t));
    if(qs == NULL){
    	fprintf(stderr, "[query_init]: allocate mem fail!\n");
    	exit(EXIT_FAILURE);
    }
    fp = gzopen(fn_fa, "r");
    if(fp == Z_NULL){
        fprintf(stderr, "[query_open]: file %s open fail!\n", fn_fa);    
        exit(1);
    }
    ks = kseq_init(fp);
    qs->kseq = ks;
    qs->fp = fp;

    return qs;
}
void query_close(queryio_t *qs)
{    
    //fprintf(stderr, "%s\n", __func__);
    
	kseq_destroy(qs->kseq);
    gzclose(qs->fp);
	free(qs);
}


/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  query_read
 *  Description:  read query 
 * =====================================================================================
 */

static inline void trim_readno(kstring_t *s)
{
	if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
		s->l -= 2, s->s[s->l] = 0;
}
#define INIT_HIT_NUM 32
#define SEED_LEN 20
int query_read_seq(queryio_t *qs, qseq_t *qseq)
{
    int i;
    kseq_t *kseq = qs->kseq;      
    int l_seq = kseq_read(kseq);   
    if(l_seq <= 0) return 0;    
    //read name and comment
    if(kseq->name.s != NULL) {
        trim_readno(&kseq->name);
        memcpy(qseq->name, kseq->name.s, kseq->name.l);
    }      
    if(kseq->comment.s != NULL) {
        memcpy(qseq->comment, kseq->comment.s, kseq->comment.l);
    }
    //read seq and generate reverse seq
    qseq->l_seq = l_seq;
    qseq->l_trim5 = 0;
    qseq->l_trim3 = 0; 
    int err_raise_pos = l_seq*3/5;
    int n_N = 1;
    int N_freq = 20;
    int n_N0 = 0, n_pos = 0, MAX_N = l_seq/N_freq *3/5;
    int flag = 0; 
    qseq->apos[0] = 0;
    for(i = 0; i < l_seq; ++i){
        unsigned nt = nst_nt4_table[(uint8_t)kseq->seq.s[i]];
        if(nt>3) { 
            qseq->apos[n_N++] = i;
            qseq->seq[i] = nt;
            qseq->rseq[l_seq-i-1] = nt;
        } else{
            qseq->seq[i] = nt; 
            qseq->rseq[l_seq-i-1] = 3-nt; 
        }
        if(i == err_raise_pos) {
            if(n_N * N_freq > i) {
                qseq->l_trim3 = l_seq; 
            } 
            n_N0 = n_N;
            n_pos = i;
        }

        if( i > err_raise_pos && nt > 3 && qseq->l_trim3 == 0){ 
            if(i - qseq->apos[n_N-1] > 2*SEED_LEN){
                n_N0 = n_N;
                n_pos = i;
            } else if( n_N > 5 && n_N - n_N0 > 1 && (i - n_pos < (n_N - n_N0) *N_freq)) {
            
                if(i - qseq->apos[n_N-1] > SEED_LEN){
                    qseq->l_trim3 = l_seq - qseq->apos[n_N]; 
                } else {
                    qseq->l_trim3 = l_seq - qseq->apos[n_N-1]; 
                }
            }
        }  
    }
    qseq->n_ambiguous = n_N - 1;
    if(kseq->qual.l > 0){
        memcpy(qseq->qual, kseq->qual.s, kseq->qual.l);
        //convert to sam base qual
    }
    return l_seq;
}		/* -----  end of function query_read  ----- */
int query_read_multiSeqs(queryio_t *qs, int n_seq, qseq_t *multiSeqs)
{
    int i=0;
    //fprintf(stderr, "%s\n", __func__);

    while(i < n_seq && query_read_seq(qs, multiSeqs+i) >0  ){
       ++i;
    }

    return i;
}

int query_read_multiPairedSeqs(queryio_t *qs[], int n_seq, qseq_t *multiSeqs)
{
    int i,j;
    i = 0, j = 1; 
    while(i < n_seq  && query_read_seq(qs[0], multiSeqs +i) >0 ){
        i += 2;
    }
    while(j < n_seq && query_read_seq(qs[1], multiSeqs+j) >0  ){
        j += 2;
    }
    if(i != j-1) {
        fprintf(stderr, "[%s]: the number of paired seqs are not same!\n", __func__);
        exit(1);
    }
    return i;

}
//[to do]: change
uint32_t gen_mapq(uint32_t b0, uint32_t b1)
{
    //uint32_t mapq = -4.343 * log(1-(double)abs(b0 - b1)/(double)b0);
    //mapq = (uint32_t)(mapq+4.99);
    if(b0 == 0) return 0;
    double a = 255.0;
    uint32_t mapq = a*((double)abs(b0-b1)/(double)b0);

     
    mapq = mapq <254?mapq:254;
    return mapq; 
}
/*
void query_gen_cigar(uint32_t l_ref, const uint32_t *mixRef, query_t *query)
{
    query->seq_start = 0;
    query->seq_end = query->l_seq-1;
    if(query->pos != 0xFFFFFFFF){
        if(query->is_gap){//have gap
           int ret = ed_diff_withcigar(mixRef, query->pos, query->l_seq+4, query->strand==0?query->seq:query->rseq, query->l_seq, query->n_diff, query->cigar->s, query->cigar->m, 1, COMPACT_CIGAR_STRING);
            if(ret == -1) fprintf(stderr, "[%s]:Erro:%s\n", __func__, query->qseq->name);
            //fprintf(stderr, "[cigar]: %s\n", query->cigar->s);           
        } else ksprintf(query->cigar, "%dM", query->l_seq);
    
    
    } //else query->cigar = NULL;
   
}*/
/*
void query_set_hits(query_t *query, int max_hits, hits_t *hits0, hits_t *hits1)
{
    int i, j;
    uint32_t primary_pos = query->pos;
    int tot_hits = 0; 
    query->b0 = query->n_diff;
    query->b1 = 100000;
    for(i = 0; i < 2; ++i){
        int n = 0;
        hits_t *hits = i==0?hits0:hits1; 
        hit_t *a = hits->a;
        uint32_t last_pos = (uint32_t)-1; 
        for(j = 0; j < hits->n; ++j){

            if(a[j].strand != i){
                fprintf(stderr, "[%s]: %s hits[%d] strand %d\n", __func__, query->qseq->name, i, a[j].strand);
            
            }

            uint32_t pos = a[j].pos;
            if (pos == last_pos) continue;
            if(a->n_diff <= query->n_diff){ 
                if(a->n_diff <= query->b1) query->b1 = a->n_diff;
                kv_push(hit_t, query->hits[i], a[j]); 
                tot_hits += 1;    
            }
            if(tot_hits == max_hits){ 
                goto end;
            }
        }  

    }
end:{
        query->mapq = gen_mapq(query->b0, query->b1);   
    
    }
}*/
