// 循环控制信息--------------------------------------------
//uint32_t queue_align_info[][2];  //开辟较大内存
//uint32_t this_stack_tree[NUM_EXT][5];  
//[i][0]是第i个节点(即，第i级索引）在queue_align_info中的起始行号
//[i][1]是第i个节点(即，第i级索引）在queue_align_info中的个数
//uint32_t this_stack_node[3];


//------------------------------------------------------------------

//void setFileName(char *in_fname, struct FileName *out_fname,int flg );
#ifndef __SEEDEXTLIB_H
#define __SEEDEXTLIB_H
#include "debug.h"
#include <stdio.h>
#include "ksw.h"
#include "seed.h"
#include "query.h"
#include "ksort.h"


typedef struct{
    int16_t seq;
    int16_t idx;
} pair_t;
#define ks_bg_pair(a, b) ((a).seq==(b).seq?(a).idx>(b).idx:(a).seq>(b).seq)
KSORT_INIT(pair_t, pair_t, ks_bg_pair)

typedef struct{
    uint32_t x[5];
}tuple5_t;
#define ks_lt_tuple5(a, b) ((a).x[0]==(b).x[0]?(a).x[1]<(b).x[1]:(a).x[0]<(b).x[0])
KSORT_INIT(tuple5_t, tuple5_t, ks_lt_tuple5)

#define NO_BWT_SUM_SIZE 64
#define ALN_THRES_SCORE 30
#define MIN_BWT_SIZE 16
/*  
#define MAX_CANDI_POS 100000
#define SW_HSH_BIT 22
#define SW_RLT_BIT 2
#define SW_DAT_BIT 8
#define SW_TO_DAT_SIZE 4 
*/
#define WOREST_PERCENT 30
#define SW_HSH_BIT 20
#define SW_RLT_BIT 4
#define SW_DAT_BIT 8
#define SW_TO_DAT_SIZE 4 
#define __MAX_CANDI_POS SW_THRES*(LEN_READ/SEED_LEN+1) 
#define MAX_CANDI_POS ( (__MAX_CANDI_POS)< 100000 ? 100000:(__MAX_CANDI_POS) )
#define THRES_ALN_OUT 20

struct SubBuf *InitIdxsArry(struct FileName *fn, struct JmpMod **out_jmp, struct CapIfo *cap, struct ExtBlck out_blck[], struct StackTree *sTree, struct SubBuf *sub);
void initStackTree(struct StackTree *sTree);

void setStackTree(struct StackTree *sTree,struct SubBuf *sub,seed_t *seed, int seed_id, struct ExtBlck*eB);

void setBlckData(struct ExtBlck *in_blck, struct SubBuf *sub, struct StackTree *sTree, uint8_t *read_seq, int seed_off, struct ExtBlck **out_blck);
void getSeed(uint8_t *seed_seq, uint8_t *read_seq, uint8_t  seed_off);
void getCapPos(struct JmpMod *jmp,uint32_t buf_algn[]);
int Algn_auto_r2l(idx_t *fm_idx, uint32_t hash_boundry[], query_t *query, struct SubBuf *sub,uint32_t best_pos[]);
//int AlgnPos(idx_t *fm_idx, uint8_t ext_cls, query_t *query, int seed_off, uint32_t pos_buf[], int pos_num, struct SubBuf *sub);

int AlgnPos(idx_t *fm_idx, query_t *query, int pos_num, struct SubBuf *sub);
int init_seed_model(seed_t *seed, uint32_t ref_len);

uint32_t aln_seed_seq(idx_t *fm_idx,uint32_t hash_boundry[], int read_len,  uint8_t *f_read_seq,  uint8_t *r_read_seq, uint8_t *qual, uint8_t* rqual, seed_t *seed );

uint32_t aln_seed_seq_1(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed );
uint32_t aln_seed_seq_2(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed, int seed_slc_buf[] );


int AlgnPos_buf(idx_t *fm_idx, struct ExtBlck *eBlck, query_t *query, seed_t *seed, int num_pair, uint8_t ext_cls,   struct SubBuf *sub);
//int AlgnPos_buf(idx_t *fm_idx, uint8_t ext_cls, query_t *query, int seed_off, int num_pair, struct ExtBlck* eBlck, struct SubBuf *sub);
uint32_t __dna2_count(uint8_t *bwt, uint32_t k , uint32_t l, uint8_t c);
void OutAlgnInfo(uint32_t (*algn_out)[5]);
int AlgnmExtSeq(struct ExtBlck *eB, struct SubBuf *sub, int flg);		
int Algn_sub_once(struct ExtBlck *eB, struct SubBuf *sub, int flg);		
uint32_t __dna2_count_small(uint8_t *bwt, uint32_t k , uint32_t l, uint8_t c);
//uint32_t getBwtSumSize(uint32_t DataNum);
void getExtSeq(struct SubBuf *sub, uint32_t cls, uint8_t *read_seq);
void gen_cnt(uint8_t *bwt, int n_data, uint8_t cnt_2nt[][17]);
int align_min(uint8_t *bwt, int n_data, uint8_t seq[16], int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4], struct SubBuf *sub);
int align_nosum(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t seq[16], int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4], struct SubBuf *sub);
int align_255(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data,uint8_t seq[16], int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4], struct SubBuf *sub);
int align_large(uint8_t *Bwt, uint8_t *cnt2,  uint32_t n_data, const uint8_t *seq , int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4], struct SubBuf *sub);
int Algn_sub_all(struct ExtBlck *eB, struct SubBuf *sub,  int flg);
int align_min_indel(uint8_t *bwt, int n_data, uint8_t seq[18], int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4], struct SubBuf *sub);
int align_nosum_indel(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t seq[18], int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4], struct SubBuf *sub);
int align_255_indel(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data,uint8_t seq[18], int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4],  struct SubBuf *sub );
int align_large_indel(uint8_t *Bwt, uint8_t *cnt2,  uint32_t n_data, const uint8_t *seq , int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4], struct SubBuf *sub);
int Algn_indel_all(struct ExtBlck *eB, struct SubBuf *sub,  int flg);
int Algn_indel_once(struct ExtBlck *eB, struct SubBuf *sub,  int flg);
int gen_aln_data(struct ExtBlck *eB, struct SubBuf *sub, struct call_t *call); 
int set_pair_data(struct SubBuf *sub, struct call_t *call); 
int get_next_aln(struct SubBuf *sub, struct call_t *call); 
int sortTreeData(struct StackTree *sTree);
int get_sw_ed_val(uint32_t pos, struct SubBuf *sub, query_t *query);

int slc_aln_pos(idx_t *fm_idx, query_t *query, int *seq_off, int pos_num, struct SubBuf *sub);
int bsearch_idx_R(idx_t *fm_idx, query_t *query, int len_buf[2], int seq_off[2], uint32_t bwt_idx[4])
{
    uint32_t bg = bwt_idx[0];
    uint32_t ed = bwt_idx[1];
    int s_off[3];
    //len_buf[0]是已经比对的长度;
    //len_buf[1]是将要比对的长度。
    //seq_off[1]是read中将要比对的开始位置
    int len_ext = len_buf[1];
    if(len_ext > 32) len_ext = 32;  
    s_off[0] = len_buf[0];
    s_off[2] = seq_off[1];
    if(seq_off[1] + len_ext > query->l_seq) {
        s_off[1] = len_buf[0] + query->l_seq - seq_off[1];
    } else {
        s_off[1] = len_buf[0] + len_ext; 
    }

    int len = s_off[0];

    uint32_t num;
    int l0 = 0, l1 = 0, n;
    l0 = get_ext_idx_R(fm_idx, query, s_off, bwt_idx); 
    len_ext = l0;
    if(l0 > 0) { 
    //if(0) { 
        bg = bwt_idx[2]; 
        ed = bwt_idx[3]; 
        num = ed + 1 - bg;
        len += l0; 
       
        //if(l0 == 32) {
        if(0) {
            s_off[0] = len;
            s_off[1] = s_off[0] + 32;
            s_off[2] = seq_off[1] + l0;
            bwt_idx[0] = bwt_idx[2];
            bwt_idx[1] = bwt_idx[3];
           
            
            int l_seq = query->l_seq;
            if(num > IS_SMLSIZ && seq_off[1]+16 < l_seq) {
                l1 = get_ext_idx_R(fm_idx, query, s_off, bwt_idx); 
            } else if(num > IS_SMLSIZ*5 && seq_off[1]+8 < l_seq){
                l1 = get_ext_idx_R(fm_idx, query, s_off, bwt_idx); 
            } else if(num > SW_THRES && seq_off[1]+4 < l_seq) {
                l1 = get_ext_idx_R(fm_idx, query, s_off, bwt_idx); 
            }
            if (l1 > 0) {
                bg = bwt_idx[2]; 
                ed = bwt_idx[3]; 
                num = ed + 1 - bg;
                len += l1; 
                //seq_off[1] += l1;
                len_ext += l1;
            }
        } 

        
    }
    bwt_idx[2] = bg;
    bwt_idx[3] = ed;
    return len_ext;
}



int find_R_exact_range(idx_t *fm_idx, query_t *query, int R_offset, uint32_t *pos_buf, int num, uint32_t out[2]) { 
    int i;
    uint32_t bg, ed, row;
    uint32_t tmp = 0, seq = 0, bg_val, ed_val, cur_val;
    uint32_t bg_row, ed_row, bot, pos;
    R_offset += 16;
    int len = query->l_seq - R_offset;
    //fprintf(stderr, "%u, R_offset = %u, query->l_seq = %u\n", __LINE__, R_offset, query->l_seq);

    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    for(i = R_offset; i < query->l_seq; ++i){
        seq <<= 2;
        seq |= read_seq[i] >3?0:read_seq[i];
    }
    //-----------------------------------
    pos = pos_buf[0] + R_offset; 
    tmp = 0;
    for(i = 0; i < 10; ++i){
        tmp <<=2;
        tmp |= __get_pac(fm_idx->pac, pos+i);
    }
    bg_val = tmp; 
    bg = 0;
    //-------------------------------------
    pos = pos_buf[num-1] + R_offset; 
    tmp = 0;
    for(i = 0; i < len; ++i){
        tmp <<=2;
        tmp |= __get_pac(fm_idx->pac, pos+i);
    }
    ed_val = tmp;    
    ed = num -1;
    bot = ed;
//fprintf(stderr, "%u, bg_val = %u, ed_val = %u, seq = %u\n", __LINE__, bg_val, ed_val, seq);
    if(bg_val == ed_val ||bg_val > seq || seq> ed_val ) {
        if(bg_val == ed_val ){
            bg_row = 0;
            ed_row = num -1;
            if(bg_val != seq) return -1;
            else return 1;
        } else {
            bg_row = num;
            ed_row = num;
            return -1; 
        }
    } 
    row = (bg+ed)/2;
    bg_row = num + 1;
    ed_row = num + 1;   
//fprintf(stderr, "%u\n", __LINE__);

    while(1){
        while(bg <= ed) { 
            pos = pos_buf[row] + R_offset; 
            tmp = 0;
            for(i = 0; i < len; ++i){
                tmp <<=2;
                tmp |= __get_pac(fm_idx->pac, pos+i);
            }
            cur_val = tmp;    
            //-------------------------------------
            if(cur_val < seq){
                bg = row + 1; 
            } else if(cur_val > seq){
                if(bg == ed) break;
                ed = row;
                bot = row; 
            } else{// cur_val == seq
                if(ed_row == num +1) ed_row = row; 
                break; 
            }
            row = (bg+ed)/2;
        } //end while(bg <= ed) 
        if(bg >= ed) break;
        ed = row;         
        row = (bg+ed)/2; 
    } //end while(1)
    if(ed_row > num) return -1;
    else bg_row = row;
    //---------------------------------------------
   
    bg = ed_row;
    ed = bot; 
    row = (bg+ed)/2;
//fprintf(stderr, "%u, num = %u, bg_row = %u, ed_row = %u\n", __LINE__, num, bg_row, ed_row); 
//fprintf(stderr, "%u, bg = %u, ed = %u, row = %u\n", __LINE__, bg, ed, row); 


    while(2){
        while(bg <= ed) { 
            pos = pos_buf[row] + R_offset; 
            tmp = 0;
            for(i = 0; i < 10; ++i){
                tmp <<=2;
                tmp |= __get_pac(fm_idx->pac, pos+i);
            }
            cur_val = tmp;    
            //-------------------------------------
            if(cur_val < seq){
                bg = row + 1; 
            } else if(cur_val > seq){
                if(bg == ed) break;
                ed = row; 
            } else{// cur_val == seq
                break; 
            }
            row = (bg+ed)/2;
        } //end while(bg <= ed) 
        if(cur_val == seq) ed_row = row; 
        if(bg >= ed) break;
        bg = row +1;
        row = (bg+ed)/2; 
    }// end while(2)+++++++++++++++
    ed_row = row;
    out[0] = bg_row;
    out[1] = ed_row;

    return ed_row + 1 - bg_row;
}
int classify_pos_0(idx_t *fm_idx, query_t *query, seed_t *seed, uint32_t pos_buf[], int pos_num, int seq_off[2], struct SubBuf *sub){
    int sn = query->seed_num;
    int sub_err = sub->sub_err;
    int l_seq = query->l_seq;
    int thres_pos_num = 0;
    
if(query->b0 < query->candi_thres + 20) {
        thres_pos_num = IS_SMLSIZ; 
    } else if(query->b0 > l_seq - sn*sub_err){
        thres_pos_num = 0; 
    } else {
        thres_pos_num = IS_SMLSIZ*(l_seq - query->b0 - sn*sub_err)/l_seq; 
    }
    int thres_sw_to = sub->thres_sw_to;
    int cur_cls = seed->cls;
    int max_cls = query->max_ext;
    uint32_t pos;
    int pos_j = 0;
    int pi;
    for(pi = 0; pi < pos_num; ++pi) {
        pos = pos_buf[pi];
        if(pos > fm_idx->bwt->seq_len) {
            printf("%u, pos = %u\n", __LINE__, pos); 
            exit(1);
        }
        int flg_ed =  get_sw_ed_val(pos, sub, query);

fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
        if(flg_ed > 0) continue;
        if(pos_num <= thres_pos_num ) {
            pos_buf[pos_j++] = pos;
            set_sw_ed_pos(pos, query, sub);
            continue;
        }
        if(cur_cls == max_cls) {
            pos_buf[pos_j++] = pos;
            set_sw_ed_pos(pos, query, sub);
            continue;
        } 

fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
        int cl_num = 0; 
        cl_num = get_set_sw_to_pos(pos, 12*2, 0, seed, query, sub);
        if(cl_num > thres_sw_to) {
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d, thres_sw_to = %d, cl_num = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls, thres_sw_to, cl_num);
            pos_buf[pos_j++] = pos;
            set_sw_ed_pos(pos, query, sub);
        } else if(cl_num > 0){
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
            int sc_num = 0;
            int olp_flg = sub->olp_flg;
            //if(olp_flg > 0) 
                sc_num = eval_seed_olp(fm_idx,  query, seed, pos, seq_off,  sub);    
            int thres_sw_olp = sub->thres_sw_olp; 
            if(sc_num >= thres_sw_olp || sc_num + cur_cls >= max_cls -1) {
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
                pos_buf[pos_j++] = pos;
                set_sw_ed_pos(pos, query, sub);
            } else{
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
                //get_set_sw_to_pos(pos, 0, 1, seed, query, sub);// set_sw_to_pos
                get_set_sw_to_pos(pos, 0, 2, seed, query, sub);// set_sw_to_pos
            }
        } else{ //cl_num == 0
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
            int n_sm;
            if(query->is_rev == 0) n_sm = sub->n_sm_f;
            else n_sm = sub->n_sm_r;
            //if(n_sm < 3 || cur_cls > 2) {
            if(n_sm < 4 || cur_cls > 2) {
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
                int sc_num = 0;
                int thres_sw_olp = sub->thres_sw_olp; 
                int olp_flg = sub->olp_flg;
                //if(olp_flg > 0) 
                sc_num = eval_seed_olp(fm_idx,  query, seed, pos, seq_off,  sub);    
fprintf(stderr, "%u, sc_num = %d, thres_sw_olp = %d, max_cls = %d, cur_cls = %d\n", __LINE__, sc_num, thres_sw_olp, max_cls , cur_cls);
                if(sc_num >= thres_sw_olp + 1|| sc_num + cur_cls > max_cls + 1) {
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d, pos_j = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls, pos_j);
                    pos_buf[pos_j++] = pos;
                    set_sw_ed_pos(pos, query, sub);
                } else{
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d, seed_id = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls, seed->id);
                    //get_set_sw_to_pos(pos, 0, 1, seed, query, sub);
                    get_set_sw_to_pos(pos, 0, 2, seed, query, sub);
                }
            } else{
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
                //get_set_sw_to_pos(pos, 0, 1, seed, query, sub);
                get_set_sw_to_pos(pos, 0, 2, seed, query, sub);
            }
        }

    }
    return pos_j;
}

int classify_pos(idx_t *fm_idx, query_t *query, seed_t *seed, uint32_t pos_buf[], int pos_num, int seq_off[2], struct SubBuf *sub){
    int sn = query->seed_num;
    int sub_err = sub->sub_err;
    int l_seq = query->l_seq;
    int thres_pos_num = sub->thres_pos_num;

    if(query->b0 < query->candi_thres + 20) {
        thres_pos_num = IS_SMLSIZ; 
    } else if(query->b0 > l_seq - sn*sub_err){
        thres_pos_num = 0; 
    } else {
        thres_pos_num = IS_SMLSIZ*(l_seq - query->b0 - sn*sub_err)/l_seq; 
    }

    int thres_sw_to = sub->thres_sw_to;
    int cur_cls = seed->cls;
    int max_cls = query->max_ext;
    uint32_t pos;
    int pos_j = 0;
    int pi;
    for(pi = 0; pi < pos_num; ++pi) {
        pos = pos_buf[pi];
        if(pos > fm_idx->bwt->seq_len) {
            printf("%u, pos = %u\n", __LINE__, pos); 
            exit(1);
        }
        int flg_ed =  get_sw_ed_val(pos, sub, query);

fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
        if(flg_ed > 0) continue;
        if(pos_num <= thres_pos_num ) {
            pos_buf[pos_j++] = pos;
            set_sw_ed_pos(pos, query, sub);
            continue;
        }
        if(cur_cls == max_cls) {
            pos_buf[pos_j++] = pos;
            set_sw_ed_pos(pos, query, sub);
            continue;
        } 

fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
        int cl_num = 0; 
        cl_num = get_set_sw_to_pos(pos, 12*2, 0, seed, query, sub);
        if(cl_num > thres_sw_to) {
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d, thres_sw_to = %d, cl_num = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls, thres_sw_to, cl_num);
            pos_buf[pos_j++] = pos;
            set_sw_ed_pos(pos, query, sub);
        } else if(cl_num > 0){
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
            int sc_num = 0;
            int olp_flg = sub->olp_flg;
            //if(olp_flg > 0) 
                sc_num = eval_seed_olp(fm_idx,  query, seed, pos, seq_off,  sub);    
            int thres_sw_olp = sub->thres_sw_olp; 
            if(sc_num >= thres_sw_olp || sc_num + cur_cls >= max_cls -1) {
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
                pos_buf[pos_j++] = pos;
                set_sw_ed_pos(pos, query, sub);
            } else{
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
                //get_set_sw_to_pos(pos, 0, 1, seed, query, sub);// set_sw_to_pos
                get_set_sw_to_pos(pos, 0, 2, seed, query, sub);// set_sw_to_pos
            }
        } else{ //cl_num == 0
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
            int n_sm;
            if(query->is_rev == 0) n_sm = sub->n_sm_f;
            else n_sm = sub->n_sm_r;
            //if(n_sm < 3 || cur_cls > 2) {
            if(n_sm < 4 || cur_cls > 2) {
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
                int sc_num = 0;
                int thres_sw_olp = sub->thres_sw_olp; 
                int olp_flg = sub->olp_flg;
                //if(olp_flg > 0) 
                sc_num = eval_seed_olp(fm_idx,  query, seed, pos, seq_off,  sub);    
fprintf(stderr, "%u, sc_num = %d, thres_sw_olp = %d, max_cls = %d, cur_cls = %d\n", __LINE__, sc_num, thres_sw_olp, max_cls , cur_cls);
                if(sc_num >= thres_sw_olp + 1|| sc_num + cur_cls > max_cls + 1) {
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d, pos_j = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls, pos_j);
                    pos_buf[pos_j++] = pos;
                    set_sw_ed_pos(pos, query, sub);
                } else{
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d, seed_id = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls, seed->id);
                    //get_set_sw_to_pos(pos, 0, 1, seed, query, sub);
                    get_set_sw_to_pos(pos, 0, 2, seed, query, sub);
                }
            } else{
fprintf(stderr, "%u, pos = %u, flg = %d, pos_num = %d, cur_cls = %d\n", __LINE__, pos, flg_ed, pos_num, cur_cls);
                //get_set_sw_to_pos(pos, 0, 1, seed, query, sub);
                get_set_sw_to_pos(pos, 0, 2, seed, query, sub);
            }
        }

    }
    return pos_j;
}

void clean_sw_pos_arry(struct SubBuf *sub)
{
    int i, j, h, k, num, n_bit, r, r1, r2;
    uint32_t pos, pos0, pos1, pos2;
     
    sub->sw_ed_rlt_max[0] = 1; 
    sub->sw_ed_rlt_max[1] = 1; 
    sub->sw_to_rlt_max[0] = 1; 
    sub->sw_to_rlt_max[1] = 1; 
    sub->sw_ed_bit_max[0] = 1;
    sub->sw_ed_bit_max[1] = 1;
    sub->sw_to_bit_max[0] = 1;
    sub->sw_to_bit_max[1] = 1;

    int to_dat_size = sub->to_dat_size;
    int dat_blk_len = sub->dat_blk_len;
    int rlt_blk_len = sub->rlt_blk_len;

    int hsh_shift = sub->hsh_shift;
    int rlt_shift = sub->rlt_shift;
    int bit_shift = sub->bit_shift; 
    for(i = 0; i < 2; ++i) {
        uint32_t *sw_ed_hsh = sub->sw_ed_hsh[i];
        uint32_t *sw_ed_rlt = sub->sw_ed_rlt[i];
        int64_t *sw_ed_bit = sub->sw_ed_bit[i];
        uint32_t *sw_ed_pos = sub->pos_ed_buf[i]; 

        num = sw_ed_pos[0];
        for(j = 1; j <= num; ++j) {
            pos = sw_ed_pos[j]; 
            pos0 = pos >>hsh_shift;
            pos1 = (pos >>SW_DAT_BIT) & (((uint32_t)1<<SW_RLT_BIT) -1);
            pos2 = pos&(((uint32_t)1 <<SW_DAT_BIT) - 1); 
            uint32_t cur_hsh = sw_ed_hsh[pos0];
            int rlt_row = cur_hsh >>rlt_shift;
            int bit_row;
            for(h = rlt_row; h < rlt_blk_len + rlt_row; ++h) {// 4个1组
                if(sw_ed_rlt[h] >0 ){
                    //bit_row = sw_ed_rlt[h] >> 10; 
                    bit_row = sw_ed_rlt[h] >> bit_shift; 
                    for(k = 0; k < dat_blk_len/64; ++k) {// 4 = 256/64
                        sw_ed_bit[bit_row + k] = 0;
                    }
                    sw_ed_rlt[h] = 0;  
                }  
            }
            sw_ed_hsh[pos0] = 0;
        } 
        sw_ed_pos[0] = 0;
    }
    for(i = 0; i < 2; ++i) {
        uint32_t *sw_to_hsh = sub->sw_to_hsh[i];
        uint32_t *sw_to_rlt = sub->sw_to_rlt[i];
        int64_t *sw_to_bit = sub->sw_to_bit[i];
        uint32_t *sw_to_pos = sub->pos_to_buf[i]; 
        num = sw_to_pos[0];
        for(j = 1; j <= num; ++j) {
            pos = sw_to_pos[j]; 
            pos0 = pos >>hsh_shift;
            //pos1 = (pos >>8) &3;
            pos1 = (pos >>SW_DAT_BIT) & (((uint32_t)1<<SW_RLT_BIT) -1);
            //pos2 = pos&0xFF; 
            pos2 = pos&(((uint32_t)1 <<SW_DAT_BIT) - 1); 
            
            uint32_t cur_hsh = sw_to_hsh[pos0];
            //int rlt_row = cur_hsh >>8;
            int rlt_row = cur_hsh >>rlt_shift;
            int bit_row;
            for(h = rlt_row; h < rlt_blk_len + rlt_row; ++h) {// 4个1组
                if(sw_to_rlt[h] >0 ){
                    bit_row = sw_to_rlt[h] >> bit_shift; 
                    for(k = 0; k < to_dat_size*dat_blk_len/64; ++k) {// 4*4 = 256*4/64
                        sw_to_bit[bit_row + k] = 0;
                    }
                    sw_to_rlt[h] = 0;  
                }  
            }
            sw_to_hsh[pos0] = 0;
        } 
        sw_to_pos[0] = 0;
    }
    //sub->pos_bk_buf[0][0][0] = 0;
    //sub->pos_bk_buf[1][0][0] = 0;
    for(i = 0; i < LEN_READ/SEED_LEN*4; ++i) {
        sub->seed_bk[i] = 0; 
    }
    sub->n_seed_bk[0] = 0;
    sub->n_seed_bk[1] = 0;

    return;
} 
int get_sw_ed_val(uint32_t pos, struct SubBuf *sub, query_t *query)
{
    int bit_row;
    uint32_t pos0, pos1, pos2, pos_r;

    int is_rev = query->is_rev; 
    uint32_t *sw_ed_hsh = sub->sw_ed_hsh[is_rev];
    uint32_t *sw_ed_rlt = sub->sw_ed_rlt[is_rev];
    int64_t *sw_ed_bit = sub->sw_ed_bit[is_rev];
    uint32_t *sw_ed_pos = sub->pos_ed_buf[is_rev]; 
    int to_dat_size = sub->to_dat_size;
    int dat_blk_len = sub->dat_blk_len;
    int rlt_blk_len = sub->rlt_blk_len;

    int hsh_shift = sub->hsh_shift;
    int rlt_shift = sub->rlt_shift;
    int bit_shift = sub->bit_shift; 

    /*
    pos0 = pos >>10;
    pos1 = (pos >>8) &3;
    pos2 = pos&0xFF;
    */
    pos0 = pos >>hsh_shift;
    pos1 = (pos >>SW_DAT_BIT) & (((uint32_t)1<<SW_RLT_BIT) -1);
    pos2 = pos&(((uint32_t)1 <<SW_DAT_BIT) - 1);

    pos_r = pos2%64; 
    uint32_t cur_hsh = sw_ed_hsh[pos0];
    if(cur_hsh == 0) { return 0;}
    //uint32_t cur_rlt = sw_ed_rlt[(cur_hsh>>8) + pos1];   
    uint32_t cur_rlt = sw_ed_rlt[(cur_hsh>>rlt_shift) + pos1];   
    if(cur_rlt == 0) { return 0;}
    //bit_row = (cur_rlt >> 10)+pos2/64;                
    bit_row = (cur_rlt >> bit_shift)+pos2/64;                
    int flg = (sw_ed_bit[bit_row] >>(63-pos_r))&1;
    
fprintf(stderr, "%u, pos = %u, bit_row = %u, cur_hsh = %u, cur_rlt = %u, sw_ed_bit= %u\n", __LINE__, pos, bit_row, cur_hsh, cur_rlt, sw_ed_bit[bit_row]); 
    
    return flg;
}

int delete_sw_ed_pos(uint32_t *pos_buf, int pos_num, struct SubBuf *sub, query_t * query)
{
    int i, j, num, n_bit, r, r1, r2, row;
    uint32_t pos, pos0, pos1, pos2, pos_r;

    int is_rev = query->is_rev; 
    uint32_t *sw_ed_hsh = sub->sw_ed_hsh[is_rev];
    uint32_t *sw_ed_rlt = sub->sw_ed_rlt[is_rev];
    int64_t *sw_ed_bit = sub->sw_ed_bit[is_rev];
    uint32_t *sw_ed_pos = sub->pos_ed_buf[is_rev]; 

    int *sw_ed_rlt_max = sub->sw_ed_rlt_max;
    int *sw_ed_bit_max = sub->sw_ed_bit_max;
    int to_dat_size = sub->to_dat_size;
    int dat_blk_len = sub->dat_blk_len;
    int rlt_blk_len = sub->rlt_blk_len;

    int hsh_shift = sub->hsh_shift;
    int rlt_shift = sub->rlt_shift;
    int bit_shift = sub->bit_shift; 


    num = pos_num;
    uint32_t pos_i = 0;
    for(j = 0; j < num; ++j) {
        pos = pos_buf[j]; 
        /*  
        pos0 = pos >>10;
        pos1 = (pos >>8) &3;
        pos2 = pos&0xFF;
        */
        pos0 = pos >>hsh_shift;
        pos1 = (pos >>SW_DAT_BIT) & (((uint32_t)1<<SW_RLT_BIT) -1);
        pos2 = pos&(((uint32_t)1 <<SW_DAT_BIT) - 1);
 
        pos_r = pos2%64; 
        uint32_t cur_hsh = sw_ed_hsh[pos0];
        if(cur_hsh == 0) {
            continue;
        }
        //uint32_t cur_rlt = sw_ed_rlt[(cur_hsh>>8) + pos1];   
        uint32_t cur_rlt = sw_ed_rlt[(cur_hsh>>rlt_shift) + pos1];   
        if(cur_rlt == 0) {//
            continue;
        }
        
        //row = (cur_rlt >> 10)+pos1+pos2/64;                
        row = (cur_rlt >> bit_shift) + pos2/64;                
        int flg = (sw_ed_bit[row] >>(63-pos_r))&1;
        if(flg > 0) {
            pos_buf[pos_i++] = pos; 
        } else {} 
       
      
    } 


    return pos_i;
}
int set_sw_ed_buf(uint32_t *pos_buf, int pos_num, query_t * query, struct SubBuf *sub)
{
    int i, j, n_bit, r, r1, r2, row;
    uint32_t pos, pos0, pos1, pos2, pos_r;

    int is_rev = query->is_rev; 
    uint32_t *sw_ed_hsh = sub->sw_ed_hsh[is_rev];
    uint32_t *sw_ed_rlt = sub->sw_ed_rlt[is_rev];
    int64_t *sw_ed_bit = sub->sw_ed_bit[is_rev];
    uint32_t *sw_ed_pos = sub->pos_ed_buf[is_rev]; 

    int *sw_ed_rlt_max = sub->sw_ed_rlt_max;
    int *sw_ed_bit_max = sub->sw_ed_bit_max;
    
    int to_dat_size = sub->to_dat_size;
    int dat_blk_len = sub->dat_blk_len;
    int rlt_blk_len = sub->rlt_blk_len;

    int hsh_shift = sub->hsh_shift;
    int rlt_shift = sub->rlt_shift;
    int bit_shift = sub->bit_shift; 

    int pos_i = 0; 
    for(j = 0; j < pos_num; ++j) {
        pos = pos_buf[j]; 
        pos0 = pos >>hsh_shift;
        pos1 = (pos >>SW_DAT_BIT) & (((uint32_t)1<<SW_RLT_BIT) -1);
        pos2 = pos&(((uint32_t)1 <<SW_DAT_BIT) - 1);
        pos_r = pos2%64; 
        
        uint32_t cur_hsh = sw_ed_hsh[pos0];
        int bit_row, rlt_row;
        if(sw_ed_rlt_max[is_rev] > MAX_CANDI_POS * rlt_blk_len) {
            printf("%u, MAX_CANDI_POS too small!!!\n", __LINE__); 

            pos_buf[pos_i] = pos; 
            pos_i++;
            continue;
        }
        if(cur_hsh == 0) {
            if(sw_ed_rlt_max[is_rev] > MAX_CANDI_POS * rlt_blk_len) {
                printf("%u, error!!!\n", __LINE__); 
                exit(1);
                //return;
            }
            cur_hsh = sw_ed_rlt_max[is_rev] << rlt_shift; 
            sw_ed_hsh[pos0] = cur_hsh;  
            sw_ed_rlt_max[is_rev] += rlt_blk_len; 
            sw_ed_pos[++sw_ed_pos[0]] = pos; 
        }
        rlt_row = (cur_hsh>>rlt_shift) + pos1;
        int cur_rlt = sw_ed_rlt[rlt_row];   
        if(sw_ed_bit_max[is_rev] > MAX_CANDI_POS * (dat_blk_len/64)) {
            printf("%u, MAX_CANDI_POS too small!!!\n", __LINE__); 
            pos_buf[pos_i] = pos; 
            pos_i++;
            continue;
        }
        if(cur_rlt == 0) {//
            if(sw_ed_bit_max[is_rev] > MAX_CANDI_POS * (dat_blk_len/64)) {
                printf("%u, error!!!\n", __LINE__);
                exit(1);
            }
            cur_rlt = sw_ed_bit_max[is_rev] << bit_shift; 
            sw_ed_rlt[rlt_row] = cur_rlt;  
            sw_ed_bit_max[is_rev] += (dat_blk_len/64);//一个区间是256bits, 64位方式访问 
        }
        if(cur_hsh == 0 ||  cur_rlt == 0) {
            printf("%u, error!!!\n", __LINE__);
            exit(1);
        }
        bit_row = (cur_rlt >>bit_shift)+pos2/64;                
        int flg = (sw_ed_bit[bit_row] >>(63-pos_r) ) & 1;
        if(flg > 0) {
            continue;
        } 
        sw_ed_bit[bit_row] |= (uint64_t)1<<(63-pos_r);
        pos_buf[pos_i] = pos; 
        pos_i++;
        int n_rlt = sw_ed_rlt[rlt_row] &(((uint32_t)1<<bit_shift)-1);
        if(n_rlt < dat_blk_len) {
            ++sw_ed_rlt[rlt_row]; 
        }   
    } 
    return pos_i;
}

int set_sw_ed_pos(uint32_t pos, query_t * query, struct SubBuf *sub)
{
    uint32_t pos0, pos1, pos2, pos_r;
    int is_rev = query->is_rev; 
    uint32_t *sw_ed_hsh = sub->sw_ed_hsh[is_rev];
    uint32_t *sw_ed_rlt = sub->sw_ed_rlt[is_rev];
    int64_t *sw_ed_bit = sub->sw_ed_bit[is_rev];
    uint32_t *sw_ed_pos = sub->pos_ed_buf[is_rev]; 
    int *sw_ed_rlt_max = sub->sw_ed_rlt_max;
    int *sw_ed_bit_max = sub->sw_ed_bit_max;
    
    int to_dat_size = sub->to_dat_size;
    int dat_blk_len = sub->dat_blk_len;
    int rlt_blk_len = sub->rlt_blk_len;

    int hsh_shift = sub->hsh_shift;
    int rlt_shift = sub->rlt_shift;
    int bit_shift = sub->bit_shift; 

    
    /*  
    pos0 = pos >>10;
    pos1 = (pos >>8) &3;
    pos2 = pos&0xFF;
    */
    pos0 = pos >>hsh_shift;
    pos1 = (pos >>SW_DAT_BIT) & (((uint32_t)1<<SW_RLT_BIT) -1);
    pos2 = pos&(((uint32_t)1 <<SW_DAT_BIT) - 1);

    pos_r = pos2%64; 
    uint32_t cur_hsh = sw_ed_hsh[pos0];
    uint32_t bit_row, rlt_row;
    if(sw_ed_rlt_max[is_rev] > MAX_CANDI_POS * rlt_blk_len) {

        printf("%u, MAX_CANDI_POS is too small!!!!, name = %s, sw_ed_rlt_max[rev] = %d\n", __LINE__, query->name, sw_ed_rlt_max[is_rev]); 
//exit(1);        
        return;
    }    
    if(cur_hsh == 0) {
        //if(sw_ed_rlt_max[is_rev] > MAX_CANDI_POS * 4*4) {
        if(sw_ed_rlt_max[is_rev] > MAX_CANDI_POS * rlt_blk_len) {
            printf("%u, error!!!\n", __LINE__); 
            exit(1);
            //return;
        }

        //cur_hsh = sw_ed_rlt_max[is_rev] << 8; 
        cur_hsh = sw_ed_rlt_max[is_rev] << rlt_shift; 
        sw_ed_hsh[pos0] = cur_hsh;  
        //sw_ed_rlt_max[is_rev] += 4*4;//1个节点4个32位组成 
//sw_ed_rlt_max[is_rev] += rlt_blk_len*sizeof(uint32_t);//1个节点4个32位组成 
        sw_ed_rlt_max[is_rev] += rlt_blk_len;//1个节点4个32位组成 
        sw_ed_pos[++sw_ed_pos[0]] = pos; 
    }
    if(sw_ed_bit_max[is_rev] > MAX_CANDI_POS * (dat_blk_len/64)) {
        printf("%u, MAX_CANDI_POS too small!!!\n", __LINE__); 
        return;
    }


    //rlt_row = (cur_hsh>>8) + pos1;
    rlt_row = (cur_hsh>>rlt_shift) + pos1;
    if(rlt_row > sub->sw_ed_rlt_size/sizeof(uint32_t)){
        fprintf(stderr, "%u, cur_hsh = %u, rlt_shift = %u, pos1 = %u\n", __LINE__, cur_hsh, rlt_shift, pos1);
        fprintf(stderr, "%u, rlt_row(%u) out of size(%u).\n", __LINE__, rlt_row, sub->sw_ed_rlt_size);
        exit(1);
    }
    
    
    uint32_t cur_rlt = sw_ed_rlt[rlt_row];   
 
    if(cur_rlt == 0) {//
        //if(sw_ed_bit_max[is_rev] > MAX_CANDI_POS * (256/64)) {
        if(sw_ed_bit_max[is_rev] > MAX_CANDI_POS * (dat_blk_len/64)){
            printf("%u, error!!!\n", __LINE__);
            exit(1);
            //return;
        }

        //cur_rlt = sw_ed_bit_max[is_rev] << 10; 
        cur_rlt = sw_ed_bit_max[is_rev] << bit_shift; 

        sw_ed_rlt[rlt_row] = cur_rlt;  
        //sw_ed_bit_max[is_rev] += 256/64;//一个区间是256bits, 64位方式访问 
        sw_ed_bit_max[is_rev] += dat_blk_len/64;//一个区间是256bits, 64位方式访问 
    }

    if(cur_hsh == 0 ||  cur_rlt == 0) {
        printf("%u, error!!!\n", __LINE__);
        exit(1);
    }
    //bit_row = (cur_rlt >>10)+pos2/64;                
    bit_row = (cur_rlt >>bit_shift)+pos2/64;                
    sw_ed_bit[bit_row] |= (uint64_t)1<<(63-pos_r);

    //int n_rlt = sw_ed_rlt[rlt_row] &0x3FF;
    int n_rlt = sw_ed_rlt[rlt_row] &(((uint32_t) 1 << bit_shift) -1);
    //if(n_rlt < 256) {
    if(n_rlt <  dat_blk_len) {
        ++sw_ed_rlt[rlt_row]; 
    }   
    return;
}
//func: get and set sw to arry
int get_set_sw_to_pos(uint32_t pos, int fnd_len, int set_flg, seed_t *seed, query_t *query, struct SubBuf *sub)
{
    int i, j, num, n_bit, r, r1, r2, row;
    uint32_t pos0, pos1, pos2, pos_r, pos_i;
   
    int is_rev = query->is_rev; 
    uint32_t *sw_to_hsh = sub->sw_to_hsh[is_rev];
    uint32_t *sw_to_rlt = sub->sw_to_rlt[is_rev];
    int64_t *sw_to_bit = sub->sw_to_bit[is_rev];
    uint32_t *sw_to_pos = sub->pos_to_buf[is_rev]; 
    uint32_t (*sw_pos_bak)[2] = sub->pos_bk_buf[is_rev]; 
    uint32_t *seed_bk = sub->seed_bk; 
    int *n_seed_bk = sub->n_seed_bk+is_rev;
    int cls = seed->cls, seed_id = seed->id;
fprintf(stderr, "%u, seed_id = %d\n", __LINE__, seed_id);
    int *sw_to_rlt_max = sub->sw_to_rlt_max;
    int *sw_to_bit_max = sub->sw_to_bit_max;
        
    int to_dat_size = sub->to_dat_size;
    int dat_blk_len = sub->dat_blk_len;
    int rlt_blk_len = sub->rlt_blk_len;

    int hsh_shift = sub->hsh_shift;
    int rlt_shift = sub->rlt_shift;
    int bit_shift = sub->bit_shift; 

    int olp_num = 0;
    if(fnd_len > 12*2 || fnd_len < 0) {
        printf("%u, error!!!\n", __LINE__);
        exit(1);
    }  
    int bg = -fnd_len/2, ed = fnd_len/2;
    if(pos < ed) ed = pos;   
    for(i = bg; i < ed; ++i) {

        pos_i = pos + i; 
        /*
        pos0 = pos_i >>10;
        pos1 = (pos_i >>8) &3;
        pos2 = pos_i&0xFF;
        */
        pos0 = pos_i >>hsh_shift;
        pos1 = (pos_i >>SW_DAT_BIT) & (((uint32_t)1<<SW_RLT_BIT) -1);
        pos2 = pos_i&(((uint32_t)1 <<SW_DAT_BIT) - 1);
 
        pos_r = pos2%64; 
        int cur_hsh = sw_to_hsh[pos0];

        uint32_t rlt_row, bit_row;
        if(cur_hsh == 0) { continue; }
        rlt_row = (cur_hsh>>rlt_shift) + pos1;
        int cur_rlt = sw_to_rlt[rlt_row];   
        if(cur_rlt == 0) { continue; }
        bit_row = (cur_rlt >>bit_shift) + (pos2*to_dat_size)/64;  
        uint64_t cur_bit = sw_to_bit[bit_row];
        cur_bit = (cur_bit >>(63- pos2*to_dat_size%64))&(((uint32_t)1 << to_dat_size) -1);  
        olp_num += cur_bit;    
    }//end for(i = bg; i < ed; ++i)++++++
fprintf(stderr, "%u\n", __LINE__);
    if(sw_to_rlt_max[is_rev] > MAX_CANDI_POS * rlt_blk_len) {
        printf("%u, MAX_CANDI_POS too small!!!\n", __LINE__); 
        return olp_num;
    }

    if(set_flg > 0){
        pos0 = pos >>hsh_shift;
        pos1 = (pos >>SW_DAT_BIT) & (((uint32_t)1<<SW_RLT_BIT) -1);
        pos2 = pos&(((uint32_t)1 <<SW_DAT_BIT) - 1);
 
        pos_r = pos2%64; 
        uint32_t cur_hsh = sw_to_hsh[pos0];
        uint32_t rlt_row, bit_row;
        if(cur_hsh == 0) {
            cur_hsh = sw_to_rlt_max[is_rev] << rlt_shift; 
            sw_to_hsh[pos0] = cur_hsh;  
            sw_to_rlt_max[is_rev] += rlt_blk_len;//1个节点4个32位组成 
            sw_to_pos[++sw_to_pos[0]] = pos;   
        }
        rlt_row = (cur_hsh>>rlt_shift) + pos1;
        uint32_t cur_rlt = sw_to_rlt[rlt_row];   
        if(sw_to_bit_max[is_rev] > MAX_CANDI_POS * (dat_blk_len/64) * to_dat_size) {
            printf("%u, MAX_CANDI_POS too small!!!\n", __LINE__); 
            return olp_num;
        }
        if(cur_rlt == 0) {//
            cur_rlt = sw_to_bit_max[is_rev] << bit_shift; 
            sw_to_rlt[rlt_row] = cur_rlt;  
            sw_to_bit_max[is_rev] += (dat_blk_len/64) * to_dat_size;
        }
        bit_row = (cur_rlt >>bit_shift) + (pos2*to_dat_size)/64;  
        uint64_t cur_bit = sw_to_bit[bit_row];
        cur_bit = (cur_bit >>(63- pos2*to_dat_size%64))& (((uint32_t)1<<to_dat_size) -1);  
        if(cur_bit < (((uint32_t)1<< to_dat_size) -1)) {
            int n_shift = 63-pos2*to_dat_size%64;
            sw_to_bit[bit_row] += ((uint64_t)1<<n_shift);
        }
        if(set_flg > 1) {
            if(cur_bit == 0 && cls == 0) {

fprintf(stderr, "%u, n_seed_bk[0] = %d, pos = %u, seed_id = %d\n", __LINE__, n_seed_bk[0], pos, seed_id);
                sw_pos_bak[n_seed_bk[0]][0] = pos;
                sw_pos_bak[n_seed_bk[0]][1] = seed_id;
                sub->seed_bk[seed_id]++;
                n_seed_bk[0]++;  
            }   
        }
        
    } 
    return olp_num;
}





#define HASH_KMER_NUM 20 
#define HASH_KMER_INTV 30 
#define HASH_SCORE_THRES 100
int clean_hash_pos(idx_t *fm_idx, uint32_t pos_buf[], int pos_num, query_t *query, struct SubBuf *sub){
    uint32_t w = 1;
    uint32_t pos, pos_r, pos_q, pos_wq, pos_wr;
    int score, i, j;

    uint16_t *hash_idx = query->hash_idx[query->is_rev];
    uint16_t *hash_pos = query->hash_pos[query->is_rev];
    int l_seq = query->l_seq;
    uint32_t seq10 = 0, pos_j = 0, pi;
    uint32_t ch, num, h_pos, l_pos = 0, r_pos = l_seq - 1, min_len = 0, match_num = 0;
   
    int stat_flg = 0, hash_sc = 0;
    pos_j = 0;
    for(pi = 0; pi < pos_num; ++pi) {
        pos = pos_buf[pi];
        seq10 = 0;
        l_pos = 0, r_pos = l_seq - 1, min_len = 0, match_num = 0;
        for(i = 0; i < 9; ++i) {
            ch = __get_pac(fm_idx->pac, pos+i);
            seq10 <<= 2;
            //seq10 |= (uint32_t)read_seq[i]&3; 
            seq10 |= ch;  
        }
  
        for(i = 9; i < l_seq; ++i) {
            ch = __get_pac(fm_idx->pac, pos+i);
            seq10 <<= 2;
            seq10 |= ch;
            seq10 &= 0xFFFFF;//20 bits 
    //fprintf(stderr, "%u, j = %d, seq10 = %x\n", __LINE__, j, seq10);
            min_len = l_seq;
            if(hash_idx[seq10] > 0) {
                h_pos = hash_idx[seq10];
                num = hash_pos[h_pos];
                for(j = h_pos+1; j <= h_pos+num; ++j) {
                    if(hash_pos[j] <= i && hash_pos[j] > l_pos) l_pos = hash_pos[j]; 
                    if(hash_pos[j] > i && hash_pos[j] < r_pos) r_pos = hash_pos[j]; 
                }
                if(i - l_pos <= r_pos - i) {
                    min_len = i - l_pos; 
                } else {
                    min_len = r_pos - i; 
                } 
                if(min_len < HASH_KMER_INTV) {
                    match_num++; 
                } 
                if(stat_flg == 0) {
                    hash_sc += 10; //seq10 长度 
                } else {
                    hash_sc += 1; 
                }
                stat_flg = 1; 
            } else {
                h_pos = 0;
                num = 0; 
                stat_flg = 0;
            } 
        
        }
        
        if(hash_sc > HASH_SCORE_THRES) {
            pos_buf[pos_j++] = pos; 
        }  
  

    }
    return pos_j;
}


/*  
int clean_swed_pos(uint32_t pos_buf[], int pos_num, query_t *query, struct SubBuf *sub){

return pos_num;
    uint64_t *sw_ed  = sub->sw_ed[query->is_rev]; 
    uint64_t *sw_to  = sub->sw_to[query->is_rev]; 
    uint32_t *pos_ed = sub->pos_ed[query->is_rev]; 
    uint32_t *pos_to = sub->pos_to[query->is_rev]; 
    uint32_t w = 1;
    uint32_t pos, pos_r, pos_q, pos_wq, pos_wr;
    int score, i, j;
    uint64_t flg_to, flg_ed;

    j = 0;
    for(i = 0; i < pos_num; ++i) {
        pos = pos_buf[i];

        pos_q = pos / 64; 
        pos_r = pos % 64;
        flg_ed = sw_ed[pos_q] & ((uint64_t)1<<(63-pos_r)); 
fprintf(stderr, "%u, pos = %u, flg_ed = %d, flg = %d\n", __LINE__, pos, flg_ed, flg_ed >>(63-pos_r));
        if(flg_ed == 0) {
            pos_buf[j++] = pos;
            if(sw_ed[pos_q] == 0 ) pos_ed[++pos_ed[0]] = pos_q;
            sw_ed[pos_q] |= ((uint64_t)1<<(63-pos_r));
        } //end if(flg_ed == 0)+++++++++++                    
    }
    return j;
}
*/
int pos_filter(uint64_t *sw_ed, uint64_t *sw_to, uint32_t *pos_ed, uint32_t *pos_to, uint32_t pos_buf[], int pos_num, int cls){
    uint32_t pos, pos_r, pos_q;
    int flg_to, flg_ed, score, i, j;
    if(cls == 0) score = 1;
    else if(cls ==1 || cls == 2) score = 2;
    else if(cls == 3 || cls == 4) score = 3;
    else score = 4;   
//fprintf(stderr, "%u, pos_num = %u\n", __LINE__, pos_num);
    j = 0;
    for(i = 0; i < pos_num; ++i) {
        pos = pos_buf[i];
        pos_q = pos/64;
        pos_r = pos%64;
        flg_ed = sw_ed[pos_q] & 1<<(63-pos_r); 
        if(flg_ed == 0) {
            pos_q = pos/64;  
            pos_r = pos%64;
            flg_to = (sw_to[pos_q] >> (62-pos_r/2*2))&3;
            flg_to += score;
            if(flg_to >= 4) { 
                if(sw_ed[pos_q] == 0) pos_ed[++pos_ed[0]] = pos;  
                sw_ed[pos_q] |= 1<<(63-pos_r);
                pos_buf[j++] = pos;
            } else if(flg_to -score == 0) {
//fprintf(stderr, "%u, %u, pos_to[0] = %u\n",__LINE__, i,  pos_to[0]);
                if(sw_to[pos_q] == 0) pos_to[++pos_to[0]] = pos;  
                sw_to[pos_q] |= flg_to <<(62-pos_r/2*2);

//fprintf(stderr, "%u, %u, pos_to[0] = %u\n",__LINE__, i,  pos_to[0]);
            } else {
                sw_to[pos_q] &= ~(3 <<(62-pos_r/2*2));
                sw_to[pos_q] |= flg_to <<(62-pos_r/2*2);
            }                         
        } //end if(flg_ed == 0)+++++++++++                    
    }
    return j;
}
int opt_sw_to_buf(uint8_t *sw_ed, uint8_t *sw_to, uint32_t *pos_ed, uint32_t *pos_to, uint32_t pos_buf[], int opt){
    uint32_t pos, pos_r, pos_q;
    int flg_to, flg_ed, score, i, j, size = sw_to[0];
    uint32_t *buf[5]; 

    for(i = 1; i < 5; ++i) {
        buf[i] = sw_to + i*size;
        buf[i][0] = 0;
    }
    for(i = 0; i < size; ++i) {
        pos = sw_to[i];
        pos_q = pos/8;
        pos_r = pos%8;
        flg_ed = sw_ed[pos_q] & 1<<(7-pos_r); 
        if(flg_ed == 0) {
            pos_q = pos/16;  //(pos－pos％2)/2/8
            pos_r = pos%16/2;//(pos－pos％2)/2%8
            flg_to = (sw_to[pos_q] >> (7-pos_r))&3;
            if(flg_to >= opt) { 
                buf[flg_to][++buf[flg_to][0]] = pos; 
            }                       
        } //end if(flg_ed == 0)+++++++++++                    
    }
    j = 0;
    for(i = opt; i < 5; ++i) {
        j += buf[i][0];
    }

    return j;
}

void gen_cnt2(const uint8_t *bwt, int n_data, uint8_t *cnt)
{
    int DataNum = n_data;
    int i, siz_cnt, rng_num, len_sum, rot, rot_off, lst_row, rot_sum, row, rng_idx, rng_off, ch_off;
    uint8_t ch, get_ch;
    if(DataNum<=254*256) siz_cnt = 2;
    if(DataNum>254*256) siz_cnt = 4;
    rng_num = ((DataNum+253) /254);
    len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;//cnt size per character per rot

    memset(cnt, 0, 16*8*len_sum);
    const uint8_t *pBwt = bwt;
    uint8_t *pSum  = cnt ;
    uint32_t glb_sum[8][16];
    for(rot=0;rot<8; rot++){ for(ch = 0; ch < 16 ; ch++)    glb_sum[rot][ch] = 0;}
    for(rot = 0; rot < 8 ; rot++ ){
        rot_off = rot*((DataNum+1)/2);
        lst_row = rot_off + DataNum/2  ;
        rot_sum = rot*16*len_sum  ;
        for(row = rot_off ;  row < lst_row ; row++){
            rng_idx = (row -rot_off)/127 ;
            rng_off = (8+ siz_cnt) * (rng_idx/8) + rng_idx%8  ;  
    
            get_ch = pBwt[row] ;
            ch = get_ch>>4;
            ch_off = rot_sum +ch*len_sum ;
            pSum[ch_off + rng_off]++;
            ch = get_ch & 0xF;
            ch_off = rot_sum +ch*len_sum ;
            pSum[ch_off + rng_off]++;
        }
            
        if(DataNum%2>0){
            rng_idx = (lst_row -rot_off)/127 ;
            rng_off = (8+ siz_cnt) * (rng_idx/8) + rng_idx%8  ;  
            get_ch = pBwt[lst_row] ;
            ch = get_ch>>4;
            ch_off = rot_sum +ch*len_sum ;
            pSum[ch_off + rng_off]++;
        }
        
        for( ch = 0 ; ch <16 ; ch++){
            ch_off = rot_sum +ch*len_sum ;
            //+++++++++++++++++++++++++++++
            rng_idx = 0 ;
            rng_off = 0 ;
            int sum_buf = 0 ;
                
            while( rng_idx < rng_num ){
                rng_off = (8+ siz_cnt) * (rng_idx/8) + rng_idx%8  ;
                sum_buf +=  pSum[ch_off + rng_off];
                rng_idx++;
                if(rng_idx%8==0){
                    int buf = sum_buf ;
                    for(i=0; i<siz_cnt; i++){
                        ////fprintf(stderr, "siz_cnt0 %u-%u\t", buf, buf & 255);
                        pSum[ch_off + rng_off+(siz_cnt-i)] = (uint8_t)(buf & 255) ;

                        buf >>= 8;
                      
                    }
                }
               
            } // End : +++++++++++++++++++while( rng_idx <= rng_num )
            glb_sum[rot][ch] = (ch==0?0:glb_sum[rot][ch-1])+sum_buf;
        } // End : +++++++++++++++++++++  for( ch = 0 ; ch <16 ; ch++)
//if(rng_num%8== 0) continue;//need this line???
        for( ch = 0 ; ch <16 ; ch++){
            ch_off = rot_sum +ch*len_sum +(len_sum -1) ;
            int buf = glb_sum[rot][ch] ;
            for(i=0; i<siz_cnt; i++){
                pSum[ch_off -i] = (uint8_t)(buf & 255) ;
                ////fprintf(stderr, "siz_cnt0 %u\t", buf & 255);
                buf >>= 8;                      
            }
        }
    }//end for(rot=0; rot <8; ++rot)
}
int PairExtSeq(struct ExtBlck *eB, struct SubBuf *sub);
uint32_t get_bwt_size(uint32_t n_data)
{
    uint32_t bwt_size = (n_data+1)/2*8; 
    int siz_cnt = n_data>254*256?4:2;
    int rng_num = (n_data+253)/254;
    int len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;
    uint32_t cnt_size = 16*8*len_sum;
    if(n_data <= MIN_BWT_SIZE) {
        return n_data*sizeof(uint32_t); 
    }
    //++++++++++++++++++++++++++++++++++++++ 
    if(n_data <= NO_BWT_SUM_SIZE) cnt_size = 0;
    else if(n_data <= 255) cnt_size = 17*8;
    return cnt_size + bwt_size; 
}

void destroyIdxsArry(struct JmpMod *jmp, struct StackTree *sTree,  struct ExtBlck eBlck[], struct SubBuf *sub)
{
    int i;
    //free(jmp->cap);
printf("%u\n", __LINE__);
    free(jmp->jmp);
printf("%u\n", __LINE__);
    for(i = 0; i < NUM_EXT; ++i) free(eBlck[i].head_cap);
printf("%u\n", __LINE__);
    free(sTree->stck_arry);
printf("%u\n", __LINE__);
	free(sTree->exted_arry);
printf("%u\n", __LINE__);
	free(sTree->smpos_arry);
printf("%u\n", __LINE__);
	free(sTree->back_arry);
printf("%u\n", __LINE__);
    free(sTree->call);
printf("%u\n", __LINE__);
    free(sub->seqRel_out);		
printf("%u\n", __LINE__);
	free(sub->seqRel_L);		
printf("%u\n", __LINE__);
	free(sub->seqRel_R);		
printf("%u\n", __LINE__);
	free(sub->pair_out->pair_arry);
printf("%u\n", __LINE__);
	free(sub->pair_out);
printf("%u\n", __LINE__);
	free(sub->algn_out);
printf("%u\n", __LINE__);
    free(sub->algnR_row);
printf("%u\n", __LINE__);
	free(sub->algnL_row);
printf("%u\n", __LINE__);
	free(sub->algnRel_row);
printf("%u\n", __LINE__);
    free(sub->algnRel_buf);
printf("%u\n", __LINE__);
    free(sub->relat_buf);
printf("%u\n", __LINE__);


printf("%u\n", __LINE__);
    free(sub->aln_out->out_buf);
printf("%u\n", __LINE__);
    free(sub->aln_out);
printf("%u\n", __LINE__);
    free(sub->err_buf);
printf("%u\n", __LINE__);
    free(sub->pos_buf);


    free(sub->kswq_f[0]);
printf("%u\n", __LINE__);
    free(sub->kswq_f[1]);
printf("%u\n", __LINE__);
    free(sub->kswq_r[0]);
printf("%u\n", __LINE__);
    free(sub->kswq_r[1]);
printf("%u\n", __LINE__);
    free(sub->kswst.qp);
printf("%u\n", __LINE__);
    free(sub->kswst.eh);
printf("%u\n", __LINE__);
    free(sub->kswgb.z);
printf("%u\n", __LINE__);
    free(sub->kswgb.qp);
printf("%u\n", __LINE__);
    free(sub->kswgb.eh);
printf("%u\n", __LINE__);
    free(sub->kswgb.cigar);
printf("%u, tgt_arry = %x\n", __LINE__, sub->tgt_arry);    
    free(sub->tgt_arry);
printf("%u\n", __LINE__);
    /*  
    free(sub->sw_ed[0]);
    free(sub->sw_ed[1]);
    free(sub->sw_to[0]);
    free(sub->sw_to[1]);
    free(sub->pos_ed[0]);
    free(sub->pos_ed[1]);
    free(sub->pos_to[0]);
    free(sub->pos_to[1]);
    free(sub->pos_bk[0]);
    free(sub->pos_bk[1]);
    */
}
struct SubBuf *InitIdxsArry(struct FileName *fn, struct JmpMod **out_jmp, struct CapIfo *cap, struct ExtBlck eBlck[], struct StackTree *sTree, struct SubBuf *sub)
{
    char in_f[128];
    FILE *fp;
	if((fp=fopen(fn->jmpmod,"rb"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}
    uint32_t buf_head[100], head_size; 
    //fread(buf_head,sizeof(uint32_t),4,fp);  
    fread(buf_head,sizeof(uint32_t),3,fp);  
    
    //int max_buf_size = buf_head[3];
    int max_buf_size = buf_head[2];
	uint32_t data_size = 0;
	uint32_t *jmp_idx;           	
    //uint32_t *cap_pos;
	uint8_t  *mod_idx; 	
    	
	//data_size = buf_head[0] + buf_head[1]+ (buf_head[2]+3)/4;
	data_size = buf_head[0] + (buf_head[1]+3)/4;
    uint32_t *jmp_mod;
    if(NULL == (jmp_mod = (uint32_t*)malloc(data_size*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    /*   
    fread(cap_pos, sizeof(uint32_t), buf_head[0], fp);
    jmp_idx = cap_pos+buf_head[0];
    */
    jmp_idx = jmp_mod;     
    fread(jmp_idx, sizeof(uint32_t), buf_head[0], fp);
    //mod_idx = (uint8_t *)(cap_pos+buf_head[0]+buf_head[1]); 
    mod_idx = (uint8_t *)(jmp_mod+buf_head[0]); 
    fread(mod_idx, sizeof(uint8_t), buf_head[1], fp);
    
    printf("[log]: read %d jmp_idx and %d mod_idx, max_buf_size = %d\n", buf_head[0], buf_head[1], max_buf_size); 
    
    fclose(fp);
	struct JmpMod *jmp;
    if(NULL == (jmp= (struct JmpMod *)malloc(sizeof(struct JmpMod)))){
	    perror("error...");
	    exit(1);
	}

    jmp->jmp = jmp_idx;
	jmp->mod = mod_idx;
	//jmp->cap = cap_pos;
	*out_jmp = jmp;
    //初始化，jmp_idx,mod_idx,cap_pos数组结束------------
    int i, j;
    uint32_t *buf_data ; 
    //以下代码段读取comfile， 
    size_t file_size[NUM_EXT];
    uint32_t *head_idx[NUM_EXT];
   
    uint32_t f_id ; 
    for(f_id = 0; f_id < NUM_EXT; ++f_id){
        file_size[f_id] = getFileSize(fn->comfile[f_id]); 
        head_idx[f_id] = (uint32_t *)malloc(file_size[f_id]); 
        FILE *fp = fopen(fn->comfile[f_id], "rb");
        
        fread(buf_head,4, NUM_FILES+1, fp); 
        fread(head_idx[f_id],1, file_size[f_id]-4*(NUM_FILES+1), fp); 
        fclose(fp);

        size_t tmp, sum=0;
        buf_head[NUM_FILES] = 0; 
        for(j = 0; j < NUM_FILES; ++j){
            tmp = sum;
            sum += buf_head[j];
            buf_head[j] = tmp;  
        } 
        buf_head[NUM_FILES] = sum;  
       
        
        memset(eBlck+f_id, 0, sizeof(struct ExtBlck)); 
        eBlck[f_id].head_cap    =  (struct CapIfo *)(head_idx[f_id]+buf_head[0]);
		eBlck[f_id].head_nxtpnt =  (uint32_t  *)(head_idx[f_id]+buf_head[1]);		   
        eBlck[f_id].head_nxtflg =  (uint8_t *)(head_idx[f_id]+buf_head[2]);
		eBlck[f_id].head_relat  =  (uint32_t *)(head_idx[f_id]+buf_head[3]);
		eBlck[f_id].head_smbwt  =  (uint8_t  *)(head_idx[f_id]+buf_head[4]);
        eBlck[f_id].head_extidx =  (uint32_t (*)[2])(head_idx[f_id]+buf_head[5]);
        eBlck[f_id].n_cap = (buf_head[1] - buf_head[0])*sizeof(uint32_t)/sizeof(struct CapIfo); 
        eBlck[f_id].n_nxtpnt = (buf_head[2] - buf_head[1])*sizeof(uint32_t)/sizeof(struct CapIfo); 
        eBlck[f_id].n_nxtflg = (buf_head[3] - buf_head[2])*sizeof(uint32_t)/sizeof(struct CapIfo); 
        eBlck[f_id].n_relat = (buf_head[4] - buf_head[3])*sizeof(uint32_t)/sizeof(struct CapIfo); 
        eBlck[f_id].n_smbwt = (buf_head[5] - buf_head[4])*sizeof(uint32_t)/sizeof(struct CapIfo); 
        eBlck[f_id].n_extidx = (file_size[f_id]-4*(NUM_FILES+1) - 4*buf_head[5])/2/sizeof(uint32_t); 

        
        printf("n_cap[%d] = %d\n", f_id, eBlck[f_id].n_cap);
        printf("n_nxtpnt[%d] = %d\n", f_id, eBlck[f_id].n_nxtpnt);
        printf("n_nxtflg[%d] = %d\n", f_id, eBlck[f_id].n_nxtflg);
        printf("n_relat[%d] = %d\n", f_id, eBlck[f_id].n_relat);
        printf("n_smbwt[%d] = %d\n", f_id, eBlck[f_id].n_smbwt);
        printf("n_extidx[%d] = %d\n", f_id, eBlck[f_id].n_extidx);
    } 
    uint32_t buf32[2];
    uint8_t buf8[2];
    struct CapIfo tmp_cap[2];

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++

    stck_arry_t *sArry = malloc( max_buf_size*sizeof(stck_arry_t)); 
	if(NULL == sArry){
	    perror("error...");
	    exit(1);
	}
    stck_arry_t *exted = malloc( max_buf_size*sizeof(stck_arry_t)); 
	if(NULL == exted){
	    perror("error...");
	    exit(1);
	}
    stck_arry_t *smpos = malloc( max_buf_size*sizeof(stck_arry_t)); 
	if(NULL == smpos){
	    perror("error...");
	    exit(1);
	}
    stck_arry_t *back_buf = malloc( max_buf_size*sizeof(stck_arry_t)); 
	if(NULL == back_buf){
	    perror("error...");
	    exit(1);
	}
    sTree->stck_arry = sArry;
	sTree->exted_arry = exted;
	sTree->smpos_arry = smpos;
	sTree->back_arry = back_buf;
    initStackTree(sTree);

    //*out_stree = sTree;
	/*  
    struct SubBuf *sub ;
    if(NULL == (sub = (struct SubBuf*)malloc(sizeof(struct SubBuf)))){
	    perror("error...");
	    exit(1);
	}
    */

    int seq_out_size = 360;
    uint32_t *seqRel_out;
	if(NULL == (seqRel_out = (uint32_t*)calloc(seq_out_size*seq_out_size, sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    uint32_t *seqRel_L;
	if(NULL == (seqRel_L = (uint32_t*)calloc(seq_out_size*seq_out_size, sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    uint32_t *seqRel_R;
	if(NULL == (seqRel_R = (uint32_t*)calloc(seq_out_size*seq_out_size, sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    pair_out_t *pair_out = malloc(sizeof(pair_out_t));
    pair_arry_t *pair_arry = malloc(max_buf_size*sizeof(pair_arry_t));
    if(pair_arry == NULL) {
    	perror("error: pair_arry...");
	    exit(1);
    }
    pair_out->pair_arry = pair_arry;
    sub->pair_out = pair_out;
    uint32_t *pos_buf = malloc(max_buf_size*sizeof(uint32_t));
    if(pos_buf == NULL) {
    	perror("error: pos_buf...");
	    exit(1);
    }
    uint32_t (*err_buf)[5] = malloc(5*max_buf_size*sizeof(uint32_t));
    if(err_buf == NULL) {
    	perror("error:err_buf...");
	    exit(1);
    }
    //sub->max_pos_buf = max_buf_size;
    sub->pos_buf = pos_buf;
    sub->err_buf = err_buf;

	uint32_t (*algn_out)[5];
	if(NULL == (algn_out = (uint32_t(*)[5])malloc(max_buf_size*5*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    uint8_t *algnR_row;
	if(NULL == (algnR_row = (uint8_t*)calloc(max_buf_size/8, sizeof(uint8_t)))){
	    perror("error...");
	    exit(1);
	}
    uint8_t *algnL_row;
	if(NULL == (algnL_row = (uint8_t*)calloc(max_buf_size/8, sizeof(uint8_t)))){
	    perror("error...");
	    exit(1);
	}
    uint8_t *algnRel_row;
	if(NULL == (algnRel_row = (uint8_t*)calloc(max_buf_size/8, sizeof(uint8_t)))){
	    perror("error...");
	    exit(1);
	}
    uint8_t *algnRel_buf;
	if(NULL == (algnRel_buf = (uint8_t*)calloc(max_buf_size/8, sizeof(uint8_t)))){
	    perror("error...");
	    exit(1);
	}
    uint32_t *relat_buf;
	if(NULL == (relat_buf = (uint32_t*)calloc(max_buf_size, sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    
    sub->seqRel_out   = seqRel_out;		
	sub->seqRel_L     = seqRel_L;		
	sub->seqRel_R     = seqRel_R;		
	sub->pair_out     = pair_out;
	sub->algn_out     = algn_out;
    sub->algnR_row    = algnR_row;
	sub->algnL_row    = algnL_row;
	sub->algnRel_row  = algnRel_row;
    sub->algnRel_buf  = algnRel_buf;
    sub->relat_buf    = relat_buf;
    sub->max_buf_size = max_buf_size; 
    uint32_t slen = (LEN_READ+15)/8; 
    sub->kswq_f[0] = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (5 + 4)); // a single block of memory
    sub->kswq_f[1] = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (5 + 4)); // a single block of memory
    sub->kswq_r[0] = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (5 + 4)); // a single block of memory
    sub->kswq_r[1] = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (5 + 4)); // a single block of memory
    for(i = 0; i < 5; ++i){
        for(j = 0; j < 5; ++j){
            if(i == j) sub->mat[i*5+j] = 1;
            else if (i == 4 || j == 4) sub->mat[i*5+j] = -1;
            else sub->mat[i*5+j] = -4;
        }
    }
    sub->kswst.qp = malloc(LEN_READ * 5);
	sub->kswst.eh = calloc(LEN_READ + 1, 8);
	

    sub->kswgb.z = malloc(LEN_READ * (LEN_READ+50));
	sub->kswgb.qp = malloc(LEN_READ * 5);
	sub->kswgb.eh = calloc(LEN_READ + 1, 8);
    sub->kswgb.cigar = calloc(LEN_READ, 4);	
    sub->kswgb.m_cigar = LEN_READ; 
    sub->tgt_arry = calloc(1, sizeof(tgt_arry_t));
printf("%u, tgt_arry = %x\n", __LINE__, sub->tgt_arry);    
    return sub;
}
void clean_sub_buf( struct ExtBlck *cB, struct SubBuf *sub)
{
    int i, j, r, r_R;
    if(sub->seqL_out[0][0] >= 360 ) {
        printf("%u, sub->seqL_out[0][0] =%u\n", __LINE__, sub->seqL_out[0][0]); 
        exit(1);
    }
    if(sub->seqR_out[0][0] >= 360 ) {
        printf("%u, sub->seqR_out[0][0] =%u\n", __LINE__, sub->seqR_out[0][0]); 
        exit(1);
    }
    /* 
    if(cB->num_relat <= 0xFF){
        uint8_t *R2rel = (uint8_t *)(cB->head_relat+cB->relat) 
                        + cB->num_relat+cB->num_seqL+1;
        for(i =1; i <=sub->seqR_out[0][0]; ++i){
            r_R = sub->seqR_out[i][0];
            for(j = R2rel[r_R]; j < R2rel[r_R+1]; j++){
                sub->algnRel_buf[j/8] = 0;
            }
        }
    } else if(cB->num_relat <= 0xFFFF){
        uint16_t *R2rel = (uint16_t *)(cB->head_relat+cB->relat) 
                          + cB->num_relat+cB->num_seqL+1;
        for(i =1; i <=sub->seqR_out[0][0]; ++i){
            r_R = sub->seqR_out[i][0];
            for(j = R2rel[r_R]; j < R2rel[r_R+1]; j++)
                sub->algnRel_buf[j/8] = 0;
        }
    } else{
        uint32_t *R2rel = cB->head_relat + cB->relat 
                          + cB->num_relat+cB->num_seqL+1;
        for(i =1; i <=sub->seqR_out[0][0]; ++i){
            r_R = sub->seqR_out[i][0];
            for(j = R2rel[r_R]; j < R2rel[r_R+1]; j++)
                sub->algnRel_buf[j/8] = 0;
        }
    }
    */
/*  
    for(i=1; i<=sub->seqL_out[0][0]; i++){
        r = sub->seqL_out[i][0];
        sub->algnL_row[r/8] = 0; 
    }
    for(i=1; i<=sub->seqR_out[0][0]; i++){
        r = sub->seqR_out[i][0];
        sub->algnR_row[r/8] = 0; 
    }
*/
    for(i=1; i<=sub->seqRel_out[0]; i++){
        r = sub->seqRel_out[i];
        sub->algnRel_row[r/8] = 0; 
    }
    /*  
    for(i=1; i<=sub->seqRel_R[0]; i++){
        r = sub->seqRel_R[i];
        sub->algnRel_buf[r/8] = 0; 
    }
    for(i=1; i<=sub->seqRel_L[0]; i++){
        r = sub->seqRel_L[i];
        sub->algnRel_buf[r/8] = 0; 
    }
    */
    sub->seqRel_out[0] = 0;
    sub->seqRel_L[0] = 0;
    sub->seqRel_R[0] = 0;
    sub->seqL_aln_old = 1;
    sub->seqR_aln_old = 1;
    sub->seqRel_aln_old = 1;
    sub->seqRel_L_old = 1;
    sub->seqRel_R_old = 1;
    //++++++++++++++++++++++++++++++++++++++++++++++++++
    for(i = 0; i < 20; ++i){
        sub->aln_L[i] = 0;
        sub->aln_R[i] = 0;
    }
    return;
}
void setBlckData_uni(struct JmpMod *jmp, struct ExtBlck *in_blck, struct ExtBlck **out_blck, query_t * query, seed_t *seed, struct StackTree *sTree,  int uni_flg, struct SubBuf *sub){
    struct ExtBlck *cB = *out_blck;      
    int i, r;
    //清零左段和右段smbwt数据 
    //int cls = query->max_ext;
fprintf(stderr, "%u\n", __LINE__);
    sub->seqL_out = sub->call_buf + sub->seqL_off[NUM_EXT]; 
    sub->seqR_out = sub->call_buf + sub->seqR_off[NUM_EXT]; 
fprintf(stderr, "%u\n", __LINE__);
    //----------------------- 
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    int l_seq = query->l_seq; 

fprintf(stderr, "%u, sTree->len_arry = %d\n", __LINE__, sTree->len_arry);
    if(sTree->len_arry < 1){
        int seed_id = seed->id;
        uint32_t bgn = seed->slc[seed_id].bgn;
        uint32_t end = seed->slc[seed_id].end;
        uint32_t num = seed->slc[seed_id].num;

        uint32_t buf_algn[4] = {};             
        buf_algn[0] = SEED_LEN;
        buf_algn[1] = bgn;
        buf_algn[2] = num;

fprintf(stderr, "%u, num = %u\n", __LINE__, num);
        if(num <= IS_SMLSIZ) {
            return; 
        } 
        getCapPos(jmp,buf_algn);
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fprintf(stderr, "%u, bgn = %u, num = %u, buf_aln[3] = %d, extidx = %u\n", __LINE__, bgn, num, buf_algn[3], in_blck->head_extidx[buf_algn[3]][0]);
        sTree->cls = 0;
        sTree->stck_arry[1].nxtpnt = buf_algn[3];
        sTree->stck_arry[1].nxtflg = num;
        sTree->stck_arry[1].err = 0;
        sTree->stck_arry[1].l_off = 0; //64
        sTree->stck_arry[1].r_off = 0; //64
        sTree->stck_arry[1].cls   = 0;
        sTree->len_arry = 1; 
        sTree->len_old = 0; 
fprintf(stderr, "%u, buf_algn[3] = %d\n", __LINE__, buf_algn[3]);
        /*  
        sub->seqL_out = sub->call_buf + sub->seqL_off[0]; 
        sub->seqR_out = sub->call_buf + sub->seqR_off[0]; 
        */
fprintf(stderr, "seqL_off[0] = %d, seqR_off[0] = %d\n", sub->seqL_off[0], sub->seqR_off[0]);
        sub->seqL_out[0][0] = 0;
        sub->seqR_out[0][0] = 0;
       
        sTree->len_exted = 0; 
        sTree->len_smpos = 0; 
      
        sub->sub_f[0] = 0; 
        sub->sub_f[1] = 1; //aln_r = 1
        sub->sub_f[2] = 1; //aln_r = 2
        sub->sub_f[3] = 1; //aln_r = 3
        sub->sub_f[4] = 0; //aln_r = 4
        sub->sub_f[5] = 0; //aln_r = 5
        sub->sub_f[6] = 0; //aln_r = 6
        sub->sub_f[7] = 0; //aln_r = 7
        sub->sub_f[8] = 0; //aln_r = 8
        sub->sub_f[9] = 0; //aln_r = 9
        sub->indel_f[0] = 0;  
        sub->indel_f[1] = 0; 
        sub->indel_f[2] = 0; 
        sub->indel_f[3] = 0; 
        sub->indel_f[4] = 0; 
        sub->indel_f[5] = 1; 
        sub->indel_f[6] = 1; 
        sub->indel_f[7] = 0; 
        sub->indel_f[8] = 0; 
        sub->indel_f[9] = 0; 

    }
fprintf(stderr, "%u, sTree->len_arry = %d\n", __LINE__, sTree->len_arry);
    sTree->cls = sTree->stck_arry[sTree->len_arry].cls;
    sTree->len_old = sTree->stck_arry[sTree->len_arry].len_p;

    uint32_t cap_row = sTree->stck_arry[sTree->len_arry].nxtpnt;
    cB = in_blck + sTree->cls;
    if(sTree->cls == 0){
        int seed_id = seed->id;
        uint32_t bgn = seed->slc[seed_id].bgn;
        uint32_t end = seed->slc[seed_id].end;
        uint32_t num = seed->slc[seed_id].num;
        cB->bg_idx = bgn;
        cB->ed_idx = end;  
    }
        
    
    struct CapIfo  *cur_cap;
fprintf(stderr, "%u\n", __LINE__);
    cur_cap = cB->head_cap + cap_row;
    cB->relat     = cur_cap->relat;	
    cB->nxtpnt    = cur_cap->nxtpnt;
    cB->smbwt     = cur_cap->smbwt;
    cB->num_seqL  = cur_cap->num_seqL;
    cB->num_seqR  = cur_cap->num_seqR;
    cB->num_relat = cur_cap->num_relat;
    cB->extidx    = cap_row;
    cB->bwtL  =  cB->smbwt ;
	uint32_t len = get_bwt_size(cB->num_seqL);

    cB->bwtR  =  cB->smbwt + len;	
    *out_blck = cB;
    sub->err_sum[0] = sTree->stck_arry[sTree->len_arry].err;
    int seed_off = seed->slc[seed->id].s_off;
    int l_off = sTree->stck_arry[sTree->len_arry].l_off;
    int r_off = sTree->stck_arry[sTree->len_arry].r_off;
    int L_offset = seed_off-sTree->cls*16-16       - l_off;
    int R_offset = SEED_LEN+seed_off+sTree->cls*16 + r_off;
    sub->l_off = l_off;
    sub->r_off = r_off;
    if(l_off > 1 || l_off < -1 || r_off > 1 || r_off < -1) {
        printf("l_off = %d, r_off =%d\n", l_off, r_off); 
        //exit(1);
    }  
    //int cut_len = 16 + sub->l_seq_delta;
    sub->l_seq_delta = 0;
    int cut_delta = sub->l_seq_delta;
    int cut_len = 16 + cut_delta;
  
    if(uni_flg == 0) { //处理左段序列N
        if(L_offset < cut_delta) {
            printf("%u, L_offset = %d\n", __LINE__, L_offset);
            exit(1); 
        }
        memcpy(sub->ext_seqL, read_seq+L_offset-sub->l_seq_delta, cut_len);
fprintf(stderr, "%u, ext_seqL, L_offset = %d\n", __LINE__, L_offset);
        for(i = 0; i < cut_len; ++i) {
            if(sub->ext_seqL[i] == 4) {
                sub->ext_seqL[i] = 0;
            }     
fprintf(stderr, "%u ", sub->ext_seqL[i]);
        }
fprintf(stderr, "\n");
        /*  
        //左序列反转 
        for(i = 0; i < 9; ++i){
            uint8_t ch_buf = sub->ext_seqL[i]; 
            sub->ext_seqL[i] = sub->ext_seqL[17-i];
            sub->ext_seqL[17-i] = ch_buf;
        }
        */
    } else if(uni_flg > 0) {   //处理右段序列N
        if(l_seq - R_offset < 16) {
            printf("%u, R_offset = %d, seed_off = %d, l_off = %d, r_off = %d\n", __LINE__, R_offset, seed_off, l_off, r_off);
            exit(1); 
        }
fprintf(stderr, "%u, R_offset = %d\n", __LINE__, R_offset);
        memcpy(sub->ext_seqR, read_seq+R_offset, cut_len);   
       
fprintf(stderr, "%u\n", __LINE__);
        for(i = 0; i < cut_len; ++i) {
            if(sub->ext_seqR[i] == 4) {
                sub->ext_seqR[i] = 0;
            }     
        }
fprintf(stderr, "%u\n", __LINE__);
    } else {
        printf("[error]: %u, uni_flg = %d\n", __LINE__); 
        exit(1);
    }
    //--------------------------	
    
fprintf(stderr, "%u\n", __LINE__);
    return;
}

void setBlckData(struct ExtBlck *in_blck, struct SubBuf *sub, struct StackTree *sTree, uint8_t *read_seq, int seed_off, struct ExtBlck **out_this){
    int i;
    struct ExtBlck *cB = *out_this;   
    if(sTree->len_arry < 1){
        printf("%u, len_arry == %d\n", __LINE__, sTree->len_arry);
        exit(1);
    }
    sTree->cls = sTree->stck_arry[sTree->len_arry].cls;
    sTree->len_old = sTree->stck_arry[sTree->len_arry].len_p;


    uint32_t cap_row = sTree->stck_arry[sTree->len_arry].nxtpnt;
    sub->err_sum[0] = sTree->stck_arry[sTree->len_arry].err;
    sub->seqL_out = sub->call_buf + sub->seqL_off[sTree->cls]; 
    sub->seqR_out = sub->call_buf + sub->seqR_off[sTree->cls]; 
    int l_off = sTree->stck_arry[sTree->len_arry].l_off;
    int r_off = sTree->stck_arry[sTree->len_arry].r_off;
    int L_offset = seed_off-sTree->cls*16-16       - l_off;
    int R_offset = SEED_LEN+seed_off+sTree->cls*16 + r_off;
    sub->l_off = l_off;
    sub->r_off = r_off;
    
    if(l_off > 1 || l_off < -1 || r_off > 1 || r_off < -1) {
        printf("l_off = %d, r_off =%d\n", l_off, r_off); 
        //exit(1);
    }  
//fprintf(stderr, "%u\n", __LINE__);
//fprintf(stderr, "L_offset = %u, l_off = %u\n", L_offset, l_off);
//fprintf(stderr, "R_offset = %u, r_off = %u\n", R_offset, r_off);
    sub->l_seq_delta = 2;
    int cut_delta = sub->l_seq_delta;
    int cut_len = 16 + cut_delta;
    memcpy(sub->ext_seqL, read_seq+L_offset-cut_delta, cut_len);
    memcpy(sub->ext_seqR, read_seq+R_offset, cut_len);

    
    //处理左段序列N
    //err->seq_L[0] = 0;
    for(i = 0; i < cut_len; ++i) {
        if(sub->ext_seqL[i] == 4) {
            //err->seq_L[++err->seq_L[0]];
            sub->ext_seqL[i] = 0;
        }     
    }
    //处理右段序列N
    //err->seq_R[0] = 0;
    for(i = 0; i < cut_len; ++i) {
        if(sub->ext_seqR[i] == 4) {
            //err->seq_R[++err->seq_R[0]];
            sub->ext_seqR[i] = 0;
        }     
    }
    /*  
    //左序列反转 
    for(i = 0; i < 9; ++i){
        uint8_t ch_buf = sub->ext_seqL[i]; 
        sub->ext_seqL[i] = sub->ext_seqL[17-i];
        sub->ext_seqL[17-i] = ch_buf;
    }
    */
    //--------------------------	
    
    cB = in_blck + sTree->cls;
    struct CapIfo  *cur_cap;
    cur_cap = cB->head_cap + cap_row;
    cB->relat     = cur_cap->relat;	
    cB->nxtpnt    = cur_cap->nxtpnt;
    cB->smbwt     = cur_cap->smbwt;
    cB->num_seqL  = cur_cap->num_seqL;
    cB->num_seqR  = cur_cap->num_seqR;
    cB->num_relat = cur_cap->num_relat;
    cB->extidx    = cap_row;
    cB->bwtL  =  cB->smbwt ;
	uint32_t len = get_bwt_size(cB->num_seqL);
    cB->bwtR  =  cB->smbwt + len;	
    *out_this = cB;
    
    return;
}

void initStackTree(struct StackTree *sTree){
	sTree->cls = 0;
    sTree->stck_arry[0].nxtpnt = 0; //其中buf_align[3]是capidx数组的行号;
	sTree->stck_arry[0].nxtflg = 0; 
	sTree->stck_arry[0].err = 0;
    sTree->len_arry = 0; 
	sTree->len_buf = 0;
    return ;
}
void set_smpos(struct StackTree *sTree, struct SubBuf *sub, seed_t *seed, int seed_id, struct ExtBlck*eB){
	uint32_t (*ext_idx)[2];
    ext_idx = (eB+sTree->cls+1)->head_extidx;
    uint32_t i= 0, j = 0;
    
    pair_out_t *pair_out = sub->pair_out;
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    //struct call_t *call = sTree->call;
    i = sTree->len_smpos;
    for (j = pair_out->s_num; j < pair_out->p_num; ++j, ++i){
        sTree->smpos_arry[i].nxtpnt =  pair_arry[j].nxtpnt;
        sTree->smpos_arry[i].nxtflg =  pair_arry[j].nxtflg;
        sTree->smpos_arry[i].err    =  pair_arry[j].err;
        sTree->smpos_arry[i].l_off  =  pair_arry[j].l_off;
        sTree->smpos_arry[i].r_off  =  pair_arry[j].r_off;
        sTree->smpos_arry[i].cls    =  sTree->cls+1;
        //sTree->smpos_arry[i].idx_bg  = ext_idx[pair_arry[j].nxtpnt][0];
        //sTree->smpos_arry[i].idx_num = ext_idx[pair_arry[j].nxtpnt][1];
        sTree->smpos_arry[i].len_p = sTree->len_arry;
        sTree->smpos_arry[i].len_g = sTree->len_old; 
        sTree->smpos_arry[i].aln_f = 0; 
        sTree->smpos_arry[i].aln_r = 0;
    } // end for
    if(sTree->cls + 1 == sTree->max_cls) {
        for(i = 1; i <= seed->flg[seed_id][0]; ++i){
            j = seed->flg[seed_id][i];//j是删除的seed 
            seed->del[j] = 1; 
        }
    }
    //sub->pair_out->p_num = 0;
    //sub->pair_out->s_num = 0;
 
////fprintf(stderr, "%u, cls = %u, len_arry = %u, sTree->len_exted = %u\n", 
//                  __LINE__,sTree->cls,  sTree->len_arry, sTree->len_exted);
    return;
}

void setStackTree(struct StackTree *sTree, struct SubBuf *sub, seed_t *seed, int seed_id, struct ExtBlck*eB){
	uint32_t (*ext_idx)[2];
    ext_idx = (eB+sTree->cls+1)->head_extidx;
    uint32_t i= 0, j = 0;
    
    pair_out_t *pair_out = sub->pair_out;
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    struct call_t *call = sTree->call;
    i = sTree->len_arry;
    if( call[sTree->cls].aln_f > 1) call[sTree->cls].aln_f = 5;//???
    sTree->stck_arry[i].aln_f = call[sTree->cls].aln_f; 
    sTree->stck_arry[i].aln_r = call[sTree->cls].aln_r;
    if(sTree->cls + 1 < sTree->max_cls) {
        i = sTree->len_arry + 1;
        for (j = pair_out->s_num; j < pair_out->p_num; ++j, ++i){
            sTree->stck_arry[i].nxtpnt =  pair_arry[j].nxtpnt;
            sTree->stck_arry[i].nxtflg =  pair_arry[j].nxtflg;
            sTree->stck_arry[i].err    =  pair_arry[j].err;
            sTree->stck_arry[i].l_off  =  pair_arry[j].l_off;
            sTree->stck_arry[i].r_off  =  pair_arry[j].r_off;
            sTree->stck_arry[i].cls    =  sTree->cls+1;
            //sTree->stck_arry[i].idx_bg  = ext_idx[pair_arry[j].nxtpnt][0];
            //sTree->stck_arry[i].idx_num = ext_idx[pair_arry[j].nxtpnt][1];
            sTree->stck_arry[i].len_p = sTree->len_arry;
            sTree->stck_arry[i].len_g = sTree->len_old; 
            sTree->stck_arry[i].aln_f = 0; 
            sTree->stck_arry[i].aln_r = 0; 
        } // end for
        sTree->len_arry = i - 1;     
    } else{ 
        i = sTree->len_exted;  
        for (j = pair_out->s_num; j < pair_out->p_num; ++j, ++i){
            sTree->exted_arry[i].nxtpnt =  pair_arry[j].nxtpnt;
            sTree->exted_arry[i].nxtflg =  pair_arry[j].nxtflg;
            sTree->exted_arry[i].err    =  pair_arry[j].err;
            sTree->exted_arry[i].l_off  =  pair_arry[j].l_off;
            sTree->exted_arry[i].r_off  =  pair_arry[j].r_off;
            sTree->exted_arry[i].cls    =  sTree->cls+1;
            //sTree->exted_arry[i].idx_bg  = ext_idx[pair_arry[j].nxtpnt][0];
            //sTree->exted_arry[i].idx_num = ext_idx[pair_arry[j].nxtpnt][1];
            sTree->exted_arry[i].len_p = sTree->len_arry; 
            sTree->exted_arry[i].len_g = sTree->len_old; 
            sTree->exted_arry[i].aln_f = 0; 
            sTree->exted_arry[i].aln_r = 0;
             
        }

        sTree->len_exted = i;
        //以下代码删除不必要的种子
        /*  
        for(i = 1; i <= seed->flg[seed_id][0]; ++i){
            j = seed->flg[seed_id][i];//j是删除的seed 
            seed->del[j] = 1; 
        }
        */
    }
    sub->pair_out->p_num = 0;
    sub->pair_out->s_num = 0;
 
fprintf(stderr, "%u, cls = %u, len_arry = %u, sTree->len_exted = %u\n", 
                  __LINE__,sTree->cls,  sTree->len_arry, sTree->len_exted);
    return;
}
    /*     
            sTree->back_buf[k][0] =  pair_out[j][0];
            sTree->back_buf[k][1] =  pair_out[j][1];
            sTree->back_buf[k][2] =  pair_out[j][2];
            sTree->back_buf[k][3] =  pair_out[j][3];
            sTree->back_buf[k][4] =  pair_out[j][4];
            sTree->back_buf[k][5] =  sTree->cls+1;
            sTree->back_buf[k][6] =  seed_id;

            if(pair_out[j][1] == 0xFF){
                sTree->back_buf[k][7] =  ext_idx[pair_out[j][1]][0];
                sTree->back_buf[k][8] =  ext_idx[pair_out[j][1]][1];
            }
    */

/*     
    int i, j, bg, ed, num, max_error, max_row;
    for (j = pair_out[0][1]+1; j <= pair_out[0][0]; ++j){
        //fprintf(stderr, "pair_out[%u] = %u  ", j, pair_out[j][0]);       
        //fprintf(stderr, "%u  ", pair_out[j][1]);       
        //fprintf(stderr, "%u  ", pair_out[j][2]);       
        //fprintf(stderr, "%u  ", pair_out[j][3]);       
        //fprintf(stderr, "%u  ", pair_out[j][4]);       
        //fprintf(stderr, "cls = %u\n", sTree->cls+1);
    }
    j = sTree->len_arry;
    bg = pair_out[0][1]+1;
    ed = pair_out[0][0]+1;

    while (bg < ed) {
        max_row = bg;
        max_error = pair_out[bg][0];
        for(i = bg+1; i < ed; ++i){
            if(pair_out[i][2] > max_error) { 
                max_error = pair_out[i][2];
                max_row = i;
            } 
        }
        sTree->stck_arry[j][0] =  pair_out[max_row][0];
        sTree->stck_arry[j][1] =  pair_out[max_row][1];
        sTree->stck_arry[j][2] =  pair_out[max_row][2];
        sTree->stck_arry[j][3] =  pair_out[max_row][3];
        sTree->stck_arry[j][4] =  pair_out[max_row][4];
        sTree->stck_arry[j][5] =  sTree->cls+1;
        ++j;
        if(max_row == bg) {
            bg++;
            continue; 
        }
        if(bg+1 < ed){
            num = 5*sizeof(uint32_t) *(pair_out[0][0]-max_row);
            memcpy(pair_out+max_row-1, pair_out+max_row, num);
        }
        ed--;
    }
    for (i = sTree->len_arry; i < j; ++i){
        //fprintf(stderr, "stck_arry[%u] = %u  ", i, sTree->stck_arry[i][0]);       
        //fprintf(stderr, "%u  ", sTree->stck_arry[i][1]);       
        //fprintf(stderr, "%u  ", sTree->stck_arry[i][2]);       
        //fprintf(stderr, "%u  ", sTree->stck_arry[i][3]);       
        //fprintf(stderr, "%u  ", sTree->stck_arry[i][4]);       
        //fprintf(stderr, "%u\n", sTree->stck_arry[i][5]);
    }
*/
int sortTreeData(struct StackTree *sTree)
{
   	int bg = sTree->len_old;
    int ed = sTree->len_arry;
    int row, max, i;
    while( bg < ed){
        max = -1;
        for(i = bg; i <= ed; ++i){
            if(sTree->stck_arry[i].err > max) {
                max = sTree->stck_arry[i].err;
                row = i; 
            } 
        }
        if(row > bg){
            memcpy(sTree->stck_arry+ed+1, sTree->stck_arry+bg,   6*sizeof(uint32_t));
    		memcpy(sTree->stck_arry+bg,   sTree->stck_arry+row,  6*sizeof(uint32_t));
		    memcpy(sTree->stck_arry+row,  sTree->stck_arry+ed+1, 6*sizeof(uint32_t));   
        }
        bg++; 
    }
    return ed-bg;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void getSeed(uint8_t *seed_seq, uint8_t *read_seq, uint8_t  seed_off){
	//在输入read_seq中获取，偏离read_seq中心seed_off(正数右偏离，负数左偏离)距离的
	//种子输出到seed_seq；

	


    int i;
    for(i=0; i < 20; ++i) seed_seq[i] = read_seq[i];
	return;
}

void AlgnSeed(uint8_t *seed_seq,uint32_t buf_algn[],struct PubParm *pParm){
	//输入的种子序列seed_seq比对到Bwt得到的结果保存到buf_algn[]中
	//其中buf_algn的[0]是输入，种子长度；即buf_algn[0] = LEN_SEED;
	//buf_algn[1]是返回值BgnIdx,buf_algn[2]是返回值NumIdx
		
	return;
}



		
void getCapPos(struct JmpMod *jmp,uint32_t buf_algn[]){
	// 其中buf_algn的[0]是输入种子长度；[1]是输入BgnIdx,[2]是输入NumIdx
	// [3]是返回值capidx或者smpos数组的行号，
	// jmp有三个分量 uint32_t *jmp，uint32_t *cap，uint8_t *mod;
	// 该函数在bgn_row=jmp[BgnIdx/256]和bgn_end=jmp[BgnIdx/256+1]范围内mod中查找，
	// 即mod[bgn_row]到mod[end_row]中查找与r_idx=BgnIdx%256的相等的行号mod_row;
	//uint32_t *cap = jmp->cap;
	uint32_t *jmp_idx = jmp->jmp;
    uint8_t *mod_idx = jmp->mod;

    uint32_t bgn = buf_algn[1];	
    uint32_t bgn_q = jmp_idx[bgn/256];
	uint32_t end_q = jmp_idx[(bgn/256)+1];
	uint32_t i_r = bgn%256, j_idx;
    for(j_idx=bgn_q; j_idx<end_q; j_idx++){
        if(i_r == mod_idx[j_idx]){
            break;
        } 
    }
    //if(j_idx == end_q){
        //fprintf(stderr, "%u, %s\n", __LINE__, __func__);
        //fprintf(stderr, "j_idx= %u, bgn_q = %u, end_q = %u\n", j_idx, bgn_q, end_q);        
        //fprintf(stderr,"bg_idx= %u, i_r = %u\n", bgn, i_r);
       // exit(1); 
   // }
    
    //buf_algn[3] = cap[j_idx];
    buf_algn[3] = j_idx;
	return;	
}



int Algn_auto_l2r(idx_t *fm_idx, uint32_t hash_boundry[], query_t *query, struct SubBuf *sub,uint32_t best_pos[])
{
    uint32_t i, k, l, num, num_r;
    uint32_t pos_buf[IS_SMLSIZ];
    uint32_t idx_buf[LEN_READ-SEED_LEN][2];
    uint32_t aln_buf[LEN_READ/SEED_LEN+2][5];
    uint32_t tmp[SEED_LEN-12][2];
    uint8_t *pseed;
    uint32_t seq12;

    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    int read_len = query->l_seq;
    uint32_t n_aln = 0, idx_num = 0, aln_len = 0;
    int flag = 0;
    int pos_c, pos_s = SEED_LEN;
    while(1){//重新选定新的种子位点
        pos_c = pos_s -12;
        pseed = read_seq + pos_c; 
        seq12 = lkt_seq2LktItem(pseed, 0, 11);
        k = fm_idx->fastmap->item[seq12];
        l = fm_idx->fastmap->item[seq12+1]-1;
        l -= get_12mer_correct(hash_boundry, l);
        n_aln = 12;
         
        while(2){//种子比对扩展
            if(n_aln <= SEED_LEN) { 
                tmp[n_aln-12][0] = k;
                tmp[n_aln-12][1] = l+1-k;
            }
            pseed = read_seq + pos_c; 
            num = bwt_match_exact_alt(fm_idx->bwt, 1, pseed, &k, &l);

            if(num > 0 ) {
                n_aln++; 
             
                if(num <= IS_SMLSIZ) {
//fprintf(stderr, "%u, num = %u\n", __LINE__, num);          
                    for(i =0; i < num; ++i) pos_buf[i] = bwt_sa(fm_idx->bwt, k+i)-pos_c; 
                    num_r = AlgnPos(fm_idx, query, num, sub);  
//fprintf(stderr, "%u, num_r = %u\n", __LINE__, num_r);          
                    if(num_r >0 && query->b0 >= query->best_thres) {
                        flag = 2;          
                        break;
                    }
                } else{
                 
                    if(aln_len >= SEED_LEN) {
//fprintf(stderr, "%u, num = %u\n", __LINE__, num);          
                        idx_buf[idx_num][0] = k;
                        idx_buf[idx_num][1] = num;
                        idx_num++;  
                    } 
                    
                }
            } else{//num = 0 的情况
                if(tmp[aln_len][1] <= 4*IS_SMLSIZ) {//上一次的比对结果, tmp[][1]是上一次的num
                    //用SW方法进行比对
                    k = tmp[aln_len][0];
                    num = tmp[aln_len][1];
//fprintf(stderr, "%u, num = %u\n", __LINE__, num);          
                    for(i =0; i < num; ++i) pos_buf[i] = bwt_sa(fm_idx->bwt, k+i)-pos_c; 
                    num_r = AlgnPos(fm_idx, query, num, sub);  

//fprintf(stderr, "%u, num_r = %u\n", __LINE__, num_r);          
                    if(num_r >0 && query->b0 >= query->best_thres) {
                        flag = 2;          
                        break;
                    }
                } 
                if(n_aln >= SEED_LEN) {
                    aln_buf[aln_len][0] = idx_num - n_aln + SEED_LEN;
                    aln_buf[aln_len][1] = n_aln - SEED_LEN;
                    aln_buf[aln_len][2] = pos_s;
                    aln_len++; 
                } 
                if(n_aln >= SEED_LEN) {
                    pos_s = pos_s + SEED_LEN; 
                } else{
                    pos_s = pos_s + SEED_LEN - n_aln; 
                }
                break; // while(2)
            }
            pos_c--;
        
            if(pos_s >= read_len || pos_c < 0){ 
                if(n_aln >= SEED_LEN) {
                    flag = 1; 
                } else{
                    flag = 0; 
                }
                break;
            }
            n_aln = 0;
        } // end while(2)
        if(flag == 2) { break; }
        if(flag == 1) {
            //如果新产生的种子序列尚未sw比对，那么进行多种子重叠比对 
        }   
    }// end while(1)
    return flag;    
}

int Algn_auto_r2l(idx_t *fm_idx, uint32_t hash_boundry[], query_t *query, struct SubBuf *sub,uint32_t best_pos[])
{
//fprintf(stderr, "%u\n", __LINE__);
    uint32_t i, k, l, num, num_r;
    uint32_t pos_buf[IS_SMLSIZ*4+1];
    uint32_t idx_buf[LEN_READ-SEED_LEN][2];
    uint32_t aln_buf[LEN_READ/SEED_LEN+2][5];
    uint32_t tmp[SEED_LEN-12][2];
    uint8_t *pseed;
    uint32_t seq12;

//fprintf(stderr, "%u\n", __LINE__);
    uint8_t read_seq[LEN_READ+50];
    if(query->is_rev) memcpy(read_seq, query->rseq, query->l_seq);
    else memcpy(read_seq, query->rseq, query->l_seq);
    for(i = 0; i < query->l_seq; ++i) read_seq[i] = read_seq[i] >3?0:read_seq[i];
    int read_len = query->l_seq;
    uint32_t n_aln = 0, idx_num = 0, aln_len = 0;
    int flag = 0;
    int pos_c, pos_s = read_len;
    while(1){
        pos_c = pos_s -12;
        pseed = read_seq + pos_c; 
        seq12 = lkt_seq2LktItem(pseed, 0, 11);
        k = fm_idx->fastmap->item[seq12];
        l = fm_idx->fastmap->item[seq12+1]-1;
        l -= get_12mer_correct(hash_boundry, l);
        n_aln = 12;
         
        while(2){
            if(n_aln <= SEED_LEN) { 
//fprintf(stderr, "%u, n_aln = %d\n", __LINE__, n_aln);
                tmp[n_aln-12][0] = k;
                tmp[n_aln-12][1] = l+1-k;
            }
            
            pseed = read_seq + pos_c; 
//fprintf(stderr, "%u, pos_c = %u, k = %u, l = %u\n", __LINE__, pos_c, k, l);
            num = bwt_match_exact_alt(fm_idx->bwt, 1, pseed, &k, &l);
//fprintf(stderr, "%u, pos_c = %u, k = %u, l = %u\n", __LINE__, pos_c, k, l);
            if(num > 0 ) {
                if(num <= IS_SMLSIZ) {
//fprintf(stderr, "%u, num = %u\n", __LINE__, num);          
                    for(i =0; i < num; ++i) pos_buf[i] = bwt_sa(fm_idx->bwt, k+i)-pos_c; 
                    num_r = AlgnPos(fm_idx, query, num, sub);  
//fprintf(stderr, "%u, num_r = %u, query->b0 = %u, best_thres = %u\n", __LINE__, num_r, query->b0, query->best_thres);          
                    if(num_r >0 && query->b0 >= query->best_thres) {
                        flag = 2;          
                    } else{ flag = 0;}
                    break;     


                } else{
                    if(aln_len >= SEED_LEN) {
//fprintf(stderr, "%u, num = %u\n", __LINE__, num);          
                        idx_buf[idx_num][0] = k;
                        idx_buf[idx_num][1] = num;
                        idx_num++;  
                    } 
                }
                n_aln++; 
                pos_c--;
            } else{// num == 0
                if(tmp[aln_len][1] <= 4*IS_SMLSIZ) {
                    //用SW方法进行比对
                    k = tmp[aln_len][0];
                    num = tmp[aln_len][1];
//fprintf(stderr, "%u, num = %u\n", __LINE__, num);          
                    for(i =0; i < num; ++i) pos_buf[i] = bwt_sa(fm_idx->bwt, k+i)-pos_c; 
                    num_r = AlgnPos(fm_idx, query, num, sub);  

//fprintf(stderr, "%u, num_r = %u\n", __LINE__, num_r);          
                    if(num_r >0 && query->b0 >= query->best_thres) {
                        flag = 2;          
                        break;
                    }
                } 
//fprintf(stderr, "%u, aln_len = %u, idx_num = %u\n", __LINE__, aln_len, idx_num);
                if(n_aln >= SEED_LEN) {
                    aln_buf[aln_len][0] = idx_num - n_aln + SEED_LEN;
                    aln_buf[aln_len][1] = n_aln - SEED_LEN;
                    aln_buf[aln_len][2] = pos_s;
                    aln_len++; 
                } 
                if(n_aln >= SEED_LEN) {
                    flag = 1; 
                } else{
                    flag = 0; 
                }
                break; //while(2)       
            }
          
             
            if(pos_c < 0){ 
                if(n_aln >= SEED_LEN) {
                    flag = 1; 
                } else{
                    flag = 0; 
                }
                break;
            }
           
        } // end while(2)
        if(flag == 2) { break; }
        if(flag == 1) {
            //如果新产生的种子序列尚未sw比对，那么进行多种子重叠比对 
        }
        pos_s = pos_s - n_aln;
//fprintf(stderr, "pos_s = %d, n_aln = %d, pos_c = %d\n", pos_s, n_aln, pos_c);
        if(pos_s < SEED_LEN || pos_c < 0){ 
            if(n_aln >= SEED_LEN) {
                flag = 1; 
            } else{
                flag = 0; 
            }
            break;
        }
    }// end while(1)

//fprintf(stderr, "%u\n", __LINE__);
    return flag;    
}

int gen_pos_exted0(idx_t *fm_idx, uint32_t hash_boundry[], query_t *query,  struct StackTree *sTree, struct SubBuf *sub, struct ExtBlck* eBlck, int SLC_CLS)
{
    uint32_t (*ext_idx)[2];
    uint32_t *pos_buf = sub->pos_buf; 
    int seed_off = sTree->seed_off; 
    uint32_t aln_len, bgn, end, num, cls, cap_row, err_sum;
    int flg = 0, best_flg = 0;
    int i, j;  
    int seq_off[2]; 
    uint8_t *seq = query->is_rev>0? query->rseq:query->seq; 
    while(sTree->len_exted > 0){
        cls = sTree->exted_arry[sTree->len_exted-1].cls;
        //if(cls < SLC_CLS) { break; }
        sTree->len_exted--;
        cap_row = sTree->exted_arry[sTree->len_exted].nxtpnt;
        err_sum = sTree->exted_arry[sTree->len_exted].err;
        sub->err_sum[0] = err_sum;
        if(err_sum > query->l_seq - query->b0 +sub->delta ) continue; 
        int l_off = sTree->exted_arry[sTree->len_exted].l_off;  
        int r_off = sTree->exted_arry[sTree->len_exted].r_off;  
        int L_offset = seed_off - cls*16  - l_off;
        int R_offset = seed_off + SEED_LEN + cls*16  + r_off;

        ext_idx = (eBlck+cls-1)->head_extidx + cap_row;
        bgn = ext_idx[0][0];
        num = ext_idx[0][1];
        end = bgn + num - 1; 
        if(L_offset >= query->l_seq - R_offset) { 
            if(L_offset > 4 && num > IS_SMLSIZ){
                uint32_t k = bgn, l = end;
                int n, len = 0;
                for(i = 1; i <= L_offset; ++i){ 
                    n = bwt_match_exact_alt(fm_idx->bwt, 1, seq+L_offset-i, &k, &l);
                    if(n > 0) {
                        bgn = k; 
                        end = l;      
                        num = n;
                        len++;
                    } else { break;}
                }
                L_offset -= len;   
            }
            if(query->l_seq - R_offset > 4 && num > IS_SMLSIZ){
                uint32_t bwt_idx[4];
                bwt_idx[0] = bgn;
                bwt_idx[1] = end;
                int len = R_offset - L_offset;
                seq_off[0] = L_offset;
                seq_off[1] = R_offset;
                int len_buf[2];
                len_buf[0] = len;
                len_buf[1] = 32;
                int l = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);

                //int l = bsearch_idx_R(fm_idx, query, len, seq_off, bwt_idx);
                if(l > 0) {
                    bgn = bwt_idx[2];
                    end = bwt_idx[3];
                    num = end+1-bgn;
                    seq_off[1] += l;
                    R_offset = seq_off[1];
                }
            }
        } else {
            if(query->l_seq - R_offset > 4 && num > IS_SMLSIZ){
                uint32_t bwt_idx[4];
                bwt_idx[0] = bgn;
                bwt_idx[1] = end;
                int len = R_offset - L_offset;
                seq_off[0] = L_offset;
                seq_off[1] = R_offset;
                int len_buf[2];
                len_buf[0] = len;
                len_buf[1] = 32;
                int l = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);

                //int l = bsearch_idx_R(fm_idx, query, len, seq_off, bwt_idx);
                if(l > 0) {
                    bgn = bwt_idx[2];
                    end = bwt_idx[3];
                    num = end+1-bgn;
                    seq_off[1] += l;
                    R_offset = seq_off[1];
                }
            }
            if(L_offset > 4 && num > IS_SMLSIZ){
                uint32_t k = bgn, l = end;
                int n, len = 0;
                for(i = 1; i <= L_offset; ++i){ 
                    n = bwt_match_exact_alt(fm_idx->bwt, 1, seq+L_offset-i, &k, &l);
                    if(n > 0) {
                        bgn = k; 
                        end = l;      
                        num = n;
                        len++;
                    } else { break;}
                }
                L_offset -= len;   
            }
        
        
        
        
        }



        //+++++++++++++++++++++++++++++++++++++++++++++
        //以下代码完成扩展比对结束的扩展种子进行边界比对
        uint32_t k, l;
        uint32_t out[2];
//fprintf(stderr, "bgn = %u, num = %u, cls = %d, cap_row = %d\n", bgn, num, cls, cap_row);
        for(i = 0; i < num; ++i) {
            uint32_t pos = bwt_sa(fm_idx->bwt, bgn+i)-L_offset;
//fprintf(stderr, "idx = %u, pos = %u\n", bgn+i, pos);
            if(pos > fm_idx->bwt->seq_len) {
                printf("error: pos out of range!!!\n");
                exit(1); 
            }
            pos_buf[i] = pos;
        }

num = set_sw_ed_buf(pos_buf, num, query, sub);                              
        //num = clean_swed_pos(pos_buf, num, query, sub);
        seq_off[0] = seed_off-cls*16-16       - l_off;
        seq_off[1] = SEED_LEN+seed_off+cls*16 + r_off;
        //int pos_n[3];
        //pos_n[0] = num;
        //num = slc_aln_pos(fm_idx,  query, seq_off, pos_n, sub);
        //pos_n[0] = num;
        if(num > 0){ 
            flg = 0; 
            if(sub->ALN_SWITCH[4] > 0) {
                for(i = 0; i < num; ++i) {
                    int sc = aln_pos_both(fm_idx,  query, pos_buf+i,  
                                        seq_off, sub);
                    rec_aln_info(query, pos_buf[i], sc, sub);
                }
            } else if(num > 0){
                flg = AlgnPos(fm_idx, query, num,  sub);
            }

          
            //flg = AlgnPos(fm_idx, query, num, sub);
            if(flg > 0){
                if(query->b0 >= query->best_thres && err_sum == 0) { 
                    best_flg = 1;
                } 
            }
        }
    } //end while(EXT_END > 0)+++++++++++++++++++++
    return;
} // end gen_pos_exted0()+++++++++++++++++++++

//int gen_pos_exted1(idx_t *fm_idx, uint32_t hash_boundry[], query_t *query,  struct StackTree *sTree, struct SubBuf *sub, struct ExtBlck* eBlck, int SLC_CLS)
int AlgnPos_exted1(idx_t *fm_idx, uint32_t hash_boundry[], query_t *query, seed_t *seed, struct StackTree *sTree, struct SubBuf *sub, struct ExtBlck* eBlck, int SLC_CLS)
{
    uint32_t (*ext_idx)[2];
    uint32_t *pos_buf = sub->pos_buf; 
    //int seed_off = sTree->seed_off; 
    int seed_off = seed->slc[seed->id].s_off; 
    uint32_t pos, aln_len, bgn, end, num, cls, cap_row, err_sum;
    int flg = 0, best_flg = 0;
    int i, j;  
    int seq_off[2]; 
////fprintf(stderr, "%u, cls = %u, slc_cls = %d\n", __LINE__, cls, SLC_CLS);
    uint8_t *seq = query->is_rev>0? query->rseq:query->seq; 
    int max_cls = 0, max_i = 0, min_err = query->l_seq;
    for(i = 0; i < sTree->len_exted; ++i){
        cls = sTree->exted_arry[i].cls;
        cap_row = sTree->exted_arry[i].nxtpnt;
        err_sum = sTree->exted_arry[i].err; 
        sub->err_sum[0] = err_sum;
        if(cls > max_cls) {
            max_i = i;
            max_cls = cls;
            min_err = err_sum;
        }
        if(cls == max_cls && err_sum < min_err){
            min_err = err_sum;
            max_i = i;
        }
    }
    int row = 0;
    for(row = 0; row < sTree->len_exted; ++row){
        cls = sTree->exted_arry[row].cls;
        cap_row = sTree->exted_arry[row].nxtpnt;
        err_sum = sTree->exted_arry[row].err; 
        sub->err_sum[0] = err_sum;
        if(err_sum > query->l_seq - query->b0 +sub->delta ) continue; 
        int l_off = sTree->exted_arry[row].l_off;  
        int r_off = sTree->exted_arry[row].r_off;  
        int L_offset = seed_off - cls*16  - l_off;
        int R_offset = seed_off + SEED_LEN + cls*16  + r_off;
        ext_idx = (eBlck+cls-1)->head_extidx + cap_row;
        bgn = ext_idx[0][0];
        num = ext_idx[0][1];
       
        end = bgn+num-1;
        if(L_offset >= query->l_seq - R_offset) { 
            if(L_offset > 4 && num > IS_SMLSIZ){
                uint32_t k = bgn, l = end;
                int n, len = 0;
                for(i = 1; i <= L_offset; ++i){ 
                    n = bwt_match_exact_alt(fm_idx->bwt, 1, seq+L_offset-i, &k, &l);
                    if(n > 0) {
                        bgn = k; 
                        end = l;      
                        num = n;
                        len++;
                    } else { break;}
                }
                L_offset -= len;   
            }
            if(query->l_seq - R_offset > 4 && num > IS_SMLSIZ){
                uint32_t bwt_idx[4];
                bwt_idx[0] = bgn;
                bwt_idx[1] = end;
                int len = R_offset - L_offset;
                seq_off[0] = L_offset;
                seq_off[1] = R_offset;
                int len_buf[2];
                len_buf[0] = len;
                len_buf[1] = 32;
                int l = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);

                //int l = bsearch_idx_R(fm_idx, query, len, seq_off, bwt_idx);
                if(l > 0) {
                    bgn = bwt_idx[2];
                    end = bwt_idx[3];
                    num = end+1-bgn;
                    seq_off[1] += l;
                    R_offset = seq_off[1];
                }
            }
        } else {
            if(query->l_seq - R_offset > 4 && num > IS_SMLSIZ){
                uint32_t bwt_idx[4];
                bwt_idx[0] = bgn;
                bwt_idx[1] = end;
                int len = R_offset - L_offset;
                seq_off[0] = L_offset;
                seq_off[1] = R_offset;
                int len_buf[2];
                len_buf[0] = len;
                len_buf[1] = 32;
                int l = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);
                //int l = bsearch_idx_R(fm_idx, query, len, seq_off, bwt_idx);
                if(l > 0) {
                    bgn = bwt_idx[2];
                    end = bwt_idx[3];
                    num = end+1-bgn;
                    seq_off[1] += l;
                    R_offset = seq_off[1];
                }
            }
            if(L_offset > 4 && num > IS_SMLSIZ){
                uint32_t k = bgn, l = end;
                int n, len = 0;
                for(i = 1; i <= L_offset; ++i){ 
                    n = bwt_match_exact_alt(fm_idx->bwt, 1, seq+L_offset-i, &k, &l);
                    if(n > 0) {
                        bgn = k; 
                        end = l;      
                        num = n;
                        len++;
                    } else { break;}
                }
                L_offset -= len;   
            }
        
        
        
        
        }
        //if(num > SW_THRES && row != max_i) {
        if(row != max_i) {
//continue;
            if(num < SW_THRES || (cls == max_cls && cls > 2 && num < SW_THRES*5)) {
                for(i = 0; i < num; ++i) {
                    pos = bwt_sa(fm_idx->bwt, bgn+i)-L_offset;
                    //get_set_sw_to_pos(pos, 0, 1, seed, query, sub);
                    get_set_sw_to_pos(pos, 0, 2, seed, query, sub);
                } 
            } else {
                continue; 
            }
        } 
        //+++++++++++++++++++++++++++++++++++++++++++++
        //以下代码完成扩展比对结束的扩展种子进行边界比对
        uint32_t k, l;
        uint32_t out[2];
        for(i = 0; i < num; ++i) {
            pos = bwt_sa(fm_idx->bwt, bgn+i)-L_offset;
            pos_buf[i] = pos;
        }
 
        //seq_off[0] = seed_off-cls*16          - l_off;
        //seq_off[1] = SEED_LEN+seed_off+cls*16 + r_off;
        seq_off[0] = L_offset;
        seq_off[1] = R_offset;  


        seed->cls = cls;
        num = classify_pos(fm_idx, query, seed, pos_buf, num, seq_off, sub);
        if(sub->SLC_SWITCH[3])
            num = slc_aln_pos(fm_idx,  query, seq_off, num, sub);
num = clean_hash_pos(fm_idx, pos_buf, num, query, sub);
        if(num > 0){ 
            flg = 0;
//sub->ALN_SWITCH[3] = 0;
            if(sub->ALN_SWITCH[3] > 0) {

                int query_b0;
                for(i = 0; i < num; ++i) {
                    int sc;
                    sc = aln_pos_both(fm_idx,  query, pos_buf+i,  
                                        seq_off, sub);
                    rec_aln_info(query, pos_buf[i], sc, sub);
                }
            } else {

                flg = AlgnPos(fm_idx, query, num,  sub);
                if(flg > 0){
                    if(query->b0 >= query->best_thres) {//????
                        best_flg = 1;//????
                    } 
                }
            
            }
        }
    } //end for(row = 0; row < sTree->len_exted; ++row)+++++++++
////fprintf(stderr, "%u\n", __LINE__); 
    return 0;
} // end aln_pos_exted()+++++++++++++++++++++



/*  
int aln_pos_exted(idx_t *fm_idx, uint32_t hash_boundry[], query_t *query,  struct StackTree *sTree, struct SubBuf *sub, struct ExtBlck* eBlck, int SLC_CLS)
{
    uint32_t (*ext_idx)[2];
    uint32_t *pos_buf = sub->pos_buf; 
    //aln_rlt_t *aln_r = sub->aln_r;
    int seed_off = sTree->seed_off; 
    uint32_t pos, aln_len, bgn, num, cls, cap_row, err_score;
    int flg = 0, best_flg = 0;
    int i, j;  
    int seq_off[2]; 
//fprintf(stderr, "%u, cls = %u, slc_cls = %d\n", __LINE__, cls, SLC_CLS);
    while(sTree->len_exted > 0){
        cls = sTree->exted_arry[sTree->len_exted-1].cls;

//fprintf(stderr, "%u, cls = %u, slc_cls = %d\n", __LINE__, cls, SLC_CLS);
        if(cls < SLC_CLS) { break; }
        sTree->len_exted--;
        cap_row = sTree->exted_arry[sTree->len_exted].nxtpnt;
        err_score = sTree->exted_arry[sTree->len_exted].err; 
//fprintf(stderr, "%u, err = %u\n", __LINE__, err_score);
//fprintf(stderr, "%u, cls = %u, sTree->cls = %d\n", __LINE__, cls, sTree->cls);
        if(err_score > query->l_seq - query->b0 +sub->delta ) continue; 
//fprintf(stderr, "%u, cls = %u, sTree->cls = %d\n", __LINE__, cls, sTree->cls);
        //int l_off = sTree->stck_arry[sTree->len_arry-1].l_off;  
        int l_off = sTree->stck_arry[sTree->len_exted].l_off;  
        int r_off = sTree->stck_arry[sTree->len_exted].r_off;  
        int L_offset = seed_off - cls*16  - l_off;
        sub->err_sum[0] = sTree->stck_arry[sTree->len_exted].err;
        //fprintf(stderr, "%u, l_off = %u, L_offset = %u\n", __LINE__, l_off, L_offset);
        ext_idx = (eBlck+cls-1)->head_extidx + cap_row;
        bgn = ext_idx[0][0];
        num = ext_idx[0][1];
        //fprintf(stderr, "%u, bg = %u, num = %u\n", __LINE__, bgn, num);
        //fprintf(stderr, "%u\n", __LINE__);

        //+++++++++++++++++++++++++++++++++++++++++++++
        //以下代码完成扩展比对结束的扩展种子进行边界比对
        uint32_t k, l,  num1 = 0, num2 = 0, num3 = 0, num4 = 0;
        uint32_t out[2];
        if( cls  == NUM_EXT ) {
            for(i = 0; i < num; ++i) {
                pos = bwt_sa(fm_idx->bwt, bgn+i)-L_offset;
                pos_buf[i] = pos;
            }
                      
            //fprintf(stderr, "%u num = %u\n", __LINE__, num);
            num = clean_swed_pos(pos_buf, num, query, sub);
            //fprintf(stderr, "%u num = %u\n", __LINE__, num);
            
            //seq_off[0] = seed_off-cls*16-16       - l_off;
            //seq_off[1] = SEED_LEN+seed_off+cls*16 + r_off;
            //int pos_n[3];
            //pos_n[0] = num;
            //num = slc_aln_pos(fm_idx,  query, seq_off, pos_n, sub);
            //pos_n[0] = num;


            //num1 = clean_swed_pos(sw_ed, sw_to, pos_ed, pos_to, pos_buf, num);
            //num1 = num;
    //fprintf(stderr, "%u\n", __LINE__);

            if(num > 0){ 
                flg = AlgnPos(fm_idx, query, num, sub);
                if(flg > 0){
                    if(query->b0 >= query->best_thres && err_score == 0) { 
                        best_flg = 1;
                    } 
                }
            }
            //break;  //while(EXT_END > 0)
            continue;  //while(EXT_END > 0)
        } //end if(MAX_SLC_SEED == NUM_EXT) ++++++++++  
//fprintf(stderr, "%u, cls = %u, max_cls = %u\n", __LINE__, cls, sTree->max_cls);
        if( cls == sTree->max_cls && sTree->max_cls < NUM_EXT ) { 
            
            //fprintf(stderr, "%u num = %u\n", __LINE__, num);
            for(i = 0; i < num; ++i) {
                pos = bwt_sa(fm_idx->bwt, bgn+i)-L_offset;
                pos_buf[i] = pos;
            }
            
            //fprintf(stderr, "%u num = %u\n", __LINE__, num);
            num = clean_swed_pos(pos_buf, num, query, sub);
            //num1 = clean_swed_pos(sw_ed, sw_to, pos_ed, pos_to, pos_buf, num);
            //num1 = num;
            seq_off[0] = seed_off-cls*16-16       - l_off;
            seq_off[1] = SEED_LEN+seed_off+cls*16 + r_off;
            int pos_n[3];
            pos_n[0] = num;
            num = slc_aln_pos(fm_idx,  query, seq_off, pos_n, sub);
            pos_n[0] = num;

            //fprintf(stderr, "%u num = %u\n", __LINE__, num);
            if(num > 0){ 
                flg = AlgnPos(fm_idx, query,  num, sub);
                if(flg > 0){
                    if(query->b0 >= query->best_thres) {//????
                        best_flg = 1;//????
                    } 
                }
            }
    //-----------------------------
 
        }//end if( cls == MAX_CLS_SEED && MAX_SLC_SEED < NUM_EXT )+++++++++++
    //fprintf(stderr, "%u\n", __LINE__);

        //if(best_flg > 0) break; //while(EXT_END >0)+++++++++++++++++++
    } //end while(EXT_END > 0)+++++++++++++++++++++
    //fprintf(stderr, "%u\n", __LINE__);
      
    return;
} // end aln_pos_exted()+++++++++++++++++++++
*/
int AlgnPos_1(idx_t *fm_idx, query_t *query, const uint32_t pos_buf[], int pos_num, struct SubBuf *sub){
    int j;

    //fprintf(stderr, "%u, pos_num = %u\n", __LINE__, pos_num);
    for(j = 0; j< pos_num; j++){
        //fprintf(stderr, "%u, pos_buf[j] = %u\n", __LINE__, pos_buf[j]);
    }

    int read_len = query->l_seq;
    int CANDI_THRES, delta = sub->delta;
    if(query->b0 > 0) CANDI_THRES = query->b0;
    else CANDI_THRES = query->candi_thres;
    uint32_t i;
    uint8_t target[LEN_READ+50] = {};
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t *aln_pos = sub->aln_r->pos; 
    int16_t (*aln_r)[4] = sub->aln_r->score;
    //sub->err_sum[1] = query->l_seq;
    //第i个pos右段比对结果，[i][0]是比对上1否则0, [i][1]比对分数
    kswq_t **qry, **kswq; 
    if(query->is_rev) {
        qry = &sub->qry_r[0];
        kswq = sub->kswq_r; 
    } else{
        qry = &sub->qry_f[0];
        kswq = sub->kswq_f; 
    }
 
    int8_t *mat= sub->mat;       
    uint32_t pos_i, pos;

    int n_aln = sub->aln_r->num;
    int n_best = sub->aln_r->best_num;
    if(n_best > n_aln) n_best = n_aln;
    for(pos_i = 0; pos_i < pos_num; ++pos_i){
        pos = pos_buf[pos_i];  
        uint32_t l_pos = pos > MAX_CLIP?pos-MAX_CLIP:0;
        uint32_t r_pos = pos+read_len+MAX_CLIP;
        for(i =0; i < r_pos - l_pos; ++i) target[i] = __get_pac(fm_idx->pac, l_pos+i);
        kswr_t r = ksw_align_R(read_len, read_seq, r_pos-l_pos, target, 
                                5, mat, 6, 1, KSW_XSTART, qry, kswq);
        int score_r = (2*read_len-r.qe-r.te-2 )/2 + (r.qb+r.tb)/2;
        int score = (read_len - r.score) -score_r/2;
//fprintf(stderr, "%u, score = %u, sw_score = %u, pos = %u\n", __LINE__, score, r.score, pos);
//fprintf(stderr, "%u, score_r= %u, r.qb = %u,  r.qe= %u\n", __LINE__, score_r, r.qb, r.qe);
//fprintf(stderr, "CANDI_THRES = %u, r.score = %u, n_aln = %u\n", CANDI_THRES, r.score, n_aln);
        //+++++++++++++++++++++++++++++++++++++++++++++++++++ 
        int nm_L =0, nm_R = 0;
        if(r.score >= CANDI_THRES ){// 
            int l = query->l_seq - r.qe;
            sub->err_sum[1] = query->l_seq - r.score; 
            if(r.tb < r.qb){
                nm_L = r.qb; 
            } else{
                for(i = 0; i < r.qb; ++i) {
                    nm_L += (read_seq[i] != target[r.tb-r.qb+i]);
                }
            }
            if(r.te+query->l_seq -r.qe> r_pos- l_pos) {
                nm_R = query->l_seq - r.qe; 
            } else {
                for(i = r.qe+1; i < query->l_seq;++i) 
                    nm_R += (read_seq[i] != target[r.te+i-r.qe]);
            }
            if(nm_L*2 > r.qb && r.qb>5 ) nm_L = (nm_L+r.qb/2)/2;
            if(nm_R*2 > l && l>5 ) nm_R = (nm_R+l/2)/2;
        }
        //-------------------------- 
        if(r.score == CANDI_THRES ){// 
            int nm = nm_L + nm_R;
            int row, bg = 1, ed = n_best;
//fprintf(stderr, "%u, n_best = %u\n", __LINE__, n_best);
            row = bg;
            while(bg <= ed){
                if(aln_r[row][1] > nm) {
                    if(bg == ed) break;
                    ed = row; 
                } else if(aln_r[row][1] < nm){
                    bg = row+1; 
                } else{
                    break;
                }
                row = (bg + ed)/2;     
            } 
//fprintf(stderr, "%u, row = %u\n", __LINE__, row);
            if(n_aln >= ALN_R_SIZE) n_aln = ALN_R_SIZE-1;
            if(n_best == 0) {
                row = 1;
            } else{
//fprintf(stderr, "%u, n_aln = %u, row = %u\n", __LINE__, n_aln, row);
                memmove(aln_r+row+1, aln_r+row, (n_aln-row)*4*sizeof(int16_t));
                memmove(aln_pos+row+1, aln_pos+row, (n_aln-row)*sizeof(uint32_t)); 
            } 
            aln_r[row][0] = r.score; 
            aln_r[row][1] = nm_L+nm_R; 
            aln_pos[row]  = pos;
            n_best++;
//fprintf(stderr, "%u, pos = %u, row = %d\n", __LINE__, pos, row);
            if(row == 1){
                query->b0 = r.score;
                query->pos = pos;
                query->strand = query->is_rev;
                //query->tlen = r.te+1-r.tb;
                query->ref_end = r.te+1;
                query->ref_start = r.tb;
                query->seq_start = r.qb;
                query->seq_end = r.qe+1;
                query->n_cigar = 0;
                CANDI_THRES = query->b0;
            }
            n_aln++;
        } else if(r.score > CANDI_THRES){
            memmove(aln_r+1, aln_r, n_aln*4*sizeof(int16_t));
            memmove(aln_pos+1, aln_pos, n_aln*sizeof(uint32_t)); 

            aln_r[0][0] = r.score; 
            aln_r[0][1] = nm_L + nm_R; 
            n_best = 1;
            aln_pos[0]  = pos;
//fprintf(stderr, "%u, pos = %u\n", __LINE__, pos);
            query->b0 = r.score;
            query->pos = pos;
            query->strand = query->is_rev;
            //query->tlen = r.te+1-r.tb;
            query->ref_end = r.te+1;
            query->ref_start = r.tb;
            query->seq_start = r.qb;
            query->seq_end = r.qe+1;
            query->n_cigar = 0;
            CANDI_THRES = query->b0;
            n_aln++;
            query->query_err = query->l_seq - query->b0;
            sub->query_err   = query->l_seq - query->b0;
        } else if(r.score + delta > CANDI_THRES && r.score >query->candi_thres+delta) { 

//fprintf(stderr, "%u, r.score = %d, CANDI_THRES = %u, n_aln = %u\n", __LINE__, r.score, CANDI_THRES, n_aln);
            
            if(n_aln >= ALN_R_SIZE) n_aln = ALN_R_SIZE-1;
            int qlen  = query->l_seq;
            int score, cur_score = r.score*qlen-qlen/2;
            int row, bg = n_best, ed = n_aln;
            row = bg;
            if(n_best + 1 < n_aln){
                while(bg <= ed){
                    score = aln_r[row][0]*qlen - aln_r[row][1];
                    if(score < cur_score) {
                        if(bg == ed) break;
                        ed = row; 
                    } else if(score > cur_score){
                        bg = row+1; 
                    } else{
                        break;
                    } 
                    row = (bg + ed)/2;
                } 
            
//fprintf(stderr, "%u, n_aln = %u, row = %u\n", __LINE__, n_aln, row);
        
                memmove(aln_r+row+1, aln_r+row, (n_aln-row)*4*sizeof(int16_t));
                memmove(aln_pos+row+1, aln_pos+row, (n_aln-row)*sizeof(uint32_t)); 
            } else if( n_best < n_aln){
                score = aln_r[row][0]*qlen - aln_r[row][1];
                if(score < cur_score ) {
                    memmove(aln_r+n_aln+1, aln_r+n_aln, 4*sizeof(int16_t));
                    memmove(aln_pos+n_aln+1, aln_pos+n_aln, sizeof(uint32_t)); 
                }
                row = n_aln; 
            } else{
                row++;
            }
            aln_r[row][0] = r.score; 
            aln_r[row][1] = qlen/10; 
            aln_pos[row]  = pos;
            n_aln++;
        }
        if(n_aln >= ALN_R_SIZE) n_aln = ALN_R_SIZE-1;
        if(n_best > n_aln) n_best = n_aln;
    }//end for(pos_i = 0; pos_i < pos_num; ++pos_i)++++++

//fprintf(stderr, "%u, n_aln = %u, aln_r[1][0] = %u\n", __LINE__, n_aln, aln_r[1][0]);
    if(aln_r[n_aln][0] < aln_r[1][0] - delta) {
        int score, best_score = aln_r[1][0] - delta;
        int row, bg = n_best+1, ed = n_aln;
        row = bg;
        while(1){
            while(bg <= ed){
                score = aln_r[row][0];
                if(score < best_score) {
                    if(bg == ed) break;
                    ed = row; 
                } else if(score > best_score){
                    bg = row+1; 
                } else{
                    break;
                } 
                row = (bg + ed)/2;
            } 
            if(bg >= ed) break;
            ed = row;         
            row = (bg+ed)/2;
        }
//fprintf(stderr, "%u, n_aln = %u\n", __LINE__, n_aln);
        n_aln = row; 
//fprintf(stderr, "%u, n_aln = %u\n", __LINE__, n_aln);
    }
//fprintf(stderr, "%u, n_aln = %u\n", __LINE__, n_aln);
    sub->aln_r->num = n_aln;
    sub->aln_r->best_num = n_best;

    return n_aln;
}

//int AlgnPos(idx_t *fm_idx, query_t *query, int pos_num, uint8_t in_flg, struct SubBuf *sub){
int AlgnPos1(idx_t *fm_idx, query_t *query, int pos_num, uint32_t *pos_buf, struct SubBuf *sub){
    kswq_t **qry, **kswq;
    int8_t *mat= sub->mat;          
    int j;
    //uint32_t (*err_buf)[5] = sub->err_buf;  
    int read_len = query->l_seq;
    int CANDI_THRES;
    int delta = sub->delta;
    if(query->b0 > 0) CANDI_THRES = query->b0;
    else CANDI_THRES = query->candi_thres;
    uint32_t i;
    uint8_t *target = query->target;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t drct_f = query->is_rev*1024;
    aln_out_t *out = sub->aln_out;
    uint32_t (*out_buf)[2] = sub->aln_out->out_buf; 
    int (*found)[4] = sub->aln_out->found; 
    int n_fnd, l_fnd; 
    
   
    if(query->is_rev) {
        qry = &sub->qry_r[0];
        kswq = sub->kswq_r; 
    } else{
        qry = &sub->qry_f[0];
        kswq = sub->kswq_f; 
    }
 
    uint32_t pos_i, pos;
    int cur_s;
    
    for(pos_i = 0; pos_i < pos_num; ++pos_i){
        pos = pos_buf[pos_i];
        //if((query->b0 >= err_buf[pos_i][0] + delta) && query->b0 >0 && in_flg) continue; 
        uint32_t l_pos = pos > MAX_CLIP?pos-MAX_CLIP:0;
        uint32_t r_pos = pos + read_len + MAX_CLIP;
        for(i =0; i < r_pos - l_pos; ++i) 
            target[i] = __get_pac(fm_idx->pac, l_pos+i);
            //kswr_t r = ksw_align_R(read_len, read_seq, r_pos-l_pos, target, 5, mat, 6, 1, KSW_XSTART, qry, kswq);
/*  
{
    int ii;
    for(ii = 0; ii < query->l_seq; ++ii) {
        //fprintf(stderr, "i = %d, t = %d, r = %d\n", ii, target[ii+pos-l_pos], read_seq[ii]);
    }
    for(ii; ii < query->l_seq + 2*MAX_CLIP; ++ii) {
        //fprintf(stderr, "i = %d, t = %d, r = --\n", ii, target[ii]);
    }
}            
*/
        
            kswq_t *qry[2] = {0, 0};       
            kswr_t r = ksw_align_R(read_len, read_seq, r_pos-l_pos, target, 5, mat, 6, 1, 0, qry, kswq);
//int score_r = (2*read_len-r.qe-r.te-2 )/2 + (r.qb+r.tb)/2;
//int score = (read_len - r.score) -score_r/2;
fprintf(stderr, "%u, pos = %u, [%u, %u], r.score = %d, r.te = %d, r.qe = %d, query->b0 = %u\n", __LINE__, pos, l_pos, r_pos, r.score, r.te, r.qe, query->b0);
        //-------------------------- 
        if(r.score > CANDI_THRES){
            cur_s = r.score - found[0][0];
            if(cur_s < 0) {
                printf("%u, found[0][0] = %d, cur_s = %d, r.score = %d\n", __LINE__, found[0][0], cur_s,r.score);     
                return; 
            } 
            if(cur_s >= delta) {
                memset(found, 0, 4*delta*sizeof(uint32_t));
            } else if(cur_s < delta) {
                memmove(found+cur_s, found, 4*cur_s*sizeof(uint32_t)); 
                memset(found, 0, 4*cur_s*sizeof(uint32_t));
            } 
            found[0][0] = r.score;
            found[0][1] = 1;
            found[0][2] = out->len;
            found[0][3] = out->len;
            
            out_buf[out->len][0] = 1;  
            out_buf[out->len][1] = 0; 
            out_buf[out->len+1][0] = pos;  
            out_buf[out->len+1][1] = drct_f + r.score;  
            out->len += 2;
            out->num++;

            query->b0 = r.score;
            query->pos = pos;
            query->strand = query->is_rev;
            //query->tlen = r.te+1-r.tb;
            query->ref_end = r.te+1;
            query->ref_start = r.tb;
            query->seq_start = r.qb;
            query->seq_end = r.qe+1;
            query->n_cigar = 0;
            CANDI_THRES = query->b0;
////fprintf(stderr, "%u, query->b0 = %u, pos = %u, ref_end = %u, ref_start = %u, seq_start = %u\n", __LINE__, query->b0, query->pos, query->ref_end, query->ref_start, query->seq_start);            
            query->query_err = query->l_seq - query->b0;
            sub->query_err   = query->l_seq - query->b0;
 
        } else if(r.score + delta > CANDI_THRES && found[0][0] >= CANDI_THRES){ 

            int flg = 0; 
            cur_s = found[0][0]-r.score;
            if(cur_s < 0) {
                printf("%u, found[0][0] = %d, cur_s = %d, r.score = %d\n", __LINE__, found[0][0], cur_s, r.score);     
                continue; 
            }
            if(found[cur_s][1] == 0){
                found[cur_s][0] = r.score;
                found[cur_s][1] = 1;
                found[cur_s][2] = out->len;
                found[cur_s][3] = out->len;
                
                out_buf[out->len][0] = 1;  
                out_buf[out->len][1] = 0; 
                out_buf[out->len+1][0] = pos;  
                out_buf[out->len+1][1] = drct_f + r.score;  
                out->len += 2;
                out->num++;
            } else{
                found[cur_s][0] = r.score;
                for(i = 0; i < delta; ++i) {
                    if( found[i][3] > found[cur_s][3]) {
                        flg = 1;
                        j = i;
                        break; 
                    }
                }
                if(flg == 0) {
                    found[cur_s][1]++;
                    out_buf[found[cur_s][3]][0]++;  
                    out_buf[out->len][0] = pos;  
                    out_buf[out->len][1] = drct_f + r.score;
                    out->len++;
                    out->num++;
                } else{
                    found[cur_s][1]++;
                    out_buf[found[cur_s][3]][1] = out->len; 
                    found[cur_s][3] = out->len;
                    out_buf[out->len][0] = 1;  
                    out_buf[out->len][1] = 0; 
                    out_buf[out->len+1][0] = pos;  
                    out_buf[out->len+1][1] = drct_f + r.score;  
                    out->len += 2;
                    out->num++;
                }
            } 
        }// end if(r.score + delta > CANDI_THRES)
    }//end for(pos_i = 0; pos_i < pos_num; ++pos_i)++++++
    return 0;
}


int AlgnPos(idx_t *fm_idx, query_t *query, int pos_num, struct SubBuf *sub){
    int rev_flg = 0;
    kswq_t **qry, **kswq;
    int8_t *mat= sub->mat;          
    int j;
    uint32_t *pos_buf = sub->pos_buf;

    //uint32_t (*err_buf)[5] = sub->err_buf;  
    int read_len = query->l_seq;
    int CANDI_THRES;
    int delta = sub->delta;
    if(query->b0 > 0) CANDI_THRES = query->b0;
    else CANDI_THRES = query->candi_thres;
    uint32_t i;
    uint8_t *target = query->target;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t drct_f = query->is_rev*1024;
    aln_out_t *out = sub->aln_out;
    uint32_t (*out_buf)[2] = sub->aln_out->out_buf; 
    int (*found)[4] = sub->aln_out->found; 
    int n_fnd, l_fnd; 
    
   
    if(query->is_rev) {
        qry = &sub->qry_r[0];
        kswq = sub->kswq_r; 
    } else{
        qry = &sub->qry_f[0];
        kswq = sub->kswq_f; 
    }
 
    uint32_t pos_i, pos;
    int cur_s;
    
    for(pos_i = 0; pos_i < pos_num; ++pos_i){
        pos = pos_buf[pos_i];
        //if((query->b0 >= err_buf[pos_i][0] + delta) && query->b0 >0 && in_flg) continue; 
        uint32_t l_pos = pos > MAX_CLIP?pos-MAX_CLIP:0;
        uint32_t r_pos = pos + read_len + MAX_CLIP;
        for(i =0; i < r_pos - l_pos; ++i) 
            target[i] = __get_pac(fm_idx->pac, l_pos+i);
            kswq_t *qry[2] = {0, 0};       
            kswr_t r = ksw_align_R(read_len, read_seq, r_pos-l_pos, target, 5, mat, 6, 1, 0, qry, kswq);
        //-------------------------- 
if(r.score >= CANDI_THRES) rev_flg = 1;
        if(r.score > CANDI_THRES){

            cur_s = r.score - found[0][0];
            if(cur_s < 0) {
                printf("%u, found[0][0] = %d, cur_s = %d, r.score = %d\n", __LINE__, found[0][0], cur_s, r.score);     
                continue; 
            }
            if(cur_s >= delta) {
                memset(found, 0, 4*delta*sizeof(uint32_t));
            } else if(cur_s < delta) {
                memmove(found+cur_s, found, 4*cur_s*sizeof(uint32_t)); 
                memset(found, 0, 4*cur_s*sizeof(uint32_t));
            } 
            found[0][0] = r.score;
            found[0][1] = 1;
            found[0][2] = out->len;
            found[0][3] = out->len;
                
            if(out->len< 0 ||  out->len > out->max_out) {
                printf("%u, out_len = %d, max_out = %d\n", __LINE__, out->len, out->max_out);
                exit(1);            
            
            }

            out_buf[out->len][0] = 1;  
            out_buf[out->len][1] = 0; 
            out_buf[out->len+1][0] = pos;  
            out_buf[out->len+1][1] = drct_f + r.score;  
            out->len += 2;
            out->num++;

            query->b0 = r.score;
            query->pos = pos;
            query->strand = query->is_rev;
            //query->tlen = r.te+1-r.tb;
            query->ref_end = r.te+1;
            query->ref_start = r.tb;
            query->seq_start = r.qb;
            query->seq_end = r.qe+1;
            query->n_cigar = 0;
            CANDI_THRES = query->b0;

            query->query_err = query->l_seq - query->b0;
            sub->query_err   = query->l_seq - query->b0;
 
        } else if(r.score + delta > CANDI_THRES && found[0][0] >= CANDI_THRES){ 

            int flg = 0; 
            cur_s = found[0][0]-r.score;
            if(cur_s < 0) {
                printf("%u, found[0][0] = %d, cur_s = %d, r.score = %d\n", __LINE__, found[0][0], cur_s, r.score);     
                continue; 
            } 
            if(found[cur_s][1] == 0){
                found[cur_s][0] = r.score;
                found[cur_s][1] = 1;
                found[cur_s][2] = out->len;
                found[cur_s][3] = out->len;
                if(out->len< 0 ||  out->len > out->max_out) {
                    printf("%u, out_len = %d, max_out = %d\n", __LINE__, out->len, out->max_out);
                    exit(1);            
            
                }
               
                out_buf[out->len][0] = 1;  
                out_buf[out->len][1] = 0; 
                out_buf[out->len+1][0] = pos;  
                out_buf[out->len+1][1] = drct_f + r.score;  
                out->len += 2;
                out->num++;
            } else{
                found[cur_s][0] = r.score;
                for(i = 0; i < delta; ++i) {
                    if( found[i][3] > found[cur_s][3]) {
                        flg = 1;
                        j = i;
                        break; 
                    }
                }
                if(flg == 0) {
                    found[cur_s][1]++;
                    
                    out_buf[found[cur_s][3]][0]++;  
                    out_buf[out->len][0] = pos;  
                    out_buf[out->len][1] = drct_f + r.score;
                    out->len++;
                    out->num++;
                } else{
                    found[cur_s][1]++;
                    out_buf[found[cur_s][3]][1] = out->len; 
                    found[cur_s][3] = out->len;
                    out_buf[out->len][0] = 1;  
                    out_buf[out->len][1] = 0; 
                    out_buf[out->len+1][0] = pos;  
                    out_buf[out->len+1][1] = drct_f + r.score;  
                    out->len += 2;
                    out->num++;
                }
            } 
        }// end if(r.score + delta > CANDI_THRES)
    }//end for(pos_i = 0; pos_i < pos_num; ++pos_i)++++++
    return rev_flg;
}
int AlgnPos_test(idx_t *fm_idx, query_t *query, int pos_i, int seq_off[2], struct SubBuf *sub){
    kswq_t **qry, **kswq;
    int8_t *mat= sub->mat;          
    int j;
    uint32_t *pos_buf = sub->pos_buf;

    //uint32_t (*err_buf)[5] = sub->err_buf;  
    int read_len = query->l_seq;
    int CANDI_THRES;
    int delta = sub->delta;
    if(query->b0 > 0) CANDI_THRES = query->b0;
    else CANDI_THRES = query->candi_thres;
    int i;
    uint8_t *target = query->target;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t drct_f = query->is_rev*1024;
    aln_out_t *out = sub->aln_out;
    uint32_t (*out_buf)[2] = sub->aln_out->out_buf; 
    int (*found)[4] = sub->aln_out->found; 
    int n_fnd, l_fnd; 
    
   
    if(query->is_rev) {
        qry = &sub->qry_r[0];
        kswq = sub->kswq_r; 
    } else{
        qry = &sub->qry_f[0];
        kswq = sub->kswq_f; 
    }
 
    uint32_t pos;
    int cur_s;
    kswr_t r;
    //for(pos_i = 0; pos_i < pos_num; ++pos_i){
    {
        pos = pos_buf[pos_i];
        //if((query->b0 >= err_buf[pos_i][0] + delta) && query->b0 >0 && in_flg) continue; 
        uint32_t l_pos = pos > MAX_CLIP?pos-MAX_CLIP:0;
        uint32_t r_pos = pos+read_len+MAX_CLIP;
        for(i =0; i<r_pos-l_pos; ++i) target[i] = __get_pac(fm_idx->pac, l_pos+i);
/*  
for(i = seq_off[0] -10; i < read_len; ++i) {
    //fprintf(stderr, "i = %d, t = %d, r = --\n", i, target[i]);
}
*/
        //for(i = seq_off[0]; i<read_len; ++i) target[i+pos-l_pos] = read_seq[i]; //??????
        //for(i = 0; i < seq_off[1]; ++i) target[i+pos-l_pos] = read_seq[i]; //??????
/*  
//fprintf(stderr, "%u, pos_i = %u, pos = %u, [%u, %u], [%d, %d]\n", __LINE__, pos_i, pos, l_pos, r_pos, seq_off[0], read_len);

for(i =pos-l_pos-1; i >=0; --i) {
    //fprintf(stderr, "i = %d, t = %d, r = --\n", i, target[i]);
}
*/

/*  
for(i = query->l_seq-1; i >= 0; --i) {
    //fprintf(stderr, "i = %d, t = %d, r = %d\n", i, target[i+pos-l_pos], read_seq[i]);
}
*/
/*  
for(i = 0; i < query->l_seq; ++i) {
    //fprintf(stderr, "i = %d, t = %d, r = %d\n", i, target[i+pos-l_pos], read_seq[i]);
}
for(i; i < query->l_seq + 2*MAX_CLIP; ++i) {
    //fprintf(stderr, "i = %d, t = %d, r = --\n", i, target[i]);
}
*/
r = ksw_align_R(read_len, read_seq, r_pos-l_pos, target, 5, mat, 6, 1, KSW_XSTART, qry, kswq);
fprintf(stderr, "%u, ref[%d, %d] -> read[%d, %d]\n", __LINE__, r.tb, r.te, r.qb, r.qe);
fprintf(stderr, "%u, cigar = \t", __LINE__);
          
        int n_cigar = 0, *cigar; 
        int score = ksw_global(r.qe+1-r.qb, read_seq+r.qb, r.te+1-r.tb, target+r.tb, 5, mat, 6, 1, 50, &n_cigar, &cigar, &sub->kswgb);
        
        return r.score;
/*  
        if(r.qb != 0) fprintf(stderr, "%uS", r.qb);
        for(i = 0; i < n_cigar; ++i){ 
            if((i == 0|| i == n_cigar-1) && (cigar[i] &3) ==2 ) {
                continue;
            }
            fprintf(stderr, "%u%c", cigar[i]>>4, "MID"[cigar[i]&3]);
        }
        if(r.qe+1 != query->l_seq) fprintf(stderr, "%uS", query->l_seq-(r.qe+1));
        fprintf(stderr, "\n");
*/
        int nm = 0, sc = 0;
        int n_exact = 0;
        int gap_o = 6, gap_e = 1, match = -4, mis = 1;
        int max_mov = MAX_CLIP;
/*backward*/
/*         
        uint8_t *x0 = target+r.te, *x1 = read_seq+r.qe;
        for(i = n_cigar-1; i >= 0; --i){
            int op = cigar[i] & 0xF;
            int len = cigar[i] >>4;
            int ti, qi;
            if(op == 0) {
                for(j =0; j < len; ++j, --x0, --x1) {
                    
                    ti = x0- target;
                    qi = x1 - read_seq;
                    if(*x0 != *x1) {
fprintf(stderr, "%3u\t%3u\t%d\t%d\tS\t%d\n", ti, qi, *x0, *x1, sc);              
                        ++nm;
                        sc += match;
                    } else {
fprintf(stderr, "%3u\t%3u\t%d\t%d\t=\t%d\n", ti, qi, *x0, *x1, sc);              
                        if( qi >= seq_off[0] && qi < seq_off[1] && ti - qi == max_mov) n_exact++;
                        sc += mis;
                    } 
                }
            } else if(op == 1) {
                for(j =len-1; j >= 0; --j, --x1) {
                    ti = x0- target;
                    qi = x1 - read_seq;
                 
fprintf(stderr, "%3u\t%3u\t%c\t%d\tI\t%d\n", ti, qi, '-', *x1, sc);              
                    if(j == len -1) sc -= gap_o;
                    sc -= gap_e;
                }
                nm += len;
                //x1 += len; 
            } else{
        
                for(j = len-1; j >= 0; --j, --x0) {
                    ti = x0 - target;
                    qi = x1 - read_seq;
                
fprintf(stderr, "%3u\t%3u\t%d\t%c\tD\t%d\n", ti, qi, *x0, '-', sc);              
                    if(j == len -1) sc -= gap_o;
                    sc -= gap_e;
                }

                nm += len;
                //x0 += len;
            }
        }
*/
/*  
        uint8_t *x0 = target+r.tb, *x1 = read_seq+r.qb;
        for(i = 0; i < n_cigar; ++i){
            int op = cigar[i] & 0xF;
            int len = cigar[i] >>4;
            if(op == 0) {
                for(j =0; j < len; ++j, ++x0, ++x1) {
                    
                    if(*x0 != *x1) {
//fprintf(stderr, "%3u\t%3u\t%d\t%d\tS\t%d\n", (x0-target), (x1-read_seq), *x0, *x1,sc);              
                        ++nm;
                        sc += match;
                    } else {
//fprintf(stderr, "%3u\t%3u\t%d\t%d\t=\t%d\n", (x0-target), (x1-read_seq), *x0, *x1, sc);              
                        sc += mis;
                    } 
                }
            } else if(op == 1) {
                for(j = 0; j < len; ++j, ++x1) {
                  
//fprintf(stderr, "%3u\t%3u\t%c\t%d\tI\t%d\n", (x0-target), (x1-read_seq), '-', *x1, sc);              
                    if(j == 0) sc -= gap_o;
                    sc -= gap_e;
                }
                nm += len;
                //x1 += len; 
            } else{
        
                for(j = 0; j < len; ++j, ++x0) {
                  
//fprintf(stderr, "%3u\t%3u\t%d\t%c\tD\t%d\n", (x0-target), (x1-read_seq), *x0, '-', sc);              
                    if(j == 0) sc -= gap_o;
                    sc -= gap_e;
                }

                nm += len;
                //x0 += len;
            }
        }
*/    
        
//fprintf(stderr, "\t nm = %d, seq_off = (%d, %d)\n", nm, seq_off[0], seq_off[1]);
        if(n_exact != SEED_LEN) {
            fprintf(stderr, "SEED not proper matched!!!!\n");
            return 0; 
        }

        //if(r.tb != 12 || r.qb != 0 || abs(r.te - r.tb - (r.qe - r.qb)) > 20) return 0;
       
        int l_t = seq_off[0] + max_mov - r.tb;
        int l_q = seq_off[0] - r.qb;
        int r_t = r.te - max_mov - seq_off[1];
        int r_q = r.qe - seq_off[1];
        int l_len, r_len; 
        if(r.tb - max_mov >= r.qb) {  //read上种子左侧有插入
            l_len = r.tb - max_mov - r.qb; 
        } else { //ref上种子左侧有插入
            l_len = r.qb - (r.tb - max_mov); 
        }
        
        if(r.te - max_mov >= r.qe) {  //read上种子左侧有插入
            r_len = r.te - max_mov - r.qe; 
        } else { //ref上种子左侧有插入
            r_len = r.qe - (r.te - max_mov); 
        }
        if(l_t < 0 || l_q < 0 || r_t < 0 || r_q < 0) return 0;
        if(r_len + l_len > 20 ) return 0;
        
        if(r.tb - max_mov >= r.qb) {              
            l_len = r.tb - max_mov - r.qb; 
            r_len = r.te - max_mov - r.qe; 
        } else {             
            l_len = r.qb - (r.tb - max_mov); 
            r_len = r.qe - (r.te - max_mov); 
        }
        if(l_len > 0 && l_len == r_len) return 0;
        return r.score; 
        //-------------------------- 
        if(r.score > CANDI_THRES){
            cur_s = r.score - found[0][0];
            if(cur_s < 0) {
                printf("%u, found[0][0] = %d, cur_s = %d, score = %d\n", __LINE__, found[0][0], cur_s, score);     
                exit(1); 
            } 
            if(cur_s >= delta) {
                memset(found, 0, 4*delta*sizeof(uint32_t));
            } else if(cur_s < delta) {
                memmove(found+cur_s, found, 4*cur_s*sizeof(uint32_t)); 
                memset(found, 0, 4*cur_s*sizeof(uint32_t));
            } 
            found[0][0] = r.score;
            found[0][1] = 1;
            found[0][2] = out->len;
            found[0][3] = out->len;
            
            out_buf[out->len][0] = 1;  
            out_buf[out->len][1] = 0; 
            out_buf[out->len+1][0] = pos;  
            out_buf[out->len+1][1] = drct_f + r.score;  
            out->len += 2;
            out->num++;

            query->b0 = r.score;
            query->pos = pos;
            query->strand = query->is_rev;
            //query->tlen = r.te+1-r.tb;
            query->ref_end = r.te+1;
            query->ref_start = r.tb;
            query->seq_start = r.qb;
            query->seq_end = r.qe+1;
            query->n_cigar = 0;
            CANDI_THRES = query->b0;
////fprintf(stderr, "%u, query->b0 = %u, pos = %u, ref_end = %u, ref_start = %u, seq_start = %u\n", __LINE__, query->b0, query->pos, query->ref_end, query->ref_start, query->seq_start);            
            query->query_err = query->l_seq - query->b0;
            sub->query_err   = query->l_seq - query->b0;
 
        } else if(r.score + delta > CANDI_THRES && found[0][0] >= CANDI_THRES){ 
            int flg = 0; 
            cur_s = found[0][0]-r.score;
            if(cur_s < 0) {
                printf("%u, found[0][0] = %d, cur_s = %d, score = %d\n", __LINE__, found[0][0], cur_s, score);     
                exit(1); 
            } 
            if(found[cur_s][1] == 0){
                found[cur_s][0] = r.score;
                found[cur_s][1] = 1;
                found[cur_s][2] = out->len;
                found[cur_s][3] = out->len;
                
                out_buf[out->len][0] = 1;  
                out_buf[out->len][1] = 0; 
                out_buf[out->len+1][0] = pos;  
                out_buf[out->len+1][1] = drct_f + r.score;  
                out->len += 2;
                out->num++;
            } else{
                found[cur_s][0] = r.score;
                for(i = 0; i < delta; ++i) {
                    if( found[i][3] > found[cur_s][3]) {
                        flg = 1;
                        j = i;
                        break; 
                    }
                }
                if(flg == 0) {
                    found[cur_s][1]++;
                    out_buf[found[cur_s][3]][0]++;  
                    out_buf[out->len][0] = pos;  
                    out_buf[out->len][1] = drct_f + r.score;
                    out->len++;
                    out->num++;
                } else{
                    found[cur_s][1]++;
                    out_buf[found[cur_s][3]][1] = out->len; 
                    found[cur_s][3] = out->len;
                    out_buf[out->len][0] = 1;  
                    out_buf[out->len][1] = 0; 
                    out_buf[out->len+1][0] = pos;  
                    out_buf[out->len+1][1] = drct_f + r.score;  
                    out->len += 2;
                    out->num++;
                }
            } 
        }// end if(r.score + delta > CANDI_THRES)
    }//end for(pos_i = 0; pos_i < pos_num; ++pos_i)++++++
    return r.score;
}

int rec_bwt_exact(idx_t *fm_idx, query_t *query, int len, struct SubBuf *sub)
{
    int CANDI_THRES;
    int delta = sub->delta;
    if(query->b0 > 0) CANDI_THRES = query->b0;
    else CANDI_THRES = query->candi_thres;

    uint32_t(*idx_buf)[2] = query->is_rev?sub->idx_rev:sub->idx_for;  
    int r_pos = query->l_seq - sub->trm_r - len;
    uint32_t st = idx_buf[r_pos][0];
    uint32_t ed = idx_buf[r_pos][1];
    int j;
    uint32_t i;
    uint32_t drct_f = query->is_rev*1024;
    aln_out_t *out = sub->aln_out;
    uint32_t (*out_buf)[2] = sub->aln_out->out_buf; 
    int (*found)[4] = sub->aln_out->found; 
    int n_fnd, l_fnd; 
    
    uint32_t pos_i;
    int cur_s;
    //-------------------------- 
    uint32_t idx, pos;
    int score = len;

    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    int trm_r = sub->trm_r, l_seq = query->l_seq; 
    uint8_t c;
    uint32_t bg_pos, ed_pos; 
    for(idx = st; idx <= ed; ++idx){
        pos = bwt_sa(fm_idx->bwt, idx) - r_pos; 
        set_sw_ed_pos(pos, query, sub);
        score = len;
        bg_pos = r_pos;
        ed_pos = query->l_seq- trm_r-1;
        for(j = r_pos-1; j >= 0; --j){
            c = __get_pac(fm_idx->pac, pos+j);
            if(c == read_seq[j]) {
                bg_pos = j;
                score++;
            } else {
                break;
            }
        }
        for(j = query->l_seq - trm_r; j < l_seq; ++j){
            c = __get_pac(fm_idx->pac, pos+j);
            if(c == read_seq[j]) {
                ed_pos = j;
                score++;
            } else {
                break;
            }
        }
         

        if(score > CANDI_THRES){
            cur_s = score - found[0][0];
            if(cur_s < 0) {
                printf("%u, found[0][0] = %d, cur_s = %d, score = %d\n", 
                        __LINE__, found[0][0], cur_s, score);     
                return; 
            }        
            if(cur_s >= delta) {
                memset(found, 0, 4*delta*sizeof(uint32_t));
            } else if(cur_s < delta) {
                if(cur_s > sub->aln_out->max_found) {
                    printf("%u, found[0][0] = %d, cur_s = %d\n", 
                            __LINE__, found[0][0], cur_s);     
                    exit(1); 
                } 
                memmove(found+cur_s, found, 4*cur_s*sizeof(uint32_t)); 
                memset(found, 0, 4*cur_s*sizeof(uint32_t));
            } 
            found[0][0] = score;
            found[0][1] = 1;
            found[0][2] = out->len;//begin idx of out_buf
            found[0][3] = out->len;
            
            out_buf[out->len][0] = 1;//number of pos  
            out_buf[out->len][1] = 0;//next pointer of pos 
            out_buf[out->len+1][0] = pos;  
            out_buf[out->len+1][1] = drct_f + score;  
            out->len += 2;
            out->num++;

            query->b0 = score;
            query->pos = pos;
            query->strand = query->is_rev;
            
            query->ref_start = bg_pos + MAX_CLIP; 
            query->ref_end   = ed_pos + 1 + MAX_CLIP; 
            query->seq_start = bg_pos; 
            query->seq_end   = ed_pos + 1; 

            query->n_cigar = 0;
            CANDI_THRES = query->b0;
            query->query_err = query->l_seq - query->b0;
            sub->query_err   = query->l_seq - query->b0;
        } else if(score + delta > CANDI_THRES && found[0][0] >= CANDI_THRES){ 
            int flg = 0; 
            cur_s = found[0][0]-score;
            if(cur_s > sub->aln_out->max_found || cur_s < 0) {
                printf("%u, found[0][0] = %d, cur_s = %d, score = %d\n", 
                        __LINE__, found[0][0], cur_s, score);     
                exit(1); 
            } 
            
            if(found[cur_s][1] == 0){
                found[cur_s][0] = score;
                found[cur_s][1] = 1;
                found[cur_s][2] = out->len;
                found[cur_s][3] = out->len;
                out_buf[out->len][0] = 1;  
                out_buf[out->len][1] = 0; 
                out_buf[out->len+1][0] = pos;  
                out_buf[out->len+1][1] = drct_f + score;  
                out->len += 2;
                out->num++;
            } else{
                found[cur_s][0] = score;
                for(i = 0; i < delta; ++i) {
                    if( found[i][3] > found[cur_s][3]) {
                        flg = 1;
                        j = i;
                        break; 
                    }
                }
                if(flg == 0) {
                    found[cur_s][1]++;
                    out_buf[found[cur_s][3]][0]++;  
                    out_buf[out->len][0] = pos;  
                    out_buf[out->len][1] = drct_f + score;
                    out->len++;
                    out->num++;
                } else{
                    found[cur_s][1]++;
                    out_buf[found[cur_s][3]][1] = out->len; 
                    found[cur_s][3] = out->len;
                    out_buf[out->len][0] = 1;  
                    out_buf[out->len][1] = 0; 
                    out_buf[out->len+1][0] = pos;  
                    out_buf[out->len+1][1] = drct_f + score;  
                    out->len += 2;
                    out->num++;
                }
            } 
        }// end if(r.score + delta > CANDI_THRES)
    }
    n_fnd  = 0;
    for(i =0; i < delta; ++i) {
        n_fnd += found[i][1];
    }
    return n_fnd;
}
/*  */
int rec_bwt_exact_0(idx_t *fm_idx, query_t *query, int len, struct SubBuf *sub)
{
    uint32_t(*idx_buf)[2] = query->is_rev?sub->idx_rev:sub->idx_for;  
    int r_pos = query->l_seq - sub->trm_r - len;
    uint32_t st = idx_buf[r_pos][0];
    uint32_t ed = idx_buf[r_pos][1];
    int j;
    uint32_t i;
    uint32_t drct_f = query->is_rev*1024;
    aln_out_t *out = sub->aln_out;
    uint32_t (*out_buf)[2] = sub->aln_out->out_buf; 
    int (*found)[4] = sub->aln_out->found; 
    int n_fnd, l_fnd; 
    
    uint32_t pos_i;
    int cur_s;
    //-------------------------- 
    uint32_t idx, pos;
    int score = len;
    int delta = 1;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    int trim = sub->trm_r, l_seq = query->l_seq; 
    uint8_t c;
    uint32_t bg_pos, ed_pos; 
    for(idx = st; idx <= ed; ++idx){
        pos = bwt_sa(fm_idx->bwt, idx) - r_pos; 
//printf("%u, idx = %u, pos = %u\n", __LINE__, idx, pos);         
        score = len;
        bg_pos = r_pos;
        ed_pos = query->l_seq- trim-1;
        for(i = r_pos-1; i >= 0; --i){
            c = __get_pac(fm_idx->pac, pos+i);
            if(c == read_seq[i]) {
                bg_pos = i;
                score++;
            } else {
                break;
            }
        }
        for(i = query->l_seq - trim; i < l_seq; ++i){
            c = __get_pac(fm_idx->pac, pos+i);
            if(c == read_seq[i]) {
                ed_pos = i;
                score++;
            } else {
                break;
            }
        }
        cur_s = score - found[0][0];
        if(cur_s < 0) {
            printf("%u, found[0][0] = %d, cur_s = %d, score = %d\n", __LINE__, found[0][0], cur_s, score);     
            return; 
        }        
        if(cur_s >= delta) {
            memset(found, 0, 4*delta*sizeof(uint32_t));

        } else if(cur_s < delta) {
            if(cur_s > sub->aln_out->max_found) {
                printf("%u, found[0][0] = %d, cur_s = %d\n", __LINE__, found[0][0], cur_s);     
                exit(1); 
            } 
            memmove(found+cur_s, found, 4*cur_s*sizeof(uint32_t)); 
            memset(found, 0, 4*cur_s*sizeof(uint32_t));
        } 
        //found[0][0] alignment score
        //found[0][1] found num
        //found[0][2] bgn address in out_buf where alignment score is found[0][0]
        //found[0][3] end address in out_buf where alignment score is found[0][0]
        if(found[0][0] < score) {
            found[0][0] = score;
            found[0][1] = 1;
            found[0][2] = out->len;
        } else {
++found[0][1]; 
            out_buf[found[0][3]][1] = out->len; 
        }
        found[0][3] = out->len;       
        
        out_buf[out->len][0] = 1;  
        out_buf[out->len][1] = 0; 
        out_buf[out->len+1][0] = pos;  
        out_buf[out->len+1][1] = drct_f + score;  
        out->len += 2;
        out->num++;
//if( idx == st) {
        if(score > query->b0) {
            query->b0 = score;
            query->pos = pos;
            query->strand = query->is_rev;
            
        /*  
            query->ref_start = query->l_seq - len + MAX_CLIP; 
            query->ref_end   = query->l_seq + MAX_CLIP; 
            query->seq_start = query->l_seq - len; 
            query->seq_end   = query->l_seq; 
        */
            query->ref_start = bg_pos + MAX_CLIP; 
            query->ref_end   = ed_pos + MAX_CLIP; 
            query->seq_start = bg_pos; 
            query->seq_end   = ed_pos; 

            
            query->n_cigar = 0;
            query->query_err = query->l_seq - query->b0;
            sub->query_err   = query->l_seq - query->b0;
        
        }
    } 
    return;
}
int rec_aln_info(query_t *query, uint32_t pos, int score, struct SubBuf *sub){
    int j;
    int CANDI_THRES;
    int delta = sub->delta;
    if(query->b0 > 0) CANDI_THRES = query->b0;
    else CANDI_THRES = query->candi_thres;
    uint32_t i;
  
   
    uint32_t drct_f = query->is_rev*1024;
    aln_out_t *out = sub->aln_out;
    uint32_t (*out_buf)[2] = sub->aln_out->out_buf; 
    int (*found)[4] = sub->aln_out->found; 
    int n_fnd, l_fnd; 
    
    uint32_t pos_i;
    int cur_s;
    //-------------------------- 
    if(score > CANDI_THRES){
        cur_s = score - found[0][0];
        if(cur_s < 0) {
            printf("%u, found[0][0] = %d, cur_s = %d, score = %d\n", 
                    __LINE__, found[0][0], cur_s, score);     
            return; 
        }        
        if(cur_s >= delta) {
            memset(found, 0, 4*delta*sizeof(uint32_t));
        } else if(cur_s < delta) {
            if(cur_s > sub->aln_out->max_found) {
                printf("%u, found[0][0] = %d, cur_s = %d\n", 
                        __LINE__, found[0][0], cur_s);     
                exit(1); 
            } 
            memmove(found+cur_s, found, 4*cur_s*sizeof(uint32_t)); 
            memset(found, 0, 4*cur_s*sizeof(uint32_t));
        } 
        found[0][0] = score;
        found[0][1] = 1;
        found[0][2] = out->len;//begin idx of out_buf
        found[0][3] = out->len;
        
        out_buf[out->len][0] = 1;//number of pos  
        out_buf[out->len][1] = 0;//next pointer of pos 
        out_buf[out->len+1][0] = pos;  
        out_buf[out->len+1][1] = drct_f + score;  
        out->len += 2;
        out->num++;

        query->b0 = score;
        query->pos = pos;
        query->strand = query->is_rev;
        
        query->ref_start = out->L_ti; 
        query->ref_end   = out->R_ti; 
        query->seq_start = out->L_qi; 
        query->seq_end   = out->R_qi; 

        query->n_cigar = 0;
        CANDI_THRES = query->b0;
        query->query_err = query->l_seq - query->b0;
        sub->query_err   = query->l_seq - query->b0;
    } else if(score + delta > CANDI_THRES && found[0][0] >= CANDI_THRES){ 
        int flg = 0; 
        cur_s = found[0][0]-score;
        if(cur_s > sub->aln_out->max_found || cur_s < 0) {
            printf("%u, found[0][0] = %d, cur_s = %d, score = %d\n", 
                    __LINE__, found[0][0], cur_s, score);     
            exit(1); 
        } 
        
        if(found[cur_s][1] == 0){
            found[cur_s][0] = score;
            found[cur_s][1] = 1;
            found[cur_s][2] = out->len;
            found[cur_s][3] = out->len;
            out_buf[out->len][0] = 1;  
            out_buf[out->len][1] = 0; 
            out_buf[out->len+1][0] = pos;  
            out_buf[out->len+1][1] = drct_f + score;  
            out->len += 2;
            out->num++;
        } else{
            found[cur_s][0] = score;
            for(i = 0; i < delta; ++i) {
                if( found[i][3] > found[cur_s][3]) {
                    flg = 1;
                    j = i;
                    break; 
                }
            }
            if(flg == 0) {
                found[cur_s][1]++;
                out_buf[found[cur_s][3]][0]++;  
                out_buf[out->len][0] = pos;  
                out_buf[out->len][1] = drct_f + score;
                out->len++;
                out->num++;
            } else{
                found[cur_s][1]++;
                out_buf[found[cur_s][3]][1] = out->len; 
                found[cur_s][3] = out->len;
                out_buf[out->len][0] = 1;  
                out_buf[out->len][1] = 0; 
                out_buf[out->len+1][0] = pos;  
                out_buf[out->len+1][1] = drct_f + score;  
                out->len += 2;
                out->num++;
            }
        } 
    }// end if(r.score + delta > CANDI_THRES)
    n_fnd  = 0;
    for(i =0; i < delta; ++i) {
        n_fnd += found[i][1];
    }
    return;
}

/*  
int AlgnPos(idx_t *fm_idx, uint8_t ext_cls, query_t *query, int seed_off, uint32_t pos_buf[], int pos_num, struct SubBuf *sub){

   	uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    int read_len = query->l_seq; 
    int MAX_ERROR, ext_err = sub->err_sum[0]*5/3;
    if(query->b0 > 0) MAX_ERROR = read_len - query->b0;
    else MAX_ERROR = query->error_thres;
    MAX_ERROR = MAX_ERROR - ext_err;
    if(MAX_ERROR < 15) MAX_ERROR = 15; 
    int CANDI_THRES;
    if(query->b0 > 0) CANDI_THRES = query->b0;
    else CANDI_THRES = query->candi_thres;

    int flg =0 ;
    int n_aln = 0;
    uint32_t *smpos;
    uint32_t i, j, j_idx, pos_s, smpos_i;
    //uint8_t *target = sub->target;
    uint8_t target[LEN_READ+50] = {};
    uint32_t *aln_pos = sub->aln_r->pos; 
    int aln_num = sub->aln_r->num;
    int16_t (*aln_r)[4] = sub->aln_r->score;//第i个pos右段比对结果，[i][0]是比对分数, [i][1]是左比对分数, [i][2]是右比对分数
    uint16_t algn_r[IS_SMLSIZ+1];
    uint32_t algn_pos[IS_SMLSIZ+1];  
    //kswq_t **swq;
    //swq = sub->swq;
    kswq_t *swq_L[2] = {0, 0};
    kswq_t *swq_R[2] = {0, 0};
    //swq[0] = 0; swq[1] = 0;
    int8_t *mat= sub->mat;       
//++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //[L_st, L_end) [R_st, R_end)
    uint32_t L_st = 0;
    uint32_t L_end = seed_off-(ext_cls*16) - (sub->aln_L[21]-64);//其中10是种子左端点中心的偏移量
    uint32_t ext_len = SEED_LEN+ext_cls*16*2;
    //uint32_t R_st = L_end + ext_len  + (sub->aln_R[21] - 64);
    uint32_t R_st = (seed_off + SEED_LEN )+(ext_cls*16) + sub->aln_R[21]-64;
    uint32_t R_end = read_len;
    uint32_t R_pos_bgn, R_pos_end;
    uint32_t L_pos_bgn, L_pos_end;
    int score;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //以下比对右段序列

    uint8_t *aln_seq; int len_seq; 
    uint32_t pos_i, pos;
    uint16_t MAX_ALN_SIZE = IS_SMLSIZ; 
    aln_seq = read_seq+R_st;
    len_seq = R_end - R_st;
//+++++++++++++++++++++++++++++++++++++++++++

if(L_end-L_st <= R_end-R_st){
    for(pos_i = 0; pos_i < pos_num; ++pos_i){
        pos = pos_buf[pos_i];  
        
        R_pos_bgn = pos+ext_len;
        R_pos_end = R_pos_bgn+len_seq;
        for(i =R_pos_bgn; i < R_pos_end; ++i) 
            target[i-R_pos_bgn] = __get_pac(fm_idx->pac, i);
        kswr_t r = ksw_align_R(len_seq, aln_seq, len_seq, target, 5, 
                                mat, 6, 1, KSW_XSTART, &swq_L[0], sub->kswq_L);

        //+++++++++++++++++++++++++++++++++++
        int score_R = len_seq - r.score;

//fprintf(stderr, "%u, score_R = %u, MAX_ERROR = %u\n",__LINE__,  score_R, MAX_ERROR);
//fprintf(stderr, "R_bg = %u, R_end = %u\n", R_st, R_end);
        if(score_R <= MAX_ERROR){
            L_pos_bgn = pos-L_end;
            uint32_t l_pos = L_pos_bgn > MAX_CLIP?L_pos_bgn-MAX_CLIP:0;
            uint32_t r_pos = L_pos_bgn+read_len+MAX_CLIP;
            for(i =0; i < r_pos - l_pos; ++i) target[i] = __get_pac(fm_idx->pac, l_pos+i);
 
            r = ksw_align_R(read_len, read_seq, r_pos-l_pos, target, 5, 
                mat, 6, 1, KSW_XSTART, &swq_R[0], sub->kswq_R);
            score = r.score;    
//fprintf(stderr, "%u, r.score = %u, CANDI_THRES = %u\n",__LINE__,  r.score, CANDI_THRES);
            if(r.score >= CANDI_THRES) {
                query->b0 = r.score;
                query->pos = L_pos_bgn;
                query->strand = query->is_rev;
                query->tlen = r.te+1-r.tb;
                query->ref_start = r.tb;
                query->seq_start = r.qb;
                query->seq_end = r.qe+1;
                query->n_cigar = 0;
                CANDI_THRES = query->b0;
            }

            if(aln_num > MAX_ALN_SIZE) { 
                int _i, min = aln_r[0][0], row = 1;
                for(_i = 1; _i < aln_num; ++_i) {
                    if(aln_r[_i][0] < min) {
                        min = aln_r[_i][0];
                        row = _i;
                    } 
                }
                if(r.score > min){
                    if(aln_r[0][0] >= r.score ) {
                        aln_r[row][0] = r.score;
                        aln_r[row][1] = score_R;
                        //aln_r[row][2] = algn_r[j];
                        aln_pos[row]  = L_pos_bgn;
                    }else{
                        aln_r[row][0] = aln_r[0][0];
                        aln_r[row][1] = aln_r[0][1];
                        aln_r[row][2] = aln_r[0][2];
                        aln_pos[row]  = aln_pos[0];
                         
                        aln_r[0][0] = r.score;
                        aln_r[0][1] = score_R;
                        //aln_r[0][2] = algn_r[j];
                        aln_pos[0]  = L_pos_bgn;
                    }
                }
                aln_num = MAX_ALN_SIZE;
            } else{
                 
                if(aln_num == 0 || aln_r[0][0] >= r.score ) {
                    aln_r[aln_num][0] = r.score;
                    aln_r[aln_num][1] = score_R;
                    //aln_r[aln_num][2] = algn_r[j];
                    aln_pos[aln_num]  = L_pos_bgn;
                }else{
                    aln_r[aln_num][0] = aln_r[0][0];
                    aln_r[aln_num][1] = aln_r[0][1];
                    aln_r[aln_num][2] = aln_r[0][2];
                    aln_pos[aln_num]  = aln_pos[0];
                    
                    aln_r[0][0] = r.score;
                    aln_r[0][1] = score_R;
                    //aln_r[0][2] = algn_r[j];
                    aln_pos[0]  = L_pos_bgn;
 
                }
              
            }
            aln_num++;
        } else{ }

    }//end for(pos_i = 0; pos_i < pos_num; ++pos_i)++++++
    
} //end if(L_end-L_st < R_end-R_st)++++++++++++++ 
if(L_end- L_st > R_end-R_st){ 
    //===============================================
    //以下代码段是左侧比对
    aln_seq = read_seq;
    len_seq = L_end - L_st;

    MAX_ALN_SIZE = 255;
    swq_L[0] = 0;
    swq_R[0] = 0;
    //for(j = 0; j < n_aln; ++j){
        //pos = algn_pos[j];
    for(pos_i = 0; pos_i < pos_num; ++pos_i){
        pos = pos_buf[pos_i];  
 
        if(pos < L_end) continue;
        L_pos_bgn = pos-L_end;
        L_pos_end = pos;
//fprintf(stderr, "%u, L_pos_bgn = %u\n", __LINE__, L_pos_bgn);
        for(i = L_pos_bgn; i < L_pos_end; ++i) 
            target[i-L_pos_bgn] = __get_pac(fm_idx->pac, i);
        kswr_t r = ksw_align_L(len_seq, aln_seq, len_seq, target, 5, 
            mat, 6, 1, KSW_XSTART, &swq_L[0], sub->kswq_L);
        int score_L = len_seq-r.score;
        score = r.score;

//fprintf(stderr, "%u, score_L = %u, MAX_ERROR = %u\n",__LINE__,  score_L, MAX_ERROR);
        if(score_L <= MAX_ERROR){
            uint32_t l_pos = L_pos_bgn > MAX_CLIP?L_pos_bgn-MAX_CLIP:0;
            uint32_t r_pos = L_pos_bgn+read_len+MAX_CLIP;
            for(i =0; i < r_pos - l_pos; ++i) target[i] = __get_pac(fm_idx->pac, l_pos+i);
            r = ksw_align_R(read_len, read_seq, r_pos-l_pos, target, 5, 
                mat, 6, 1, KSW_XSTART, &swq_R[0], sub->kswq_R);
//fprintf(stderr, "r.score = %u, CANDI_THRES = %u\n", r.score, CANDI_THRES);
            if(r.score > CANDI_THRES) {
                query->b0 = r.score;
                query->pos = L_pos_bgn;
                query->strand = query->is_rev;
                query->tlen = r.te+1-r.tb;
                query->ref_start = r.tb;
                query->seq_start = r.qb;
                query->seq_end = r.qe+1;
                query->n_cigar = 0;


                CANDI_THRES = query->b0;
            }
            if(aln_num > MAX_ALN_SIZE) { 
                int _i, max = aln_r[0][0], row = 1;
                for(_i = 1; _i < aln_num; ++_i) {
                    if(aln_r[_i][0] > max) {
                        max = aln_r[_i][0];
                        row = _i;
                    } 
                }
                if(score < max){
                    if(aln_r[0][0] <= score ) {
                        aln_r[row][0] = score;
                        aln_r[row][1] = score_L;
                        //aln_r[row][2] = algn_r[j];
                        aln_pos[row]  = L_pos_bgn;
                    }else{
                        aln_r[row][0] = aln_r[0][0];
                        aln_r[row][1] = aln_r[0][1];
                        aln_r[row][2] = aln_r[0][2];
                        aln_pos[row]  = aln_pos[0];
                        
                        aln_r[0][0] = score;
                        aln_r[0][1] = score_L;
                        //aln_r[0][2] = algn_r[j];
                        aln_pos[0]  = L_pos_bgn;
                    }
                }
                aln_num = MAX_ALN_SIZE;
            } else{
                 
                if(aln_num == 0 || aln_r[0][0] <= score ) {
                    aln_r[aln_num][0] = score;
                    aln_r[aln_num][1] = score_L;
                    //aln_r[aln_num][2] = algn_r[j];
                    aln_pos[aln_num]  = L_pos_bgn;
                }else{
                    aln_r[aln_num][0] = aln_r[0][0];
                    aln_r[aln_num][1] = aln_r[0][1];
                    aln_r[aln_num][2] = aln_r[0][2];
                    aln_pos[aln_num]  = aln_pos[0];
                    
                    aln_r[0][0] = score;
                    aln_r[0][1] = score_L;
                    //aln_r[0][2] = algn_r[j];
                    aln_pos[0]  = L_pos_bgn;
 
                }
              
            }
            aln_num++;
        } 
    }//end  for(j = 0; j < n_aln; ++j)+++++++++++++
} // end if(L_end - L_st >= R_end-R_st)+++++++++++++
    sub->aln_r->num = aln_num; 
    flg = aln_num;
    return flg;
}
*/

int init_seed_model(seed_t *seed, uint32_t ref_len)
{
    //seed->MAX_SIZE = 18;
    int SEED_SLC_SIZE = seed->MAX_SIZE;
    seed->slc   = calloc(SEED_SLC_SIZE, sizeof(seed_slc_t)); 
    if(seed->slc == NULL) {
        printf("seed->slc == NULL!!!\n"); 
        exit(1);
    }

    //seed->slc_i = calloc(SEED_SLC_SIZE, sizeof(seed_slc_t)); 
    seed->err   = calloc(SEED_SLC_SIZE, sizeof(struct ext_err_t)); 
    if(seed->err == NULL) {
        printf("seed->err == NULL!!!\n"); 
        exit(1);
    }
    
    //seed->pos_tmp = calloc(MAX_SEED_NUM, sizeof(uint32_t));
    //seed->m_pos   = calloc(MAX_SEED_NUM*2, sizeof(uint32_t));
    //seed->m_flg   = calloc(ref_len/64+1, sizeof(uint8_t));//64 = 8*8, 8个间隔存储数据
    //uint32_t *seq_pos = calloc(seed->MAX_SIZE, sizeof(uint32_t));
    //seed->L_len = seq_pos; 
    /*  
    int i;
    for(i = 0; i < seed->MAX_SIZE/2;  ++i) {
        seed->slc_i[i].drct   = 0;
        seed->slc_i[i].s_off  = 0;
    }
    for(i = seed->MAX_SIZE/2; i < seed->MAX_SIZE;  ++i) {
        seed->slc_i[i].drct   = 1;
        seed->slc_i[i].s_off  = 0;
    }
    */
    return;
}
/*  
uint32_t aln_seed_seq(idx_t *fm_idx,uint32_t hash_boundry[], int read_len,  uint8_t *f_read_seq,  uint8_t *r_read_seq, uint8_t *qual, uint8_t* r_qual, seed_t *seed )
{
            int i, j;
            int seed_slc_size = seed->slc_size;
         
            //确定第seed_id新的种子序列           
            seed_slc_t *slc   = seed->slc; 
            seed_slc_t *slc_i = seed->slc_i; 
            seed_slc_t *slc_b = seed->slc_b; 
            aln_mer8_t *aln_8 = seed->aln_8;
            int find_num = seed->find_num;
            int cur_row  = seed->cur_row; 
            int seed_id  = seed->id;
            uint32_t seq12;

            uint8_t *pseed, *seed_buf= seed->seed_buf; 
            int seq_off = read_len/2-slc_i[seed_id].s_off;

            int seed_flg = seed->seed_flg;
            uint32_t bgn, end, num = 0, k, l;


            int MIN_Q = 0;
            if(seed_id == 0) {
                 //1）获取read信息,获取read反向互补序列；
                //2）生成正向序列与反向序列上的种子序列；
                int row_b = 0, row = 0;
                for(i = 0; i < seed->MAX_SIZE; ++i){
                    uint32_t seed_bg, seed_ed;
                    seed_bg = read_len/2-(slc_i[i].s_off-slc_i[i].h_off);
                    seed_ed = seed_bg+12;
                    int sum = 0;
                    if(slc_i[i].drct == 0) {
                        for(j = seed_bg; j < seed_ed; ++j) {
                             if( qual[j] < MIN_Q ){ 
                                 //sum++; 
                             }
                         }
                    } else{
                        for(j = seed_bg; j < seed_ed; ++j) {
                            if( r_qual[j] < MIN_Q ){
                                 //sum++; 
                            }
                        }
                    }
                    if(sum > 2) { continue;}
                    else if(sum > 0){ 
                        memcpy(slc_b+row_b, slc_i+i, sizeof(seed_slc_t));
                        row_b++;
                    } else{
                        memcpy(slc+row, slc_i+i, sizeof(seed_slc_t));
                        row++;
                    }
                } 
                memcpy(slc+row, slc_b, row_b*sizeof(seed_slc_t));
                //3）赋值当前read_id的seed_slc_size;
                seed_slc_size = row + row_b;
                seed->slc_size = seed_slc_size;
                find_num = 0; 
            } //end if(seed_id == 0)  
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            //种子序列指针赋值
            if(find_num==0) {
                if(slc[seed_id].drct == 0){
                    pseed = f_read_seq+seq_off;
                } else{//slc[seed_id].drct ==1
                    pseed = r_read_seq+seq_off;
                }     
            }             
        while(seed_id < seed_slc_size){ 
            if(find_num==0) {
                seq_off = read_len/2-slc_i[seed_id].s_off;
                if(slc[seed_id].drct == 0){
                    pseed = f_read_seq+seq_off;
                } else{//slc[seed_id].drct ==1
                    pseed = r_read_seq+seq_off;
                }
            }   
            if(find_num == 0) {//选定新的种子序列  
                if(slc[seed_id].flag == 0) {//通过BWT精确比对选定种子序列
                    seq12 = lkt_seq2LktItem(pseed, 8, 19);
                    k = fm_idx->fastmap->item[seq12];
                    l = fm_idx->fastmap->item[seq12+1]-1;
                    l -= get_12mer_correct(hash_boundry, l);
                    num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, pseed, &k, &l);
                    bgn = k;
                    end = l;
                    seed_flg = 0;
                } else if(slc[seed_id].flag ==1){//通过BWT选定近似种子序列集合
                    cur_row = 0;
                    seed->hash_12mer = lkt_seq2LktItem(pseed, 8, 19);
                    seed->seed_8mer_L = pseed;
                    seed->seed_seq = pseed;
                    find_num = get_seed_8mer_bwt_L(fm_idx, hash_boundry, seed);
                    if(find_num == 0) {
                        seed_id++;
                        continue;
                    }
         

                    bgn = seed->aln_buf[cur_row][2]; 
                    end = seed->aln_buf[cur_row][3]; 
                    num = end+1- bgn;
                    seed_flg = 1;
                    if(find_num == 1) { 
                        find_num = 0; 
                    }
                } else{//slc[seed_id] ==2
                    //调用右侧模糊比对查找函数,获得模糊比对序列集合数组
                    int bg_i = 0, ed_i = 8, num_seq;
                    num_seq = (ed_i-bg_i)*3; 
                    seed->hash_12mer = lkt_seq2LktItem(pseed, 0, 11);
                    seed->seed_8mer_R = pseed+12;
                    seed->seed_seq = pseed;
                    int var_seq_num = gen_8mer_var_buf(seed, bg_i, ed_i, 1); 
                    find_num = get_seed_8mer_seq_R(fm_idx, seed, num_seq);
                    if(find_num == 0) {
                        seed_id++;
                        continue;
                    }
                    memcpy(seed_buf, pseed, SEED_LEN); 
                    cur_row = 0;
                    uint8_t f_i = seed->find_buf_R[cur_row];  
                    uint8_t var_pos  = seed->var_info[f_i][0];
                    uint8_t var_ch   = seed->var_info[f_i][1];
                    uint8_t var_ch_o = seed_buf[var_pos+12]; 
                    seed_buf[var_pos+12] = var_ch; 
                    seq12 = lkt_seq2LktItem(seed_buf, 8, 19);
                    k = fm_idx->fastmap->item[seq12];
                    l = fm_idx->fastmap->item[seq12+1]-1;
                    l -= get_12mer_correct(hash_boundry, l);
                    num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, seed_buf, &k, &l);
                    bgn = k; end = l;
                    seed_flg = 2;
                    if(find_num == 1) { 
                        find_num = 0; 
                    }
                    seed_buf[var_pos+12] = var_ch_o; 
                }
            } else{//if(find_num >0)
                if(seed_flg ==1) {//第一类近似比对，即BWT近似比对, 种子序列继续
                    seed_id--;
                    cur_row++; 
                    bgn = seed->aln_buf[cur_row][2]; 
                    end = seed->aln_buf[cur_row][3]; 
                    num = end+1- bgn;
                    //seed_flg = 1;
                } else if(seed_flg ==2) {//第二类近似比对继续
                    seed_id--;
                    cur_row++; 
                    
                    uint8_t f_i = seed->find_buf_R[cur_row];  
                    uint8_t var_pos  = seed->var_info[f_i][0];
                    uint8_t var_ch   = seed->var_info[f_i][1];
                    uint8_t var_ch_o = seed_buf[var_pos+12]; 
                    seed_buf[var_pos+12] = var_ch; 
                    seq12 = lkt_seq2LktItem(seed_buf, 8, 19);
                    k = fm_idx->fastmap->item[seq12];
                    l = fm_idx->fastmap->item[seq12+1]-1;
                    l -= get_12mer_correct(hash_boundry, l);
                    num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, seed_buf, &k, &l);
                    bgn = k; 
                    end = l;

                    seed_buf[var_pos+12] = var_ch_o; 
                }

                if(find_num - cur_row == 1) {//cur_row最后一个行的情况 
                    find_num = 0; 
                }
            }

            seed->find_num = find_num ;
            seed->cur_row = cur_row ; 
            seed->id = seed_id ;
            seed->seed_flg = seed_flg;
            seed->bgn = bgn;
            seed->end = end;
            seed->num = num;

            if(num >0){//种子序列比对失败
                break; 
            }
            seed->id = ++seed_id;
        }
        
        seed->id = seed_id ;
        return num;
    }
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*  
uint32_t aln_seed_seq_1(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed )
{
    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq;
    uint8_t *qual = query->qual;
//uint8_t r_qual[LEN_READ] = {};
    uint8_t r_qual[LEN_READ] = {};

    int i, j;
    int seed_slc_size = seed->slc_size;
 
    //确定第seed_id新的种子序列           
    seed_slc_t *slc   = seed->slc; 
    seed_slc_t *slc_i = seed->slc_i; 
    seed_slc_t *slc_b = seed->slc_b; 
    aln_mer8_t *aln_8 = seed->aln_8;
    int find_num = seed->find_num;
    int cur_row  = seed->cur_row; 
    int seed_id  = seed->id;
    uint32_t seq12;

    uint8_t *pseed, *seed_buf= seed->seed_buf; 
    //int seq_off = read_len/2-slc_i[seed_id].s_off;
    int seq_off[2];
    seq_off[0] = query->center_off[0];
    seq_off[1] = query->center_off[1];

    int seed_flg = seed->seed_flg;
    uint32_t bgn, end, num = 0, k, l;


    int MIN_Q = 0;
    if(seed_id == 0) {
         //1）获取read信息,获取read反向互补序列；
        //2）生成正向序列与反向序列上的种子序列；
        int row_b = 0, row = 0;
        
        for(i = 0; i < seed->MAX_SIZE; ++i){
            uint32_t seed_bg, seed_ed;
            //seed_bg = read_len/2-(slc_i[i].s_off-slc_i[i].h_off);
            seed_bg = seq_off[slc_i[i].drct] -(slc_i[i].s_off-slc_i[i].h_off);
            seed_ed = seed_bg+12;
            int sum = 0;
            if(slc_i[i].drct == 0) {
                for(j = seed_bg; j < seed_ed; ++j) {
                     if( qual[j] < MIN_Q ){ 
                         //sum++; 
                     }
                 }
            } else{
                for(j = seed_bg; j < seed_ed; ++j) {
                    if( r_qual[j] < MIN_Q ){
                         //sum++; 
                    }
                }
            }
            if(sum > 2) { continue;}
            else if(sum > 0){ 
                memcpy(slc_b+row_b, slc_i+i, sizeof(seed_slc_t));
                row_b++;
            } else{
                memcpy(slc+row, slc_i+i, sizeof(seed_slc_t));
                row++;
            }
        } 
        memcpy(slc+row, slc_b, row_b*sizeof(seed_slc_t));
        //3）赋值当前read_id的seed_slc_size;
        seed_slc_size = row + row_b;
        seed->slc_size = seed_slc_size;
        find_num = 0; 
    } //end if(seed_id == 0)  
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //种子序列指针赋值
    if(find_num==0) {
        if(slc[seed_id].drct == 0){
            //pseed = f_read_seq+seq_off[0] - slc[seed_id].s_off;
            pseed = f_read_seq+seq_off[0] + slc[seed_id].s_off;
//fprintf(stderr, "%u, seed_id = %u, L_offset = %u\n", __LINE__, seed_id, seq_off[0]+slc[seed_id].s_off); 
        } else{//slc[seed_id].drct ==1
            //pseed = r_read_seq+seq_off[1] - slc[seed_id].s_off;
            pseed = r_read_seq+seq_off[1] + slc[seed_id].s_off;
//fprintf(stderr, "%u, seed_id = %u, L_offset = %u\n", __LINE__, seed_id, seq_off[1]+slc[seed_id].s_off); 
        }     
    }             
    
    while(seed_id < seed_slc_size){ 
        if(find_num==0) {
            //seq_off = read_len/2-slc_i[seed_id].s_off;
            if(slc[seed_id].drct == 0){
                //pseed = f_read_seq+seq_off[0] - slc[seed_id].s_off;
                pseed = f_read_seq+seq_off[0] + slc[seed_id].s_off;

//fprintf(stderr, "%u, seed_id = %u, L_offset = %u\n", __LINE__, seed_id, seq_off[0]+slc[seed_id].s_off); 
            } else{//slc[seed_id].drct ==1
                //pseed = r_read_seq+seq_off[1] - slc[seed_id].s_off;
                pseed = r_read_seq+seq_off[1] + slc[seed_id].s_off;
//fprintf(stderr, "%u, seed_id = %u, L_offset = %u\n", __LINE__, seed_id, seq_off[1]+slc[seed_id].s_off); 
            }
        }   
        if(find_num == 0) {//选定新的种子序列  
            if(slc[seed_id].flag == 0) {//通过BWT精确比对选定种子序列
                seq12 = lkt_seq2LktItem(pseed, 8, 19);
                k = fm_idx->fastmap->item[seq12];
                l = fm_idx->fastmap->item[seq12+1]-1;
                l -= get_12mer_correct(hash_boundry, l);
                num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, pseed, &k, &l);
                bgn = k;
                end = l;
                seed_flg = 0;
            } else if(slc[seed_id].flag ==1){//通过BWT选定近似种子序列集合
                cur_row = 0;
                seed->hash_12mer = lkt_seq2LktItem(pseed, 8, 19);
                seed->seed_8mer_L = pseed;
                seed->seed_seq = pseed;
                find_num = get_seed_8mer_bwt_L(fm_idx, hash_boundry, seed);
                if(find_num == 0) {
                    seed_id++;
                    continue;
                }
//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
for(i = 0; i < find_num; ++i){
//fprintf(stderr, "err_pos = %u,alt=%u,  bgn=%u, end = %u\n", seed->aln_buf[i][0], seed->aln_buf[i][1], seed->aln_buf[i][2],                     seed->aln_buf[i][3]); 
}
                bgn = seed->aln_buf[cur_row][2]; 
                end = seed->aln_buf[cur_row][3]; 
                num = end+1- bgn;
                seed_flg = 1;
                if(find_num == 1) { 
                    find_num = 0; 
                }
            } else{//slc[seed_id] ==2
                //调用右侧模糊比对查找函数,获得模糊比对序列集合数组
                int bg_i = 0, ed_i = 8, num_seq;
                num_seq = (ed_i-bg_i)*3; 
                seed->hash_12mer = lkt_seq2LktItem(pseed, 0, 11);
                seed->seed_8mer_R = pseed+12;
                seed->seed_seq = pseed;
                int var_seq_num = gen_8mer_var_buf(seed, bg_i, ed_i, 1); 
                find_num = get_seed_8mer_seq_R(fm_idx, seed, num_seq);
                if(find_num == 0) {
                    seed_id++;
                    continue;
                }
                memcpy(seed_buf, pseed, SEED_LEN); 
                cur_row = 0;
                uint8_t f_i = seed->find_buf_R[cur_row];  
                uint8_t var_pos  = seed->var_info[f_i][0];
                uint8_t var_ch   = seed->var_info[f_i][1];
                uint8_t var_ch_o = seed_buf[var_pos+12]; 
                seed_buf[var_pos+12] = var_ch; 
                seq12 = lkt_seq2LktItem(seed_buf, 8, 19);
                k = fm_idx->fastmap->item[seq12];
                l = fm_idx->fastmap->item[seq12+1]-1;
                l -= get_12mer_correct(hash_boundry, l);
                num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, seed_buf, &k, &l);
                bgn = k; end = l;
                seed_flg = 2;
                if(find_num == 1) { 
                    find_num = 0; 
                }
                seed_buf[var_pos+12] = var_ch_o; 
            }
        } else{//if(find_num >0)
            if(seed_flg ==1) {//第一类近似比对，即BWT近似比对, 种子序列继续
                seed_id--;
                cur_row++; 
                bgn = seed->aln_buf[cur_row][2]; 
                end = seed->aln_buf[cur_row][3]; 
                num = end+1- bgn;
                //seed_flg = 1;
            } else if(seed_flg ==2) {//第二类近似比对继续
                seed_id--;
                cur_row++; 
                
                uint8_t f_i = seed->find_buf_R[cur_row];  
                uint8_t var_pos  = seed->var_info[f_i][0];
                uint8_t var_ch   = seed->var_info[f_i][1];
                uint8_t var_ch_o = seed_buf[var_pos+12]; 
                seed_buf[var_pos+12] = var_ch; 
                seq12 = lkt_seq2LktItem(seed_buf, 8, 19);
                k = fm_idx->fastmap->item[seq12];
                l = fm_idx->fastmap->item[seq12+1]-1;
                l -= get_12mer_correct(hash_boundry, l);
                num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, seed_buf, &k, &l);
                bgn = k; 
                end = l;
                seed_buf[var_pos+12] = var_ch_o; 
            }

            if(find_num - cur_row == 1) {//cur_row最后一个行的情况 
                find_num = 0; 
            }
        }

        seed->find_num = find_num ;
        seed->cur_row = cur_row ; 
        seed->id = seed_id ;
        seed->seed_flg = seed_flg;
        seed->bgn = bgn;
        seed->end = end;
        seed->num = num;
        if(num >0){//种子序列比对失败
            break; 
        }
        seed->id = ++seed_id;
    }
    seed->id = seed_id;
    return num;
}
*/
/*  
int get_seed_seq(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed)
{
    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq;

    seed_slc_t *slc = seed->slc; 
    int seed_id  = seed->id;
    
    uint32_t seq12;
    uint8_t *pseed; 
    uint32_t k, l, bgn, end, num = 0;
    int i, j;
    //------------------------------------------------------------ 
    if(slc[seed_id].drct == 0){
        pseed = f_read_seq+query->m_off[0] + slc[seed_id].s_off;
    } else{//slc[seed_id].drct ==1
        pseed = r_read_seq+query->m_off[1] + slc[seed_id].s_off;
    }

    seq12 = lkt_seq2LktItem(pseed, 8, 19);
    k = fm_idx->fastmap->item[seq12];
    l = fm_idx->fastmap->item[seq12+1]-1;
    l -= get_12mer_correct(hash_boundry, l);
    num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, pseed, &k, &l);
    bgn = k;
    end = l;
    seed->id = seed_id;
    seed->bgn = bgn;
    seed->end = end;
    seed->num = num;
    if(num >0){//种子序列产生成功
        slc[seed_id].bgn = bgn; 
        slc[seed_id].end = end; 
        slc[seed_id].num = num; 
    } 
    return num;
}
*/

void init_seed_pos_r(query_t *query, seed_t *seed, struct StackTree *sTree)
{

    int read_len = query->l_seq;
    int i, j, len;
    seed_slc_t *slc   = seed->slc; 
    //确定第seed_id新的种子序列           
   
    int sn = query->seed_num;  
    int l_seed = query->l_seq/sn;   
    //int m_off = read_len/3 + l_seed - (read_len/3)%l_seed;
    int m_off = read_len/3;
    int s_off = m_off - SEED_LEN/2;
    s_off += l_seed - s_off%l_seed; 

    int l_len, r_len, num;
    
    //slc[0].s_off  = m_off - SEED_LEN/2 ;
    slc[0].s_off  = s_off;

    slc[0].ext_num = slc[0].s_off/16;
    if(slc[0].ext_num > query->max_ext) 
        slc[0].ext_num = query->max_ext;


    l_len = slc[0].s_off;
    r_len = slc[0].s_off;

    //for(i = 1; i < query->seed_num; ++i) {
    i = 1;
    int seed_num = 1;

    //while(1){
    for(i = 1; i < sn; ++i) {

        if(i % 2 == 1) {
            if(l_len < l_seed) break; 
            l_len -= l_seed;
            slc[i].s_off = l_len; 
            slc[i].ext_num = slc[i].s_off/16;  //17 = 16+1
        } else{
            r_len += l_seed;
            //if(r_len > read_len*2/3) break; 
            slc[i].s_off = r_len;
            slc[i].ext_num = (2*m_off - slc[i].s_off)/16;
        }
fprintf(stderr, "%u, seed_id = %d, l_len = %d, r_len = %d, s_off = %d\n", __LINE__, i, l_len, r_len, slc[i].s_off);
        if(slc[i].ext_num > query->max_ext) 
            slc[i].ext_num = query->max_ext;
        seed_num++;
    }
    for(i = seed_num; i < sn; ++i){
        r_len += l_seed;
        if(r_len  > read_len - SEED_LEN) {
            printf("%s, %u, r_len = %d\n", __func__, __LINE__, r_len);
            exit(1); 
        }

        slc[i].s_off = r_len;
fprintf(stderr, "%u, seed_id = %d, l_len = %d, r_len = %d, s_off = %d\n", __LINE__, i, l_len, r_len, slc[i].s_off);
        if(2*m_off > slc[i].s_off) {
            slc[i].ext_num = (2*m_off - slc[i].s_off)/16; 
        } else {
            slc[i].ext_num = 0;
        }
        if(slc[i].ext_num > query->max_ext) 
            slc[i].ext_num = query->max_ext;
        //seed_num++;
    }
    //m_off = read_len*2/3;
    s_off = m_off - SEED_LEN;
    s_off = read_len - (s_off - s_off%l_seed + SEED_LEN);


    //slc[sn].s_off  = 2*m_off - SEED_LEN/2 ;
    slc[sn].s_off  = s_off;
    slc[sn].ext_num = (read_len - slc[sn].s_off - SEED_LEN)/16;
    if(slc[sn].ext_num > query->max_ext) 
        slc[sn].ext_num = query->max_ext;

    l_len = slc[sn].s_off; 
    r_len = l_len;
    seed_num = sn+1;
    for(i = sn + 1; i < sn*2; ++i) {

fprintf(stderr, "seed_id = %d, l_len = %d, r_len = %d\n", i, l_len, r_len);
        if(i % 2 == 1) {
            if(r_len + l_seed > read_len - SEED_LEN) break; 
            r_len += l_seed;
            slc[i].s_off = r_len;
            slc[i].ext_num = (read_len - slc[i].s_off-SEED_LEN)/16;
        } else{
            l_len -= l_seed;
            slc[i].s_off = l_len;     
            slc[i].ext_num = (slc[i].s_off - m_off + SEED_LEN)/16;
        }
        if(slc[i].ext_num > query->max_ext) 
            slc[i].ext_num = query->max_ext;

        seed_num++;
    }

    for(i = seed_num; i < 2*sn; ++i){
        l_len -= l_seed;
fprintf(stderr, "seed_id = %d, l_len = %d, r_len = %d, slc[]\n", seed_num, l_len, r_len);
        if(l_len < 0) {
            printf("%s, %u, l_len = %d\n", __func__, __LINE__, l_len);
            exit(1); 
        }
        slc[i].s_off = l_len;
        if(slc[i].s_off < m_off) {
            slc[i].ext_num = (slc[i].s_off - m_off + SEED_LEN)/16; 
        } else {
            slc[i].ext_num = 0;
        }
        if(slc[i].ext_num > query->max_ext) 
            slc[i].ext_num = query->max_ext;
        //seed_num++;
    }

   
    len = read_len/2-SEED_LEN/2;
    for(i = 0; i < sn*2; ++i) { 

fprintf(stderr, "s_off[%d] = %d\n", i, slc[i].s_off);
        slc[i].bgn = 0;
        slc[i].num = 0;
        slc[i].end = 0;
    }
    /*  
    for(i = 0; i < seed->m_num; ++i){
        seed->m_flg[seed->m_pos[i]/8] = 0;
    }
    */
    slc[0].sTree_bg  = 0;  
    slc[0].sTree_num = 0;
    sTree->bg_buf    = 0;
    sTree->len_buf   = 0;

    return;
}
void init_seed_pos1(query_t *query, seed_t *seed, struct StackTree *sTree)
{
    int buf[LEN_READ/16][2];
    int sort[LEN_READ/16];
     
    int read_len = query->l_seq;
    int i, j, len;
    seed_slc_t *slc   = seed->slc; 
    //确定第seed_id新的种子序列           

    int l_seed, l_len, r_len, num;
    int sn0 = query->seed_add;
    int sn1 = query->seed_num - sn0;
    int sn2 = 2*(sn0 + sn1);
    //slc[0].s_off  = query->m_off - SEED_LEN/2 ;
    int flg;
    int m_off = query->m_off;
    l_seed = read_len/sn1;
    int l_bg = (read_len%l_seed)/2;
    for(i = 0; i < sn1; ++i) {
        buf[i][0] = slc[i].s_off;
        buf[i][1] = slc[i].ext_num;
    }
    
    l_len = (m_off - l_bg)/l_seed;
    r_len = l_len;

    if(m_off <= read_len/2) {
        flg = 0; 
    } else {
        flg = 1;
    }  
    sort[0] = l_len; 
//fprintf(stderr, "l_len = %d, i = 0\n", l_len);
    for(i = 1; i < sn1; ++i) {
//fprintf(stderr, "l_len = %d, r_len = %d, i = %d\n", l_len, r_len, i);
        if(flg > 0) {
            if(i%2 == 1) {
                sort[i] = --l_len;
                //break; 
            } else {
                sort[i] = ++r_len; 
            } 
        } else {
            if(i%2 == 1) {
                sort[i] = ++r_len; 
            } else {
                sort[i] = --l_len;
                //break; 
            }
        } 
        if(l_len < 0 || r_len*l_seed + SEED_LEN > read_len) {
            break; 
        }
    }
    int ni = i;
    for(i = 0; i < ni; ++i) {
        slc[i].s_off = buf[sort[i]][0];
        slc[i].ext_num = buf[sort[i]][1];
    }
/*  
    for(i = 0; i < sn1; ++i) {
        fprintf(stderr, "%u, sort[%d] = %d\n", __LINE__, i, sort[i]); 
    }
*/
    
    m_off = read_len - m_off;   
    //for(i = sn1; i < sn2; ++i) {
    for(i = sn0+sn1; i < sn2-sn0; ++i) {
        buf[i-(sn0+sn1)][0] = slc[i].s_off;
        buf[i-(sn0+sn1)][1] = slc[i].ext_num;
    }
    
    l_len = (m_off - l_bg)/l_seed;
    r_len = l_len;
    
    if(read_len - m_off <= read_len/2) {
        flg = 0; 
    } else {
        flg = 1;
    }  

    
    sort[0] = l_len; 
//fprintf(stderr, "l_len = %d, i = 0\n", l_len);
    for(i = 1; i < sn1; ++i) {
//fprintf(stderr, "l_len = %d, r_len = %d, i = %d\n", l_len, r_len, i);
        if(flg > 0) {
            if(i%2 == 1) {
                sort[i] = --l_len;
                //break; 
            } else {
                sort[i] = ++r_len; 
            } 
        } else {
            if(i%2 == 1) {
                sort[i] = ++r_len; 
            } else {
                sort[i] = --l_len;
                //break; 
            }
        } 
        if(l_len < 0 || r_len*l_seed + SEED_LEN > read_len) {
            break; 
        }
    }
/*  
    for(i = 0; i < sn1; ++i) {
        fprintf(stderr, "%u, sort[%d] = %d\n", __LINE__, i, sort[i]); 
    }
*/
    ni = i;
    for(i = ni; i < sn1; ++i) {
        sort[i] = sn1 - i - 1;
    }
    for(i = 0; i < sn1; ++i) {
        slc[i+sn1].s_off = buf[sort[i]][0];
        slc[i+sn1].ext_num = buf[sort[i]][1];
    }
    /*  
    for(i = ni; i < sn1; ++i) {
        //slc[sn2-i].s_off = buf[i][0];
        slc[i+sn1].s_off = buf[i][0];
        //slc[sn2 -i].ext_num = buf[i][1];
        slc[i+sn1].ext_num = buf[i][1];
    }
    */
    int flg_val;
    if(flg > 0) {
        flg_val = 1;
    } else {
        flg_val = -1; 
    }
    
    int sf, delta = 0;
    int bg, ed;
    bg = 0, ed = sn1;
    j = 0;
    for(i = bg; i < ed; ++i){
        sf =  slc[i].s_off + flg_val * SEED_LEN/2;
        if(sf > 0 && sf + SEED_LEN < read_len) {
            slc[ed+j].s_off = sf; 
            if(slc[i].s_off <= sf) {
                if(slc[i].ext_num*(16+delta) + sf+ SEED_LEN > read_len) {
                    slc[ed+j].ext_num = (read_len - SEED_LEN - sf)/(16+delta); 
                } else {
                    slc[ed+j].ext_num = slc[i].ext_num; 
                } 
            } else {
                if(sf - slc[i].ext_num*(16+delta) < 0) {
                    slc[ed+j].ext_num = sf/(16+delta); 
                } else {
                    slc[ed+j].ext_num = slc[i].ext_num; 
                } 
            } 
            ++j;
            if(j >= sn0) break;
        }
                
    } 
    j = 0;
    bg = sn0 + sn1, ed = sn2 - sn0;
    for(i = bg; i < ed; ++i){
        sf =  slc[i].s_off - flg_val * SEED_LEN/2;
        if(sf > 0 && sf + SEED_LEN < read_len) {
            slc[ed+j].s_off = sf; 
            if(slc[i].s_off <= sf) {
                if(slc[i].ext_num*(16+delta) + sf+ SEED_LEN > read_len) {
                    slc[ed+j].ext_num = (read_len - SEED_LEN - sf)/(16+delta); 
                } else {
                    slc[ed+j].ext_num = slc[i].ext_num; 
                } 
            } else {
                if(sf - slc[i].ext_num*(16+delta) < 0) {
                    slc[ed+j].ext_num = sf/(16+delta); 
                } else {
                    slc[ed+j].ext_num = slc[i].ext_num; 
                } 
            } 
            ++j;
            if(j >= sn0) break;
        }


    }

fprintf(stderr, "%u  -----------\n\n", __LINE__);
    for(i = 0; i < sn2; ++i) { 
        slc[i].bgn = 0;
        slc[i].num = 0;
        slc[i].end = 0;
      
fprintf(stderr, "%u, slc[%u].s_off = %d, ext_num = %d, l_len = %d, r_len = %d\n", __LINE__, i, slc[i].s_off, slc[i].ext_num, slc[i].s_off - 16*slc[i].ext_num, query->l_seq - (slc[i].s_off + SEED_LEN + 16*slc[i].ext_num));
    }
fprintf(stderr, "%u  -----------\n\n", __LINE__);
    return;
}

void init_seed_pos(query_t *query, seed_t *seed, struct StackTree *sTree)
{
    int read_len = query->l_seq;
    int i, j, len;
    seed_slc_t *slc   = seed->slc; 
    int l_seed, l_len, r_len, num;
    int sn0 = query->seed_add;
    int sn1 = query->seed_num - sn0;
    int sn2 = 2*(sn0 + sn1);
    l_seed = read_len/sn1;
    int l_bg = (read_len - l_seed*sn1 + l_seed - SEED_LEN)/2;
    if(l_bg < 0) l_bg = 0;
    l_len = l_bg; 
    for(i = 0; i < sn1; ++i) {
        slc[i].s_off = l_len; 
        l_len += l_seed;
    }
    if(slc[sn1-1].s_off > read_len - SEED_LEN) {
        slc[sn1-1].s_off = read_len - SEED_LEN; 
    } 
    
     
    
    int buf[LEN_READ/16][2];
    int sort[LEN_READ/16];
    //确定第seed_id新的种子序列           
    int flg;
    int m_off = query->m_off;
    l_seed = read_len/sn1;
    for(i = 0; i < sn1; ++i) {
        buf[i][0] = slc[i].s_off;
    }
    
    l_len = (m_off - l_bg)/l_seed;
    r_len = l_len;

    if(m_off <= read_len/2) {
        flg = 0; 
    } else {
        flg = 1;
    }  
    sort[0] = l_len; 
    for(i = 1; i < sn1; ++i) {
        if(flg > 0) {
            if(i%2 == 1) {
                sort[i] = --l_len;
                //break; 
            } else {
                sort[i] = ++r_len; 
            } 
        } else {
            if(i%2 == 1) {
                sort[i] = ++r_len; 
            } else {
                sort[i] = --l_len;
                //break; 
            }
        } 
        if(l_len < 0 || r_len*l_seed + SEED_LEN > read_len) {
            break; 
        }
    }
    int ni = i;
    for(i = 0; i < ni; ++i) {
        slc[i].s_off = buf[sort[i]][0];
    }
    /*  
    m_off = read_len - m_off;   
    for(i = sn0+sn1; i < sn2-sn0; ++i) {
        buf[i-(sn0+sn1)][0] = slc[i].s_off;
    }
    
    l_len = (m_off - l_bg)/l_seed;
    r_len = l_len;
    
    if(read_len - m_off <= read_len/2) {
        flg = 0; 
    } else {
        flg = 1;
    }  
    sort[0] = l_len; 
    for(i = 1; i < sn1; ++i) {
        if(flg > 0) {
            if(i%2 == 1) {
                sort[i] = --l_len;
                //break; 
            } else {
                sort[i] = ++r_len; 
            } 
        } else {
            if(i%2 == 1) {
                sort[i] = ++r_len; 
            } else {
                sort[i] = --l_len;
                //break; 
            }
        } 
        if(l_len < 0 || r_len*l_seed + SEED_LEN > read_len) {
            break; 
        }
    }


    ni = i;
    for(i = ni; i < sn1; ++i) {
        sort[i] = sn1 - i - 1;
    }
    for(i = 0; i < sn1; ++i) {
        slc[i+sn1].s_off = buf[sort[i]][0];
    }

    int flg_val;
    if(flg > 0) {
        flg_val = 1;
    } else {
        flg_val = -1; 
    }
*/    
    int sf, delta = 0;
    int bg, ed, mid_off;
    bg = sn1, ed = sn1 + sn0;
    j = 0;
    mid_off = slc[0].s_off + SEED_LEN/2;
    int sign = 1;
    int l_d = 0;
    for(i = 0; i < sn0; ++i){
        slc[sn1 + i].s_off = mid_off + sign*l_d; 
        if(i%2 == 0) l_d += l_seed;     
        sign *= -1;
    }
    
    for(i = 0; i < sn1+sn0; ++i){
        int l_off = query->l_seq - (slc[i].s_off + SEED_LEN);
        slc[sn0+sn1+i].s_off = l_off;
    } 
    
    
    
    for(i =0 ; i < sn2; ++i) {
        int l_off = slc[i].s_off; 
        int r_off = read_len - slc[i].s_off - SEED_LEN; 
        int min_off = l_off < r_off?l_off:r_off; 
        slc[i].ext_num = min_off/16; 
        if(slc[i].ext_num > query->max_ext) slc[i].ext_num = query->max_ext;
    }


    
    

fprintf(stderr, "%u  -----------\n\n", __LINE__);
    for(i = 0; i < sn2; ++i) { 
        slc[i].bgn = 0;
        slc[i].num = 0;
        slc[i].end = 0;

fprintf(stderr, "%u, slc[%u].s_off = %d, ext_num = %d, l_len = %d, r_len = %d\n", __LINE__, i, slc[i].s_off, slc[i].ext_num, slc[i].s_off - 16*slc[i].ext_num, query->l_seq - (slc[i].s_off + SEED_LEN + 16*slc[i].ext_num));

    }
fprintf(stderr, "%u  -----------\n\n", __LINE__);
    return;
}

void init_seed_pos0(query_t *query, seed_t *seed, struct StackTree *sTree)
{
//printf("max_ext = %d\n", query->max_ext);
    int read_len = query->l_seq;
    int i, j, len;
    seed_slc_t *slc   = seed->slc; 
    //确定第seed_id新的种子序列           

    int l_seed, l_len, r_len, num;
    int sn0 = query->seed_add;
    int sn1 = query->seed_num - sn0;
    int sn2 = 2*(sn0 + sn1);
    //slc[0].s_off  = query->m_off - SEED_LEN/2 ;
    l_seed = read_len/sn1;
    //l_len = slc[0].s_off;
    //r_len = slc[0].s_off;  

    //int l_bg = (read_len%l_seed)/2;
    int l_bg = (read_len - l_seed*sn1 + l_seed - SEED_LEN)/2;
    l_len = l_bg; 
    for(i = 0; i < sn1; ++i) {
        
        slc[i].s_off = l_len; 
        if(l_len < read_len - SEED_LEN - l_len) {
            len = l_len; 
        } else {
            len = read_len - SEED_LEN - l_len; 
        }
        slc[i].ext_num = len/16;  //17 = 16+1
        if(slc[i].ext_num > query->max_ext) 
            slc[i].ext_num = query->max_ext; 
        l_len += l_seed;
    }
    if(slc[sn1-1].s_off > read_len - SEED_LEN) {
        slc[sn1-1].s_off = read_len - SEED_LEN; 
    } 
    l_len = l_bg; 
    //len = read_len/2-SEED_LEN/2;
    for(i = sn1+sn0; i < sn2-sn0; ++i) {
        slc[i].s_off = l_len; 
        if(l_len < read_len - SEED_LEN - l_len) {
            len = l_len; 
        } else {
            len = read_len - SEED_LEN - l_len; 
        }
        slc[i].ext_num = len/16;  //17 = 16+1
        if(slc[i].ext_num > query->max_ext) 
            slc[i].ext_num = query->max_ext; 
        l_len += l_seed;
    }
    if(slc[sn2-sn0-1].s_off > read_len - SEED_LEN) {
        slc[sn2-sn0-1].s_off = read_len - SEED_LEN; 
    } 
    int bg = sn1+sn0, ed = bg + sn1/2+1;
    int sf, ext_num; 
    for(i = bg; i < ed; ++i) {
        sf = slc[i].s_off;
        ext_num = slc[i].ext_num;
        j = sn2 - sn0 -1 - (i - bg); 
        slc[i].s_off = slc[j].s_off;
        slc[i].ext_num = slc[j].ext_num;
        slc[j].s_off = sf;
        slc[j].ext_num = ext_num; 
    }        
        
//fprintf(stderr, "%u  -----------\n\n", __LINE__);
//for(i = 0; i < sn2; ++i) fprintf(stderr, "%u, slc[%u].s_off = %d\n", __LINE__, i, slc[i].s_off);
    return;
}

void init_seed_pos2(query_t *query, seed_t *seed, struct StackTree *sTree)
{

    int read_len = query->l_seq;
    int i, j, len;
    seed_slc_t *slc   = seed->slc; 
    //确定第seed_id新的种子序列           

    int l_seed, l_len, r_len, num;
    slc[0].s_off  = query->m_off - SEED_LEN/2 ;
    l_seed = query->l_seq/query->seed_num;
    l_len = slc[0].s_off;
    r_len = slc[0].s_off;  
    for(i = 1; i < query->seed_num; ++i) {
        if(i % 2 == 1) {
            l_len -= l_seed;
            slc[i].s_off = l_len; 
        } else{
            r_len += l_seed;
            slc[i].s_off = r_len;
        }
    }
   
    slc[query->seed_num].s_off  = query->m_off - SEED_LEN/2 ;
    l_len = slc[query->seed_num].s_off; 
    r_len = l_len;
    for(i = query->seed_num + 1; i < query->seed_num*2; ++i) {
        if(i % 2 == 1) {
            l_len -= l_seed;
            slc[i].s_off = l_len; 
        } else{
            r_len += l_seed;
            slc[i].s_off = r_len;
        }
    }
    
    
    for(i = 0; i < query->seed_num*2; ++i) fprintf(stderr, "%u, slc[%u].s_off = %d\n", __LINE__, i, slc[i].s_off);

    len = read_len/2-SEED_LEN/2;
    for(i = 0; i < query->seed_num*2; ++i) { 
        if(slc[i].s_off < len) {
            slc[i].ext_num = slc[i].s_off/16;  //17 = 16+1
        } else{
            slc[i].ext_num = (read_len-(slc[i].s_off+SEED_LEN))/16;
        }
        

        if(slc[i].ext_num > query->max_ext) slc[i].ext_num = query->max_ext;
        slc[i].bgn = 0;
        slc[i].num = 0;
        slc[i].end = 0;
    }
    /*  
    for(i = 0; i < seed->m_num; ++i){
        seed->m_flg[seed->m_pos[i]/8] = 0;
    }
    */
    slc[0].sTree_bg  = 0;  
    slc[0].sTree_num = 0;
    sTree->bg_buf    = 0;
    sTree->len_buf   = 0;

    return;
}
int bwt_exact_aln1(idx_t *fm_idx, uint32_t hash_boundry[],  query_t *query, uint32_t idx[2], uint32_t seq_off[2], struct SubBuf *sub)
{
    int i, j;
    int read_len = query->l_seq;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t (*idx_buf)[2] = query->is_rev?sub->idx_rev:sub->idx_for;
    int bg = 0, ed = read_len - 1 ;

    //经过BWT比对，产生种子序列   
    uint32_t k, l, len = 0;
    /*  
    uint32_t seq12 = lkt_seq2LktItem(read_seq, ed-11, ed);
    k = fm_idx->fastmap->item[seq12];
    l = fm_idx->fastmap->item[seq12+1]-1;
    l -= get_12mer_correct(hash_boundry, l);
    ed = ed - 11;    
    */
    k = idx[0], l = idx[1];
    ed = seq_off[0]; 
    idx_buf[ed][0] = k;
    idx_buf[ed][1] = l;
    ed--;
    if(l +1 - k > 0){ 
        //len = 12;
        len = seq_off[1] - seq_off[0];
        for(i = ed; i >= bg; --i) {
            idx_buf[i][0] = 1;
            idx_buf[i][1] = 0;
            uint32_t num = bwt_match_exact_alt(fm_idx->bwt, 1, read_seq+i, &k, &l);
            if(num == 0) break;
            idx_buf[i][0] = k;
            idx_buf[i][1] = l;
        }
        len += ed - i; 
    }
    idx_buf[LEN_READ][0] = len;
     
   
    return len;
}

int bwt_exact_aln(idx_t *fm_idx, uint32_t hash_boundry[],  query_t *query, struct SubBuf *sub)
{
 
    int i, j;
    int read_len = query->l_seq;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t (*idx_buf)[2] = query->is_rev?sub->idx_rev:sub->idx_for;
    int trim = sub->trm_r;
    int bg = 0 + sub->trm_l, ed = read_len -1 - sub->trm_r;
    //经过BWT比对，产生种子序列   
    uint32_t k, l, len = 0;
    uint32_t seq12 = lkt_seq2LktItem(read_seq, ed-11, ed);
    k = fm_idx->fastmap->item[seq12];
    l = fm_idx->fastmap->item[seq12+1]-1;
    l -= get_12mer_correct(hash_boundry, l);
    ed = ed - 11;    
    idx_buf[ed][0] = k;
    idx_buf[ed][1] = l;
    ed--;
    if(l +1 - k > 0){ 
        len = 12;
        uint32_t num;
        for(i = ed; i >= bg; --i) {
            num = bwt_match_exact_alt(fm_idx->bwt, 1, read_seq+i, &k, &l);
            if(num == 0) break;
            idx_buf[i][0] = k;
            idx_buf[i][1] = l;
        }
        len += ed - i; 
    }
    idx_buf[LEN_READ][0] = len;
     

   
    return len;
}


//int gen_seed_kmer(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed, struct StackTree *sTree)

int gen_best_qual_kmer(idx_t *fm_idx, uint32_t hash_boundry, query_t *query, seed_t *seed)
{


}
int gen_seed_kmer(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed)
{
    uint32_t k, l, bgn, end, num = 0;

    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq;
    uint32_t (*idx_buf)[17][3] = seed->idx_buf; 

    seed_slc_t *slc = seed->slc; 
    int seed_id  = seed->id;
    
    uint32_t seq12;
    uint8_t pseed[SEED_LEN+8], *seed_p; 

    int i, j;
    //选定候选种子序列 
////fprintf(stderr, "%u, seed_id = %u, s_off = %u\n", __LINE__, seed_id, slc[seed_id].s_off);
    /*  
    if(slc[seed_id].drct == 0){
        seed_p = f_read_seq+slc[seed_id].s_off - 8;
    } else{
        seed_p = r_read_seq+slc[seed_id].s_off - 8;
    }
    */
    int s_off, ed, bg;
    if(slc[seed_id].s_off < 8) {
        s_off = 0;
        ed = slc[seed_id].s_off+SEED_LEN;
    } else {
        s_off = slc[seed_id].s_off - 8;
        ed = SEED_LEN + 8; 
    } 
   
////fprintf(stderr, "%u, seed_id = %u, s_off = %u\n", __LINE__, seed_id, slc[seed_id].s_off);
    //if(seed_id < seed->MAX_SIZE/2){
    if(seed_id < query->seed_num){
        seed_p = f_read_seq + s_off;
    } else{
        seed_p = r_read_seq + s_off;
    }
//fprintf(stderr, "%u, seed_id = %u, s_off = %u, ed = %u\n", __LINE__, seed_id, slc[seed_id].s_off, ed);

    if(seed->err == NULL){
        printf("%u, seed->err == NULL!!!\n", __LINE__);
        exit(1);
    }
    //struct ext_err_t *err = seed->err+seed_id;


    //err->sd_err[0] = 0;


    int N_pos = 0;
    for(i = ed-1; i >=0; --i) {
        if(seed_p[i] == 4) {
            //err->sd_err[++err->sd_err[0]] = i;
            pseed[i] = 0;
          
            if(ed-i < 16) {
                slc[seed_id].len = 0;  
                slc[seed_id].bgn = 0; 
                slc[seed_id].end = 0; 
                slc[seed_id].num = 0;
                return 0; 
            } else {
                N_pos = i+1;
                break; 
            } 
         
        } else{
            pseed[i] = seed_p[i];
        } 
    }
    //经过BWT比对，产生种子序列   
    //seq12 = lkt_seq2LktItem(pseed, 16, 27);
    seq12 = lkt_seq2LktItem(pseed, ed-12, ed-1);
    k = fm_idx->fastmap->item[seq12];
    l = fm_idx->fastmap->item[seq12+1]-1;
    l -= get_12mer_correct(hash_boundry, l);
    idx_buf[seed_id][0][0] = k;
    idx_buf[seed_id][0][1] = l;
    idx_buf[seed_id][0][2] = l+1-k;
    
    ed = (SEED_LEN+8-12)-1;
    bg = N_pos;
    if(slc[seed_id].s_off < 8) {
        ed = (SEED_LEN + slc[seed_id].s_off - 12) - 1;
    }
    slc[seed_id].len = 0;
    for(i = ed; i >= bg; --i) {
        num = bwt_match_exact_alt(fm_idx->bwt, 1, pseed+i, &k, &l);
        if(num == 0) break;
        idx_buf[seed_id][ed-i+1][0] = k;
        idx_buf[seed_id][ed-i+1][1] = l;
        idx_buf[seed_id][ed-i+1][2] = num; 
//fprintf(stderr, "%u, seed_id = %u, i = %u, k = %u, l = %u, num = %u\n", __LINE__, seed_id, i, k, l, num);
        //old_num = num;
        slc[seed_id].len++;  
        if(slc[seed_id].len == 8) {
            slc[seed_id].bgn = k; 
            slc[seed_id].end = l; 
            slc[seed_id].num = l+1-k; 
        }
    }
    
    if(i < 8){ 
        num = slc[seed_id].num; 
    }
    /*  
    for(seed_id = 0; seed_id < seed->MAX_SIZE; ++seed_id){
        fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
        for(i = 0; i < 17; ++i) {
            fprintf(stderr, "idx_buf[%u]: [0] = %u, ", i, idx_buf[seed_id][i][0]);
            fprintf(stderr, "[1] = %u, ", idx_buf[seed_id][i][1]);
            fprintf(stderr, "[2] = %u\n", idx_buf[seed_id][i][2]);
        }
        fprintf(stderr, "----------\n");

    }
    */
   
    return num;
}

int gen_seed_kmer_2(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed)
{
    uint32_t k, l, bgn, end, num = 0;
    int read_len = query->l_seq;
    uint8_t *seq = query->is_rev? query->rseq:query->seq;    
    seed_slc_t *slc = seed->slc; 
    int sid  = seed->id;
    uint32_t seq12;
    uint8_t *pseed = seq; 
    int i, j;
    int r_pos, ed, bg;
    r_pos = seed->slc[sid].s_off;     
    
    ed = r_pos; 
    int N_pos = 0;
    slc[sid].len = 0;  
    slc[sid].bgn = 0; 
    slc[sid].end = 0; 
    slc[sid].num = 0;

    for(i = ed-1; i >= ed-12; --i) {
        if(pseed[i] == 4) {
         
            if(ed-i < 16) {
                slc[sid].s_off = i;
                return 0; 
            } else {
                break; 
            } 
         
        }     
    }
    //经过BWT比对，产生种子序列   
    seq12 = lkt_seq2LktItem(pseed, ed-12, ed-1);
    k = fm_idx->fastmap->item[seq12];
    l = fm_idx->fastmap->item[seq12+1]-1;
    l -= get_12mer_correct(hash_boundry, l);
   
    ed -= 12;
    bg = 0;
    slc[sid].len = 12;
    slc[sid].bgn = k; 
    slc[sid].end = l; 
    slc[sid].num = l+1-k; 
    for(i = ed-1; i >= bg; --i) {
        if(pseed[i] == 4) break;
        num = bwt_match_exact_alt(fm_idx->bwt, 1, pseed+i, &k, &l);
        if(num == 0) { break; }
        slc[sid].len++;  
        slc[sid].bgn = k; 
        slc[sid].end = l; 
        slc[sid].num = l+1-k; 
        slc[sid].s_off = i; 
        if(num <= IS_SMLSIZ/4) { break; }
    }
    num = slc[sid].num; 
    if(slc[sid].len <= 16) { num = 0;}
    return num;
}
int set_12mer_qual(query_t *query)
{
    int read_len = query->l_seq;
    uint8_t *q_flg = query->q_flg;
    uint8_t *qual = query->qual;
    int thres = 3*4*5;
    int th_buf[5] = {1+33, 10+33, 15+33, 20+33, 25+33};  
    int val[5];
    int q_val[2], q_err[2];  
    int o_val = 0, c_val = 0;
    int i, j, k;
    val[0] = thres + 1;
    val[1] = 2*thres/3;
    val[2] = thres/2;
    val[3] = thres/3;
    val[4] = thres/4 ;

/*  
    val[0] = thres + 1;
    val[1] = thres + 1;
    val[2] = 2*thres/3;
    val[3] = thres/2;
    val[4] = thres/3;
*/
    for(i = 0; i < read_len; ++i) {
        q_flg[i] = 0;  
    }
    for(i = 0; i < 12; ++i) {
        q_err[0] = 0;
        for(k = 0; k < 3; k++) {
            if(qual[i] < th_buf[k]) {
                q_err[0] = val[k];  
                break;
            } 
        }
        o_val += q_err[0];
    }    
    if(o_val > thres) {
        q_flg[0] += 16;
        q_flg[11] += 1; 
    }

    for(i = 12; i < read_len; ++i) {
        q_val[0] = qual[i];
        q_val[1] = qual[i-12];
        for(j = 0; j < 2; ++j){
            q_err[j] = 0;
            for(k = 0; k < 3; k++) {
                if(q_val[j] < th_buf[k]) {
                    q_err[j] = val[k];  
                    break;
                } 
            }
        }
        c_val = o_val - q_err[1] + q_err[0];
        o_val = c_val;
        
        if(c_val > thres) {
            q_flg[i-11] += 16;
            q_flg[i] += 1; 
        }
    }
    return 0;
}

int gen_seed_kmer_1(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed)
{
    uint32_t k, l, bgn, end, num = 0;
    int read_len = query->l_seq;
    //uint8_t *f_read_seq = query->seq; 
    //uint8_t *r_read_seq = query->rseq;
    seed_slc_t *slc = seed->slc; 
    int seed_id  = seed->id;
    uint8_t *seq;
    uint32_t seq12;
    uint8_t pseed[SEED_LEN+8], *seed_p; 

    int i, j;
    int s_off, ed, bg, len;
    //int thres_qual = sub->thres_qual;
    int thres_qual = 30-1 + 33;
    //int min_len_qual = sub->min_len_qual;
    int min_len_qual = 16;
    
    uint8_t *qual = query->qual;
    int bg_buf[5]= {}, ed_buf[5] = {}, len_buf[5] = {};
    int thres_buf[5] = {30-1+33, 25-1+33, 20-1+33, 15-1+33, 0};
    //bg = 0, ed = 0, len = 0;
    int qual_cls = 1; 
    int qual_val;
    int rev = query->is_rev; 
    for(i = 0; i < read_len; ++i) {
        qual_val = rev? qual[i]:qual[read_len-1-i];
        for(j = 0; j < qual_cls; ++j) {
            if(qual_val < thres_buf[j]) {
                if(len_buf[j] > min_len_qual){
                    if(len_buf[j] >= ed_buf[j] - bg_buf[j]){
                        ed_buf[j] = i;
                        bg_buf[j] = ed_buf[j] - len_buf[j];
                    }
                }
                len_buf[j] = 0;      
            } else {
                ++len_buf[j];         
            }
        }
        /*  
        if(qual[i] < thres_qual) {
            if(len > min_len_qual){
                if(len >= ed - bg){
                    ed = i;
                    bg = ed - len;
                }
            }
            len = 0;      
        } else {
            ++len;         
        }
        */
    
    
    }
    for(j = 0; j < qual_cls; ++j) {
        if(len_buf[j] > min_len_qual){
            if(len_buf[j] >= ed_buf[j] - bg_buf[j]){
                ed_buf[j] = i;
                bg_buf[j] = ed_buf[j] - len_buf[j];
            }
        } 
        len_buf[j] = ed_buf[j] - bg_buf[j];  
    }
    int best_j = -1;
    for(j = 0; j < qual_cls; ++j) {
        if(len_buf[j] >= SEED_LEN + 2) {
            best_j = j;
if(j != 0) printf("%u, best_j = %d\n", __LINE__, best_j);
            break;
        }  else if(len_buf[j] >= min_len_qual && best_j < 0) {
            best_j = j; 
        } 
    
    }
    if(best_j < 0) {
        return 0; 
    }
    bg = bg_buf[best_j];
    ed = ed_buf[best_j];
    len = len_buf[best_j];
//printf("%u, bg = %d, ed = %d, len = %d\n", __LINE__, bg, ed, len);
/*  
    if(ed_buf[j] - bg_buf[j] < min_len_qual) return 0;  
    else len_buf[j] = ed_buf[j] - bg_buf[j];  
*/

fprintf(stderr, "%u, bg = %d, ed = %d\n", __LINE__, bg, ed); 

    int find_num = 0;
    int for_ed = ed;
    int rev_ed = query->l_seq - bg; 
    //int rev;
    int si, off;
    for(rev= 0; rev < 2; ++rev) { 
        seed_id = 2*query->seed_num+rev; 
        seq = rev?query->rseq:query->seq;
        bg = rev?(rev_ed-len): (for_ed-len); 
        ed = rev?rev_ed:for_ed; 
    
        if(rev == 0) {
            for(si = 0; si < query->seed_num; ++si) {
                off = seed->slc[si].s_off + SEED_LEN;
                if(off == ed) {
                    ed = ed - 2;
                    break;
                }     
            }    
        } else {
            for(si = query->seed_num; si < 2*query->seed_num; ++si) {
                off = seed->slc[si].s_off + SEED_LEN;
                if(off == ed) {
                    ed = ed - 2;
                    break;
                }     
            }
        } 

fprintf(stderr, "%u, bg = %d, ed = %d\n", __LINE__, bg, ed); 
        seq12 = lkt_seq2LktItem(seq, ed-12, ed-1);
fprintf(stderr, "seq12 = %x\n", seq12);
        if(seq12 = 0xFFFFFFFF) {return 0;}
        k = fm_idx->fastmap->item[seq12];
        l = fm_idx->fastmap->item[seq12+1]-1;
        l -= get_12mer_correct(hash_boundry, l);
        ed -= (12 + 1);
        
fprintf(stderr, "%u, bg = %d, ed = %d\n", __LINE__, bg, ed); 
        slc[seed_id].len = 0;
        for(i = ed; i >= bg; --i) {
            num = bwt_match_exact_alt(fm_idx->bwt, 1, seq+i, &k, &l);
            if(num == 0) break;
            slc[seed_id].len++;  
            slc[seed_id].bgn = k; 
            slc[seed_id].end = l; 
            slc[seed_id].num = l+1-k; 
            slc[seed_id].s_off = i;
            if(slc[seed_id].len == 8) break;      
        } 
        if(slc[seed_id].len + 12 >= 16) find_num++;
    
fprintf(stderr, "%u, bg = %d, ed = %d\n", __LINE__, bg, ed); 
    }
fprintf(stderr, "%u, bg = %d, ed = %d\n", __LINE__, bg, ed); 
    return find_num;
}

int gen_seed_kmer1(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed, struct StackTree *sTree)
{
    uint32_t k, l, bgn, end, num = 0;

    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq;

    seed_slc_t *slc = seed->slc; 
    int seed_id  = seed->id;
    
    uint32_t seq12;
    uint8_t pseed[SEED_LEN], *seed_p; 

    int i, j;
    //------------------------------------------------------------ 
    //保存扩展比对结果 
    if(seed_id > 0) {
        slc[seed_id-1].sTree_bg  = sTree->bg_buf;        
        if(slc[seed_id-1].num >0 && sTree->len_buf >0) {
            slc[seed_id-1].sTree_num = sTree->len_buf;
        } else{
            slc[seed_id-1].sTree_num = 0;
        }
    }
    //选定候选种子序列 
    if(slc[seed_id].drct == 0){
        seed_p = f_read_seq+slc[seed_id].s_off;
    } else{
        seed_p = r_read_seq+slc[seed_id].s_off;
    }
    //struct ext_err_t *err = seed->err+seed_id;
    //err->sd_err[0] = 0;
    for(i = 0; i < SEED_LEN; ++i) {
        if(seed_p[i] == 4) {
            //err->sd_err[++err->sd_err[0]];
            pseed[i] = 0;
        } else{
            pseed[i] = seed_p[i];
        } 
    }
    //经过BWT比对，产生种子序列   
    seq12 = lkt_seq2LktItem(pseed, 8, 19);
    k = fm_idx->fastmap->item[seq12];
    l = fm_idx->fastmap->item[seq12+1]-1;
    l -= get_12mer_correct(hash_boundry, l);
    num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, pseed, &k, &l);

    bgn = k;
    end = l;
    if(num >0){//种子序列产生成功
        slc[seed_id].bgn = bgn; 
        slc[seed_id].end = end; 
        slc[seed_id].num = num; 
    }         
    seed->bgn = bgn;
    seed->end = end;
    seed->num = num;
    return num;
}
/*  
void set_aln_path(idx_t *fm_idx, query_t *query, seed_t *seed, struct StackTree *sTree)
{
   
    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq;
    uint8_t *qual = query->qual;

    uint8_t r_qual[LEN_READ] = {};
    uint8_t *read_seq;
    uint32_t pos_buf[2];
    uint32_t best_pos[10];

    int i, j, len;
    //确定第seed_id新的种子序列           
    seed_slc_t *slc   = seed->slc; 
    seed_slc_t *slc_i  = seed->slc_i; 
    int seed_id  = seed->id;
    uint32_t bgn, end, num = 0, k, l;
    int bg, ed, max_cls, max_row, slc_drct;
    int drct_score[2] = {};
    uint8_t *m_flg = seed->m_flg;
    uint32_t *m_pos = seed->m_pos;
    uint32_t bg_idx, idx_num; 
    int seed_off;

    for(i = 0; i < 6; ++i) {
        if(slc[i].num > 0) { 
            bg = slc[i].sTree_bg;
            ed = bg + slc[i].sTree_num; 
            max_cls = 0;
            max_row = bg;

            slc[i].max_cls = max_cls;
            slc[i].max_row = max_row; 
            //read进行打分，正向read和反向read分别打分 
            drct_score[slc[i].drct] += max_cls+1;
        } else{
            drct_score[slc[i].drct] = 0; 
            slc[i].max_cls = 0;
            slc[i].max_row = 0; 
        }
    }
    //2. 选定后续比对的路径.
//fprintf(stderr, "%u\n", __LINE__);
//for(j = 0; j < 19; ++j) //fprintf(stderr, "slc[%u].s_off = %u\n", j, slc[j].s_off);
//+++++++++++++++++++++++++++++++++++++++++++
    for(i = 6; i < 18; ++i) {
        len = read_len/2-SEED_LEN/2; 
        if(slc[i].s_off < len) {
            slc[i].ext_num = (slc[i].s_off-1)/16;  //17 = 16+1
        } else{
            slc[i].ext_num = (read_len-(slc[i].s_off+SEED_LEN)-1)/16;
        }
        if(slc[i].ext_num > query->max_ext) slc[i].ext_num = query->max_ext;
    }
    //3. mid_flag数组当中标记出现过的pos
    uint32_t pos, pos_8, pos_q, pos_r;
    for(i = 0; i < 6; ++i){
        if(slc[i].drct != slc_drct) continue;
        max_cls = slc[i].max_cls;
        seed_off = slc[i].s_off + query->m_off[slc[i].drct];
        if(max_cls > 0){
            bg = slc[i].sTree_bg;
            ed = bg + slc[i].sTree_num; 
            for(j = bg; j < ed; ++j){
 
            } 
        } else if (slc[i].num >0){
            
            bg_idx = slc[i].bgn;    
            idx_num = slc[i].num;    
          
            for(k = 0; k < num; ++k){

                pos = bwt_sa(fm_idx->bwt, bg_idx+k); 
                pos = pos + max_cls*16 - seed_off;
                pos_8 = pos / 8;//8间隔存储数据
                pos_q = pos_8 / 8;
                pos_r = pos_8 % 8;
                if((m_flg[pos_q] && (1<<pos_r)) == 0){
                    m_flg[pos_q] |= (1<<pos_r);
                    m_pos[seed->m_num++] = pos_8;
                } 
            }
        }
    }//end  for(i = 0; i < 6; ++i)+++++++++++
    seed->back_num = 0;

    return;
}
*/
/*  
int aln_overlap_pos(idx_t *fm_idx,query_t *query, seed_t *seed, struct SubBuf *sub,  struct StackTree *sTree)
{
    uint32_t bg, ed, num = 0, k, l;
    int i, j, flg, sw_flg = 0;

    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq;
    uint8_t *qual = query->qual;
    uint8_t *read_seq;
    uint32_t pos_buf[2];
    uint32_t best_pos[10];
    uint32_t pos, pos_q, pos_r, pos_8;
    int max_cls, max_row; 
    seed_slc_t *slc   = seed->slc; 
    int seed_id  = seed->id;
    
    uint32_t back_num  = seed->back_num;
    uint32_t *pos_back = seed->pos_tmp;
    uint8_t  *m_flg    = seed->m_flg;
    uint32_t *m_pos    = seed->m_pos;
    uint32_t bg_idx, idx_num;   
    int seed_off;    
    i = seed_id -1; 
    slc[i].sTree_bg  = sTree->bg_buf;        
    if(slc[i].num >0 && sTree->len_buf >0) {
        slc[i].sTree_num = sTree->len_buf;
    } else{
        slc[i].sTree_num = 0;
    }
    //++++++++++++++++++++++++++++++++++++++++ 
    if(slc[i].num > 0) { 
        bg = slc[i].sTree_bg;
        ed = bg + slc[i].sTree_num; 
        max_cls = 0;
        max_row = bg;

//fprintf(stderr, "bg = %u, ed = %u\n", bg, ed);
        slc[i].max_cls = max_cls;
        slc[i].max_row = max_row; 
        //read进行打分，正向read和反向read分别打分 
        //drct_score[slc[i].drct] += max_cls+1;
    } else{
        //drct_score[slc[i].drct] = 0; 
        slc[i].max_cls = 0;
        slc[i].max_row = 0; 
    }
    //++++++++++++++++++++++++++++++++++++++
    //以下代码查找重叠的pos
    max_cls = slc[i].max_cls;
    if(seed->slc_drct == 0) {
        seed_off = slc[i].s_off + query->l_off[0]; 
        read_seq = query->seq;
    } else{
        seed_off = slc[i].s_off + query->r_off[1]; 
        read_seq = query->rseq;
    }
    read_len = query->l_seq;
    if(max_cls > 0){
        bg = slc[i].sTree_bg;
        ed = bg + slc[i].sTree_num; 
        for(j = bg; j < ed; ++j){
 
        }//end for(j = bg; j < ed; ++j) +++++++++++++++++++++++++++++
    } else if (slc[i].num >0){// if(max_cls > 0)
        bg_idx = slc[i].bgn;    
        idx_num = slc[i].num;    

//fprintf(stderr, "%u, bg_idx = %u, idx_num = %u\n",__LINE__,  bg_idx, idx_num);
        for(k = 0; k < idx_num; ++k){
            pos = bwt_sa(fm_idx->bwt, bg_idx+k); 
            pos = pos + max_cls*16 - seed_off;
            pos_8 = pos / 8;//8间隔存储数据
            pos_q = pos_8 / 8;
            pos_r = pos_8 % 8;
            //+++++++++++++++++++++++++++++++++++++++++++++ 
            //判断是否存在重叠的pos,如果有用sw方法进行比对。
            l = 0; 
            if((m_flg[pos_q] && (1<<pos_r)) == 1){
                pos_buf[l++] = pos; 
            }
            if(pos%8 <4) {
                pos_q = (pos_8-1)/8; 
                pos_r = (pos_8-1)%8;
            } else{
                pos_q = (pos_8+1)/8; 
                pos_r = (pos_8+1)%8;
            }
            if((m_flg[pos_q] && (1<<pos_r)) == 1){
                pos_buf[l++] = pos; 
            }
            if(l >0) {
              //sw_flg = AlgnPos_1(fm_idx, read_seq, read_len, pos_buf, l, sub, best_pos); 
              sw_flg = AlgnPos(fm_idx, query, l, 1, sub); 
            } else{
                pos_back[back_num++] = pos; 
            }
            if(sw_flg > 0){
                //如果比对精度高跳出，否则继续。 
            } //sw方法比对代码段结束。
            //---------------------------------------------
        }//end for(k = 0; k < num; ++k) ++++++++++++++++++++++++
    } //end else if (slc[i].num >0)   ++++++++++++++++++++++++++

    return sw_flg;
}
*/
/*  
void sup_pos_flg(seed_t *seed)
{

    uint32_t back_num  = seed->back_num;
    uint32_t *pos_back = seed->pos_tmp;
    uint32_t *m_pos = seed->m_pos; 
    int i;
    //3. mid_flag数组当中标记出现过的pos
    uint8_t *m_flg = seed->m_flg;
    uint32_t bg_idx, ed_idx, num; 
    int seed_off;
    uint32_t pos, pos_8, pos_q, pos_r;
    for(i = 0; i < back_num; ++i){
        pos = pos_back[i];
        pos_8 = pos / 8;//8间隔存储数据
        pos_q = pos_8 / 8;
        pos_r = pos_8 % 8;
        if((m_flg[pos_q] && (1<<pos_r)) == 0){
            m_flg[pos_q] |= (1<<pos_r);
            m_pos[seed->m_num++] = pos_8;
        } 
    }

    return;
}
*/
/*  
//以下临时保存废弃的代码+++++++++++++++++++++++++++++
for(i = 6; i < 9; ++i){
    max_cls = slc[i].max_cls;
    if(slc_drct == 0) { 
        seed_off = slc[i].s_off + query->l_off[0];
    } else{
        seed_off = slc[i].s_off + query->r_off[1];
    }
    if(max_cls > 0){
        bg = slc[i].sTree_bg;
        ed = bg + slc[i].sTree_num; 
        for(j = bg; j < ed; ++j){
        
            if(max_cls == sTree->back_buf[j][5]) {
                bg_idx = sTree->back_buf[j][7];    
                num = sTree->back_buf[j][8];    
              
                for(k = 0; k < num; ++k){
                    pos = bwt_sa(fm_idx->bwt, bg_idx+k); 
                    pos = pos + max_cls*16 - seed_off;
                    pos_8 = pos / 8;//8间隔存储数据
                    pos_q = pos_8 / 8;
                    pos_r = pos_8 % 8;
                    if((m_flg[pos_q] && (1<<pos_r)) == 0){
                        m_flg[pos_q] |= (1<<pos_r);
                        m_pos[seed->m_num++] = pos_8;
                    } 
                }
            }
        } 
    } else if (slc[i].num >0){
        
        bg_idx = slc[i].bgn;    
        num = slc[i].num;    
      
        for(k = 0; k < num; ++k){
            pos = bwt_sa(fm_idx->bwt, bg_idx+k); 
            pos = pos + max_cls*16 - seed_off;
            pos_8 = pos / 8;//8间隔存储数据
            pos_q = pos_8 / 8;
            pos_r = pos_8 % 8;
            if((m_flg[pos_q] && (1<<pos_r)) == 0){
                m_flg[pos_q] |= (1<<pos_r);
                m_pos[seed->m_num++] = pos_8;
            } 
        }
    }
}//end  for(i = 6; i < 9; ++i)+++++++++++
*/
//-----------------------------------------------------
/*  
uint32_t gen_seed_seq(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed, struct SubBuf *sub,  struct StackTree *sTree)
{
    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq;
    uint8_t *qual = query->qual;
    uint8_t r_qual[LEN_READ] = {};
    uint8_t *read_seq;
    int read_len;
    uint32_t pos_buf[2];
    uint32_t best_pos[10];
    int i, j, k, flg, sw_flg;
    int seed_slc_size = seed->slc_size;
    seed_slc_t *slc   = seed->slc; 
    seed_slc_t *slc_i  = seed->slc_i; 
    int seed_id  = seed->id;
    uint32_t seq12;
    uint8_t *pseed, *seed_buf= seed->seed_buf; 
    uint32_t bgn, end, num = 0, k, l;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    while(seed_id < seed_slc_size){ 
        //以下代码段对于给定的seed_id，通过BWT比对产生初始种子序列
        num = get_seed_seq(fm_idx, hash_boundry, query, seed);
        if(num == 0) {
            seed_id = ++seed->id;
            continue; 
        } else{
            seed_id = seed->id; 
            break; 
        } 
    }
    //+++++++++++++++++++++++++++++++++++++++++ 
    //初始化变量
    if(seed_id == 0) {
        //以下代码段初始化种子信息
        init_seed_pos(query, seed, sTree);
        seed_slc_size = 18;
    }

    while(seed_id < seed_slc_size){ 
        //1. 评估已经比对的结果; 
        //2. 选定优先比对序列，即选定正反向序列 
        //3. 选定后续比对的路径.
        if(seed_id == 6) {
            //以下代码对上述比对结果进行评估，选定比对路径
            set_aln_path(fm_idx, query, seed, sTree);
        }//end if(seed_id == 6)++++++++ 
        //查找重叠pos,并进行SW比对 
        if(seed_id >= 6 && seed_id < 9) {
            num = aln_overlap_pos(fm_idx,query, seed, sub,  sTree);
            if(num > 0) {
                //？？？
                break;
            }

        }//end  if(seed_id >= 6 && seed_id < 9) ++++++++++++++++++++++++
        //评估已经比对的结果
        //判断是否可以进行多种子的重叠pos计算
        //如果满足条件查找重叠pos
        //如果不满足条件或者比对失败，确定后续比对路径
        if(seed_id ==9){
            sup_pos_flg(seed);
        } 
        //生成新的种子序列
        if(seed_id >=9 && seed_id < 12){
            num = aln_overlap_pos(fm_idx,query, seed, sub,  sTree);
            if(num > 0) {
                //？？？
                break;
            }
        } 
        //评估已经比对的结果
        //判断是否可以进行多种子的重叠pos计算
        //如果满足条件查找重叠pos
        //如果不满足条件或者比对失败，确定后续比对路径
        if(seed_id ==12){
            num = seed->m_num;      
            uint32_t *m_pos = seed->m_pos;
            uint8_t  *m_flg = seed->m_flg;
            for(i = 0; i < num; ++i){
                pos = m_pos[i]/8; 
                m_flg[pos] = 0;
            }
        }
        //生成新的种子序列 
        if(seed_id >=12 && seed_id < 15){
        
        }
        if(seed_id ==15){
             
        }
        //生成新的种子序列 
        if(seed_id >=15 && seed_id < 18){
        
        }
        if(seed_id ==18){
             
        }
        //产生种子新的序列 
        //while(seed_id < seed_slc_size){
        while(seed_id < 12){
            num = gen_seed_kmer(fm_idx, hash_boundry, query, seed, sTree);
            if(num == 0) {
                seed_id = ++seed->id;
                continue; 
            } else{
                seed_id = seed->id; 
                break; 
            }
        }//end while(seed_id < seed_slc_size)+++++++++++++
    }// end while(seed_id < seed_slc_size)+++++++++++ 


    return num;
}

uint32_t aln_seed_seq_2(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed, int seed_slc_buf[] )
{
    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq;
    uint8_t *qual = query->qual;
//uint8_t r_qual[LEN_READ] = {};
    uint8_t r_qual[LEN_READ] = {};

    int i, j;
    //int seed_slc_size = seed->slc_size;
    int seed_slc_size = 4;
 
    //确定第seed_id新的种子序列           
    seed_slc_t *slc   = seed->slc; 
    seed_slc_t *slc_i = seed->slc_i; 
    seed_slc_t *slc_b = seed->slc_b; 
    aln_mer8_t *aln_8 = seed->aln_8;
    int find_num = seed->find_num;
    int cur_row  = seed->cur_row; 
    int seed_id  = seed->id;
    uint32_t seq12;

    uint8_t *pseed, *seed_buf= seed->seed_buf; 
    //int seq_off = read_len/2-slc_i[seed_id].s_off;
    int seq_off = read_len/2+seed_slc_buf[seed_slc_buf[0]];

    int seed_flg = seed->seed_flg;
    uint32_t bgn, end, num = 0, k, l;

//fprintf(stderr, "%u, seed_id = %u, seq_off = %d, seed_slc_buf[0] = %d\n", __LINE__, seed_id, seq_off, seed_slc_buf[0]);
    int MIN_Q = 0;
    if(seed_id == 0) {
        int row_b = 0, row = 0;
        //for(i = 0; i < seed->MAX_SIZE; ++i){
        for(i = 0; i < seed_slc_size; ++i){
            uint32_t seed_bg, seed_ed;
            //seed_bg = read_len/2-(slc_i[i].s_off-slc_i[i].h_off);
            seed_bg = seq_off+slc_i[i].h_off;
            seed_ed = seed_bg+12;
            int sum = 0;
            if(slc_i[i].drct == 0) {
                for(j = seed_bg; j < seed_ed; ++j) {
                     if( qual[j] < MIN_Q ){ 
                         //sum++; 
                     }
                 }
            } else{
                for(j = seed_bg; j < seed_ed; ++j) {
                    if( r_qual[j] < MIN_Q ){
                         //sum++; 
                    }
                }
            }
            if(sum > 2) { continue;}
            else if(sum > 0){ 
                memcpy(slc_b+row_b, slc_i+i, sizeof(seed_slc_t));
                row_b++;
            } else{
                memcpy(slc+row, slc_i+i, sizeof(seed_slc_t));
                row++;
            }
        } 
        memcpy(slc+row, slc_b, row_b*sizeof(seed_slc_t));
        //3）赋值当前read_id的seed_slc_size;
        seed_slc_size = row + row_b;
        seed->slc_size = seed_slc_size;
        find_num = 0; 

//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
    } //end if(seed_id == 0)  

//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //种子序列指针赋值
    if(find_num==0) {
        if(slc[seed_id].drct == 0){
            pseed = f_read_seq+seq_off;
        } else{//slc[seed_id].drct ==1
            pseed = r_read_seq+seq_off;
        }     
    }             

//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
    while(seed_id < seed_slc_size){ 
        if(find_num==0) {
            //seq_off = read_len/2-slc_i[seed_id].s_off;
            if(slc[seed_id].drct == 0){
                pseed = f_read_seq+seq_off;
            } else{//slc[seed_id].drct ==1
                pseed = r_read_seq+seq_off;
            }
        }   

//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);

        if(find_num == 0) {//选定新的种子序列  

//fprintf(stderr, "%u, find_num = %u\n", __LINE__, find_num);
            if(slc[seed_id].flag == 0) {//通过BWT精确比对选定种子序列

//fprintf(stderr, "%u, find_num = %u\n", __LINE__, find_num);
                seq12 = lkt_seq2LktItem(pseed, 8, 19);

//fprintf(stderr, "%u, seq12 = %u\n", __LINE__, seq12);
                k = fm_idx->fastmap->item[seq12];
                l = fm_idx->fastmap->item[seq12+1]-1;
//fprintf(stderr, "%u, find_num = %u\n", __LINE__, find_num);
                l -= get_12mer_correct(hash_boundry, l);
//fprintf(stderr, "%u, find_num = %u\n", __LINE__, find_num);
                num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, pseed, &k, &l);
                bgn = k;
                end = l;
                seed_flg = 0;
//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
            } else if(slc[seed_id].flag ==1){//通过BWT选定近似种子序列集合
                cur_row = 0;

//fprintf(stderr, "%u, find_num = %u\n", __LINE__, find_num);
                seed->hash_12mer = lkt_seq2LktItem(pseed, 8, 19);
//fprintf(stderr, "%u, hash_12mer = %u\n", __LINE__, seed->hash_12mer);
                seed->seed_8mer_L = pseed;
                seed->seed_seq = pseed;                
                find_num = get_seed_8mer_bwt_L(fm_idx, hash_boundry, seed);
//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
                if(find_num == 0) {
                    seed_id++;
                    continue;
                }
                bgn = seed->aln_buf[cur_row][2]; 
                end = seed->aln_buf[cur_row][3]; 
                num = end+1- bgn;
                seed_flg = 1;

//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
                if(find_num == 1) { 
                    find_num = 0; 
                }
            } 
        } else{//if(find_num >0)

//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
            if(seed_flg ==1) {//第一类近似比对，即BWT近似比对, 种子序列继续
                seed_id--;
                cur_row++; 
                bgn = seed->aln_buf[cur_row][2]; 
                end = seed->aln_buf[cur_row][3]; 
                num = end+1- bgn;
                //seed_flg = 1;
            } 

            if(find_num - cur_row == 1) {//cur_row最后一个行的情况 
                find_num = 0; 
            }
        }

//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
        seed->find_num = find_num ;
        seed->cur_row = cur_row ; 
        seed->id = seed_id ;
        seed->seed_flg = seed_flg;
        seed->bgn = bgn;
        seed->end = end;
        seed->num = num;

        if(num >0){//种子序列比对失败
            break; 
        }
//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
        seed->id = ++seed_id;
    }
    seed->id = seed_id;
    return num;
}
*/
/*  
int AlgnPos_smpos(idx_t *fm_idx, uint8_t ext_cls, query_t *query, int seed_off, int num_pair, struct ExtBlck *eBlck,  struct SubBuf *sub)
{
    uint64_t *sw_ed = sub->sw_ed[query->is_rev]; 
    uint64_t *sw_to = sub->sw_to[query->is_rev]; 
    uint32_t *pos_ed = sub->pos_ed[query->is_rev]; 
    uint32_t *pos_to = sub->pos_to[query->is_rev]; 

    
    uint32_t nxt_pnt, nxt_end, pos, bgn, end, num;
    uint8_t  nxt_flg;
    uint32_t i, j;
    int Flg_Algn = 0;
//fprintf(stderr, "%u, num_pair = %u\n", __LINE__, num_pair);
    uint32_t *pos_buf = sub->pos_buf;
    //int *err_buf = (int *)sub->err_buf;
    int seq_off[2];
   

    for (i = 0; i < num_pair; ++i){
        nxt_pnt = sub->pair_out->pair_arry[i].nxtpnt;
        nxt_flg = sub->pair_out->pair_arry[i].nxtflg;

//fprintf(stderr,"%u, i = %u, nxt_flg = %u, nxt_pnt = %u\n", __LINE__, i, nxt_flg, nxt_pnt);
if(nxt_flg > sub->NEXT_EXT_NUM) {
    printf("%u, nxt_flg = %u\n", __LINE__, nxt_flg);
    exit(1);
}
        nxt_end = nxt_pnt + nxt_flg;
        int l_off = sub->pair_out->pair_arry[i].l_off; 
        int r_off = sub->pair_out->pair_arry[i].r_off;
        sub->err_sum[0] = sub->pair_out->pair_arry[i].err;
        uint32_t pos_i = 0;
////fprintf(stderr,"%u, ------------------------------\n", __LINE__);
//fprintf(stderr,"%u, i = %u, nxt_flg = %u, nxt_pnt = %u\n", __LINE__, i, nxt_flg, nxt_pnt);
//fprintf(stderr,"%u\n", __LINE__);

        //uint32_t L_end = seed_off-(ext_cls+1)*16 - (sub->aln_L[21]-64);//其中10是种子左端点中心的偏移量
        uint32_t L_end = seed_off-(ext_cls+1)*16 - l_off;//其中10是种子左端点中心的偏移量
//fprintf(stderr,"%u, nxtflg = %u\n", __LINE__, nxt_flg);
        if(nxt_flg > IS_SMLSIZ) {
//++++++++++++++++++++++++++++++++++
            uint32_t (*ext_idx)[2] = (eBlck+ext_cls)->head_extidx + nxt_pnt;
//fprintf(stderr, "%u, ext_cls = %d\n", __LINE__, ext_cls);
            bgn = ext_idx[0][0];
            num = ext_idx[0][1];
//fprintf(stderr, "%u, bgn = %d, num = %u\n", __LINE__, bgn, num);
            end = bgn + num;
            for(j = bgn; j < end; j++){
                pos = bwt_sa(fm_idx->bwt, j) - L_end;
                pos_buf[pos_i++] = pos;
//fprintf(stderr, "j = %u, pos = %u, pos_i = %u, pos_buf[pos_i] = %u\n", j, pos, pos_i, pos_buf[pos_i-1]);
            }

//----------------------------------         
        } else if(nxt_flg>1){
//fprintf(stderr,"%u\n", __LINE__);
            for(j=nxt_pnt; j<nxt_end; j++){
                pos = bwt_sa(fm_idx->bwt, j) - L_end;
                pos_buf[pos_i++] = pos;
//fprintf(stderr, "j = %u, pos = %u, pos_i = %u, pos_buf[pos_i] = %u\n", j, pos, pos_i, pos_buf[pos_i-1]);
            }
        } else{
//fprintf(stderr,"%u, pos_i = %u\n", __LINE__, pos_i);
            pos_buf[pos_i++] = nxt_pnt - L_end; 
//fprintf(stderr, "pos[%u] = %u\n", 0, nxt_pnt);
        }
//fprintf(stderr,"%u\n", __LINE__);
        uint32_t best_pos[10];
        for(j = 0; j< pos_i; j++){
            //fprintf(stderr, "pos_buf[j] = %u\n", pos_buf[j]);
        }
//fprintf(stderr,"%u\n", __LINE__);
        //fprintf(stderr, "%u pos_i = %u\n", __LINE__, pos_i);
        pos_i = clean_swed_pos(pos_buf, pos_i, query, sub);
        //fprintf(stderr, "%u pos_i = %u\n", __LINE__, pos_i);
        //seq_off[0] = seed_off-(ext_cls+1)*16-16       - l_off;
        seq_off[0] = seed_off-(ext_cls+1)*16          - l_off;
        seq_off[1] = SEED_LEN+seed_off+(ext_cls+1)*16 + r_off;
        int pos_n[3];
        pos_n[0] = pos_i;
        pos_i = slc_aln_pos(fm_idx,  query, seq_off, pos_n, sub);
        pos_n[0] = pos_i;

        Flg_Algn = 0;
        //if(pos_i > 0) Flg_Algn =  AlgnPos(fm_idx, query, pos_i, 1, sub);
        if(pos_i > 0) Flg_Algn =  AlgnPos(fm_idx, query, pos_i, sub);
//fprintf(stderr,"%u\n", __LINE__);
        if(sub->aln_r->num > 0) {
//if(sub->aln_r->score[0][0] <= ALN_THRES_SCORE)  break;
        }


    }
//fprintf(stderr,"%u\n", __LINE__);   
    return Flg_Algn;
}
*/

//int AlgnPos_buf(idx_t *fm_idx, uint8_t ext_cls, query_t *query, int seed_off, int num_pair, struct ExtBlck *eBlck,  struct SubBuf *sub)
int AlgnPos_buf(idx_t *fm_idx, struct ExtBlck *eBlck, query_t *query, seed_t *seed, int num_pair, uint8_t ext_cls, struct SubBuf *sub)
{
    /*  
    uint64_t *sw_ed  = sub->sw_ed[query->is_rev]; 
    uint64_t *sw_to  = sub->sw_to[query->is_rev]; 
    uint32_t *pos_ed = sub->pos_ed[query->is_rev]; 
    uint32_t *pos_to = sub->pos_to[query->is_rev]; 
    uint32_t *pos_bk = sub->pos_bk[query->is_rev]; 
    */
    
    uint32_t nxt_pnt, nxt_end, pos, bgn, end, num;
    uint8_t  nxt_flg;
    uint32_t i, j;
    int Flg_Algn = 0;

    uint32_t *pos_buf = sub->pos_buf;
    int seq_off[2];
    int seed_off = seed->slc[seed->id].s_off;

    for (i = 0; i < num_pair; ++i){
        nxt_pnt = sub->pair_out->pair_arry[i].nxtpnt;
        nxt_flg = sub->pair_out->pair_arry[i].nxtflg;
        if(nxt_flg > sub->NEXT_EXT_NUM) {
            printf("%u, nxt_flg = %u\n", __LINE__, nxt_flg);
            exit(1);
        }
        nxt_end = nxt_pnt + nxt_flg;
        int l_off = sub->pair_out->pair_arry[i].l_off; 
        int r_off = sub->pair_out->pair_arry[i].r_off;
        sub->err_sum[0] = sub->pair_out->pair_arry[i].err;
        uint32_t pos_i = 0;
        uint32_t L_end = seed_off-(ext_cls+1)*16 - l_off;
        if(nxt_flg > IS_SMLSIZ) {
            uint32_t (*ext_idx)[2] = (eBlck+ext_cls)->head_extidx + nxt_pnt;
            bgn = ext_idx[0][0];
            num = ext_idx[0][1];
            end = bgn + num;
            for(j = bgn; j < end; j++){
                pos = bwt_sa(fm_idx->bwt, j) - L_end;
                pos_buf[pos_i++] = pos;
            }
        } else if(nxt_flg>1){
            for(j=nxt_pnt; j<nxt_end; j++){
                pos = bwt_sa(fm_idx->bwt, j) - L_end;
                pos_buf[pos_i++] = pos;
            }
        } else{
            pos_buf[pos_i++] = nxt_pnt - L_end; 
        }
fprintf(stderr, "%u,  pos_i = %u\n", __LINE__, pos_i);
//print_pos(fm_idx, pos_i, pos_buf, query->l_seq);               
 

        uint32_t best_pos[10];

////fprintf(stderr, "%u, sub->thres_sw_to = %u, sub->thres_pos_num= %u, pos_i = %u\n", __LINE__,sub->thres_sw_to, sub->thres_pos_num, pos_i);       
        //pos_i = classify_pos(pos_buf, pos_i, query, sub);

        
        seq_off[0] = seed_off-(ext_cls+1)*16          - l_off;
        seq_off[1] = SEED_LEN+seed_off+(ext_cls+1)*16 + r_off;
        seed->cls = ext_cls+1;
        
////fprintf(stderr, "%u, pos_i = %u, seq_off[0] = %u, seq_off[1] = %u\n", __LINE__, pos_i, seq_off[0], seq_off[1]);       
fprintf(stderr, "%u, pos_i =%u\n", __LINE__, pos_i);
//print_pos(fm_idx, pos_i, pos_buf, query->l_seq);               
        pos_i = classify_pos(fm_idx, query, seed, pos_buf, pos_i, seq_off, sub);
////fprintf(stderr, "%u, pos_i =%u\n", __LINE__, pos_i);
//print_pos(fm_idx, pos_i, pos_buf, query->l_seq);               
         
     
       
      
        if(sub->SLC_SWITCH[2]) pos_i = slc_aln_pos(fm_idx,  query, seq_off, pos_i, sub);
        /*  
        if(seed->cls == query->max_ext) {
            //pos_i = clean_swed_pos(pos_buf, pos_i, query, sub);
        } else {
            if(sub->SLC_SWITCH[2]) pos_i = slc_aln_pos(fm_idx,  query, seq_off, pos_n, sub);
        } 
        */
        //pos_i = clean_swed_pos(pos_buf, pos_i, query, sub);
        //if(sub->SLC_SWITCH[2]) pos_i = slc_aln_pos(fm_idx,  query, seq_off, pos_n, sub);
    
////fprintf(stderr, "%u, pos_i =%u\n", __LINE__, pos_i);
//print_pos(fm_idx, pos_i, pos_buf, query->l_seq);               
        //if((query->b0 > 0 && query->b0 > sub->err_buf[0][0] +sub->delta) || sub->err_buf[0][0] < query->candi_thres + 40) pos_i = 0; //????? 

////fprintf(stderr, "%u, pos_i = %u, query->b0 = %u, sub->err_buf[0][0] = %d, sub->delta = %d, query->candi_thres = %d\n", __LINE__, pos_i, query->b0, sub->err_buf[0][0], sub->delta, query->candi_thres);       
/*  
for(j = 0; j < pos_i; ++j) {
    int sc1 = 0; 
    query->b1 = sc1;
    int sc0 = aln_pos_both(fm_idx,  query, pos_buf+j,  seq_off, sub);
    uint32_t pos = pos_buf[j];
    rec_aln_info(query, pos, sc0, sub);
}
*/


        Flg_Algn = 0;
        if(pos_i > 0) {
            if(sub->ALN_SWITCH[2] > 0) {
                
                for(i = 0; i < pos_i; ++i) {
                    int sc = aln_pos_both(fm_idx,  query, pos_buf+i,  
                                        seq_off, sub);
                    rec_aln_info(query, pos_buf[i], sc, sub);
                }
            } else {
                Flg_Algn = AlgnPos(fm_idx, query, pos_i,  sub);
            }
        }
        
    }
    return Flg_Algn;
}
void OutAlgnInfo(uint32_t (*algn_out)[5]){

	return;
}
int get_ch_bwt(int r_R, struct ExtBlck *eB, struct SubBuf *sub, int flg, uint8_t *ch){
	uint8_t*  smbwt;
	uint8_t*  sumbuf;
    uint8_t    *seq;
	uint32_t   alg_in[3] = {0,0,0};
    uint32_t (*alg_out)[4];
    uint32_t DataNum = 0;
	uint32_t num_seqL  = eB->num_seqL;
	uint32_t num_seqR  = eB->num_seqR;
    if(flg<1){
		smbwt = eB->head_smbwt + eB->bwtL;
		DataNum = num_seqL;
		alg_out = sub->seqL_out;
		seq = sub->ext_seqL;
	}else{
		smbwt = eB->head_smbwt + eB->bwtR;
		DataNum = num_seqR;
		alg_out = sub->seqR_out;
		seq = sub->ext_seqR;
	}
    uint32_t n = DataNum;
    uint8_t *cnt = smbwt+((n+1)/2)*8;

    if(n<=MIN_BWT_SIZE){
        //get_ch_min(r_R, smbwt, seq, 7, n, alg_out); 
        //break; 
    } else if(n<=NO_BWT_SUM_SIZE){ 
        uint8_t cnt1[8][17];
        get_ch_nosum(r_R, smbwt, seq, 7, cnt1, n, alg_out); 
    } else if(n <= 255) {
        get_ch_255(r_R, smbwt, seq, 7, (uint8_t (*)[17])cnt, n, alg_out);     
    } else{//>=256
        uint32_t InitInfo[10];
        InitInfo[0] = n;
        get_ch_large_alt(r_R, smbwt, cnt, seq, InitInfo, alg_out);
    }
    if(sub->seqR_out[0][0] == 2 && sub->seqL_out[0][0] == 2) //fprintf(stderr, "dataNum = %u\n", n);
    return alg_out[0][0];

}
int Algn_indel_once(struct ExtBlck *eB, struct SubBuf *sub,  int flg){
	uint8_t  *smbwt;
	uint8_t  *sumbuf;
	uint8_t  *seq;
    int8_t  *aln_in;	
	uint32_t (*alg_out)[4];
	uint32_t n = 0;
	uint32_t num_seqL  = eB->num_seqL;
	uint32_t num_seqR  = eB->num_seqR;
    uint8_t *algn_row;
    int r_once = 0; 

    sub->flg_off = flg;
    if(flg<1){
		smbwt = eB->head_smbwt + eB->bwtL;
		n = num_seqL;
		alg_out = sub->seqL_out;
		aln_in = sub->aln_L;
        seq = sub->ext_seqL + sub->l_seq_delta - aln_in[13];
        sub->seqL_aln_old = alg_out[0][0];
        algn_row = sub->algnL_row;
    }else{
        smbwt = eB->head_smbwt + eB->bwtR;
		n = num_seqR;
		alg_out = sub->seqR_out;
		aln_in = sub->aln_R;
        seq = sub->ext_seqR;
        sub->seqR_aln_old = alg_out[0][0];
        algn_row = sub->algnR_row;
    }
    int r, i;
    if(aln_in[13] == -1 ){
        aln_in[16] = alg_out[0][0];        
    }else if(aln_in[13] == 1){
        aln_in[17] = alg_out[0][0];        
    }else if(aln_in[13] == -2){
        aln_in[18] = alg_out[0][0];   
    }else if(aln_in[13] == 2){
        aln_in[19] = alg_out[0][0];        
    }else{
        printf("%s, %u, Error: InDel = %d!\n",__func__, __LINE__, aln_in[13]);
        exit(1);
    }
    uint8_t *cnt = smbwt+((n+1)/2)*8;
    aln_in[10] = 8; 
    aln_in[11] = 7;
    aln_in[12] = 0;//once : 0, all:1 
    //aln_in[15] = 3;
    aln_in[15] = sub->ins_err;
    if(n<=MIN_BWT_SIZE){
        aln_in[10] = 0;
        aln_in[11] = 1;
        r_once = align_min_indel(smbwt, n,seq, aln_in, algn_row, alg_out, sub); 
     } else if(n<=NO_BWT_SUM_SIZE){ 
        uint8_t cnt1[8][17];
        gen_cnt(smbwt, n, cnt1);
        r_once = align_nosum_indel(smbwt, cnt1, n, seq, aln_in, algn_row, alg_out, sub); 
    } else if(n <= 255) {
        r_once = align_255_indel(smbwt, (uint8_t (*)[17])cnt, n, seq, aln_in, algn_row, alg_out, sub);     
    } else if(n > 255){//>=256
        r_once = align_large_indel(smbwt, cnt, n, seq, aln_in, algn_row, alg_out, sub);    
    }
    for(i = 1; i <= alg_out[0][0]; ++i) { 
        r  = alg_out[i][0];
        algn_row[r/8] = 0; 
    } 
    return r_once;
}


int Algn_indel_all(struct ExtBlck *eB, struct SubBuf *sub,  int flg){
	uint8_t  *smbwt;
	uint8_t  *sumbuf;
	uint8_t  *seq;
    int8_t  *aln_in;	
	uint32_t (*alg_out)[4];
	uint32_t n = 0;
	uint32_t num_seqL  = eB->num_seqL;
	uint32_t num_seqR  = eB->num_seqR;
    uint8_t *algn_row;
    int r_all = 0, r_once;  
    if(flg<1){
		smbwt = eB->head_smbwt + eB->bwtL;
		n = num_seqL;
		alg_out = sub->seqL_out;
        aln_in = sub->aln_L;
        seq = sub->ext_seqL + sub->l_seq_delta - aln_in[13];
        sub->seqL_aln_old = alg_out[0][0];
        algn_row = sub->algnL_row;
    }else{
		smbwt = eB->head_smbwt + eB->bwtR;
		n = num_seqR;
		alg_out = sub->seqR_out;
        aln_in = sub->aln_R;
        seq = sub->ext_seqR;
        sub->seqR_aln_old = alg_out[0][0];
        algn_row = sub->algnR_row;
    }
    int st_row;
    if(aln_in[13] == -1 ){
        st_row = aln_in[16];        
    }else if(aln_in[13] == 1){
        st_row = aln_in[17];        
    }else if(aln_in[13] == -2){
        st_row = aln_in[18];   
    }else if(aln_in[13] == 2){
        st_row = aln_in[19];        
    }else{
        printf("%s, %u, Error: InDel = %d!!!!!\n",__func__, __LINE__, aln_in[13]);
        exit(1);
    }
    int r, i;
    for(i = st_row+1; i <= alg_out[0][0]; ++i) { 
        r  = alg_out[i][0];
        algn_row[r/8] |= 1<<(7-(r)%8); 
    } 
    uint8_t *cnt = smbwt+((n+1)/2)*8;
    aln_in[10] = 8; 
    aln_in[12] = 1; 
    //aln_in[15] = 3;
    aln_in[15] = sub->ins_err;
    if(n<=MIN_BWT_SIZE){
        aln_in[10] = 0;
        aln_in[11] = 1;
        r_all = align_min_indel(smbwt, n,seq, aln_in, algn_row, alg_out, sub); 
     } else if(n<=NO_BWT_SUM_SIZE){ 
        uint8_t cnt1[8][17];
        gen_cnt(smbwt, n, cnt1);
        for(i = 0; i < 8; ++i){
            aln_in[11] = i;
            r_once = align_nosum_indel(smbwt, cnt1, n, seq, aln_in, algn_row, alg_out, sub); 
            if(r_once >0) r_all += 2;
        }
    } else if(n <= 255) {
        for(i = 0; i < 8; ++i){
            aln_in[11] = i;
            r_once = align_255_indel(smbwt, (uint8_t (*)[17])cnt, n, seq, aln_in, algn_row, alg_out, sub);     
            if(r_once >0) r_all += 2;
        }
    } else if(n > 255){//>=256
        for(i = 0; i < 8; ++i){//++++++++++++++++++++++++++++++++?????????????????
            aln_in[11] = i;
            r_once = align_large_indel(smbwt, cnt, n, seq, aln_in, algn_row, alg_out, sub);    
            if(r_once >0) r_all += 2;
        }
    }
    
    for(i = 1; i <= alg_out[0][0]; ++i) { 
        r  = alg_out[i][0];
        algn_row[r/8] = 0; 
    } 
    if(aln_in[13] == -1 ){
        aln_in[16] = alg_out[0][0];        
    }else if(aln_in[13] == 1){
        aln_in[17] = alg_out[0][0];        
    }else if(aln_in[13] == -2){
        aln_in[18] = alg_out[0][0];   
    }else if(aln_in[13] == 2){
        aln_in[19] = alg_out[0][0];        
    }else{
        printf("%s, %u, Error: InDel = %d!!!!!\n",__func__, __LINE__, aln_in[13]);
        exit(1);
    }
    return r_all;
}

int Algn_sub_all(struct ExtBlck *eB, struct SubBuf *sub,  int flg){
	uint8_t  *smbwt;
	uint8_t  *sumbuf;
	uint8_t  *seq;
    int8_t  *aln_in;	
	uint32_t (*alg_out)[4];
	uint32_t n = 0;
	uint32_t num_seqL  = eB->num_seqL;
	uint32_t num_seqR  = eB->num_seqR;
    uint8_t *algn_row;
    int r_all = 0, r_once; 
    if(flg<1){
		smbwt = eB->head_smbwt + eB->bwtL;
		n = num_seqL;
		alg_out = sub->seqL_out;
		seq = sub->ext_seqL + sub->l_seq_delta;
        aln_in = sub->aln_L;
        sub->seqL_aln_old = alg_out[0][0];
        algn_row = sub->algnL_row;
    }else{
		smbwt = eB->head_smbwt + eB->bwtR;
		n = num_seqR;
		alg_out = sub->seqR_out;
		seq = sub->ext_seqR;
        aln_in = sub->aln_R;
        sub->seqR_aln_old = alg_out[0][0];
        algn_row = sub->algnR_row;
    }
//fprintf(stderr, "%u, n = %u\n", __LINE__, n);
    int i, r;  
    uint8_t *cnt = smbwt+((n+1)/2)*8;
    aln_in[10] = 8;  
    aln_in[12] = 1;  
    //aln_in[15] = 2;
    aln_in[15] = sub->sub_err;
    for(i = 1; i <= alg_out[0][0]; ++i) { 
        r  = alg_out[i][0];
        algn_row[r/8] |= 1<<(7-(r)%8); 
    }
    if(n<=MIN_BWT_SIZE){
        aln_in[10] = 1;
        aln_in[11] = 1;
        r_all = align_min(smbwt, n,seq, aln_in, algn_row, alg_out, sub); 
     } else if(n<=NO_BWT_SUM_SIZE){ 
        uint8_t cnt1[8][17];
        gen_cnt(smbwt, n, cnt1);
        for(i = 0; i < 8; ++i){
            aln_in[11] = i;
            r_once = align_nosum(smbwt, cnt1, n, seq, aln_in, algn_row, alg_out, sub); 
            if(r_once >0) r_all += 2;
        }
    } else if(n <= 255) {
        for(i = 0; i < 8; ++i){
            aln_in[11] = i;
            r_once = align_255(smbwt, (uint8_t (*)[17])cnt, n, seq, aln_in, algn_row, alg_out, sub);     
            if(r_once >0) r_all += 2;
        }
    } else if(n > 255){//>=256
        for(i = 0; i < 8; ++i){
            aln_in[11] = i;
            r_once = align_large(smbwt, cnt, n, seq, aln_in, algn_row, alg_out, sub);    
            if(r_once >0) r_all += 2;
        }
    }
    for(i = 1; i <= alg_out[0][0]; ++i) { 
        r  = alg_out[i][0];
        algn_row[r/8] = 0; 
    }
    return r_all;
}

int Algn_sub_once(struct ExtBlck *eB, struct SubBuf *sub, int flg){
	uint8_t  *smbwt;
	uint8_t  *sumbuf;
	uint8_t  *seq;
	int8_t  *aln_in;
	uint32_t (*alg_out)[4];
	uint32_t n = 0;
	uint32_t num_seqL  = eB->num_seqL;
	uint32_t num_seqR  = eB->num_seqR;
    uint8_t *algn_row;

    if(flg<1){
		smbwt = eB->head_smbwt + eB->bwtL;
		n = num_seqL;
		alg_out = sub->seqL_out;
		seq = sub->ext_seqL+sub->l_seq_delta;
	    aln_in = sub->aln_L;
        algn_row = sub->algnL_row;
    }else{
		smbwt = eB->head_smbwt + eB->bwtR;
		n = num_seqR;
		alg_out = sub->seqR_out;
		seq = sub->ext_seqR;
	    aln_in = sub->aln_R;
        algn_row= sub->algnR_row;
    }

    int r, i; 

    aln_in[10] = 8;
    aln_in[11] = 7;
    aln_in[12] = 0;
    aln_in[15] = sub->sub_err;
    int r_once = 0;
    uint8_t *cnt = smbwt+((n+1)/2)*8;
   
    for(i = 1; i <= alg_out[0][0]; ++i) { 
        r  = alg_out[i][0];
        algn_row[r/8] |= 1<<(7-(r)%8); 
    } 

fprintf(stderr, "%u, alg_out = %d, n = %d\n", __LINE__, alg_out[0][0], n);    
    if(n<=MIN_BWT_SIZE){
        aln_in[10] = 1;
        aln_in[11] = 1;
        r_once = align_min(smbwt, n,seq, aln_in,algn_row, alg_out, sub); 
    } else if(n<=NO_BWT_SUM_SIZE){ 
        uint8_t cnt1[8][17];
        gen_cnt(smbwt, n, cnt1);
        r_once = align_nosum(smbwt, cnt1, n, seq, aln_in,algn_row, alg_out, sub); 
    } else if(n <= 255) {
        r_once = align_255(smbwt, (uint8_t (*)[17])cnt, n, seq, aln_in, algn_row, alg_out, sub);     
    } else{//>=256
        r_once = align_large(smbwt, cnt, n, seq, aln_in, algn_row, alg_out, sub);
    }
     
fprintf(stderr, "%u, alg_out = %d\n", __LINE__, alg_out[0][0]);    

    for(i = 1; i <=alg_out[0][0]; ++i ){
        r = alg_out[i][0];
        int k = algn_row[r/8]&(1<<(7-r%8));
    }
    for(i = 1; i <= alg_out[0][0]; ++i) { 
        r = alg_out[i][0];
        algn_row[r/8] = 0; 
    } 
   
    return r_once;

}
/*  
int AlgnmExtSeq(struct ExtBlck *eB, struct SubBuf *sub, int flg){
	uint8_t  *smbwt;
	uint8_t  *sumbuf;
	uint8_t  *seq;
	int8_t  *aln_in;
	uint32_t (*alg_out)[4];
	uint32_t n = 0;
	uint32_t num_seqL  = eB->num_seqL;
	uint32_t num_seqR  = eB->num_seqR;
    uint8_t *algn_row;

    if(flg<1){
		smbwt = eB->head_smbwt + eB->bwtL;
		n = num_seqL;
		alg_out = sub->seqL_out;
		seq = sub->ext_seqL;
	    aln_in = sub->aln_L;
        algn_row = sub->algnL_row;
    }else{
		smbwt = eB->head_smbwt + eB->bwtR;
		n = num_seqR;
		alg_out = sub->seqR_out;
		seq = sub->ext_seqR;
	    aln_in = sub->aln_R;
        algn_row= sub->algnR_row;
    }
    aln_in[10] = 8;
    aln_in[11] = 7;
    aln_in[12] = 0;

    uint8_t *cnt = smbwt+((n+1)/2)*8;
//fprintf(stderr, "%u n_data = %u, aln_in[14] = %u\n", __LINE__, n, aln_in[14]);

    if(n<=MIN_BWT_SIZE){
        aln_in[10] = 2;
        aln_in[11] = 1;
 
        align_min(smbwt, n,seq, aln_in,algn_row, alg_out); 
    } else if(n<=NO_BWT_SUM_SIZE){ 
        uint8_t cnt1[8][17];
        gen_cnt(smbwt, n, cnt1);
        align_nosum(smbwt, cnt1, n, seq, aln_in,algn_row, alg_out); 
    } else if(n <= 255) {
        align_255(smbwt, (uint8_t (*)[17])cnt, n, seq, aln_in, algn_row, alg_out);     
    } else{//>=256
        align_large(smbwt, cnt, n, seq, aln_in, algn_row, alg_out);
    }
    if(sub->seqR_out[0][0] == 2 && sub->seqL_out[0][0] == 2) fprintf(stderr, "dataNum = %u\n", n);
    int i;

    //fprintf(stderr, "%u, %d alg_out[0][0] =%u, n_data = %u\n", __LINE__, flg, alg_out[0][0], n);
    for(i = 1; i <=alg_out[0][0]; ++i ){
        int r = alg_out[i][0];
        int k = algn_row[r/8]&(1<<(7-r%8));
        //fprintf(stderr, "seq_out[%u] = %u, algn_row = %u\n", i, alg_out[i][0], k);
    }
    for(i = 1; i <=n; ++i ){
       
        int k = algn_row[i/8]&(1<<(7-i%8));
        if(k >0){ 
            //fprintf(stderr, "seq_row = %u\n", i); 
        }
    }
    
    return alg_out[0][0];

}
*/

void test_smbwt_retire_seq(uint8_t *bwt, uint32_t DataNum, int i, uint8_t seq[16])
{
    int rot;
    rot = 0; 
 
    uint8_t cnt[8][17];
    gen_cnt(bwt, DataNum, cnt);
    uint32_t k = i;
    for(rot=0; rot<8; ++rot){
        uint32_t x = rot*((DataNum+1)/2*2);
        //uint8_t ch = __get_bwt2(bwt, x+k);
        uint8_t ch = bwt[(x+k)/2]; 
        ch = (k%2==0)?ch>>4:ch&15;
        seq[15-rot*2] = ch&3;
        seq[14-rot*2] = ch>>2;
        uint32_t occ = __dna2_count(bwt, x, x+k, ch);
        uint32_t C = cnt[rot][ch]; 
        k = C+occ;
    }
}
void test_AlgnBwtSml(uint8_t*bwt, uint32_t DataNum, uint8_t *seq_buf)
{

    uint8_t cnt[8][17];
    gen_cnt(bwt, DataNum, cnt);
    //恢复序列算法
    int i, rot;
    rot = 0; 
    uint8_t seq[16];
    for(i = 0; i < DataNum; ++i){
        uint32_t k = i;

        for(rot=0; rot<8; ++rot){
            uint32_t x = rot*((DataNum+1)/2*2);
            //uint8_t ch = __get_bwt2(bwt, x+k);
            uint8_t ch = bwt[(x+k)/2]; 
            ch = (k%2==0)?ch>>4:ch&15;
            seq[15-rot*2] = ch&3;
            seq[14-rot*2] = ch>>2;
            uint32_t occ = __dna2_count(bwt, x, x+k, ch);
            uint32_t C = cnt[rot][ch]; 
            k = C+occ;
        }
        memcpy(seq_buf+i*16, seq, 16); 
    }    
    
    
    return;
}
void test_smbwt_retire_seq_large(uint8_t*bwt, uint32_t DataNum, int i, uint8_t seq[16])
{

    uint8_t cnt[8][17];
    gen_cnt(bwt, DataNum, cnt);
    //恢复序列算法
    int rot;
    rot = 0; 
    
 
    uint32_t k = i;

    for(rot=0; rot<8; ++rot){
        uint32_t x = rot*((DataNum+1)/2*2);
        //uint8_t ch = __get_bwt2(bwt, x+k);
        uint8_t ch = bwt[(x+k)/2]; 
        ch = (k%2==0)?ch>>4:ch&15;
        seq[15-rot*2] = ch&3;
        seq[14-rot*2] = ch>>2;
        uint32_t occ = __dna2_count(bwt, x, x+k, ch);
        //uint32_t C = cnt[rot][ch]; 
        uint32_t C = __dna2_count_small(bwt, x, x+DataNum, ch); 
        k = C+occ;
    }
  

    
    
    return;
}
void test_AlgnBwtSml_large(uint8_t*bwt, uint32_t DataNum, uint8_t *seq_buf)
{

    uint8_t cnt[8][17];
    gen_cnt(bwt, DataNum, cnt);
    //恢复序列算法
    int i, rot;
    rot = 0; 
    uint8_t seq[16];
    for(i = 0; i < DataNum; ++i){
        uint32_t k = i;

        for(rot=0; rot<8; ++rot){
            uint32_t x = rot*((DataNum+1)/2*2);
            //uint8_t ch = __get_bwt2(bwt, x+k);
            uint8_t ch = bwt[(x+k)/2]; 
            ch = (k%2==0)?ch>>4:ch&15;
            seq[15-rot*2] = ch&3;
            seq[14-rot*2] = ch>>2;
            uint32_t occ = __dna2_count(bwt, x, x+k, ch);
            //uint32_t C = cnt[rot][ch]; 
            uint32_t C = __dna2_count_small(bwt, x, x+DataNum, ch); 
            k = C+occ;
        }
        memcpy(seq_buf+i*16, seq, 16); 
    }    
    
    
    return;
}
void test_smbwt(struct ExtBlck *eB, int flg, uint8_t seq[16])
{
    uint8_t *seq_buf, *seq_buf0, *bwt;
    uint32_t num, *sort_seq;
    
    if(flg<1){
		bwt = eB->head_smbwt + eB->bwtL;
		num = eB->num_seqL;
        //seq = sub->ext_seqL;
	}else{
		bwt = eB->head_smbwt + eB->bwtR;
		num = eB->num_seqR;
        //seq = sub->ext_seqR;
	}


    seq_buf = malloc(16*num*2);  

    test_AlgnBwtSml(bwt, num, seq_buf);


    int i, j;
   
    for(i = 0; i < num; ++i){
        int count = 0;;
        for(j = 0; j < 16; ++j){
            if(seq_buf[i*16+j] == seq[j]) ++count;
            else break;  
        }
        if(count == 16){ 
            break; 
        } 
    }
    if(i < num){
        //fprintf(stderr, "find cur_seq i = %u\n", i);
    } else{
    
        //fprintf(stderr, "can't find cur_seq i = %u\n", i);

        //fprintf(stderr, "Seq:\t");
        for(j = 0; j < 16; ++j){
            //fprintf(stderr, "%u", seq[j]);
        }
        //fprintf(stderr, "\n");
        for(i = 0; i < num; ++i){
            //fprintf(stderr, "%u\t", i);
            for(j = 0; j < 16; ++j){
                //fprintf(stderr, "%u", seq_buf[i*16+j]);
            }
            //fprintf(stderr, "\n");
        }
 
    }

    
/*
    seq_buf0 = seq_buf+num*16;
  
    for(i = 0; i < num; ++i){
        uint8_t ch;
        for(j = 0; j < 16; ++j){
            seq_buf0[i*16+j] = (sort_seq[i]>>((15-j)*2))&3;      
        }
    }

    for(i = 0; i < num; ++i){
        for(j = 0; j < 16; ++j){
            ////fprintf(stderr, "[%u, %u]: %u %u\n", i, j, seq_buf[i*16+j], seq_buf0[i*16+j]);
            if(seq_buf[i*16+j] != seq_buf0[i*16+j]){
                //fprintf(stderr, "test_smbwt error!!!!\n"); 
                exit(1);
            
            }  
        }
       

    
    }
*/
    free(seq_buf);
    return;
}
/*
void AlgnBwtBig(uint8_t*Bwt, uint8_t*Seq, uint32_t DataNum, uint32_t AlgIn[],uint32_t(*AlgOut)[3]){
	uint32_t len_sum = getBwtSumSize(DataNum);
	uint8_t* seq ;
	uint8_t* bwt ;
	uint8_t (*sum_buf)[17][len_sum];
	uint32_t(*alg_info)[3];
	seq = Seq;
	bwt = Bwt;
//+++++++++++++++++++++++++++++++++++++++++++
//临时屏蔽	
//sum_buf = (uint8_t (*)[17][len_sum])(bwt+8*((DataNum+1)/2));
//+++++++++++++++++++++++++++++++++++++++++++
	alg_info = AlgOut;

	return;
}

uint32_t getBwtSumSize(uint32_t DataNum){
	uint32_t siz_cnt, rng_num, len_sum;
	len_sum = 0;
	if(DataNum<=254*256) siz_cnt = 2 ;
	if(DataNum>256*256) siz_cnt = 4 ;

	rng_num = ((DataNum + 255) /256) ;
	len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;

	return len_sum;
}
*/
/*  
void getExtSeq(struct SubBuf *sub, uint32_t cls, uint8_t *read_seq){ 
	uint32_t i,bgnL,bgnR;
	uint8_t ch;
	sub->off_L = LEN_READ/2 - LEN_SEED/2 - LEN_EXT*(cls+1);
	sub->off_R = LEN_READ/2 + LEN_SEED/2 + LEN_EXT*cls;

	bgnL = sub->off_L;
	bgnR = sub->off_R;
	
	for(i=0;i<16;i++){
		sub->ext_seqL[i] = read_seq[bgnL+i];
		sub->ext_seqR[i] = read_seq[bgnR+i];
	}

    return;
}
*/
#include "ksort.h"
typedef struct{
    uint32_t x, y, z;
} turple_t;
#define ks_turple(a, b) ((a).z >(b).z)
KSORT_INIT(turple_t, turple_t, ks_turple)
#define set_align(aln, x) aln|=
int PairExtSeq_1(struct ExtBlck *eB, struct SubBuf *sub){

	uint32_t (*seqL)[3];
	uint32_t (*seqR)[3];
	uint32_t *seqRel;
    	
    	    
    uint32_t (*pair)[3];

	uint8_t *algnR_row;	
	uint8_t *algnRel_row;	
	
    uint32_t *relat;
	uint32_t *L2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	
    seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
	pair = sub->pair_out;
    algnR_row = sub->algnR_row; 
	algnRel_row = sub->algnRel_row;
    relat = eB->head_relat + eB->relat; 
	L2rel = relat + eB->num_relat; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int numL,numR,i,j,k,r_L,r_R,p_num;
	numL = seqL[0][0]; 
	numR = seqR[0][0];
//fprintf(stderr, "seqL[0][1] = %u, seqR[0][1] = %u\n", seqL[0][1], seqR[0][1]);
//fprintf(stderr, "%s, 1129, numL, numR = %u, %u\n", __func__, numL, numR);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//测试代码
/*  
    for(i =1; i <=eB->num_seqR; ++i){
        if(algnR_row[i/8] > 0) {
            printf("algnR_row[] init error!!\n"); 
            exit(1);
        }
    }
*/ 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int st_i = sub->seqR_aln_old;
/*  
    for(i=st_i; i<=numR; i++){
		r_R = seqR[i][0];
        algnR_row[r_R/8] = algnR_row[r_R/8] | (1<<(7-r_R%8));
//fprintf(stderr, "i, r_R = %u %u\n", i, r_R);
    }
*/

    if(pair[0][0]>0){
        p_num = pair[0][0]-pair[0][1]; 
        memcpy(pair+1, pair+1+pair[0][1], p_num*(3*sizeof(uint32_t)));
    } else{
        p_num = 0;
    }
//+++++++++++++++++++++++++++++++++++++++++++++++++++++    
st_i = sub->seqL_aln_old;
//st_i = 1;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++
    for(i=st_i; i<=numL; i++){
		r_L = seqL[i][0];
////fprintf(stderr, "i, r_L = %u %u\n", i, r_L);

////fprintf(stderr, "i = %u, L2rel[%u] = %u, %u\n", i, r_L, L2rel[r_L],L2rel[r_L+1]);
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
    	    int flag_Rel = algnRel_row[j/8] & (1<<(7-j%8)); 
            if(flag_Rel > 0) continue;	
            r_R = relat[j];



    		int flag_R = algnR_row[r_R/8] & (1<<(7-r_R%8));
            if(flag_R>0){
				for(k = 1; k <= numR; k++){
                    if(seqR[k][0] == r_R) break;
                }
                if(k > numR){
                    printf("%u, %s algnR_row[] error!, k =%u, numR = %u\n", __LINE__, __func__, k, numR);
                    exit(1);
                }
                p_num++;
                pair[p_num][0] = nxtpnt[j];
				pair[p_num][1] = nxtflg[j];
				pair[p_num][2] = seqL[i][1] + seqR[k][1];
////fprintf(stderr, "pair = %u, %u, %u\n", pair[p_num][0], pair[p_num][1], pair[p_num][2]); 
////fprintf(stderr, "%u, p_num = %u, j = %u, r_L = %u, r_R= %u\n", __LINE__, p_num, j, r_L, r_R);   
    	        algnRel_row[j/8] |= (1<<(7-j%8));
                seqRel[0]++;
                seqRel[seqRel[0]] = j;

            }
		}  
    }
 
////fprintf(stderr, "%u, %s p_num=%u, numR = %u, numL= %u\n", __LINE__, __func__, p_num, numR, numL);
    int s_num = 1;
    int b_num = p_num;
    if(p_num == 1){
        if(pair[1][1]<=IS_SMLSIZ){ 
            s_num = 2;
        } 
    }

	while(p_num>1){
		while(pair[s_num][1]<=IS_SMLSIZ){
			s_num++;
			if(s_num>b_num)break;
		}
		if(s_num>b_num)break;
		while(pair[b_num][1]>IS_SMLSIZ){
			b_num--;
			if(s_num>b_num)break;
		}
		if(s_num>b_num)break;
		if(s_num==b_num){
            if(pair[s_num][1]<=IS_SMLSIZ)s_num++; 
            else b_num--;
            break;
		}

		memcpy(pair+p_num+1, pair+s_num,3*sizeof(uint32_t));
		memcpy(pair+s_num, pair+b_num,3*sizeof(uint32_t));
		memcpy(pair+b_num, pair+p_num+1,3*sizeof(uint32_t));
		s_num++; 
		b_num--;
		if(s_num>b_num){
          break;
		}
		if(s_num == b_num){
          if(pair[s_num][1]<=IS_SMLSIZ)s_num++; 
          else b_num--;
          break;
		}	
	}
	pair[0][0] = p_num;
	pair[0][1] = s_num-1;
	pair[0][2] = 0;
	return (int) p_num;
}
/*  
int PairExtSeq_32_uni(int flag, struct ExtBlck *eB, struct SubBuf *sub, uint32_t *relat){
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint32_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + eB->num_relat; 
    R2rel = L2rel + eB->num_seqL+1; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
   
    int st_L = 0;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 0;
    int ed_R = sub->seqR_out[0][0];
 
    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, rror!!!\n", __LINE__, __func__);
        exit(1); 
    }
    if(flag == 0) {
        for(i=st_L; i<=numL; i++){
            r_L = seqL[i][0];
            for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                r_Rel = relat[j];
                int err_val = sub->err_sum[0];
                if(seqL[i][1] > 0) err_val += seqL[i][1];
                int max_err = sub->query_err + sub->delta;
                if(err_val < max_err) {
                    pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
                    pair_arry[p_num].nxtflg = nxtflg[r_Rel];
                    pair_arry[p_num].err = err_val;
                    pair_arry[p_num].l_off = seqL[i][3];
                    p_num++;                    
                }
             
            }  
        }
    }
    if(flag == 1) {
        for(i=st_R; i<=numR; i++){
            r_R = seqR[i][0];
            for(r_Rel=R2rel[r_R]; r_Rel<R2rel[r_R+1]; r_Rel++) {
                int err_val = sub->err_sum[0];
                if(seqR[i][1] > 0) err_val += seqR[i][1];
                int max_err = sub->query_err + sub->delta;
                if(err_val < max_err) {
                    pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
                    pair_arry[p_num].nxtflg = nxtflg[r_Rel];
                    pair_arry[p_num].err = err_val;
                    pair_arry[p_num].r_off = seqR[i][3];
                    p_num++;                   
                }
            }
        }
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}
*/
int PairExtSeq_8_uni(idx_t *fm_idx, struct ExtBlck *cB, uint8_t *relat, struct StackTree *sTree, int flag, struct SubBuf *sub){


    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint8_t *L2rel, *R2rel, *R2idx;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + cB->num_relat; 
    R2rel = L2rel + cB->num_seqL+1; 
    R2idx = R2rel + cB->num_seqR+1; 
    nxtpnt= cB->head_nxtpnt + cB->nxtpnt; 
	nxtflg= cB->head_nxtflg + cB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
    int st_L = 1;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 1;
    int ed_R = sub->seqR_out[0][0];
 
    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, rror!!!\n", __LINE__, __func__);
        exit(1); 
    }
fprintf(stderr, "%u, flag = %d\n", __LINE__, flag);
    int len_arry = sTree->len_arry;
    uint32_t (*ext_idx)[2];
    uint32_t cls = sTree->stck_arry[len_arry].cls;
    uint32_t cap_row, cap_flg;
    int max_err = sub->query_err + sub->delta;
    int err_val = sub->err_sum[0];
    uint32_t bgn, end, num, pos;
    if(flag == 0) {
        uint32_t cur_bg, cur_ed, bg_cap, bg_nxt, cap_num, top_num = 0, bot_num = 0, bg_idx, ed_idx;
fprintf(stderr, "%u, numL = %d\n", __LINE__, numL);
        for(i=st_L; i<=numL; i++){
            max_err = sub->query_err + sub->delta;
            err_val = sub->err_sum[0] + seqL[i][1];
            if(err_val > max_err) { continue; } 
            r_L = seqL[i][0];
            bgn = (uint32_t)-1;
            end = 0;
            bg_cap = 0, bg_nxt = 0, cap_num = 0;
            
            //r_Rel = relat[L2rel[r_L]];
            r_Rel = L2rel[r_L];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                bgn = bwt_get_idx(fm_idx->bwt, cap_row); 
                num = 1;
            } else if(cap_flg <= IS_SMLSIZ) {
                bgn = cap_row;
                num = cap_flg; 
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                pos = bwt_sa(fm_idx->bwt, bg_idx);
                num = ext_idx[0][1];
            }
            
            r_Rel = L2rel[r_L+1]-1;
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                end = bwt_get_idx(fm_idx->bwt, cap_row); 
                num = 1;
            } else if(cap_flg <= IS_SMLSIZ) {
                num = cap_flg; 
                end = cap_row + num - 1;
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                num = ext_idx[0][1];
                end = ext_idx[0][0] + ext_idx[0][1] -1; 
            }
            if(bgn > end) continue;
            num = end + 1 - bgn;
            if(num > 1) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;  
            } else {
                pair_arry[p_num].idx_bg = 0;
                pair_arry[p_num].idx_num = 0; 
                pair_arry[p_num].nxtpnt = nxtpnt[seqR[i][0]]+16;
                pair_arry[p_num].nxtflg = 1;  
            }
            pair_arry[p_num].err = err_val + seqR[i][1];
            pair_arry[p_num].l_off = sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3];
            p_num++;  
            
            
            /*   
            for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg == 1) continue;
                else break;
            }
            if(j == L2rel[r_L+1]) {
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
                    pair_arry[p_num].nxtpnt = cap_row;
                    pair_arry[p_num].nxtflg = cap_flg;
                    pair_arry[p_num].idx_bg = 0;
                    pair_arry[p_num].idx_num = 0;
                    pair_arry[p_num].err = err_val + seqL[i][1];
                    pair_arry[p_num].l_off = seqL[i][3];
                    pair_arry[p_num].r_off = sub->r_off;
                    p_num++;  
                } 
                continue;
            }
            top_num = j;
            for(j=L2rel[r_L+1]-1; j >= top_num; --j){
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg == 1) continue;
                else break;
            }
            if(j < top_num) {
                printf("%u, error, j = %d, top_num = %d\n", __LINE__, j, top_num);
                exit(1);
            } else {
               bot_num = L2rel[r_L+1] -1 -j;  
            }  
            r_Rel = relat[top_num];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg <= IS_SMLSIZ) {
                bg_idx = cap_row; 
            } else {
                ext_idx = cB->head_extidx + cap_row;
                bg_idx = ext_idx[0][0]; 
            }
            bg_idx -= top_num - L2rel[r_L]; 
            j = L2rel[r_L+1]-1-bot_num; 
            r_Rel = relat[j];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            uint32_t cap_num;
            if(cap_flg <= IS_SMLSIZ) {
                ed_idx = cap_row; 
                cap_num = cap_flg;
            } else {
                ext_idx = cB->head_extidx + cap_row;
                ed_idx = ext_idx[0][0];
                cap_num = ext_idx[0][1]; 
            }
            ed_idx += cap_num - 1 + bot_num; 
            pair_arry[p_num].nxtpnt = 0;
            pair_arry[p_num].nxtflg = 0;
            pair_arry[p_num].idx_bg = bg_idx;
            pair_arry[p_num].idx_num = ed_idx + 1 - bg_idx;
            pair_arry[p_num].err = err_val + seqL[i][1];
            pair_arry[p_num].l_off = seqL[i][3];
            pair_arry[p_num].r_off = sub->r_off;
            p_num++;
            */
        } // end for(i=st_L; i<=numL; i++)++++++++++++++++
    }
    if(flag == 1) {
        /*
        cap_row = sTree->stck_arry[len_arry].nxtpnt;
        if(sTree->cls == 0) {
            bg_idx = cB->bg_idx;
        } else {
            ext_idx = cB->head_extidx + cap_row; 
            bg_idx = ext_idx[0][0];
        }
        */
        uint32_t bg_idx;
fprintf(stderr, "%u, st_R = %d, numR = %d\n", __LINE__, st_R, numR);
//exit(1);
        for(i=st_R; i<=numR; i++){
            max_err = sub->query_err + sub->delta;
            err_val = sub->err_sum[0] + seqR[i][1];
//printf("%u, st_R = %d, numR = %d, err_val = %d, max_err = %d\n", __LINE__, st_R, numR, err_val, max_err);
            if(err_val > max_err) { continue; } 
            
            r_R = R2rel[seqR[i][0]];
            r_Rel = relat[r_R]; 
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                pos = cap_row;
                num = cap_flg; 
            } else if(cap_flg <= IS_SMLSIZ) {
                pos = bwt_sa(fm_idx->bwt, cap_row); 
                num = cap_flg;
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                pos = bwt_sa(fm_idx->bwt, bg_idx);
                num = ext_idx[0][1];
            }
            bgn = bwt_get_idx(fm_idx->bwt, pos + 16); 
            //----- 
            r_R = R2rel[seqR[i][0]+1]-1;
            r_Rel = relat[r_R]; 
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                pos = cap_row;
                num = cap_flg; 
            } else if(cap_flg <= IS_SMLSIZ) {
                num = cap_flg;
                pos = bwt_sa(fm_idx->bwt, cap_row+num-1); 
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                num = ext_idx[0][1];
                pos = bwt_sa(fm_idx->bwt, bg_idx+num-1);
            }
            end = bwt_get_idx(fm_idx->bwt, pos + 16); 
            if(bgn > end) { continue; }
            num = end + 1 -bgn; 
            if(num > 1) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;  
            } else {
                pair_arry[p_num].idx_bg = 0;
                pair_arry[p_num].idx_num = 0; 
                pair_arry[p_num].nxtpnt = nxtpnt[seqR[i][0]]+16;
                pair_arry[p_num].nxtflg = 1;  
            }
            pair_arry[p_num].err = err_val + seqR[i][1];
            pair_arry[p_num].l_off = sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3];
            p_num++;                   
        }
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}
int PairExtSeq_16_uni(idx_t *fm_idx, struct ExtBlck *cB, uint16_t *relat, struct StackTree *sTree, int flag, struct SubBuf *sub){


    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint16_t *L2rel, *R2rel, *R2idx;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + cB->num_relat; 
    R2rel = L2rel + cB->num_seqL+1; 
    R2idx = R2rel + cB->num_seqR+1; 
    nxtpnt= cB->head_nxtpnt + cB->nxtpnt; 
	nxtflg= cB->head_nxtflg + cB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
    int st_L = 1;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 1;
    int ed_R = sub->seqR_out[0][0];
 
    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, rror!!!\n", __LINE__, __func__);
        exit(1); 
    }
fprintf(stderr, "%u, flag = %d\n", __LINE__, flag);
    int len_arry = sTree->len_arry;
    uint32_t (*ext_idx)[2];
    uint32_t cls = sTree->stck_arry[len_arry].cls;
    uint32_t cap_row, cap_flg;
    int max_err = sub->query_err + sub->delta;
    int err_val = sub->err_sum[0];
    uint32_t bgn, end, num, pos;
    if(flag == 0) {
        uint32_t cur_bg, cur_ed, bg_cap, bg_nxt, cap_num, top_num = 0, bot_num = 0, bg_idx, ed_idx;
fprintf(stderr, "%u, numL = %d\n", __LINE__, numL);
        for(i=st_L; i<=numL; i++){
            max_err = sub->query_err + sub->delta;
            err_val = sub->err_sum[0] + seqL[i][1];
            if(err_val > max_err) { continue; } 
            r_L = seqL[i][0];
            bgn = (uint32_t)-1;
            end = 0;
            bg_cap = 0, bg_nxt = 0, cap_num = 0;
            
            //r_Rel = relat[L2rel[r_L]];
            r_Rel = L2rel[r_L];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                bgn = bwt_get_idx(fm_idx->bwt, cap_row); 
                num = 1;
            } else if(cap_flg <= IS_SMLSIZ) {
                bgn = cap_row;
                num = cap_flg; 
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                pos = bwt_sa(fm_idx->bwt, bg_idx);
                num = ext_idx[0][1];
            }
            
            r_Rel = L2rel[r_L+1]-1;
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                end = bwt_get_idx(fm_idx->bwt, cap_row); 
                num = 1;
            } else if(cap_flg <= IS_SMLSIZ) {
                num = cap_flg; 
                end = cap_row + num - 1;
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                num = ext_idx[0][1];
                end = ext_idx[0][0] + ext_idx[0][1] -1; 
            }
            if(bgn > end) continue;
            num = end + 1 - bgn;
            if(num > 1) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;  
            } else {
                pair_arry[p_num].idx_bg = 0;
                pair_arry[p_num].idx_num = 0; 
                pair_arry[p_num].nxtpnt = nxtpnt[seqR[i][0]]+16;
                pair_arry[p_num].nxtflg = 1;  
            }
            pair_arry[p_num].err = err_val + seqR[i][1];
            pair_arry[p_num].l_off = sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3];
            p_num++;  
            
            
            /*   
            for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg == 1) continue;
                else break;
            }
            if(j == L2rel[r_L+1]) {
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
                    pair_arry[p_num].nxtpnt = cap_row;
                    pair_arry[p_num].nxtflg = cap_flg;
                    pair_arry[p_num].idx_bg = 0;
                    pair_arry[p_num].idx_num = 0;
                    pair_arry[p_num].err = err_val + seqL[i][1];
                    pair_arry[p_num].l_off = seqL[i][3];
                    pair_arry[p_num].r_off = sub->r_off;
                    p_num++;  
                } 
                continue;
            }
            top_num = j;
            for(j=L2rel[r_L+1]-1; j >= top_num; --j){
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg == 1) continue;
                else break;
            }
            if(j < top_num) {
                printf("%u, error, j = %d, top_num = %d\n", __LINE__, j, top_num);
                exit(1);
            } else {
               bot_num = L2rel[r_L+1] -1 -j;  
            }  
            r_Rel = relat[top_num];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg <= IS_SMLSIZ) {
                bg_idx = cap_row; 
            } else {
                ext_idx = cB->head_extidx + cap_row;
                bg_idx = ext_idx[0][0]; 
            }
            bg_idx -= top_num - L2rel[r_L]; 
            j = L2rel[r_L+1]-1-bot_num; 
            r_Rel = relat[j];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            uint32_t cap_num;
            if(cap_flg <= IS_SMLSIZ) {
                ed_idx = cap_row; 
                cap_num = cap_flg;
            } else {
                ext_idx = cB->head_extidx + cap_row;
                ed_idx = ext_idx[0][0];
                cap_num = ext_idx[0][1]; 
            }
            ed_idx += cap_num - 1 + bot_num; 
            pair_arry[p_num].nxtpnt = 0;
            pair_arry[p_num].nxtflg = 0;
            pair_arry[p_num].idx_bg = bg_idx;
            pair_arry[p_num].idx_num = ed_idx + 1 - bg_idx;
            pair_arry[p_num].err = err_val + seqL[i][1];
            pair_arry[p_num].l_off = seqL[i][3];
            pair_arry[p_num].r_off = sub->r_off;
            p_num++;
            */
        } // end for(i=st_L; i<=numL; i++)++++++++++++++++
    }
    if(flag == 1) {
        /*  
        cap_row = sTree->stck_arry[len_arry].nxtpnt;
        
        if(sTree->cls == 0) {
            bg_idx = cB->bg_idx;
        } else {
            ext_idx = cB->head_extidx + cap_row; 
            bg_idx = ext_idx[0][0];
        }
        */
        uint32_t bg_idx;

fprintf(stderr, "%u, st_R = %d, numR = %d\n", __LINE__, st_R, numR);
//exit(1);
        /*  
        for(i = 0; i < cB->num_relat; ++i) {
            uint32_t pnt = nxtpnt[i]; 
            uint8_t flg = nxtflg[i];
            if(flg == 1) {
                pnt = bwt_get_idx(fm_idx->bwt, pnt);
            } 
fprintf(stderr, "i = %d, nxtpnt = %u, nxtflg = %u\n", i, pnt, flg);
        }
        */

        for(i=st_R; i<=numR; i++){
            max_err = sub->query_err + sub->delta;
            err_val = sub->err_sum[0] + seqR[i][1];
//printf("%u, st_R = %d, numR = %d, err_val = %d, max_err = %d\n", __LINE__, st_R, numR, err_val, max_err);
            if(err_val > max_err) { continue; } 
            r_R = R2rel[seqR[i][0]];
            r_Rel = relat[r_R]; 
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                pos = cap_row;
                num = cap_flg; 
            } else if(cap_flg <= IS_SMLSIZ) {
                pos = bwt_sa(fm_idx->bwt, cap_row); 
                num = cap_flg;
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                pos = bwt_sa(fm_idx->bwt, bg_idx);
                num = ext_idx[0][1];
            }
            bgn = bwt_get_idx(fm_idx->bwt, pos + 16); 
            //----- 
            r_R = R2rel[seqR[i][0]+1]-1;
            r_Rel = relat[r_R]; 
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                pos = cap_row;
                num = cap_flg; 
            } else if(cap_flg <= IS_SMLSIZ) {
                num = cap_flg;
                pos = bwt_sa(fm_idx->bwt, cap_row+num-1); 
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                num = ext_idx[0][1];
                pos = bwt_sa(fm_idx->bwt, bg_idx+num-1);
            }
            end = bwt_get_idx(fm_idx->bwt, pos + 16); 
fprintf(stderr, "%u, bgn = %u, end = %u\n", __LINE__, bgn, end);
            if(bgn > end) { continue; }
            num = end + 1 -bgn; 
            if(num > 1) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;  
            } else {
                pair_arry[p_num].idx_bg = 0;
                pair_arry[p_num].idx_num = 0; 
                pair_arry[p_num].nxtpnt = nxtpnt[seqR[i][0]]+16;
                pair_arry[p_num].nxtflg = 1;  
            }
            pair_arry[p_num].err = err_val + seqR[i][1];
            pair_arry[p_num].l_off = sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3];
            p_num++;                   
        }
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}
int PairExtSeq_32_uni(idx_t *fm_idx, struct ExtBlck *cB, uint32_t *relat, struct StackTree *sTree, int flag, struct SubBuf *sub){


    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint32_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + cB->num_relat; 
    R2rel = L2rel + cB->num_seqL+1; 
    //R2idx = R2rel + cB->num_seqR+1; 
    nxtpnt= cB->head_nxtpnt + cB->nxtpnt; 
	nxtflg= cB->head_nxtflg + cB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
    int st_L = 1;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 1;
    int ed_R = sub->seqR_out[0][0];
 
    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, rror!!!\n", __LINE__, __func__);
        exit(1); 
    }
fprintf(stderr, "%u, flag = %d\n", __LINE__, flag);
    int len_arry = sTree->len_arry;
    uint32_t (*ext_idx)[2];
    uint32_t cls = sTree->stck_arry[len_arry].cls;
    uint32_t cap_row, cap_flg;
    int max_err = sub->query_err + sub->delta;
    int err_val = sub->err_sum[0];
    uint32_t bgn, end, num, pos;
    if(flag == 0) {
        uint32_t cur_bg, cur_ed, bg_cap, bg_nxt, cap_num, top_num = 0, bot_num = 0, bg_idx, ed_idx;
fprintf(stderr, "%u, numL = %d\n", __LINE__, numL);
        for(i=st_L; i<=numL; i++){
            max_err = sub->query_err + sub->delta;
            err_val = sub->err_sum[0] + seqL[i][1];
            if(err_val > max_err) { continue; } 
            r_L = seqL[i][0];
            bgn = (uint32_t)-1;
            end = 0;
            bg_cap = 0, bg_nxt = 0, cap_num = 0;
            
            //r_Rel = relat[L2rel[r_L]];
            r_Rel = L2rel[r_L];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                bgn = bwt_get_idx(fm_idx->bwt, cap_row); 
                num = 1;
            } else if(cap_flg <= IS_SMLSIZ) {
                bgn = cap_row;
                num = cap_flg; 
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                pos = bwt_sa(fm_idx->bwt, bg_idx);
                num = ext_idx[0][1];
            }
            
            r_Rel = L2rel[r_L+1]-1;
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                end = bwt_get_idx(fm_idx->bwt, cap_row); 
                num = 1;
            } else if(cap_flg <= IS_SMLSIZ) {
                num = cap_flg; 
                end = cap_row + num - 1;
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                num = ext_idx[0][1];
                end = ext_idx[0][0] + ext_idx[0][1] -1; 
            }
            if(bgn > end) continue;
            num = end + 1 - bgn;
            if(num > 1) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;  
            } else {
                pair_arry[p_num].idx_bg = 0;
                pair_arry[p_num].idx_num = 0; 
                pair_arry[p_num].nxtpnt = nxtpnt[seqR[i][0]]+16;
                pair_arry[p_num].nxtflg = 1;  
            }
            pair_arry[p_num].err = err_val + seqR[i][1];
            pair_arry[p_num].l_off = sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3];
            p_num++;  
            
            
            /*   
            for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg == 1) continue;
                else break;
            }
            if(j == L2rel[r_L+1]) {
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
                    pair_arry[p_num].nxtpnt = cap_row;
                    pair_arry[p_num].nxtflg = cap_flg;
                    pair_arry[p_num].idx_bg = 0;
                    pair_arry[p_num].idx_num = 0;
                    pair_arry[p_num].err = err_val + seqL[i][1];
                    pair_arry[p_num].l_off = seqL[i][3];
                    pair_arry[p_num].r_off = sub->r_off;
                    p_num++;  
                } 
                continue;
            }
            top_num = j;
            for(j=L2rel[r_L+1]-1; j >= top_num; --j){
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg == 1) continue;
                else break;
            }
            if(j < top_num) {
                printf("%u, error, j = %d, top_num = %d\n", __LINE__, j, top_num);
                exit(1);
            } else {
               bot_num = L2rel[r_L+1] -1 -j;  
            }  
            r_Rel = relat[top_num];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg <= IS_SMLSIZ) {
                bg_idx = cap_row; 
            } else {
                ext_idx = cB->head_extidx + cap_row;
                bg_idx = ext_idx[0][0]; 
            }
            bg_idx -= top_num - L2rel[r_L]; 
            j = L2rel[r_L+1]-1-bot_num; 
            r_Rel = relat[j];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            uint32_t cap_num;
            if(cap_flg <= IS_SMLSIZ) {
                ed_idx = cap_row; 
                cap_num = cap_flg;
            } else {
                ext_idx = cB->head_extidx + cap_row;
                ed_idx = ext_idx[0][0];
                cap_num = ext_idx[0][1]; 
            }
            ed_idx += cap_num - 1 + bot_num; 
            pair_arry[p_num].nxtpnt = 0;
            pair_arry[p_num].nxtflg = 0;
            pair_arry[p_num].idx_bg = bg_idx;
            pair_arry[p_num].idx_num = ed_idx + 1 - bg_idx;
            pair_arry[p_num].err = err_val + seqL[i][1];
            pair_arry[p_num].l_off = seqL[i][3];
            pair_arry[p_num].r_off = sub->r_off;
            p_num++;
            */
        } // end for(i=st_L; i<=numL; i++)++++++++++++++++
    }
    if(flag == 1) {
        /*  
        cap_row = sTree->stck_arry[len_arry].nxtpnt;
       
        if(sTree->cls == 0) {
            bg_idx = cB->bg_idx;
        } else {
            ext_idx = cB->head_extidx + cap_row; 
            bg_idx = ext_idx[0][0];
        }
        */
        uint32_t bg_idx;
fprintf(stderr, "%u, st_R = %d, numR = %d\n", __LINE__, st_R, numR);
//exit(1);
        for(i=st_R; i<=numR; i++){
            max_err = sub->query_err + sub->delta;
            err_val = sub->err_sum[0] + seqR[i][1];
//printf("%u, st_R = %d, numR = %d, err_val = %d, max_err = %d\n", __LINE__, st_R, numR, err_val, max_err);
            if(err_val > max_err) { continue; } 
            
            r_R = R2rel[seqR[i][0]];
            r_Rel = relat[r_R]; 
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
fprintf(stderr, "%u, p_num = %d, cap_row = %d, cap_flg = %d, r_Rel = %d\n", __LINE__, p_num, cap_row, cap_flg, r_Rel);
            if(cap_flg == 1) {
                pos = cap_row;
                num = cap_flg; 
            } else if(cap_flg <= IS_SMLSIZ) {
                pos = bwt_sa(fm_idx->bwt, cap_row); 
                num = cap_flg;
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                pos = bwt_sa(fm_idx->bwt, bg_idx);
                num = ext_idx[0][1];
            }
            bgn = bwt_get_idx(fm_idx->bwt, pos + 16); 
            //----- 
            r_R = R2rel[seqR[i][0]+1]-1;
            r_Rel = relat[r_R]; 
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
fprintf(stderr, "%u, p_num = %d, cap_row = %d, cap_flg = %d, r_Rel = %d\n", __LINE__, p_num, cap_row, cap_flg, r_Rel);
            if(cap_flg == 1) {
                pos = cap_row;
                num = cap_flg; 
            } else if(cap_flg <= IS_SMLSIZ) {
                num = cap_flg;
                pos = bwt_sa(fm_idx->bwt, cap_row+num-1); 
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                num = ext_idx[0][1];
                pos = bwt_sa(fm_idx->bwt, bg_idx+num-1);
            }
            end = bwt_get_idx(fm_idx->bwt, pos + 16); 

            if(bgn > end) { continue; }
            num = end + 1 -bgn; 
fprintf(stderr, "%u, bgn = %u, end = %u, num = %u, p_num = %d\n", __LINE__, bgn, end, num, p_num); 
            if(num > 1) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;  
            } else {
                pair_arry[p_num].idx_bg = 0;
                pair_arry[p_num].idx_num = 0; 
                pair_arry[p_num].nxtpnt = nxtpnt[seqR[i][0]]+16;
                pair_arry[p_num].nxtflg = 1;  
            }
            pair_arry[p_num].err = err_val + seqR[i][1];
            pair_arry[p_num].l_off = sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3];
            p_num++;                   
        }
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}

int PairExtSeq_32_uni_2(idx_t *fm_idx, struct ExtBlck *cB, uint32_t *relat, struct StackTree *sTree, int flag, struct SubBuf *sub){


    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint32_t *L2rel, *R2rel, *R2idx;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + cB->num_relat; 
    R2rel = L2rel + cB->num_seqL+1; 
    R2idx = R2rel + cB->num_seqR+1; 
    nxtpnt= cB->head_nxtpnt + cB->nxtpnt; 
	nxtflg= cB->head_nxtflg + cB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
    int st_L = 1;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 1;
    int ed_R = sub->seqR_out[0][0];
 
    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, rror!!!\n", __LINE__, __func__);
        exit(1); 
    }
fprintf(stderr, "%u, flag = %d\n", __LINE__, flag);
    int len_arry = sTree->len_arry;
    uint32_t (*ext_idx)[2];
    uint32_t cls = sTree->stck_arry[len_arry].cls;
    uint32_t cap_row, cap_flg;
    int max_err = sub->query_err + sub->delta;
    int err_val = sub->err_sum[0];
    uint32_t bgn, end, num, pos;
    if(flag == 0) {
        uint32_t cur_bg, cur_ed, bg_cap, bg_nxt, cap_num, top_num = 0, bot_num = 0, bg_idx, ed_idx;
fprintf(stderr, "%u, numL = %d\n", __LINE__, numL);
        for(i=st_L; i<=numL; i++){
            max_err = sub->query_err + sub->delta;
            err_val = sub->err_sum[0] + seqL[i][1];
            if(err_val < max_err) { continue; } 
            r_L = seqL[i][0];
            bgn = (uint32_t)-1;
            end = 0;
            bg_cap = 0, bg_nxt = 0, cap_num = 0;
            
            //r_Rel = relat[L2rel[r_L]];
            r_Rel = L2rel[r_L];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                bgn = bwt_get_idx(fm_idx->bwt, cap_row); 
                num = 1;
            } else if(cap_flg <= IS_SMLSIZ) {
                bgn = cap_row;
                num = cap_flg; 
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                pos = bwt_sa(fm_idx->bwt, bg_idx);
                num = ext_idx[0][1];
            }
            
            //r_Rel = relat[L2rel[r_L+1]-1];
            r_Rel = L2rel[r_L+1] -1;
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                end = bwt_get_idx(fm_idx->bwt, cap_row); 
                num = 1;
            } else if(cap_flg <= IS_SMLSIZ) {
                num = cap_flg; 
                end = cap_row + num - 1;
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                num = ext_idx[0][1];
                end = ext_idx[0][0] + ext_idx[0][1] -1; 
            }
            if(bgn > end) continue;
            num = end + 1 - bgn;
            if(num > 1) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;  
            } else {
                pair_arry[p_num].idx_bg = 0;
                pair_arry[p_num].idx_num = 0; 
                pair_arry[p_num].nxtpnt = nxtpnt[seqR[i][0]]+16;
                pair_arry[p_num].nxtflg = 1;  
            }
            pair_arry[p_num].err = err_val + seqR[i][1];
            pair_arry[p_num].l_off = sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3];
            p_num++;  
            
            
            /*   
            for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg == 1) continue;
                else break;
            }
            if(j == L2rel[r_L+1]) {
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
                    pair_arry[p_num].nxtpnt = cap_row;
                    pair_arry[p_num].nxtflg = cap_flg;
                    pair_arry[p_num].idx_bg = 0;
                    pair_arry[p_num].idx_num = 0;
                    pair_arry[p_num].err = err_val + seqL[i][1];
                    pair_arry[p_num].l_off = seqL[i][3];
                    pair_arry[p_num].r_off = sub->r_off;
                    p_num++;  
                } 
                continue;
            }
            top_num = j;
            for(j=L2rel[r_L+1]-1; j >= top_num; --j){
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg == 1) continue;
                else break;
            }
            if(j < top_num) {
                printf("%u, error, j = %d, top_num = %d\n", __LINE__, j, top_num);
                exit(1);
            } else {
               bot_num = L2rel[r_L+1] -1 -j;  
            }  
            r_Rel = relat[top_num];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg <= IS_SMLSIZ) {
                bg_idx = cap_row; 
            } else {
                ext_idx = cB->head_extidx + cap_row;
                bg_idx = ext_idx[0][0]; 
            }
            bg_idx -= top_num - L2rel[r_L]; 
            j = L2rel[r_L+1]-1-bot_num; 
            r_Rel = relat[j];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            uint32_t cap_num;
            if(cap_flg <= IS_SMLSIZ) {
                ed_idx = cap_row; 
                cap_num = cap_flg;
            } else {
                ext_idx = cB->head_extidx + cap_row;
                ed_idx = ext_idx[0][0];
                cap_num = ext_idx[0][1]; 
            }
            ed_idx += cap_num - 1 + bot_num; 
            pair_arry[p_num].nxtpnt = 0;
            pair_arry[p_num].nxtflg = 0;
            pair_arry[p_num].idx_bg = bg_idx;
            pair_arry[p_num].idx_num = ed_idx + 1 - bg_idx;
            pair_arry[p_num].err = err_val + seqL[i][1];
            pair_arry[p_num].l_off = seqL[i][3];
            pair_arry[p_num].r_off = sub->r_off;
            p_num++;
            */
        } // end for(i=st_L; i<=numL; i++)++++++++++++++++
    }
    if(flag == 1) {
        cap_row = sTree->stck_arry[len_arry].nxtpnt;
        uint32_t bg_idx;
        if(sTree->cls == 0) {
            bg_idx = cB->bg_idx;
        } else {
            ext_idx = cB->head_extidx + cap_row; 
            bg_idx = ext_idx[0][0];
        }
        for(i=st_R; i<=numR; i++){
            max_err = sub->query_err + sub->delta;
            err_val = sub->err_sum[0] + seqR[i][1];
            if(err_val < max_err) { continue; } 
            r_Rel = seqR[i][0];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                pos = cap_row;
                num = cap_flg; 
            } else if(cap_flg <= IS_SMLSIZ) {
                pos = bwt_sa(fm_idx->bwt, cap_row); 
                num = cap_flg;
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                pos = bwt_sa(fm_idx->bwt, bg_idx);
                num = ext_idx[0][1];
            }
            bgn = bwt_get_idx(fm_idx->bwt, pos + 16); 
            
            r_Rel = seqR[i][0] + 1;
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                pos = cap_row;
                num = cap_flg; 
            } else if(cap_flg <= IS_SMLSIZ) {
                num = cap_flg;
                pos = bwt_sa(fm_idx->bwt, cap_row+num-1); 
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                num = ext_idx[0][1];
                pos = bwt_sa(fm_idx->bwt, bg_idx+num-1);
            }
            end = bwt_get_idx(fm_idx->bwt, pos + 16); 
            if(bgn > end) { continue; }
            num = end + 1 -bgn; 
            if(num > 1) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;  
            } else {
                pair_arry[p_num].idx_bg = 0;
                pair_arry[p_num].idx_num = 0; 
                pair_arry[p_num].nxtpnt = nxtpnt[seqR[i][0]]+16;
                pair_arry[p_num].nxtflg = 1;  
            }
            pair_arry[p_num].err = err_val + seqR[i][1];
            pair_arry[p_num].l_off = sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3];
            p_num++;                   
        }
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}


int PairExtSeq_16_uni_2(idx_t *fm_idx, struct ExtBlck *cB, uint16_t *relat, struct StackTree *sTree, int flag, struct SubBuf *sub){


    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint16_t *L2rel, *R2rel, *R2idx;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + cB->num_relat; 
    R2rel = L2rel + cB->num_seqL+1; 
    R2idx = R2rel + cB->num_seqR+1; 
    nxtpnt= cB->head_nxtpnt + cB->nxtpnt; 
	nxtflg= cB->head_nxtflg + cB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
    int st_L = 1;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 1;
    int ed_R = sub->seqR_out[0][0];
 
    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, rror!!!\n", __LINE__, __func__);
        exit(1); 
    }
fprintf(stderr, "%u, flag = %d\n", __LINE__, flag);
    int len_arry = sTree->len_arry;
    uint32_t (*ext_idx)[2];
    uint32_t cls = sTree->stck_arry[len_arry].cls;
    uint32_t cap_row, cap_flg;
    int max_err = sub->query_err + sub->delta;
    int err_val = sub->err_sum[0];
    uint32_t bgn, end, num, pos;
    if(flag == 0) {
        uint32_t cur_bg, cur_ed, bg_cap, bg_nxt, cap_num, top_num = 0, bot_num = 0, bg_idx, ed_idx;
fprintf(stderr, "%u, numL = %d\n", __LINE__, numL);
        for(i=st_L; i<=numL; i++){
            max_err = sub->query_err + sub->delta;
            err_val = sub->err_sum[0] + seqL[i][1];
            if(err_val < max_err) { continue; } 
            r_L = seqL[i][0];
            bgn = (uint32_t)-1;
            end = 0;
            bg_cap = 0, bg_nxt = 0, cap_num = 0;
            
            r_Rel = relat[L2rel[r_L]];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                bgn = bwt_get_idx(fm_idx->bwt, cap_row); 
                num = 1;
            } else if(cap_flg <= IS_SMLSIZ) {
                bgn = cap_row;
                num = cap_flg; 
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                pos = bwt_sa(fm_idx->bwt, bg_idx);
                num = ext_idx[0][1];
            }
            
            r_Rel = relat[L2rel[r_L+1]-1];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                end = bwt_get_idx(fm_idx->bwt, cap_row); 
                num = 1;
            } else if(cap_flg <= IS_SMLSIZ) {
                num = cap_flg; 
                end = cap_row + num - 1;
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                num = ext_idx[0][1];
                end = ext_idx[0][0] + ext_idx[0][1] -1; 
            }
            if(bgn > end) continue;
            num = end + 1 - bgn;
            if(num > 1) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;  
            } else {
                pair_arry[p_num].idx_bg = 0;
                pair_arry[p_num].idx_num = 0; 
                pair_arry[p_num].nxtpnt = nxtpnt[seqR[i][0]]+16;
                pair_arry[p_num].nxtflg = 1;  
            }
            pair_arry[p_num].err = err_val + seqR[i][1];
            pair_arry[p_num].l_off = sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3];
            p_num++;  
            
            
            /*   
            for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg == 1) continue;
                else break;
            }
            if(j == L2rel[r_L+1]) {
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
                    pair_arry[p_num].nxtpnt = cap_row;
                    pair_arry[p_num].nxtflg = cap_flg;
                    pair_arry[p_num].idx_bg = 0;
                    pair_arry[p_num].idx_num = 0;
                    pair_arry[p_num].err = err_val + seqL[i][1];
                    pair_arry[p_num].l_off = seqL[i][3];
                    pair_arry[p_num].r_off = sub->r_off;
                    p_num++;  
                } 
                continue;
            }
            top_num = j;
            for(j=L2rel[r_L+1]-1; j >= top_num; --j){
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg == 1) continue;
                else break;
            }
            if(j < top_num) {
                printf("%u, error, j = %d, top_num = %d\n", __LINE__, j, top_num);
                exit(1);
            } else {
               bot_num = L2rel[r_L+1] -1 -j;  
            }  
            r_Rel = relat[top_num];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg <= IS_SMLSIZ) {
                bg_idx = cap_row; 
            } else {
                ext_idx = cB->head_extidx + cap_row;
                bg_idx = ext_idx[0][0]; 
            }
            bg_idx -= top_num - L2rel[r_L]; 
            j = L2rel[r_L+1]-1-bot_num; 
            r_Rel = relat[j];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            uint32_t cap_num;
            if(cap_flg <= IS_SMLSIZ) {
                ed_idx = cap_row; 
                cap_num = cap_flg;
            } else {
                ext_idx = cB->head_extidx + cap_row;
                ed_idx = ext_idx[0][0];
                cap_num = ext_idx[0][1]; 
            }
            ed_idx += cap_num - 1 + bot_num; 
            pair_arry[p_num].nxtpnt = 0;
            pair_arry[p_num].nxtflg = 0;
            pair_arry[p_num].idx_bg = bg_idx;
            pair_arry[p_num].idx_num = ed_idx + 1 - bg_idx;
            pair_arry[p_num].err = err_val + seqL[i][1];
            pair_arry[p_num].l_off = seqL[i][3];
            pair_arry[p_num].r_off = sub->r_off;
            p_num++;
            */
        } // end for(i=st_L; i<=numL; i++)++++++++++++++++
    }
    if(flag == 1) {
        cap_row = sTree->stck_arry[len_arry].nxtpnt;
        uint32_t bg_idx;
        if(sTree->cls == 0) {
            bg_idx = cB->bg_idx;
        } else {
            ext_idx = cB->head_extidx + cap_row; 
            bg_idx = ext_idx[0][0];
        }
        for(i=st_R; i<=numR; i++){
            max_err = sub->query_err + sub->delta;
            err_val = sub->err_sum[0] + seqR[i][1];
            if(err_val < max_err) { continue; } 
            r_Rel = seqR[i][0];
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                pos = cap_row;
                num = cap_flg; 
            } else if(cap_flg <= IS_SMLSIZ) {
                pos = bwt_sa(fm_idx->bwt, cap_row); 
                num = cap_flg;
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                pos = bwt_sa(fm_idx->bwt, bg_idx);
                num = ext_idx[0][1];
            }
            bgn = bwt_get_idx(fm_idx->bwt, pos + 16); 
            
            r_Rel = seqR[i][0] + 1;
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
            if(cap_flg == 1) {
                pos = cap_row;
                num = cap_flg; 
            } else if(cap_flg <= IS_SMLSIZ) {
                num = cap_flg;
                pos = bwt_sa(fm_idx->bwt, cap_row+num-1); 
            } else {
                ext_idx = cB->head_extidx + cap_row; 
                bg_idx = ext_idx[0][0];
                num = ext_idx[0][1];
                pos = bwt_sa(fm_idx->bwt, bg_idx+num-1);
            }
            end = bwt_get_idx(fm_idx->bwt, pos + 16); 
            if(bgn > end) { continue; }
            num = end + 1 -bgn; 
            if(num > 1) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;  
            } else {
                pair_arry[p_num].idx_bg = 0;
                pair_arry[p_num].idx_num = 0; 
                pair_arry[p_num].nxtpnt = nxtpnt[seqR[i][0]]+16;
                pair_arry[p_num].nxtflg = 1;  
            }
            pair_arry[p_num].err = err_val + seqR[i][1];
            pair_arry[p_num].l_off = sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3];
            p_num++;                   
        }
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}

int PairExtSeq_8_uni_1( struct ExtBlck *cB, uint8_t *relat, struct StackTree *sTree, int flag, struct SubBuf *sub){
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint8_t *L2rel, *R2rel, *R2idx;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + cB->num_relat; 
    R2rel = L2rel + cB->num_seqL+1; 
    R2idx = R2rel + cB->num_seqR+1; 
    nxtpnt= cB->head_nxtpnt + cB->nxtpnt; 
	nxtflg= cB->head_nxtflg + cB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
    /*  
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    */
    int st_L = 1;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 1;
    int ed_R = sub->seqR_out[0][0];
 
    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, rror!!!\n", __LINE__, __func__);
        exit(1); 
    }
fprintf(stderr, "%u, flag = %d\n", __LINE__, flag);
    int len_arry = sTree->len_arry;
    uint32_t (*ext_idx)[2];
    uint32_t cls = sTree->stck_arry[len_arry].cls;
    uint32_t cap_row, cap_flg;
    int max_err = sub->query_err + sub->delta;
    int err_val = sub->err_sum[0];
    uint32_t bgn, end, num;
    if(flag == 0) {

        for(r_Rel=0; r_Rel<cB->num_relat; r_Rel++){
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
if(cap_flg > 1) fprintf(stderr, "%u, r_Rel = %d, cap_row = %u, cap_flg = %d\n", __LINE__, r_Rel, cap_row, cap_flg);
        }

        uint32_t cur_bg, cur_ed, bg_cap, bg_nxt, cap_num, top_num = 0, bot_num = 0, bg_idx, ed_idx;
fprintf(stderr, "%u, numL = %d\n", __LINE__, numL);
        for(i=st_L; i<=numL; i++){
            r_L = seqL[i][0];
            //if(seqL[i][1] > 0) err_val += seqL[i][1];
            bgn = (uint32_t)-1;
            end = 0;
fprintf(stderr, "%u, i = %d, L2rel[%d] = %u, L2rel[%d+1] = %u\n", __LINE__, i, r_L, L2rel[r_L], r_L, L2rel[r_L+1]);
            if(err_val < max_err) {
                bg_cap = 0, bg_nxt = 0, cap_num = 0;
                
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
if(cap_flg > 1) fprintf(stderr, "%u, r_Rel = %d, cap_row = %u, cap_flg = %d\n", __LINE__, r_Rel, cap_row, cap_flg);
                }

                
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
                    if(cap_flg == 1) continue;
                    else break;
                }
                if(j == L2rel[r_L+1]) {
                    for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                        r_Rel = relat[j];
                        cap_row = nxtpnt[r_Rel];
                        cap_flg = nxtflg[r_Rel];

                        pair_arry[p_num].nxtpnt = cap_row;
                        pair_arry[p_num].nxtflg = cap_flg;
                        pair_arry[p_num].idx_bg = 0;
                        pair_arry[p_num].idx_num = 0;
                        pair_arry[p_num].err = err_val + seqL[i][1];
                        pair_arry[p_num].l_off = seqL[i][3];
                        pair_arry[p_num].r_off = sub->r_off;
                        p_num++;  
                    } 
                    continue;
                }
                top_num = j;
                for(j=L2rel[r_L+1]-1; j >= top_num; --j){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
                    if(cap_flg == 1) continue;
                    else break;
                }
                if(j < top_num) {
                    printf("%u, error, j = %d, top_num = %d\n", __LINE__, j, top_num);
                    exit(1);
                } else {
                   bot_num = L2rel[r_L+1] -1 -j;  
                }  
                r_Rel = relat[top_num];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg <= IS_SMLSIZ) {
                    bg_idx = cap_row; 
                } else {
                    ext_idx = cB->head_extidx + cap_row;
                    bg_idx = ext_idx[0][0]; 
                }
                bg_idx -= top_num - L2rel[r_L]; 
                
                j = L2rel[r_L+1]-1-bot_num; 
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                uint32_t cap_num;
                if(cap_flg <= IS_SMLSIZ) {
                    ed_idx = cap_row; 
                    cap_num = cap_flg;
                } else {
                    ext_idx = cB->head_extidx + cap_row;
                    ed_idx = ext_idx[0][0];
                    cap_num = ext_idx[0][1]; 
                }
                ed_idx += cap_num - 1 + bot_num; 
fprintf(stderr, "%u, cap_num = %u, top_num = %u, bot_num = %u, bg_idx = %u, ed_idx = %u, cap_row = %u, cap_flg = %u\n", __LINE__, cap_num, top_num, bot_num, bg_idx, ed_idx, cap_row, cap_flg);
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;
                pair_arry[p_num].idx_bg = bg_idx;
                pair_arry[p_num].idx_num = ed_idx + 1 - bg_idx;
                pair_arry[p_num].err = err_val + seqL[i][1];
                pair_arry[p_num].l_off = seqL[i][3];
                pair_arry[p_num].r_off = sub->r_off;
                p_num++;
            } //if(err_val < max_err)+++++++++++++++++++++
        } // end for(i=st_L; i<=numL; i++)++++++++++++++++
    }
    if(flag == 1) {
        cap_row = sTree->stck_arry[len_arry].nxtpnt;
        //ext_idx = (cB+cls)->head_extidx + cap_row;
        uint32_t bg_idx;
        if(sTree->cls == 0) {
            bg_idx = cB->bg_idx;
        } else {
            ext_idx = cB->head_extidx + cap_row; 
            bg_idx = ext_idx[0][0];
        }
        for(i=st_R; i<=numR; i++){
            r_R = seqR[i][0];
            bgn = bg_idx + R2idx[r_R];
            num = R2idx[r_R+1] - R2idx[r_R];
            if(err_val < max_err) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].err = err_val + seqR[i][1];
                pair_arry[p_num].l_off = sub->l_off;
                pair_arry[p_num].r_off = seqR[i][3];
                p_num++;                   
            }
        }
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}


int PairExtSeq_32_uni_1( struct ExtBlck *cB, uint32_t *relat, struct StackTree *sTree, int flag, struct SubBuf *sub){
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint32_t *L2rel, *R2rel, *R2idx;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + cB->num_relat; 
    R2rel = L2rel + cB->num_seqL+1; 
    R2idx = R2rel + cB->num_seqR+1; 
    nxtpnt= cB->head_nxtpnt + cB->nxtpnt; 
	nxtflg= cB->head_nxtflg + cB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
    /*  
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    */
    int st_L = 1;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 1;
    int ed_R = sub->seqR_out[0][0];
 
    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, rror!!!\n", __LINE__, __func__);
        exit(1); 
    }
fprintf(stderr, "%u, flag = %d\n", __LINE__, flag);
    int len_arry = sTree->len_arry;
    uint32_t (*ext_idx)[2];
    uint32_t cls = sTree->stck_arry[len_arry].cls;
    uint32_t cap_row, cap_flg;
    int max_err = sub->query_err + sub->delta;
    int err_val = sub->err_sum[0];
    uint32_t bgn, end, num;
    if(flag == 0) {

        for(r_Rel=0; r_Rel<cB->num_relat; r_Rel++){
            cap_row = nxtpnt[r_Rel];
            cap_flg = nxtflg[r_Rel];
if(cap_flg > 1) fprintf(stderr, "%u, r_Rel = %d, cap_row = %u, cap_flg = %d\n", __LINE__, r_Rel, cap_row, cap_flg);
        }

        uint32_t cur_bg, cur_ed, bg_cap, bg_nxt, cap_num, top_num = 0, bot_num = 0, bg_idx, ed_idx;
fprintf(stderr, "%u, numL = %d\n", __LINE__, numL);
        for(i=st_L; i<=numL; i++){
            r_L = seqL[i][0];
            //if(seqL[i][1] > 0) err_val += seqL[i][1];
            bgn = (uint32_t)-1;
            end = 0;
fprintf(stderr, "%u, i = %d, L2rel[%d] = %u, L2rel[%d+1] = %u\n", __LINE__, i, r_L, L2rel[r_L], r_L, L2rel[r_L+1]);
            if(err_val < max_err) {
                bg_cap = 0, bg_nxt = 0, cap_num = 0;
                
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
if(cap_flg > 1) fprintf(stderr, "%u, r_Rel = %d, cap_row = %u, cap_flg = %d\n", __LINE__, r_Rel, cap_row, cap_flg);
                }

                
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
                    if(cap_flg == 1) continue;
                    else break;
                }
                if(j == L2rel[r_L+1]) {
                    for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                        r_Rel = relat[j];
                        cap_row = nxtpnt[r_Rel];
                        cap_flg = nxtflg[r_Rel];

                        pair_arry[p_num].nxtpnt = cap_row;
                        pair_arry[p_num].nxtflg = cap_flg;
                        pair_arry[p_num].idx_bg = 0;
                        pair_arry[p_num].idx_num = 0;
                        pair_arry[p_num].err = err_val + seqL[i][1];
                        pair_arry[p_num].l_off = seqL[i][3];
                        pair_arry[p_num].r_off = sub->r_off;
                        p_num++;  
                    } 
                    continue;
                }
                top_num = j;
                for(j=L2rel[r_L+1]-1; j >= top_num; --j){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    cap_flg = nxtflg[r_Rel];
                    if(cap_flg == 1) continue;
                    else break;
                }
                if(j < top_num) {
                    printf("%u, error, j = %d, top_num = %d\n", __LINE__, j, top_num);
                    exit(1);
                } else {
                   bot_num = L2rel[r_L+1] -1 -j;  
                }  
                r_Rel = relat[top_num];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                if(cap_flg <= IS_SMLSIZ) {
                    bg_idx = cap_row; 
                } else {
                    ext_idx = cB->head_extidx + cap_row;
                    bg_idx = ext_idx[0][0]; 
                }
                bg_idx -= top_num - L2rel[r_L]; 
                
                j = L2rel[r_L+1]-1-bot_num; 
                r_Rel = relat[j];
                cap_row = nxtpnt[r_Rel];
                cap_flg = nxtflg[r_Rel];
                uint32_t cap_num;
                if(cap_flg <= IS_SMLSIZ) {
                    ed_idx = cap_row; 
                    cap_num = cap_flg;
                } else {
                    ext_idx = cB->head_extidx + cap_row;
                    ed_idx = ext_idx[0][0];
                    cap_num = ext_idx[0][1]; 
                }
                ed_idx += cap_num - 1 + bot_num; 
fprintf(stderr, "%u, cap_num = %u, top_num = %u, bot_num = %u, bg_idx = %u, ed_idx = %u, cap_row = %u, cap_flg = %u\n", __LINE__, cap_num, top_num, bot_num, bg_idx, ed_idx, cap_row, cap_flg);
                pair_arry[p_num].nxtpnt = 0;
                pair_arry[p_num].nxtflg = 0;
                pair_arry[p_num].idx_bg = bg_idx;
                pair_arry[p_num].idx_num = ed_idx + 1 - bg_idx;
                pair_arry[p_num].err = err_val + seqL[i][1];
                pair_arry[p_num].l_off = seqL[i][3];
                pair_arry[p_num].r_off = sub->r_off;
                p_num++;
            } //if(err_val < max_err)+++++++++++++++++++++
        } // end for(i=st_L; i<=numL; i++)++++++++++++++++
    }
    if(flag == 1) {
        cap_row = sTree->stck_arry[len_arry].nxtpnt;
        //ext_idx = (cB+cls)->head_extidx + cap_row;
        uint32_t bg_idx;
        if(sTree->cls == 0) {
            bg_idx = cB->bg_idx;
        } else {
            ext_idx = cB->head_extidx + cap_row; 
            bg_idx = ext_idx[0][0];
        }
        for(i=st_R; i<=numR; i++){
            r_R = seqR[i][0];
            bgn = bg_idx + R2idx[r_R];
            num = R2idx[r_R+1] - R2idx[r_R];
            if(err_val < max_err) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].err = err_val + seqR[i][1];
                pair_arry[p_num].l_off = sub->l_off;
                pair_arry[p_num].r_off = seqR[i][3];
                p_num++;                   
            }
        }
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}


/*  
int PairExtSeq_16_uni(int flag, struct ExtBlck *eB, struct SubBuf *sub, uint16_t *relat){
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint16_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + eB->num_relat; 
    R2rel = L2rel + eB->num_seqL+1; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
    
    int st_L = 0;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 0;
    int ed_R = sub->seqR_out[0][0];
 

    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, flag = %d, error!!!\n", __LINE__, __func__, flag);
        exit(1); 
    }
    if(flag == 0) {
        for(i=st_L; i<=numL; i++){
            r_L = seqL[i][0];
            for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                r_Rel = relat[j];
                int err_val = sub->err_sum[0];
                if(seqL[i][1] > 0) err_val += seqL[i][1];
                int max_err = sub->query_err + sub->delta;
                if(err_val < max_err) {
                    pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
                    pair_arry[p_num].nxtflg = nxtflg[r_Rel];
                    pair_arry[p_num].err = err_val;
                    pair_arry[p_num].l_off = seqL[i][3];
                    p_num++;                    
                }
             
            }  
        }
    }
    if(flag == 1) {
        for(i=st_R; i<=numR; i++){
            r_R = seqR[i][0];
            for(r_Rel=R2rel[r_R]; r_Rel<R2rel[r_R+1]; r_Rel++) {
                int err_val = sub->err_sum[0];
                if(seqR[i][1] > 0) err_val += seqR[i][1];
                int max_err = sub->query_err + sub->delta;
                if(err_val < max_err) {
                    pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
                    pair_arry[p_num].nxtflg = nxtflg[r_Rel];
                    pair_arry[p_num].err = err_val;
                    pair_arry[p_num].r_off = seqR[i][3];
                    p_num++;                   
                }
            }
        }
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}
*/
int PairExtSeq_8_uni0( struct ExtBlck *cB, uint8_t *relat, struct StackTree *sTree, int flag, struct SubBuf *sub){
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint8_t *L2rel, *R2rel, *R2idx;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + cB->num_relat; 
    R2rel = L2rel + cB->num_seqL+1; 
    R2idx = R2rel + cB->num_seqR+1; 
	nxtpnt= cB->head_nxtpnt + cB->nxtpnt; 
	nxtflg= cB->head_nxtflg + cB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
    /*  
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    */
    int st_L = 1;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 1;
    int ed_R = sub->seqR_out[0][0];
 
    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, rror!!!\n", __LINE__, __func__);
        exit(1); 
    }
    int len_arry = sTree->len_arry;
    uint32_t (*ext_idx)[2];
    uint32_t cls = sTree->stck_arry[len_arry].cls;
    uint32_t cap_row;
    int max_err = sub->query_err + sub->delta;
    int err_val = sub->err_sum[0];
    uint32_t bgn, end, num;
    if(flag == 0) {
        uint32_t cur_bg, cur_ed;
        for(i=st_L; i<=numL; i++){
            r_L = seqL[i][0];
            //if(seqL[i][1] > 0) err_val += seqL[i][1];
            bgn = (uint32_t)-1;
            end = 0;
            if(err_val < max_err) {
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    ext_idx = (cB+1)->head_extidx + cap_row;
                    //ext_idx = cB->head_extidx + cap_row;
                    cur_bg = ext_idx[0][0];
                    cur_ed = ext_idx[0][0] + ext_idx[0][1] -1;
                    if(bgn > cur_bg) {
                        bgn = cur_bg;
                    }
                    if(end < cur_ed) {
                        end = cur_ed;
                    }
                } 
            
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = end+1-bgn;
                pair_arry[p_num].err = err_val + seqL[i][1];
                pair_arry[p_num].l_off = seqL[i][3];
                pair_arry[p_num].r_off = sub->r_off;
                p_num++;                    
            }
        }
    }
    if(flag == 1) {
        cap_row = sTree->stck_arry[len_arry].nxtpnt;
        ext_idx = cB->head_extidx + cap_row;
        //bgn = ext_idx[0][0];
        //num = ext_idx[0][1];
        for(i=st_R; i<=numR; i++){
            r_R = seqR[i][0];
            bgn = ext_idx[0][0] + R2idx[r_R];
            num = R2idx[r_R+1] - R2idx[r_R];
            //if(seqR[i][1] > 0) err_val += seqR[i][1];
   
            if(err_val < max_err) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].err = err_val + seqR[i][1];
                pair_arry[p_num].l_off = sub->l_off;
                pair_arry[p_num].r_off = seqR[i][3];
                p_num++;                   
            }
        }
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}
int PairExtSeq_32_3(struct ExtBlck *eB, struct SubBuf *sub, uint32_t *relat){
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
	uint32_t *seqRel, *seqRel_L, *seqRel_R;
    //uint32_t (*pair)[5];
    uint8_t *algnRel_row, *algnRel_buf;	
	uint32_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
    seqRel_L = sub->seqRel_L;
    seqRel_R = sub->seqRel_R;
    algnRel_row = sub->algnRel_row;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + eB->num_relat; 
    R2rel = L2rel + eB->num_seqL+1; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int i, j, k; 
    int r_L, r_R, r_Rel;
	int p_num;
    
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    int numL = ed_L; 
	int numR = ed_R;
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++) {
            algnRel_buf[j/8] |= (1<<(7-j%8));
        }
    } 
    
    p_num = sub->pair_out->p_num; 
    for(i=st_R; i<=numR; i++){
        r_R = seqR[i][0];
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++){
            r_Rel = relat[j];
            int flag_Rel = algnRel_row[r_Rel/8] & (1<<(7-r_Rel%8)); 
            if(flag_Rel > 0) continue;	
            int flag_R = algnRel_buf[r_Rel/8] & (1<<(7-r_Rel%8));
            if(flag_R>0){
                for(k = st_L; k <= numL; k++){
                    r_L = seqL[k][0];
                    if(r_Rel >= L2rel[r_L] &&r_Rel < L2rel[r_L+1]) break; 
                }
                if(k > numL){
                    printf("%u, %s algnR_row[] error!, k =%u, numL = %u\n", __LINE__, __func__, k, numL);
exit(1);
                    //continue;
                }
                int err_val = sub->err_sum[0];
                if(seqR[i][1] > 0) err_val += seqR[i][1];
                if(seqL[k][1] > 0) err_val += seqL[k][1];
                int max_err = sub->query_err + sub->delta;
                if(err_val < max_err) {
                    pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
                    pair_arry[p_num].nxtflg = nxtflg[r_Rel];
                    pair_arry[p_num].err = err_val;
                    /*  
                    pair_arry[p_num].l_off = seqL[i][3];
                    pair_arry[p_num].r_off = seqR[k][3];
                    */
                    pair_arry[p_num].l_off = seqL[k][3] + sub->l_off;
                    pair_arry[p_num].r_off = seqR[i][3] + sub->r_off;

                    
                    p_num++;                   
                    algnRel_row[r_Rel/8] |= (1<<(7-r_Rel%8));
                    seqRel[0]++;
                    seqRel[seqRel[0]] = r_Rel;
                }
            }
		}  
    }
    sub->pair_out->p_num = p_num; 
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++) {
            algnRel_buf[j/8] = 0;
        }
    } 
    return p_num;
}

int PairExtSeq_16_3(struct ExtBlck *eB, struct SubBuf *sub, uint16_t *relat){
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
	uint32_t *seqRel, *seqRel_L, *seqRel_R;
    //uint32_t (*pair)[5];
    uint8_t *algnRel_row, *algnRel_buf;	
	uint16_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
    seqRel_L = sub->seqRel_L;
    seqRel_R = sub->seqRel_R;
    algnRel_row = sub->algnRel_row;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + eB->num_relat; 
    R2rel = L2rel + eB->num_seqL+1; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int i, j, k; 
    int r_L, r_R, r_Rel;
	int p_num;
    
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    int numL = ed_L; 
	int numR = ed_R;
    /*  
    for(i=st_R; i<=numR; i++){
        r_R = seqR[i][0];
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++) {
            algnRel_buf[j/8] |= (1<<(7-j%8));
        }
    } 
    */
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++) {
            algnRel_buf[j/8] |= (1<<(7-j%8));
        }
    } 
    
    p_num = sub->pair_out->p_num; 
    for(i=st_R; i<=numR; i++){
        r_R = seqR[i][0];
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++){
            r_Rel = relat[j];
            int flag_Rel = algnRel_row[r_Rel/8] & (1<<(7-r_Rel%8)); 
            if(flag_Rel > 0) continue;	
            int flag_R = algnRel_buf[r_Rel/8] & (1<<(7-r_Rel%8));
            if(flag_R>0){
                for(k = st_L; k <= numL; k++){
                    r_L = seqL[k][0];
                    if(r_Rel >= L2rel[r_L] &&r_Rel < L2rel[r_L+1]) break; 
                }
                if(k > numL){
                    printf("%u, %s algnR_row[] error!, k =%u, numL = %u\n", __LINE__, __func__, k, numL);
exit(1);
                    //continue;
                }
                int err_val = sub->err_sum[0];
                if(seqR[i][1] > 0) err_val += seqR[i][1];
                if(seqL[k][1] > 0) err_val += seqL[k][1];
                int max_err = sub->query_err + sub->delta;
                if(err_val < max_err) {
                    pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
                    pair_arry[p_num].nxtflg = nxtflg[r_Rel];
                    pair_arry[p_num].err = err_val;
                    /*  
                    pair_arry[p_num].l_off = seqL[i][3];
                    pair_arry[p_num].r_off = seqR[k][3];
                    */
                    pair_arry[p_num].l_off = seqL[k][3] + sub->l_off;
                    pair_arry[p_num].r_off = seqR[i][3] + sub->r_off;

                    
                    p_num++;                   
                    algnRel_row[r_Rel/8] |= (1<<(7-r_Rel%8));
                    seqRel[0]++;
                    seqRel[seqRel[0]] = r_Rel;
                }
            }
		}  
    }
    sub->pair_out->p_num = p_num; 
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++) {
            algnRel_buf[j/8] = 0;
        }
    } 
    return p_num;
}
int PairExtSeq_32(struct ExtBlck *eB, struct SubBuf *sub, uint32_t *relat)
{
    int err_val = sub->err_sum[0];
    int max_err = sub->query_err + sub->delta;
    if(err_val >= max_err) return 0;
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
	uint32_t *seqRel, *seqRel_L, *seqRel_R;
    uint8_t *algnRel_row;	
	uint32_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
    seqRel_L = sub->seqRel_L;
    seqRel_R = sub->seqRel_R;
    algnRel_row = sub->algnRel_row;
    L2rel = relat + eB->num_relat; 
    R2rel = L2rel + eB->num_seqL+1; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int i, j, k; 
    int r_L, r_R, r_Rel;
	int p_num;
    
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    int numL = ed_L; 
	int numR = ed_R;
    uint32_t *relat_buf = sub->relat_buf;

    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++) {
            relat_buf[j] = i + 1;
        }
    }

    p_num = sub->pair_out->p_num; 
    for(i=st_R; i<=numR; i++){
        r_R = seqR[i][0];
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++){
            r_Rel = relat[j];
            int flag_Rel = algnRel_row[r_Rel/8] & (1<<(7-r_Rel%8)); 
            if(flag_Rel > 0) continue;	
            k = relat_buf[r_Rel]; 
            if(k == 0){
                continue;
            } else {
                k--;
            }
            pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
            pair_arry[p_num].nxtflg = nxtflg[r_Rel];
            pair_arry[p_num].err = err_val + seqR[i][1] + seqL[k][1];
            pair_arry[p_num].l_off = seqL[k][3] + sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3] + sub->r_off;
            p_num++;                   
            algnRel_row[r_Rel/8] |= (1<<(7-r_Rel%8));
            seqRel[0]++;
            seqRel[seqRel[0]] = r_Rel;
        }  
    }
    sub->pair_out->p_num = p_num; 
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++) {
            relat_buf[j] = 0;
        }
    }
    return p_num;
}

int PairExtSeq_16(struct ExtBlck *eB, struct SubBuf *sub, uint16_t *relat)
{
    int err_val = sub->err_sum[0];
    int max_err = sub->query_err + sub->delta;
    if(err_val >= max_err) return 0;
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
	uint32_t *seqRel, *seqRel_L, *seqRel_R;
    uint8_t *algnRel_row;	
	uint16_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
    seqRel_L = sub->seqRel_L;
    seqRel_R = sub->seqRel_R;
    algnRel_row = sub->algnRel_row;
    L2rel = relat + eB->num_relat; 
    R2rel = L2rel + eB->num_seqL+1; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int i, j, k; 
    int r_L, r_R, r_Rel;
	int p_num;
    
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    int numL = ed_L; 
	int numR = ed_R;
    uint32_t *relat_buf = sub->relat_buf;
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++) {
            //relat_buf[j] = r_L + 1;
            relat_buf[j] = i + 1;
        }
    }
    p_num = sub->pair_out->p_num; 
    for(i=st_R; i<=numR; i++){
        r_R = seqR[i][0];
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++){
            r_Rel = relat[j];
            int flag_Rel = algnRel_row[r_Rel/8] & (1<<(7-r_Rel%8)); 
            if(flag_Rel > 0) continue;	
            k = relat_buf[r_Rel]; 
            if(k == 0){
                continue;
            } else {
                k--;
            }
            pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
            pair_arry[p_num].nxtflg = nxtflg[r_Rel];
            pair_arry[p_num].err = err_val + seqR[i][1] + seqL[k][1];
            pair_arry[p_num].l_off = seqL[k][3] + sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3] + sub->r_off;
            p_num++;                   
            algnRel_row[r_Rel/8] |= (1<<(7-r_Rel%8));
            seqRel[0]++;
            seqRel[seqRel[0]] = r_Rel;
        }  
    }
    sub->pair_out->p_num = p_num; 
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++) {
            relat_buf[j] = 0;
        }
    }
    return p_num;
}

int PairExtSeq_8(struct ExtBlck *eB, struct SubBuf *sub, uint8_t *relat)
{
    int err_val = sub->err_sum[0];
    int max_err = sub->query_err + sub->delta;
    if(err_val >= max_err) return 0;
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
	uint32_t *seqRel, *seqRel_L, *seqRel_R;
    uint8_t *algnRel_row;	
	uint8_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
    seqRel_L = sub->seqRel_L;
    seqRel_R = sub->seqRel_R;
    algnRel_row = sub->algnRel_row;
    L2rel = relat + eB->num_relat; 
    R2rel = L2rel + eB->num_seqL+1; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int i, j, k; 
    int r_L, r_R, r_Rel;
	int p_num;
    
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    int numL = ed_L; 
	int numR = ed_R;
    uint32_t *relat_buf = sub->relat_buf;
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++) {
            //relat_buf[j] = r_L + 1;
            relat_buf[j] = i + 1;
        }
    }
    p_num = sub->pair_out->p_num; 
    for(i=st_R; i<=numR; i++){
        r_R = seqR[i][0];
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++){
            r_Rel = relat[j];
            int flag_Rel = algnRel_row[r_Rel/8] & (1<<(7-r_Rel%8)); 
            if(flag_Rel > 0) continue;	
            k = relat_buf[r_Rel]; 
            if(k == 0){
                continue;
            } else {
                k--;
            }
            pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
            pair_arry[p_num].nxtflg = nxtflg[r_Rel];
            pair_arry[p_num].err = err_val + seqR[i][1] + seqL[k][1];
            pair_arry[p_num].l_off = seqL[k][3] + sub->l_off;
            pair_arry[p_num].r_off = seqR[i][3] + sub->r_off;
            p_num++;                   
            algnRel_row[r_Rel/8] |= (1<<(7-r_Rel%8));
            seqRel[0]++;
            seqRel[seqRel[0]] = r_Rel;
        }  
    }
    sub->pair_out->p_num = p_num; 
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++) {
            relat_buf[j] = 0;
        }
    }
    return p_num;
}
int PairExtSeq_16_uni0( struct ExtBlck *cB, uint16_t *relat, struct StackTree *sTree, int flag, struct SubBuf *sub){
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
    uint8_t *algnRel_buf;	
	uint16_t *L2rel, *R2rel, *R2idx;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + cB->num_relat; 
    R2rel = L2rel + cB->num_seqL+1; 
    R2idx = R2rel + cB->num_seqR+1; 
    nxtpnt= cB->head_nxtpnt + cB->nxtpnt; 
	nxtflg= cB->head_nxtflg + cB->nxtpnt; 
	int i, j; 
    int r_L, r_R, r_Rel;
    /*  
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    */
    int st_L = 1;
    int ed_L = sub->seqL_out[0][0]; 
    int st_R = 1;
    int ed_R = sub->seqR_out[0][0];
 
    int numL = ed_L; 
	int numR = ed_R;
	int p_num = 0;
    if(flag != 0 && flag != 1) {
        printf("%u, %s, rror!!!\n", __LINE__, __func__);
        exit(1); 
    }

    int len_arry = sTree->len_arry;
    uint32_t (*ext_idx)[2];
    uint32_t cls = sTree->stck_arry[len_arry].cls;
    uint32_t cap_row;
    int max_err = sub->query_err + sub->delta;
    int err_val = sub->err_sum[0];
    uint32_t bgn, end, num;
fprintf(stderr, "%u\n", __LINE__);
    if(flag == 0) {
        uint32_t cur_bg, cur_ed;
        for(i=st_L; i<=numL; i++){
            r_L = seqL[i][0];
            bgn = (uint32_t)-1;
            end = 0;
            if(err_val < max_err) {
fprintf(stderr, "%u\n", __LINE__);
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                    r_Rel = relat[j];
                    cap_row = nxtpnt[r_Rel];
                    //ext_idx = (eB+cls)->head_extidx + cap_row;

fprintf(stderr, "%u, i = %d, j = %d, r_Rel = %d, n_relat = %d, p_num = %d, cls = %d, cap_row = %d, n_cap = %d\n", __LINE__, i, j, r_Rel, cB->n_relat, p_num, cls, cap_row, cB->n_cap);
                    ext_idx = (cB+1)->head_extidx + cap_row;
                    //ext_idx = cB->head_extidx + r_Rel;
                    cur_bg = ext_idx[0][0];
                    cur_ed = ext_idx[0][0] + ext_idx[0][1] -1;
                    if(bgn > cur_bg) {
                        bgn = cur_bg;
                    }
                    if(end < cur_ed) {
                        end = cur_ed;
                    }
                } 
            
fprintf(stderr, "%u\n", __LINE__);
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = end+1-bgn;
                pair_arry[p_num].err = err_val + seqL[i][1];
                pair_arry[p_num].l_off = seqL[i][3];
                pair_arry[p_num].r_off = sub->r_off;
                p_num++;                    
            }
        }
    }
fprintf(stderr, "%u\n", __LINE__);
    if(flag == 1) {
       
        cap_row = sTree->stck_arry[len_arry].nxtpnt;
        //ext_idx = (eB+cls)->head_extidx + cap_row;
        ext_idx = cB->head_extidx + cap_row;
        //bgn = ext_idx[0][0];
        //num = ext_idx[0][1];
fprintf(stderr, "%u, st_R = %d, numR = %d, cls = %d, cap_row = %d\n", __LINE__, st_R, numR, cls, cap_row);
        for(i=st_R; i<=numR; i++){
            r_R = seqR[i][0];
            bgn = ext_idx[0][0] + R2idx[r_R];
            num = R2idx[r_R+1] - R2idx[r_R];
            //if(seqR[i][1] > 0) err_val += seqR[i][1];
            if(err_val < max_err) {
                pair_arry[p_num].idx_bg = bgn;
                pair_arry[p_num].idx_num = num; 
                pair_arry[p_num].err = err_val + seqR[i][1];
                pair_arry[p_num].l_off = sub->l_off;
                pair_arry[p_num].r_off = seqR[i][3];
                p_num++;                   
            }
        }
fprintf(stderr, "%u\n", __LINE__);
    }
    sub->pair_out->p_num = p_num; 
    return p_num;
}

int PairExtSeq_16_2(struct ExtBlck *eB, struct SubBuf *sub, uint16_t *relat){
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
	uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
	uint32_t *seqRel, *seqRel_L, *seqRel_R;
    //uint32_t (*pair)[5];
    uint8_t *algnRel_row, *algnRel_buf;	
	uint16_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
    seqRel_L = sub->seqRel_L;
    seqRel_R = sub->seqRel_R;
    algnRel_row = sub->algnRel_row;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + eB->num_relat; 
    R2rel = L2rel + eB->num_seqL+1; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int i, j, k; 
    int r_L, r_R, r_Rel;
	int p_num;
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    int numL = ed_L; 
	int numR = ed_R;
 
    for(i=st_R; i<=numR; i++){
        r_R = seqR[i][0];
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++) {
            algnRel_buf[j/8] |= (1<<(7-j%8));  
        }
    } 
    p_num = sub->pair_out->p_num; 
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
    	    r_Rel = relat[j];
            int flag_Rel = algnRel_row[r_Rel/8] & (1<<(7-r_Rel%8)); 
            if(flag_Rel > 0) continue;	
            int flag_R = algnRel_buf[r_Rel/8] & (1<<(7-r_Rel%8));
            if(flag_R>0){
                for(k = st_R; k <= numR; k++){
                    r_R = seqR[k][0];
                    if(r_Rel >= R2rel[r_R] &&r_Rel < R2rel[r_R+1]) break; 
                }
                if(k > numR){
                    printf("%u, %s algnR_row[] error!, k =%u, numR = %u\n", __LINE__, __func__, k, numR);
                    exit(1);
                }
                /*  
                int penalty = sub->err_sum[0];
                if(seqL[i][1] > 0) penalty += seqL[i][1]*2+1;
                if(seqR[k][1] > 0) penalty += seqR[k][1]*2+1;
                */
                int err_val = sub->err_sum[0];
                if(seqL[i][1] > 0) err_val += seqL[i][1];
                if(seqR[k][1] > 0) err_val += seqR[k][1];
                int max_err = sub->query_err + sub->delta;
                if(err_val < max_err) {
                    pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
                    pair_arry[p_num].nxtflg = nxtflg[r_Rel];
                    pair_arry[p_num].err = err_val;
                    //pair_arry[p_num].l_off = seqL[i][3];
                    //pair_arry[p_num].r_off = seqR[k][3];
                
                    pair_arry[p_num].l_off = seqL[i][3] + sub->l_off;
                    pair_arry[p_num].r_off = seqR[k][3] + sub->r_off;


                    
                    p_num++;                   
                    algnRel_row[r_Rel/8] |= (1<<(7-r_Rel%8));
                    seqRel[0]++;
                    seqRel[seqRel[0]] = r_Rel;
                }
            }
		}  
    }
    sub->pair_out->p_num = p_num;
    for(i=st_R; i<=numR; i++){
        r_R = seqR[i][0];
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++) {
            algnRel_buf[j/8] = 0;
        }
    }
    return p_num;
}

int PairExtSeq_32_2(struct ExtBlck *eB, struct SubBuf *sub, uint32_t *relat){

    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
	uint32_t (*seqL)[4];
	uint32_t (*seqR)[4];
	uint32_t *seqRel, *seqRel_L, *seqRel_R;
	uint8_t *algnRel_row, *algnRel_buf;	
	uint32_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
    seqRel_L = sub->seqRel_L;
    seqRel_R = sub->seqRel_R;
    algnRel_row = sub->algnRel_row;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + eB->num_relat; 
    R2rel = relat+eB->num_relat + eB->num_seqL+1;

	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int i, j, k; 
    int r_L, r_R, r_Rel;
	int p_num;
    int p_id = sub->pair_b2e_id;
    int st_L = sub->pair_b2e[p_id][0];
    int ed_L = sub->pair_b2e[p_id][1]; 
    int st_R = sub->pair_b2e[p_id][2];
    int ed_R = sub->pair_b2e[p_id][3];
    int numL = ed_L; 
	int numR = ed_R;
    for(i=st_R; i<=numR; i++){
        r_R = seqR[i][0];

        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++) {
            algnRel_buf[j/8] |= (1<<(7-j%8));  
        }
    } 
    p_num = sub->pair_out->p_num; 
    for(i=st_L; i<=numL; i++){
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
    	    r_Rel = relat[j];
            int flag_Rel = algnRel_row[r_Rel/8] & (1<<(7-r_Rel%8)); 
            if(flag_Rel > 0) continue;	
            int flag_R = algnRel_buf[r_Rel/8] & (1<<(7-r_Rel%8));
            if(flag_R>0){
                for(k = st_R; k <= numR; k++){
                    r_R = seqR[k][0];
                    if(r_Rel >= R2rel[r_R] &&r_Rel < R2rel[r_R+1]) break; 
                }
                if(k > numR){
                    printf("%u, %s algnR_row[] error!, k =%u, numR = %u\n", __LINE__, __func__, k, numR);
                    exit(1);
                    //continue;
                }
                /*  
                int penalty = sub->err_sum[0];
                if(seqL[i][1] > 0) penalty += seqL[i][1]*2+1;
                if(seqR[k][1] > 0) penalty += seqR[k][1]*2+1;
                */
                int err_val = sub->err_sum[0];
                if(seqL[i][1] > 0) err_val += seqL[i][1];
                if(seqR[k][1] > 0) err_val += seqR[k][1];
                int max_err = sub->query_err + sub->delta;
                if(err_val < max_err) {
                    pair_arry[p_num].nxtpnt = nxtpnt[r_Rel];
                    pair_arry[p_num].nxtflg = nxtflg[r_Rel];
                    pair_arry[p_num].err = err_val;
                    //pair_arry[p_num].l_off = seqL[i][3];
                    //pair_arry[p_num].r_off = seqR[k][3];
                    pair_arry[p_num].l_off = seqL[i][3] + sub->l_off;
                    pair_arry[p_num].r_off = seqR[k][3] + sub->r_off;
                    
                    p_num++;
                    algnRel_row[r_Rel/8] |= (1<<(7-r_Rel%8));
                    seqRel[0]++;
                    seqRel[seqRel[0]] = r_Rel;
                }
            }
		}  
    }
    sub->pair_out->p_num = p_num;
    for(i=st_R; i<=numR; i++){
        r_R = seqR[i][0];
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++) {
            algnRel_buf[j/8] = 0;
        }
    }
    return p_num;
}
int PairExtSeq_all(struct ExtBlck *eB, struct SubBuf *sub)
{
    uint32_t *relat = eB->head_relat + eB->relat; 
    int flag = 0;
    int i;
    if(eB->num_relat<=0xFF) {
        uint8_t *relat_8 = (uint8_t *)relat;
        flag = PairExtSeq_8(eB, sub, relat_8);
    } else if(eB->num_relat<=0xFFFF){
        uint16_t *relat_16 = (uint16_t *)relat;
        flag = PairExtSeq_16(eB, sub, relat_16);
    } else{
        flag = PairExtSeq_32(eB, sub, relat);
    }    
    int p_num = sub->pair_out->p_num;
    return p_num;
}
int PairExtSeq_all_uni(idx_t *fm_idx, struct ExtBlck *cB, struct StackTree *sTree, int flg,  struct SubBuf *sub)
{
    uint32_t *relat = cB->head_relat + cB->relat; 
    int ret_flg = 0;
    int i;

    if(cB->num_relat<=0xFF) {
        uint8_t *relat_8 = (uint8_t *)relat;
fprintf(stderr, "%u\n", __LINE__);
        ret_flg = PairExtSeq_8_uni(fm_idx, cB, relat_8, sTree, flg, sub);
    } else if(cB->num_relat<=0xFFFF){
        uint16_t *relat_16 = (uint16_t *)relat;
fprintf(stderr, "%u\n", __LINE__);
        ret_flg = PairExtSeq_16_uni(fm_idx, cB, relat_16, sTree, flg, sub);
    } else{
fprintf(stderr, "%u\n", __LINE__);
        ret_flg = PairExtSeq_32_uni(fm_idx, cB, relat, sTree, flg, sub);
    }    
fprintf(stderr, "%u\n", __LINE__);
    int p_num = sub->pair_out->p_num;
    return p_num;
}
int aln_smbwt_uni_L(idx_t *fm_idx, struct JmpMod *jmp, struct ExtBlck *eBlck, struct ExtBlck *cB, query_t *query, seed_t *seed, struct StackTree *sTree,  int gen_flg, struct SubBuf *sub)
{
    int p_num = 0;
    //设定read的左侧16个长度的子序列
    if(gen_flg > 0) { 
fprintf(stderr, "%u, seed_id = %d\n", __LINE__, seed->id);
        sTree->len_arry = 0;
        setBlckData_uni(jmp, eBlck, &cB, query, seed, sTree, 0, sub);
        //调用
        int max_err = get_max_err(sub);
        int sub_err = sub->sub_err;
        int ret_flg = Algn_sub_once(cB,sub,0);
fprintf(stderr, "%u, ret_flg = %d\n", __LINE__, ret_flg);
        if(ret_flg > 0) {
            if(max_err > sub_err) {
                ret_flg = Algn_sub_all(cB, sub, 0);
            }         
        } else if(ret_flg == 0) {
            return p_num; 
        }
    } else {
        if(sub->seqL_out[0][0] <= 0) {
            printf("[error]: %u, num_seqL <= 0!\n", __LINE__);
            exit(1); 
        } 
    }
    if(sub->seqL_out[0][0] > 0) {
        p_num = PairExtSeq_all_uni(fm_idx, cB, sTree, 0, sub);
    }
    //sub->seqL_out = sub->call_buf + sub->seqL_off[NUM_EXT]; 
    //sub->seqR_out = sub->call_buf + sub->seqR_off[NUM_EXT]; 
    /*  
    int i, r;
    for(i=1; i<= sub->seqL_out[0][0]; i++){
        r = sub->seqL_out[i][0];
        sub->algnL_row[r/8] = 0; 
    }
    for(i=1; i<=sub->seqR_out[0][0]; i++){
        r = sub->seqR_out[i][0];
        sub->algnR_row[r/8] = 0; 
    }
    */
    
    return p_num;
}
int aln_smbwt_uni_R(idx_t *fm_idx, struct JmpMod *jmp, struct ExtBlck *eBlck, query_t *query, seed_t *seed, struct StackTree *sTree,  int gen_flg, struct SubBuf *sub)
{
    int p_num = 0;
    struct ExtBlck *cB; 
    //设定read的左侧16个长度的子序列
    if(gen_flg > 0) { 
fprintf(stderr, "%u, seed_id = %d\n", __LINE__, seed->id);
        sTree->len_arry = 0;
        setBlckData_uni(jmp, eBlck, &cB, query, seed, sTree, 1, sub);
fprintf(stderr, "%u\n", __LINE__);
        //调用
        int max_err = get_max_err(sub);
        int sub_err = sub->sub_err;
        int ret_flg = Algn_sub_once(cB,sub,1);
        if(ret_flg > 0) {
            //if(ret_flg > 1) sub->seqR_out[0][0] = 0;
            if(max_err > sub_err) {
                ret_flg = Algn_sub_all(cB, sub, 1);
            }         
        } else {
            return p_num; 
        }
    } else {
        if(sub->seqR_out[0][0] <= 0) {
            printf("[error]: %u, num_seqR <= 0!\n", __LINE__);
            exit(1); 
        } 
    }
fprintf(stderr, "%u \n", __LINE__);
    if(sub->seqR_out[0][0] > 0) {
fprintf(stderr, "%u sub->seq_R[0][0] = %d\n", __LINE__, sub->seqR_out[0][0]);
        p_num = PairExtSeq_all_uni(fm_idx, cB, sTree, 1,sub);
fprintf(stderr, "%u\n", __LINE__);
    }
    /*  
    sub->seqL_out = sub->call_buf + sub->seqL_off[NUM_EXT]; 
    sub->seqR_out = sub->call_buf + sub->seqR_off[NUM_EXT]; 
    int i, r;
    for(i=1; i<= sub->seqL_out[0][0]; i++){
        r = sub->seqL_out[i][0];
        sub->algnL_row[r/8] = 0; 
    }
    for(i=1; i<=sub->seqR_out[0][0]; i++){
        r = sub->seqR_out[i][0];
        sub->algnR_row[r/8] = 0; 
    }
    */
    return p_num;
}



int sortPairExtSeq(struct SubBuf *sub)
{
    int i, j, k;
    pair_arry_t *pair_arry = sub->pair_out->pair_arry;
    int p_num = sub->pair_out->p_num;
    int s_num = 0;
    int b_num = p_num-1;
    int NEXT_EXT_NUM = sub->NEXT_EXT_NUM;    
    if(p_num == 1){
        //if(pair_arry[0].nxtflg <= IS_SMLSIZ){ 
        if(pair_arry[0].nxtflg <= NEXT_EXT_NUM){ 
            sub->pair_out->s_num = 1;
        } else{
            sub->pair_out->s_num = 0; 
        }
        return 1;  
    }
    //+++++++++++++++++++++++++ 
    j = 0;
    k = p_num;
    for(i = 0; i< p_num; ++i){
        if( pair_arry[i].nxtflg > NEXT_EXT_NUM) {
            memcpy(pair_arry+k,   pair_arry+i,   sizeof(pair_arry_t));
            k++;
        } else{
            memcpy(pair_arry+j,   pair_arry+i,   sizeof(pair_arry_t));
            j++;
        }
    }

    memcpy(pair_arry+j,   pair_arry+p_num,   (k-p_num)*sizeof(pair_arry_t));
    s_num = j;
/*  
    while(p_num>1){
		//while(pair_arry[s_num].nxtflg<=IS_SMLSIZ){
		while(pair_arry[s_num].nxtflg<=NEXT_EXT_NUM){
			s_num++;
			if(s_num>b_num)break;
		}
		if(s_num>b_num)break;
		while(pair_arry[b_num].nxtflg>NEXT_EXT_NUM){
			b_num--;
			if(s_num>b_num)break;
		}
		if(s_num>b_num)break;
		if(s_num==b_num){
            if(pair_arry[s_num].nxtflg<=NEXT_EXT_NUM)s_num++; 
            else b_num--;
            break;
		}
        memcpy(pair_arry+p_num,   pair_arry+s_num,   sizeof(pair_arry_t));
		memcpy(pair_arry+s_num,   pair_arry+b_num,   sizeof(pair_arry_t));
		memcpy(pair_arry+b_num,   pair_arry+p_num,   sizeof(pair_arry_t));
        s_num++; 
		b_num--;
		if(s_num>b_num){
            break;
		}
		if(s_num == b_num){
          if(pair_arry[s_num].nxtflg<=NEXT_EXT_NUM)s_num++; 
          else b_num--;
          break;
		}	
	}
*/
    //+++++++++++++++++++++++++++++++++++
    sub->pair_out->s_num = s_num;
    int bg = s_num, ed = p_num-1;
    int max, row;
    while( bg < ed){
        max = pair_arry[bg].err*256 + 256 - pair_arry[bg].nxtflg;
        row = bg;
        for(i = bg+1; i <= ed; ++i){
            if(pair_arry[i].err*256 + 256 - pair_arry[i].nxtflg > max) {
                max = pair_arry[i].err*256 + 256 - pair_arry[i].nxtflg;
                row = i; 
            } 
        }
        if(row > bg){
            memcpy(pair_arry+p_num,   pair_arry+bg,      sizeof(pair_arry_t));
    		memcpy(pair_arry+bg,      pair_arry+row,     sizeof(pair_arry_t));
		    memcpy(pair_arry+row,     pair_arry+p_num,   sizeof(pair_arry_t));   
        }
        bg++; 
    }
/*  
    //fprintf(stderr, "%u, p_num = %u, s_num = %u\n", 
            __LINE__, p_num, s_num);
    for(i = 0; i< p_num; ++i){
         //fprintf(stderr, "err = %u, nxtflg = %u\n", 
                 pair_arry[i].err, pair_arry[i].nxtflg);
    }

    int err_flg = 0;
    for(i = 0; i< s_num; ++i){
         if( pair_arry[i].nxtflg > NEXT_EXT_NUM) {
            printf("%u, error nxtflg > IS_SMLSIZ!\n", __LINE__); 
            //fprintf(stderr, "nxtflg = %u\n", pair_arry[i].nxtflg);
            //exit(1);
            //err_flg = -1;
         }
    } 
    for(i = s_num; i< p_num; ++i){
        if( pair_arry[i].nxtflg <= NEXT_EXT_NUM) {
            printf("%u, error nxtflg > IS_SMLSIZ!\n", __LINE__); 
            //fprintf(stderr, "nxtflg = %u\n", pair_arry[i].nxtflg);
            exit(1);
            err_flg = -1;
        }
        int penalty0 = pair_arry[i].err*256+256-pair_arry[i].nxtflg;
        int penalty1 = pair_arry[i+1].err*256+256-pair_arry[i+1].nxtflg;
        if(i + 1< p_num) {
            if(penalty0 < penalty1) {
                printf("%u, %s, error sort\n", __LINE__, __func__); 
                //fprintf(stderr, "%u, %s, error sort\n", __LINE__, __func__);
                //exit(1);
                //err_flg = -1;
            } 
        }
    }

    if(err_flg < 0) exit(1);
*/
    return p_num;
}


/*
  
for(i=0; i<=p_num; ++i){
//fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
}

    //fprintf(stderr, "pair_num = %u\n", p_num); 
    for(i=0; i<=p_num; ++i){
        //fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
    
    }
*/

//++++++++++++++++++++++++++++++++++++
    //int s_num = 1;
/*
    //fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
    for(i=0; i<=p_num; ++i){
        //fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
    
    }

*/
    //+++++++++++++++++++++++++++++++++++++++++++
    //sort(pair, s_num, p_num, key=pair[i][2])
    //从1开始s_num升序排序， 以后需要测试
//ks_introsort(turple_t, p_num-s_num+1, pair[0]+s_num); 
   
//以下代码处理左扩展序列查找与配对++++++++++++++++++++++++++++++++++
/*  
            for(i = 0; i < numR; ++i){
     		    r_R = seqR[i][0];
                for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
    		        r_R = relat[j];
                    int err_num = 0;
                    int k;
                    for(k = 7; k>=0; --k){
                        seq_ch = seq[2*k]<<2|seq[2*k+1];      
                        //get_ch_bwt(r_R, cB, sub, 1, &ch);
                        if(seq_ch != ch){
                            err_num++; 
                        } 
                        if(err_num > 1){
                            break; 
                        }
                    }
                    if(k <0){
                        p_num++;
                        pair[p_num][0] = nxtpnt[j];
                        pair[p_num][1] = nxtflg[j];
                        pair[p_num][2] = seqL[i][1] + err_num;
                    } 
                }
            }//end  for(i = 0; i < numL; ++i)++++++++++++++++++++++++++++++
*/       






    //清理缓冲区++++++++++++++++++++++++++++++++++++++
/*  
    numR = seqR[0][0];
	for(i=1; i<=numR; i++){
		r_R = seqR[i][0];
		algnR_row[r_R] = 0;
	}
*/
/*
//fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
for(i=0; i<=p_num; ++i){
    //fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 

}
*/

/*
for(i = 0; i < eB->num_seqL; i++){
    //fprintf(stderr, "L2rel[%u] = %u\n", i, L2rel[i]);
}
for(i = 0; i < eB->num_relat; ++i){
    //if(nxtpnt[i]>17484874&& nxtpnt[i]<17484874+200)
    //fprintf(stderr, "relat[%u] = %u, nxtpnt = %u, nxtflg = %u\n", i, relat[i], nxtpnt[i], nxtflg[i]);
}
*/
////fprintf(stderr, "%s, 1130, seqL, seqR = %u, %u\n", __func__, seqL[1][0], seqR[1][0]);
/*
	for(i=0; i<=eB->num_seqR; i++){
      


        if(algnR_row[i] >0){
            //fprintf(stderr, " algnR_row[%u] = %u\n", i, algnR_row[i]); 
            exit(1);
        }	
	
	}
*/


int PairExtSeq_R(struct ExtBlck *eB, struct SubBuf *sub){

	uint32_t (*seqL)[3];
	uint32_t (*seqR)[3];
	uint32_t (*pair)[3];
	uint8_t *algnR_row;	
	uint32_t *relat;
	uint32_t *L2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	pair = sub->pair_out;
    algnR_row = sub->algnR_row; 
	relat = eB->head_relat + eB->relat; 
	L2rel = relat + eB->num_relat; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int numL,numR,i,j,k,r_L,r_R,p_num=0;
	numL = seqL[0][0]; 
	numR = seqR[0][0];
//fprintf(stderr, "seqL[0][1] = %u, seqR[0][1] = %u\n", seqL[0][1], seqR[0][1]);
//fprintf(stderr, "%s, 1129, numL, numR = %u, %u\n", __func__, numL, numR);
	
   
    numR = seqR[0][0];
    numL = seqL[0][0];
    uint8_t *seq;    
    uint8_t seq_ch, ch;
    int row = 0;
//fprintf(stderr, "%u, %s p_num=%u, numR = %u, numL= %u\n", __LINE__, __func__, p_num, numR, numL);
//fprintf(stderr, "%u, %s p_num=%u, numR = %u, numL= %u\n", __LINE__, __func__, p_num, numR, numL);

    if(numL >0 && numR > 0){
        //以下代码处理右扩展序列查找与配对++++++++++++++++++++++++++++++++++
     
        for(i=1; i <= numR; ++i){
            r_R = seqR[i][0];
            algnR_row[r_R/8] = algnR_row[r_R/8] | (1<<(7-r_R%8)) ;
        } 
    
        uint8_t *smbwt = eB->head_smbwt + eB->bwtR;
        uint32_t n = eB->num_seqR; 
        uint8_t *cnt = smbwt+((n+1)/2)*8;
        uint8_t cnt1[8][17];
        if(n <= NO_BWT_SUM_SIZE && n > MIN_BWT_SIZE) gen_cnt(smbwt, n, cnt1);
        seq = sub->ext_seqR;    

//fprintf(stderr, "%u, %s n_data=%u\n", __LINE__, __func__, n);          
        for(i = 1; i <= numL; ++i){
            r_L = seqL[i][0];
            for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                r_R = relat[j];
                int flag_R = algnR_row[r_R/8] & (1<<(7-r_R%8));
                if( flag_R > 0) continue;
                int err_num = 0;
                //int k;

//fprintf(stderr, "%u, %s r_R=%u\n", __LINE__, __func__, r_R);          

                    for(k = 7; k>=0; --k){
                        seq_ch = seq[2*k]<<2|seq[2*k+1];      
                        if(n<=MIN_BWT_SIZE){
                            //r_R = get_ch_min(r_R, k, smbwt, n, &ch); 
                            break;
                        } else if(n<=NO_BWT_SUM_SIZE){ 
                            r_R = get_ch_nosum(r_R, k, smbwt, cnt1, n, &ch); 
                        } else if(n <= 255) {
                            r_R = get_ch_255(r_R, k, smbwt, (uint8_t (*)[17])cnt, n, &ch);     
    //fprintf(stderr, "%u, %s r_R=%u, k = %u, ch=%u\n", __LINE__, __func__, r_R, k, ch);          
                        } else{//>=256
                            r_R = get_ch_large_alt(r_R, k, smbwt, cnt, n, &ch);
                        }
                        if(seq_ch != ch){
                            err_num++; 
                        } 
                        if(err_num > 1){
                            break; 
                        }
                    }
                    if(k <0){
                        p_num++;
                        pair[p_num][0] = nxtpnt[j];
                        pair[p_num][1] = nxtflg[j];
                        pair[p_num][2] = seqL[i][1] + err_num;
                    }
/*  
                        seqR[0][0]++;                    
                        algnR_row[r_R/8] = algnR_row[r_R/8] | (1<<(7-r_R%8));
                        seqR[seqR[0][0]][0] = r_R;
                        seqR[seqR[0][0]][1] = err_num;
                       

                    }
                } else{//if(algnR_row[r_R] == 0) ++++++++++++++++++++++++
                    p_num++;
                    pair[p_num][0] = nxtpnt[j];
                    pair[p_num][1] = nxtflg[j];
                    pair[p_num][2] = seqL[i][1] + err_num;
                } 
*/
            }//end for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++)+++++++++++++++++++
        }//end  for(i = 0; i < numL; ++i)+++++++++++++++++++++++++++++++++++

//fprintf(stderr, "%u, %s p_num=%u\n", __LINE__, __func__, p_num);          
        //以下代码处理左扩展序列查找与配对++++++++++++++++++++++++++++++++++
/*  
        for(i = 0; i < numR; ++i){
            r_R = seqR[i][0];
            for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                r_R = relat[j];
                int err_num = 0;
                int k;
                for(k = 7; k>=0; --k){
                    seq_ch = seq[2*k]<<2|seq[2*k+1];      
                    //get_ch_bwt(r_R, cB, sub, 1, &ch);
                    if(seq_ch != ch){
                        err_num++; 
                    } 
                    if(err_num > 1){
                        break; 
                    }
                }
                if(k <0){
                    p_num++;
                    pair[p_num][0] = nxtpnt[j];
                    pair[p_num][1] = nxtflg[j];
                    pair[p_num][2] = seqL[i][1] + err_num;
                } 
            }
        }//end  for(i = 0; i < numL; ++i)++++++++++++++++++++++++++++++
*/       
    }
    int s_num = 1;
    int b_num = p_num;
    if(p_num == 1){
        if(pair[1][1]<=IS_SMLSIZ){ 
            s_num = 2;
        } 
    }

	while(p_num>1){
		while(pair[s_num][1]<=IS_SMLSIZ){
			s_num++;

			if(s_num>b_num)break;
		}
		if(s_num>b_num)break;
		while(pair[b_num][1]>IS_SMLSIZ){
			b_num--;
			if(s_num>b_num)break;
		}
		if(s_num>b_num)break;
		if(s_num==b_num){
          if(pair[s_num][1]<=IS_SMLSIZ)s_num++; 
          else b_num--;
          break;
		}

		memcpy(pair+p_num+1, pair+s_num,3*sizeof(uint32_t));
		memcpy(pair+s_num, pair+b_num,3*sizeof(uint32_t));
		memcpy(pair+b_num, pair+p_num+1,3*sizeof(uint32_t));
		s_num++; 
		b_num--;
		if(s_num>b_num){
          break;
		}
		if(s_num == b_num){
          if(pair[s_num][1]<=IS_SMLSIZ)s_num++; 
          else b_num--;
          break;
		}	
	}

	pair[0][0] = p_num;
	pair[0][1] = s_num-1;
	pair[0][2] = 0;

   

    //清理缓冲区++++++++++++++++++++++++++++++++++++++
/*  
    numR = seqR[0][0];
	for(i=1; i<=numR; i++){
		r_R = seqR[i][0];
		algnR_row[r_R] = 0;
	}
*/
/*
//fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
for(i=0; i<=p_num; ++i){
    //fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 

}
*/

	return (int) p_num;
}

int PairExtSeq_L(struct ExtBlck *eB, struct SubBuf *sub){

	uint32_t (*seqL)[3];
	uint32_t (*seqR)[3];
	uint32_t (*pair)[3];

	uint8_t *algnR_row;	
	uint32_t *relat;
	uint32_t *L2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	pair = sub->pair_out;
    algnR_row = sub->algnR_row; 
	relat = eB->head_relat + eB->relat; 
	L2rel = relat + eB->num_relat; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int numL,numR,i,j,k,r_L,r_R,p_num=0;
	numL = seqL[0][0]; 
	numR = seqR[0][0];
//fprintf(stderr, "seqL[0][1] = %u, seqR[0][1] = %u\n", seqL[0][1], seqR[0][1]);
//fprintf(stderr, "%s, 1129, numL, numR = %u, %u\n", __func__, numL, numR);
	
    uint8_t *seq;    
    uint8_t seq_ch, ch;
    int row = 0;
//fprintf(stderr, "%u, %s p_num=%u, numR = %u, numL= %u\n", __LINE__, __func__, p_num, numR, numL);
//fprintf(stderr, "%u, %s p_num=%u, numR = %u, numL= %u\n", __LINE__, __func__, p_num, numR, numL);

    if(numR >0 && numL >0){
        //以下代码处理右扩展序列查找与配对++++++++++++++++++++++++++++++++++
/*  
        for(i=1; i <= numR; ++i){
            r_R = seqR[i][0];
            algnR_row[r_R/8] = algnR_row[r_R/8] | (1<<(7-r_R%8));
        } 
*/
        uint8_t *smbwt = eB->head_smbwt + eB->bwtL;
        uint32_t n_L = eB->num_seqL; 
        uint32_t n_relat = eB->num_relat; 
        uint8_t *cnt = smbwt+((n_L+1)/2)*8;
        uint8_t cnt1[8][17];
        if(n_L <= NO_BWT_SUM_SIZE && n_L > MIN_BWT_SIZE) gen_cnt(smbwt, n_L, cnt1);
        seq = sub->ext_seqL;    
//fprintf(stderr, "%u, num_seqL = %u, num_relat = %u, num_seqR =%u\n", __LINE__, n_L, n_relat, eB->num_seqR);        
        //if(n_relat < 256) 
        {
            for(i = 0; i < n_L; ++i){
                int find_flag = 0;
                int r_L = i;
                for(j =  L2rel[i]; j < L2rel[i+1]; ++j){
                    r_R = relat[j];
                    int flag_R = algnR_row[r_R/8] & (1<<(7-r_R%8));
                    if(flag_R > 0){
                        int err_num = 0;
                        //int k;

//fprintf(stderr, "%u\n", __LINE__);        
                        if(find_flag == 0){

//fprintf(stderr, "%u\n", __LINE__);        
                            for(k = 7; k>=0; --k){
                                seq_ch = seq[2*k]<<2|seq[2*k+1];      
                                if(n_L<=MIN_BWT_SIZE){
                                    //r_R = get_ch_min(r_R, k, smbwt, n, &ch); 
                                    break;
                                } else if(n_L<=NO_BWT_SUM_SIZE){ 
                                    r_L = get_ch_nosum(r_L, k, smbwt, cnt1, n_L, &ch); 
                                } else if(n_L <= 255) {
                                    r_L = get_ch_255(r_L, k, smbwt, (uint8_t (*)[17])cnt, n_L, &ch);     
            //fprintf(stderr, "%u, %s r_R=%u, k = %u, ch=%u\n", __LINE__, __func__, r_R, k, ch);          
                                } else{//>=256
                                    r_L = get_ch_large_alt(r_L, k, smbwt, cnt, n_L, &ch);
                                }
                                if(seq_ch != ch){
                                    err_num++; 
                                } 
                                if(err_num > 1){
                                    break; 
                                }
                            }

//fprintf(stderr, "%u\n", __LINE__);        
                            if(k <0){

//fprintf(stderr, "%u\n", __LINE__);        
//fprintf(stderr, "p_num =%u. i =%u, j = %u\n", p_num, i, j);
                                for(k = 1; k <= numR; k++){
                                    if(seqR[k][0] == r_R) break;
                                }
                                if(k > numR){
                                    printf("%u, %s algnR_row[] error!, k =%u, numR = %u\n", __LINE__, __func__, k, numR);
                                    exit(1);
                                }

                                p_num++;
                                pair[p_num][0] = nxtpnt[j];
                                pair[p_num][1] = nxtflg[j];
                                pair[p_num][2] = seqR[k][1] + err_num;
                                find_flag = 1; 

//fprintf(stderr, "%u\n", __LINE__);        
                            } else{

//fprintf(stderr, "%u\n", __LINE__);        
                                break;//for(j =  L2rel[i]; j < L2rel[i+1]; ++j)++++
                            }

//fprintf(stderr, "%u\n", __LINE__);        
                        } else{
//fprintf(stderr, "%u\n", __LINE__);        
                            p_num++;
                            pair[p_num][0] = nxtpnt[j];
                            pair[p_num][1] = nxtflg[j];
                          
                            for(k = 1; k <= numR; k++){
                                if(seqR[k][0] == r_R) break;
                            }
                            if(k > numR){
                                printf("%u, %s algnR_row[] error!, k =%u, numR = %u\n", __LINE__, __func__, k, numR);
                                exit(1);
                            }
                            pair[p_num][2] =err_num+seqR[k][0];
                            
                            find_flag = 1; 
                        }

//fprintf(stderr, "%u\n", __LINE__);        
                    }// end if(algnR_row[r_R] > 0)+++++++++++++++++++

                } // end for(j =  L2rel[i]; j < L2rel[i+1]; ++j)++++++++
            }// end for(i = 0; i < n_L; ++i)++++++++++++++++++++++++++++
        } // end if(n_relat < 256) +++++++++++++++++++++++++++++++++++++
    } // end if(numR >0 && numL >0)+++++++++++++++++++++++++++++++++++++
    
    int s_num = 1;
    int b_num = p_num;
    if(p_num == 1){
        if(pair[1][1]<=IS_SMLSIZ){ 
            s_num = 2;
        } 
    }

	while(p_num>1){
		while(pair[s_num][1]<=IS_SMLSIZ){
			s_num++;
            if(s_num>b_num)break;
		}
		if(s_num>b_num)break;
		while(pair[b_num][1]>IS_SMLSIZ){
			b_num--;
			if(s_num>b_num)break;
		}
		if(s_num>b_num)break;
		if(s_num==b_num){
          if(pair[s_num][1]<=IS_SMLSIZ)s_num++; 
          else b_num--;
          break;
		}

		memcpy(pair+p_num+1, pair+s_num,3*sizeof(uint32_t));
		memcpy(pair+s_num, pair+b_num,3*sizeof(uint32_t));
		memcpy(pair+b_num, pair+p_num+1,3*sizeof(uint32_t));
		s_num++; 
		b_num--;
		if(s_num>b_num){
          break;
		}
		if(s_num == b_num){
          if(pair[s_num][1]<=IS_SMLSIZ)s_num++; 
          else b_num--;
          break;
		}	
	}

	pair[0][0] = p_num;
	pair[0][1] = s_num-1;
	pair[0][2] = 0;

   

    //清理缓冲区++++++++++++++++++++++++++++++++++++++
/*  
    numR = seqR[0][0];
	for(i=1; i<=numR; i++){
		r_R = seqR[i][0];
		algnR_row[r_R] = 0;
	}
*/
/*
//fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
for(i=0; i<=p_num; ++i){
    //fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 

}
*/

	return (int) p_num;
}

int PairExtSeq(struct ExtBlck *eB, struct SubBuf *sub){

	uint32_t (*seqL)[3];
	uint32_t (*seqR)[3];
	uint32_t (*pair)[3];

	uint8_t *algnR_row;	
	uint32_t *relat;
	uint32_t *L2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
	pair = sub->pair_out;
    algnR_row = sub->algnR_row; 
	relat = eB->head_relat + eB->relat; 
	L2rel = relat + eB->num_relat; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	uint32_t numL,numR,i,j,k, r_L,r_R,p_num;
	numL = seqL[0][0]; 
	numR = seqR[0][0];
//fprintf(stderr, "seqL[0][1] = %u, seqR[0][1] = %u\n", seqL[0][1], seqR[0][1]);
//fprintf(stderr, "%s, 1129, numL, numR = %u, %u\n", __func__, numL, numR);
/*
for(i = 0; i < eB->num_seqL; i++){
    //fprintf(stderr, "L2rel[%u] = %u\n", i, L2rel[i]);
}
for(i = 0; i < eB->num_relat; ++i){
    //if(nxtpnt[i]>17484874&& nxtpnt[i]<17484874+200)
    //fprintf(stderr, "relat[%u] = %u, nxtpnt = %u, nxtflg = %u\n", i, relat[i], nxtpnt[i], nxtflg[i]);
}
*/
////fprintf(stderr, "%s, 1130, seqL, seqR = %u, %u\n", __func__, seqL[1][0], seqR[1][0]);
/*
	for(i=0; i<=eB->num_seqR; i++){
        //algnR_row[i] = 0;


        if(algnR_row[i] >0){
            //fprintf(stderr, " algnR_row[%u] = %u\n", i, algnR_row[i]); 
            exit(1);
        }	
	
	}
*/
	for(i=1; i<=numR; i++){
		r_R = seqR[i][0];
        algnR_row[r_R/8] = algnR_row[r_R/8] | (1<<(7-r_R%8));
//fprintf(stderr, "i, r_R = %u %u\n", i, r_R);
    }
    p_num = 0;
	for(i=1; i<=numL; i++){
		r_L = seqL[i][0];

//fprintf(stderr, "i, r_L = %u %u\n", i, r_L);
//fprintf(stderr, "i = %u, L2rel[%u] = %u, %u\n", i, r_L, L2rel[r_L],L2rel[r_L+1]);

        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
    		r_R = relat[j];

//fprintf(stderr, "r_L = %u, r_R= %u\n", r_L, r_R);
////fprintf(stderr, "algnR_row[%u]= %u\n", r_R, algnR_row[r_R]);
    		int flag_R = algnR_row[r_R/8] & (1<<(7-r_R%8));
            if(flag_R>0){
				for(k = 1; k <= numR; k++){
                    if(seqR[k][0] == r_R) break;
                }
                if(k > numR){
                    printf("%u, %s algnR_row[] error!\n", __LINE__, __func__);
                    exit(1);
                }
                p_num++;
                pair[p_num][0] = nxtpnt[j];
				pair[p_num][1] = nxtflg[j];
				pair[p_num][2] = seqL[i][1] + seqR[k][1];

////fprintf(stderr, "pair = %u, %u, %u\n", pair[p_num][0], pair[p_num][1], pair[p_num][2]); 
    		}
		}  
    	
    }


/*
  
for(i=0; i<=p_num; ++i){
//fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
}

    //fprintf(stderr, "pair_num = %u\n", p_num); 
    for(i=0; i<=p_num; ++i){
        //fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
    
    }
*/

//++++++++++++++++++++++++++++++++++++
    //int s_num = 1;
    int s_num = 1;
    int b_num = p_num;
    if(p_num == 1){
        if(pair[1][1]<=IS_SMLSIZ){ 
            s_num = 2;
        } 
    }

	while(p_num>1){
		while(pair[s_num][1]<=IS_SMLSIZ){
			s_num++;

			if(s_num>b_num)break;
		}
		if(s_num>b_num)break;
		while(pair[b_num][1]>IS_SMLSIZ){
			b_num--;
			if(s_num>b_num)break;
		}
		if(s_num>b_num)break;
		if(s_num==b_num){
          if(pair[s_num][1]<=IS_SMLSIZ)s_num++; 
          else b_num--;
          break;
		}

		memcpy(pair+p_num+1, pair+s_num,3*sizeof(uint32_t));
		memcpy(pair+s_num, pair+b_num,3*sizeof(uint32_t));
		memcpy(pair+b_num, pair+p_num+1,3*sizeof(uint32_t));
		s_num++; 
		b_num--;
		if(s_num>b_num){
          break;
		}
		if(s_num == b_num){
          if(pair[s_num][1]<=IS_SMLSIZ)s_num++; 
          else b_num--;
          break;
		}	
	}

	pair[0][0] = p_num;
	pair[0][1] = s_num-1;
	pair[0][2] = 0;
/*
    //fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
    for(i=0; i<=p_num; ++i){
        //fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
    
    }

*/
    //+++++++++++++++++++++++++++++++++++++++++++
    //sort(pair, s_num, p_num, key=pair[i][2])
    //从1开始s_num升序排序， 以后需要测试
//ks_introsort(turple_t, p_num-s_num+1, pair[0]+s_num); 
    //algnR_row[0] = p_num;
    //清理缓冲区++++++++++++++++++++++++++++++++++++++
    numR = seqR[0][0];
	for(i=1; i<=numR; i++){
		r_R = seqR[i][0];
        algnR_row[r_R/8] = 0;
	}
/*
//fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
for(i=0; i<=p_num; ++i){
    //fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 

}
*/

	return (int) p_num;
}
int align_255(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t seq[16], int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4] , struct SubBuf *sub)
{
    int8_t seq_len = aln_in[10];
    int8_t seq_st = aln_in[11];
    int8_t err_len = aln_in[12];     
    int8_t penalty = aln_in[15]; 
    int n_err; 
    int r_flag;
    if(err_len >0) r_flag = -2;
    else r_flag =-1;

    int stat_flag, mid_num, top_num, bot_num, bg_row, ed_row, row, rot_st, rot_ed;
    uint32_t bg_idx, ed_idx, pos, j;//bg_idx,ed_idx begin and end of 2nt-encoded bwt 
    uint8_t seq_ch, ch;
    int n_aln=0, n_exact_aln=0;
    pos = (seq_st+1)%8; //alignment sequence pos [0, 8)
    int rot = (7-seq_st)%8; //bw rotation
    int cur_len = 0; // current alignment length

    bg_idx = 0; ed_idx = n_data-1;//bwt interval [bg_idx, ed_idx]
    int aln_num=0, st_len = seq_len, err_num = 0;
    n_aln = aln_out[0][0]; 
    
    while(1){
        pos = (pos+7)%8; 
        seq_ch = seq[pos*2]<<2| seq[pos*2+1];
        rot = 7-pos;
        if(bg_idx == 0 && ed_idx==n_data-1){
            if(cnt_2nt[rot][seq_ch] == cnt_2nt[rot][seq_ch+1]){
                if(err_num >=1) {break;}
                err_num++; 
                st_len = (cur_len+1)%8;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
                continue;
            }
            bg_idx = cnt_2nt[rot][seq_ch]; 
            ed_idx = cnt_2nt[rot][seq_ch+1]-1;
            ++cur_len;
            if(cur_len == 8) {
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[n_aln][3] = 0;
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
                    r_flag = 1; 
                }
                return r_flag;
            } 
            continue; 
        }

        mid_num =0; top_num = 0; bot_num = 0;        
        int n_ch = cnt_2nt[rot][seq_ch+1]-cnt_2nt[rot][seq_ch];
        rot_st = rot*((n_data+1)/2);
        bg_row = rot_st+(bg_idx+1)/2;
        ed_row = rot_st+(ed_idx+1)/2;
        rot_ed = rot_st+n_data/2;//if n_data%2==0, row rot_ed is not in current rotation bwt
     
        if(bg_idx <= n_data-ed_idx){
            stat_flag = 1;
            for(row = rot_st; row<bg_row; ++row){
                ch = bwt[row] >>4;
                if(seq_ch==ch) ++top_num;
                ch = bwt[row] &0xF;
                if(seq_ch==ch) ++top_num; 
                if(n_ch == top_num) break; 
            }
        }
        if(bg_idx%2>0){
            ch = bwt[bg_row-1]&0xF;
            if(seq_ch == ch) {
                if(top_num >0)--top_num;
                ++mid_num;
            }
        }
        for(row = bg_row; row < ed_row; ++row){
            ch = bwt[row] >>4;
            if(seq_ch == ch) ++mid_num;
            ch = bwt[row]&0xF;
            if(seq_ch == ch) ++mid_num;
            if(n_ch == mid_num) break;
        }
        if(ed_idx%2==0){
            ch = bwt[ed_row]>>4;
            if(seq_ch==ch) ++mid_num;
            if(ed_idx < n_data-1){
                ch = bwt[ed_row]&0xF;
                if(seq_ch == ch) ++bot_num;
                ++ed_row;
            }
        }
        if(bg_idx > n_data-ed_idx){ 
            for(row=ed_row; row< rot_ed; ++row){
                ch = bwt[row]>>4;
                if(seq_ch==ch) ++bot_num; 
                ch = bwt[row]&0xF;
                if(seq_ch==ch) ++bot_num;
            }
            if(n_data%2>0 && ed_idx < n_data-1){
                ch = bwt[rot_ed]>>4;
                if(seq_ch==ch) {++bot_num; }
            }
        }
        
        if(bg_idx >n_data-ed_idx)  top_num = n_ch - bot_num-mid_num;
        if(err_len>0){
            ++cur_len;
            if(mid_num < 1 && cur_len <= 7 ) { return -2;}
            if(mid_num >=1 && cur_len < 7){
                if(cur_len <8-err_len){
                    bg_idx = top_num+cnt_2nt[rot][seq_ch];
                    ed_idx = bg_idx+mid_num-1;
                    //rot = (rot+1)%8;
                    continue;
                }
            }
        }
        if(err_len >0 && cur_len == 7 && mid_num >=1) {
            uint32_t i_idx, j_idx, j_row, j_rot, j_pos;
            int err_pos = (pos+7)%8;
            bg_idx = top_num+cnt_2nt[rot][seq_ch];
            ed_idx = bg_idx+mid_num-1;
            uint8_t err_ch = (seq[err_pos*2]<<2) | (seq[err_pos*2+1]);
            for(i_idx = bg_idx; i_idx<=ed_idx; ++i_idx){
                j_rot = (rot+1)%8;
                j_idx = i_idx;
                rot_st =j_rot*((n_data+1)/2);
                j_row = rot_st+(j_idx+1)/2;
                rot_ed = rot_st+n_data/2;
                if(j_idx%2>0) {
                    ch = bwt[j_row-1]&0xF; 
                } else{
                    ch = bwt[j_row]>>4;
                }
                
                n_err = 0; 
                //比较后两个bits
                if((ch & 3) != (err_ch &3)){ ++n_err; }
                //比较后3、4位置bit
                if((ch & 12) != (err_ch &12)){ ++n_err; }
                //仅保留单碱基匹配错误
                //if(i == 0 || i == 2) continue; 
                j_pos = pos;
                while(j_pos != 0){//[fixme?] rot != 0?
                    rot_st =j_rot*((n_data+1)/2);
                    j_row = rot_st+(j_idx+1)/2;
                    rot_ed = rot_st+n_data/2;
                    if(j_idx%2>0) {
                        seq_ch = bwt[j_row-1]&0xF; 
                    } else{
                        seq_ch = bwt[j_row]>>4;
                    }
                    n_ch = cnt_2nt[j_rot][seq_ch+1] - cnt_2nt[j_rot][seq_ch]; 
                    int num=0, k; 
                    if(j_idx <=(n_data+3)/4){
                        for(k = rot_st; k < j_row; ++k){
                            ch = bwt[k]>>4;
                            if(seq_ch == ch) ++num;
                            ch = bwt[k] &0xF;
                            if(seq_ch==ch) ++num;
                        }
                        if(j_idx%2==0){
                            ch = bwt[j_row]>>4;
                            if(seq_ch == ch) ++num;
                        }
                    }
                    
                    if(j_idx > (n_data+3)/4){
                        for(k = j_row; k < rot_ed; ++k){
                            ch = bwt[k]>>4;
                            if(seq_ch == ch){
                                if(k > j_row) ++num;
                                else if(j_idx%2>0) ++num;
                            }
                            ch = bwt[k]&0xF;
                            if(seq_ch==ch) ++num;
                        
                        } 
                        if(n_data%2>0){
                            if(j_idx < n_data -1){
                                ch = bwt[rot_ed]>>4;
                                if(seq_ch == ch) ++num;
                            }
                        } 
                        num = n_ch -num; 
                    }
                    j_idx = num + cnt_2nt[j_rot][seq_ch]-1;
                    --j_pos; 
                    j_rot  = (j_rot+1)%8; 
                }//end while(j_pos != 0)+++++++++++++++++++++++++ 
                int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = j_idx;
                    aln_out[n_aln][1] = penalty*n_err;
                    aln_out[n_aln][2] = err_pos;
                    aln_out[n_aln][3] = 0;
                    algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                    r_flag = 2;
                   

                }
            }//end for+++++++++++++++++++++++++++
            aln_out[0][0] = n_aln;
            aln_out[0][1] = n_exact_aln;
            return r_flag;
        }
        if(mid_num >0){
            ++cur_len;
            if(cur_len <8){
                bg_idx = top_num+cnt_2nt[rot][seq_ch];
                ed_idx = bg_idx+mid_num-1;
               // rot = (rot+1)%8;

            } else if(cur_len == 8){
                bg_idx = top_num+cnt_2nt[rot][seq_ch];
                ed_idx = bg_idx+mid_num-1;
if(mid_num > 1) {
    //fprintf(stderr, "line = %u, func = %s\n", __LINE__, __func__);
    printf("line = %u, func = %s\n", __LINE__, __func__);
    exit(1);
}
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    ++n_aln; ++n_exact_aln;
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[n_aln][3] = 0;
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] = n_exact_aln;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
                }
                return 1;
            }
        } else if(mid_num ==0){
            if(cur_len == 7){
                uint32_t i_idx, j_idx, j_row, j_rot, j_pos;
                int err_pos = pos;            
                uint8_t err_ch = (seq[err_pos*2]<<2) | (seq[err_pos*2+1]);
                for(i_idx = bg_idx; i_idx<=ed_idx; ++i_idx){
                    j_idx = i_idx;
                    j_pos = (pos+1)%8;
                    int j_rot = rot;
                    //++++++++++++++++++++++++++++++++++++
                    rot_st =j_rot*((n_data+1)/2);
                    j_row = rot_st+(j_idx+1)/2;
                    rot_ed = rot_st+n_data/2;
                    if(j_idx%2>0) {
                        ch = bwt[j_row-1]&0xF; 
                    } else{
                        ch = bwt[j_row]>>4;
                    }
                    n_err = 0; 
                    //比较后两个bits
                    if((ch & 3) != (err_ch &3)){ ++n_err; }
                    //比较后3、4位置bit
                    if((ch & 12) != (err_ch &12)){ ++n_err; }
                    //仅保留单碱基匹配错误
                    //------------------------------------    
                    while(j_pos != 0){//[fixme?] rot != 0?
                        rot_st =j_rot*((n_data+1)/2);
                        j_row = rot_st+(j_idx+1)/2;
                        rot_ed = rot_st+n_data/2;
                        if(j_idx%2>0) {
                            seq_ch = bwt[j_row-1]&0xF; 
                        } else{
                            seq_ch = bwt[j_row]>>4;
                        }
                        n_ch = cnt_2nt[j_rot][seq_ch+1] - cnt_2nt[j_rot][seq_ch]; 
                        int num=0, k; 
                        if(j_idx <=(n_data+3)/4){
                            for(k = rot_st; k < j_row; ++k){
                                ch = bwt[k]>>4;
                                if(seq_ch == ch) ++num;
                                ch = bwt[k] &0xF;
                                if(seq_ch==ch) ++num;
                            }
                            if(j_idx%2==0){
                                ch = bwt[j_row]>>4;
                                if(seq_ch == ch) ++num;
                            }
                        }
                        if(j_idx > (n_data+3)/4){
                            for(k = j_row; k < rot_ed; ++k){
                                ch = bwt[k]>>4;
                                if(seq_ch == ch){
                                    if(k > j_row) ++num;
                                    else if(j_idx%2>0) ++num;
                                }
                                ch = bwt[k]&0xF;
                                if(seq_ch==ch) ++num;
                            
                            } 
                            if(n_data%2>0){
                                if(j_idx < n_data -1){
                                    ch = bwt[rot_ed]>>4;
                                    if(seq_ch == ch) ++num;
                                }
                            } 
                            num = n_ch -num; 
                        }
                        j_idx = num + cnt_2nt[j_rot][seq_ch]-1;
                        --j_pos; 
                        j_rot  = (j_rot+1)%8; 
                    }//end while 
                    int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                    if(flag_aln == 0){
                        ++n_aln; 
                        aln_out[n_aln][0] = j_idx;
                        aln_out[n_aln][1] = penalty*n_err;
                        aln_out[n_aln][2] = err_pos;
                        aln_out[n_aln][3] = 0;                         
                        algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                        r_flag = 2;
                       

                    } 
                }//end for
                aln_out[0][0] = n_aln;
                aln_out[0][1] = n_exact_aln;
                return r_flag;
            } else if(cur_len <7){
                aln_num++;
                if(aln_num >7) break;
                if(pos == 0) break;
                if(cur_len + st_len <7) break; 
                err_num = 0;
                st_len = cur_len+st_len-7;
                aln_in[(pos+1)%8] =1;
                //rot=(rot+1)%8;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
           } 
        }        
    }//end while(1) 
    return r_flag;
}

int align_255_indel(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t seq[18], int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4] , struct SubBuf *sub )
{
    int8_t seq_len = aln_in[10];
    int8_t seq_st  = aln_in[11];
    int8_t err_len = aln_in[12];     
    int8_t InDel   = aln_in[13];
    int8_t penalty = aln_in[15]; 
    int stat_flag, mid_num, top_num, bot_num, bg_row, ed_row, row, rot_st, rot_ed;
    uint32_t bg_idx, ed_idx, pos, j;//bg_idx,ed_idx begin and end of 2nt-encoded bwt 
    uint8_t seq_ch, ch;
    int n_aln=0, n_exact_aln=0;
    pos = (seq_st+1)%8; //alignment sequence pos [0, 8)
    int rot = (7-seq_st)%8; //bw rotation
    int cur_len = 0; // current alignment length
    if(InDel == 0 || InDel < -2 || InDel > 2 ) {
        printf("%u, %s, InDel == %d , error!!\n", __LINE__, __func__, InDel);
        exit(1); 
    }
    if(n_data <= 1){
        printf("%u, %s, n_data == %u , error!!\n", __LINE__, __func__, n_data);
        exit(1); 
    }
    int r_flag;
    if(err_len > 0) r_flag = -2;
    else r_flag = -1; 

    bg_idx = 0; ed_idx = n_data-1;//bwt interval [bg_idx, ed_idx]
    int aln_num=0, st_len = seq_len, err_num = 0;
    n_aln = aln_out[0][0]; 
   
   
    int edit_indel[4];
    int seq_pos = 0;
    if(err_len ==0 ) {
        pos = 0;
        seq_pos = 0; 
    } else{
        seq_pos = pos*2; 
    }

    int len_seq = 16 + InDel;
    while(1){
        pos = (pos+7)%8; 
        seq_pos = (seq_pos+len_seq-2) % len_seq;
        seq_ch = seq[seq_pos]<<2| seq[seq_pos+1];
        rot = 7-pos;
        if(bg_idx == 0 && ed_idx==n_data-1){
            if(cnt_2nt[rot][seq_ch] == cnt_2nt[rot][seq_ch+1]){
                if(err_num >=1) { break;}//goto return
                err_num++; 
                seq_pos = 2 * pos;
                if(cur_len > 0) st_len = cur_len+st_len-7;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
                continue;//while(1)
            }
            bg_idx = cnt_2nt[rot][seq_ch]; 
            ed_idx = cnt_2nt[rot][seq_ch+1]-1;
            ++cur_len;
            if(cur_len == 8) {
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[n_aln][3] = InDel;
                     
                    
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
                    r_flag = 1;
                }
                return r_flag;
            }             
            continue;//while(1) 
        }
        mid_num =0; top_num = 0; bot_num = 0;        
        int n_ch = cnt_2nt[rot][seq_ch+1]-cnt_2nt[rot][seq_ch];
        rot_st = rot*((n_data+1)/2);
        bg_row = rot_st+(bg_idx+1)/2;
        ed_row = rot_st+(ed_idx+1)/2;
        rot_ed = rot_st+n_data/2;//if n_data%2==0, row rot_ed is not in current rotation bwt
        if(bg_idx <= n_data-ed_idx){
             stat_flag = 1;
            for(row = rot_st; row<bg_row; ++row){
                ch = bwt[row] >>4;
                if(seq_ch==ch) ++top_num;
                ch = bwt[row] &0xF;
                if(seq_ch==ch) ++top_num; 
                if(n_ch == top_num) break;//for(rpw = rot_st; row<bg_row; ++row) 
            }
        }
        if(bg_idx%2>0){
            ch = bwt[bg_row-1]&0xF;
            if(seq_ch == ch) {
                if(top_num >0)--top_num;
                ++mid_num;
            }
        }
        for(row = bg_row; row < ed_row; ++row){

            ch = bwt[row] >>4;
            if(seq_ch == ch) ++mid_num;
            ch = bwt[row]&0xF;
            if(seq_ch == ch) ++mid_num;
            if(n_ch == mid_num) break;  //for(row = bg_row; row < ed_row; ++row)
        }
        if(ed_idx%2==0){
            ch = bwt[ed_row]>>4;
            if(seq_ch==ch) ++mid_num;
            if(ed_idx < n_data-1){
                ch = bwt[ed_row]&0xF;
                if(seq_ch == ch) ++bot_num;
                ++ed_row;
            }
        }
        if(bg_idx > n_data-ed_idx){ 
            for(row=ed_row; row< rot_ed; ++row){
                ch = bwt[row]>>4;
                if(seq_ch==ch) ++bot_num; 
                ch = bwt[row]&0xF;
                if(seq_ch==ch) ++bot_num;
            }
            if(n_data%2>0 && ed_idx < n_data-1){
                ch = bwt[rot_ed]>>4;
                if(seq_ch==ch) {++bot_num; }
            }
        }
        if(bg_idx >n_data-ed_idx)  top_num = n_ch - bot_num-mid_num;
        ++cur_len;
////fprintf(stderr, "%u,err_len = %u, pos = %u, seq_pos = %u, cur_len = %u, mid_num = %u\n", __LINE__, err_len, pos, seq_pos, cur_len, mid_num);
        if(mid_num >=1 && cur_len < 7){
            bg_idx = top_num+cnt_2nt[rot][seq_ch];
            ed_idx = bg_idx+mid_num-1;
            continue;//while(1)
        }
        if(err_len>0 && mid_num < 1 && cur_len <= 7 ) { return -2;}
        if(cur_len == 7 && mid_num >=1) {
            uint32_t i_idx, j_idx, j_row, j_pos, j_rot;
            int err_pos = (pos+7)%8;
            bg_idx = top_num+cnt_2nt[rot][seq_ch];
            ed_idx = bg_idx+mid_num-1;
            uint8_t err_ch = (seq[err_pos*2]<<2) | (seq[err_pos*2+1]);
            if(InDel == -1) seq_pos = (seq_pos+len_seq-1)%len_seq;
            else if(InDel == 1) seq_pos = (seq_pos+len_seq-3)%len_seq;
            else if(InDel == 2) seq_pos = (seq_pos+len_seq-4)%len_seq;
////fprintf(stderr, "%u, bg_idx = %u, ed_idx = %u, InDel = %d\n", __LINE__, bg_idx, ed_idx, InDel);
            for(i_idx = bg_idx; i_idx<=ed_idx; ++i_idx){
                j_rot = (rot+1)%8;
                j_idx = i_idx;
                rot_st =j_rot*((n_data+1)/2);
                j_row = rot_st+(j_idx+1)/2;
                rot_ed = rot_st+n_data/2;
                if(j_idx%2>0) {
                    ch = bwt[j_row-1]&0xF; 
                } else{
                    ch = bwt[j_row]>>4;
                }
                edit_indel[0] = InDel;
                if(InDel < 0){ 
                    if(InDel == -2){
                        edit_indel[1] = seq_pos;
                        edit_indel[3] = ch;

                    } else if(InDel == -1){
                        if(seq[seq_pos] == (ch&3) ) {
                            edit_indel[1] = (seq_pos+len_seq-1)%len_seq;
                            edit_indel[3] = ((ch>>2)&3);
                        } else if(seq[seq_pos] == (ch>>2) ) {
                            edit_indel[1] = seq_pos;
                            edit_indel[3] = (ch&3);
                        } else{
                            continue;//for 
                        }
                    } 
                } //end if(InDel < 0)  
                int edit_buf[4];
                if(InDel >0) {
                    if(InDel == 1){
                        edit_buf[0] = (seq[seq_pos]<<2)|(seq[seq_pos+1]);
                        edit_buf[1] = (seq[seq_pos]<<2)|(seq[seq_pos+2]);
                        edit_buf[2] = (seq[seq_pos+1]<<2)|(seq[seq_pos+2]);
                        int i;
                        for(i = 0; i < 3; ++i){
                            if(ch == edit_buf[i]) break;//for(i) 
                        }

                        if(i == 3) {continue;}//while(1) 
                        else{ 
                            if(i == 0) {                            
                                edit_indel[1] = (seq_pos+2)%len_seq;
                            } else if(i == 1) {                            
                                edit_indel[1] = (seq_pos+1)%len_seq;
                            } else if(i == 2) {                            
                                edit_indel[1] = (seq_pos+0)%len_seq;
                            }
                            edit_indel[3] = seq[edit_indel[1]]; 
                        }
                    } else if(InDel == 2){
                        edit_buf[0] = (seq[seq_pos]<<2)|(seq[seq_pos+1]);
                        edit_buf[1] = (seq[seq_pos]<<2)|(seq[seq_pos+3]);
                        edit_buf[2] = (seq[seq_pos+2]<<2)|(seq[seq_pos+3]);
                        int i;
                        for(i = 0; i < 3; ++i){
                            if(ch == edit_buf[i]) break; 
                        }

                        if(i == 3) {continue;}//while(1) 
                        else{
                            if(i == 0) {                            
                                edit_indel[1] = (seq_pos+2)%len_seq;
                            } else if(i == 1) {                            
                                edit_indel[1] = (seq_pos+1)%len_seq;
                            } else if(i == 2) {                            
                                edit_indel[1] = (seq_pos+0)%len_seq;
                            }
                            edit_indel[3] = (seq[edit_indel[1]]<<2) | (seq[edit_indel[1]+1]);
                        }
                    } //end if(InDel == 1)
                } //end if(InDel >0)
                j_pos = pos;
                while(j_pos != 0){//[fixme?] rot != 0?
                    rot_st =j_rot*((n_data+1)/2);
                    j_row = rot_st+(j_idx+1)/2;
                    rot_ed = rot_st+n_data/2;
                    if(j_idx%2>0) {
                        seq_ch = bwt[j_row-1]&0xF; 
                    } else{
                        seq_ch = bwt[j_row]>>4;
                    }
                    n_ch = cnt_2nt[j_rot][seq_ch+1] - cnt_2nt[j_rot][seq_ch]; 
                    int num=0, k; 
                    if(j_idx <=(n_data+3)/4){
                        for(k = rot_st; k < j_row; ++k){
                            ch = bwt[k]>>4;
                            if(seq_ch == ch) ++num;
                            ch = bwt[k] &0xF;
                            if(seq_ch==ch) ++num;
                        }
                        if(j_idx%2==0){
                            ch = bwt[j_row]>>4;
                            if(seq_ch == ch) ++num;
                        }
                    }
                    
                    if(j_idx > (n_data+3)/4){
                        for(k = j_row; k < rot_ed; ++k){
                            ch = bwt[k]>>4;
                            if(seq_ch == ch){
                                if(k > j_row) ++num;
                                else if(j_idx%2>0) ++num;
                            }
                            ch = bwt[k]&0xF;
                            if(seq_ch==ch) ++num;
                        
                        } 
                        if(n_data%2>0){
                            if(j_idx < n_data -1){
                                ch = bwt[rot_ed]>>4;
                                if(seq_ch == ch) ++num;
                            }
                        } 
                        num = n_ch -num; 
                    }
                    j_idx = num + cnt_2nt[j_rot][seq_ch]-1;
                    --j_pos; 
                    j_rot  = (j_rot+1)%8; 
                }//end while(j_pos != 0)+++++++++++++++++++++++++ 
                int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = j_idx;
                    aln_out[n_aln][1] = penalty;
                    aln_out[n_aln][2] = err_pos;
                    aln_out[n_aln][3] = InDel;
                    algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                    if(err_len > 0) r_flag = 2;
                    else r_flag = 1; 
                }
            }//end for+++++++++++++++++++++++++++
            if(aln_out[0][0] == n_aln && err_len == 0) {
                mid_num = 0;
                cur_len--;
            }else{
                aln_out[0][0] = n_aln;
                aln_out[0][1] = n_exact_aln;
                break; 
            }
        } //if(cur_len == 7 && mid_num >=1)+++++++++++
        if(err_len == 0 && mid_num ==0 && cur_len <=7){
            aln_num++;
            if(aln_num >7) break; //goto return 
            if(cur_len + st_len <7) break; //goto return 
            if(pos == 0) {
                if(cur_len == 7 ) pos++;
                else {
                    pos = 7;
                    aln_num = 7; 
                }
            }
            seq_pos = pos*2;
            st_len = cur_len+st_len-7;
            aln_in[(pos+1)%8] =1;
            err_num = 0;
            cur_len = 0;
            bg_idx = 0;
            ed_idx = n_data-1; 
            continue;//while(1)
        } 
    }//end while(1) 
    return r_flag;
}
//+++++++++++++++++++++++++++++++++++++++++++++++
/*
uint32_t __i, __j;
for(__i = 0; __i < 8; ++__i){
    //fprintf(stderr, "%u \t", __i);
    for(__j = 0; __j < 17; ++__j)
        //fprintf(stderr, "%u \t", cnt_2nt[__i][__j]);
    //fprintf(stderr, "\n");

}
*/
//++++++++++++++++++++++++++++++++++++++++++++++++

static inline void accumulate_cnt_uint8(int n, uint8_t *cnt_2nt)
{
    int i;
    for(i = 1; i < n-1; ++i) cnt_2nt[i] += cnt_2nt[i-1]; 
    for(i = n-1; i >0; --i) cnt_2nt[i] = cnt_2nt[i-1];
    cnt_2nt[0] = 0;
 
}
void gen_cnt(uint8_t *bwt, int n_data, uint8_t cnt_2nt[][17])
{
    
    uint8_t rot, ch;
    uint32_t bg_row, ed_row, i;
    for(rot=0; rot<8; ++rot){
        for(ch=0; ch<17; ++ch) {cnt_2nt[rot][ch] = 0;} 
    }
    for(rot =0; rot<8; ++rot){
        bg_row = rot*((n_data+1)/2);
        ed_row = bg_row +n_data/2; 
        for(i = bg_row; i < ed_row; ++i){
            ch = bwt[i]>>4;
            ++cnt_2nt[rot][ch];
            ch = bwt[i]&0xF;
            ++cnt_2nt[rot][ch]; 

        }
        if(n_data%2 != 0){
            ch = bwt[ed_row]>>4;
            ++cnt_2nt[rot][ch];
        }
        uint32_t __j;
        accumulate_cnt_uint8(17, cnt_2nt[rot]);
   
    }
}


uint32_t __dna2_count(uint8_t *bwt, uint32_t k , uint32_t l, uint8_t c){
    uint32_t i, n;
    for(n=0, i =k; i <l; ++i){
        n += (__get_bwt2(bwt, i)==c)?1:0;
    }
    return n;

}
uint32_t __dna2_count_small(uint8_t *bwt, uint32_t k , uint32_t l, uint8_t c){
    uint32_t i, n;
    for(n=0, i =k; i <l; ++i){
        n += (__get_bwt2(bwt, i)<c)?1:0;
    }
    return n;

}
int align_min(uint8_t *bwt, int n_data, uint8_t seq[16],int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4], struct SubBuf *sub )
{
    int8_t err_num = aln_in[10];
    int8_t flag    = aln_in[11]; 
    int8_t err_len = aln_in[12];
    int8_t penalty = aln_in[15]; 
    int i, j;
    int n_aln = 0, n_exact_aln = 0;
    n_aln = aln_out[0][0];
    uint8_t *sort_seq = (uint8_t *)bwt; 
    uint32_t row;
    uint8_t seq_buf;
    int n_diff, err_pos = 0;
    int st_row = 0;
    int r_flag;
    if(err_len > 0) r_flag = -2;
    else r_flag = -1; 
    if(err_len == 0){
        uint8_t in_seq[4];
        for(i =0; i < 4; ++i){
            in_seq[i] = seq[i*4]; 
            in_seq[i] = (in_seq[i] <<2)| (seq[i*4+1]); 
            in_seq[i] = (in_seq[i] <<2)| (seq[i*4+2]); 
            in_seq[i] = (in_seq[i] <<2)| (seq[i*4+3]); 
        } 
        int bg = 0, end = n_data-1;
     
        while(1){
            if(bg > end){break; }
            row = (bg+end)/2;        
            for(i =0;i < 4; ++i){
                if(in_seq[i] > sort_seq[4*row+i]){
                    bg = row+1;
                    break; 
                } else if(in_seq[i] < sort_seq[4*row+i]){
                    end = row-1; 
                    break; 
                }             
            } 
            if(i == 4){//find exact match
                n_aln++;
                aln_out[n_aln][0] = row;
                aln_out[n_aln][1] = 0;
                aln_out[n_aln][2] = err_pos;
                aln_out[n_aln][3] = 0; 
                 
                aln_out[0][0] = n_aln;
                algn_row[row/8] |= 1<<(7-row%8);
                r_flag = 1;
                break;//while(1)
            }
        } // while(1) +++++++
        aln_out[n_data+1][0] = 0; 
        if(r_flag >0) return r_flag;
        //精确比对段结束-----------------------------------------------------
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //精确比对失败时，进一步判断是否存在近似比对,存在时找到一个
        for(row = 0; row < n_data; row++){ 
            n_diff = 0; 
            for(i = 0; i < 4; ++i){
                for(j = 0; j < 4; ++j){
                    seq_buf = (sort_seq[row*4+i]>>((3-j)*2))&3;
                    if(seq_buf != seq[i*4+j] ) {
                        ++n_diff;
                        err_pos = i*4+j;
                    }
                    if(n_diff > err_num) break;
                }
                if(n_diff > err_num) break;
            }
            if(n_diff <= err_num){
                int flag_aln = algn_row[row/8] & (1<<(7-row%8));
                if(flag_aln == 0){
                    n_aln++;
                    aln_out[n_aln][0] = row;
                    aln_out[n_aln][1] = penalty*n_diff;
                    aln_out[n_aln][2] = err_pos;
                    aln_out[n_aln][3] = 0; 
                     
                    
                    aln_out[0][0] = n_aln;
                    algn_row[row/8] |= 1<<(7-row%8);
                    r_flag = 2;
                    if(n_diff >0) aln_in[err_pos/2] = ++aln_in[8];

                    if(n_diff == 0) {
                        aln_out[0][1] = 1;
                        for(i =0; i < 3; ++i){
                            uint32_t tmp;
                            tmp = aln_out[n_aln][i];
                            aln_out[n_aln][i] = aln_out[1][i];
                            aln_out[1][i] = tmp;    
                        }

                        if(flag > 0) { break;} 
                    } 
                }// end if(flag_aln == 0)+++++++++++++++++++++++++++++++
            }
        } // for(row = 0; row < n_data; row++) ++++++++++++++
        aln_out[0][0] = n_aln; 
        aln_out[0][2] = n_aln - aln_out[0][1];
        aln_out[n_data+1][0] = row+1; 
        return r_flag;   
    }//end if(err_len == 0)-------------------------------------------------

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //全局近似比对查找全部
    if(aln_out[0][0] >0){
        st_row = aln_out[n_data+1][0]; 
    }     
    for(row = st_row; row < n_data; row++){ 
        n_diff = 0; 
        for(i = 0; i < 4; ++i){
            for(j = 0; j < 4; ++j){
                seq_buf = (sort_seq[row*4+i]>>((3-j)*2))&3;
                if(seq_buf != seq[i*4+j] ) {
                    ++n_diff;
                    err_pos = i*4+j;
                }
                if(n_diff > err_num) break;
            }
            if(n_diff > err_num) break;
        }
        if(n_diff <= err_num){
            int flag_aln = algn_row[row/8] & (1<<(7-row%8));
            if(flag_aln == 0){
                n_aln++;
                aln_out[n_aln][0] = row;
                aln_out[n_aln][2] = err_pos;
                aln_out[n_aln][3] = 0;
                 
                
                aln_out[0][0] = n_aln;
                algn_row[row/8] |= 1<<(7-row%8);
                if(err_len >0) {
                    aln_out[n_aln][1] = penalty*n_diff;
                    r_flag = 2;
                }else {
                    aln_out[n_aln][1] = 0;
                    r_flag = 1;
                }
                    
                if(n_diff >0) aln_in[err_pos/2] = ++aln_in[8];
                if(n_diff == 0) {
                    aln_out[0][1] = 1;
                    for(i =0; i < 3; ++i){
                        uint32_t tmp;
                        tmp = aln_out[n_aln][i];
                        aln_out[n_aln][i] = aln_out[1][i];
                        aln_out[1][i] = tmp;    
                    }

                    if(flag > 0) { break;} 
                } 
            }// end if(flag_aln == 0)+++++++++++++++++++++++++++++++
        }
    }
    aln_out[0][0] = n_aln; 
    aln_out[0][2] = n_aln - aln_out[0][1];
    aln_out[n_data+1][0] = row+1; 
    return r_flag;   
    
}

int align_min_indel(uint8_t *bwt, int n_data, uint8_t seq[18],int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4] , struct SubBuf *sub)
{
    int8_t err_num = aln_in[10];
    int8_t flag    = aln_in[11];
    int8_t err_len = aln_in[12];
    int8_t InDel   = aln_in[13];
    int8_t penalty = aln_in[15]; 
    int r_flag; 
    if(err_len > 0) r_flag = -2;
    else r_flag = -1; 

    int i, j;
    int n_aln = 0, n_exact_aln = 0;
    n_aln = aln_out[0][0];
    uint8_t *sort_seq = (uint8_t *)bwt; 
    uint32_t row;
    uint8_t seq_buf;
    int n_diff, err_pos = 0;
    int st_row = 0;
     //判断是否存在插入删除比对,存在时只找一个答案。
    int len_L = 0, len_R = 0;
    if(err_len == 0) {//运行indel_once
        for(row = 0; row < n_data; row++){ 
            len_L = 0, len_R = 0;
            n_diff = 0; 
            for(i = 0; i < 4; ++i){
                for(j = 0; j < 4; ++j){
                    seq_buf = (sort_seq[row*4+i]>>((3-j)*2))&3;
                    if(seq_buf != seq[i*4+j] ) {
                        ++n_diff;
                        break;//for(j =0; j <4; ++j)
                    }
                    len_L++;              
                }
                if(n_diff > err_num) break;// for(i = 0; i < 4; ++i)
            }
            n_diff = 0;
            for(i = 3; i >= 0; --i){
                for(j = 3; j >= 0; --j){
                    seq_buf = (sort_seq[row*4+i]>>((3-j)*2))&3;
if(i*4+j+InDel < 0){
    //printf("i*4+j+InDel = %d\n", i*4+j+InDel);
    break;
    //exit(1);
}
                    if(seq_buf != seq[(i*4+j)+InDel] ) {
                        ++n_diff;
                        break;//for(j =0; j <4; ++j)
                    }
                    len_R++;
                }
                if(n_diff > err_num) break;// for(i = 0; i < 4; ++i)
            }
            if((InDel < 0 && len_L+len_R >= 16+InDel) || (InDel > 0 && len_L+len_R >= 16)){
                n_aln++;
                aln_out[n_aln][0] = row;
                aln_out[n_aln][1] = penalty;
                aln_out[n_aln][3] = InDel;  
                algn_row[row/8] |= 1<<(7-row%8);
                r_flag = 1;
                break;
            } 
        }
        aln_out[0][0] = n_aln; 
        aln_out[0][2] = n_aln - aln_out[0][1];
        if(aln_in[13] == -1 ){
            aln_in[16] = row+1;        
        }else if(aln_in[13] == 1){
            aln_in[17] = row+1;        
        }else if(aln_in[13] == -2){
            aln_in[18] = row+1;   
        }else if(aln_in[13] == 2){
            aln_in[19] = row+1;        
        }else{
            printf("%s, %u, Error: InDel = %d!\n",__func__, __LINE__, aln_in[13]);
            exit(1);
        }
        return r_flag; 
    }
    //存在插入删除时找出全部答案
    if(aln_out[0][0] >0){
        if(aln_in[13] == -1 ){
            st_row = aln_in[16];        
        }else if(aln_in[13] == 1){
            st_row = aln_in[17];        
        }else if(aln_in[13] == -2){
            st_row = aln_in[18];   
        }else if(aln_in[13] == 2){
            st_row = aln_in[19];        
        }else{
            printf("%s, %u, Error: InDel = %d!\n",__func__, __LINE__, aln_in[13]);
            exit(1);
        }
    }     
    for(row = st_row; row < n_data; row++){ 
        len_L = 0, len_R = 0;
        n_diff = 0; 
        for(i = 0; i < 4; ++i){
            for(j = 0; j < 4; ++j){
                seq_buf = (sort_seq[row*4+i]>>((3-j)*2))&3;
                if(seq_buf != seq[i*4+j] ) {
                    ++n_diff;
                    break;//for(j =0; j <4; ++j)
                }
                len_L++;
            }
            if(n_diff > err_num) break;// for(i = 0; i < 4; ++i)
        }
        n_diff = 0;
        for(i = 3; i >= 0; --i){
            for(j = 3; j >= 0; --j){
                seq_buf = (sort_seq[row*4+i]>>((3-j)*2))&3;
if(i*4+j+InDel < 0){
    //printf("i*4+j+InDel = %d\n", i*4+j+InDel);
    break;
    //exit(1);
}
 
                if(seq_buf != seq[(i*4+j)+InDel] ) {
                    ++n_diff;
                    break;//for(j =0; j <4; ++j)
                }
                len_R++;
            }
            if(n_diff > err_num) break;// for(i = 0; i < 4; ++i)
        }
        
        if((InDel < 0 && len_L+len_R >= 16+InDel) || (InDel > 0 && len_L+len_R >= 16)){
            n_aln++;
            aln_out[n_aln][0] = row;
            aln_out[n_aln][1] = penalty;
            aln_out[n_aln][2] = err_pos;
            aln_out[n_aln][3] = InDel; 
            aln_out[0][0] = n_aln;
            algn_row[row/8] |= 1<<(7-row%8);
            r_flag = 2;
        } 
    }
    
    aln_out[0][0] = n_aln; 
    aln_out[0][2] = n_aln - aln_out[0][1];
    if(aln_in[13] == -1 ){
        aln_in[16] = row+1;        
    }else if(aln_in[13] == 1){
        aln_in[17] = row+1;        
    }else if(aln_in[13] == -2){
        aln_in[18] = row+1;   
    }else if(aln_in[13] == 2){
        aln_in[19] = row+1;        
    }else{
        printf("%s, %u, Error: InDel = %d!\n",__func__, __LINE__, aln_in[13]);
        exit(1);
    }
    return r_flag;   
}

int get_ch_nosum( int idx, int pos, uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t *ret_ch)
{
    int rot_st, rot_ed, i;
    uint8_t seq_ch, ch;
    uint32_t j_idx, j_row;
    j_idx = idx;
    int cur_rot = 7-pos;
  
    rot_st =cur_rot*((n_data+1)/2);
    j_row = rot_st+(j_idx+1)/2;
    rot_ed = rot_st+n_data/2;
    if(j_idx%2>0) {
        seq_ch = bwt[j_row-1]&0xF; 
    } else{
        seq_ch = bwt[j_row]>>4;
    }
    *ret_ch = seq_ch; 
    uint32_t n_ch = cnt_2nt[cur_rot][seq_ch+1] - cnt_2nt[cur_rot][seq_ch]; 
    int num=0, k; 
    if(j_idx <=(n_data+3)/4){
        for(k = rot_st; k < j_row; ++k){
            ch = bwt[k]>>4;
            if(seq_ch == ch) ++num;
            ch = bwt[k] &0xF;
            if(seq_ch==ch) ++num;
        }
        if(j_idx%2==0){
            ch = bwt[j_row]>>4;
            if(seq_ch == ch) ++num;
        }
    }
    
    if(j_idx > (n_data+3)/4){
        for(k = j_row; k < rot_ed; ++k){
            ch = bwt[k]>>4;
            if(seq_ch == ch){
                if(k > j_row) ++num;
                else if(j_idx%2>0) ++num;
            }
            ch = bwt[k]&0xF;
            if(seq_ch==ch) ++num;
        } 
        if(n_data%2>0){
            if(j_idx < n_data -1){
                ch = bwt[rot_ed]>>4;
                if(seq_ch == ch) ++num;
            }
        } 
        num = n_ch -num; 
    }
    j_idx = num + cnt_2nt[cur_rot][seq_ch]-1;
    
    return j_idx; 
}
int get_ch_255(int idx, int pos, uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t *ret_ch)
{
    int rot_st, rot_ed;
    uint32_t j;//bg_idx,ed_idx begin and end of 2nt-encoded bwt 
    uint8_t seq_ch, ch;
    uint32_t j_idx, j_row;
    j_idx = idx;
    int cur_rot = (7-pos);

    rot_st =cur_rot*((n_data+1)/2);
    j_row = rot_st+(j_idx+1)/2;
    rot_ed = rot_st+n_data/2;
    if(j_idx%2>0) {
        seq_ch = bwt[j_row-1]&0xF; 
    } else{
        seq_ch = bwt[j_row]>>4;
    }
    *ret_ch = seq_ch;
    int n_ch = cnt_2nt[cur_rot][seq_ch+1] - cnt_2nt[cur_rot][seq_ch]; 
    
    //fprintf(stderr, "rot_off=%u, seq_ch = %u\n",rot_st, seq_ch);
    //fprintf(stderr, "j_idx=%u, j_row = %u\n",j_idx, j_row);
    int num=0, k; 
    if(j_idx <=(n_data+3)/4){
        for(k = rot_st; k < j_row; ++k){
            ch = bwt[k]>>4;
            if(seq_ch == ch) ++num;
            ch = bwt[k] &0xF;
            if(seq_ch==ch) ++num;
        }
        if(j_idx%2==0){
            ch = bwt[j_row]>>4;
            if(seq_ch == ch) ++num;
        }
    }
    
    if(j_idx > (n_data+3)/4){
        for(k = j_row; k < rot_ed; ++k){
            ch = bwt[k]>>4;
            if(seq_ch == ch){
                if(k > j_row) ++num;
                else if(j_idx%2>0) ++num;
            }
            ch = bwt[k]&0xF;
            if(seq_ch==ch) ++num;
        
        } 
        if(n_data%2>0){
            if(j_idx < n_data -1){
                ch = bwt[rot_ed]>>4;
                if(seq_ch == ch) ++num;
            }
        } 
        num = n_ch -num; 
    }
    j_idx = num + cnt_2nt[cur_rot][seq_ch]-1;
    return j_idx;
}
int get_ch_large_alt(int idx, int pos, uint8_t *Bwt, uint8_t *cnt2, uint32_t n_data, uint8_t *ret_ch){
    uint8_t  *pBwt, *pSum;
    uint32_t sum_buf;
    uint32_t idx_i[2], idx_q[2] , idx_r[2] , idx_row[2], rot_off;
    uint32_t glb_cnt[2] = {}; 
    uint32_t lst_row, ch_off, rng_cnt[2], rng_end[2], ch_cnt[2], rng_len, seq_num[2];
    
    int i, j, k;
    pBwt =Bwt ;
    pSum = cnt2 ;
    int siz_cnt;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    if( n_data <= 254*256 ) siz_cnt = 2 ;
    if( n_data > 254*256 ) siz_cnt = 4 ;
    int rng_num = ((n_data+253) /254) ;
    int len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;
    int data_q = (n_data)/254;
    int data_r = (n_data)%254;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // bg_idx 和ed_idx是下一个轮换内部的相对位置+++++


    uint32_t j_idx = idx;    
    uint32_t rot= 7-pos;
    uint8_t seq_ch, ch, get_ch;
    rot_off = rot*((n_data+1)/2);
    uint32_t j_row = rot_off +( j_idx+1)/2;
    lst_row = rot_off  + n_data/2 ;
    if(j_idx%2>0)  {
        get_ch =pBwt[j_row-1];
        seq_ch =  get_ch&0xF;
    } else if(j_idx%2==0)  {
        get_ch =pBwt[j_row];
        seq_ch  =  get_ch>>4;
    }
    *ret_ch = seq_ch;
    rot_off = rot*((n_data+1)/2);
    ch_off = 16*rot*len_sum + len_sum*seq_ch ;
    idx_i[0] = j_idx;
    idx_q[0] = idx_i[0]/254;
    idx_r[0] = idx_i[0]%254;
    idx_row[0] = rot_off +( idx_i[0]+1)/2;
    rng_cnt[0] = (8+siz_cnt)*(idx_q[0] / 8);
    rng_end[0] = rng_cnt[0] + idx_q[0] % 8;
  

    
    ch_cnt[0] = 0;
    if(rng_cnt[0]>0){
        for(j = 0; j< siz_cnt; j++){
            ch_cnt[0]=(ch_cnt[0]<<8)|(uint32_t)(*(pSum+ch_off+rng_cnt[0]- siz_cnt +j)); 
        }
    }
    for(k= rng_cnt[0]; k< rng_end[0]; k++){ 
        ch_cnt[0] += (uint32_t)(*(pSum+ ch_off + k ));
    }       
    seq_num[0]  = (uint32_t)(*(pSum+ ch_off + rng_end[0] ));
    glb_cnt[0] = 0 ;
    if(seq_ch > 0){
        for(j = 0; j<siz_cnt; j++) 
            glb_cnt[0] = (glb_cnt[0]<<8) | (uint32_t)(*(pSum+ch_off - siz_cnt+j)) ; 
    }
    ch_cnt[0] += glb_cnt[0];
    uint32_t bot_num =(uint32_t)pSum[ch_off +rng_end[0]];
    rng_len = 254 ;
    if(idx_q[0]==data_q) rng_len = data_r;
    rng_end[0] = rot_off+ (idx_i[0] - idx_r[0])/2 + rng_len/2;
    for(k= idx_row[0]; k< rng_end[0] ; k++){
        get_ch = pBwt[k] ;
        ch = get_ch>>4;
        if(seq_ch ==ch) { 
            if(k > idx_row[0] || idx_r[0]%2 >0 ){--bot_num;}
        }
        ch = get_ch&0xF ;
        if(seq_ch ==ch)  --bot_num;
    }
    if((rng_len %2 >0 ) &&(rng_len-1 > idx_r[0])){  
        get_ch = pBwt[rng_end[0]];
        ch = get_ch>>4;
        if(seq_ch == ch) --bot_num; 
    }
    bot_num += ch_cnt[0];
    j_idx = bot_num-1;


    return j_idx;
} // End : Alignm_big( )-
int align_nosum(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t seq[16], int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4], struct SubBuf *sub)
{
    int i;
    int8_t seq_len = aln_in[10];
    int8_t seq_st = aln_in[11]; 
    int8_t err_len = aln_in[12]; 
    int8_t penalty = aln_in[15]; 
    int n_err;
    int r_flag; 
    if(err_len > 0) r_flag = -2;
    else r_flag = -1; 
    if(n_data < 2) {
        printf("%u, %s, call func error: n_data < 2!\n", __LINE__, __func__);
        exit(1);
    } 
    int stat_flag, mid_num, top_num, bot_num, bg_row, ed_row, row, rot_st, rot_ed;
    uint32_t bg_idx, ed_idx, j;//bg_idx,ed_idx begin and end of 2nt-encoded bwt 
    uint8_t seq_ch, ch;
    int n_aln = 0, n_exact_aln = 0;
    n_aln = aln_out[0][0];
    uint32_t pos = (seq_st+1)%8; 
    int rot = (7-seq_st)%8; //bw rotation
    int cur_len = 0; // current alignment length
    int st_len = seq_len;
 
    bg_idx = 0; ed_idx = n_data-1;//bwt interval [bg_idx, ed_idx]
    n_aln = aln_out[0][0];
    int aln_num = 0, err_num = 0;
    while(1){
        pos = (pos+7)%8; 
        seq_ch = seq[pos*2]<<2| seq[pos*2+1];
        rot = 7-pos;
        if(bg_idx == 0 && ed_idx==n_data-1){
            if(cnt_2nt[rot][seq_ch] == cnt_2nt[rot][seq_ch+1]){
                if(err_num >=1) {break;}
                err_num++; 
                st_len = (cur_len+1)%8;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
                continue;
            }
            bg_idx = cnt_2nt[rot][seq_ch]; 
            ed_idx = cnt_2nt[rot][seq_ch+1]-1;
            ++cur_len;
            if(cur_len == 8) {
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    ++n_aln;
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[n_aln][3] = 0;                     
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
                    r_flag = 1;
                }
                return r_flag;
            } 
            continue; 
        } //end else if(bg_idx == 0 && ed_idx==n_data-1)+++++++++++++++++
        mid_num =0; top_num = 0; bot_num = 0;        
        int n_ch = cnt_2nt[rot][seq_ch+1]-cnt_2nt[rot][seq_ch];
        rot_st = rot*((n_data+1)/2);
        bg_row = rot_st+(bg_idx+1)/2;
        ed_row = rot_st+(ed_idx+1)/2;
        rot_ed = rot_st+n_data/2;//if n_data%2==0, row rot_ed is not in current rotation bwt
     
        if(bg_idx <= n_data-ed_idx){
            stat_flag = 1;
            for(row = rot_st; row<bg_row; ++row){
                ch = bwt[row] >>4;
                if(seq_ch==ch) ++top_num;
                ch = bwt[row] &0xF;
                if(seq_ch==ch) ++top_num; 
                if(n_ch == top_num) break; 
            }
        }
        if(bg_idx%2>0){
            ch = bwt[bg_row-1]&0xF;
            if(seq_ch == ch) {
                if(top_num >0)--top_num;
                ++mid_num;
            }
        }
        for(row = bg_row; row < ed_row; ++row){
            ch = bwt[row] >>4;
            if(seq_ch == ch) ++mid_num;
            ch = bwt[row]&0xF;
            if(seq_ch == ch) ++mid_num;
            if(n_ch == mid_num) break;
        }
        if(ed_idx%2==0){
            ch = bwt[ed_row]>>4;
            if(seq_ch==ch) ++mid_num;
            if(ed_idx < n_data-1){
                ch = bwt[ed_row]&0xF;
                if(seq_ch == ch) ++bot_num;
                ++ed_row;
            }
        }
        if(bg_idx > n_data-ed_idx){ 
            for(row=ed_row; row< rot_ed; ++row){
                ch = bwt[row]>>4;
                if(seq_ch==ch) ++bot_num; 
                ch = bwt[row]&0xF;
                if(seq_ch==ch) ++bot_num;
            }
            if(n_data%2>0 && ed_idx < n_data-1){
                ch = bwt[rot_ed]>>4;
                if(seq_ch==ch) {++bot_num; }
            }
        }
        
       
        if(bg_idx >n_data-ed_idx)  top_num = n_ch - bot_num-mid_num;
        if(err_len>0){
            ++cur_len;     
            if(cur_len >8) {
                //fprintf(stderr, "%s, %u, ERR: cur_len = %u\n", __func__, __LINE__, cur_len);
                exit(1);
            }
            if(mid_num < 1 && cur_len <= 7 ) { return -2;}
            if(mid_num >=1 && cur_len < 7){
                //----------------------------------------------
                if(cur_len <8-err_len){
                    bg_idx = top_num+cnt_2nt[rot][seq_ch];
                    ed_idx = bg_idx+mid_num-1;
                    continue;
                }
            }
        }
        if(err_len >0 && ((cur_len == 7 && mid_num >=1)||(cur_len ==8 && mid_num <=1)) ){
            uint32_t i_idx, j_idx, j_row, j_pos, j_rot;
            int err_pos = (pos+7)%8;
            if(mid_num > 0) {
                bg_idx = top_num+cnt_2nt[rot][seq_ch];
                ed_idx = bg_idx+mid_num-1;
            }
            uint8_t err_ch = (seq[err_pos*2]<<2) | (seq[err_pos*2+1]);  
            if(bg_idx > ed_idx) { 
                //fprintf(stderr, "%u, %s, Error bg_idx > ed_idx!!!!!\n", __LINE__, __func__);
                exit(1); 
            }
            for(i_idx = bg_idx; i_idx<=ed_idx; ++i_idx){
                j_idx = i_idx;
                j_rot = (rot+1)%8;
                rot_st =j_rot*((n_data+1)/2);
                j_row = rot_st+(j_idx+1)/2;
 
                if(j_idx%2>0) {
                    ch = bwt[j_row-1]&0xF; 
                } else{
                    ch = bwt[j_row]>>4;
                }
                n_err = 0; 
                //比较后两个bits
                if((ch & 3) != (err_ch &3)){ ++n_err; }
                //比较后3、4位置bit
                if((ch & 12) != (err_ch &12)){ ++n_err; }
                //仅保留单碱基匹配错误
                //if(i == 0 || i == 2) continue; 
                j_pos = pos;
                while(j_pos != 0){//[fixme?] rot != 0?
                    rot_st =j_rot*((n_data+1)/2);
                    j_row = rot_st+(j_idx+1)/2;
                    rot_ed = rot_st+n_data/2;
                    if(j_idx%2>0) {
                        seq_ch = bwt[j_row-1]&0xF; 
                    } else{
                        seq_ch = bwt[j_row]>>4;
                    }
                    n_ch = cnt_2nt[j_rot][seq_ch+1] - cnt_2nt[j_rot][seq_ch]; 
                    int num=0, k; 
                    if(j_idx <=(n_data+3)/4){
                        for(k = rot_st; k < j_row; ++k){
                            ch = bwt[k]>>4;
                            if(seq_ch == ch) ++num;
                            ch = bwt[k] &0xF;
                            if(seq_ch==ch) ++num;
                        }
                        if(j_idx%2==0){
                            ch = bwt[j_row]>>4;
                            if(seq_ch == ch) ++num;
                        }
                    }
                    
                    if(j_idx > (n_data+3)/4){
                        for(k = j_row; k < rot_ed; ++k){
                            ch = bwt[k]>>4;
                            if(seq_ch == ch){
                                if(k > j_row) ++num;
                                else if(j_idx%2>0) ++num;
                            }
                            ch = bwt[k]&0xF;
                            if(seq_ch==ch) ++num;
                        
                        } 
                        if(n_data%2>0){
                            if(j_idx < n_data -1){
                                ch = bwt[rot_ed]>>4;
                                if(seq_ch == ch) ++num;
                            }
                        } 
                        num = n_ch -num; 
                    }
                    j_idx = num + cnt_2nt[j_rot][seq_ch]-1;
                    --j_pos; 
                    j_rot  = (j_rot+1)%8; 
                }//end while(j_pos != 0)+++++++++++++++++++++++++ 
                int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = j_idx;
                    aln_out[n_aln][1] = penalty * n_err;
                    aln_out[n_aln][2] = err_pos;
                    aln_out[n_aln][3] = 0;                     
                    algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                    r_flag =2; 
                }
            }//end for+++++++++++++++++++++++++++
            aln_out[0][0] = n_aln;
            aln_out[0][1] = n_exact_aln;
            return r_flag;
        }
        if(mid_num >0){
            ++cur_len;
            if(cur_len <8){
                bg_idx = top_num+cnt_2nt[rot][seq_ch];
                ed_idx = bg_idx+mid_num-1;
                //rot = (rot+1)%8;

            } else if(cur_len == 8){ 
if(mid_num > 1) {
    //fprintf(stderr, "mid_num = %u, cur_len == 8", mid_num);
    //fprintf(stderr, "line = %u, func = %s\n", __LINE__, __func__);
    exit(1);
}
                bg_idx = top_num+cnt_2nt[rot][seq_ch];
                ed_idx = bg_idx+mid_num-1;
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    n_aln++; n_exact_aln++;
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[n_aln][3] = 0;  
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
                    r_flag = 1;
                }
                return r_flag;
            }
        } else if(mid_num == 0){
            if(cur_len == 7){
                uint32_t i_idx, j_idx, j_row, j_pos, j_rot;
                int err_pos = pos;
                uint8_t err_ch = (seq[err_pos*2]<<2) | (seq[err_pos*2+1]);  
                for(i_idx = bg_idx; i_idx<=ed_idx; ++i_idx){
                    j_idx = i_idx;
                    j_pos = (pos+1)%8;
                    j_rot = rot;
                    //++++++++++++++++++++++++++++++++++++++
                    rot_st =j_rot*((n_data+1)/2);
                    j_row = rot_st+(j_idx+1)/2;
     
                    if(j_idx%2>0) {
                        ch = bwt[j_row-1]&0xF; 
                    } else{
                        ch = bwt[j_row]>>4;
                    }
                    n_err = 0; 
                    //比较后两个bits
                    if((ch & 3) != (err_ch &3)){ ++n_err; }
                    //比较后3、4位置bit
                    if((ch & 12) != (err_ch &12)){ ++n_err; }
                    //仅保留单碱基匹配错误
                    //if(i == 0 || i == 2) continue; 
                    //--------------------------------------
                    while(j_pos != 0){//[fixme?] rot != 0?
                        rot_st =j_rot*((n_data+1)/2);
                        j_row = rot_st+(j_idx+1)/2;
                        rot_ed = rot_st+n_data/2;
                        if(j_idx%2>0) {
                            seq_ch = bwt[j_row-1]&0xF; 
                        } else{
                            seq_ch = bwt[j_row]>>4;
                        }
                        n_ch = cnt_2nt[j_rot][seq_ch+1] - cnt_2nt[j_rot][seq_ch]; 
                        int num=0, k; 
                        if(j_idx <=(n_data+3)/4){
                            for(k = rot_st; k < j_row; ++k){
                                ch = bwt[k]>>4;
                                if(seq_ch == ch) ++num;
                                ch = bwt[k] &0xF;
                                if(seq_ch==ch) ++num;
                            }
                            if(j_idx%2==0){
                                ch = bwt[j_row]>>4;
                                if(seq_ch == ch) ++num;
                            }
                        }
                        
                        if(j_idx > (n_data+3)/4){
                            for(k = j_row; k < rot_ed; ++k){
                                ch = bwt[k]>>4;
                                if(seq_ch == ch){
                                    if(k > j_row) ++num;
                                    else if(j_idx%2>0) ++num;
                                }
                                ch = bwt[k]&0xF;
                                if(seq_ch==ch) ++num;
                            
                            } 
                            if(n_data%2>0){
                                if(j_idx < n_data -1){
                                    ch = bwt[rot_ed]>>4;
                                    if(seq_ch == ch) ++num;
                                }
                            } 
                            num = n_ch -num; 
                        }
                        j_idx = num + cnt_2nt[j_rot][seq_ch]-1;
                        --j_pos; 
                        j_rot  = (j_rot+1)%8; 
                    }//end while(j_pos != 0)+++++++++++++++++++++++++ 
                    int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                    if(flag_aln == 0){
                        n_aln++; 
                        aln_out[n_aln][0] = j_idx;
                        aln_out[n_aln][1] = penalty*n_err;
                        aln_out[n_aln][2] = err_pos;
                        aln_out[n_aln][3] = 0;                         
                        algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                        r_flag = 2;
                    }
                }//end for+++++++++++++++++++++++++++
                aln_out[0][0] = n_aln;
                aln_out[0][1] = n_exact_aln;
                return r_flag;
            } else if(cur_len <7){
                aln_num++;
                if(aln_num >7) break;
                if(pos == 0) break;
                if(cur_len + st_len <7) break; 
                err_num = 0;
                st_len = cur_len+st_len-7;
                aln_in[(pos+1)%8] = 1;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
            } 
        }// if(mid_num == 0)+++++++++++++++++++++++++++++++++
    }//end while(1)++++++++++++++++++++++++++++++++++++++++++
   return r_flag; 
}
int align_nosum_indel(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t seq[18], int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4] , struct SubBuf *sub)
{
    int8_t seq_len = aln_in[10];
    int8_t seq_st  = aln_in[11];
    int8_t err_len = aln_in[12]; 
    int8_t InDel   = aln_in[13];
    int8_t penalty = aln_in[15]; 
    //gen_cnt(bwt, n_data, cnt_2nt);
    int stat_flag, mid_num, top_num, bot_num, bg_row, ed_row, row, rot_st, rot_ed;
    uint32_t bg_idx, ed_idx, j;//bg_idx,ed_idx begin and end of 2nt-encoded bwt 
    uint8_t seq_ch, ch;
    int n_aln = 0, n_exact_aln = 0;
    n_aln = aln_out[0][0];
    uint32_t pos = (seq_st+1)%8; 
    int rot = (7-seq_st)%8; //bw rotation
    int cur_len = 0; // current alignment length
    int st_len = seq_len;
    if(InDel == 0 || InDel < -2 || InDel > 2 ) {
        printf("%u, %s, InDel == %u , error!!\n", __LINE__, __func__, InDel);
        exit(1); 
    }
    if(n_data <= 1){
        printf("%u, %s, n_data == %u , error!!\n", __LINE__, __func__, n_data);
        exit(1); 
    } 
    int r_flag;
    if(err_len > 0) r_flag = -2;
    else r_flag = -1; 
    bg_idx = 0; ed_idx = n_data-1;//bwt interval [bg_idx, ed_idx]
    n_aln = aln_out[0][0];
    int aln_num = 0, err_num = 0;
    int edit_indel[4];
    int i;
    int seq_pos = 0;
    if(err_len ==0 ) {
        pos = 0;
        seq_pos = 0; 
    } else{
        seq_pos = pos*2; 
    }
    
    int len_seq = 16 + InDel;
    while(1){
        pos = (pos+7)%8; 
        seq_pos = (seq_pos+len_seq-2) % len_seq;
        seq_ch = seq[seq_pos]<<2| seq[seq_pos+1];
        rot = 7 - pos;
        if(bg_idx == 0 && ed_idx==n_data-1){
            if(cnt_2nt[rot][seq_ch] == cnt_2nt[rot][seq_ch+1]){
                if(err_num >=1) {break; }
                err_num++; 
                seq_pos = 2 * pos;
                if(cur_len > 0) st_len = cur_len+st_len-7;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
                continue;
            }
            bg_idx = cnt_2nt[rot][seq_ch]; 
            ed_idx = cnt_2nt[rot][seq_ch+1]-1;
            ++cur_len;
            if(cur_len == 8) {
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    ++n_aln;
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[n_aln][3] = InDel; 
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
                    r_flag = 1;
                }
                return r_flag;
            } 
            continue; 
        } //end else if(bg_idx == 0 && ed_idx==n_data-1)+++++++++++++++++
        mid_num =0; top_num = 0; bot_num = 0;        
        int n_ch = cnt_2nt[rot][seq_ch+1]-cnt_2nt[rot][seq_ch];
        rot_st = rot*((n_data+1)/2);
        bg_row = rot_st+(bg_idx+1)/2;
        ed_row = rot_st+(ed_idx+1)/2;
        rot_ed = rot_st+n_data/2;//if n_data%2==0, row rot_ed is not in current rotation bwt
        if(bg_idx <= n_data-ed_idx){
            stat_flag = 1;
            for(row = rot_st; row<bg_row; ++row){
                ch = bwt[row] >>4;
                if(seq_ch==ch) ++top_num;
                ch = bwt[row] &0xF;
                if(seq_ch==ch) ++top_num; 
                if(n_ch == top_num) break; 
            }
        }
        if(bg_idx%2>0){
            ch = bwt[bg_row-1]&0xF;
            if(seq_ch == ch) {
                if(top_num >0)--top_num;
                ++mid_num;
            }
        }
        for(row = bg_row; row < ed_row; ++row){
            ch = bwt[row] >>4;
            if(seq_ch == ch) ++mid_num;
            ch = bwt[row]&0xF;
            if(seq_ch == ch) ++mid_num;
            if(n_ch == mid_num) break;
        }
        if(ed_idx%2==0){
            ch = bwt[ed_row]>>4;
            if(seq_ch==ch) ++mid_num;
            if(ed_idx < n_data-1){
                ch = bwt[ed_row]&0xF;
                if(seq_ch == ch) ++bot_num;
                ++ed_row;
            }
        }
        if(bg_idx > n_data-ed_idx){ 
            for(row=ed_row; row< rot_ed; ++row){
                ch = bwt[row]>>4;
                if(seq_ch==ch) ++bot_num; 
                ch = bwt[row]&0xF;
                if(seq_ch==ch) ++bot_num;
            }
            if(n_data%2>0 && ed_idx < n_data-1){
                ch = bwt[rot_ed]>>4;
                if(seq_ch==ch) {++bot_num; }
            }
        }
        
        if(bg_idx >n_data-ed_idx)  top_num = n_ch - bot_num-mid_num;
        ++cur_len;      
        if(mid_num >=1 && cur_len < 7){
            bg_idx = top_num+cnt_2nt[rot][seq_ch];
            ed_idx = bg_idx+mid_num-1;
            continue;
        }
        if(err_len>0 && mid_num < 1 && cur_len <= 7 ) { 
            r_flag = -2;
            return r_flag;
        }
        r_flag = -2; 
        if((cur_len == 7 && mid_num >=1) || (cur_len == 8 && mid_num == 0)){
            uint32_t i_idx, j_idx, j_row, j_pos, j_rot;
            int err_pos = (pos+7)%8;
            bg_idx = top_num+cnt_2nt[rot][seq_ch];
            ed_idx = bg_idx+mid_num-1;
            uint8_t err_ch = (seq[err_pos*2]<<2) | (seq[err_pos*2+1]);  
            if(InDel == -1) seq_pos = (seq_pos+len_seq-1)%len_seq; //if(InDel == -2) seq_pos = seq_pos;
            else if(InDel == 1) seq_pos = (seq_pos+len_seq-3)%len_seq;
            else if(InDel == 2) seq_pos = (seq_pos+len_seq-4)%len_seq;
            for(i_idx = bg_idx; i_idx<=ed_idx; ++i_idx){
                j_idx = i_idx;
                j_rot = (rot+1)%8;
                rot_st =j_rot*((n_data+1)/2);
                j_row = rot_st+(j_idx+1)/2;
 
                if(j_idx%2>0) {
                    ch = bwt[j_row-1]&0xF; 
                } else{
                    ch = bwt[j_row]>>4;
                }
                edit_indel[0] = InDel;
                if(InDel < 0){ 
                    if(InDel == -2){
                        edit_indel[1] = seq_pos;
                        edit_indel[3] = ch;
                    } else if(InDel == -1){
                        if(seq[seq_pos] == (ch&3) ) {
                            edit_indel[1] = (seq_pos+len_seq-1)%len_seq;
                            edit_indel[3] = ((ch>>2)&3);
                        } else if(seq[seq_pos] == (ch>>2) ) {
                            edit_indel[1] = seq_pos;
                            edit_indel[3] = (ch&3);
                        } else{
                            continue; 
                        }
                    } 
                } 
                int edit_buf[4];
                if(InDel >0) {
                    if(InDel == 1){
                        edit_buf[0] = (seq[seq_pos]<<2)|(seq[seq_pos+1]);
                        edit_buf[1] = (seq[seq_pos]<<2)|(seq[seq_pos+2]);
                        edit_buf[2] = (seq[seq_pos+1]<<2)|(seq[seq_pos+2]);
                        for(i = 0; i < 3; ++i){
                            if(ch == edit_buf[i]) break; 
                        }
                        if(i == 3) {continue;} 
                        else{ 
                            if(i == 0) {                            
                                edit_indel[1] = (seq_pos+2)%len_seq;
                            } else if(i == 1) {                            
                                edit_indel[1] = (seq_pos+1)%len_seq;
                            } else if(i == 2) {                            
                                edit_indel[1] = (seq_pos+0)%len_seq;
                            }
                            edit_indel[3] = seq[edit_indel[1]]; 
                        }
                    } else if(InDel == 2){
                        edit_buf[0] = (seq[seq_pos]<<2)|(seq[seq_pos+1]);
                        edit_buf[1] = (seq[seq_pos]<<2)|(seq[seq_pos+3]);
                        edit_buf[2] = (seq[seq_pos+2]<<2)|(seq[seq_pos+3]);
                        for(i = 0; i < 3; ++i){
                            if(ch == edit_buf[i]) break; 
                        }

                        if(i == 3) {continue;} 
                        else{
                            if(i == 0) {                            
                                edit_indel[1] = (seq_pos+2)%len_seq;
                            } else if(i == 1) {                            
                                edit_indel[1] = (seq_pos+1)%len_seq;
                            } else if(i == 2) {                            
                                edit_indel[1] = (seq_pos+0)%len_seq;
                            }
                            edit_indel[3] = (seq[edit_indel[1]]<<2) | (seq[edit_indel[1]+1]);
                        }
                    } 
                }
                j_pos = pos;
                while(j_pos != 0){//[fixme?] rot != 0?
                    rot_st =j_rot*((n_data+1)/2);
                    j_row = rot_st+(j_idx+1)/2;
                    rot_ed = rot_st+n_data/2;
                    if(j_idx%2>0) {
                        seq_ch = bwt[j_row-1]&0xF; 
                    } else{
                        seq_ch = bwt[j_row]>>4;
                    }
                    n_ch = cnt_2nt[j_rot][seq_ch+1] - cnt_2nt[j_rot][seq_ch]; 
                    int num=0, k; 
                    if(j_idx <=(n_data+3)/4){
                        for(k = rot_st; k < j_row; ++k){
                            ch = bwt[k]>>4;
                            if(seq_ch == ch) ++num;
                            ch = bwt[k] &0xF;
                            if(seq_ch==ch) ++num;
                        }
                        if(j_idx%2==0){
                            ch = bwt[j_row]>>4;
                            if(seq_ch == ch) ++num;
                        }
                    }
                    if(j_idx > (n_data+3)/4){
                        for(k = j_row; k < rot_ed; ++k){
                            ch = bwt[k]>>4;
                            if(seq_ch == ch){
                                if(k > j_row) ++num;
                                else if(j_idx%2>0) ++num;
                            }
                            ch = bwt[k]&0xF;
                            if(seq_ch==ch) ++num;
                        } 
                        if(n_data%2>0){
                            if(j_idx < n_data -1){
                                ch = bwt[rot_ed]>>4;
                                if(seq_ch == ch) ++num;
                            }
                        } 
                        num = n_ch -num; 
                    }
                    j_idx = num + cnt_2nt[j_rot][seq_ch]-1;
                    --j_pos; 
                    j_rot  = (j_rot+1)%8;
                    
                }//end while(pos != 0)+++++++++++++++++++++++++ 
                int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = j_idx;
                    aln_out[n_aln][1] = penalty;
                    aln_out[n_aln][2] = err_pos;
                    aln_out[n_aln][3] = InDel; 
                    algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                    if(err_len > 0 ) r_flag = 2;
                    else r_flag = 1;
                }
            }//end for+++++++++++++++++++++++++++
            if(aln_out[0][0] == n_aln && err_len == 0) {
                mid_num = 0;
                cur_len--;
            }else{
                aln_out[0][0] = n_aln;
                aln_out[0][1] = n_exact_aln;
                break; 
            }
        }
        if(err_len == 0 && mid_num== 0 && cur_len <=7){
            aln_num++;
            if(aln_num >7) break;
            if(cur_len + st_len <7) break; 
            if(pos == 0) {
                if(cur_len == 7 ) pos++;
                else {
                    pos = 7;
                    aln_num = 7; 
                }
            }
            seq_pos = pos*2;
            st_len = cur_len+st_len-7;
            aln_in[(pos+1)%8] = 1;
            err_num = 0;
            cur_len = 0;
            bg_idx = 0;
            ed_idx = n_data-1;
            continue;
        } 
    }//end while(1)++++++++++++++++++++++++++++++++++++++++++
    return r_flag; 
}

/*
//fprintf(stderr, "\n========================================================\n");
//fprintf(stderr, "(===bg_idx, ed_idx, n_data) = (%u, %u, %u)\n", bg_idx, ed_idx, n_data);
//fprintf(stderr, "(===cur_len, cur_rot) = (%u, %u)\n", cur_len, cur_rot);
int __n0 = __dna2_count(Bwt+rot_off, 0, bg_idx, seq_ch);
int __tot = __dna2_count_small(Bwt+rot_off, 0, n_data, seq_ch);
int __n1 = __dna2_count(Bwt+rot_off, bg_idx, ed_idx+1, seq_ch);
int __n2 = __dna2_count(Bwt+rot_off, ed_idx+1, n_data, seq_ch);
//fprintf(stderr, "seq_ch, tot, n0, n1, n2= %u %u %u %u %u\n", seq_ch, __tot, __n0, __n1, __n2);
//fprintf(stderr, "bg_idx, ed_idx= %u %u\n", __n0+__tot, __n0+__tot+__n1-1);
//fprintf(stderr, "\n--------------------------------------------------------\n");
*/  
/*
{
uint8_t *seq_buf = calloc(100000, 1);
test_AlgnBwtSml_large(Bwt, n_data, seq_buf);
uint32_t __j;
for(__j= 0; __j<n_data; __j++){
    uint8_t *cur_seq = seq_buf+16*__j;
    //fprintf(stderr, "%2u ", __j);
    for(i = 0; i < 8; ++i){
        //fprintf(stderr, "%3u ", (cur_seq[i*2]<<2)|(cur_seq[i*2+1]));
    }
    //fprintf(stderr, "\n");

}
for(i = 0; i < 8; ++i) //fprintf(stderr, "%u\t", (seq[i*2]<<2)|(seq[i*2+1]));
//fprintf(stderr, "\n");

free(seq_buf);
}
*/
        
int align_large(uint8_t *Bwt, uint8_t *cnt2, uint32_t n_data, const uint8_t *seq , int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4] , struct SubBuf *sub ){

    int i, j, k;
    int8_t seq_len = aln_in[10]; 
    int8_t seq_st = aln_in[11]; 
    int8_t err_len = aln_in[12]; 
    int8_t penalty = aln_in[15]; 
    int8_t n_err; 
    
    int r_flag;
    if(err_len >0) r_flag = -2;
    else r_flag =-1;

    int n_aln=0, n_exact_aln=0;
    n_aln = aln_out[0][0]; 
 
    uint8_t  *pBwt, *pSum ;
    uint32_t sum_buf;
    uint32_t idx_i[2], idx_q[2] , idx_r[2] , idx_row[2], rot_off;
    uint32_t upp_num2[2];//upp_num[0] = number of ch_seq in [0, bg_idx);
                         //upp_num[1] = number of ch_seq in [bg_idx, ed_idx];
    
    pBwt =Bwt ;
    pSum = cnt2 ;
    int bgn_pos = 0;        
   
    int bg_idx = 0, ed_idx = n_data-1;
    int alg_num = 0, bgn_rot = 0, siz_cnt, bgn_stat;
    int pos = (seq_st+1)%8;
    int rot = (7-seq_st)%8;   
    int cur_len = 0, st_len = seq_len, err_num = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    if( n_data <= 254*256 ) siz_cnt = 2 ;
    if( n_data > 254*256 ) siz_cnt = 4 ;
    int rng_num = ((n_data+253) /254) ;
    int len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;
    int rot_sum =rot*16*len_sum;
    int data_q = (n_data)/254;
    int data_r = (n_data)%254;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // bg_idx 和ed_idx是下一个轮换内部的相对位置+++++
    uint32_t upp_num=0, mid_num=0, seq_num[2];
    uint32_t rng_cnt[2], rng_end[2], ch_cnt[2], rng_len, fst_row, bgn_row, end_row, lst_row, row;
    int ch_off, idx_num, seq_ch;
    uint8_t get_ch, ch;

    while(1){
        pos = (pos+7)%8;
        seq_ch = (seq[2*pos]<<2) | seq[2*pos+1] ;
        rot = 7 -pos; 
        rot_sum =rot*16*len_sum;       
        rot_off = rot*((n_data+1)/2);
/*   
//fprintf(stderr, "\n========================================================\n");
//fprintf(stderr, "(===bg_idx, ed_idx, n_data) = (%u, %u, %u)\n", bg_idx, ed_idx, n_data);
//fprintf(stderr, "(===cur_len, cur_rot) = (%u, %u)\n", cur_len, rot);
int __n0 = __dna2_count(Bwt+rot_off, 0, bg_idx, seq_ch);
int __tot = __dna2_count_small(Bwt+rot_off, 0, n_data, seq_ch);
int __n1 = __dna2_count(Bwt+rot_off, bg_idx, ed_idx+1, seq_ch);
int __n2 = __dna2_count(Bwt+rot_off, ed_idx+1, n_data, seq_ch);
//fprintf(stderr, "seq_ch, tot, n0, n1, n2= %u %u %u %u %u\n", seq_ch, __tot, __n0, __n1, __n2);
//fprintf(stderr, "bg_idx, ed_idx= %u %u\n", __n0+__tot, __n0+__tot+__n1-1);
//fprintf(stderr, "\n--------------------------------------------------------\n");
*/


        if((bg_idx== 0)&&( ed_idx == n_data-1)) {
            mid_num = 0; upp_num =0; 
            if(seq_ch ==0){
                ch_off = rot_sum + (len_sum -siz_cnt);
                for(i=0; i<siz_cnt; i++){
                    mid_num <<=8 ;
                    mid_num |= pSum[ch_off+i];
                }
                bg_idx = 0;
                ed_idx = mid_num - 1;
            }else{
                ch_off = rot_sum +seq_ch*len_sum + (len_sum -siz_cnt) ;
                for(i=0; i<siz_cnt; i++){
                    upp_num <<=8 ;
                    upp_num |= pSum[ch_off -len_sum+i];
                    mid_num <<=8;
                    mid_num |= pSum[ch_off +i];                   
                }
                bg_idx = upp_num ;
                ed_idx = mid_num - 1;
            }
            if((mid_num == 0) || (bg_idx> ed_idx)){
                if(err_num >=1) {break;}
                err_num++; 
                st_len = (cur_len+1)%8;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
                continue;
            }
            cur_len = (cur_len+1)%8;  
            if(cur_len == 8) {
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[n_aln][3] = 0; 
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
fprintf(stderr, "%u\n", __LINE__);
                }
fprintf(stderr, "%u\n", __LINE__);
                return 1;
            } 
            continue;
        }
        mid_num = 0; upp_num = 0;
        rot_off = rot*((n_data+1)/2);
        ch_off = 16*rot*len_sum + len_sum*seq_ch ;
        idx_i[0] = (bg_idx>0) ? (bg_idx-1) : 0;
        idx_i[1] = ed_idx;
        for(i=0; i<2; i++){
            idx_q[i] = idx_i[i]/254;
            idx_r[i] = idx_i[i]%254;
            idx_row[i] = rot_off +(idx_i[i]+1)/2;
        }
        idx_num=2;
        uint32_t glb_cnt[2];
        ch_cnt[0] = 0; ch_cnt[1] = 0;
        for(i=0; i< idx_num ; i++){
            if(bg_idx == 0 && 0 == i){
                continue;
            }
            rng_cnt[i] = (8+siz_cnt)*(idx_q[i] / 8);
            rng_end[i] = rng_cnt[i] + idx_q[i] % 8;
            if(rng_cnt[i]>0){
                for(j = 0; j< siz_cnt; j++){
                    sum_buf = (uint32_t)pSum[ch_off+rng_cnt[i] - siz_cnt +j];
                    ch_cnt[i] = ((ch_cnt[i]<<8)|sum_buf); 
                }
            }
            for(k= rng_cnt[i]; k< rng_end[i]; k++){
                sum_buf = (uint32_t)pSum[ch_off + k];
                ch_cnt[i] += sum_buf;
            }      
            seq_num[i]  = (uint32_t)pSum[ch_off + rng_end[i]];
            glb_cnt[i] = 0 ;
            if(seq_ch > 0){
                for(j = 0; j<siz_cnt; j++) {
                    sum_buf =  (uint32_t)pSum[ch_off -siz_cnt+j]; 
                    glb_cnt[i] = ((glb_cnt[i]<<8)|sum_buf);
                }
            }
            ch_cnt[i] += glb_cnt[i];
        }// end for
        if(bg_idx==0) ch_cnt[0] += glb_cnt[1];
        upp_num2[0] = 0;upp_num2[1] = 0;              
        if(idx_q[0] < idx_q[1] ){ 
            for(i = 0 ;  i<2 ; i++ ){
                if(bg_idx == 0 && 0 ==i){ 
                    upp_num2[0] += ch_cnt[0];
                    continue;
                }
                if(idx_r[i] <= 134){
                    for( k = rot_off+(idx_i[i]-idx_r[i])/2;  k < idx_row[i];  k++){
                        get_ch = pBwt[k];
                        ch = get_ch>>4;
                        if(seq_ch ==ch) upp_num2[i]++;
                        ch = get_ch &0xF;
                        if(seq_ch == ch) upp_num2[i]++;
                    }
                    if(idx_r[i] % 2 == 0 ){// 该点已经计算过，但需要考虑bg_idx点是否包含在其中++++ 
                        get_ch = pBwt[idx_row[i]];
                        ch = get_ch>>4;
                        if(seq_ch == ch) ++upp_num2[i];
                    }
                }
                if(idx_r[i] > 134){  //离下面端点更近的情况运行代码+++
                    upp_num2[i]=(uint32_t)pSum[ch_off +rng_end[i]];
                    rng_len = (idx_q[i] == data_q) ? data_r : 254;
                    rng_end[i] = rot_off+ (idx_i[i] - idx_r[i])/2 + rng_len /2;
                    for(k= idx_row[i]; k< rng_end[i]; k++){
                        get_ch = pBwt[k] ;
                        if(k > idx_row[i] || idx_r[i]%2 >0 ){
                            ch = get_ch >> 4;
                            if(seq_ch ==ch)  --upp_num2[i];
                        }
                        ch = get_ch & 0xF ;
                        if(seq_ch ==ch)  --upp_num2[i];
                    }
                    if((rng_len %2 >0)&&(idx_r[i] < rng_len-1)){
                        get_ch = pBwt[rng_end[i]];
                        ch = get_ch >>4;
                        if(seq_ch == ch){ 
                            --upp_num2[i];
                        } 
                    }
                }
                upp_num2[i] += ch_cnt[i] ;
            } // End: ++++++++++++++ for(i = 0 ;  i<2 ; i++ )
        }// End : +++++++ if(bg_idx_q < ed_idx_q )+++++++++++++++++++
        // 以下是bg_idx和ed_idx在同一个间隔内++++++++++++++++++++++++
        uint32_t bot_num = 0, mid_num = 0;
        if(idx_q[0] == idx_q[1] ) {
            rng_len = (data_q == idx_q[1]) ? data_r : 254;  
            bgn_row = rot_off+(idx_i[0]+1)/2;
            end_row = rot_off+(idx_i[1]+1)/2;
            fst_row = rot_off +(idx_i[0]-idx_r[0])/2;
            lst_row = fst_row + rng_len/2;
            if( (idx_r[0] <= rng_len - idx_r[1]) && bg_idx >0 ){  
                //上段距离小于等于下段距离++++++++
                for(row= fst_row ; row < bgn_row; row++){
                    get_ch = pBwt[row];
                    ch = get_ch>>4;
                    if(seq_ch ==ch) upp_num2[0]++;
                    ch = get_ch &0xF;
                    if(seq_ch ==ch) upp_num2[0]++;
                }
            }
            //+++++++++++++++++++++++++++++++++++++++++++
            if(bg_idx%2 >0){
                get_ch = pBwt[bgn_row]; 
                ch = get_ch>>4;
                if(seq_ch ==ch) upp_num2[0]++;
                ch = get_ch&0xF;
                if(seq_ch==ch)  upp_num2[1]++;
                bgn_row++;
            } else if(ed_idx == 0){
                get_ch = pBwt[bgn_row]; 
                ch = get_ch>>4;
                if(seq_ch==ch)  upp_num2[1]++;
                bgn_row++;
            }
            for(row= bgn_row ; row < end_row ; row++){
                get_ch = pBwt[row];
                ch = get_ch>>4;
                if(seq_ch ==ch) upp_num2[1]++;
                ch = get_ch &0xF;
                if(seq_ch ==ch) upp_num2[1]++;
            }
            if(ed_idx%2==0 ) { // ed_idx端边界处理++++++++++++++++++++++++
                get_ch = pBwt[end_row] ;              
                ch = get_ch >>4;
                if(seq_ch ==ch && ed_idx >0) upp_num2[1]++;
                if(idx_r[1] < rng_len-1){  // idx_r[1]是ed_idx_r
                    ch = get_ch & 0xF;
                    if(seq_ch ==ch) {
                        bot_num++;
                    }
                    end_row++;
                } 
            }
            if( idx_r[0] > rng_len - idx_r[1] ) { //上段距离大于等于下段距离++++++++
                for(row= end_row; row < lst_row; row++){
                    get_ch = pBwt[row]; 
                    ch = get_ch>>4;        
                    if(seq_ch ==ch){ 
                        bot_num++;
                    }
                    ch = get_ch &0xF;    
                    if(seq_ch ==ch) {
                        bot_num++;
                    }
                }
                if(rng_len%2 > 0 && idx_r[1] < rng_len-1){ // 最末端边界+++
                    get_ch = (*(pBwt + lst_row));              
                    ch = get_ch >>4;   
                    if(seq_ch ==ch) bot_num++;
                }
                // 计算bg_idx点开始第一个出现seq_ch的为止，所有seq_ch出现的个数++
                upp_num2[0] = seq_num[1] - bot_num - upp_num2[1];
            } // End: +++++++++++++if( idx_r[0] > rng_len – idx_r[1] )
            upp_num2[0] += ch_cnt[1];
            upp_num2[1] += upp_num2[0];
        }  //  End: +++++++++++++ if(idx_q[0] == idx_q[1] )
        mid_num = upp_num2[1] - upp_num2[0];
        int alg_err_num = 0;
        if(err_len>0){
            ++cur_len;
            if(mid_num < 1 && cur_len <= 7 ) { 
                return -2;}
            if(mid_num >=1 && cur_len < 7){
                if(cur_len <8-err_len){
                    bg_idx =  upp_num2[0]; 
                    ed_idx =  upp_num2[1]-1;
                    continue; 
                }
            }
        }
        if(err_len >0 && cur_len == 7 && mid_num >=1) {
            uint32_t i_idx, j_idx, j_rot, j_pos;
            int err_pos= (pos+7)%8;
            bg_idx =  upp_num2[0]; 
            ed_idx =  upp_num2[1]-1;
            uint8_t err_ch = (seq[err_pos*2]<<2) | (seq[err_pos*2+1]);
            for(i_idx = bg_idx ;  i_idx <= ed_idx ;  i_idx++){
                j_idx = i_idx;    
                j_rot = (rot+1)%8;
                rot_off =j_rot*((n_data+1)/2);
                uint32_t j_row = rot_off+(j_idx+1)/2;

                if(j_idx%2>0) {
                    ch = pBwt[j_row-1]&0xF; 
                } else{
                    ch = pBwt[j_row]>>4;
                }
                n_err = 0; 
                //比较后两个bits
                if((ch & 3) != (err_ch &3)){ ++n_err; }
                //比较后3、4位置bit
                if((ch & 12) != (err_ch &12)){ ++n_err; }
                //仅保留单碱基匹配错误
                j_pos = pos;
                // 查找某一近似比对成功序列的行号 ++++++++++++++++++
                /*计算编辑距离*/
                while(j_pos != 0){
                    rot_off = j_rot*((n_data+1)/2);
                    uint32_t j_row = rot_off +( j_idx+1)/2;
                    lst_row = rot_off  + n_data/2 ;
                    if(j_idx%2>0)  {
                        get_ch =pBwt[j_row-1];
                        seq_ch =  get_ch&0xF;
                    } else if(j_idx%2==0)  {
                        get_ch =pBwt[j_row];
                        seq_ch  =  get_ch>>4;
                    }
                    ch_off = 16*j_rot*len_sum + len_sum*seq_ch ;
                    idx_i[0] = j_idx;
                    idx_q[0] = idx_i[0]/254;
                    idx_r[0] = idx_i[0]%254;
                    idx_row[0] = rot_off +( idx_i[0]+1)/2;
                    rng_cnt[0] = (8+siz_cnt)*(idx_q[0] / 8);
                    rng_end[0] = rng_cnt[0] + idx_q[0] % 8;
                    mid_num = 0;
                    upp_num = 0;
                    bot_num = 0;
                    ch_cnt[0] = 0;
                    if(rng_cnt[0]>0){
                        for(j = 0; j< siz_cnt; j++){
                            ch_cnt[0]=(ch_cnt[0]<<8)|(uint32_t)(*(pSum+ch_off+rng_cnt[0]- siz_cnt +j)); 
                        }
                    }
                    for(k= rng_cnt[0]; k< rng_end[0]; k++){ 
                        ch_cnt[0] += (uint32_t)(*(pSum+ ch_off + k ));
                    }       
                    seq_num[0]  = (uint32_t)(*(pSum+ ch_off + rng_end[0] ));
                    glb_cnt[0] = 0 ;
                    if(seq_ch > 0){
                        for(j = 0; j<siz_cnt; j++) 
                            glb_cnt[0] = (glb_cnt[0]<<8) | (uint32_t)(*(pSum+ch_off - siz_cnt+j)) ; 
                    }
                    ch_cnt[0] += glb_cnt[0];
                    bot_num = 0 ;
                    bot_num=(uint32_t)pSum[ch_off +rng_end[0]];
                    rng_len = 254 ;
                    if(idx_q[0]==data_q) rng_len = data_r;
                    rng_end[0] = rot_off+ (idx_i[0] - idx_r[0])/2 + rng_len/2;
                    for(k= idx_row[0]; k< rng_end[0] ; k++){
                        get_ch = pBwt[k] ;
                        ch = get_ch>>4;
                        if(seq_ch ==ch) { 
                            if(k > idx_row[0] || idx_r[0]%2 >0 ){--bot_num;}
                        }
                        ch = get_ch&0xF ;
                        if(seq_ch ==ch)  --bot_num;
                    }
                    if((rng_len %2 >0 ) &&(rng_len-1 > idx_r[0])){  
                        get_ch = pBwt[rng_end[0]];
                        ch = get_ch>>4;
                        if(seq_ch == ch) --bot_num; 
                    }
                    bot_num += ch_cnt[0];
                    j_idx = bot_num-1;
                    j_pos = (j_pos +7)%8;
                    j_rot = (j_rot +1)%8;
                }// End : while(j_pos != 0) ----------------------------------
                int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                if(flag_aln == 0){
                    n_aln++;
                    aln_out[n_aln][0] = j_idx;
                    aln_out[n_aln][1] = penalty*n_err;
                    aln_out[n_aln][2] = err_pos;
                    aln_out[n_aln][3] = 0;  
                    algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                    r_flag = 2;

                }
            }// for(i_idx = bg_idx;i_idx <= ed_idx; i_idx ++)--------------------------
            aln_out[0][0] = n_aln; 
            aln_out[0][1] = n_exact_aln; 
fprintf(stderr, "%u\n", __LINE__);
            return r_flag;
        }
        if( mid_num>0) { // 比对尚未结束++++++++++
            cur_len += 1;  
            if( cur_len < 8){
                bg_idx =  upp_num2[0]; 
                ed_idx =  upp_num2[1]-1;
                continue;
            }
            if( cur_len == 8){
                bg_idx =  upp_num2[0]; 
                ed_idx = upp_num2[1]-1;
if(mid_num > 1) {
    //fprintf(stderr, "mid_num = %u, cur_len == 8\n", mid_num);
    //fprintf(stderr, "%u, func = %s\n", __LINE__, __func__);
    printf("%u, func = %s\n", __LINE__, __func__);
    exit(1);
}
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    alg_err_num  = 0 ;
                    n_aln++; n_exact_aln++;
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[n_aln][3] = 0; 
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] = n_exact_aln;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
fprintf(stderr, "%u\n", __LINE__);
                }
                return 1;     
            }
        } // End : if(mid_num>0) // 比对尚未结束------
        if(mid_num==0){//本次 比对结束 +++++++++
            if(cur_len == 7){ //输出近似比对成功结果+++++++++++++++
                uint32_t i_idx, j_idx, j_rot, j_pos;
                int err_pos= pos;
                uint8_t err_ch = (seq[err_pos*2]<<2) | (seq[err_pos*2+1]);
                for(i_idx = bg_idx ;  i_idx <= ed_idx ;  i_idx++){
                    uint32_t j_row;
                    j_idx = i_idx;    
                    j_pos = (pos+1)%8;
                    j_rot= rot;
                    
                    //++++++++++++++++++++++++++++++++ 
                    rot_off =j_rot*((n_data+1)/2);
                    j_row = rot_off+(j_idx+1)/2;
                    if(j_idx%2>0) {
                        ch = pBwt[j_row-1]&0xF; 
                    } else{
                        ch = pBwt[j_row]>>4;
                    }
                    n_err = 0; 
                    //比较后两个bits
                    if((ch & 3) != (err_ch &3)){ ++n_err; }
                    //比较后3、4位置bit
                    if((ch & 12) != (err_ch &12)){ ++n_err; }
                    //仅保留单碱基匹配错误
                    //---------------------------------- 
                    /*计算编辑距离*/
                    while(j_pos != 0){
                        rot_off = j_rot*((n_data+1)/2);
                        j_row = rot_off +( j_idx+1)/2;
                        lst_row = rot_off  + n_data/2 ;
                        if(j_idx%2>0)  {
                            get_ch =pBwt[j_row-1];
                            seq_ch =  get_ch&0xF;
                        } else if(j_idx%2==0)  {
                            get_ch =pBwt[j_row];
                            seq_ch  =  get_ch>>4;
                        }
                        ch_off = 16*j_rot*len_sum + len_sum*seq_ch ;
                        idx_i[0] = j_idx;
                        idx_q[0] = idx_i[0]/254;
                        idx_r[0] = idx_i[0]%254;
                        idx_row[0] = rot_off +( idx_i[0]+1)/2;
                        rng_cnt[0] = (8+siz_cnt)*(idx_q[0] / 8);
                        rng_end[0] = rng_cnt[0] + idx_q[0] % 8;
                        mid_num = 0;
                        upp_num = 0;
                        bot_num = 0;
                        ch_cnt[0] = 0;
                        if(rng_cnt[0]>0){
                            for(j = 0; j< siz_cnt; j++){
                                ch_cnt[0]=(ch_cnt[0]<<8)|(uint32_t)(*(pSum+ch_off+rng_cnt[0]- siz_cnt +j)); 
                            }
                        }
                        for(k= rng_cnt[0]; k< rng_end[0]; k++){ 
                            ch_cnt[0] += (uint32_t)(*(pSum+ ch_off + k ));
                        }       
                        seq_num[0]  = (uint32_t)(*(pSum+ ch_off + rng_end[0] ));
                        glb_cnt[0] = 0 ;
                        if(seq_ch > 0){
                            for(j = 0; j<siz_cnt; j++) 
                                glb_cnt[0] = (glb_cnt[0]<<8) | (uint32_t)(*(pSum+ch_off - siz_cnt+j)) ; 
                        }
                        ch_cnt[0] += glb_cnt[0];
                        bot_num = 0 ;
                        bot_num=(uint32_t)pSum[ch_off +rng_end[0]];
                        rng_len = 254 ;
                        if(idx_q[0]==data_q) rng_len = data_r;
                        rng_end[0] = rot_off+ (idx_i[0] - idx_r[0])/2 + rng_len/2;
                        for(k= idx_row[0]; k< rng_end[0] ; k++){
                            get_ch = pBwt[k] ;
                            ch = get_ch>>4;
                            if(seq_ch ==ch) { 
                                if(k > idx_row[0] || idx_r[0]%2 >0 ){--bot_num;}
                            }
                            ch = get_ch&0xF ;
                            if(seq_ch ==ch)  --bot_num;
                        }
                        if((rng_len %2 >0 ) &&(rng_len-1 > idx_r[0])){  
                            get_ch = pBwt[rng_end[0]];
                            ch = get_ch>>4;
                            if(seq_ch == ch) --bot_num; 
                        }
                        bot_num += ch_cnt[0];
                        j_idx = bot_num-1;
                        j_pos = (j_pos +7)%8;
                        j_rot = (j_rot +1)%8;
                    }// End : while(j_pos != 0) ----------------------------------
                    int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                    if(flag_aln == 0){
                        n_aln++;
                        aln_out[n_aln][0] = j_idx;
                        aln_out[n_aln][1] = penalty * n_err;
                        aln_out[n_aln][2] = err_pos;
                        aln_out[n_aln][3] = 0;                         
                        algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                        r_flag = 2;
                        
fprintf(stderr, "%u\n", __LINE__);
                    }
                }// for(i_idx = bg_idx;i_idx <= ed_idx; i_idx ++)--------------------------
                aln_out[0][0] = n_aln; 
                aln_out[0][1] = n_exact_aln; 
fprintf(stderr, "%u\n", __LINE__);
                return r_flag;
            } //  if(cur_len == 7)  //近似比对成功结果处理结束------------------
            if(cur_len < 7){ // 重新选定比对起始位置进行比对 ++++++++++++++++
                alg_num++;
                if(alg_num >7) {break;}
                if(pos ==0) {break;}
                if(cur_len +st_len <7) { break;} 
                err_num = 0;
                st_len = cur_len+st_len-7;
                aln_in[(pos+1)%8] =1;
                cur_len = 0 ;
                bg_idx = 0; 
                ed_idx = n_data-1;
            }  // End:  ++++++++++++++++++if(cur_len < 7) 
        }// End : +++++++++++++++++++if((mid_num==0)
    }//++++++++++++++++++++++++++++ while(1)
fprintf(stderr, "%u\n", __LINE__);
    return r_flag;
} // End : Alignm_big( )------------------------------------------------------

int align_large_indel(uint8_t *Bwt, uint8_t *cnt2, uint32_t n_data, const uint8_t *seq , int8_t *aln_in, uint8_t *algn_row, uint32_t (*aln_out)[4] , struct SubBuf *sub){
    
    int8_t seq_len = aln_in[10];
    int8_t seq_st  = aln_in[11]; 
    int8_t err_len = aln_in[12];
    int8_t InDel   = aln_in[13];
    int8_t penalty = aln_in[15]; 
    int n_aln=0, n_exact_aln=0;
    n_aln = aln_out[0][0]; 
 
    uint8_t  *pBwt, *pSum ;
    uint32_t sum_buf;
    uint32_t idx_i[2], idx_q[2] , idx_r[2] , idx_row[2], rot_off;
    uint32_t upp_num2[2];//upp_num[0] = number of ch_seq in [0, bg_idx);
                         //upp_num[1] = number of ch_seq in [bg_idx, ed_idx];
    int i, j, k;
    pBwt =Bwt ;
    pSum = cnt2 ;
    int bgn_pos = 0;        
    if(InDel == 0 || InDel < -2 || InDel > 2 ) {
        printf("%u, %s, InDel == %u , error!!\n", __LINE__, __func__, InDel);
        exit(1); 
    }
    if(n_data <= 1){
        printf("%u, %s, n_data == %u , error!!\n", __LINE__, __func__, n_data);
        exit(1); 
    }
    int r_flag;
    if(err_len > 0) r_flag = -2;
    else r_flag = -1; 

    int bg_idx = 0, ed_idx = n_data-1;
    int alg_num = 0, bgn_rot = 0, siz_cnt, bgn_stat;
    int pos = (seq_st+1)%8;
    int rot = (7-seq_st)%8;   
    int cur_len = 0, st_len = seq_len, err_num = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    if( n_data <= 254*256 ) siz_cnt = 2 ;
    if( n_data > 254*256 ) siz_cnt = 4 ;
    int rng_num = ((n_data+253) /254) ;
    int len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;
    int rot_sum =rot*16*len_sum  ;
    int data_q = (n_data)/254;
    int data_r = (n_data)%254;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // bg_idx 和ed_idx是下一个轮换内部的相对位置+++++
    int edit_indel[4];
    uint32_t upp_num=0, mid_num=0, seq_num[2];
    uint32_t rng_cnt[2], rng_end[2], ch_cnt[2], rng_len, fst_row, bgn_row, end_row, lst_row, row;
    int seq_ch;   
    
    int seq_pos = 0;
    if(err_len ==0 ) {
        pos = 0;
        seq_pos = 0; 
    } else{
        seq_pos = pos*2; 
    }
    int len_seq = 16 + InDel;
    
    while(1){
        pos = (pos+7)%8; 
        seq_pos = (seq_pos+len_seq-2) % len_seq;
        seq_ch = seq[seq_pos]<<2| seq[seq_pos+1];
        rot = 7-pos;

        int rot_sum =rot*16*len_sum  ;
        int ch_off, idx_num;
        uint8_t get_ch, ch;
        rot_off = rot*((n_data+1)/2);

        if((bg_idx== 0)&&( ed_idx == n_data-1)) {
            mid_num = 0; upp_num =0; 
            if(seq_ch ==0){
                ch_off = rot_sum + (len_sum -siz_cnt);
                for(i=0; i<siz_cnt; i++){
                    mid_num <<=8 ;
                    mid_num |= pSum[ch_off+i];
                }
                bg_idx = 0;
                ed_idx = mid_num - 1;
            }else{
                ch_off = rot_sum +seq_ch*len_sum + (len_sum -siz_cnt) ;
                for(i=0; i<siz_cnt; i++){
                    upp_num <<=8 ;
                    upp_num |= pSum[ch_off -len_sum+i];
                    mid_num <<=8;
                    mid_num |= pSum[ch_off +i];                   
                }
                bg_idx = upp_num ;
                ed_idx = mid_num - 1;
            }

           
           if((mid_num == 0) || (bg_idx> ed_idx)){
                if(err_num >=1) {break;}
                err_num++; 
                //st_len = (cur_len+1)%8;
                if(cur_len > 0) st_len = cur_len+st_len-7;
                seq_pos = 2 * pos;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
                continue;
            }
            cur_len = (cur_len+1)%8;  
            if(cur_len == 8) {
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[n_aln][3] = InDel;                      
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
                    r_flag = 1;
                }
                return r_flag;
            } 
            continue;
        }
        mid_num = 0; upp_num = 0;
        rot_off = rot*((n_data+1)/2);
        ch_off = 16*rot*len_sum + len_sum*seq_ch ;
        idx_i[0] = (bg_idx>0) ? (bg_idx-1) : 0;
        idx_i[1] = ed_idx;
        for(i=0; i<2; i++){
            idx_q[i] = idx_i[i]/254;
            idx_r[i] = idx_i[i]%254;
            idx_row[i] = rot_off +(idx_i[i]+1)/2;
        }
        idx_num=2;
        uint32_t glb_cnt[2];
        ch_cnt[0] = 0; ch_cnt[1] = 0;
        for(i=0; i< idx_num ; i++){
            if(bg_idx == 0 && 0 == i){
                continue;
            }
            rng_cnt[i] = (8+siz_cnt)*(idx_q[i] / 8);
            rng_end[i] = rng_cnt[i] + idx_q[i] % 8;
            if(rng_cnt[i]>0){
                for(j = 0; j< siz_cnt; j++){
                    sum_buf = (uint32_t)pSum[ch_off+rng_cnt[i] - siz_cnt +j];
                    ch_cnt[i] = ((ch_cnt[i]<<8)|sum_buf); 
                }
            }
    
            for(k= rng_cnt[i]; k< rng_end[i]; k++){
                sum_buf = (uint32_t)pSum[ch_off + k];
                ch_cnt[i] += sum_buf;
            }      
            seq_num[i]  = (uint32_t)pSum[ch_off + rng_end[i]];


            glb_cnt[i] = 0 ;
            if(seq_ch > 0){
                for(j = 0; j<siz_cnt; j++) {
                    sum_buf =  (uint32_t)pSum[ch_off -siz_cnt+j]; 
                    glb_cnt[i] = ((glb_cnt[i]<<8)|sum_buf);
                }
            }
            ch_cnt[i] += glb_cnt[i];
        }// end for
        if(bg_idx==0) ch_cnt[0] += glb_cnt[1];
        upp_num2[0] = 0;upp_num2[1] = 0;              

        if(idx_q[0] < idx_q[1] ){ 
            for(i = 0 ;  i<2 ; i++ ){
                if(bg_idx == 0 && 0 ==i){ 
                    upp_num2[0] += ch_cnt[0];
                    continue;
                }
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                if(idx_r[i] <= 134){
                    for( k = rot_off+(idx_i[i]-idx_r[i])/2;  k < idx_row[i];  k++){
                        get_ch = pBwt[k];
                        ch = get_ch>>4;
                        if(seq_ch ==ch) upp_num2[i]++;
                        ch = get_ch &0xF;
                        if(seq_ch == ch) upp_num2[i]++;
                    }
                    // 该点已经计算过，但需要考虑bg_idx点是否包含在其中++++ 
                    if(idx_r[i] % 2 == 0 ){
                        get_ch = pBwt[idx_row[i]];
                        ch = get_ch>>4;
                        if(seq_ch == ch) ++upp_num2[i];
                    }
                }
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                if(idx_r[i] > 134){  //离下面端点更近的情况运行代码+++
                    upp_num2[i]=(uint32_t)pSum[ch_off +rng_end[i]];
                    rng_len = (idx_q[i] == data_q) ? data_r : 254;
                    rng_end[i] = rot_off+ (idx_i[i] - idx_r[i])/2 + rng_len /2;
                    for(k= idx_row[i]; k< rng_end[i]; k++){
                        get_ch = pBwt[k] ;
                        if(k > idx_row[i] || idx_r[i]%2 >0 ){
                            ch = get_ch >> 4;
                            if(seq_ch ==ch)  --upp_num2[i];
                        }
                        ch = get_ch & 0xF ;
                        if(seq_ch ==ch)  --upp_num2[i];
                    }
                    if((rng_len %2 >0)&&(idx_r[i] < rng_len-1)){
                        get_ch = pBwt[rng_end[i]];
                        ch = get_ch >>4;
                        if(seq_ch == ch){ 
                            --upp_num2[i];
                        } 
                    }
                }
                upp_num2[i] += ch_cnt[i] ;
            } // End: ++++++++++++++ for(i = 0 ;  i<2 ; i++ )
        }// End : +++++++ if(bg_idx_q < ed_idx_q )+++++++++++++++++++
        // 以下是bg_idx和ed_idx在同一个间隔内++++++++++++++++++++++++
        uint32_t bot_num = 0, mid_num = 0;
        if(idx_q[0] == idx_q[1] ) {
            rng_len = (data_q == idx_q[1]) ? data_r : 254;  
            bgn_row = rot_off+(idx_i[0]+1)/2;
            end_row = rot_off+(idx_i[1]+1)/2;
            fst_row = rot_off +(idx_i[0]-idx_r[0])/2;
            lst_row = fst_row + rng_len/2;
            if( (idx_r[0] <= rng_len - idx_r[1]) && bg_idx >0 ){  
                //上段距离小于等于下段距离++++++++
                for(row= fst_row ; row < bgn_row; row++){
                    get_ch = pBwt[row];
                    ch = get_ch>>4;
                    if(seq_ch ==ch) upp_num2[0]++;
                    ch = get_ch &0xF;
                    if(seq_ch ==ch) upp_num2[0]++;
                }
            }
            //+++++++++++++++++++++++++++++++++++++++++++
            if(bg_idx%2 >0){
                get_ch = pBwt[bgn_row]; 
                ch = get_ch>>4;
                if(seq_ch ==ch) upp_num2[0]++;
                ch = get_ch&0xF;
                if(seq_ch==ch)  upp_num2[1]++;
                bgn_row++;
            } else if(ed_idx == 0){
                get_ch = pBwt[bgn_row]; 
                ch = get_ch>>4;
                if(seq_ch==ch)  upp_num2[1]++;
                bgn_row++;
            }
            for(row= bgn_row ; row < end_row ; row++){
                get_ch = pBwt[row];
                ch = get_ch>>4;
                if(seq_ch ==ch) upp_num2[1]++;
                ch = get_ch &0xF;
                if(seq_ch ==ch) upp_num2[1]++;
            }
            if(ed_idx%2==0 ) { // ed_idx端边界处理++++++++++++++++++++++++
                get_ch = pBwt[end_row] ;              
                ch = get_ch >>4;
                if(seq_ch ==ch && ed_idx >0) upp_num2[1]++;
                if(idx_r[1] < rng_len-1){  // idx_r[1]是ed_idx_r
                    ch = get_ch & 0xF;
                    if(seq_ch ==ch) {
                        bot_num++;
                    }
                    end_row++;
                } 
            }
            if( idx_r[0] > rng_len - idx_r[1] ) { //上段距离大于等于下段距离++++++++
                for(row= end_row; row < lst_row; row++){
                    get_ch = pBwt[row]; 
                    ch = get_ch>>4;        
                    if(seq_ch ==ch){ 
                        bot_num++;
                    }
                    ch = get_ch &0xF;    
                    if(seq_ch ==ch) {
                        bot_num++;
                    }
                }
                if(rng_len%2 > 0 && idx_r[1] < rng_len-1){ // 最末端边界+++
                    get_ch = (*(pBwt + lst_row));              
                    ch = get_ch >>4;   
                    if(seq_ch ==ch) bot_num++;
                }
                // 计算bg_idx点开始第一个出现seq_ch的为止，所有seq_ch出现的个数++
                upp_num2[0] = seq_num[1] - bot_num - upp_num2[1];
            } // End: +++++++++++++if( idx_r[0] > rng_len – idx_r[1] ) 
            upp_num2[0] += ch_cnt[1];
            upp_num2[1] += upp_num2[0]; 
        }  //  End: +++++++++++++ if(idx_q[0] == idx_q[1] )
        mid_num = upp_num2[1] - upp_num2[0];
        int alg_err_num = 0;
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ++cur_len;
        if(mid_num >=1 && cur_len < 7){
            bg_idx =  upp_num2[0]; 
            ed_idx =  upp_num2[1]-1;            
            continue; 
        }
        if(err_len>0 && mid_num < 1 && cur_len <= 7 ) { return -2;}
        if(cur_len == 7 && mid_num >=1) {
            uint32_t i_idx, j_idx, j_rot, j_pos;
            int err_pos= (pos+7)%8;
            bg_idx =  upp_num2[0]; 
            ed_idx =  upp_num2[1]-1;
            if(InDel == -1) seq_pos = (seq_pos+len_seq-1)%len_seq;
            else if(InDel == 1) seq_pos = (seq_pos+len_seq-3)%len_seq;
            else if(InDel == 2) seq_pos = (seq_pos+len_seq-4)%len_seq;
            for(i_idx = bg_idx ;  i_idx <= ed_idx ;  i_idx++){
                j_idx = i_idx;    
                j_rot= (rot+1)%8;
               
                rot_off =j_rot*((n_data+1)/2);
                uint32_t j_row = rot_off+(j_idx+1)/2;
                if(j_idx%2>0) {
                    ch = pBwt[j_row-1]&0xF; 
                } else{
                    ch = pBwt[j_row]>>4;
                }

                edit_indel[0] = InDel;
                if(InDel < 0){ 
                    if(InDel == -2){
                        edit_indel[1] = seq_pos;
                        edit_indel[3] = ch;
                    } else if(InDel == -1){
                        if(seq[seq_pos] == (ch&3) ) {
                            edit_indel[1] = (seq_pos+len_seq-1)%len_seq;
                            edit_indel[3] = ((ch>>2)&3);
                        } else if(seq[seq_pos] == (ch>>2) ) {
                            edit_indel[1] = seq_pos;
                            edit_indel[3] = (ch&3);
                        } else{
                            continue; 
                        }
                    } 
                } 
                int edit_buf[4];
                if(InDel >0) {
                    if(InDel == 1){
                        edit_buf[0] = (seq[seq_pos]<<2)|(seq[seq_pos+1]);
                        edit_buf[1] = (seq[seq_pos]<<2)|(seq[seq_pos+2]);
                        edit_buf[2] = (seq[seq_pos+1]<<2)|(seq[seq_pos+2]);
                        for(i = 0; i < 3; ++i){
                            if(ch == edit_buf[i]) break; 
                        }

                        if(i == 3) {continue;} 
                        else{ 
                            if(i == 0) {                            
                                edit_indel[1] = (seq_pos+2)%len_seq;
                            } else if(i == 1) {                            
                                edit_indel[1] = (seq_pos+1)%len_seq;
                            } else if(i == 2) {                            
                                edit_indel[1] = (seq_pos+0)%len_seq;
                            }
                            edit_indel[3] = seq[edit_indel[1]]; 
                        }
                    } else if(InDel == 2){
                        edit_buf[0] = (seq[seq_pos]<<2)|(seq[seq_pos+1]);
                        edit_buf[1] = (seq[seq_pos]<<2)|(seq[seq_pos+3]);
                        edit_buf[2] = (seq[seq_pos+2]<<2)|(seq[seq_pos+3]);
                        for(i = 0; i < 3; ++i){
                            if(ch == edit_buf[i]) break; 
                        }

                        if(i == 3) {continue;}
                        else{
                            if(i == 0) {                            
                                edit_indel[1] = (seq_pos+2)%len_seq;
                            } else if(i == 1) {                            
                                edit_indel[1] = (seq_pos+1)%len_seq;
                            } else if(i == 2) {                            
                                edit_indel[1] = (seq_pos+0)%len_seq;
                            }
                            edit_indel[3] = (seq[edit_indel[1]]<<2) | (seq[edit_indel[1]+1]);
                        }
                    } 
                }
    /*  
                i = 0; 
                //比较后两个bits
                if((ch & 3) == (err_ch &3)){ ++i; }
                //比较后3、4位置bit
                if((ch & 12) == (err_ch &12)){ ++i; }
                //仅保留单碱基匹配错误
                if(i == 0 || i == 2) continue; 
    */
                j_pos = pos;
                // 查找某一近似比对成功序列的行号 ++++++++++++++++++
                /*计算编辑距离*/
                while(j_pos != 0){
                    rot_off = j_rot*((n_data+1)/2);
                    uint32_t j_row = rot_off +( j_idx+1)/2;
                    lst_row = rot_off  + n_data/2 ;
                    if(j_idx%2>0)  {
                        get_ch =pBwt[j_row-1];
                        seq_ch =  get_ch&0xF;
                    } else if(j_idx%2==0)  {
                        get_ch =pBwt[j_row];
                        seq_ch  =  get_ch>>4;
                    }
                    ch_off = 16*j_rot*len_sum + len_sum*seq_ch ;
                    idx_i[0] = j_idx;
                    idx_q[0] = idx_i[0]/254;
                    idx_r[0] = idx_i[0]%254;
                    idx_row[0] = rot_off +( idx_i[0]+1)/2;
                    rng_cnt[0] = (8+siz_cnt)*(idx_q[0] / 8);
                    rng_end[0] = rng_cnt[0] + idx_q[0] % 8;
                    mid_num = 0;
                    upp_num = 0;
                    bot_num = 0;
                    ch_cnt[0] = 0;
                    if(rng_cnt[0]>0){
                        for(j = 0; j< siz_cnt; j++){
                            ch_cnt[0]=(ch_cnt[0]<<8)|(uint32_t)(*(pSum+ch_off+rng_cnt[0]- siz_cnt +j)); 
                        }
                    }
                    for(k= rng_cnt[0]; k< rng_end[0]; k++){ 
                        ch_cnt[0] += (uint32_t)(*(pSum+ ch_off + k ));
                    }       
                    seq_num[0]  = (uint32_t)(*(pSum+ ch_off + rng_end[0] ));
                    glb_cnt[0] = 0 ;
                    if(seq_ch > 0){
                        for(j = 0; j<siz_cnt; j++) 
                            glb_cnt[0] = (glb_cnt[0]<<8) | (uint32_t)(*(pSum+ch_off - siz_cnt+j)) ; 
                    }
                    ch_cnt[0] += glb_cnt[0];
                    bot_num = 0 ;
                    bot_num=(uint32_t)pSum[ch_off +rng_end[0]];
                    rng_len = 254 ;
                    if(idx_q[0]==data_q) rng_len = data_r;
                    rng_end[0] = rot_off+ (idx_i[0] - idx_r[0])/2 + rng_len/2;
                    for(k= idx_row[0]; k< rng_end[0] ; k++){
                        get_ch = pBwt[k] ;
                        ch = get_ch>>4;
                        if(seq_ch ==ch) { 
                            if(k > idx_row[0] || idx_r[0]%2 >0 ){--bot_num;}
                        }
                        ch = get_ch&0xF ;
                        if(seq_ch ==ch)  --bot_num;
                    }
                    if((rng_len %2 >0 ) &&(rng_len-1 > idx_r[0])){  
                        get_ch = pBwt[rng_end[0]];
                        ch = get_ch>>4;
                        if(seq_ch == ch) --bot_num; 
                    }
                    bot_num += ch_cnt[0];
                    j_idx = bot_num-1;
                    j_pos = (j_pos +7)%8;
                    j_rot = (j_rot +1)%8;
                }// End : while(j_pos != 0) ----------------------------------
                int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));

                if(flag_aln == 0){
                    n_aln++;
                    aln_out[n_aln][0] = j_idx;
                    aln_out[n_aln][1] = penalty;
                    aln_out[n_aln][2] = err_pos;
                    aln_out[n_aln][3] = InDel; 
                    algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                    if(err_len > 0) r_flag = 2;
                    else r_flag = 1; 
/*  
{
uint8_t _seq[16];
test_smbwt_retire_seq_large(Bwt, n_data, j_idx, _seq);
//fprintf(stderr, "\n\n%u, aln_seq[%u]: ", __LINE__, j_idx);
for(i = 0; i < 16; ++i) //fprintf(stderr, "%u ", _seq[i]);
//fprintf(stderr, "\n\n");
}
*/


                }

            }// for(i_idx = bg_idx;i_idx <= ed_idx; i_idx ++)--------------------------
            if(aln_out[0][0] == n_aln && err_len == 0) {
                mid_num = 0;
                cur_len--;
            }else{
                aln_out[0][0] = n_aln; 
                aln_out[0][1] = n_exact_aln; 
                break; 
            }
        }
        // 重新选定比对起始位置进行比对 ++++++++++++++++
        if(err_len == 0 && mid_num==0 && cur_len <= 7){             
            alg_num++;
            if(alg_num >7) {break;}
            if(cur_len +st_len <7) { break;} 
            if(pos == 0) {
                if(cur_len == 7 ) pos++;
                else {
                    pos = 7;
                    alg_num = 7; 
                }
            }
            seq_pos = pos*2; 
            st_len = cur_len+st_len-7;
            aln_in[(pos+1)%8] =1;
            err_num = 0;    
            cur_len = 0 ;
            bg_idx = 0; 
            ed_idx = n_data-1;
            continue;
        }// End : +++++++++++++++++++if((mid_num==0)
    }//++++++++++++++++++++++++++++ while(1)
    return r_flag;
} // End : Alignm_big( )------------------------------------------------------
//++++++++++++++++++++++++++++
int gen_bwt_idx(idx_t *fm_idx,  uint32_t hash_boundry[], query_t *query, int seed_off, int len, uint32_t out[2]) 
{   
    int i; 
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint8_t pseed[LEN_READ];
    for(i = 0; i < LEN_READ; ++i) {
        if(read_seq[i] > 3) pseed[i] = 0;
        else pseed[i] = read_seq[i];
    }  
    uint32_t k, l, num;
    k = 0;
    l = fm_idx->bns->l_pac;

    num = bwt_match_exact_alt(fm_idx->bwt, len, pseed+(seed_off-len), &k, &l);
    out[0] = k;
    out[1] = num;
    return num;
}
int gen_ext_idx(idx_t *fm_idx,  uint32_t hash_boundry[], query_t *query, int seed_off, int max_cls, struct ExtBlck *in_blck, struct SubBuf *sub, struct JmpMod *jmp) 
{   
    int i, len = max_cls; 
    struct ExtBlck *cB;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    //uint8_t pseed[SEED_LEN];   
//fprintf(stderr, "%u, query->is_rev = %d\n", __LINE__, query->is_rev); 
    uint8_t pseed[LEN_READ];
    for(i = 0; i < LEN_READ; ++i) {
        if(read_seq[i] > 3) pseed[i] = 0;
        else pseed[i] = read_seq[i];
    }  
    uint32_t k, l, num;
/*  
    for(i = seed_off; i < seed_off+SEED_LEN; ++i) {
        pseed[i-seed_off] = read_seq[i] >3?0:read_seq[i]; 
    }

    
    uint32_t seq12 = lkt_seq2LktItem(pseed, 8, 19);
    k = fm_idx->fastmap->item[seq12];
    l = fm_idx->fastmap->item[seq12+1]-1;
    l -= get_12mer_correct(hash_boundry, l);
    num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, pseed, &k, &l);
*/
/*  
    uint32_t seq12 = lkt_seq2LktItem(pseed+seed_off-12, 0, 11);
    k = fm_idx->fastmap->item[seq12];
    l = fm_idx->fastmap->item[seq12+1]-1;
    l -= get_12mer_correct(hash_boundry, l);
//fprintf(stderr, "%u, len = %u, seed_off = %u, k = %u, l = %u\n", __LINE__, len, seed_off, k, l);
    num = bwt_match_exact_alt(fm_idx->bwt, len -12, pseed+seed_off-len, &k, &l);
*/
    //fprintf(stderr, "%u,seed_off = %d, len = %d\n", __LINE__, seed_off, len);
    k = 0;
    l = fm_idx->bns->l_pac;
//fprintf(stderr, "k = %u, l = %u\n", k, l);
    num = bwt_match_exact_alt(fm_idx->bwt, len-2, pseed + (seed_off-len)+1, &k, &l);
    sub->pair_out->pair_arry->idx_bg = k;
    sub->pair_out->pair_arry->idx_num = num;
    return num;

    //++++++++++++++++++++++++++++++++++++ 
    int flg = 0, cls = 0;

    uint32_t buf_algn[20];
    buf_algn[0] = SEED_LEN;
    buf_algn[1] = k;
    buf_algn[2] = num;
    getCapPos(jmp,buf_algn);
    uint32_t cap_row = buf_algn[3];
    cls = 1;
    while(cls <= max_cls) {
        cB = in_blck + cls-1;  
        clean_sub_buf(cB, sub);
        int L_offset = seed_off - cls*16;
        int R_offset = seed_off + SEED_LEN +  cls*16 - 16;
        memcpy(sub->ext_seqL, read_seq+L_offset-2, 18);
        memcpy(sub->ext_seqR, read_seq+R_offset, 18);
        //处理左段序列N
        //err->seq_L[0] = 0;
        for(i = 0; i < 18; ++i) {
            if(sub->ext_seqL[i] == 4) {
                sub->ext_seqL[i] = 0;
            }     
        }
        //处理右段序列N
        //err->seq_R[0] = 0;
        for(i = 0; i < 18; ++i) {
            if(sub->ext_seqR[i] == 4) {
                //err->seq_R[++err->seq_R[0]];
                sub->ext_seqR[i] = 0;
            }     
        }
        //左序列反转 
        for(i = 0; i < 9; ++i){
            uint8_t ch_buf = sub->ext_seqL[i]; 
            sub->ext_seqL[i] = sub->ext_seqL[17-i];
            sub->ext_seqL[17-i] = ch_buf;
        }

        struct CapIfo  *cur_cap;
        cur_cap = cB->head_cap + cap_row;
        cB->relat     = cur_cap->relat;	
        cB->nxtpnt    = cur_cap->nxtpnt;
        cB->smbwt     = cur_cap->smbwt;
        cB->num_seqL  = cur_cap->num_seqL;
        cB->num_seqR  = cur_cap->num_seqR;
        cB->num_relat = cur_cap->num_relat;
        cB->extidx    = cap_row;
        cB->bwtL  =  cB->smbwt ;
        uint32_t len = get_bwt_size(cB->num_seqL);
        cB->bwtR  =  cB->smbwt + len;	
//fprintf(stderr, "%u, num_seqL = %u, num_seqR = %u, relat = %u\n",__LINE__, cB->num_seqL, cB->num_seqR, cB->num_relat);       
        //-----------------------------
        flg = Algn_sub_once(cB,sub,0);
//fprintf(stderr, "%u, cls = %u\n", __LINE__, cls);
        if(flg != 1) {
            printf("%u, Error, left exact match fail!!!\n", __LINE__);
            exit(1); 
        }         
        flg = Algn_sub_once(cB,sub,1);
//fprintf(stderr, "%u, cls = %u\n", __LINE__, cls);
        if(flg != 1) {
            printf("%u, Error, right exact match fail!!!\n", __LINE__);
            exit(1); 
        }         

        sub->seqL_aln_old = 0;
        sub->seqR_aln_old = 0;
        //set_pair_data(sub, call+sTree->cls); //设定配对数据 
        sub->query_len = query->l_seq;
        sub->query_b0 = query->b0;
        sub->query_err = query->l_seq - query->b0;

//fprintf(stderr, "%u, seqL = %u\n", __LINE__, sub->seqL_out[0][0]);
//fprintf(stderr, "%u, seqR = %u\n", __LINE__, sub->seqR_out[0][0]);
        flg = PairExtSeq_all(cB,sub); //完成配对过程
//fprintf(stderr, "%u, cls = %u, p_num = %u\n", __LINE__, cls, flg);
        if(flg != 1) {
            printf("%u, pair_num != 1\n", __LINE__);
            exit(1); 
        }
        cap_row = sub->pair_out->pair_arry[0].nxtpnt;
        cls++;

    }

//fprintf(stderr, "%u, cls = %u\n", __LINE__, cls);
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    return flg;
}//end gen_aln_data()++++++++++++++++++++++
int get_max_err(struct SubBuf *sub) 
{
    //int sub_err = 5;
    //int ins_err = 6; 
    //int del_err = 7;
    //int delta = sub->aln_r->delta;
    //int cur_err = sub->err_sum[0];
    //int query_err = sub->query_err;
    //int *path_err = sub->path_err;
    //int cls = sub->cls;
    int max_err;
    max_err = sub->query_err + sub->delta - sub->err_sum[0];
    
    //max_err = max_err - sub->err_sum[0];
    //sub->err_sum[1] = sub_err;
    //sub->err_sum[2] = ins_err;
    //sub->err_sum[3] = del_err;
    return max_err;
    
}
    /*   
    if(query_err < sub_err) {
        max_err = delta; 
    } else if(query_err < 2*sub_err){
        max_err = delta+sub_err; 
    } else {
          
        max_err = path_err[cls];
        if(sub->err_sum[0] > path_err[cls]) {
            max_err = 1; 
        } else{
            max_err = delta; 
        }
        
        if(query_err <= sub->err_sum[0]){
            max_err = 1; 
        } else if(  query_err <= sub->err_sum[0]+delta ){
            max_err = delta; 
        } else{
            max_err = query_err - sub->err_sum[0] +delta/2;
        } 
    }
    */
   


int gen_aln_data(struct ExtBlck *cB, struct SubBuf *sub, struct call_t *call) 
{
    int *indel_f = sub->indel_f; 
    int max_err = get_max_err(sub);
    int sub_err = sub->sub_err;
    int ins_err = sub->ins_err; 
    int del_err = sub->del_err;
    int r_flag = 0;
    int old_num_L = 0;
    int old_num_R = 0;
 
    int aln_f  = call->aln_f;
    int aln_r  = call->aln_r;
    int pair_f = call->pair_f;
//fprintf(stderr, "%u, aln_f = %d\n", __LINE__, aln_f);
    if(aln_f == 0) {//首次调用时完成查找左右两段的精确比对或者置换比对
        r_flag = -5; 
        call->sub_O_L = Algn_sub_once(cB,sub,0);
        if(call->sub_O_L >1) {
            if(max_err > sub_err) {
                call->sub_A_L = Algn_sub_all(cB, sub, 0);
                call->num_L[0] = 0;
                call->num_L[1] = sub->seqL_out[0][0];
            } else{
                call->num_L[0] = 0;
            }
        } else if(call->sub_O_L == 1){
            call->num_L[0] = 1;
            call->num_L[1] = 0;
        } else{
            call->num_L[0] = 0;
            call->num_L[1] = 0;
        }
        call->sub_O_R = Algn_sub_once(cB,sub,1);
        if(call->sub_O_R >1) {
            if(max_err > sub_err) {
                call->sub_A_R = Algn_sub_all(cB, sub, 1);
                call->num_R[0] = 0;
                call->num_R[1] = sub->seqR_out[0][0];
            } else{
                call->num_R[0] = 0;
            }
        } else if(call->sub_O_R == 1){
            call->num_R[0] = 1;
            call->num_R[1] = 0;
        } else{
            call->num_R[0] = 0;
            call->num_R[1] = 0;
        }
        if(call->sub_O_L ==1 && call->sub_O_R == 1) aln_r = 1;
        if(call->sub_O_L ==1 && call->sub_O_R >  1) aln_r = 2;
        if(call->sub_O_L > 1 && call->sub_O_R == 1) aln_r = 3; 
        if(call->sub_O_L > 1 && call->sub_O_R >  1) aln_r = 4; 
        if(call->sub_O_L ==1 && call->sub_O_R <= 0) aln_r = 5;
        if(call->sub_O_L <=0 && call->sub_O_R == 1) aln_r = 6;   
        if(call->sub_O_L > 1 && call->sub_O_R <= 0) aln_r = 7;
        if(call->sub_O_L <=0 && call->sub_O_R  > 1) aln_r = 8;
        if(call->sub_O_L <=0 && call->sub_O_R <= 0) aln_r = 9;
//fprintf(stderr, "%u, aln_r = %d\n", __LINE__, aln_r); 
        call->aln_f  = aln_f;
        call->aln_r  = aln_r;
        call->pair_f = pair_f;
        if(aln_r > 6 ) {//修改日期:2018-07-07
            return -1;
        } else if(aln_r == 6) {//修改日:2018-07-08
            if(max_err > del_err && indel_f[6] > 0){
                aln_f = 2;//左段插入删除比对                          
            } else{
                call->aln_f = 5;
                return -1; 
            }
        } else if(aln_r == 5)  {//修改日: 2018-07-08
            if(max_err > del_err && indel_f[5] > 0){
                aln_f = 3;//右段插入删除比对 
            } else{
                call->aln_f = 5;
                return -1;  
            }
        } else{
            if(aln_r == 1) {
                pair_f = 0;   
            } 
            if(aln_r == 2 || aln_r == 3) {
                if(max_err > sub_err) {
                    pair_f = 0;     
                } else{
                    call->aln_f = 5;
                    return -1;   
                }
            }
            if(aln_r == 4) {
                if(max_err > 2*sub_err) {
                    pair_f = 0;     
                } else{
                    call->aln_f = 5;
                    return -1;   
                }
            }
        } 
        call->aln_f  = aln_f;
        call->pair_f = pair_f;
        r_flag = 5; 
    }//end  if(aln_f == 0) ++++++++++++++++++++++++++++++++++++++++++              
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(aln_f == 1){//再次（非首次）调用时完成左右两段置换近似比对
        r_flag = -1;
        if(max_err <= sub_err) { 
            call->aln_f = 5;
            return -1;  
        }
        int r_L = 0;
        int r_R = 0;
        if(aln_r == 1 || aln_r == 2) { 
            if(call->sub_O_L == 1) {//aln_r = 1, 2的情况
                old_num_L = sub->seqL_out[0][0];
                call->sub_A_L = Algn_sub_all(cB, sub, 0);
                call->num_L[1] = sub->seqL_out[0][0]-old_num_L;
                r_L = call->num_L[1];
            }
        }
        if( aln_r == 1 || aln_r == 3 ){
            if(call->sub_O_R == 1) {//aln_r = 1, 3的情况
                old_num_R = sub->seqR_out[0][0];
                call->sub_A_R = Algn_sub_all(cB, sub, 1);
                call->num_R[1] = sub->seqR_out[0][0]-old_num_R;
                r_R = call->num_R[1];
            }
        }
        if(r_L + r_R == 0) {
            if(max_err <= del_err) { 
                call->aln_f = 5;
                return -1; 
            }
            if(aln_r == 1 || aln_r == 3) aln_f = 2; 
            else if(aln_r == 2) aln_f = 3;
            else {
                printf("%u, %s error!!\n", __LINE__, __func__);
                exit(1); 
            }
        } else {
        /*  
            if(r_L >0 && r_R >0) {
                pair_f = 1;//两侧都有增量
            } else if(r_L >0){
                pair_f = 1;//只有左侧有增量
            } else if(r_R >0){
                pair_f = 1;//只有右侧有增量 
            }
        */
            pair_f = 1;
        }
        call->aln_f  = aln_f;
        call->pair_f = pair_f;
        r_flag = 1;
    }
    
    if(aln_f  == 2){//首次或再次调用时完成左侧插入删除查找比对
    /*  
        if(pair_f == 0 && aln_r == 6 && call->sub_O_R ==1) {
            call->sub_A_R = Algn_sub_all(cB, sub, 1);
            call->num_R[1] = sub->seqR_out[0][0] - call->num_R[0];
        }
    */
        r_flag = -2;
        //int i;
        if(aln_r == 1 || aln_r == 3 || aln_r == 6) { 
            sub->aln_L[13] = -1;//删除一个碱基的比对
            old_num_L = sub->seqL_out[0][0];
            call->num_L[2] = 0;
            call->ind_O_L[0] = Algn_indel_once(cB, sub, 0);
            if(call->ind_O_L[0] >0) {
                call->ind_A_L[0] = Algn_indel_all(cB, sub, 0);
                call->num_L[2] = sub->seqL_out[0][0]-old_num_L;
            }
          
            sub->aln_L[13] = 1;//插入一个碱基的比对
            old_num_L = sub->seqL_out[0][0];
            call->num_L[3] = 0;
            call->ind_O_L[1] = Algn_indel_once(cB, sub, 0);

            if(call->ind_O_L[1] >0) {
                call->ind_A_L[1] = Algn_indel_all(cB, sub, 0);
                call->num_L[3] = sub->seqL_out[0][0]-old_num_L;
            }
        
            call->aln_f  = aln_f;
            call->pair_f = pair_f;
            if(call->num_L[2] + call->num_L[3] == 0) {
            /*  
                if(aln_r == 6 || aln_r == 8) {
                    return r_flag; 
                }
            */
                if(aln_r == 1) aln_f = 3;//考虑到aln_r = 1, 3, 4保留的出口
                else {
                    call->aln_f = 5;
                    return -1;
                }
            } else {
                if(aln_r == 6) pair_f = 0;
                else pair_f = 2; 
            }
        }
        call->aln_f  = aln_f;
        call->pair_f = pair_f;
        r_flag = 2;
    }
    if(aln_f  == 3){//首次或再次调用时完成右侧插入删除比对
    /*  
        if(aln_r == 5 && call->sub_O_L ==1) {
            call->sub_A_L = Algn_sub_all(cB, sub, 0);
            call->num_L[1] = sub->seqL_out[0][0] - call->num_L[0];
        }
    */
    
        r_flag = -3;
        if(aln_r == 1 || aln_r == 2 || aln_r == 5) { 
            sub->aln_R[13] = -1;
            old_num_R = sub->seqR_out[0][0];
            call->num_R[2] = 0;
            call->ind_O_R[0] = Algn_indel_once(cB, sub, 1);
            if(call->ind_O_R[0] >0) {
                call->ind_A_R[0] = Algn_indel_all(cB, sub, 1);
                call->num_R[2] = sub->seqR_out[0][0]-old_num_R;
            }                     
            sub->aln_R[13] = 1; 
            old_num_R = sub->seqR_out[0][0];
            call->num_R[3] = 0;
            call->ind_O_R[1] = Algn_indel_once(cB, sub, 1);
            if(call->ind_O_R[1] >0) {
                call->ind_A_R[1] = Algn_indel_all(cB, sub, 1);
                call->num_R[3] = sub->seqR_out[0][0]-old_num_R;
            }
            if(call->num_R[2] + call->num_R[3] == 0) {
                //call->aln_f  = aln_f;
                //call->pair_f = pair_f;
                call->aln_f  = 5;
                return -1;
            }else{
                if(aln_r == 5) pair_f = 0; 
                else pair_f = 3; 
            }
        }
        call->aln_f  = aln_f;
        call->pair_f = pair_f;
        r_flag = 3;  
    }
    return r_flag;
    //以下代码无效
    if(aln_f == 4) {//再次调用时左右两侧都进行插入删除比对
        r_flag = -4;
        int r_L = 0;
        int r_R = 0;
        sub->aln_L[13] = -1;
        call->ind_O_L[0] = Algn_indel_once(cB, sub, 0);
        if(call->ind_O_L[0] > 0){
            r_L = 1;
        }else{
            sub->aln_L[13] = 1;
            call->ind_O_L[1] = Algn_indel_once(cB, sub, 0);
            if(call->ind_O_L[1] >0) {
                r_L = 2;
            } 
        } 
        if(r_L == 0) {
            aln_f = 0;
            call->aln_f  = aln_f;
            return r_flag; 
        }
        sub->aln_R[13] = -1;
        call->ind_O_R[0] = Algn_indel_once(cB, sub, 1);
        if(call->ind_O_R[0] > 0){
            r_R = 1;
        }else{
            sub->aln_R[13] = 1;
            call->ind_O_R[1] = Algn_indel_once(cB, sub, 1);
            if(call->ind_O_R[1] >0) {
                r_R = 2;
            } 
        } 
        if(r_R == 0) {
            aln_f = 0;
            call->aln_f  = aln_f;
            return r_flag; 
        }

        if(r_L == 1) {
            sub->aln_L[13] = -1;
            call->ind_A_L[0] = Algn_indel_all(cB, sub, 0);
            call->num_L[2] = sub->seqL_out[0][0];
            sub->aln_L[13] = 1;
            old_num_L = sub->seqL_out[0][0];
            call->ind_O_L[1] = Algn_indel_once(cB, sub, 0);
            if(call->ind_O_L[1] >0) {
                call->num_L[3] = 0;
                call->ind_A_L[1] = Algn_indel_all(cB, sub, 0);
                call->num_L[3] = sub->seqL_out[0][0]-old_num_L;
            } 
        } else{
            sub->aln_L[13] = 1;
            call->ind_A_L[1] = Algn_indel_all(cB, sub, 0);
            call->num_L[2] = sub->seqL_out[0][0];
        }
        if(r_R == 1) {
            sub->aln_R[13] = -1;
            call->ind_A_R[0] = Algn_indel_all(cB, sub, 1);
            call->num_R[2] = sub->seqR_out[0][0];
            sub->aln_R[13] = 1;
            old_num_R = sub->seqR_out[0][0];
            call->ind_O_R[1] = Algn_indel_once(cB, sub, 1);
            if(call->ind_O_R[1] >0) {
                call->num_R[3] = 0;
                call->ind_A_R[1] = Algn_indel_all(cB, sub, 1);
                call->num_R[3] = sub->seqR_out[0][0]-old_num_R;
            } 
        } else{
            sub->aln_R[13] = 1;
            call->ind_A_R[1] = Algn_indel_all(cB, sub, 1);
            call->num_R[2] = sub->seqR_out[0][0];
        }
        pair_f = 0;
        call->aln_f  = aln_f;
        call->pair_f = pair_f;
        r_flag = 4;
    }
    return r_flag;
}//end gen_aln_data()++++++++++++++++++++++

int set_pair_data(struct SubBuf *sub, struct call_t *call) 
{
    int aln_f = call->aln_f;
    int aln_r = call->aln_r;
    int pair_f = call->pair_f;
    int rep_i = call->rep_i;
    int r_flag = 0;

    if(pair_f ==0) {
        sub->pair_b2e[0][0] = 1;  
        sub->pair_b2e[0][1] = sub->seqL_out[0][0]; 
        sub->pair_b2e[0][2] = 1; 
        sub->pair_b2e[0][3] = sub->seqR_out[0][0];
        sub->pair_b2e_num = 1;
    }

if(pair_f  == 1) {
    if(aln_r == 1) {
        if(call->num_L[1] == 0 && call->num_R[1] == 0) { 
            sub->pair_b2e_num = 0;
        }
        if(call->num_L[1] > 0 && call->num_R[1] == 0) {
            sub->pair_b2e[0][0] = sub->seqL_out[0][0]-call->num_L[1]+1;  
            sub->pair_b2e[0][1] = sub->seqL_out[0][0]; 
            sub->pair_b2e[0][2] = 1; 
            sub->pair_b2e[0][3] = 1;
            sub->pair_b2e_num = 1;
        }
        if(call->num_L[1] == 0 && call->num_R[1] > 0) {
            sub->pair_b2e[0][0] = 1;  
            sub->pair_b2e[0][1] = 1; 
            sub->pair_b2e[0][2] = sub->seqR_out[0][0]-call->num_R[1]+1; 
            sub->pair_b2e[0][3] = sub->seqR_out[0][0];
            sub->pair_b2e_num = 1;
        }
   
        if(call->num_L[1] > 0 && call->num_R[1] > 0) {
            sub->pair_b2e[0][0] = sub->seqL_out[0][0]-call->num_L[1]+1;  
            sub->pair_b2e[0][1] = sub->seqL_out[0][0]; 
            sub->pair_b2e[0][2] = 1; 
            sub->pair_b2e[0][3] = 1;
            
            sub->pair_b2e[1][0] = 1;  
            sub->pair_b2e[1][1] = 1; 
            sub->pair_b2e[1][2] = sub->seqR_out[0][0]-call->num_R[1]+1; 
            sub->pair_b2e[1][3] = sub->seqR_out[0][0];

            int max_err = get_max_err(sub);
            int sub_err = sub->sub_err;
 
            if(max_err > 2*sub_err) {
                sub->pair_b2e[2][0] = sub->seqL_out[0][0]-call->num_L[1]+1;  
                sub->pair_b2e[2][1] = sub->seqL_out[0][0]; 
                sub->pair_b2e[2][2] = sub->seqR_out[0][0]-call->num_R[1]+1; 
                sub->pair_b2e[2][3] = sub->seqR_out[0][0];
                sub->pair_b2e_num = 3;
            } else{
                sub->pair_b2e_num = 2;
            }
        } 
    }
    if(aln_r == 2) {
        if(call->num_L[1] > 0) {
            sub->pair_b2e[0][0] = sub->seqL_out[0][0]-call->num_L[1]+1;  
            sub->pair_b2e[0][1] = sub->seqL_out[0][0]; 
            sub->pair_b2e[0][2] = 1; 
            sub->pair_b2e[0][3] = sub->seqR_out[0][0];
            sub->pair_b2e_num = 1;
        } else{
            sub->pair_b2e_num = 0;
        } 
    } 
    if(aln_r == 3) {
        if(call->num_R[1] > 0) {
            sub->pair_b2e[0][0] = 1;  
            sub->pair_b2e[0][1] = sub->seqL_out[0][0]; 
            sub->pair_b2e[0][2] = sub->seqR_out[0][0]-call->num_R[1] + 1; 
            sub->pair_b2e[0][3] = sub->seqR_out[0][0];
            sub->pair_b2e_num = 1; 
        } else{
            sub->pair_b2e_num = 0;
        } 
    }
} 
    if(pair_f == 2) {//处理aln_r = 1, 3, 只有左侧有增量
        /*
        sub->seqL_aln_old = sub->seqL_out[0][0]-call->num_L[1]; 
        sub->seqR_aln_old = 0;       
        sub->seqR_out[0][0] = 1;
        */
        int num_L =  call->num_L[2] + call->num_L[3];
        sub->pair_b2e[0][0] = sub->seqL_out[0][0]-num_L+1;  
        sub->pair_b2e[0][1] = sub->seqL_out[0][0]; 
        sub->pair_b2e[0][2] = 1; 
        sub->pair_b2e[0][3] = 1;
        sub->pair_b2e_num = 1;

    /*   
        if(aln_f == 2) {
            if(rep_i == 0) {
                sub->seqL_aln_old = sub->seqL_out[0][0]-call->num_L[2]-call->num_L[3];
                sub->seqL_out[0][0] = sub->seqL_out[0][0]-call->num_L[3]; 
            }
            if(rep_i == 1){
                sub->seqL_aln_old = sub->seqL_out[0][0];
                sub->seqL_out[0][0] = sub->seqL_out[0][0]+call->num_L[3]; 
            }
        }
    */ 
    } 
    if(pair_f == 3) {//处理aln_r = 1, 2, 只有右侧有增量 
        /*  
        sub->seqL_aln_old = 0;
        sub->seqL_out[0][0] = 1;
        sub->seqR_aln_old = sub->seqR_out[0][0]-call->num_R[1];
        */
       
        int num_R =  call->num_R[2] + call->num_R[3];
        sub->pair_b2e[0][0] = 1;  
        sub->pair_b2e[0][1] = 1; 
        sub->pair_b2e[0][2] = sub->seqR_out[0][0]-num_R +1; 
        sub->pair_b2e[0][3] = sub->seqR_out[0][0] ;
        sub->pair_b2e_num = 1;


    /*   
        if( aln_f == 3)  {
            //fprintf(stderr, "%u, rep_i = %u, call->num_R[2] = %u, call->num_R[3] = %u, sub->seqR_out[0][0] = %u\n", __LINE__, rep_i, call->num_R[2], call->num_R[3], sub->seqR_out[0][0]);
            if(rep_i == 0) {
                sub->seqR_aln_old = sub->seqR_out[0][0]-call->num_R[2]-call->num_R[3];
                sub->seqR_out[0][0] = sub->seqR_out[0][0]-call->num_R[3]; 
//fprintf(stderr, "%u, seqR_aln_old = %u, seqR_out[0][0] = %u\n", __LINE__, sub->seqR_aln_old, sub->seqR_out[0][0]);
            }
            if(rep_i == 1){
                sub->seqR_aln_old = sub->seqR_out[0][0];
                sub->seqR_out[0][0] = sub->seqR_out[0][0]+call->num_R[3]; 
//fprintf(stderr, "%u, seqR_aln_old = %u, seqR_out[0][0] = %u\n", __LINE__, sub->seqR_aln_old, sub->seqR_out[0][0]);
            }
        }
    */
    } 
    r_flag = pair_f;
    return r_flag;
}// end set_pair_data()++++++++++++++ 
int get_next_aln(struct SubBuf *sub, struct call_t *call) 
{ 
    
    int *sub_f = sub->sub_f;
    int *indel_f = sub->indel_f;
    int aln_f = call->aln_f;
    int aln_r = call->aln_r;
    int pair_f = call->pair_f;
    int r_flag = 1;
   
    if(aln_r == 4 || aln_r == 5 || aln_r == 6) {//pair_f == 0
        if(pair_f > 0) {
            printf("%u, aln_r  = %d, pair_f = %d, aln_f = %u\n", 
                    __LINE__, aln_r, pair_f, aln_f);
            printf("%u, %s error!!\n", __LINE__, __func__);
            exit(1); 
        } 
        call->aln_f = 5;
        return r_flag;  
    }  
    if(pair_f == 2) {
        if(aln_r == 3) call->aln_f = 5;
        else if(aln_r == 1){ 
            if(indel_f[1] > 0) {
                call->aln_f = 3; 
            } else{
                call->aln_f = 5;
            }
        }
        return r_flag+1; 
    } 
    if(pair_f == 3) {
        call->aln_f = 5;
        return r_flag+2; 
    }
    if(pair_f == 1) {
        //if(aln_r == 1 || aln_r == 3) {
        if(aln_r == 1) {
            if(indel_f[1] > 0) call->aln_f = 2;
            else call->aln_f = 5;
        } else if(aln_r == 2){
            if(indel_f[2] > 0) call->aln_f = 3;
            else call->aln_f = 5;
        } else if(aln_r == 3) { 
            if(indel_f[3] > 0) call->aln_f = 2;
            else call->aln_f = 5;
        } else{
            printf("ERROR, %u, %s\n", __LINE__, __func__);
            exit(1); 
        } 
        return -1; 
    }
    if(pair_f == 0){//aln_r = 1, 2, 3
        if(aln_r == 1) {
            if(sub_f[1] > 0) call->aln_f = 1;
            else call->aln_f = 5;
        } else if(aln_r == 2){
            if(sub_f[2] > 0) call->aln_f = 1;
            else call->aln_f = 5;
         } else if(aln_r == 3){
            if(sub_f[3] > 0) call->aln_f = 1;
            else call->aln_f = 5;
        } else{
            printf("ERROR, %u, %s\n", __LINE__, __func__);
            exit(1); 
        } 
        return -1;       
    }
    printf("ERROR, %u, %s\n", __LINE__, __func__);
    exit(1); 
/*  
    if(pair_f == 0 ){  
        if(aln_f == 0) {
            if(aln_r ==1) { 
                aln_f = 1; 
            }else if(aln_r <  4) { 
                aln_f = 1; 
            }else if(aln_r == 4) { 
                aln_f = 2; 
            }else {
                printf("%u, pair_f = %u, aln_f = %u, aln_r = %u\n", 
                        __LINE__, pair_f, aln_f, aln_r );
                //fprintf(stderr, "%u, pair_f = %u, aln_f = %u, aln_r = %u\n", 
                        __LINE__, pair_f, aln_f, aln_r );
                exit(1);
            }
        }else if(aln_f == 2) {

            //if(aln_r == 5 || aln_r == 6){ 
            if(aln_r == 6 || aln_r == 8){ 
                aln_f = 5; 
            } else if(aln_r==4){
                aln_f = 3;
            } else{
                printf("%u, pair_f = %u, aln_f = %u, aln_r = %u\n", __LINE__, pair_f, aln_f, aln_r );
                //fprintf(stderr, "%u, pair_f = %u, aln_f = %u, aln_r = %u\n", __LINE__, pair_f, aln_f, aln_r );
                exit(1);
            }
        }else if(aln_f == 3) {

            //if(aln_r == 7 || aln_r == 8){ 
            if(aln_r == 5 || aln_r == 7){ 
                aln_f = 5; 
            }else if(aln_r==4){
                aln_f = 5;
            }else{
                printf("%u, pair_f = %u, aln_f = %u, aln_r = %u\n", 
                        __LINE__, pair_f, aln_f, aln_r );
                //fprintf(stderr,"%u, pair_f = %u, aln_f = %u, aln_r = %u\n",
                        __LINE__, pair_f, aln_f, aln_r );
                exit(1);
            }
        }
        call->aln_f = aln_f;
        call->pair_f = pair_f;
        return r_flag;
    }
    if(pair_f == 1) {
        r_flag = pair_f;
        if(aln_f != 1 || aln_r != 1) {
            printf("%u, pair_f = %u, aln_f = %u, aln_r = %u\n", 
                    __LINE__, pair_f, aln_f, aln_r );
            //fprintf(stderr, "%u, pair_f = %u, aln_f = %u, aln_r = %u\n", 
                    __LINE__, pair_f, aln_f, aln_r );
            exit(1);
        } 
        sub->seqL_out[0][0] = call->num_L[5];
        sub->seqR_out[0][0] = call->num_R[5]; 
        call->aln_f = aln_f;
        call->pair_f = pair_f;
        return r_flag; 
    }
    if(pair_f == 2) {
        r_flag = pair_f;
        if(aln_f == 1 && (aln_r ==1 || aln_r ==2)){ 
            aln_f = 2;
        } else if (aln_f==2 && rep_i == 1) {
            aln_f = 3;
        } 
        call->aln_f = aln_f;
        call->pair_f = pair_f;
        return r_flag;
    }
    if(pair_f == 3) {
        r_flag = pair_f;
        if(aln_f == 1 && (aln_r ==1 || aln_r ==3)){ 
            aln_f = 2;
        } else if (aln_f==3 && rep_i == 1) {
            aln_f = 5;
        }
        call->aln_f = aln_f;
        call->pair_f = pair_f;
 
        return r_flag;
    }
    */
} //end get_next_aln()++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++
/*
if(n_data == 419){
uint32_t __i, __j, __k, __l;
//fprintf(stderr, "\n----------------------\n");
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        //fprintf(stderr, "seq_ch = %u\t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            //fprintf(stderr, "%u\t", pSum[__j+__l]);
        }
        //fprintf(stderr, "\n");
    }
    //fprintf(stderr, "-------------------\n\n");
}
//fprintf(stderr, "\n+++++++++++++++++\n");
}*/
//++++++++++++++++++++++++++++++++++

/*
uint32_t __i, __j, __k, __l;
//fprintf(stderr, "\n----------------------\n");
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        //fprintf(stderr, "seq_ch = %u\t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            //fprintf(stderr, "%u\t", pSum[__j+__l]);
        }
        //fprintf(stderr, "\n");
    }
    //fprintf(stderr, "-------------------\n\n");
}
//fprintf(stderr, "\n+++++++++++++++++\n");

//++++++++++++++++++++++++++++++++++
uint8_t *cnt = calloc(800000, sizeof(uint8_t));
gen_cnt2(Bwt, n_data, cnt);
//fprintf(stderr, "\n----------------------\n");
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        //fprintf(stderr, "seq_ch = %u \t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            //fprintf(stderr, "%u\t", cnt[__j+__l]);
        }
        //fprintf(stderr, "\n");
    }
    //fprintf(stderr, "-------------------\n\n");
}
//fprintf(stderr, "\n+++++++++++++++++\n");




free(cnt);
exit(1);
*/
//++++++++++++++++++++++++++++++++++
//传入变量：种子的index区间; 最长的延伸长度, 延伸序列的起始位置
//返回：返回每次的比对区间序列。
//比对过程中，比对失败终止。

//int bwt_ext_kmer(idx_t *fm_idx, uint8_t *read_seq, uint32_t idx[2], uint16_t s_off[2], uint32_t (*out)[2])
int bwt_ext_kmer(idx_t *fm_idx, uint8_t *read_seq, uint32_t idx[2], int *s_off, uint32_t (*out)[2])
{
    if(s_off[0] == 0) return 0;
    int L_off = s_off[0]-1;
    //int len = s_off[1];
    int len = s_off[1]<L_off? s_off[1]:L_off;
    uint8_t *seq = read_seq;

    uint32_t i, j, k, l, num;

    k = idx[0];
    l = idx[1];

    out[0][0] = k;
    out[0][1] = l;  

    for(i = 0, j = 1; i < len; ++i){

        num = bwt_match_exact_alt(fm_idx->bwt, 1,  seq+L_off-i, &k, &l);

        if(num > 1) {
            out[j][0] = k;
            out[j][1] = l;
            ++j; 
        } else if(num ==1){
            out[j][0] = k;
            out[j][1] = l;
            ++j;
            break;
        } else {
            break;
        }
        if(num < IS_SMLSIZ) break;
    }
    return j-1;
}
//传入变量：种子的index区间; 比对目标序列的起始位置与长度数组
//返回：返回候选SW比对pos集合。
//比对过程中，比对失败终止。

int init_idx_target2(idx_t *fm_idx, uint32_t hash_boundry[], uint8_t *read_seq, uint16_t (*seq_pos)[2], int seq_num, uint32_t (*target_idx)[3])
{

    int bg = seq_pos[0];
    int len = seq_pos[1];
    uint8_t *seq = read_seq+bg, *pseq, seq_buf[12];
    uint32_t i, j, k, l, num, pos, pos_i, i_bwt, seq12;
//fprintf(stderr, "%u, seq_num = %u\n", __LINE__, seq_num);
    for(i = 0; i < seq_num; ++i){
        pseq = read_seq+seq_pos[i][0];
//fprintf(stderr, "seq_pos[%u][0] = %u\n", i, seq_pos[i][0]);
//fprintf(stderr, "seq_pos[%u][1] = %u\n", i, seq_pos[i][1]);
        len = seq_pos[i][1];
        for(k = 0; k < len; ++k){
            seq_buf[k] = pseq[k]; 
        }
        for(k = len; k < 12; ++k){
            seq_buf[k] = 0;
        }
        seq12 = lkt_seq2LktItem(seq_buf, 0, 11);
        target_idx[i][0] = fm_idx->fastmap->item[seq12];

        for(k = len; k < 12; ++k){
            seq_buf[k] = 3;
        }
        seq12 = lkt_seq2LktItem(seq_buf, 0, 11);
        l = fm_idx->fastmap->item[seq12+1]-1;
        l -= get_12mer_correct(hash_boundry, l);
        target_idx[i][1] = l;
        target_idx[i][2] = seq_pos[i][0];

    } 
    return i;
}
int init_idx_target1(idx_t *fm_idx, uint32_t hash_boundry[], query_t *query,  uint32_t (*target_idx)[3])
{
    uint32_t i, j, k, l,  seq12;
    uint8_t seq_buf[12];
    //uint8_t *read_seq = query->read_seq;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    //int bg = query->l_seq % 12 / 2;
    //int ed = bg + query->l_seq / 12 * 12;
    int bg = 0;
    int ed = query->l_seq-11;
    int num = query->l_seq/12;
    int q, r;
//fprintf(stderr, "%u\n", __LINE__);
    for(i = bg; i < ed; ++i) {
        r = (i - bg)%12;
        q = (i - bg)/12;
        int bg_k = i;
        int ed_k = i + 12; 
//fprintf(stderr, "%u, i = %u\n", __LINE__, i);
        for(k = bg_k; k < ed_k; ++k){
            seq_buf[k-bg_k] = read_seq[k]; 
        }
//fprintf(stderr, "%u, k = %u\n", __LINE__, k);
        seq12 = lkt_seq2LktItem(seq_buf, 0, 11);
//fprintf(stderr, "%u\n", __LINE__);
        l = fm_idx->fastmap->item[seq12+1]-1;
//fprintf(stderr, "%u\n", __LINE__);
        l -= get_12mer_correct(hash_boundry, l);
//fprintf(stderr, "%u, r*num+q = %u\n", __LINE__, r*num+q);
        target_idx[r*num+q][0] = fm_idx->fastmap->item[seq12];
//fprintf(stderr, "%u\n", __LINE__);
        target_idx[r*num+q][1] = l;
//fprintf(stderr, "%u\n", __LINE__);
        target_idx[r*num+q][2] = i;
//fprintf(stderr, "%u\n", __LINE__);
    } 
       
//fprintf(stderr, "%u, i = %u\n", __LINE__, i);
   return i;
}


int init_idx_target(idx_t *fm_idx, uint32_t hash_boundry[], query_t *query,  uint32_t (*target_idx)[2])
{
    uint32_t i, j, k, l,  seq12 = 0;
    //uint8_t *read_seq = query->read_seq;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    int len = query->l_seq;
    int flg = 0;
    for(i = 0; i < 11; ++i) {
        seq12 <<= 2;
        seq12 |= (uint32_t)read_seq[i]&3; 
        if(read_seq[i] > 3) flg = 1;
    }
    for(i = 11; i < len; ++i) {
        seq12 <<= 2;
        seq12 |= (uint32_t)read_seq[i]&3;
        
        if(read_seq[i] > 3) flg = 1;
        seq12 &= 0xFFFFFF; 
        l = fm_idx->fastmap->item[seq12+1]-1;
        l -= get_12mer_correct(hash_boundry, l);
        target_idx[i-11][0] = fm_idx->fastmap->item[seq12];
        target_idx[i-11][1] = l;
        if(flg > 0) {
            target_idx[i-11][0] = (uint32_t)-1;
            target_idx[i-11][1] = 0;
            flg++;
            flg = flg%13;
        }
    }
  
    uint32_t seq12_0 = seq12; 
 
    for(i = len - 11; i < len; ++i) {
        seq12_0 <<= 2;
        seq12_0 &= 0xFFFFFF; 
        l = fm_idx->fastmap->item[seq12_0+1]-1;
        l -= get_12mer_correct(hash_boundry, l);
        target_idx[i][0] = fm_idx->fastmap->item[seq12_0];
    } 
    uint32_t seq12_3 = seq12; 
    for(i = len - 11; i < len; ++i) {
        seq12_3 <<= 2;
        seq12_3 |= (uint32_t)3; 
        seq12_3 &= 0xFFFFFF; 
        l = fm_idx->fastmap->item[seq12_3+1]-1;
        l -= get_12mer_correct(hash_boundry, l);
        target_idx[i][1] = l;
    }      
    //+++++++++++++++++++++++++++++++++++

    int16_t *hash_idx = query->hash_idx[query->is_rev];
    int16_t *hash_pos = query->hash_pos[query->is_rev]; 
    int32_t *hash_bak = query->hash_bak[query->is_rev]; 
    int16_t buf[LEN_READ];
    uint32_t seq10 = 0;
    for(i = 0; i < 9; ++i) {
        seq10 <<= 2;
        seq10 |= (uint32_t)read_seq[i]&3; 
    }
    j = 1;
    for(i = 9; i < len; ++i) {

        seq10 <<= 2;
        seq10 |= (uint32_t)read_seq[i]&3;
        seq10 &= 0xFFFFF;//20 bits 

        if(hash_idx[seq10] == 0) {
            hash_bak[j++] = seq10;  
        }

        hash_idx[seq10]++;

    }
    hash_bak[0] = j - 1;
  
    int tot = 1;
    for(i = 1; i <= hash_bak[0]; ++i) {
        uint16_t cnt = hash_idx[hash_bak[i]];
        hash_idx[hash_bak[i]] = tot; 

        hash_pos[tot] = cnt;
//fprintf(stderr, "hash_pos[tot = %d] = %d\n", tot, cnt);
//fprintf(stderr, "hash_idx[%x] = %d\n", hash_bak[i], tot);
        tot += cnt + 1; 
    }
    seq10 = 0;
    for(i = 0; i < 9; ++i) {
        seq10 <<= 2;
        seq10 |= (uint32_t)read_seq[i]&3; 
    }
//fprintf(stderr, "%u, hash_idx[2] = %d\n", __LINE__, hash_idx[2]);
    for(i = 9; i < len; ++i) {
        seq10 <<= 2;
        seq10 |= (uint32_t)read_seq[i]&3;
        seq10 &= 0xFFFFF;//20 bits 
        j = ++hash_idx[seq10];
        hash_pos[j] = i - 9;
//fprintf(stderr, "hash_pos[j = %d] = %d\n", j, i-9);
    }
//fprintf(stderr, "%u, hash_idx[2] = %d\n", __LINE__, hash_idx[2]);
    for(i = 1; i <= hash_bak[0]; ++i) {
        uint16_t idx = hash_idx[hash_bak[i]];
//fprintf(stderr, "i = %d, hash_idx[%x] = %d, num = %d\n", i, hash_bak[i], idx, hash_pos[idx]);
    
    }



    for(i = hash_bak[0]; i > 1; --i) {
        hash_idx[hash_bak[i]] = hash_idx[hash_bak[i-1]] + 1;
    }   
//fprintf(stderr, "%u, hash_idx[2] = %d\n", __LINE__, hash_idx[2]);
    hash_idx[hash_bak[1]] = 1;

    //-----------------------------------
//fprintf(stderr, "%u, hash_bak[0] = %d\n", __LINE__, hash_bak[0]);
/*  
    for(i = 1; i <= hash_bak[0]; ++i) {
        uint16_t idx = hash_idx[hash_bak[i]];
fprintf(stderr, "i = %d, hash_idx[%x] = %d, num = %d\n", i, hash_bak[i], idx, hash_pos[idx]);
        for(j = 1; j <= hash_pos[idx]; ++j) {
fprintf(stderr, "hash_pos[%d] = %d\n", j, hash_pos[idx+j]);
        }
    
    }
*/
    for(i = 0; i < len; ++i) {
    
fprintf(stderr, "%u, i = %d, target_idx = (%u, %u)\n", __LINE__, i, target_idx[i][0], target_idx[i][1]);
    }
    return len;
}

int aln_idx_map2(idx_t *fm_idx,  uint8_t *read_seq, uint32_t idx[2], uint32_t (*target_idx)[3], int target_num, int L_offset, uint32_t *pos_buf, struct SubBuf *sub)
{
    uint32_t i, j, k, l, pos, pos_i, i_bwt, seq12;
    uint32_t *isa = fm_idx->bwt->isa;
    //---------------------------
    int err_num = sub->err_sum[1];
    err_num = 1;
//fprintf(stderr, "%u\n", __LINE__);
    if(sub->query_err <= 10) err_num = 1;
    else err_num = 2;
    k = idx[0];
    l = idx[1];
    int row = 0, max_row = 0, max_score = 0, score;
//fprintf(stderr, "%u, k= %u, l = %u\n", __LINE__, k, l );
    for(i = k; i <= l; ++i){
        pos = bwt_sa(fm_idx->bwt, i) - L_offset;  
        score = 0;  
//fprintf(stderr, "%u, i = %u, pos = %u\n", __LINE__, i-k, pos);
        //i_bwt = isa[pos];
        //pos = bwt_sa(fm_idx->bwt, i_bwt);  
//fprintf(stderr, "%u, i = %u, pos = %u\n", __LINE__, i-k, pos);
        for(j = 0; j < target_num; ++j){
            pos_i = pos + target_idx[j][2];
            //i_bwt = isa[pos_i];
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
//fprintf(stderr, "%u, j = %u, pos_i = %u, target_idx[2] = %u\n", __LINE__, j, pos_i, target_idx[j][2]);
//fprintf(stderr, "%u, ibwt = %u, target_idx[0] = %u, target_idx[1] = %u\n", __LINE__,i_bwt,  target_idx[j][0], target_idx[j][1]);
            if(i_bwt >= target_idx[j][0] && i_bwt < target_idx[j][1] ){
                score++; 
            }
        }
//fprintf(stderr, "%u, score = %d, err_num = %d, target_num = %d\n", __LINE__, score, err_num, target_num);
        if(score > target_num - err_num) {
            pos_buf[row] = pos;
            if(score > max_score) {
                max_score = score;
                max_row = row; 
            }
            row++;
        }
    } 
    if(max_row > 0) {
        pos_buf[row] = pos_buf[0];
        pos_buf[0] = pos_buf[max_row];
        pos_buf[max_row] = pos_buf[row]; 
    }
//fprintf(stderr, "%u\n", __LINE__);
    return row;
}



int aln_idx_map1(idx_t *fm_idx,  query_t *query, uint32_t idx[], uint32_t (*target_idx)[3], int L_offset, int8_t *score, int err_num[3], struct SubBuf *sub)
{
    //uint8_t *read_seq = query->read_seq;

    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t i, j, l, pos, pos_i, i_bwt, seq12;
    uint32_t *isa = fm_idx->bwt->isa;
    int k;
    //---------------------------
    //int8_t score[LEN_READ/12];
    int bg = query->l_seq%12/2;
    int ed = (query->l_seq-bg)/12;
    int num = query->l_seq/12;
    int q, r;
    err_num[0] = ed;
    err_num[1] = 0;
    err_num[2] = 0;
    int r_flg = 0;
    int l_flg = 0;
    pos = bwt_sa(fm_idx->bwt, idx[0]) - L_offset;
    idx[1] = pos;
    for(i = 0; i < ed; ++i) {
        //i_bwt = isa[pos + i * 12 + bg];
        pos_i = pos+ i*12 +bg;
        i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
        score[i] = 0;
        j = bg * num + i; 
        if(i_bwt >= target_idx[j][0] && i_bwt < target_idx[j][1] ){
            score[i]++;
            err_num[0]--;
            l_flg = 1;
            err_num[2] = 0;
        } else{
            if(l_flg == 0) err_num[1] = i;
            err_num[2]++;  
        }
    }
   

    
    //-------------------------------- 
    //fprintf(stderr, "%u\n", __LINE__);
    return err_num[0];
}
int slc_aln_pos(idx_t *fm_idx, query_t *query, int *seq_off, int pos_num, struct SubBuf *sub)
{
    uint32_t *pos_buf = sub->pos_buf;
    uint32_t (*err_buf)[5] = sub->err_buf;
    int pos_n[3];
    pos_n[0] = pos_num;
    

    if(pos_num == 0) return pos_num;    
    int thres_score, max_score;
    uint32_t max_pos;
    //uint32_t err_out[5];
    uint32_t *err_out = (uint32_t *)(sub->err_buf + pos_num + 1);
    err_out[0] = pos_num;
    if(query->b0 > 0) thres_score = query->b0;
    else thres_score = query->candi_thres;
    max_score = thres_score;

    ////fprintf(stderr, "%u, query->b0 = %u, max_score = %u\n", __LINE__, query->b0, max_score);    
    uint32_t pos;  
    int pos_j = -1;                                     
    uint32_t pos_i = 0;
    int i, j;
    int err_sum = sub->err_sum[0]/5;
    int cur_err;
    for(j=0; j < pos_num; j++){
        

        eval_pos(fm_idx,  query, pos_buf+j,  seq_off, err_out, sub);
        pos = pos_buf[j];
fprintf(stderr, "%u, pos = %u, err_out[3] = %d, j = %d, pos_num = %d\n", __LINE__, pos, err_out[3], j, pos_num);
        if(pos > fm_idx->bwt->seq_len) {
            printf("%u, pos = %u\n", __LINE__, pos); 
            exit(1);
        }
//fprintf(stderr, "%u, err_out[3] = %d, pos = %u\n", __LINE__, err_out[3], pos);
        if(err_out[3] == 0) { 
            continue;
        } else if(err_out[3] == 1){
            pos_buf[pos_i] = pos; 
            memcpy(err_buf+pos_i, err_out, 5*sizeof(uint32_t));
            pos_i++;
        } else if(err_out[3] == 2) {
            pos_buf[pos_i] = pos; 
            memcpy(err_buf+pos_i, err_out, 5*sizeof(uint32_t));
            if(err_out[0] > max_score){
                max_score = err_out[0];
                max_pos = pos;
                pos_j = pos_i;
            }
            pos_i++; 
        } else if(err_out[3] == 3) {
            pos_buf[pos_i] = pos; 
            //err_out[0] += 8; 
            memcpy(err_buf+pos_i, err_out, 5*sizeof(uint32_t));
            pos_i++;
        } else {
            printf("error: %s, %u\n", __func__, __LINE__); 
            exit(1);
        }
    } // end for(j=0; j < pos_num; j++)++++
    pos_n[1] = pos_i;

//fprintf(stderr, "%u, pos_i = %d\n", __LINE__, pos_i); 
    if(pos_i == 0) {
        pos_n[2] = 0;
        return pos_i;
    } else if(pos_i == 1){
        pos_n[2] = 1;
        return pos_i;
    } else {
        pos_n[2] = 0; 
    }

//fprintf(stderr, "%u, pos_i = %d\n", __LINE__, pos_i); 
    if(pos_i > 1 && max_score > thres_score) {
        pos_buf[pos_i] = pos_buf[0];
        pos_buf[0] = pos_buf[pos_j]; 
        pos_buf[pos_j] = pos_buf[pos_i]; 
        for(i = 0; i < 5; ++i)         
            err_buf[pos_i][i] = err_buf[0][i];
        for(i = 0; i < 5; ++i)         
            err_buf[0][i] = err_buf[pos_j][i];
        for(i = 0; i < 5; ++i)         
            err_buf[pos_j][i] = err_buf[pos_i][i];
    }     

//fprintf(stderr, "%u, pos_i = %d\n", __LINE__, pos_i); 
    //++++++++++++++++++++++++++++++++ 
    for(i = 1, j = 2;  j < pos_i;) {
        if(err_buf[i][0] >= max_score){
            ++i;
            j = i +1;
        } else if( err_buf[j][0] < max_score ){
            ++j; 
        } else if(err_buf[j][0] >= max_score ){
            int k;
            pos_buf[pos_i] = pos_buf[i];
            pos_buf[i] = pos_buf[j]; 
            pos_buf[j] = pos_buf[pos_i];

            for(k = 0; k < 5; ++k)         
                err_buf[pos_i][k] = err_buf[i][k];
            for(k = 0; k < 5; ++k)         
                err_buf[i][k] = err_buf[j][k];
            for(k = 0; k < 5; ++k)         
                err_buf[j][k] = err_buf[pos_i][k];
            i++;
            j++;
        } else{
        }
    }
    pos_n[2] = i;

//fprintf(stderr, "%u, pos_i = %d\n", __LINE__, pos_i); 
    return pos_i; 
}

int slc_aln_pos1(idx_t *fm_idx, query_t *query, int *seq_off, int pos_n[], struct SubBuf *sub)
{
    uint32_t *pos_buf = sub->pos_buf;
    int *err_buf = (int *)sub->err_buf;
    int pos_num = pos_n[0];
    int indel_num = 2;
    int indel_ed = indel_num*2 + 1; 
    int err[3] = {16, 16, 16};
    int err_sum = sub->err_sum[0]/5;   
    uint32_t min_pos, min_err,MIN_ERR, err_num[32][3];
    MIN_ERR = (query->l_seq - query->b0)/5 + 1;
    if(MIN_ERR > 12) {
        MIN_ERR = 12;
    }
    min_err = MIN_ERR;

    //fprintf(stderr, "%u, query->b0 = %u, min_err = %u\n", __LINE__, query->b0, min_err);    
    uint32_t pos;  
    int pos_j = -1;                                     
    uint32_t pos_i = 0;
    int i, j;
    for(j=0; j < pos_num; j++){
        aln_pos_map(fm_idx,  query, pos_buf+j,  seq_off, indel_num, err_num, sub);
        //eval_pos(fm_idx,  query, pos_buf+j,  seq_off, err_out, sub);
        
        pos = pos_buf[j];                     
//fprintf(stderr, "%u, pos = %u, indel_ed = %u\n", __LINE__, pos_buf[j], indel_ed);
    
        for(i = 0; i <= indel_ed; ++i) {
            if(err_num[i][0] < err[0]) err[0] = err_num[i][0]; 
            if(err_num[i][1] < err[1]) err[1] = err_num[i][1]; 
            if(err_num[i][2] < err[2]) err[2] = err_num[i][2]; 
        }
//fprintf(stderr, "%u, err[0] = %u\n", __LINE__, err[0]); 
        if(err[0] > err_sum + err[1] + err[2] + 1){
            err[0] = err_sum + err[1] + err[2] + 1;
        }
        if(pos_num < IS_SMLSIZ ) {
            int len = query->l_seq;
            int L = seq_off[0];
            int R = (len - seq_off[1]);
            int err0; 
            if(L*3 > len && err[1] < 2){ 
                err0 = err[1] + 1;
            } else if (R*3 > len && err[2] < 2) {
                err0 = err[2] + 1;
            } 
            if(err0 < err[0]) err[0] = err0; 
        }
        int L = query->l_seq;
        if( seq_off[0]*3 > L && (L - seq_off[1])*3 > L ) {//???
            int err0;
            if(err_num[indel_num*2+1][1] <= err_num[indel_num*2+1][2]){
                err0 = err_num[indel_num*2+1][1];
            } else{
                err0 = err_num[indel_num*2+1][2];
            }
            if(err0 + 2 < err[0]) {
                err[0] = err0 + 2; 
            }
        }  
//fprintf(stderr, "%u, err[0] = %u\n", __LINE__, err[0]); 
        if(err[0] < MIN_ERR +1) {
            pos_buf[pos_i] = pos; 
            err_buf[pos_i] = err[0];
//fprintf(stderr, "%u, pos = %u\n", __LINE__, pos);
            if(err[0] < min_err) {
                min_err = err[0];
                min_pos = pos;
                pos_j = pos_i;
//fprintf(stderr, "%u, min_pos = %u, min_err = %u, pos_j = %u\n", __LINE__, min_pos, min_err,  pos_j);
            }
            pos_i++;
//fprintf(stderr, "%u, pos = %u, pos_i = %u\n", __LINE__, pos, pos_i);
        }
//fprintf(stderr, "%u, pos = %u, pos_i = %u, pos_j = %u\n", __LINE__, pos, pos_i, pos_j);
        
    }
    pos_n[1] = pos_i;
    pos_n[2] = 0;
    if(pos_i > 1 && min_err < MIN_ERR) {
        pos_buf[pos_j] = pos_buf[0];
        err_buf[pos_j] = err_buf[0];
        
        pos_buf[0] = min_pos;
        err_buf[0] = min_err;
        if(pos_i < 2) {
            pos_n[2] = 1;
            return pos_i;
        }
        //++++++++++++++++++++++++++++++++ 
        j = 1;
        for(i = 1; i < pos_i; ++i) {
            if(err_buf[i] <= min_err){
                if(i > j) {
                    pos_buf[pos_i] = pos_buf[j];
                    err_buf[pos_i] = err_buf[j];
                    
                    pos_buf[j] = pos_buf[i];
                    err_buf[j] = err_buf[i];
                    
                    pos_buf[i] = pos_buf[pos_i];
                    err_buf[i] = err_buf[pos_i];
                
                }
                j++;
            }
        }
        pos_n[2] = j;
    }
    return pos_i; 
}

int aln_out( query_t *query, uint32_t pos, int score, int len_L_cut, int len_R_cut, int in_flg,  struct SubBuf *sub)
{
//return 0;
    int delta = sub->delta;
    uint32_t drct_f = query->is_rev*1024;
    aln_out_t *out = sub->aln_out;
    uint32_t (*out_buf)[2] = sub->aln_out->out_buf; 
    int (*found)[4] = sub->aln_out->found; 
    int cur_s;
////fprintf(stderr, "%u, score = %d, found[0][0] = %u\n", __LINE__, score, found[0][0]);   
    //---------------------------------------- 
    if(in_flg == 2) { 
        cur_s = score - found[0][0];
        if(cur_s >= delta) {
            memset(found, 0, 4*delta*sizeof(uint32_t));
        } else if(cur_s < delta) {
            memmove(found+cur_s, found, 4*cur_s*sizeof(uint32_t)); 
            memset(found, 0, 4*cur_s*sizeof(uint32_t));
        } 
        found[0][0] = score;
        found[0][1] = 1;
        found[0][2] = out->len;
        found[0][3] = out->len;
        
        out_buf[out->len][0] = 1;  
        out_buf[out->len][1] = 0; 
        out_buf[out->len+1][0] = pos;  
        out_buf[out->len+1][1] = drct_f + score;  
        out->len += 2;
        out->num++;

        query->b0 = score;
        query->pos = pos;
        query->strand = query->is_rev;
        //query->tlen = query->l_seq - len_L_cut - len_R_cut;//r.te+1-r.tb;
        query->ref_end = query->l_seq - len_R_cut;//r.te+1-r.tb;
        query->ref_start = len_L_cut+MAX_CLIP; //r.tb;
        query->seq_start = len_L_cut; //r.qb;
        query->seq_end = query->l_seq - len_R_cut;//r.qe+1;
        query->n_cigar = 0;
        query->query_err = query->l_seq - query->b0;
        sub->query_err   = query->l_seq - query->b0;

////fprintf(stderr, "%u, query->b0 = %u, pos = %u, tlen = %u, ref_start = %u, seq_start = %u\n", __LINE__, query->b0, query->pos, query->tlen, query->ref_start, query->seq_start);            
    } else if (in_flg == 1){
        int flg = 0, i; 
        cur_s = found[0][0]-score;
        if(found[cur_s][1] == 0){
            found[cur_s][0] = score;
            found[cur_s][1] = 1;
            found[cur_s][2] = out->len;
            found[cur_s][3] = out->len;
            
            out_buf[out->len][0] = 1;  
            out_buf[out->len][1] = 0; 
            out_buf[out->len+1][0] = pos;  
            out_buf[out->len+1][1] = drct_f + score;  
            out->len += 2;
            out->num++;
        } else{
            found[cur_s][0] = score;
            for(i = 0; i < delta; ++i) {
                if( found[i][3] > found[cur_s][3]) {
                    flg = 1;
                    break; 
                }
            }
            if(flg == 0) {
                found[cur_s][1]++;
                out_buf[found[cur_s][3]][0]++;  
                out_buf[out->len][0] = pos;  
                out_buf[out->len][1] = drct_f + score;
                out->len++;
                out->num++;
            } else{
                found[cur_s][1]++;
                out_buf[found[cur_s][3]][1] = out->len; 
                found[cur_s][3] = out->len;
                out_buf[out->len][0] = 1;  
                out_buf[out->len][1] = 0; 
                out_buf[out->len+1][0] = pos;  
                out_buf[out->len+1][1] = drct_f + score;  
                out->len += 2;
                out->num++;
            }
        } 
    }
    return 0;
}
int calcu_aln_score(idx_t *fm_idx,  query_t *query, uint32_t pos_buf[],  int seq_off[],  uint32_t err_out[], struct SubBuf *sub)
{
    int max_err = query->l_seq - query->b0 + sub->delta;
    if(max_err > query->error_thres) max_err = query->error_thres;  
    int sub_flg   = sub->calcu_flg[0];
    int indel_flg = sub->calcu_flg[1];
    int bg_flg    = sub->calcu_flg[2];

    uint32_t *NT_sum = sub->NT_sum; 
    uint32_t (*target_idx)[2] = query->target_idx; 
    //uint8_t *read_seq = query->read_seq;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t pos, pos_i, i_bwt,l_flg[33], seq12;
    uint32_t *isa = fm_idx->bwt->isa;
    int i, j, k;
    int l_len = seq_off[0];
    int r_len = query->l_seq - seq_off[1];
    int l_num = l_len/12;
    int r_num = r_len/12;
    int l_bg = l_len-l_num*12;
    int r_bg = seq_off[1];
    int l_sum[LEN_READ/12], r_sum[LEN_READ/12], 
        l_score[LEN_READ/12], r_score[LEN_READ/12]; 
    int cur_score, l_max_score, r_max_score, l_max_i, r_max_i;
    pos = pos_buf[0];
    //+++++++++++++++++++++++++++++++++++++++++++
    if(err_out[3] > 0) {
        for(i = 0; i < l_num; ++i){
            if( (err_out[4] & (uint32_t)1 << (31-i)) > 0){
                l_score[i] = 1;  
            } else {
                l_score[i] = 0; 
            }
        } 
        l_sum[0] = l_score[0];
        for(i = 1; i < l_num; ++i) {
            l_sum[i] = l_sum[i-1] + l_score[i]; 
        }
        for(i = 0; i < r_num; ++i){
            if( (err_out[4] & (uint32_t)1 << (15-i)) > 0){
                r_score[i] = 1;  
                r_sum[i] = r_sum[i-1] + 1; 
            } else {
                r_score[i] = 0; 
            }
        }
        r_sum[0] = r_score[0];
        for(i = 1; i < r_num; ++i) {
            r_sum[i] = r_sum[i-1] + r_score[i]; 
        }
        if(err_out[1] == 1) {
            l_max_i = err_out[2]%256; 
        } else if(err_out[1] == 2){
            r_max_i = err_out[2]%256;
        }
    } else{
        printf("%s, %u\n", __func__, __LINE__); 
        exit(1);
    } 
    uint8_t ch;
    uint32_t bg, ed;
    int err_num = 0;
    int len_L_cut = 0, len_R_cut = 0;
    if( (sub_flg > 0 && err_out[1] == 0) || (indel_flg > 0 && err_out[1] == 2)) {
        for(i = 0; i < l_num; ++i){
            if(l_score[i] == 0) {
                for(j = 0; j < 12; ++j){
                    k = l_bg + 12*i + j; 
                    ch = read_seq[k];
                    bg = NT_sum[ch];
                    ed = NT_sum[ch+1]; 
                    pos_i = pos + k; 
                    i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
                    if(i_bwt < bg || i_bwt >= ed ){
                        err_num += 5;  
                    } 
                }
            } 
        }
        for(i = 5; i < l_bg; ++i) {
            ch = read_seq[i];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos + i; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                err_num += 5;  
            } 
        }
        for(i = 4; i >= 0; --i) {
            ch = read_seq[i];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos + i; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                err_num += i+1;
                len_L_cut = i; 
                break;
            } 
        }
        //++++++++++++++++++++++++++++++++++
    }
    if( (sub_flg > 0 && err_out[1] == 0) || (indel_flg > 0 && err_out[1] == 1)) {
        for(i = 0; i < r_num; ++i){
            if(r_score[i] == 0) {
                for(j = 0; j < 12; ++j){
                    k = r_bg + 12*i + j; 
                    ch = read_seq[k];
                    bg = NT_sum[ch];
                    ed = NT_sum[ch+1]; 
                    pos_i = pos + k;
                    i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
                    if(i_bwt < bg || i_bwt >= ed ){
                        err_num += 5;  
                    } 
                }
            } 
        }
        for(i = r_num*12; i < r_len-5; ++i) {
            k = r_bg + i;
            ch = read_seq[k];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos + k; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                err_num += 5;  
            } 
        }
        for(i = r_len-5; i < r_len; ++i) {
            k = r_bg + i;
            ch = read_seq[k];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos + k;
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                err_num += r_len -i;
                len_R_cut = r_len - i;
                break;  
            } 
        }
    } 
    //+++++++++++++++++++++++++++++
    if(sub_flg > 0 && err_out[1] == 0) {
        int old_score = err_out[0];
        int sub_score = query->l_seq - err_num - sub->err_sum[0];

        for(i = 0; i < 4; ++i) err_out[i] = 0; 
        if(sub_score < old_score) {
            err_out[3] = 3;
        }
        err_out[0] = sub_score;
        if(query->b0 < sub_score) {
            err_out[2] = 2;
            aln_out(query, pos, sub_score, len_L_cut, len_R_cut, 2, sub);
        } else  if(query->b0 - sub->delta < sub_score){
            err_out[2] = 1;
            aln_out(query, pos, sub_score, len_L_cut, len_R_cut, 1, sub);
        } else{
            err_out[2] = 0;
        } 
        return l_sum[l_num-1] + r_sum[r_num-1];   
    }
    //-----------------------------
    if(indel_flg == 0) {
        return l_sum[l_num-1] + r_sum[r_num-1];   
    }
    int ins_score, del_score, indel_i, len_indel;
    if(bg_flg > 0) { 
        len_indel = bg_flg; 
    } else if(indel_flg > 0){
        len_indel = indel_flg;
    } else {
        len_indel = 0; 
    }
    int ins_num, del_num;
    int ins_buf[LEN_READ/12], del_buf[LEN_READ/12];

    if(indel_flg > 0 && err_out[1] == 1) {//左侧插入删除
        //发生插入删除区间
        int err_num_R = err_num;
        int r_ed;
        int l_ed;
        for(j = 11; j >= 0; j--){
            i = l_max_i*12 + l_bg + j; 
            ch = read_seq[i];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos + i; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                r_ed = j;
                break;
            }        
        } 
        indel_i = err_out[2];
        indel_i >>= 8;
        indel_i -= 32;
        uint32_t pos_0 = pos + indel_i; 
        for(k = 0; k < 12; ++k){
            i = l_max_i*12 + l_bg + k;
            ch = read_seq[i];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos_0 + i;
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                l_ed = k;
                break;
            }        
        }
        int len;
        int l;
        if(indel_i > 0) l = indel_i;
        else l = -indel_i;
        if(indel_i > 0) {
            len = r_ed -l_ed + 1;
        } else { 
            len = r_ed - l_ed -l +1;
        }

        if(len <0) len = 0; 
        if(indel_i > 0) {
            err_num += 6 + l*2 + 5*len;  
        } else{
            err_num += 6 + l + 5*len; 
        }
////fprintf(stderr, "%u, indel_i = %u, l_ed = %u, r_ed = %u, len = %u, l = %u, err_num = %u\n", __LINE__, indel_i, l_ed, r_ed, len, l, err_num);
        //插入删除的右侧
        int *indel_buf;
        if(indel_i > 0) indel_buf = ins_buf;
        else indel_buf = del_buf;
        for(i = l_max_i+1; i < l_num; ++i){
            if(l_score[i] == 0) {
                for(j = 0; j < 12; ++j){
                    k = l_bg + 12*i + j; 
                    ch = read_seq[k];
                    bg = NT_sum[ch];
                    ed = NT_sum[ch+1]; 
                    pos_i = pos_0 + k; 
                   
                    i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
                    if(i_bwt < bg || i_bwt >= ed ){
                        err_num += 5;  
                    } 
                }
            } 
        }
        //插入删除的左侧
        int l_err = 0;
        for(i = 0; i < l_max_i; ++i) {
            j = l_bg + 12*i; 
            pos_i = pos_0 + j; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt >= target_idx[j][0] && i_bwt < target_idx[j][1] ){
                l_score[i] = 1;  
            } else{
                l_score[i] = 0;
                l_err++; 
            }
        }
        if(l_err >0) {
            for(i = 0; i < l_max_i; ++i){
                if(l_score[i] == 0) {
                    for(j = 0; j < 12; ++j){
                        k = l_bg + 12*i + j; 
                        ch = read_seq[k];
                        bg = NT_sum[ch];
                        ed = NT_sum[ch+1]; 
                        pos_i = pos_0 + k; 
                        i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
                        if(i_bwt < bg || i_bwt >= ed ){
                            err_num += 5;  
                        } 
                    }
                } 
            }
        }
        //左侧边界 
        for(i = 5; i < l_bg; ++i) {
            ch = read_seq[i];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos_0 + i; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                err_num += 5;  
            } 
        }
        for(i = 4; i >= 0; --i) {
            ch = read_seq[i];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos_0 + i; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                err_num += i+1;
                len_L_cut = i; 
                break;  
            } 
        }
        int l_indel_score = query->l_seq - err_num - sub->err_sum[0];
        err_out[0] = l_indel_score;
        err_out[1] = 1;
        err_out[3] = 0;
        if(query->b0 < l_indel_score) {
            err_out[2] = 2;
            aln_out(query, pos, l_indel_score, len_L_cut, len_R_cut, 2, sub);
        } else  if(query->b0 - sub->delta < l_indel_score){
            err_out[2] = 1;
            aln_out(query, pos, l_indel_score, len_L_cut, len_R_cut, 1, sub);
        } else{
            err_out[2] = 0;
            l_indel_score = query->l_seq - err_num_R - sub->err_sum[0] - 7 - 5;
            if(query->b0 - sub->delta < l_indel_score) {
                err_out[3] = 1; 
            }         
        } 
        return l_sum[l_num-1] + r_sum[r_num-1]; 
    }
    if(indel_flg > 0 && err_out[1] == 2 ) {
        //发生插入删除区间
        int err_num_L = err_num;
        int r_ed;
        int l_ed;
        int p = r_max_i*12 + r_bg;        
        for(j = 0; j < 12; j++){
            i = p + j; 
            ch = read_seq[i];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos + i; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                l_ed = j;
                break;
            }        
        } 
        indel_i = err_out[2];
        indel_i >>= 8;
        indel_i -= 32;
        uint32_t pos_0 = pos + indel_i; 
        for(k = 11; k >= 0; --k){
            i = p + k;
////fprintf(stderr, "%u, i = %u, r_max_i = %d, r_bg = %d, k = %d\n", __LINE__, i, r_max_i, r_bg, k);
            ch = read_seq[i];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos_0 + i;
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                r_ed = k;
                break;
            }        
        }
        int len, l;
        if(indel_i > 0) l = indel_i;
        else l = -indel_i;
        if(indel_i > 0) {
            len = r_ed -l_ed + 1;
        } else { 
            len = r_ed - l_ed -l +1;
        }
        if(len <0) len = 0; 
        if(indel_i > 0) {
            err_num += 6 + indel_i + 5*len;  
        } else{
            err_num += 6 + indel_i*2 + 5*len; 
        }
////fprintf(stderr, "%u, indel_i = %u, l_ed = %u, r_ed = %u, len = %u, l = %u, err_num = %u\n", __LINE__, indel_i, l_ed, r_ed, len, l, err_num);
        //插入删除的右侧
        int *indel_buf;
        if(indel_i > 0) indel_buf = ins_buf;
        else indel_buf = del_buf;
        for(i = r_max_i +1; i < r_num; ++i){
////fprintf(stderr, "%u, i = %u, r_max_i = %u, r_num = %u\n", __LINE__, i, r_max_i, r_num);
            if(indel_buf[i] == 0) {
                for(j = 0; j < 12; ++j){
                    k = r_bg + 12*i + j; 
                    ch = read_seq[k];
                    bg = NT_sum[ch];
                    ed = NT_sum[ch+1]; 
                    pos_i = pos_0 + k; 
                    i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
                    if(i_bwt < bg || i_bwt >= ed ){
                        err_num += 5;  
                    } 
                }
            } 
        }
        //插入删除的左侧
        for(i = 0; i < r_max_i; ++i){
            if(r_score[i] == 0) {
                for(j = 0; j < 12; ++j){
                    k = r_bg + 12*i + j; 
                    ch = read_seq[k];
                    bg = NT_sum[ch];
                    ed = NT_sum[ch+1]; 
                    pos_i = pos_0 + k; 
                    i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
                    if(i_bwt < bg || i_bwt >= ed ){
                        err_num += 5;  
                    } 
                }
            } 
        }
        //右侧边界
        for(i = r_num*12; i < r_len-5; ++i) {
            k = r_bg + i;
            ch = read_seq[k];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos_0 + k; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                err_num += 5;  
            } 
        }
        for(i = r_len-5; i < r_len; ++i) {
            k = r_bg +i;
            ch = read_seq[k];
            bg = NT_sum[ch];
            ed = NT_sum[ch+1]; 
            pos_i = pos_0 + k; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt < bg || i_bwt >= ed ){
                err_num += r_len -i;
                len_R_cut = r_len - i;
                break;  
            } 
        }
        int r_err = 0;
        for(i = r_max_i+1; i < r_num; ++i) {
            j = r_bg + 12*i; 
            pos_i = pos_0 + j; 
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
            if(i_bwt >= target_idx[j][0] && i_bwt < target_idx[j][1] ){
                r_score[i] = 1;  
            } else{
                r_score[i] = 0;
                r_err++;  
            }
        }
        if(r_err > 0) {
            for(i = r_max_i+1; i < r_num; ++i) {
                if(r_score[i] == 0) {
                    for(j = 0; j < 12; ++j){
                        k = r_bg + 12*i + j; 
                        ch = read_seq[k];
                        bg = NT_sum[ch];
                        ed = NT_sum[ch+1]; 
                        pos_i = pos_0 + k; 
                        i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
                        if(i_bwt < bg || i_bwt >= ed ){
                            err_num += 5;  
                        } 
                    }
                } 
            }

        } 
        int r_indel_score = query->l_seq - err_num - sub->err_sum[0];
////fprintf(stderr, "%u, b0 = %d, err_num = %u, r_indel_score = %u, err_sum[0] = %u\n", __LINE__, query->b0, err_num, r_indel_score, sub->err_sum[0]);
        err_out[0] = r_indel_score;
        err_out[1] = 1;
        err_out[3] = 0;
        if(query->b0 < r_indel_score) {
            err_out[2] = 2;
            aln_out(query, pos, r_indel_score, len_L_cut, len_R_cut, 2, sub);
        } else  if(query->b0 - sub->delta < r_indel_score){
            err_out[2] = 1;
            aln_out(query, pos, r_indel_score, len_L_cut, len_R_cut, 1, sub);
        } else{
            err_out[2] = 0;
            r_indel_score = query->l_seq - err_num_L - sub->err_sum[0] - 7 - 5;
            if(query->b0 - sub->delta < r_indel_score) {
                err_out[3] = 1; 
            }         
        } 
        return l_sum[l_num-1] + r_sum[r_num-1]; 
    }
//fprintf(stderr, "%u\n", __LINE__);
printf("error: %u %s\n", __LINE__, __func__);
exit(1);
    return l_sum[l_num-1] + r_sum[r_num-1];   
}

int eval_12mer_olp(idx_t *fm_idx,  query_t *query, seed_t *seed, uint32_t pos,  int seq_off[],   struct SubBuf *sub)
{
    int rev = query->is_rev;
    uint32_t (*tidx)[2]; 
    if(rev == 0) tidx = query->target_idx_f; 
    else tidx = query->target_idx_r; 
    uint8_t *q_flg = query->q_flg;
    int l_seq = query->l_seq;
    int l_off = seq_off[0];
    int r_off = seq_off[1];
    int l_kmer[LEN_READ], r_kmer[LEN_READ];
    int l_num = 0, r_num = 0;
    int q_val = 0;
    int i, j, k;
    int len = 0;
    j = 1, k = 1;
    for(i = 0; i < l_seq; ++i){
        if(rev == 0) {
            q_val = q_flg[i]/16;
        } else {
            q_val = q_flg[l_seq-i-1]%16;
        }
        if(q_val > 0){
            len = 0;
            continue; 
        } else {
            len++;
            if(len == 12) {
                if(i < l_off) { 
                    l_kmer[j++] = i - len+1;  
                } else if(i > r_off){
                    r_kmer[k++] = i - len+1;  
                }
                len = 0;  
            }                 
           
        } 
    }
    l_kmer[0] = j;
    r_kmer[0] = k;

fprintf(stderr, "%u, num_lkmer = %d, num_rkmer = %d, name = %s\n", __LINE__, l_kmer[0], r_kmer[0], query->name); 

    int sd_len;
    uint32_t idx;
    int find_flg, l_score_sum = 0, r_score_sum = 0;
    uint32_t find_i, pos_idx; 
    int kmer_i;
    int ins_w = seq_off[2];
    for(i = 1; i < l_kmer[0]; ++i) {
        kmer_i = l_kmer[i];
        pos_idx = pos + kmer_i;
        for(j = 0; j < ins_w; ++j){
            k = kmer_i - ins_w/2 + j; 
            if(k == kmer_i || k < 0 || k > l_seq - 12) continue; 
            if(idx >= tidx[k][0] && idx < tidx[k][1] ){
                l_num++;
                break; 
            } 
        }
    }
    for(i = 1; i < r_kmer[0]; ++i) {
        kmer_i = r_kmer[i];
        pos_idx = pos + kmer_i;
        for(j = 0; j < ins_w; ++j){
            k = kmer_i - ins_w/2 + j; 
            if(k == kmer_i || k < 0 || k > l_seq - 12) continue; 
            if(idx >= tidx[k][0] && idx < tidx[k][1] ){
                r_num++;
                break; 
            } 
        }
    }
//return 0;
    return l_num + r_num;   
}

int eval_seed_olp(idx_t *fm_idx,  query_t *query, seed_t *seed, uint32_t pos,  int seq_off[],   struct SubBuf *sub)
{
    //uint32_t pos_num = err_out[0];  
    //uint8_t *read_seq = query->read_seq;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t pos_i, i_bwt, seq12;
    int i, j, k;
    int l_len = seq_off[0];
    //int r_len = query->l_seq - seq_off[1];
    int r_len = seq_off[1];
    int cur_score, l_max_score, r_max_score, l_max_i, r_max_i;
    int l_seed_buf[20], r_seed_buf[20];

//++++++++++++++++++++++++++
    seed_t *sd;
    int s_id;
    int bg_seed, ed_seed;
    if(query->is_rev == 0) {
        bg_seed = 0;
        ed_seed = query->seed_num; 
    } else{
        bg_seed = query->seed_num;
        ed_seed = query->seed_num*2; 
    }
    j = 1;
    k = 1;
    for(i = bg_seed; i < ed_seed; ++i){
        if(seed->slc[i].s_off < l_len){
            l_seed_buf[j++] = i;  
        }
        if(seed->slc[i].s_off + SEED_LEN > r_len){
            r_seed_buf[k++] = i;  
        } 
    }

    l_seed_buf[0] = j;
    r_seed_buf[0] = k;
    int sd_len, mid;
    uint32_t idx;
    int find_flg, l_score_sum = 0, r_score_sum = 0;
    uint32_t find_i, pos_idx; 
    for(i = 1; i < l_seed_buf[0]; ++i) {
        
        s_id = l_seed_buf[i];
        sd_len = seed->slc[s_id].len;
        if(sd_len  == 0) continue;
        mid = sd_len/2 + 1;
        int pos_off = seed->slc[s_id].s_off + SEED_LEN - 12 -mid;
        pos_idx = pos + pos_off;
        idx = bwt_get_idx(fm_idx->bwt, pos_idx); 
        fprintf(stderr, "%u, sid = %d, pos_off = %d, idx = %u\n", __LINE__, s_id, pos_off, idx);
        find_flg = -1;
        if(idx >= seed->idx_buf[s_id][mid][0] && idx <= seed->idx_buf[s_id][mid][1]) {
            l_score_sum++;
            find_flg = 0;
            find_i = mid;
            continue;            
        } else{
        
        }
        for(j = 1; j < mid; ++j){
            k = mid - j; 
            if(idx >= seed->idx_buf[s_id][k][0] && idx <= seed->idx_buf[s_id][k][1]) {
                l_score_sum++;
                find_i = k;
                find_flg = 1;
                break; 
            } 

            k = mid + j; 
            if(idx >= seed->idx_buf[s_id][k][0] && idx <= seed->idx_buf[s_id][k][1]) {
                l_score_sum++;
                find_i = k;
                find_flg = 1;
                break; 
            }
        }
    }
//--------------------------
    for(i = 1; i < r_seed_buf[0]; ++i) {
        s_id = r_seed_buf[i];
        sd_len = seed->slc[s_id].len;
        mid = sd_len/2 + 1;
        int pos_off = seed->slc[s_id].s_off + SEED_LEN - 12 -mid;
        pos_idx = pos + pos_off;
        idx = bwt_get_idx(fm_idx->bwt, pos_idx); 
        find_flg = -1; 
        if(idx >= seed->idx_buf[s_id][mid][0] && idx <= seed->idx_buf[s_id][mid][1]) {
            r_score_sum++;
            find_flg = 0;
            find_i = mid;
            continue;            
        } else{
        
        }
        for(j = 1; j < mid; ++j){
            k = mid - j; 
            if(idx >= seed->idx_buf[s_id][k][0] && idx <= seed->idx_buf[s_id][k][1]) {
                r_score_sum++;
                find_i = k;
                find_flg = 1;
                break; 
            } 
            k = mid + j; 
            if(idx >= seed->idx_buf[s_id][k][0] && idx <= seed->idx_buf[s_id][k][1]) {
                r_score_sum++;
                find_i = k;
                find_flg = 1;
                break; 
            }
        }
    }
    return l_score_sum + r_score_sum;   
}


int eval_pos_sub(idx_t *fm_idx,  query_t *query, int pos_n, uint32_t pos_buf[],  int seq_off[], struct SubBuf *sub)
{

    int max_mov = MAX_CLIP; 
    tgt_arry_t *tgt = sub->tgt_arry;
    int min_err_n = LEN_READ;
    uint32_t min_pos;
    uint8_t *min_tgt; 
    uint32_t *min_err_pos;
    int pi;
    int i, j, k;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    for(pi = 0; pi < pos_n; ++pi){
        uint32_t pos = pos_buf[pi]; 
        uint8_t *target = tgt->tgt[pi].target; 
        tgt->tgt[pi].seq_off[0] = seq_off[0]; 
        tgt->tgt[pi].seq_off[1] = seq_off[1]; 
        int bg, ed;
        if(pos < max_mov) { bg = max_mov-pos; }
        else { bg = 0;}
        for(i = bg; i < max_mov; ++i){
            target[i] = __get_pac(fm_idx->pac, pos-i);
        } 

        if(pos >= max_mov) { bg = 0; }
        else { bg = max_mov - pos;}
        ed = max_mov;

        uint8_t ch;
        for(i = bg; i < ed; ++i) {
            target[i] = __get_pac(fm_idx->pac, pos + i - max_mov);
        }
        uint32_t *err_pos = tgt->tgt[pi].err_pos;
        int cur_err_n = 0; 
        for(i = 0; i < LEN_READ; ++i) { 
            ch = __get_pac(fm_idx->pac, pos+i);
            target[i+max_mov] = ch;
            if(ch != read_seq[i]) {
                err_pos[cur_err_n++] = i; 
            }
        }
        tgt->tgt[pi].err_n = cur_err_n;
        if(cur_err_n < min_err_n) {
            min_err_n = cur_err_n;
            min_pos = pos;
            min_tgt = target; 
        }
        bg = query->l_seq;
        ed = bg + max_mov;
        if(pos + ed > fm_idx->bns->l_pac ) ed = fm_idx->bns->l_pac - pos;
        for(i = bg; i < ed; ++i){
            ch = __get_pac(fm_idx->pac, pos+i);
            target[i+max_mov] = ch;

        }
    }
    
    return min_err_n;
}
int aln_pos_both(idx_t *fm_idx,  query_t *query, uint32_t pos_buf[],  int seq_off[], struct SubBuf *sub)
{
fprintf(stderr, "%u, seq_off[0] = %d, seq_off[1] = %d\n", __LINE__, seq_off[0], seq_off[1]);
    int sc_L = 0, sc_R = 0;
    aln_out_t *out = sub->aln_out; 
    query->query_err = query->l_seq - query->b0;
    int delta = sub->delta;
    int l_num = seq_off[0], r_num = query->l_seq - seq_off[1];
    if(l_num >= r_num) {
        sc_L = aln_pos_L(fm_idx, query, pos_buf, seq_off, sub);
        if(query->query_err != query->l_seq - query->b0) {
            printf("query->query_err != query->l_seq - query->b0!\n");
            exit(1); 
        }
        if(sc_L < 1) return 0;
        if(sc_L +delta > query->b0) {
            query->query_err -= query->l_seq - sc_L; 
        } else return 0;
        sc_R = aln_pos_R(fm_idx, query, pos_buf, seq_off,  sub);           
        if(sc_R < 1) return 0;
    } else {
        sc_R = aln_pos_R(fm_idx, query, pos_buf, seq_off, sub);
        if(query->query_err != query->l_seq - query->b0) {
            printf("query->query_err != query->l_seq - query->b0!\n");
            exit(1); 
        }
        if(sc_R < 1) return 0;
        if(sc_R +delta > query->b0){ 
            query->query_err -= query->l_seq - sc_R; 
        } else return 0;
        sc_L = aln_pos_L(fm_idx, query, pos_buf, seq_off,  sub);
        if(sc_L < 1) return 0;
    }
  
    fprintf(stderr, "%u, %s end\n", __LINE__, __func__);
    fprintf(stderr, "L_sc = %d, L_ti = %d, L_qi = %d\n", 
                        out->L_sc, out->L_ti, out->L_qi);
    fprintf(stderr, "R_sc = %d, R_ti = %d, R_qi = %d\n", 
                        out->R_sc, out->R_ti, out->R_qi);
    fprintf(stderr, "score = %d, pos = %u\n", 
                    sc_L + sc_R - query->l_seq - sub->err_sum[0], pos_buf[0]);
    
    return (sc_L + sc_R - query->l_seq - sub->err_sum[0]);
}
/*  
int eval_pos_0(idx_t *fm_idx,  query_t *query, uint32_t pos_buf[],  int seq_off[], struct SubBuf *sub)
{
    int sc_L = 0, sc_R = 0;
    aln_out_t *out = sub->aln_out; 
    query->query_err = query->l_seq - query->b0;
    int delta = sub->delta;
    int l_num = seq_off[0], r_num = query->l_seq - seq_off[1];
    if(l_num >= r_num) {
        sc_L = eval_pos_L(fm_idx, query, pos_buf, seq_off, sub);
        if(query->query_err != query->l_seq - query->b0) {
            printf("query->query_err != query->l_seq - query->b0!\n");
            exit(1); 
        }
        if(sc_L < 1) return 0;
        if(sc_L +delta > query->b0) {
            query->query_err -= query->l_seq - sc_L; 
        } else return 0;
        sc_R = eval_pos_R(fm_idx, query, pos_buf, seq_off,  sub);           
        if(sc_R < 1) return 0;
    } else {
        sc_R = eval_pos_R(fm_idx, query, pos_buf, seq_off, sub);
        if(query->query_err != query->l_seq - query->b0) {
            printf("query->query_err != query->l_seq - query->b0!\n");
            exit(1); 
        }
        if(sc_R < 1) return 0;
        if(sc_R +delta > query->b0){ 
            query->query_err -= query->l_seq - sc_R; 
        } else return 0;
        sc_L = eval_pos_L(fm_idx, query, pos_buf, seq_off,  sub);
        if(sc_L < 1) return 0;
    }
  
    fprintf(stderr, "%u, %s end\n", __LINE__, __func__);
    fprintf(stderr, "L_sc = %d, L_ti = %d, L_qi = %d\n", 
                        out->L_sc, out->L_ti, out->L_qi);
    fprintf(stderr, "R_sc = %d, R_ti = %d, R_qi = %d\n", 
                        out->R_sc, out->R_ti, out->R_qi);
    fprintf(stderr, "score = %d, pos = %u\n", 
                    sc_L + sc_R - query->l_seq - sub->err_sum[0], pos_buf[0]);
    
    return (sc_L + sc_R - query->l_seq - sub->err_sum[0]);
}
*/

int aln_pos_R_0(idx_t *fm_idx,  query_t *query, uint32_t pos_buf[],  int seq_off[], struct SubBuf *sub)
{
    int max_mov = MAX_CLIP; 
    int i, j, k;
    aln_out_t *out = sub->aln_out; 
    uint32_t pos, pos_i;
    pos = pos_buf[0];
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint8_t *target_e = query->target + max_mov;     
    uint8_t *target = target_e + max_mov;// 20 ????
    int *err_pos     = sub->eval_pos;
    int l_seq = query->l_seq;
    int seq_sc = -1;
    int best_sc = 0, best_i = -1;
    int tar_mov = 0;
    int R_bg = seq_off[1] -1 - max_mov; 
    for(j = R_bg; j < l_seq + 2*max_mov; ++j){ 
        pos_i = pos - max_mov + j;// 12 or 24 ???? 
        target_e[j] = __get_pac(fm_idx->pac, pos_i);
    }
    R_bg += max_mov; 
    int sc_flg = 0;
    int cur_sc = 0;
    int pd;
    int err_n = 0, err_n0 = 0; 
    err_pos[0] = -1; 
    err_n0 = 0;
    for(i = R_bg; i < l_seq; ++i) {
        if(target[i] != read_seq[i])  {
            err_pos[0] = i;
            err_n0++;
            break; 
        } 
    }
    if(err_pos[0] != -1) seq_sc = err_pos[0];
    else seq_sc = l_seq;
    best_sc = seq_sc;
    best_i = err_pos[0] - 1;
    pd = err_pos[0];
    if(err_n0 == 0) {
        out->R_ti = max_mov + tar_mov + l_seq;
        out->R_qi = l_seq;
        out->R_sc = best_sc;
        return best_sc; 
    }
    if(l_seq - pd <= sub->sub_err) {
        out->R_ti = max_mov + tar_mov + l_seq;
        out->R_qi = l_seq;
        out->R_sc = best_sc; 
        return best_sc; 
    }
    
    k = pd + 1;

    err_n = err_n0;
    cur_sc = -5; 
/*  
    int ct_sc = 0;
    int ct_i = k;
    int sm_sc = 0;
    int sm_i = k;
    int sm_n = 0; 
*/
    while(k < l_seq){ 
        if(read_seq[k] < 4) { //N = 4
            if(target[k] == read_seq[k]) {
                cur_sc++; 
            } else{
                cur_sc -= 4;
                err_pos[err_n++] = k; 
            } 
        } else{ // N得分
            cur_sc -= 1;
            err_pos[err_n++] = k; 
        }
/*      
        if(sm_sc > cur_sc) {
            sm_sc = cur_sc;
            sm_i = k;
            sm_n = err_n; 
        } 
        if(ct_sc < cur_sc) {
            ct_sc = cur_sc;
            ct_i = k; 
        } 
*/

        if( err_n > 1) {// 12 or 24 ???
            ++k; 
            break; 
        }
        ++k;
    } //  end while(k >= 0) ++++++++++
    
         
    pd = k - 1; 
    //评估当前区间插入删除比对的可能性
    if(err_n == 0) {
        printf("%u, error : err_n = %d, pd = %d!!!\n", __LINE__, err_n, pd);    
        exit(1);
    }
    if(err_n <= 2) {
        //cur_sc = l_seq - err_pos[0] - 5;
        if(read_seq[err_pos[0]] < 4) {
            cur_sc = l_seq - err_pos[0] - 5;
        } else {
            cur_sc = l_seq - err_pos[0] - 2;
        }           
        
        seq_sc += cur_sc;
        best_i = l_seq; 
    } // end if(err_n == 1) +++
    if(err_n == 1) {
fprintf(stderr, "%u, cur_sc = %d, seq_sc = %d, err_pos[1] = %d\n", __LINE__, cur_sc, seq_sc, err_pos[1]);
        best_sc = seq_sc;
        out->R_ti = max_mov + tar_mov + best_i;
        out->R_qi = best_i;
        out->R_sc = best_sc;
        return best_sc; 
    }
    if(err_n == 2 && l_seq - err_pos[1] -1 <= sub->sub_err) {
        cur_sc = -(l_seq - err_pos[1] -1); 
fprintf(stderr, "%u, cur_sc = %d, seq_sc = %d, err_pos[1] = %d\n", __LINE__, cur_sc, seq_sc, err_pos[1]);
        seq_sc += cur_sc;
        best_i = err_pos[1]; 
        best_sc = seq_sc;
        out->L_ti = max_mov + tar_mov + best_i;
        out->L_qi = best_i;
        out->L_sc = best_sc; 
        return best_sc;
    } 
    return 0;
}

int aln_pos_R(idx_t *fm_idx,  query_t *query, uint32_t pos_buf[],  int seq_off[], struct SubBuf *sub)
{
    int max_mov = MAX_CLIP, len_km = 12; 
    int err = 0;
    int max_err = query->l_seq - query->b0 + sub->delta;
    int i, j, k, ret = 0;
    aln_out_t *out = sub->aln_out; 
    uint32_t pos, pos_i;
    pos = pos_buf[0];
    uint32_t (*target_idx)[2] = query->target_idx; 
//uint8_t *read_seq = query->read_seq;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint8_t *target_e = query->target + max_mov;     
    uint8_t *target = target_e + max_mov;// 20 ????
    kswst_t *kswst = &sub->kswst;
    int *err_pos     = sub->eval_pos;
    int *err_pos0    = sub->eval_pos + LEN_READ;
    int *err_pos_buf = sub->eval_pos + 2*LEN_READ; 

    int qlen = query->l_seq;
    int read_len = query->l_seq + max_mov;
    int l_len = seq_off[0] + len_km - 1;
    int l_num = l_len/len_km;
    int l_bg = l_len-l_num*len_km;
    int r_len = query->l_seq - seq_off[1] + len_km -1;
    int r_num = r_len/len_km;
    int r_bg = seq_off[1] - 1;

    int *l_sc  = sub->eval_sc;  
    int *l_sum = sub->eval_sc + LEN_READ/12;
    int *r_sc  = sub->eval_sc + 2*(LEN_READ/12); 
    int *r_sum = sub->eval_sc + 3*(LEN_READ/12); 

    int l_max_sc, r_max_sc, l_max_i, r_max_i;
    int L_ed = l_bg;
    int R_bg = r_bg;
   
    int seq_sc = -1;
    int best_sc = 0, best_i = -1;
    int tar_mov = 0;
    R_bg = seq_off[1] -1 - max_mov; 
    for(j = R_bg; j < read_len + max_mov; ++j){ 
        pos_i = pos - max_mov + j;// 12 or 24 ???? 
        if(pos_i > fm_idx->bns->l_pac ) {
            printf("%u, pos = %u, out of rer len　%u\n", __LINE__, pos, fm_idx->bns->l_pac);
            exit(1);
        }
        target_e[j] = __get_pac(fm_idx->pac, pos_i);
    }
    R_bg += max_mov; 
    int sc_flg = 0;
    int cur_sc = 0, cut_sc = 0;
    int err_num = 0;
    int err_num0 = 0, rep;
    int pd;
    int pe;
    int err_n = 0, err_n0 = 0, err_n1 = 0, N_sum1 = 0, N_sum = 0, pos1 = 0; 
    int max_gap = max_mov;
    int len;
    int gap_o = 6, gap_e = 1;
    
    err_pos[0] = -1; 
    err_n0 = 0;
    for(i = R_bg; i < qlen; ++i) {
        if(target[i] != read_seq[i])  {
            err_pos[0] = i;
            err_n0++;
            break; 
        } 
    }
    if(err_pos[0] != -1) seq_sc = err_pos[0];
    else seq_sc = qlen;
    best_sc = seq_sc;
    best_i = err_pos[0] - 1;
    pd = err_pos[0];
    if(err_n0 == 0) {
        out->R_ti = max_mov + tar_mov + query->l_seq;
        out->R_qi = query->l_seq;
        out->R_sc = best_sc;
        return best_sc; 
    }
    int back_n0 = 0; 
    while(pd < qlen) {
        int k = pd + 1;
        err_n = err_n0;
        cur_sc = -5; 
        int ct_sc = 0;
        int ct_i = k;
        int sm_sc = 0;
        int sm_i = k;
        int sm_n = 0; 
        while(k < qlen){ 
            if(read_seq[k] < 4) { //N = 4
                if(target[k] == read_seq[k]) {
                    cur_sc++; 
                } else{
                    cur_sc -= 4;
                    err_pos[err_n++] = k; 
                } 
            } else{ // N得分
                cur_sc -= 1;
                err_pos[err_n++] = k; 
            }
            
            if(sm_sc > cur_sc) {
                sm_sc = cur_sc;
                sm_i = k;
                sm_n = err_n; 
            } 
            if(ct_sc < cur_sc) {
                ct_sc = cur_sc;
                ct_i = k; 
            } 


if( err_n > 1 && k - err_pos[0] >= 2*max_mov && qlen - k > max_mov) {// 12 or 24 ???
    ++k; 
    break; 
}
            ++k;
        } //  end while(k >= 0) ++++++++++
        pd = k - 1; 
        if(sm_sc < -5*3 && sm_i > max_mov){
            pd = sm_i;
            err_n = sm_n;
        }        

        //评估当前区间插入删除比对的可能性
        if(err_n == 0) {
            printf("%u, error : err_n = %d, pd = %d!!!\n", __LINE__, err_n, pd);    
            exit(1);
        }
        if(err_n == 1) {
            if(pd + 1 < qlen) { 
                printf("%u, error : err_n = %d, pd = %d!!!\n", __LINE__, err_n, pd);    
                exit(1);
            }
            
            cur_sc = qlen - err_pos[0] - 5;
            if(read_seq[err_pos[0]] < 4) {
                cur_sc = qlen - err_pos[0] - 5;
            } else {
                cur_sc = qlen - err_pos[0] - 2;
            }           
            cut_sc = 0;
            if(cur_sc < cut_sc) {
                seq_sc += cut_sc;
                best_i = err_pos[0]; 
                sc_flg = 1;
            } else {
                seq_sc += cur_sc;
                best_i = qlen; 
            } 
        
            best_sc = seq_sc;
            out->R_ti = max_mov + tar_mov + best_i;
            out->R_qi = best_i;
            out->R_sc = best_sc;
            return best_sc; 
        } // end if(err_n == 1) +++
    {
        err = err_pos[0] - seq_off[1] + 1 - seq_sc;
        int err0 = 0;
        for(i = 0; i < err_n; ++i) {
            if(read_seq[err_pos[i]] < 4) {
                err0 += 4; 
            } else {
                err0 += 1; 
            } 
        }
        if(err + err0 > max_err) {
            if(err0 < 8 || err + 7 > max_err) { 
                best_i = 0;
                out->R_ti = max_mov + tar_mov + best_i;
                out->R_qi = best_i;
                out->R_sc = -(err + err0);
                return -(err + err0);          
            }         
        }
    }  
        int buf_fnd_flg, buf_fnd_i, buf_gap_n, buf_gap_sc;
        int buf_cur_sc, buf_ext_L, buf_ext_R, buf_err_n;
        int buf_m_i, buf_m_sc, buf_ext_b, buf_ext_bi;
        int buf_pe;
        int ins_sc, del_sc, ins_i = 0, del_i = 0, ins_p, del_p;
        int ins_flg, del_flg;
        int PLS_SC = 5;
       
        int fnd_flg = 0;
        int gap_n = 0;
        int s_len = pd - err_pos[0] + 1;
        int min_len = max_mov/2; 
        len = err_n*5 - gap_o;
        if(len > max_mov) len = max_mov; 
        if(s_len < max_mov) {
            len = s_len/3;
            min_len = s_len/2; 
        }
        
        int ext_L = 0, ext_R = 0, ext_L0 = 0;
        int ext_b = 0, ext_bi = 0, ext_bn = 0;
        int ext_s = 0, ext_si = 0, ext_sn = 0;
        int cur_b = 0, cur_bi = 0;
        
        int sub_b = 0, sub_bi = err_pos[0]; 
        int sub_sc = pd - err_pos[0] + 1;
        N_sum = 0;
        for(j = 0; j < err_n; ++j) {
            int sc = sub_sc - (pd - err_pos[j] + 1);
            if(sc > sub_b) {
                sub_b = sc;
                sub_bi = err_pos[j]; 
            }
            if(read_seq[err_pos[j]] < 4){
                sub_sc -= 5;
            } else {
                sub_sc -= 2;
                ++N_sum;  
            } 
        }
        int bg = pd, end;
        int init_flg = 0;
        int gap_sc;
        int m_sc;
        int m_i, m_sc0;

        buf_fnd_flg = 0;
        buf_cur_sc = sub_sc;
        if(err_n > 3) len = 0;

        for(i = 1; i <= len; ++i){
            if(init_flg == 0) {
                ins_p   = err_n - 1;
                del_p   = err_n - 1;
                ins_sc  = 0;
                del_sc  = 0;
                gap_n   = 0;
                ins_flg = 0;
                del_flg = 0;
                fnd_flg = 0;
                int dp;
                for(dp = pd; dp <= pd+i; ++dp) {
                    if(err_pos[del_p] == dp) {
                        --del_p; 
                    } 
                }
                bg = pd;
                end = err_pos[0] + 1;
            }
            for(j = bg; j > end; --j){
                init_flg = 0;
                if(ins_flg == 0) {
                    k = j + i;
                    if(target[k] == read_seq[j]) ins_sc++;
                    else ins_sc -= 4;
                    if(ins_sc + (j - err_pos[0]) < sub_sc) { 
                        ins_flg = -1;
                    }
                    if(j == err_pos[ins_p]) {
                        if(pd - j > min_len) {
                            int eval_sc;
                            eval_sc = ins_sc + (j - err_pos[0]) - 5*ins_p;
                            eval_sc -= gap_o + gap_e*i; 
                            int g_sc = gap_o + i*gap_e;   
                            if( (pd >= qlen && eval_sc + qlen - err_pos[0] > 0)|| (pd < qlen && eval_sc + qlen - pd > 0 )) {
                                if(eval_sc > sub_sc + PLS_SC && j > g_sc) {//找到的read当中可插入比对位点
                                    //插入位置赋值
                                    gap_n = -i;
                                    fnd_flg = 1;
                                    --ins_p;
                                    bg = j;
                                    break; 
                                }
                            }
                        }
                        ins_p--; 
                    }
                }
                if(del_flg == 0 && j - i > end) {
                    k = j - i;
                    if(k > end) {
                        if(target[k] == read_seq[j]) del_sc++;
                        else del_sc -= 4;
                    }
                    if(del_sc + (k - err_pos[0]) < sub_sc) { 
                        del_flg = -1;
                    } 
                    if(k == err_pos[del_p]) {
                        if(pd - j > min_len) {
                            int eval_sc;
                            eval_sc = del_sc + (k - err_pos[0]) - 5*del_p;
                            eval_sc -= gap_o + gap_e*i; 
                            int g_sc = gap_o + 2*i*gap_e; 
                            if( (pd >= qlen && eval_sc + qlen - err_pos[0] > 0) || (pd < qlen && eval_sc + qlen - pd > 0 )) {
                                if(eval_sc > sub_sc + PLS_SC && j > g_sc) {//找到的read当中可删除比对位点
                                    //删除位置赋值
                                    gap_n = i;
                                    fnd_flg = 2;
                                    --del_p;
                                    bg = j; 
                                    break; 
                                }
                            }
                        }
                        del_p--; 
                    }
                }
                if(fnd_flg > 0) break;
                if(ins_flg != 0 && del_flg != 0) break;
            }// end for(j = bg; j < end; ++j)+++++++++
            if(j == end) {// j是计算过的 
                j = end + 1; 
            }
            ++j; // j是未计算的
            if(fnd_flg > 0) {
                int j0 = j, j1;
                int cur_p;
                if(fnd_flg == 1) {
                    cur_p = ins_p;
                    cur_sc = ins_sc;
                    m_sc = ins_sc;
                    m_i = err_pos[ins_p+1];
                } else {
                    cur_p = del_p;
                    cur_sc = del_sc;
                    m_sc = del_sc;
                    m_i = err_pos[del_p+1];
                }
                m_sc0 = m_sc + (j - err_pos[0]) - 5*cur_p;                
                ext_L = 0; 
                ext_R = 0;
                ext_L0 = 0;
                err_n0 = 0;
                int g_n = 0;
                if(gap_n > 0) g_n = gap_n;
                int back_n = back_n0;
                int max_i = j, max_n = 0, max_sum = 0;  
                int max_L = 0, max_L0 = 0;
                while(j > R_bg){
                    k = j + gap_n;
                    int h = j - g_n;
                    if(h < R_bg) break;
                    if(j >= err_pos[0]) {
                        if(target[k] == read_seq[j]) { 
                            cur_sc++;
                        } else{
                            err_pos0[err_n0++] = j;
                            cur_sc-=4;                            
                        }
                        if(h == err_pos[cur_p]) {
                            cur_p--; 
                        }
                        int c_sc = j - err_pos[0] - 5*(cur_p + 1); 
                        if(h >= err_pos[0]) {
                            if(m_sc0 <= cur_sc + c_sc) {
                                m_sc0 = cur_sc + c_sc;
                                m_sc = cur_sc;
                                m_i = j; 
                            }         
                        }
                        if(h < err_pos[0]) {
                            if(target[k] == read_seq[j]) { 
                                if(target[h] != read_seq[h]){ 
                                    ext_L++;
                                }
                            } else {
                                if(err_n0 == 0) {
                                    err_pos0[0] = j; 
                                    err_n0 = 1;
                                }   
                                if(back_n0 == 0) { 
                                    break;
                                } else if( target[h] == read_seq[h] ){
                                    ext_L0++;     
                                } 
                            }
                            if(back_n0 > 0) {
                                if(ext_L - ext_L0 >= max_sum) {
                                    max_sum = ext_L - ext_L0-1;
                                    max_L = ext_L;
                                    max_L0 = ext_L0; 
                                    max_i = j;
                                    max_n = err_n0;  
                                }
                                if(back_n <= 0) { break;}
                                if(ext_L - ext_L0 + 1< 0) break;
                                if(target[h] != read_seq[h]){ 
                                    back_n--;
                                    if(back_n <= 0) break; 
                                }
                            }
                        }
                    } else { //j > err_pos[0]+++
                     
                        if(target_e[k] == read_seq[j]) { 
                            if(target[h] != read_seq[h]){ 
                                ext_L++;
                            }
                        } else {
                            if(err_n0 == 0) {
                                err_pos0[0] = j; 
                                err_n0 = 1;
                            }   
                            if(back_n0 == 0) { 
                                break;
                            } else if( target[h] == read_seq[h] ){
                                ext_L0++;     
                            }
                        }
                        if(back_n0 > 0) {
                            if(ext_L - ext_L0 >= max_sum) {
                                max_sum = ext_L - ext_L0;
                                max_i = j;
                                max_n = err_n0; 
                                max_L = ext_L;
                                max_L0 = ext_L0;
                            }
                            if(back_n <= 0) { break;}
                            if(ext_L - ext_L0 + 1< 0) break;
                            if(target[h] != read_seq[h]){ 
                                back_n--;
                                if(back_n <= 0) break; 
                            }
                        }
                    } 
                    ++j;
                }// end while(j < L_ed) +++++
                if(back_n0 > 0) {
                    ext_L = max_L;
                    ext_L0 = max_L0;
                }
                if(err_n0 == 0) {
                    err_pos0[0] = j;
                    err_n0 = 1; 
                }
                cur_sc = m_sc;               
                int h;
                for(h = err_n0-1; h >= 0; --h) {
                    if(err_pos0[h] <= m_i) {
                        err_pos0[h] = err_pos0[0];     
                        for(k = 1; k <= h/2; ++k) {
                            int tmp = err_pos0[k];
                            err_pos0[k] = err_pos0[h-k]; 
                            err_pos0[h-k] = tmp; 
                        } 
                        err_pos0[0] = m_i;
                        break; 
                    }
                    if(h == 0) break; 
                } 
                err_n0 = h + 1;
                cur_b = cur_sc;
                cur_bi = pd;
                for(h = err_n0-1; h >= 0; --h) {
                    int sc = err_pos0[0] - err_pos0[h] - 5*h;
                    if(sc > cur_b) {
                        cur_b = sc;
                        cur_bi = err_pos0[h];
                    }                 
                }
                int un_n = 0;
                for(h = 0; h < err_n; ++h) {
                    if(err_pos[h] > err_pos0[0] + g_n) {
                        un_n++; 
                    } else { 
                        break;
                    }
                }
                int errp_var = err_pos0[0] - err_pos[0];
                if(errp_var < 0) errp_var = 0;
                cur_sc += errp_var + 1 - 5*un_n + 5*ext_L - 5*ext_L0;
                h = j0-1;
                ext_R = 0;
                ext_b = 0, ext_bi = 0, ext_bn = 0;
                ext_s = 0, ext_si = 0, ext_sn = 0;
                while(h <= qlen){
                    k = h + gap_n;
                    if(h <= pd) {//pb是计算过的点
                        if(target[k] != read_seq[h]) { 
                            err_pos0[err_n0++] = h;    
                        }
                        ++h;
                        continue;
                    } 
                    if(target[k] == read_seq[h]) { 
                        ext_R++;   
                    } else {
                        ext_R -= 4;
                        err_pos0[err_n0++] = h;
                        if(ext_R + qlen - h < ext_b && qlen - h >= max_mov) {//????
                            ext_sn = err_n0;
                            break; //while(h > = 0)
                        } else {
                        
                        } 
                    }
                    if(ext_R > ext_b) { 
                        ext_b = ext_R;
                        ext_bi = h;
                        ext_bn = err_n0;
                    }
                    ++h;
                } // end  while(h >= 0) +++
                pe = h;
                if(ext_b > 0) {
                    if(ext_bi > max_mov) {
                        //这时考虑左端存在插入删除的情况
                        cur_sc += ext_b;
                        err_n0 = ext_bn;
                        pe = ext_bi; 
                    } else {
                        cur_sc += ext_b;
                    } 
                } else {
                    pe = pd;
                    ext_R = 0;
                }         
                gap_sc = 0;
                if(fnd_flg == 1) {
                    gap_sc = gap_o + (-gap_n)*gap_e;
                } else if(fnd_flg == 2){
                    gap_sc = gap_o + 2*gap_n*gap_e;
                } else {
                    exit(1); 
                }
                cur_sc -= gap_sc; 
                if(pd <= 0) {
                    int bs = qlen - err_pos0[0] - 5*err_n0;            
                    int bi = 0, bn = err_n0; 
                    for(j = err_n0 -1; j >= 0; --j) {
                        int sc = err_pos0[j] - err_pos0[0] - 5*j;
                        if(sc > bs) {
                            bs = sc;
                            bi = err_pos0[j];
                            bn = j; 
                        }
                    }
                    cur_sc -= bi - 5*(err_n0 - bn); 
                } 
                if(cur_sc > sub_sc + PLS_SC && cur_sc >= 0) { //???
                    break;
                } else{ //if(cur_sc > sub_sc && cur_sc >= 0) +++
                    if(qlen - pd >0 && cur_sc > buf_cur_sc && cur_sc - sub_sc > gap_o + PLS_SC  && cur_sc + qlen - pd > 0) {
                        buf_fnd_flg = fnd_flg;
                        buf_fnd_i = i;
                        buf_gap_n = gap_n;
                        buf_gap_sc = gap_sc;
                        buf_cur_sc = cur_sc;
                        buf_ext_L = ext_L;
                        buf_ext_R = ext_R;
                        buf_err_n = err_n;
                        buf_m_i = m_i;
                        buf_m_sc = m_sc;
                        buf_ext_b = ext_b;
                        buf_ext_bi = ext_bi;
                        buf_pe = pe;         
                        for(j = 0; j < err_n0; ++j){
                            err_pos_buf[j] = err_pos0[j];
                        }
                        buf_err_n = err_n0; 
                    }
                    cur_sc = sub_sc;
                    
                    if(fnd_flg == 2) {
                        bg++;
                        del_flg = -1;
                    } else if(fnd_flg == 1){
                        ins_flg = -1; 
                    }
                   
                    fnd_flg = 0;
                    init_flg = 1;
                    --i;
                }
            } // if(fnd_flg > 0)                
        }// end for(i = 1; i < len; ++i) +++++++++
        if(buf_fnd_flg > 0 && i > len) {
            fnd_flg = buf_fnd_flg;
            i = buf_fnd_i;
            gap_n = buf_gap_n;
            gap_sc = buf_gap_sc;
            cur_sc = buf_cur_sc;
            ext_L = buf_ext_L;
            ext_R = buf_ext_R;
            err_n = buf_err_n;
            m_i = buf_m_i;
            m_sc = buf_m_sc;
            ext_b = buf_ext_b;
            ext_bi = buf_ext_bi;
            pe = buf_pe;
            for(j = 0; j < buf_err_n; ++j){
                err_pos0[j] = err_pos_buf[j];
            }
            err_n0 = buf_err_n;
        }
        if(len > 0) {
            if(fnd_flg > 0) {
                if(cur_sc > sub_b && seq_sc + cur_sc > best_sc) {
                    best_sc = seq_sc + cur_sc;
                    if(ext_b > 0) {
                        best_i = ext_bi;  
                    } else {
                        best_i = pd;  
                    }
                } else if(sub_b >= cur_sc && seq_sc + sub_b > best_sc) {
                    best_sc = seq_sc + sub_b;
                    best_i = sub_bi; 
                } 
            } else {
                if( seq_sc + sub_b > best_sc) {
                    best_sc = seq_sc + sub_b;
                    best_i = sub_bi;
                } 
            }
            if(fnd_flg > 0) {
                for(i = err_n0-1; i >= 0; --i ) {
                    if(err_pos0[i] >= pe) {
                        err_n0 = i + 1;
                        break; 
                    } 
                }
                int m_sc = pe - err_pos0[0] - 5*err_n0;
                int m_i = err_pos0[0];
                int c_sc, p_i = 0;
                int m_sc0, m_sc1;
                int sc0 = seq_sc + cur_sc;
                int sc1 = seq_sc + sub_sc;
                if(sc0 > best_sc) {
                    best_sc = sc0;
                    if(ext_bi == 0) {
                        best_i = pd; 
                    } else {
                        best_i = ext_bi; 
                    }
                } else if(sc1 > best_sc){
                    best_sc = sc1;
                    best_i = pd;
                }
                if(pe > pd) pd = pe;
                else if (pe < pd){
                    printf("%u, error!!!!, pe = %d, pd = %d, err_n = %d\n", __LINE__, pe, pd, err_n);
                    exit(1); 
                } 
            } // end if(fnd_flg > 0) ++++ 
        }
        if(len == 0) {
            int i, j, k, h;
            pd = qlen - 1; 
            if(err_n1 > 0) {
                seq_sc -= err_pos[0] - pos1 - 5*err_n1 + 3*N_sum1;
                err_pos[0] = pos1;
            }         
            kswst->gapo = 6;
            kswst->gape = 1;

if(err_pos[0] < 0) {
    exit(1);
}
            kswst->qlen = pd - err_pos[0] + 1;
            kswst->tlen = kswst->qlen + MAX_CLIP;
            if(pd < seq_sc) kswst->h0 = pd - err_pos[0] + 1;
            else kswst->h0 = seq_sc;
            kswst->h0 = LEN_READ;
            kswst->w = KSW_EXT_W;
            int bg = err_pos[0];
            for(i = bg; i <= pd; ++i) {
                j = i - bg;
                kswst->target[j] = target[i];
                kswst->query[j] = read_seq[i];
            }
            for(; i <= pd + max_mov; ++i) {
                kswst->target[++j] = target[i];
            }
            kswst->sc[0]  = 0;
            kswst->qle[0] = 0; 
            kswst->tle[0] = 0; 
            kswst->sc[1]  = 0;
            kswst->qle[1] = 0;
            kswst->tle[1] = 0;

            ksw_ext_short(5, sub->mat, kswst);

            int cut_sc = kswst->sc[0];
            int cut_i = kswst->qle[0]; 
            int cut_j = kswst->tle[0]; 
            int cur_sc = kswst->sc[1];
            int cur_i = kswst->qle[1];
            int cur_j = kswst->tle[1];
            int ret_sc = cut_sc + seq_sc - LEN_READ;

            out->R_ti = max_mov + bg + cut_j;
            out->R_qi = bg + cut_i;
            out->R_sc = ret_sc;
            return ret_sc;
        }
        int cut_sc = 0, cut_i = -1;  
        if(pd + 1 >= qlen || fnd_flg == 0) {
            for(i = 1; i < err_n; ++i) {
                int c_sc = err_pos[0] - err_pos[i] - 5*i;
                if(c_sc > cut_sc) {
                    cut_sc = c_sc; 
                    cut_i = err_pos[i];
                }     
            }
            if(sub_sc < cut_sc) {
                if(cut_sc + seq_sc > best_sc) {
                    best_sc = cut_sc + seq_sc;
                    best_i = cut_i+1;
                }
                if(pd + 1 < qlen) {
                    if(fnd_flg > 0) {
                        seq_sc += cur_sc;
                    } else {
                        seq_sc += sub_sc;
                    }  
                } else {
                    seq_sc = best_sc; 
                }
            } else { // if(sub_sc >= cut_sc)+++
                if(fnd_flg > 0) {
                    if(cur_sc + seq_sc > best_sc) {
                        best_sc = cur_sc + seq_sc;
                        best_i = pd; 
                    }  
                    seq_sc += cur_sc;
                } else {
                    if(sub_sc + seq_sc > best_sc) {
                        best_sc = sub_sc + seq_sc;
                        best_i = pd; 
                    }
                    seq_sc += sub_sc;
                }
            }
        } 
        if(fnd_flg > 0) {
            if(pd + 1 < qlen) seq_sc += cur_sc;
            if(err_n0 > 0) { 
                for(i = 0; i < err_n0; ++i) {
                    err_pos[i] = err_pos0[i];
                }    
                err_n = err_n0;
            } else {
                err_n = 1; 
            } 
            pos += gap_n;// ??????
            target_e += gap_n;
            target += gap_n;
            tar_mov += gap_n;
            if(gap_n >= 0) { max_gap = 2 + gap_n; }
            else { max_gap = 2 - gap_n; }
        }
        if(err_n <= 1 && fnd_flg == 0) {
            printf("%u, error: err_n == 0, pd = %d\n", __LINE__, pd);    
            exit(1);
        }
        if(pd + 1 >= qlen){
            break;
        }
       
        err = err_pos[0] - seq_off[1] + 1 - seq_sc;
        if(read_seq[err_pos[0]] < 4) {
            err += 4;
        } else {
            err += 1; 
        }
       
        if(max_err < err) { // ????
            best_i = 0;
            out->R_ti = max_mov + tar_mov + best_i;
            out->R_qi = best_i;
            out->R_sc = -err;
            return -err;
        } 

        err_n1 = 0;
        pos1 = 0;

        if(err_n > 0) {
            j = err_n - 1;
            N_sum1 = 0;
            if(j > 0 && pd > err_pos[j] && pd - err_pos[j] <= max_mov) {
                --j;
                if(read_seq[err_pos[j]] == 4) {
                    N_sum1++; 
                }
            }

            err_n1 = err_n - 1 - j;
            pos1 = err_pos[j]; 
            int sc = seq_sc - (pd - err_pos[err_n-1]);
            if(read_seq[err_pos[err_n-1]] < 4 ) seq_sc = sc + 4;
            else seq_sc = sc + 1;
            err_pos[0] = err_pos[err_n-1];

        } else {
            printf("%u, error : errn == 0!!!\n", __LINE__);
            exit(1);
        }
        err_n0 = 1; 
        if(err_n < 3) back_n0 = 0;
        else back_n0 = err_n - 2;
    } //end while(1) +++++
    if(seq_sc >= best_sc) { //是否有必要 ???
        best_sc = seq_sc;
        best_i = query->l_seq;
    }

    out->R_ti = max_mov + tar_mov + best_i;
    out->R_qi = best_i;
    out->R_sc = best_sc;
    return best_sc; 
}
int eval_pos_L0(idx_t *fm_idx,  query_t *query, uint32_t pos,  int seq_off[], struct SubBuf *sub)
{
    int max_mov = MAX_CLIP, len_km = 12; 
    int err = 0;
    int max_err = query->l_seq - query->b0 + sub->delta - sub->err_sum[0];
    aln_out_t *out = sub->aln_out; 
    int i, j, k;
    uint32_t pos_i;
  
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint8_t *target_e = query->target + max_mov;     
    uint8_t *target = target_e + max_mov;// 20 ????
    kswst_t *kswst = &sub->kswst;
    int *err_pos     = sub->eval_pos;
    int *err_pos0    = sub->eval_pos + LEN_READ;
    int *err_pos_buf = sub->eval_pos + 2*LEN_READ; 
    int qlen = query->l_seq;
    int read_len = query->l_seq + max_mov; 
    
    int seq_sc = -1;
    int best_sc = 0, best_i = -1;
    int tar_mov = 0;
    int L_ed = seq_off[0] + max_mov;
    for(j = 0; j < L_ed; ++j){ 
        pos_i = pos - max_mov + j;// 12 or 24 ???? 
        target_e[j] = __get_pac(fm_idx->pac, pos_i);
    }
    L_ed -= max_mov; 
    int sc_flg = 0;
    int cur_sc = 0;
    int cut_sc;
    int pb;
    int pe;
    int err_n = 0, err_n0 = 0, err_n1 = 0, N_sum1 = 0, N_sum = 0, pos1 = 0; 
    int max_gap = max_mov;
    int len;
    int gap_o = 6, gap_e = 1;
    
    err_pos[0] = -1; 
    err_n0 = 0;
    for(i = L_ed-1; i >= 0; --i) {
        if(target[i] != read_seq[i])  {
            err_pos[0] = i;
            err_n0++;
            break; 
        } 
    }

    if(err_pos[0] != -1) seq_sc = qlen - 1 - err_pos[0];
    else seq_sc = qlen;
    best_sc = seq_sc;
    best_i = err_pos[0] + 1;
    pb = err_pos[0];
    if(err_n0 == 0) {
        out->L_ti = max_mov + tar_mov + best_i;
        out->L_qi = best_i;
        out->L_sc = best_sc;
        return best_sc; 
    }
    int back_n0 = 0; 
    while(pb > 0) {
        int k = pb - 1;
        err_n = err_n0;
        cur_sc = -5; 
        int ct_sc = 0;
        int ct_i = k;
        int sm_sc = 0;
        int sm_i = k;
        int sm_n = 0; 
        while(k >= 0){ 
            if(read_seq[k] < 4) { //N = 4
                if(target[k] == read_seq[k]) {
                    cur_sc++; 
                } else{
                    cur_sc -= 4;
                    err_pos[err_n++] = k; 
                } 
            } else{ // N得分
                cur_sc -= 1;
                err_pos[err_n++] = k; 
            }
            if(sm_sc > cur_sc) {
                sm_sc = cur_sc;
                sm_i = k;
                sm_n = err_n; 
            } 
            if(ct_sc < cur_sc) {
                ct_sc = cur_sc;
                ct_i = k; 
            } 
//+++++++++++++++++++++++++++++++++++++++

if( err_n > 1 && err_pos[0] - k >= 2*max_mov && k > max_mov) {// 12 or 24 ???
    --k;
    break; 
}

//--------------------------------------
            --k;
        } //  end while(k >= 0) ++++++++++

fprintf(stderr, "%u, seq_sc = %d, k = %d, err_n = %d, cur_sc = %d, sm_i = %d, sm_sc = %d\n", __LINE__, seq_sc, k, err_n, cur_sc, sm_i, sm_sc);    
        pb = k + 1; 
        if(sm_sc < -5*3 && sm_i > max_mov){
            pb = sm_i;
            err_n = sm_n;
        }        
fprintf(stderr, "%u, seq_sc = %d, err_pos[0] = %d, err_n = %d, pb = %d\n", __LINE__, seq_sc, err_pos[0], err_n, pb);    
        //评估当前区间插入删除比对的可能性
        
        if(err_n == 0) {
            printf("%u, error : err_n = %d, pb = %d!!!\n", __LINE__, err_n, pb);    
            exit(1);
        }
        
fprintf(stderr, "%u, seq_sc = %d\n", __LINE__, seq_sc);    

        if(err_n == 1) {
            if(pb > 0) { 
                printf("%u, error : err_n = %d, pb = %d!!!\n", __LINE__, err_n, pb);    
                exit(1);
            } 
            if(read_seq[err_pos[0]] < 4) {
                cur_sc = err_pos[0] + 1 -5;
            } else {
                cur_sc = err_pos[0] + 1 -2;
            }

            cut_sc = 0;
            if(cur_sc < cut_sc) {
                seq_sc += cut_sc;
                best_i = err_pos[0];     
                sc_flg = 1;
            } else {
                seq_sc += cur_sc;
                best_i = 0;
            } 
            best_sc = seq_sc;
            out->L_ti = max_mov + tar_mov + best_i;
            out->L_qi = best_i;
            out->L_sc = best_sc;

            
            return best_sc; 
        } // end if(err_n == 1) +++
fprintf(stderr, "%u, seq_sc = %d, cur_sc = %d\n", __LINE__, seq_sc, cur_sc);    
//+++++++++++++++++++++++++++++++++
    {
        err = err_pos[0] - seq_off[1] + 1 - seq_sc;
        int err0 = 0;
        for(i = 0; i < err_n; ++i) {
            if(read_seq[err_pos[i]] < 4) {
                err0 += 4; 
            } else {
                err0 += 1; 
            } 
        }
        if(err + err0 > max_err) {
            if(err0 < 8 || err + 7 > max_err) { 
                best_i = 0;
                out->L_ti = max_mov + tar_mov + best_i;
                out->L_qi = best_i;
                out->L_sc = -(err + err0);
                return -(err + err0);          
            }         
        }
    }  
//---------------------------------

        int buf_fnd_flg, buf_fnd_i, buf_gap_n, buf_gap_sc;
        int buf_cur_sc, buf_ext_L, buf_ext_R, buf_err_n;
        int buf_m_i, buf_m_sc, buf_ext_b, buf_ext_bi;
        int buf_pe;
        int ins_sc, del_sc, ins_i = 0, del_i = 0, ins_p, del_p;
        int ins_flg, del_flg;
        int PLS_SC = 5;
       
        int fnd_flg = 0;
        int gap_n = 0;
        int s_len = err_pos[0] - pb + 1;
        int min_len = max_mov/2; 
        len = err_n*5 - gap_o;
        if(len > max_mov) len = max_mov; 
        if(s_len < max_mov) {
            len = s_len/3;
            min_len = s_len/2; 
        }
//+++++++++++++++++++++
//len = 0; 
//--------------------
fprintf(stderr, "%u, seq_sc = %d, max_gap = %d, min_len = %d, s_len = %d, pb = %d\n", __LINE__, seq_sc, max_gap, min_len, s_len, pb);    
        cur_sc = 0;

       
        int ext_L = 0, ext_R = 0, ext_R0 = 0;
        int ext_b = 0, ext_bi = 0, ext_bn = 0;
        int ext_s = 0, ext_si = 0, ext_sn = 0;
        int cur_b = 0, cur_bi = 0;
        int sub_sc = err_pos[0] - pb + 1;            
        int sub_b = 0, sub_bi = err_pos[0]; 
        for(j = 0; j < err_n; ++j) {
            int sc = sub_sc - (err_pos[j] - pb + 1);
            if(sc > sub_b) {
                sub_b = sc;
                sub_bi = err_pos[j] + 1; 
            }
            if(read_seq[err_pos[j]] < 4) {
                sub_sc -= 5;
            } else {
                sub_sc -= 2;
                ++N_sum;
            }
        
        }
        
        
        int bg = pb, end;
        int init_flg = 0;
        int gap_sc;
        int m_sc;
        int m_i, m_sc0;

        buf_fnd_flg = 0;
        buf_cur_sc = sub_sc;
fprintf(stderr, "%u, len = %d, sub_sc = %d\n", __LINE__, len, sub_sc);       
        if(err_n > 3) len = 0;
        for(i = 1; i <= len; ++i){
            if(init_flg == 0) {
                ins_p   = err_n - 1;
                del_p   = err_n - 1;
                ins_sc  = 0;
                del_sc  = 0;
                gap_n   = 0;
                ins_flg = 0;
                del_flg = 0;
                fnd_flg = 0;
                int dp;

                for(dp = pb; dp <= pb+i; ++dp) {
                    if(err_pos[del_p] == dp) {
                        --del_p; 
                    } 
                }

                bg = pb;
                end = err_pos[0] + 1;
            }
fprintf(stderr, "%u, i = %d, pb = %d, end = %d\n", __LINE__, i, pb, end);    
            
            for(j = bg; j < end; ++j){
                init_flg = 0;
                if(ins_flg == 0) {
                    k = j - i + max_mov;
                    if(target_e[k] == read_seq[j]) ins_sc++;
                    else ins_sc -= 4;
                    if(ins_sc + (err_pos[0] - j) < sub_sc) { 
                        ins_flg = -1;
                    }
                    if(j == err_pos[ins_p]) {
                        if(j - pb > min_len) {
                            int eval_sc;
                            eval_sc = ins_sc + (err_pos[0] - j) - 5*ins_p;
                            eval_sc -= gap_o + gap_e*i; 
                            int g_sc = gap_o + i*gap_e;   
                            if( (pb <= 0 && eval_sc + err_pos[0] > 0)|| (pb > 0 && eval_sc + pb > 0 )) {
                                if(eval_sc > sub_sc + PLS_SC && j > g_sc) {//找到的read当中可插入比对位点
                                    //插入位置赋值
                                    gap_n = -i;
                                    fnd_flg = 1;
                                    --ins_p;
                                    bg = j;
                                    break; 
                                }
                            }
                        }
                        ins_p--; 
                    }
                }

                if(del_flg == 0 && j + i < end) {
                    k = j + i;
                    if(k < end) {
                        if(target[k] == read_seq[j]) del_sc++;
                        else del_sc -= 4;
                    }
                    if(del_sc + (err_pos[0] - k) < sub_sc) { 
                        del_flg = -1;
                    } 
                    if(j+i == err_pos[del_p]) {
                        if(j - pb > min_len) {
                            int eval_sc;
                            eval_sc = del_sc + (err_pos[0] - k) - 5*del_p;
                            eval_sc -= gap_o + gap_e*i; 
                            int g_sc = gap_o + 2*i*gap_e; 
                            if( (pb <= 0 && eval_sc + err_pos[0] > 0) || (pb > 0 && eval_sc + pb > 0 )) {
                                if(eval_sc > sub_sc + PLS_SC && j > g_sc) {//找到的read当中可删除比对位点
                                    //删除位置赋值
                                    gap_n = i;
                                    fnd_flg = 2;
                                    --del_p;
                                    bg = j; 
                                    break; 
                                }
                            }
                        }
                        del_p--; 
                    }
                }
                if(fnd_flg > 0) break;
                if(ins_flg != 0 && del_flg != 0) break;
            }// end for(j = bg; j < end; ++j)+++++++++
            if(j == end) {// j是计算过的 
                j = end - 1; 
            }
            ++j; // j是未计算的
            //+++++++++++++++++++++++++++
            //---------------------------
            if(fnd_flg > 0) {
                int j0 = j, j1;
                int cur_p;
                if(fnd_flg == 1) {
                    cur_p = ins_p;
                    cur_sc = ins_sc;
                    m_sc = ins_sc;
                    m_i = err_pos[ins_p+1];
                } else {
                    cur_p = del_p;
                    cur_sc = del_sc;
                    m_sc = del_sc;
                    m_i = err_pos[del_p+1];
                }

                m_sc0 = m_sc - err_pos[0] - j - 5*cur_p;
                ext_L = 0; 
                ext_R = 0;
                ext_R0 = 0;
                err_n0 = 0;
                int g_n = 0;
                if(gap_n > 0) g_n = gap_n;
                int back_n = back_n0;
                int max_i = j, max_n = 0, max_sum = 0;  
                int max_j = k;
                int max_R = 0, max_R0 = 0;
                while(j < L_ed){
                    k = j + gap_n + max_mov;
                    int h = j + g_n;
                    if(h >= L_ed) break;
                    //if(j + g_n <= err_pos[0]) {
                    if(j <= err_pos[0]) {
                        if(target_e[k] == read_seq[j]) { 
                            cur_sc++;
                        } else{
                            err_pos0[err_n0++] = j;
                            cur_sc-=4;                            
                        }
                        if(h == err_pos[cur_p]) {
                            cur_p--; 
                        }
                        int c_sc = err_pos[0] - j - 5*(cur_p + 1); 
                        if(h <= err_pos[0]) {
                            if(m_sc0 <= cur_sc + c_sc) {
                                m_sc0 = cur_sc + c_sc;
                                m_sc = cur_sc;
                                m_i = j; 
                            }         
                        }
                        if(h > err_pos[0]) {
                            if(target_e[k] == read_seq[j]) { 
                                if(target[h] != read_seq[h]){ 
                                    ext_R++;
                                }
                            } else {
                                if(err_n0 == 0) {
                                    err_pos0[0] = j; 
                                    err_n0 = 1;
                                }   
                                if(back_n0 == 0) { 
                                    break;
                                } else if( target[h] == read_seq[h] ){
                                    ext_R0++;     
                                } 
                            }
                            
                            if(back_n0 > 0) {
                                if(ext_R - ext_R0 >= max_sum) {
                                    max_sum = ext_R - ext_R0;
                                    max_i = j;
                                    max_j = k;
                                    max_n = err_n0; 
                                    max_R = ext_R;
                                    max_R0 = ext_R0;
                                }
                                if(back_n <= 0) { break;}
                                if(ext_R - ext_R0 + 1< 0) break;

                                if(target[h] != read_seq[h]){ 
                                    back_n--;
                                    if(back_n <= 0) break; 
                                }
                            }
                        }
                    } else { //j > err_pos[0]+++
                     
                        if(target_e[k] == read_seq[j]) { 
                            if(target[h] != read_seq[h]){ 
                                ext_R++;
                            }
                        } else {
                            if(err_n0 == 0) {
                                err_pos0[0] = j; 
                                err_n0 = 1;
                            }   
                            if(back_n0 == 0) { 
                                break;
                            } else if( target[h] == read_seq[h] ){
                                ext_R0++;     
                            }
                        }
                        if(back_n0 > 0) {
                            if(ext_R - ext_R0 >= max_sum) {
                                max_sum = ext_R - ext_R0;
                                max_i = j;
                                max_j = k;
                                max_n = err_n0; 
                                max_R = ext_R;
                                max_R0 = ext_R0;
                            }
                            if(back_n <= 0) { break;}
                            if(ext_R - ext_R0 + 1< 0) break;

                            if(target[h] != read_seq[h]){ 
                                back_n--;
                                if(back_n <= 0) break; 
                            }
                        }
                    } 
                    ++j;
                }// end while(j < L_ed) +++++
fprintf(stderr, "%u, j = %d, err_n0 = %d, m_sc = %d, m_i = %d, cur_sc = %d, ext_R= %d, ext_R0 = %d, max_i = %d, max_n = %d, max_sum = %d, back_n0 = %d\n", __LINE__, j, err_n0, m_sc, m_i, cur_sc, ext_R, ext_R0, max_i, max_n, max_sum, back_n0);
                if(back_n0 > 0) {
                    ext_R = max_R;
                    ext_R0 = max_R0;
                }
                if(err_n0 == 0) {
                    err_pos0[0] = j;
                    err_n0 = 1; 
                }
                cur_sc = m_sc;               
                int h;
fprintf(stderr, "%u, err_n0 = %d, err_pos0[0] = %d\n", __LINE__, err_n0, err_pos0[0]);
                for(h = err_n0-1; h >= 0; --h) {
                    if(err_pos0[h] <= m_i) {
                        err_pos0[h] = err_pos0[0];     
                        for(k = 1; k <= h/2; ++k) {
                            int tmp = err_pos0[k];
                            err_pos0[k] = err_pos0[h-k]; 
                            err_pos0[h-k] = tmp; 
                        } 
                        err_pos0[0] = m_i;
                        break; 
                    }
                    if(h == 0) break; 
                } 
                err_n0 = h + 1;
                cur_b = cur_sc;
                cur_bi = pb;
                for(h = err_n0-1; h >= 0; --h) {
                    int sc = err_pos0[0] - err_pos0[h] - 5*h;
                    if(sc > cur_b) {
                        cur_b = sc;
                        cur_bi = err_pos0[h] + 1;
                    }                 
                }
                 
                
                int un_n = 0;
                for(h = 0; h < err_n; ++h) {
                    if(err_pos[h] > err_pos0[0] + g_n) {
                        un_n++; 
                    } else { 
                        break;
                    }
                }
fprintf(stderr, "%u, un_n = %d, ext_R = %d, ext_R0 = %d, cur_sc = %d, err_pos[0] = %d, err_pos0[0] = %d\n", __LINE__, un_n, ext_R, ext_R0, cur_sc, err_pos[0], err_pos0[0]);
                int errp_var = err_pos[0] - err_pos0[0];
                if(errp_var < 0) errp_var = 0;
                cur_sc += errp_var + 1 - 5*un_n + 5*ext_R - 5*ext_R0;

                h = j0-1;
                ext_L = 0;
                ext_b = 0, ext_bi = 0, ext_bn = 0;
                ext_s = 0, ext_si = 0, ext_sn = 0;
                  
fprintf(stderr, "%d, j0 = %d, pb = %d\n", __LINE__, j0, pb); 
                while(h >= 0){
                    k = h + gap_n + max_mov;
                    if(h >= pb) {//pb是计算过的点
                        if(target_e[k] != read_seq[h]) { 
                            err_pos0[err_n0++] = h;    
                        }
                        --h;
                        continue;
                    } 
                    if(target_e[k] == read_seq[h]) { 
                        ext_L++;   
                    } else {
                        ext_L -= 4;
                        err_pos0[err_n0++] = h;
                        if(ext_L + h < ext_b && h >= max_mov) {//????
                            ext_sn = err_n0;
                            break; //while(h > = 0)
                        } else {
                        
                        } 
                    }
                    if(ext_L > ext_b) { 
                        ext_b = ext_L;
                        ext_bi = h;
                        ext_bn = err_n0;
                    }
                    --h;
                } // end  while(h >= 0) +++
                pe = h;
                if(ext_b > 0) {
                    if(ext_bi > max_mov) {
                        //这时考虑左端存在插入删除的情况
                        cur_sc += ext_b;
                        err_n0 = ext_bn;
                        pe = ext_bi; 
                    } else {
                        cur_sc += ext_b;
                    } 
                } else {
                    pe = pb;
                    ext_L = 0;
                }         
                gap_sc = 0;
                if(fnd_flg == 1) {
                    gap_sc = gap_o + (-gap_n)*gap_e;
                } else if(fnd_flg == 2){
                    gap_sc = gap_o + 2*gap_n*gap_e;
                } else {
                    exit(1); 
                }
                cur_sc -= gap_sc; 
                if(pb <= 0) {
                    int bs = err_pos0[0] + 1 - 5*err_n0;            
                    int bi = 0, bn = err_n0; 
                    for(j = err_n0 -1; j >= 0; --j) {
                        int sc = err_pos0[0] - err_pos0[j] - 5*j;
                        if(sc > bs) {
                            bs = sc;
                            bi = err_pos0[j] + 1;
                            bn = j; 
                        }
                    }
                    cur_sc -= bi - 5*(err_n0 - bn); 
                } 
fprintf(stderr, "%u, seq_sc = %d, sub_sc = %d, cur_sc = %d, gap_sc = %d, gap_n = %d, cur_b = %d, cur_bi = %d\n", __LINE__, seq_sc, sub_sc, cur_sc, gap_sc, gap_n, cur_b, cur_bi); 
                if(cur_sc > sub_sc + PLS_SC && cur_sc >= 0) { //???
                    break;
                } else{ //if(cur_sc > sub_sc && cur_sc >= 0) +++
                    if(pb >0 && cur_sc > buf_cur_sc && cur_sc - sub_sc > gap_o + PLS_SC  && cur_sc + pb > 0) {
                        buf_fnd_flg = fnd_flg;
                        buf_fnd_i = i;
                        buf_gap_n = gap_n;
                        buf_gap_sc = gap_sc;
                        buf_cur_sc = cur_sc;
                        buf_ext_L = ext_L;
                        buf_ext_R = ext_R;
                        buf_err_n = err_n;
                        buf_m_i = m_i;
                        buf_m_sc = m_sc;
                        buf_ext_b = ext_b;
                        buf_ext_bi = ext_bi;
                        buf_pe = pe;         
                        for(j = 0; j < err_n0; ++j){
                            err_pos_buf[j] = err_pos0[j];
                        }
                        buf_err_n = err_n0; 
                    }
                    cur_sc = sub_sc;
                    
                    if(fnd_flg == 2) {
                        bg++;
                        del_flg = -1;
                    } else if(fnd_flg == 1){
                        ins_flg = -1; 
                    }
                   
                    fnd_flg = 0;
                    init_flg = 1;
                    --i;
                }
            } // if(fnd_flg > 0)                
        }// end for(i = 1; i < len; ++i) +++++++++
        if(buf_fnd_flg > 0 && i > len) {
            fnd_flg = buf_fnd_flg;
            i = buf_fnd_i;
            gap_n = buf_gap_n;
            gap_sc = buf_gap_sc;
            cur_sc = buf_cur_sc;
            ext_L = buf_ext_L;
            ext_R = buf_ext_R;
            err_n = buf_err_n;
            m_i = buf_m_i;
            m_sc = buf_m_sc;
            ext_b = buf_ext_b;
            ext_bi = buf_ext_bi;
            pe = buf_pe;
            for(j = 0; j < buf_err_n; ++j){
                err_pos0[j] = err_pos_buf[j];
            }
            err_n0 = buf_err_n;
        }
//fprintf(stderr, "%u, err_n = %d, err_pos[0] = %d, pb = %d, fnd_flg = %d, ext_b = %d, ext_bi = %d, m_sc = %d, m_i = %d, gap_n = %d\n", __LINE__, err_n, err_pos[0], pb, fnd_flg, ext_b, ext_bi, m_sc, m_i, gap_n);
        if(len > 0) {
            if(fnd_flg > 0) {
                if(cur_sc > sub_b && seq_sc + cur_sc > best_sc) {
                    best_sc = seq_sc + cur_sc;
                    if(ext_b > 0) {
                        best_i = ext_bi;  
                    } else {
                        best_i = pb;  
                    }
                } else if(sub_b >= cur_sc && seq_sc + sub_b > best_sc) {
                    best_sc = seq_sc + sub_b;
                    best_i = sub_bi; 
                } 
            } else {
                if( seq_sc + sub_b > best_sc) {
                    best_sc = seq_sc + sub_b;
                    best_i = sub_bi; 
                } 
            }
            if(fnd_flg > 0) {
                for(i = err_n0-1; i >= 0; --i ) {
                    if(err_pos0[i] >= pe) {
                        err_n0 = i + 1;
                        break; 
                    } 
                }
                int m_sc = err_pos0[0] - pe - 5*err_n0;
                int m_i = err_pos0[0];
                int c_sc, p_i = 0;
                int m_sc0, m_sc1;
                int sc0 = seq_sc + cur_sc;
                int sc1 = seq_sc + sub_sc;
                if(sc0 > best_sc) {
                    best_sc = sc0;
                    if(ext_bi == 0) {
                        best_i = pb; 
                    } else {
                        best_i = ext_bi; 
                    }
                } else if(sc1 > best_sc){
                    best_sc = sc1;
                    best_i = pb;
                }
                if(pe < pb) pb = pe;
                else if (pe > pb){
                    printf("%u, error!!!!, pe = %d, pb = %d, err_n = %d\n", __LINE__, pe, pb, err_n);
                    exit(1); 
                } 
            } // end if(fnd_flg > 0) ++++ 
fprintf(stderr, "%u, fnd_flg = %d, seq_sc = %d, cur_sc = %d, best_sc = %d, err_pos[0] = %d, err_n = %d\n", __LINE__, fnd_flg, seq_sc, cur_sc, best_sc, err_pos[0], err_n);           
        }
fprintf(stderr, "%u, len = %d, err_n = %d, err_n1 = %d, err_pos[0] = %d, pos1 = %d, N_sum1 = %d\n", __LINE__, len, err_n, err_n1, err_pos[0], pos1, N_sum1);        
        if(len == 0) {
            int i, j, k, h;
            pb = 0; 
            if(err_n1 > 0) {
                //seq_sc -= err_pos[0] - pos1 - 5*err_n1 + 3*N_sum1;
                seq_sc -= pos1 - err_pos[0] - 5*err_n1 + 3*N_sum1;
                err_pos[0] = pos1;
            }
            
            
            kswst->gapo = 6;
            kswst->gape = 1;
            kswst->qlen = err_pos[0] - pb + 1;
            kswst->tlen = kswst->qlen + MAX_CLIP;
            if(pb < seq_sc) kswst->h0 = err_pos[0] - pb;
            else kswst->h0 = seq_sc;
            kswst->h0 = LEN_READ;
            kswst->w = KSW_EXT_W;
            for(i = err_pos[0]; i >= pb; --i) {
                j = err_pos[0] - i;
                kswst->target[j] = target[i];
                kswst->query[j] = read_seq[i];
            }
fprintf(stderr, "%u, len = %d, err_n = %d, err_pos[0] = %d, pb = %d, i = %d, j = %d\n", __LINE__, len, err_n, err_pos[0], pb, i, j);        
            for(i = pb+MAX_CLIP-1; i >= pb; --i) {
                kswst->target[++j] = target_e[i];
            }
fprintf(stderr, "%u, len = %d, err_n = %d, err_pos[0] = %d, pb = %d, i = %d, j = %d\n", __LINE__, len, err_n, err_pos[0], pb, i, j);        
            kswst->sc[0]  = 0;
            kswst->qle[0] = 0; 
            kswst->tle[0] = 0; 
            kswst->sc[1]  = 0;
            kswst->qle[1] = 0;
            kswst->tle[1] = 0;
            ksw_ext_short(5, sub->mat, kswst);
            int cut_sc = kswst->sc[0];
            int cut_i = kswst->qle[0]; 
            int cut_j = kswst->tle[0]; 
            int cur_sc = kswst->sc[1];
            int cur_i = kswst->qle[1];
            int cur_j = kswst->tle[1];
            int ret_sc = cut_sc + seq_sc - LEN_READ;
            out->L_ti = max_mov + err_pos[0] + 1 - cut_j;
            out->L_qi = err_pos[0] + 1 - cut_i;
            out->L_sc = ret_sc;
            return ret_sc;
        }
        int cut_sc = 0, cut_i = -1;  
fprintf(stderr, "%u, pb = %d, err_pos[0] = %d, err_n = %d, cur_sc = %d, sub_sc = %d\n", __LINE__, pb, err_pos[0], err_n, cur_sc, sub_sc);              
        if(pb <= 0 || fnd_flg == 0) {
            for(i = 1; i < err_n; ++i) {
                int c_sc = err_pos[0] - err_pos[i] - 5*i;
                if(c_sc > cut_sc) {
                    cut_sc = c_sc; 
                    cut_i = err_pos[i];
                }     
            }
fprintf(stderr, "%u, cut_i = %d, cut_sc = %d, sub_sc = %d, cur_sc = %d\n", __LINE__, cut_i, cut_sc, sub_sc, cur_sc);              
            if(sub_sc < cut_sc) {
                if(cut_sc + seq_sc > best_sc) {
                    best_sc = cut_sc + seq_sc;
                    best_i = cut_i+1;
fprintf(stderr, "%u, best_sc = %d, best_i = %d\n", __LINE__, best_sc, best_i);    
                }
                if(pb > 0) {
                    if(fnd_flg > 0) {
                        seq_sc += cur_sc;
                    } else {
fprintf(stderr, "%u, seq_sc = %d, sub_sc = %d\n", __LINE__, seq_sc, sub_sc);    
                        seq_sc += sub_sc;
                    }  

                } else {
                    seq_sc = best_sc; 
                }
            } else { // if(sub_sc >= cut_sc)+++
                if(fnd_flg > 0) {
                    if(cur_sc + seq_sc > best_sc) {
                        best_sc = cur_sc + seq_sc;
                        best_i = pb; 
                    }  
                    seq_sc += cur_sc;
                } else {
                    if(sub_sc + seq_sc > best_sc) {
                        best_sc = sub_sc + seq_sc;
                        best_i = pb; 
                    }
                    seq_sc += sub_sc;
                }
            }
        } 
        if(fnd_flg > 0) {
            if(pb > 0) seq_sc += cur_sc;
            if(err_n0 > 0) { 
                for(i = 0; i < err_n0; ++i) {
                    err_pos[i] = err_pos0[i];
                }    
                err_n = err_n0;
            } else {
                err_n = 1; 
            } 
fprintf(stderr, "%u, gap_n = %d\n", __LINE__, gap_n);    
            pos += gap_n;// ??????
            target_e += gap_n;
            target += gap_n;
            tar_mov += gap_n;
            if(gap_n >= 0) { max_gap = 2 + gap_n; }
            else { max_gap = 2 - gap_n; }
        }
        if(err_n <= 1 && fnd_flg == 0) {
            printf("%u, error: err_n == 0, pb = %d\n", __LINE__, pb);    
            exit(1);
        }
        if(pb <= 0){ break; }
        int len_aln = 0;
        if(pb == 0) {
            len_aln = seq_off[1];
        } else {
            len_aln = seq_off[1] - err_pos[0];
        }
        err = len_aln - seq_sc;
        if(max_err < err) { // ????
            best_i = 0;
            out->L_ti = max_mov + tar_mov + best_i;
            out->L_qi = best_i;
            out->L_sc = -err;
            return -err;
        } 
        err_n1 = 0;
        pos1 = 0;
        if(err_n > 0) {
            j = err_n - 1;
            N_sum1 = 0;
            if(j > 0 && pb < err_pos[j] && err_pos[j] - pb <= max_mov) {
                --j;
                if(read_seq[err_pos[j]] == 4){
                    N_sum1++;
                }
            } 
            err_n1 = err_n - 1 - j;
            pos1 = err_pos[j]; 
            int sc = seq_sc - (err_pos[err_n-1] - pb);
            if (read_seq[err_pos[err_n -1]] < 4) seq_sc = sc + 4;
            else seq_sc = sc + 1;
            err_pos[0] = err_pos[err_n-1];
fprintf(stderr, "%u, seq_sc = %d, err_pos[0] = %d, sc = %d, err_n = %d, pb = %d, err_pos[err_n-1] = %d, read_seq[i] = %d\n", __LINE__, seq_sc, err_pos[0], sc, err_n, pb, err_pos[err_n-1], read_seq[err_pos[err_n-1]]);            
        } else {
            printf("%u, error : errn == 0!!!\n", __LINE__);
            exit(1);
        }
        err_n0 = 1; 
        if(err_n < 3) back_n0 = 0;
        else back_n0 = err_n - 2;
    } //end while(1) +++++
    if(seq_sc > best_sc) { //是否有必要 ???
        best_sc = seq_sc;
        best_i = 0;
    }
    out->L_ti = max_mov + tar_mov + best_i;
    out->L_qi = best_i;
    out->L_sc = best_sc;
    return best_sc; 
}
int aln_pos_L_0(idx_t *fm_idx,  query_t *query, uint32_t pos_buf[],  int seq_off[], struct SubBuf *sub)
{
    aln_out_t *out = sub->aln_out; 
    int i, j, k;
    uint32_t pos, pos_i;
    pos = pos_buf[0];
    int max_mov = MAX_CLIP;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint8_t *target_e = query->target + max_mov;     
    uint8_t *target = target_e + max_mov;// 20 ????
    int *err_pos     = sub->eval_pos;
    int seq_sc = -1;
    int best_sc = 0, best_i = -1;
    int tar_mov = 0;

    int L_ed = seq_off[0] + max_mov; 
    for(j = 0; j < L_ed; ++j){ 
        pos_i = pos - max_mov + j;// 12 or 24 ???? 
        target_e[j] = __get_pac(fm_idx->pac, pos_i);
    }
    L_ed -= max_mov; 
    int cur_sc = 0;
    int pb;
    int err_n = 0, err_n0 = 0; 
    err_pos[0] = -1; 
    err_n0 = 0;

    for(i = L_ed-1; i >= 0; --i) {
        if(target[i] != read_seq[i])  {
            err_pos[0] = i;
            err_n0++;
            break; 
        } 
    }
fprintf(stderr, "%u, seq_sc = %d, i = %d\n", __LINE__, seq_sc, i);
    if(err_pos[0] != -1) seq_sc = seq_off[1] - 1 - err_pos[0];
    else seq_sc = seq_off[1];
fprintf(stderr, "%u, seq_sc = %d\n", __LINE__, seq_sc);
    best_sc = seq_sc;
    best_i = err_pos[0] + 1;
    pb = err_pos[0];

    if(err_n0 == 0) {
        out->L_ti = max_mov + tar_mov + best_i;
        out->L_qi = best_i;
        out->L_sc = best_sc;
        return best_sc; 
    }
    if(pb <= sub->sub_err) {
        out->L_ti = max_mov + tar_mov + best_i;
        out->L_qi = best_i;
        out->L_sc = best_sc;
        return best_sc; 
    }
    k = pb - 1;
    err_n = err_n0;
    cur_sc = -5; 
    while(k >= 0){ 
        if(read_seq[k] < 4) { //N = 4
            if(target[k] == read_seq[k]) {
                cur_sc++; 
            } else{
                cur_sc -= 4;
                err_pos[err_n++] = k; 
            } 
        } else{ // N得分
            cur_sc -= 1;
            err_pos[err_n++] = k; 
        }

        if( err_n > 1) {
            --k;
            break; 
        }            
        --k;
    } //  end while(k >= 0) ++++++++++
    pb = k + 1; 
    //评估当前区间插入删除比对的可能性
    if(err_n == 0) {
        printf("%u, error : err_n = %d, pb = %d!!!\n", __LINE__, err_n, pb);    
        exit(1);
    }
    if(err_n <= 2) {
        if(read_seq[err_pos[0]] < 4) {
            cur_sc = err_pos[0] + 1 -5;
        } else {
            cur_sc = err_pos[0] + 1 -2;
        }
        seq_sc += cur_sc;
    } // end if(err_n <= 2) +++
    if(err_n == 1) {
        best_i = 0; 
        best_sc = seq_sc;
        out->L_ti = max_mov + tar_mov + best_i;
        out->L_qi = best_i;
        out->L_sc = best_sc; 
        return best_sc;
    } 
    
     
    if(err_n == 2 && err_pos[1] < sub->sub_err) {
      
        cur_sc = -(err_pos[1] + 1); 
fprintf(stderr, "%u, cur_sc = %d, seq_sc = %d, err_pos[1] = %d\n", __LINE__, cur_sc, seq_sc, err_pos[1]);
     
        seq_sc += cur_sc;
        best_i = err_pos[1] + 1; 
        best_sc = seq_sc;
        out->L_ti = max_mov + tar_mov + best_i;
        out->L_qi = best_i;
        out->L_sc = best_sc; 
        return best_sc;
    } 
            
    return 0; 
}
int aln_pos_L(idx_t *fm_idx,  query_t *query, uint32_t pos_buf[],  int seq_off[], struct SubBuf *sub)
{
    int max_mov = MAX_CLIP, len_km = 12; 
    int err = 0;
    int max_err = query->l_seq - query->b0 + sub->delta - sub->err_sum[0];
    //int max_err = query->query_err + sub->delta;

    aln_out_t *out = sub->aln_out; 
    
    int i, j, k, ret = 0;
    uint32_t pos, pos_i;
    pos = pos_buf[0];
    uint32_t (*target_idx)[2] = query->target_idx; 
    //uint8_t *read_seq = query->read_seq;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint8_t *target_e = query->target + max_mov;     
    uint8_t *target = target_e + max_mov;// 20 ????
    kswst_t *kswst = &sub->kswst;
    int *err_pos     = sub->eval_pos;
    int *err_pos0    = sub->eval_pos + LEN_READ;
    int *err_pos_buf = sub->eval_pos + 2*LEN_READ; 
    //int eval_pos[2*LEN_READ];
 
    int qlen = query->l_seq;
    int read_len = query->l_seq + max_mov; 
    int l_len = seq_off[0] + len_km - 1;
    int l_num = l_len/len_km;
    int l_bg = l_len-l_num*len_km;
    int r_len = query->l_seq - seq_off[1] + len_km -1;
    int r_num = r_len/len_km;
    int r_bg = seq_off[1] - 1;
    
    int *l_sc  = sub->eval_sc;  
    int *l_sum = sub->eval_sc + LEN_READ/12;
    int *r_sc  = sub->eval_sc + 2*(LEN_READ/12); 
    int *r_sum = sub->eval_sc + 3*(LEN_READ/12); 
    
    int l_max_sc, r_max_sc, l_max_i, r_max_i;
    int L_ed = l_bg;

    int seq_sc = -1;
    int best_sc = 0, best_i = -1;
    int tar_mov = 0;

    L_ed = seq_off[0] + max_mov; 
    for(j = 0; j < L_ed; ++j){ 
        pos_i = pos - max_mov + j;// 12 or 24 ???? 
        target_e[j] = __get_pac(fm_idx->pac, pos_i);
    }
    L_ed -= max_mov; 

    int sc_flg = 0;
    int cur_sc = 0,  cut_sc = 0;
    int err_num = 0;
    int err_num0 = 0, rep;

    int pb;
    int pe;
    int err_n = 0, err_n0 = 0, err_n1 = 0, N_sum1 = 0, N_sum = 0, pos1 = 0; 
    int max_gap = max_mov;
    int len;
    int gap_o = 6, gap_e = 1;
    
    err_pos[0] = -1; 
    err_n0 = 0;

    for(i = L_ed-1; i >= 0; --i) {
        if(target[i] != read_seq[i])  {
            err_pos[0] = i;
            err_n0++;
            break; 
        } 
    }

    if(err_pos[0] != -1) seq_sc = qlen - 1 - err_pos[0];
    else seq_sc = qlen;
    best_sc = seq_sc;
    best_i = err_pos[0] + 1;
    pb = err_pos[0];

    if(err_n0 == 0) {
        out->L_ti = max_mov + tar_mov + best_i;
        out->L_qi = best_i;
        out->L_sc = best_sc;
        return best_sc; 
    }

    int back_n0 = 0; 
    while(pb > 0) {
        int k = pb - 1;
        err_n = err_n0;
        cur_sc = -5; 
        int ct_sc = 0;
        int ct_i = k;
        int sm_sc = 0;
        int sm_i = k;
        int sm_n = 0; 
        while(k >= 0){ 
            if(read_seq[k] < 4) { //N = 4
                if(target[k] == read_seq[k]) {
                    cur_sc++; 
                } else{
                    cur_sc -= 4;
                    //if(err_n == 0) 
                    err_pos[err_n++] = k; 
                } 
            } else{ // N得分
                cur_sc -= 1;
                err_pos[err_n++] = k; 
            }
            if(sm_sc > cur_sc) {
                sm_sc = cur_sc;
                sm_i = k;
                sm_n = err_n; 
            } 
            if(ct_sc < cur_sc) {
                ct_sc = cur_sc;
                ct_i = k; 
            } 
//+++++++++++++++++++++++++++++++++++++++

if( err_n > 1 && err_pos[0] - k >= 2*max_mov && k > max_mov) {// 12 or 24 ???
    --k;
    break; 
}

//--------------------------------------
            --k;
        } //  end while(k >= 0) ++++++++++

        pb = k + 1; 
        if(sm_sc < -5*3 && sm_i > max_mov){
            pb = sm_i;
            err_n = sm_n;
        }        
        //if(seq_sc < 0) seq_sc = query->l_seq - err_pos[0];         
        //评估当前区间插入删除比对的可能性
        if(err_n == 0) {
            printf("%u, error : err_n = %d, pb = %d!!!\n", __LINE__, err_n, pb);    
            exit(1);
        }
        if(err_n == 1) {
            if(pb > 0) { 
                printf("%u, error : err_n = %d, pb = %d!!!\n", __LINE__, err_n, pb);    
                exit(1);
            } 
            if(read_seq[err_pos[0]] < 4) {
                cur_sc = err_pos[0] + 1 -5;
            } else {
                cur_sc = err_pos[0] + 1 -2;
            }
            cut_sc = 0;
            if(cur_sc < cut_sc) {
                seq_sc += cut_sc;
                best_i = err_pos[0];     
                sc_flg = 1;
            } else {
                seq_sc += cur_sc;
                best_i = 0;
            } 
            best_sc = seq_sc;
            out->L_ti = max_mov + tar_mov + best_i;
            out->L_qi = best_i;
            out->L_sc = best_sc;

            
            return best_sc; 
        } // end if(err_n == 1) +++
    {
        err = err_pos[0] - seq_off[1] + 1 - seq_sc;
        int err0 = 0;
        for(i = 0; i < err_n; ++i) {
            if(read_seq[err_pos[i]] < 4) {
                err0 += 4; 
            } else {
                err0 += 1; 
            } 
        }
        if(err + err0 > max_err) {
            if(err0 < 8 || err + 7 > max_err) { 
                best_i = 0;
                out->L_ti = max_mov + tar_mov + best_i;
                out->L_qi = best_i;
                out->L_sc = -(err + err0);
                return -(err + err0);          
            }         
        }
    }  
//---------------------------------
        //if(err_n > 1) {
        int buf_fnd_flg, buf_fnd_i, buf_gap_n, buf_gap_sc;
        int buf_cur_sc, buf_ext_L, buf_ext_R, buf_err_n;
        int buf_m_i, buf_m_sc, buf_ext_b, buf_ext_bi;
        int buf_pe;
        int ins_sc, del_sc, ins_i = 0, del_i = 0, ins_p, del_p;
        int ins_flg, del_flg;
        int PLS_SC = 5;
       
        int fnd_flg = 0;
        int gap_n = 0;
        int s_len = err_pos[0] - pb + 1;
        int min_len = max_mov/2; 
        len = err_n*5 - gap_o;
        if(len > max_mov) len = max_mov; 
        if(s_len < max_mov) {
            len = s_len/3;
            min_len = s_len/2; 
        }
//+++++++++++++++++++++
//len = 0; 
//--------------------
fprintf(stderr, "%u, seq_sc = %d, max_gap = %d, min_len = %d, s_len = %d, pb = %d\n", __LINE__, seq_sc, max_gap, min_len, s_len, pb);    
        cur_sc = 0;

       
        int ext_L = 0, ext_R = 0, ext_R0 = 0;
        int ext_b = 0, ext_bi = 0, ext_bn = 0;
        int ext_s = 0, ext_si = 0, ext_sn = 0;
        int cur_b = 0, cur_bi = 0;

        int sub_sc = err_pos[0] - pb + 1;            
        int sub_b = 0, sub_bi = err_pos[0]; 
        for(j = 0; j < err_n; ++j) {
            int sc = sub_sc - (err_pos[j] - pb + 1);
            if(sc > sub_b) {
                sub_b = sc;
                sub_bi = err_pos[j] + 1; 
            }
            if(read_seq[err_pos[j]] < 4) {
                sub_sc -= 5;
            } else {
                sub_sc -= 2;
                ++N_sum;
            }
        }
        
        
        int bg = pb, end;
        int init_flg = 0;
        int gap_sc;
        int m_sc;
        int m_i, m_sc0;

        buf_fnd_flg = 0;
        buf_cur_sc = sub_sc;

        if(err_n > 3) len = 0;
        for(i = 1; i <= len; ++i){
            if(init_flg == 0) {
                ins_p   = err_n - 1;
                del_p   = err_n - 1;
                ins_sc  = 0;
                del_sc  = 0;
                gap_n   = 0;
                ins_flg = 0;
                del_flg = 0;
                fnd_flg = 0;
                int dp;

                for(dp = pb; dp <= pb+i; ++dp) {
                    if(err_pos[del_p] == dp) {
                        --del_p; 
                    } 
                }

                bg = pb;
                end = err_pos[0] + 1;
            }

            for(j = bg; j < end; ++j){
                init_flg = 0;
                if(ins_flg == 0) {
                    k = j - i + max_mov;
                    if(target_e[k] == read_seq[j]) ins_sc++;
                    else ins_sc -= 4;
                    if(ins_sc + (err_pos[0] - j) < sub_sc) { 
                        ins_flg = -1;
                    }
                    if(j == err_pos[ins_p]) {
                        if(j - pb > min_len) {
                            int eval_sc;
                            eval_sc = ins_sc + (err_pos[0] - j) - 5*ins_p;
                            eval_sc -= gap_o + gap_e*i; 
                            int g_sc = gap_o + i*gap_e;   
                            if( (pb <= 0 && eval_sc + err_pos[0] > 0)|| (pb > 0 && eval_sc + pb > 0 )) {
                                if(eval_sc > sub_sc + PLS_SC && j > g_sc) {//找到的read当中可插入比对位点
                                    //插入位置赋值
                                    gap_n = -i;
                                    fnd_flg = 1;
                                    --ins_p;
                                    bg = j;
                                    break; 
                                }
                            }
                        }
                        ins_p--; 
                    }
                }

                if(del_flg == 0 && j + i < end) {
                    k = j + i;
                    if(k < end) {
                        if(target[k] == read_seq[j]) del_sc++;
                        else del_sc -= 4;
                    }
                    if(del_sc + (err_pos[0] - k) < sub_sc) { 
                        del_flg = -1;
                    } 
                    if(j+i == err_pos[del_p]) {
                        if(j - pb > min_len) {
                            int eval_sc;
                            eval_sc = del_sc + (err_pos[0] - k) - 5*del_p;
                            eval_sc -= gap_o + gap_e*i; 
                            int g_sc = gap_o + 2*i*gap_e; 
                            if( (pb <= 0 && eval_sc + err_pos[0] > 0) || (pb > 0 && eval_sc + pb > 0 )) {
                                if(eval_sc > sub_sc + PLS_SC && j > g_sc) {//找到的read当中可删除比对位点
                                    //删除位置赋值
                                    gap_n = i;
                                    fnd_flg = 2;
                                    --del_p;
                                    bg = j; 
                                    break; 
                                }
                            }
                        }
                        del_p--; 
                    }
                }
                if(fnd_flg > 0) break;
                if(ins_flg != 0 && del_flg != 0) break;
            }// end for(j = bg; j < end; ++j)+++++++++
            if(j == end) {// j是计算过的 
                j = end - 1; 
            }
            ++j; // j是未计算的
            //+++++++++++++++++++++++++++
            //---------------------------
            if(fnd_flg > 0) {
                int j0 = j, j1;
                int cur_p;
                if(fnd_flg == 1) {
                    cur_p = ins_p;
                    cur_sc = ins_sc;
                    m_sc = ins_sc;
                    m_i = err_pos[ins_p+1];
                } else {
                    cur_p = del_p;
                    cur_sc = del_sc;
                    m_sc = del_sc;
                    m_i = err_pos[del_p+1];
                }

                m_sc0 = m_sc - err_pos[0] - j - 5*cur_p;
                ext_L = 0; 
                ext_R = 0;
                ext_R0 = 0;
                err_n0 = 0;
                int g_n = 0;
                if(gap_n > 0) g_n = gap_n;
                int back_n = back_n0;
                int max_i = j, max_n = 0, max_sum = 0;  
                int max_j = k;
                int max_R = 0, max_R0 = 0;
                while(j < L_ed){
                    k = j + gap_n + max_mov;
                    int h = j + g_n;
                    if(h >= L_ed) break;
                    //if(j + g_n <= err_pos[0]) {
                    if(j <= err_pos[0]) {
                        if(target_e[k] == read_seq[j]) { 
                            cur_sc++;
                        } else{
                            err_pos0[err_n0++] = j;
                            cur_sc-=4;                            
                        }
                        if(h == err_pos[cur_p]) {
                            cur_p--; 
                        }
                        int c_sc = err_pos[0] - j - 5*(cur_p + 1); 
                        if(h <= err_pos[0]) {
                            if(m_sc0 <= cur_sc + c_sc) {
                                m_sc0 = cur_sc + c_sc;
                                m_sc = cur_sc;
                                m_i = j; 
                            }         
                        }
                        if(h > err_pos[0]) {
                            if(target_e[k] == read_seq[j]) { 
                                if(target[h] != read_seq[h]){ 
                                    ext_R++;
                                }
                            } else {
                                if(err_n0 == 0) {
                                    err_pos0[0] = j; 
                                    err_n0 = 1;
                                }   
                                if(back_n0 == 0) { 
                                    break;
                                } else if( target[h] == read_seq[h] ){
                                    ext_R0++;     
                                } 
                            }
                            
                            if(back_n0 > 0) {
                                if(ext_R - ext_R0 >= max_sum) {
                                    max_sum = ext_R - ext_R0;
                                    max_i = j;
                                    max_j = k;
                                    max_n = err_n0; 
                                    max_R = ext_R;
                                    max_R0 = ext_R0;
                                }
                                if(back_n <= 0) { break;}
                                if(ext_R - ext_R0 + 1< 0) break;

                                if(target[h] != read_seq[h]){ 
                                    back_n--;
                                    if(back_n <= 0) break; 
                                }
                            }
                        }
                    } else { //j > err_pos[0]+++
                        if(target_e[k] == read_seq[j]) { 
                            if(target[h] != read_seq[h]){ 
                                ext_R++;
                            }
                        } else {
                            if(err_n0 == 0) {
                                err_pos0[0] = j; 
                                err_n0 = 1;
                            }   
                            if(back_n0 == 0) { 
                                break;
                            } else if( target[h] == read_seq[h] ){
                                ext_R0++;     
                            }
                        }
                        if(back_n0 > 0) {
                            if(ext_R - ext_R0 >= max_sum) {
                                max_sum = ext_R - ext_R0;
                                max_i = j;
                                max_j = k;
                                max_n = err_n0; 
                                max_R = ext_R;
                                max_R0 = ext_R0;
                            }
                            if(back_n <= 0) { break;}
                            if(ext_R - ext_R0 + 1< 0) break;

                            if(target[h] != read_seq[h]){ 
                                back_n--;
                                if(back_n <= 0) break; 
                            }
                        }
                    } 
                    ++j;
                }// end while(j < L_ed) +++++

                if(back_n0 > 0) {
                    ext_R = max_R;
                    ext_R0 = max_R0;
                }
                if(err_n0 == 0) {
                    err_pos0[0] = j;
                    err_n0 = 1; 
                }
                cur_sc = m_sc;               
                int h;
                for(h = err_n0-1; h >= 0; --h) {
                    if(err_pos0[h] <= m_i) {
                        err_pos0[h] = err_pos0[0];     
                        for(k = 1; k <= h/2; ++k) {
                            int tmp = err_pos0[k];
                            err_pos0[k] = err_pos0[h-k]; 
                            err_pos0[h-k] = tmp; 
                        } 
                        err_pos0[0] = m_i;
                        break; 
                    }
                    if(h == 0) break; 
                } 
                err_n0 = h + 1;
                cur_b = cur_sc;
                cur_bi = pb;
                for(h = err_n0-1; h >= 0; --h) {
                    int sc = err_pos0[0] - err_pos0[h] - 5*h;
                    if(sc > cur_b) {
                        cur_b = sc;
                        cur_bi = err_pos0[h] + 1;
                    }                 
                }
                 
                
                int un_n = 0;
                for(h = 0; h < err_n; ++h) {
                    if(err_pos[h] > err_pos0[0] + g_n) {
                        un_n++; 
                    } else { 
                        break;
                    }
                }
                int errp_var = err_pos[0] - err_pos0[0];
                if(errp_var < 0) errp_var = 0;
                cur_sc += errp_var + 1 - 5*un_n + 5*ext_R - 5*ext_R0;

                h = j0-1;
                ext_L = 0;
                ext_b = 0, ext_bi = 0, ext_bn = 0;
                ext_s = 0, ext_si = 0, ext_sn = 0;
                  
                while(h >= 0){
                    k = h + gap_n + max_mov;
                    if(h >= pb) {//pb是计算过的点
                        if(target_e[k] != read_seq[h]) { 
                            err_pos0[err_n0++] = h;    
                        }
                        --h;
                        continue;
                    } 
                    if(target_e[k] == read_seq[h]) { 
                        ext_L++;   
                    } else {
                        ext_L -= 4;
                        err_pos0[err_n0++] = h;
                        if(ext_L + h < ext_b && h >= max_mov) {//????
                            ext_sn = err_n0;
                            break; //while(h > = 0)
                        } else {
                        
                        } 
                    }
                    if(ext_L > ext_b) { 
                        ext_b = ext_L;
                        ext_bi = h;
                        ext_bn = err_n0;
                    }
                    --h;
                } // end  while(h >= 0) +++
                pe = h;
                if(ext_b > 0) {
                    if(ext_bi > max_mov) {
                        //这时考虑左端存在插入删除的情况
                        cur_sc += ext_b;
                        err_n0 = ext_bn;
                        pe = ext_bi; 
                    } else {
                        cur_sc += ext_b;
                    } 
                } else {
                    pe = pb;
                    ext_L = 0;
                }         
                gap_sc = 0;
                if(fnd_flg == 1) {
                    gap_sc = gap_o + (-gap_n)*gap_e;
                } else if(fnd_flg == 2){
                    gap_sc = gap_o + 2*gap_n*gap_e;
                } else {
                    exit(1); 
                }
                cur_sc -= gap_sc; 
                if(pb <= 0) {
                    int bs = err_pos0[0] + 1 - 5*err_n0;            
                    int bi = 0, bn = err_n0; 
                    for(j = err_n0 -1; j >= 0; --j) {
                        int sc = err_pos0[0] - err_pos0[j] - 5*j;
                        if(sc > bs) {
                            bs = sc;
                            bi = err_pos0[j] + 1;
                            bn = j; 
                        }
                    }
                    cur_sc -= bi - 5*(err_n0 - bn); 
                } 
                if(cur_sc > sub_sc + PLS_SC && cur_sc >= 0) { //???
                    break;
                } else{ //if(cur_sc > sub_sc && cur_sc >= 0) +++
                    if(pb >0 && cur_sc > buf_cur_sc && cur_sc - sub_sc > gap_o + PLS_SC  && cur_sc + pb > 0) {
                        buf_fnd_flg = fnd_flg;
                        buf_fnd_i = i;
                        buf_gap_n = gap_n;
                        buf_gap_sc = gap_sc;
                        buf_cur_sc = cur_sc;
                        buf_ext_L = ext_L;
                        buf_ext_R = ext_R;
                        buf_err_n = err_n;
                        buf_m_i = m_i;
                        buf_m_sc = m_sc;
                        buf_ext_b = ext_b;
                        buf_ext_bi = ext_bi;
                        buf_pe = pe;         
                        for(j = 0; j < err_n0; ++j){
                            err_pos_buf[j] = err_pos0[j];
                        }
                        buf_err_n = err_n0; 
                    }
                    cur_sc = sub_sc;
                    
                    if(fnd_flg == 2) {
                        bg++;
                        del_flg = -1;
                    } else if(fnd_flg == 1){
                        ins_flg = -1; 
                    }
                   
                    fnd_flg = 0;
                    init_flg = 1;
                    --i;
                }
            } // if(fnd_flg > 0)                
        }// end for(i = 1; i < len; ++i) +++++++++
        if(buf_fnd_flg > 0 && i > len) {
            fnd_flg = buf_fnd_flg;
            i = buf_fnd_i;
            gap_n = buf_gap_n;
            gap_sc = buf_gap_sc;
            cur_sc = buf_cur_sc;
            ext_L = buf_ext_L;
            ext_R = buf_ext_R;
            err_n = buf_err_n;
            m_i = buf_m_i;
            m_sc = buf_m_sc;
            ext_b = buf_ext_b;
            ext_bi = buf_ext_bi;
            pe = buf_pe;
            for(j = 0; j < buf_err_n; ++j){
                err_pos0[j] = err_pos_buf[j];
            }
            err_n0 = buf_err_n;
        }
        if(len > 0) {
            if(fnd_flg > 0) {
                if(cur_sc > sub_b && seq_sc + cur_sc > best_sc) {
                    best_sc = seq_sc + cur_sc;
                    if(ext_b > 0) {
                        best_i = ext_bi;  
                    } else {
                        best_i = pb;  
                    }
                } else if(sub_b >= cur_sc && seq_sc + sub_b > best_sc) {
                    best_sc = seq_sc + sub_b;
                    best_i = sub_bi; 
                } 
            } else {
                if( seq_sc + sub_b > best_sc) {
                    best_sc = seq_sc + sub_b;
                    best_i = sub_bi; 
                } 
            }
            if(fnd_flg > 0) {
                for(i = err_n0-1; i >= 0; --i ) {
                    if(err_pos0[i] >= pe) {
                        err_n0 = i + 1;
                        break; 
                    } 
                }
                int m_sc = err_pos0[0] - pe - 5*err_n0;
                int m_i = err_pos0[0];
                int c_sc, p_i = 0;
                int m_sc0, m_sc1;
               
                int sc0 = seq_sc + cur_sc;
                int sc1 = seq_sc + sub_sc;
                if(sc0 > best_sc) {
                    best_sc = sc0;
                    if(ext_bi == 0) {
                        best_i = pb; 
                    } else {
                        best_i = ext_bi; 
                    }
                } else if(sc1 > best_sc){
                    best_sc = sc1;
                    best_i = pb;
                }
                if(pe < pb) pb = pe;
                else if (pe > pb){
                    printf("%u, error!!!!, pe = %d, pb = %d, err_n = %d\n", __LINE__, pe, pb, err_n);
                    exit(1); 
                } 
            } // end if(fnd_flg > 0) ++++ 
        }

        if(len == 0) {
            int i, j, k, h;
            pb = 0; 
            if(err_n1 > 0) {
                seq_sc -= pos1 - err_pos[0] - 5*err_n1 + 3*N_sum1;
                err_pos[0] = pos1;
            }
            
            kswst->gapo = 6;
            kswst->gape = 1;
            kswst->qlen = err_pos[0] - pb + 1;
            kswst->tlen = kswst->qlen + MAX_CLIP;

            if(pb < seq_sc) kswst->h0 = err_pos[0] - pb;
            else kswst->h0 = seq_sc;
            kswst->h0 = LEN_READ;
            kswst->w = KSW_EXT_W;

            for(i = err_pos[0]; i >= pb; --i) {
                j = err_pos[0] - i;
                kswst->target[j] = target[i];
                kswst->query[j] = read_seq[i];
            }
            for(i = pb+MAX_CLIP-1; i >= pb; --i) {
                kswst->target[++j] = target_e[i];
            }

            kswst->sc[0]  = 0;
            kswst->qle[0] = 0; 
            kswst->tle[0] = 0; 
            kswst->sc[1]  = 0;
            kswst->qle[1] = 0;
            kswst->tle[1] = 0;

            ksw_ext_short(5, sub->mat, kswst);

            int cut_sc = kswst->sc[0];
            int cut_i = kswst->qle[0]; 
            int cut_j = kswst->tle[0]; 
            int cur_sc = kswst->sc[1];
            int cur_i = kswst->qle[1];
            int cur_j = kswst->tle[1];

            int ret_sc = cut_sc + seq_sc - LEN_READ;
            out->L_ti = max_mov + err_pos[0] + 1 - cut_j;
            out->L_qi = err_pos[0] + 1 - cut_i;
            out->L_sc = ret_sc;


            return ret_sc;
        }
        int cut_sc = 0, cut_i = -1;  
        if(pb <= 0 || fnd_flg == 0) {
            for(i = 1; i < err_n; ++i) {
                int c_sc = err_pos[0] - err_pos[i] - 5*i;
                if(c_sc > cut_sc) {
                    cut_sc = c_sc; 
                    cut_i = err_pos[i];
                }     
            }
            if(sub_sc < cut_sc) {
                if(cut_sc + seq_sc > best_sc) {
                    best_sc = cut_sc + seq_sc;
                    best_i = cut_i+1;
                }
                if(pb > 0) {
                    if(fnd_flg > 0) {
                        seq_sc += cur_sc;
                    } else {
                        seq_sc += sub_sc;
                    }  

                } else {
                    seq_sc = best_sc; 
                }
            } else { // if(sub_sc >= cut_sc)+++
                if(fnd_flg > 0) {
                    if(cur_sc + seq_sc > best_sc) {
                        best_sc = cur_sc + seq_sc;
                        best_i = pb; 
                    }  
                    seq_sc += cur_sc;
                } else {
                    if(sub_sc + seq_sc > best_sc) {
                        best_sc = sub_sc + seq_sc;
                        best_i = pb; 
                    }
                    seq_sc += sub_sc;
                }
            }
        } 
        if(fnd_flg > 0) {
            if(pb > 0) seq_sc += cur_sc;
            if(err_n0 > 0) { 
                for(i = 0; i < err_n0; ++i) {
                    err_pos[i] = err_pos0[i];
                }    
                err_n = err_n0;
            } else {
                err_n = 1; 
            } 
            pos += gap_n;// ??????
            target_e += gap_n;
            target += gap_n;
            tar_mov += gap_n;
            if(gap_n >= 0) { max_gap = 2 + gap_n; }
            else { max_gap = 2 - gap_n; }
        }
        if(err_n <= 1 && fnd_flg == 0) {
            printf("%u, error: err_n == 0, pb = %d\n", __LINE__, pb);    
            exit(1);
        }        
        if(pb <= 0){
            break;
        }
        int len_aln = 0;
        if(pb == 0) {
            len_aln = seq_off[1];
        } else {
            len_aln = seq_off[1] - err_pos[0];
        }
        err = len_aln - seq_sc;
        if(max_err < err) { // ????
            best_i = 0;
            out->L_ti = max_mov + tar_mov + best_i;
            out->L_qi = best_i;
            out->L_sc = -err;
            return -err;
        } 
        err_n1 = 0;
        pos1 = 0;
        if(err_n > 0) {
            j = err_n - 1;
            N_sum1 = 0;
            if(j > 0 && pb < err_pos[j] && err_pos[j] - pb <= max_mov) {
                --j;
                if(read_seq[err_pos[j]] == 4){
                    N_sum1++;
                }
            
            } 
            err_n1 = err_n - 1 - j;
            pos1 = err_pos[j]; 
            
            int sc = seq_sc - (err_pos[err_n-1] - pb);
            if (read_seq[err_pos[err_n -1]] < 4) seq_sc = sc + 4;
            else seq_sc = sc + 1;
            err_pos[0] = err_pos[err_n-1];

        } else {
            printf("%u, error : errn == 0!!!\n", __LINE__);
            exit(1);
        }
        err_n0 = 1; 
        if(err_n < 3) back_n0 = 0;
        else back_n0 = err_n - 2;
    } //end while(1) +++++
    if(seq_sc > best_sc) { //是否有必要 ???
        best_sc = seq_sc;
        best_i = 0;
    }
    out->L_ti = max_mov + tar_mov + best_i;
    out->L_qi = best_i;
    out->L_sc = best_sc;
    return best_sc; 
}


int eval_pos(idx_t *fm_idx,  query_t *query, uint32_t pos_buf[],  int seq_off[], uint32_t err_out[], struct SubBuf *sub)
{
    uint32_t pos_num = err_out[0];  
    sub->calcu_flg[0] = 1;
    sub->calcu_flg[1] = 2;
    sub->calcu_flg[2] = 12;
    int sub_flg   = sub->calcu_flg[0];
    int indel_flg = sub->calcu_flg[1];
    int bg_flg    = sub->calcu_flg[2];
    int ins_buf[LEN_READ/12], del_buf[LEN_READ/12];
    int i, j, k;
    uint32_t (*target_idx)[2] = query->target_idx; 
    //uint8_t *read_seq = query->read_seq;
    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t pos, pos_i, i_bwt, seq12;
    //uint32_t l_flg[33];
    uint32_t *isa = fm_idx->bwt->isa;

    int l_len = seq_off[0];
    int r_len = query->l_seq - seq_off[1];
    int l_num = l_len/12;
    int r_num = r_len/12;
    int l_bg = l_len-l_num*12;
    int r_bg = seq_off[1];
    int *l_sc  = sub->eval_sc;  
    int *l_sum = sub->eval_sc + LEN_READ/12;
    int *r_sc  = sub->eval_sc + 2*(LEN_READ/12); 
    int *r_sum = sub->eval_sc + 3*(LEN_READ/12); 
    int cur_sc, l_max_sc, r_max_sc, l_max_i, r_max_i;
    pos = pos_buf[0];
    for(i = 0; i < l_num; ++i) {
        j = l_bg + 12*i; 
        pos_i = pos + j; 
        i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
        if(i_bwt >= target_idx[j][0] && i_bwt < target_idx[j][1] ){
            l_sc[i] = 1;  
        } else{
            l_sc[i] = 0; 
        }
//fprintf(stderr, "%u, l_sc[%d] = %d\n", __LINE__, i, l_sc[i]);
    }
    for(i = 0; i < r_num; ++i) {
        j = r_bg + 12*i; 
        pos_i = pos + j; 
        i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
        if(i_bwt >= target_idx[j][0] && i_bwt < target_idx[j][1] ){
            r_sc[i] = 1;  
        } else{
            r_sc[i] = 0; 
        }
    }
    l_sum[0] = l_sc[0];
    for(i = 1; i < l_num; ++i) {
        l_sum[i] = l_sum[i-1] + l_sc[i]; 
    }
    r_sum[0] = r_sc[0];
    for(i = 1; i < r_num; ++i) {
        r_sum[i] = r_sum[i-1] + r_sc[i]; 
    }

    l_max_sc = -7;
    if(l_num > 0) {
        l_max_sc = 5*(l_sum[l_num-1] - l_sum[1])  + 4*(1-l_sum[0]) -7;
        l_max_i = 1;
        for(i = 2; i < l_num; ++i){
           cur_sc = 5*(l_sum[l_num-1] - l_sum[i])  + 4*(i - l_sum[i-1]) -7;
           if(l_max_sc < cur_sc) {
                l_max_sc = cur_sc;
                l_max_i = i;
           }
        }
    }
    r_max_sc = -7;
    if(r_num > 0) {
        r_max_sc =  4*((r_num -1) - (r_sum[r_num-1]-r_sum[0])) -7; 
        r_max_i = 0;
        for(i = 1; i < r_num -2; ++i){
            cur_sc = 5*r_sum[i-1] + 4*((r_num -1 - i) - (r_sum[r_num-1] - r_sum[i])) - 7;
            if(r_max_sc < cur_sc) {
                r_max_sc = cur_sc;
                r_max_i = i;
            }
        }
    }
//fprintf(stderr, "%u, l_num = %u, l_sum[l_num-1] = %u\n", __LINE__, l_num, l_sum[l_num-1]);
//fprintf(stderr, "%u, r_num = %u, r_sum[r_num-1] = %u\n", __LINE__, r_num, r_sum[r_num-1]);
    int l_sub_err = 0, r_sub_err = 0;
    if(l_num > 0) l_sub_err = (l_num - l_sum[l_num-1])*5;
    //else l_sub_err = 0;
    if(r_num > 0) r_sub_err = (r_num - r_sum[r_num-1])*5;
    //else r_sub_err = 0;
    int sub_err = l_sub_err + r_sub_err;
    int l_indel = 5*l_num - l_max_sc  + r_sub_err;
    int r_indel = 5*r_num - r_max_sc  + l_sub_err;

    
////fprintf(stderr, "%u, l_num = %u, l_sum[l_num-1] = %u\n", __LINE__, l_num, l_sum[l_num-1]); 
////fprintf(stderr, "%u, r_num = %u, r_sum[r_num-1] = %u\n", __LINE__, r_num, r_sum[r_num-1]); 
////fprintf(stderr, "%u, l_max_i = %u, l_max_sc = %u\n", __LINE__, l_max_i, l_max_sc); 
////fprintf(stderr, "%u, r_max_i = %u, r_max_sc = %u\n", __LINE__, r_max_i, r_max_sc); 
     
////fprintf(stderr, "%u, sub_err = %d, l_indel = %d, r_indel = %u\n", __LINE__, sub_err, l_indel, r_indel); 
    //++++++++++++++++++++++++++++++++
    int delta = sub->delta;
    int err_sum = sub->err_sum[0];
    int thres_score;
    if(query->b0 >0) {
        thres_score = query->b0; 
    } else{
        thres_score = query->candi_thres; 
    }

    int ins_num, del_num, len_indel;
    if(bg_flg > 0) { 
        len_indel = bg_flg; 
    } else if(indel_flg > 0){
        len_indel = indel_flg;
    } else {
        len_indel = 0; 
    }

    int err_num;
    int l_indel_sc = query->l_seq - l_indel;
    int r_indel_sc = query->l_seq - r_indel;
    int l_indel_q = 0, l_indel_qi = 0; 
    if(l_indel < sub_err && l_indel <= r_indel){
        err_num = 0;
        for(i = 1; i < len_indel; ++i) {
            ins_num = 0;    
            del_num = 0;
            for(j = 0; j < l_max_i; ++j){
                k =  j*12 + l_bg;
                pos_i = pos + k + i; 
                i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
                if(i_bwt >= target_idx[k][0] && i_bwt < target_idx[k][1] ){
                    ins_num++;
                    //ins_buf[i] = 1;  
                } else{
                    //ins_buf[i] = 0; 
                }
                //--------------------------------
                k =  j*12 + l_bg;
                pos_i = pos + k - i; 
                i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
 
                if(i_bwt >= target_idx[k][0] && i_bwt < target_idx[k][1] ){
                    del_num++;  
                    //del_buf[i] = 1;  
                } else{
                    //del_buf[i] = 0; 
                }
            }
            l_indel_qi = 0;
            int n_lsum = l_sum[l_max_i-1];
            if(del_num > n_lsum|| ins_num > n_lsum) {
                if(ins_num >= del_num) {
                    l_indel_q = ins_num; 
                    l_indel_qi = i;
                } else{
                    l_indel_q = del_num; 
                    l_indel_qi = -i;
                }
                break; 
            }         
        }//end  for(i = 1; i < len_indel; ++i) +++ 
        int l;
        if(l_indel_qi >0) l = 2*l_indel_qi;
        else l = -l_indel_qi;
        cur_sc = 6 + l + 5*(l_max_i - l_indel_q ); 
        l_indel_sc = query->l_seq - cur_sc - err_sum; 
    }
    int r_indel_q = 0, r_indel_qi = 0; 
    if(r_indel < sub_err && r_indel <= l_indel){
        err_num = 0;
        for(i = 1; i < len_indel; ++i) {
            ins_num = 0;
            del_num = 0;
            for(j = r_max_i+1; j < r_num; ++j){ //????
                k =  j*12 + r_bg;
                pos_i = pos + k + i; 
                i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
                if(i_bwt >= target_idx[k][0] && i_bwt < target_idx[k][1] ){
                    ins_num++;
                    ins_buf[i] = 1;  
                } else{
                    ins_buf[i] = 0; 
                }
                //--------------------------------
                k =  j*12 + r_bg;
                pos_i = pos + k - i; 
                i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);
                if(i_bwt >= target_idx[k][0] && i_bwt < target_idx[k][1] ){
                    del_num++;  
                    del_buf[i] = 1;  
                } else{
                    del_buf[i] = 0; 
                }
            }
////fprintf(stderr, "%u, del_num = %u, ins_num = %u\n", __LINE__, del_num, ins_num);
            int n_rsum = r_sum[r_num -1] - r_sum[r_max_i];  
            r_indel_qi = 0;
            if(del_num > n_rsum || ins_num > n_rsum) {//????
                if(ins_num >= del_num) {
                    r_indel_q = ins_num; 
                    r_indel_qi = i;
                } else{
                    r_indel_q = del_num;
                    r_indel_qi = -i;
                }
                break; 
            }         
        }//end  for(i = 1; i < len_indel; ++i) +++
        int l;
        if(r_indel_qi > 0) l = r_indel_qi;
        else l = -2*r_indel_qi; 
        cur_sc = 6 + l + 5*(r_num - r_max_i -1 - r_indel_q); 
        r_indel_sc = query->l_seq - cur_sc - err_sum;

////fprintf(stderr, "%u, r_indel_qi = %d, r_indel_q = %d\n", __LINE__, r_indel_qi, r_indel_q); 
    }


    int sub_score = query->l_seq - sub_err - err_sum;
////fprintf(stderr, "%u, sub_score = %d, l_indel_sc = %u, r_indel_sc = %d\n", __LINE__, sub_score, l_indel_sc, r_indel_sc); 
//fprintf(stderr, "%u, sub_err = %d, l_indel = %u, r_indel = %d\n", __LINE__, sub_err, l_indel, r_indel); 
//fprintf(stderr, "%u, l_indel_qi = %d, r_indel_qi = %d\n", __LINE__, l_indel_qi, r_indel_qi); 
    err_out[3] = 0;
    err_out[4] = 0;
    if(l_indel < sub_err && l_indel <= r_indel){
        if(l_indel < r_indel) {
            if(l_indel_qi != 0){ 
                err_out[0] = l_indel_sc;
                err_out[1] = 1;
                err_out[2] = l_max_i;
                err_out[2] += (l_indel_qi+32) * 256; 
                cur_sc = l_indel_sc;
            } else{
                err_out[0] = query->l_seq - err_sum - l_indel; 
                err_out[1] = 0; 
                err_out[2] = 0; 
                err_out[3] = 3; 
            }
        } else if(l_indel == r_indel){
            if(l_indel_qi != 0 && r_indel_qi != 0){
                if(l_indel_sc >= r_indel_sc) {
                    err_out[0] = l_indel_sc; 
                } else {
                    err_out[0] = r_indel_sc; 
                }
                err_out[1] = 0; 
                err_out[2] = 0; 
                err_out[3] = 3; 
            } else if(l_indel_qi != 0){ 
                err_out[0] = l_indel_sc;
                err_out[1] = 1;
                err_out[2] = l_max_i;
                err_out[2] += (l_indel_qi+32) * 256; 
                cur_sc = l_indel_sc;
            } else if(r_indel_qi != 0){ 
                err_out[0] = r_indel_sc;
                err_out[1] = 1;
                err_out[2] = r_max_i;
                err_out[2] += (r_indel_qi+32) * 256; 
                cur_sc = r_indel_sc;
            } else {
                err_out[0] = query->l_seq - err_sum - l_indel; 
                err_out[1] = 0; 
                err_out[2] = 0; 
                err_out[3] = 3; 
            }
        } 
    } else if(r_indel < sub_err && r_indel < l_indel) {
        if(r_indel_qi != 0) {
            err_out[0] = r_indel_sc;           
            err_out[1] = 2;
            err_out[2] = r_max_i; 
            err_out[2] += (r_indel_qi+32) * 256; 
            cur_sc = r_indel_sc;
        } else{
            err_out[0] = query->l_seq - err_sum - r_indel; 
            err_out[1] = 0; 
            err_out[2] = 0; 
            err_out[3] = 3; 
        }

    } else {
        err_out[0] = sub_score;
        err_out[1] = 0;
        err_out[2] = 0;
        cur_sc = sub_score;
    }

    //---------------------------
    if(cur_sc > query->b0 && err_out[3] == 0) {
        err_out[3] = 2;
    } else if(cur_sc + delta > query->b0 && err_out[3] == 0) {
        err_out[3] = 1;
    } 
    if(err_out[3] > 0) {
        for(i = 0; i < l_num; ++i){
            if(l_sc[i] > 0) {
                err_out[4] |= (uint32_t)1 << (31-i); 
            }
        } 
        for(i = 0; i < r_num; ++i){
            if(r_sc[i] > 0) {
                err_out[4] |= (uint32_t)1 << (15-i); 
            }
        } 
    } 
fprintf(stderr, "%u, l_num = %d, r_num = %d\n", __LINE__, l_num, r_num);
    return l_sum[l_num-1] + r_sum[r_num-1];   
}


int aln_pos_map(idx_t *fm_idx,  query_t *query, uint32_t pos_buf[],  int seq_off[], int indel_num, int (*err_num)[3], struct SubBuf *sub)
{
   
    uint32_t (*target_idx)[3] = query->target_idx; 
    //uint8_t *read_seq = query->read_seq;

    uint8_t *read_seq = query->is_rev?query->rseq:query->seq;
    uint32_t i, j, pos, pos_i, i_bwt,l_flg[33], seq12;
    uint32_t *isa = fm_idx->bwt->isa;
    int k;
    //---------------------------
    int seq_bg = seq_off[0];
    int seq_ed = seq_off[1];
    int bg = query->l_seq%12/2;
    int ed = query->l_seq/12;
    int num = query->l_seq/12;
    int j_ed = indel_num*2+1;

    for(i = 0; i <= j_ed; ++i){
        err_num[i][0] = 0;
        err_num[i][1] = 0;
        err_num[i][2] = 0;
        l_flg[i] = 0;
    }

    //pos = bwt_sa(fm_idx->bwt, idx[0]) - seq_bg;
    pos = pos_buf[0];
    //idx[1] = pos;
    for(i = 0; i < ed; ++i) {
        for(j = 0; j < j_ed; j++){
            pos_i = pos + i * 12 + bg + j -indel_num; 
          
            /*  
            if(pos_i >= 0 && pos_i < fm_idx->bns->l_pac ){
                i_bwt = isa[pos_i];
            } else i_bwt = 0;
            */
            i_bwt = bwt_get_idx(fm_idx->bwt, pos_i);


            if(i_bwt >= target_idx[i][0] && i_bwt < target_idx[i][1] ){
                l_flg[j] = 1;
                err_num[j][2] = 0;
            } else{
                err_num[j][0]++;
                if(l_flg[j] == 0) err_num[j][1] = i;
                err_num[j][2]++;
                k = i * 12 + bg; 
                if(j == indel_num) ++err_num[indel_num*2+1][0];
                if(k <= seq_bg && j == indel_num) ++err_num[indel_num*2+1][1];
                if(k >= seq_ed && j == indel_num) ++err_num[indel_num*2+1][2];
            }
        }
    }

    return err_num[indel_num][0];
}
int aln_long_seed(idx_t *fm_idx, uint32_t hash_boundry[], struct ExtBlck*eB, query_t *query, seed_t *seed, struct SubBuf *sub)
{
    int i, j, aln_ret = 0;
    uint32_t idx[2][3] = {{0, 0, 0}, {0, 0, 0}};
    
    int st_p[2] = {-1, -1};
    int aln_len = 0;
    int len[2], l_off, r_off, s_off, sid, s_id, s_len, ext_num; 
    len[0] = sub->idx_for[LEN_READ][0]; 
    len[1] = sub->idx_rev[LEN_READ][0]; 
    uint32_t seq_off[2], bg, ed, num[2];
    int sn = query->seed_num, l_seq = query->l_seq;
    uint32_t pos;
    int trm_r[2];
    trm_r[0] = sub->trm_r;
    trm_r[1] = sub->trm_r;
    int rev;
    if(len[0] >= l_seq - 2*sub->sub_err || 
            len[1] >= l_seq - 2*sub->sub_err){
        printf("%u, query->name = %s, len[0] = %d, len[1] = %d, trm_r = %d, trm_l = %d\n", 
                __LINE__, query->name, len[0], len[1], sub->trm_r, sub->trm_l);
        exit(1);  
    }
    uint32_t(*idx_tmp)[2];
    for(rev = 0; rev < 2; ++rev) {        
        idx_tmp = rev?sub->idx_rev:sub->idx_for;
        if(len[rev] <= SEED_LEN) continue;
        query->is_rev = rev;
        seq_off[1] = l_seq - trm_r[rev];   
        seq_off[0] = seq_off[1] - len[rev];
        bg = idx_tmp[seq_off[0]][0]; 
        ed = idx_tmp[seq_off[0]][1]; 
        num[rev] = ed + 1 - bg;
        if(num[rev] == 0) {
            printf("%u, error, num == 0\n", __LINE__);
            exit(1); 
        } 
fprintf(stderr, "%u, len[rev] = %d, num[rev] = %d\n", __LINE__, len[rev], num[rev]);
         
        if(len[rev] >= 2*l_seq/3 && num[rev] <= SW_THRES){
            int bg_r = l_seq - trm_r[rev] - SEED_LEN;
            int ed_r = l_seq - trm_r[rev] - len[rev];
            uint32_t n = 0;
            for(i = bg_r; i > ed_r; --i){
                bg = idx_tmp[i][0]; 
                ed = idx_tmp[i][1];  
                n = ed + 1 - bg;
                if(n <= SW_THRES) {
                    st_p[rev] = i;
                    idx[rev][0] = bg;
                    idx[rev][1] = ed; 
                    break;
                }
            }   
            i = SEED_LEN + 5;
            if(st_p[rev] < i) {
                st_p[rev] = i; 
                bg = idx_tmp[i][0]; 
                ed = idx_tmp[i][1];  
            } 
            len[rev] = l_seq - trm_r[rev] - st_p[rev];   
        } else {
            i = l_seq - trm_r[rev] - len[rev];
            st_p[rev] = i;
            idx[rev][0] = idx_tmp[i][0];
            idx[rev][1] = idx_tmp[i][1];  
        }
        seq_off[1] = l_seq - trm_r[rev];   
        seq_off[0] = seq_off[1] - len[rev];
        bg = idx_tmp[seq_off[0]][0]; 
        ed = idx_tmp[seq_off[0]][1]; 
        num[rev] = ed + 1 - bg;
fprintf(stderr, "%u, rev = %d, len[rev] = %d, trm_r[rev] = %d, num[rev] = %d\n", __LINE__, rev, len[rev], trm_r[rev], num[rev]);
        if(num[rev] > IS_SMLSIZ){
            uint32_t bwt_idx[4];
            int len_buf[2];    
            uint32_t k, l, n;
            bwt_idx[0] = bg;
            bwt_idx[1] = ed;
            len_buf[0] = len[rev];
            len_buf[1] = trm_r[rev];
            int l0 = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);
            if(l0 > 0) {
                k = bwt_idx[2];
                l = bwt_idx[3];
                n = l+1-k;
                seq_off[1] += l0;
                bg = k;
                ed = l;
                num[rev] = l + 1 - k;
                len[rev] += l0;
                trm_r[rev] -= l0; 
            }

fprintf(stderr, "%u, rev = %d, len[rev] = %d, trm_r[rev] = %d\n", __LINE__, rev, len[rev], trm_r[rev]);
        }
        if(len[rev] >= WOREST_PERCENT){
fprintf(stderr, "%u, len[rev] = %d, num[rev] = %d\n", __LINE__, len[rev], num[rev]);
fprintf(stderr, "%u, bg = %u, ed = %u, rev = %d\n", __LINE__, bg, ed, rev);
            uint32_t ix, pos;
            uint32_t *pos_buf = sub->pos_buf;
            //trm_r = sub->trm_r;
            uint32_t pos_i = 0, pi;

            if(ed +1 - bg > SW_THRES) ed = bg + SW_THRES;
            for(ix=bg; ix <= ed; ix++){
                pos = bwt_sa(fm_idx->bwt, ix);
                pos_buf[pos_i] = pos>seq_off[0]?pos-seq_off[0]:0;
                pos_i++;  
            }
           
fprintf(stderr, "%u, pos_i = %u\n", __LINE__, pos_i);
            if(pos_i > THRES_ALN_OUT*2) {

                l_off = l_seq - trm_r[rev] - len[rev];
                seed->id = 2*sn;
                for(sid = rev*sn; sid < sn+rev*sn; ++sid){
fprintf(stderr, "%u, sid = %d, num = %u, s_off = %d\n", __LINE__, sid, seed->slc[sid].num, seed->slc[sid].s_off);
                    if(seed->slc[sid].num == 0) continue;
                    if(seed->slc[sid].s_off > l_off) {
                        l_off = seed->slc[sid].s_off;
                        seed->id = sid;
                        break; 
                    }    
                }
                //seed->id = s_id;
fprintf(stderr, "%u, l_off = %u, seed->id = %d, len[rev] = %d, trm_r[rev] = %d\n", __LINE__, l_off, seed->id, len[rev], trm_r[rev]);
                pos_i = classify_pos_0(fm_idx, query, seed, pos_buf, pos_i, seq_off, sub);
            }
fprintf(stderr, "%u, pos_i = %u, seq_off[0] = %d, seq_off[1] = %d\n", __LINE__, pos_i, seq_off[0], seq_off[1]);
            query->l_seq -= trm_r[rev]; 
            uint8_t c;
            uint8_t *read_seq = query->is_rev?query->rseq:query->seq; 
            for(pi = 0; pi < pos_i; ++pi){   
                int sc = aln_pos_L(fm_idx, query, pos_buf+pi, seq_off, sub);
                
fprintf(stderr, "%u, pos = %u, sc = %d, trm_r = %d, l_seq = %d, query->l_seq = %d\n", __LINE__, pos_buf[pi], sc, trm_r[rev], l_seq, query->l_seq);
                int ed_pos = l_seq - trm_r[rev];
                for(j = l_seq - trm_r[rev]; j < l_seq; ++j){
                    c = __get_pac(fm_idx->pac, pos_buf[pi]+j);
                    if(c == read_seq[j]) {
                        ed_pos = j;
                        sc++; 
                    } else {
                        break;
                    }
                }
                if(sc >= l_seq) {
                    printf("%u, pos = %u, sc = %d\n", __LINE__, pos_buf[pi], sc);
                    exit(1); 
                }
                
                sub->aln_out->R_ti = MAX_CLIP + ed_pos;
                sub->aln_out->R_qi = ed_pos;
                rec_aln_info(query, pos_buf[pi], sc, sub);
            }
            query->l_seq += trm_r[rev]; 
fprintf(stderr, "%u, b0 = %d\n", __LINE__, query->b0);
            continue;
        }// end if(len[rev] >= 2*query->l_seq/3 || num[rev] <= SW_THRES)+++++
    }

    fprintf(stderr, "%u, query->b0 = %d\n", __LINE__, query->b0);
    aln_len = len[0]>len[1]?len[0]:len[1];
    if(query->b0 > l_seq - sub->sub_err) {
        return aln_len;
    }
    int ext_cls = 2;
    int ext_flg = SEED_LEN + 16*2*ext_cls; 
    if(ext_cls){
        for(rev = 0; rev < 2; ++rev) { 
            if(len[rev] > ext_flg && num[rev] < SW_THRES) {
                for(sid = rev*sn ; sid < sn + rev*sn; ++sid) {
                    s_off = seed->slc[sid].s_off;
                    ext_num = seed->slc[sid].ext_num;
                    l_off = s_off - 16 * ext_num;
                    r_off = s_off + SEED_LEN + 16 * ext_num;
                    if(l_off >= l_seq - len[rev] && r_off <= l_seq) {
                        seed->slc[sid].ext_num = 0; 
                    }
                } 
            }
        }
    }// end if(ext_cls)++++
    //-------------------------------------------- 
    int m_off, m_off0 = l_seq - l_seq/3 - SEED_LEN;          
    int flg[2] = {0, 0};
    int l, ms_off[2] = {0,0}, ms_i[2] = {0, 0};
    for(rev = 0; rev < 2; ++rev) { 
        m_off = st_p[rev]; 
        if(m_off > m_off0) m_off = m_off0;
        for(sid = rev*sn; sid < sn+rev*sn; ++sid){
            if(seed->slc[sid].s_off < m_off) {
                if(seed->slc[sid].s_off > ms_off[rev]) {
                    ms_off[rev] = seed->slc[sid].s_off;
                    ms_i[rev] = sid; 
                }    
                if(seed->slc[sid].num < 1) { 
                    flg[rev] = 1;
                    break;  
                } else {
                    l = seed->slc[sid].len; 
                    if(l < 16 && l < seed->slc[sid].s_off) {
                        flg[rev] = 1;
                        break; 
                    }
                }
            }                     
        }
    } 
    uint32_t (*idx_buf)[17][3] = seed->idx_buf; 
    int st_pos[2];

    for(rev = 0; rev < 2; ++rev) { 
        query->fr_flg[rev] = 0;
        st_pos[rev] = seed->slc[ms_i[rev]].s_off + SEED_LEN; 
        if(flg[rev] == 0) {
            query->is_rev = rev;
            seed->id = ms_i[rev];
            s_off = seed->slc[ms_i[rev]].s_off;
            s_len = seed->slc[ms_i[rev]].len;
            seq_off[0] = s_off + 8-s_len;
            seq_off[1] = s_off + SEED_LEN;
            uint32_t idx[2], bg, ed;
            idx[0] = idx_buf[ms_i[rev]][s_len][0];
            idx[1] = idx_buf[ms_i[rev]][s_len][1];
fprintf(stderr, "%u, seq_off[0] = %u, seq_off[1] = %u\n", __LINE__, seq_off[0], seq_off[1]);
fprintf(stderr, "%u, idx[0] = %u, idx[1] = %u\n", __LINE__, idx[0], idx[1]);
            len[rev] = bwt_exact_aln1(fm_idx, hash_boundry, query, 
                                        idx, seq_off, sub);

fprintf(stderr, "%u, idx[0] = %u, idx[1] = %u, len[rev] = %d, rev = %d\n", __LINE__, idx[0], idx[1], len[rev], rev);
            uint32_t(*idx_tmp)[2] = rev?sub->idx_rev:sub->idx_for; 
            bg = idx_tmp[seq_off[1] - len[rev]][0];
            ed = idx_tmp[seq_off[1] - len[rev]][1];
            /*  
            if(rev == 0) {
                bg = sub->idx_for[seq_off[1] - len[rev]][0]; 
                ed = sub->idx_for[seq_off[1] - len[rev]][1]; 
            } else {
                bg = sub->idx_rev[seq_off[1] - len[rev]][0]; 
                ed = sub->idx_rev[seq_off[1] - len[rev]][1]; 
            }
            */
            
fprintf(stderr, "%u, bg = %u, ed = %u\n", __LINE__,bg, ed);
            num[rev] = ed + 1 - bg;
            if(num[rev] > IS_SMLSIZ) {
                uint32_t bwt_idx[4];
                bwt_idx[0] = bg;
                bwt_idx[1] = ed;
                int len_buf[2];
                len_buf[0] = len[rev];
                len_buf[1] = 32;
                l = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);
                if(l > 0) {
                    bg = bwt_idx[2];
                    ed = bwt_idx[3];
                    num[rev] = ed+1-bg;
                    seq_off[1] += l;
                    len[rev] += l;
                }
            }

            if(num[rev] == 0 ) {
                printf("%u, error!!!!\n", __LINE__);
                exit(1); 
            } 
            //if(aln_len < len[rev]) { aln_len = len[rev];}
            if(len[rev] > 2*SEED_LEN && num[rev] < SW_THRES) { 
                if(num[rev] < IS_SMLSIZ || len[rev] == seq_off[1]) {
                    for(i =0; i < num[rev]; ++i){
                        pos = bwt_sa(fm_idx->bwt, bg+i) 
                                       - (seq_off[1] - len[rev]);
                        sub->pos_buf[0] = pos;
                        int is_new = AlgnPos(fm_idx, query, 1, sub);
                        set_sw_ed_pos(pos, query, sub);
                        if(query->b0 > l_seq - sub->delta && is_new > 0) {
                            query->fr_flg[rev] = query->b0; 
                        }

                    } 
                } else {
                    for(i =0; i < num[rev]; ++i){
                        pos = bwt_sa(fm_idx->bwt, bg+i) 
                                        - (seq_off[1] - len[rev]);
                        get_set_sw_to_pos(pos, 0, 2, seed, query, sub);
                    }                               
                }
            }                     
        }
    } 
    aln_len = len[0] > len[1]?len[0]:len[1]; 
/*  
    if(query->b0 > l_seq - sub->delta) {
        return aln_len;
    } 
*/
    if(ext_cls){
        for(rev = 0; rev < 2; ++rev) {
            if(len[rev] > ext_flg && num[rev] < SW_THRES) {
                for(sid = rev*sn ; sid < sn + rev*sn; ++sid) {
                    s_off = seed->slc[sid].s_off;
                    ext_num = seed->slc[sid].ext_num;
                    l_off = s_off - 16 * ext_num;
                    r_off = s_off + SEED_LEN + 16 * ext_num;
                    if(l_off >= st_pos[rev] - len[rev] && r_off <= st_pos[rev]) {
                        seed->slc[sid].ext_num = 0; 
                    }
                } 
            }
        }
    }

    return aln_len;
}

int aln_long_exact(idx_t *fm_idx, uint32_t hash_boundry[], struct ExtBlck*eB, query_t *query, seed_t *seed, struct SubBuf *sub)
{
    int aln_flg = 0, n_fnd = 0;
    int l_seq = query->l_seq; 
    int sub_err = sub->sub_err;
    int trm_r = sub->trm_r, trm_l = sub->trm_l;
    int len[2] = {0, 0};  
    int rev;
    uint32_t bgn, end, num;
    uint32_t (*idx_buf)[2];
    int i, j;
    for(rev = 0; rev < 2; ++rev){
        query->is_rev = rev;       
        uint8_t *read_seq = rev?query->rseq:query->seq;
        len[rev] = bwt_exact_aln(fm_idx, hash_boundry, query, sub);

        
fprintf(stderr, "%u, len[rev] = %d\n", __LINE__, len[rev]);
        if(len[rev] >= l_seq - trm_l - trm_r){
            query->fr_flg[rev] = len[rev];
            uint32_t bg_p = l_seq - trm_r -len[rev];
            idx_buf = rev?sub->idx_rev:sub->idx_for;
            bgn = idx_buf[bg_p][0];    
            end = idx_buf[bg_p][1];  
            num = end + 1 - bgn;  
fprintf(stderr, "%u, rev = %d, num = %u\n", __LINE__, rev, num);
            if(num >= THRES_ALN_OUT) {
                aln_flg = 2;
            } else {
                if(aln_flg == 0) aln_flg = 1;
            }
        } else {
            query->fr_flg[rev] = 0;
        } 
    }
    if(aln_flg == 0) return aln_flg;
fprintf(stderr, "%u, aln_flg = %d\n", __LINE__, aln_flg); 
    if(aln_flg == 2) {
        for(rev = 0; rev < 2; ++rev){
            query->is_rev = rev;
            uint8_t *read_seq = rev?query->rseq:query->seq;
            if(query->fr_flg[rev] == 0) continue;
            idx_buf = rev?sub->idx_rev:sub->idx_for;


            int bg_p = l_seq - trm_r - len[rev];
            bgn = idx_buf[bg_p][0]; 
            end = idx_buf[bg_p][1];
            num = end + 1 - bgn; 
fprintf(stderr, "%u, bg_p = %d, len[rev] = %d, trm_r\n", __LINE__, len[rev]); 
            //if(num > THRES_ALN_OUT*10){
            if(num > SW_THRES){
                uint32_t k = bgn, l = end, n = num;
                //int l0 = 0;
                for(i = bg_p-1; i >=0; --i) {
                    n = bwt_match_exact_alt(fm_idx->bwt, 1, read_seq+i, &k, &l);
                    if(n <= THRES_ALN_OUT) break;
                    idx_buf[i][0] = k;
                    idx_buf[i][1] = l;
                    len[rev]++;
                    //l0++;  
                }  
           
                bgn = k;
                end = l;
                num = l + 1 - k;
                //len[rev] += l0; 
            
                if(num > SW_THRES) {
                //if(num > THRES_ALN_OUT*10){
                    uint32_t bwt_idx[4];
                    bwt_idx[0] = bgn;
                    bwt_idx[1] = end;
                    int len_buf[2], seq_off[2];
                    len_buf[0] = len[rev];
                    len_buf[1] = trm_r;
                   
                    seq_off[1] = query->l_seq - trm_r;
                    seq_off[0] = seq_off[1] - len[rev];
                    int l0 = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);
                    if(l0 > 0) {
                        k = bwt_idx[2];
                        l = bwt_idx[3];
                        n = l+1-k;
                        seq_off[1] += l0;
                        if(n > THRES_ALN_OUT) {
                            bgn = k;
                            end = l;
                            num = l + 1 - k;
                            len[rev] += l0;
                            trm_r -= l0;
                            sub->trm_r = trm_r;
                           
                        }
                    }
                }
            } // end if(num > THRES_ALN_OUT*10)++++++
            
            n_fnd = rec_bwt_exact(fm_idx, query, len[rev], sub);    

            aln_flg = 3;
        }
    }
    if(aln_flg == 3) {
        return aln_flg; 
    } 
    int end_p[2] = {0, 0};  
    int pi;
    int n_f[2] = {0, 0};
    
fprintf(stderr, "%u, b0 = %d\n", __LINE__, query->b0);
    int rev_i, rev_flg = 0;
    if(query->fr_flg[0] == 0) {
        rev_flg = 1; 
    }
    for(rev_i = 0; rev_i < 2; ++rev_i){
        rev = (rev_i + rev_flg)%2;
        query->is_rev = rev;
        uint8_t *read_seq = rev?query->rseq:query->seq;
        uint32_t b = 1, e = 0, n = 0;
        uint32_t p_R = l_seq - sub->sub_err;
        uint32_t p_L = l_seq - trm_r - len[rev];

        idx_buf = rev?sub->idx_rev:sub->idx_for;
        if(query->fr_flg[rev] == 0){ 
            int b0 = query->b0;
            int n_f = 0; 
            int (*found)[4] = sub->aln_out->found; 
            int n_fnd  = 0;
            for(i =0; i < sub->delta; ++i) {
                n_fnd += found[i][1];
                if(found[i][0] >= l_seq - sub->sub_err && n_fnd >= THRES_ALN_OUT){
                    break; 
                } else if(found[i][0] > l_seq - sub->sub_err) {
                    break;  
                } 
            }
            if(n_fnd >= THRES_ALN_OUT) {
                continue;
            }
            if(p_L >= SEED_LEN) {
                query->l_seq = p_L;   
                int tm_r = sub->trm_r, tm_l = sub->trm_l;
                sub->trm_l = 0;
                sub->trm_r = 0;
                int l0 = bwt_exact_aln(fm_idx, hash_boundry, query, sub);
                sub->trm_l = tm_l;
                sub->trm_r = tm_r;
                query->l_seq = l_seq;
                if(len[rev] + l0 < l_seq -tm_r - tm_l){
                    continue; 
                } else {
                    query->fr_flg[rev] = len[rev]; 
                }
            }   
        }
        
      
   
        //int idx = l_seq - trm_r - len[rev];
   
        for(pi = p_R - SEED_LEN; pi > p_L; --pi){
            n = idx_buf[pi][1] + 1 - idx_buf[pi][0]; 
            if(n <= 10*THRES_ALN_OUT) {
                break;
            }
        }
fprintf(stderr, "%u, pi = %d, p_L = %d, rev = %d, n = %u\n", __LINE__, pi, p_L, rev, n);
        if(pi > p_L) {
            b = idx_buf[pi][0]; 
            e = idx_buf[pi][1]; 
            n = e + 1 -b;
            if(n < THRES_ALN_OUT && pi < l_seq - sub->sub_err) {
                b = idx_buf[pi+1][0]; 
                e = idx_buf[pi+1][1];
                pi += 1; 
            }                 
        } else {
            pi = p_L;
            b = idx_buf[pi][0]; 
            e = idx_buf[pi][1]; 
            n = e + 1 -b;

        }
        end_p[rev] = 0; 
        if(pi >= SEED_LEN) {
        
            end_p[rev] = pi - 1; 
//end_p[rev] = SEED_LEN + sub->trm_l;
        }        
        uint32_t idx;
        uint32_t seq_off[2];
  
        //seq_off[0] = p_L;
        //seq_off[1] = pi; 
        
        seq_off[0] = pi;
        seq_off[1] = p_R; 


fprintf(stderr, "%u, bg = %u, ed = %u, seq_off[0] = %d, seq_off[1] = %d, pi = %d, p_R = %d\n", __LINE__, b, e, seq_off[0], seq_off[1], pi, p_R);
        for(idx = b; idx <=e; ++idx){
            uint32_t pos = bwt_sa(fm_idx->bwt, idx) - seq_off[0]; 
            uint32_t pos_buf[2];
            pos_buf[0] = pos; 
            
            int flg_ed =  get_sw_ed_val(pos, sub, query);
            if(flg_ed > 0) continue;
            int sc = aln_pos_L_0(fm_idx, query, pos_buf, seq_off, sub);
            
            fprintf(stderr, "%u, pos = %u, sc = %d, trm_r = %d\n", __LINE__, pos, sc, trm_r);
            
            set_sw_ed_pos(pos, query, sub);
            int ed_pos = query->l_seq - trm_r;
            for(j = l_seq - trm_r; j < l_seq; ++j){
                uint8_t c = __get_pac(fm_idx->pac, pos+j);
                if(c == read_seq[j]) {
                    ed_pos = j;
                    sc++; 
                } else {
                    break;
                }
            }
            sub->aln_out->R_ti = MAX_CLIP + ed_pos+1;
            sub->aln_out->R_qi = ed_pos+1;
            n_f[rev] = rec_aln_info(query, pos, sc, sub);  

        } 
    }

fprintf(stderr, "%u, b0 = %d\n", __LINE__, query->b0);
    int n_flg = n_f[0] > n_f[1] ?n_f[0]:n_f[1];
    if(n_flg >= THRES_ALN_OUT) {
        aln_flg = 2; 
    }
               
fprintf(stderr, "%u, n_flg = %d\n", __LINE__, n_flg);
    
    rev_flg = 0;
    if(end_p[0] < l_seq - sub->trm_r - sub->trm_l) {
        rev_flg = 1; 
    }

    for(rev_i = 0; rev_i < 2; ++rev_i){
fprintf(stderr, "%u, end_p[rev] = %d\n", __LINE__, end_p[rev]);
        rev = (rev_i + rev_flg)%2;
        if(end_p[rev] == 0) continue;
        
        if(query->fr_flg[rev] < l_seq - sub->trm_r - sub->trm_l){ 
            int b0 = query->b0;
            int n_f = 0; 
            int (*found)[4] = sub->aln_out->found; 
            int n_fnd  = 0;
            for(i =0; i < sub->delta; ++i) {
                n_fnd += found[i][1];
                if(found[i][0] >= l_seq - sub->sub_err && n_fnd >= THRES_ALN_OUT){
                    break; 
                } else if(found[i][0] > l_seq - sub->sub_err) {
                    break;  
                } 
            }
            if(n_fnd >= THRES_ALN_OUT) {
                continue;
            }
        }

        query->is_rev = rev;
        idx_buf = rev?sub->idx_rev:sub->idx_for;
       
        uint8_t *read_seq = rev?query->rseq:query->seq;
        uint32_t b, e, n;
        int l_seq = query->l_seq; 
        query->l_seq = end_p[rev];   
fprintf(stderr, "%u\n", __LINE__);
        int tm_r = sub->trm_r, tm_l = sub->trm_l;
        sub->trm_l = 0;
        sub->trm_r = 0;
        len[rev] = bwt_exact_aln(fm_idx, hash_boundry, query, sub);
fprintf(stderr, "%u, rev = %d, len[rev] = %d, end_p[rev] = %d\n", __LINE__, rev, len[rev], end_p[rev]);
        sub->trm_l = tm_l;
        sub->trm_r = tm_r;
        query->l_seq = l_seq;
        
        //if(len[rev] == query->l_seq) {
        if( end_p[rev] - len[rev] <= trm_l) {
            b = idx_buf[trm_l-1][0];    
            e = idx_buf[trm_l-1][1];  
            /*   
            if(rev == 0) {
                b = sub->idx_for[trm_l-1][0]; 
                e = sub->idx_for[trm_l-1][1]; 
            } else {
                b = sub->idx_rev[trm_l-1][0]; 
                e = sub->idx_rev[trm_l-1][1]; 
            }
            */
            len[rev] = end_p[rev] - trm_l;
        } else {
            continue; 
        }            
fprintf(stderr, "%u, b = %u, e = %u\n", __LINE__, b, e);

        uint32_t idx;
        uint32_t seq_off[2];
        seq_off[0] = end_p[rev]-len[rev];
        seq_off[1] = end_p[rev];
        for(idx = b; idx <=e; ++idx){
            uint32_t pos = bwt_sa(fm_idx->bwt, idx) - seq_off[0]; 
            uint32_t pos_buf[2];
            pos_buf[0] = pos; 
            int flg_ed =  get_sw_ed_val(pos, sub, query);
            if(flg_ed > 0) continue;


            int sc = aln_pos_R_0(fm_idx, query, pos_buf, seq_off, sub);
fprintf(stderr, "%u, pos = %u, sc = %d\n", __LINE__, pos, sc);
            set_sw_ed_pos(pos, query, sub);
            
            int ed_pos = query->l_seq - trm_r;
           
            sub->aln_out->R_ti = MAX_CLIP + ed_pos;
            sub->aln_out->R_qi = ed_pos;
            n_f[rev] = rec_aln_info(query, pos, sc, sub);  

        } 
fprintf(stderr, "%u\n", __LINE__);
    }
    
fprintf(stderr, "%u\n", __LINE__);
    if(aln_flg < 2) {
        aln_flg = 2;      
    }
return aln_flg;
}

int aln_sw_seed_all(idx_t *fm_idx, uint32_t hash_boundry[], struct ExtBlck*eB, query_t *query, seed_t *seed, struct SubBuf *sub)
{
    int bg, ed; 
    int s_id;
    uint32_t bgn, end, num;
    int s_off[3], len = 0, L_off;
    uint32_t (*out)[2] = sub->out_buf;
    uint32_t pos_i = 0;
    uint32_t *pos_buf = sub->pos_buf;
    uint32_t idx_buf[2];
    int seed_size = query->seed_num*2;
    int *SLC_SWITCH = sub->SLC_SWITCH;
    int *ALN_SWITCH = sub->ALN_SWITCH;
    //+++++++++++++++++++++++++++++++++++++ 
    int overlap_num = 0;
    int SUCEED_FLAG = 0;
    seed->cls = 0;
    if(seed->n_for >= seed->n_rev) {
        bg = 0;
        ed = seed_size/2;
        query->is_rev = 0;
        overlap_num = seed->n_for; 
    } else{
        bg = seed_size/2;
        ed = seed_size; 
        query->is_rev = 1;
        overlap_num = seed->n_rev; 
    }
    if(overlap_num > 5) {
      sub->thres_sw_to = overlap_num/3;
    } else {
      sub->thres_sw_to = overlap_num/2; 
    } 
    sub->olp_flg = 1;
    //------------------------------
    int rep = 0;
    for(rep = 0; rep < 2; ++rep) { 
        for(s_id = bg; s_id < ed; s_id++) {
            bgn = seed->slc[s_id].bgn;
            num = seed->slc[s_id].num;
            end = seed->slc[s_id].end;  
            if(num == 0) continue;
            s_off[0] = seed->slc[s_id].s_off;
            s_off[1] = 16;
            L_off = s_off[0];
            pos_i = 0;
            if(s_id < query->seed_num) {
                query->is_rev = 0;
                query->read_seq = query->seq;
                query->target_idx = query->target_idx_f;
            } else {
                query->is_rev = 1;
                query->read_seq = query->rseq;
                query->target_idx = query->target_idx_r;
            }
            len = 0;
            //if(num > 1) {
            if(num > IS_SMLSIZ) {
                idx_buf[0] = bgn;
                idx_buf[1] = end;
                len = bwt_ext_kmer(fm_idx, query->read_seq, 
                                    idx_buf, s_off, out);
                if(len > 0) {
                    if(len > s_off[1]) {
                        printf("%u, error: len < s_off[1] !!\n", __LINE__); 
                        exit(1);
                    }
                     
                    bgn = out[len][0];
                    end = out[len][1];
                    num = end + 1 - bgn; 
                    L_off -= len;
                }
            }
            int seq_off[2];
            seq_off[0] = L_off;
            seq_off[1] = L_off + SEED_LEN + len;
            //if(num > IS_SMLSIZ) {
            if(0) {
                uint32_t bwt_idx[4];
                bwt_idx[0] = bgn;
                bwt_idx[1] = end;
                int len_buf[2];
                len_buf[0] = len+SEED_LEN;
                len_buf[1] = 32;
                int l = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);
                if(l > 0) {
                    bgn = bwt_idx[2];
                    end = bwt_idx[3];
                    num = end+1-bgn;
                    seq_off[1] += l;
                }
            }
            uint32_t ix, pos;
            pos_i = 0;
            for(ix=bgn; ix <= end; ix++){
                pos = bwt_sa(fm_idx->bwt, ix);
                pos_buf[pos_i] = pos>L_off?pos-L_off:0;
                pos_i++;  
            }
            seed->id = s_id;
            pos_i = classify_pos(fm_idx, query, seed, 
                                pos_buf, pos_i, seq_off, sub);
            if(pos_i == 0) continue;                       
            if(SLC_SWITCH[0]) pos_i = slc_aln_pos(fm_idx,  query, 
                                                    seq_off, pos_i, sub);
            if(pos_i > 0){
                AlgnPos(fm_idx, query, pos_i,  sub);
                /*  
                sub->err_sum[0] = 0;
                for(i = 0; i < pos_i; ++i) {
                    int sc = aln_pos_both(fm_idx,  query, pos_buf+i,  
                                        seq_off, sub);
                    rec_aln_info(query, pos_buf[i], sc, sub);
                }
                */
            }

            int max_err = 0;
            int query_err = query->l_seq - query->b0 + sub->delta;
            int seed_err = query->seed_num*sub->sub_err;
            int cur_err = s_id % query->seed_num; 
            if(query->is_rev == 0) seed_err -= seed->n_for*sub->sub_err;
            else seed_err -= seed->n_rev*sub->sub_err;
            max_err = query_err - seed_err - cur_err;
            if(max_err <= 0) {
                 if(query->rep_num > 0) SUCEED_FLAG =1;
                 break;
            }
        }
        if(seed->n_for >= seed->n_rev) {
            bg = seed_size/2;
            ed = seed_size;  
            query->is_rev = 1; 
        } else{
            bg = 0;
            ed = seed_size/2;  
            query->is_rev = 0; 
        }
        
    }// end for(rep =0; rep < 2; ++rep)++++++
    //seed_id = query->seed_num*2;
    return SUCEED_FLAG;
}//


int aln_sw_seed_some(idx_t *fm_idx, uint32_t hash_boundry[], struct ExtBlck*eB, query_t *query, seed_t *seed, struct SubBuf *sub)
{    
    int rep_num, sd, i, j, seed_i, L_off;
    //int rev_flg[2] = {0, 0};
    uint32_t bgn, end, num, pos_i, pos, idx_buf[2];
    int s_off[3], len = 0;
    uint32_t (* out)[2] = sub->out_buf;
    uint32_t seed_buf[LEN_READ];
    seed->cls = 0;
    j = 0; 
    int seed_id = seed->id;
    int seed_num = query->seed_num;
    uint32_t *pos_buf = sub->pos_buf;
   
    int seed_size = query->seed_num*2;
    int *SLC_SWITCH = sub->SLC_SWITCH;
    int *ALN_SWITCH = sub->ALN_SWITCH;
    //+++++++++++++++++++++++++++++++++++++ 
    int overlap_num = 0;
    int SUCEED_FLAG = 0;
 
    
    
    if(seed_id < seed_num) {
        for(i = 0; i < seed_num; ++i){
            if(seed->slc[i].num > 0 && seed->slc[i].num <= SW_THRES){
                seed_buf[j++] = i;
                seed->del[i] = 1;
            } else{
            
            } 
        }
        query->is_rev = 0;
        sub->thres_sw_to = query->n_olp_for;
    } else {
        for(i = seed_num; i < seed_num*2; ++i){
            if(seed->slc[i].num > 0 && seed->slc[i].num <= SW_THRES ){
                seed_buf[j++] = i;
                seed->del[i] = 1;
            } else{ } 
        }
        query->is_rev = 1; 
        sub->thres_sw_to = query->n_olp_rev;
    }
    rep_num = j; 
//cur_time = outUsedTime();
    //++++++++++++++++++++++++++++++++++
    sub->thres_sw_to = overlap_num/3;
    sub->thres_pos_num = IS_SMLSIZ/2; // num == 1
    sub->thres_sw_olp = 2;
    sub->olp_flg = 1;
    sub->err_sum[0] = 0;
    //----------------------------------
    for(sd = rep_num -1; sd >= 0; --sd) {
        seed_i = seed_buf[sd];
        //++rev_flg[seed_i/seed_num];
        num = seed->slc[seed_i].num; 
        bgn = seed->slc[seed_i].bgn; 
        end = seed->slc[seed_i].end; 

        pos_i = 0;
        idx_buf[0] = bgn;
        idx_buf[1] = end;
        s_off[0] = seed->slc[seed_i].s_off;
        s_off[1] = 16;
        L_off = s_off[0];
        if(seed_i < query->seed_num) {
            query->is_rev = 0;
            query->read_seq = query->seq;
            query->target_idx = query->target_idx_f;
        } else {
            query->is_rev = 1;
            query->read_seq = query->rseq;
            query->target_idx = query->target_idx_r;
        }
        len = 0; 
        //if(num > 1){
        if(num > IS_SMLSIZ){
            len = bwt_ext_kmer(fm_idx, query->read_seq, 
                                idx_buf, s_off, out);
            if(len > s_off[1]) {
                printf("%u, error: len < seed_off_L[1] !!\n", __LINE__); 
                exit(1);
            }
             
            bgn = out[len][0];
            end = out[len][1];
            num = end+1 - bgn;
            L_off -= len;
        }
        int seq_off[2];
        seq_off[0] = L_off;
        seq_off[1] = L_off + SEED_LEN + len;
        //if(num > IS_SMLSIZ) {
        if(0) {
            uint32_t bwt_idx[4];
            bwt_idx[0] = bgn;
            bwt_idx[1] = end;
            int len_buf[2];
            len_buf[0] = len+SEED_LEN;
            len_buf[1] = 32;
            int l = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);
            if(l > 0) {
                bgn = bwt_idx[2];
                end = bwt_idx[3];
                num = end+1-bgn;
                seq_off[1] += l;
            }
        }
        
        uint32_t ix, pos;
        pos_i = 0;
        for(ix=bgn; ix <= end; ix++){
            pos = bwt_sa(fm_idx->bwt, ix);
            pos_buf[pos_i] = pos>L_off ? pos-L_off:0;
            pos_i++; 
        }
        /*  
        int seq_off[2];
        seq_off[0] = L_off;
        seq_off[1] = s_off[0]+SEED_LEN;
        */
        seed->id = seed_i;
        pos_i = classify_pos(fm_idx, query, seed, 
                                pos_buf, pos_i, seq_off, sub);

        if(SLC_SWITCH[1]) pos_i = slc_aln_pos(fm_idx,  query, 
                                                seq_off, pos_i, sub);
        if(pos_i > 0) {
            if(ALN_SWITCH[1] > 0) {
                sub->err_sum[0] = 0;
                for(i = 0; i < pos_i; ++i) {
                    int sc = aln_pos_both(fm_idx,  query, pos_buf+i,  
                                        seq_off, sub);
                    rec_aln_info(query, pos_buf[i], sc, sub);
                }
            } else if(pos_i > 0){
                AlgnPos(fm_idx, query, pos_i,  sub);
            }
        }
    } // end for(sd = rep_num -1; sd >= 0; --sd)++++
    if(seed_id <  seed_num) seed->n_sm_f = 0;
    if(seed_id >= seed_num) seed->n_sm_r = 0;
    if(seed->slc[seed_id].num <= SW_THRES) {
        seed_id++;
        //continue;
    } 
    return SUCEED_FLAG;
}  
//pos_bk_buf数据处理
int aln_sw_bak_buf(idx_t *fm_idx, uint32_t hash_boundry[], struct ExtBlck*eB, query_t *query, seed_t *seed, struct SubBuf *sub)
{
    int aln_flg = 0;
    int THRES_BAK = sub->THRES_BAK;
    int sid, i, rev, sn = query->seed_num;
    uint32_t pos, pos_i = 0, j;
    for(rev = 0; rev < 2; ++rev){
        if(query->b0 > THRES_BAK) {return aln_flg;}
        query->is_rev = rev;
        uint8_t *seq = rev? query->rseq:query->seq;
        pos_i = 0;
        for(i = 0; i < sub->n_seed_bk[rev]; ++i) {
            pos = sub->pos_bk_buf[rev][i][0];
            int flg_ed =  get_sw_ed_val(pos, sub, query);
            if(flg_ed > 0) continue;
            sid = sub->pos_bk_buf[rev][i][1];
            if(sid == 2*sn) {
                query->is_rev = rev;
                sub->pos_buf[0] = pos;
                pos_i = 1; 
                pos_i = set_sw_ed_buf(sub->pos_buf, pos_i, query, sub);    
                AlgnPos(fm_idx, query, pos_i, sub); 
                pos_i = 0;
                continue;
            }
fprintf(stderr, "%u\n", __LINE__);           
            if(sub->seed_bk[sid] > IS_SMLSIZ) {
            //if(0) {
                uint32_t bwt_idx[4], seq_off[2], num;
                uint32_t bg = seed->slc[sid].bgn; 
                uint32_t ed = seed->slc[sid].end;
                seq_off[0] = seed->slc[sid].s_off;  
                seq_off[1] = seed->slc[sid].s_off+SEED_LEN;
                int len = SEED_LEN;  
                bwt_idx[0] = bg;
                bwt_idx[1] = ed;
                int len_buf[2];
                len_buf[0] = len;
                len_buf[1] = 32;
                int l0 = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);
                if(l0 > 0) {
                    seq_off[1] += l0;
                    num = bwt_idx[3] + 1 - bwt_idx[2];
                   
                    bg = bwt_idx[2];
                    ed = bwt_idx[3];
                   
                    int num0 = num;
                    if(num > IS_SMLSIZ) {
                        uint32_t k = bg, l = ed;
                        int n, len1 = 0;
                        for(i = 1; i <= seq_off[0]; ++i){ 
                            n = bwt_match_exact_alt(fm_idx->bwt, 1, 
                                                    seq+seq_off[0]-i, &k, &l);
                            if(n > 0) {
                                bg = k; 
                                ed = l;      
                                num = n;
                                len1++;
                            } else { break;}
                        }
                        seq_off[0] -= len1; 
                    }        
                    if(num < IS_SMLSIZ) {
                        for(j = bg; j <= ed; ++j) {
                            pos = bwt_sa(fm_idx->bwt, j) - seq_off[0];
                            flg_ed =  get_sw_ed_val(pos, sub, query);
                            if(flg_ed > 0) continue;
                            sub->pos_buf[pos_i++] = pos; 
                        } 
                        sub->seed_bk[sid] = 0;
                    } else if(seed->slc[sid].num < SW_THRES){
                        sub->seed_bk[sid] = -1; 
                    } else{
                        sub->seed_bk[sid] = 0; 
                    }
                } // end if(l > 0) +++++++++
                if(sub->seed_bk[sid] > 0) sub->seed_bk[sid] = -1; 
            } else if(sub->seed_bk[sid] != 0){
                sub->pos_buf[pos_i++] = pos;
                sub->seed_bk[sid]--;   
            } 
            if(pos_i > 0){
                pos_i = set_sw_ed_buf(sub->pos_buf, pos_i, query, sub);    
                AlgnPos1(fm_idx, query, pos_i, sub->pos_buf, sub);
                pos_i = 0;
            }
fprintf(stderr, "%u\n", __LINE__);           
       
        }// end for(i = 0; i < sub->n_seed_bk[rev]; ++i)+++
        
    }// end for(rev = 0; rev <2; ++rev)+++++
    
    return aln_flg;     
    //if(query->b0 < bk_thres){ 

   
    if(0){
        int THRES_DELTA = (query->l_seq - query->candi_thres)/2;  
        uint32_t (*pos_bk)[2] = sub->pos_bk_buf[0]; 
        uint32_t (*pos_bk1)[2] = sub->pos_bk_buf[1]; 
        int n_pos_bk0 = sub->n_seed_bk[0];
        int n_pos_bk1 = sub->n_seed_bk[1];
        int bk_thres = query->candi_thres + THRES_DELTA;
        uint32_t *pos_buf = sub->pos_buf;
        if(n_pos_bk0 > 0){
            query->is_rev = 0; 
            int pos_i = 0;
            for(pos_i =0; pos_i < n_pos_bk0; ++pos_i) {
                pos_buf[pos_i] = pos_bk[pos_i][0];
            } 
            pos_i = clean_hash_pos(fm_idx, pos_buf, pos_i, query, sub);
            pos_i = set_sw_ed_buf(pos_buf, pos_i, query, sub);      
            AlgnPos1(fm_idx, query, pos_i, pos_buf, sub);
        }                            
        if(n_pos_bk1 > 0){
            query->is_rev = 1; 
            int pos_i = 0;
            for(pos_i =0; pos_i < n_pos_bk1; ++pos_i) {
                pos_buf[pos_i] = pos_bk1[pos_i][0];
            } 
            pos_i = clean_hash_pos(fm_idx, pos_buf, pos_i, query, sub);
            pos_i = set_sw_ed_buf(pos_buf, pos_i, query, sub);      
            AlgnPos1(fm_idx, query, pos_i, pos_buf, sub);
        }
     
    }

}

//种子向右精确扩展 
int aln_seed_exact_R(idx_t *fm_idx, uint32_t hash_boundry[], struct ExtBlck*eB, query_t *query, seed_t *seed, struct SubBuf *sub)
{                        
    int ret_flg = 0;
    uint32_t (*idx_buf)[17][3] = seed->idx_buf; 
    uint32_t *pos_buf = sub->pos_buf;
    uint32_t ix, bg, ed, num, len, L_off, si, pos_i;
    for(si = 0; si < query->seed_num*2; ++si){
        if(si < query->seed_num) {
            query->is_rev = 0;
        } else {
            query->is_rev = 1;
        }

        len = seed->slc[si].len; 
        if(len <= 2) {
            continue;
        } else if(len > 8) {
            num = idx_buf[si][SEED_LEN-12][2];
            if(num > SW_THRES && seed->slc[si].ext_num > 0) {
                continue; 
            }   
            if(num <= SW_THRES) {
                len = 7; 
            } 
        }
        bg = idx_buf[si][len][0];
        ed = idx_buf[si][len][1];
        num = idx_buf[si][len][2];
        L_off = seed->slc[si].s_off + SEED_LEN -12 -len;
        int seq_off[2], l = 0;
        seq_off[0] = L_off;
        seq_off[1] = seed->slc[si].s_off+SEED_LEN;
        
        if(num > IS_SMLSIZ) {
            uint32_t bwt_idx[4];
            bwt_idx[0] = bg;
            bwt_idx[1] = ed;
            int len_buf[2];
            len_buf[0] = len + 12;
            len_buf[1] = 10;
            l = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);
            if(l > 0) {
                bg = bwt_idx[2];
                ed = bwt_idx[3];
                num = ed+1-bg;
                seq_off[1] += l;
            }
        }
        if(num <= SW_THRES) {//???
            pos_i = 0;
            for(ix=bg; ix <= ed; ix++){
                pos_buf[pos_i] = bwt_sa(fm_idx->bwt, ix);
                
                if(pos_buf[pos_i] > L_off){ 
                    pos_buf[pos_i] = pos_buf[pos_i] - L_off;
                    pos_i++;
                } 
            }
            sub->thres_sw_olp = 1;
            sub->olp_flg = 1;
            sub->err_sum[0] = 0;
            if(query->b0 < query->l_seq/3) {
                sub->thres_pos_num = SW_THRES;
            } else if(query->b0 < query->l_seq/2) {
                sub->thres_pos_num = SW_THRES/5;
            }  else if(query->b0 < query->l_seq*3/5) {
                sub->thres_pos_num = SW_THRES/10;
            } else {
                sub->thres_pos_num = IS_SMLSIZ;
            }
            if(sub->thres_pos_num < IS_SMLSIZ) 
                sub->thres_pos_num = IS_SMLSIZ;
            //pos_i = set_sw_ed_buf(pos_buf, pos_i, query, sub);                              
fprintf(stderr, "%u, pos_i = %d\n", __LINE__, pos_i);
            seed->id = si;
            pos_i = classify_pos(fm_idx, query, seed, 
                        pos_buf, pos_i, seq_off, sub);
fprintf(stderr, "%u, pos_i = %d\n", __LINE__, pos_i);
            if(pos_i > 0) {
                
                //pos_i = clean_hash_pos(fm_idx, pos_buf, pos_i, query, sub);
                AlgnPos(fm_idx, query, pos_i, sub); 
            }
        } 
       
    } //for(si = 0; si < query->seed_num*2; ++si)++++++++
    return ret_flg;
}
int aln_apro_R(idx_t *fm_idx, uint32_t hash_boundry[], struct JmpMod *jmp, struct ExtBlck*eBlck, struct ExtBlck *cB, query_t *query, seed_t *seed, struct StackTree *sTree, struct SubBuf *sub)
{                        
    int i;
    int aln_flg = 0;
    uint32_t ix, bg, ed, num, len, L_off, si, pos_i;
    uint32_t bgn, end, pos, j;
    int l_seed = 0, l_off = query->l_seq;
    for(i = 0; i < query->seed_num; ++i){
        if(seed->slc[i].s_off < l_off) {
            l_off = seed->slc[i].s_off;
            l_seed = i; 
        }
    } 
    if(seed->slc[l_seed].num > SW_THRES) {
        query->is_rev = 0;
        seed->id = l_seed; 
        int b0 = query->b0;
        int p_num = aln_smbwt_uni_R(fm_idx, jmp, eBlck, query, seed, sTree, 1, sub);
        uint32_t k;
        pair_arry_t *pair = sub->pair_out->pair_arry;
        int pi, flg_ed;
        for(pi = 0; pi < p_num; ++pi) {
            bgn = pair[pi].idx_bg;
            num = pair[pi].idx_num;
            if(num == 0) {
                uint8_t nxtflg = pair[pi].nxtflg; 
                if(nxtflg != 1) {
                    printf("%u, read_id=%s, bgn=%u, num=%u, pi=%d\n", 
                            __LINE__, query->name, bgn, num, pi); 
                    exit(1);     
                }                                       
                pos = pair[pi].nxtpnt;
                if(pos < l_off) continue;
                else pos -= l_off;
                flg_ed =  get_sw_ed_val(pos, sub, query);
                if(flg_ed > 0) continue;
                sub->pos_buf[0] = pos ;
                AlgnPos(fm_idx, query, 1, sub);
                set_sw_ed_pos(pos, query, sub);
                continue; 
                
            }
            //if(num > SW_THRES) {
            if(0) {
                uint32_t bwt_idx[4];
                bwt_idx[0] = bgn;
                bwt_idx[1] = bgn+num-1;
                int len_buf[2], seq_off[2];
                len_buf[0] = SEED_LEN+16;
                len_buf[1] = 32;
                seq_off[0] = seed->slc[l_seed].s_off;
                seq_off[1] = seq_off[0] + len_buf[0];
                int l = bsearch_idx_R(fm_idx, query, 
                            len_buf, seq_off, bwt_idx);
                if(l > 16) {
                    printf("read_id=%s, l=%d, bwt_idx[2]=%u, bwt_idx[3]=%u\n", 
                            query->name, l, bwt_idx[2], bwt_idx[3]);
                    exit(1);
                }
                if(l > 0) {
                    bgn = bwt_idx[2];
                    end = bwt_idx[3];
                    num = end+1-bgn;
                    seq_off[1] += l;
                }
            }
            for(j =0; j < num; ++j){
                uint32_t pos = bwt_sa(fm_idx->bwt, bgn+j) 
                                - l_off ;
                flg_ed =  get_sw_ed_val(pos, sub, query);
                if(flg_ed > 0) continue;
                sub->pos_buf[0] = pos ;
                AlgnPos(fm_idx, query, 1, sub);
                set_sw_ed_pos(pos, query, sub);
            }
            
        } 
    }//end if(seed->slc[l_seed].num > SW_THRES) +++++++
    return aln_flg;
}

int aln_apro_L(idx_t *fm_idx, uint32_t hash_boundry[], struct JmpMod *jmp, struct ExtBlck*eBlck, struct ExtBlck *cB, query_t *query, seed_t *seed, struct StackTree *sTree, struct SubBuf *sub)
{                        
    int i;
    int aln_flg = 0;
    uint32_t ix, bg, ed, num, len, L_off, si, pos_i;
    int r_seed = 0, r_off = 0;
    uint32_t bgn, end, pos, j;
    for(i = query->seed_num; i < 2*query->seed_num; ++i){
        if(seed->slc[i].s_off > r_off) {
            r_off = seed->slc[i].s_off;
            r_seed = i; 
        }
    }
    uint32_t k = seed->slc[r_seed].bgn, l = seed->slc[r_seed].end;
    uint32_t n, len1 = 0;
    for(i = r_off-1; i >= r_off-16; --i){ 
        n = bwt_match_exact_alt(fm_idx->bwt, 1, query->rseq+i, &k, &l);
        if(n > 0) {
            len1++;
        } else { break;}
    }
   
    if(seed->slc[r_seed].num > SW_THRES) {
    //if(0)
        query->is_rev = 1;
        seed->id = r_seed; 
        int p_num = aln_smbwt_uni_L(fm_idx, jmp, eBlck, cB, query, 
                                    seed, sTree, 1, sub);
        uint32_t k;
        pair_arry_t *pair = sub->pair_out->pair_arry;
        int pi;
        for(pi = 0; pi < p_num; ++pi) {
            bgn = pair[pi].idx_bg;
            num = pair[pi].idx_num;
            if(num == 0) {
                uint8_t nxtflg = pair[pi].nxtflg; 
                if(nxtflg != 1) {
                    printf("%u, read_id=%s, bgn=%u, num=%u, pi=%d\n", 
                                __LINE__, query->name, bgn, num, pi); 
                    exit(1);     
                }                                       
                uint32_t pos = pair[pi].nxtpnt ;
                if(pos < r_off-16) continue;
                else pos -= r_off-16;
                int flg_ed =  get_sw_ed_val(pos, sub, query);
                if(flg_ed > 0) continue;
                sub->pos_buf[0] = pos ;
                AlgnPos(fm_idx, query, 1, sub);
                set_sw_ed_pos(pos, query, sub);
                continue; 
            }
            int flg_ed = 0;
            int MAX_NUM = IS_SMLSIZ*10;
            if(num > MAX_NUM) {
                r_off = r_off - 16;
                int r_ext = r_off;           
                n = num;
                k = bgn, l = bgn + num-1;
                while( n > MAX_NUM && r_ext >= 0) {  
                    n = bwt_match_exact_alt(fm_idx->bwt, 
                            1, query->rseq+r_ext, &k, &l);
                    if(n > 0) {
                        bgn = k, end = l;
                        num = end+1-bgn;
                        r_off = r_ext; 
                    }
                    r_ext--;
                }
            }
            if(num == 0) {
                if(pair[pi].nxtflg != 1) {
                    printf("%u, error!!!\n", __LINE__);
                    printf("nxtflg=%u, nxtpnt=%u\n", pair[pi].nxtflg, pair[pi].nxtpnt);
                    exit(1); 
                }
                pos = pair[pi].nxtpnt;
                flg_ed =  get_sw_ed_val(pos, sub, query);
                if(flg_ed > 0) continue;
                sub->pos_buf[0] = pos ;
                AlgnPos(fm_idx, query, 1, sub);
                //set_sw_ed_pos(pos, query, sub);
                continue;
            } 
            flg_ed = 0;
            for(j =0; j < num; ++j){
                uint32_t pos = bwt_sa(fm_idx->bwt, bgn+j) 
                                - r_off;
                flg_ed =  get_sw_ed_val(pos, sub, query);
                if(flg_ed > 0) continue;
                sub->pos_buf[0] = pos ;
                AlgnPos(fm_idx, query, 1, sub);
                //set_sw_ed_pos(pos, query, sub);
            }
        }
    } //end if(seed->slc[r_seed].num > SW_THRES) ++++++
    return aln_flg; 
}
int aln_qual_best_one(idx_t *fm_idx, uint32_t hash_boundry[], struct JmpMod *jmp, struct ExtBlck*eBlck, struct ExtBlck *cB, query_t *query, seed_t *seed, struct StackTree *sTree, struct SubBuf *sub)
{
    int i;
    int MAX_NUM = IS_SMLSIZ*10;
    int seed_i = 0, l_off = query->l_seq;
    int s_find = gen_seed_kmer_1(fm_idx,hash_boundry, query, seed);
    int rev;
    uint32_t bgn, end, num, j, k, l, pos;
    for(rev= 0; rev < 2; ++rev) {
        if(s_find == 0) break;
        query->is_rev = rev; 
        int sid = 2*query->seed_num+rev; 
        uint8_t *seq = rev?query->rseq:query->seq;
        int len = seed->slc[sid].len + 12;
        if(len < 16) continue;
        bgn = seed->slc[sid].bgn;
        end = seed->slc[sid].end;
        num = seed->slc[sid].num;
        int s_off = seed->slc[sid].s_off;
        if(num <= MAX_NUM) {
            for(j =0; j < num; ++j){
                uint32_t pos = bwt_sa(fm_idx->bwt, bgn+j) 
                                - s_off ;
                int flg_ed =  get_sw_ed_val(pos, sub, query);
                if(flg_ed > 0) continue;
                sub->pos_buf[0] = pos ;
                AlgnPos(fm_idx, query, 1, sub);
                set_sw_ed_pos(pos, query, sub);
            }
            continue;
        }   
        //---------------------- 
        //以下处理num > MAX_NUM 
        int end_len;
        if(len < SEED_LEN) {
            int flg;
            if(rev == 0) {
                end_len = query->l_seq - (s_off+len);
                if(end_len > 16) flg = 1;
                else flg = 0; 
            } else {
                end_len = s_off;
                if(end_len > 16) flg = 0;
                else flg = 1; 
            }
            if(flg == 1) {  
                uint32_t bwt_idx[4];
                bwt_idx[0] = bgn;
                bwt_idx[1] = end;
                int len_buf[2], seq_off[2];
                len_buf[0] = len;
                len_buf[1] = SEED_LEN-len;
                seq_off[0] = seed->slc[sid].s_off;
                seq_off[1] = seq_off[0] + len_buf[0];
                int l = bsearch_idx_R(fm_idx, query, 
                    len_buf, seq_off, bwt_idx);
                bgn = bwt_idx[2];
                end = bwt_idx[3];
                num = end+1-bgn;
                len += l; 
            } else {
                int l0 = len; 
                for(i = 1; i <= SEED_LEN-len; ++i){ 
                    uint32_t n = bwt_match_exact_alt(fm_idx->bwt, 1, 
                                            seq + s_off - i, &k, &l);
                    if(n > 0) {
                        bgn = k, end = l;
                        num = end+1-bgn;
                        --s_off;
                        l0++;
                    } else {
                        break;
                    }
                }
                len = l0; 
            }
        } //  if(len < SEED_LEN) +++++
        if(num < MAX_NUM) {
            for(j =0; j < num; ++j){
                uint32_t pos = bwt_sa(fm_idx->bwt, bgn+j) 
                                - s_off ;
                int flg_ed =  get_sw_ed_val(pos, sub, query);
                if(flg_ed > 0) continue;
                sub->pos_buf[0] = pos ;
                AlgnPos(fm_idx, query, 1, sub);
                set_sw_ed_pos(pos, query, sub);
            }
            continue;
        }
        if(len < SEED_LEN && num > MAX_NUM ) {
            printf("%u, len = %d, num = %d, query->b0 = %d\n", 
                    __LINE__, len, num, query->b0); 
        }
        if(len == SEED_LEN && num > MAX_NUM ) {
            seed->slc[sid].bgn = bgn;
            seed->slc[sid].end = end;
            seed->slc[sid].num = num;
            seed->slc[sid].len = SEED_LEN;
            end_len = query->l_seq - (s_off+len);
            query->is_rev = rev;
            seed->id = sid; 
            int p_num;
            if(s_off > 16 ) {
                p_num = aln_smbwt_uni_L(fm_idx, jmp, eBlck, cB, 
                                    query, seed, sTree, 1, sub);
            } else if(end_len > 16){
                p_num = aln_smbwt_uni_R(fm_idx, jmp, eBlck, query, 
                                        seed, sTree, 1, sub);
            }
            if(p_num == 0) {
                continue; // for( rev; )    
            }
            pair_arry_t *pair = sub->pair_out->pair_arry;
            int pi;
            for(pi = 0; pi < p_num; ++pi) {
                bgn = pair[pi].idx_bg;
                num = pair[pi].idx_num;
                int r_off = s_off;
                if(num == 0) {
                    uint8_t nxtflg = pair[pi].nxtflg; 
                    if(nxtflg != 1) {
                        printf("%u, read_id = %s\n", 
                                __LINE__, query->name); 
                        exit(1);     
                    }                                       
                    uint32_t pos = pair[pi].nxtpnt ;
                    if(pos < r_off-16) continue;
                    else pos -= r_off-16;
                    int flg_ed =  get_sw_ed_val(pos, sub, query);
                    if(flg_ed > 0) continue;
                    sub->pos_buf[0] = pos ;
                    AlgnPos(fm_idx, query, 1, sub);
                    set_sw_ed_pos(pos, query, sub);
                    continue; 
                }
                int flg_ed = 0;
                if(num > MAX_NUM) {
                    r_off = r_off - 16;
                    int r_ext = r_off;           
                    uint32_t n = num;
                    k = bgn, l = bgn + num-1;
                    while( n > MAX_NUM && r_ext >= 0) {  
                        n = bwt_match_exact_alt(fm_idx->bwt, 
                                1, query->rseq+r_ext, &k, &l);
                        if(n > 0) {
                            bgn = k, end = l;
                            num = end+1-bgn;
                            r_off = r_ext; 
                        }
                        r_ext--;
                    }
                }
                if(num == 0) {
                    if(pair[pi].nxtflg != 1) {
                        printf("%u, error!!!\n", __LINE__);
                        exit(1); 
                    }
                    pos = pair[pi].nxtpnt;
                    flg_ed =  get_sw_ed_val(pos, sub, query);
                    if(flg_ed > 0) continue;
                    sub->pos_buf[0] = pos ;
                    AlgnPos(fm_idx, query, 1, sub);
                    //set_sw_ed_pos(pos, query, sub);
                    continue;
                } 
                flg_ed = 0;
                for(j =0; j < num; ++j){
                    uint32_t pos = bwt_sa(fm_idx->bwt, bgn+j) 
                                    - r_off;
                    flg_ed =  get_sw_ed_val(pos, sub, query);
                    if(flg_ed > 0) continue;
                    sub->pos_buf[0] = pos ;
                    AlgnPos(fm_idx, query, 1, sub);
                    //set_sw_ed_pos(pos, query, sub);
                }
            } //end for(pi = 0; pi < p_num; ++pi)++++
        }//end  if(len == SEED_LEN && num < MAX_NUM ) ++++ 
    } //end for(rev)
    //break;
}// end if(query->b0 < query->l_seq*WOREST_PERCENT*5/3/100)+++
int aln_qual_good_all(idx_t *fm_idx, uint32_t hash_boundry[], struct JmpMod *jmp, struct ExtBlck*eBlck, struct ExtBlck *cB, query_t *query, seed_t *seed, struct StackTree *sTree, struct SubBuf *sub)
{
    int LEN_JUMP = 8;
    int THRES_SUCC = 60;
    int THRES_MIN = 60;
    int suc_flg = 0; 
    set_12mer_qual(query);
    
    int pi, pk, q_val;
    int b0 = query->b0;
    int sn = query->seed_num;
    int sid = 2*sn + 1;
    int sid_bak = seed->id;
    int rev, l_seq = query->l_seq; 
    uint8_t *q_flg = query->q_flg;
    int16_t s_pos[2][LEN_READ/4];
    uint32_t idx_buf[LEN_READ][3];
    int16_t pq_buf[LEN_READ][2];
    int pq_num = 0;
    int idx_num = 0, i, row; 
    uint32_t bgn, end, num, j;
    int len; 
    int s;     
    for(rev = 0; rev < 2; ++rev) {
        for(s = 0; s < sn; ++s) {
            s_pos[rev][s] = seed->slc[s+rev*sn].s_off + SEED_LEN;  
        } 
    }
    //精确扩展中，l_seq-1已处理
    for(rev = 0; rev < 2; ++rev) {
        for(pi = l_seq - 2; pi > 16; --pi) {
            for(s = 0; s < sn; ++s) {
                if(pi == s_pos[rev][s]) {
                    break; 
                } 
            }
            if(s < sn) continue;
            if(rev == 0) {
                q_val = q_flg[pi]%16;
            } else {
                pk = l_seq - pi -1;
                q_val = q_flg[pk]/16; 
            } 
            if(q_val > 0) continue; 
            for(j = pi; j > 16; --j) {
                if(rev == 0) {
                    q_val = q_flg[j]%16;
                } else {
                    pk = l_seq - j -1;
                    q_val = q_flg[pk]/16; 
                }
                if(q_val > 0) break;
            }
            pq_buf[pq_num][0] = pi - j;
            pq_buf[pq_num][1] = pi*16+rev; 
            pq_num++;

            if(pi - j > LEN_JUMP) {
                pi = pi - LEN_JUMP;
            } else {
                pi = j;
            }  
        }
    }
    ks_introsort(pair_t, pq_num, (pair_t *)pq_buf); 
/*  
int rev0, pi0, len0;
if(pq_num > LEN_READ) {
    printf("%u, pq_num = %d\n", __LINE__, pq_num);
    exit(1);
}
for(row = 0; row < pq_num; ++row) {
    len0 = pq_buf[row][0];
    rev0 = pq_buf[row][1]%16;
    pi0 = pq_buf[row][1]/16;
    if(len0 > l_seq-16 || len0 < 0 || 
        rev0 < 0 || rev0 > 1 || 
        pi0> l_seq-2 ||pi0 <16) 
        printf("%u, pq_buf[%d][%d], pq_buf[%d][%d]\n", __LINE__, row, 0, pq_buf[row][0], row, 1, pq_buf[row][1]);
}

   

    for(row = 0; row < pq_num; ++row) {
    fprintf(stderr, "%u, pq_buf[%d][%d], pq_buf[%d][%d]\n", __LINE__, row, 0, pq_buf[row][0], row, 1, pq_buf[row][1]);
}
*/
    idx_num = 0;
    for(row = 0; row < pq_num; ++row) {
        if(suc_flg > 0) break;
        rev = pq_buf[row][1]%16;
        pi = pq_buf[row][1]/16;
        query->is_rev = rev;
        const uint8_t *seq = rev?query->rseq:query->seq;
        int s;     
        seed->id = sid;
        seed->slc[sid].s_off = pi; 
        num = gen_seed_kmer_2(fm_idx,hash_boundry, query, seed);
        if(num == 0) { continue; } 

        bgn = seed->slc[sid].bgn;
        end = seed->slc[sid].end;
        num = seed->slc[sid].num;
        len = seed->slc[sid].len;
        int s_off = seed->slc[sid].s_off;
        if(num > IS_SMLSIZ*10) {
//continue;
            idx_buf[idx_num][0] = num;
            idx_buf[idx_num][1] = pi*1024 + len*16 + rev;
            idx_buf[idx_num][2] = bgn;
            idx_num++;
            continue;
        }         
        for(j =0; j < num; ++j){
            uint32_t pos = bwt_sa(fm_idx->bwt, bgn+j) - s_off ;
            int flg_ed =  get_sw_ed_val(pos, sub, query);
fprintf(stderr, "%u, pos = %u, flg = %d\n", __LINE__, pos, flg_ed);
            if(flg_ed > 0) continue;
            sub->pos_buf[0] = pos;
            AlgnPos(fm_idx, query, 1, sub);
            if(query->b0 >= THRES_SUCC) {
                suc_flg = 1;
                break;
            }
        } 
    }
    seed->id = sid_bak; 
//return suc_flg;
    if(query->b0  > THRES_MIN) return suc_flg; 
    uint32_t bg, ed, pos_i, ix;
    int l;
    uint32_t *pos_buf = sub->pos_buf;
    int seq_off[2];
    uint32_t bwt_idx[4];
    int len_buf[2];
//idx_num = 0;
    for(row = 0; row < idx_num; ++row){
        if(query->b0 > THRES_MIN) break;
        num = idx_buf[row][0];
        bgn = idx_buf[row][2];       
        
        pi = idx_buf[row][1]/1024; 
        len = (idx_buf[row][1]%1024)/16; 
        rev = (idx_buf[row][1]%1024)%16;
        query->is_rev = rev;
        seq_off[1] = pi+1;
        seq_off[0] = seq_off[1] - len;
        //if(num > IS_SMLSIZ) {
        if(0) {
            bwt_idx[0] = bgn;
            bwt_idx[1] = bgn+num-1;
            len_buf[0] = len;
            len_buf[1] = 32;
printf("%u, len_buf[0] = %d, len_buf[1] = %d, bwt_idx[0] = %u, bwt_idx[1] = %u\n", __LINE__, len_buf[0], len_buf[1], bwt_idx[0], bwt_idx[1]);
            l = bsearch_idx_R(fm_idx, query, len_buf, seq_off, bwt_idx);
            if(l > 0) {
                bgn= bwt_idx[2];
                end = bwt_idx[3];
                num = ed+1-bg;
                seq_off[1] += l;
                len += l;
            }
        }
continue;
        if(len >= THRES_MIN || num <= IS_SMLSIZ) {
            if(num > IS_SMLSIZ) {
                num = IS_SMLSIZ; 
            }
            pos_i = 0;
            for(ix=bgn; ix <= end; ix++){
                pos_buf[pos_i] = bwt_sa(fm_idx->bwt, ix);
                if(pos_buf[pos_i] > seq_off[0]){ 
                    pos_buf[pos_i] = pos_buf[pos_i] - seq_off[0];
                    pos_i++;
                } 
            }
            if(pos_i > 0) {
                AlgnPos(fm_idx, query, pos_i, sub); 
            }
        }
    }  

    return suc_flg; 
}

 
#endif
