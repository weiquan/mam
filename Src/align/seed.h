/*
 * =====================================================================================
 *
 *       Filename:  seed.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2018年03月14日 15时05分08秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wei Quan (wq), wquanhit@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef __SEED_H
#define __SEED_H
#define SEED_LEN 20
typedef struct{
    int drct;//read方向,正向序列0， 反向序列1 
    int s_off;//种子序列左端相对于read中心的偏离的碱基个数
    //int h_off;//12kmer hash序列的左端对种子序列左端的偏离碱基个数
    //int m8_off;//8mer序列左端对种子序列左端的相对偏离
    //int aln_p;//指向aln_mer8比对结果的行号
    //int aln_n;//在aln_mer8中比对上的序列数目
    //int flag;//表示操作近似方案,0:不做近似比对，1:bwt比对, 2:二分查找比对
    uint32_t bgn, end, num, pos, len, sum;
    uint32_t max_row, sTree_bg, sTree_num;
    uint8_t ext_num, max_cls; 
} seed_slc_t;

typedef struct{
    int slc_row;//该mer8中在aln_mer8中的行号
    int err_pos;//错误匹配出现的位置；
    int err_val;//错误匹配的碱基
} aln_mer8_t;
#define MAX_SEED LEN_READ/16*3
struct ext_err_t{
    int16_t score;
    int16_t aln_score;
    uint8_t num_N;
    uint8_t sd_err[SEED_LEN+8];//［0］是N出现的个数，［1］开始是种子当中N出现的位置。
    uint8_t seq_L[SEED_LEN+8];//[0]是左扩展中出现的N的个数, [1]开始是N出现的相对位置.
    uint8_t seq_R[SEED_LEN+8];//[0]是右扩展中出现的N的个数, [1]开始是N出现的相对位置.
};

typedef struct{
   
    uint32_t hash_12mer; 
    uint8_t *seed_seq;
    uint8_t seed_buf[20];//SEED_LEN=20 
/*
    uint16_t *seed_8mer_R; 
    uint16_t sort_buf_R[25];
    uint8_t find_buf_R[24];
    int num_find_R;
    
    uint16_t *seed_8mer_L; 
    uint16_t sort_buf_L[25];
    uint8_t find_buf_L[24];
    int num_find_L;
  
    uint8_t var_info[24][3];   
    uint32_t aln_buf[24][6];
    int num_aln; 
*/
    //++++++++++++++++++++++++++++
    

    uint32_t seq12;
    //uint8_t *m_flg;//基因组长度len_ref_seq/8 bits
    int m_num;
    uint32_t *m_pos;//MAX_SEED_NUM 个
    int back_num;
    uint32_t *pos_tmp;
/*      
    uint8_t *l_flg;
    uint8_t *r_flg;
    uint32_t *l_pos;
    uint32_t *r_pos;
*/
    int slc_drct;
    int MAX_SIZE;
    int slc_size;
    seed_slc_t *slc;  
    seed_slc_t *slc_i; 
    //seed_slc_t *slc_b;     
    //aln_mer8_t *aln_8;
    uint32_t *L_len, *R_len;
    int id, cls;

    uint32_t bgn, end, num, max_r, max_f, n_sm_r, n_sm_f, n_rev, n_for;

    //uint32_t idx_arry[2*LEN_READ/SEED_LEN][2];
    //uint32_t idx_buf[(LEN_READ/16 +1)*2][17][3];
    uint32_t idx_buf[MAX_SEED][17][3];

    struct ext_err_t *err;//种子id可以确定当前read的正反向和相对位置
    //int flg[LEN_READ/16*2][LEN_READ/16*2];
    int flg[MAX_SEED][MAX_SEED];
    int del[MAX_SEED];
    
} seed_t;

int get_12mer_correct(uint32_t hash_boundry[], uint32_t end_kmer);

int gen_8mer_sort_varseq(seed_t *seed, int flag);

int get_seed_8mer_bwt_L(idx_t *fm_idx, uint32_t *hash_boundry, seed_t *seed);
#endif
