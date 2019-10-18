
//基于多级索引的Read比对系统框架
//种子序列双向扩展索引的数据结构如下：

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <time.h>
#include <sys/stat.h>
#include <stdint.h>
#include <errno.h>
#include <dirent.h>
#include <stdint.h>
#include <assert.h>
#include <zlib.h>
#include "lookup.h"
#include "ksw.h"
#define LEN_SEED  20 
#define OFF_SEED  10 
#define LEN_EXT   16 
//#define LEN_READ  250
#define  __NUM_EXT  (((LEN_READ- LEN_SEED )/2))/LEN_EXT  
#define NUM_EXT ((__NUM_EXT) < 5 ?(__NUM_EXT):5) 
#define  MAX_SEED_NUM 600000
//#define  MAX_BLCK_SIZE 600000
#define  IS_SMLSIZ 16
//#define LEN_READ 200
//#define SW_THRES (IS_SMLSIZ*20)
#define SW_THRES (IS_SMLSIZ*50)
//#define SW_THRES (IS_SMLSIZ*100)

#define  LEN_FILE_NAME 100
#define  MAX_NAME 128
#define  NUM_FILES 6

#define  ONE_INPUT_SIZE 1000000
#define  READ_DATA_FILE ".\\read_data\\read_data_file.txt" 
#define  REFSEQ_SIZE   1000
#define  BWT_SUM_SIZE   1000/4
#define  OCC_INTERVAL 0x80

#define ALN_R_SIZE 100000

//#define MAX_CLIP 12
#define MAX_CLIP 30 


typedef unsigned char ubyte_t;

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
    int sa_intv, isa_intv;
    bwtint_t n_sa, n_isa;
    bwtint_t *sa;
    bwtint_t *isa; //inverse suffix array
    
    uint32_t n_mod, n_flg, n_pos;     
    bwtint_t *sa_pos;
    bwtint_t *sa_mod;
    uint64_t *sa_flg;    


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

    lkt_t *fastmap;/*fastmap 12bp seq*/
    
    bwt_t *bwt;/*FM index*/
    bntseq_t *bns;
    uint8_t *pac;/* packed reference seq */
    /* snp-aware seed index */
    int n_pmap;
    uint32_t *pmap;/* preseq map */ 
    int n_var;
    uint8_t *refseq;/*2bp refseq at variant site*/
    uint32_t *sai;/* sa index */
    
    /* approximate match */
    bp_t *is_multiseeds;/*is_multiseeds[i] = 1 if Reference[sa[i]:sa[i]+20] is multi location seeds */
    uint32_t n_isa2seq16;
    uint32_t* isa2seq16; /* isa2seq16[i] = start of seq16s, i is occ of 1 in is_multiseeds before isa*/
    uint32_t n_seq16s;
    uint32_t *seq16s;

    int n_tot;
    uint32_t *lext_idx;
    uint32_t *lext1_idx;
    uint32_t *rext_idx; 
    
    int n_rext, n_lext0, n_lext1;
    ext_t *rext0, *lext0;
    uint32_t *lext1;
    /* 
    bp_t *is_var;// variants map
    uint32_t n_var;
    uint32_t *var_seqs;// var_seqs[2i] = variant pre seq; var_seqs[2i+1] = variant index 
    all_var_t *var;*/
} idx_t;
// 循环控制信息--------------------------------------------
//uint32_t queue_align_info[][2];  //开辟较大内存
//uint32_t this_stack_tree[NUM_EXT][5];  
//[i][0]是第i个节点(即，第i级索引）在queue_align_info中的起始行号
//[i][1]是第i个节点(即，第i级索引）在queue_align_info中的个数
//uint32_t this_stack_node[3];



  
typedef struct{//比对结果结构体
    int num, len_buf, best_num;
    int16_t score[ALN_R_SIZE+2][4];
    uint32_t pos[ALN_R_SIZE+2];
    uint32_t *pos_buf;
    int delta;
    int size;
} aln_rlt_t;
typedef struct{//比对结果结构体
    int best[10];//最好的比对信息
    int len;//比对上数据保存的个数
    int num;//比对上的分数名称的个数
    int max_found;
    int found[10][4];//已经比对的当中保存的信息
    //[i][0] 是第i个的分数名称, 如195;
    //[i][1] 是第i个的保存数目;
    //[i][2] 是第i个的第一个数据区的指针；
    //[i][3] 是第i个的尾数据区的指针。
    int max_out;
    uint32_t (*out_buf)[2]; 
    int L_sc;
    int L_ti, L_qi;
    int R_sc;
    int R_ti, R_qi;
} aln_out_t;

typedef struct {
    uint32_t nxtpnt;
    uint32_t idx_bg;
    uint32_t idx_num;
    uint8_t nxtflg;
    uint8_t err;
    int8_t l_off;
    int8_t r_off;
} pair_arry_t;
typedef struct {
    pair_arry_t *pair_arry;
    int p_num;
    int s_num;
    int p_buf[3];
} pair_out_t;
typedef struct{
    uint8_t target[LEN_READ+40];
    uint32_t pos;
    uint32_t err_pos[LEN_READ/4+10];
    int seq_off[2];
    int err_n;
} tgt_t;
typedef struct{
    tgt_t tgt[SW_THRES+10];
    int tgt_n;
    int cur_err_n;
    uint32_t cur_pos;
    uint8_t *cur_tgt; 
    uint32_t *cur_err_pos;
    int min_err_n;
    uint32_t min_pos;
    uint8_t *min_tgt; 
    uint32_t *min_err_pos;
} tgt_arry_t; 

#define SEED_LEN 20
struct SubBuf{  

    int max_buf_size;
    //uint32_t (*seqL_out)[4];//[0]: 行号； [1]:错误分数；[2]:错误位置；[3]:偏移量
	//uint32_t (*seqR_out)[4];	
    int32_t (*call_buf)[4];
    int32_t (*seqL_out)[4];//[0]: 行号； [1]:错误分数；[2]:错误位置；[3]:偏移量
	int32_t (*seqR_out)[4];	
    int32_t seqL_off[NUM_EXT+2]; 
    int32_t seqR_off[NUM_EXT+2];

    uint32_t *seqRel_out;		
	uint32_t *seqRel_L;		
	uint32_t *seqRel_R;		
	
    uint32_t seqL_aln_old; 
    uint32_t seqR_aln_old; 
    uint32_t seqRel_aln_old; 
    uint32_t seqRel_L_old; 
    uint32_t seqRel_R_old; 
    //uint32_t (*pair_out)[5];
    pair_out_t *pair_out;
    uint32_t *pos_buf;
    uint32_t (*err_buf)[5];
    uint32_t out_buf[12*LEN_READ];// 
    //1. [0][0]是匹配上的个数, [0][1]是精确个数，[0][2]是最大错配个数
	//2. 左坐标1开始数据，[1][0]是nxtpnt,[1][1]是nxtflg，[1][2]是错配分数,
    //   [1][3]是左偏移量,[1][4]是右偏移量 
    int pair_b2e[3][4];
    int pair_b2e_num, pair_b2e_id;
    uint32_t (*algn_out)[5];
    uint32_t (*sort_buf)[2]; 
	// [i][0]是第i个比对成功的左序列行号，[i][1]是比对的分数
	uint8_t *algnR_row;
	uint8_t *algnL_row;
	uint8_t *algnRel_row;
	uint8_t *algnRel_buf;
    uint32_t *relat_buf;
	// [i][0]是1表示第i个右序列比对上，0表示第i个右序列没有比对，
	//uint32_t *algnR_buf;
	// [0]是右序列比对上的个数，[1]开始是比对的行号，该数组用于清理algnR_row[][2]；
	
	uint32_t  off_L ; //LEN_READ/2 - LEN_SEED/2 - LEN_EXT*(this_cls+1);
	uint32_t  off_R ; //LEN_READ/2 + LEN_SEED/2 + LEN_EXT*this_cls;
	uint8_t   ext_seqL[18];
	uint8_t   ext_seqR[18];
    int8_t   aln_L[30];
    //0-7起始pos，10是len_seq, 11是起始seq位置, 12是err_len，13是InDel,14是 err_mov, 
    int8_t   aln_R[30];//
    uint8_t seqL[3];
	uint8_t seqR[3]; 
    kswq_t *qry_f[2], *qry_r[2];
    kswq_t *kswq_f[2], *kswq_r[2]; 
    int8_t mat[25];
    kswgb_t kswgb;
    kswst_t kswst;
    aln_rlt_t *aln_r;
    aln_out_t *aln_out;
    uint32_t *buf;//用于SW函数的pos内存空间缓冲区

    int err_sum[5];
    int trm_l, trm_r;
    int sub_err, ins_err, del_err, delta;
    //int64_t *sw_ed[2];//［0］是正向read，［1］是反向read
    //int64_t *sw_to[2];
    int max_pos;
    //uint32_t *pos_ed[2];
    //uint32_t *pos_to[2];
    //uint32_t *pos_bk[2];


    //++++++++++++++++++++++++++++++
    int to_dat_size;
    int dat_blk_len;
    int rlt_blk_len;
    int hsh_shift;
    int rlt_shift; 
    int bit_shift; 

    int sw_ed_bit_size;
    int sw_to_bit_size; 
    int sw_ed_hsh_size;
    int sw_to_hsh_size;
    int sw_ed_rlt_size;
    int sw_to_rlt_size;
    int pos_ed_buf_size;
    int pos_to_buf_size;
    int pos_bk_buf_size;

    int sw_ed_rlt_max[2], sw_to_rlt_max[2], sw_ed_bit_max[2], sw_to_bit_max[2];
    uint64_t *sw_ed_bit[2];//［0］是正向read，［1］是反向read
    uint64_t *sw_to_bit[2];
    
    uint32_t *sw_ed_hsh[2];//［0］是正向read，［1］是反向read
    uint32_t *sw_to_hsh[2];
 
    uint32_t *sw_ed_rlt[2];//［0］是正向read，［1］是反向read
    uint32_t *sw_to_rlt[2];
 
    uint32_t *pos_ed_buf[2];
    uint32_t *pos_to_buf[2];
    uint32_t (*pos_bk_buf[2])[2];
    int seed_bk[LEN_READ/SEED_LEN*4]; 
    int seed_id, n_seed_bk[2];
//------------------------------
    int query_b0, query_len, query_err;
    int path_err[NUM_EXT+1];
    int cls;

    int indel_f[10];
    int sub_f[10];
    int NEXT_EXT_NUM;
    uint32_t NT_sum[5];
    int calcu_flg[4];
    
    int thres_sw_to;
    int thres_pos_num;
    int thres_sw_olp;
    int olp_flg;
    int n_sm_f, n_sm_r;  
    tgt_arry_t *tgt_arry;
    uint32_t idx_for[LEN_READ+1][2];
    uint32_t idx_rev[LEN_READ+1][2];
    int l_off, r_off, flg_off;
    int l_seq_delta;

   
    int eval_sc[2*2*(LEN_READ/12)];
    int eval_pos[3*LEN_READ];
    //int ref_bg, ref_ed, read_bg, read_ed; 
    int SLC_SWITCH[10];
    int ALN_SWITCH[10];
	int THRES_BAK;
    
 
};  
//LEN_EXT_IDX 是所有ExtIdx[][3]数组长度中最大的，可以用文件大小来计算。
struct JmpMod{  
   uint32_t *jmp;
   //uint32_t *cap;   
   uint8_t *mod;
};
/*
struct FileName{
	int  numidxfiles;//init 6

	char (*capidx)[LEN_FILE_NAME]  ;
	char (*relat)[LEN_FILE_NAME]   ;
	char (*smbwt)[LEN_FILE_NAME]   ;
	char (*nxtpnt)[LEN_FILE_NAME]  ;
	char (*nxtflg)[LEN_FILE_NAME]  ;
	char (*smpos)[LEN_FILE_NAME]   ;

	char (*comfile)[LEN_FILE_NAME] ;
	char (*seedidx)[LEN_FILE_NAME];

	//char (*jmpmod)[LEN_FILE_NAME];
    char *jmpmod;
};
*/

typedef struct {
    uint32_t nxtpnt;
    uint32_t idx_bg;
    uint32_t idx_num;
    uint32_t len_p, len_g;
    uint8_t nxtflg;
    uint8_t err; 
    int8_t l_off;
    int8_t r_off;
    int8_t cls;
    int8_t m_cls;
    int8_t aln_f;
    int8_t aln_r;
} stck_arry_t;
struct StackTree {
    //uint32_t (*stck_arry)[6];  
    stck_arry_t *stck_arry; 
    //uint32_t (*back_buf)[9];  
    //uint32_t (*exted_arry)[6];  
    stck_arry_t *exted_arry;
    stck_arry_t *smpos_arry;
    stck_arry_t *back_arry;

    int len_cls;
    uint32_t cls;
    int max_cls, max_ext;
    int len_arry;
    int len_exted;
    int len_smpos;
    int add_num;
    int len_old, len_old_back;
    int len_buf;
    int bg_buf;
    int len_back;
    //uint32_t tree_node[NUM_EXT];
    //uint32_t num_trv;
    //uint32_t cap_row, cap_num;
    //int err;
    int seed_off, seed_id;    
    struct call_t *call;

}StackTree;

struct ExtBlck{
	uint32_t  *head_All ;
	//uint32_t* head_capidx ;
	uint32_t  *head_relat  ;	
	uint8_t   *head_smbwt  ;	
	uint32_t  *head_nxtpnt ;
	uint8_t   *head_nxtflg ;
	uint32_t  (*head_extidx)[2];
	
    int  n_relat  ;	
	int  n_smbwt  ;	
	int  n_nxtpnt ;
	int  n_nxtflg ;
	int  n_extidx;
    int  n_cap;

    
    struct    CapIfo *head_cap;

	//+++++++++++++++++++++++++++++++++++++++++++++++
	//下面是指针偏移量
	uint32_t relat  ; //head_relat  + head_capidx[bgn_row].relat;	
	uint32_t smbwt  ; //head_smbwt  + head_capidx[bgn_row].smbwt;	
	uint32_t nxtpnt ; //head_nxtpnt + head_capidx[bgn_row].nxtpnt;
	//uint32_t nxtflg ; //head_nxtflg + head_capidx[bgn_row].nxtflg;
	uint32_t extidx;
    uint32_t bwtL   ; //smbwt ;
	uint32_t len    ;  // getBwtSumSize(uint32_t num);
	uint32_t bwtR   ;  //smbwt + (bwt_L+1)/2 + len;

	
	uint32_t num_seqL ;
	uint32_t num_seqR ;
	uint32_t num_relat;

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//seqL[0] = this_num_seq_L; // DataNum
	//seqL[1] = LEN_EXT-1;      // begining pos of seq[];
	//seqL[2] = 0;      		 // begining rot of  Bwt

	//seqR[0] = this_num_seq_R; // DataNum
	//seqR[1] = LEN_EXT-1;      // begining pos of seq[];
	//seqR[2] = 0;  
    uint32_t bg_idx, ed_idx;
};

struct CapIfo{
	uint32_t relat  ;	
	uint32_t smbwt  ;	
	uint32_t nxtpnt ;
	//uint32_t nxtflg ;
	uint32_t num_seqL ;
	uint32_t num_seqR ;
	uint32_t num_relat ;	

};
struct call_t{
    int8_t sub_O_L;
    int8_t sub_O_R;
    int8_t sub_A_L; 
    int8_t sub_A_R; 
    int8_t ind_O_L[4];
    int8_t ind_O_R[4];
    int8_t ind_A_L[4]; 
    int8_t ind_A_R[4];

    int num_L[8];
    int num_R[8];
    int aln_f;
    int aln_r;
    int pair_f;
    int rep_i;
    uint32_t (*buf)[4];
    uint32_t off_L;
    uint32_t off_R;
    uint32_t cap_row, cap_num;
    int err, len_old;
};
//LEN_EXT_IDX 是所有ExtIdx[][3]数组长度中最大的，可以用文件大小来计算。

