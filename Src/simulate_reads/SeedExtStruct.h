
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
#define  LEN_SEED  20 
#define  OFF_SEED  10 
#define  LEN_EXT   16 
#define  LEN_READ  150
#define  NUM_EXT  (((LEN_READ- LEN_SEED )/2)+OFF_SEED+ LEN_EXT-1)/LEN_EXT  
#define  MAX_SEED_NUM 600000
//#define  MAX_BLCK_SIZE 600000
#define  IS_SMLSIZ 16

#define  LEN_FILE_NAME 100
#define  MAX_NAME 128
#define  NUM_FILES 6

#define  ONE_INPUT_SIZE 1000000
#define  READ_DATA_FILE ".\\read_data\\read_data_file.txt" 
#define  REFSEQ_SIZE   1000
#define  BWT_SUM_SIZE   1000/4
#define  OCC_INTERVAL 0x80





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




struct SubBuf{  
	int max_buf_size;
    uint32_t (*seqL_out)[3];
	uint32_t (*seqR_out)[3];		
	uint32_t *seqRel_out;		
	uint32_t *seqRel_L;		
	uint32_t *seqRel_R;		
	
    uint32_t seqL_aln_old; 
    uint32_t seqR_aln_old; 
    uint32_t seqRel_aln_old; 
    uint32_t seqRel_L_old; 
    uint32_t seqRel_R_old; 
    //seqL_out[0][0]  是比对成功的个数
	//seqL_out[0][1]  是精确比对的个数
	//seqL_out[0][2]  最大错配的个数
	//seqL_out[1][ ]  1开始是，比对上的SeqL[]的行号和质量,错配位置	

	uint32_t (*pair_out)[3];
	//[0][1]是精确个数，[0][2]是最大错配个数
	//1是开始是nxtpnt和nxtflg，错配分数, 
    //[0][0]是匹配上的个数???
	uint32_t (*algn_out)[5];
//uint32_t (*algnL)[2]; 
uint32_t (*sort_buf)[2]; 
	// [i][0]是第i个比对成功的左序列行号，[i][1]是比对的分数
	uint8_t *algnR_row;
	uint8_t *algnL_row;
	uint8_t *algnRel_row;
	uint8_t *algnRel_buf;

	// [i][0]是1表示第i个右序列比对上，0表示第i个右序列没有比对，
	//uint32_t *algnR_buf;
	// [0]是右序列比对上的个数，[1]开始是比对的行号，该数组用于清理algnR_row[][2]；
	
	uint32_t  off_L ; //LEN_READ/2 - LEN_SEED/2 - LEN_EXT*(this_cls+1);
	uint32_t  off_R ; //LEN_READ/2 + LEN_SEED/2 + LEN_EXT*this_cls;
	uint8_t   ext_seqL[16];
	uint8_t   ext_seqR[16];
    uint8_t   st_posL[10];
    uint8_t   st_posR[10];

	uint8_t seqL[3];
	uint8_t seqR[3]; 
//用于AlnPos函数
    uint8_t target[128];
    uint32_t aln_pos[128]; 
    uint8_t aln_r[128][2];//第i个pos右段比对结果，[i][0]是比对上1否则0, [i][1]比对分数
    
    //kswq_t *swq[2];
    kswq_t *kswq_L[2], *kswq_R[2]; 
    //q = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (m + 4)); // a single block of memory
    int8_t mat[25];
};  


//LEN_EXT_IDX 是所有ExtIdx[][3]数组长度中最大的，可以用文件大小来计算。
struct JmpMod{  
   uint32_t *jmp;
   uint32_t *cap;   
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
struct StackTree {
 
    
    uint32_t (*stck_arry)[4];  
    uint32_t (*back_buf)[7];  
    uint32_t cls;
    int len_arry;
    int len_buf;
    uint32_t tree_node[NUM_EXT];
    uint32_t num_trv;

}StackTree;

struct ExtBlck{
	uint32_t  *head_All ;
	//uint32_t* head_capidx ;
	uint32_t  *head_relat  ;	
	uint8_t   *head_smbwt  ;	
	uint32_t  *head_nxtpnt ;
	uint8_t   *head_nxtflg ;
	uint32_t  (*head_extidx)[2];
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

//LEN_EXT_IDX 是所有ExtIdx[][3]数组长度中最大的，可以用文件大小来计算。

