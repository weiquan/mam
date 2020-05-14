
//种子扩展多级索引层次数据结构
//种子序列双向扩展索引的数据结构如下：
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
//#include <dirent.h>
//#include <sys/time.h>
//#include <time.h>
////#include <sys/stat.h>
#include <stdint.h>
#include "ksort.h"
#include "utils.h"
#include "setFileName.h"
#include "index.h"

#include "malloc_wrap.h"
#define USE_MALLOC_WRAPPERS


typedef struct{
    uint32_t seq;
    uint32_t idx;
} pair_t;
#define ks_lt_pair(a, b) ((a).seq==(b).seq?(a).idx<(b).idx:(a).seq<(b).seq)
KSORT_INIT(pair_t, pair_t, ks_lt_pair)
#define SWAP(type_t, a, b) do{type_t x=(a); (a)=(b); (b)=x;} while(0)

//KSORT_INIT_GENERIC(uint64_t)

#define  LEN_SEED  20 
#define  OFF_SEED  10 
#define  LEN_EXT   16 
#define  LEN_READ  150
#define  NUM_EXT  (((LEN_READ- LEN_SEED )/2)+OFF_SEED+ LEN_EXT-1)/LEN_EXT  
#define  MAX_SEED_NUM 600000
#define  LEN_FILE_NAME 100 
#define  IS_SMLSIZ 16
#define  NUM_SCALE 50
#define  MIN_BWT_SIZE 16
#define  NO_BWT_SUM_SIZE 64

    //LEN_EXT_IDX 是所有ExtIdx[][3]数组长度中最大的，可以用文件大小来计算。

//--------------------------------------------------------------
//扩展种子Index区间生成过程：
//主BWT索引已经创建，包括:
//Bwt后缀前置序列数组SA[]，
//累加数组Rank[]，
//Index到Pos的映射数组ItoPos[]数组，
//参考基因组序列RefSeq[]数组，
//-----------------------------------------------------------------

struct CapIfo {
	uint32_t relat;	
	uint32_t smbwt;
  uint32_t nxtcap;	
	uint32_t num_seqL ;  // 当前块中seqL数据个数；
	uint32_t num_seqR ;  // 当前块中seqR数据个数；
	uint32_t num_relat;  // 当前块中relat数据个数；
     
};

struct SubBuf{  

  int max_buf_size;	
  int num_ext;
  uint32_t (*seqL_buf)[2];  //截取的数据缓冲区
	uint32_t (*seqR_buf)[2];
//+++++++++++++++++++++++++++++++++++++++++
	uint32_t *sortL; //排序后的唯一序列
	uint32_t *sortR;
	uint32_t *relat;
    //uint32_t *R_relat;	
	uint32_t *L2rel ;
  uint32_t *R2rel ;
  uint32_t (*idxR)[2]; // sortR[]的每一个序列的(bgnIdx,ednIdx)
  uint8_t  *buf_8;
  uint16_t *buf_16;
  uint32_t *buf_32;

	//注：idxR[][2]的数据存放在;seqL_buf[][2]的物理空间。
//++++++++++++++++++++++++++++++++++++++++++++++++++++
	uint32_t (*extIdx)[2]; // sortR[]的每一个序列的(bgnIdx,ednIdx)
	//注：extIdx[][2]的数据存放在;seqR_buf[][2]的物理空间。	
	uint32_t (*seqL_idxR)[3];  // relat[]的每个行上左序列和右序列对应的idx,即(seqL,bgnIdxR,endIdxR)	

	uint32_t *pos_buf;  // 保存当前快中的pos数据； 	
	uint32_t *nxt_idx;  // 保存当前快中的nxt_idx数据；
	uint8_t  *nxt_flg;  // 保存当前快中的nxt_flg数据；


	uint8_t  *bwt_seqL;   
	uint8_t  *bwt_sumL; 
	uint8_t  *bwt_seqR;   
	uint8_t  *bwt_sumR;
	uint8_t  *smbwt  ;

	uint32_t num_seqL ;  // 当前块中seqL数据个数；
	uint32_t num_seqR ;  // 当前块中seqR数据个数；
	uint32_t num_relat;  // 当前块中relat数据个数；	
	uint32_t num_pos ; 	 // 当前块中pos数据个数；

	uint8_t  num_bwtseqL;   
	uint8_t  num_bwtsumL; 
	uint8_t  num_bwtseqR;   
	uint8_t  num_bwtsumR;

	uint32_t num_bwtL;
	uint32_t num_bwtR;
	uint32_t len_capidx ;  // 当前块前为止capidx数组中的全局capidx数据长度；
	uint32_t len_relat  ;  // 当前块前为止relat数组中的全局relat数据长度；	
	uint32_t len_smbwt  ;  // 当前块前为止smbwt数组中的全局smbwt数据长度；	
	uint32_t len_nxtpnt ;  // 当前块前为止nxtpnt数组中的全局nxtpnt数据长度；
	uint32_t len_nxtflg ;  // 当前块前为止nxtflg数组中的全局flg数据长度；
	uint32_t len_extidx  ;  // 当前块前为止extidxs数组中的全局pos数据长度；

	struct CapIfo cap[1];
	//uint32_t off_pos ; 
    //uint32_t blck_id;

};  


struct SeedIdx {
	uint32_t (*cur_seedidx)[2];
	uint32_t (*ord_seedidx)[2];
	uint32_t (*nxt_seedidx)[2];
	uint32_t *jmp_idx;
	uint32_t *cap_pos;
	uint8_t  *mod_idx; 
};
void InitSeedIdx(struct SeedIdx *sidx,uint32_t MaxLen, uint32_t ref_len);

//void setFileName(char *in_fname, struct FileName *out_fname);
uint32_t getlen(uint64_t LenExtIdx[],struct FileName *inf);
void InitSubBuf(struct SubBuf *sub_buf);
//int CrtJmpMod(int n,uint32_t cur_seedidx[][2], uint32_t* jmp_idx, uint8_t *mod_idx, uint32_t *cap_pos);

int CrtJmpMod(int n, uint32_t cur_seedidx[][2], uint32_t* jump_idx, uint8_t *mod_idx, uint32_t *out_buf, int n_ext, idx_t *fm_idx);
void OpenSeedIdxfiles(struct  FileName *f, uint32_t (*sidx)[2], uint64_t LenExtIdx[],uint32_t repNum,FILE *fpw[]);

void relatNxtCapIdx(struct SubBuf *subBuf,uint32_t *jmp_idx, uint8_t *jmp_mod,uint32_t *cap_pos,idx_t *fm_idx);
void comFile(struct FileName *rf, struct FileName *wf);
uint64_t getFileSize(const char *file);
void getBlckData(uint32_t sidx[2],uint32_t num_ext, struct SubBuf *sub_buf, idx_t *idx);

void gen_cnt2(const uint8_t *bwt, int n_data, uint8_t *cnt);
void gen_cnt_all(const uint8_t* bwt, int n_data, uint8_t *cnt);
void buldSmBwt(struct SubBuf *sub,int flg);
void AlgnSeq(PubParm  *pParm, char seq[], uint32_t buf_idx[]);
void test_AlgnBwtSml(uint8_t*bwt, uint32_t DataNum, uint8_t *seq_buf);

uint32_t get_bwt_size(uint32_t n_data);



void log_sub(struct SubBuf *sub)
{

    int i, j;
   


    printf("============================\n");
    for(i=0; i< sub->num_seqL; ++i){
        uint32_t seq = sub->sortL[i];
        for(j=0;  j<16; ++j){
            printf("%u", (seq>>(30-j*2))&3); 
        }
        printf("\t%u\n", sub->L2rel[i]);
    }
    printf("\n---------------\n");
    for(i=0; i< sub->num_relat; ++i){
        printf("%u\n", sub->relat[i]); 

    }
    
    printf("\n---------------\n");
    for(i=0; i< sub->num_seqR; ++i){
        uint32_t seq = sub->sortR[i];
        for(j=0;  j<16; ++j){
            printf("%u", (seq>>(30-j*2))&3); 
        }
    
        printf("\n");
    }
   
    printf("\n============================\n");




}
void log_seqs(idx_t *idx, uint32_t k, uint32_t l, int n_ext)
{
    uint32_t i, j;
    for(i = k; i <l; ++i){
        uint32_t pos = bwt_sa(idx->bwt, i); 
        printf("%u\t", i);      
        if(pos < 16 || pos +20+32*n_ext >= idx->bwt->seq_len) {
            

            if(pos < 16){
                for(j = 0; j < 16; ++j){
                    printf("*"); 

                }   
            } else{
                
                for(j = pos - 16; j < pos; ++j){
                    printf("%u", __get_pac(idx->pac, j)); 

                }
   
            }
            printf("\t");
            
            for(j = pos; j < pos+20+32*n_ext; ++j){
                printf("%u", __get_pac(idx->pac, j)); 
            }

            printf("\t");
            if(pos+20+32*n_ext >= idx->bwt->seq_len){
                for(j =0; j < 16; ++j){
                    printf("*"); 
                }
            } else{
            uint32_t ed = pos+20+32*n_ext;
            for(j = ed; j < ed+16; ++j){
                printf("%u", __get_pac(idx->pac, j)); 
            }
   
            
            }
            printf("\n");
   
            continue; 
        } 
        uint32_t bg = pos -16;
        uint32_t ed = pos +20+32*n_ext; 
        for(j = bg; j < bg+16; ++j){
            printf("%u", __get_pac(idx->pac, j)); 

        }
        printf("\t");
        for(j = bg+16; j < ed; ++j){
            printf("%u", __get_pac(idx->pac, j)); 
        }
        printf("\t");
        for(j = ed; j < ed+16; ++j){
            printf("%u", __get_pac(idx->pac, j)); 
        }
        printf("\n");
    
    }

}



#define FILE_PREFIX "../Index/GRCh38"
//int main(int argc, char *argv[]){
//build hierarchical index
//h1 index for 20bp seeds
//h2 index for 52bp seeds
//h3 index for 84bp reads
//...
int build_hier_idx(const char* prefix){
  //+++++++++++++++++++++++++++++++++++++++++++
  //参考基因组，Bwt，I2Pos数组初始化
  //idx_t *fm_idx = (idx_t *)calloc(1, sizeof(idx_t)); 
	fprintf(stderr, "restore fm_index %s\n", prefix);
  idx_t* fm_idx = fmidx_restore(prefix);

  //+++++++++++++++++++++++++++++++++++++++++++
	//初始化
  FileName  *f = (FileName*)malloc(sizeof(FileName));
	setFileName(f, prefix);
  uint64_t LenExtIdx[NUM_EXT+1] ;
	uint32_t MaxLen; 
	uint32_t i = getlen(LenExtIdx, f);  
  //LenExtIdx是所有ExtIdx[][3]数组长度中最大的，
	//可以用文件大小来计算,返回最大值序号。
	MaxLen = LenExtIdx[i];
    //内存初始化------------------------------------------------------
	//uint32_t (*cur_seedidx)[3], (*ord_seedidx)[3], (*nxt_seedidx)[3];
	uint32_t (*cur_seedidx)[2], (*ord_seedidx)[2], (*nxt_seedidx)[2];
	uint32_t *jmp_idx, *cap_pos;
	uint8_t  *mod_idx; 
	struct SeedIdx *sIdx = (struct SeedIdx*)malloc(sizeof(struct SeedIdx));

	InitSeedIdx(sIdx, MaxLen, fm_idx->bwt->seq_len);
	cur_seedidx = sIdx->cur_seedidx;
	nxt_seedidx = sIdx->nxt_seedidx;
	jmp_idx     = sIdx->jmp_idx;
	mod_idx     = sIdx->mod_idx;
	cap_pos     = sIdx->cap_pos;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//数组初始化----------------------------------------------------
	FILE *fp = xopen(f->seedidx[0],"rb");
  fread(cur_seedidx, 2*sizeof(uint32_t), LenExtIdx[0], fp);
  err_fclose(fp);
	
    
    //================================================================
	//计算cur_seedidx的jmp_idx,mod_idx, cap_pos数组，并输出到FlgIdxFile_1
	uint32_t out_buf[2];    
  int max_buf_size = CrtJmpMod(LenExtIdx[0],cur_seedidx, jmp_idx, mod_idx,out_buf, 0, fm_idx);
  fprintf(stderr, "%u, out_buf[0] = %d, max_buf_size = %u\n", __LINE__, out_buf[0], max_buf_size);

//++++++++++++++++++++++++++++++++++++++++++++++++++
	struct SubBuf *sub;
    if(NULL == (sub = (struct SubBuf*)malloc(sizeof(struct SubBuf)))){
	    perror("error:[malloc(NUM_EXT*sizeof(SubBuf))]");
	    exit(1);
	}
  sub->max_buf_size = max_buf_size;
  InitSubBuf(sub);
  struct CapIfo  *cap;
	cap = sub->cap;

//++++++++++++++++++++++++++++++++++++++++++++++++++
    // 权威
  if((fp=fopen(f->jmpmod,"wb"))==NULL){//判断是否打开文件
    printf("can't open file = %s\n",f->jmpmod);
    //exit(0);
	}

  uint32_t l_jmpmod;

  uint32_t n_jmp = (fm_idx->bwt->seq_len+512-1)/256;
  fwrite(&n_jmp, sizeof(uint32_t), 1, fp);
  uint32_t n_mod = (out_buf[0]+3)/4*4;
  fwrite(&n_mod, sizeof(uint32_t), 1, fp);
  fwrite(&sub->max_buf_size, sizeof(uint32_t), 1, fp);
	fwrite(jmp_idx,sizeof(uint32_t), n_jmp, fp);  
	fwrite(mod_idx,sizeof(uint8_t) , n_mod, fp);  

  fclose(fp); 




	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//开始循环处理----------------------------------------------------
	uint8_t repNum = 0;
	FILE *fpw[6];
	FILE *fp_capidx;
  FILE *fp_relat;
	FILE *fp_smbwt;
	FILE *fp_nxtpnt;
	FILE *fp_nxtflg;
	FILE *fp_extidx;

	int NumExt = NUM_EXT;
  repNum = 0;
  uint32_t MAX_RELAT = 0;
  while(repNum < NumExt){
    sub->len_capidx = 0;	
    cap[0].relat = 0;
    cap[0].smbwt = 0;
    cap[0].nxtcap = 0; 
     
    OpenSeedIdxfiles(f, nxt_seedidx,LenExtIdx,repNum,fpw);
    
    fp_relat  = fpw[0];
		fp_smbwt  = fpw[1];
		fp_nxtpnt = fpw[2];
		fp_nxtflg = fpw[3];
		fp_capidx = fpw[4];
//fp_extidx = fpw[5];
    fprintf(stderr, "OpenSeedIdxFiles finish!\n");
		uint32_t bgn,end,num,j,idx, len_seedidx;
    len_seedidx = out_buf[0];       
    fprintf(stderr, "%u, out_buf[0] = %d, max_buf_size = %u\n", __LINE__, out_buf[0], max_buf_size);
    CrtJmpMod(LenExtIdx[repNum+1],nxt_seedidx, jmp_idx, mod_idx, out_buf, repNum+1, fm_idx);
    fprintf(stderr, "%u, rep_Num = %d, out_buf[0] = %d, max_buf_size = %u\n", __LINE__, repNum, out_buf[0], max_buf_size);
    idx = 0;
    while( idx < len_seedidx){
      bgn = cur_seedidx[idx][0];
      num = cur_seedidx[idx][1];
      end = bgn + num-1; 
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //以下代码段 repNum = NUM_EXT , 即最后一个三元数据时运行
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //以下代码段 repNum < NUM_EXT 时运行

      if(num <= IS_SMLSIZ){
        idx++;
        continue;
        //小规模数据不做任何处理
        //比对过程中可以通过自身的jmp_idx,mod_idx信息，
        //可以找到相应的SmPos[]数组的行号
      } 
      //fprintf(stderr, "382, if(num <= IS_SMLSIZ :  end\n");
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //以下是当num_idx>IS_SMLSIZ时进行------------------------- 
      //生成左扩展序列集合SeqL[]和右扩展序列集合SeqR[]，
      //计算关联关系矩阵数组Relat[ ]和映射数组LtoRel[ ], 

      uint32_t sidx[2];
      sidx[0] = bgn;
      sidx[1] = num;
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //预编译期间临时屏蔽
      getBlckData(sidx,repNum,sub,fm_idx);  

      //log_sub(sub);
      //log_seqs(fm_idx, sidx[0], sidx[0]+sidx[1], repNum);
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Relat[  ]和LtoRel[ ]分段输出到同一个文件RelatFile_i。
     
      int bg, ed, relat_len;           
      if(sub->num_relat >MAX_RELAT) MAX_RELAT = sub->num_relat;
      if(sub->num_relat <= 0xFF)	{
        uint8_t *buf_8 = sub->buf_8;
        bg = 0; 
        relat_len = sub->num_relat;
        for(i = bg; i < relat_len; ++i){
            buf_8[i] = (uint8_t)sub->relat[i]; 
        }
        bg = sub->num_relat;
        relat_len += sub->num_seqL+1;
        for(i = bg, j =0; i < relat_len; ++i, ++j){
            buf_8[i] = (uint8_t)sub->L2rel[j]; 
        }
        bg = relat_len;
        relat_len += sub->num_seqR+1;
        for(i = bg, j=0; i < relat_len; ++i, ++j){
            buf_8[i] = (uint8_t)sub->R2rel[j]; 
        }
        /*  
        bg = relat_len; 
        relat_len += sub->num_seqR+1;
        for(i = bg, j=0; i < relat_len; ++i, ++j){
            buf_8[i] = (uint8_t)(sub->idxR[j][0] - sub->idxR[0][0]); 
        }
        */
        relat_len = (relat_len+3)/4*4;
        fwrite(buf_8,sizeof(uint8_t),relat_len,fp_relat);
        relat_len = relat_len/4;
      } else if(sub->num_relat <=0xFFFF ){
        uint16_t *buf_16 = sub->buf_16; 
        bg = 0;
        relat_len = sub->num_relat;
        for(i = bg; i < relat_len; ++i){
            buf_16[i] = (uint16_t)sub->relat[i]; 
        }
        
        bg = relat_len; 
        relat_len += sub->num_seqL+1;
        for(i = bg; i < relat_len; ++i){
            j = i-sub->num_relat;
            buf_16[i] = (uint16_t)sub->L2rel[j]; 
        }
        
        bg = relat_len; 
        relat_len += sub->num_seqR+1;
        for(i = bg, j=0; i < relat_len; ++i, ++j){
            buf_16[i] = (uint16_t)sub->R2rel[j]; 
        }
        /*   
        bg = relat_len; 
        relat_len += sub->num_seqR+1;
        for(i = bg, j=0; i < relat_len; ++i, ++j){
            buf_16[i] = (uint16_t)(sub->idxR[j][0] - sub->idxR[0][0]); 
        }
        */
        relat_len = (relat_len+1)/2*2;
        fwrite(buf_16,sizeof(uint16_t),relat_len,fp_relat);
        relat_len = relat_len/2;
      } else{
        uint32_t *buf_32 = sub->buf_32; 
        fwrite(sub->relat,sizeof(uint32_t),sub->num_relat,fp_relat);
        // buf_relat是当前Relat[]数据的缓冲区
        // len_relat是当前buf_relat数据的个数
        fwrite(sub->L2rel,sizeof(uint32_t),sub->num_seqL+1,fp_relat);
        fwrite(sub->R2rel,sizeof(uint32_t),sub->num_seqR+1,fp_relat);
        /*  
        for(i = 0; i < sub->num_seqR+1; ++i){
            buf_32[i] = (uint32_t)(sub->idxR[i][0] - sub->idxR[0][0]); 
        }
        fwrite(buf_32,sizeof(uint32_t),sub->num_seqR+1,fp_relat);
        */
        // buf_L_Rel是当前LtoRel[ ]数据的缓冲区
        // num_seq_L是当前SeqL数据的个数
        //num_seq_L, num_seq_R，len_relat, 分别保存在cap_idx的相应分量中；
        relat_len = sub->num_relat+sub->num_seqL+sub->num_seqR+2;
          //relat_len += sub->num_seqR + 1;
      }
      
      //---------------------------------------------------------	
      //计算SeqL[]的SmBwt,生成buf_bwt_seq_L和buf_bwt_sum_L;
      //len_bwtseq_L = (num_seq_L+1)/2;
      //len_bwtsum_L = getbwtsumsize(num_seq_L);
      sub->num_ext = repNum;	
      fprintf(stderr, "%u, rep_Num = %d, out_buf[0] = %d, max_buf_size = %u\n", __LINE__, repNum, out_buf[0], max_buf_size);
      relatNxtCapIdx(sub,jmp_idx,mod_idx,cap_pos,fm_idx);
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //预编译期间临时屏蔽
      //0生成左Bwt，1生成右Bwt
      //fprintf(stderr, "444, buldSmBwt(sub,1) end\n"); 
      // 前面已经出现+++++++++++++++++++++++++
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
      //uint8_t test_buf0[5], test_buf1[5];
      uint8_t ch_buf[5], row;  
      if(sub->num_seqL > MIN_BWT_SIZE) {
        buldSmBwt(sub,0); //0生成左Bwt，1生成右Bwt
        //test_smbwt(sub, 0);
        uint32_t len_bwtL = (sub->num_seqL+1)/2*8; 
        fwrite(sub->bwt_seqL,sizeof(uint8_t),len_bwtL,fp_smbwt);
        uint32_t len_sumL = get_bwt_size(sub->num_seqL)-len_bwtL; 		
//fprintf(stderr, "[len_sumL]:num_seqL, len_sumL = %u %u\n", sub->num_seqL, len_sumL);
        if(len_sumL > 0){
          //以下当满足累加索引条件时调用＋＋＋＋＋＋＋＋＋＋＋＋＋＋＋
          gen_cnt_all(sub->bwt_seqL, sub->num_seqL, sub->bwt_sumL);
          fwrite(sub->bwt_sumL,sizeof(uint8_t), len_sumL,fp_smbwt);
        }
      } else{
        for(row = 0; row < sub->num_seqL; ++row){
          uint32_t buf = sub->sortL[row]; 
          ch_buf[0] = (uint8_t )((buf>>24)&0xFF);
          ch_buf[1] = (uint8_t )((buf>>16)&0xFF);
          ch_buf[2] = (uint8_t )((buf>>8)&0xFF);
          ch_buf[3] = (uint8_t )((buf)&0xFF);
          fwrite(ch_buf,sizeof(uint8_t),4,fp_smbwt);
        } 
      //log_seqs(fm_idx, sidx[0], sidx[0]+sidx[1], repNum);
      }
      //fwrite(sub->sortL,sizeof(uint8_t),sub->num_seqL*4,fp_smbwt);
      

      if(sub->num_seqR > MIN_BWT_SIZE) {
        buldSmBwt(sub,1); //0生成左Bwt，1生成右Bwtt
        //test_smbwt(sub, 1);
        uint32_t len_bwtR = (sub->num_seqR+1)/2*8; 
        uint32_t len_sumR = get_bwt_size(sub->num_seqR)- len_bwtR;
//fprintf(stderr, "[len_sumR]: num_seqR, len_sumR = %u, %u\n", sub->num_seqR, len_sumR);
        fwrite(sub->bwt_seqR,sizeof(uint8_t),len_bwtR,fp_smbwt);
        if(len_sumR > 0){
            //以下当满足累加索引条件时调用＋＋＋＋＋＋＋＋＋＋＋＋＋＋＋
          gen_cnt_all(sub->bwt_seqR, sub->num_seqR, sub->bwt_sumR);
          fwrite(sub->bwt_sumR,sizeof(uint8_t), len_sumR,fp_smbwt);
        }
      }else{
        for(row = 0; row < sub->num_seqR; ++row){
          uint32_t buf = sub->sortR[row]; 
          ch_buf[0] = (uint8_t )((buf>>24)&0xFF);
          ch_buf[1] = (uint8_t )((buf>>16)&0xFF);
          ch_buf[2] = (uint8_t )((buf>>8)&0xFF);
          ch_buf[3] = (uint8_t )((buf)&0xFF);
          fwrite(ch_buf,sizeof(uint8_t),4,fp_smbwt);
        } 
      }

    //-------------------------------------------------------
  //配对左序列seq_L[]与相关联的右序列seq_L[]的Indx区间
  //关联数据保存在buf_rel_idx[][]数组	
        //输出nxt_idx[i],nxt_idx[i],和apIdx[]数据
      fwrite(sub->nxt_idx,sizeof(uint32_t),sub->num_relat,fp_nxtpnt);
      fwrite(sub->nxt_flg,sizeof(uint8_t), sub->num_relat,fp_nxtflg);
    
//fwrite(sidx,sizeof(uint32_t)*2, 1,fp_extidx);
  
        
      cap[0].num_seqL  = sub->num_seqL;
      cap[0].num_seqR  = sub->num_seqR;
      cap[0].num_relat = sub->num_relat;
      fwrite(cap,sizeof(struct CapIfo),1,fp_capidx);	
      cap[0].nxtcap+= sub->num_relat;
      //cap[0].relat += sub->num_relat + sub->num_seqL+1;	 
      cap[0].relat += relat_len;	 
      cap[0].smbwt += get_bwt_size(sub->num_seqL) + get_bwt_size(sub->num_seqR);
      //计算relat数组的指针增量 
      sub->len_capidx++;
      //printf("repNum = %u, capidx = %u\n", repNum, sub->len_capidx); 
			idx++;
 		  //fprintf(stderr, "idx=%u\n", idx);	
		} // End: while(idx<len_seedidx) +++++++++++++++++++++++++++
    
    fclose(fp_nxtpnt);
		fclose(fp_nxtflg);
		fclose(fp_capidx);
		fclose(fp_relat);
		fclose(fp_smbwt);
//fclose(fp_extidx);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    FILE *fp_extidx;
    if((fp_extidx = fopen(f->extidx[repNum],"w"))== NULL){
        printf("can't open file = %s\n",f->extidx[repNum]);
        exit(1);
    }
    idx = 0;
    //uint32_t len_ext_idx = LenExtIdx[repNum+1];
    uint32_t len_ext_idx = out_buf[0];

    /*  
    while(idx < len_ext_idx){
        if(nxt_seedidx[idx][1] <= IS_SMLSIZ){
    idx++;
    continue;
      //小规模数据不做任何处理
      //比对过程中可以通过自身的jmp_idx,mod_idx信息，
      //可以找到相应的SmPos[]数组的行号
  }
        fwrite(nxt_seedidx+idx,sizeof(uint32_t)*2, 1,fp_extidx);
        ++idx;
    }
    */

    fwrite(nxt_seedidx,sizeof(uint32_t)*2, len_ext_idx, fp_extidx);
    fclose(fp_extidx);
//-------------------------------------------------------
    ord_seedidx = cur_seedidx;
		cur_seedidx = nxt_seedidx ;
		nxt_seedidx = ord_seedidx; 
    repNum++;
	}
/*  
    repNum = 0;
    FILE *fp;
    while(repNum < NUM_EXT) {

	    if((fp=fopen(f->seedidx[repNum+1],"rb"))==NULL){//判断是否打开文件
	        printf("can't open file = %s\n",f->seedidx[repNum+1]);
	        exit(0);
	    }
        uint32_t  numread=fread(nxt_seedidx,2*sizeof(uint32_t),LenExtIdx[repNum+1],fp);
	    fclose(fp);
        
        
        FILE *fp_extidx;
        if((fp_extidx = fopen(f->extidx[repNum],"w"))== NULL){
            printf("can't open file = %s\n",f->extidx[repNum]);
            exit(0);
        }
        uint32_t idx = 0;
        while(idx < LenExtIdx[repNum+1]){
        
            if(nxt_seedidx[idx][1] <= IS_SMLSIZ){
				idx++;
				continue;
					//小规模数据不做任何处理
					//比对过程中可以通过自身的jmp_idx,mod_idx信息，
					//可以找到相应的SmPos[]数组的行号
			}
            fwrite(nxt_seedidx+idx,sizeof(uint32_t)*2, 1,fp_extidx);
        }
        fclose(fp_extidx);
    }
*/
    set_config(MAX_RELAT);
    fprintf(stderr, "out of loop!\n");
    comFile(f,f);
	return 0;
}

int set_config(int MAX_RELAT)
{
    FILE *fp = fopen("config_file.txt", "w");
    fprintf(fp, "MAX_RELAT = %u\n", MAX_RELAT);
    fclose(fp);
}
void comFile(struct FileName *rf, struct FileName *wf){
	//+++++++++++++++++++++++++++++++++++++++++++++
printf("comFile() = beginning OK! ");
	FILE *fp_com;
	FILE *fp;
	int f_id;
	uint32_t num = 0;
	char  buf[2];
    int i,j;
    int fsize[NUM_FILES][NUM_EXT];
    int head[NUM_FILES][NUM_EXT];
    char *fName[NUM_FILES][NUM_EXT] ;

    for(i=0;i<NUM_EXT;i++){
    	fName[0][i]	=  rf->capidx[i];
	    fName[1][i]	=  rf->nxtpnt[i];
	    fName[2][i]	=  rf->nxtflg[i];
	    fName[3][i]	=  rf->relat[i];
	    fName[4][i]	=  rf->smbwt[i];
	    fName[5][i]	=  rf->extidx[i];
    }
    printf("comFile() = %s\n","461 row OK! ");
   
	for(i=0;i<NUM_EXT;i++){
		fsize[0][i] =  getFileSize(rf->capidx[i]);
		fsize[1][i] =  getFileSize(rf->nxtpnt[i]);
		fsize[2][i] =  getFileSize(rf->nxtflg[i]);
		fsize[3][i] =  getFileSize(rf->relat[i]);
		fsize[4][i] =  getFileSize(rf->smbwt[i]);		
		fsize[5][i] =  getFileSize(rf->extidx[i]);
	}

    printf("comFile() = %s\n","472 row OK! ");
	for(i=0;i<NUM_EXT;i++){
		//head[0][i] = 5;
		head[0][i] = ((fsize[0][i]+3)/4) ; 
		head[1][i] = ((fsize[1][i]+3)/4) ; 
		head[2][i] = ((fsize[2][i]+3)/4) ; 
		head[3][i] = ((fsize[3][i]+3)/4) ; 
		head[4][i] = ((fsize[4][i]+3)/4) ; 
		head[5][i] = ((fsize[5][i]+3)/4) ; 
	}

    printf("comFile() = %s\n","483  row OK! ");

	int num_ext = 0 ;
	while(num_ext<NUM_EXT){
		if((fp_com=fopen(wf->comfile[num_ext],"wb"))==NULL){
		    printf("can't open file : %s\n", wf->comfile[num_ext]);
		    exit(0);
		}
        uint32_t head_buf[NUM_FILES+1] = {};
       
        //for(i =0; i < NUM_FILES-1; ++i) {
        for(i =0; i < NUM_FILES; ++i) {
            head_buf[i] = head[i][num_ext];
            head_buf[NUM_FILES] += head_buf[i];
        }
        fwrite(head_buf, sizeof(uint32_t), NUM_FILES+1, fp_com);
        f_id = 0;
		//while(f_id<NUM_FILES-1){
		while(f_id<NUM_FILES){
   			//file_capidx数据输出到file_extidx
			if((fp=fopen(fName[f_id][num_ext],"rb"))==NULL){
			    printf("can't open file : %s\n",fName[f_id][num_ext]);
			    exit(0);
			}
			num = 0;
			while(fread(buf, sizeof(uint8_t), 1, fp) >0){
                num++;
				fwrite(buf,sizeof(uint8_t),1,fp_com);
			}
			if(num != fsize[f_id][num_ext]){
				fprintf(stderr, "file read Err\n");
				fprintf(stderr, "%u != %u, %u, %u, %s\n", num, fsize[f_id][num_ext], f_id, num_ext, fName[f_id][num_ext]);
                exit(0);
			}
            for(j=0; j<(head[f_id][num_ext]*4-fsize[f_id][num_ext]);j++ ){
		        buf[0] = 0;
				fwrite(buf,sizeof(uint8_t),1,fp_com);
			}
			fclose(fp);
			f_id++;
		}//End:  while(num_file<NumFiles) +++++++++++++++++++
		fclose(fp_com);
        num_ext++;
	}
	return;
} 
// End： void comDataFile(FileName *file) +++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void OpenSeedIdxfiles( struct FileName *f, uint32_t (*sidx)[2], uint64_t LenExtIdx[],uint32_t repNum,FILE *fpw[]){
	
printf("OpenSeedIdxfiles() %s\n","beginning ! ");

	FILE *fp;
	if((fp=fopen(f->seedidx[repNum+1],"rb"))==NULL){//判断是否打开文件
	    printf("can't open file = %s\n",f->seedidx[repNum+1]);
	    exit(0);
	}
	//uint32_t  numread=fread(sidx,3*sizeof(uint32_t),LenExtIdx[repNum+1],fp);
	uint32_t  numread=fread(sidx,2*sizeof(uint32_t),LenExtIdx[repNum+1],fp);
	fclose(fp);
	if((fpw[0]=fopen(f->relat[repNum],"w"))==NULL){
	    printf("can't open file = %s\n", f->relat[repNum] );
	    exit(0);
	}
	
	if((fpw[1]=fopen(f->smbwt[repNum],"w"))==NULL){
	    printf("can't open file = %s\n", f->smbwt[repNum]);
	    exit(0);
	}

	if((fpw[2]=fopen(f->nxtpnt[repNum],"w"))==NULL){
	    printf("can't open file = %s\n", f->nxtpnt[repNum]);
	    exit(0);
	}

	if((fpw[3]=fopen(f->nxtflg[repNum],"w"))==NULL){
	    printf("can't open file = %s\n", f->nxtflg[repNum]);
	    exit(0);
	}

	if((fpw[4]=fopen(f->capidx[repNum],"w"))==NULL){
	    printf("can't open file = %s\n", f->capidx[repNum]);
	    exit(0);
	}
/*    
	if((fpw[5] = fopen(f->extidx[repNum],"w"))== NULL){
	    printf("can't open file = %s\n",f->extidx[repNum]);
	    exit(0);
	}
*/
	return;
} // End: void Openfiles( ) ++++++++++++++++++++++++++++++++

//*************************************************
uint32_t getlen(uint64_t LenExtIdx[],struct FileName *f){
    int i,m_i;
    uint64_t max;
    m_i=0;
    max = 0 ;
    for(i=0;i<NUM_EXT+1;i++){
        //LenExtIdx[i] = getFileSize(f->seedidx[i])/12;
    	LenExtIdx[i] = getFileSize(f->seedidx[i])/8;
       if(max < LenExtIdx[i]){
       	 max = LenExtIdx[i];
       	 m_i = i;
       } 	
    }
	return m_i;
}
//*************************************************

/*
void CrtJmpMod(uint32_t cur_seedidx[][3], uint32_t* jmp_idx, 
	uint8_t *mod_idx, uint32_t *cap_pos){

	printf("CrtJmpMod() = %s\n","OK");

	return;
}
*/
int CrtJmpMod0(int n, uint32_t cur_seedidx[][2], uint32_t* jump_idx, uint8_t *mod_idx, uint32_t *cap_pos){

    int max_buf_size = 0;
    //uint32_t row_jump = 0;
    uint32_t last_idx_q = 0;
    jump_idx[0] = 0;
    //uint32_t jump_idx[LEN_NEXT_IDX+1]; 
    uint32_t row, i;
    uint32_t num, idx;
    uint32_t cap_buf[2]={};
    for(row=0; row<n; ++row){
        if(cur_seedidx[row][1] > max_buf_size) max_buf_size = cur_seedidx[row][1];
        //cap_pos[row] = cur_seedidx[row][2]; 
        idx = cur_seedidx[row][0];
        num = cur_seedidx[row][1];
         
        if(num <= IS_SMLSIZ){
            cap_pos[row] = cap_buf[0];
            cap_buf[0] += num;
        } else{
            cap_pos[row] = cap_buf[1]++;
        }
        uint32_t idx_q = idx/256;
        mod_idx[row] = idx%256;

        if(idx_q == last_idx_q){
            continue;
        } else{
            for(i= last_idx_q+1; i <= idx_q; ++i){
                jump_idx[i] = row; 
            }
            last_idx_q = idx_q;
        }

    } 
    //jump_idx[(n+255)/256] = n; 
    jump_idx[i] = n; 

//fprintf(stderr, "[%s, %u]: idx != 844453254 \n", __func__,__LINE__); 
    return max_buf_size;
}

int CrtJmpMod(int n, uint32_t cur_seedidx[][2], uint32_t* jump_idx, uint8_t *mod_idx, uint32_t *out_buf, int n_ext, idx_t *fm_idx){

    int max_buf_size = 0;
    //uint32_t row_jump = 0;
    uint32_t last_idx_q = 0, idx_q = 0;
    jump_idx[0] = 0;
    //uint32_t jump_idx[LEN_NEXT_IDX+1]; 
    uint32_t row, i;
    uint32_t num, idx, pos, k, l, num0 = 0;
    uint8_t seq[LEN_SEED];
    int row_j = 0;
    for(row=0; row<n; ++row){
        if(cur_seedidx[row][1] > max_buf_size) 
            max_buf_size = cur_seedidx[row][1];
        idx = cur_seedidx[row][0];
        num = cur_seedidx[row][1];
        //+++++++++++++++++++++++++++++++++
        if(num <= IS_SMLSIZ){ continue;}  
        if(n_ext == 0) {
            if(num > IS_SMLSIZ * NUM_SCALE){
                cur_seedidx[row_j][0] = cur_seedidx[row][0];
                cur_seedidx[row_j][1] = cur_seedidx[row][1];
                row_j++; 
            } else{
                continue;
            }
        } else if(n_ext > 0) {
            pos = bwt_sa(fm_idx->bwt, idx) + 16*n_ext;
            for(i = 0; i < LEN_SEED; ++i) seq[i] = __get_pac(fm_idx->pac, pos+i); 
            k = 0, l = fm_idx->bwt->seq_len;
            num0 = bwt_match_exact_alt(fm_idx->bwt, LEN_SEED, seq, &k, &l);
            if(num0 < num) {
                printf("num0 = %u, num = %u\n", num0, num);
                exit(1);  
            }
            if(num0 <= IS_SMLSIZ * NUM_SCALE) { continue;}
            cur_seedidx[row_j][0] = cur_seedidx[row][0];
            cur_seedidx[row_j][1] = cur_seedidx[row][1];
            row_j++;
        }  
        
        idx_q = idx/256;
        mod_idx[row_j-1] = idx%256;
        if(idx_q == last_idx_q){
            continue;
        } else{
            for(i= last_idx_q+1; i <= idx_q; ++i){
                jump_idx[i] = row_j-1; 
            }
            last_idx_q = idx_q;
        }
    } 
    //jump_idx[i] = row_j; 
    for(i= last_idx_q+1; i <= idx_q; ++i){
        jump_idx[i] = row_j-1; 
    }
    jump_idx[idx_q+1] = row_j;
    
    out_buf[0] = row_j;
    out_buf[1] = idx_q+1;
    return max_buf_size;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++
/*  */
int hier_cut_ext_seq(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, uint32_t (*seq_buf)[2])
{
    bwtint_t i, pos, st, ed;
    uint32_t seq16;
    uint32_t *ItoPos = bwt->sa;
    uint32_t c;
    uint32_t ii;
    uint32_t r = 0, r0;
    uint8_t ch, jj, kk;
    //fprintf(stderr, "\n"); 
    /* fetch 16 mer preseq in [k, l) and rm duplicate preseq */
    uint32_t max_pos, min_pos;
    int j = 0;
    if(offset<0){
        min_pos = 16;
        max_pos = bwt->seq_len;
    } else{
        min_pos = 0;
        max_pos = bwt->seq_len-offset-16;
    }
    seq16 =1;
    for(i = k; i < l; ++i){
        pos = bwt_sa(bwt, i);
        st = pos+offset; //ed = st+16;
        int flag;
        if(pos >=min_pos && pos<=max_pos){
            r =0;
            ch = pac[st/4];
            jj = (st%4)*2;
            kk = (8-jj); 
            r = (uint32_t )((0xFF>>jj)&ch);
            for(ii = st/4+1; ii <st/4+4 ; ++ii){
                r<<=8;
                r0 = (uint32_t)pac[ii];  
                r|= r0;
            }
            if(st%4>0){
                r <<=jj;
                r0 = (uint32_t)pac[ii];
                r0>>=kk;
                r|=r0;
            }
            seq16 = r;
            flag = 0;
        } else if(pos < min_pos){//if(st<0) 
            flag =1;
        } else if(pos > max_pos){//if ed > bwt->seq_len
            flag=1;
        }
        /*  
        if(offset < 0){
        	
            uint32_t __i;
        	uint32_t buf0, buf;
			buf = seq16;
        	buf0 =0;
        	for(__i = 0; __i < 16; ++__i){
        		buf0 <<=2;

        		buf0 = buf0|(buf&3);
        		buf >>=2;

        	}
        	seq_buf[j][0] = buf0;

        } else{
        	seq_buf[j][0] = seq16;
		}
        */
        seq_buf[j][0] = seq16;
        seq_buf[j++][1] = flag;       
        //printf("pos, st, seq16= %u, %u, %u\n", pos, st, seq16);
    }  
    return j;
}
/*  
int rem_euqal(vec_uint_t *seq_buf, bwtint_t k, vec_ext_t *rext)
{
    int i, n_uniq;
    uint32_t *b = seq_buf->a;
    
    ext_t *p = kv_pushp(ext_t, *rext);
    p->seq = b[0];
    p->idx = k;
 
    n_uniq = 1; 
    for(i = 1; i < seq_buf->n; ++i){
        if(b[i] > b[i-1]) {
            assert(rext->m >n_uniq*2+1);
            ext_t *p = kv_pushp(ext_t, *rext);
            p->seq = b[i];
            p->idx = k+i;
            //a[n_uniq*2] = b[i];
            //a[n_uniq*2+1] = k+i;          
            ++n_uniq;
        } 
    }
    return n_uniq;
}*/

/*
int filter_ext_seq(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, vec_uint_t *seq16s, bp_t *is_repeat)
{
    
    int seq16_st = seq16s->n;
    bwtint_t i, pos;
    uint32_t seq16;
    //fprintf(stderr, "\n"); 
  
    for(i = k; i < l; ++i){
        pos = bwt_sa(bwt, i);
        if(pos <16) continue; 
        seq16 = bns_extract_seq16(pac, pos-16, pos);//???
        //fprintf(stderr, "%u\t", seq16); 
        if(bp_get(is_repeat, seq16) == 1) continue;//skip if occured before
        bp_set1(is_repeat, seq16);
        kv_push(uint32_t, *seq16s, seq16);
    }
    ks_introsort(uint32_t, seq16s->n-seq16_st, seq16s->a+seq16_st);    
    //fprintf(stderr, "\n");
    //for(i = seq16_st; i <seq16s->n; ++i) fprintf(stderr, "%u\t", seq16s->a[i]);
    //fprintf(stderr, "\n");
    for(i=seq16_st; i < seq16s->n; ++i) bp_set0(is_repeat, seq16s->a[i]);    
    //++global_stat_20seed0[l-k+1]; 
    //++global_stat_36seed[seq16s->n-seq16_st]; 
    //global_stat_20seed1[l-k+1] += seq16s->n -seq16_st+1; 
    return seq16s->n;
}*/
void log_seqbuf(uint32_t n, uint32_t (*seqbuf)[2])
{
    uint32_t i, j;
    for(i=0; i< n; ++i){
        uint32_t seq = seqbuf[i][0];
        uint32_t idx = seqbuf[i][1];

        printf("%u\t", i);
        for(j=0;  j<16; ++j){
            printf("%u", (seq>>(30-j*2))&3); 
        }
        printf("\t%u\n", idx);
    }


}
void getBlckData(uint32_t sidx[2],uint32_t num_ext, struct SubBuf *sub, idx_t *idx){

    uint8_t *pac = idx->pac;
    bwt_t *bwt = idx->bwt;    

    uint32_t k = sidx[0], l = k+sidx[1];
    int i, j;
    uint32_t bg_idx, ed_idx;
    int n = sidx[1];


    int off_pos = LEN_SEED+2*16*num_ext; 
    int n_r = hier_cut_ext_seq(pac, bwt, k, l,  off_pos, sub->seqR_buf);//right seq buf
    //log_seqbuf(n, sub->seqR_buf);
    off_pos = -16; 
    int n_l = hier_cut_ext_seq(pac, bwt, k, l,  off_pos, sub->seqL_buf);//left seq buf
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//下面代码段完成以下几个功能：
//1. 过滤左段和右段不合理的数据导致的扩展。
//2. 生成右段的排序后的唯一序列。
//3. 生成右段唯一序列对应的index区间。
//4. 过滤后的左段序列集合中的第二列输入了相应的右段排序后序列的行号以便后面生成关联数组

    uint32_t flag=0, old_flag=1; 
    uint32_t row_l = 0, row_r = 0, seq; 
    for(i = 0; i < n; ++i ){
        flag = sub->seqR_buf[i][1];
        if(flag == 0 && old_flag == 0) {
            if(sub->seqR_buf[i][0] > seq ){
                ed_idx = i-1; 
                sub->sortR[row_r] = seq;
                sub->idxR[row_r][0] = bg_idx+sidx[0];
                sub->idxR[row_r][1] = ed_idx+sidx[0];
                for(j = bg_idx; j <= ed_idx; ++j){
                    if(sub->seqL_buf[j][1] == 0) {     
                        sub->seqL_buf[row_l][0] = sub->seqL_buf[j][0]; 
                        sub->seqL_buf[row_l][1] = row_r;
                        row_l++;
                    } 
                }
                row_r++; 
                seq = sub->seqR_buf[i][0];
                bg_idx = i;
            }         
        } else if(flag == 0 && old_flag ==1 ){
            seq = sub->seqR_buf[i][0];
            bg_idx = i;
        } else if(flag == 1 && old_flag == 0){
            ed_idx = i-1; 
            sub->sortR[row_r] = seq;
            sub->idxR[row_r][0] = bg_idx+sidx[0];
            sub->idxR[row_r][1] = ed_idx+sidx[0];
            for(j = bg_idx; j <= ed_idx; ++j){
                if(sub->seqL_buf[j][1] == 0) {     
                    sub->seqL_buf[row_l][0] = sub->seqL_buf[j][0]; 
                    sub->seqL_buf[row_l][1] = row_r;
                    row_l++; 
                }
            }
            row_r++; 
        }         
        old_flag = flag;
    } 
    if(flag==0){
        ed_idx = n-1; 
        sub->sortR[row_r] = seq;
        sub->idxR[row_r][0] = bg_idx+sidx[0];
        sub->idxR[row_r][1] = ed_idx+sidx[0];
        for(j = bg_idx; j <=ed_idx; ++j){
            if(sub->seqL_buf[j][1] == 0) {     
                sub->seqL_buf[row_l][0] = sub->seqL_buf[j][0]; 
                sub->seqL_buf[row_l][1] = row_r;
                row_l++; 
            }
        }
        row_r++; 
    } 
    sub->idxR[row_r][0] = ed_idx+sidx[0];
    sub->idxR[row_r][1] = ed_idx+sidx[0];

    uint32_t num_seqL = row_l;
    sub->num_seqR = row_r;
    //log_sub(sub);
	uint32_t (*seqL)[2] = sub->seqL_buf;
    ks_introsort(pair_t, num_seqL, (pair_t *)seqL); 
    uint32_t last_seq = seqL[0][0]; 
	sub->sortL[0] = last_seq;
	sub->L2rel[0] = 0;
    sub->relat[0] = seqL[0][1];	
	int last_relat = seqL[0][1];
    for(j =0,k=0,  i = 1; i < num_seqL; ++i){
		if(seqL[i][0] > last_seq){
            k++;
            j++;

            sub->sortL[j] = seqL[i][0];
            last_seq = seqL[i][0];
            
            sub->L2rel[j] = k;
            sub->relat[k] = seqL[i][1];
			last_relat = seqL[i][1];
        
        } else if(seqL[i][1]>last_relat){
			k++;
            sub->relat[k] = seqL[i][1];
			last_relat = seqL[i][1];
			
		}
	}

    ++j; ++k;
	sub->num_seqL = j;
	sub->num_relat = k;
    sub->L2rel[j] = k;

    uint32_t (*relat_buf)[2] = sub->seqL_buf;
    uint32_t *relat = sub->relat;
    uint32_t *L2rel = sub->L2rel;
    uint32_t *R2rel = sub->R2rel;
    //uint32_t *R_relat = sub->R_relat;
    for(i = 0; i < sub->num_relat; ++i){
        relat_buf[i][0] = relat[i]; 
        relat_buf[i][1] = i; 
    }
    
    ks_introsort(pair_t, sub->num_relat, (pair_t *)relat_buf); 
    j = 1;
    R2rel[0] = 0;
    int row = 0, num = 0;
    for(i = 0; i < sub->num_relat; ++i){
        //R_relat[i] = relat_buf[i][1];
        relat[i] = relat_buf[i][1];
        if(relat_buf[i][0] > row){
            R2rel[j] = num;
            row = relat_buf[i][0];
            j++;
        }             
        num++; 
    } 
    R2rel[j] = sub->num_relat;
    /*  
    for(i = 0; i < sub->num_relat; ++i){
        relat_buf[i][0] = relat_buf[i][1]; 
        relat_buf[i][1] = i; 
    }
    ks_introsort(pair_t, sub->num_relat, (pair_t *)relat_buf); 
    for(i = 0; i < sub->num_relat; ++i){
        relat[i] = relat_buf[i][1];  
    }
    */
    
    return;
}
/*
void build_local_bwt_alt(uint32_t *sort_seq, int l, int r, uint32_t *sa0, uint32_t *sa1, uint8_t *bwt)
{
    int i, j;
    uint32_t cnt_2nt[17];
    uint32_t *last_sa, *cur_sa;
    sort_seq += l;
    int n = r-l;
    for(j=0; j< 17; ++j) cnt_2nt[j] = 0;
    for(j = 0; j < n; ++j){
        uint32_t x = sort_seq[j]&0xF;
        __set_bwt(bwt, j*2+1, x&3);
        __set_bwt(bwt, j*2, (x>>2)&3);
        //__set_bwt2(bwt, j, x);
        ++cnt_2nt[x];
    }
    for(j =0; j < n; ++j) sa0[j] = j;
    last_sa = sa0; cur_sa = sa1;
    for(i = 1; i < 8; ++i){
        //fprintf(stderr, "\n==Iter = %u\n", i-1);

        //log_array(17, cnt_2nt);
        accumulate_cnt(17, cnt_2nt);
        //log_array(17, cnt_2nt);
        //log_array(r-l, last_sa);
        //log_array(r-l, cur_sa);
        //sort cur_sa 
        for(j = 0; j < n; ++j){
            //uint32_t x = (lext[last_sa[j]].seq>>((i-1)*4))&0xF;
            uint32_t x = __get_col(sort_seq[last_sa[j]], i-1);
            //fprintf(stderr, "%u->%u\n", x, cnt_2nt[x]);
            cur_sa[cnt_2nt[x]++]  = last_sa[j];
        }
        //generate bwt from cur_sa        
        for(j=0; j< 17; ++j) cnt_2nt[j] = 0;
        for(j = 0; j < n; ++j){
            //uint32_t x = (lext[cur_sa[j]].seq>>(i*4))&0xF;
            uint32_t x = __get_col(sort_seq[cur_sa[j]], i);
            //__set_bwt(bwt, i*n*2+j*2+1, x&3);
            //__set_bwt(bwt, i*n*2+j*2, (x>>2)&3);
            int n_2nt = (n+1)/2*2;
            __set_bwt2(bwt, i*n_2nt+j, x);
            ++cnt_2nt[x];
        }
        SWAP(uint32_t *, last_sa, cur_sa); 
    }
    fprintf(stderr, "\n==Iter = %u\n", i-1);
    //log_seq16array(n, last_sa, lext, i-1);
}
*/


//+++++++++++++++++++++++++++++++++++++++++++++++++++
void InitSeedIdx(struct SeedIdx *sIdx,uint32_t MaxLen, uint32_t ref_len){
	uint32_t (*cur_seedidx)[2];
	uint32_t (*nxt_seedidx)[2];
	uint32_t *jmp_idx;
	uint32_t *cap_pos;
	uint8_t  *mod_idx; 

	//if(NULL == (cur_seedidx = (uint32_t(*)[3])malloc(MaxLen*3*sizeof(uint32_t)))){
	if(NULL == (cur_seedidx = (uint32_t(*)[2])malloc(MaxLen*2*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    //if(NULL == (nxt_seedidx = (uint32_t(*)[3])malloc(MaxLen*3*sizeof(uint32_t)))){
    if(NULL == (nxt_seedidx = (uint32_t(*)[2])malloc(MaxLen*2*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    fprintf(stderr, "jmp_size = %u\n", (((ref_len+511)/256)*sizeof(uint32_t)));
    if(NULL == (jmp_idx = (uint32_t*)malloc(((ref_len+511)/256)*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    if(NULL == (mod_idx = (uint8_t*)malloc(MaxLen*sizeof(uint8_t)))){
	    perror("error...");
	    exit(1);
	}
    if(NULL == (cap_pos = (uint32_t*)malloc(MaxLen*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
	sIdx->cur_seedidx = cur_seedidx;
	sIdx->nxt_seedidx = nxt_seedidx;
	sIdx->jmp_idx     = jmp_idx;
	sIdx->mod_idx     = mod_idx;
	sIdx->cap_pos     = cap_pos;
	return;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//#define MAX_BLCK_SIZE 500000
#define MAX_SEED_SIZE 500000 
#define ALL_BUFS_SIZE 13  // All bufs size in SubBuf



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void InitSubBuf(struct SubBuf *sub){
	uint32_t *buf;
	int max_buf_size = sub->max_buf_size;
    if(NULL == (buf = (uint32_t*)malloc(max_buf_size*ALL_BUFS_SIZE*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
	
    sub->seqR_buf = (uint32_t(*)[2])buf;
	sub->seqL_buf = (uint32_t(*)[2])(buf + max_buf_size*2);

	sub->sortL    = buf + max_buf_size*4;
	sub->sortR    = buf + max_buf_size*5;

	sub->relat    = buf + max_buf_size*6;
	sub->L2rel    = buf + max_buf_size*7;
	sub->R2rel    = buf + max_buf_size*8;
    	
    sub->seqL_idxR= (uint32_t(*)[3])(buf + max_buf_size*9);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	sub->idxR     = (uint32_t(*)[2])buf; //注：idxR[][2]的数据存放在;seqL_buf[][2]的物理空间。
	sub->extIdx   = (uint32_t(*)[2])(buf + max_buf_size*2); //注：extIdx[][2]的数据存放在;seqR_buf[][2]的物理空间。
  	
    sub->pos_buf = buf + max_buf_size*10 ; 	
	sub->nxt_idx = buf + max_buf_size*11;
	sub->nxt_flg = (uint8_t*)(buf + max_buf_size*12);
    //sub->R_relat = buf + max_buf_size*13;
    uint8_t *buf_8;
	if(NULL == (buf_8 = (uint8_t*)malloc(max_buf_size*3*sizeof(uint8_t)))){
	    perror("error...");
	    exit(1);
	}
    uint16_t *buf_16;
	if(NULL == (buf_16 = (uint16_t*)malloc(max_buf_size*3*sizeof(uint16_t)))){
	    perror("error...");
	    exit(1);
	}
    uint32_t *buf_32;
	if(NULL == (buf_32 = (uint32_t*)malloc(max_buf_size*3*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    
    sub->buf_8 = buf_8;
    sub->buf_16 = buf_16;
    sub->buf_32 = buf_32;

//+++++++++++++++++++++++++++++++++++
	uint8_t *bwt;
	if(NULL == (bwt = (uint8_t*)malloc(12*max_buf_size*sizeof(uint8_t)))){
	    perror("error...");
	    exit(1);
	}
	sub->smbwt    = bwt;
	sub->bwt_seqL = bwt ;   
	sub->bwt_sumL = bwt + 4*max_buf_size ;
	sub->bwt_seqR = bwt + 6*max_buf_size ;    
	sub->bwt_sumR = bwt + 10*max_buf_size ;  

//+++++++++++++++++++++++++++++++++++++++++++++++++
	sub->num_seqL = 0;
	sub->num_seqR = 0;
	sub->num_relat= 0;	
	//sub->num_ext  = 0;

	//sub->off_pos  = 0; 
    //sub->blck_id  = 0; 

	sub->len_capidx = 0;
	sub->len_relat  = 0;	
	sub->len_smbwt  = 0;	
	sub->len_nxtpnt = 0;
	sub->len_nxtflg = 0;
	sub->len_extidx = 0;
	sub->num_pos    = 0;
	sub->num_bwtL   = 0;
	sub->num_bwtR   = 0;


	sub->cap[0].relat    = 0;	
	sub->cap[0].smbwt    = 0;
    sub->cap[0].nxtcap   = 0;	
	sub->cap[0].num_seqL = 0;
	sub->cap[0].num_seqR = 0;
	sub->cap[0].num_relat= 0;	
	//sub->cap[0].num_pos  = 0;

//+++++++++++++++++++++++++++++++++++++++++++++++++
	return ;
}
#define __set_bwt(bwt, l, c) ((bwt)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define __get_bwt(bwt, l) ((bwt)[(l)>>2]>>((~(l)&3)<<1)&3)
#define __set_bwt2(bwt, l, c) ((bwt)[(l)>>1] |= (c)<<((~(l)&1)<<2))
#define __get_bwt2(bwt, l) ((bwt)[(l)>>1] >> ((~(l)&1)<<2)&15)


#define __get_col(seq, i) (((seq)>>((i)*4))&0xF)
static inline void accumulate_cnt(int n, int *cnt_2nt)
{
    int i;
    for(i = 1; i < n-1; ++i) cnt_2nt[i] += cnt_2nt[i-1]; 
    for(i = n-1; i >0; --i) cnt_2nt[i] = cnt_2nt[i-1];
    cnt_2nt[0] = 0;
 
}
static inline void accumulate_cnt_uint8(int n, uint8_t *cnt_2nt)
{
    int i;
    for(i = 1; i < n-1; ++i) cnt_2nt[i] += cnt_2nt[i-1]; 
    for(i = n-1; i >0; --i) cnt_2nt[i] = cnt_2nt[i-1];
    cnt_2nt[0] = 0;


}


void build_local_bwt_alt(uint32_t *sort_seq, int l, int r, uint32_t *sa0, uint32_t *sa1, uint8_t *bwt)
{
    int i, j;
    uint32_t cnt_2nt[17];
    uint32_t *last_sa, *cur_sa;
    sort_seq += l;
    int n = r-l;
    for(j=0; j< 17; ++j) cnt_2nt[j] = 0;
    for(j = 0; j < n; ++j){
        uint32_t x = sort_seq[j]&0xF;
        __set_bwt(bwt, j*2+1, x&3);
        __set_bwt(bwt, j*2, (x>>2)&3);
        //__set_bwt2(bwt, j, x);
        ++cnt_2nt[x];
    }
    for(j =0; j < n; ++j) sa0[j] = j;
    last_sa = sa0; cur_sa = sa1;
    for(i = 1; i < 8; ++i){
        //fprintf(stderr, "\n==Iter = %u\n", i-1);

        //log_array(17, cnt_2nt);
        accumulate_cnt(17, cnt_2nt);
        //log_array(17, cnt_2nt);
        //log_array(r-l, last_sa);
        //log_array(r-l, cur_sa);
        //sort cur_sa 
        for(j = 0; j < n; ++j){
            //uint32_t x = (lext[last_sa[j]].seq>>((i-1)*4))&0xF;
            uint32_t x = __get_col(sort_seq[last_sa[j]], i-1);
            //fprintf(stderr, "%u->%u\n", x, cnt_2nt[x]);
            cur_sa[cnt_2nt[x]++]  = last_sa[j];
        }
        //generate bwt from cur_sa        
        for(j=0; j< 17; ++j) cnt_2nt[j] = 0;
        for(j = 0; j < n; ++j){
            //uint32_t x = (lext[cur_sa[j]].seq>>(i*4))&0xF;
            uint32_t x = __get_col(sort_seq[cur_sa[j]], i);
            //__set_bwt(bwt, i*n*2+j*2+1, x&3);
            //__set_bwt(bwt, i*n*2+j*2, (x>>2)&3);
            int n_2nt = (n+1)/2*2;
            __set_bwt2(bwt, i*n_2nt+j, x);
            ++cnt_2nt[x];
        }
        SWAP(uint32_t *, last_sa, cur_sa); 
    }
    //fprintf(stderr, "\n==Iter = %u\n", i-1);
    //log_seq16array(n, last_sa, lext, i-1);
}
uint32_t get_bwt_size(uint32_t n_data)
{
    uint32_t bwt_size = (n_data+1)/2*8; 
    int siz_cnt = n_data>254*256?4:2;
    int rng_num = (n_data+253)/254;
    int len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;
    uint32_t cnt_size = 16*8*len_sum;
    //uint32_t size = 16*8*len_sum+seq_size;
    if(n_data <= MIN_BWT_SIZE) {
        return n_data*sizeof(uint32_t); 
    }
    //++++++++++++++++++++++++++++++++++++++++++
    if(n_data <= NO_BWT_SUM_SIZE) cnt_size = 0;
    else if(n_data <= 255) cnt_size = 17*8;
 
    return cnt_size + bwt_size; 
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
        //log_array2(17, rot, cnt_2nt);
        accumulate_cnt_uint8(17, cnt_2nt[rot]);
        //log_array2(17, rot, cnt_2nt);

    }
}


void gen_cnt_all(const uint8_t* bwt, int n_data, uint8_t *cnt)
{
    //fprintf(stderr, "[gen_cnt_all]\n");
    if(n_data<=NO_BWT_SUM_SIZE) {
        fprintf(stderr, "Error in gen_cnt_all!!\n"); 
        exit(1);
    } else if(n_data <= 255) gen_cnt(bwt, n_data, (uint8_t (*)[17])cnt);
    else gen_cnt2(bwt, n_data, cnt);
    return;
}
//+++++++++++++++++++++++++++++++
/*
if(n_data==134){
	int i, rot;
	fprintf(stderr, "\n");

	for(rot =0 ; rot < 8; ++rot){
		for(i = 0; i < 17; ++i) {
			fprintf(stderr, "%u\t", cnt[rot*17+i]);
		}	
		fprintf(stderr, "\n");
	
	}
}*/
/*
if(n_data == 419){
int DataNum = n_data;
int siz_cnt; 
uint32_t __i, __j, __k, __l;
fprintf(stderr, "\n----------------------\n");
if( DataNum <= 254*256 ) siz_cnt = 2 ;
if( DataNum > 254*256 ) siz_cnt = 4 ;
int rng_num = ((DataNum+253) /254) ;
int len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        fprintf(stderr, "gen_cnt_all = %u   \t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            fprintf(stderr, "%u\t", cnt[__j+__l]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "-------------------\n\n");
}
fprintf(stderr, "\n+++++++++++++++++\n");
}*/
//+++++++++++++++++++++++++++++++
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
                
            //while( rng_idx <= rng_num ){
            while( rng_idx < rng_num ){
                rng_off = (8+ siz_cnt) * (rng_idx/8) + rng_idx%8  ;
                sum_buf +=  pSum[ch_off + rng_off];
                rng_idx++;
                if(rng_idx%8==0){
                    int buf = sum_buf ;
                    for(i=0; i<siz_cnt; i++){
                        //fprintf(stderr, "siz_cnt0 %u-%u\t", buf, buf & 255);
						pSum[ch_off + rng_off+(siz_cnt-i)] = (uint8_t)(buf & 255) ;

                        buf >>= 8;
                      
                    }
                }
               
            } // End : +++++++++++++++++++while( rng_idx <= rng_num )
            glb_sum[rot][ch] = (ch==0?0:glb_sum[rot][ch-1])+sum_buf;
        } // End : +++++++++++++++++++++  for( ch = 0 ; ch <16 ; ch++)
        for( ch = 0 ; ch <16 ; ch++){
            ch_off = rot_sum +ch*len_sum +(len_sum -1) ;
            int buf = glb_sum[rot][ch] ;
            for(i=0; i<siz_cnt; i++){
				pSum[ch_off -i] = (uint8_t)(buf & 255) ;
                //fprintf(stderr, "siz_cnt0 %u\t", buf & 255);
                buf >>= 8;                      
            }
        }
    }//end for(rot=0; rot <8; ++rot)
//++++++++++++++++++++++++++++++++++++++
/*
if(n_data == 3701){
int DataNum = n_data;
int siz_cnt; 
uint32_t __i, __j, __k, __l;
fprintf(stderr, "\n----------------------\n");
if( DataNum <= 254*256 ) siz_cnt = 2 ;
if( DataNum > 254*256 ) siz_cnt = 4 ;
int rng_num = ((DataNum+253) /254) ;
int len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        fprintf(stderr, "gen_cnt2 =  %3u ", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            fprintf(stderr, "%3u ", cnt[__j+__l]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "-------------------\n\n");
}
fprintf(stderr, "\n+++++++++++++++++\n");
}
*/
	return;
}

//++++++++++++++++++++++++++++++++++++++
/*
if(n_data == 419){
int DataNum = n_data;
int siz_cnt; 
uint32_t __i, __j, __k, __l;
fprintf(stderr, "\n----------------------\n");
if( DataNum <= 254*256 ) siz_cnt = 2 ;
if( DataNum > 254*256 ) siz_cnt = 4 ;
int rng_num = ((DataNum+253) /254) ;
int len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        fprintf(stderr, "gen_cnt2 =  %4u   \t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            fprintf(stderr, "%u\t", cnt[__j+__l]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "-------------------\n\n");
}
fprintf(stderr, "\n+++++++++++++++++\n");
}
*/
//++++++++++++++++++++++++++++++++++++++
//生成序列
uint32_t __dna2_count(uint8_t *bwt, uint32_t k , uint32_t l, uint8_t c){
    uint32_t i, n;
    for(n=0, i =k; i <l; ++i){
        n += (__get_bwt2(bwt, i)==c)?1:0;
    }
    return n;

}

void test_smbwt(struct SubBuf *sub, int flg)
{
    uint8_t *seq_buf, *seq_buf0, *bwt;
    uint32_t num, *sort_seq;
    if(flg == 0){
        bwt = sub->bwt_seqL;
        num = sub->num_seqL;
        sort_seq = sub->sortL;
    } else if(flg == 1){
        bwt = sub->bwt_seqR;
        num = sub->num_seqR;
        sort_seq = sub->sortR;
    }
    seq_buf = malloc(16*num*2);  
    test_AlgnBwtSml(bwt, num, seq_buf);
    int i, j;
    seq_buf0 = seq_buf+num*16;
    for(i = 0; i < num; ++i){
        uint8_t ch;
        for(j = 0; j < 16; ++j){
            seq_buf0[i*16+j] = (sort_seq[i]>>((15-j)*2))&3;      
        }
    }

    for(i = 0; i < num; ++i){
        for(j = 0; j < 16; ++j){
            //printf("[%u, %u]: %u %u\n", i, j, seq_buf[i*16+j], seq_buf0[i*16+j]);
            if(seq_buf[i*16+j] != seq_buf0[i*16+j]){
                fprintf(stderr, "test_smbwt error!!!!\n"); 
                exit(1);
            
            }  
        }
       
    
    
    }
    free(seq_buf);
    return;
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


void buldSmBwt(struct SubBuf *sub,int flg){
    uint8_t *bwt; uint32_t *sa0, *sa1;
    //bwt = (uint8_t *)calloc((n+1)/2*8, sizeof(uint8_t));  
    int n;
    uint32_t *sort_seq;
    //bwt = flg==0?sub->bwt_seqL:sub->bwt_seqR;
    if(flg == 0 ){
        bwt = sub->bwt_seqL; 
        n = sub->num_seqL; 
        sort_seq = sub->sortL;
    } else{
        bwt = sub->bwt_seqR; 
        n = sub->num_seqR;
        sort_seq = sub->sortR;
    }
    //sa0 = (uint32_t *)calloc(n, sizeof(uint32_t));
    //sa1 = (uint32_t *)calloc(n, sizeof(uint32_t));
    sa0 = (uint32_t *)sub->seqL_buf;
    sa1 = (uint32_t *)sub->seqR_buf;
    
    memset(bwt, 0, (n+1)/2*8);
    build_local_bwt_alt(sort_seq, 0, n, sa0, sa1, bwt); 
    //free(sa0);free(sa1);
    //build_occ(bwt, n); 
    //dump bwt
        
       
 
        

   

  
	return;
}


//**************************************************
//*************************************************************

void relatNxtCapIdx(struct SubBuf *subBuf,uint32_t *jmp_idx, uint8_t *mod_idx,uint32_t *cap_pos, idx_t *fm_idx)
{
	uint32_t  i, j,	bgn,end, num,  bgn_q, end_q, i_r, seq, pos, row, j_idx;
  uint32_t buf_idx[5];
  struct SubBuf *sub;
	sub = subBuf;
  sub->num_pos = 0;
	for(i=0;i<sub->num_seqL;i++){
		bgn = sub->L2rel[i];
		end = sub->L2rel[i+1];
		for(j=bgn; j<end; j++){
			//row = sub->relat[j];
			//sub->seqL_idxR[row][0] = sub->sortL[i];
			sub->seqL_idxR[j][0] = sub->sortL[i];
		}				
	}
	for(i=0;i<sub->num_seqR;i++){
		bgn = sub->R2rel[i];
		end = sub->R2rel[i+1];
		for(row=bgn; row<end; row++){
			j = sub->relat[row];
            //sub->seqL_idxR[row][1] = sub->idxR[i][0]; 
			//sub->seqL_idxR[row][2] = sub->idxR[i][1];						
		  sub->seqL_idxR[j][1] = sub->idxR[i][0]; 
			sub->seqL_idxR[j][2] = sub->idxR[i][1];		
    }				
	}

  fprintf(stderr, "%u, num_seqL = %u, num_seqR = %d, num_relat = %d\n", __LINE__, sub->num_seqL, sub->num_seqR, sub->num_relat);	
  for(i=0;i<sub->num_relat;i++){
    seq = sub->seqL_idxR[i][0];
    buf_idx[0] = sub->seqL_idxR[i][1];
    buf_idx[1] = sub->seqL_idxR[i][2];
    fprintf(stderr, "i = %u, seq = %x, seqL_idxR[1] = %u, seqL_idxR[2] = %u\n", i, seq, sub->seqL_idxR[i][1], sub->seqL_idxR[i][2]);	
    uint8_t ch[200];
    /*  
    for(j=0;j<16;j++){
        ch[j] = seq&3;
        seq >>= 2;
    }
    */
    for(j=0;j<16;j++){
      ch[15-j] = seq&3;
      seq >>= 2;
    }
    int is_aln = bwt_match_exact_alt(fm_idx->bwt, 16, ch, &buf_idx[0], &buf_idx[1]);
    if(is_aln <1 ){ 
      printf("[%u, error]:is_aln = %d, k = %u, l = %u\n", __LINE__, is_aln, buf_idx[0], buf_idx[1]);
        //exit(1);
    }
    sub->extIdx[i][0]= buf_idx[0]; //返回的bgnIdx
    sub->extIdx[i][1]= buf_idx[1]; //返回的endIdx
    fprintf(stderr, "i = %u, bg_idx = %u, ed_idx = %u, pos = %u\n", i, sub->extIdx[i][0], sub->extIdx[i][1], bwt_sa(fm_idx->bwt, sub->extIdx[i][0]));	
  //test code测试代码
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++
   
    uint32_t tmp[3];
    uint32_t pos = bwt_sa(fm_idx->bwt, buf_idx[0]);
    uint32_t __l =  LEN_SEED+(sub->num_ext+1)*32;
    if(pos >fm_idx->bwt->seq_len-__l) {
      fprintf(stderr, "pos out range!\n");
      exit(1);
    }
    int __i;
    for(__i= 0; __i < __l; ++__i){
      ch[__i] = __get_pac(fm_idx->pac, pos+__i); 
    } 
  /*
    printf("lseq1: ");
    for(j = 0; j < 16; ++j) printf("%2u ", ch[j]);
    printf("\n");
  if(is_aln <1 ){ 
        //printf("[%u, error]:is_aln = %d, k = %u, l = %u\n", __LINE__, is_aln, buf_idx[0], buf_idx[1]);
        exit(1);
    }
  */
    tmp[0] = 0; tmp[1] = fm_idx->bwt->seq_len; 
    tmp[2] = bwt_match_exact_alt(fm_idx->bwt, __l, ch, &tmp[0], &tmp[1]);

    //fprintf(stderr, "[Test all]: buf=(%u, %u, %u), tmp = (%u, %u, %u)\n", buf_idx[0], buf_idx[1], is_aln, tmp[0], tmp[1], tmp[2]); 
  /*  
    if(tmp[2] < 1){
        fprintf(stderr, "[Test ]: not align buf=(%u, %u)\n", buf_idx[0], buf_idx[1]); 
        exit(1);
    } else{
        if(tmp[0] != buf_idx[0] || tmp[1] != buf_idx[1]){
            fprintf(stderr, "[Test %u]: buf=(%u, %u, %u), tmp = (%u, %u, %u)\n", __l, buf_idx[0], buf_idx[1], is_aln, tmp[0], tmp[1], tmp[2]); 
            
        } 
    }
    */
    if(tmp[0] != buf_idx[0] || tmp[1] != buf_idx[1]){

      fprintf(stderr, "[Test %u]: buf=(%u, %u, %u), tmp = (%u, %u, %u)\n", __l, buf_idx[0], buf_idx[1], is_aln, tmp[0], tmp[1], tmp[2]); 
      uint32_t __k0 = 0, __l0 = fm_idx->bwt->seq_len; 
      tmp[2] = bwt_match_exact_alt(fm_idx->bwt, __l-32, ch+16, &__k0, &__l0);
      fprintf(stderr, "[Test]: %u, %u\n", __k0, __l0); 
      for(j = __k0; j <= __l0; ++j) {
        pos = bwt_sa(fm_idx->bwt, j);
        fprintf(stderr, "pos[%u] out of range, %u\n", j, pos); 
        if(pos <= 16 || pos >= fm_idx->bwt->seq_len-16){
          fprintf(stderr, "pos out of range, %u\n", pos); 
        } 
      }
      exit(1); 
    } 
  }
	sub->num_pos = 0;
  for(i=0;i<sub->num_relat;i++){
		bgn = sub->extIdx[i][0];
		end = sub->extIdx[i][1];
		num = end - bgn + 1;
    if(num==1){
			pos = bwt_sa(fm_idx->bwt, bgn);
      sub->nxt_idx[i] = pos;
			sub->nxt_flg[i] = num;
			continue;
		} else if(num <= IS_SMLSIZ){
      sub->nxt_idx[i] = bgn;           
      sub->nxt_flg[i] = num;
      continue;
    }
		bgn_q = jmp_idx[bgn/256];
		end_q = jmp_idx[(bgn/256)+1];
    i_r = bgn%256;
		for(j_idx=bgn_q; j_idx<end_q; j_idx++){
			if(i_r == mod_idx[j_idx]){
				break;
			} 
		}
    if(j_idx == end_q){
        fprintf(stderr,"[repeat_num = %u]: no j_idx!\n", sub->num_ext);
        fprintf(stderr,"bg_idx, ed_idx= %u, %u!\n", bgn, end);
        fprintf(stderr,"bgn_q, end_q= %u, %u!\n", jmp_idx[bgn/256], jmp_idx[bgn/256+1]);
        exit(1); 
    }
    sub->nxt_idx[i] = j_idx; //指向下一级CapIdx[]的行号
		if(num < 255){
			sub->nxt_flg[i] = (uint8_t)num ; //指向CapIdx[]的类型
 		}else{
			sub->nxt_flg[i] = 0xFF ;
		}
  } // End : for(i=0;i<len_relat;i++)
	return;
}

