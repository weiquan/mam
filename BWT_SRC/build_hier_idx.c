
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
#include "rbwt.h"

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
//#define  MAX_EXT  (((LEN_READ- LEN_SEED )/2)+OFF_SEED+ LEN_EXT-1)/LEN_EXT  
#define  MAX_SEED_NUM 600000
#define  LEN_FILE_NAME 100 

#define  NUM_SCALE 1
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

  int32_t (*relat_sai)[2]; // sortR[]的每一个序列的(bgnIdx,ednIdx)
	//注：extIdx[][2]的数据存放在;seqR_buf[][2]的物理空间。	
	uint32_t (*seqL_idxR)[3];  // relat[]的每个行上左序列和右序列对应的idx,即(seqL,bgnIdxR,endIdxR)	

	uint32_t *pos_buf;  // 保存当前快中的pos数据； 	
	uint32_t *nxt_idx;  // 保存当前快中的nxt_idx数据；
	uint8_t  *nxt_flg;  // 保存当前快中的nxt_flg数据；


	rbwt_t  *bwt_seqL;   
	uint8_t  *bwt_sumL; 
	rbwt_t  *bwt_seqR;   
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
	uint32_t (*last_seedidx)[2];
	uint32_t (*nxt_seedidx)[2];
	int n_jmp;
	int n_mod;
	int n_cap;

  uint32_t *jmp_idx;
	uint32_t *cap_pos;
	uint8_t  *mod_idx; 
};


struct SeedIdx* seedidx_init(uint32_t max_seeds, uint32_t ref_len);
//void setFileName(char *in_fname, struct FileName *out_fname);
uint32_t getlen(uint64_t len_seedidx[],struct FileName *inf);

struct SubBuf *SubBuf_init(int max_buf_size);
int gen_jmpidx(int n, uint32_t cur_seedidx[][2], uint32_t* jump_idx, uint8_t *mod_idx, uint32_t *out_buf, int n_ext, idx_t *fm_idx);
void seedidx_open(struct  FileName *f, uint32_t (*sidx)[2], uint64_t len_seedidx[],uint32_t exti,FILE *fpw[]);
void gen_relat(uint32_t bgn,  uint32_t num, uint32_t num_ext, struct SubBuf *sub_buf, idx_t *idx);
void gen_nxt_cap(struct SubBuf *subBuf,int n_jmp, uint32_t *jmp_idx, int n_mod, uint8_t *jmp_mod,uint32_t *cap_pos,idx_t *fm_idx);
void comFile(struct FileName *rf, struct FileName *wf);
uint64_t getFileSize(const char *file);








uint32_t get_bwt_size(uint32_t n_data);



int __test_sort(int n, const uint32_t *sort_seq)
{
  int i;
  for(i = 1; i <n; ++i ){
    if(sort_seq[i] <= sort_seq[i-1]) {
      return i;
    } 
  }
  return -1;
}
void log_sub(struct SubBuf *sub)
{ 
  int i, j;
  printf("==========n_seqL = %u==================\n", sub->num_seqL);
  for(i=0; i< sub->num_seqL; ++i){
    uint32_t seq = sub->sortL[i];
    for(j=0;  j<16; ++j){
      printf("%u", (seq>>(30-j*2))&3); 
    }
    printf("\t%u\n", sub->L2rel[i]);
  }
  printf("\n-------n_relat = %u------\n", sub->num_relat);
  for(i=0; i< sub->num_relat; ++i){
    printf("%u\n", sub->relat[i]); 
  } 
  printf("\n-------n_seqR = %u--------\n", sub->num_seqR);
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

// build hierarchical index
// h1 index for 20bp seeds
// h2 index for 52bp seeds
// h3 index for 84bp reads
//...
int build_hier_idx(const char* prefix)
{
  bwtint_t i;
  fprintf(stderr, "[%s]:  restore fm_index %s\n", __func__, prefix);
  idx_t *idx = fbwt_fmidx_restore(prefix);
  
  //+++++++++++++++++++++++++++++++++++++++++++
	//初始化
  FileName *f = setFileName(prefix);
  uint64_t len_seedidx[MAX_EXT+1] ;
	uint32_t max_seeds; 
	i = getlen(len_seedidx, f);  
  max_seeds = len_seedidx[i];
  //内存初始化------------------------------------------------------
  uint32_t (*cur_seedidx)[2], (*last_seedidx)[2], (*nxt_seedidx)[2];
	uint32_t *jmp_idx, *cap_pos;
	uint8_t  *mod_idx = NULL; 
	
  fprintf(stderr, "[%s]:  generate seed index %s\n", __func__);
  struct SeedIdx* sIdx = seedidx_init(max_seeds, idx->bwt->seq_len);
	cur_seedidx = sIdx->cur_seedidx;
	nxt_seedidx = sIdx->nxt_seedidx;
	jmp_idx     = sIdx->jmp_idx;
	mod_idx     = sIdx->mod_idx;
	cap_pos     = sIdx->cap_pos;
  //数组初始化----------------------------------------------------
	FILE *fp = xopen(f->seedidx[0],"rb");
  fread(cur_seedidx, 2*sizeof(uint32_t), len_seedidx[0], fp);
  err_fclose(fp);
	//计算cur_seedidx的jmp_idx,mod_idx, cap_pos数组，并输出到FlgIdxFile_1
	uint32_t out_buf[2];    
  int max_buf_size = gen_jmpidx(len_seedidx[0], cur_seedidx, jmp_idx, mod_idx,out_buf, 0, idx);

  //struct SubBuf *sub= (struct SubBuf*)malloc(sizeof(struct SubBuf));
  //sub->max_buf_size = max_buf_size;
  struct SubBuf *sub = SubBuf_init(max_buf_size);
  struct CapIfo  *cap = sub->cap;
  fp = xopen(f->jmpmod, "wb");
  
  fwrite(&sIdx->n_jmp, sizeof(uint32_t), 1, fp);
  uint32_t n_mod = (out_buf[0]+3)/4*4;
  fwrite(&n_mod, sizeof(uint32_t), 1, fp);
  fwrite(&sub->max_buf_size, sizeof(uint32_t), 1, fp);
	fwrite(jmp_idx, sizeof(uint32_t), sIdx->n_jmp, fp);  
	fwrite(mod_idx, sizeof(uint8_t) , n_mod, fp);  
  err_fclose(fp); 
  int num_ext = MAX_EXT;
  bwtint_t max_relat = 0;
  int exti;
  for(exti= 0; exti < num_ext; ++exti){
    int tot_relat = 0;
    sub->len_capidx = 0;	
    cap[0].relat = 0; 
    cap[0].smbwt = 0; 
    cap[0].nxtcap = 0; 
    FILE *fp_seedidx = xopen(f->seedidx[exti+1], "rb");
    fread(nxt_seedidx, 2*sizeof(uint32_t), len_seedidx[exti+1], fp_seedidx);
    err_fclose(fp_seedidx);

    FILE *fp_relat = xopen(f->relat[exti],"w");
    FILE *fp_smbwt = xopen(f->smbwt[exti],"w");
    FILE *fp_nxtpnt = xopen(f->nxtpnt[exti],"w");
    FILE *fp_nxtflg = xopen(f->nxtflg[exti],"w");
    FILE *fp_capidx = xopen(f->capidx[exti],"w");

    uint32_t l_seedidx = out_buf[0];
    fprintf(stderr, "[%s]:  %u, fn_relat[%d] = %s, l_seedidx = %d\n", __func__, __LINE__, exti, f->relat[exti], l_seedidx);
    fwrite(&l_seedidx, sizeof(uint32_t), 1, fp_relat);
    fprintf(stderr, "[%s]:  open seedidx finish!\n", __func__);
    bwtint_t bgn, end, num;
    //fprintf(stderr, "[%s]:  %u, out_buf[0] = %d, max_buf_size = %u\n", __func__, __LINE__, out_buf[0], max_buf_size);
    gen_jmpidx(len_seedidx[exti+1], nxt_seedidx, jmp_idx, mod_idx, out_buf, exti+1, idx);
    //fprintf(stderr, "[%s]:  %u, ext_num = %d, out_buf[0] = %d, max_buf_size = %u\n", __func__, __LINE__, exti, out_buf[0], max_buf_size);
    bwtint_t si;
    for(si = 0; si < l_seedidx; ++si){
      bgn = cur_seedidx[si][0];
      num = cur_seedidx[si][1];
      end = bgn + num-1;
      if(num <= IS_SMLSIZ){ continue; } 
      //fprintf(stderr, "[%s]:  %u, generate relation array for left-seqs and right-seqs info.\n",__func__,  __LINE__);
      //fprintf(stderr, "[%s]:  %u, bgn = %u, end = %u.\n",__func__,  __LINE__, bgn, end);
      gen_relat(bgn, num, exti, sub, idx);  
      bwtint_t bg, ed, relat_len;           
      if(sub->num_relat >max_relat) max_relat = sub->num_relat;
      fwrite(&bgn, sizeof(uint32_t), 1, fp_relat);
      fwrite(&num, sizeof(uint32_t), 1, fp_relat);
      int x = sub->num_seqL +1;
      fwrite(&x, sizeof(uint32_t), 1,fp_relat);
      fwrite(sub->L2rel,sizeof(uint32_t), x, fp_relat);
      x = sub->num_seqR +1;
      fwrite(&x, sizeof(uint32_t), 1,fp_relat);
      //fwrite(sub->R2rel,sizeof(uint32_t), x, fp_relat);
      fwrite(&tot_relat,sizeof(uint32_t), 1,fp_relat);
      fwrite(&sub->num_relat,sizeof(uint32_t), 1,fp_relat);
      fwrite(sub->relat,sizeof(uint32_t),sub->num_relat,fp_relat);
      tot_relat += sub->num_relat;
      relat_len = sub->num_relat+sub->num_seqL+sub->num_seqR+2;
      sub->num_ext = exti;	
      //fprintf(stderr, "[%s]:  %u, generate next generation info.\n", __func__, __LINE__);
      gen_nxt_cap(sub, sIdx->n_jmp, jmp_idx, sIdx->n_mod, mod_idx, cap_pos, idx);
      fwrite(sub->nxt_idx,sizeof(uint32_t),sub->num_relat, fp_relat);
      fwrite(sub->nxt_flg,sizeof(uint8_t), sub->num_relat, fp_relat);
      //fprintf(stderr, "[%s]:  %u, generate rbwt.\n", __func__, __LINE__);
      if(sub->num_seqL > MIN_BWT_SIZE) {
#ifdef DEBUG
        if(__test_sort(sub->num_seqL , sub->sortL) > 0){
          fprintf(stderr, "[%s:%u]: seq %u wrong!\n", __func__, __LINE__, i);
          int ii;
          for(ii = 0; ii < sub->num_seqL; ++ii) {
            fprintf(stderr, "[%u]: %u!\n", ii, sub->sortL[ii]); 
          } 
          exit(1);
        }
#endif
        rbwt_t *rbwt = rbwt_bwt_build(sub->num_seqL, 16, sub->sortL, sub->seqL_buf, sub->seqR_buf);
        rbwt_bwt_dump(rbwt, fp_smbwt);
        rbwt_bwt_destroy(rbwt);
      } else{
        //log_seqs(fm_idx, sidx[0], sidx[0]+sidx[1], exti);
      }
      if(sub->num_seqR > MIN_BWT_SIZE) {
#ifdef DEBUG
        if(__test_sort(sub->num_seqR , sub->sortR) > 0){
          fprintf(stderr, "[%s:%u]: seq %u wrong!\n", __func__, __LINE__, i);
          int ii;
          for(ii = 0; ii < sub->num_seqL; ++ii) {
            fprintf(stderr, "[%u]: %u!\n", ii, sub->sortL[ii]); 
          }
          exit(1);
        }
#endif
        rbwt_t *rbwt = rbwt_bwt_build(sub->num_seqR, 16, sub->sortR, sub->seqL_buf, sub->seqR_buf);
        rbwt_bwt_dump(rbwt, fp_smbwt);
        rbwt_bwt_destroy(rbwt);
      } else{
      
      }

      cap[0].num_seqL  = sub->num_seqL;
      cap[0].num_seqR  = sub->num_seqR;
      cap[0].num_relat = sub->num_relat;
      //fwrite(cap, sizeof(struct CapIfo), 1, fp_capidx);	
      cap[0].nxtcap += sub->num_relat;
      cap[0].relat += relat_len;	 
      sub->len_capidx++;
      if(si % 100000 == 0) fprintf(stderr, "[%s]:  %u, %u rbwt finish!.\n", __func__, __LINE__, si);
    } // End: while(idx<len_seedidx) +++++++++++++++++++++++++++
    
    fclose(fp_nxtpnt);
    fclose(fp_nxtflg);
    fclose(fp_capidx);
    fclose(fp_relat);
    fclose(fp_smbwt);
    /*   
    FILE *fp_extidx = xopen(f->extidx[exti], "w");
    uint32_t len_ext_idx = out_buf[0];
    fwrite(nxt_seedidx, sizeof(uint32_t)*2, len_ext_idx, fp_extidx);
    fclose(fp_extidx);
    */
    last_seedidx = cur_seedidx;
		cur_seedidx = nxt_seedidx;
		nxt_seedidx = last_seedidx; 
    fprintf(stderr, "[%s]:  %u, hier%u idx finish!\n", __func__, __LINE__, exti);
  }
  fp = xopen("max_relat", "wb");
  fprintf(fp, "max_relat = %u\n", max_relat);
  err_fclose(fp);
  
  SubBuf_destroy(sub);
  destroyFileName(f);
  seedidx_destroy(sIdx); 
  fbwt_fmidx_destroy(idx); 
  //comFile(f,f);
	return 0;
}

int set_config(int max_relat)
{
    FILE *fp = fopen("config_file.txt", "w");
    fprintf(fp, "max_relat = %u\n", max_relat);
    fclose(fp);
}
void comFile(struct FileName *rf, struct FileName *wf)
{
	printf("comFile() = beginning OK! ");
	FILE *fp_com;
	FILE *fp;
	int f_id;
	uint32_t num = 0;
	char  buf[2];
  int i,j;
  int fsize[NUM_FILES][MAX_EXT];
  int head[NUM_FILES][MAX_EXT];
  char *fName[NUM_FILES][MAX_EXT] ;

  for(i=0;i<MAX_EXT;i++){
    fName[0][i]	=  rf->capidx[i];
    fName[1][i]	=  rf->nxtpnt[i];
    fName[2][i]	=  rf->nxtflg[i];
    fName[3][i]	=  rf->relat[i];
    fName[4][i]	=  rf->smbwt[i];
    fName[5][i]	=  rf->extidx[i];
  }
  printf("comFile() = %s\n","461 row OK! ");
   
	for(i=0;i<MAX_EXT;i++){
		fsize[0][i] =  getFileSize(rf->capidx[i]);
		fsize[1][i] =  getFileSize(rf->nxtpnt[i]);
		fsize[2][i] =  getFileSize(rf->nxtflg[i]);
		fsize[3][i] =  getFileSize(rf->relat[i]);
		fsize[4][i] =  getFileSize(rf->smbwt[i]);		
		fsize[5][i] =  getFileSize(rf->extidx[i]);
	}

  printf("comFile() = %s\n","472 row OK! ");
	for(i=0;i<MAX_EXT;i++){
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
	while(num_ext<MAX_EXT){
		if((fp_com=fopen(wf->comfile[num_ext],"wb"))==NULL){
      printf("can't open file : %s\n", wf->comfile[num_ext]);
      exit(0);
		}
    uint32_t head_buf[NUM_FILES+1] = {};
    for(i =0; i < NUM_FILES; ++i) {
      head_buf[i] = head[i][num_ext];
      head_buf[NUM_FILES] += head_buf[i];
    }
    fwrite(head_buf, sizeof(uint32_t), NUM_FILES+1, fp_com);
    f_id = 0;
    while(f_id<NUM_FILES){
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

void seedidx_open( struct FileName *f, uint32_t (*sidx)[2], uint64_t len_seedidx[],uint32_t exti,FILE *fpw[])
{
	FILE *fp;
	if((fp=fopen(f->seedidx[exti+1],"rb"))==NULL){//判断是否打开文件
	    printf("can't open file = %s\n",f->seedidx[exti+1]);
	    exit(0);
	}
	//uint32_t  numread=fread(sidx,3*sizeof(uint32_t),len_seedidx[exti+1],fp);
	uint32_t  numread=fread(sidx,2*sizeof(uint32_t),len_seedidx[exti+1],fp);
	fclose(fp);
	if((fpw[0]=fopen(f->relat[exti],"w"))==NULL){
	    printf("can't open file = %s\n", f->relat[exti] );
	    exit(0);
	}
	
	if((fpw[1]=fopen(f->smbwt[exti],"w"))==NULL){
	    printf("can't open file = %s\n", f->smbwt[exti]);
	    exit(0);
	}

	if((fpw[2]=fopen(f->nxtpnt[exti],"w"))==NULL){
	    printf("can't open file = %s\n", f->nxtpnt[exti]);
	    exit(0);
	}

	if((fpw[3]=fopen(f->nxtflg[exti],"w"))==NULL){
	    printf("can't open file = %s\n", f->nxtflg[exti]);
	    exit(0);
	}

	if((fpw[4]=fopen(f->capidx[exti],"w"))==NULL){
	    printf("can't open file = %s\n", f->capidx[exti]);
	    exit(0);
	}
/*    
	if((fpw[5] = fopen(f->extidx[exti],"w"))== NULL){
	    printf("can't open file = %s\n",f->extidx[exti]);
	    exit(0);
	}
*/
	return;
} // End: void Openfiles( ) ++++++++++++++++++++++++++++++++

//*************************************************
uint32_t getlen(uint64_t len_seedidx[],struct FileName *f){
    int i,m_i;
    uint64_t max;
    m_i=0;
    max = 0 ;
    for(i=0;i<MAX_EXT+1;i++){
        //len_seedidx[i] = getFileSize(f->seedidx[i])/12;
    	len_seedidx[i] = getFileSize(f->seedidx[i])/8;
       if(max < len_seedidx[i]){
       	 max = len_seedidx[i];
       	 m_i = i;
       } 	
    }
	return m_i;
}
//*************************************************

/*
void gen_jmpidx(uint32_t cur_seedidx[][3], uint32_t* jmp_idx, 
	uint8_t *mod_idx, uint32_t *cap_pos){

	printf("gen_jmpidx() = %s\n","OK");

	return;
}
*/
int gen_jmpidx0(int n, uint32_t cur_seedidx[][2], uint32_t* jump_idx, uint8_t *mod_idx, uint32_t *cap_pos){

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

int gen_jmpidx(int n_seeds, uint32_t seedidx[][2], uint32_t* jump_idx, uint8_t *mod_idx, uint32_t *out_buf, int exti, idx_t *fm_idx)
{
    int max_rep = 0;
    //uint32_t row_jump = 0;
    uint32_t last_idx_q = 0, idx_q = 0;
    jump_idx[0] = 0;
    //uint32_t jump_idx[LEN_NEXT_IDX+1]; 
    uint32_t row, i;
    uint32_t num, bg, pos, k, l, num0 = 0;
    uint8_t seq[LEN_SEED];
    int row_j = 0;
    for(row=0; row < n_seeds; ++row){
        if(seedidx[row][1] > max_rep) max_rep = seedidx[row][1];
        bg = seedidx[row][0];
        num = seedidx[row][1];
        //__check_seedidx(bg, num, 20 + 32* exti,fm_idx);
        if(num <= IS_SMLSIZ){ continue;}  
        if(exti == 0) {
            if(num <= IS_SMLSIZ * NUM_SCALE){ continue;}  
            seedidx[row_j][0] = seedidx[row][0];
            seedidx[row_j][1] = seedidx[row][1];
            row_j++; 
        } else if(exti > 0) {
            pos = bwt_sa(fm_idx->bwt, bg) + 16*exti;
            for(i = 0; i < LEN_SEED; ++i) seq[i] = __get_pac(fm_idx->pac, pos+i); 
            k = 0, l = fm_idx->bwt->seq_len;
            num0 = bwt_match_exact_alt(fm_idx->bwt, LEN_SEED, seq, &k, &l);
            if(num0 < num) {
                fprintf(stderr, "[Error]: num0 = %u, num = %u, %s, %u\n", num0, num, __func__, __LINE__);
                exit(1);  
            }
            if(num0 <= IS_SMLSIZ * NUM_SCALE) { continue;}
            seedidx[row_j][0] = seedidx[row][0];
            seedidx[row_j][1] = seedidx[row][1];
            row_j++;
        }  
        
        idx_q = bg/256;
        mod_idx[row_j-1] = bg%256;
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
    return max_rep;
}
int __check_seedidx(uint32_t bg, uint32_t num, int len, idx_t *idx)
{
  uint32_t i, k, l;
  uint32_t pos = bwt_sa(idx->bwt, bg);
  uint8_t *seq = calloc(len, sizeof(uint8_t));
  for(i = 0; i < len; ++i) seq[i] = __get_pac(idx->pac, pos+i);
  k =0; l = idx->bwt->seq_len;
  uint32_t num0 = bwt_match_exact_alt(idx->bwt, len, seq, &k, &l);

  if(num0 != num) {
    fprintf(stderr, "[%s:%u]: bg = %u, num = %u, num0 = %u\n", __func__, __LINE__,  bg, num, num0);
  }
  free(seq);
  return num0 == num;

}
//++++++++++++++++++++++++++++++++++++++++++++++++++
/*  */
#define IS_SEQ 0
int hier_cut_ext_seq(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, uint32_t (*seq_buf)[2])
{
  bwtint_t i, j, ii, pos, st, ed;
  //fprintf(stderr, "\n"); 
  /* fetch 16 mer preseq in [k, l) and rm duplicate preseq */
  uint32_t max_pos, min_pos;
  if(offset<0){
    min_pos = 16;
    max_pos = bwt->seq_len;
  } else{
    min_pos = 0;
    max_pos = bwt->seq_len-offset-16;
  }

  for(j = 0,i = k; i < l; ++i){
    int flag;
    uint32_t seq16;
    pos = bwt_sa(bwt, i);


    if(pos >=min_pos && pos<=max_pos){
      st = pos+offset; //ed = st+16;
      seq16 = 0;
      for(ii = st; ii < st+16; ++ii) {
        seq16 <<= 2;
        seq16 |= __get_pac(pac, ii);
      }
      flag = IS_SEQ;
    } else if(pos < min_pos){//if(st<0) 
      flag = 1;
    } else if(pos > max_pos){//if ed > bwt->seq_len
      flag = 1;
    }
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
/*  
void gen_relat(uint32_t bg, uint32_t num, uint32_t num_ext, struct SubBuf *sub, idx_t *idx)
{
  uint8_t *pac = idx->pac;
  bwt_t *bwt = idx->bwt;   

  bwtint_t i, j, k = bg, l = bg + num, bg_idx, ed_idx;
  int n = num;
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
  uint32_t flag=0, last_flag=1; 
  uint32_t row_l = 0, row_r = 0, seq; 
  for(i = 0; i < n; ++i ){
    flag = sub->seqR_buf[i][1];
    if(flag == IS_SEQ && last_flag == IS_SEQ) {
      if(sub->seqR_buf[i][0] > seq ){
        ed_idx = i-1; 
        sub->sortR[row_r] = seq;
        sub->idxR[row_r][0] = bg_idx+bg;
        sub->idxR[row_r][1] = ed_idx+bg;
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
    } else if(flag == IS_SEQ && last_flag != IS_SEQ ){
      seq = sub->seqR_buf[i][0];
      bg_idx = i;
    } else if(flag != IS_SEQ && last_flag == IS_SEQ){
      ed_idx = i-1; 
      sub->sortR[row_r] = seq;
      sub->idxR[row_r][0] = bg_idx+bg;
      sub->idxR[row_r][1] = ed_idx+bg;
      for(j = bg_idx; j <= ed_idx; ++j){
        if(sub->seqL_buf[j][1] == 0) {     
          sub->seqL_buf[row_l][0] = sub->seqL_buf[j][0]; 
          sub->seqL_buf[row_l][1] = row_r;
          row_l++; 
        }
      }
      row_r++; 
    }         
    last_flag = flag;
  } 
  if(flag==IS_SEQ){
    ed_idx = n-1; 
    sub->sortR[row_r] = seq;
    sub->idxR[row_r][0] = bg_idx+bg;
    sub->idxR[row_r][1] = ed_idx+bg;
    for(j = bg_idx; j <=ed_idx; ++j){
      if(sub->seqL_buf[j][1] == 0) {     
        sub->seqL_buf[row_l][0] = sub->seqL_buf[j][0]; 
        sub->seqL_buf[row_l][1] = row_r;
        row_l++; 
      }
    }
    row_r++; 
  } 
  sub->idxR[row_r][0] = ed_idx+bg;
  sub->idxR[row_r][1] = ed_idx+bg;

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
  int row = 0; num = 0;
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
   
  return;
}
*/
int __check_idx(idx_t *idx, bwtint_t k, bwtint_t l, int len)
{
  int i;
  uint8_t seq[1024];
  uint32_t pos = bwt_sa(idx->bwt, k);
  for(i =0; i < len; ++i) {
    seq[i] = __get_pac(idx->pac, pos + i); 
  } 
  bwtint_t k0 = 0, l0 = idx->bwt->seq_len;
  bwt_match_exact_alt(idx->bwt, len, seq, &k0, &l0);
  if(k != k0 || l != l0) {
    fprintf(stderr, "(k, l) = (%u, %u) (k0, l0) = (%u, %u)\n", k, l, k0, l0); 
    return 1;
  } else {
    return 0;
  }
}
int __log_sortL(idx_t *idx, bwtint_t k, bwtint_t l)
{
  uint32_t (*sortL)[2] = calloc(500000, 2*sizeof(uint32_t));
  int n_sortL = hier_cut_ext_seq(idx->pac, idx->bwt, k, l, -16, sortL);
  int64_t i, j,si, sseq= sortL[0][0];
  for(si = 1, i = 0; i < n_sortL; ++i){
    if(sortL[i][0] == sseq) continue;
    sortL[si][0] = sortL[i][0];
    sseq = sortL[i][0];
    si++; 
  }
  for(i = 0; i < si; ++i) {
    for(j = 0; j < 16; ++j) {
      fprintf(stderr, "%u", (sortL[i][0]>>(30-j))&3); 
    }
    fprintf(stderr, "\t%u\n", sortL[i][1]);
  }
  for(i = k; i <= l; ++i) {
    uint32_t pos = bwt_sa(idx->bwt, i) - 16;
     for(j = 0; j < 16; ++j) {
      fprintf(stderr, "%u", __get_pac(idx->pac, pos+j)); 
    }
    fprintf(stderr, "\n");
  } 
  
  free(sortL);
  return 0;
}
int __log_sortR(idx_t *idx, bwtint_t k, bwtint_t l, int exti)
{
  uint8_t seq[1024];
  int len = 20 + 32*exti;

  int64_t i, j,si, k0 =0, l0 = idx->bwt->seq_len;
  for(i = k; i < l; ++i) {
    uint32_t pos = bwt_sa(idx->bwt, i);
    
    for(j = 0; j < len+16; ++j) {
      seq[j] = __get_pac(idx->pac, pos+j);
      if(j >= len) fprintf(stderr, "%u", __get_pac(idx->pac, pos+j)); 
    }
    k0 =0, l0 = idx->bwt->seq_len;
    bwt_match_exact_alt(idx->bwt, len+16, seq, &k0, &l0);
    fprintf(stderr, "\t%u, %u\n", k0, l0);
  } 
  

  return 0;
}

void log_seq16(uint32_t seq16)
{
  int i;
  for(i= 0; i < 16; ++i){
    fprintf(stderr, "%u", ((seq16>>(30-i*2)) &3));
  }
}

void gen_relat(uint32_t bg, uint32_t num, uint32_t num_ext, struct SubBuf *sub, idx_t *idx)
{
  //log_seqs(idx, bg, bg+num, num_ext);
  uint8_t *pac = idx->pac;
  bwt_t *bwt = idx->bwt;   

  bwtint_t i, j, k = bg, l = bg + num, bg_idx, ed_idx;
  int n = num;
  int off_pos = LEN_SEED+2*16*num_ext; 
  int n_r = hier_cut_ext_seq(pac, bwt, k, l,  off_pos, sub->seqR_buf);//right seq buf
  //log_seqbuf(n_r, sub->seqR_buf);
  off_pos = -16; 
  int n_l = hier_cut_ext_seq(pac, bwt, k, l,  off_pos, sub->seqL_buf);//left seq buf
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //下面代码段完成以下几个功能：
  //1. 过滤左段和右段不合理的数据导致的扩展。
  //2. 生成右段的排序后的唯一序列。
  //3. 生成右段唯一序列对应的index区间。
  //4. 过滤后的左段序列集合中的第二列输入了相应的右段排序后序列的行号以便后面生成关联数组
  uint32_t flag=0; 
  uint32_t row_l = 0, row_r = 0, last_seq; 
  bwtint_t ri = 0, rii, li=0;
  uint32_t bg_ri, ed_ri;
  while(ri < n_r){
    //last_seq = sub->seqR_buf[ri][0]; 
    flag = sub->seqR_buf[ri][1];
    while(ri < n_r && flag != IS_SEQ) {
      ++ri; 
      flag = sub->seqR_buf[ri][1];
    }
    bg_idx = ri;
    while(ri < n_r && flag == IS_SEQ) {
      ++ri;      
      flag = sub->seqR_buf[ri][1]; 
    }
    ed_idx = ri;
    
    rii = bg_idx;
    while(rii < ed_idx){
      bg_ri = rii;
      /*  
      while(bg_ri == rii|| cur_seq == last_seq){
        last_seq = cur_seq;
        ++rii; 
        cur_seq = sub->seqR_buf[rii][0];
      }
      */
      while(rii < ed_idx && sub->seqR_buf[rii][0] == sub->seqR_buf[bg_ri][0]){
        ++rii; 
      }

      ed_ri = rii -1; 
      for(li = bg_ri; li <= ed_ri; ++li) {
        if(sub->seqL_buf[li][1] == IS_SEQ) {     
          sub->seqL_buf[row_l][0] = sub->seqL_buf[li][0]; 
          sub->seqL_buf[row_l][1] = row_r;
          row_l++;
        }
      
      }
      sub->sortR[row_r] = sub->seqR_buf[bg_ri][0];      
      sub->idxR[row_r][0] = bg + bg_ri;
      sub->idxR[row_r][1] = bg + ed_ri;
      //fprintf(stderr, "[RSEQ]: seq = %u, %u, %u\n",sub->seqR_buf[bg_ri][0], sub->idxR[row_r][0], sub->idxR[row_r][1]);
      //__log_sortL(idx, bg+bg_ri, bg+ed_ri); 
      ++row_r;
      
      //if(0) {
      if(__check_idx(idx, bg+bg_ri, bg+ed_ri, 20 + 32*num_ext+16) != 0) {
        fprintf(stderr, "[%s]:  error in line %u!\n", __func__, __LINE__); 
        __log_sortR(idx, k, l, num_ext);
        exit(1); 
      }
    }
  } 
  
  sub->sortR[row_r] = (uint32_t)-1;      
  sub->idxR[row_r][0] = ed_ri+bg;
  sub->idxR[row_r][1] = ed_ri+bg;

  uint32_t num_seqL = row_l;
  sub->num_seqR = row_r;

  uint32_t (*seqL)[2] = sub->seqL_buf;
  //sort(key= (uniq_lid, uniq_rid))
  ks_introsort(pair_t, num_seqL, (pair_t *)seqL); 
 
  
  last_seq = seqL[0][0]; 
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



  /*  
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
  int row = 0; num = 0;
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
  */
  //log_sub(sub); 
  return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++
struct SeedIdx* seedidx_init(uint32_t max_seeds, uint32_t ref_len)
{
	struct SeedIdx *sIdx = (struct SeedIdx*)malloc(sizeof(struct SeedIdx));
  uint32_t (*cur_seedidx)[2];
	uint32_t (*nxt_seedidx)[2];
	uint32_t *jmp_idx;
	uint32_t *cap_pos;
	uint8_t  *mod_idx; 


	if(NULL == (cur_seedidx = (uint32_t(*)[2])calloc(max_seeds, 2*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}

  if(NULL == (nxt_seedidx = (uint32_t(*)[2])calloc(max_seeds, 2*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}

	sIdx->cur_seedidx = cur_seedidx;
	sIdx->nxt_seedidx = nxt_seedidx;
	sIdx->n_jmp = (ref_len+511)/256;
	sIdx->n_mod = (ref_len+511)/256;
	sIdx->n_cap = max_seeds;
  sIdx->jmp_idx     = (uint32_t*)calloc(((ref_len+511)/256), sizeof(uint32_t));
	sIdx->mod_idx     = (uint32_t*)calloc(((ref_len+511)/256), sizeof(uint32_t));
	sIdx->cap_pos     = (uint32_t*)calloc(max_seeds, sizeof(uint32_t));
	
  
  return sIdx;
}
void seedidx_destroy(struct SeedIdx *sIdx)
{
  free(sIdx->jmp_idx);
  free(sIdx->mod_idx);
  free(sIdx->cap_pos);  
  free(sIdx->cur_seedidx);
  free(sIdx->nxt_seedidx);
  free(sIdx);  
  return;
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//#define MAX_BLCK_SIZE 500000
#define MAX_SEED_SIZE 500000 
#define ALL_BUFS_SIZE 13  // All bufs size in SubBuf



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct SubBuf *SubBuf_init(int max_buf_size)
{
	
  struct SubBuf *sub= (struct SubBuf*)malloc(sizeof(struct SubBuf));
 	sub->max_buf_size = max_buf_size; 
  uint32_t *buf;
  if(NULL == (buf = (uint32_t*)calloc(max_buf_size*ALL_BUFS_SIZE, sizeof(uint32_t)))){
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
  sub->relat_sai   = (uint32_t(*)[2])(buf + max_buf_size*2); //注：extIdx[][2]的数据存放在;seqR_buf[][2]的物理空间。
  sub->idxR     = calloc(max_buf_size, 2*sizeof(uint32_t)); //注：idxR[][2]的数据存放在;seqL_buf[][2]的物理空间。
  sub->pos_buf = buf + max_buf_size*10 ;
  sub->nxt_idx = buf + max_buf_size*11;
  sub->nxt_flg = (uint8_t*)(buf + max_buf_size*12);
 
 
  //sub->R_relat = buf + max_buf_size*13;
  
  /*
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

  */
	sub->num_seqL = 0;
	sub->num_seqR = 0;
	sub->num_relat= 0;	

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
  
  return sub;
}

void SubBuf_destroy(struct SubBuf *sub)
{
  free(sub->seqR_buf);
  free(sub->idxR);
  free(sub);
  return;
}


/*  
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
*/
//uint32_t get_jmpidx(int n_jmp_idx, uint32_t *jmp_idx, int n_mod_idx, uint8_t *mod_idx, uint32_t bgn, uint32_t end)

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  sai2local
 *  Description: convert suffix array interval [bgn, end] to local index id 
 * =====================================================================================
 */
uint32_t sai_to_local_id(int n_jmp, uint32_t *jmp_idx, int n_mod, uint8_t *mod_idx, uint32_t bgn, uint32_t end)
{
  bwtint_t mod_bg, mod_ed, sai_r, i; 
  if(bgn/256 > n_jmp) {
    fprintf(stderr, "[%s:%u]: bgn/256(%u) > n_jmp(%u)!\n", __func__, __LINE__, bgn/256, n_jmp);
    exit(1); 
  }
  mod_bg = jmp_idx[bgn/256];
  mod_ed = jmp_idx[(bgn/256)+1];
  sai_r = bgn%256;
  if(mod_ed > n_mod) {
    fprintf(stderr, "[%s:%u]: mod_ed(%u) > n_mod(%u)!\n", __func__, __LINE__, mod_ed, n_mod);
    exit(1); 
 
  
  }
  for(i = mod_bg; i < mod_ed; i++){
    if(sai_r == mod_idx[i]){
      break;
    }
  }
  if(i == mod_ed){ return -1; }
  return i;
}
uint32_t aln_mam(int l_seq, uint8_t *seq, int s_bg, int s_ed, bwtint_t *bg, bwtint_t *ed, idx_t *idx, int max_loc)
{
  int64_t i, j;
  uint32_t k = *bg;
  uint32_t l = *ed;
  if(k + max_loc > l) {
    return 1; 
  }
  uint32_t lid = sai_to_local_id(idx->n_jmp, idx->jmp_idx, idx->n_mod, idx->mod_idx, k, l);
  int exti = 0;
  hidx_t *lidx = idx->hier_idx[exti]+lid;
  while(k + max_loc < l && s_bg > 16 && s_ed+16 <l_seq ) { 
    uint32_t k0 = 0, l0 = lidx->bwt0->n_seqs, n0;
    n0 = rbwt_exact_match(lidx->bwt0, lidx->bwt0->n_rot, idx->cnt_table, 16, seq+s_bg-16, &k0, &l0);
    uint32_t k1 = 0, l1 = lidx->bwt0->n_seqs, n1;
    n1 = rbwt_exact_match(lidx->bwt1, lidx->bwt1->n_rot, idx->cnt_table, 16, seq+s_ed, &k1, &l1);
    if(n0 !=1 || n1 != 1) {
      fprintf(stderr, "[%s:%lu]: n0 = %d, n1 = %d\n", __func__, __LINE__, n0, n1); 
      exit(1); 
    }
    for(i = k0; i < l0; ++i){
      uint32_t bg, ed;
      bg = lidx->L2rel[i];
      ed = lidx->L2rel[i+1];
      for(j = bg; j < ed; ++j) {
        if(lidx->relat[j] == k1) {
          break;
        } 
      }
      if(j == ed) {
        lid = (uint32_t)-1; 
      } else {
        lid = lidx->nxt_idx[j]; 
        lidx = idx->hier_idx[exti]+lid;
        k = lidx->bgn;
        l = k + lidx->num;
        s_bg -= 16;
        s_ed += 16;
      } 
    }
    exti++; 
  }
}
uint32_t aln_mam_alt(idx_t *idx ,int l_seq, uint8_t *seq, int *s_k, int *s_l, bwtint_t *bg, bwtint_t *ed,  int max_loc)
{
  int64_t i, j;
  uint32_t k = *bg;
  uint32_t l = *ed;
  uint32_t s_bg = *s_k;
  uint32_t s_ed = *s_l;
  if(k + max_loc > l) {
    return l - k; 
  }
  uint32_t lid = sai_to_local_id(idx->n_jmp, idx->jmp_idx, idx->n_mod, idx->mod_idx, k, l);
  int exti = 0;
  int nxt_flg = 0;
  while(s_bg >= 16 && s_ed+16 <l_seq ) { 
    hidx_t *lidx = idx->hier_idx[exti]+lid;
    k = lidx->bgn;
    l = k + lidx->num-1;
    if(exti >= 4 || k + max_loc > l) { break; }
#ifdef DEBUG 
        uint32_t kk = 0, ll = idx->bwt->seq_len;
        bwt_match_exact_alt(idx->bwt, s_ed - s_bg, seq+s_bg, &kk, &ll);
        if(k != kk || l != ll) {
          fprintf(stderr, "[%s:%u]: exti = %u, lid = %d, (%u, %u) != (%u, %u)\n", __func__, __LINE__, exti, lid, k, l, kk, ll); 
          exit(1); 
        }
#endif
    uint32_t k0 = 0, l0 = lidx->bwt0->n_seqs, n0 = 0, rot0;
    n0 = rbwt_appr_match(lidx->bwt0, lidx->bwt0->n_rot, idx->cnt_table, 16, seq+s_bg-16, &rot0, &k0, &l0);
    if(n0 <= 0) goto end;
    uint32_t k1 = 0, l1 = lidx->bwt1->n_seqs, n1 = 0, rot1;
    n1 = rbwt_appr_match(lidx->bwt1, lidx->bwt1->n_rot, idx->cnt_table, 16, seq+s_ed, &rot1, &k1, &l1);
    if(n1 <= 0) goto end;
    int ri, li; 
    for(k = k1; k < l1; ++k) {
      ri = rbwt_seq_rank(lidx->bwt1, lidx->bwt1->n_rot, idx->cnt_table, rot1, k);
      for(i = k0; i < l0; ++i){
        li = rbwt_seq_rank(lidx->bwt0, lidx->bwt0->n_rot, idx->cnt_table, rot0, i);
        uint32_t bg, ed;
        bg = lidx->L2rel[li];
        ed = lidx->L2rel[li+1];
        //pairing lseq and rseq [todo: use binary search]
        for(j = bg; j < ed; ++j) {
          if(lidx->relat[j] == ri) {
            break;
          } 
        }
        if(j == ed) { // pairing fail
          break;
          lid = (uint32_t)-1;
          fprintf(stderr, "[%s:%u]: error!!!\n", __func__, __LINE__);
          exit(1); 
        }
        lid = lidx->nxt_idx[j];
        nxt_flg = lidx->nxt_flg[j];
        s_bg -= 16;
        s_ed += 16;
        
        if(nxt_flg <= IS_SMLSIZ) {
          k = lid;
          l = lid + nxt_flg-1;
          goto end;
        } 
      }
    }
    exti++; 
  }
end:
  *bg = k;
  *ed = l;
  *s_k = s_bg;
  *s_l = s_ed;
  return l - k;  
}


static inline int __cmp(const void *x, const void *y)
{
  uint32_t *a = (uint32_t *)x;
  uint32_t *b = (uint32_t *)y;
  if(*a < *b) return -1;
  else if(*a == *b) return 0;
  else return 1;
}
uint32_t aln_mem_alt(idx_t *idx ,int l_seq, uint8_t *seq, int *s_k, int *s_l, bwtint_t *bg, bwtint_t *ed,  int max_loc)
{
  int64_t i, j;
  uint32_t k = *bg;
  uint32_t l = *ed;
  uint32_t s_bg = *s_k;
  uint32_t s_ed = *s_l;
  if(k + max_loc > l) {
    return l - k; 
  }
  uint32_t lid = sai_to_local_id(idx->n_jmp, idx->jmp_idx, idx->n_mod, idx->mod_idx, k, l);
  int exti = 0;
  int nxt_flg = 0;
  //while(s_bg >= 16 && s_ed+16 <l_seq ) { 
  while(exti <= 4) { 
    hidx_t *lidx = idx->hier_idx[exti]+lid;
    k = lidx->bgn;
    l = k + lidx->num-1;
    if(exti >= 4 || k + max_loc > l|| s_bg<16||s_ed+16>l_seq) { goto end;}
#ifdef DEBUG 
        uint32_t kk = 0, ll = idx->bwt->seq_len;
        bwt_match_exact_alt(idx->bwt, s_ed - s_bg, seq+s_bg, &kk, &ll);
        if(k != kk || l != ll) {
          fprintf(stderr, "[%s:%u]: exti = %u, lid = %d, (%u, %u) != (%u, %u)\n", __func__, __LINE__, exti, lid, k, l, kk, ll); 
          exit(1); 
        }
#endif


    uint32_t k0 = 0, l0 = lidx->bwt0->n_seqs, n0 = 0;
    n0 = rbwt_exact_match(lidx->bwt0, lidx->bwt0->n_rot, idx->cnt_table, 16, seq+s_bg-16, &k0, &l0);
    if(n0 != 1) goto end;
    uint32_t k1 = 0, l1 = lidx->bwt1->n_seqs, n1 = 0;
    n1 = rbwt_exact_match(lidx->bwt1, lidx->bwt1->n_rot, idx->cnt_table, 16, seq+s_ed, &k1, &l1);
    if(n1 != 1) goto end;
    
    for(i = k0; i < l0; ++i){
      uint32_t bg, ed;
      bg = lidx->L2rel[i];
      ed = lidx->L2rel[i+1];
      //pairing lseq and rseq [todo: use binary search]
      /*  
      for(j = bg; j < ed; ++j) {
        if(lidx->relat[j] == k1) {
          break;
        } 
      }
      */
  
      uint32_t *is_find = bsearch(&k1, lidx->relat+bg, ed - bg, sizeof(uint32_t), __cmp);
      /*  
      if(is_find == NULL) {
        if(j != ed) {
          fprintf(stderr, "[%s:%u]: key = %u, j = %u, ed = %u\n", __func__, __LINE__, is_find, j, ed);
          for(j = bg; j < ed; ++j) {
            fprintf(stderr, "[%s:%u]: j = %u, %u, k1 = %u\n", __func__, __LINE__, j, lidx->relat[j], k1);
          }
          exit(1);
        } 
      } else if(is_find != lidx->relat+j) {
        fprintf(stderr, "[%s:%u]: key = %u, j = %u\n", __func__, __LINE__, (is_find-lidx->relat- bg)/4, j);
        exit(1);
      }
      */
      //if(j == ed) { // pairing fail
      if(is_find == NULL) { // pairing fail
        goto end;
        lid = (uint32_t)-1;
        fprintf(stderr, "[%s:%u]: error!!!\n", __func__, __LINE__);
        exit(1); 
      }
      j = (is_find-lidx->relat)/sizeof(uint32_t);
      lid = lidx->nxt_idx[j];
      nxt_flg = lidx->nxt_flg[j];
      s_bg -= 16;
      s_ed += 16;
      
      if(nxt_flg <= IS_SMLSIZ) {
        k = lid;
        l = lid + nxt_flg-1;
        goto end;
      } 
    }
    exti++; 
  }
end:
  *bg = k;
  *ed = l;
  *s_k = s_bg;
  *s_l = s_ed;
  return l - k;  
}
uint32_t aln_mem(idx_t *idx ,int l_seq, uint8_t *seq, int *s_k, int *s_l, bwtint_t *bg, bwtint_t *ed,  int max_loc)
{
  int64_t i, j;
  uint32_t k = *bg;
  uint32_t l = *ed;
  uint32_t s_bg = *s_k;
  uint32_t s_ed = *s_l;
  if(k + max_loc > l) {
    return l - k; 
  }
  uint32_t lid = sai_to_local_id(idx->n_jmp, idx->jmp_idx, idx->n_mod, idx->mod_idx, k, l);
  int exti = 0;
  int nxt_flg = 0;
  while(s_bg >= 16 && s_ed+16 <l_seq ) { 
    hidx_t *lidx = idx->hier_idx[exti]+lid;
    k = lidx->bgn;
    l = k + lidx->num-1;
    if(exti >= 4 || k + max_loc > l) { break; }
#ifdef DEBUG 
        uint32_t kk = 0, ll = idx->bwt->seq_len;
        bwt_match_exact_alt(idx->bwt, s_ed - s_bg, seq+s_bg, &kk, &ll);
        if(k != kk || l != ll) {
          fprintf(stderr, "[%s:%u]: exti = %u, lid = %d, (%u, %u) != (%u, %u)\n", __func__, __LINE__, exti, lid, k, l, kk, ll); 
          exit(1); 
        }
#endif


    uint32_t k0 = 0, l0 = lidx->bwt0->n_seqs, n0 = 0;
    n0 = rbwt_exact_match(lidx->bwt0, lidx->bwt0->n_rot, idx->cnt_table, 16, seq+s_bg-16, &k0, &l0);
    uint32_t k1 = 0, l1 = lidx->bwt1->n_seqs, n1 = 0;
    n1 = rbwt_exact_match(lidx->bwt1, lidx->bwt1->n_rot, idx->cnt_table, 16, seq+s_ed, &k1, &l1);
    if(n0 !=1 || n1 != 1) {
      if(n0 != 1) {
        uint8_t *seq_buf = calloc(lidx->bwt0->n_seqs, 16*sizeof(uint8_t));
        rbwt_bwt_to_seqs(lidx->bwt0, lidx->bwt0->n_rot, seq_buf, idx->cnt_table);
        int ii, jj;
        for(jj = 0; jj < 16; ++jj) {
          fprintf(stderr, "%u\t", seq[s_bg-16+jj]);
        }

        fprintf(stderr, "\n");
        for(ii = 0; ii < lidx->bwt0->n_seqs; ++ii) {
          for(jj = 0; jj < 16; ++jj) {
            fprintf(stderr, "%u\t", seq_buf[ii*16+jj]);
          }
          
          fprintf(stderr, "\n");
        }
        __log_sortL(idx, k, l);
        
        free(seq_buf);
      } else if(n1 != 1) {
        uint8_t *seq_buf = calloc(lidx->bwt1->n_seqs, 16*sizeof(uint8_t));
        rbwt_bwt_to_seqs(lidx->bwt1, lidx->bwt1->n_rot, seq_buf, idx->cnt_table);
        int ii, jj;
        for(jj = 0; jj < 16; ++jj) {
          fprintf(stderr, "%u\t", seq[s_ed+jj]);
        }

        fprintf(stderr, "\n");
        for(ii = 0; ii < lidx->bwt1->n_seqs; ++ii) {
          for(jj = 0; jj < 16; ++jj) {
            fprintf(stderr, "%u\t", seq_buf[ii*16+jj]);
          }
          
          fprintf(stderr, "\n");
        }
         
        __log_sortR(idx, k, l, exti);
      
        free(seq_buf);
      }
      
      
      fprintf(stderr, "[%s:%lu]: n0 = %d, n1 = %d\n", __func__, __LINE__, n0, n1); 
      exit(1); 
    }
    for(i = k0; i < l0; ++i){
      uint32_t bg, ed;
      bg = lidx->L2rel[i];
      ed = lidx->L2rel[i+1];
      //pairing lseq and rseq [todo: use binary search]
      for(j = bg; j < ed; ++j) {
        if(lidx->relat[j] == k1) {
          break;
        } 
      }
      if(j == ed) {
        lid = (uint32_t)-1;
        fprintf(stderr, "[%s:%u]: error!!!\n", __func__, __LINE__);
        exit(1); 
      }
      lid = lidx->nxt_idx[j];
      nxt_flg = lidx->nxt_flg[j];
      s_bg -= 16;
      s_ed += 16;
      
      if(nxt_flg <= IS_SMLSIZ) {
        k = lid;
        l = lid + nxt_flg-1;
        goto end;
      } 
    }
    exti++; 
  }
end:
  *bg = k;
  *ed = l;
  *s_k = s_bg;
  *s_l = s_ed;
  return l - k;  
}

void gen_nxt_cap(struct SubBuf *subBuf,int n_jmp, uint32_t *jmp_idx, int n_mod, uint8_t *mod_idx,uint32_t *cap_pos, idx_t *fm_idx)
{
  uint32_t  i, j,	k, bgn, end, num, i_r, seq, pos, row, local_id;
  uint32_t buf_idx[5];
  struct SubBuf *sub = subBuf;
  sub->num_pos = 0;
  for(i=0; i<sub->num_seqL; i++){
    bgn = sub->L2rel[i];
    end = sub->L2rel[i+1];

    for(j=bgn; j<end; j++){
      sub->seqL_idxR[j][0] = sub->sortL[i];
      k = sub->relat[j];
      sub->seqL_idxR[j][1] = sub->idxR[k][0];
      sub->seqL_idxR[j][2] = sub->idxR[k][1];

    }
  }
  /*  
  for(i=0; i<sub->num_seqR; i++){
    bgn = sub->R2rel[i];
    end = sub->R2rel[i+1];
    for(row=bgn; row<end; row++){
      j = sub->relat[row];
      sub->seqL_idxR[j][1] = sub->idxR[i][0];
      sub->seqL_idxR[j][2] = sub->idxR[i][1];
    }
  }
  */
  //fprintf(stderr, "[%s]:  %u, num_seqL = %u, num_seqR = %d, num_relat = %d\n",__func__,  __LINE__, sub->num_seqL, sub->num_seqR, sub->num_relat);	
  for(i=0;i<sub->num_relat;i++){
    seq = sub->seqL_idxR[i][0]; // seq of lseq
    buf_idx[0] = sub->seqL_idxR[i][1]; // begin sai of mseq + rseq
    buf_idx[1] = sub->seqL_idxR[i][2]; // end sai of mseq + rseq

    //fprintf(stderr, "[%s]:  i = %u, seq = %x, seqL_idxR[1] = %u, seqL_idxR[2] = %u\n", __func__, i, seq, sub->seqL_idxR[i][1], sub->seqL_idxR[i][2]);	
    uint8_t seq_nt[256];
/*  
    for(j=0;j<16;j++){
        ch[j] = seq&3;
        seq >>= 2;
    }
*/
    for(j=0;j<16;j++){
      seq_nt[15-j] = seq&3;
      seq >>= 2;
    }

    int is_aln = bwt_match_exact_alt(fm_idx->bwt, 16, seq_nt, &buf_idx[0], &buf_idx[1]);
    if(is_aln <1 ){ 
      fprintf(stderr, "[%u, error]:is_aln = %d, k = %u, l = %u\n", __LINE__, is_aln, buf_idx[0], buf_idx[1]);
      uint8_t sq[256], sq1[256];
      uint32_t pos = bwt_sa(fm_idx->bwt, buf_idx[0])-16;
      uint32_t pos1 = bwt_sa(fm_idx->bwt, buf_idx[1])-16;
      for(j = 0; j < 16; ++j) sq[j] = __get_pac(fm_idx->pac, pos+j);
      for(j = 0; j < 16; ++j) sq1[j] = __get_pac(fm_idx->pac, pos1+j);
      for(j = 0; j < 16; ++j) fprintf(stderr, "%u\t%u\t%u\n", seq_nt[j], sq[j], sq1[j]); 
      exit(1);

    }
    sub->relat_sai[i][0]= buf_idx[0]; //返回的bgnIdx
    sub->relat_sai[i][1]= buf_idx[1]; //返回的endIdx
#ifdef DEBUG    
    if(__check_idx(fm_idx, sub->relat_sai[i][0], sub->relat_sai[i][1], (sub->num_ext+1)*32+LEN_SEED) != 0) {
      fprintf(stderr, "i = %u, relat_sai = (%u, %u), pos = %u\n", i, sub->relat_sai[i][0], sub->relat_sai[i][1], bwt_sa(fm_idx->bwt, sub->relat_sai[i][0]));	
      exit(1);
    
    }
#endif
    
    
    
  /*   
  //test code测试代码
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++
   
    uint32_t tmp[3];
    uint32_t pos = bwt_sa(fm_idx->bwt, buf_idx[0]);
    uint32_t __l =  LEN_SEED+(sub->num_ext+1)*32;
    if(pos >fm_idx->bwt->seq_len-__l) {
      fprintf(stderr, "pos out range!\n");
      exit(1);
    }
    bwtint_t __i;
    for(__i= 0; __i < __l; ++__i){
      seq_nt[__i] = __get_pac(fm_idx->pac, pos+__i); 
    } 
     tmp[0] = 0; tmp[1] = fm_idx->bwt->seq_len; 
    tmp[2] = bwt_match_exact_alt(fm_idx->bwt, __l, seq_nt, &tmp[0], &tmp[1]);

     if(tmp[0] != buf_idx[0] || tmp[1] != buf_idx[1]){

      fprintf(stderr, "[Test %u]: buf=(%u, %u, %u), tmp = (%u, %u, %u)\n", __l, buf_idx[0], buf_idx[1], is_aln, tmp[0], tmp[1], tmp[2]); 
      uint32_t __k0 = 0, __l0 = fm_idx->bwt->seq_len; 
      tmp[2] = bwt_match_exact_alt(fm_idx->bwt, __l-32, seq_nt+16, &__k0, &__l0);
      fprintf(stderr, "[Test]: %u, %u\n", __k0, __l0); 
      for(j = __k0; j <= __l0; ++j) {
        pos = bwt_sa(fm_idx->bwt, j);
        //fprintf(stderr, "pos[%u] out of range, %u\n", j, pos); 
        if(pos <= 16 || pos >= fm_idx->bwt->seq_len-16){
          fprintf(stderr, "pos out of range, %u\n", pos); 
        } 
      }
      exit(1); 
    }
   */ 
  }
  sub->num_pos = 0;
  for(i=0;i<sub->num_relat;i++){
    bgn = sub->relat_sai[i][0];
    end = sub->relat_sai[i][1];
    num = end+1 - bgn;
    /*  
    if(num==1){
      pos = bwt_sa(fm_idx->bwt, bgn);
      sub->nxt_idx[i] = pos;
      sub->nxt_flg[i] = num;
   
    } else 
    */
    if(num <= IS_SMLSIZ){
      sub->nxt_idx[i] = bgn;           
      sub->nxt_flg[i] = num;
  
    } else{
      if(sub->num_ext == 4){
        sub->nxt_idx[i] = bgn;           
      } else { 
        local_id = sai_to_local_id(n_jmp, jmp_idx, n_mod, mod_idx, bgn, end);
        /*  
        bgn_q = jmp_idx[bgn/256];
        end_q = jmp_idx[(bgn/256)+1];
        i_r = bgn%256;
        for(local_id=bgn_q; local_id<end_q; local_id++){
          if(i_r == mod_idx[local_id]){
            break;
          }
        }
        */
        if(local_id == end){
            fprintf(stderr,"[repeat_num = %u]: no local_id!\n", sub->num_ext);
            fprintf(stderr,"bg_idx, ed_idx= %u, %u!\n", bgn, end);
            fprintf(stderr,"bgn_q, end_q= %u, %u!\n", jmp_idx[bgn/256], jmp_idx[bgn/256+1]);
            exit(1); 
        }
        sub->nxt_idx[i] = local_id; //指向下一级CapIdx[]的行号
      }
      if(num < 255){
        sub->nxt_flg[i] = (uint8_t)num ; //指向CapIdx[]的类型
      }else{
        sub->nxt_flg[i] = 0xFF ;
      }
    }
  } // End : for(i=0;i<len_relat;i++)
	return;
}

