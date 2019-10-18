
//种子扩展多级索引层次数据结构
//种子序列双向扩展索引的数据结构如下：
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <dirent.h>
#include <sys/time.h>
#include <time.h>
#include <sys/stat.h>
#include <stdint.h>
#include "setFileName.h"
#include "index.h"


#include <zlib.h>
#include "utils.h"
#include "khash.h"
#include "index.h"
#include "kvec.h"
#include "ksort.h"
#include "malloc_wrap.h"
#define USE_MALLOC_WRAPPERS

#define MAX_NAME 128
KHASH_MAP_INIT_STR(str, int)
typedef kvec_t(uint32_t) vec_uint_t;
typedef kvec_t(ext_t) vec_ext_t;
#define SWAP(type_t, a, b) do{type_t x=(a); (a)=(b); (b)=x;} while(0)
#define __set_bwt(bwt, l, c) ((bwt)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define __get_bwt(bwt, l) ((bwt)[(l)>>2]>>((~(l)&3)<<1)&3)


//#define __sort_lt_ext(a, b)((a).seq == (a).seq?(a).idx<(b).idx:(a).seq < (b).seq ) 
#define __sort_lt_ext(a, b)((a).seq < (b).seq||((a).seq==(b).seq&&(a).idx < (b).idx )) 

KSORT_INIT_GENERIC(uint32_t)
KSORT_INIT(ext, ext_t, __sort_lt_ext)



#define  LEN_SEED  20 
#define  OFF_SEED  10 
#define  LEN_EXT   16 
#define  LEN_READ  150
#define  NUM_EXT  (((LEN_READ- LEN_SEED )/2)+OFF_SEED+ LEN_EXT-1)/LEN_EXT  
#define  MAX_SEED_NUM 600000
#define  LEN_FILE_NAME 100 
#define  IS_SMLSIZ 8
#define  LEN_EXT_IDX  
//LEN_EXT_IDX 是所有ExtIdx[][3]数组长度中最大的，可以用文件大小来计算。

//--------------------------------------------------------------
//扩展种子Index区间生成过程：
//主BWT索引已经创建，包括:
//Bwt后缀前置序列数组SA[]，
//累加数组Rank[]，
//Index到Pos的映射数组ItoPos[]数组，
//参考基因组序列RefSeq[]数组，
//-----------------------------------------------------------------

typedef struct{
	uint32_t relat;	
	uint8_t  smbwt;
	//uint32_t smpos;
    uint32_t nxtcap;	

	uint32_t num_seqL ;  // 当前块中seqL数据个数；
	uint32_t num_seqR ;  // 当前块中seqR数据个数；
	uint32_t num_relat;  // 当前块中relat数据个数；	
	uint32_t num_pos ; 	 // 当前块中pos数据个数；

} CapInfo;

typedef struct{  

	uint32_t (*seqL_buf)[2];  //截取的数据缓冲区
	uint32_t (*seqR_buf)[2];
//+++++++++++++++++++++++++++++++++++++++++
	uint32_t *sortL; //排序后的唯一序列
	uint32_t *sortR;
	uint32_t *relat ;	
	uint32_t *L2rel ;
	uint32_t (*idxR)[2]; // sortR[]的每一个序列的(bgnIdx,ednIdx)

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
	uint32_t len_smpos  ;  // 当前块前为止smpos数组中的全局pos数据长度；

	CapInfo cap[1];
	//uint32_t off_pos ; 
    //uint32_t blck_id;

} SubBuf;  


void setFileName(char *in_fname, FileName *out_fname,int flg );
uint32_t getlen(uint32_t LenExtIdx[],FileName *f);
void InitSubBuf(SubBuf *sub_buf);
void CrtJmpMod(int n, uint32_t cur_seedidx[][3], uint32_t* jmp_idx, uint8_t *mod_idx, uint32_t *cap_pos);

void OpenSeedIdxfiles(FileName *f, uint32_t (*sidx)[3], uint32_t LenExtIdx[],uint32_t repNum,FILE *fpw[6]);
void relatNxtCapIdx(SubBuf *subBuf,uint32_t *jmp_idx, uint8_t *jmp_mod,uint32_t *cap_pos,uint32_t *ItoPos);
void comFile(FileName *rf, FileName *wf);
long getFileSize(char* file);
void getBlckData(uint32_t sidx[2],uint32_t num_ext, 
        	SubBuf *sub_buf, idx_t *idx/*pubParm pParm*/ );

void buldSmBwt(SubBuf *sub,int flg);
int build_smbwt_idx(idx_t *idx){

	//+++++++++++++++++++++++++++++++++++++++++++
	//初始化
	char in_f[1024];
	FileName  f;
	strcpy(in_f,"seedidx");
	setFileName(in_f,&f,0);
	
	uint32_t LenExtIdx[NUM_EXT];
	uint32_t MaxLen; 
	uint32_t i = getlen(LenExtIdx,&f);  
	//LenExtIdx是所有ExtIdx[][3]数组长度中最大的，
	//可以用文件大小来计算,返回最大值序号。
	MaxLen = LenExtIdx[i];

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//内存初始化------------------------------------------------------
	uint32_t (*cur_seedidx)[3], (*ord_seedidx)[3], (*nxt_seedidx)[3];
	uint32_t *jmp_idx, *cap_pos;
	uint8_t  *mod_idx; 

	if(NULL == (cur_seedidx = (uint32_t(*)[3])malloc(MaxLen*3*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}

	if(NULL == (nxt_seedidx = (uint32_t(*)[3])malloc(MaxLen*3*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}

	if(NULL == (jmp_idx = (uint32_t*)malloc(((MaxLen+255)/256)*sizeof(uint32_t)))){
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

	//==============================================================
	//
	//
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//数组初始化----------------------------------------------------
	SubBuf sub;
	InitSubBuf(&sub);

	CapInfo  *cap;
	cap = sub.cap;

	FILE *fp;
	if((fp=fopen64(f.seedidx[0],"r"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}
	uint32_t numread=fread(cur_seedidx,3*sizeof(uint32_t),LenExtIdx[0],fp);
	//printf("Number of items read=%d\n",numread);
	fclose(fp);
	//================================================================
	//计算cur_seedidx的jmp_idx,mod_idx, cap_pos数组，并输出到FlgIdxFile_1
	CrtJmpMod(LenExtIdx[0], cur_seedidx, jmp_idx, mod_idx,cap_pos);
	//FileName  jf;
	strcpy(in_f,"jmpmod");
	setFileName(in_f,&f,1);

	if((fp=fopen64(*(f.jmpmod),"w"))==NULL){//判断是否打开文件
	    printf("can't open file\n");
	    exit(0);
	}
	fwrite(cap_pos,sizeof(uint32_t),LenExtIdx[0],fp);  
	fwrite(jmp_idx,sizeof(uint32_t),(LenExtIdx[0]+512-1)/256,fp);  
	fwrite(mod_idx,sizeof(uint32_t),(LenExtIdx[0]+3)/4,fp);  
	fclose(fp); 

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//开始循环处理----------------------------------------------------
	uint8_t repNum = 0 ;
	strcpy(in_f,"idxfile");
	setFileName(in_f,&f,1);
	
	FILE *fpw[6];
	FILE *fp_capidx;
	FILE *fp_relat;
	FILE *fp_smbwt;
	FILE *fp_nxtpnt;
	FILE *fp_nxtflg;
	FILE *fp_smpos;
	while(repNum < NUM_EXT){
		OpenSeedIdxfiles( &f, nxt_seedidx,LenExtIdx,repNum, fpw);
		fp_relat  = fpw[0];
		fp_smbwt  = fpw[1];
		fp_nxtpnt = fpw[2];
		fp_nxtflg = fpw[3];
		fp_capidx = fpw[4];
		fp_smpos  = fpw[5];

		CrtJmpMod(LenExtIdx[repNum], nxt_seedidx, jmp_idx, mod_idx,cap_pos);

		uint32_t bgn,end,num,j,idx,len_seedidx ;


		len_seedidx = LenExtIdx[repNum];
		idx = 0;
		while(idx<len_seedidx){
			bgn = cur_seedidx[idx][0];
			num = cur_seedidx[idx][1];
			end = bgn + num-1; 

			if(num <= IS_SMLSIZ){
				sub.num_pos = 0 ;
				for(j=bgn; j<=end; j++) 
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//预编译期间临时屏蔽
					//sub.pos_buf[sub.num_pos++]=ItoPos[j]; 
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++				
				fwrite(sub.pos_buf,sizeof(uint32_t), sub.num_pos,fp_smpos);
				continue;
				//小规模数据不做任何处理
				//比对过程中可以通过自身的jmp_idx,mod_idx信息，
				//可以找到相应的SmPos[]数组的行号
			} 
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			//以下是当num_idx>IS_SMLSIZ时进行------------------------- 
			//生成左扩展序列集合SeqL[]和右扩展序列集合SeqR[]，
			//计算关联关系矩阵数组Relat[ ]和映射数组LtoRel[ ], 

			uint32_t sidx[2];
			sidx[0] = bgn;
			sidx[1] = num;
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//预编译期间临时屏蔽
	getBlckData(sidx,(uint32_t)repNum,&sub,idx);  
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	        

			//Relat[  ]和LtoRel[ ]分段输出到同一个文件RelatFile_i。

			fwrite(sub.relat,sizeof(uint32_t),sub.num_relat,fp_relat);
			// buf_relat是当前Relat[]数据的缓冲区
			// len_relat是当前buf_relat数据的个数
			fwrite(sub.L2rel,sizeof(uint32_t),sub.num_seqL,fp_relat);
			// buf_L_Rel是当前LtoRel[ ]数据的缓冲区
			// num_seq_L是当前SeqL数据的个数
			//num_seq_L, num_seq_R，len_relat, 分别保存在cap_idx的相应分量中；
		
			//---------------------------------------------------------	
			//计算SeqL[]的SmBwt,生成buf_bwt_seq_L和buf_bwt_sum_L;
			//len_bwtseq_L = (num_seq_L+1)/2;
			//len_bwtsum_L = getbwtsumsize(num_seq_L);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//预编译期间临时屏蔽
	        buldSmBwt(&sub,0); //0生成左Bwt，1生成右Bwt
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++		

			

			fwrite(sub.bwt_seqL,sizeof(uint8_t),sub.num_bwtseqL,fp_smbwt);
			fwrite(sub.bwt_sumL,sizeof(uint8_t),sub.num_bwtsumL,fp_smbwt);

			//计算SeqR[]的SmBwt,生成buf_bwt_seq_R和buf_bwt_sum_R;
			//len_bwtseq_R = (num_seq_R+1)/2;
			//len_bwtsum_R = getbwtsumsize(num_seq_R);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	//预编译期间临时屏蔽
	        buldSmBwt(&sub,1); //0生成左Bwt，1生成右Bwtt
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


			fwrite(sub.bwt_seqR,sizeof(uint8_t),sub.num_bwtseqR,fp_smbwt);
			fwrite(sub.bwt_sumR,sizeof(uint8_t),sub.num_bwtsumR,fp_smbwt);
			
			//计算smbwt数组的指针增量  
			//cur_num_smbwt = len_bwtseq_L+ len_bwtsum_L + len_bwtseq_R+ len_bwtsum_R
				
			//-------------------------------------------------------
			//配对左序列seq_L[]与相关联的右序列seq_L[]的Indx区间
			//关联数据保存在buf_rel_idx[][]数组	
			uint32_t *ItoPos;
			relatNxtCapIdx(&sub,jmp_idx,mod_idx,cap_pos,ItoPos);

			//输出nxt_idx[i],nxt_idx[i],和apIdx[]数据
			fwrite(sub.pos_buf,sizeof(uint32_t), sub.num_pos,fp_smpos);

			fwrite(sub.nxt_idx,sizeof(uint32_t),sub.len_relat,fp_nxtpnt);
			fwrite(sub.nxt_flg,sizeof(uint8_t), sub.len_relat,fp_nxtflg);

			cap[0].num_seqL  = sub.num_seqL;
			cap[0].num_seqR  = sub.num_seqR;
			cap[0].num_relat = sub.num_relat;

			fwrite(cap,sizeof(CapInfo),1,fp_capidx);	

	        cap[0].relat += sub.num_relat + sub.num_seqL ;	 
			cap[0].smbwt += sub.num_bwtL + sub.num_bwtR ; 
			cap[0].nxtcap+= sub.num_relat;
  
	        //计算relat数组的指针增量 
			sub.len_capidx++;

			idx++;
		} // End: while(idx<len_seedidx) +++++++++++++++++++++++++++
			//---------------------------------------------------------
		fclose(fp_nxtpnt);
		fclose(fp_nxtflg);
		fclose(fp_capidx);
		fclose(fp_relat);
		fclose(fp_smbwt);
		fclose(fp_smpos);

		ord_seedidx = cur_seedidx;
		cur_seedidx = nxt_seedidx ;
		nxt_seedidx = NULL;
	}
	strcpy(in_f,"comfile");
	setFileName(in_f,&f,1);
    comFile(&f,&f);
	return 0;
}


//==================================================================
void comFile(FileName *rf, FileName *wf){
	//+++++++++++++++++++++++++++++++++++++++++++++
	FILE *fp_com;
	FILE *fp;
	int f_id;
	uint32_t num = 0;
	char  buf[2];

    int fsize[NUM_FILES+1][NUM_EXT];
    int head[NUM_FILES+1][NUM_EXT];

    char *fName[NUM_FILES][NUM_EXT] ;

	fName[0][0]	=  rf->capidx[0];
	fName[1][0]	=  rf->nxtpnt[0];
	fName[2][0]	=  rf->nxtflg[0];
	fName[3][0]	=  rf->relat[0];
	fName[4][0]	=  rf->smbwt[0];
	fName[5][0]	=  rf->smpos[0];

    int i,j;
	for(i=0;i<NUM_EXT;i++){
		fsize[1][i] =  getFileSize(rf->capidx[i]);
		fsize[2][i] =  getFileSize(rf->nxtpnt[i]);
		fsize[3][i] =  getFileSize(rf->nxtflg[i]);
		fsize[4][i] =  getFileSize(rf->relat[i]);
		fsize[5][i] =  getFileSize(rf->smbwt[i]);		
		fsize[6][i] =  getFileSize(rf->smpos[i]);
	}
	for(i=0;i<NUM_EXT;i++){
		head[0][i] = 5;
		head[1][i] = ((fsize[1][i]+3)/4)*4 ; 
		head[2][i] = ((fsize[2][i]+3)/4)*4 ; 
		head[3][i] = ((fsize[3][i]+3)/4)*4 ; 
		head[4][i] = ((fsize[4][i]+3)/4)*4 ; 
		head[5][i] = ((fsize[5][i]+3)/4)*4 ; 
		head[6][i] = ((fsize[6][i]+3)/4)*4 ; 
	}
	int num_ext = 0 ;
	while(num_ext<NUM_EXT){
		if((fp_com=fopen64(wf->comfile[num_ext],"w"))==NULL){
		    printf("can't open file\n");
		    exit(0);
		}
		f_id = 0;
		while(f_id<NUM_FILES){
			//++++++++++++++++++++++++++++++++++++++++++
			//file_capidx数据输出到file_extidx
			if((fp=fopen64(fName[f_id][num_ext],"r"))==NULL){
			    printf("can't open file\n");
			    exit(0);
			}
			num = 0;
			//whiel(num < head[i][1]){
			
			while((buf[0] = fgetc(fp)) != EOF){
				num++;
				fwrite(buf,sizeof(uint8_t),1,fp_com);
			}
			if(num != fsize[num_ext][1]){
				    printf("file read Err\n");
				    exit(0);
			}

			for(j=0; j<(head[num_ext][1]-fsize[num_ext][1]);j++ ){
				buf[0] = 0;
				fwrite(buf,sizeof(uint8_t),1,fp_com);
			}
			fclose(fp);
			f_id++;
			//++++++++++++++++++++++++++++++++++++++++++++
		}//End:  while(num_file<NumFiles) +++++++++++++++++++
		fclose(fp_com);
		num_ext++;
	}
	return;
} 
// End： void comDataFile(FileName *file) +++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void OpenSeedIdxfiles(FileName *f, uint32_t (*sidx)[3], uint32_t LenExtIdx[],uint32_t repNum,FILE *fpw[]){
	FILE *fp;
	if((fp=fopen64(f->seedidx[repNum+1],"r"))==NULL){//判断是否打开文件
	    printf("can't open file\n");
	    exit(0);
	}
	uint32_t  numread=fread(sidx,3*sizeof(uint32_t),LenExtIdx[repNum+1],fp);
	//printf("Number of items read=%d\n",numread);
	fclose(fp);
	
	if((fpw[0]=fopen64(f->relat[repNum+1],"w"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}
	
	if((fpw[1]=fopen64(f->smbwt[repNum+1],"w"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}

	if((fpw[2]=fopen64(f->nxtpnt[repNum+1],"w"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}

	if((fpw[3]=fopen64(f->nxtflg[repNum+1],"w"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}

	if((fpw[4]=fopen64(f->capidx[repNum+1],"w"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}

	if((fpw[5] = fopen64(f->smpos[repNum+1],"w"))== NULL){
	    printf("can't open file\n");
	    exit(0);
	}
	return;
} // End: void Openfiles( ) ++++++++++++++++++++++++++++++++

//*************************************************
//*************************************************


void CrtJmpMod(int n, uint32_t cur_seedidx[][3], uint32_t* jump_idx, 
	uint8_t *mod_idx, uint32_t *cap_pos){

    //uint32_t row_jump = 0;
    uint32_t last_idx_q = 0;
    //uint32_t jump_idx[LEN_NEXT_IDX+1]; 
    uint32_t row, i;
    for(row=0; row<n; ++row){
        cap_pos[row] = cur_seedidx[row][2]; 
        uint32_t idx = cur_seedidx[row][0];
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
    jump_idx[(n+255)/256] = n; 
    return;
}
/*
int cut_ext_seq(uint8_t *pac, bwt_t *bwt, bwtint_t k, bwtint_t l, int offset, vec_uint_t *seq_buf)
{
    bwtint_t i, pos, st, ed;
    uint32_t seq16;
    //fprintf(stderr, "\n"); 
    //fetch 16 mer preseq in [k, l) and rm duplicate preseq
    for(i = k; i < l; ++i){
        pos = bwt_sa(bwt, i);
        st = pos+offset; ed = st+16;
        if(offset < 0 && pos < (bwtint_t)-offset) st = 0;
        if(offset > 0 && ed > bwt->seq_len) ed = bwt->seq_len; 
        seq16 = bns_extract_seq16(pac, st, ed);//???
        //fprintf(stderr, "%u\t", seq16); 
        kv_push(uint32_t, *seq_buf, seq16);
    }
    //fprintf(stderr, "\n");
    //for(i = seq16_st; i <seq16s->n; ++i) fprintf(stderr, "%u\t", seq16s->a[i]);
    //fprintf(stderr, "\n");
    return seq_buf->n;
}*/
int rm_euqal(vec_uint_t *seq_buf, bwtint_t k, vec_ext_t *rext)
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
}


void getBlckData(uint32_t sidx[2],uint32_t num_ext, 
        	SubBuf *sub_buf, idx_t *idx/*pubParm pParm*/ ){
//
#define MAX_COUNT 500000
    vec_uint_t seq_buf; 
    vec_ext_t rext, lext;
    kv_init(seq_buf); kv_init(rext);kv_init(lext);
    kv_resize(uint32_t, seq_buf, MAX_COUNT);
    kv_resize(ext_t, rext, MAX_COUNT);
    kv_resize(ext_t, lext, MAX_COUNT);
    ext_t tmp; 
    uint32_t *lext_uniq_idx = (uint32_t *)calloc(MAX_COUNT, sizeof(ext_t));       
    int tot_lext0 = 0, tot_lext1 = 0, tot_rext = 0;
 
    uint8_t *pac = idx->pac;
    bwt_t *bwt = idx->bwt;    
    //rext.n = 0; lext.n =0; seq_buf.n = 0;
    uint32_t k = sidx[0], l = k+sidx[1];
    cut_ext_seq(pac, bwt, k, l,  20, &seq_buf);//right seq buf
    rem_euqal(&seq_buf, k, &rext);//gen right extend seq 
    uint32_t i;
    for(i =0; i < rext.n; ++i) { 
        seq_buf.n = 0;
        int st, ed;
        st = rext.a[i].idx;
        ed = i+1 >= rext.n?l:rext.a[i+1].idx;
        cut_ext_seq(pac, bwt, st, ed, -16, &seq_buf);//left seq buf 
        ks_introsort(uint32_t, seq_buf.n, seq_buf.a); 
        rem_euqal(&seq_buf, st, &lext);
    }
    int row_r, row_l;
    /*
        fprintf(stderr, "----------------------\n");
 
        for(row_l = 0; row_l < lext.n; ++row_l){
            fprintf(stderr, "[lext]:\t%u\t", lext.a[row_l].idx);
            log_seq162nt(lext.a[row_l].seq); 
        }
        for(row_r = 0; row_r < rext.n; ++row_r){
            fprintf(stderr, "[rext]:\t%u\t", rext.a[row_r].idx);
            log_seq162nt(rext.a[row_r].seq); 
        }
    */  
    row_r = 1, row_l = 0; 
    while(row_l < lext.n && row_r < rext.n){
        if(lext.a[row_l].idx  < rext.a[row_r].idx) lext.a[row_l++].idx = row_r-1;
        else row_r++; 
    } 
    if(row_l < lext.n) for(i = row_l; i < lext.n; ++i) lext.a[i].idx = rext.n-1; 
    ks_introsort(ext, lext.n, lext.a);

    int n_uniq = 0;
    uint32_t last_seq16 = lext.a[0].seq;
    lext_uniq_idx[n_uniq++] = 0;
    for(i = 1; i < lext.n; ++i){
        if(lext.a[i].seq!= last_seq16) {
            //lext_uniq_idx[n_uniq++] = i+tot_lext0; 
            lext_uniq_idx[n_uniq++] = i; 
            last_seq16 = lext.a[i].seq;
        }
    }
    //dump_cnt(fp_totext, tot_lext0, tot_lext1, tot_rext);
    fprintf(stderr, "%u\t%u\t%u\n", tot_lext0, tot_lext1, tot_rext); 
    tot_lext1 += lext.n; 
    tot_lext0 += n_uniq;
    tot_rext += rext.n; 
    //fwrite(rext.a, sizeof(ext_t), rext.n, fp_rext0);         
    
    //for(i = 0; i < rext.n; ++i) log_seq162nt(rext.a[i].seq);

    for(i = 0; i < rext.n; ++i) {
        sub_buf->sortR[i] = rext.a[i].seq;
        sub_buf->idxR[i][0] = rext.a[i].idx;
        sub_buf->idxR[i][1] = i == rext.n?l:rext.a[i].idx-1;
    }
    for(i = 0; i < n_uniq; ++i) {

        sub_buf->sortL[i] = lext.a[lext_uniq_idx[i]].seq;
        sub_buf->L2rel[i] = lext_uniq_idx[i]; 
        //fwrite(&tmp, sizeof(ext_t), 1, fp_lext0);
        //fwrite(&lext.a[lext_uniq_idx[i]].seq, sizeof(uint32_t), 1, fp_lext0); 
        //fwrite(lext_uniq_idx+i, sizeof(uint32_t), 1, fp_lext0);
        fprintf(stderr, "%u\t", lext_uniq_idx[i]);
        log_seq162nt(lext.a[lext_uniq_idx[i]].seq);

    } 
    for(i = 0; i < lext.n; ++i) sub_buf->relat[i] = lext.a[i].idx;
    
    kv_destroy(seq_buf);
    kv_destroy(lext);
    kv_destroy(rext);
    return;
}

void build_local_bwt(uint32_t *sort_seq, int l, int r, uint32_t *sa0, uint32_t *sa1, uint8_t *bwt)
{
    int i, j;
    uint32_t cnt_2nt[17];
    uint32_t *last_sa, *cur_sa;
    sort_seq += l;
    int n = r-l;
    for(j=0; j< 17; ++j) cnt_2nt[j] = 0;
    for(j = 0; j < n; ++j){
        uint32_t x = sort_seq[j]&0xF;
        //__set_bwt(bwt, j*2+1, x&3);
        //__set_bwt(bwt, j*2, (x>>2)&3);
        __set_bwt2(bwt, j, x);
        ++cnt_2nt[x];
    }
    for(j =0; j < n; ++j) sa0[j] = j;
    last_sa = sa0; cur_sa = sa1;
    for(i = 1; i < 8; ++i){
        fprintf(stderr, "\n==Iter = %u\n", i-1);

        //log_array(17, cnt_2nt);
        accumulate_cnt(17, cnt_2nt);
        log_array(17, cnt_2nt);
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

uint8_t *build_seq16_bwt (int n, uint32_t *sort_seq)
{
   
    uint8_t *bwt; uint32_t *sa0, *sa1;
    bwt = (uint8_t *)calloc((n+1)/2*8, sizeof(uint8_t));  
    sa0 = (uint32_t *)calloc(n, sizeof(uint32_t));
    sa1 = (uint32_t *)calloc(n, sizeof(uint32_t));
    build_local_bwt(sort_seq, 0, n, sa0, sa1, bwt); 
    free(sa0);free(sa1);
    build_occ(bwt, n); 
        
        
       //dump bwt
        
       
 
        

    return bwt;
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define MAX_BLCK_SIZE 500000
#define MAX_SEED_SIZE 500000 
#define ALL_BUFS_SIZE 12  // All bufs size in SubBuf



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void InitSubBuf(SubBuf *sub_buf){

	SubBuf *sub ;

	if(NULL == (sub = (SubBuf*)malloc(sizeof(SubBuf)))){
	    perror("error...");
	    exit(1);
	}
	uint32_t *buf;
	if(NULL == (buf = (uint32_t*)malloc(MAX_BLCK_SIZE*ALL_BUFS_SIZE*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
	sub->seqL_buf = (uint32_t(*)[2])buf;
	sub->seqR_buf = (uint32_t(*)[2])(buf + MAX_BLCK_SIZE*2);

	sub->sortL    = buf + MAX_BLCK_SIZE*4;
	sub->sortR    = buf + MAX_BLCK_SIZE*5;

	sub->relat    = buf + MAX_BLCK_SIZE*6;
	sub->L2rel    = buf + MAX_BLCK_SIZE*7;
	sub->seqL_idxR= (uint32_t(*)[3])(buf + MAX_BLCK_SIZE*8);
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	sub->idxR     = (uint32_t(*)[2])buf; //注：idxR[][2]的数据存放在;seqL_buf[][2]的物理空间。
	sub->extIdx   = (uint32_t(*)[2])(buf + MAX_BLCK_SIZE*2); //注：extIdx[][2]的数据存放在;seqR_buf[][2]的物理空间。

	

	sub->pos_buf = buf + MAX_BLCK_SIZE*9 ; 	
	sub->nxt_idx = buf + MAX_BLCK_SIZE*10;
	sub->nxt_flg = (uint8_t*)(buf + MAX_BLCK_SIZE*11);


//+++++++++++++++++++++++++++++++++++
	uint8_t *bwt;
	if(NULL == (bwt = (uint8_t*)malloc(12*MAX_BLCK_SIZE*sizeof(uint8_t)))){
	    perror("error...");
	    exit(1);
	}
	sub->smbwt    = bwt;
	sub->bwt_seqL = bwt ;   
	sub->bwt_sumL = bwt + 4*MAX_BLCK_SIZE ;
	sub->bwt_seqR = bwt + 6*MAX_BLCK_SIZE ;    
	sub->bwt_sumR = bwt + 10*MAX_BLCK_SIZE ;  

//+++++++++++++++++++++++++++++++++++++++++++++++++
	sub->num_seqL = 0;
	sub->num_seqR = 0;
	sub->num_relat= 0;	
	//sub->num_ext  = 0;

	//sub->off_pos  = 0; 
    //sub->blck_id  = 0; 

	sub->len_capidx= 0;
	sub->len_relat = 0;	
	sub->len_smbwt = 0;	
	sub->len_nxtpnt= 0;
	sub->len_nxtflg= 0;
	sub->len_smpos = 0;
	sub->num_pos   = 0;
	sub->num_bwtL  = 0;
	sub->num_bwtR  = 0;


	sub->cap[0].relat    = 0;	
	sub->cap[0].smbwt    = 0;
    sub->cap[0].nxtcap   = 0;	
	sub->cap[0].num_seqL = 0;
	sub->cap[0].num_seqR = 0;
	sub->cap[0].num_relat= 0;	
	sub->cap[0].num_pos  = 0;

//+++++++++++++++++++++++++++++++++++++++++++++++++
	sub_buf = sub;

	return ;
}


void buldSmBwt(SubBuf *sub,int flg){
    uint8_t *bwt; uint32_t *sa0, *sa1;
    //bwt = (uint8_t *)calloc((n+1)/2*8, sizeof(uint8_t));  
    int n;
    uint32_t *sort_seq;
    bwt = flg==0?sub->bwt_seqL:sub->bwt_seqR;
    if(flg == 0 ){
        bwt = sub->bwt_seqL; 
        n = sub->num_seqL; 
        sort_seq = sub->sortL;
    } else{
        bwt = sub->bwt_seqR; 
        n = sub->num_seqR;
        sort_seq = sub->sortR;
    }
    sa0 = (uint32_t *)calloc(n, sizeof(uint32_t));
    sa1 = (uint32_t *)calloc(n, sizeof(uint32_t));
    build_local_bwt(sort_seq, 0, n, sa0, sa1, bwt); 
    free(sa0);free(sa1);
    build_occ(bwt, n); 

        
       //dump bwt
        
       
 
        

   

  
	return;
}

//**************************************************
//*************************************************************
void relatNxtCapIdx(SubBuf *subBuf,uint32_t *jmp_idx, uint8_t *mod_idx,uint32_t *cap_pos,uint32_t *ItoPos){
	uint32_t  i, j,	bgn,end, num,  bgn_q, end_q, i_r, 
			  seq, pos, row, j_idx;

	uint32_t buf_idx[5];

	SubBuf *sub;
	sub = subBuf;

	sub->num_pos = 0;

	for(i=0;i<sub->num_seqL;i++){
		bgn = sub->L2rel[i];
		end = sub->L2rel[i+1];
		for(j=bgn; j<end; j++){
			row = sub->relat[j];
			sub->seqL_idxR[j][0] = sub->sortL[i];
			sub->seqL_idxR[j][1] = sub->idxR[row][0]; 
			sub->seqL_idxR[j][2] = sub->idxR[row][1];										
		}				
	}
	//-------------------------------------------------------
	//buf_rel_idx[][]中的每一条数据，进行Bwt比对获得(bgnIdx,endIdx)的数据，
	//保存在sub.extIdx[][]
	for(i=0;i<sub->num_relat;i++){
		seq = sub->seqL_idxR[i][0];
		buf_idx[0] = sub->seqL_idxR[i][1];
		buf_idx[1] = sub->seqL_idxR[i][2];
//++++++++++++++++++++++++++++++++++++++++++++++
//预编译期间屏蔽
		//AlgnSeq(pParm,seq,buf_idx); //seq 是整数
//++++++++++++++++++++++++++++++++++++++++++++++++
		sub->extIdx[i][0]= buf_idx[3]; //返回的bgnIdx
		sub->extIdx[i][1]= buf_idx[4]; //返回的endIdx
	}
	//---------------------------------------------------------
	//buf_rel_idx[][]中的每一条数据，进行Bwt比对获得(bgnIdx,endIdx)的数据，
	//保存在sub.extIdx[][]
	sub->num_pos = 0;
	for(i=0;i<sub->num_relat;i++){
		bgn = sub->extIdx[i][0];
		end = sub->extIdx[i][1];
		num = end - bgn + 1;
		if(num==1){
			pos = ItoPos[bgn];
			sub->nxt_idx[i] = pos;
			sub->nxt_flg[i] = num;
			continue;
		}
		//if(num_idx>1)+++++++++++++++++++++++++++++++++++++++++++++++++++
		bgn_q = jmp_idx[i/256];
		end_q = jmp_idx[(i/256)+1];
		i_r = i%256;
		for(j_idx=bgn_q; j_idx<end_q; j_idx++){
			if(i_r == mod_idx[j_idx]){
				break;
			} 
		}		
		sub->nxt_idx[i] = cap_pos[j_idx]; //指向下一级CapIdx[]的行号
		if(num < 255){
			sub->nxt_flg[i] = (uint8_t)num ; //指向CapIdx[]的类型
		}else{
			sub->nxt_flg[i] = 0xFF ;
		}
		if(num<=IS_SMLSIZ){
			for(j=bgn; j<end; j++){
				pos = ItoPos[j]; 
				sub->pos_buf[sub->num_pos++] = pos; 
			}
		}
	} // End : for(i=0;i<len_relat;i++)
	return;
}

uint32_t getlen(uint32_t l_extidx[], FileName *f)
{
	int i;
	long max_len = -1, max_idx = -1;
	for(i =0; i < NUM_EXT; ++i){
        long len = getFileSize(f->seedidx[i]);
        l_extidx[i] = (uint32_t)len;
		if(len > max_len) {
			max_len = len;
			max_idx = i; 
		}

	}
	return max_idx;

}
long getFileSize(char *fn)
{
	FILE *fp = fopen(fn, "r");
	fseek(fp, 0, SEEK_END);
	long len = ftell(fp);
	fclose(fp);
    return len; 


}
