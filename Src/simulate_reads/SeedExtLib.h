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
//#include "debug.h"
#include <stdio.h>
#include "ksw.h"
#include "seed.h"
#include "query.h"
#define SEED_LEN 20
#define NO_BWT_SUM_SIZE 64
#define ALN_THRES_SCORE 30
#define MIN_BWT_SIZE 16

//struct SubBuf *InitIdxsArry(struct FileName *fn, struct SubBuf **out_sub, struct ExtBlck **out_blck, struct StackTree **out_stree, struct CapIfo *cap, struct JmpMod **out_jmp);
struct SubBuf *InitIdxsArry(struct FileName *fn, struct SubBuf **out_sub, struct ExtBlck out_blck[], struct StackTree **out_stree, struct CapIfo *cap, struct JmpMod **out_jmp);
void initStackTree(struct StackTree *sTree);


void setStackTree(struct StackTree *sTree,uint32_t (*pair_out)[3], int seed_id, uint32_t (*ext_idx)[2]);
void setBlckData(struct ExtBlck *in_blck, uint8_t cls, uint32_t nxt_cap, struct SubBuf *sub, struct ExtBlck **out_this);
//uint32_t getFileSize(char *file);
void getSeed(uint8_t *seed_seq, uint8_t *read_seq, uint8_t  seed_off);
//void AlgnSeed(uint8_t *seed_seq,uint32_t buf_algn[], idx_t *idx);
void getCapPos(struct JmpMod *jmp,uint32_t buf_algn[]);
//int AlgnPos(struct PubParm *pParm, struct ExtBlck *eB,struct SubBuf *sub, uint8_t *read_seq,uint32_t pos);

int AlgnPos(idx_t *fm_idx, uint8_t ext_cls, uint8_t *read_seq,int read_len, int seed_off, uint32_t pos_buf[], int pos_num, struct SubBuf *sub);

int init_seed_model(seed_t *seed );

uint32_t aln_seed_seq(idx_t *fm_idx,uint32_t hash_boundry[], int read_len,  uint8_t *f_read_seq,  uint8_t *r_read_seq, uint8_t *qual, uint8_t* rqual, seed_t *seed );

uint32_t aln_seed_seq_1(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed );
int AlgnPos_buf(idx_t *fm_idx, uint8_t ext_cls, uint8_t *read_seq,int read_len, int seed_off, int num_pair, uint32_t pos_buf[], struct SubBuf *sub);
uint32_t __dna2_count(uint8_t *bwt, uint32_t k , uint32_t l, uint8_t c);
void OutAlgnInfo(uint32_t (*algn_out)[5]);

int AlgnmExtSeq(struct ExtBlck *eB, struct SubBuf *sub, int flg);		

//void AlgnBwtSml(uint8_t*Bwt,uint8_t*Seq, uint32_t AlgIn[],uint32_t(*AlgOut)[3]);
//void AlgnBwtBig(uint8_t*Bwt, uint8_t*SumBuf,uint8_t*Seq, uint32_t DataNum, uint32_t AlgIn[],uint32_t(*AlgOut)[3]);
uint32_t __dna2_count_small(uint8_t *bwt, uint32_t k , uint32_t l, uint8_t c);

//void AlgnBwtSml(uint8_t*Bwt, uint8_t*Seq, uint32_t DataNum, uint32_t AlgIn[], uint32_t(*AlgOut)[3]);
//void AlgnBwtBig(uint8_t*Bwt, uint8_t*Seq, uint32_t DataNum, uint32_t AlgIn[], uint32_t(*AlgOut)[3]);
//uint32_t getBwtSumSize(uint32_t DataNum);
void getExtSeq(struct SubBuf *sub, uint32_t cls, uint8_t *read_seq);
void gen_cnt(uint8_t *bwt, int n_data, uint8_t cnt_2nt[][17]);
void align_min(uint8_t *bwt, int n_data, uint8_t seq[16], int err_num, int flag, uint8_t *st_pos, uint8_t *algn_row, uint32_t (*aln_out)[3]);
void align_nosum(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t seq[16], int seq_len, int seq_st,  uint8_t *st_pos, uint8_t *algn_row, uint32_t (*aln_out)[3]);
void align_255(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data,uint8_t seq[16], int seq_len, int seq_st,  uint8_t *st_pos, uint8_t *algn_row, uint32_t (*aln_out)[3]);
void align_large_alt(uint8_t *Bwt, uint8_t *cnt2,  uint32_t n_data, const uint8_t *seq , int seq_len, int seq_st, uint8_t *st_pos, uint8_t *algn_row, uint32_t (*aln_out)[3]);

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
                        //fprintf(stderr, "siz_cnt0 %u-%u\t", buf, buf & 255);
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
                //fprintf(stderr, "siz_cnt0 %u\t", buf & 255);
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

//struct SubBuf *InitIdxsArry(struct FileName *fn, struct SubBuf **out_sub, struct ExtBlck **out_blck, struct StackTree **out_stree, struct CapIfo *cap, struct JmpMod **out_jmp)
struct SubBuf *InitIdxsArry(struct FileName *fn, struct SubBuf **out_sub, struct ExtBlck eBlck[], struct StackTree **out_stree, struct CapIfo *cap, struct JmpMod **out_jmp)
{
    char in_f[128];
    //内存初始化------------------------------------------------------
	//初始化，jmp_idx,mod_idx,cap_pos数组开始
	FILE *fp;
	if((fp=fopen(fn->jmpmod,"rb"))==NULL){
	    printf("can't open file\n");
	    exit(0);
	}
    uint32_t buf_head[100], head_size; 
    fread(buf_head,sizeof(uint32_t),4,fp);  
    int max_buf_size = buf_head[3];

fprintf(stderr, "n_cap, n_jmp, n_mod = %u, %u, %u\n", buf_head[0], buf_head[1], buf_head[2]);
	uint32_t data_size = 0;
	uint32_t *jmp_idx;           	
    uint32_t *cap_pos;
	uint8_t  *mod_idx; 	
    	
	data_size = buf_head[0] + buf_head[1]+ (buf_head[2]+3)/4;
	if(NULL == (cap_pos = (uint32_t*)malloc(data_size*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    
    fread(cap_pos, sizeof(uint32_t), buf_head[0], fp);
    jmp_idx = cap_pos+buf_head[0];
    fread(jmp_idx, sizeof(uint32_t), buf_head[1], fp);
    mod_idx = (uint8_t *)(cap_pos+buf_head[0]+buf_head[1]); 
    fread(mod_idx, sizeof(uint8_t), buf_head[2], fp);
    fclose(fp);
	struct JmpMod *jmp;
    if(NULL == (jmp= (struct JmpMod *)malloc(sizeof(struct JmpMod)))){
	    perror("error...");
	    exit(1);
	}

    jmp->jmp = jmp_idx;
	jmp->mod = mod_idx;
	jmp->cap = cap_pos;
	*out_jmp = jmp;
    //初始化，jmp_idx,mod_idx,cap_pos数组结束------------
    int i;
    //for(i = 0; i < buf_head[0]; ++i) fprintf(stderr, "cap_pos[%u] = %u\n", i, cap_pos[i]);
    //for(i = 0; i < buf_head[1]; ++i) fprintf(stderr, "jmp_idx[%u] = %u\n", i, jmp_idx[i]);
    //for(i = 0; i < buf_head[2]; ++i) fprintf(stderr, "mod_idx[%u] = %u\n", i, mod_idx[i]);
    uint32_t f_id ; 
	uint32_t *buf_data ; 
    //以下代码段读取comfile， 
    size_t file_size[NUM_EXT];
    uint32_t *head_idx[NUM_EXT];
   
    int j;  
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
            //printf("buf_head[%u] = %u \n", j, buf_head[j]);
        } 
        buf_head[NUM_FILES] = sum;  
        memset(eBlck+f_id, 0, sizeof(struct ExtBlck)); 
        eBlck[f_id].head_cap    =  (struct CapIfo *)(head_idx[f_id]+buf_head[0]);
		eBlck[f_id].head_nxtpnt =  (uint32_t  *)(head_idx[f_id]+buf_head[1]);		   
        eBlck[f_id].head_nxtflg =  (uint8_t *)(head_idx[f_id]+buf_head[2]);
		eBlck[f_id].head_relat  =  (uint32_t *)(head_idx[f_id]+buf_head[3]);
		eBlck[f_id].head_smbwt  =  (uint8_t  *)(head_idx[f_id]+buf_head[4]);
        eBlck[f_id].head_extidx =  (uint32_t (*)[2])(head_idx[f_id]+buf_head[5]);
        
    } 
    uint32_t buf32[2];
    uint8_t buf8[2];
    struct CapIfo tmp_cap[2];
/*      
    for(f_id = 0; f_id < NUM_EXT; ++f_id){
        file_size[f_id] = getFileSize(fn->capidx[f_id])/sizeof(struct CapIfo); 
        FILE *fp = fopen(fn->capidx[f_id], "rb");
        fprintf(stderr, "num_ext = %u, %s\n", f_id, fn->capidx[f_id]);
        fprintf(stderr, "file_size = %u\n", file_size[f_id]);
        for(i = 0; i < file_size[f_id]; ++i){ 
            fread(tmp_cap, sizeof(struct CapIfo), 1, fp); 
            if(tmp_cap[0].relat != eBlck[f_id].head_cap[i].relat ){
                fprintf(stderr, "relat,num_ext = %u,  i = %u\n", f_id, i);
                fprintf(stderr, "tmp = %u,  eBlck = %u\n", tmp_cap[0].relat , eBlck[f_id].head_cap[i].relat);
            }  
            if(tmp_cap[0].smbwt != eBlck[f_id].head_cap[i].smbwt ){
                fprintf(stderr, "smbwt, num_ext = %u,  i = %u\n", f_id,  i);
                fprintf(stderr, "tmp = %u,  eBlck = %u\n", tmp_cap[0].smbwt, eBlck[f_id].head_cap[i].smbwt);
            }  
            if(tmp_cap[0].nxtpnt != eBlck[f_id].head_cap[i].nxtpnt ){
                fprintf(stderr, "nxtpnt, num_ext = %u,  i = %u\n", f_id, i);
                fprintf(stderr, "tmp = %u,  eBlck = %u\n", tmp_cap[0].nxtpnt, eBlck[f_id].head_cap[i].nxtpnt);
            }  
            if(tmp_cap[0].num_seqL != eBlck[f_id].head_cap[i].num_seqL ){
                fprintf(stderr, "num_seqL, num_ext = %u,  i = %u\n", f_id, i);
                fprintf(stderr, "tmp = %u,  eBlck = %u\n", tmp_cap[0].num_seqL, eBlck[f_id].head_cap[i].num_seqL);
            }
            if(tmp_cap[0].num_seqR != eBlck[f_id].head_cap[i].num_seqR ){
                fprintf(stderr, "num_seqR, num_ext = %u,  i = %u\n", f_id, i);
                fprintf(stderr, "tmp = %u,  eBlck = %u\n", tmp_cap[0].num_seqR , eBlck[f_id].head_cap[i].num_seqR);
            } 
            if(tmp_cap[0].num_relat != eBlck[f_id].head_cap[i].num_relat ){
                fprintf(stderr, "num_relat, num_ext = %u, i = %u\n", f_id, i); 
                fprintf(stderr, "tmp = %u,  eBlck = %u\n", tmp_cap[0].num_relat, eBlck[f_id].head_cap[i].num_relat);
            }
        }
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
        for(i = 0; i < file_size[f_id]; ++i){ 
            struct CapIfo *cap_buf = eBlck[f_id].head_cap+i;            
            fread(tmp_cap, sizeof(struct CapIfo), 1, fp); 
            fprintf(stderr, "num_ext = %u,  i = %u\n", f_id, i);
            fprintf(stderr, "tmp   = %u\t", tmp_cap[0].relat);
            fprintf(stderr,         "%u\t", tmp_cap[0].smbwt);
            fprintf(stderr,         "%u\t", tmp_cap[0].nxtpnt);
            fprintf(stderr,         "%u\t", tmp_cap[0].num_seqL);
            fprintf(stderr,         "%u\t", tmp_cap[0].num_seqR);
            fprintf(stderr,         "%u\n", tmp_cap[0].num_relat);
            fprintf(stderr, "eBlck = %u\t", cap_buf->relat);
            fprintf(stderr,         "%u\t", cap_buf->smbwt);
            fprintf(stderr,         "%u\t", cap_buf->nxtpnt);
            fprintf(stderr,         "%u\t", cap_buf->num_seqL);
            fprintf(stderr,         "%u\t", cap_buf->num_seqR);
            fprintf(stderr,         "%u\n", cap_buf->num_relat);
        }
*/
/*  
    //Test nxtflg and nxtpnt    
    for(f_id = 0; f_id < NUM_EXT; ++f_id){
        uint32_t tmp[2];
        uint8_t tmp0[2];
        file_size[f_id] = getFileSize(fn->nxtpnt[f_id])/sizeof(uint32_t); 
        FILE *fp = fopen(fn->nxtpnt[f_id], "rb");
        FILE *fp1 = fopen(fn->nxtflg[f_id], "rb");
        fprintf(stderr, "num_ext = %u, %s\n", f_id, fn->nxtpnt[f_id]);
        fprintf(stderr, "file_size = %lu\n", file_size[f_id]);
        for(i = 0; i < file_size[f_id]; ++i){ 
            uint32_t *nxtpnt = eBlck[f_id].head_nxtpnt;            
            uint8_t *nxtflg = eBlck[f_id].head_nxtflg;            
            fread(tmp, sizeof(uint32_t), 1, fp); 
            fread(tmp0, sizeof(uint8_t), 1, fp1); 
            if(tmp[0] != nxtpnt[i] || tmp0[0] != nxtflg[i]){ 
                fprintf(stderr, "ERROR!!!!\n");
                fprintf(stderr, "num_ext = %u,  i = %u\t", f_id, i);
                fprintf(stderr, "tmp = %u \t", tmp[0]);
                fprintf(stderr, "pnt = %u \t", nxtpnt[i]);
                fprintf(stderr, "tmp = %u \t", tmp0[0]);
                fprintf(stderr, "flg = %u \n", nxtflg[i]);
                exit(1);
            }
        }
        fprintf(stderr, "nxtpnt, nxtflg ok!\n"); 
        fclose(fp);
        fclose(fp1);
    } 
    for(f_id = 0; f_id < NUM_EXT; ++f_id){
        uint32_t tmp[2];
        file_size[f_id] = getFileSize(fn->relat[f_id])/sizeof(uint32_t); 
        FILE *fp = fopen(fn->relat[f_id], "rb");
        for(i = 0; i < file_size[f_id]; ++i){ 
            uint32_t *relat = eBlck[f_id].head_relat;            
            fread(tmp, sizeof(uint32_t), 1, fp); 
            if(tmp[0] != relat[i]){ 
                fprintf(stderr, "num_ext = %u,  i = %u\t", f_id, i);
                fprintf(stderr, "tmp = %u \t", tmp[0]);
                fprintf(stderr, "relat = %u \n", relat[i]);
                exit(1);
            }
        }
        fprintf(stderr, "relat ok!\n"); 
        fclose(fp);
    } 
    for(f_id = 0; f_id < NUM_EXT; ++f_id){
        uint8_t tmp[2];
        file_size[f_id] = getFileSize(fn->smbwt[f_id])/sizeof(uint8_t); 
printf("file size = %lu\n", file_size[f_id]);
        FILE *fp = fopen(fn->smbwt[f_id], "rb");
        if(fp ==NULL) {
            fprintf(stderr, "FILE %s open fail!\n"); 
            exit(1);
        }
        for(i = 0; i < file_size[f_id]; ++i){ 
            uint8_t *smbwt = eBlck[f_id].head_smbwt;            
            fread(tmp, sizeof(uint8_t), 1, fp); 
            if(tmp[0] != smbwt[i]){ 
                fprintf(stderr, "num_ext = %u,  i = %u\t", f_id, i);
                fprintf(stderr, "tmp = %u \t", tmp[0]);
                fprintf(stderr, "smbwt = %u \n", smbwt[i]);
                exit(1);
            }
        }
        fprintf(stderr, "smbwt ok!\n"); 
        fclose(fp);
    } 
*/
/* 
    for(f_id = 0; f_id < NUM_EXT; ++f_id){
        uint32_t tmp[2];
        file_size[f_id] = getFileSize(fn->extidx[f_id])/sizeof(uint32_t); 
        FILE *fp = fopen(fn->extidx[f_id], "rb");
        for(i = 0; i < file_size[f_id]; ++i){ 
            uint32_t *extidx = eBlck[f_id].head_extidx;            
            fread(tmp, sizeof(uint32_t), 1, fp); 
            if(tmp[0] != extidx[i]){ 
                fprintf(stderr, "num_ext = %u,  i = %u\t", f_id, i);
                fprintf(stderr, "tmp = %u \t", tmp[0]);
                fprintf(stderr, "extidx = %u \n", extidx[i]);
                exit(1);
            }
        }
        fprintf(stderr, "extidx ok!\n"); 
        fclose(fp);
    } 
*/  
/*   
    for(f_id=0;f_id<NUM_EXT; f_id++){
		if((fp=fopen(fn->comfile[f_id],"rb"))==NULL){
		    printf("can't open file_extidx[f_id]\n");
		    exit(0);
		}
        int i, j;
        //文件名初始化-----------------------------------------------------		
  
        //FILE *fp_com=fopen(fn->comfile[j],"rb");
 
        uint32_t tmp[7];
        fread(tmp, sizeof(uint32_t), 7, fp);
        for(i = 0; i < 7; ++i) fprintf(stderr, "%u\n", tmp[i]);
        fclose(fp);
        fread(buf_head,sizeof(uint32_t),7,fp);  
		data_size = buf_head[0]-6;
        uint32_t buf0=0,buf1=0;
		uint32_t i; 
		for(i=0;i<head_size;i++){
			buf0 = buf_head[i];
			buf_head[i] = buf1;
			buf1+=buf0 ;
		}
        data_size = buf1;
		buf_head[i] = data_size;
        if(NULL == (buf_data = (uint32_t*)malloc(data_size*sizeof(uint32_t)))){
		    perror("error...");
		    exit(1);
		}
        if(data_size != fread(buf_data,sizeof(uint32_t),data_size,fp)){
		    printf("can't fread file\n");
		    exit(0);		
		}
        fclose(fp);
        uint32_t sum= 0, buf;
        for(i = 0; i < 6; ++i){
            buf = buf_head[i]; 
            buf_head[i] = sum;
            sum += buf;
        }
		eBlck[f_id].head_cap    =  (uint32_t *)(buf_data+buf_head[0]);
		eBlck[f_id].head_nxtpnt =  (uint8_t  *)(buf_data+buf_head[1]);		   
        eBlck[f_id].head_nxtflg =  (uint32_t *)(buf_data+buf_head[2]);
		eBlck[f_id].head_relat  =  (uint32_t *)(buf_data+buf_head[3]);
		eBlck[f_id].head_smbwt  =  (uint8_t  *)(buf_data+buf_head[4]);
		eBlck[f_id].head_extidx =  (uint32_t *)(buf_data+buf_head[5]);

	}  //End:+++++++++++++++++++++++++for(f_id=0;f_id<NUM_EXT; f_id++)
    *out_blck = eBlck;
*/

	
	struct StackTree *sTree ;
	if(NULL == (sTree = (struct StackTree*)malloc(sizeof(struct StackTree)))){
	    perror("error...");
	    exit(1);
	}

	uint32_t (*sArry)[4];

	if(NULL == (sArry = (uint32_t(*)[4])malloc(max_buf_size*4*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    uint32_t (*sBuf)[7];
    if(NULL == (sBuf = (uint32_t(*)[7])malloc(max_buf_size*7*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
	sTree->stck_arry = sArry;
	sTree->back_buf = sBuf;

	initStackTree(sTree);

    *out_stree = sTree;
	struct SubBuf *sub ;

	if(NULL == (sub = (struct SubBuf*)malloc(sizeof(struct SubBuf)))){
	    perror("error...");
	    exit(1);
	}
    int seq_out_size = 128;
	uint32_t (*seqL_out)[3];
	if(NULL == (seqL_out = (uint32_t(*)[3])calloc(seq_out_size, 3*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
    uint32_t (*seqR_out)[3];
	if(NULL == (seqR_out = (uint32_t(*)[3])calloc(seq_out_size,3*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
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
	uint32_t (*pair_out)[3];
	if(NULL == (pair_out = (uint32_t(*)[3])malloc(max_buf_size*3*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}


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

/*
	uint32_t *algnR_buf;
	if(NULL == (algnR_buf = (uint32_t*)malloc(max_buf_size*sizeof(uint32_t)))){
	    perror("error...");
	    exit(1);
	}
*/
	sub->seqL_out     = seqL_out;
	sub->seqR_out     = seqR_out;		
	sub->seqRel_out   = seqRel_out;		
	sub->seqRel_L     = seqRel_L;		
	sub->seqRel_R     = seqRel_R;		
	sub->pair_out     = pair_out;
	sub->algn_out     = algn_out;
    sub->sort_buf     = (uint32_t (*)[2])seqL_out; 
	sub->algnR_row    = algnR_row;
	sub->algnL_row    = algnL_row;
	sub->algnRel_row  = algnRel_row;
    sub->algnRel_buf  = algnRel_buf;
    sub->max_buf_size = max_buf_size; 
  
    uint32_t slen = (LEN_READ+15)/8; 
    sub->kswq_L[0] = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (5 + 4)); // a single block of memory
    sub->kswq_L[1] = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (5 + 4)); // a single block of memory
    sub->kswq_R[0] = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (5 + 4)); // a single block of memory
    sub->kswq_R[1] = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (5 + 4)); // a single block of memory

    for(i = 0; i < 5; ++i){
        
        for(j = 0; j < 5; ++j){
            if(i == j) sub->mat[i*5+j] = 1;
            else sub->mat[i*5+j] = -1;
        }
    }
  
 


  
        *out_sub = sub; 
	
   
    
    
    
    
    
    return sub;
}
void setBlckData(struct ExtBlck *in_blck, uint8_t cls, uint32_t off_cap,struct SubBuf *sub, struct ExtBlck **out_this){
    struct ExtBlck *cB = *out_this;   
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
    //以下代码段清理sub工作空间 
     //以下代码是sub初始化	
   
    int i, j, r, r_R;

fprintf(stderr, "%u, ", __LINE__);
fprintf(stderr, "num_seqL, num_seqR, num_relat = %u, %u, %u\n", cB->num_seqL, cB->num_seqR, cB->num_relat);  


    for(i = 0; i < 8; ++i){
        uint8_t ch_buf = sub->ext_seqL[i]; 
        sub->ext_seqL[i] = sub->ext_seqL[15-i];
        sub->ext_seqL[15-i] = ch_buf;
    }
    for(i = 0; i < 10; ++i){
        sub->st_posL[i] = 0;
        sub->st_posR[i] = 0;
    }
    if(cB->num_relat <= 0xFF){
        uint8_t *R2rel = (uint8_t *)(cB->head_relat+cB->relat) +cB->num_relat+cB->num_seqL+1;
       
        for(i =1; i <=sub->seqR_out[0][0]; ++i){
            r_R = sub->seqR_out[i][0];
            for(j = R2rel[r_R]; j < R2rel[r_R+1]; j++)
                sub->algnRel_buf[j/8] = 0;
        }
    } else if(cB->num_relat <= 0xFFFF){
        uint16_t *R2rel = (uint16_t *)(cB->head_relat+cB->relat) +cB->num_relat+cB->num_seqL+1;
        for(i =1; i <=sub->seqR_out[0][0]; ++i){
            r_R = sub->seqR_out[i][0];
            for(j = R2rel[r_R]; j < R2rel[r_R+1]; j++)
                sub->algnRel_buf[j/8] = 0;
        }
    } else{
        uint32_t *R2rel = cB->head_relat + cB->relat +cB->num_relat+cB->num_seqL+1;
        for(i =1; i <=sub->seqR_out[0][0]; ++i){
            r_R = sub->seqR_out[i][0];
            for(j = R2rel[r_R]; j < R2rel[r_R+1]; j++)
                sub->algnRel_buf[j/8] = 0;
        }

    }
    
    
    for(i=1; i<=sub->seqL_out[0][0]; i++){
        r = sub->seqL_out[i][0];
        sub->algnL_row[r/8] = 0; 
    }
    
    for(i=1; i<=sub->seqR_out[0][0]; i++){
        r = sub->seqR_out[i][0];
        sub->algnR_row[r/8] = 0; 
    }

    for(i=1; i<=sub->seqRel_out[0]; i++){
        r = sub->seqRel_out[i];
        sub->algnRel_row[r/8] = 0; 
    }
  
    for(i=1; i<=sub->seqRel_R[0]; i++){
        r = sub->seqRel_R[i];
        sub->algnRel_buf[r/8] = 0; 
    }
    for(i=1; i<=sub->seqRel_L[0]; i++){
        r = sub->seqRel_L[i];
        sub->algnRel_buf[r/8] = 0; 
    }

    for(i =0; i <3; ++i){
        sub->seqL_out[0][i] = 0;             
        sub->seqR_out[0][i] = 0;
    }

    sub->seqRel_out[0] = 0;
    sub->seqRel_L[0] = 0;
    sub->seqRel_R[0] = 0;
    sub->seqL_aln_old = 1;
    sub->seqR_aln_old = 1;
    sub->seqRel_aln_old = 1;
    sub->seqRel_L_old = 1;
    sub->seqRel_R_old = 1;

    sub->pair_out[0][0] = 0;
    sub->pair_out[0][1] = 0;
    sub->pair_out[0][2] = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++
//以下代码是临时检测代码

for(i =1; i <=cB->num_seqL; ++i){
    if((sub->algnL_row[i/8]&(1<<(7-i%8))) > 0) {
        printf("%s, %u, algnL_row[%u/8] init error!!\n", __func__,__LINE__, i); 
        exit(1);
    }
}
fprintf(stderr, "%u\n", __LINE__);
  
for(i =1; i <=cB->num_seqR; ++i){
    if((sub->algnR_row[i/8] &(1<<(7-i%8))) > 0) {
        printf("%s, %u, algnR_row[%u/8] init error!!\n", __func__,__LINE__, i); 
        exit(1);
    }
}

fprintf(stderr, "%u\n", __LINE__);
for(i =1; i <=cB->num_relat; ++i){
    if((sub->algnRel_row[i/8] &(1<<(7-i%8)))>0) {
        printf("%s, %u, algnRel_row[%u/8] init error!!\n", __func__,__LINE__, i); 
        exit(1);
    }
}
fprintf(stderr, "%u\n", __LINE__);
for(i =1; i <=cB->num_relat; ++i){
    if((sub->algnRel_buf[i/8] &(1<<(7-i%8)))>0) {
        printf("%s, %u, algnRel_buf[%u/8] init error!!\n", __func__,__LINE__, i); 
        exit(1);
    }
}

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
    //以下代码段更新cB数据块
    
	cB = in_blck + cls;
    struct CapIfo  *cur_cap;
	cur_cap = cB->head_cap + off_cap;
    cB->relat     = cur_cap->relat;	
	cB->nxtpnt    = cur_cap->nxtpnt;
	//cB->nxtflg    = cur_cap->nxtflg;
	cB->smbwt     = cur_cap->smbwt;
	cB->num_seqL  = cur_cap->num_seqL;
	cB->num_seqR  = cur_cap->num_seqR;
	cB->num_relat = cur_cap->num_relat;
    //cB->extidx    = off_cap;
    cB->extidx = cB->nxtpnt;
fprintf(stderr, "%u, off_cap = %u, ", __LINE__, off_cap);
fprintf(stderr, "num_seqL, num_seqR, num_relat = %u, %u, %u\n", cB->num_seqL, cB->num_seqR, cB->num_relat);  
fprintf(stderr, "cap: relat, nxtpnt, smbwt = %u, %u, %u\n", cur_cap->relat, cur_cap->nxtpnt, cur_cap->smbwt);  
	cB->bwtL  =  cB->smbwt ;
	uint32_t len = get_bwt_size(cB->num_seqL);
    cB->bwtR  =  cB->smbwt + len;	
    *out_this = cB; 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    return;
}

void initStackTree(struct StackTree *sTree){
	sTree->cls = 0;
    sTree->stck_arry[0][0] = 0; //其中buf_align[3]是capidx数组的行号;
	sTree->stck_arry[0][1] = 0; 
	sTree->stck_arry[0][2] = 0;
    sTree->len_arry = 0; 
	sTree->len_buf = 0;
    return ;
}


void setStackTree(struct StackTree *sTree,uint32_t (*pair_out)[3], int seed_id, uint32_t (*ext_idx)[2]){
	uint32_t bgn, end, i, j, k; 
    i = sTree->len_arry;
    k = sTree->len_buf; 
    for (j = pair_out[0][1]+1; j <= pair_out[0][0]; ++j, ++i, ++k){
	    sTree->stck_arry[i][0] =  pair_out[j][0];
		sTree->stck_arry[i][1] =  pair_out[j][1];
		sTree->stck_arry[i][2] =  pair_out[j][2];
		sTree->stck_arry[i][3] =  sTree->cls+1;
   
        sTree->back_buf[k][0] =  pair_out[j][0];
		sTree->back_buf[k][1] =  pair_out[j][1];
		sTree->back_buf[k][2] =  pair_out[j][2];
		sTree->back_buf[k][3] =  sTree->cls+1;
		sTree->back_buf[k][4] =  seed_id;
        if(pair_out[j][1] == 0xFF){
    	    sTree->back_buf[k][5] =  ext_idx[pair_out[j][1]][0];
		    sTree->back_buf[k][6] =  ext_idx[pair_out[j][1]][1];
        }

	}
    sTree->len_arry = i;
    sTree->len_buf = k;
//sTree->len_buf = 0;

    pair_out[0][0] = 0;
    pair_out[0][1] = 0;
    pair_out[0][2] = 0;
    return;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	uint8_t  *read_seq;
	uint8_t	 *seed_seq;

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
	uint32_t *cap = jmp->cap;
	uint32_t *jmp_idx = jmp->jmp;
    uint8_t *mod_idx = jmp->mod;

    uint32_t bgn = buf_algn[1];	
    uint32_t bgn_q = jmp_idx[bgn/256];
	uint32_t end_q = jmp_idx[(bgn/256)+1];
	uint32_t i_r = bgn%256, j_idx;
    for(j_idx=bgn_q; j_idx<end_q; j_idx++){
        
//fprintf(stderr, "mod_idx = %u\n", mod_idx[j_idx]);
        if(i_r == mod_idx[j_idx]){
            break;
        } 
    }
    if(j_idx == end_q){
//fprintf(stderr, "idx = 844453254, num = 315, row = 32174737\n");
//fprintf(stderr, "mod_idx = %u\n", mod_idx[32174737]);
        fprintf(stderr, "%u, %s\n", __LINE__, __func__);
        fprintf(stderr, "j_idx= %u, bgn_q = %u, end_q = %u\n", j_idx, bgn_q, end_q);        
        fprintf(stderr,"bg_idx= %u, i_r = %u\n", bgn, i_r);
        exit(1); 
    }
    
    buf_algn[3] = cap[j_idx];
	return;	
}

//int AlgnPos(struct PubParm *pParm, uint32_t pos, uint8_t *read_seq,uint32_t (*algn_out)[5])

int AlgnPos(idx_t *fm_idx, uint8_t ext_cls, uint8_t *read_seq,int read_len, int seed_off, uint32_t pos_buf[], int pos_num, struct SubBuf *sub){

   	int flg =0 ;
    int n_aln = 0;
    uint32_t *smpos;
    uint32_t i, j, j_idx, pos_s, smpos_i;
    uint8_t *target = sub->target;
    uint32_t *aln_pos = sub->aln_pos; 
    uint8_t (*aln_r)[2] = sub->aln_r;//第i个pos右段比对结果，[i][0]是比对上1否则0, [i][1]比对分数
    
    //kswq_t **swq;
    //swq = sub->swq;
    kswq_t *swq_L[2] = {0, 0};
    kswq_t *swq_R[2] = {0, 0};
    //swq[0] = 0; swq[1] = 0;
    int8_t *mat= sub->mat;       
//++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //[L_st, L_end) [R_st, R_end)
    uint32_t L_st = 0;
    uint32_t L_end = (read_len/2+seed_off)-(ext_cls*16);//其中10是种子左端点中心的偏移量
    uint32_t ext_len = SEED_LEN+ext_cls*16*2;
    uint32_t R_st = L_end+ext_len;
    uint32_t R_end = read_len;
    uint32_t R_pos_bgn, R_pos_end;
    uint32_t L_pos_bgn, L_pos_end;
   //++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //以下比对右段序列
    uint8_t *aln_seq; int len_seq; 
    uint32_t pos_i, pos;
    swq_R[0] = 0;
    aln_seq = read_seq+R_st;
    len_seq = R_end - R_st;
    for(pos_i = 0; pos_i < pos_num; ++pos_i){
        pos = pos_buf[pos_i];  
        R_pos_bgn = pos+ext_len;
        R_pos_end = R_pos_bgn+len_seq;
        for(i =R_pos_bgn; i < R_pos_end; ++i) target[i-R_pos_bgn] = __get_pac(fm_idx->pac, i);
/*  
for(i = 0; i < len_seq; ++i){
    fprintf(stderr, "%u\t%u\n", aln_seq[i], target[i]);
}
*/
        kswr_t r = ksw_align_R(len_seq, aln_seq, len_seq, target, 5, mat, 5, 2, KSW_XSTART, &swq_R[0], sub->kswq_R);
        int score_r = (2*len_seq-r.qe-r.te-2)/2+(r.qb+r.tb);
        int score = (len_seq - r.score + score_r)*100/len_seq;
//fprintf(stderr, "score = %u, ref_bg = %u, ref_end = %u\n", r.score, R_pos_bgn+r.tb, R_pos_bgn+r.te);
//fprintf(stderr, "%u, score_r = %d, score = %d\n", __LINE__, score_r, score);
        if(score < ALN_THRES_SCORE){ 
            aln_r[n_aln][0] = 0; 
            aln_r[n_aln][1] = score;
            aln_pos[n_aln] = pos;
            n_aln++;
//fprintf(stderr, "%u, n_aln = %u, pos = %u\n", __LINE__, n_aln, pos);
        } else{
        }
    }//end for(pos_i = 0; pos_i < pos_num; ++pos_i)++++++
    //===============================================
    //以下代码段是左侧比对
    aln_seq = read_seq;
    len_seq = L_end - L_st;
    uint8_t aln_num = 0;
    swq_L[0] = 0;
    for(j = 0; j < n_aln; ++j){
        pos = aln_pos[j];
        if(pos < L_end) continue;
        L_pos_bgn = pos-L_end;
        L_pos_end = pos;
        for(i = L_pos_bgn; i < L_pos_end; ++i) target[i-L_pos_bgn] = __get_pac(fm_idx->pac, i);
        kswr_t r = ksw_align_L(len_seq, aln_seq, len_seq, target, 5, mat, 5, 2, KSW_XSTART, &swq_L[0], sub->kswq_L);
        int score_r = (2*len_seq-r.qe-r.te-2)/2+(r.qb+r.tb);
        int score = (len_seq - r.score + score_r)*100/len_seq;

//fprintf(stderr, "score = %u, ref_bg = %u, ref_end = %u\n", r.score, L_pos_bgn+r.tb, L_pos_bgn+r.te);
//fprintf(stderr, "%u, score_r = %d, score = %d\n", __LINE__, score_r, score);
        if(score< ALN_THRES_SCORE){
            aln_r[aln_num][0] = 0;
            aln_r[aln_num][1] = aln_r[j][1]+score;
            aln_pos[aln_num] = L_pos_bgn;
            aln_num++; 
//fprintf(stderr, "%u, n_aln = %u, pos = %u\n", __LINE__, aln_num, pos);
        } else{ }
//fprintf(stderr, "L_pos_end = %u, score=%u, r.score=%u, l_target = %u, st, se, qt, qe = %u, %u,%u,%u\n", L_pos_end, score, r.score, len_seq, r.tb, r.te, r.qb, r.qe); 
    }//end  for(j = 0; j < n_aln; ++j)+++++++++++++
    flg = aln_num;
//fprintf(stderr, "%u, aln_num = %u\n", __LINE__, aln_num);
    return flg;
}
int init_seed_model(seed_t *seed )
{
    seed->MAX_SIZE = 16;
    int SEED_SLC_SIZE = seed->MAX_SIZE;
    seed->slc   = calloc(SEED_SLC_SIZE, sizeof(seed_slc_t)); 
    seed->slc_i = calloc(SEED_SLC_SIZE, sizeof(seed_slc_t)); 
    seed->slc_b = calloc(SEED_SLC_SIZE, sizeof(seed_slc_t)); 
    seed->aln_8 = calloc(8*3, sizeof(aln_mer8_t));
    int a[16][5] =  {  //drct,    s_off,            h_off,   m8_off,  flag;
                    {    0,       SEED_LEN/2,       8,       0,       0    }, //0
                    {    1,       SEED_LEN/2,       8,       0,       0    }, //1
                    {    0,       SEED_LEN/2-8,     8,       0,       0    }, //2
                    {    1,       SEED_LEN/2-8,     8,       0,       0    }, //3
                    {    0,       SEED_LEN/2+8,     8,       0,       0    }, //4
                    {    1,       SEED_LEN/2+8,     8,       0,       0    }, //5
                    {    0,       SEED_LEN/2,       8,       0,       1    }, //6
                    {    1,       SEED_LEN/2,       8,       0,       1    }, //7
                    {    0,       SEED_LEN/2-8,     8,       0,       1    }, //8
                    {    1,       SEED_LEN/2-8,     8,       0,       1    }, //9
                    {    0,       SEED_LEN/2,       0,       12,      2    }, //10
                    {    1,       SEED_LEN/2,       0,       12,      2    }, //11
                    {    0,       SEED_LEN/2+8,     8,       0,       1    }, //12
                    {    1,       SEED_LEN/2+8,     8,       0,       1    }, //13
                    {    0,       SEED_LEN/2+4,     8,       0,       1    }, //14
                    {    1,       SEED_LEN/2+4,     8,       0,       1    }  //15
                    };
    int i;
    for(i = 0; i < 16; ++i) {
        seed->slc_i[i].drct   = a[i][0];
        seed->slc_i[i].s_off  = a[i][1];
        seed->slc_i[i].h_off  = a[i][2];
        seed->slc_i[i].m8_off = a[i][3];
        seed->slc_i[i].flag   = a[i][4]; 
    }
    return;
}

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
/*  
                for(i= 0; i < read_len; ++i){ 
                    r_read_seq[read_len-i-1] = 3-f_read_seq[i];
                    r_qual[read_len-i-1] = qual[i];
                }
*/
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
/*
for(i = 0; i < find_num; ++i){
    fprintf(stderr, "err_pos = %u,alt=%u,  bgn=%u, end = %u\n", seed->aln_buf[i][0], seed->aln_buf[i][1], seed->aln_buf[i][2],                     seed->aln_buf[i][3]); 
}
*/          

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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
uint32_t aln_seed_seq_1(idx_t *fm_idx,uint32_t hash_boundry[], query_t *query, seed_t *seed )
{
    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq;
    uint8_t *qual = query->qual;
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
    int seq_off = read_len/2-slc_i[seed_id].s_off;

    int seed_flg = seed->seed_flg;
    uint32_t bgn, end, num = 0, k, l;


    int MIN_Q = 0;
    if(seed_id == 0) {
         //1）获取read信息,获取read反向互补序列；
/*  
        for(i= 0; i < read_len; ++i){ 
            r_read_seq[read_len-i-1] = 3-f_read_seq[i];
            r_qual[read_len-i-1] = qual[i];
        }
*/
        //2）生成正向序列与反向序列上的种子序列；
        int row_b = 0, row = 0;
fprintf(stderr, "%u\n", __LINE__);
        for(i = 0; i < seed->MAX_SIZE; ++i){
            uint32_t seed_bg, seed_ed;
            seed_bg = read_len/2-(slc_i[i].s_off-slc_i[i].h_off);
            seed_ed = seed_bg+12;
            int sum = 0;

fprintf(stderr, "%u\n", __LINE__);
            if(slc_i[i].drct == 0) {

fprintf(stderr, "%u\n", __LINE__);
                for(j = seed_bg; j < seed_ed; ++j) {
                     if( qual[j] < MIN_Q ){ 
                         //sum++; 
                     }
                 }
            } else{

fprintf(stderr, "%u\n", __LINE__);
                for(j = seed_bg; j < seed_ed; ++j) {
                    if( r_qual[j] < MIN_Q ){
                         //sum++; 
                    }
                }
            }

fprintf(stderr, "%u\n", __LINE__);
            if(sum > 2) { continue;}
            else if(sum > 0){ 

fprintf(stderr, "%u\n", __LINE__);
                memcpy(slc_b+row_b, slc_i+i, sizeof(seed_slc_t));
                row_b++;
            } else{
                memcpy(slc+row, slc_i+i, sizeof(seed_slc_t));
                row++;
            }

fprintf(stderr, "%u\n", __LINE__);
        } 

fprintf(stderr, "%u\n", __LINE__);
        memcpy(slc+row, slc_b, row_b*sizeof(seed_slc_t));
        //3）赋值当前read_id的seed_slc_size;
        seed_slc_size = row + row_b;
        seed->slc_size = seed_slc_size;
        find_num = 0; 
    } //end if(seed_id == 0)  
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //种子序列指针赋值

fprintf(stderr, "%u\n", __LINE__);
    if(find_num==0) {
        if(slc[seed_id].drct == 0){
            pseed = f_read_seq+seq_off;
        } else{//slc[seed_id].drct ==1
            pseed = r_read_seq+seq_off;
        }     
    }             

fprintf(stderr, "%u\n", __LINE__);
    while(seed_id < seed_slc_size){ 

fprintf(stderr, "%u\n", __LINE__);
        if(find_num==0) {
            seq_off = read_len/2-slc_i[seed_id].s_off;
            if(slc[seed_id].drct == 0){
                pseed = f_read_seq+seq_off;
            } else{//slc[seed_id].drct ==1
                pseed = r_read_seq+seq_off;
            }
        }   

fprintf(stderr, "%u\n", __LINE__);
        if(find_num == 0) {//选定新的种子序列  

fprintf(stderr, "%u\n", __LINE__);
            if(slc[seed_id].flag == 0) {//通过BWT精确比对选定种子序列

fprintf(stderr, "%u\n", __LINE__);
                seq12 = lkt_seq2LktItem(pseed, 8, 19);
fprintf(stderr, "seq12 = %u, seq_off = %u\n", seq12, seq_off);
                k = fm_idx->fastmap->item[seq12];
                l = fm_idx->fastmap->item[seq12+1]-1;

fprintf(stderr, "%u\n", __LINE__);
                l -= get_12mer_correct(hash_boundry, l);

fprintf(stderr, "%u\n", __LINE__);
                num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, pseed, &k, &l);

fprintf(stderr, "%u\n", __LINE__);
                bgn = k;
                end = l;
                seed_flg = 0;
            } else if(slc[seed_id].flag ==1){//通过BWT选定近似种子序列集合
                cur_row = 0;

fprintf(stderr, "%u\n", __LINE__);
                seed->hash_12mer = lkt_seq2LktItem(pseed, 8, 19);

fprintf(stderr, "%u\n", __LINE__);
                seed->seed_8mer_L = pseed;
                seed->seed_seq = pseed;

fprintf(stderr, "%u\n", __LINE__);
                find_num = get_seed_8mer_bwt_L(fm_idx, hash_boundry, seed);

fprintf(stderr, "%u\n", __LINE__);
                if(find_num == 0) {
                    seed_id++;
                    continue;
                }
/*
for(i = 0; i < find_num; ++i){
fprintf(stderr, "err_pos = %u,alt=%u,  bgn=%u, end = %u\n", seed->aln_buf[i][0], seed->aln_buf[i][1], seed->aln_buf[i][2],                     seed->aln_buf[i][3]); 
}
*/          

                bgn = seed->aln_buf[cur_row][2]; 
                end = seed->aln_buf[cur_row][3]; 
                num = end+1- bgn;
                seed_flg = 1;
                if(find_num == 1) { 
                    find_num = 0; 
                }

fprintf(stderr, "%u\n", __LINE__);
            } else{//slc[seed_id] ==2
                //调用右侧模糊比对查找函数,获得模糊比对序列集合数组
                int bg_i = 0, ed_i = 8, num_seq;
                num_seq = (ed_i-bg_i)*3; 

fprintf(stderr, "%u\n", __LINE__);
                seed->hash_12mer = lkt_seq2LktItem(pseed, 0, 11);
                seed->seed_8mer_R = pseed+12;
                seed->seed_seq = pseed;
           
fprintf(stderr, "%u\n", __LINE__);
                int var_seq_num = gen_8mer_var_buf(seed, bg_i, ed_i, 1); 

fprintf(stderr, "%u\n", __LINE__);
                find_num = get_seed_8mer_seq_R(fm_idx, seed, num_seq);

fprintf(stderr, "%u\n", __LINE__);
                if(find_num == 0) {
                    seed_id++;
                    continue;
                }

fprintf(stderr, "%u\n", __LINE__);
                memcpy(seed_buf, pseed, SEED_LEN); 
                cur_row = 0;
                uint8_t f_i = seed->find_buf_R[cur_row];  
                uint8_t var_pos  = seed->var_info[f_i][0];
                uint8_t var_ch   = seed->var_info[f_i][1];
                uint8_t var_ch_o = seed_buf[var_pos+12]; 
                seed_buf[var_pos+12] = var_ch; 

fprintf(stderr, "%u\n", __LINE__);
                seq12 = lkt_seq2LktItem(seed_buf, 8, 19);
                k = fm_idx->fastmap->item[seq12];
                l = fm_idx->fastmap->item[seq12+1]-1;

fprintf(stderr, "%u\n", __LINE__);
                l -= get_12mer_correct(hash_boundry, l);

fprintf(stderr, "%u\n", __LINE__);
                num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, seed_buf, &k, &l);

fprintf(stderr, "%u\n", __LINE__);
                bgn = k; end = l;
                seed_flg = 2;
                if(find_num == 1) { 
                    find_num = 0; 
                }
                seed_buf[var_pos+12] = var_ch_o; 
            }
        } else{//if(find_num >0)

fprintf(stderr, "%u\n", __LINE__);
            if(seed_flg ==1) {//第一类近似比对，即BWT近似比对, 种子序列继续
                seed_id--;
                cur_row++; 

fprintf(stderr, "%u\n", __LINE__);
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

fprintf(stderr, "%u\n", __LINE__);
                seq12 = lkt_seq2LktItem(seed_buf, 8, 19);
                k = fm_idx->fastmap->item[seq12];
                l = fm_idx->fastmap->item[seq12+1]-1;

fprintf(stderr, "%u\n", __LINE__);
                l -= get_12mer_correct(hash_boundry, l);

fprintf(stderr, "%u\n", __LINE__);
                num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, seed_buf, &k, &l);
                bgn = k; 
                end = l;

fprintf(stderr, "%u\n", __LINE__);
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

fprintf(stderr, "%u\n", __LINE__);
        if(num >0){//种子序列比对失败
            break; 
        }

fprintf(stderr, "%u\n", __LINE__);
        seed->id = ++seed_id;
    }
    
fprintf(stderr, "%u\n", __LINE__);
    seed->id = seed_id;
    return num;
}

int AlgnPos_buf(idx_t *fm_idx, uint8_t ext_cls, uint8_t *read_seq,int read_len, int seed_off, int num_pair, uint32_t pos_buf[], struct SubBuf *sub)
{
    uint32_t nxt_pnt, nxt_end, pos;
    uint8_t  nxt_flg;
    uint32_t i, j;
    int Flg_Algn = 0;
    for (i = 1; i <= num_pair; ++i){
        nxt_pnt = sub->pair_out[i][0];
        nxt_flg = sub->pair_out[i][1];
        nxt_end = nxt_pnt + nxt_flg;
        uint32_t pos_i = 0;
        if(nxt_flg>1){
            for(j=nxt_pnt; j<nxt_end; j++){
                pos = bwt_sa(fm_idx->bwt, j);
                pos_buf[pos_i++] = pos;
            }
        } else{
            pos_buf[pos_i++] = nxt_pnt; 
        }
for(j = 0; j < pos_i; ++j) {
    fprintf(stderr, "pos_buf: ");
    fprintf(stderr, "%u  ", pos_buf[j]);
    fprintf(stderr, "\n");
}
        //利用SW方法比对进行比对，返回成功与否标记，成功返回比对位置个数，失败返回0；
        Flg_Algn = AlgnPos(fm_idx, ext_cls+1, read_seq, read_len, seed_off, pos_buf, pos_i, sub);  
        if(Flg_Algn>0){
            OutAlgnInfo(sub->algn_out);
            break;
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
    uint32_t (*alg_out)[3];
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
    if(sub->seqR_out[0][0] == 2 && sub->seqL_out[0][0] == 2) fprintf(stderr, "dataNum = %u\n", n);
    return alg_out[0][0];

}

int Algn_apro_all(struct ExtBlck *eB, struct SubBuf *sub,  int flg){
	uint8_t  *smbwt;
	uint8_t  *sumbuf;
	uint8_t  *seq;
    uint8_t  *st_pos;	
	uint32_t (*alg_out)[3];
	uint32_t n = 0;
	uint32_t num_seqL  = eB->num_seqL;
	uint32_t num_seqR  = eB->num_seqR;
    uint8_t *algn_row;
    
    if(flg<1){
		smbwt = eB->head_smbwt + eB->bwtL;
		n = num_seqL;
		alg_out = sub->seqL_out;
		seq = sub->ext_seqL;
        st_pos = sub->st_posL;
        sub->seqL_aln_old = alg_out[0][0];
        algn_row = sub->algnL_row;
    }else{
		smbwt = eB->head_smbwt + eB->bwtR;
		n = num_seqR;
		alg_out = sub->seqR_out;
		seq = sub->ext_seqR;
        st_pos = sub->st_posR;
        sub->seqR_aln_old = alg_out[0][0];
        algn_row = sub->algnR_row;
    }
fprintf(stderr, "%u, n_data = %u, flag = %d\n", __LINE__, n, flg);
    
    st_pos[9] = 1;
    uint8_t *cnt = smbwt+((n+1)/2)*8;

    if(n<=MIN_BWT_SIZE){
        align_min(smbwt, n,seq, 2, 1,st_pos, algn_row, alg_out); 
        //return alg_out[0][0];
    } 
    int i;
   
    if(n<=NO_BWT_SUM_SIZE){ 
        uint8_t cnt1[8][17];
        gen_cnt(smbwt, n, cnt1);
        for(i = 0; i < 8; ++i){
//if(st_pos[i] > 0) continue;

fprintf(stderr, "%u\n", __LINE__);        
            align_nosum(smbwt, cnt1, n, seq, 8, i,st_pos, algn_row, alg_out); 

fprintf(stderr, "%u\n", __LINE__);        
        }
        //return alg_out[0][0];
    } else if(n <= 255) {
        for(i = 0; i < 8; ++i){

            

//if(st_pos[i] > 0) continue;

fprintf(stderr, "%u\n", __LINE__);        
            align_255(smbwt, (uint8_t (*)[17])cnt, n, seq, 8, i, st_pos, algn_row, alg_out);     

fprintf(stderr, "%u\n", __LINE__);        
        }
        //return alg_out[0][0];
    } else if(n > 255){//>=256
        for(i = 0; i < 8; ++i){

//if(st_pos[i] > 0) continue;
            
fprintf(stderr, "%u\n", __LINE__);        
            align_large_alt(smbwt, cnt, n, seq, 8, i, st_pos, algn_row, alg_out);    
fprintf(stderr, "%u\n", __LINE__);        
        }
fprintf(stderr, "%u, %s, %u\n", __LINE__, __func__, alg_out[0][0]); 
        //return alg_out[0][0];
    }
    
//if(sub->seqR_out[0][0] == 2 && sub->seqL_out[0][0] == 2) fprintf(stderr, "dataNum = %u\n", n);
    
    fprintf(stderr, "%u, %d alg_out[0][0] =%u, n_data = %u\n", __LINE__, flg, alg_out[0][0], n);
    
    return alg_out[0][0];

}


int AlgnmExtSeq(struct ExtBlck *eB, struct SubBuf *sub, int flg){
	uint8_t  *smbwt;
	uint8_t  *sumbuf;
	uint8_t  *seq;
	uint8_t  *st_pos;
	uint32_t (*alg_out)[3];
	uint32_t n = 0;
	uint32_t num_seqL  = eB->num_seqL;
	uint32_t num_seqR  = eB->num_seqR;
    uint8_t *algn_row;

    if(flg<1){
		smbwt = eB->head_smbwt + eB->bwtL;
		n = num_seqL;
		alg_out = sub->seqL_out;
		seq = sub->ext_seqL;
	    st_pos = sub->st_posL;
        algn_row = sub->algnL_row;
    }else{
		smbwt = eB->head_smbwt + eB->bwtR;
		n = num_seqR;
		alg_out = sub->seqR_out;
		seq = sub->ext_seqR;
	    st_pos = sub->st_posR;
        algn_row= sub->algnR_row;
    }

    uint8_t *cnt = smbwt+((n+1)/2)*8;
fprintf(stderr, "%u n_data = %u\n", __LINE__, n);
    if(n<=MIN_BWT_SIZE){

        align_min(smbwt, n,seq, 1, 1,  st_pos,algn_row, alg_out); 
    } else if(n<=NO_BWT_SUM_SIZE){ 
        uint8_t cnt1[8][17];
        gen_cnt(smbwt, n, cnt1);
        align_nosum(smbwt, cnt1, n, seq, 8, 7,st_pos,algn_row, alg_out); 
    } else if(n <= 255) {
        align_255(smbwt, (uint8_t (*)[17])cnt, n, seq, 8, 7, st_pos, algn_row, alg_out);     
    } else{//>=256
        align_large_alt(smbwt, cnt, n, seq, 8, 7, st_pos, algn_row, alg_out);
    }
    if(sub->seqR_out[0][0] == 2 && sub->seqL_out[0][0] == 2) fprintf(stderr, "dataNum = %u\n", n);
    int i;

    fprintf(stderr, "%u, %d alg_out[0][0] =%u, n_data = %u\n", __LINE__, flg, alg_out[0][0], n);
    for(i = 1; i <=alg_out[0][0]; ++i ){
        int r = alg_out[i][0];
        int k = algn_row[r/8]&(1<<(7-r%8));
        fprintf(stderr, "seq_out[%u] = %u, algn_row = %u\n", i, alg_out[i][0], k);
    }
    for(i = 1; i <=n; ++i ){
       
        int k = algn_row[i/8]&(1<<(7-i%8));
        if(k >0){ 
            fprintf(stderr, "seq_row = %u\n", i); 
        }
    }
    
    return alg_out[0][0];

}


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
        fprintf(stderr, "find cur_seq i = %u\n", i);
    } else{
    
        fprintf(stderr, "can't find cur_seq i = %u\n", i);

        fprintf(stderr, "Seq:\t");
        for(j = 0; j < 16; ++j){
            fprintf(stderr, "%u", seq[j]);
        }
        fprintf(stderr, "\n");
        for(i = 0; i < num; ++i){
            fprintf(stderr, "%u\t", i);
            for(j = 0; j < 16; ++j){
                fprintf(stderr, "%u", seq_buf[i*16+j]);
            }
            fprintf(stderr, "\n");
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
            //fprintf(stderr, "[%u, %u]: %u %u\n", i, j, seq_buf[i*16+j], seq_buf0[i*16+j]);
            if(seq_buf[i*16+j] != seq_buf0[i*16+j]){
                fprintf(stderr, "test_smbwt error!!!!\n"); 
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
fprintf(stderr, "seqL[0][1] = %u, seqR[0][1] = %u\n", seqL[0][1], seqR[0][1]);
fprintf(stderr, "%s, 1129, numL, numR = %u, %u\n", __func__, numL, numR);
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
fprintf(stderr, "i, r_R = %u %u\n", i, r_R);
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
//fprintf(stderr, "i, r_L = %u %u\n", i, r_L);

//fprintf(stderr, "i = %u, L2rel[%u] = %u, %u\n", i, r_L, L2rel[r_L],L2rel[r_L+1]);
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
//fprintf(stderr, "pair = %u, %u, %u\n", pair[p_num][0], pair[p_num][1], pair[p_num][2]); 
//fprintf(stderr, "%u, p_num = %u, j = %u, r_L = %u, r_R= %u\n", __LINE__, p_num, j, r_L, r_R);   
    	        algnRel_row[j/8] |= (1<<(7-j%8));
                seqRel[0]++;
                seqRel[seqRel[0]] = j;

            }
		}  
    }
 
//fprintf(stderr, "%u, %s p_num=%u, numR = %u, numL= %u\n", __LINE__, __func__, p_num, numR, numL);
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

int PairExtSeq_8(struct ExtBlck *eB, struct SubBuf *sub, uint8_t *relat){

	uint32_t (*seqL)[3];
	uint32_t (*seqR)[3];
	uint32_t *seqRel, *seqRel_L, *seqRel_R;
    uint32_t (*pair)[3];
    uint8_t *algnRel_row, *algnRel_buf;	
	uint8_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
    seqRel_L = sub->seqRel_L;
    seqRel_R = sub->seqRel_R;
    pair = sub->pair_out;
    algnRel_row = sub->algnRel_row;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + eB->num_relat; 
    R2rel = L2rel + eB->num_seqL+1; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int i, j, k; 
    int r_L, r_R, r_Rel;
	int p_num;
    
    int numL = seqL[0][0]; 
	int numR = seqR[0][0];
    
    int st_i = sub->seqR_aln_old;
fprintf(stderr, "%u, st_i = %u, numR = %u, numL = %u\n", __LINE__, st_i, numR, numL); 
    for(i=st_i; i<=numR; i++){
        r_R = seqR[i][0];
fprintf(stderr, "seqR[%u] -= %u, %u, %u\n", i, seqR[i][0],  seqR[i][1], seqR[i][2]);
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++) {
            algnRel_buf[j/8] |= (1<<(7-j%8));  
            //seqRel_R[0]++;
            //seqRel_R[seqRel_R[0]] = j;
        }
    } 
    
    if(pair[0][0]>0){
        p_num = pair[0][0]-pair[0][1]; 
        memcpy(pair+1, pair+1+pair[0][1], p_num*(3*sizeof(uint32_t)));
    } else{
        p_num = 0;
    }
st_i = sub->seqL_aln_old;
    for(i=st_i; i<=numL; i++){

fprintf(stderr, "seqL[%u] -= %u, %u, %u\n", i, seqL[i][0],  seqL[i][1], seqL[i][2]);
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
    	    int flag_Rel = algnRel_row[j/8] & (1<<(7-j%8)); 
fprintf(stderr, "%u, st_i = %u, numR = %u, numL = %u, flag_Rel = %u\n", __LINE__, st_i, numR, numL, flag_Rel); 
            if(flag_Rel > 0) continue;	
            
            r_Rel = relat[j];
            int flag_R = algnRel_buf[r_Rel/8] & (1<<(7-r_Rel%8));

fprintf(stderr, "%u, r_Rel = %u, flag_R = %u\n", __LINE__, r_Rel, flag_R); 
            if(flag_R>0){
				for(k = 1; k <= numR; k++){
                    if(seqR[k][0] == r_R) break;
                }
                if(k > numR){
                    printf("%u, %s algnR_row[] error!, k =%u, numR = %u\n", __LINE__, __func__, k, numR);
                    exit(1);
                }
                p_num++;
                pair[p_num][0] = nxtpnt[r_Rel];
				pair[p_num][1] = nxtflg[r_Rel];
				pair[p_num][2] = seqL[i][1] + seqR[k][1];
    	        algnRel_row[r_Rel/8] |= (1<<(7-r_Rel%8));
                seqRel[0]++;
                seqRel[seqRel[0]] = r_Rel;
            }
		}  
    }
fprintf(stderr, "%u, p_num = %u\n", __LINE__, p_num); 
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

int PairExtSeq_16(struct ExtBlck *eB, struct SubBuf *sub, uint16_t *relat){

	uint32_t (*seqL)[3];
	uint32_t (*seqR)[3];
	uint32_t *seqRel, *seqRel_L, *seqRel_R;
    uint32_t (*pair)[3];
    uint8_t *algnRel_row, *algnRel_buf;	
	uint16_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
    seqRel_L = sub->seqRel_L;
    seqRel_R = sub->seqRel_R;
	pair = sub->pair_out;
    algnRel_row = sub->algnRel_row;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + eB->num_relat; 
    R2rel = L2rel + eB->num_seqL+1; 
	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int i, j, k; 
    int r_L, r_R, r_Rel;
	int p_num;
    
    int numL = seqL[0][0]; 
	int numR = seqR[0][0];
    
    int st_i = sub->seqR_aln_old;
fprintf(stderr, "%u, st_i = %u, numR = %u, numL = %u\n", __LINE__, st_i, numR, numL); 
    for(i=st_i; i<=numR; i++){
        r_R = seqR[i][0];
        
fprintf(stderr, "seqR[%u] -= %u, %u, %u\n", i, seqR[i][0],  seqR[i][1], seqR[i][2]);
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++) {
            algnRel_buf[j/8] |= (1<<(7-j%8));  
            //seqRel_R[0]++;
            //seqRel_R[seqRel_R[0]] = j;
        }
    } 
    if(pair[0][0]>0){
        p_num = pair[0][0]-pair[0][1]; 
        memcpy(pair+1, pair+1+pair[0][1], p_num*(3*sizeof(uint32_t)));
    } else{
        p_num = 0;
    }
st_i = sub->seqL_aln_old;
    for(i=st_i; i<=numL; i++){

fprintf(stderr, "seqL[%u] -= %u, %u, %u\n", i, seqL[i][0],  seqL[i][1], seqL[i][2]);
        r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
    	    int flag_Rel = algnRel_row[j/8] & (1<<(7-j%8)); 
            if(flag_Rel > 0) continue;	
            
            r_Rel = relat[j];
    		int flag_R = algnRel_buf[r_Rel/8] & (1<<(7-r_Rel%8));
            if(flag_R>0){
				for(k = 1; k <= numR; k++){
                    if(seqR[k][0] == r_R) break;
                }
                if(k > numR){
                    printf("%u, %s algnR_row[] error!, k =%u, numR = %u\n", __LINE__, __func__, k, numR);
                    exit(1);
                }
                p_num++;
                pair[p_num][0] = nxtpnt[r_Rel];
				pair[p_num][1] = nxtflg[r_Rel];
				pair[p_num][2] = seqL[i][1] + seqR[k][1];
    	        algnRel_row[r_Rel/8] |= (1<<(7-r_Rel%8));
                seqRel[0]++;
                seqRel[seqRel[0]] = r_Rel;
            }
		}  
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
	return (int) p_num;
}

int PairExtSeq_32(struct ExtBlck *eB, struct SubBuf *sub, uint32_t *relat){

	uint32_t (*seqL)[3];
	uint32_t (*seqR)[3];
	uint32_t *seqRel, *seqRel_L, *seqRel_R;
    uint32_t (*pair)[3];
	uint8_t *algnRel_row, *algnRel_buf;	
	uint32_t *L2rel, *R2rel;
	uint32_t *nxtpnt ; 
	uint8_t  *nxtflg ; 
	seqL = sub->seqL_out;
	seqR = sub->seqR_out;
    seqRel = sub->seqRel_out;
    seqRel_L = sub->seqRel_L;
    seqRel_R = sub->seqRel_R;
	pair = sub->pair_out;
	algnRel_row = sub->algnRel_row;
	algnRel_buf = sub->algnRel_buf;
    L2rel = relat + eB->num_relat; 
    R2rel = relat+eB->num_relat + eB->num_seqL+1;

	nxtpnt= eB->head_nxtpnt + eB->nxtpnt; 
	nxtflg= eB->head_nxtflg + eB->nxtpnt; 
	
	int i, j, k; 
    int r_L, r_R, r_Rel;
	int p_num;
    
    int numL = seqL[0][0]; 
	int numR = seqR[0][0];
/*  
{
fprintf(stderr, "%u\n:===========================================\n", __LINE__); 
fprintf(stderr, "relat_num= %u\n", eB->num_relat); 
     for(i = 1; i <= eB->num_seqL; ++i ){
        int k = sub->algnL_row[i/8]&(1<<(7-i%8));
        if(k >0){ 
            fprintf(stderr, "algnL_row[%u/8] >0\n", i); 
        }
    }
fprintf(stderr, "\n---------------------------\n"); 
}
*/

    int st_i = sub->seqR_aln_old;
fprintf(stderr, "numR = %u, st_i = %u\n", numR, st_i);
    for(i=st_i; i<=numR; i++){
        r_R = seqR[i][0];
        for(j=R2rel[r_R]; j<R2rel[r_R+1]; j++) {
            algnRel_buf[j/8] |= (1<<(7-j%8));  
        }
    } 
    if(pair[0][0]>0){
        p_num = pair[0][0]-pair[0][1]; 
        memcpy(pair+1, pair+1+pair[0][1], p_num*(3*sizeof(uint32_t)));
    } else{
        p_num = 0;
    }

st_i = sub->seqL_aln_old;
    for(i=st_i; i<=numL; i++){
		r_L = seqL[i][0];
        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
    	    int flag_Rel = algnRel_row[j/8] & (1<<(7-j%8)); 
            if(flag_Rel > 0) continue;	
            
            r_Rel = relat[j];
    		int flag_R = algnRel_buf[r_Rel/8] & (1<<(7-r_Rel%8));
            if(flag_R>0){
				for(k = 1; k <= numR; k++){
                    if(seqR[k][0] == r_R) break;
                }
                if(k > numR){
                    printf("%u, %s algnR_row[] error!, k =%u, numR = %u\n", __LINE__, __func__, k, numR);
                    exit(1);
                }
                p_num++;
                pair[p_num][0] = nxtpnt[r_Rel];
				pair[p_num][1] = nxtflg[r_Rel];
				pair[p_num][2] = seqL[i][1] + seqR[k][1];
    	        algnRel_row[r_Rel/8] |= (1<<(7-r_Rel%8));
                seqRel[0]++;
                seqRel[seqRel[0]] = r_Rel;
            }
		}  
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
    return (int) p_num;
}
int PairExtSeq_all(struct ExtBlck *eB, struct SubBuf *sub)
{
    uint32_t *relat = eB->head_relat + eB->relat; 
    int flag = 0;
    int i;

    if(eB->num_relat<=0xFF) {
        uint8_t *relat_8 = (uint8_t *)relat;
        flag = PairExtSeq_8(eB, sub, relat_8);
fprintf(stderr, "%u, %s, flag = %u\n", __LINE__, __func__, flag);
    } else if(eB->num_relat<=0xFFFF){
        uint16_t *relat_16 = (uint16_t *)relat;
        flag = PairExtSeq_16(eB, sub, relat_16);
fprintf(stderr, "%u, %s, flag = %u\n", __LINE__, __func__, flag);
    } else{
        flag = PairExtSeq_32(eB, sub, relat);
fprintf(stderr, "%u, %s, flag = %u\n", __LINE__, __func__, flag);
    }    
fprintf(stderr, "%u, %s\n", __LINE__, __func__);
    return flag;
}
 

/*
  
for(i=0; i<=p_num; ++i){
fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
}

    fprintf(stderr, "pair_num = %u\n", p_num); 
    for(i=0; i<=p_num; ++i){
        fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
    
    }
*/

//++++++++++++++++++++++++++++++++++++
    //int s_num = 1;
/*
    fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
    for(i=0; i<=p_num; ++i){
        fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
    
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
fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
for(i=0; i<=p_num; ++i){
    fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 

}
*/

/*
for(i = 0; i < eB->num_seqL; i++){
    fprintf(stderr, "L2rel[%u] = %u\n", i, L2rel[i]);
}
for(i = 0; i < eB->num_relat; ++i){
    //if(nxtpnt[i]>17484874&& nxtpnt[i]<17484874+200)
    fprintf(stderr, "relat[%u] = %u, nxtpnt = %u, nxtflg = %u\n", i, relat[i], nxtpnt[i], nxtflg[i]);
}
*/
//fprintf(stderr, "%s, 1130, seqL, seqR = %u, %u\n", __func__, seqL[1][0], seqR[1][0]);
/*
	for(i=0; i<=eB->num_seqR; i++){
      


        if(algnR_row[i] >0){
            fprintf(stderr, " algnR_row[%u] = %u\n", i, algnR_row[i]); 
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
fprintf(stderr, "seqL[0][1] = %u, seqR[0][1] = %u\n", seqL[0][1], seqR[0][1]);
fprintf(stderr, "%s, 1129, numL, numR = %u, %u\n", __func__, numL, numR);
	
   
    numR = seqR[0][0];
    numL = seqL[0][0];
    uint8_t *seq;    
    uint8_t seq_ch, ch;
    int row = 0;
fprintf(stderr, "%u, %s p_num=%u, numR = %u, numL= %u\n", __LINE__, __func__, p_num, numR, numL);
fprintf(stderr, "%u, %s p_num=%u, numR = %u, numL= %u\n", __LINE__, __func__, p_num, numR, numL);

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

fprintf(stderr, "%u, %s n_data=%u\n", __LINE__, __func__, n);          
        for(i = 1; i <= numL; ++i){
            r_L = seqL[i][0];
            for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
                r_R = relat[j];
                int flag_R = algnR_row[r_R/8] & (1<<(7-r_R%8));
                if( flag_R > 0) continue;
                int err_num = 0;
                //int k;

fprintf(stderr, "%u, %s r_R=%u\n", __LINE__, __func__, r_R);          

                    for(k = 7; k>=0; --k){
                        seq_ch = seq[2*k]<<2|seq[2*k+1];      
                        if(n<=MIN_BWT_SIZE){
                            //r_R = get_ch_min(r_R, k, smbwt, n, &ch); 
                            break;
                        } else if(n<=NO_BWT_SUM_SIZE){ 
                            r_R = get_ch_nosum(r_R, k, smbwt, cnt1, n, &ch); 
                        } else if(n <= 255) {
                            r_R = get_ch_255(r_R, k, smbwt, (uint8_t (*)[17])cnt, n, &ch);     
    fprintf(stderr, "%u, %s r_R=%u, k = %u, ch=%u\n", __LINE__, __func__, r_R, k, ch);          
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

fprintf(stderr, "%u, %s p_num=%u\n", __LINE__, __func__, p_num);          
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
fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
for(i=0; i<=p_num; ++i){
    fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 

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
fprintf(stderr, "seqL[0][1] = %u, seqR[0][1] = %u\n", seqL[0][1], seqR[0][1]);
fprintf(stderr, "%s, 1129, numL, numR = %u, %u\n", __func__, numL, numR);
	
    uint8_t *seq;    
    uint8_t seq_ch, ch;
    int row = 0;
fprintf(stderr, "%u, %s p_num=%u, numR = %u, numL= %u\n", __LINE__, __func__, p_num, numR, numL);
fprintf(stderr, "%u, %s p_num=%u, numR = %u, numL= %u\n", __LINE__, __func__, p_num, numR, numL);

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
fprintf(stderr, "%u, num_seqL = %u, num_relat = %u, num_seqR =%u\n", __LINE__, n_L, n_relat, eB->num_seqR);        
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

fprintf(stderr, "%u\n", __LINE__);        
                        if(find_flag == 0){

fprintf(stderr, "%u\n", __LINE__);        
                            for(k = 7; k>=0; --k){
                                seq_ch = seq[2*k]<<2|seq[2*k+1];      
                                if(n_L<=MIN_BWT_SIZE){
                                    //r_R = get_ch_min(r_R, k, smbwt, n, &ch); 
                                    break;
                                } else if(n_L<=NO_BWT_SUM_SIZE){ 
                                    r_L = get_ch_nosum(r_L, k, smbwt, cnt1, n_L, &ch); 
                                } else if(n_L <= 255) {
                                    r_L = get_ch_255(r_L, k, smbwt, (uint8_t (*)[17])cnt, n_L, &ch);     
            fprintf(stderr, "%u, %s r_R=%u, k = %u, ch=%u\n", __LINE__, __func__, r_R, k, ch);          
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

fprintf(stderr, "%u\n", __LINE__);        
                            if(k <0){

fprintf(stderr, "%u\n", __LINE__);        
fprintf(stderr, "p_num =%u. i =%u, j = %u\n", p_num, i, j);
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

fprintf(stderr, "%u\n", __LINE__);        
                            } else{

fprintf(stderr, "%u\n", __LINE__);        
                                break;//for(j =  L2rel[i]; j < L2rel[i+1]; ++j)++++
                            }

fprintf(stderr, "%u\n", __LINE__);        
                        } else{
fprintf(stderr, "%u\n", __LINE__);        
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

fprintf(stderr, "%u\n", __LINE__);        
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
fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
for(i=0; i<=p_num; ++i){
    fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 

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
fprintf(stderr, "seqL[0][1] = %u, seqR[0][1] = %u\n", seqL[0][1], seqR[0][1]);
fprintf(stderr, "%s, 1129, numL, numR = %u, %u\n", __func__, numL, numR);
/*
for(i = 0; i < eB->num_seqL; i++){
    fprintf(stderr, "L2rel[%u] = %u\n", i, L2rel[i]);
}
for(i = 0; i < eB->num_relat; ++i){
    //if(nxtpnt[i]>17484874&& nxtpnt[i]<17484874+200)
    fprintf(stderr, "relat[%u] = %u, nxtpnt = %u, nxtflg = %u\n", i, relat[i], nxtpnt[i], nxtflg[i]);
}
*/
//fprintf(stderr, "%s, 1130, seqL, seqR = %u, %u\n", __func__, seqL[1][0], seqR[1][0]);
/*
	for(i=0; i<=eB->num_seqR; i++){
        //algnR_row[i] = 0;


        if(algnR_row[i] >0){
            fprintf(stderr, " algnR_row[%u] = %u\n", i, algnR_row[i]); 
            exit(1);
        }	
	
	}
*/
	for(i=1; i<=numR; i++){
		r_R = seqR[i][0];
        algnR_row[r_R/8] = algnR_row[r_R/8] | (1<<(7-r_R%8));
fprintf(stderr, "i, r_R = %u %u\n", i, r_R);
    }
    p_num = 0;
	for(i=1; i<=numL; i++){
		r_L = seqL[i][0];

fprintf(stderr, "i, r_L = %u %u\n", i, r_L);
fprintf(stderr, "i = %u, L2rel[%u] = %u, %u\n", i, r_L, L2rel[r_L],L2rel[r_L+1]);

        for(j=L2rel[r_L]; j<L2rel[r_L+1]; j++){
    		r_R = relat[j];

fprintf(stderr, "r_L = %u, r_R= %u\n", r_L, r_R);
//fprintf(stderr, "algnR_row[%u]= %u\n", r_R, algnR_row[r_R]);
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

//fprintf(stderr, "pair = %u, %u, %u\n", pair[p_num][0], pair[p_num][1], pair[p_num][2]); 
    		}
		}  
    	
    }


/*
  
for(i=0; i<=p_num; ++i){
fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
}

    fprintf(stderr, "pair_num = %u\n", p_num); 
    for(i=0; i<=p_num; ++i){
        fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
    
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
    fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
    for(i=0; i<=p_num; ++i){
        fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 
    
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
fprintf(stderr, "pair_num = %u, s_num = %u\n", p_num, s_num); 
for(i=0; i<=p_num; ++i){
    fprintf(stderr, "pair[%u] = %u, %u, %u\n", i, pair[i][0], pair[i][1], pair[i][2]); 

}
*/

	return (int) p_num;
}
void align_255(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t seq[16], int seq_len, int seq_st, uint8_t *st_pos, uint8_t *algn_row, uint32_t (*aln_out)[3])
{
    int stat_flag, mid_num, top_num, bot_num, bg_row, ed_row, row, rot_st, rot_ed;
    uint32_t bg_idx, ed_idx, cur_pos, j;//bg_idx,ed_idx begin and end of 2nt-encoded bwt 
    uint8_t seq_ch, ch;
    int n_aln=0, n_exact_aln=0;
    cur_pos = (seq_st+1)%8; //alignment sequence pos [0, 8)
    int rot = (7-seq_st)%8; //bw rotation
    int cur_len = 0; // current alignment length
    uint8_t err_len = st_pos[9];     
    bg_idx = 0; ed_idx = n_data-1;//bwt interval [bg_idx, ed_idx]
    int aln_num=0, st_len = seq_len, err_num = 0;
    n_aln = aln_out[0][0]; 
    
    while(1){
        cur_pos = (cur_pos+7)%8; 
        seq_ch = seq[cur_pos*2]<<2| seq[cur_pos*2+1];
        if(bg_idx == 0 && ed_idx==n_data-1){
          //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
            if(cnt_2nt[rot][seq_ch] == cnt_2nt[rot][seq_ch+1]){
                if(err_num >=1) {break;}
                err_num++; 
                st_len = (cur_len+1)%8;
                rot=(rot+1)%8;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
//fprintf(stderr, "%u, cur_len = %u, rot = %u\n", __LINE__, cur_len, rot);
                continue;
            }
 
            
            bg_idx = cnt_2nt[rot][seq_ch]; 
            ed_idx = cnt_2nt[rot][seq_ch+1]-1;
            rot=(rot+1)%8;
            ++cur_len;
            if(cur_len == 8) {
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
//fprintf(stderr, "%u, n_aln = %u, idx = %u\n", __LINE__, n_aln, bg_idx);
                }
                return;
            } 
       
//
//fprintf(stderr, "[top_num, mid_num, bot_num] = %u\t%u\t%u\n", top_num, mid_num, bot_num);
//fprintf(stderr, "align seq %u\t%u\n",  cur_pos, seq_ch);
//fprintf(stderr, "[bg_idx, ed_idx] = %u\t%u\n",  bg_idx, ed_idx);
           
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
/*
        fprintf(stderr,"\n+++++++++++++++++++++++++++++++++++++\n");
        if(rot == 1 && n_data == 134 && n_ch == 2) {
            uint8_t __cnt[8][17];
            gen_cnt(bwt, n_data, __cnt);

            uint32_t __i;
            for(__i = rot_st*2; __i < rot_st*2+bg_idx; ++__i){
                fprintf(stderr, "%u\t", __get_bwt2(bwt, __i));
            }
            fprintf(stderr, "\n-------------------\n");
            for(__i = rot_st*2+bg_idx; __i <= rot_st*2+ed_idx; ++__i){
                fprintf(stderr, "%u\t", __get_bwt2(bwt, __i));
            }
            fprintf(stderr, "\n-------------------\n");
            for(__i = rot_st*2+ed_idx; __i <= rot_st*2+n_data; ++__i){
                fprintf(stderr, "%u\t", __get_bwt2(bwt, __i));
            }

            fprintf(stderr, "\n==========================\n");
        }

        fprintf(stderr, "[n_data, n_ch, rot, seq_ch] = %u\t%u\t%u\t%u\n",  n_data, n_ch, rot, seq_ch);
        fprintf(stderr, "[cnt+1, cnt] = %u\t%u\n",  cnt_2nt[rot][seq_ch+1], cnt_2nt[rot][seq_ch]);
        fprintf(stderr, "[top_num, mid_num, bot_num] = %u\t%u\t%u\n", top_num, mid_num, bot_num);
        fprintf(stderr, "align seq %u\t%u\n",  cur_pos, seq_ch);
        fprintf(stderr, "[bg_idx, ed_idx] = %u\t%u\n",  bg_idx, ed_idx);
        fprintf(stderr,"-------------------------------------\n");
 */       
        //for(j = rot_st; j < rot_ed; ++j){
            //fprintf(stderr, "%u\t%u\n",  bwt[j]>>4, bwt[j]&0xF);
        //}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
//fprintf(stderr, "%u, err_len = %u, mid_num = %u, cur_len = %u\n", __LINE__, err_len, mid_num, cur_len);
        if(err_len>0){
//fprintf(stderr, "%u, err_len = %u, mid_num = %u, cur_len = %u\n", __LINE__, err_len, mid_num, cur_len);
            if(mid_num < 1 && cur_len <= 6 ) { return;}
//fprintf(stderr, "%u, err_len = %u, mid_num = %u, cur_len = %u\n", __LINE__, err_len, mid_num, cur_len);
            if(mid_num >=1 && cur_len < 6){
//fprintf(stderr, "%u, err_len = %u, mid_num = %u, cur_len = %u\n", __LINE__, err_len, mid_num, cur_len);
                ++cur_len;
                if(cur_len <8-err_len){
                    bg_idx = top_num+cnt_2nt[rot][seq_ch];
                    ed_idx = bg_idx+mid_num-1;
                    rot = (rot+1)%8;
                    continue;
                }
            }
        }

//fprintf(stderr, "%u, err_len = %u, mid_num = %u, cur_len = %u\n", __LINE__, err_len, mid_num, cur_len);
        if(err_len >0 && cur_len == 6 && mid_num >=1) {
            uint32_t i_idx, j_idx, j_row;
            int err_pos = cur_pos, pos;
            bg_idx = top_num+cnt_2nt[rot][seq_ch];
            ed_idx = bg_idx+mid_num-1;
 
//fprintf(stderr, "%u, err_len = %u, mid_num = %u, cur_len = %u\n", __LINE__, err_len, mid_num, cur_len);
            for(i_idx = bg_idx; i_idx<=ed_idx; ++i_idx){

//fprintf(stderr, "%u, err_len = %u, mid_num = %u, cur_len = %u\n", __LINE__, err_len, mid_num, cur_len);
                j_idx = i_idx;
                pos = cur_pos;
                int cur_rot = (rot+1)%8;
                while(pos != 0){//[fixme?] rot != 0?
                    rot_st =cur_rot*((n_data+1)/2);
                    j_row = rot_st+(j_idx+1)/2;
                    rot_ed = rot_st+n_data/2;
                    if(j_idx%2>0) {
                        seq_ch = bwt[j_row-1]&0xF; 
                    } else{
                        seq_ch = bwt[j_row]>>4;
                    }
                    n_ch = cnt_2nt[cur_rot][seq_ch+1] - cnt_2nt[cur_rot][seq_ch]; 
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
                    --pos; 
                    cur_rot  = (cur_rot+1)%8; 
                }//end while(pos != 0)+++++++++++++++++++++++++ 
                int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = j_idx;
                    aln_out[n_aln][1] = 1;
                    aln_out[n_aln][2] = err_pos;
                    algn_row[j_idx/8] |= 1<<(7-j_idx%8);
//fprintf(stderr, "%u, n_aln = %u, idx = %u\n", __LINE__, n_aln, j_idx);
                }
            }//end for+++++++++++++++++++++++++++
            aln_out[0][0] = n_aln;
            aln_out[0][1] = n_exact_aln;

//fprintf(stderr, "%u, err_len = %u, mid_num = %u, cur_len = %u\n", __LINE__, err_len, mid_num, cur_len);
            return ;
        }

//fprintf(stderr, "%u, err_len = %u, mid_num = %u, cur_len = %u\n", __LINE__, err_len, mid_num, cur_len);




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        
        
        if(mid_num >0){
            ++cur_len;
            if(cur_len <8){
                bg_idx = top_num+cnt_2nt[rot][seq_ch];
                ed_idx = bg_idx+mid_num-1;
                rot = (rot+1)%8;

            } else if(cur_len == 8){
                bg_idx = top_num+cnt_2nt[rot][seq_ch];
                ed_idx = bg_idx+mid_num-1;
//++++++++++++++++++++++++++++++++++++++++++++++++++++

if(mid_num > 1) {
    fprintf(stderr, "mid_num >0, len ==8\n"); 
    fprintf(stderr, "mid_num >0,[%u,%u]\n\n\n", bg_idx, ed_idx); 
    fprintf(stderr, "line = %u, func = %s\n", __LINE__, __func__);
    exit(1);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    ++n_aln; ++n_exact_aln;
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] = n_exact_aln;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
//fprintf(stderr, "%u, n_aln = %u, idx = %u\n", __LINE__, n_aln, bg_idx);
                }
                //return;
goto end;
            }
        } else if(mid_num ==0){
    //fprintf(stderr, "cur_len = %u, %u\t%u\n", cur_len, bg_idx, ed_idx);
//+++++++++++++++++++++++++++++++++++++++++++++++++++
//fprintf(stderr, "%s, n_data= %u\n", __func__, n_data);
//fprintf(stderr, "line = %u, func = %s\n", __LINE__, __func__);
//fprintf(stderr, "if(mid_num == 0)\n");
//return;
//exit(1);
//+++++++++++++++++++++++++++++++++++++++++++++++++++
            //cur_len = 7;
            if(cur_len == 7){
                //fprintf(stderr, "mid_num ==0 , len == 7\n"); 
                uint32_t i_idx, j_idx, j_row;
                int err_pos = cur_pos, pos;            
                for(i_idx = bg_idx; i_idx<=ed_idx; ++i_idx){
                    j_idx = i_idx;
                    pos = (cur_pos+1)%8;
                    //fprintf(stderr, "[pos]:%u\t%u\n", pos, cur_pos);
                    //fprintf(stderr, "%u\n", __LINE__);
                    int cur_rot = rot;
                    while(pos != 0){//[fixme?] rot != 0?
                        rot_st =cur_rot*((n_data+1)/2);
                        j_row = rot_st+(j_idx+1)/2;
                        rot_ed = rot_st+n_data/2;
                        if(j_idx%2>0) {
                            seq_ch = bwt[j_row-1]&0xF; 
                        } else{
                            seq_ch = bwt[j_row]>>4;
                        }
                       //n_ch = cnt_2nt[rot][seq_ch+1] - cnt_2nt[rot][seq_ch]; 
                        n_ch = cnt_2nt[cur_rot][seq_ch+1] - cnt_2nt[cur_rot][seq_ch]; 
                        
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
                        
//fprintf(stderr, "[j_idx, num, cnt]: %u\t%u\t%u\n", j_idx, num, cnt_2nt[cur_rot][seq_ch]);                    
                        --pos; 
                        cur_rot  = (cur_rot+1)%8; 
                        //fprintf(stderr, "%u\t%u\n", pos, __LINE__);
                    }//end while 
//fprintf(stderr, "[aln_idx]: %u\n", j_idx);                    
                    int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                    if(flag_aln == 0){
                        ++n_aln; 
                        aln_out[n_aln][0] = j_idx;
                        aln_out[n_aln][1] = 1;
                        aln_out[n_aln][2] = err_pos;
                        algn_row[j_idx/8] |= 1<<(7-j_idx%8);
//fprintf(stderr, "%u, n_aln = %u, idx = %u\n", __LINE__, n_aln, j_idx);
                    } 
                    //fprintf(stderr, "%u\t%u\n", pos, __LINE__);
                }//end for
                
//fprintf(stderr, "%u\t%u\n", pos,__LINE__);
                
//fprintf(stderr, "[pos]:%u\t%u\n", pos, cur_pos);
                aln_out[0][0] = n_aln;
                aln_out[0][1] = n_exact_aln;
//fprintf(stderr, "n_aln = %u, n_exact_aln = %u\n", n_aln, n_exact_aln);
                //return ;
goto end;
            } else if(cur_len <7){
//fprintf(stderr, "mid_num ==0 , len < 7\n"); 
//fprintf(stderr, "%u\t%u\t%u\n", aln_num, cur_pos, cur_len);
                //aln_num += 8-cur_len;
                aln_num ++;
                if(aln_num >7) break;
                if(cur_pos == 0) break;
                if(cur_len + st_len <7) break; 
                //if(err_num >= 1) break;
                err_num = 0;
                st_len = cur_len+st_len-7;
                st_pos[(cur_pos+1)%8] =1;
                rot=(rot+1)%8;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
           } 
        }        
    }//end while(1) 
end:
    {
        int i; 
        fprintf(stderr, "\n[align_255]: query_seq\n");
        for(i =0; i <16; ++i) fprintf(stderr, "%u", seq[i]);
        fprintf(stderr, "\nalignment seq\n");
        for(j =1; j <= n_aln; ++j){
            uint8_t aln_seq[16];
            test_smbwt_retire_seq(bwt, n_data, aln_out[j][0], aln_seq); 
            for(i =0; i <16; ++i) fprintf(stderr, "%u", aln_seq[i]);
            fprintf(stderr, "\n\n");
        }
    }
    return;
}
//+++++++++++++++++++++++++++++++++++++++++++++++
/*
uint32_t __i, __j;
for(__i = 0; __i < 8; ++__i){
    fprintf(stderr, "%u \t", __i);
    for(__j = 0; __j < 17; ++__j)
        fprintf(stderr, "%u \t", cnt_2nt[__i][__j]);
    fprintf(stderr, "\n");

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
        //log_array2(17, rot, cnt_2nt);
        accumulate_cnt_uint8(17, cnt_2nt[rot]);
        //log_array2(17, rot, cnt_2nt);
//++++++++++++++++++++++++++++++++++++
/*
if(rot==1 && n_data == 134) {
    fprintf(stderr, "cnt++++++++++++++++++++\n");
    for(__j =0; __j < 17; ++__j){
        fprintf(stderr, "%u\t", cnt_2nt[rot][__j]);

    }
    fprintf(stderr, "\n");
}*/
//++++++++++++++++++++++++++++++++++++
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
void align_min(uint8_t *bwt, int n_data, uint8_t seq[16],int err_num,  int flag, uint8_t *st_pos, uint8_t *algn_row, uint32_t (*aln_out)[3])
{
    int i, j;
    int n_aln = 0, n_exact_aln = 0;
    //aln_out[0][0] = 0;
    //aln_out[0][1] = 0;
    //aln_out[0][2] = 0;
    n_aln = aln_out[0][0];
    uint8_t *sort_seq = (uint8_t *)bwt; 
    uint32_t row;
    uint8_t seq_buf;
    int n_diff, err_pos = 0;
    int st_row = 0;
    if(aln_out[0][0] >0){
//fprintf(stderr, "%s %u, n_data = %u\n", __func__, __LINE__, n_data);
//exit(1);

        st_row = aln_out[n_data+1][0]; 
    }     

//fprintf(stderr, "%s %u, n_data = %u\n", __func__, __LINE__, n_data);
    for(row = st_row; row < n_data; row++){ 
        n_diff = 0; 
        for(i = 0; i < 4; ++i){
            for(j = 0; j < 4; ++j){
                seq_buf = (sort_seq[row*4+i]>>((3-j)*2))&3;
                if(seq_buf != seq[i*4+j] ) {
                    ++n_diff;
                    err_pos = i*4+j;
fprintf(stderr, "%s %u, n_data = %u, seq_buf = %u\n", __func__, __LINE__, n_data, seq_buf);
                }


                if(n_diff > err_num) break; 
            }
            if(n_diff > err_num) break;

//fprintf(stderr, "%s %u, n_data = %u\n", __func__, __LINE__, n_data);
        }
       
//fprintf(stderr, "%s %u, n_data = %u\n", __func__, __LINE__, n_data);
        
        if(n_diff <= err_num){

//fprintf(stderr, "%s %u, n_data = %u\n", __func__, __LINE__, n_data);
            int flag_aln = algn_row[row/8] & (1<<(7-row%8));
            if(flag_aln == 0){
            
                n_aln++;
                aln_out[n_aln][0] = row;
                aln_out[n_aln][1] = n_diff;
                aln_out[n_aln][2] = err_pos;
                aln_out[0][0] = n_aln;
                algn_row[row/8] |= 1<<(7-row%8);
                if(n_diff >0) st_pos[err_pos/2] = ++st_pos[8];
//fprintf(stderr, "%s %u, n_data = %u\n", __func__, __LINE__, n_data);
                if(n_diff == 0) {
                    aln_out[0][1] = 1;
                    for(i =0; i < 3; ++i){
                        uint32_t tmp;
                        tmp = aln_out[n_aln][i];
                        aln_out[n_aln][i] = aln_out[1][i];
                        aln_out[1][i] = tmp;    
                    }
//fprintf(stderr, "%s %u, n_data = %u\n", __func__, __LINE__, n_data);
                    if(flag > 0) { break;} 
                } 
            }// end if(flag_aln == 0)+++++++++++++++++++++++++++++++
        }

    }

//fprintf(stderr, "%s %u, n_data = %u\n", __func__, __LINE__, n_data);
    aln_out[0][0] = n_aln; 
    aln_out[0][2] = n_aln - aln_out[0][1];
    aln_out[n_data+1][0] = row+1; 

//fprintf(stderr, "%s %u, n_data = %u\n", __func__, __LINE__, n_data);
    return;   
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
    
    fprintf(stderr, "rot_off=%u, seq_ch = %u\n",rot_st, seq_ch);
    fprintf(stderr, "j_idx=%u, j_row = %u\n",j_idx, j_row);
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
void align_nosum(uint8_t *bwt, uint8_t cnt_2nt[][17], int n_data, uint8_t seq[16], int seq_len ,int seq_st, uint8_t *st_pos, uint8_t *algn_row, uint32_t (*aln_out)[3])
{
    //gen_cnt(bwt, n_data, cnt_2nt);
    int stat_flag, mid_num, top_num, bot_num, bg_row, ed_row, row, rot_st, rot_ed;
    uint32_t bg_idx, ed_idx, j;//bg_idx,ed_idx begin and end of 2nt-encoded bwt 
    uint8_t seq_ch, ch;
    int n_aln = 0, n_exact_aln = 0;
    n_aln = aln_out[0][0];
    uint32_t cur_pos = (seq_st+1)%8; 
    int rot = (7-seq_st)%8; //bw rotation
    int cur_len = 0; // current alignment length
    int st_len = seq_len;
    uint8_t err_len = st_pos[9]; 
    bg_idx = 0; ed_idx = n_data-1;//bwt interval [bg_idx, ed_idx]
    n_aln = aln_out[0][0];
    int aln_num = 0, err_num = 0;
    int i;
/*  
fprintf(stderr, "align_seq\n");

fprintf(stderr, "%2u ", 0);
for(i = 0; i < 8; ++i){
    fprintf(stderr, "%2u ", (seq[i*2]<<2)|(seq[i*2+1]));
}
fprintf(stderr, "\n");
uint8_t *seq_buf = calloc(100000, 1);
test_AlgnBwtSml(bwt, n_data, seq_buf);
uint32_t __j;
for(__j= 0; __j<n_data; __j++){
    uint8_t *cur_seq = seq_buf+16*__j;
    fprintf(stderr, "%2u ", __j);
    for(i = 0; i < 8; ++i){
        fprintf(stderr, "%2u ", (cur_seq[i*2]<<2)|(cur_seq[i*2+1]));
    }
    fprintf(stderr, "\n");

}
for(i = 0; i < 8; ++i) fprintf(stderr, "%u\t", (seq[i*2]<<2)|(seq[i*2+1]));
fprintf(stderr, "\n");

free(seq_buf);
*/
    while(1){
        cur_pos = (cur_pos+7)%8; 
        seq_ch = seq[cur_pos*2]<<2| seq[cur_pos*2+1];


/*  
{
fprintf(stderr, "\n========================================================\n");
fprintf(stderr, "(===bg_idx, ed_idx, n_data) = (%u, %u, %u)\n", bg_idx, ed_idx, n_data);
fprintf(stderr, "(===cur_len, cur_rot) = (%u, %u)\n", cur_len, rot);


int rot_off = rot*((n_data+1)/2);
int __n0 = __dna2_count(bwt+rot_off, 0, bg_idx, seq_ch);
int __tot = __dna2_count_small(bwt+rot_off, 0, n_data, seq_ch);
int __n1 = __dna2_count(bwt+rot_off, bg_idx, ed_idx+1, seq_ch);
int __n2 = __dna2_count(bwt+rot_off, ed_idx+1, n_data, seq_ch);
fprintf(stderr, "seq_ch, tot, n0, n1, n2= %u %u %u %u %u\n", seq_ch, __tot, __n0, __n1, __n2);
fprintf(stderr, "bg_idx, ed_idx= %u %u\n", __n0+__tot, __n0+__tot+__n1-1);
fprintf(stderr, "\n--------------------------------------------------------\n");
}
*/
        if(n_data == 1){
            //恢复序列算法
            int i;
            rot = 0; 
            uint8_t bwt_seq[16];
            uint32_t k = 0;

            for(rot=0; rot<8; ++rot){
                uint32_t x = rot*((n_data+1)/2*2);
                //uint8_t ch = __get_bwt2(bwt, x+k);
                uint8_t ch = bwt[(x+k)/2]; 
                ch = (k%2==0)?ch>>4:ch&15;
                bwt_seq[15-rot*2] = ch&3;
                bwt_seq[14-rot*2] = ch>>2;
                uint32_t occ = __dna2_count(bwt, x, x+k, ch);
                uint32_t C = cnt_2nt[rot][ch]; 
                k = C+occ;
            }
            int n_diff = 0;
            for(i = 0; i < 16; ++i){
                if(bwt_seq[i]!= seq[i]) n_diff++; 
            }
            if(n_diff <=2){
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    n_aln++;
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = n_diff;
                    aln_out[n_aln][2] = 0;
                    aln_out[0][0] = n_aln;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
                    if(n_diff == 0) aln_out[0][1] = 1;
                }
                return;
            }
        } else if(bg_idx == 0 && ed_idx==n_data-1){
/*
fprintf(stderr, "align seq %u\t%u\n",  cur_pos, seq_ch);
for(j = 0; j <= (n_data+1)/2; ++j){
    fprintf(stderr, "%u\t%u\n",  bwt[j]>>4, bwt[j]&0xF);
}
*/
            if(cnt_2nt[rot][seq_ch] == cnt_2nt[rot][seq_ch+1]){
                if(err_num >=1) {break;}
                err_num++; 
                st_len = (cur_len+1)%8;
                rot=(rot+1)%8;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
//fprintf(stderr, "%u, cur_len = %u, rot = %u\n", __LINE__, cur_len, rot);
                continue;
            }
            
            bg_idx = cnt_2nt[rot][seq_ch]; 
            ed_idx = cnt_2nt[rot][seq_ch+1]-1;
            rot=(rot+1)%8;
            ++cur_len;
 
//fprintf(stderr, "%u, cur_len = %u, rot = %u\n", __LINE__, cur_len, rot);
            if(cur_len == 8) {
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    ++n_aln;
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
                }
                return;
            } 
            
//fprintf(stderr, "[top_num, mid_num, bot_num] = %u\t%u\t%u\n", top_num, mid_num, bot_num);
//fprintf(stderr, "align seq %u\t%u\n",  cur_pos, seq_ch);
//fprintf(stderr, "[bg_idx, ed_idx] = %u\t%u\n",  bg_idx, ed_idx);
           
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
/*  
fprintf(stderr, "[top_num, mid_num, bot_num] = %u\t%u\t%u\n", top_num, mid_num, bot_num);
fprintf(stderr, "align seq %u\t%u\n",  cur_pos, seq_ch);
fprintf(stderr, "[bg_idx, ed_idx] = %u\t%u\n",  bg_idx, ed_idx);
*/
/*  
        for(j = rot_st; j <= rot_ed; ++j){
            fprintf(stderr, "%u\t%u\n",  bwt[j]>>4, bwt[j]&0xF);
        }
*/
/*
*/
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        if(err_len>0){
            if(mid_num < 1 && cur_len <= 6 ) { return;}
            if(mid_num >=1 && cur_len < 6){
                ++cur_len;
                if(cur_len <8-err_len){
                    bg_idx = top_num+cnt_2nt[rot][seq_ch];
                    ed_idx = bg_idx+mid_num-1;
                    rot = (rot+1)%8;
                    continue;
                }
            }
        }

        if(err_len >0 && cur_len == 6 && mid_num >=1) {
            uint32_t i_idx, j_idx, j_row;
            int err_pos = cur_pos, pos;
            bg_idx = top_num+cnt_2nt[rot][seq_ch];
            ed_idx = bg_idx+mid_num-1;
             
            for(i_idx = bg_idx; i_idx<=ed_idx; ++i_idx){
                j_idx = i_idx;
                pos = cur_pos;
                int cur_rot = (rot+1)%8;
                while(pos != 0){//[fixme?] rot != 0?
                    rot_st =cur_rot*((n_data+1)/2);
                    j_row = rot_st+(j_idx+1)/2;
                    rot_ed = rot_st+n_data/2;
                    if(j_idx%2>0) {
                        seq_ch = bwt[j_row-1]&0xF; 
                    } else{
                        seq_ch = bwt[j_row]>>4;
                    }
                    n_ch = cnt_2nt[cur_rot][seq_ch+1] - cnt_2nt[cur_rot][seq_ch]; 
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
                    --pos; 
                    cur_rot  = (cur_rot+1)%8; 
                }//end while(pos != 0)+++++++++++++++++++++++++ 
                int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = j_idx;
                    aln_out[n_aln][1] = 1;
                    aln_out[n_aln][2] = err_pos;
                    algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                }
            }//end for+++++++++++++++++++++++++++
            aln_out[0][0] = n_aln;
            aln_out[0][1] = n_exact_aln;
            return ;
        }


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++ 




        if(mid_num >0){
            ++cur_len;
            if(cur_len <8){
                bg_idx = top_num+cnt_2nt[rot][seq_ch];
                ed_idx = bg_idx+mid_num-1;
                rot = (rot+1)%8;

            } else if(cur_len == 8){ 
//++++++++++++++++++++++++++++++++++++++++++++++++++

if(mid_num > 1) {
    fprintf(stderr, "mid_num = %u, cur_len == 8", mid_num);
    fprintf(stderr, "line = %u, func = %s\n", __LINE__, __func__);
    exit(1);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++
                bg_idx = top_num+cnt_2nt[rot][seq_ch];
                ed_idx = bg_idx+mid_num-1;
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    n_aln++; n_exact_aln++;
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
                }
                return;
            }
        } else if(mid_num == 0){
//++++++++++++++++++++++++++++++++++++++++++++++++++++
//fprintf(stderr, "%s, n_data= %u\n", __func__, n_data);
//fprintf(stderr, "cur_len = %u, rot = %u\n",cur_len, rot);
//fprintf(stderr, "line = %u, func = %s\n", __LINE__, __func__);
//exit(1);
//++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(cur_len == 7){
                uint32_t i_idx, j_idx, j_row;
                int err_pos = cur_pos, pos;
                for(i_idx = bg_idx; i_idx<=ed_idx; ++i_idx){
                    j_idx = i_idx;
                    pos = (cur_pos+1)%8;
                    int cur_rot = rot;
                    while(pos != 0){//[fixme?] rot != 0?
                        rot_st =cur_rot*((n_data+1)/2);
                        j_row = rot_st+(j_idx+1)/2;
                        rot_ed = rot_st+n_data/2;
                        if(j_idx%2>0) {
                            seq_ch = bwt[j_row-1]&0xF; 
                        } else{
                            seq_ch = bwt[j_row]>>4;
                        }
                        n_ch = cnt_2nt[cur_rot][seq_ch+1] - cnt_2nt[cur_rot][seq_ch]; 
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
                        --pos; 
                        cur_rot  = (cur_rot+1)%8; 
                    }//end while(pos != 0)+++++++++++++++++++++++++ 
                    int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                    if(flag_aln == 0){
                        n_aln++; 
                        aln_out[n_aln][0] = j_idx;
                        aln_out[n_aln][1] = 1;
                        aln_out[n_aln][2] = err_pos;
                        algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                    }
                }//end for+++++++++++++++++++++++++++
                aln_out[0][0] = n_aln;
                aln_out[0][1] = n_exact_aln;
                return ;
            } else if(cur_len <7){
    
       //fprintf(stderr, "aln_num = %u, cur_pos = %u, cur_len+st_len =%u, err_num = %u\n", aln_num, cur_pos, cur_len+st_len, err_num); 
                //aln_num += 8-cur_len;
                aln_num++;
                //if(aln_num >8) break;
                if(aln_num >7) break;
                if(cur_pos == 0) break;
                if(cur_len + st_len <7) break; 
                //if(cur_len + st_len <7) break; 
                //if(err_num >= 1) break;
                err_num = 0;
                st_len = cur_len+st_len-7;
                rot=(rot+1)%8;
                st_pos[(cur_pos+1)%8] = 1;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
            } 
        }// if(mid_num == 0)+++++++++++++++++++++++++++++++++
    }//end while(1)++++++++++++++++++++++++++++++++++++++++++
   return; 
}
/*
fprintf(stderr, "\n========================================================\n");
fprintf(stderr, "(===bg_idx, ed_idx, n_data) = (%u, %u, %u)\n", bg_idx, ed_idx, n_data);
fprintf(stderr, "(===cur_len, cur_rot) = (%u, %u)\n", cur_len, cur_rot);
int __n0 = __dna2_count(Bwt+rot_off, 0, bg_idx, seq_ch);
int __tot = __dna2_count_small(Bwt+rot_off, 0, n_data, seq_ch);
int __n1 = __dna2_count(Bwt+rot_off, bg_idx, ed_idx+1, seq_ch);
int __n2 = __dna2_count(Bwt+rot_off, ed_idx+1, n_data, seq_ch);
fprintf(stderr, "seq_ch, tot, n0, n1, n2= %u %u %u %u %u\n", seq_ch, __tot, __n0, __n1, __n2);
fprintf(stderr, "bg_idx, ed_idx= %u %u\n", __n0+__tot, __n0+__tot+__n1-1);
fprintf(stderr, "\n--------------------------------------------------------\n");
*/  
/*
{
uint8_t *seq_buf = calloc(100000, 1);
test_AlgnBwtSml_large(Bwt, n_data, seq_buf);
uint32_t __j;
for(__j= 0; __j<n_data; __j++){
    uint8_t *cur_seq = seq_buf+16*__j;
    fprintf(stderr, "%2u ", __j);
    for(i = 0; i < 8; ++i){
        fprintf(stderr, "%3u ", (cur_seq[i*2]<<2)|(cur_seq[i*2+1]));
    }
    fprintf(stderr, "\n");

}
for(i = 0; i < 8; ++i) fprintf(stderr, "%u\t", (seq[i*2]<<2)|(seq[i*2+1]));
fprintf(stderr, "\n");

free(seq_buf);
}
*/
/*
if(n_data == 571){
uint32_t __i, __j, __k, __l;
fprintf(stderr, "\n----------------------\n");
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        fprintf(stderr, "seq_ch = %u\t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            fprintf(stderr, "%u\t", pSum[__j+__l]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "-------------------\n\n");
}
fprintf(stderr, "\n+++++++++++++++++\n");
}
*/
//+++++++++++++++++++++++++++++++++++++++++++++++++++++
    /*
            if(cnt_2nt[rot][seq_ch] == cnt_2nt[rot][seq_ch+1]){
                ++cur_len;
                if(cur_len >= 2) break;
                cur_len %= 8;
                len = cur_len;
                bg_idx = 0;
                ed_idx = n_data-1; 
                ++rot;
                rot%=8;
                continue;
            }
            bg_idx = cnt_2nt[rot][seq_ch]; 
            ed_idx = cnt_2nt[rot][seq_ch+1]-1;
            rot=(rot+1)%8;
            ++cur_len;
            if(cur_len == 8) {
                aln_out[1][0] = bg_idx;
                aln_out[1][1] = 0;
                aln_out[1][2] = 0;
                aln_out[0][0] = aln_out[0][1] =1;
                return;
            } 
    */
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         
//fprintf(stderr, "[%u]: error %s\n", __LINE__, __func__);
                //exit(1); 
         /*  
               ++cur_len;
                if(cur_len >= 2) break;
                cur_len %= 8;
                len = cur_len;
                bg_idx = 0;
                ed_idx = n_data-1; 
                ++cur_rot;
                cur_rot%=8;
        */


void align_large_alt(uint8_t *Bwt, uint8_t *cnt2, uint32_t n_data, const uint8_t *seq , int seq_len, int seq_st, uint8_t *st_pos, uint8_t *algn_row, uint32_t (*aln_out)[3]){
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
    uint8_t err_len = st_pos[9];
    int bg_idx = 0, ed_idx = n_data-1;
    int alg_num = 0, bgn_rot = 0, siz_cnt, bgn_stat;
    int cur_pos = (seq_st+1)%8;
    int cur_rot = (7-seq_st)%8;   
    int cur_len = 0, st_len = seq_len, err_num = 0;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    if( n_data <= 254*256 ) siz_cnt = 2 ;
    if( n_data > 254*256 ) siz_cnt = 4 ;
    int rng_num = ((n_data+253) /254) ;
    int len_sum = (8+ siz_cnt) *( rng_num /8) + rng_num %8 + siz_cnt ;
    int rot_sum =cur_rot*16*len_sum  ;
    int data_q = (n_data)/254;
    int data_r = (n_data)%254;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // bg_idx 和ed_idx是下一个轮换内部的相对位置+++++
    uint32_t upp_num=0, mid_num=0, seq_num[2];
    uint32_t rng_cnt[2], rng_end[2], ch_cnt[2], rng_len, fst_row, bgn_row, end_row, lst_row, row;
    while(1){
        int rot_sum =cur_rot*16*len_sum  ;
        
        int ch_off, idx_num, seq_ch;
        uint8_t get_ch, ch;
        cur_pos = (cur_pos+7)%8;
        seq_ch = (seq[2*cur_pos]<<2) | seq[2*cur_pos+1] ;
        rot_off = cur_rot*((n_data+1)/2);

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
                cur_rot=(cur_rot+1)%8;
                cur_len = 0;
                bg_idx = 0;
                ed_idx = n_data-1; 
                continue;
            }
            cur_len = (cur_len+1)%8;  
            cur_rot = (cur_rot+1)%8;
             

         
            //rot=(rot+1)%8;
            //++cur_len;
            if(cur_len == 8) {
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
                    n_aln++; 
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] =1;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);
//fprintf(stderr, "%s, %u, idx = %u, err_pos = %u\n", __func__,__LINE__, bg_idx, 0); 

                }
                return;
            } 
            
       
           //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
            continue;
        }
        mid_num = 0; upp_num = 0;
        rot_off = cur_rot*((n_data+1)/2);
        ch_off = 16*cur_rot*len_sum + len_sum*seq_ch ;
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
//fprintf(stderr, "sum_buf = %u\t", sum_buf);
                    ch_cnt[i] = ((ch_cnt[i]<<8)|sum_buf); 
                }

//fprintf(stderr, "\n");
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
                    if(idx_r[i] % 2 == 0 ){// 该点已经计算过，但需要考虑bg_idx点是否包含在其中++++ 
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
//fprintf(stderr, "i=%u, rng_len = %u, rng_end = %u, rot_off = %u\n", i, rng_len, rng_end[i], rot_off);
//fprintf(stderr, "idx_row = %u, idx_i[i] = %u, idx_r = %u\n", idx_row[i], idx_i[i], idx_r[i]);
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
                            //if(i== 0) ++upp_num2[i]; 
                            --upp_num2[i];
                        } 
                   
                        //if(seq_ch == ch) --upp_num2[i]; 
                    }


                }

                upp_num2[i] += ch_cnt[i] ;
            } // End: ++++++++++++++ for(i = 0 ;  i<2 ; i++ )
        }// End : +++++++ if(bg_idx_q < ed_idx_q )+++++++++++++++++++
        
        //*****************************************************************
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
        if(err_len>0){
            if(mid_num < 1 && cur_len <= 6 ) { return;}
            if(mid_num >=1 && cur_len < 6){
                ++cur_len;
                if(cur_len <8-err_len){
                    bg_idx =  upp_num2[0]; 
                    ed_idx =  upp_num2[1]-1;
                    cur_rot = (cur_rot+1)% 8;
                    continue; 
                }
            }
        }

        if(err_len >0 && cur_len == 6 && mid_num >=1) {
            uint32_t i_idx, j_idx;
            int pos, err_pos= (cur_pos+7)%8;
            bg_idx =  upp_num2[0]; 
            ed_idx =  upp_num2[1]-1;

//fprintf(stderr, "%s, %u, bg_idx = %u, ed_idx = %u\n", __func__,__LINE__, bg_idx, ed_idx); 
            for(i_idx = bg_idx ;  i_idx <= ed_idx ;  i_idx++){
                j_idx = i_idx;    
                pos = cur_pos;
                // 查找某一近似比对成功序列的行号 ++++++++++++++++++
                uint32_t rot= (cur_rot+1)%8;
                /*计算编辑距离*/
                while(pos != 0){
//fprintf(stderr, "%u, pos = %u, rot = %u\n", __LINE__, pos, rot); 
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
                    rot_off = rot*((n_data+1)/2);
                    ch_off = 16*rot*len_sum + len_sum*seq_ch ;
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
                    pos = (pos +7)%8;
                    rot = (rot +1)%8;
                }// End : while(pos != 0) ----------------------------------
                int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));

//fprintf(stderr, "%s, %u, idx = %u, err_pos = %u\n", __func__,__LINE__, j_idx, err_pos); 
                if(flag_aln == 0){
                    n_aln++;
                    aln_out[n_aln][0] = j_idx;
                    aln_out[n_aln][1] = 1;
                    aln_out[n_aln][2] = err_pos;
                    algn_row[j_idx/8] |= 1<<(7-j_idx%8);
                
//fprintf(stderr, "%s, %u, idx = %u, err_pos = %u\n", __func__,__LINE__, j_idx, err_pos); 
                }
            }// for(i_idx = bg_idx;i_idx <= ed_idx; i_idx ++)--------------------------
            aln_out[0][0] = n_aln; 
            aln_out[0][1] = n_exact_aln; 
            return;
        }


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        
        
        
        if( mid_num>0) { // 比对尚未结束++++++++++
            cur_len += 1;  
            if( cur_len < 8){
                bg_idx =  upp_num2[0]; 
                ed_idx =  upp_num2[1]-1;
                cur_rot = (cur_rot+1)% 8;
                continue;
            }
            if( cur_len == 8){

                bg_idx =  upp_num2[0]; 
                ed_idx = upp_num2[1]-1;
//++++++++++++++++++++++++++++++++++++++++++++++++++
if(mid_num > 1) {
    fprintf(stderr, "mid_num = %u, cur_len == 8", mid_num);
    fprintf(stderr, "line = %u, func = %s\n", __LINE__, __func__);
    exit(1);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++
//fprintf(stderr, "cur_len == 8, cur_rot = %u,  bg_idx = %u, ed_idx = %u\n",cur_rot, bg_idx, ed_idx);
//fprintf(stderr, "\n========================================================\n");
                int flag_aln = algn_row[bg_idx/8] & (1<<(7-bg_idx%8));
                if(flag_aln == 0){
            
                    alg_err_num  = 0 ;
                    n_aln++; n_exact_aln++;
                    aln_out[n_aln][0] = bg_idx;
                    aln_out[n_aln][1] = 0;
                    aln_out[n_aln][2] = 0;
                    aln_out[0][0] = n_aln;
                    aln_out[0][1] = n_exact_aln;
                    algn_row[bg_idx/8] |= 1<<(7-bg_idx%8);


                }
    //+++++++++++++++++++++++++++++
/*
if(n_data==305){
    fprintf(stderr, "2227, n_aln = %u, bg_idx = %u\n", n_aln, bg_idx);
}*/
//+++++++++++++++++++++++++++++
                return;     
            }
        } // End : if(mid_num>0) // 比对尚未结束------
        if(mid_num==0){//本次 比对结束 +++++++++
//++++++++++++++++++++++++++++++++++++++++++++

//fprintf(stderr, "cur_len = %u, bg_idx = %u, ed_idx = %u\n", cur_len, bg_idx, ed_idx);
//fprintf(stderr, "cur_pos = %u, cur_rot = %u, seq_ch = %u\n", cur_pos, cur_rot, seq_ch);
//fprintf(stderr, "line = %u, func = %s\n", __LINE__, __func__);
/*  
uint32_t __i, __j, __k, __l;
fprintf(stderr, "\n----------------------\n");
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        fprintf(stderr, "seq_ch = %3u ", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            fprintf(stderr, "%3u ", pSum[__j+__l]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "-------------------\n\n");
}
*/
//fprintf(stderr, "if(mid_num ==0) \n+++++++++++++++++\n");
//exit(1);
//++++++++++++++++++++++++++++++++++++++++++++++

            if(cur_len == 7){ //输出近似比对成功结果+++++++++++++++
                //比对上的最多有15个，都不相等，有序，可以计算汉明距离
                //比对上的序列号 ；   //编辑位置;   //编辑值;   //alg_ err_num = 1或者2 ;  
                // 查找每个近似比对成功序列的行号 ++++++++++++++++++
                uint32_t i_idx, j_idx;
                int pos, err_pos= cur_pos;
                for(i_idx = bg_idx ;  i_idx <= ed_idx ;  i_idx++){
                    j_idx = i_idx;    
                    pos = (cur_pos+1)%8;
                    // 查找某一近似比对成功序列的行号 ++++++++++++++++++
                    uint32_t rot= cur_rot;
                    /*计算编辑距离*/
                    while(pos != 0){
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
                        rot_off = rot*((n_data+1)/2);
                        ch_off = 16*rot*len_sum + len_sum*seq_ch ;
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
                        pos = (pos +7)%8;
                        rot = (rot +1)%8;
                    }// End : while(pos != 0) ----------------------------------
                    int flag_aln = algn_row[j_idx/8] & (1<<(7-j_idx%8));
                    if(flag_aln == 0){
            
                        n_aln++;
                        aln_out[n_aln][0] = j_idx;
                        aln_out[n_aln][1] = 1;
                        aln_out[n_aln][2] = err_pos;
                        algn_row[j_idx/8] |= 1<<(7-j_idx%8);


                    }
                }// for(i_idx = bg_idx;i_idx <= ed_idx; i_idx ++)--------------------------
                aln_out[0][0] = n_aln; 
                aln_out[0][1] = n_exact_aln; 
                return;
            } //  if(cur_len == 7)  //近似比对成功结果处理结束------------------

            if(cur_len < 7){ // 重新选定比对起始位置进行比对 ++++++++++++++++
                //保存本次比对的相关信息
                //包括[bg_idx, ed_idx]， len ，pos, cur_rot, seq_ch  
                //alg_num  += (8-cur_len);
                alg_num++;
                if(alg_num >7) {break;}
                if(cur_pos ==0) {break;}
                if(cur_len +st_len <7) { break;} 
                //if(err_num >=1) break;
                err_num = 0;
                st_len = cur_len+st_len-7;
                st_pos[(cur_pos+1)%8] =1;
                cur_rot = (cur_rot+1)%8;
                cur_len = 0 ;
                bg_idx = 0; 
                ed_idx = n_data-1;
            }  // End:  ++++++++++++++++++if(cur_len < 7) 
        }// End : +++++++++++++++++++if((mid_num==0)
    }//++++++++++++++++++++++++++++ while(1)
    return;
} // End : Alignm_big( )------------------------------------------------------
//++++++++++++++++++++++++++++


//++++++++++++++++++++++++++++
/*
if(n_data == 419){
uint32_t __i, __j, __k, __l;
fprintf(stderr, "\n----------------------\n");
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        fprintf(stderr, "seq_ch = %u\t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            fprintf(stderr, "%u\t", pSum[__j+__l]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "-------------------\n\n");
}
fprintf(stderr, "\n+++++++++++++++++\n");
}*/
//++++++++++++++++++++++++++++++++++

/*
uint32_t __i, __j, __k, __l;
fprintf(stderr, "\n----------------------\n");
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        fprintf(stderr, "seq_ch = %u\t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            fprintf(stderr, "%u\t", pSum[__j+__l]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "-------------------\n\n");
}
fprintf(stderr, "\n+++++++++++++++++\n");

//++++++++++++++++++++++++++++++++++
uint8_t *cnt = calloc(800000, sizeof(uint8_t));
gen_cnt2(Bwt, n_data, cnt);
fprintf(stderr, "\n----------------------\n");
for(__k =0; __k < 8; __k++){
    for(__i = 0; __i < 16; ++__i){
        fprintf(stderr, "seq_ch = %u \t", __i);
        __j = 16*__k*len_sum + len_sum*__i;
        for(__l = 0; __l < len_sum; ++__l){
            fprintf(stderr, "%u\t", cnt[__j+__l]);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "-------------------\n\n");
}
fprintf(stderr, "\n+++++++++++++++++\n");




free(cnt);
exit(1);
*/
//++++++++++++++++++++++++++++++++++
#endif
