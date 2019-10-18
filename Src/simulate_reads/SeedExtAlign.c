// 循环控制信息--------------------------------------------
//uint32_t queue_align_info[][2];  //开辟较大内存
//uint32_t this_stack_tree[NUM_EXT_SEED][5];  
//[i][0]是第i个节点(即，第i级索引）在queue_align_info中的起始行号
//[i][1]是第i个节点(即，第i级索引）在queue_align_info中的个数
//uint32_t this_stack_node[3];
#include <stdio.h>
#include "SeedExtStruct.h"
#include "ksort.h"
#include "kvec.h"
#include "bwt.h"
#include "setFileName.h"
#include "SeedExtLib.h"
#include "query.h"
#include "ksw.h"
#include "lookup.h"
#include "seed.h"
//#include "kseq.h"
//KSEQ_INIT(gzFile, gzread)
#define SEED_LEN 20

#define FILE_PREFIX "../../BWT_INDEX/prefix"
const char *read_file = "../../Reads/sample.fq";
#define MAX_READS 100000
#define MIN_Q 255

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define READ_LEN 200
//#define NUM_EXT_SEED 2
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void set_err_mdl(uint8_t err_mdl[]);
void set_read_data(idx_t *fm_idx, uint32_t read_id, uint8_t *err_mdl, int read_len, uint8_t *f_read_seq, uint8_t *r_read_seq, uint8_t *qual, uint8_t *r_qual);

void set_read_data_1(idx_t *fm_idx, uint32_t read_id, uint8_t *err_mdl, query_t *query);

typedef struct{
    uint8_t name[128];
    uint8_t comment[128];
    uint8_t seq[READ_LEN];
    uint8_t rseq[READ_LEN];
    uint8_t qual[READ_LEN];
    uint8_t apos[READ_LEN];
} read_t; 
typedef struct{
    int mov_L;
    int mov_R;
    int num_L;
    int num_R;
    int num_T;
    uint8_t *apos;
} err_t;

int main(){	
    int NUM_EXT_SEED = 2;
    int REAL_FLAG = 1;
    //+++++++++++++++++++++++++++++++++++++++++++
	//系统参数初始化，包括：
	//1.全局变量数组头指针，
	//2.全局控制变量，
	//3.文件名字符数组，以及全局文件缓冲区变量，

	//-----------------------------------------
	//主Bwt索引初始化
	//文件中读入主BWT的序列数组SA[],累加数组Rank[]
	//文件中读入参考序列数组RefSeq[]
	//文件中读入Index到Pos映射数组ItoPos[];
	//------------------------------------------
	idx_t *fm_idx = (idx_t *)calloc(1, sizeof(idx_t)); 
	idx_restore(FILE_PREFIX, fm_idx);
	uint32_t *I2Pos= fm_idx->bwt->sa;
	//+++++++++++++++++++++++++++++++++++++++++++
	//以下是多级索引初始化
    struct JmpMod *jmp;
	struct FileName  *fn;	
    struct ExtBlck   eBlck[NUM_EXT+1];
	struct StackTree *sTree = ( struct StackTree *)malloc(sizeof(struct StackTree));
	struct SubBuf    *sub;
	struct CapIfo    *cap;
	if(NULL == (fn = (struct FileName*)malloc(sizeof(struct FileName)))){
	    perror("error:[malloc(NUM_EXT_SEED*sizeof(FileName))]");
	    exit(1);
	}
    setFileName(fn);
   
    InitIdxsArry(fn, &sub, eBlck,&sTree,cap,&jmp);
    //struct ExtBlck *iB = calloc(1, sizeof(ExtBlck));
    struct ExtBlck *cB = eBlck;
    


            
    int i, j;
    uint32_t hash_boundry[12] = {};
    uint8_t kmer[12] = {};
    for(i =0; i < 12; ++i){
        kmer[i] = __get_pac(fm_idx->pac, fm_idx->bwt->seq_len-12+i);
    }
    uint32_t k = 0, l = 0;
    for(i =0; i < 11; ++i){
        int n = bwt_match_exact_alt(fm_idx->bwt, 1, kmer+12-1-i, &k, &l);
        hash_boundry[i] = k;
    }
    for(i = 0; i < 11; ++i){
        for(j = i+1; j < 11; ++j){
            if(hash_boundry[i] > hash_boundry[j]){
                uint32_t tmp;
                tmp = hash_boundry[i];
                hash_boundry[i] = hash_boundry[j];
                hash_boundry[j] = tmp;         
            }
        }
    }
 
    fprintf(stderr, "NUM_EXT_SEED = %u\n", NUM_EXT_SEED);
    fprintf(stderr, "IS_SMLSIZ = %u\n", IS_SMLSIZ);
    fprintf(stderr, "READ_LEN = %u\n", READ_LEN);
    fprintf(stderr, "SEED_LEN = %u\n", SEED_LEN);
    
    //+++++++++++++++++++++++++++++++++++++++++++
	//1.readArry数组开辟空间，
	//2.Read,Seed等的状态信息
	//3.输入输出格式信息
	//4.其它
	//
	//+++++++++++++++++++++++++++++++++++++++++++
    
    uint8_t  *f_read_seq, *r_read_seq, *read_seq;
	uint8_t	 seed_seq[32];
    char     *buf_read;	
    int      seed_off;
	uint32_t get_num;
	uint32_t read_id;
    uint32_t buf_algn[20]; 	
	int      Flg_Algn  = 0 ;
	uint32_t cur_read_num = 0;

    queryio_t *qs = query_open(read_file);;
    read_t *read_buf = (read_t *)calloc(MAX_READS, sizeof(read_t)); 
    query_t *query_buf = (query_t *)calloc(MAX_READS, sizeof(query_t));
    for(i = 0; i < MAX_READS; ++i){
        query_buf[i].name = read_buf[i].name; 
        query_buf[i].comment = read_buf[i].comment; 
        query_buf[i].seq = read_buf[i].seq; 
        query_buf[i].rseq = read_buf[i].rseq; 
        query_buf[i].qual = read_buf[i].qual; 
        query_buf[i].apos = read_buf[i].apos; 
    } 
    int read_num = 0;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//以下代码是对参考基因组中获取read数据进行变异错误的模版
  

    //以下代码段是种子相关变量初始化 
    seed_t *seed = calloc(1, sizeof(seed_t)); 
    init_seed_model(seed);
    int SEED_SLC_SIZE = seed->MAX_SIZE;//应该满足SEED_SLC_SIZE>seed_slc_size
    int MAX_SLC = SEED_SLC_SIZE+2;
    
    uint32_t pos_buf[IS_SMLSIZ+1];
    int seed_slc_size, seed_slc_size_i;
    int seed_id;
  
uint8_t err_mdl[READ_LEN+20] = {}; 
set_err_mdl(err_mdl);  
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //以下代码段是read相关变量定义及其常用变量    
    uint32_t read_len = READ_LEN;
    uint8_t  ch[READ_LEN], r_ch[READ_LEN];
    uint8_t  qual[READ_LEN] = {}, r_qual[READ_LEN] = {}; 
    f_read_seq = ch;
    r_read_seq = r_ch;
    query_t *query;
            
    int l_seed_st = 0, r_seed_ed=READ_LEN;
    for(i = 0; i < SEED_SLC_SIZE; ++i){
        if(seed->slc_i[i].s_off>l_seed_st) l_seed_st = seed->slc_i[i].s_off; 
        if(seed->slc_i[i].s_off<r_seed_ed) r_seed_ed = seed->slc_i[i].s_off; 
    } 
    //while( (n_seqs = query_read_query_buf(qs,MAX_READS, query_buf)) >0){ //主循环++++++++++++++++++++++++++++++++
    
    uint32_t bgn,end,num,bgn_row,end_row,cap_row, old_seed_id;
    int tot_reads = 0, n_seqs;
    uint32_t aln_result_sum[100] = {};   
    int seed_slc_num = 0; 
    int max_stre_cls = 0;
    while(1){ //主循环++++++++++++++++++++++++++++++++
//if(tot_reads >= 10) break;
        if(REAL_FLAG == 1){ 
            n_seqs = query_read_multiSeqs(qs,MAX_READS, query_buf);
            if(n_seqs == 0) break;
//n_seqs = 11;

        }
        // 循环过程中每次读入规定个数的reads
        tot_reads += n_seqs;

//++++++++++++++++++++++++++++++++
//测试代码
        //seed_off = 0;
        read_id = 0;

        seed_id = 0;
        seed_slc_size = 1;
        //++++++++++++++++++++++++++++++++
        //一次读入的read数组循环处理
       
        //set_read_data(fm_idx, 0, err_mdl, read_len, f_read_seq, r_read_seq, qual, r_qual);
        query = query_buf;
if(REAL_FLAG == 0){    
            //read_id = 140000000;
            //n_seqs = 240000000;
            read_id = 0;
            n_seqs = 100000;
           
            query->l_seq = READ_LEN;
            query->seq = ch; 
            query->rseq = r_ch; 
            query->qual = qual; 
            set_read_data_1(fm_idx, read_id, err_mdl, query);
            //for(i = 0; i <query_buf[0].l_seq; ++i ) fprintf(stderr, "%u", ch[i]);
            //fprintf(stderr, "----------\n");
        } 



   

/*  
           
*/
   
    int read_flag = 0; 
    while(2) {
        if(read_id >= n_seqs) break;      
        printf("read_id = %u\n", read_id);
        fprintf(stderr, "read_id = %u\n", read_id);
            
 
if(REAL_FLAG == 1){  
    seed_slc_num = 0;
    max_stre_cls = 0;
    sTree->cls = 0;
    query = query_buf+read_id;
    sTree->len_arry = 0;
    sTree->len_buf = 0;
/*  
    if(query->n_ambiguous >1) {
        read_id++;
        continue;
    }
*/
} else{

    query = query_buf+0;
    query->l_seq = READ_LEN;
    query->seq = ch; 
    query->rseq = r_ch; 
    query->qual = qual;
    query->n_ambiguous = 0; 
    set_read_data_1(fm_idx, read_id, err_mdl, query);

}
++read_id;


        int left_len, right_len, buf[200];
        if(query->n_ambiguous == 0){
            read_flag = 0; 
        } else {
//if(query->n_ambiguous >2) aln_result_sum[7]++;
//if(query->n_ambiguous ==2) aln_result_sum[6]++;
            left_len = query->l_seq -l_seed_st;
            right_len = query->l_seq -r_seed_ed+SEED_LEN; 
            for(i = 0; i < query->n_ambiguous; ++i){
                int j = 0;
                if(query->apos[i] >= left_len && query->apos[i] <= right_len){
                    buf[j++] = query->apos[i]; 
                }     
            }
            if(j < 2) {read_flag = j;}
            else {
                for(i = 0; i <j; ++i){
                    if(buf[i+1]-buf[i] <= SEED_LEN) read_flag = 10; 
                }
                if(read_flag < 10) read_flag = 1;
            } 
        }
        //以下代码处理非种子扩展比对过程
        //if(read_flag > 1)  aln_no_seed_ext(fm_idx, hash_boundry, query);       
            //aln_no_seed_ext(idx_t *fm_idx, uint32_t *hash_boundry, query_t* query);
        
        seed_id = 0;
        seed_slc_size = 1;
        while(3) {//种子扩展比对过程+++++++++++++++++++++++++++++++++++++++++++++++++

            
fprintf(stderr, "seed_id = %u, seed_slc_size = %u\n", seed_id, seed_slc_size);
if(seed_id >= seed_slc_size && seed_id != MAX_SLC){
    fprintf(stderr, "446, read_id = %u unmapped\n", read_id);
    fprintf(stderr, "sTree->len_buf = %u\n", sTree->len_buf);
    fprintf(stderr, "num = %u\n", sTree->back_buf[0][1]);
    //return;
}    
            if( seed_id >= seed_slc_size){
                //输出read_id比对结果的信息
                //seed_id的大小判定，比对失败与成功
                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                if(query->n_ambiguous == 1) aln_result_sum[5]++;
                if(query->n_ambiguous == 2) aln_result_sum[6]++;
                if(query->n_ambiguous >2) aln_result_sum[7]++;
                if(seed_id== MAX_SLC) {
                    aln_result_sum[0]++;
                } else  {
   
                    if(query->n_ambiguous == 0) aln_result_sum[1]++;
                    else if( query->n_ambiguous < 2) aln_result_sum[2]++;
                    else aln_result_sum[3]++;
                    printf("read_id = %u, cls = %u, seed_slc_num = %u\n", read_id, sTree->cls, seed_slc_num);               
                 
                    aln_result_sum[8+max_stre_cls]++;
                    if(seed_slc_num +13> 98) aln_result_sum[99]++;
                    else aln_result_sum[13+seed_slc_num]++;
                    

                }
                break; 
            }
            
        
//++++++++++++++++++++++++++++++++++++++++++++++
//以下临时注释
/*
            query_t *q = query_buf+read_id;
            f_read_seq = q->seq;
            getSeed(seed_seq,f_read_seq,seed_off);
*/          
//++++++++++++++++++++++++++++++++++++++++++++
//生成带变异的reads

            seed->id = seed_id;
            //num = aln_seed_seq(fm_idx, hash_boundry, read_len,  f_read_seq, r_read_seq, qual, r_qual, seed );

//fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
            num = aln_seed_seq_1(fm_idx, hash_boundry, query,  seed );

            seed_id = seed->id;
            seed_slc_size = seed->slc_size;                   

            if(num ==0) continue;
fprintf(stderr, "%u, seed_id = %u, num = %u, find_num = %u\n", __LINE__, seed_id, num, seed->find_num);
            seed_slc_num++;
            bgn = seed->bgn;
            end = seed->end;
            num = seed->num;
  

            //选定的read序列指针赋值给read_seq
            if(seed->slc[seed_id].drct == 0){
                //read_seq = f_read_seq;
                read_seq = query->seq;
            } else{//slc[seed_id].drct ==1
                //read_seq = r_read_seq;
                read_seq = query->rseq;
            }
            read_len = query->l_seq;
            seed_off = -seed->slc[seed_id].s_off;
            old_seed_id = seed_id;
            
            int L_offset = read_len/2+seed_off;
            int R_offset = read_len/2+SEED_LEN+seed_off; 

            if(L_offset >(read_len-R_offset)){
                NUM_EXT_SEED = (read_len-R_offset)/16; 
            } else{
                NUM_EXT_SEED = L_offset/16;
            }

            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            uint32_t j_idx = bgn;
            uint32_t pos_i = 0;
            buf_algn[0] = SEED_LEN;
            buf_algn[1] = bgn;
            buf_algn[2] = num;
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            //小规模idx区间直接通过sw方法比对
            if(num <= IS_SMLSIZ){
                j_idx = bgn;
                pos_i = 0;
                for(j_idx=bgn;j_idx<=end;j_idx++){
                    uint32_t pos = bwt_sa(fm_idx->bwt, j_idx);
                    pos_buf[pos_i++] = pos;
                }
                Flg_Algn = AlgnPos(fm_idx, 0, read_seq, read_len, seed_off, pos_buf, pos_i, sub);  
                if(Flg_Algn>0) {
                    //OutAlgnInfo();
        //printf("qname = %s, alned to pos %u\n", query->name, pos_buf[pos_i++]);
                    seed_id = MAX_SLC;//比对成功
                }else{
                    seed_id++;
                }
                continue;
            }// end if(num<=IS_SMLSIZ)++++++++++++++++++++++++++++++++++++++++++++
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            //以下利用多级索引进行比对
            // 函数定义：int getCapPos(JmpMod *jmp，uint32_t buf_algn[]);
            // 其中buf_algn的[0]是输入，种子长度；[1]是输入BgnIdx,[2]是输入Num
            // [3]是返回值capidx或者smpos数组的行号，
            getCapPos(jmp,buf_algn);
            cap_row = buf_algn[3];
            sub->pair_out[0][0] = 1;
            sub->pair_out[0][1] = 1;
            sub->pair_out[1][0] = buf_algn[3];
            sub->pair_out[1][1] = num;
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            sTree->cls = 0;
            sTree->stck_arry[0][0] = buf_algn[3];
            sTree->stck_arry[0][1] = num;
            sTree->stck_arry[0][2] = 0;
            sTree->stck_arry[0][3] = sTree->cls;
            //sTree->len_arry = 0; 
            sTree->len_arry = 1; 
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
            sTree->cls = 0;
            sTree->back_buf[0][0] = buf_algn[3];
            sTree->back_buf[0][1] = num;
            sTree->back_buf[0][2] = 0;
            sTree->back_buf[0][3] = sTree->cls;
            //sTree->len_arry = 0; 
            sTree->len_buf = 1; 
  
 
 
            //以下循环是多级索引中比对			
/*     
fprintf(stderr, "------------------------------\n");
int ii;
for(ii = 0; ii < sTree->len_arry; ++ii){
    fprintf(stderr, "sTree->stck_arry[%u] = ", ii);
    fprintf(stderr, "  %2u  ", sTree->stck_arry[ii][0]);
    fprintf(stderr, "  %2u  ", sTree->stck_arry[ii][1]);
    fprintf(stderr, "  %2u  ", sTree->stck_arry[ii][2]);
    fprintf(stderr, "  %2u  ", sTree->stck_arry[ii][3]);
    fprintf(stderr, "  %2u  ", sTree->stck_arry[ii][4]);
    fprintf(stderr, "\n");
}
fprintf(stderr, "sTree->cls = %u\n", sTree->cls);
fprintf(stderr, "cap_row = %u, sTree->cls = %u, sTree->len_arry=%u\n", cap_row, sTree->cls, sTree->len_arry);
*/
 

            int aln_err_flg = 0; 
            int aln_ext_flg = 0; 
            uint32_t (*ext_idx)[2];



            while(4){


/*  
                if(aln_ext_flg == 0){
                    memcpy(sTree->back_buf+sTree->len_buf, sTree->stck_arry[sTree->len_arry], sizeof(uint32_t)*4); 
                    sTree->back_buf[sTree->len_buf][4] = seed_id;
                    sTree->back_buf[sTree->len_buf][5] = aln_err_flg;
                }
*/
                aln_ext_flg = 0;
                if(sTree->len_arry == 0){ break;} //while(4)+++ 
                
                sTree->cls = sTree->stck_arry[sTree->len_arry-1][3];
            
                if(sTree->cls >max_stre_cls) max_stre_cls = sTree->cls;
                if(sTree->cls == NUM_EXT_SEED-1){
                    seed_id = MAX_SLC;
                    break; //while(4)+++ 
                }
                cap_row = sTree->stck_arry[sTree->len_arry-1][0];
                --sTree->len_arry;
                //下面的函数赋值当前块索引数据
                //以下扩展序列段数据赋值
                L_offset = read_len/2+seed_off-sTree->cls*16-16;
                R_offset = read_len/2+SEED_LEN+seed_off+sTree->cls*16; 
fprintf(stderr, "%u, sTree->cls = %u, L_offset = %u, R_offset = %u\n", __LINE__, sTree->cls, L_offset, R_offset);
fprintf(stderr, "sTree->len_arry = %u, stck_arry[len][1] = %u\n", sTree->len_arry, sTree->stck_arry[sTree->len_arry][1]); 
                memcpy(sub->ext_seqL, read_seq+L_offset, 16);
                memcpy(sub->ext_seqR, read_seq+R_offset, 16);

           
                setBlckData(eBlck,sTree->cls,cap_row, sub, &cB);
                //以下smBWT左段比对+++++++++++++++++++++++++++++++++++++++++++++++++++
                Flg_Algn = AlgnmExtSeq(cB,sub,0);
                if(Flg_Algn == 0){
                    //aln_ext_flg = 0;
                    aln_err_flg = 1;
                    continue; //while(4)+++ 
                }
                //以下smBWT右段比对+++++++++++++++++++++++++++++++++++++++++++++++++++
                Flg_Algn = AlgnmExtSeq(cB,sub,1);		
                if(Flg_Algn == 0){ 
                    //aln_ext_flg = 0;
                    aln_err_flg = 2;
                    continue; //while(4)+++  
                }
//test_smbwt(cB, 0, sub->ext_seqL); 
//test_smbwt(cB, 1, sub->ext_seqR); 
                uint32_t nxt_pnt,nxt_end,pos,i,j ;
                uint8_t  nxt_flg ;
                uint32_t num_smpos;
                uint32_t *nxt_smpos;
                int Flg_Appr = 0;
                //以下Flg_Algn >0时，即左右都存在bwt的比对序列时运行
fprintf(stderr, "%u\n", __LINE__);        
                Flg_Algn = PairExtSeq_all(cB,sub);

fprintf(stderr, "%u\n", __LINE__);        
                if(Flg_Algn == 0){
                    Algn_apro_all(cB, sub, 1);  
                    Flg_Appr++;
                    Flg_Algn = PairExtSeq_all(cB,sub);
                }

fprintf(stderr, "%u\n", __LINE__);        
                if(Flg_Algn == 0){
                    Algn_apro_all(cB, sub, 0);  
                    Flg_Appr++;
                    Flg_Algn = PairExtSeq_all(cB,sub);
                }

fprintf(stderr, "%u\n", __LINE__);        
                if(Flg_Algn == 0) { 
                    //aln_ext_flg = 0;
                    aln_err_flg = 3;
                    continue; //while(4)+++ 
                } 

fprintf(stderr, "%u\n", __LINE__);        
                if(Flg_Algn > 0){
                    Flg_Algn = 0;         
                    num_smpos =  sub->pair_out[0][1];
                    if(num_smpos >0){
fprintf(stderr, "%u, num_smpos = %u\n", __LINE__, num_smpos);
                        Flg_Algn = AlgnPos_buf(fm_idx, sTree->cls, read_seq, read_len, seed_off, num_smpos, pos_buf, sub);  
                        
fprintf(stderr, "%u Flg_Algn = %u\n", __LINE__, Flg_Algn);
                        if(Flg_Algn>0){
                            seed_id = MAX_SLC;
                
        //printf("qname = %s, alned to pos %u\n", query->name, pos_buf[pos_i++]);
                            break;   //while(4){  //多级索引中比对
                        } 
                    }//end if(num_smpos>0)+++++++
                    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    //多级索引之间转换控制核心控制块信息更新
                    //控制块堆栈信息更新开始			
                    if(sub->pair_out[0][0]- sub->pair_out[0][1] >0){
                        ext_idx = cB->head_extidx+cB->extidx;
                        setStackTree(sTree,sub->pair_out, seed_id, ext_idx);
                        //if(sub->pair_out[0][0]- sub->pair_out[0][1] >1 || Flg_Appr >1 ) {
                        //if(sTree->len_arry >1) {
                        //    continue;
                        //} 
                        aln_ext_flg = 1;
                    }                    
                    /*  
                    if(sTree->len_arry >1) {
                        continue;//while(4)
                    }
                    */
                }// 1081:if(Flg_Algn > 0)++++++++++++++++++++++++++++++++++++++++++++++++++

fprintf(stderr, "%u\n", __LINE__);        
                if(Flg_Appr >1) { 
                    //aln_ext_flg = 0:1
                    aln_err_flg = 4;
                    continue; //while(4)+++ 
                }

fprintf(stderr, "%u\n", __LINE__);        
                if(Flg_Appr == 0) {
                    Algn_apro_all(cB, sub, 1);  
                    Flg_Appr++; 
                    Flg_Algn = PairExtSeq_all(cB,sub);

fprintf(stderr, "%u\n", __LINE__);        
                    if(Flg_Algn == 0){
                        Algn_apro_all(cB, sub, 0);  
fprintf(stderr, "%u\n", __LINE__);        
                        Flg_Appr++; 
                        Flg_Algn = PairExtSeq_all(cB,sub); 
fprintf(stderr, "%u\n", __LINE__);        
                    }

fprintf(stderr, "%u\n", __LINE__);        
                    if(Flg_Algn == 0){ 
                        aln_err_flg = 5; 
                        continue; //while(4)+++ 
                    } 
                }else {//if(Flg_App==1)
                    Algn_apro_all(cB, sub, 0);  
                    Flg_Appr++; 
                    Flg_Algn = PairExtSeq_all(cB,sub); 

fprintf(stderr, "%u\n", __LINE__);        
                    if(Flg_Algn == 0){ 
                        aln_err_flg = 6;
                        continue; //while(4)+++ 
                    } 
                }
                if(Flg_Algn > 0){
                    num_smpos =  sub->pair_out[0][1];
                    Flg_Algn = 0;         
                    if(num_smpos >0){
                        Flg_Algn = AlgnPos_buf(fm_idx, sTree->cls, read_seq, read_len, seed_off, num_smpos, pos_buf, sub);  
                        if(Flg_Algn>0){
                            seed_id = MAX_SLC;
                    
        //printf("qname = %s, alned to pos %u\n", query->name, pos_buf[pos_i++]);
                            break;   //while(4){  //多级索引中比对
                        } 
                    }//end if(num_smpos>0)+++++++
                    //多级索引之间转换控制核心控制块信息更新
                    //控制块堆栈信息更新开始			

fprintf(stderr, "%u\n", __LINE__);        
                    if(sub->pair_out[0][0]- sub->pair_out[0][1] >0 ){
                        ext_idx = cB->head_extidx+cB->extidx;
                        setStackTree(sTree,sub->pair_out, seed_id, ext_idx);
                        //setStackTree(sTree,sub->pair_out);
                        //if(sub->pair_out[0][0]- sub->pair_out[0][1] >1 || Flg_Appr >1 ) {
                        aln_ext_flg = 1;
                    }
                    /*  
                    if(sTree->len_arry >1) {
                        continue; //while(4)
                    }
                    */
                }                
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(Flg_Appr == 0) {
    printf("%u\n", __LINE__);
    exit(1);
}

fprintf(stderr, "%u\n", __LINE__);        
                if(Flg_Appr >1) { 
                    aln_err_flg = 7; 
                    Flg_Algn = 0;
                    continue; //while(4)+++ 
                } //while(4)+++ 

fprintf(stderr, "%u\n", __LINE__);        
                if(Flg_Appr == 1) {
                    Algn_apro_all(cB, sub, 0);  
                    Flg_Appr++; 
                    Flg_Algn = PairExtSeq_all(cB,sub); 
                    if(Flg_Algn == 0){ 
                        aln_err_flg = 8; 
                        continue; //while(4)+++  
                    } 
                }
                num_smpos =  sub->pair_out[0][1];
                Flg_Algn = 0;         

fprintf(stderr, "%u\n", __LINE__);        
                if(num_smpos >0){

fprintf(stderr, "%u read_len = %u\n", __LINE__, read_len);
                    Flg_Algn = AlgnPos_buf(fm_idx, sTree->cls, read_seq, read_len, seed_off, num_smpos, pos_buf, sub);  

fprintf(stderr, "%u read_id = %u, seed_id = %u\n", __LINE__, read_id, seed_id);
    if(Flg_Algn>0){

fprintf(stderr, "%u read_id = %u, seed_id = %u\n", __LINE__, read_id, seed_id);
                    seed_id = MAX_SLC;
                        
        //printf("qname = %s, alned to pos %u\n", query->name, pos_buf[pos_i++]);
                        break;   //while(4){  //多级索引中比对
                    } 
                   
                }//end if(num_smpos>0)+++++++
                //多级索引之间转换控制核心控制块信息更新
                //控制块堆栈信息更新开始			

fprintf(stderr, "%u\n", __LINE__);        
                if(sub->pair_out[0][0]- sub->pair_out[0][1] >0 ){
                    ext_idx = cB->head_extidx+cB->extidx;
                    setStackTree(sTree,sub->pair_out, seed_id, ext_idx);
                      
                    //setStackTree(sTree,sub->pair_out);
                    aln_ext_flg = 1;
                    continue;
                }

fprintf(stderr, "%u\n", __LINE__);        
            } // End: while(4) //一个read的扩展比对结束+++++++++

fprintf(stderr, "%u\n", __LINE__);        
            if(seed_id < MAX_SLC) seed_id++;
        }// End : while(3) 一次读入read序列集合的比对结束+++++++
        if(seed_id == MAX_SLC) { 
            //read边界扩展处理 
        } else{//if(seed_id == seed_slc_size)
            //非种子扩展比对方法处理 
        }
    } // End: while(2)+++++++++++++++++++++++++++++++++++++++++++     
    //for(i= 0; i < n_seqs; ++i) query_destroy(query_buf+i); 
if(REAL_FLAG==0)  break;//while(1)++++++++++
    }//End while(1) 主循环++++++++++++++++++++++++++++++++++++++++
printf("aln_sum =  ");
    for(i= 0;i < 100; ++i) {
        printf("%u :  %u  \n", i, aln_result_sum[i]);
        printf("\n");
        if(i == 7) printf("------------------------------\n");
        if(i == 13) printf("------------------------------\n");

    }
printf("\n");
printf("----------------------\n");


    free(query_buf);

  
    //kseq_destroy(seq);
    //gzclose(fp);
    return 1;

} //End ：main() +++++++++++++++++++++++++++++++++++++++++++ 
void set_read_data(idx_t *fm_idx, uint32_t read_id, uint8_t *err_mdl, int read_len, uint8_t *f_read_seq, uint8_t *r_read_seq, uint8_t *qual, uint8_t *r_qual)
{
    uint8_t *cur_err_mdl = err_mdl+read_id%20;
    uint32_t pos = read_id;
    int i; 
    for(i= 0; i < read_len; ++i){
        f_read_seq[i] = __get_pac(fm_idx->pac, pos+i); 
    }
    fprintf(stderr, "origin:\n");
    for(i=0; i <read_len; ++i) fprintf(stderr, "%u", f_read_seq[i]);
    fprintf(stderr, "\n");
    for(i= 0; i < read_len; ++i){ 
        r_read_seq[read_len-i-1] = 3-f_read_seq[i];
        qual[read_len-i-1] = r_qual[i];
    }

    for(i = 0; i < read_len; ++i ){
        uint8_t c  = f_read_seq[i];
        if(cur_err_mdl[i] > 0) {
            f_read_seq[i] = (f_read_seq[i]+cur_err_mdl[i])%4;
            fprintf(stderr, "err[%u] = %u->%u\n", i, c, f_read_seq[i]);
        }
       
    }
    fprintf(stderr, "var:\n");
    for(i=0; i <read_len; ++i) fprintf(stderr, "%u", r_read_seq[i]);
    fprintf(stderr, "\n");
    return;
}
void set_read_data_1(idx_t *fm_idx, uint32_t read_id, uint8_t *err_mdl, query_t *query)
{
    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq; 
    uint8_t *qual = query->qual; 
    //uint8_t *r_qual = query->qual; 
    uint8_t *cur_err_mdl = err_mdl+read_id%20;
    uint32_t pos = read_id;
    int i; 
    for(i= 0; i < read_len; ++i){
        f_read_seq[i] = __get_pac(fm_idx->pac, pos+i); 
    }
/*  
    fprintf(stderr, "origin:\n");
    for(i=0; i <read_len; ++i) fprintf(stderr, "%u", f_read_seq[i]);
    fprintf(stderr, "\n");
*/
    for(i= 0; i < read_len; ++i){ 
        r_read_seq[read_len-i-1] = 3-f_read_seq[i];
        //r_qual[read_len-i-1] = qual[i];
    }

    for(i = 0; i < read_len; ++i ){
        uint8_t c  = f_read_seq[i];
        if(cur_err_mdl[i] > 0) {
            f_read_seq[i] = (f_read_seq[i]+cur_err_mdl[i])%4;
            fprintf(stderr, "err[%u] = %u->%u\n", i, c, f_read_seq[i]);
        }
       
    }
/*  
    fprintf(stderr, "var:\n");
    for(i=0; i <read_len; ++i) fprintf(stderr, "%u", r_read_seq[i]);
    fprintf(stderr, "\n");
*/
    return;
}

void set_err_mdl(uint8_t err_mdl[]){
    int i;
    for(i=0; i < READ_LEN; ++i) {
        int j = ((i+19)/20)%4;
        if(i%20 ==5) err_mdl[i+j] =((i+19)/20)%4; 
    }   
    /*  
    for(i=89; i < 111; ++i) {
        err_mdl[i] = 0; 
    }
    */
    /*  
    fprintf(stderr, "0-90\n");
    for(i=0; i < 90; ++i) {
        if(err_mdl[i] > 0) fprintf(stderr, "%u, %u\n", i, err_mdl[i]);
    }
    fprintf(stderr, "90-130\n");
    for(i=90; i < 130; ++i) {
        if(err_mdl[i] > 0)  fprintf(stderr, "%u, %u\n", i, err_mdl[i]);
    }
    fprintf(stderr, "130-200\n");
    for(i=130; i < 200; ++i) {
        if(err_mdl[i] > 0)  fprintf(stderr, "%u, %u\n", i, err_mdl[i]);
    }
    */
    return;
}
