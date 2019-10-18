// 循环控制信息--------------------------------------------
//uint32_t queue_align_info[][2];  //开辟较大内存
//uint32_t this_stack_tree[NUM_EXT][5];  
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
#define READ_LEN 200
#define FILE_PREFIX "../../BWT_INDEX/prefix"
const char *read_file = "../../Reads/Read.fq";
#define MAX_READS 5
#define MIN_Q 255
typedef struct{
    int slc_row;//该mer8中在aln_mer8中的行号
    int err_pos;//错误匹配出现的位置；
    int err_val;//错误匹配的碱基
} aln_mer8_t;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(){	
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
	    perror("error:[malloc(NUM_EXT*sizeof(FileName))]");
	    exit(1);
	}
    setFileName(fn);
    InitIdxsArry(fn, &sub, eBlck,&sTree,cap,&jmp);

    int i, j;
    uint32_t hash_boundry[12] = {};
    uint8_t seed[12] = {};
    for(i =0; i < 12; ++i){
        seed[i] = __get_pac(fm_idx->pac, fm_idx->bwt->seq_len-12+i);
    }
    uint32_t k = 0, l = 0;
    for(i =0; i < 11; ++i){
        int n = bwt_match_exact_alt(fm_idx->bwt, 1, seed+12-1-i, &k, &l);
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
    int  seed_off;
	uint32_t get_num;
	uint32_t read_id;
    uint32_t buf_algn[20]; 	
	int      Flg_Algn  = 0 ;
	uint32_t cur_read_num = 0;

    queryio_t *qs = query_open(read_file);;
    query_t *multiSeqs = (query_t *)calloc(MAX_READS, sizeof(query_t));
    int read_num = 0;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//以下代码是对参考基因组中获取read数据进行变异错误的模版
uint8_t err_mdl[READ_LEN+20] = {};
for(i=0; i < READ_LEN; ++i) {
    int j = ((i+19)/20)%4;
    if(i%20 ==5) err_mdl[i+j] =((i+19)/20)%4; 
}
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
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    fprintf(stderr, "NUM_EXT = %u\n", NUM_EXT);
    fprintf(stderr, "IS_SMLSIZ = %u\n", IS_SMLSIZ);
    fprintf(stderr, "READ_LEN = %u\n", READ_LEN);
    fprintf(stderr, "SEED_LEN = %u\n", SEED_LEN);
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //以下代码段是种子相关变量初始化 
    uint32_t pos_buf[IS_SMLSIZ+1];
    int SEED_SLC_SIZE = 16;//应该满足SEED_SLC_SIZE>seed_slc_size
    int MAX_SLC = SEED_SLC_SIZE+2;
    int seed_slc_size, seed_slc_size_i;
    int seed_id;
    int find_num = 0;
    seed_slc_t *slc = calloc(SEED_SLC_SIZE, sizeof(seed_slc_t));
    seed_slc_t *slc_i = calloc(SEED_SLC_SIZE, sizeof(seed_slc_t));
    seed_slc_t *slc_b = calloc(SEED_SLC_SIZE, sizeof(seed_slc_t));
    aln_mer8_t *aln_8 = calloc(8*3, sizeof(aln_mer8_t));
    seed_init(slc_i);
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //以下代码段是read相关变量定义及其常用变量    
    uint32_t read_len = READ_LEN;
    uint8_t ch[READ_LEN], r_ch[READ_LEN];
    uint8_t qual[READ_LEN] = {}, r_qual[READ_LEN] = {}; 
    uint32_t bgn,end,num,bgn_row,end_row,cap_row, old_seed_id;
    uint8_t aln_rst_flg = 0;
    uint32_t aln_result_sum[3] = {};
    
    while(1){ //主循环++++++++++++++++++++++++++++++++++++++++
        // 循环过程中每次读入规定个数的reads
        uint32_t n_seqs = query_read_multiSeqs(qs,MAX_READS, multiSeqs);
        if(n_seqs <1) break;
//++++++++++++++++++++++++++++++++
//测试代码
        //seed_off = 0;
        read_id = 0;
        n_seqs = 100000;
        seed_id = 0;
//++++++++++++++++++++++++++++++++
        //一次读入的read数组循环处理
        while(read_id < n_seqs) {

            
        
//++++++++++++++++++++++++++++++++++++++++++++++
//以下临时注释
/*
            query_t *q = multiSeqs+read_id;
            f_read_seq = q->seq;
            getSeed(seed_seq,f_read_seq,seed_off);
*/          
//++++++++++++++++++++++++++++++++++++++++++++
//生成带变异的reads
if(seed_id ==0)
{
    uint8_t *cur_err_mdl = err_mdl+read_id%20;
    uint32_t pos = read_id;
    for(i= 0; i < read_len; ++i){
        ch[i] = __get_pac(fm_idx->pac, pos+i); 
    }
    fprintf(stderr, "origin:\n");
    for(i=0; i <read_len; ++i) fprintf(stderr, "%u", ch[i]);
    fprintf(stderr, "\n");


    f_read_seq = ch;
    r_read_seq = r_ch;
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
    /*  
    //j = read_len/2-SEED_LEN/2+read_id%SEED_LEN;
    j = read_id%READ_LEN;
    uint8_t err = read_id%3+1;
    fprintf(stderr, "\nedit pos = %u, err = %u, ch = %u\n", j, (ch[j]+err)%4, ch[j]);
    ch[j] = (ch[j] +err)%4;
    //r_read_seq = ch;
    */
    fprintf(stderr, "var:\n");
    for(i=0; i <read_len; ++i) fprintf(stderr, "%u", r_read_seq[i]);
    fprintf(stderr, "\n");
}
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            fprintf(stderr, "seed_id = %u\n", seed_id);
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
                for(i = 0; i < SEED_SLC_SIZE; ++i){
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
                find_num = 0; 
            } //end if(seed_id == 0)  
           //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(seed_id == seed_slc_size){
    fprintf(stderr, "446, read_id = %u unmapped\n", read_id);
    return;
}    
            if( seed_id >= seed_slc_size){
printf("472, read_id = %u\n", read_id);
                //输出read_id比对结果的信息
                //seed_id的大小判定，比对失败与成功
                //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                if(seed_id== MAX_SLC+1) {
                    aln_result_sum[2]++;
                } else if(seed_id == MAX_SLC) {
                    aln_result_sum[1]++;
                } else{
                    printf("485, read_id = %u, seed_id= %u, seed_slc_size = %u\n", read_id, seed_id, seed_slc_size);                
                    aln_result_sum[0]++;
                }
                read_id++; 
                aln_rst_flg = 0;
                seed_id = 0;
fprintf(stderr, "===================================\n", read_id);
fprintf(stderr, "read_id = %u\n", read_id);
                continue;
            }
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            //确定第seed_id新的种子序列           
            int cur_row; 
            uint32_t seq12;
            uint8_t *pseed, seed_buf[SEED_LEN]; 
            int seq_off = read_len/2-slc_i[seed_id].s_off;
            seed_t seed;
            int seed_flg;
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            //种子序列指针赋值
            if(slc[seed_id].drct == 0){
                pseed = f_read_seq+seq_off;
            } else{//slc[seed_id].drct ==1
                pseed = r_read_seq+seq_off;
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
                    seed.hash_12mer = lkt_seq2LktItem(pseed, 8, 19);
                    seed.seed_8mer_L = pseed;
                    seed.seed_seq = pseed;
                    find_num = get_seed_8mer_bwt_L(fm_idx, hash_boundry, &seed);
                    if(find_num == 0) {
                        seed_id++;
                        continue;
                    }
for(i = 0; i < find_num; ++i){
    fprintf(stderr, "err_pos = %u,alt=%u,  bgn=%u, end = %u\n", seed.aln_buf[i][0], seed.aln_buf[i][1], seed.aln_buf[i][2],                     seed.aln_buf[i][3]); 
}
            
                    bgn = seed.aln_buf[cur_row][2]; 
                    end = seed.aln_buf[cur_row][3]; 
                    num = end+1- bgn;
                    seed_flg = 1;
                    if(find_num == 1) { 
                        find_num = 0; 
                    }

fprintf(stderr, "%u, seed_id=%u, bgn = %u, end=%u\n", __LINE__, seed_id, bgn, end);
                } else{//slc[seed_id] ==2
                    //调用右侧模糊比对查找函数,获得模糊比对序列集合数组
                    int bg_i = 0, ed_i = 8, num_seq;
                    num_seq = (ed_i-bg_i)*3; 
                    seed.hash_12mer = lkt_seq2LktItem(pseed, 0, 11);
                    seed.seed_8mer_R = pseed+12;
                    seed.seed_seq = pseed;
                    int var_seq_num = gen_8mer_var_buf(&seed, bg_i, ed_i, 1); 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
{
    fprintf(stderr, "%u, var_seq_num = %d\n", __LINE__, var_seq_num);
    for(k = 0; k <SEED_LEN; k++){
        fprintf(stderr, "%u ", pseed[k]); 
    } 
    fprintf(stderr, "-------\n"); 
    for(i = 0; i < num_seq; ++i){
        uint32_t seq = seed.sort_buf_R[i]; 
        int k;
        for(k = 14; k >=0; k-=2){
            fprintf(stderr, "%u ", (seq>>k)&3); 
        } 
        fprintf(stderr, "\t%u\n", seq); 
    }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    find_num = get_seed_8mer_seq_R(fm_idx, &seed, num_seq);
fprintf(stderr, "%u, find_num = %d\n", __LINE__, find_num);
                    if(find_num == 0) {
                        seed_id++;
                        continue;
                    }
                    memcpy(seed_buf, pseed, SEED_LEN); 
{
    for(i =0; i< find_num; ++i){
        uint8_t f_i = seed.find_buf_R[i];  
        uint8_t var_pos  = seed.var_info[f_i][0];
        uint8_t var_ch   = seed.var_info[f_i][1];
        uint8_t var_ch_o = seed_buf[var_pos+12]; 
        fprintf(stderr, "%u, f_i = %u, var_pos = %u, var_ch = %u, var_ch_o = %u\n", __LINE__, f_i, var_pos, var_ch, var_ch_o);
    }

}
                    cur_row = 0;
                    uint8_t f_i = seed.find_buf_R[cur_row];  
                    uint8_t var_pos  = seed.var_info[f_i][0];
                    uint8_t var_ch   = seed.var_info[f_i][1];
                    uint8_t var_ch_o = seed_buf[var_pos+12]; 
fprintf(stderr, "%u, f_i = %u, var_pos = %u, var_ch = %u, var_ch_o = %u\n", __LINE__, f_i, var_pos, var_ch, var_ch_o);
                    seed_buf[var_pos+12] = var_ch; 
                    seq12 = lkt_seq2LktItem(seed_buf, 8, 19);
                    k = fm_idx->fastmap->item[seq12];
                    l = fm_idx->fastmap->item[seq12+1]-1;
                    l -= get_12mer_correct(hash_boundry, l);
                    num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, seed_buf, &k, &l);
                    bgn = k; end = l;
fprintf(stderr, "%u, num = %u, bgn = %u, end = %u\n", __LINE__, num, bgn, end);
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
                    bgn = seed.aln_buf[cur_row][2]; 
                    end = seed.aln_buf[cur_row][3]; 
                    num = end+1- bgn;
                    //seed_flg = 1;
                } else if(seed_flg ==2) {//第二类近似比对继续
                    seed_id--;
                    cur_row++; 
                    
                    uint8_t f_i = seed.find_buf_R[cur_row];  
                    uint8_t var_pos  = seed.var_info[f_i][0];
                    uint8_t var_ch   = seed.var_info[f_i][1];
                    uint8_t var_ch_o = seed_buf[var_pos+12]; 
                    seed_buf[var_pos+12] = var_ch; 
                    seq12 = lkt_seq2LktItem(seed_buf, 8, 19);
                    k = fm_idx->fastmap->item[seq12];
                    l = fm_idx->fastmap->item[seq12+1]-1;
                    l -= get_12mer_correct(hash_boundry, l);
                    num = bwt_match_exact_alt(fm_idx->bwt, SEED_LEN-12, seed_buf, &k, &l);
                    bgn = k; 
                    end = l;
fprintf(stderr, "%u, num = %u, bgn = %u, end = %u\n", __LINE__, num, bgn, end);
                    seed_buf[var_pos+12] = var_ch_o; 
                }
                if(find_num - cur_row == 1) {//cur_row最后一个行的情况 
                    find_num = 0; 
                }
            }
    

            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            //选定的read序列指针赋值给read_seq
            if(slc[seed_id].drct == 0){
                read_seq = f_read_seq;
            } else{//slc[seed_id].drct ==1
                read_seq = r_read_seq;
            }
            seed_off = -slc[seed_id].s_off;
            old_seed_id = seed_id;

fprintf(stderr, "------------------------------------------------\n");
fprintf(stderr, "read_id = %u, seed_id = %u, strand = %u\n", read_id, seed_id, slc[seed_id].drct);
fprintf(stderr, "find_num = %u, cur_row = %u\n", find_num, cur_row);
fprintf(stderr, "bgn = %u, end =%u, num = %u\n", bgn, end, num); 
fprintf(stderr, "seed_flag = %u\n", seed_flg);
fprintf(stderr, "seed_off = %u, hash_off = %u\n", slc[seed_id].s_off, slc[seed_id].h_off);
fprintf(stderr, "bgn = %u, end =%u, num = %u\n", bgn, end, num); 
/*  
fprintf(stderr, "forward:\n");
for(i=74; i <126; ++i) fprintf(stderr, "%u", f_read_seq[i]);
fprintf(stderr, "\n");
fprintf(stderr, "reverse:\n");
for(i=74; i <126; ++i) fprintf(stderr, "%u", r_read_seq[i]);
fprintf(stderr, "\n");
fprintf(stderr, "read_seq:\n");
for(i=74; i <126; ++i) fprintf(stderr, "%u", read_seq[i]);
fprintf(stderr, "\n");
fprintf(stderr, "seq_off = %d\n", seq_off); 
fprintf(stderr, "seed_off = %d\n", seed_off); 
*/

            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            uint32_t j_idx = bgn;
            uint32_t pos_i = 0;
            buf_algn[0] = SEED_LEN;
            buf_algn[1] = bgn;
            buf_algn[2] = num;
            if(num < 1){//种子序列比对失败
                seed_id++;
                continue; 
            }
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
 

            struct ExtBlck *cB;
            while(3){
                if(sTree->len_arry == 0){ break;}
                sTree->cls = sTree->stck_arry[sTree->len_arry-1][3];
                if(sTree->cls == NUM_EXT-1){
                    seed_id = MAX_SLC;
                    break; 
                }
                cap_row = sTree->stck_arry[sTree->len_arry-1][0];
                --sTree->len_arry;
                //下面的函数赋值当前块索引数据
                //以下扩展序列段数据赋值
                int L_offset = read_len/2+seed_off-sTree->cls*16-16;
                int R_offset = read_len/2+SEED_LEN+seed_off+sTree->cls*16; 
                memcpy(sub->ext_seqL, read_seq+L_offset, 16);
                memcpy(sub->ext_seqR, read_seq+R_offset, 16);
                setBlckData(eBlck,sTree->cls,cap_row, sub, &cB);
                //以下smBWT左段比对+++++++++++++++++++++++++++++++++++++++++++++++++++
fprintf(stderr, "%u\n", __LINE__);
                Flg_Algn = AlgnmExtSeq(cB,sub,0);
fprintf(stderr, "%u\n", __LINE__);
                if(Flg_Algn == 0){ continue;} 
fprintf(stderr, "%u\n", __LINE__);
                //以下smBWT右段比对+++++++++++++++++++++++++++++++++++++++++++++++++++
                Flg_Algn = AlgnmExtSeq(cB,sub,1);		
                if(Flg_Algn == 0){ continue; } 
//test_smbwt(cB, 0, sub->ext_seqL); 
//test_smbwt(cB, 1, sub->ext_seqR); 
                uint32_t nxt_pnt,nxt_end,pos,i,j ;
                uint8_t  nxt_flg ;
                uint32_t num_smpos;
                uint32_t *nxt_smpos;
                int Flg_Appr = 0;
                //以下Flg_Algn >0时，即左右都存在bwt的比对序列时运行
fprintf(stderr, "1099, seqL_out = %u, %u, %u\n", sub->seqL_out[0][0], sub->seqL_out[0][1], sub->seqL_out[0][2]);
for(i = 1; i <= sub->seqL_out[0][0]; ++i){
    fprintf(stderr, "bg_idx = %u, err_pos= %u\n", sub->seqL_out[i][0], sub->seqL_out[i][1]);
}
fprintf(stderr, "1104, seqL_out = %u, %u, %u\n", sub->seqR_out[0][0], sub->seqR_out[0][1], sub->seqR_out[0][2]);
for(i = 1; i <= sub->seqR_out[0][0]; ++i){
    fprintf(stderr, "bg_idx = %u, err_pos = %u\n", sub->seqR_out[i][0], sub->seqR_out[i][1]);
}
                
                Flg_Algn = PairExtSeq_1(cB,sub);
                if(Flg_Algn == 0){
                    Algn_apro_all(cB, sub, 1);  
                    Flg_Appr++;
                    Flg_Algn = PairExtSeq_1(cB,sub);
                }
                if(Flg_Algn == 0){
                    Algn_apro_all(cB, sub, 0);  
                    Flg_Appr++;
                    Flg_Algn = PairExtSeq_1(cB,sub);
                }
                if(Flg_Algn == 0) { continue;}
                if(Flg_Algn > 0){
                    Flg_Algn = 0;         
                    num_smpos =  sub->pair_out[0][1];
                    if(num_smpos >0){
                        Flg_Algn = AlgnPos_buf(fm_idx, sTree->cls, read_seq, read_len, seed_off, num_smpos, pos_buf, sub);  
                        if(Flg_Algn>0){
                            seed_id = MAX_SLC;
                            break;   //while(3){  //多级索引中比对
                        } 
                    }//end if(num_smpos>0)+++++++
                    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    //多级索引之间转换控制核心控制块信息更新
                    //控制块堆栈信息更新开始			
                    if(sub->pair_out[0][0]- sub->pair_out[0][1] >0){
                        setStackTree(sTree,sub->pair_out);
                        //if(sub->pair_out[0][0]- sub->pair_out[0][1] >1 || Flg_Appr >1 ) {
                        //if(sTree->len_arry >1) {
                        //    continue;
                        //} 
                    } 
                }// 1081:if(Flg_Algn > 0)++++++++++++++++++++++++++++++++++++++++++++++++++
                if(Flg_Appr >1) { continue;}
                if(Flg_Appr == 0) {
                    Algn_apro_all(cB, sub, 1);  
                    Flg_Appr++; 
                    Flg_Algn = PairExtSeq_1(cB,sub);
                    if(Flg_Algn == 0){
                        Algn_apro_all(cB, sub, 0);  
                        Flg_Appr++; 
                        Flg_Algn = PairExtSeq_1(cB,sub); 
                    }
                    if(Flg_Algn == 0){ continue;}
                }else {//if(Flg_App==1)
                    Algn_apro_all(cB, sub, 0);  
                    Flg_Appr++; 
                    Flg_Algn = PairExtSeq_1(cB,sub); 
                    if(Flg_Algn == 0){ continue;}
                }
                if(Flg_Algn > 0){
                    num_smpos =  sub->pair_out[0][1];
                    Flg_Algn = 0;         
                    if(num_smpos >0){
                        Flg_Algn = AlgnPos_buf(fm_idx, sTree->cls, read_seq, read_len, seed_off, num_smpos, pos_buf, sub);  
                        if(Flg_Algn>0){
                            seed_id = MAX_SLC;
                            break;   //while(3){  //多级索引中比对
                        } 
                    }//end if(num_smpos>0)+++++++
                    //多级索引之间转换控制核心控制块信息更新
                    //控制块堆栈信息更新开始			
                    if(sub->pair_out[0][0]- sub->pair_out[0][1] >0 ){
                        setStackTree(sTree,sub->pair_out);
                        //if(sub->pair_out[0][0]- sub->pair_out[0][1] >1 || Flg_Appr >1 ) {
                        //if(sTree->len_arry >1) {
                            //continue;
                        //}
                    }
                }                
                //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(Flg_Appr == 0) {
    printf("%u\n", __LINE__);
    exit(1);
}
                if(Flg_Appr >1) { continue;}
                if(Flg_Appr == 1) {
                    Algn_apro_all(cB, sub, 0);  
                    Flg_Appr++; 
                    Flg_Algn = PairExtSeq_1(cB,sub); 
                    if(Flg_Algn == 0){ continue; }
                }
                num_smpos =  sub->pair_out[0][1];
                Flg_Algn = 0;         
                if(num_smpos >0){
                    Flg_Algn = AlgnPos_buf(fm_idx, sTree->cls, read_seq, read_len, seed_off, num_smpos, pos_buf, sub);  
                    if(Flg_Algn>0){
                        seed_id = MAX_SLC;
                        break;   //while(3){  //多级索引中比对
                    } 
                   
                }//end if(num_smpos>0)+++++++
                //多级索引之间转换控制核心控制块信息更新
                //控制块堆栈信息更新开始			
                if(sub->pair_out[0][0]- sub->pair_out[0][1] >0 ){
                    setStackTree(sTree,sub->pair_out);
                    continue;
                }
            } // End: while(3) //一个read的扩展比对结束+++++++++
            if(seed_id < MAX_SLC) seed_id++;
        }// End : while(read_id < ONE_INPUT_SIZE) 一次读入read序列集合的比对结束 ++++++
        //for(i= 0; i < n_seqs; ++i) query_destroy(multiSeqs+i); 
        break;
    }//End while(1) 主循环++++++++++++++++++++++++++++++++++++++++


    printf("aln_sum = %u %u %u\n", aln_result_sum[0], aln_result_sum[1], aln_result_sum[2]);
    free(multiSeqs);
	free(slc);
    free(aln_8);
    //kseq_destroy(seq);
    //gzclose(fp);
    return 1;

} //End ：main() +++++++++++++++++++++++++++++++++++++++++++ 

