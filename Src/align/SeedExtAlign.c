    // /环控制信息--------------------------------------------
//uint32_t queue_align_info[][2];  //开辟较大内存
//uint32_t this_stack_tree[NUM_EXT_SEED][5];  
//[i][0]是第i个节点(即，第i级索引）在queue_align_info中的起始行号
//[i][1]是第i个节点(即，第i级索引）在queue_align_info中的个数1
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

extern void query_aln2sam1(query_t *query, idx_t *idx, struct SubBuf *sub, FILE *fp);
extern int query_aln2sam2(query_t *query, idx_t *idx, struct SubBuf *sub, char *str_buf);

#define SEED_LEN 20
#define SEED_FLT 0



#define FILE_PREFIX "./GRCh38"
//const char *read_file = "../../Reads/sample.fq";
//const char *read_file = "../../Reads/grch38_mason2.1.fq";
//sim data
//const char *read_file = "../../../Genome/Grch38/s1_1.fq";
//const char *read_file = "../../../Genome/Grch38/s1_1_200.fq";
//const char *read_file = "../../../Genome/Grch38/s1_1_250_1M.fq";
//const char *read_file = "../../../Genome/Grch38/s1_1_150_1M.fq";
//const char *read_file = "../../../Genome/Grch38/s1_1_100_1M.fq";
//const char *read_file = "../../../Genome/Grch38/s1_1_200_10M.fq";

//real data
//const char *read_file = "../../../Reads/downloading/real_250_100K2.fastq";
//const char *read_file = "../../../Reads/downloading/SRR826460_1.fastq";
//const char *read_file = "../../../Reads/downloading/ERR1009351_1.fastq";
const char *read_file = "./reads.fq";
//debug data
//const char *read_file = "./wrong.fa";

//const char *read_file = "./repeat_200.fa";
//const char *read_file = "./repeat_reads_100000_x0.5M.fa";


//const char *sam_file = "aln_100000_255.sam";
//const char *sam_file = "s1_1_200.sam";
const char *sam_file = "wrong.sam";
//const char *read_file = "./fault.fa";
//const char *sam_file = "grch38_mason2_100k";
#define MAX_READS 10000
#define MIN_Q 255
#define LEFT_TRIM 0
#define RIGHT_TRIM 0

#define BEST_PERCENT 95
#define ALN_R_DELTA 6
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//#define LEN_READ 200
//#define NUM_EXT_SEED 2
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double outUsedTime();
void set_err_mdl(uint8_t err_mdl[]);
void set_read_data(idx_t *fm_idx, uint32_t read_id, uint8_t *err_mdl, int read_len, uint8_t *f_read_seq, uint8_t *r_read_seq, uint8_t *qual, uint8_t *r_qual);

void set_read_data_1(idx_t *fm_idx, uint32_t read_id, uint8_t *err_mdl, query_t *query);

typedef struct{
    uint8_t name[128];
    uint8_t comment[128];
    uint8_t seq[LEN_READ];
    uint8_t rseq[LEN_READ];
    uint8_t qual[LEN_READ];
    uint8_t apos[LEN_READ];
} read_t; 
typedef struct{
    int mov_L;
    int mov_R;
    int num_L;
    int num_R;
    int num_T;
    uint8_t *apos;
} err_t;

int main()
{	
    double time_of_tot_run = 0.0;
    double time_of_init_sys = 0.0;
    double time_of_read_data = 0.0;  
    double time_of_bwt_aln = 0.0;  
    double time_of_aln_pos = 0.0;
    double time_of_aln_pos_buf = 0.0;
    double time_of_gen_aln = 0.0;
    double time_of_get_seed = 0.0;
    double time_of_smpos = 0.0;
    double time_of_pair = 0.0;
    double time_of_exted0 = 0.0;
    double time_of_exted1 = 0.0;
    double time_of_setBlck = 0.0;
    double time_of_intswed = 0.0;
    double time_of_out_pos = 0.0;
    double time_of_last_sw = 0.0;
    double time_of_LONG = 0.0;
    double cur_time; 
    double start_time, end_time;
    time_of_tot_run = outUsedTime();
    //time_of_init_sys = outUsedTime();
    start_time = time_of_tot_run;

   
    uint32_t START_READ_ID, END_READ_ID; 
    //int NUM_EXT_SEED;
    int MAX_CLS_SEED;
    int SEED_ADD = 4;
    //int SEED_ADD = 0;
    int MAX_CLS_READ;
    int MAX_REAL_NUM = 0;

    //++++++++++++++++++++++++++++++++++++++++++++++++++ 
    int PRINT_INTV = 10000;//打印间隔 
    int PRINT_READS = 0;
    int REAL_FLAG = 1;
    FILE *fp_out_reads;
    if(PRINT_READS == 1) {
        fp_out_reads = fopen("./repeat_reads_100000_x1M.fa", "w");
    }
    //--------------------------------------------------
    if(REAL_FLAG == 1){
        START_READ_ID = 0;
        //NUM_EXT_SEED = 5;
        //NUM_EXT_SEED = NUM_EXT;
        //MAX_CLS_READ = NUM_EXT;
        MAX_REAL_NUM = 1000000;//真实数据数量
        printf("\nIt's REAL data alignment!\n\n");
    }else{
        printf("\nIt's SIM data alignment!\n\n");
        //NUM_EXT_SEED = 5;
        //MAX_CLS_READ = 5;
        //START_READ_ID   = 0; 
        //START_READ_ID   = 801456194; 
        //START_READ_ID   = 804429481; 
        //START_READ_ID = 804429471; 
        START_READ_ID = 1983156844; 
        //START_READ_ID =950000000;
        //END_READ_ID   = 1150000000; 
        END_READ_ID   = 1983156845; 
    }
    //-----------------------------------------
	//主Bwt索引初始化
	//------------------------------------------
	idx_t *fm_idx = (idx_t *)calloc(1, sizeof(idx_t)); 
    end_time = outUsedTime();
    printf("fm_idx mem calloc() end:   %lf\n", end_time-start_time);
    start_time = end_time;

    idx_restore(FILE_PREFIX, fm_idx);
	
    end_time = outUsedTime();
    printf("idx_restore() end:   %lf\n", end_time-start_time);
    start_time = end_time;
    
    printf("Ref_len = %lu\n", fm_idx->bns->l_pac);
    uint32_t *I2Pos= fm_idx->bwt->sa;
	//+++++++++++++++++++++++++++++++++++++++++++
	//以下是多级索引初始化
    struct JmpMod *jmp;
	struct FileName  *fn;	
    struct ExtBlck   eBlck[NUM_EXT+1];
    struct CapIfo    *cap;
	if(NULL == (fn = (struct FileName*)malloc(sizeof(struct FileName)))){
	    perror("error:[malloc(NUM_EXT_SEED*sizeof(FileName))]");
	    exit(1);
    }
    struct SubBuf *sub;
    if(NULL == (sub = (struct SubBuf*)calloc(1, sizeof(struct SubBuf)))){
	    perror("error...");
	    exit(1);
	}
	struct StackTree *sTree ;
	if(NULL == (sTree = (struct StackTree*)malloc(sizeof(struct StackTree)))){
	    perror("error...");
	    exit(1);
	}
    setFileName(fn);
    end_time = outUsedTime();
    printf("setFileName() end:   %lf\n", end_time-start_time);
    start_time = end_time;
 
    InitIdxsArry(fn, &jmp, cap, eBlck, sTree, sub);
    
    end_time = outUsedTime();
    printf("InitIdxsArry() end:   %lf\n", end_time-start_time);
    start_time = end_time;
 
    int i, j;
    for(i = 0; i < 5; ++i) sub->NT_sum[i] = fm_idx->bwt->L2[i]; 
    struct ExtBlck *cB = eBlck;
    //+++++++++++++++++++++++++++++++++++++++++++ 

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
    end_time = outUsedTime();
    printf("hash_boundry() end:   %lf\n", end_time-start_time);
    start_time = end_time;
 
    uint8_t  *f_read_seq, *r_read_seq, *read_seq;
	uint8_t	 seed_seq[32];
    char     *buf_read;	
    int      seed_off;
	uint32_t get_num;
	uint32_t read_id;
    uint32_t buf_algn[20]; 	
	int      Flg_Algn  = 0 ;
	uint32_t cur_read_num = 0;
    FILE *fp_sam = fopen(sam_file, "w");
    if(fp_sam == NULL) {
        printf("open file %s fail!!\n", sam_file);
        exit(1);
    }
    queryio_t *qs = query_open(read_file);;
    read_t *read_buf = (read_t *)calloc(MAX_READS, sizeof(read_t)); 
    qseq_t *query_buf = (qseq_t *)calloc(MAX_READS, sizeof(qseq_t));
    for(i = 0; i < MAX_READS; ++i){
        query_buf[i].name = read_buf[i].name; 
        query_buf[i].comment = read_buf[i].comment; 
        query_buf[i].seq = read_buf[i].seq; 
        query_buf[i].rseq = read_buf[i].rseq; 
        query_buf[i].qual = read_buf[i].qual; 
        query_buf[i].apos = read_buf[i].apos; 
    }
    end_time = outUsedTime();
    printf("query_open() end:   %lf\n", end_time-start_time);
    start_time = end_time;
 
    int read_num = 0;
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int32_t (*call_buf)[4] = (int32_t(*)[4])calloc(360*2*(NUM_EXT+1), 4*sizeof(int32_t)); 
    sub->call_buf = call_buf;
    sub->seqL_off[0] = 0; 
    sub->seqR_off[0] = 360; 
    for(i = 1; i < NUM_EXT+2; i++){
        sub->seqL_off[i] = sub->seqL_off[i-1] + 360*2;
        sub->seqR_off[i] = sub->seqR_off[i-1] + 360*2; 
    }
    sub->sub_err = 5;
    sub->ins_err = 7;
    sub->del_err = 7;
    sub->delta = 6;
    sub->l_seq_delta = 2;
    struct call_t *call = (struct call_t *)calloc(NUM_EXT, sizeof(struct call_t));
    sTree->call = call; 
    sub->aln_out = (aln_out_t *)calloc(1, sizeof(aln_out_t));
    sub->aln_out->max_out = MAX_SEED_NUM;
    sub->aln_out->max_found = 10;
    sub->aln_out->out_buf = (uint32_t (*)[2])calloc(MAX_SEED_NUM, 2*sizeof(uint32_t));

    //以下代码段是种子相关变量初始化 
    seed_t *seed = (seed_t *)calloc(1, sizeof(seed_t)); 

    seed->MAX_SIZE = ((LEN_READ - SEED_LEN)/16 + SEED_ADD)*2;//min_seed = 16
    init_seed_model(seed, fm_idx->bns->l_pac);
    int SEED_SLC_SIZE = seed->MAX_SIZE;
    int MAX_SLC = LEN_READ;
    
    end_time = outUsedTime();
    printf("init_seed_model() end:   %lf\n", end_time-start_time);
    start_time = end_time;
 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    int MAX_POS = MAX_SEED_NUM*10;
    {
     
    uint32_t pos_width = 1; 
    uint32_t len = fm_idx->bns->l_pac/pos_width/64+1; 
        

    sub->NEXT_EXT_NUM = IS_SMLSIZ*1;
    //-------------------------------------------------------
     if(SW_HSH_BIT + SW_RLT_BIT + SW_DAT_BIT != 32) {
        printf("%u, sw_bit_size  error!!!!\n", __LINE__);
        exit(1);
    }

    int to_dat_size = SW_TO_DAT_SIZE;
    int dat_blk_len = (uint32_t)1<<SW_DAT_BIT;
    int rlt_blk_len = (uint32_t)1<<SW_RLT_BIT;
    sub->to_dat_size = to_dat_size;
    sub->dat_blk_len = dat_blk_len;
    sub->rlt_blk_len = rlt_blk_len;
    sub->hsh_shift = 32 - SW_HSH_BIT;
    //sub->rlt_shift = SW_RLT_BIT + 6; // 8 -> 4 ??? 
    sub->rlt_shift = 32 - SW_HSH_BIT; // 8 -> 4 ??? 
    sub->bit_shift = SW_DAT_BIT + 2; //10 ???
    //sub->sw_ed_hsh_size = (((uint32_t)1<<22)+2) * 4 ; // 16M个8bit的字节
    sub->sw_ed_hsh_size = (((uint32_t)1<<SW_HSH_BIT)+2) * sizeof(uint32_t) ; // 16M个8bit的字节
    //sub->sw_to_hsh_size = (((uint32_t)1<<22)+2) * 4; // 16M个8bit的字节 
    sub->sw_to_hsh_size = (((uint32_t)1<<SW_HSH_BIT)+2) * sizeof(uint32_t); // 16M个8bit的字节 
    //sub->sw_ed_rlt_size = 4*4*(MAX_CANDI_POS+2);//1.6M
    sub->sw_ed_rlt_size = sizeof(uint32_t)*rlt_blk_len* (MAX_CANDI_POS+2);//1.6M
    //sub->sw_to_rlt_size = 4*4*(MAX_CANDI_POS+2);//1.6M
    sub->sw_to_rlt_size = sizeof(uint32_t)*rlt_blk_len*(MAX_CANDI_POS+2);//1.6M
    //sub->sw_ed_bit_size = (256/8*(MAX_CANDI_POS + 2)); //3M个8bit的字节,1个数据1bit保存
    sub->sw_ed_bit_size = (dat_blk_len/8*(MAX_CANDI_POS + 2)); //3M个8bit的字节,1个数据1bit保存
    //sub->sw_to_bit_size = (4*256/8*(MAX_CANDI_POS + 2)); //12M个8bit的字节,1个数据4bit保存 
    sub->sw_to_bit_size = (to_dat_size*dat_blk_len/8*(MAX_CANDI_POS + 2)); //12M个8bit的字节,1个数据4bit保存 



    sub->pos_ed_buf_size = sizeof(uint32_t)*(MAX_CANDI_POS+2);
    sub->pos_to_buf_size = sizeof(uint32_t)*(MAX_CANDI_POS+2);
    sub->pos_bk_buf_size = sizeof(uint32_t)*(MAX_CANDI_POS+2);
   
    sub->sw_ed_bit[0] = (uint64_t *)calloc(sub->sw_ed_bit_size/sizeof(uint64_t), sizeof(uint64_t));
    sub->sw_ed_bit[1] = (uint64_t *)calloc(sub->sw_ed_bit_size/sizeof(uint64_t), sizeof(uint64_t));
    sub->sw_to_bit[0] = (uint64_t *)calloc(sub->sw_to_bit_size/sizeof(uint64_t), sizeof(uint64_t));
    sub->sw_to_bit[1] = (uint64_t *)calloc(sub->sw_to_bit_size/sizeof(uint64_t), sizeof(uint64_t));

    sub->sw_ed_hsh[0] = (uint32_t *)calloc(sub->sw_ed_hsh_size/sizeof(uint32_t), sizeof(uint32_t));
    sub->sw_ed_hsh[1] = (uint32_t *)calloc(sub->sw_ed_hsh_size/sizeof(uint32_t), sizeof(uint32_t));
    sub->sw_to_hsh[0] = (uint32_t *)calloc(sub->sw_to_hsh_size/sizeof(uint32_t), sizeof(uint32_t));
    sub->sw_to_hsh[1] = (uint32_t *)calloc(sub->sw_to_hsh_size/sizeof(uint32_t), sizeof(uint32_t));
     
    sub->sw_ed_rlt[0] = (uint32_t *)calloc(sub->sw_ed_rlt_size/sizeof(uint32_t), sizeof(uint32_t));
    sub->sw_ed_rlt[1] = (uint32_t *)calloc(sub->sw_ed_rlt_size/sizeof(uint32_t), sizeof(uint32_t));
    sub->sw_to_rlt[0] = (uint32_t *)calloc(sub->sw_to_rlt_size/sizeof(uint32_t), sizeof(uint32_t));
    sub->sw_to_rlt[1] = (uint32_t *)calloc(sub->sw_to_rlt_size/sizeof(uint32_t), sizeof(uint32_t));
     

   
    sub->pos_ed_buf[0] = malloc(sub->pos_ed_buf_size);
    sub->pos_ed_buf[1] = malloc(sub->pos_ed_buf_size);
    sub->pos_to_buf[0] = malloc(sub->pos_to_buf_size);
    sub->pos_to_buf[1] = malloc(sub->pos_to_buf_size);
    sub->pos_bk_buf[0] = malloc(sub->pos_bk_buf_size*2);
    sub->pos_bk_buf[1] = malloc(sub->pos_bk_buf_size*2);

    sub->pos_ed_buf[0][0] = 0;
    sub->pos_ed_buf[1][0] = 0;
    sub->pos_to_buf[0][0] = 0;
    sub->pos_to_buf[1][0] = 0;
    sub->pos_bk_buf[0][0][0] = 0;
    sub->pos_bk_buf[1][0][0] = 0;
    
     
 

    }
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    uint32_t *pos_buf = sub->pos_buf;
    int seed_slc_size, seed_slc_size_i;
    int seed_id;
  
    uint8_t err_mdl[LEN_READ+20] = {}; 

    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //以下代码段是read相关变量定义及其常用变量    
    uint32_t read_len = LEN_READ;
    uint8_t  ch[LEN_READ], r_ch[LEN_READ];
    uint8_t  qual[LEN_READ] = {}, r_qual[LEN_READ] = {}; 
    f_read_seq = ch;
    r_read_seq = r_ch;
    query_t *query = calloc(1, sizeof(query_t));       
    query->target_idx_r = calloc(LEN_READ, 2*sizeof(uint32_t));
    query->target_idx_f = calloc(LEN_READ, 2*sizeof(uint32_t));
    query->target = calloc(LEN_READ+100, sizeof(uint8_t));
    query->hash_idx[0] = calloc(0xFFFFF+1, sizeof(uint16_t));
    query->hash_idx[1] = calloc(0xFFFFF+1, sizeof(uint16_t));
    
    int l_seed_st = 0, r_seed_ed=LEN_READ;
    uint32_t bgn,end,num,bgn_row,end_row,cap_row, old_seed_id;
    int tot_reads = 0, n_seqs;
    int aln_result_sum[100] = {};   
    int seed_slc_num = 0; 
    int max_stre_cls = 0;
    int MAX_ERR_RATE = 5;
    int MAX_ERR_NUM = LEN_READ * MAX_ERR_RATE/100; 
    int DEF_SCORE = 5;   
    int BEST_SCORE = 0;
    //int ERROR_SCORE = 0;
    //int best_mdl[LEN_READ/16*2];
    int best_mdl[32];
    //int error_aln_score[100];
    float SUB_RATE = 0.02, GAP_RATE=0.01;

    //time_of_init_sys = outUsedTime() - time_of_init_sys;
 
    end_time = outUsedTime();
    printf("init_sys() end:   %lf\n", end_time-start_time);
    start_time = end_time;

    char *str_buf = malloc(MAX_READS * 4096); 
    int n_exact_aln = 0;
//++++++++++++++++++++++++++++ 
int sw_ed_rlt_max[2] = {}; 
int sw_to_rlt_max[2] = {}; 
int sw_ed_bit_max[2] = {};
int sw_to_bit_max[2] = {};
//----------------------------

 
    while(1){ //主循环++++++++++++++++++++++++++++++++
    //cur_time = outUsedTime();        
        if(REAL_FLAG == 1){ 
            if(tot_reads >= MAX_REAL_NUM) break;
            for(i = 0; i < MAX_READS; ++i){
                memset(read_buf+i, 0, sizeof(read_t));  
            }           
            n_seqs = query_read_multiSeqs(qs, MAX_READS, query_buf);
            if(n_seqs == 0) break;
        }
        // 循环过程中每次读入规定个数的reads
        //++++++++++++++++++++++++++++++++
        read_id = 0;
        seed_id = 0;
        seed_slc_size = 1;
        if(REAL_FLAG == 0){    
            tot_reads = 0; 
            //read_id = START_READ_ID;
            n_seqs  = END_READ_ID; 
            query->l_seq = LEN_READ;
            query->seq = ch; 
            query->rseq = r_ch; 
            query->qual = qual; 
            //set_read_data_1(fm_idx, read_id, err_mdl, query);
        } 
        int read_flag = 0; 
        //time_of_read_data +=  outUsedTime() - cur_time;
        char *str_sam = str_buf;
        int *SLC_SWITCH = sub->SLC_SWITCH;
        int *ALN_SWITCH = sub->ALN_SWITCH;
        { 
            int i;
            int sn, seed_num = LEN_READ/SEED_LEN; 
            for(i = 0; i < 10; ++i) {
                SLC_SWITCH[i] = 1;
                ALN_SWITCH[i] = 1;
            }
            if(seed_num < 7) {
               
            } else if(seed_num < 9) { //150
                ALN_SWITCH[2] = 0;  
                //ALN_SWITCH[3] = 0;  
                ALN_SWITCH[4] = 0;  
                //SLC_SWITCH[2] = 0;  
                SLC_SWITCH[3] = 0;  
            } else if(seed_num < 11) {
           
            } else if(seed_num < 13) {
          
            } else {
         
            }
      
        } 
                


        for(read_id = START_READ_ID; read_id < n_seqs; ++read_id, ++tot_reads){ 

//if(tot_reads >= 10290) break;
            /*   
if(tot_reads < 450000) continue;
if(tot_reads == 494892) exit(1);
*/

//if(tot_reads < 900000) continue;
            if(REAL_FLAG&& tot_reads >= MAX_REAL_NUM) break;

            
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             
            if(REAL_FLAG == 1){  
                query->qseq = query_buf+read_id;
                query->name = query->qseq->name;
fprintf(stderr, "\n\nread_id = %u, name = %s\n", read_id, query->name);
                if(tot_reads % PRINT_INTV  == 0) {
                    printf("id = %u, SEQ_NAME = %s\n", tot_reads, query->name);
                }
                query->comment = query->qseq->comment;
                query->l_trim3 = query->qseq->l_trim3; 
                query->l_trim5 = query->qseq->l_trim5; 
                query->seq = query->qseq->seq;
                query->rseq = query->qseq->rseq;
                query->l_seq = query->qseq->l_seq;
                query->qual = query->qseq->qual;
                query->n_ambiguous = query->qseq->n_ambiguous;
                query->apos = query->qseq->apos;
                query->pos = (uint32_t)-1;
                query->cigar = NULL;
                query->seq_start = 0;
                query->seq_end = query->l_seq-1; 
                query->is_gap = -1;
                query->pos = 0xFFFFFFFF;
                query->n_diff = -1;
                query->strand = 3;
                query->b0 = -1;
                query->b1 = -1;
                query->rep_num = 0;
                //query->candi_thres = query->l_seq * WOREST_PERCENT/100;
                
                query->fr_flg[0] = 0;
                query->fr_flg[1] = 0;

                query->candi_thres = 29;
                { 
                    int sn;
                    int sn0 = query->l_seq/SEED_LEN;
                    
                    if(sn0 < 7) {
                        sn = query->l_seq/(SEED_LEN-4); 
                    } else if(sn0 < 9) {
                        sn = query->l_seq/(SEED_LEN-3); 
                    } else if(sn0 < 11) {
                        sn = query->l_seq/(SEED_LEN-2); 
                    } else if(sn0 < 13) {
                        sn = query->l_seq/(SEED_LEN-1); 
                    } else {
                        sn = query->l_seq/SEED_LEN;
                    }
                    if(sn%2 == 0) { sn -= 1;}
                    if(SEED_ADD < 0) SEED_ADD = 0;  
                    if(SEED_ADD*2 > sn) {
                        query->seed_add = (sn + 1)/2;
                        query->seed_num = sn + query->seed_add;
                    } else {
                        query->seed_add = SEED_ADD;
                        query->seed_num = sn + query->seed_add;
                    }
                } 
                int l_seq = query->l_seq;
                sub->qry_f[0] = 0;
                sub->qry_f[1] = 0;
                sub->qry_r[0] = 0;
                sub->qry_r[1] = 0;
                
                for(i = 0; i < sub->aln_out->max_found; ++i) {
                    sub->aln_out->found[i][0] = 0; 
                    sub->aln_out->found[i][1] = 0; 
                    sub->aln_out->found[i][2] = 0; 
                    sub->aln_out->found[i][3] = 0; 
                }  
                sub->aln_out->num = 0; 
                sub->aln_out->len = 0; 
                
                sub->n_seed_bk[0] = 0;
                sub->n_seed_bk[1] = 0;
                for(i = 0; i < query->seed_num*2; ++i) {
                    sub->seed_bk[i] = 0; 
                }
                sub->trm_r = sub->sub_err;
                sub->trm_l = sub->sub_err;
                sub->idx_for[LEN_READ][0] = 0;
                sub->idx_rev[LEN_READ][0] = 0;
                clean_sw_pos_arry(sub);
                //++++++++++++++++++++++++++++++++++++++++++
                if( query->n_ambiguous == 0 ){
cur_time = outUsedTime();
fprintf(stderr, "%u, pos = %u, b0 = %d\n", __LINE__, query->pos, query->b0);
                    int aln_flg = aln_long_exact(fm_idx, hash_boundry, eBlck, query, seed, sub);


                    if(aln_flg == 1) {
printf("%u, pos = %u, b0 = %d\n", __LINE__, query->pos, query->b0);
                        exit(1); 
                    }
                    if(aln_flg > 1) {

fprintf(stderr, "%u, pos = %u, b0 = %d\n", __LINE__, query->pos, query->b0);
                        str_sam += query_aln2sam2(query, fm_idx,sub, str_sam);
                        if(str_sam - str_buf > MAX_READS * 4096) {
                            printf("%u, too much data!!!!\n", __LINE__);
                            printf("%u, read_id = %s, len(str_buf) = %u\n", 
                                    __LINE__, query->name, strlen(str_buf));
                            exit(1);
                        
                        }
                        if(read_id == n_seqs-1 ||tot_reads == MAX_REAL_NUM-1) {
                            fputs(str_buf, fp_sam);
                            str_sam = str_buf; 
                        }
fprintf(stderr, "%u\n", __LINE__);
time_of_bwt_aln += outUsedTime() - cur_time;
                        n_exact_aln++;
                        continue; // while(2)
 
                    }
time_of_bwt_aln += outUsedTime() - cur_time;
                }
                
                fprintf(stderr, "%u, read_id = %s\n", __LINE__, query->name );
                //------------------------------------------ 
                num = query->hash_bak[0][0];
                for(i = 1; i <= num; ++i ){
                    uint32_t seq10 = query->hash_bak[0][i];
                    query->hash_idx[0][seq10] = 0;
                } 
                query->hash_bak[0][0] = 0; 
                
                num = query->hash_bak[1][0];
                for(i = 1; i <= num; ++i ){
                    uint32_t seq10 = query->hash_bak[1][i];
                    query->hash_idx[1][seq10] = 0;
                }
                query->hash_bak[1][0] = 0;
                
                query->is_rev = 0;
                query->read_seq = query->seq; 
                init_idx_target(fm_idx, hash_boundry, 
                                query,  query->target_idx_f);

                query->is_rev = 1;
                query->read_seq = query->rseq; 
                init_idx_target(fm_idx, hash_boundry, 
                                query,  query->target_idx_r);
                
                int trim_len = l_seq - LEFT_TRIM - RIGHT_TRIM;
                query->m_off = trim_len/2;
                query->max_ext = (query->m_off-SEED_LEN/2)/16;
                if(query->max_ext > NUM_EXT) query->max_ext = NUM_EXT;
                MAX_CLS_READ = query->max_ext;
            } else{//if(REAL_FLAG == 0)
                query->qseq = query_buf+0;
                char name[128] = {};
                sprintf(name, "%u", read_id);
                query->qseq->name = name;
                query->l_seq = LEN_READ;
                query->seq = ch; 
                query->rseq = r_ch; 
                query->qual = qual;
                query->n_ambiguous = 0; 
                query->cigar = NULL;
                query->is_rev = 0;
                int trim_len = l - LEFT_TRIM - RIGHT_TRIM;
                query->m_off = trim_len/2;
                query->max_ext = (query->m_off-SEED_LEN/2)/16;
                if(query->max_ext > NUM_EXT) query->max_ext = NUM_EXT;
                MAX_CLS_READ = query->max_ext; 
            }
            
            seed_slc_num = 0;
            max_stre_cls = 0;
            sTree->cls = 0;
            sTree->max_cls = NUM_EXT;
            sTree->max_ext = query->max_ext;

            sTree->len_arry = 0;
            sTree->len_exted = 0;
            sTree->len_buf = 0;
            int left_len, right_len, buf[LEN_READ];
            if(query->n_ambiguous == 0){
                read_flag = 0; 
            } else {
                left_len  = (query->l_seq/2 -l_seed_st);
                right_len = (query->l_seq/2 -r_seed_ed)+SEED_LEN; 
                int j = 0;
                for(i = 0; i < query->n_ambiguous; ++i){
                    if(query->apos[i] >= left_len && query->apos[i] <= right_len){
                        buf[j++] = query->apos[i]; 
                    }     
                }
                if(j < 2) {read_flag = j;}
                else {
                    for(i = 0; i <j-1; ++i){
                        if(buf[i+1]-buf[i] <= SEED_LEN) read_flag = 10; 
                    }
                    if(read_flag < 10) read_flag = 1;
                } 
            }
            //++++++++++++++++++++++++++++++++++++++++++++
            uint32_t pos;
/*  
{
    int __i;
    for(__i = 0; __i < 2; ++__i) {
        if(sub->sw_ed_rlt_max[__i] > sw_ed_rlt_max[__i]) {
            sw_ed_rlt_max[__i] = sub->sw_ed_rlt_max[__i];
        } 
        if(sub->sw_to_rlt_max[__i] > sw_to_rlt_max[__i]) {
            sw_to_rlt_max[__i] = sub->sw_to_rlt_max[__i];
        }
        if(sub->sw_ed_bit_max[__i] > sw_ed_bit_max[__i]) {
            sw_ed_bit_max[__i] = sub->sw_ed_bit_max[__i];
        }
        if(sub->sw_to_bit_max[__i] > sw_to_bit_max[__i]) {
            sw_to_bit_max[__i] = sub->sw_to_bit_max[__i];
        }
    
    }
}
*/
//----------------------------
            //clean_sw_pos_arry(sub);
        
            //--------------------------------------------  
            //以下代码处理非种子扩展比对过程
            seed_id = 0;
            seed_slc_size = 2;
            int SUCEED_FLAG = 0;
            //int seed_num[MAX_SEED];
            int seed_flg[MAX_SEED];
            int seed_bgn[MAX_SEED];
            int seed_end[MAX_SEED];
            int SEED_INIT_ED = 0;
            int LONG_SEED_ED = 0; 
            int THRES_DELTA = (query->l_seq - query->candi_thres)/2;  
            {
                best_mdl[0] = 0;
                best_mdl[1] = -5;
                best_mdl[2] = -5;
                best_mdl[3] = -5;
                best_mdl[4] = -5;
                //for(i = 3; i < seed->MAX_SIZE/2; ++i) {
                for(i = 5; i < query->seed_num; ++i) {
                    best_mdl[i] = -11;
                }
                best_mdl[i] = 0;
                best_mdl[i+1] = -5;
                best_mdl[i+2] = -5;
                best_mdl[i+3] = -5;
                best_mdl[i+4] = -5;
                for(i = query->seed_num+5; i < query->seed_num*2; ++i) {
                    best_mdl[i] = -11;
                }
            } 
            int overlap_num = 0;           
if((seed_id == query->seed_num || seed_id >= query->seed_num*2) ){
    printf("%u, seed_id = %u, num = %d\n", __LINE__, seed_id, query->seed_num);
    exit(1);
}
            while(3) {//种子扩展比对过程++++++++++++++++++++++++++++++++++ 
                fprintf(stderr, "%u, seed_id = %u, num = %d\n", 
                                __LINE__, seed_id, query->seed_num);
                fprintf(stderr, "query->rep_num = %u, b0 = %d\n",
                                query->rep_num, query->b0);
                if((seed_id == query->seed_num || seed_id >= query->seed_num*2) && query->rep_num == 1){
cur_time = outUsedTime();
                    //pos_bk_buf数据处理                   
                    sub->THRES_BAK = 100;
                    aln_sw_bak_buf(fm_idx, hash_boundry, eBlck, query, seed, sub); 

                    //int THRES_SC = 50, THRES_MIN = 70;
                    int THRES_MIN = 70;

                    if(query->b0 < THRES_MIN) {
                        aln_seed_exact_R(fm_idx, hash_boundry, eBlck, query, seed, sub);
                    }                    

                    THRES_MIN = 60;
                    //左右两端近似扩展                   
                    if(query->b0 < THRES_MIN) {
                        aln_apro_R(fm_idx, hash_boundry, jmp, eBlck, 
                                    cB, query, seed, sTree, sub);
                    }
                   
                    if(query->b0 < THRES_MIN) {
                        aln_apro_L(fm_idx, hash_boundry, jmp, eBlck,
                                    cB, query, seed, sTree, sub);
                    }    

         //找质量最好的区间作为比对区间
        if(0){
            aln_qual_best_one(fm_idx, hash_boundry, jmp, eBlck, cB, query, 
                                seed, sTree, sub);
        }// end if(query->b0 < query->l_seq*WOREST_PERCENT*5/3/100)+++
        //所有符合条件质量区域进行比对 
THRES_MIN = 60; 
        if(query->b0 < THRES_MIN){
            int suc = aln_qual_good_all(fm_idx, hash_boundry, jmp, eBlck, cB, query, 
                                seed, sTree, sub);
fprintf(stderr, "%u, suc = %d\n", __LINE__, suc);
        }
        aln_result_sum[19]++;
        aln_result_sum[50+seed_id]++;
time_of_last_sw += outUsedTime() - cur_time;                  
cur_time = outUsedTime(); 
fprintf(stderr, "%u\n", __LINE__);
                    if(REAL_FLAG) {
fprintf(stderr, "%u\n", __LINE__);
                        str_sam += query_aln2sam2(query, fm_idx,sub, str_sam);
                        if(str_sam - str_buf > MAX_READS * 4096) {
                            printf("%u, too much data!!!!\n", __LINE__);
                            printf("%u, read_id = %s, len(str_buf) = %u\n", 
                                    __LINE__, query->name, strlen(str_buf));
                            //exit(1);
                            //fputs(str_buf, fp_sam);
                        } else if(read_id == n_seqs-1 ||tot_reads == MAX_REAL_NUM-1) {
                            fputs(str_buf, fp_sam);
                            str_sam = str_buf; 
                        }
fprintf(stderr, "%u\n", __LINE__);
                    }
                    time_of_out_pos += outUsedTime() - cur_time;
                    break; //while(3)
                }  //if(seed_id == query->seed_num || seed_id >= query->seed_num*2)++++++
                if((seed_id == query->seed_num || seed_id >= query->seed_num*2)
                     && query->rep_num < 1) {
                    seed_id = seed_id % (query->seed_num*2); 
                    query->rep_num++;
                }
                int best_score = best_mdl[seed_id] + query->l_seq;
                fprintf(stderr, "%u, rep_num = %d\n", __LINE__, query->rep_num);
                fprintf(stderr, "best_sc = %d, best_mdl = %d, seed_id = %d\n", 
                                best_score, best_mdl[seed_id], seed_id);
                fprintf(stderr, "b0 = %d, SUCEED_FLAG = %d\n", query->b0, SUCEED_FLAG);
                if(SUCEED_FLAG > 0 || query->b0 > best_score){
                    if(query->rep_num == 0) {
                        if(query->is_rev == 0){
                            seed_id = query->seed_num; 
                        } else{
                            seed_id = 0;
                        }
                        query->rep_num++; 
                    } else{
                        cur_time = outUsedTime(); 
                        if(REAL_FLAG){ 
                            fprintf(stderr, "%u, best_sc = %d, b0 = %d\n", 
                                            __LINE__, best_score, query->b0);
                            str_sam += query_aln2sam2(query, fm_idx,sub, str_sam);
                            if(str_sam - str_buf > MAX_READS * 4096) {
                                printf("%u, too much data!!!!\n", __LINE__);
                                printf("read_id = %s, len(str_buf) = %u\n", 
                                        query->name, strlen(str_buf));
                                exit(1);
                            
                            }

                            if(read_id == n_seqs-1 ||tot_reads == MAX_REAL_NUM -1) {
                                fputs(str_buf, fp_sam);
                                str_sam = str_buf; 
                            }
                        }
                        time_of_out_pos += outUsedTime() - cur_time;
                        int score_row = (query->l_seq - query->b0)/5;
                        if(query->b0 < 0) aln_result_sum[20]++;
                        else if(score_row > 16) aln_result_sum[19]++; 
                        else  aln_result_sum[score_row]++;
                        aln_result_sum[30+seed_id]++;
                        break;//while(3) 
                    }
                    SUCEED_FLAG = 0;                
                }                
                //初始化种子变量
                fprintf(stderr, "%u, seed_id = %u\n", __LINE__, seed_id);
                seed->id = seed_id;

                if(SEED_INIT_ED == 0) {
cur_time = outUsedTime(); 
                    int slc_buf_size = seed_slc_size*sizeof(seed_slc_t);
                    memset(seed->slc, 0, slc_buf_size);
                    if(SEED_FLT) {
                        query->m_off = query->l_seq/3; 
                    } else {
                        query->m_off = query->l_seq/2; 
                    }
                    init_seed_pos(query, seed, sTree);
                    seed_slc_size = query->seed_num*2; 
                    for(i = 0; i < seed_slc_size; ++i) { 
                        seed->del[i] = 0;
                        //seed_num[i] = 0;
                    }
                    seed->n_for = 0;
                    seed->n_rev = 0;
                    seed->max_f = 0;
                    seed->max_r = 0;
                    seed->n_sm_f = 0;
                    seed->n_sm_r = 0;
                    int max_num_for = 0;
                    int max_num_rev = 0;
                    for(i = 0; i < seed_slc_size; i++) {
                        seed->id = i; 
                        num = gen_seed_kmer(fm_idx, hash_boundry, query, seed);
                        fprintf(stderr, "%u, seed_id=%d, num=%d, len = %d\n", __LINE__, i, num, seed->slc[i].len); 
                        if(num > 0) {
                            if(i < query->seed_num){
                                seed->n_for++;
                                if(num > seed->max_f) seed->max_f = num;
                                if(num <= SW_THRES) {
                                    seed->n_sm_f++;
                                } else{
                                    max_num_for++; 
                                    seed->slc[i].sum = max_num_for;
                                }
                            } else{
                                seed->n_rev++; 
                                seed->slc[i].sum = seed->n_rev;
                                if(num > seed->max_r) seed->max_r = num;
                                if(num <= SW_THRES) seed->n_sm_r++;
                            }
                        }
                    }// end for(i = 0; i < seed_slc_size; i++)++++++         
                    if(seed->n_for >= seed->n_rev) {
                        seed_id = 0;
                    } else {
                        seed_id = query->seed_num;
                    }
                    int n_bg_r = seed->n_rev - seed->n_sm_r;
                    int n_bg_f = seed->n_for - seed->n_sm_f;
                    if(n_bg_f <= n_bg_r){
                        query->is_rev = 0; 
                    } else{
                        query->is_rev = 1;
                    }
                    sub->n_sm_f = seed->n_sm_f;
                    sub->n_sm_r = seed->n_sm_r;
                    query->n_olp_for = seed->n_for -1;
                    query->n_olp_rev = seed->n_rev -1;
                    seed->cls = 0;
                    SEED_INIT_ED = 1;
                    LONG_SEED_ED = 0; 
time_of_get_seed += outUsedTime() - cur_time;
                }
                
                                //+++++++++++++++++++++++++++++++++
                if(LONG_SEED_ED == 0){
cur_time = outUsedTime();
                    LONG_SEED_ED = 1; 
                    
//printf("%u, read_id  = %s, b0 = %d\n", __LINE__, query->name, query->b0); 
fprintf(stderr, "%u, read_id  = %s, b0 = %d\n", __LINE__, query->name, query->b0); 
                    sub->thres_sw_to = 2;
                    sub->thres_pos_num = IS_SMLSIZ; // num = 1
                    sub->thres_sw_olp = 1;
                    sub->olp_flg = 1;
                    sub->err_sum[0] = 0;
                

                    
                    
                    if(query->b0 >= query->l_seq - sub->sub_err) {
                        printf("%u, read_id  = %s, b0 = %d\n", __LINE__, query->name, query->b0); 
                        exit(1);
                       
                    }  
                    fprintf(stderr, "%u, read_id  = %s, b0 = %d\n", __LINE__, query->name, query->b0);                    
                    int aln_len = aln_long_seed(fm_idx, hash_boundry, eBlck, query, seed, sub);
if(0) {
                    //if(query->b0 > query->l_seq - sub->delta) {
                       
                        if(query->b0 > query->l_seq - sub->sub_err) {
                            printf("%u, read_id  = %s, b0 = %d\n", __LINE__, query->name, query->b0); 
                            exit(1);
                        }
                        
                        
                        int rev, sid, sn = query->seed_num;
                        int n_sm[2], n_mx[2];
                        int s_off0, s_off1, d;
                        int min0, min1, m_i0, m_i1;    
                        n_sm[0] = seed->n_sm_f;
                        n_sm[1] = seed->n_sm_r;
                        n_mx[0] = seed->max_f;
                        n_mx[1] = seed->max_r;
                        for(rev = 0; rev < 2; ++rev){
                            if(query->fr_flg[rev] == 0) {
                                for(sid = 0; sid < sn; ++sid){
                                    seed->slc[sid+rev*sn].num = 0; 
                                } 
                            } else {
                                if(n_sm[rev] == 0 && n_mx[rev] > SW_THRES) {
                                    m_i0 = rev*sn;
                                    m_i1 = rev*sn;  
                                    s_off0 = seed->slc[m_i0].s_off;
                                    min1 = 0; 
                                    for(sid = rev*sn; sid < sn+rev*sn; ++sid){
                                        if(sid == m_i0) continue;
                                        if(seed->slc[sid].ext_num > min1){
                                            s_off1 = seed->slc[sid].s_off;
                                            d = s_off0 - s_off1;
                                            if(d < 0) d = -d;
                                            if(d < SEED_LEN) continue;
                                            min1 = seed->slc[sid].num;
                                            m_i1 = sid;
                                        } 
                                    }
                                } else {
                                    min0 = 0xFFFFFF, min1 = 0xFFFFFF;    
                                    for(sid = rev*sn; sid < sn+rev*sn; ++sid){
                                        if(seed->slc[sid].num >0 && 
                                            seed->slc[sid].num < min0){
                                            min0 = seed->slc[sid].num;
                                            m_i0 = sid; 
                                        }
                                    } 
                                    s_off0 = seed->slc[m_i0].s_off;
                                    for(sid = rev*sn; sid < sn+rev*sn; ++sid){
                                        if(sid == m_i0) continue;
                                        if(seed->slc[sid].num < min1){
                                            s_off1 = seed->slc[sid].s_off;
                                            d = s_off0 - s_off1;
                                            if(d < 0) d = -d;
                                            if(d < SEED_LEN) continue;
                                            min1 = seed->slc[sid].num;
                                            m_i1 = sid; 
                                        }
                                    }
                                    if(seed->slc[m_i1].num > SW_THRES) {
                                        s_off0 = seed->slc[m_i0].s_off; 
                                        m_i1 = rev*sn;  
                                        min1 = 0; 
                                        for(sid = rev*sn; sid < sn+rev*sn; ++sid){
                                            if(sid == m_i0) continue;
                                            if(seed->slc[sid].ext_num > min1){
                                                s_off1 = seed->slc[sid].s_off;
                                                d = s_off0 - s_off1;
                                                if(d < 0) d = -d;
                                                if(d < SEED_LEN) continue;
                                                min1 = seed->slc[sid].num;
                                                m_i1 = sid;
                                            } 
                                        }
                                    }
                                } 
                                for(sid = rev*sn; sid < rev*sn + sn; ++sid) {
                                    if(sid == m_i0 || sid == m_i1) continue;
                                    seed->slc[sid].num = 0; 
                                }
                            } 
                        }     

/*   
                        int (*found)[4] = sub->aln_out->found; 
                        int tot = 0;
                        int fi;
                        for(fi = 0; fi < sub->delta; ++fi){
                            tot += found[fi][1]; 
                        }

                        if(tot > THRES_ALN_OUT) {
                       
                            int qlen = query->l_seq;
                            str_sam += query_aln2sam2(query, fm_idx,sub, str_sam);
                            if(str_sam - str_buf > MAX_READS * 4096) {
                                printf("%u, too much data!!!!\n", __LINE__);
                                printf("%u, read_id = %s, len(str_buf) = %u\n", 
                                        __LINE__, query->name, strlen(str_buf));
                                exit(1);
                            
                            }
                            if(read_id == n_seqs-1 ||tot_reads == MAX_REAL_NUM-1) {
                                fputs(str_buf, fp_sam);
                                str_sam = str_buf; 
                            }
time_of_LONG +=  outUsedTime() - cur_time;
                            break;// while(3)
                        }
*/
                    }
time_of_LONG +=  outUsedTime() - cur_time;
                }
                
                //---------------------------------
                sub->thres_sw_to = 2;
                sub->thres_pos_num = IS_SMLSIZ; // num = 1
                sub->thres_sw_olp = 1;
                sub->olp_flg = 1;
                sub->err_sum[0] = 0;
                
fprintf(stderr, "%u, seed->max_r = %d, max_f = %d\n", __LINE__, seed->max_r, seed->max_f);
                //所有的种子都进行SW比对
                if(seed->max_r <= SW_THRES && seed->max_f <= SW_THRES ){
cur_time = outUsedTime();
                    SUCEED_FLAG = aln_sw_seed_all(fm_idx, hash_boundry, eBlck, query, seed, sub);
                    seed_id = query->seed_num*2;
time_of_smpos += outUsedTime() - cur_time;
                    continue;
                }
                
                //部分种子进行SW比对
                int sm_f = seed->n_sm_f, sm_r = seed->n_sm_r;
                int seed_num = query->seed_num;
fprintf(stderr, "%u, b0 = %d\n", __LINE__, query->b0);
                if((seed_id < seed_num && sm_f > 0) || (seed_id >= seed_num && sm_r > 0 )){
cur_time = outUsedTime();
                    seed->id = seed_id;               
                    SUCEED_FLAG = aln_sw_seed_some(fm_idx, hash_boundry, eBlck, 
                                                    query, seed, sub);
time_of_aln_pos += outUsedTime() - cur_time; 
                }
                //产生种子新的序列 
fprintf(stderr, "%u, b0 = %d, query->pos = %u\n", __LINE__, query->b0, query->pos);
                while(seed_id < seed_slc_size){
                    if(seed->del[seed_id] > 0) {
                        seed_id++;
                        continue; 
                    }
                    num = seed->slc[seed_id].num;
                    if(num == 0) {
                        seed_id++;
                        continue; 
                    } else{// num > 0
                        MAX_CLS_SEED = seed->slc[seed_id].ext_num;
                        if(MAX_CLS_SEED == 0) {
                            seed_id++;
                            continue;
                        }
                        sTree->max_cls = MAX_CLS_SEED;
                        seed->bgn = seed->slc[seed_id].bgn;
                        seed->num = seed->slc[seed_id].num;
                        seed->end = seed->slc[seed_id].end;
                        break; 
                    }
                }//end while(seed_id < seed_slc_size)+++++++++++++
                if(SUCEED_FLAG > 0) continue;
                if(seed_id <  query->seed_num && seed->n_sm_f > 0) continue;
                if(seed_id >= query->seed_num && seed->n_sm_r > 0) continue;
                if(seed_id == seed_slc_size) continue;
                
                BEST_SCORE = best_mdl[seed_id] + query->l_seq;
                query->best_thres  = BEST_SCORE; 
                query->error_thres = MAX_ERR_NUM*DEF_SCORE; 
                seed_slc_num++;
                bgn = seed->slc[seed_id].bgn;
                end = seed->slc[seed_id].end;
                num = seed->slc[seed_id].num;
                //seed_num[seed_id] = num;
                if(num ==0) continue;    

cur_time = outUsedTime(); 
                //选定的read序列指针赋值给read_seq
                if(seed_id < query->seed_num){
                    read_seq = query->seq;
                    query->is_rev = 0;
                    sub->thres_sw_to = query->n_olp_for;
                    query->target_idx = query->target_idx_f;
                } else{//slc[seed_id].drct ==1
                    read_seq = query->rseq;
                    query->is_rev = 1;
                    sub->thres_sw_to = query->n_olp_rev;
                    query->target_idx = query->target_idx_r;
                }
                sub->thres_pos_num = 0;
                //sub->thres_pos_num = IS_SMLSIZ/4;
                sub->thres_sw_to = sub->thres_sw_to/3;
                query->read_seq = read_seq;
                read_len = query->l_seq;
                seed_off = seed->slc[seed_id].s_off;
                seed->cls = 0; 
                old_seed_id = seed_id;
                int L_offset = seed_off;
                int R_offset = SEED_LEN + seed_off;
                //+++++++++++++++++++++++++++++++++++++++++
                uint32_t pos_i = 0;
               
                buf_algn[0] = SEED_LEN;

                buf_algn[1] = bgn;
                buf_algn[2] = num;

                if(num > IS_SMLSIZ) {
                    sTree->len_buf++;
                } 
                getCapPos(jmp,buf_algn);


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    uint32_t SEED_NXT_PNT = buf_algn[3];
    uint32_t SEED_NUM = num;

    int cls = 0, max_cls = 0; 
    int err_score = 0;
    uint32_t (*ext_idx)[2];
    int test_rep = 0;
    int sTree_finish = 0;
    struct call_t *cal; 
    uint32_t cap_row;
    memset(call, 0, NUM_EXT*sizeof(struct call_t));
    sub->query_err = query->l_seq-query->b0;
    for(i = 0; i <= NUM_EXT; ++i){
        sub->path_err[i] = sub->query_err;
    } 
   
        sTree->cls = 0;
        sTree->stck_arry[1].nxtpnt = SEED_NXT_PNT;
        sTree->stck_arry[1].nxtflg = SEED_NUM;
        sTree->stck_arry[1].err = 0;
        sTree->stck_arry[1].l_off = 0; //64
        sTree->stck_arry[1].r_off = 0; //64
        sTree->stck_arry[1].cls   = 0;
        sTree->stck_arry[1].len_p = 0;
        sTree->stck_arry[1].len_g = 0;
        sTree->stck_arry[1].aln_f = 0;
        sTree->stck_arry[1].aln_r = 0;
 
        sTree->len_arry = 1; 
        //sTree->len_old = 1; 
        sTree->len_old = 0; 
        sTree->len_old_back = 0; 
        sub->seqL_out = sub->call_buf + sub->seqL_off[0]; 
        sub->seqR_out = sub->call_buf + sub->seqR_off[0]; 
        sub->seqL_out[0][0] = 0;
        sub->seqR_out[0][0] = 0;
       
        sTree->len_exted = 0; 
        sTree->len_smpos = 0; 
        sTree->seed_off = seed_off;
        sTree->seed_id = seed_id;
        test_rep = 0;
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
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
        //以下循环是多级索引中比对			
        int cal_flag = 0;
        int EXT_END = 0;
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
        //一个种子扩展比对过程开始
        //----------------------------------------------------------------
        if(sTree->max_cls > 3 ) { 
            sTree->max_cls = 3;
        }
        while(5){
            fprintf(stderr, "\n%u test_rep = %u  =========\n", 
                                __LINE__, test_rep++);
            fprintf(stderr, "sTree->max_ext = %d, max_cls = %d, seed_id = %d\n", 
                                sTree->max_ext, sTree->max_cls, sTree->seed_id);
            //++++++++++++++++++++++++++++++++
            int len_old_flg = sTree->stck_arry[sTree->len_arry].len_p - sTree->len_old;
            if(sTree->len_arry == 0) break; 
            if(len_old_flg == 0) { 
                if(call[sTree->cls].aln_f > 4) { 
                    sTree->len_arry--;
                    if(sTree->len_arry > 0) {
                        sTree->cls = sTree->stck_arry[sTree->len_arry].cls;
                        sTree->len_old = sTree->stck_arry[sTree->len_arry].len_p;
                        call[sTree->cls].aln_f = sTree->stck_arry[sTree->len_arry].aln_f; 
                        call[sTree->cls].aln_r = sTree->stck_arry[sTree->len_arry].aln_r; 
                    }
                    continue;
                }             
            } else if(len_old_flg > 0){
                sTree->len_old = sTree->stck_arry[sTree->len_arry].len_p;
                sTree->cls = sTree->stck_arry[sTree->len_arry].cls;
                call[sTree->cls].aln_f = sTree->stck_arry[sTree->len_arry].aln_f; 
                call[sTree->cls].aln_r = sTree->stck_arry[sTree->len_arry].aln_r; 
                if(sTree->len_old >= sTree->len_arry) {
                    printf("%u, %s, error!!!\n", __LINE__, __func__);
                    exit(1);                    
                }             
            } else{
                printf("%u, %s, error!!!\n", __LINE__, __func__);
                exit(1);
            } 
            
            fprintf(stderr, "%u\n", __LINE__);
            fprintf(stderr, "len_arry = %u, old_arry = %u, cls = %u\n", 
                                sTree->len_arry, sTree->len_old, sTree->cls);
            fprintf(stderr, "call.aln_f = %u, aln_r = %u, pair_f = %u\n", 
                                call[sTree->cls].aln_f, call[sTree->cls].aln_r, 
                                call[sTree->cls].pair_f);
//cur_time = outUsedTime(); 
            clean_sub_buf(cB, sub);
            setBlckData(eBlck,sub,  sTree, read_seq,  seed_off, &cB);  
/*  
fprintf(stderr, "%u, num_realt = %d\n", __LINE__, cB->num_relat);
for(i = 0; i < cB->num_relat; ++i) {
    uint32_t *nxtpnt = cB->head_nxtpnt + cB->nxtpnt;
    uint8_t *nxtflg = cB->head_nxtflg + cB->nxtpnt;
    uint32_t cap_row = nxtpnt[i];
    uint8_t cap_flg = nxtflg[i];
    if(cap_flg == 1) cap_row = bwt_get_idx(fm_idx->bwt, cap_row+16);
    fprintf(stderr, "%u, cap_row = %u, cap_flg = %d\n", __LINE__, cap_row, cap_flg);

}     
exit(1);
*/
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            //以下代码段完成比对过程
            //---------------------------------------------------------------- 
            cls = 0;
            max_cls = 0;
            int r_flag; 
            int max_err, sub_err, del_err;
            if(call[sTree->cls].aln_f == 0) {
                for(i =0; i <3; ++i){
                    sub->seqL_out[0][i] = 0;             
                    sub->seqR_out[0][i] = 0;
                }
            }
            while(6) {
                fprintf(stderr, "%u, seed_id = %u, aln_f = %u\n", 
                                    __LINE__, seed_id, call[sTree->cls].aln_f);

                if(call[sTree->cls].aln_f >= 5) { 
                    break; //while(6)
                }
                //以下函数生成比对数据 
                sub->cls = sTree->cls;
                sub->query_err = query->l_seq - query->b0;
                r_flag = gen_aln_data(cB, sub, call+sTree->cls);
fprintf(stderr, "%u, seed_id = %d\n",__LINE__, seed_id); 
                if(r_flag <=0 ) {
                    uint32_t nxtflg = sTree->stck_arry[sTree->len_arry].nxtflg;
                    call[sTree->cls].aln_f = 5;    
                    break;//while(6)
                }
                if(call[sTree->cls].aln_f  >=5) {
                    break;  //while(6)
                }
                set_pair_data(sub, call+sTree->cls); //设定配对数据 
                int flg = get_next_aln(sub, call+sTree->cls );//设定后续比对过程参数 
                sub->pair_out->p_num = 0;
                sub->pair_out->s_num = 0;
                Flg_Algn = 0;
                for(i = 0; i < sub->pair_b2e_num; ++i) { 
                    sub->pair_b2e_id = i;
                    Flg_Algn += PairExtSeq_all(cB,sub); //完成配对过程
                    sub->pair_out->p_buf[i] = sub->pair_out->p_num;
                    if(i == 1) {
                        if(sub->pair_out->p_num > 10){ 
                            sub->pair_b2e_num = 2;
                        }
                    } 
                }
                if(Flg_Algn == 0){ 
                    if(call[sTree->cls].aln_f  >=5) {
                        break;  //while(6)
                    }
                    continue;
                }
                if(sub->pair_out->p_num >0) { 
                    Flg_Algn = sortPairExtSeq(sub);//配对结果分类与排序 
//time_of_pair += outUsedTime() - cur_time;
                    cls =  sTree->stck_arry[sTree->len_arry].cls;
                    int err = sTree->stck_arry[sTree->len_arry].err;
                    if(sub->path_err[cls] > err) sub->path_err[cls] = err;    
                    uint32_t num_smpos = sub->pair_out->s_num;
                    uint32_t num_bgpos = sub->pair_out->p_num - num_smpos;

                    if(num_smpos >0){
                        seed->id = seed_id; 
                        Flg_Algn = AlgnPos_buf(fm_idx, eBlck, query, seed, 
                                                num_smpos, sTree->cls, sub);  
                    }//end if(num_smpos>0)+++++++
                    fprintf(stderr, "%u, cls = %u, aln_f = %u, len_exted = %u\n", 
                     __LINE__, sTree->cls, call[sTree->cls].aln_f, sTree->len_exted); 
                    if(sub->pair_out->p_num- sub->pair_out->s_num >0){
                        setStackTree(sTree, sub, seed, seed_id, eBlck);
                        sub->err_sum[0] = sTree->stck_arry[sTree->len_arry].err;
                        if(sTree->len_exted > 0 ){ 
                            if(sTree->exted_arry[sTree->len_exted-1].cls 
                                == sTree->max_ext ) {
                                    Flg_Algn = gen_pos_exted0(fm_idx, hash_boundry, query, 
                                        sTree, sub, eBlck, sTree->max_ext);                 
                            }
                        }
                    }//end if(sub->pair_out[0][0]- sub->pair_out[0][1] >0)+++
                }//end if(sub->pair_out[0][0] >0)++++++++++++ 
                break; //while(6) 深度优先搜索
            }//end while(6)++++++++++++++++++++++++++++++++++++++++++++
            if(SUCEED_FLAG > 0 || EXT_END > 0) break; //while(5)
            if(sTree->len_arry == 1 && call[sTree->cls].aln_f > 4) { 
                break;
            }
            if(sTree->len_arry == 0) {
                printf("%u, error: len_arry = 0!\n", __LINE__);
                exit(1);
            }
        }// End: while(5) //一个read的扩展比对结束+++++++++
time_of_aln_pos_buf += outUsedTime() - cur_time;
cur_time = outUsedTime(); 
        if(sTree->len_exted > 0){ 
            for( i = 0; i < sTree->len_exted; ++i) {
                int flg_num = sTree->exted_arry[i].nxtflg;
                if(flg_num >= 255 ){
                    aln_result_sum[60+1]++; 
                } else{
                    aln_result_sum[60+0]++; 
                }
                cls = sTree->exted_arry[i].cls;
                aln_result_sum[60+1+cls]++;               
            }
            aln_result_sum[70+seed_id]++;               
            AlgnPos_exted1(fm_idx, hash_boundry, query, seed, sTree, 
                            sub, eBlck, sTree->max_cls);
        }
time_of_exted1 += outUsedTime() - cur_time; 
//++++++++++++++++++++++++++++++
        seed_id++; 
    }// End : while(3) 一次读入read序列集合的比对结束+++++++
    } // End: while(2)+++++++++++++++++++++++++++++++++++++++++++     
        if(REAL_FLAG==0)  break;//while(1)++++++++++
    }//End while(1) 主循环++++++++++++++++++++++++++++++++++++++++
    

/*  
    printf("aln_sum =  ");
    for(i= 0;i < 100; ++i) {
        printf("%u :  %d  \n", i, aln_result_sum[i]);
        printf("\n");
        if(i == 9) printf("------------------------------\n");

    }
    printf("\n");
    printf("----------------------\n");
    fprintf(stderr, "aln_sum =  ");
    for(i= 0;i < 100; ++i) {
        fprintf(stderr, "%u :  %d  \n", i, aln_result_sum[i]);
        fprintf(stderr, "\n");
        if(i == 9) fprintf(stderr, "------------------------------\n");
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "----------------------\n");
*/
printf("%u\n", __LINE__);
    if(PRINT_READS == 1) { fclose(fp_out_reads);}
printf("%u\n", __LINE__);
    free(fn->jmpmod);
printf("%u\n", __LINE__);
    free(fn->seedidx);
printf("%u\n", __LINE__);
    free(fn->capidx);
printf("%u\n", __LINE__);
    free(fn->comfile);
printf("%u\n", __LINE__);
    free(fn);
printf("%u\n", __LINE__);
    free(str_buf);
printf("%u\n", __LINE__);
    free(call_buf);
printf("%u\n", __LINE__);
    free(seed->slc);
printf("%u\n", __LINE__);
    free(seed->err);
printf("%u\n", __LINE__);
    free(seed);
printf("%u\n", __LINE__);
    free(query->target_idx_r);
printf("%u\n", __LINE__);
    free(query->target_idx_f);
printf("%u\n", __LINE__);
    free(query->target);
printf("%u\n", __LINE__);
    free(query->hash_idx[0]);
printf("%u\n", __LINE__);
    free(query->hash_idx[1]);
printf("%u\n", __LINE__);
    free(query);
printf("%u\n", __LINE__);
    query_close(qs);
printf("%u\n", __LINE__);

    idx_destroy(fm_idx);
printf("%u\n", __LINE__);
    free(fm_idx);
printf("%u\n", __LINE__);
    destroyIdxsArry(jmp, sTree, eBlck, sub);
printf("%u\n", __LINE__);
    
    free(sub->sw_ed_bit[0]);
printf("%u\n", __LINE__);
    free(sub->sw_ed_bit[1]);
printf("%u\n", __LINE__);
    free(sub->sw_to_bit[0]);
printf("%u\n", __LINE__);
    free(sub->sw_to_bit[1]);
printf("%u\n", __LINE__);
    free(sub->sw_ed_hsh[0]);
printf("%u\n", __LINE__);
    free(sub->sw_ed_hsh[1]);
printf("%u\n", __LINE__);
    free(sub->sw_to_hsh[0]);
printf("%u\n", __LINE__);
    free(sub->sw_to_hsh[1]);
printf("%u\n", __LINE__);
    free(sub->sw_ed_rlt[0]);
printf("%u\n", __LINE__);
    free(sub->sw_ed_rlt[1]);
printf("%u\n", __LINE__);
    free(sub->sw_to_rlt[0]);
printf("%u\n", __LINE__);
    free(sub->sw_to_rlt[1]);
printf("%u\n", __LINE__);
    free(sub->pos_ed_buf[0]);
printf("%u\n", __LINE__);
    free(sub->pos_ed_buf[1]);
printf("%u\n", __LINE__);
    free(sub->pos_to_buf[0]);
printf("%u\n", __LINE__);
    free(sub->pos_to_buf[1]);
printf("%u\n", __LINE__);
    free(sub->pos_bk_buf[0]);
printf("%u\n", __LINE__);
    free(sub->pos_bk_buf[1]);
printf("%u\n", __LINE__);


    free(sub);
printf("%u\n", __LINE__);
    free(jmp);
printf("%u\n", __LINE__);
    free(sTree);
printf("%u\n", __LINE__);
    free(query_buf);
printf("%u\n", __LINE__);
    free(read_buf);
printf("%u\n", __LINE__);
    fclose(fp_sam);
printf("%u\n", __LINE__);
  
    //kseq_destroy(seq);
    //gzclose(fp);
    
    //++++++++++++++++++++++++++++++++++++++++++++++
    time_of_tot_run = outUsedTime()-time_of_tot_run;
printf("%u\n", __LINE__);
    int time, h, m, s;

    //-------------------------------------- 
    printf("time_tot_run: ");
    fprintf(stderr, "time_of_tot_run:     ");
    time = time_of_tot_run;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_of_sys : ");
    fprintf(stderr, "time_of_sys:    ");
    time = time_of_init_sys;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_of_read: ");
    fprintf(stderr, "time_of_read_data:   ");
    time =  time_of_read_data;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 

//-------------------------------------- 
    printf("time_of_bwt : ");
    fprintf(stderr, "time_of_bwt_aln:   ");
    time = time_of_bwt_aln;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 





    printf("time_get_sed: ");
    fprintf(stderr, "time_of_get_seed:    ");
    time = time_of_get_seed;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_smpos  : ");
    fprintf(stderr, "time_of_smpos:    ");
    time = time_of_smpos;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_sw_pos : ");
    fprintf(stderr, "time_of_aln_pos:     ");
    time = time_of_aln_pos;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_pos_buf: ");
    fprintf(stderr, "time_of_aln_pos_buf: ");
    time = time_of_aln_pos_buf;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_exted0 : ");
    fprintf(stderr, "time_of_exted0:      ");
    time = time_of_exted0;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_exted1 : ");
    fprintf(stderr, "time_of_exted1:      ");
    time = time_of_exted1;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_gen_aln: ");
    fprintf(stderr, "time_of_gen_aln:     ");
    time = time_of_gen_aln;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_of_pair: ");
    fprintf(stderr, "time_of_pair:        ");
    time = time_of_pair;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_setBlck: ");
    fprintf(stderr, "time_of_setBlck:     ");
    time = time_of_setBlck;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_intswed: ");
    fprintf(stderr, "time_of_intswed:     ");
    time = time_of_intswed;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u\n", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_of_last: ");
    fprintf(stderr, "time_of_last:     ");
    time = time_of_last_sw;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u ", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_out_pos: ");
    fprintf(stderr, "time_of_out_pos:     ");
    time = time_of_out_pos;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u ", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
    //-------------------------------------- 
    printf("time_out_LNG: ");
    fprintf(stderr, "time_of_out_LNG:     ");
    time = time_of_LONG;
    h = time/3600;
    m = (time%3600)/60;
    s = (time%3600)%60;
    printf("%6u ", time);
    fprintf(stderr, "%u ", time);
    printf("%3d h %2d m %2d s\n", h, m, s);
    fprintf(stderr, "%d h %d m %d s\n", h, m, s);
/*  
    int __i;
    for(__i = 0; __i < 2; ++__i) {
        printf("sw_ed_rlt_max[%d] = %d\n", __i, sw_ed_rlt_max[__i]);
        printf("sw_to_rlt_max[%d] = %d\n", __i, sw_to_rlt_max[__i]);
        printf("sw_ed_bit_max[%d] = %d\n", __i, sw_ed_bit_max[__i]);
        printf("sw_to_bit_max[%d] = %d\n", __i, sw_to_bit_max[__i]);    
    }
*/
    return 0;
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
void print_read_data_1(idx_t *fm_idx, uint32_t read_id, uint8_t *err_mdl, query_t *query)
{

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
/* ins  
    int L_seg = 73 - 0;
    int R_seg = 112 + 0; 
    int indel_pos = (74 + read_id%14 + 1) - 0 ; 
*/    
//----------------------------------------------------------------------   
    int L_seg = 88 - 0;
    int R_seg = 128 + 0; 
    int indel_pos = (110 + read_id%14 + 1) + 0 ;   
//----------------------------------------------------------------------   
 
    int read_len = query->l_seq;
    uint8_t *f_read_seq = query->seq; 
    uint8_t *r_read_seq = query->rseq; 
    uint8_t *qual = query->qual; 
    //uint8_t *r_qual = query->qual; 
    uint8_t *cur_err_mdl = err_mdl+read_id%20;
    uint32_t pos = read_id;
    int i; 
/*  
    fprintf(stderr, "origin:\n");
    for(i=0; i <read_len; ++i) fprintf(stderr, "%u", f_read_seq[i]);
    fprintf(stderr, "\n");
*/
    fprintf(stderr,"%u, ------------------\n", __LINE__);
    for(i=74; i <90; ++i) fprintf(stderr, "%u ", __get_pac(fm_idx->pac, atol(query->qseq->name)+i));
    fprintf(stderr, "\n");
    for(i=74; i <90; ++i) fprintf(stderr, "%u ", __get_pac(fm_idx->pac, read_id+i));
    fprintf(stderr, "\n");
    for(i=74; i <90; ++i) fprintf(stderr, "%u ", f_read_seq[i]);
    fprintf(stderr, "\n");



    for(i=110; i <=126; ++i){ 
        if(i == indel_pos) fprintf(stderr, "- ");
        fprintf(stderr, "%u ", __get_pac(fm_idx->pac, read_id+i));
    }
    fprintf(stderr, "\n");

    for(i = 0; i < read_len; ++i ){
        uint8_t c  = __get_pac(fm_idx->pac, read_id+i);
        if(cur_err_mdl[i] > 0) {
            if(i<L_seg||i>R_seg){
                fprintf(stderr, "err[%u] = %u->%u\n", i, c, (c+cur_err_mdl[i])%4);
            }
        }
    }
    //f_read_seq[112] = (f_read_seq[112]+1)%4;  
    //for(i = 112; i < read_len-1; ++i ){
/*  
    for(i = indel_pos; i < read_len-1; ++i ){
        //f_read_seq[88] = (f_read_seq[88]+1)%4;
        f_read_seq[i] = f_read_seq[i+1];
    }
*/
    fprintf(stderr, "indel_pos = %u, indel = %d\n", indel_pos, 1);
  
/*  
    for(i=L_seg; i <90; ++i) {
        if(i == indel_pos) fprintf(stderr, "- ");
        fprintf(stderr, "%u ", f_read_seq[i]);
    }
*/
    for(i=110; i <=126; ++i) {
        fprintf(stderr, "%u ", f_read_seq[i]);
    }
    
    fprintf(stderr, "\n");
    fprintf(stderr,"%u, ------------------\n", __LINE__);
/*  
    fprintf(stderr, "var:\n");
    for(i=0; i <read_len; ++i) fprintf(stderr, "%u", r_read_seq[i]);
    fprintf(stderr, "\n");
*/
 

    return;
}

void set_read_data_1(idx_t *fm_idx, uint32_t read_id, uint8_t *err_mdl, query_t *query)
{

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
/* ins  
    int L_seg = 73 - 0;
    int R_seg = 112 + 0; 
    int indel_pos = (74 + read_id%14 + 1) - 0 ; 
*/    
//----------------------------------------------------------------------   
    int L_seg = 88 -  0;
    int R_seg = 128 + 0; 
    int indel_pos = (110 + read_id%14 + 1) + 0 ;   
//----------------------------------------------------------------------   
 
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
/*  
    fprintf(stderr,"%u, ------------------\n", __LINE__);
    //for(i=L_seg; i <=90; ++i) fprintf(stderr, "%u ", f_read_seq[i]);
    for(i=110; i <=126; ++i){ 
        if(i == indel_pos) fprintf(stderr, "- ");
        fprintf(stderr, "%u ", f_read_seq[i]);
    }
    fprintf(stderr, "\n");

    for(i = 0; i < read_len; ++i ){
        uint8_t c  = f_read_seq[i];
        if(cur_err_mdl[i] > 0) {
            if(i<L_seg||i>R_seg){
                f_read_seq[i] = (f_read_seq[i]+cur_err_mdl[i])%4;
                fprintf(stderr, "err[%u] = %u->%u\n", i, c, f_read_seq[i]);
            }
        }
    }
*/
    //f_read_seq[112] = (f_read_seq[112]+1)%4;  
    //for(i = 112; i < read_len-1; ++i ){
/*  
    for(i = indel_pos; i < read_len-1; ++i ){
        //f_read_seq[88] = (f_read_seq[88]+1)%4;
        f_read_seq[i] = f_read_seq[i+1];
    }
*/
    fprintf(stderr, "indel_pos = %u, indel = %d\n", indel_pos, 1);
/*    
    for(i = read_len-1; i >=indel_pos; --i ){
        //f_read_seq[88] = (f_read_seq[88]+1)%4;
        //f_read_seq[112] = (f_read_seq[112]+1)%4;
        f_read_seq[i] = f_read_seq[i-1];
    }
*/
/*  
    for(i=L_seg; i <90; ++i) {
        if(i == indel_pos) fprintf(stderr, "- ");
        fprintf(stderr, "%u ", f_read_seq[i]);
    }
*/
    for(i=110; i <=126; ++i) {
        fprintf(stderr, "%u ", f_read_seq[i]);
    }
    
    fprintf(stderr, "\n");
    fprintf(stderr,"%u, ------------------\n", __LINE__);
    
/*  
    fprintf(stderr, "var:\n");
    for(i=0; i <read_len; ++i) fprintf(stderr, "%u", r_read_seq[i]);
    fprintf(stderr, "\n");
*/
    return;
}

void set_err_mdl(uint8_t err_mdl[]){
    int i;
    for(i=0; i < LEN_READ; ++i) {
        int j = ((i+19)/20)%4;
        if(i%20 ==5) err_mdl[i+j] =((i+19)/20)%4; 
    }   
    return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void outUsedTime0(int flag, double *add_time)
{
    static struct timeval tp_start,tp_end;
    double time_used;
    static unsigned long start_cpu, end_cpu; // CPU time
    if(flag ==0){
        gettimeofday(&tp_start,NULL);
        //设置开始时间: CPU time
        start_cpu = clock();
    } else {
        //设置结束时间: CPU time
        end_cpu = clock();
        //print the system time
        gettimeofday(&tp_end,NULL);
        //print the CPU time
        //printf("Total CPU time used: %lf seconds.\n", (double)(end_cpu-start_cpu)/CLOCKS_PER_SEC);
        time_used = tp_end.tv_sec-tp_start.tv_sec + 
                (double)(tp_end.tv_usec-tp_start.tv_usec)/1000000;
        //printf("Total Used Time By gettimeofday(): %lf Seconds.\n", time_used);
        *add_time += time_used; 
        gettimeofday(&tp_start,NULL);
        //设置开始时间: CPU time
        start_cpu = clock();
    }
    return ; 
}
double outUsedTime()
{
//return 0.0;
 
    struct timeval tp;
    gettimeofday(&tp,NULL);
    double cur_time = tp.tv_sec + (double)(tp.tv_usec)/1000000;
   
    return cur_time; 

}
void outTime(int flag, double *add_time)
{
    static struct timeval tp_start,tp_end;
    double time_used;
    static unsigned long start_cpu, end_cpu; // CPU time
    if(flag ==0){
        gettimeofday(&tp_start,NULL);
        //设置开始时间: CPU time
        start_cpu = clock();

    } else {
        //设置结束时间: CPU time
        end_cpu = clock();
        //print the system time
        gettimeofday(&tp_end,NULL);
        //print the CPU time
        //printf("Total CPU time used: %lf seconds.\n", (double)(end_cpu-start_cpu)/CLOCKS_PER_SEC);
        time_used = tp_end.tv_sec-tp_start.tv_sec + (double)(tp_end.tv_usec-tp_start.tv_usec)/1000000;
        //printf("Total Used Time By gettimeofday(): %lf Seconds.\n", time_used);
        *add_time += time_used; 
        gettimeofday(&tp_start,NULL);
        //设置开始时间: CPU time
        start_cpu = clock();
    }
    return ; 
}
