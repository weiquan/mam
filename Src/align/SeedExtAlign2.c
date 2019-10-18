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
//#include "bwt.h"
//#include "setFileName.h"
//#include "SeedExtLib.h"
#include "query.h"
#include "ksw.h"
#include "lookup.h"
#include "seed.h"
#include "debug.h"
//#include "kseq.h"
//KSEQ_INIT(gzFile, gzread)
#define SEED_LEN 20
#define READ_LEN 20
#define FILE_PREFIX "../../BWT_INDEX/prefix"
//const char *read_file = "../../Reads/Read.fq";
#define MAX_READS 5


#define __get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)
int get_12mer_correct(uint32_t hash_boundry[], uint32_t end_kmer)
{
    int i;
    //uint32_t top = hash_boundry[0], bot = hash_boundry[LEN_SEED-1];
    uint32_t top = 0, bot = 11;
    uint32_t mid;
    int flag = 0;
    while(top < bot) {
        mid = (top+bot)/2;
        if(end_kmer == hash_boundry[mid]){
            flag = 1;
            break;
        } else if(end_kmer < hash_boundry[mid]){
            bot = mid;
        } else{
            top = mid+1;
        }
    }
    return flag;

}
int get_ext_idx_R(idx_t *fm_idx, query_t *query, int s_off[3], uint32_t idx[4]){ 
    uint8_t *seq = query->is_rev? query->rseq: query->seq;
    /*  
    if(s_off[2] + 32 > query->l_seq) {
        s_off[1] = s_off[0] + query->l_seq - s_off[2];
        //s_off[2] = query->l_seq;
    }
    if(s_off[1] < s_off[0] || s_off[0] < 0 
        || s_off[1] > query->l_seq) {
        printf("%u, error!!!\n", __LINE__);
        exit(1);
    }
    */
    int ed = s_off[1] - s_off[0];
    //if(ed > 32) ed = 32;
    uint32_t bg_p = s_off[0];
    uint32_t ed_p = s_off[1];
    uint32_t read_p = s_off[2];
   
    uint32_t t_i = idx[0];
    uint32_t b_i = idx[1];
    uint32_t m_i;
    uint32_t i;
    uint32_t pos = bwt_sa(fm_idx->bwt, t_i) + bg_p; 
    if(pos + ed > fm_idx->bwt->seq_len) {
        printf("%u, pos +ed > ref length!!!, pos = %u, ed = %d, ref_len = %u\n", __LINE__, pos, ed, fm_idx->bwt->seq_len);
        exit(1); 
    }
    uint64_t seq64, t_seq, b_seq, m_seq;
    seq64 = 0;
    for(i =0; i < ed; ++i){
        seq64 <<=2;
        seq64 |= seq[read_p+i];
    }
    t_seq = 0;
    for(i =0; i < ed; ++i){
        t_seq <<=2;
        t_seq |= __get_pac(fm_idx->pac, pos+i);
    }
    pos = bwt_sa(fm_idx->bwt, b_i) + bg_p;
    b_seq = 0;
    for(i =0; i < ed; ++i){
        b_seq <<=2;
        b_seq |= __get_pac(fm_idx->pac, pos+i);
    }
    int len = 0;
    if( t_seq == b_seq ) {
        for(i = 0; i < ed; ++i) {
            if(t_seq == seq64) {
                len = ed - i; 
                break;
            } 
            t_seq >>= 2;
            seq64 >>= 2;
        }
        idx[2] = t_i;
        idx[3] = b_i;  
        return len;
    }  
    int stat_t, stat_b;
    uint64_t s_buf = seq64, t_buf = t_seq, b_buf = b_seq;  

    uint64_t top[32][2], bot[32][2];
    uint32_t top_sch[2] = {}, bot_sch[2] = {};
    uint32_t top_n = 1, bot_n = 1;
    top[0][0] = t_i;
    top[0][1] = t_seq;
    bot[0][0] = b_i;
    bot[0][1] = b_seq;

    top_sch[0] = t_i;
    top_sch[1] = b_i;
    bot_sch[0] = t_i;
    bot_sch[1] = b_i;
    if(seq64 <= t_seq) {
        for(i = 0; i < ed; ++i) {
            if(t_buf == s_buf) {
                len = ed - i; 
                break;
            } 
            t_buf >>= 2;
            s_buf >>= 2;
        }
        stat_t = 1;
        stat_b = 0;
        if(len == 0) {
            idx[2] = idx[0];
            idx[3] = idx[0] -1; 
            return 0; 
        } 
    } 
    if(seq64 >= b_seq) {
        for(i = 0; i < ed; ++i) {
            if(b_buf == s_buf) {
                len = ed - i; 
                break;
            } 
            b_buf >>= 2;
            s_buf >>= 2;
        }
        stat_t = 0;
        stat_b = 1;
        if(len == 0) {
            idx[2] = idx[1]+1;
            idx[3] = idx[1]; 
            return 0; 
        }
    }
    uint64_t o_i, m_val;
    int find_flg;
    if(seq64 > t_seq && seq64 < b_seq) {
        o_i = b_i+1;
        find_flg = 0;
        while(t_i <= b_i){
            m_i = t_i+(b_i-t_i)/2;
            pos = bwt_sa(fm_idx->bwt, m_i) + bg_p;
            m_val = 0;
            for(i =0; i < ed; ++i){
                m_val <<=2;
                m_val |= __get_pac(fm_idx->pac, pos + i);
            }
            if(m_val == seq64){
                find_flg = 1;
                break;
            } else if(m_val < seq64){
                t_i = m_i+1;
                top[top_n][0] = m_i; 
                top[top_n][1] = m_val;
                top_n++; 
            } else{
                b_i = m_i-1;
                bot[bot_n][0] = m_i; 
                bot[bot_n][1] = m_val;
                bot_n++; 
            }
            o_i = m_i;
        }    

        int j;
        uint32_t buf;
        if(find_flg == 1) {
            s_buf = seq64;
            len = ed;
            top_sch[0] = m_i; 
            for(j = bot_n -1; j >=0; j--) {
                buf = bot[j][1];
                if(s_buf < buf)  {
                    top_sch[1] = bot[j][0]; 
                    break; 
                }
            }
            bot_sch[1] = m_i; 
            for(j = top_n -1; j >=0; j--) {
                buf = top[j][1];
                if(s_buf < buf)  {
                    bot_sch[0] = top[j][0]; 
                    break; 
                }
            }
            stat_t = 1;
            stat_b = 1; 
        } else {
            int len_t = 0, len_b = 0;
            t_buf = top[top_n-1][1];
            s_buf = seq64;
            for(i = 0; i < ed; ++i) {
                if(t_buf == s_buf) {
                    len_t = ed - i; 
                    break;
                } 
                t_buf >>= 2;
                s_buf >>= 2;
            }
            b_buf = bot[bot_n-1][1];
            s_buf = seq64;
            for(i = 0; i < ed; ++i) {
                if(b_buf == s_buf) {
                    len_b = ed - i; 
                    break;
                } 
                b_buf >>= 2;
                s_buf >>= 2;
            }
            if(len_t == len_b){
                len = len_t;
                s_buf = seq64 >>(2*(ed -len));
                top_sch[0] = top[top_n-1][0]; 
                for(j = bot_n -1; j >=0; j--) {
                    buf = bot[j][1] >> (2*(ed-len));
                    if(s_buf < buf) {
                        top_sch[1] = bot[j][0]; 
                        break; 
                    }
                }
        
                bot_sch[1] = bot[bot_n-1][0]; 
                for(j = top_n -1; j >=0; j--) {
                    buf = top[j][1] >> (2*(ed-len));
                    if(s_buf < buf)  {
                        bot_sch[0] = top[j][0]; 
                        break; 
                    }
                }
                stat_t = 1;
                stat_b = 1;
            } else if(len_t > len_b) {
                len = len_t; 
                s_buf = seq64 >>(2*(ed -len));
                
                top_sch[0] = top[top_n-1][0]; 
                top_sch[1] = bot[bot_n-1][0]; 
                bot_sch[1] = top[top_n-1][0]; 
                for(j = top_n -1; j >=0; j--) {
                    buf = top[j][1] >> (2*(ed-len));
                    if(s_buf < buf)  {
                        bot_sch[0] = top[j][0]; 
                        break; 
                    }
                }
                stat_t = 0;
                stat_b = 1;
            } else { // len_t < len_b
                len = len_b; 
                s_buf = seq64 >>(2*(ed -len));
                bot_sch[0] = top[top_n-1][0]; 
                bot_sch[1] = bot[bot_n-1][0]; 
                top_sch[0] = bot[bot_n-1][0]; 
                for(j = bot_n -1; j >=0; j--) {
                    buf = bot[j][1] >> (2*(ed-len));
                    if(s_buf < buf) {
                        top_sch[1] = bot[j][0]; 
                        break; 
                    }
                }
                stat_t = 1;
                stat_b = 0;
            }
        } 
        if(len == 0) {
            idx[2] = idx[1]+1;
            idx[3] = idx[1]; 
            return 0; 
        }
    } 
    if(stat_t == 1) {
        t_i = top_sch[0];
        b_i = top_sch[1];
        s_buf = seq64 >> ((ed-len)*2);
        o_i = b_i+1;
        while(t_i <= b_i){
            m_i = t_i+(b_i-t_i+1)/2;
            if(m_i == o_i) { break; }
            pos = bwt_sa(fm_idx->bwt, m_i) + bg_p;
            m_val = 0;
            for(i =0; i < len; ++i){
                m_val <<=2;
                m_val |= __get_pac(fm_idx->pac, pos + i);
            } 
            if(s_buf == m_val){
                t_i = m_i;
            } else if(s_buf < m_val){
                b_i = m_i;
            } else{
                printf("%u, error!!!!\n", __LINE__); 
                exit(1);
            }
            o_i = m_i;
        }
        if(stat_b == 0) {
            idx[2] = top_sch[0];
            idx[3] = t_i;
            return len; 
        } else {// stat_b == 1
            idx[3] = t_i; 
        } 
    }// end if(stat_t == 1) ++++++
    if(stat_b == 1) {
        t_i = bot_sch[0];
        b_i = bot_sch[1];
        s_buf = seq64 >> ((ed-len)*2);
        o_i = b_i+1;
        while(t_i <= b_i){
            m_i = t_i+(b_i-t_i)/2;
            if(m_i == o_i) break; 
            pos = bwt_sa(fm_idx->bwt, m_i) + bg_p;
            m_val = 0;
            for(i =0; i < len; ++i){
                m_val <<=2;
                m_val |= __get_pac(fm_idx->pac, pos + i);
            }
            if(s_buf == m_val ){
                b_i = m_i;
            } else if(s_buf > m_val){
                t_i = m_i+1;
            } else{
                printf("%u, error!!!!\n", __LINE__); 
                exit(1);
            }
            o_i = m_i;
        }
        if(stat_t == 0) {
            idx[2] = b_i;
            idx[3] = bot_sch[1];
            return len; 
        } else {
            idx[2] = b_i; 
        } 
    }
return len; 
//-----------------------------------
    
    if(idx[2] >= idx[0] && idx[2] <= idx[1]) {
    } else {
        printf("%u, error!!!\n", __LINE__);
        exit(1); 
    } 
    if(idx[3] >= idx[0] && idx[3] <= idx[1]) {
    } else {
        printf("%u, error!!!\n", __LINE__);
        exit(1); 
    }
    seq64 = 0;
    for(i =0; i < len; ++i){
        seq64 <<=2;
        seq64 |= seq[read_p+i];
    }
    uint32_t test_idx0 = idx[2] - 1;
    uint32_t test_idx1 = idx[3] + 1;
fprintf(stderr, "%u, idx[0] = %u, idx[1] = %u, idx[2] = %u, idx[3] = %u\n", __LINE__, idx[0], idx[1], idx[2], idx[3]);
    if(test_idx0 >= idx[0]){
        pos = bwt_sa(fm_idx->bwt, test_idx0) + bg_p;
        t_seq = 0;
        for(i =0; i < len; ++i){
            t_seq <<=2;
            t_seq |= __get_pac(fm_idx->pac, pos+i);
        }
        if(t_seq >= seq64) {
            printf("%u, t_seq = %lx, seq64 = %lx, error!!!\n", __LINE__, t_seq, seq64);
            exit(1); 
        }
    } else {
        test_idx0 = 0; 
    } 
        
    if(test_idx1 <= idx[1]){
        pos = bwt_sa(fm_idx->bwt, test_idx1) + bg_p;
        b_seq = 0;
        for(i =0; i < len; ++i){
            b_seq <<=2;
            b_seq |= __get_pac(fm_idx->pac, pos+i);
        }
        if(b_seq <= seq64) {
            printf("%u, error!!!\n", __LINE__);
            exit(1); 
        }
    } else {
        test_idx1 = 0; 
    } 
    
fprintf(stderr, "%u, test_idx1 = %u\n", __LINE__, test_idx1); 
    if(len < 32) {
        uint8_t seq_c;
        uint32_t k; 
        seq64 = 0;
        for(i = 0; i < len; ++i) {
            seq_c = seq[read_p+i];
            for(k = idx[2]; k <= idx[3]; ++k){
                pos = bwt_sa(fm_idx->bwt, k) + bg_p; 
                uint8_t c = __get_pac(fm_idx->pac, pos+i);
                
        fprintf(stderr, "%u, idx = %u, c = %u, seq_c = %u\n", __LINE__, k, c, seq_c);
                if(c != seq_c) {
                    printf("%u, error!!!\n", __LINE__);
                    exit(1);
                }
            }
        }
        seq_c = seq[read_p+len];
       
        for(k = idx[2]; k <= idx[3]; ++k){
            pos = bwt_sa(fm_idx->bwt, k) + bg_p; 

            uint8_t c = __get_pac(fm_idx->pac, pos+len);
fprintf(stderr, "%u, idx = %u, c = %u, seq_c = %u\n", __LINE__, k, c, seq_c);
            if(c == seq_c) {
                
fprintf(stderr, "%u, idx = %u, c = %u, seq_c = %u\n", __LINE__, k, c, seq_c);
                uint64_t tmp = 0;
                for(i = 0; i < len+1; ++i) {
                    c = __get_pac(fm_idx->pac, pos+i);
                    fprintf(stderr, "%u ", c);  
                }
                fprintf(stderr, "\n");  
                for(i = 0; i < len+1; ++i) {
                    c = seq[read_p+i];
                    fprintf(stderr, "%u ", c);  
                }
                fprintf(stderr, "\n");  
                printf("%u, error!!!\n", __LINE__);
                exit(1);
            }
        }
    }
    
   
    fprintf(stderr, "%u\n", __LINE__);
    return len;
}

/*  
int get_8mer_R_idx(idx_t *fm_idx, uint8_t find[],  uint16_t ext_R[24], uint32_t buf[6]){ 
    uint32_t i;
    uint32_t bgn_row = buf[0];
    uint32_t end_row = buf[1];
    uint32_t t_i = buf[2];
    uint32_t b_i = buf[3];
    int find_row = buf[4]; 
    int flag = 0;

    uint32_t o_i = b_i+1;
    while(t_i <= b_i){
        uint32_t m_i = t_i+(b_i-t_i)/2;
        uint32_t pos = bwt_sa(fm_idx->bwt, m_i);
        uint32_t m_val = 0;
        for(i =0; i <8; ++i){
            m_val <<=2;
            m_val |= __get_pac(fm_idx->pac, pos+12+i);
        }
        if(m_val == ext_R[end_row]){
            find[find_row++] = end_row;
            break;
        } else if(m_val < ext_R[end_row]){
            t_i = m_i+1; 
        } else{
            b_i = m_i-1;
        }
        o_i = m_i;
    }

    flag = find_row - buf[4];
    buf[4] = find_row;
    return flag;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
//查找右段8mer当中近似序列函数-----------------------------------------------------
//num_var_seq是变异序列数目

int get_seed_8mer_seq_R(idx_t *fm_idx, seed_t *seed, int num_var_seq)
{
    uint32_t ref_seq_len = fm_idx->bwt->seq_len;
    int i, j;
    //24种子变异序列已经产生，排序
    uint32_t buf[6];
    uint32_t stk_array[24][6];
    uint32_t seq12 = seed->hash_12mer;
    uint16_t *ext_R = seed->sort_buf_R;
    uint8_t *find = seed->find_buf_R;
    uint32_t find_row = 0;
    seed->num_find_R = find_row;
     
    uint32_t bgn_idx = fm_idx->fastmap->item[seq12];
    uint32_t pos = bwt_sa(fm_idx->bwt, bgn_idx);
    while(pos +20 > ref_seq_len){
        ++bgn_idx;
        pos = bwt_sa(fm_idx->bwt, bgn_idx);
    }

    uint16_t seq16 = 0;
    for(i =0; i <8; ++i){
        seq16 <<=2;
        seq16 |= __get_pac(fm_idx->pac, pos+12+i);
    }
    uint16_t bgn_seq = seq16;
    uint32_t end_idx = fm_idx->fastmap->item[seq12+1]-1;
    pos = bwt_sa(fm_idx->bwt, end_idx);
    while(pos +20 > ref_seq_len){
        end_idx--;
        pos = bwt_sa(fm_idx->bwt, end_idx);
    }
    seq16 = 0;
    for(i =0; i <8; ++i){
        seq16 <<=2;
        seq16 |= __get_pac(fm_idx->pac, pos+12+i);
    }
    uint16_t end_seq = seq16;
    uint32_t bgn_row = 0, end_row=0; 

    for(i =0; i<num_var_seq; ++i) {
        if(bgn_seq > ext_R[i]) {
        } else if(bgn_seq == ext_R[i]){
            find[find_row++] = i;
            seed->num_find_R = find_row;
            buf[4] = find_row;
            break; 
        } else{
            break;      
        }
    }
    if(i == num_var_seq) return 0;
    bgn_row = i;
 
    for(i = num_var_seq-1; i >= 0; --i) {
        if(end_seq < ext_R[i]) {
        } else if(end_seq == ext_R[i]){
            find[find_row++] = i;
            seed->num_find_R = find_row;
            buf[4] = find_row;
            break; 
        } else{
            break;      
        }
    }

    if(i < 0) return 0;
    end_row = i;
    

    buf[0] = bgn_row;
    buf[1] = end_row;
    buf[4] = find_row; 

    uint32_t t_i, m_i, b_i, o_i;
    int stk_num = 0;
    if(end_row < bgn_row) { 
        return 0;
    } else if(bgn_row == end_row){//只有一个, 进行二分查找
        buf[2] = bgn_idx;
        buf[3] = end_idx;
        get_8mer_R_idx(fm_idx, find, ext_R, buf); 
        find_row = buf[4]; 
    }
    //以下代码end_row - bgn_row >=2时运行
    stk_array[stk_num][0] = bgn_row;
    stk_array[stk_num][1] = bgn_row;
    stk_array[stk_num][2] = bgn_idx;
    stk_array[stk_num][3] = bgn_seq;
    stk_array[stk_num][4] = bgn_idx;
    stk_array[stk_num][5] = bgn_seq;
    stk_num++;
    stk_array[stk_num][0] = end_row;
    stk_array[stk_num][1] = end_row;
    stk_array[stk_num][2] = end_idx;
    stk_array[stk_num][3] = end_seq;
    stk_array[stk_num][4] = end_idx;
    stk_array[stk_num][5] = end_seq;

    uint32_t old_idx = end_idx+1; 
    
    uint32_t mid_idx = bgn_idx/2+end_idx/2+(bgn_idx%2+end_idx%2)/2;
    while(1){
        if(stk_num == 0){break;}
        bgn_idx = stk_array[stk_num-1][4];
        end_idx = stk_array[stk_num][2];   
     
        if(bgn_idx+1>= end_idx) {

            if(old_idx == mid_idx) {
                --stk_num;
                continue; 
            }  
        }

        bgn_row = stk_array[stk_num-1][1];
        end_row = stk_array[stk_num][0];
        bgn_seq = stk_array[stk_num-1][5];
        end_seq = stk_array[stk_num][3];
 
        mid_idx = bgn_idx/2+end_idx/2+(bgn_idx%2+end_idx%2)/2;
        old_idx = mid_idx;
        pos = bwt_sa(fm_idx->bwt, mid_idx);
        seq16 = 0;
        for(i =0; i <8; ++i){
            seq16 <<=2;
            seq16 |= __get_pac(fm_idx->pac, pos+12+i);
        }
        uint16_t mid_seq = seq16;
        if(mid_seq == bgn_seq){
            stk_array[stk_num-1][4] = mid_idx;
            stk_array[stk_num-1][5] = mid_seq;
            continue;
        } 
        if(mid_seq == end_seq){
            stk_array[stk_num][2] = mid_idx;
            stk_array[stk_num][3] = mid_seq;
            continue; 
        }
        if(mid_seq < ext_R[bgn_row]){
            stk_array[stk_num-1][4] = mid_idx;
            stk_array[stk_num-1][5] = mid_seq;
            continue;
        }
        if(mid_seq > ext_R[end_row]){
            stk_array[stk_num][2] = mid_idx;
            stk_array[stk_num][3] = mid_seq;
            continue; 
        } 
        if(mid_seq == ext_R[bgn_row]){
            find[find_row] = bgn_row;  
            ++find_row; 
            seed->num_find_R = find_row;
            buf[4] = find_row;
            bgn_row++; 
            stk_array[stk_num-1][4] = mid_idx;
            stk_array[stk_num-1][5] = mid_seq;
            stk_array[stk_num-1][1] = bgn_row;
            continue; 
        }
        if(mid_seq == ext_R[end_row]){
            find[find_row] = end_row;  
            ++find_row; 
            seed->num_find_R = find_row;
            buf[4] = find_row;
            end_row--;
            stk_array[stk_num][2] = mid_idx;
            stk_array[stk_num][3] = mid_seq;
            stk_array[stk_num][0] = end_row;
            continue; 
        }

        int flag = 0;
        int t_row = bgn_row;
        int b_row = end_row; 
        int m_row;

//fprintf(stderr, "mid_idx = %u, mid_seq = %u\n", mid_idx,  mid_seq); 
        int m_row_bg, m_row_ed;        
        
        while(1) {
            m_row = (t_row+b_row)/2;
            if(mid_seq == ext_R[m_row]){
                flag = 1;
                if(m_row>0 ) m_row_bg = m_row-1;
                else{ 
                    m_row_bg = m_row; 
                }
                if(m_row < num_var_seq-1 ) m_row_ed = m_row+1;
                else{ 
                    m_row_ed = m_row; 
                } 
                break;
            } else if(mid_seq < ext_R[m_row]){
                b_row = m_row;
            } else{
                t_row = m_row;
            }
            if(t_row+1 == b_row) {
                m_row_bg = t_row;
                m_row_ed = b_row;
                break;   
            }
        }
        if(flag == 1){
            find[find_row] = m_row;  
            ++find_row; 
            seed->num_find_R = find_row;
            buf[4] = find_row;
        } else{        }
        bgn_idx = stk_array[stk_num-1][4];
        end_idx = stk_array[stk_num][2];   
 
        int top_flag = 0; 
        if(m_row_bg > bgn_row){
            top_flag = 1;
        }else{//可能需要二分查找代码 
            buf[0] = bgn_row;
            buf[1] = m_row_bg;
            buf[2] = bgn_idx;
            buf[3] = mid_idx;
            get_8mer_R_idx(fm_idx, find, ext_R, buf); 
            find_row = buf[4]; 
        }  
        int bot_flag = 0;
        if(m_row_ed< end_row){
           bot_flag = 1;
        } else{//可能需要二分查找代码 
            buf[0] = m_row_ed;
            buf[1] = end_row;
            buf[2] = mid_idx;
            buf[3] = end_idx;
            get_8mer_R_idx(fm_idx, find, ext_R, buf);
            find_row = buf[4]; 
        }    
       if(top_flag+bot_flag == 2){//上下两断信息不变，插入中间分段信息
            for(i = 0; i < 6; ++i) stk_array[stk_num+1][i] = stk_array[stk_num][i];
            stk_array[stk_num][0] = m_row_bg;
            stk_array[stk_num][1] = m_row_ed;
            stk_array[stk_num][2] = mid_idx;
            stk_array[stk_num][3] = mid_seq;
            stk_array[stk_num][4] = mid_idx;
            stk_array[stk_num][5] = mid_seq;
            stk_num++;
        } else if(top_flag +bot_flag == 1){
            if(top_flag == 1) {//修正下段信息
                stk_array[stk_num][0] = m_row_bg;                
                stk_array[stk_num][2] = mid_idx;
                stk_array[stk_num][3] = mid_seq;
            } 
            if(bot_flag == 1) {//修正上段信息
                stk_array[stk_num-1][1] = m_row_ed;                
                stk_array[stk_num-1][4] = mid_idx;
                stk_array[stk_num-1][5] = mid_seq;
            } 
        } else{
            --stk_num;
        }
    }//end while(1)++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    seed->num_find_R = find_row;
    return find_row;
}//end get_seed_8mer_R()----------------------------------------------------


//查找右段8mer当中近似序列函数-----------------------------------------------------
int get_seed_8mer_R(idx_t *fm_idx, seed_t *seed)
{
    uint32_t ref_seq_len = fm_idx->bwt->seq_len;
//fprintf(stderr, "\n%u, ", __LINE__);
    int i, j;
    //24种子变异序列已经产生，排序
    uint32_t buf[6];
    uint32_t stk_array[24][6];
    uint32_t seq12 = seed->hash_12mer;
    uint16_t *ext_R = seed->sort_buf_R;
    uint8_t *find = seed->find_buf_R;
    uint32_t find_row = 0;
    seed->num_find_R = find_row;
    

//fprintf(stderr, "\n%u, ", __LINE__);
    uint32_t bgn_idx = fm_idx->fastmap->item[seq12];
    uint32_t pos = bwt_sa(fm_idx->bwt, bgn_idx);
    while(pos +20 > ref_seq_len){
        ++bgn_idx;
        pos = bwt_sa(fm_idx->bwt, bgn_idx);
    }

    uint16_t seq16 = 0;
    for(i =0; i <8; ++i){
        seq16 <<=2;
        seq16 |= __get_pac(fm_idx->pac, pos+12+i);
    }
    uint16_t bgn_seq = seq16;


//fprintf(stderr, "\n%u, ", __LINE__);
    uint32_t end_idx = fm_idx->fastmap->item[seq12+1]-1;
    pos = bwt_sa(fm_idx->bwt, end_idx);
//fprintf(stderr, "bgn_idx = %u, end_idx = %u\n", bgn_idx, end_idx);
    while(pos +20 > ref_seq_len){
        end_idx--;
        pos = bwt_sa(fm_idx->bwt, end_idx);
    }
//fprintf(stderr, "pos = %u, ref_seq_len = %u", pos, ref_seq_len);
    seq16 = 0;
    for(i =0; i <8; ++i){
        seq16 <<=2;
        seq16 |= __get_pac(fm_idx->pac, pos+12+i);
    }
    uint16_t end_seq = seq16;
    uint32_t bgn_row = 0, end_row=0; 

    for(i =0; i<24; ++i) {
        //if(bgn_seq < ext_R[i]) {
        if(bgn_seq > ext_R[i]) {
        
        } else if(bgn_seq == ext_R[i]){
            find[find_row++] = i;

//fprintf(stderr, "\n%u, i = %d\n ", __LINE__, i);
//printf("%u, find_row = %u, i = %u\n", __LINE__, find_row, i);
            seed->num_find_R = find_row;
            buf[4] = find_row;
            break; 
        } else{
            break;      
        }
    }
    if(i == 24) return 0;
    bgn_row = i;
 
//fprintf(stderr, "\n%u, ", __LINE__);
//fprintf(stderr, "end_seq = %u, ext_R[23] = %u\n", end_seq, ext_R[23]);
    for(i =23; i>=0; --i) {
        if(end_seq < ext_R[i]) {
        
        } else if(end_seq == ext_R[i]){
            find[find_row++] = i;
            seed->num_find_R = find_row;
            buf[4] = find_row;
            break; 
        } else{
            break;      
        }
    }

//fprintf(stderr, "\n%u, i = %d\n ", __LINE__, i);
    if(i < 0) return 0;
    end_row = i;
    
//fprintf(stderr, "\n%u, ", __LINE__);
    buf[0] = bgn_row;
    buf[1] = end_row;
    
    buf[4] = find_row; 

    uint32_t t_i, m_i, b_i, o_i;
    int stk_num = 0;
    if(end_row < bgn_row) { 
        return 0;
    } else if(bgn_row == end_row){//只有一个, 进行二分查找
        buf[2] = bgn_idx;
        buf[3] = end_idx;
        get_8mer_R_idx(fm_idx, find, ext_R, buf); 
        find_row = buf[4]; 
    }


    //以下代码end_row - bgn_row >=2时运行
    stk_array[stk_num][0] = bgn_row;
    stk_array[stk_num][1] = bgn_row;
    stk_array[stk_num][2] = bgn_idx;
    stk_array[stk_num][3] = bgn_seq;
    stk_array[stk_num][4] = bgn_idx;
    stk_array[stk_num][5] = bgn_seq;
    stk_num++;
    stk_array[stk_num][0] = end_row;
    stk_array[stk_num][1] = end_row;
    stk_array[stk_num][2] = end_idx;
    stk_array[stk_num][3] = end_seq;
    stk_array[stk_num][4] = end_idx;
    stk_array[stk_num][5] = end_seq;
    uint32_t old_idx = end_idx+1; 
    
    uint32_t mid_idx = bgn_idx/2+end_idx/2+(bgn_idx%2+end_idx%2)/2;
    while(1){
        if(stk_num == 0){break;}
        bgn_idx = stk_array[stk_num-1][4];
        end_idx = stk_array[stk_num][2];   
        if(bgn_idx+1>= end_idx) {
        if(old_idx == mid_idx) {
                --stk_num;
                continue; 
            }  
        };

        bgn_row = stk_array[stk_num-1][1];
        end_row = stk_array[stk_num][0];
        bgn_seq = stk_array[stk_num-1][5];
        end_seq = stk_array[stk_num][3];
 
        mid_idx = bgn_idx/2+end_idx/2+(bgn_idx%2+end_idx%2)/2;
        old_idx = mid_idx;
        pos = bwt_sa(fm_idx->bwt, mid_idx);
        seq16 = 0;
        for(i =0; i <8; ++i){
            seq16 <<=2;
            seq16 |= __get_pac(fm_idx->pac, pos+12+i);
        }
        uint16_t mid_seq = seq16;
        
//fprintf(stderr, "mid_seq = %u, mid_idx = %u\n",  mid_seq, mid_idx); 
        if(mid_seq == bgn_seq){
            stk_array[stk_num-1][4] = mid_idx;
            stk_array[stk_num-1][5] = mid_seq;

            continue;
        } else if(mid_seq == end_seq){
            stk_array[stk_num][2] = mid_idx;
            stk_array[stk_num][3] = mid_seq;

            continue; 
        }
        if(mid_seq < ext_R[bgn_row]){
            stk_array[stk_num-1][4] = mid_idx;
            stk_array[stk_num-1][5] = mid_seq;

//fprintf(stderr, "\n%u, ", __LINE__);
            continue;
        } else if(mid_seq > ext_R[end_row]){
            stk_array[stk_num][2] = mid_idx;
            stk_array[stk_num][3] = mid_seq;


//fprintf(stderr, "\n%u, ", __LINE__);
            continue; 
        }

//fprintf(stderr, "\n%u, ", __LINE__);
        int flag = 0;
        int t_row = bgn_row;
        //int b_row = end_row+1; 
        int b_row = end_row; 
        int m_row;

//fprintf(stderr, "mid_idx = %u\n",  mid_seq); 
        int m_row_bg, m_row_ed;        
        while(1) {
            m_row = (t_row+b_row)/2;
            
            if(mid_seq == ext_R[m_row]){
                flag = 1;
                if(m_row>0 ) m_row_bg = m_row-1;
                else{ 
                    m_row_bg = m_row; 
                }
                if(m_row <23 ) m_row_ed = m_row+1;
                else{ 
                    m_row_ed = m_row; 
                } 
                break;
            } else if(mid_seq < ext_R[m_row]){
                b_row = m_row;
            } else{
                t_row = m_row;
            }
            if(t_row +1== b_row) {
                m_row_bg = t_row;
                m_row_ed = b_row;
                break;   
            }
        }
//fprintf(stderr, "m_row = %u\n", m_row);        
        
      
        if(flag == 1){
            find[find_row] = m_row;  
            
            ++find_row; 
            seed->num_find_R = find_row;
            buf[4] = find_row;
            
//printf("%u, find_row = %u, m_row = %u\n", __LINE__, find_row, m_row);
            //m_row_bg = m_row-1;
            //m_row_ed = m_row+1;
        } else{
            //m_row_bg = m_row;
            //m_row_ed = m_row+1; 
        
//printf("%u, find_row = %u, m_row = %u\n", __LINE__, find_row, m_row);
        }


//fprintf(stderr, "\n%u, ", __LINE__);
        bgn_idx = stk_array[stk_num-1][4];
        end_idx = stk_array[stk_num][2];   
 
        int top_flag = 0; 
        //if(m_row-1> bgn_row){
        if(m_row_bg > bgn_row){

//fprintf(stderr, "\n%u, ", __LINE__);
            top_flag = 1;
        }else{//可能需要二分查找代码 
      
//fprintf(stderr, "\n%u, ", __LINE__);
            buf[0] = bgn_row;
            buf[1] = m_row_bg;
            buf[2] = bgn_idx;
            buf[3] = mid_idx;
            get_8mer_R_idx(fm_idx, find, ext_R, buf); 
            find_row = buf[4]; 
        }  
        int bot_flag = 0;
        if(m_row_ed< end_row){
       
//fprintf(stderr, "\n%u, ", __LINE__);
           bot_flag = 1;
        } else{//可能需要二分查找代码 

//fprintf(stderr, "\n%u, ", __LINE__);
//fprintf(stderr, "m_row_ed = %u, end_row = %u, mid_idx = %u, end_idx = %u\n", m_row_ed, end_row, mid_idx, end_idx);
            buf[0] = m_row_ed;
            buf[1] = end_row;
            buf[2] = mid_idx;
            buf[3] = end_idx;
            get_8mer_R_idx(fm_idx, find, ext_R, buf); 
            find_row = buf[4]; 
        }

//fprintf(stderr, "\n%u, ", __LINE__);
       if(top_flag+bot_flag == 2){//上下两断信息不变，插入中间分段信息
            for(i = 0; i < 6; ++i) stk_array[stk_num+1][i] = stk_array[stk_num][i];
            stk_array[stk_num][0] = m_row_bg;
            stk_array[stk_num][1] = m_row_ed;
            stk_array[stk_num][2] = mid_idx;
            stk_array[stk_num][3] = mid_seq;
            stk_array[stk_num][4] = mid_idx;
            stk_array[stk_num][5] = mid_seq;
            stk_num++;

//fprintf(stderr, "\n%u, ", __LINE__);

        } else if(top_flag +bot_flag == 1){
            if(top_flag == 1) {//修正下段信息
                stk_array[stk_num][0] = m_row_bg;                
                //stk_array[stk_num][1] = m_row_ed;                
                stk_array[stk_num][2] = mid_idx;
                stk_array[stk_num][3] = mid_seq;

//printf("%u, find_row = %u, m_row = %u\n", __LINE__, find_row, m_row);
            } 
            if(bot_flag == 1) {//修正上段信息

                stk_array[stk_num-1][1] = m_row_ed;                

                stk_array[stk_num-1][4] = mid_idx;
                stk_array[stk_num-1][5] = mid_seq;
//printf("%u, find_row = %u, m_row = %u\n", __LINE__, find_row, m_row);
            } 
        } else{
            --stk_num;
//printf("%u, find_row = %u, m_row = %u\n", __LINE__, find_row, m_row);
        }
    }//end while(1)++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //seed->num_find_R = find_row;
    seed->num_find_R = buf[4];
    return find_row;
//if (stk_num >0)fprintf(stderr, "[%s, %u]: stk_num = %u, bg_row = %u, ed_row = %u\n", __func__, __LINE__, stk_num, stk_array[stk_num-1][1], stk_array[stk_num][0]);
//fprintf(stderr, "[%s, %u]: stk[0][4] = %u, stk[1][2] = %u\n", __func__, __LINE__, stk_array[0][4], stk_array[1][2]);
}

int get_seed_8mer_L(idx_t *fm_idx, seed_t *seed)
{
    int64_t i, j;
    //24种子变异序列已经产生，排序

    uint32_t seq12 = seed->hash_12mer;
    uint16_t *ext_R = seed->sort_buf_L;
    uint8_t *find = seed->find_buf_L;
    uint32_t find_row = 0;
    seed->num_find_L = find_row;
    
    
    uint32_t bgn_idx = fm_idx->fastmap->item[seq12];
    uint32_t pos = bwt_sa(fm_idx->bwt, bgn_idx);
    uint16_t seq16 = 0;
    for(i =0; i <8; ++i){
        seq16 <<=2;
        seq16 |= __get_pac(fm_idx->pac, pos+12+i);
    }
    uint16_t bgn_seq = seq16;

    uint32_t end_idx = fm_idx->fastmap->item[seq12+1]-1;
    for(i = bgn_idx; i <= end_idx; ++i) {
        pos = bwt_sa(fm_idx->bwt, i);
        seq16 = 0;
        if(pos <=8) continue;  
        for(j =0; j <8; ++j){
            seq16 <<=2;
            seq16 |= __get_pac(fm_idx->pac, pos-8+j);
        }
        //用seq16在ext_R中进行二分查找,如果有追加find

        int t_i = 0, b_i = 23, m_i, is_aln=0;
        while(t_i <= b_i){
            m_i = (t_i+b_i)/2;
            if(seq16 == ext_R[m_i]) {


                printf("seq16 = %u\n");
                for(j = 0; j < 24; ++j) printf("%u\t", ext_R[j]);
                printf("\n");
                int j;
                for(j=0;j < find_row;++j) {
                    if(find[j] == m_i) break;
                }
                if(find_row == 0 ||j < find_row) find[find_row++] = m_i;
                is_aln = 1;
                break;
            } else if(seq16 < ext_R[m_i]){
               b_i = m_i-1;
            } else{
                t_i = m_i+1;
            }
        }
    }
    seed->num_find_L = find_row;
    return find_row;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//生成8mer变异序列函数----------------------------------------------------------
//seed是种子序列，bg_i变异的起始位置，ed_i是变异的结束边界，flag:0和1分别是左右侧。
int gen_8mer_var_buf(seed_t *seed, int bg_i, int ed_i, int flag)
{
    int i, j;
    uint8_t *seq;     
    uint16_t *ext_16;
    uint8_t (*var_info)[3] = seed->var_info;
    if(bg_i<0 || ed_i>8) {
        printf("bg_i or ed_i out of range!\n");
        exit(1);
    }
    if(flag == 0) {
        seq = seed->seed_8mer_L; 
        ext_16 = seed->sort_buf_L;
    } else{
        seq = seed->seed_8mer_R; 
        ext_16 = seed->sort_buf_R;
    } 
    
    int bg_row =0, ed_row =(ed_i-bg_i)*3;
    ext_16[ed_row] = 0;
    for(i = 0; i < 8; ++i){
        ext_16[ed_row] <<= 2;
        ext_16[ed_row] |= seq[i];
    }
    for(i = 0; i < ed_row; ++i){
        ext_16[i] = ext_16[ed_row];
    }
    ed_row--;
    for(i = bg_i; i < ed_i; ++i){
        uint8_t ch = seq[i];
        int offset = 14-i*2; 
        for(j = 0; j < ch; ++j) {
            //ext_16[bg_row+j][i] = j;
            ext_16[bg_row+j] &= ~(3<<offset);
            ext_16[bg_row+j] |= j<<offset;
            var_info[bg_row+j][0] = i;
            var_info[bg_row+j][1] = j;
            var_info[bg_row+j][2] = ch;
        }
        for(j = ch+1; j< 4; ++j){
            //ext_16[ed_row-(3-j)][i] = j; 
            ext_16[ed_row-(3-j)] &= ~(3<<offset); 
            ext_16[ed_row-(3-j)] |= j<<offset; 
            var_info[ed_row-(3-j)][0] = i;
            var_info[ed_row-(3-j)][1] = j;
            var_info[ed_row-(3-j)][2] = ch;
        }
        bg_row += ch;  
        ed_row -= 3 - ch;
    }
    return ext_16[(ed_i-bg_i)*3]; 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
int gen_8mer_sort_varseq(seed_t *seed, int flag)
{
    int i, j;
    uint8_t *seq;     
    uint16_t *ext_16;
    if(flag == 0) {
        seq = seed->seed_8mer_L; 
        ext_16 = seed->sort_buf_L;
    } else{
        seq = seed->seed_8mer_R; 
        ext_16 = seed->sort_buf_R;
    } 
    ext_16[24] = 0;
    for(i = 0; i < 8; ++i){
        ext_16[24] <<= 2;
        ext_16[24] |= seq[i];
    }
    for(i = 0; i < 24; ++i){
        ext_16[i] = ext_16[24];
    }
    int bg_row =0, ed_row =23;
    for(i = 0; i < 8; ++i){
        uint8_t ch = seq[i];
        int offset = 14-i*2; 
        for(j = 0; j < ch; ++j) {
            //ext_16[bg_row+j][i] = j;
            ext_16[bg_row+j] &= ~(3<<offset);
            ext_16[bg_row+j] |= j<<offset;
        }
        for(j = ch+1; j< 4; ++j){
            //ext_16[ed_row-(3-j)][i] = j; 
            ext_16[ed_row-(3-j)] &= ~(3<<offset); 
            ext_16[ed_row-(3-j)] |= j<<offset; 
        }
        bg_row +=ch;  
        ed_row -=3-ch;
    }
    return 0; 
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
int get_seed_8mer_bwt_L(idx_t *fm_idx, uint32_t *hash_boundry, seed_t *seed)
{ 
    int i, j;
    uint32_t (*buf)[6] = seed->aln_buf;
    uint8_t *cur_seed = seed->seed_seq;
    uint32_t seq12 = lkt_seq2LktItem(cur_seed+8, 0, 11);
    uint32_t k, l;
    k = fm_idx->fastmap->item[seq12];
    l = fm_idx->fastmap->item[seq12+1]-1;
    l -= get_12mer_correct(hash_boundry, l);

    int row = 0; 
    for(i = 0; i < 8; ++i){//
        uint8_t cur_ch = cur_seed[7-i];

        for(j=0; j<3; ++j){
            //cur_seed[i] = (cur_seed[i]+1)%4; 
            cur_seed[7-i] = (cur_seed[7-i]+1)%4; 
           
            uint32_t kk=k, ll=l;
            //int num = bwt_match_exact_alt(fm_idx->bwt, 8-i, cur_seed+7-i, &kk, &ll);                 
            int num = bwt_match_exact_alt(fm_idx->bwt, 8-i, cur_seed, &kk, &ll);                 
            if(num > 0){ 
                buf[row][0] = 7-i; 
                buf[row][1] = (cur_ch+j+1)%4; 
                buf[row][2] = kk; 
                buf[row][3] = ll; 
                buf[row][4] = cur_ch;
                ++row;
            }
        } 
        cur_seed[7-i] = (cur_seed[7-i]+1)%4; 
    
        int num = bwt_match_exact_alt(fm_idx->bwt, 1, cur_seed+7-i, &k, &l);
        

        if(num == 0) break;        
    }
    seed->num_aln = row;
    return row;
}


#ifdef __TEST_SEED
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

if(0)
{
    int buf_256[256] = {};
    int buf_1[8] = {};
    int buf_2[8] = {};
    int buf_3[8] = {};
    int buf_4[8] = {};

    for(i=0; i < 0xFFFFFF; ++i){
        int num = fm_idx->fastmap->item[i+1] - fm_idx->fastmap->item[i];
        if(num <= 0xFF) { 
            ++buf_256[num];
            ++buf_1[num/32]; 
        } else if (num <= 0xFFFF){ 
        
            ++buf_2[num/(32*256)]; 
        } else if(num <= 0xFFFFFF){
        
            ++buf_3[num/(32*256*256)]; 
        } else{
        
            ++buf_4[0]; 
        }
    }
    for(i = 0; i < 256; ++i) printf("[ %3u , %8u ]\n", i, buf_256[i]);
    for(i=0; i < 8; ++i) printf("[%u, %4uk ] ", i, buf_1[i]/1000);
    printf("\n");
    for(i=0; i < 8; ++i) printf("[%u, %10u ] ", i, buf_2[i]);
    printf("\n");
    for(i=0; i < 8; ++i) printf("[%u, %10u ] ", i, buf_3[i]);
    printf("\n");
    for(i=0; i < 8; ++i) printf("[%u, %10u ] ", i, buf_4[i]);
    printf("\n");
    
    
    return 0;
}
//++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++

    //+++++++++++++++++++++++++++++++++++++++++++
	//1.readArry数组开辟空间，
	//2.Read,Seed等的状态信息
	//3.输入输出格式信息
	//4.其它
	//
	//+++++++++++++++++++++++++++++++++++++++++++


	uint8_t  *read_seq;
	uint8_t	 seed_seq[32];
    char     *buf_read;	


	uint32_t get_num;
	uint32_t read_id;

    uint32_t buf_algn[20]; 	
	int      Flg_Algn  = 0 ;
	uint32_t cur_read_num = 0;
    int read_num = 0;

    read_id = 0;
    //uint32_t n_seqs = fm_idx->bwt->seq_len-READ_LEN; 
uint32_t n_seqs = 1000000; 
 
  
    int n_rev_aln = 0;
    
    uint32_t sum_aln[128]= {};
    int seed_slc_num = 0;
    while(read_id < n_seqs) {
        //if(read_id %1000000 == 0) printf("read_id = %u, n_rev_aln= %u, percent = %f\n", read_id, n_rev_aln, (double)n_rev_aln/(double)read_id);
//++++++++++++++++++++++++++++++++++++++++++++++++++
        uint32_t __l =  READ_LEN;
        uint8_t ch[READ_LEN];
        int __i;
        uint32_t pos = read_id;
        for(__i= 0; __i < __l; ++__i){
            ch[__i] = __get_pac(fm_idx->pac, pos+__i); 
        } 
        uint32_t read_len = READ_LEN;
        read_seq = ch;
        uint8_t rev_seed[20];

        uint32_t k=0, l=fm_idx->bwt->seq_len;
        int i, j = 0, n_aln=0;
        uint8_t cur_seed[20];
 
        int is_var = 0;



    uint32_t err_idx_8 = read_id%4;
    uint32_t err_val = read_id%3+1; 
    ch[12+err_idx_8+4] = (ch[12+err_idx_8+4]+err_val)%4;  

    is_var = 1; 
if(is_var)
{
    seed_t *seed = calloc(1, sizeof(seed_t));
     
//printf("%u\n", __LINE__);
    seed->seed_seq = ch;
    
//printf("has_12mer = %u, %u\n",seed->hash_12mer, __LINE__);

   

    seed->hash_12mer = lkt_seq2LktItem(ch, 0, 11);
//int n_12mer = fm_idx->fastmap->item[seed->hash_12mer+1] - fm_idx->fastmap->item[seed->hash_12mer];
    seed->seed_8mer_R = ch+12;
//if(n_12mer < 256*4) {
        //gen_8mer_sort_varseq(seed, 1);
        //get_seed_8mer_R(fm_idx, seed);
       
        int bg_i = 5, ed_i = 8, num_seq;
        num_seq = (ed_i-bg_i)*3; 
        gen_8mer_var_buf(seed, bg_i, ed_i, 1); 
        get_seed_8mer_seq_R(fm_idx, seed, num_seq);
//seed_slc_num++;

        if(seed->num_find_R <= 0) {
            printf("read_id = %u not mapped\n", read_id);
         
            for(i = 12; i <20; i++){printf("%u", ch[i]); } 
            printf("\n"); 

            for(i = 0; i < 8; ++i){
                for(j = 0; j <3; ++j) {  
                    //printf("sort_buf_L[%u, %u] = %u\n", i, j, seed->sort_buf_L[i*3+j]); 
                    uint32_t seq = seed->sort_buf_R[i*3+j]; 
                    int k;
                    printf("%u\t", i*3+j);
                    for(k = 14; k >=0; k-=2){
                        printf("%u", (seq>>k)&3); 
                    } 
                    printf("\t%u\n", seq); 
                }
            }
            for(i=0; i < seed->num_find_R;++i){
                printf("find_row = %u\n", seed->find_buf_R[i]); 
            } 
            printf("------------------------\n\n");      
            exit(1);
        }
     

//} 
    //printf("12mer = %u\n", seed->hash_12mer);
    //printf("8mer = %u\n", lkt_seq2LktItem(ch, 12, 19));
    

if(read_id %1000000 == 0) printf("read_id = %u, seed_slc_num = %u, percent = %u\n", read_id, seed_slc_num, seed_slc_num/(1+read_id/100));
    free(seed);

}






        read_id++;   
    }        
//++++++++++++++++++++++++++++++++++++++

//++++++++++++++++++++++++++++++++++++++



    


  
//printf("read_id = %u, n_rev_aln= %u\n", read_id, n_rev_aln);
    return 1;

} //End ：main() +++++++++++++++++++++++++++++++++++++++++++ 

#endif
*/
