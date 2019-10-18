
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

#define  LEN_SEED  20 
#define  OFF_SEED  10 
#define  LEN_EXT   16 
#define  LEN_READ  150
#define  NUM_EXT  (((LEN_READ- LEN_SEED )/2)+OFF_SEED+ LEN_EXT-1)/LEN_EXT  
#define  MAX_SEED_NUM 600000
#define   LEN_FILE_NAME 100 

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
//种子处理标记数组FlgPos[]数组，
//NUM_EXT个扩展类别的种子序列区间标记数组FlgIdx[][NUM_EXT]数组;
//------------------------------------------------------------------

typedef struct {
	//int  numidxfiles = {6} ;
	int  numidxfiles;

	char (*capidx)[LEN_FILE_NAME]  ;
	char (*relat)[LEN_FILE_NAME]   ;
	char (*smbwt)[LEN_FILE_NAME]   ;
	char (*nxtpnt)[LEN_FILE_NAME]  ;
	char (*nxtflg)[LEN_FILE_NAME]  ;
	char (*smpos)[LEN_FILE_NAME]   ;
	char (*extidx)[LEN_FILE_NAME]   ;

	char (*comfile)[LEN_FILE_NAME] ;
	char (*seedidx)[LEN_FILE_NAME];

	char *jmpmod;
} FileName;

#define NUM_FILES 6
//NUM_FILE_NAME 
void setFileName(char *in_fname, FileName *out_fname,int flg );
