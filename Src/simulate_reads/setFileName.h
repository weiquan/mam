
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
#define  LEN_FILE_NAME 100 

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
#define NUM_FILES 6
#define REFSEQ_SIZE   1000
#define BWT_SUM_SIZE   1000/4



struct FileName{
	int  numidxfiles;//init 6

	char (*capidx)[LEN_FILE_NAME]  ;
	char (*relat)[LEN_FILE_NAME]   ;
	char (*smbwt)[LEN_FILE_NAME]   ;
	char (*nxtpnt)[LEN_FILE_NAME]  ;
	char (*nxtflg)[LEN_FILE_NAME]  ;
	char (*extidx)[LEN_FILE_NAME]   ;

	char (*comfile)[LEN_FILE_NAME] ;
	char (*seedidx)[LEN_FILE_NAME];

	//char (*jmpmod)[LEN_FILE_NAME];
    char *jmpmod;
};


/*需要传递的全局数据变量*/
typedef struct PubParmType
{
	uint32_t	*RefSeq;
	uint32_t	*I2Pos;
	uint32_t 	*Bwt ;
	uint32_t 	*Sum ;
}PubParm;
//NUM_FILE_NAME 
void setFileName(struct FileName *out_fname);

void  InitRefBwt(PubParm *pParm);


/*
void setFileName(char *in_fname, struct FileName *out_fname){
	//文件名初始化-----------------------------------------------------
printf("%s\n","setFileName(in_f,&f,0); beging Ok" );
	char jmpmod[LEN_FILE_NAME];
	char seedidx[LEN_FILE_NAME];
	char idxfile[LEN_FILE_NAME];
	char comfile[LEN_FILE_NAME];

	strcpy(jmpmod, "jmpmod" );
	strcpy(seedidx,"seedidx");
	strcpy(idxfile,"idxfile");	
	strcpy(comfile,"comfile");
	
	char (*pf)[LEN_FILE_NAME];
	struct FileName *f;
	
	int i;
	char ch[5];

	if(0==flg){
		if(NULL == (f = (struct FileName*)malloc(sizeof(struct FileName)))){
		    perror("error:[malloc(NUM_EXT*sizeof(FileName))]");
		    exit(1);
		}
	}

	if(NULL == f){
	    perror("error:[NULL == f]");
	    exit(1);
	}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (strcmp(in_fname, jmpmod) == 0) {
		if(NULL == (pf = (char(*)[LEN_FILE_NAME])malloc(LEN_FILE_NAME * sizeof(char)))){
		    perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
		    exit(1);
		}
		f->jmpmod  = pf; 
		strcpy (*(f->jmpmod),"./jmpmod/jmpmod");
		strcat (*(f->jmpmod),".txt");	
		out_fname = f;
		return;
    } 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (strcmp(in_fname, seedidx) == 0) {

		if(NULL == (pf = (char(*)[LEN_FILE_NAME])malloc(NUM_EXT*LEN_FILE_NAME * sizeof(char)))){
		    perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
		    exit(1);
		}
		f->seedidx = pf; 
		//NUM_EXT个输出文件seedidx_i.txt


		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->seedidx+i),"./seedidx/seedidx_");

			itoa(i,ch,10);
			strcat (*(f->seedidx+i),ch);	
			strcat (*(f->seedidx+i),".txt");


		}
		out_fname = f;
		return;
    } 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (strcmp(in_fname, comfile) == 0) {

		if(NULL == (pf = (char(*)[LEN_FILE_NAME])malloc(NUM_EXT*LEN_FILE_NAME * sizeof(char)))){
		    perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
		    exit(1);
		}
		f->comfile = pf; 
		//NUM_EXT个输出文件seedidx_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->comfile+i),"./seedidx/seedidx_");
			itoa(i,ch,10);
			strcat (*(f->comfile+i),ch);	
			strcat (*(f->comfile+i),".txt");	
		}
		out_fname = f;
		return;
    } 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (strcmp(in_fname, idxfile) == 0) {
		if(NULL == (pf = (char(*)[LEN_FILE_NAME])malloc(NUM_FILES * NUM_EXT * LEN_FILE_NAME * sizeof(char)))){
		    perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
		    exit(1);
		}
		int spc = (NUM_EXT*LEN_FILE_NAME);
		f->capidx  = pf + 0*spc; 
		f->relat   = pf + 1*spc; 
		f->smbwt   = pf + 2*spc; 
		f->nxtpnt  = pf + 3*spc; 
		f->nxtflg  = pf + 4*spc; 
		f->smpos   = pf + 5*spc; 

		//NUM_EXT个输出文件CapIdxFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->capidx+i),"./capidx/capidx_");
			itoa(i,ch,10);
			strcat (*(f->capidx+i),ch);	
			strcat (*(f->capidx+i),".txt");	
		}

		//NUM_EXT个输出文件RelatFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->relat+i),"./relat/relat_");
			itoa(i,ch,10);
			strcat (*(f->relat+i),ch);	
			strcat (*(f->relat+i),".txt");	
		}
		//NUM_EXT个输出文件SmBwtFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->smbwt+i),"./smbwt/smbwt_");
			itoa(i,ch,10);
			strcat (*(f->smbwt+i),ch);	
			strcat (*(f->smbwt+i),".txt");	
		}

		//NUM_EXT个输出文件NxtPntFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->nxtpnt+i),"./nxtpnt/nxtpnt_");
			itoa(i,ch,10);
			strcat (*(f->nxtpnt+i),ch);	
			strcat (*(f->nxtpnt+i),".txt");	
		}

		//NUM_EXT个输出文件NxtFlgFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->nxtflg+i),"./nxtflg/nxtflg_");
			itoa(i,ch,10);
			strcat (*(f->nxtflg+i),ch);	
			strcat (*(f->nxtflg+i),".txt");	
		}

		//NUM_EXT个输出文件SmPosFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->smpos+i),"./smpos/smpos_");
			itoa(i,ch,10);
			strcat (*(f->smpos+i),ch);	
			strcat (*(f->smpos+i),".txt");	
		}
		//NUM_EXT个输入文件ExtIdxFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->seedidx+i),"./extidx/extidx_");
			itoa(i,ch,10);
			strcat (*(f->seedidx+i),ch);	
			strcat (*(f->seedidx+i),".txt");	
		}
	
		out_fname = f;
		return;
    } 
    perror("error：there are not [strcmp(in_fname, idxfile) == 0] ");

	return;
}

*/

#define ROOT_PATH "../../Index/"

void setFileName(struct FileName *f){
    char jmpmod[LEN_FILE_NAME];
	char seedidx[LEN_FILE_NAME];
	char idxfile[LEN_FILE_NAME];
	char comfile[LEN_FILE_NAME];

	strcpy(jmpmod, "jmpmod" );
	strcpy(seedidx,"seedidx");
	strcpy(idxfile,"idxfile");	
	strcpy(comfile,"comfile");
	
	char (*pf)[LEN_FILE_NAME];
	
	int i;
	char ch[5];

    f->numidxfiles = 6;

	if(NULL == f){
	    perror("error:[NULL == f]");
	    exit(1);
	}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    //if (strcmp(in_fname, jmpmod) == 0) {
    if(NULL == (f->jmpmod  = (char*)malloc(LEN_FILE_NAME * sizeof(char)))){
        //perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
        

        exit(1);
    }

    strcpy (f->jmpmod,ROOT_PATH);
    strcat (f->jmpmod,"data/jmpmod/jmpmod");
    strcat (f->jmpmod,".txt");	
    

        //return;
    //}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if(NULL == (pf = (char(*)[LEN_FILE_NAME])malloc((NUM_EXT+1)*LEN_FILE_NAME * sizeof(char)))){
        perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
        exit(1);
    }

    f->seedidx = pf; 
    //NUM_EXT个输出文件seedidx_i.txt


    for(i=0;i<NUM_EXT+1; i++){
        
        strcpy (f->seedidx[i],ROOT_PATH);
        strcat (*(f->seedidx+i),"data/seedidx/seedidx_");
        sprintf(ch, "%d", i);  
        //itoa(i,ch,10);
        strcat (*(f->seedidx+i),ch);	
        strcat (*(f->seedidx+i),".txt");	
    }
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if(NULL == (f->comfile = (char(*)[LEN_FILE_NAME])malloc(NUM_EXT*LEN_FILE_NAME * sizeof(char)))){
        perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
        exit(1);
    }
    //f->comfile = pf; 
    //NUM_EXT个输出文件seedidx_i.txt
    for(i=0;i<NUM_EXT; i++){
        strcpy(f->comfile[i], ROOT_PATH);
        strcat (*(f->comfile+i),"data/comfile/comfile_");
        sprintf(ch, "%d", i);
        strcat (*(f->comfile+i),ch);	
        strcat (*(f->comfile+i),".txt");	
    }

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    if(NULL == (pf = (char(*)[LEN_FILE_NAME])malloc(NUM_FILES * (NUM_EXT+1) * LEN_FILE_NAME * sizeof(char)))){
        perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
        exit(1);
    }
    f->capidx  = pf + 0*NUM_EXT; 
    f->relat   = pf + 1*NUM_EXT; 
    f->smbwt   = pf + 2*NUM_EXT; 
    f->nxtpnt  = pf + 3*NUM_EXT; 
    f->nxtflg  = pf + 4*NUM_EXT; 
    f->extidx  = pf + 5*NUM_EXT; 

    //NUM_EXT个输出文件CapIdxFile_i.txt
    for(i=0;i<NUM_EXT+1; i++){
        strcpy(f->capidx[i], ROOT_PATH);
        strcat (f->capidx[i],"data/capidx/capidx_");
        sprintf(ch, "%d", i);
        strcat (f->capidx[i],ch);	
        strcat (f->capidx[i],".txt");	
    }

    //NUM_EXT个输出文件RelatFile_i.txt
    for(i=0;i<NUM_EXT+1; i++){
        strcpy(f->relat[i], ROOT_PATH);
        strcat (f->relat[i],"data/relat/relat_");
        sprintf(ch, "%d", i);
        strcat (f->relat[i],ch);	
        strcat (f->relat[i],".txt");	
    }
//printf("f->idxfile = %s\n",f->relat[0]);			

    //NUM_EXT个输出文件SmBwtFile_i.txt
    for(i=0;i<NUM_EXT+1; i++){
        strcpy(f->smbwt[i], ROOT_PATH);
        strcat (*(f->smbwt+i),"data/smbwt/smbwt_");
        sprintf(ch, "%d", i);
        strcat (*(f->smbwt+i),ch);	
        strcat (*(f->smbwt+i),".txt");	
    }

    //NUM_EXT个输出文件NxtPntFile_i.txt
    for(i=0;i<NUM_EXT+1; i++){
        strcpy(f->nxtpnt[i], ROOT_PATH);
        strcat (*(f->nxtpnt+i),"data/nxtpnt/nxtpnt_");
        sprintf(ch, "%d", i);
        strcat (*(f->nxtpnt+i),ch);	
        strcat (*(f->nxtpnt+i),".txt");	
    }

    //NUM_EXT个输出文件NxtFlgFile_i.txt
    for(i=0;i<NUM_EXT+1; i++){
        strcpy(f->nxtflg[i], ROOT_PATH);
        strcat (*(f->nxtflg+i),"data/nxtflg/nxtflg_");
        sprintf(ch, "%d", i);
        strcat (*(f->nxtflg+i),ch);	
        strcat (*(f->nxtflg+i),".txt");	
    }

    //NUM_EXT个输出文件SmPosFile_i.txt
    for(i=0;i<NUM_EXT+1; i++){
        strcpy(f->extidx[i], ROOT_PATH);
        strcat (*(f->extidx+i),"data/extidx/extidx_");
        sprintf(ch, "%d", i);
        strcat (*(f->extidx+i),ch);	
        strcat (*(f->extidx+i),".txt");	
    }
   
    return;
}
void getFileName(char *in_fname, struct FileName *f){
	//文件名初始化-----------------------------------------------------
printf("setFileName(in_f,f); 243: Ok\n" );

	char jmpmod[LEN_FILE_NAME];
	char seedidx[LEN_FILE_NAME];
	char idxfile[LEN_FILE_NAME];
	char comfile[LEN_FILE_NAME];

	strcpy(jmpmod, "jmpmod" );
	strcpy(seedidx,"seedidx");
	strcpy(idxfile,"idxfile");	
	strcpy(comfile,"comfile");
	
	char (*pf)[LEN_FILE_NAME];
	
	int i;
	char ch[5];

	if(NULL == f){
	    perror("error:[NULL == f]");
	    exit(1);
	}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    if (strcmp(in_fname, jmpmod) == 0) {
   
        strcpy(f->jmpmod, ROOT_PATH); 
        strcat (f->jmpmod,"data/jmpmod/jmpmod");
        strcat (f->jmpmod,".txt");	
    

        return;
    }


printf("jmpmod out3 ok\n");

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (strcmp(in_fname, seedidx) == 0) {
		if(NULL == (pf = (char(*)[LEN_FILE_NAME])malloc(NUM_EXT*LEN_FILE_NAME * sizeof(char)))){
		    perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
		    exit(1);
		}

		f->seedidx = pf; 
		//NUM_EXT个输出文件seedidx_i.txt


		for(i=0;i<NUM_EXT+1; i++){
			strcpy(f->seedidx[i], ROOT_PATH);
            strcat (*(f->seedidx+i),"data/seedidx/seedidx_");
      		sprintf(ch, "%d", i);  
			//itoa(i,ch,10);
			strcat (*(f->seedidx+i),ch);	
			strcat (*(f->seedidx+i),".txt");	
		}
		return;
    } 

printf("jmpmod out4 ok\n");
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (strcmp(in_fname, comfile) == 0) {

		if(NULL == (pf = (char(*)[LEN_FILE_NAME])malloc(NUM_EXT*LEN_FILE_NAME * sizeof(char)))){
		    perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
		    exit(1);
		}
		f->comfile = pf; 
		//NUM_EXT个输出文件seedidx_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy(f->comfile[i], ROOT_PATH);
            strcat (*(f->comfile+i),"data/comfile/comfile_");
			sprintf(ch, "%d", i);
			strcat (*(f->comfile+i),ch);	
			strcat (*(f->comfile+i),".txt");	
		}
		return;
   } 
 
printf("jmpmod out5 ok\n");
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (strcmp(in_fname, idxfile) == 0) {
		if(NULL == (pf = (char(*)[LEN_FILE_NAME])malloc(NUM_FILES * NUM_EXT * LEN_FILE_NAME * sizeof(char)))){
		    perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
		    exit(1);
		}
		f->capidx  = pf + 0*NUM_EXT; 
		f->relat   = pf + 1*NUM_EXT; 
		f->smbwt   = pf + 2*NUM_EXT; 
		f->nxtpnt  = pf + 3*NUM_EXT; 
		f->nxtflg  = pf + 4*NUM_EXT; 
		f->extidx  = pf + 5*NUM_EXT; 

		//NUM_EXT个输出文件CapIdxFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (f->capidx[i],"data/capidx/capidx_");
			sprintf(ch, "%d", i);
			strcat (f->capidx[i],ch);	
			strcat (f->capidx[i],".txt");	
		}

		//NUM_EXT个输出文件RelatFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (f->relat[i],"data/relat/relat_");
			sprintf(ch, "%d", i);
			strcat (f->relat[i],ch);	
			strcat (f->relat[i],".txt");	
		}
//printf("f->idxfile = %s\n",f->relat[0]);			
	
		//NUM_EXT个输出文件SmBwtFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->smbwt+i),"data/smbwt/smbwt_");
			sprintf(ch, "%d", i);
			strcat (*(f->smbwt+i),ch);	
			strcat (*(f->smbwt+i),".txt");	
		}

		//NUM_EXT个输出文件NxtPntFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->nxtpnt+i),"data/nxtpnt/nxtpnt_");
			sprintf(ch, "%d", i);
			strcat (*(f->nxtpnt+i),ch);	
			strcat (*(f->nxtpnt+i),".txt");	
		}

		//NUM_EXT个输出文件NxtFlgFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->nxtflg+i),"data/nxtflg/nxtflg_");
			sprintf(ch, "%d", i);
			strcat (*(f->nxtflg+i),ch);	
			strcat (*(f->nxtflg+i),".txt");	
		}

		//NUM_EXT个输出文件SmPosFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			strcpy (*(f->extidx+i),"data/extidx/extidx_");
			sprintf(ch, "%d", i);
			strcat (*(f->extidx+i),ch);	
			strcat (*(f->extidx+i),".txt");	
		}
		//NUM_EXT个输入文件ExtIdxFile_i.txt
		return;
    } 

printf("jmpmod out6 ok\n");
    perror("error：there are not [strcmp(in_fname, idxfile) == 0] ");
	return;

}


#include <sys/stat.h>  
 
uint64_t getFileSize(const char *file)  
{  
    uint64_t fsize = 0;      
    struct stat st;  
    if(stat(file, &st) == -1){  
        return fsize;  
    }else{  
        fsize =(uint64_t) st.st_size;  
    }  
    return fsize;  
} 
/*
 int checkDataFile(PubParm *pParm){
	// Check for  all genome files existence
	FILE *fp;
	//char strPath[FILENAME_LEN];
	//int path_len;
	struct stat st;
	char fileName[FILENAME_LEN];
	int path_len;
    char ch;
	char *strGenPath =	pParm->strGenPath;
	char *strTempPath =	pParm->strTemPath;
	char *strFileInfo =	pParm->strFileInfo;

	int i;
	int row;
	int fileID;

	strcpy(strGenPath, INPUT_GENOMES_PATH);
	pParm->genPathLen = strlen(strGenPath);
	if(stat(strGenPath, &st) == 0)
	{
		printf("the directory [ %s ] has already exist.\n", strGenomesPath);
	}else
	{
		printf("the directory [ %s ] has not exist.\n", strGenomesPath);
		return 0;
	}

	strcpy(strTemPath, TEMPER_PATH);
	pParm->temPathLen = strlen(strTemPath);
	if(stat(strTemPath, &st) == 0)
	{
		printf("the directory [ %s ] has already exist.\n", strTemPath);
	}else
	{
		//mkdir(strTemperPath, 0755); //make the folder 
	}

	strcpy(strGenFileInfo,GENOMES_FILE_INFO);
	pParm->genFileInfoLen = strlen(strGenFileInfo);
	if(stat(strGenFileInfo, &st) == 0)
	{
		printf("the [ %s ] has already exist.\n", strGenFileInfo);
	}else
	{
		printf("the [ %s ] has not exist.\n", strGenFileInfo);
		return 0;
	}

	fp = fopen(GENOMES_FILE_INFO,"r");

	while(( ch = fgetc(fp))!='\n');

	int genFilesNum=0;
	int curFileNameLen=0;

	while(feof(fp)==0)
	{
	    if((ch=fgetc( fp))!='\n' && ch!=' ')
	    {
              fileName[curFileNameLen]=ch;
              //fileName[GenomesFilesNum][curFileNameLen]=ch;
              curFileNameLen++;
		}else if(ch=='\n' && curFileNameLen>0)
		{
              for( i=0;i<curFileNameLen;i++)
              {
					pParm->fileName[genFilesNum*FILENAME_LEN+i]=fileName[i];
			  }
              pParm->fileNameLen[genFilesNum] = curFileNameLen;
              curFileNameLen=0;
              genFilesNum++;
		}
	}


	if(curFileNameLen>0)
	{
		for(i=0;i<curFileNameLen;i++)
        {
			row=genFilesNum*FILENAME_LEN;
			pParm->fileName[row+i]=fileName[i];
		}
		genFilesNum++;
	}

	pParm->genFilesNum=genFilesNum;
    fclose(fp);

	strcpy(fileName, INPUT_GENOMES_PATH);
	path_len = strlen(fileName);

	for(fileID=0; fileID<genFilesNum; fileID++)
	{
		row=fileID*FILENAME_LEN;
		curFileNameLen=pParm->fileNameLen[fileID];
        for(i=0;i<curFileNameLen;i++)
        {
			fileName[path_len+i]=pParm->fileName[row+i];
		}
        fileName[path_len+curFileNameLen]='\0';

		if(stat(fileName, &st) == 0)
		{
			printf("the directory [ %s ] has already exist.\n", fileName);
		}else
		{
			printf("the directory [ %s ] has not exist.\n", fileName);
			return 0;
		}
	}
	//------------------------------------------------------------------------------
	// set files Name
	//char  fileName[256];
    char  pathFile[256];
	char* tmpPath;
	char* sID;
	int name_len;
	int sid_len;
	strcpy(tmpPath, TEMPER_PATH);
	path_len = strlen(tmpPath) ;

	// set posFileW string value
	char posFileW[256] = {0};
	strcpy(fileName, "/posfile_w.bin");
	strcpy(pathFile, tmpPath);
	strcpy(pathFile+path_len, fileName);
	name_len = strlen(fileName) ;
	pathFile[path_len+name_len-1] ='\0';
	strcpy(posFileW,pathFile);
	pParm->posFileW=posFileW;

	// set posFileR string value
	char posFileR[256] = {0};
	strcpy(fileName, "/posfile_r.bin");
	strcpy(pathFile, tmpPath);
	strcpy(pathFile+path_len, fileName);
	name_len = strlen(fileName) ;
	pathFile[path_len+name_len-1] ='\0';
	strcpy(posFileR,pathFile);
	pParm->posFileR=posFileR;

	return 1;
}

*/

//#define REFSEQ_SIZE 10000;

void  InitRefBwt(PubParm *pParm){

//printf("%s\n", "InitRefBwt(PubParm *pParm):  begin OK!" ); // ++++++++++++++++++++++++
	uint32_t	*RefSeq;
	uint32_t	*I2Pos;
	uint32_t 	*Bwt ;
	uint32_t  	*Sum ;
/*	
	if(NULL == (pParm = (PubParm*)malloc(sizeof(PubParm)))){
	    perror("error...PubParm");
	    exit(1);
	}
*/
	if(NULL == (RefSeq = (uint32_t*)malloc(((REFSEQ_SIZE+15)/16)*sizeof(uint32_t)))){
	    perror("error...RefSeq");
	    exit(1);
	}

	if(NULL == (I2Pos = (uint32_t*)malloc(REFSEQ_SIZE*sizeof(uint32_t)))){
	    perror("error...I2Pos");
	    exit(1);
	}

	if(NULL == (Bwt = (uint32_t*)malloc(((REFSEQ_SIZE+15)/16)*sizeof(uint32_t)))){
	    perror("error...Bwt");
	    exit(1);
	}
	
	if(NULL == (Sum = (uint32_t*)malloc(BWT_SUM_SIZE*sizeof(uint32_t)))){
	    perror("error...Sum");
	    exit(1);
	}
    
    // 打开数据文件，分别读入数据到相应的数组

	pParm->RefSeq = RefSeq;
	pParm->I2Pos  = I2Pos;
	pParm->Bwt    = Bwt;
	pParm->Sum    = Sum;

	I2Pos[5] = 1000;

//printf("I2Pos[5] = %d\n", I2Pos[5] ); 
//printf("%s\n", "InitRefBwt(PubParm *pParm):  end OK!" ); // ++++++++++++++++++++++++
	return;
}




