#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "setFileName.h"
//void setFileName(char *in_fname, FileName *out_fname,int flg ){

#define ROOT_PATH "../Index/"
void setFileName(FileName *f, const char* prefix)
{
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
    sprintf(f->jmpmod, "%s_jmpmod.txt", prefix);

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
        sprintf(f->seedidx[i], "%s_seedidx_%d.txt", prefix, i);
    }
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if(NULL == (f->comfile = (char(*)[LEN_FILE_NAME])malloc(NUM_EXT*LEN_FILE_NAME * sizeof(char)))){
        perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
        exit(1);
    }
    //f->comfile = pf; 
    //NUM_EXT个输出文件seedidx_i.txt
    for(i=0;i<NUM_EXT; i++){
        sprintf(f->comfile[i], "%s_comfile_%d.txt", prefix, i);
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
    f->extidx   = pf + 5*NUM_EXT; 

    //NUM_EXT个输出文件CapIdxFile_i.txt
    for(i=0;i<NUM_EXT+1; i++){
      sprintf(f->capidx[i], "%s_capidx_%d.txt", prefix, i);
      sprintf(f->relat[i],  "%s_relat_%d.txt", prefix, i);
      sprintf(f->smbwt[i],  "%s_smbwt_%d.txt", prefix, i);
      sprintf(f->nxtpnt[i], "%s_nxtpnt_%d.txt", prefix, i);
      sprintf(f->nxtflg[i], "%s_nxtflg_%d.txt", prefix, i);
      sprintf(f->extidx[i], "%s_extidx_%d.txt", prefix, i);
    }


    return;

}

void getFileName(char *in_fname, FileName *f){
	//文件名初始化-----------------------------------------------------
printf("%s\n","setFileName(in_f,f); 243: Ok" );

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
   
  
        //strcpy (f->jmpmod,"data/jmpmod/jmpmod");
        strcpy (f->jmpmod,"jmpmod");
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
			//strcpy (*(f->seedidx+i),"data/seedidx/seedidx_");
			strcpy (*(f->seedidx+i),"seedidx_");
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
			//strcpy (*(f->comfile+i),"data/comfile/comfile_");
			strcpy (*(f->comfile+i),"comfile_");
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
			//strcpy (f->capidx[i],"data/capidx/capidx_");
			strcpy (f->capidx[i],"capidx_");
			sprintf(ch, "%d", i);
			strcat (f->capidx[i],ch);	
			strcat (f->capidx[i],".txt");	
		}

		//NUM_EXT个输出文件RelatFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (f->relat[i],"data/relat/relat_");
			strcpy (f->relat[i],"relat_");
			sprintf(ch, "%d", i);
			strcat (f->relat[i],ch);	
			strcat (f->relat[i],".txt");	
		}
//printf("f->idxfile = %s\n",f->relat[0]);			
	
		//NUM_EXT个输出文件SmBwtFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->smbwt+i),"data/smbwt/smbwt_");
			strcpy (*(f->smbwt+i),"smbwt_");
			sprintf(ch, "%d", i);
			strcat (*(f->smbwt+i),ch);	
			strcat (*(f->smbwt+i),".txt");	
		}

		//NUM_EXT个输出文件NxtPntFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->nxtpnt+i),"data/nxtpnt/nxtpnt_");
			strcpy (*(f->nxtpnt+i),"nxtpnt_");
			sprintf(ch, "%d", i);
			strcat (*(f->nxtpnt+i),ch);	
			strcat (*(f->nxtpnt+i),".txt");	
		}

		//NUM_EXT个输出文件NxtFlgFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->nxtflg+i),"data/nxtflg/nxtflg_");
			strcpy (*(f->nxtflg+i),"nxtflg_");
			sprintf(ch, "%d", i);
			strcat (*(f->nxtflg+i),ch);	
			strcat (*(f->nxtflg+i),".txt");	
		}

		//NUM_EXT个输出文件SmPosFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->extidx+i),"data/extidx/extidx_");
			strcpy (*(f->extidx+i),"extidx_");
			sprintf(ch, "%d", i);
			strcat (*(f->extidx+i),ch);	
			strcat (*(f->extidx+i),".txt");	
		}
		return;
    } 

printf("jmpmod out6 ok\n");
    perror("error：there are not [strcmp(in_fname, idxfile) == 0] ");
	return;

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



 
