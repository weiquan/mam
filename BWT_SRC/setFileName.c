#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "setFileName.h"
//void setFileName(char *in_fname, FileName *out_fname,int flg ){
const char* ROOT=".";
void setFileName(char *in_fname, FileName *f,int flg ){
	//文件名初始化-----------------------------------------------------

	char jmpmod[LEN_FILE_NAME];
	char seedidx[LEN_FILE_NAME];
	char idxfile[LEN_FILE_NAME];
	char comfile[LEN_FILE_NAME];

	strcpy(jmpmod, "jmpmod" );
	strcpy(seedidx,"seedidx");
	strcpy(idxfile,"idxfile");	
	strcpy(comfile,"comfile");
	
	char (*pf)[LEN_FILE_NAME];
	//FileName *f;
	
	int i;
	char ch[5];
    /*	
	if(0==flg){
		if(NULL == (f = (FileName*)malloc(sizeof(FileName)))){
		    perror("error:[malloc(NUM_EXT*sizeof(FileName))]");
		    exit(1);
		}
	}
	if(NULL == f){
	    perror("error:[NULL == f]");
	    exit(1);
	}*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (strcmp(in_fname, jmpmod) == 0) {
		char *s = (char*)malloc(LEN_FILE_NAME * sizeof(char));
        if(NULL == s){
		    perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
		    exit(1);
		}
		f->jmpmod  = s; 
		//strcpy (*(f->jmpmod),"./jmpmod/jmpmod");
		//strcat (*(f->jmpmod),".txt");	
		//strcpy (*(f->seedidx+i),"./seedidx/seedidx_");
        //itoa(i,ch,10);
        //strcat (*(f->seedidx+i),ch);	
        //strcat (*(f->seedidx+i),".txt");	
        //sprintf(f->jmpmod, "%s/jmpmod/jmpmod.txt", ROOT);
        sprintf(f->jmpmod, "%s/jmpmod.txt", ROOT);
        //printf("%s\n", f->seedidx[i]);
        //out_fname = f;
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
			//strcpy (*(f->seedidx+i),"./seedidx/seedidx_");
			//itoa(i,ch,10);
			//strcat (*(f->seedidx+i),ch);	
			//strcat (*(f->seedidx+i),".txt");	
		    //sprintf(f->seedidx[i], "%s/seedidx/seedidx_%u.txt", ROOT, i);
		    sprintf(f->seedidx[i], "%s/seedidx_%u.txt", ROOT, i);
            //printf("%s\n", f->seedidx[i]);
        }
		//out_fname = f;
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
			//strcpy (*(f->comfile+i),"./seedidx/seedidx_");
			//itoa(i,ch,10);
			//strcat (*(f->comfile+i),ch);	
			//strcat (*(f->comfile+i),".txt");	
		    //sprintf(f->comfile[i], "%s/comfile/comfile_%u.txt", ROOT, i);
		    sprintf(f->comfile[i], "%s/comfile_%u.txt", ROOT, i);
        }
		//out_fname = f;
		return;
    } 

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (strcmp(in_fname, idxfile) == 0) {
		if(NULL == (pf = (char(*)[LEN_FILE_NAME])malloc(NUM_FILES * NUM_EXT * LEN_FILE_NAME * sizeof(char)))){
		    perror("error：[malloc(NUM_FILE_NAME*NUM_EXT*LEN_FILE_NAME*sizeof(char))]");
		    exit(1);
		}
		//int spc = (NUM_EXT*LEN_FILE_NAME);
		int spc = NUM_EXT;
        f->capidx  = pf + 0*spc; 
		f->relat   = pf + 1*spc; 
		f->smbwt   = pf + 2*spc; 
		f->nxtpnt  = pf + 3*spc; 
		f->nxtflg  = pf + 4*spc; 
		f->smpos   = pf + 5*spc; 

		//NUM_EXT个输出文件CapIdxFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->capidx+i),"./capidx/capidx_");
			//itoa(i,ch,10);
			//strcat (*(f->capidx+i),ch);	
			//strcat (*(f->capidx+i),".txt");	
		    //sprintf(f->capidx[i], "%s/capidx/capidx_%u.txt", ROOT, i);
		    sprintf(f->capidx[i], "%s/capidx_%u.txt", ROOT, i);
        }

		//NUM_EXT个输出文件RelatFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->relat+i),"./relat/relat_");
			//itoa(i,ch,10);
			//strcat (*(f->relat+i),ch);	
			//strcat (*(f->relat+i),".txt");	
		    //sprintf(f->relat[i], "%s/relat/relat_%u.txt", ROOT, i);
		    sprintf(f->relat[i], "%s/relat_%u.txt", ROOT, i);
        }
		//NUM_EXT个输出文件SmBwtFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->smbwt+i),"./smbwt/smbwt_");
			//itoa(i,ch,10);
			//strcat (*(f->smbwt+i),ch);	
			//strcat (*(f->smbwt+i),".txt");	
		    //sprintf(f->smbwt[i], "%s/smbwt/smbwt_%u.txt", ROOT, i);
		    sprintf(f->smbwt[i], "%s/smbwt_%u.txt", ROOT, i);
        }

		//NUM_EXT个输出文件NxtPntFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->nxtpnt+i),"./nxtpnt/nxtpnt_");
			//itoa(i,ch,10);
			//strcat (*(f->nxtpnt+i),ch);	
			//strcat (*(f->nxtpnt+i),".txt");	
		    //sprintf(f->nxtpnt[i], "%s/nxtpnt/nxtpnt_%u.txt", ROOT, i);
		    sprintf(f->nxtpnt[i], "%s/nxtpnt_%u.txt", ROOT, i);
        }

		//NUM_EXT个输出文件NxtFlgFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->nxtflg+i),"./nxtflg/nxtflg_");
			//itoa(i,ch,10);
			//strcat (*(f->nxtflg+i),ch);	
			//strcat (*(f->nxtflg+i),".txt");	
            //sprintf(f->nxtflg[i], "%s/nxtflg/nxtflg_%u.txt", ROOT, i);
            sprintf(f->nxtflg[i], "%s/nxtflg_%u.txt", ROOT, i);
        }

		//NUM_EXT个输出文件SmPosFile_i.txt
		for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->smpos+i),"./smpos/smpos_");
			//itoa(i,ch,10);
			//strcat (*(f->smpos+i),ch);	
			//strcat (*(f->smpos+i),".txt");	
		    //sprintf(f->smpos[i], "%s/smpos/smpos_%u.txt", ROOT, i);
		    sprintf(f->smpos[i], "%s/smpos_%u.txt", ROOT, i);
        }
		//NUM_EXT个输入文件ExtIdxFile_i.txt
        /*
        for(i=0;i<NUM_EXT; i++){
			//strcpy (*(f->seedidx+i),"./extidx/extidx_");
			//itoa(i,ch,10);
			//strcat (*(f->seedidx+i),ch);	
			//strcat (*(f->seedidx+i),".txt");	
		    sprintf(f->seedidx[i], "./extidx/extidx_%u.txt", i);
        }*/
	
		//out_fname = f;
		return;
    } 
    
    perror("error：there are not [strcmp(in_fname, idxfile) == 0] ");
	return;
}


 
