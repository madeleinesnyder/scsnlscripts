/*
  (c) Copyright Jonas Larsson 2008

    This file is part of mlpipe: a pipe I/O library for Matlab and Octave.

    mlpipe is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mlpipe is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with mlpipe.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <mex.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#define MAXPIPES 1000

static FILE* pSystemFIDs[MAXPIPES];     /* Global system file IDs */
static int   pPositions[MAXPIPES];      /* Current pipe positions */
static char  pPermissions[MAXPIPES][2]; /* Pipe permissions ('r' or 'w') */
static int   pByteSwap[MAXPIPES];       /* Whether to byte swap or not; 0: no, 1: yes */
static char* pFilenames[MAXPIPES];      /* Filenames associated with each pipe (including path) */
static char* pCommands[MAXPIPES];       /* Commands used to open each pipe */
static int   nPipes;

void validatePipe( pid ) {
  if (pid<0 || pid>=MAXPIPES)
    mexErrMsgTxt("(mlpipe) Invalid pipe ID\n");
  if (pSystemFIDs[pid]==0) 
    mexErrMsgTxt("(pclose) Pipe not opened\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  static int firsttime=1;
  int i,j,pid;
  if (firsttime) {
    /* initialize all static variables to zero  */
    nPipes=0;
    memset(pSystemFIDs,0,MAXPIPES);
    memset(pPositions,0,MAXPIPES);
    memset(pPermissions,0,MAXPIPES);
    memset(pByteSwap,0,MAXPIPES);
    memset(pFilenames,0,MAXPIPES);
    memset(pCommands,0,MAXPIPES);
    firsttime=0;
  }

  /* Parse inputs. Always at least 2. */
  if (nrhs<2)
    mexErrMsgTxt("mlpipe(): not enough input arguments.\n");
 
  /* First input is numeric (pid) */
  pid=(int) *mxGetPr( prhs[0] );
  /* Pipe IDs are negative to distinguish them from file IDs, so change sign */
  pid=-pid;

  /* Second input is command to run */
  int buflen = mxGetN(prhs[1])*mxGetM( prhs[1] )+1;
  char *command= (char*)malloc(buflen);
  mxGetString( prhs[1], command, buflen);
  if (command == NULL) {
    mexErrMsgTxt("(mlpipe) Invalid command!\n");
  }

  if (!strcmp(command,"open")) {
    /*
      popen():
      [pid,message]=mlpipe(0,'open',filename,permission,byteswap,pipecommand) : open pipe with command
      [filename,permission,byteswap,pipecommand]=mlpipe(pid,'open')           : return pipe status
      pid=mlpipe(0,'open',[])                                                 : return vector with open pipes
    */

    if (nrhs>2) {
      /* open file */
      /* 1. Find first available empty pipe ID */
      pid=-1;
      for (i=0; i<MAXPIPES; i++) {
	if (pSystemFIDs[i]==0) {
	  pid=i;
	  nPipes++;
	  break;
	}
      }
      
      if (pid<0)
	mexErrMsgTxt("(mlpipe) Maximum number of open pipes reached; close one first.\n");
      
      /* Determine pipe command: concatenate command and filename */
      buflen = mxGetN(prhs[2])*mxGetM( prhs[2] )+1;
      char* filename=(char*) malloc(buflen);
      mxGetString( prhs[2], filename, buflen);
      
      char permission[2];
      mxGetString( prhs[3], permission, 2);

      int byteswap =(int) *mxGetPr( prhs[4] );

      buflen = mxGetN(prhs[5])*mxGetM( prhs[5] )+1;
      char* pipecommand=(char*) malloc(buflen+strlen(filename));
      mxGetString( prhs[5], pipecommand, buflen);
      
      strcat(pipecommand, filename);
      
      /* Open pipe */
      pSystemFIDs[pid]=popen( pipecommand, permission);

      /* Check for errors */
      if (pSystemFIDs[pid]==NULL)
	mexErrMsgTxt("(mlpipe) Error opening pipe\n");
      
      /* Save parameters */
      strcpy(pPermissions[pid],permission);
      pPositions[pid]=0;
      pByteSwap[pid]=byteswap; 
      pFilenames[pid]=filename;
      pCommands[pid]=pipecommand;
  
      /* Return pipe ID (pid+1) - note negative sign */
      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(plhs[0]) = (double) -(pid+1);
      
      
    } else {
      if (pid==0) {
	/* return open pipe IDs */
	plhs[0] = mxCreateDoubleMatrix(1,nPipes,mxREAL);
	double* dp=mxGetPr(plhs[0]);
	for (i=0; i<MAXPIPES; i++) {
	  if (pSystemFIDs[i]>0) {
	    *dp=(double) -(i+1);
	    dp++;
	  }
	}
	
      } else {
	/* return info about specified pipe ID */
	pid--;
	validatePipe( pid );
	if (nlhs<4) 
	  mexErrMsgTxt("(mlpipe) Invalid syntax|n");
	/* Everything OK, create output arrays */
	plhs[0] = mxCreateString( pFilenames[pid] );
	plhs[1] = mxCreateString( pPermissions[pid] );
	plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
	*mxGetPr(plhs[2]) = (double) pByteSwap[pid];
	plhs[3] = mxCreateString( pCommands[pid] );
      }
      
    }
    
  } else if (!strcmp(command,"close")) {
    /*
      pclose():
      mlpipe(pid,'close')    : close pipe pid
      mlpipe(0,'close')      : close all pipes
    */
    if (pid==0) {
      /* close all pipes */
      for (i=0; i<MAXPIPES; i++) {
	if (pSystemFIDs[i]>0) {
	  pclose(pSystemFIDs[i]);
	  free(pCommands[i]);
	  free(pFilenames[i]);
	}
      }
      nPipes=0;
      memset(pSystemFIDs,0,MAXPIPES);
      memset(pPositions,0,MAXPIPES);
      memset(pPermissions,0,MAXPIPES*2);
      memset(pByteSwap,0,MAXPIPES);
      memset(pFilenames,0,MAXPIPES);
      memset(pCommands,0,MAXPIPES);


    } else {
      /* close this pipe */
      pid--; /* convert to zero-offset */
      validatePipe( pid );
      memset(pPermissions[pid],0,2);
      pPositions[pid]=0;
      pByteSwap[pid]=0;
      free(pFilenames[pid]);
      pFilenames[pid]=NULL;
      free(pCommands[pid]);
      pCommands[pid]=NULL;
      pclose(pSystemFIDs[pid]);
      pSystemFIDs[pid]=0;
      nPipes--;
    }

  } else if (!strcmp(command,"read")) {
    /*
      pread():
      [data,count]=mlpipe(pid,'read',size,nbytes,byteswap)
     */
    pid--;
    validatePipe( pid );
    int rows=1,cols=1;
    rows = (int) *mxGetPr( prhs[2]);
    if (mxGetM(prhs[2])>1)
      cols = (int) *(mxGetPr( prhs[2] )+1);
    int nbytes = (int) *mxGetPr( prhs[3]);

    buflen = mxGetN(prhs[4])*mxGetM( prhs[4] )+1;
    char* datatype=(char*) malloc(buflen);
    mxGetString( prhs[4], datatype, buflen);

    int byteswap = (int) *mxGetPr( prhs[5]);
    if (byteswap<0)
      byteswap=pByteSwap[pid];

    char* buf;
    int count;
    if (rows<1) {
      /* read to eof */
      cols=1;
      int blocksize=100000;
      buf=(char*)malloc(blocksize*nbytes);
      count=fread( buf, nbytes, blocksize, pSystemFIDs[i] );
      if (count==0)
	mexErrMsgTxt("(pread) Error reading file\n");
      while (count==blocksize) {
	/* reallocate buffer and continue reading */
	buf=(char*)realloc(buf,(count+blocksize)*nbytes);
	count+=fread( buf, nbytes, blocksize, pSystemFIDs[i] );
      } 
      rows=count;

    } else {
      /* allocate buffer */
      buf=(char*)malloc(rows*cols*nbytes);
      
      /* read data */ 
      count=fread( buf, nbytes, rows*cols, pSystemFIDs[i] );
      if (count==0)
	mexErrMsgTxt("(pread) Error reading file\n");
      if (count<rows*cols)
	mexPrintf("(pread) Warning: could only read %i elements out of %i expected.\n",count,rows*cols);
    }

    /* Update file position */
    pPositions[pid]+=count*nbytes;

    /* byte swap */
    if (byteswap) {
      char* tmp = malloc(nbytes);
      char* to=buf;
      char* from;
      for (i=0; i<count; i++) {
	memcpy( tmp, to, nbytes );
	from = tmp + nbytes - 1;
	for (j=0; j<nbytes; j++, from--, to++)
	  *to=*from;
      }
      free(tmp);      
    }

    /* return result */
    plhs[0]=mxCreateDoubleMatrix(rows,cols,mxREAL);
    double* dp=mxGetPr( plhs[0] );

    if (!strcmp(datatype,"uchar")) {
      unsigned char* bp=(unsigned char*)buf;
      for (i=0; i<count; i++, dp++,bp++)
	*dp = (double) *bp;

    } else if (!strcmp(datatype,"char")) {
      char* bp=(char*)buf;
      for (i=0; i<count; i++, dp++,bp++)
	*dp = (double) *bp;

    } else if (!strcmp(datatype,"ushort")) {
      unsigned short* bp=(unsigned short*)buf;
      for (i=0; i<count; i++, dp++,bp++)
	*dp = (double) *bp;

    } else if (!strcmp(datatype,"short")) {
      short* bp=(short*)buf;
      for (i=0; i<count; i++, dp++,bp++)
	*dp = (double) *bp;

    } else if (!strcmp(datatype,"int")) {
      int* bp=(int*)buf;
      for (i=0; i<count; i++, dp++,bp++)
	*dp = (double) *bp;

    } else if (!strcmp(datatype,"float")) {
      float* bp=(float*)buf;
      for (i=0; i<count; i++, dp++,bp++)
	*dp = (double) *bp;

    } else if (!strcmp(datatype,"double")) {
      double* bp=(double*)buf;
      for (i=0; i<count; i++, dp++,bp++)
	*dp = (double) *bp;
    } else 
      mexErrMsgTxt("(mlpipe) Unknown datatype\n");

    plhs[1]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[1]) = (double) count;

    free(datatype);

  } else if (!strcmp(command,"write")) {
    /*
      pwrite():
      count=mlpipe(pid,'read',data,nbytes,datatype,byteswap)
     */
    pid--; /* convert to zero-offset */
    validatePipe(pid);

    /* parse arguments */
    double* data = mxGetPr( prhs[2]);
    /* get size of data */
    int datasize = mxGetN(prhs[2])*mxGetM( prhs[2] );

    int nbytes = (int) *mxGetPr( prhs[3]);

    buflen = mxGetN(prhs[4])*mxGetM( prhs[4] )+1;
    char* datatype=(char*) malloc(buflen);
    mxGetString( prhs[4], datatype, buflen);

    int byteswap = (int) *mxGetPr( prhs[5]);
    if (byteswap<0)
      byteswap=pByteSwap[pid];

    /* convert to byte buffer */
    char* buf=malloc(datasize*nbytes);
    double* dp=data;

    if (!strcmp(datatype,"uchar")) {
      unsigned char* bp=(unsigned char*)buf;
      for (i=0; i<datasize; i++, dp++,bp++)
	*bp = (unsigned char) *dp;

    } else if (!strcmp(datatype,"char")) {
      char* bp=(char*)buf;
      for (i=0; i<datasize; i++, dp++,bp++)
	*bp = (char) *dp;

    } else if (!strcmp(datatype,"ushort")) {
      unsigned short* bp=(unsigned short*)buf;
      for (i=0; i<datasize; i++, dp++,bp++)
	*bp = (unsigned short) *dp;

    } else if (!strcmp(datatype,"short")) {
      short* bp=(short*)buf;
      for (i=0; i<datasize; i++, dp++,bp++)
	*bp = (short) *dp;

    } else if (!strcmp(datatype,"int")) {
      int* bp=(int*)buf;
      for (i=0; i<datasize; i++, dp++,bp++)
	*bp = (int) *dp;

    } else if (!strcmp(datatype,"float")) {
      float* bp=(float*)buf;
      for (i=0; i<datasize; i++, dp++,bp++)
	*bp = (float) *dp;

    } else if (!strcmp(datatype,"double")) {
      double* bp=(double*)buf;
      for (i=0; i<datasize; i++, dp++,bp++)
	*bp = (double) *dp;

    } else 
       mexErrMsgTxt("(mlpipe) Unknown datatype\n");

    /* byte swap */
    if (byteswap) {
      char* tmp = malloc(nbytes);
      char* to=buf;
      char* from;
      for (i=0; i<datasize; i++) {
	memcpy( tmp, to, nbytes );
	from = tmp + nbytes - 1;
	for (j=0; j<nbytes; j++, from--, to++)
	  *to=*from;
      }
      free(tmp);      
    }

    /* write buffer to file */ 
    int count=fwrite( buf, nbytes, datasize, pSystemFIDs[pid] );
    free(buf);
    
    /* update file position */
    pPositions[pid]+=count*nbytes;

    /* return count */
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0]) = (double) count;

  } else if (!strcmp(command,"scanf")) {
    /*
      pscanf():
      [data,count]=mlpipe(pid,'scanf',format,size)
     */

    mexErrMsgTxt("(mlpipe) scanf not yet supported.\n");

    pid--; /* convert to zero-offset */
    validatePipe(pid);

    /* get format string */
    buflen = mxGetN(prhs[2])*mxGetM( prhs[2] )+1;
    char *format= (char*)malloc(buflen);
    mxGetString( prhs[2], format, buflen);

    /* get size argument */
    int rows=1,cols=1;
    rows = (int) *mxGetPr( prhs[3]);
    if (mxGetM(prhs[3])>1)
      cols = (int) *(mxGetPr( prhs[3] )+1);

    /* initialize output array */

    /* read from file */
    /*    int count=fscanf(pSystemFIDs[pid], format, &data ); */

    /* update file position */ 

    /* Return result */

    free(format);

  } else if (!strcmp(command,"printf")) {
    /* Print formatted data to a file */
    mexErrMsgTxt("(mlpipe) printf not yet supported.\n");
  } else if (!strcmp(command,"seek")) {
    /*
      pseek():
      status=mlpipe(pid,'seek',offset,origin)
     */
    pid--; /* convert to zero-offset */
    validatePipe(pid);

    int offset =(int) *mxGetPr( prhs[2] );
    int origin =(int) *mxGetPr( prhs[3] );

    int count=offset;

    if (origin<0) 
      count=offset-pPositions[pid];

    if (count<0 || origin>0) {
      /* close and reopen file */
      if (origin>0) {
	pclose(pSystemFIDs[pid]);
	pSystemFIDs[pid]=popen( pCommands[pid], "r" );
	/* read ahead to eof to determine file size */
	int blocksize=10000;
	char* buf=malloc(blocksize);
	int fs=0,c=0;
	do { 
	  c=fread( buf, 1, blocksize, pSystemFIDs[pid] );
	  fs+=c;
	} while (c>0);
	free(buf);
    
	/* compute count from bof */
	count=fs+offset; /* offset must be negative for offset=1 */

	/* close and reopen file */
	pclose(pSystemFIDs[pid]);
	pSystemFIDs[pid]=popen( pCommands[pid], "r" );
	pPositions[pid]=0;

      } else {
	/* close and reopen file */
	pclose(pSystemFIDs[pid]);
	pSystemFIDs[pid]=popen( pCommands[pid], "r" );
	count=pPositions[pid]+count;
	pPositions[pid]=0;
      }
    } 

    /* Read dummy data */
    char* buf=malloc(count);
    count=fread( buf, 1, count, pSystemFIDs[pid] );
    free(buf);
    
    /* Update file position */
    pPositions[pid]+=count;    

    plhs[0]= mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0]) = (double) pPositions[pid];

  } else if (!strcmp(command,"tell")) {
    /*
      ptell():
      pos=mlpipe(pid,'tell') : return position of pipe pid
     */
    pid--; /* convert to zero-offset */
    validatePipe(pid);
    plhs[0]= mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0]) = (double) pPositions[pid];

  } else {
    mexErrMsgTxt("(mlpipe) Invalid command\n");
  }

  free(command);


}
