#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"/usr/include/cfitsio/fitsio.h"
#include<glob.h>
#include<time.h>
#include<float.h>
#include"rw_header.h"
#include<omp.h>
#define NUM_PSTOKES  4
#define NUM_CHANS  25
#define NUM_BEAMS  7
#define NUM_SAMPS_PER_BLOCK  100
#define TOT_CHANS  500
#define DATA_SIZE NUM_PSTOKES*NUM_CHANS*NUM_BEAMS*NUM_SAMPS_PER_BLOCK


/* C-program to convert bfpulsarfits to filterbank */

/* Kaustubh Rajwade, WVU, 03 November 2015 */

/* requantization sub-routine */

/*
void Float2Byte(float *pfbuf, int ilen,float *pcbuf){
 
 int i = 0;
 float fMax = -FLT_MAX;
 float fMin = FLT_MAX;
 float frange;
 float fIntMax = (float) (powf(2.0 , 8) - 1.0);

 for (i = 0; i < ilen; ++i){
  if (pfbuf[i] > fMax)
           fMax = pfbuf[i];
    if (pfbuf[i] < fMin)
           fMin = pfbuf[i];

 }

 frange = fMax - fMin;

// Change number of threads here
// Devansh Agarwal, 05 July 2017
// Parallelized requantization


#pragma omp parallel num_threads(6)
{
#pragma omp for
 for (i = 0; i < ilen; ++i){

   pcbuf[i] = (float) roundf(((pfbuf[i]-fMin)/frange)*fIntMax);
 }
}
 return;
}
*/

/* Main code */
int main(int argc, char *argv[])
{

 int numFiles = 20,e,i, j, k,l, m, n,o,status=0,step_size,ind;
 long nrows;
 char filename[128];
 FILE *f1=NULL,*f2=NULL,*f3=NULL,*f4=NULL,*f5=NULL,*f6=NULL,*f7=NULL;
 fitsfile *fptr,*fptr1[20];
 float XX[DATA_SIZE/4],YY[DATA_SIZE/4],XY[DATA_SIZE/4],YX[DATA_SIZE/4], StokesI[DATA_SIZE/4];
 float *beam0, *beam1, *beam2, *beam3, *beam4, *beam5, *beam6;
 glob_t globbuf;
 int colnum=4;
 char timestamp[1024];
 char bank[20];
 char alphabet[1];
 float *chan=NULL;

/* Removing Previous files */

 char remcmd[128];

 strcpy(remcmd,"rm -rf BF_Pulsar_0.fil");
 system(remcmd);
 strcpy(remcmd,"rm -rf BF_Pulsar_1.fil");
 system(remcmd);
 strcpy(remcmd,"rm -rf BF_Pulsar_2.fil");
 system(remcmd);
 strcpy(remcmd,"rm -rf BF_Pulsar_3.fil");
 system(remcmd);
 strcpy(remcmd,"rm -rf BF_Pulsar_4.fil");
 system(remcmd);
 strcpy(remcmd,"rm -rf BF_Pulsar_5.fil");
 system(remcmd);
 strcpy(remcmd,"rm -rf BF_Pulsar_6.fil");
 system(remcmd);

  
/* parsing all the fits files */

 sprintf(timestamp,"%s",argv[1]);
 strcat(timestamp,"*.fits*");
 glob(timestamp,GLOB_ERR,NULL,&globbuf);

 if(globbuf.gl_pathv == NULL)
 {
      fprintf(stderr,"Error: Fits files not found. The code should be run from the same directory as the fits files.\n");
      exit(1);
 }
 strcpy(filename,globbuf.gl_pathv[0]);

/* defining arrays */

 fits_open_file(&fptr, filename, READONLY, &status);
 fits_movabs_hdu(fptr, 2, NULL, &status);
 fits_get_num_rows(fptr, &nrows, &status);
 fits_close_file(fptr,&status);
 printf("nrows=%ld\n",nrows);
 strcpy(bank,"ABCDEFGHIJKLMNOPQRST");

/* Opening all files together*/
 int flag[20];
 for (e=0; e<numFiles;e++)
 {
    sprintf(timestamp,"%s",argv[1]);
    sprintf(alphabet,"%c",bank[e]);
    strcat(timestamp,alphabet);
    strcat(timestamp,".fits");
    fits_open_file(&fptr1[e], timestamp, READONLY, &status);
  
    if (status)
      flag[e] = 1;
    else
      flag[e] = 0;
    status=0;
 }

 nrows=500;
 beam0 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
 beam1 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
 beam2 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
 beam3 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
 beam4 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
 beam5 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
 beam6 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
 chan = (float*)malloc(sizeof(float)*DATA_SIZE);


/* Converting the spliced files to filterbank file for each beam */
 
 step_size=0;
 clock_t begin = clock();
 for(k=1; k<nrows;k++)
 {
    //printf("nrow=%d\n",k);
    for (i=0; i<NUM_SAMPS_PER_BLOCK;i++)
    {
       for (o=0;o<5;o++)
       {

          for(j=0;j<numFiles;j++)
          {

             if (flag[j] == 0)
             {
                 fits_movabs_hdu(fptr1[j], 2, NULL, &status);
                 fits_report_error(stderr,status);
                 fits_read_col(fptr1[j], TFLOAT,colnum, k, 1, DATA_SIZE ,NULL, chan,NULL, &status);
                 fits_report_error(stderr,status);
             }
            else
            {
                 memset(chan,0,sizeof(float)*DATA_SIZE);  
            }

/* Converting the chunks to contiguous frequency */
            int index=7,count=0,idx;
            for (m = 0;m < 5; m++ )
            {
            
               idx = 28*m + 140*o + 700*i;
               for (n=0;n<NUM_BEAMS;n++)
               {
           
                   XX[count] = chan[idx];
                   YY[count] = chan[idx+index];
                   XY[count] = chan[idx + 2*index];
                   YX[count] = chan[idx + 3*index];

                   StokesI[count] = XX[count] + YY[count];
           
                   count++;
                   idx++;
               }

            }


/* Reading data from each fits file to an array */
           for(l=0;l<5;l++)
           {
 
               ind = NUM_BEAMS*l;
               beam0[step_size] = StokesI[ind];
               beam1[step_size] = StokesI[ind+1];
               beam2[step_size] = StokesI[ind+2];
               beam3[step_size] = StokesI[ind+3];
               beam4[step_size] = StokesI[ind+4];
               beam5[step_size] = StokesI[ind+5]; 
               beam6[step_size] = StokesI[ind+6];
               step_size++ ;   
           
           }
 
         }
 
       }

    }

 }
 
 clock_t end = clock();
/* converting the filterbank to contiguous frequency channels */

 
  f1 = fopen("BF_pulsar_0.fil","ab");
  if (f1 == NULL) 
  {
      printf("Error opening File\n");
      exit(1);
  }   
  f2 = fopen("BF_pulsar_1.fil","ab");
  if (f2 == NULL) 
  {
      printf("Error opening File\n");
      exit(1);
  }
  f3 = fopen("BF_pulsar_2.fil","ab");
  if (f3 == NULL) 
  {

      printf("Error opening File\n");
      exit(1);
  }
  f4 = fopen("BF_pulsar_3.fil","ab");
  if (f4 == NULL) 
  {

      printf("Error opening File\n");
      exit(1);
  }
  f5 = fopen("BF_pulsar_4.fil","ab");
  if (f5 == NULL) 
  {

      printf("Error opening File\n");
      exit(1);
  }
  f6 = fopen("BF_pulsar_5.fil","ab");
  if (f6 == NULL) 
  {

      printf("Error opening File\n");
      exit(1);
  }
  f7 = fopen("BF_pulsar_6.fil","ab");
  if (f7 == NULL) 
  {
     printf("Error opening File\n");
     exit(1);
  }

  double time_spent = (double)(end - begin)/ CLOCKS_PER_SEC;
 
  printf("Total time taken:%lf seconds\n",time_spent);

/* Writing the header */

  strcpy(filename,globbuf.gl_pathv[0]);

  fits_open_file(&fptr, filename, READONLY, &status);

  rw_header(fptr,f1);
  rw_header(fptr,f2);
  rw_header(fptr,f3);
  rw_header(fptr,f4);
  rw_header(fptr,f5);
  rw_header(fptr,f6);
  rw_header(fptr,f7);

 fits_close_file(fptr,&status);

/* Writing data to each filterbank file */

  fwrite(beam0,NUM_SAMPS_PER_BLOCK*nrows*TOT_CHANS,sizeof(float),f1);
  fwrite(beam1,NUM_SAMPS_PER_BLOCK*nrows*TOT_CHANS,sizeof(float),f2);
  fwrite(beam2,NUM_SAMPS_PER_BLOCK*nrows*TOT_CHANS,sizeof(float),f3);
  fwrite(beam3,NUM_SAMPS_PER_BLOCK*nrows*TOT_CHANS,sizeof(float),f4);
  fwrite(beam4,NUM_SAMPS_PER_BLOCK*nrows*TOT_CHANS,sizeof(float),f5);
  fwrite(beam5,NUM_SAMPS_PER_BLOCK*nrows*TOT_CHANS,sizeof(float),f6);
  fwrite(beam6,NUM_SAMPS_PER_BLOCK*nrows*TOT_CHANS,sizeof(float),f7);
  


 /* Closing the files */

  for (j = 0; j<numFiles;j++)
  {

     fits_close_file(fptr1[j],&status);

  }  
 
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
  fclose(f5);
  fclose(f6);
  fclose(f7);

  free(beam0);
  free(beam1);
  free(beam2);
  free(beam3);
  free(beam4);
  free(beam5);
  free(beam6);
  free(chan);
 
  return(0);

}
