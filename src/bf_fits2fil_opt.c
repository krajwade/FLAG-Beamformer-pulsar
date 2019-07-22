#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<glob.h>
#include<float.h>
#include<time.h>
#include"rw_header.h"
#include "cuda_utils.h"
#include "madfilter_small.h"
#include<omp.h>
#include<getopt.h>

/* Defining macros */
#define NUM_PSTOKES  4
#define NUM_CHANS  25
#define NUM_BEAMS  7
#define NUM_SAMPS_PER_BLOCK  100
#define TOT_CHANS  500
#define DATA_SIZE NUM_PSTOKES*NUM_CHANS*NUM_BEAMS*NUM_SAMPS_PER_BLOCK


/* C-program to convert bfpulsarfits to filterbank */
/* This version is more optimized and reduces cfitsio overheads by clever placing for nested for-loops */

/* Kaustubh Rajwade, JBCA, August 2018 */

/* requantization sub-routine */

void Float2Byte(float *pfbuf, int ilen,unsigned char *pcbuf){
 
 int i,j;
 float fMax = -FLT_MAX;
 float fMin = FLT_MAX;
 float frange;
 float fIntMax = (float) (powf(2.0 , 8) - 1.0);

 /* Get rid of the DC component */
 /* Extremely 'ad-hoc' way of doing things. Works for now!! */
 for (i = 0; i < ilen; ++i)
 {
     if (abs(pfbuf[i] - pfbuf[200]) > 100.0)
     {
         pfbuf[i] = pfbuf[i - 1];	
     }
 }

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


//#pragma omp parallel num_threads(6)
//#pragma omp for
 for (j = 0; j < ilen; ++j){

   pcbuf[j] = (unsigned char) roundf(((pfbuf[j]-fMin)/frange)*fIntMax);
 }
 return;
}

void print_usage() 
{
    printf("Usage: bf_fits2fil  -b <filebase>  -p <project id> <Optional flags>\n");
    printf("Options are:\n");
    printf("-h            prints the help\n");
    printf("-b            complete path to the basefilename\n");
    printf("-p            project ID\n");
    printf("-q [0 or 1]   toggles 8-bit quantization of the data\n");
    printf("-f [0 or 1]   toggles bandpass flip\n");
    printf("-m [0 or 1]   toggle MAD-based filtering of the bandpass\n");
}

/* Main code */
int main(int argc, char *argv[])
{

 int iNextOpt = 0;
 int numFiles = 20,e,i, j, k ,m, n,o,status=0,step_size;
 long nrows;
 char filename[128];
 FILE *f1=NULL,*f2=NULL,*f3=NULL,*f4=NULL,*f5=NULL,*f6=NULL,*f7=NULL;
 fitsfile *fptr,*fptr1[20];
 float *beam0, *beam1, *beam2, *beam3, *beam4, *beam5, *beam6;
 glob_t globbuf;
 int colnum=4;
 char timestamp[1024];
 char basetimestamp[1024];
 char proj_id[128];
 char bank[20];
 char alphabet[1];
 long mcnt, mcnt_int[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
 float *chan=NULL;
 unsigned char *beam0_8bit, *beam1_8bit, *beam2_8bit, *beam3_8bit, *beam4_8bit, *beam5_8bit, *beam6_8bit;
 int quant_flag=0, flip_flag=0, filt_flag=0;
 
/* Removing Previous files */

remove("BF_pulsar_0.fil");
remove("BF_pulsar_1.fil");
remove("BF_pulsar_2.fil");
remove("BF_pulsar_3.fil");
remove("BF_pulsar_4.fil");
remove("BF_pulsar_5.fil");
remove("BF_pulsar_6.fil");

 
/* Options */
  while ((iNextOpt = getopt(argc,argv, "hq:f:m:b:p:")) != -1)
  { 
      switch (iNextOpt)
        {


          case 'h':
          print_usage();
          return(0);
          break;

          case 'q':
          quant_flag = atoi(optarg);
          break;

          case 'f':
          flip_flag = atoi(optarg);
          break;

          case 'm':
          filt_flag = atoi(optarg);
          break;
  	
  	  case 'b':
	  sprintf(timestamp,"%s",optarg); 

	  case 'p':
	  sprintf(proj_id,"%s",optarg);

        }
  
   }

  if (quant_flag != 0 && quant_flag != 1)
       {
          fprintf(stderr,"Incorrect value for the quantization flag. Can be either 0 or 1.\n");
	  print_usage();
	  exit(0);	
       }

     if (flip_flag != 0 && flip_flag != 1)
       {
	  fprintf(stderr,"Incorrect value for the channel flip flag. Can be either 0 or 1.\n");
	  print_usage();
	  exit(0);	 
       }
 
     if (filt_flag !=0 && filt_flag != 1 )
     {
	  fprintf(stderr,"Incorrect value for filtering flag. Can be either 0 or 1.\n");
	  print_usage();
	  exit(0);
     }

     if (strcmp(timestamp,"") == 0)
     {
         printf("No time stamp given!\n");
	 print_usage();
         exit(0);
     } 

     if (strcmp(proj_id,"") == 0)
     {
         printf("No project id given!\n");
	 print_usage();
         exit(0);
     }
 
/* parsing all the fits files */
 
 strcpy(basetimestamp,timestamp);
 strcat(basetimestamp,"*.fits*");
 glob(basetimestamp,GLOB_ERR,NULL,&globbuf);

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
    sprintf(alphabet,"%c",bank[e]);
    strcpy(basetimestamp, timestamp);
    strcat(basetimestamp,alphabet);
    strcat(basetimestamp,".fits");
    fits_open_file(&fptr1[e], basetimestamp, READONLY, &status);

/* Check the size */ 
    if (status)
      flag[e] = 1;
    else
      flag[e] = 0;

    status=0;
 }

 //nrows=1000;
 beam0 =(float*)malloc(sizeof(float)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam1 =(float*)malloc(sizeof(float)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam2 =(float*)malloc(sizeof(float)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam3 =(float*)malloc(sizeof(float)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam4 =(float*)malloc(sizeof(float)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam5 =(float*)malloc(sizeof(float)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam6 =(float*)malloc(sizeof(float)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 chan = (float*)malloc(sizeof(float)*DATA_SIZE*numFiles);

 beam0_8bit =(unsigned char*)malloc(sizeof(unsigned char)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam1_8bit =(unsigned char*)malloc(sizeof(unsigned char)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam2_8bit =(unsigned char*)malloc(sizeof(unsigned char)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam3_8bit =(unsigned char*)malloc(sizeof(unsigned char)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam4_8bit =(unsigned char*)malloc(sizeof(unsigned char)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam5_8bit =(unsigned char*)malloc(sizeof(unsigned char)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 beam6_8bit =(unsigned char*)malloc(sizeof(unsigned char)*TOT_CHANS*NUM_SAMPS_PER_BLOCK);
 

/* open filterbank */

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

/* Writing the header */

 strcpy(filename,globbuf.gl_pathv[0]);

 fits_open_file(&fptr, filename, READONLY, &status);

 rw_header(fptr,f1,proj_id, timestamp,quant_flag, flip_flag);
 rw_header(fptr,f2,proj_id, timestamp, quant_flag, flip_flag);
 rw_header(fptr,f3,proj_id, timestamp, quant_flag, flip_flag);
 rw_header(fptr,f4,proj_id, timestamp,quant_flag, flip_flag);
 rw_header(fptr,f5,proj_id, timestamp,quant_flag, flip_flag);
 rw_header(fptr,f6,proj_id, timestamp, quant_flag, flip_flag);
 rw_header(fptr,f7,proj_id, timestamp,quant_flag, flip_flag);


/* Converting the spliced files to filterbank file for each beam */

 clock_t begin = clock(); 
 step_size=0;
 for(k=1; k<nrows-10;k++)
 {
    for(j=0;j<numFiles;j++)
    {
	    if (flag[j] == 0)
	    {
		    fits_movabs_hdu(fptr1[j], 2, NULL, &status);
		    fits_report_error(stderr,status);
		    fits_read_col(fptr1[j], TFLOAT,colnum, k, 1, DATA_SIZE ,NULL, chan + (j*DATA_SIZE),NULL, &status);
		    fits_read_col(fptr1[j], TULONG,2, k, 1, 1 ,NULL, &mcnt ,NULL, &status);

	    /* Logic to ignore smaller fits files and have internal counters to fill missing data */
		    if (status)
		    {
			    fits_report_error(stderr,status);
			    memset(chan+(j*DATA_SIZE),0,sizeof(float)*DATA_SIZE);
		    }
		    else if (mcnt - mcnt_int[j] != 0)
		    {
			    memset(chan+(j*DATA_SIZE),0,sizeof(float)*DATA_SIZE);
			    mcnt_int[j] += 200;
		    }
		    else
		    {
			    mcnt_int[j] += 200;
		    }
	    }
	    else
	    {
		    memset(chan+(j*DATA_SIZE),0,sizeof(float)*DATA_SIZE);  
	    }

    }
    step_size=0;
    for (i=0; i<NUM_SAMPS_PER_BLOCK;i++)
    {
	    /* Converting the chunks to contiguous frequency */

	    for (m = 0;m < 5; m++ )
	    {
                for (o=0; o < numFiles; o++ )
		{
		    for (n=0; n < 5 ;n++)
		    {

			    beam0[step_size] = chan[28*n + DATA_SIZE*o + 140*m + 700*i ] + chan[28*n +DATA_SIZE*o + 140*m + 700*i + NUM_BEAMS];
			    beam1[step_size] = chan[28*n + DATA_SIZE*o + 140*m + 700*i + 1] + chan[28*n +DATA_SIZE*o + 140*m + 700*i + NUM_BEAMS + 1];
			    beam2[step_size] = chan[28*n + DATA_SIZE*o + 140*m + 700*i + 2] + chan[28*n +DATA_SIZE*o + 140*m + 700*i + NUM_BEAMS + 2];
			    beam3[step_size] = chan[28*n + DATA_SIZE*o + 140*m + 700*i + 3] + chan[28*n +DATA_SIZE*o + 140*m + 700*i + NUM_BEAMS + 3];
			    beam4[step_size] = chan[28*n + DATA_SIZE*o + 140*m + 700*i + 4] + chan[28*n +DATA_SIZE*o + 140*m + 700*i + NUM_BEAMS + 4];
			    beam5[step_size] = chan[28*n + DATA_SIZE*o + 140*m + 700*i + 5] + chan[28*n +DATA_SIZE*o + 140*m + 700*i + NUM_BEAMS + 5];
			    beam6[step_size] = chan[28*n + DATA_SIZE*o + 140*m + 700*i + 6] + chan[28*n +DATA_SIZE*o + 140*m + 700*i + NUM_BEAMS + 6];
                            step_size++;
		    }
                
                }
	    }

      }
      if (flip_flag == 1)
	    {
		channel_flip_float(beam0, TOT_CHANS*NUM_SAMPS_PER_BLOCK , TOT_CHANS);
		channel_flip_float(beam1, TOT_CHANS*NUM_SAMPS_PER_BLOCK , TOT_CHANS);
		channel_flip_float(beam2, TOT_CHANS*NUM_SAMPS_PER_BLOCK , TOT_CHANS);
		channel_flip_float(beam3, TOT_CHANS*NUM_SAMPS_PER_BLOCK , TOT_CHANS);
		channel_flip_float(beam4, TOT_CHANS*NUM_SAMPS_PER_BLOCK , TOT_CHANS);
		channel_flip_float(beam5, TOT_CHANS*NUM_SAMPS_PER_BLOCK , TOT_CHANS);
		channel_flip_float(beam6, TOT_CHANS*NUM_SAMPS_PER_BLOCK , TOT_CHANS);
	    }

      /* Running MAD filter */
      if (quant_flag == 1)
      {
	   Float2Byte(beam0,TOT_CHANS*NUM_SAMPS_PER_BLOCK, beam0_8bit);
	   Float2Byte(beam1,TOT_CHANS*NUM_SAMPS_PER_BLOCK, beam1_8bit);
	   Float2Byte(beam2,TOT_CHANS*NUM_SAMPS_PER_BLOCK, beam2_8bit);
	   Float2Byte(beam3,TOT_CHANS*NUM_SAMPS_PER_BLOCK, beam3_8bit);
	   Float2Byte(beam4,TOT_CHANS*NUM_SAMPS_PER_BLOCK, beam4_8bit);
	   Float2Byte(beam5,TOT_CHANS*NUM_SAMPS_PER_BLOCK, beam5_8bit);
	   Float2Byte(beam6,TOT_CHANS*NUM_SAMPS_PER_BLOCK, beam6_8bit);
      }
      if (filt_flag == 1 && quant_flag != 1)
      {
	  fprintf(stderr, "Error: MAD filter only runs on 8-bit data\n");
	  exit(0);
      }
      if (filt_flag == 1 && quant_flag == 1)
      {
          run_madfilter(beam0_8bit, TOT_CHANS*NUM_SAMPS_PER_BLOCK, TOT_CHANS);
          run_madfilter(beam1_8bit, TOT_CHANS*NUM_SAMPS_PER_BLOCK, TOT_CHANS);
          run_madfilter(beam2_8bit, TOT_CHANS*NUM_SAMPS_PER_BLOCK, TOT_CHANS);
          run_madfilter(beam3_8bit, TOT_CHANS*NUM_SAMPS_PER_BLOCK, TOT_CHANS);
          run_madfilter(beam4_8bit, TOT_CHANS*NUM_SAMPS_PER_BLOCK, TOT_CHANS);
          run_madfilter(beam5_8bit, TOT_CHANS*NUM_SAMPS_PER_BLOCK, TOT_CHANS);
          run_madfilter(beam6_8bit, TOT_CHANS*NUM_SAMPS_PER_BLOCK, TOT_CHANS);
      }
      
      /* Writing data to each filterbank file */

      if (quant_flag != 1)
      {
	      fwrite(beam0,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(float),f1);
	      fwrite(beam1,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(float),f2);
	      fwrite(beam2,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(float),f3);
	      fwrite(beam3,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(float),f4);
	      fwrite(beam4,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(float),f5);
	      fwrite(beam5,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(float),f6);
	      fwrite(beam6,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(float),f7);
      }

      else
      {
	      fwrite(beam0_8bit,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(unsigned char),f1);
	      fwrite(beam1_8bit,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(unsigned char),f2);
	      fwrite(beam2_8bit,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(unsigned char),f3);
	      fwrite(beam3_8bit,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(unsigned char),f4);
	      fwrite(beam4_8bit,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(unsigned char),f5);
	      fwrite(beam5_8bit,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(unsigned char),f6);
	      fwrite(beam6_8bit,NUM_SAMPS_PER_BLOCK*TOT_CHANS,sizeof(unsigned char),f7);

      } 
}

 
/* converting the filterbank to contiguous frequency channels */

  clock_t end = clock();

 

  double time_spent = (double) ((end - begin)/ CLOCKS_PER_SEC) ;
  printf("Total time:%lf seconds \n",time_spent);


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
  free(beam0_8bit);
  free(beam1_8bit);
  free(beam2_8bit);
  free(beam3_8bit);
  free(beam4_8bit);
  free(beam5_8bit);
  free(beam6_8bit);
  free(chan);
 
  return(0);

}
