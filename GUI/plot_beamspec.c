#include<stdio.h>
#include<stdlib.h>
#include<float.h>
#include<math.h>
#include<string.h>
#include"/usr/include/cfitsio/fitsio.h"
#include<glob.h>
#include<cpgplot.h>
#include<getopt.h>
#define TOT_CHANS 500
#define NUM_PSTOKES  4
#define NUM_CHANS  25
#define NUM_BEAMS  7
#define NUM_SAMPS_PER_BLOCK  10
#define DATA_SIZE NUM_PSTOKES*NUM_CHANS*NUM_BEAMS*NUM_SAMPS_PER_BLOCK


/* C-program to plot spectra of 7 beams */

/* Kaustubh Rajwade, WVU, 03 November 2015 */


void print_usage() {
    printf("Usage: plot_beamspec -f <filebase> -l <numlow> -g <numhigh>\n");
}



int main(int argc, char *argv[]){


 int iNextOpt = 0 ;
 int numFiles = 20,x,y, z,i, j, k,l, m, n,o,status=0,step_size,ind;
 long nrows;
 char filename[128];
 fitsfile *fptr,*fptr1[20];
 float XX[DATA_SIZE/4],YY[DATA_SIZE/4],XY[DATA_SIZE/4],YX[DATA_SIZE/4], StokesI[175],*chan=NULL,nosig=0.0;
 float *beam1, *beam2, *beam3, *beam4, *beam5, *beam6, *beam7;
 glob_t globbuf;
 int colnum=4,idx=0;
 char timestamp[1024];
 float beam1_avg[TOT_CHANS],beam2_avg[TOT_CHANS],beam3_avg[TOT_CHANS],beam4_avg[TOT_CHANS],beam5_avg[TOT_CHANS],beam6_avg[TOT_CHANS],beam7_avg[TOT_CHANS],chidx[TOT_CHANS];
 char bank[20],alphabet[1];  
 int flow,fhigh, index;
 char basename[1024];


/* Options */
  while ((iNextOpt = getopt(argc,argv, "hf:l:g:")) != -1)
  { 
      switch (iNextOpt)
        {


          case 'h':
          print_usage();
          return(0);
          break;

          case 'f':
          sprintf(timestamp,"%s",optarg);;
          break;

          case 'l':
          flow = atoi(optarg);
          break;

          case 'g':
          fhigh = atoi(optarg);
          break;

        }
  
   } 


   if (0 == flow)
       {
          flow=0;
       }

     if (0  == fhigh)
       {
          fhigh = 500;
       }
 
     if (strcmp(timestamp,"") == 0)
     {
         printf("No time stamp given!\n");
         exit(0);
     }
/* parsing all the fits files */

 strcpy(bank,"AABCDEFGHIJKLMNOPQRST");
 sprintf(basename,"%s",timestamp);
 strcat(basename,"*.fits*");
 glob(basename,GLOB_ERR,NULL,&globbuf);

 if(globbuf.gl_pathv == NULL){
  fprintf(stderr,"Error: Fits files not found. The code should be run from the same directory as the fits files.\n");
  exit(1);
 }
 strcpy(filename,globbuf.gl_pathv[0]);

 fits_open_file(&fptr, filename, READONLY, &status);
 fits_movabs_hdu(fptr, 2, NULL, &status);
 fits_get_num_rows(fptr, &nrows, &status);
 fits_close_file(fptr,&status);
 printf("nrows=%ld\n",nrows);


/* opening all files */

 int flag[20],e;
for (e=1; e<=numFiles;e++){
  sprintf(basename,"%s",timestamp);
  sprintf(alphabet,"%c",bank[e]);
  strcat(basename,alphabet);
  strcat(basename,".fits");
  fits_open_file(&fptr1[e], basename, READONLY, &status);
  if (status)
   flag[e] = 1;
  else
   flag[e] = 0;
  status=0;
}


/* defining arrays */


beam1 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
beam2 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
beam3 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
beam4 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
beam5 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
beam6 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
beam7 =(float*)malloc(sizeof(float)*TOT_CHANS*nrows*NUM_SAMPS_PER_BLOCK);
chan = (float*)malloc(sizeof(float)*DATA_SIZE);



/* Converting the spliced files to filterbank file for each beam */
 
if (cpgopen("?") < 1)
     return 1;

 for(k=1; k<=nrows;k++){
    step_size=0;

    for (i=0; i<NUM_SAMPS_PER_BLOCK;i++){

      for (o=0;o<5;o++){

       for(j=1;j<=numFiles;j++){

        if (flag[j] == 0){
          fits_movabs_hdu(fptr1[j], 2, NULL, &status);
          fits_report_error(stderr,status);
          fits_read_col(fptr1[j], TFLOAT,colnum, k, 1, DATA_SIZE ,NULL, chan,NULL, &status);
          fits_report_error(stderr,status);
         }
        else{

          memset(chan,0,sizeof(float)*DATA_SIZE); 
 
        }
/* Converting the chunks to contiguous frequency */
        int index=7,count=0;
        for (m = 0;m < 5; m++ ){

           idx = 28*m + 140*o + 700*i;
         for (n=0;n<NUM_BEAMS;n++){

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
        for(l=0;l<5;l++){
 
           ind = NUM_BEAMS*l;
           beam1[step_size] = StokesI[ind];
           beam2[step_size] = StokesI[ind+1];
           beam3[step_size] = StokesI[ind+2];
           beam4[step_size] = StokesI[ind+3];
           beam5[step_size] = StokesI[ind+4];
           beam6[step_size] = StokesI[ind+5]; 
           beam7[step_size] = StokesI[ind+6];
           step_size++ ;   
           
         }
 
       }
 
      }

     }

/* Averaging the spectra */
    for (z=flow;z < fhigh;z++){
      beam1_avg[z] = (beam1[z] + beam1[500 + z ] + beam1[1000 + z] + beam1[1500 + z] + beam1[2000 + z] + beam1[2500 + z]+beam1[3000 + z] + beam1[3500 + z]+beam1[4000 + z] + beam1[4500 + z])/10;
      beam2_avg[z] = (beam2[z] + beam2[500 + z] + beam2[1000 + z] + beam2[1500 + z] + beam2[2000 + z] + beam2[2500 + z]+beam2[3000 + z] + beam2[3500 + z]+beam2[4000 + z] + beam2[4500 + z])/10;
      beam3_avg[z] = (beam3[z] + beam3[500 + z] + beam3[1000 + z] + beam3[1500 + z] + beam3[2000 + z] + beam3[2500 + z]+beam3[3000 + z] + beam3[3500 + z]+beam3[4000 + z] + beam3[4500 + z])/10;
      beam4_avg[z] = (beam4[z] + beam4[500 + z] + beam4[1000 + z] + beam4[1500 + z] + beam4[2000 + z] + beam4[2500 + z]+beam4[3000 + z] + beam4[3500 + z]+beam4[4000 + z] + beam4[4500 + z])/10;
      beam5_avg[z] = (beam5[z] + beam5[500 + z] + beam5[1000 + z] + beam5[1500 + z] + beam5[2000 + z] + beam5[2500 + z]+beam5[3000 + z] + beam5[3500 + z]+beam5[4000 + z] + beam5[4500 + z])/10;
      beam6_avg[z] = (beam6[z] + beam6[500 + z] + beam6[1000 + z] + beam6[1500 + z] + beam6[2000 + z] + beam6[2500 + z]+beam6[3000 + z] + beam6[3500 + z]+beam6[4000 + z] + beam6[4500 + z])/10;
      beam7_avg[z] = (beam7[z] + beam7[500 + z ] + beam7[1000 + z] + beam7[1500 + z] + beam7[2000 + z] + beam7[2500 + z]+beam7[3000 + z] + beam7[3500 + z]+beam7[4000 + z] + beam7[4500 + z])/10;

     }    
    printf("Plotting Bandpass for all beams ...\n");
/* Plotting the spectra */

    for (m=0;m<500;m++){
      chidx[m] = m;
    } 

    cpgsubp(3,3); 
/* Beam 1 */
    cpgpanl(2,2);
    cpgsch(1);
    cpgscr(2,0.1,0.3,0.5);
    float dataMin = FLT_MAX;
    float dataMax = -FLT_MAX;
    for (y= 0; y < 500; y++){

    if (beam1_avg[y] > dataMax)
           dataMax = beam1_avg[y];
    if (beam1_avg[y] < dataMin)
           dataMin = beam1_avg[y];
    }
    cpgsvp(0.1, 0.9 ,0.1, 0.9);
    cpgswin(chidx[0], chidx[499], dataMin, dataMax);
    cpgsci(3);
    cpgline(500, chidx, beam1_avg);
    cpglab("Channel", "Amplitude", "Beam 1");
    cpgsci(1);
    cpgbox("BCNT", 0.0, 100, "BCNT", 0.0, 100);
   
/* Beam2*/

    cpgpanl(1,3);
    cpgsch(0.5);
    cpgscr(2,0.1,0.3,0.5);
    for (y=0; y < 500; y++){

    if (beam2_avg[y] > dataMax)
           dataMax = beam2_avg[y];
    if (beam2_avg[y] < dataMin)
           dataMin = beam2_avg[y];
    }

    cpgsvp(0.1, 0.9, 0.1, 0.9);

    cpgswin(chidx[0], chidx[499], dataMin, dataMax);
    cpglab("Channel", "Amplitude", "Beam 2");
    cpgsci(3);
    cpgline(500, chidx, beam2_avg);
    cpgsci(1);
    cpgbox("BCNT", 0.0, 100, "BCNT", 0.0, 100);

/*Beam 3 */
    cpgpanl(1,2);
    cpgsch(0.5);
    cpgscr(2,0.1,0.3,0.5);
    for (y=0; y < 500; y++){

    if (beam3_avg[y] > dataMax)
           dataMax = beam3_avg[y];
    if (beam3_avg[y] < dataMin)
           dataMin = beam3_avg[y];
    }

    cpgsvp(0.1, 0.9, 0.1, 0.9);

    cpgswin(chidx[0], chidx[499], dataMin, dataMax);
    cpglab("Channel", "Amplitude", "Beam 3");
    cpgsci(3);
    cpgline(500, chidx, beam3_avg);
    cpgsci(1);
    cpgbox("BCNT", 0.0, 100, "BCNT", 0.0, 100);
/* Beam 4 */
    
    cpgpanl(1,1);
    cpgsch(0.5);
    cpgscr(2,0.1,0.3,0.5);
    for (y=0; y < 500; y++){

    if (beam4_avg[y] > dataMax)
           dataMax = beam4_avg[y];
    if (beam4_avg[y] < dataMin)
           dataMin = beam4_avg[y];
    }

    cpgsvp(0.1, 0.9, 0.1, 0.90);

    cpgswin(chidx[0], chidx[499], dataMin, dataMax);
    cpglab("Channel", "Amplitude", "Beam 4");
    cpgsci(3);
    cpgline(500, chidx, beam4_avg);
    cpgsci(1);
    cpgbox("BCNT", 0.0, 100, "BCNT", 0.0, 100);
    
 
/* Beam 5 */

    cpgpanl(3,1);
    cpgsch(0.5);
    cpgscr(2,0.1,0.3,0.5);
    for (y=0; y < 500; y++){

    if (beam5_avg[y] > dataMax)
           dataMax = beam5_avg[y];
    if (beam5_avg[y] < dataMin)
           dataMin = beam5_avg[y];
    }

    cpgsvp(0.1, 0.9, 0.1, 0.90);

    cpgswin(chidx[0], chidx[499], dataMin, dataMax);
    cpglab("Channel", "Amplitude", "Beam 5");
    cpgsci(3);
    cpgline(500, chidx, beam5_avg);
    cpgsci(1);
    cpgbox("BCNT", 0.0, 100, "BCNT", 0.0, 100);

/* Beam 6 */
   
    cpgpanl(3,2);
    cpgsch(0.5);
    cpgscr(2,0.1,0.3,0.5);
    for (y=0; y < 500; y++){

    if (beam6_avg[y] > dataMax)
           dataMax = beam6_avg[y];
    if (beam6_avg[y] < dataMin)
           dataMin = beam6_avg[y];
    }

    cpgsvp(0.1, 0.9, 0.1, 0.90);

    cpgswin(chidx[0], chidx[499], dataMin, dataMax);
    cpglab("Channel", "Amplitude", "Beam 6");
    cpgsci(3);
    cpgline(500, chidx, beam6_avg);
    cpgsci(1);
    cpgbox("BCNT", 0.0, 100, "BCNT", 0.0, 100);

    cpgpage();
/* Beam 7 */
    cpgpanl(3,3);
    cpgsch(1);
    cpgscr(2,0.1,0.3,0.5);
    for (y=0; y < 500; y++){

    if (beam7_avg[y] > dataMax)
           dataMax = beam7_avg[y];
    if (beam7_avg[y] < dataMin)
           dataMin = beam7_avg[y];
    }

    cpgsvp(0.1, 0.9, 0.1, 0.9);
    cpgsci(1);
    cpgswin(chidx[0], chidx[499], dataMin, dataMax);
    cpglab("Channel", "Amplitude", "Beam 7");
    cpgsci(3);
    cpgline(500, chidx, beam7_avg);
    cpgsci(1);
    cpgbox("BCNT", 0.0, 100, "BCNT", 0.0,100);
    cpgpage();

 }
 cpgclos();

 free(beam1);
 free(beam2);
 free(beam3);
 free(beam4);
 free(beam5);
 free(beam6);
 free(beam7);
 free(chan);
 
 return(0);

}
