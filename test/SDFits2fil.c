// Program to convert beamformed SDFITS data to filterbank


#include"/usr/include/cfitsio/fitsio.h"
#include<math.h>
#include<stdio.h>
#include<glob.h>
#include<string.h>


#define NUMROWS 1812
#define NUMCHANS 16384

// Main program

int main(int argc, char *argv[]){

// Checking for arguments

 if (argc<2) {
    printf("usage: %s sdfits filename\n",argv[0]);
    exit(0);
  }

// declaring variables
 
 fitsfile *fptr;
 FILE *fptr1=NULL;
 float XX[NUMCHANS],YY[NUMCHANS],I[NUMCHANS];
 int status=0,rownum;
 char filename[128];
 
 strcpy(filename,"sample.dat");
// strcat(filename,argv[2]);

 fptr1 = fopen(filename,"wb");
 if(fptr1 == NULL){
 printf("Error opening file\n");
 exit(1);
 }

 int i,j,k;
// Parsing data from the table

 fits_open_file(&fptr,argv[1],READONLY,&status);
 fits_movabs_hdu(fptr,2,NULL,&status);

 for (i=0;i<NUMROWS;i++){

     rownum = 2*i + 1;
     fits_read_col(fptr,TFLOAT,7,rownum,1,NUMCHANS,NULL,XX,NULL,&status);
     fits_read_col(fptr,TFLOAT,7,rownum+1,1,NUMCHANS,NULL,YY,NULL,&status);     

// Creating Stokes I
     for (j=0;j<NUMCHANS;j++){

       I[j] = sqrt(pow(XX[(NUMCHANS-1)-j],2) + pow(YY[(NUMCHANS-1)-j],2));

     }

     fwrite(I,sizeof(float),NUMCHANS,fptr1);

 }


 fclose(fptr1);
 fits_close_file(fptr,&status);

 return(0);
}
