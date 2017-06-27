#include<stdio.h>
#include<stdlib.h>
#include"/usr/include/cfitsio/fitsio.h"
#include<complex.h>

/* A test to simultaneously read and write from a fits file */

int main(int argc, char *argv[])
{

  fitsfile *fptr1;
  FILE *fptr=NULL;
  int i,j,iter,status=0,ii,nkeys,a,keyval;
  char card[FLEN_CARD];
  long nrows;
  float elem[70000];
 
  keyval = atoi(argv[2]);
 
  fits_open_file(&fptr1,argv[1], READWRITE, &status);
  //fits_movabs_hdu(fptr1, 2, NULL, &status); 
  //fits_update_key(fptr1,TINT,"NAXIS2",&keyval,NULL,&status);
  //fits_close_file(fptr1, &status);

  //fits_open_file(&fptr1,argv[1], READWRITE, &status);  

  fptr = fopen("test.dat","wb");
  if (fptr == NULL){
    fprintf(stderr,"Error opening file\n");
    exit(1);
  }

  fits_movabs_hdu(fptr1, 1, NULL, &status);
  fits_get_hdrspace(fptr1, &nkeys, NULL, &status);
/*  for (ii = 1; ii <= nkeys; ii++)  {
          fits_read_record(fptr1, ii, card, &status);
          printf("%s\n", card);
        }*/

  fits_movabs_hdu(fptr1, 2, NULL, &status);
  for (i=0; i < 5; i++){

          fits_read_col(fptr1,TFLOAT,3,i+1,1,70000, NULL, elem,NULL,&status);
          fwrite(elem,sizeof(float),70000,fptr);
    }

  if (status)          /* print any error messages */
            fits_report_error(stderr, status);
          return(status);
 
  for (i=0;i< 5;i++){
          printf("%f\n",elem[i]);
  }

  fits_close_file(fptr1, &status); 
  fclose(fptr);
  if (status)          /* print any error messages */
            fits_report_error(stderr, status);
          return(status);
  return(0);  

}

