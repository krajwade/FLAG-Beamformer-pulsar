//Routine to read and write header parameters to the filterbank file //

#define _GNU_SOURCE 
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include"/usr/include/cfitsio/fitsio.h"
#include<math.h>
#include"rw_header.h"

/* Remove a single character from the  string */

char *del(char str[], char str1[128]) {
   int i;
   int size;
   char ch1;

   size = strlen(str);

   for (i = 0; i < size-1; i++) {
         ch1 = str[i+1];
         str1[i] = ch1;
   }
   str1[size] = '\0';
   strcpy(str,str1);
   return str;

   }


int rw_header(fitsfile *fptr,FILE *fptr1)
{
 

  int  status=0,nkeys,ilen=0;
  char record[FLEN_CARD],headstring[1024],strval[1024];
  double keyval;
  char pathname[128];
  int iTemp;



const char g_aacSP_ObsNames[NUMOBS][LEN_GENSTRING] = {
    OBS_FAKE,
    OBS_AO,
    OBS_ORT,
    OBS_NANCAY,
    OBS_PARKES,
    OBS_JB,
    OBS_GBT,
    OBS_GMRT,
    OBS_EFFELSBERG,
};

const char g_aacSP_MacNames[NUMMAC][LEN_GENSTRING] = {
    Mac_WAPP,
    Mac_PSPM,
    Mac_MOCK,
    Mac_AOFTM,
    Mac_BCPM,
    Mac_OOTY,
    Mac_SCAMP,
    Mac_SPIGOT,
    Mac_PUPPI,
    Mac_GUPPI,
    Mac_VEGAS,
    Mac_FLAGBF,
};



// Get observatory ID from Name

int GetObsIDFromName(char *pcObs)
{
    int i = 0;

    for (i = 0; i < NUMOBS; ++i)
    {
        if (strcasestr(pcObs, g_aacSP_ObsNames[i]) != NULL)
        {
            return i;
        }
    }

  return(0);
}

// Getting Backend from ID

int GetMacIDFromName(char *pcObs)
{
    int i = 0;

    for (i = 0; i < NUMMAC; ++i)
    {
        if (strcasestr(pcObs, g_aacSP_MacNames[i]) != NULL)
        {
            return i;
        }
    }

    return(0);
}

// Inserting the header parameters with the header keywords that header program in Sigproc uses to output the header.

  fits_movabs_hdu(fptr, 1, NULL, &status);
  fits_get_hdrspace(fptr, &nkeys, NULL, &status);

// Start
  strcpy(headstring,"HEADER_START");
  ilen = strlen(headstring);
  fwrite(&ilen,sizeof(int),1,fptr1);
  fwrite(headstring,sizeof(char),ilen,fptr1);
// writing to file
 

// Telescope 
  ilen=strlen("telescope_id");
  fwrite(&ilen,sizeof(ilen),1,fptr1); 
  strcpy(headstring,"telescope_id");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  iTemp = GetObsIDFromName(OBS_GBT);
  fwrite(&iTemp,sizeof(iTemp),1,fptr1);
// writing to file

  
//Backend
  ilen=strlen("machine_id");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"machine_id");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  iTemp = GetMacIDFromName(Mac_FLAGBF);
  fwrite(&iTemp,sizeof(iTemp),1,fptr1);
  
//writing to file

// Data type  
  ilen=strlen("data_type");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"data_type");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  iTemp=1;
  fwrite(&iTemp,sizeof(iTemp),1,fptr1);

//writing to file


// Write the barycentric flag
  iTemp = 0;
  ilen=strlen("barycentric");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"barycentric");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&iTemp,sizeof(iTemp),1,fptr1);
//

//Source
/* Opening the GO Fits file */

/* fits_open_file(&fptr2, go_file, READONLY, &status);*/
  fits_read_card(fptr, "OBJECT", record, &status); /* read keyval */
  fits_read_key(fptr,TSTRING, "OBJECT", &strval, NULL, &status);
  ilen = strlen("source_name");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"source_name");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  ilen = strlen(strval);
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  fwrite(strval,sizeof(char),ilen,fptr1);
// writing to file
  fits_close_file(fptr, &status);
  if (status){          /* print any error messages */
     fits_report_error(stderr, status);
     return(0);
   }



//RA
  keyval = 12345.0;                         // Still not sure where to get the information for RA
  ilen = strlen("src_raj");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"src_raj");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&keyval,sizeof(keyval),1,fptr1);
// writing to file


//DEC
  keyval = 12345.0;                        // Same for DEC
  ilen = strlen("src_dej");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"src_dej");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&keyval,sizeof(keyval),1,fptr1);
// writing to file

// nchans
  iTemp = 500;
  ilen = strlen("nchans");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"nchans");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&iTemp,sizeof(int),1,fptr1);


// Number of bits per sample
  iTemp = 8;
  ilen = strlen("nbits");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"nbits");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&iTemp,sizeof(iTemp),1,fptr1);

//writing to file


// Number of Ifs
  iTemp = 1;
  ilen = strlen("nifs");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"nifs");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&iTemp,sizeof(iTemp),1,fptr1);
// Writing to file  

 
// Frequency channels
// A routine to write frequency channels out
/*  strcpy(headstring,"FREQUENCY_START");
  ilen = strlen(headstring);
  fwrite(&ilen,sizeof(int),1,fptr1);
  fwrite(headstring,sizeof(char),ilen,fptr1);
  strcpy(headstring,"nchans");
  ilen = strlen(headstring);
  fwrite(&ilen,sizeof(int),1,fptr1);
  fwrite(headstring,sizeof(char),ilen,fptr1);
  iTemp = nchans;
  fwrite(&iTemp,sizeof(iTemp),1,fptr1);

  int i,j,k=0;
  double *frmhz=NULL;
  frmhz = (double *) malloc(sizeof(double)*nchans);
  for (i=0; i<5; i++){

    for  (j = 0; j<90;j++){
       ilen = strlen("fchannel");
       fwrite(&ilen,sizeof(int),1,fptr1);
       strcpy(headstring,"fchannel");
       fwrite(headstring,sizeof(char),ilen,fptr1);
       frmhz[k] = fch1[i] + j*foff;
       fwrite(&frmhz[k],sizeof(double),1,fptr1);
       k++;
    }
  
  }


  strcpy(headstring,"FREQUENCY_END");
  ilen = strlen(headstring);
  fwrite(&ilen,sizeof(int),1,fptr1);
  fwrite(headstring,sizeof(char),ilen,fptr1);
// Writing to file
*/

// Simple frequency write 
  keyval = 1525;
  ilen = strlen("fch1");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"fch1");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&keyval,sizeof(keyval),1,fptr1);

  keyval = -0.30318;
  ilen = strlen("foff");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"foff");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&keyval,sizeof(keyval),1,fptr1);

// MJD
  fits_read_card(fptr, "STRTDMJD", record, &status); /* read keyval */
  fits_read_key(fptr,TDOUBLE, "STRTDMJD", &keyval, NULL, &status);
  ilen = strlen("tstart");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"tstart");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&keyval,sizeof(keyval),1,fptr1);

//Writing to file


// Sampling time
//  fits_read_card(fptr, "DURATION", record, &status); /* read keyval */
//  fits_read_key(fptr,TDOUBLE, "DURATION", &keyval, NULL, &status);
  keyval = 0.00013;
  ilen = strlen("tsamp");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"tsamp");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&keyval,sizeof(keyval),1,fptr1);

//Writing to file

// End header
  strcpy(headstring,"HEADER_END");
  ilen = strlen(headstring);
  fwrite(&ilen,sizeof(int),1,fptr1);
  fwrite(headstring,sizeof(char),ilen,fptr1);  


  return(0);
}
