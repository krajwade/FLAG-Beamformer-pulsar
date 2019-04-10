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

double radeg2hms(double ra)
{
    float ra_m, ra_s;
    char rastr[1024], ra_hstr[128], ra_mstr[128], ra_sstr[128];
    sprintf(ra_hstr, "%lf",floor(ra/15.0));
    ra_m = (ra/15 - floor(ra/15.0))*60.0;
    sprintf(ra_mstr, "%lf",ra_m - floor(ra_m));
    ra_s = (ra_m - floor(ra_m))*60.0;
    sprintf(ra_sstr, "%lf",ra_s);
    strcat(rastr,ra_hstr);
    strcat(rastr,ra_mstr);
    strcat(rastr,ra_sstr);
    return atof(rastr);
}

double decdeg2dms(double dec)
{
    float dec_m, dec_s;
    char decstr[1024], dec_dstr[128], dec_mstr[128], dec_sstr[128];
    if (dec < 0.0)
    {
      float val = abs(dec); 
      sprintf(dec_dstr, "%f",floor(val));
      dec_m = (val - floor(val))*60.0;
      sprintf(dec_mstr, "%lf",dec_m - floor(dec_m));
      dec_s = (dec_m - floor(dec_m))*60.0;
      sprintf(dec_sstr, "%lf",dec_s);
      strcat(decstr,dec_dstr);
      strcat(decstr,dec_mstr);
      strcat(decstr,dec_sstr);
      return -1.0*atof(decstr);
    }

    else
    {
      sprintf(dec_dstr, "%lf",floor(dec));
      dec_m = (dec - floor(dec))*60.0;
      sprintf(dec_mstr, "%lf",dec_m - floor(dec_m));
      dec_s = (dec_m - floor(dec_m))*60.0;
      sprintf(dec_sstr, "%lf",dec_s);
      strcat(decstr,dec_dstr);
      strcat(decstr,dec_mstr);
      strcat(decstr,dec_sstr);
      return atof(decstr);
    }
}

double getRA(fitsfile *fp)
{

    int status=0;
    //FILE* fptr=NULL;
    //char command[1024];
    static double ra;
    double ra_fin;
    fits_read_key(fp,TDOUBLE, "RA", &ra, NULL, &status);
    if (status)
    {          /* print any error messages */
	    fits_report_error(stderr, status);
	    return(status);
    }

    /*fits_read_key(fp,TDOUBLE, "DEC", &dec, NULL, &status);
    if (status)
    {           print any error messages */
    /*	    fits_report_error(stderr, status);
	    return(status);
    }*/

    /*sprintf(command,"python /users/krajwade/bf/FLAG-Beamformer-pulsar/src/ra_conv.py %lf %lf",ra, dec);
    fptr = popen(command,"r");

    if (fptr == NULL){
      printf("Error!!\n") ;
      exit(1);
    }
    fread(ra_fin, sizeof(double),1,fptr);*/
    ra_fin = radeg2hms(ra);
    return ra_fin;
}


double getDec(fitsfile* fp)
{

    int status=0;
   // FILE* fptr=NULL;
 // char command[1024];
    static double dec;
    double dec_fin;
    fits_read_key(fp,TDOUBLE, "DEC", &dec, NULL, &status);

    if (status)
    {         /*  print any error messages */
    	    fits_report_error(stderr, status);
	    return(status);
    }

    /*fits_read_key(fp,TDOUBLE, "RA", &ra, NULL, &status);

    if (status)
    {           print any error messages */
    /*        fits_report_error(stderr, status);
	    return(status);
    }*/


    /*sprintf(command,"python /users/krajwade/bf/FLAG-Beamformer-pulsar/src/dec_conv.py %lf %lf",ra, dec);
    fptr = popen(command,"r");

    if (fptr == NULL){
      printf("Error!!\n") ;
      exit(1);
    }
    fread(dec_fin, sizeof(double),1,fptr);*/
    dec_fin = decdeg2dms(dec);
    return dec_fin;

}

double getFreq(fitsfile* fp)
{

    int status=0;

    static double freq;
    fits_read_key(fp,TDOUBLE, "RESTFRQ", &freq, NULL, &status);

    if (status)
    {          /* print any error messages */
	    fits_report_error(stderr, status);
	    return(status);
    }

    return freq/1000000.0;

}
int rw_header(fitsfile *fptr,FILE *fptr1, char* proj_id, char* timestamp, int qflag, int fflag)
{
 

  int  status=0,nkeys,ilen=0;
  char record[FLEN_CARD],headstring[1024],strval[1024];
  double keyval;
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
  if (status){          /* print any error messages */
     fits_report_error(stderr, status);
     return(0);
   }


// Opening the Corresponding GO fits file
  fitsfile *gofile;
  char path[1024];
  char timestampcpy[1024];
  strcpy(timestampcpy, timestamp);
  strcpy(path,"/home/gbtdata/");
  strcat(path,proj_id);
  strcat(path, "/GO/");
  strcat(timestampcpy,".fits");
  strcat(path,timestampcpy);

// DoI need to add the null char here?

  fits_open_file(&gofile, path, READONLY, &status);

  if (status) 
  {         /* print any error messages */
     fits_report_error(stderr, status);
     exit(1);
  }


//RA
  keyval = getRA(gofile);                         // Still not sure where to get the information for RA
  ilen = strlen("src_raj");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"src_raj");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&keyval,sizeof(keyval),1,fptr1);
// writing to file

//DEC
  keyval = getDec(gofile);                        // Same for DEC
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
  if (qflag == 1)
      iTemp = 8;
  else
      iTemp = 32;
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
  if (fflag == 1)
  	keyval = getFreq(gofile) + (250*0.30318) ; // Getting to the highest frequency channel
  else
	keyval = getFreq(gofile) - (250*0.30318);
  ilen = strlen("fch1");
  fwrite(&ilen,sizeof(ilen),1,fptr1);
  strcpy(headstring,"fch1");
  fwrite(headstring,sizeof(char),ilen,fptr1);
  fwrite(&keyval,sizeof(keyval),1,fptr1);

  if (fflag == 1)
  	keyval = -0.30318;
  else
	keyval = 0.30318;
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
