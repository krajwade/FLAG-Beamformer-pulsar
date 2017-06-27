#ifndef __RW_HEADER_H__
#define __RW_HEADER_H__

#include"/usr/include/cfitsio/fitsio.h"
#include<string.h>
#include<unistd.h>

int rw_header(fitsfile *fp, FILE *fptr);

enum tagObservatory
{
    OBSID_FAKE = 0,
    OBSID_AO,
    OBSID_ORT,
    OBSID_NANCAY,
    OBSID_PARKES,
    OBSID_JB,
    OBSID_GBT,
    OBSID_GMRT,
    OBSID_EFFELSBERG
};

#define NUMOBS          9   /* number of supported sites */

#define OBS_FAKE        "Fake"
#define OBS_AO          "Arecibo"
#define OBS_ORT         "Ooty"
#define OBS_NANCAY      "Nancay"
#define OBS_PARKES      "Parkes"
#define OBS_JB          "Jodrell"
#define OBS_GBT         "GBT"
#define OBS_GMRT        "GMRT"
#define OBS_EFFELSBERG  "Effelsberg"

#define LEN_GENSTRING   256   /* Length of generic string */


enum tagMachineID{
     
   MacID_WAPP=0,
   MacID_PSPM,
   MacID_MOCK,
   MacID_AOFTM,
   MacID_BCPM,
   MacID_OOTY,
   MacID_SCAMP,
   MacID_SPIGOT,
   MacID_PUPPI,
   MacID_GUPPI,
   MacID_VEGAS,
   MacID_FLAGBF
};
        
#define Mac_WAPP   "WAPP"
#define Mac_PSPM   "PSPM"
#define Mac_MOCK   "MOCK"
#define Mac_AOFTM  "AOFTM"
#define Mac_BCPM   "BCPM"
#define Mac_OOTY   "OOTY"
#define Mac_SCAMP  "SCAMP"
#define Mac_SPIGOT "SPIGOT"
#define Mac_PUPPI  "PUPPI"
#define Mac_GUPPI  "GUPPI"
#define Mac_VEGAS  "VEGAS"
#define Mac_FLAGBF "FLAGBF"

#define NUMMAC      12    /* NUmber of backends */

enum tagDataType{         /* Data Types */
    Format_Fil=0,
    Format_tim
};

#define iFormat_Fil "filterbank"
#define iFormat_tim "timeseries"

#define NUMTYP       2

#define foff         0.30318     /*Channel width in MHz*/
#endif /* __RW_HEADER_H__ */

