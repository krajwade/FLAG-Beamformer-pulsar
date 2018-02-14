#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import sys
import glob
import argparse
# change presto path here
sys.path.append("/opt/pulsar/src/PRESTOv2.7.13/lib/python/")
import filterbank as fb 
import sigproc

# written by Devansh Agarwal
# Feb 12, 2018
# mail: devanshkv@gmail.com

# A function to add a SigProc Header

def getsigprocheader(file):

  telescope_ids = {"Fake": 0, "Arecibo": 1, "Ooty": 2, "Nancay": 3,
                 "Parkes": 4, "Jodrell": 5, "GBT": 6, "GMRT": 7,
                 "Effelsberg": 8}

  machine_ids = {"WAPP": 0, "PSPM": 1, "Wapp": 2,"AOFTM": 3,
               "BCPM1": 4, "FLAGBF": 5, "SCAMP": 6,
               "GBT Pulsar Spigot": 7, "SPIGOT": 7, "PUPPI": 8 ,
               "GUPPI": 9,"PA":10,"VEGAS": 11}


  hdulist = fits.open(file)


  def prep_string(string):
    return struct.pack('i', len(string))+string

  def prep_double(name, value):
    return prep_string(name)+struct.pack('d', float(value))

  def prep_int(name, value):
    return prep_string(name)+struct.pack('i', int(value))

  hdr = prep_string("HEADER_START")
  hdr += prep_int("telescope_id", 6)
  hdr += prep_int("machine_id",5)
  hdr += prep_int("data_type", 1) # 1 = filterbank, 2 = timeseries

  source =hdulist[0].header['OBJECT']

#TODO: parse RA DEC from the ancillary fits files

  hdr += prep_string("source_name")
  hdr += prep_string(source)
  hdr += prep_int("barycentric", 1)
  hdr += prep_int("pulsarcentric", 0)
  hdr += prep_double("src_raj",12345)
  hdr += prep_double("src_dej", 12345)
  hdr += prep_int("nbits", 32)
  hdr += prep_int("nifs", 1)

  hdr += prep_int("nchans", 500)
  hdr += prep_double("fch1", 1325.0)
  hdr += prep_double("foff",-0.30318)

  MJD = hdulist[5].header['UTDSTART']

  hdr += prep_double("tsamp",0.000130)

  return hdr


# Main code begins here


parser = argparse.ArgumentParser()
parser.add_argument('-f','--all_files', metavar="path", type=str,\
        help="Path to files to be merged; enclose in quotes, accepts * as wildcard for directories or filenames")
parser.add_argument('-m', '--mem', type=float, help = "% RAM you want to use [0-1], def = 1", default = 1.0)
parser.add.argument('-fb' '--flipband', type=int,help = "Flag to flip the band pass while storing to filterbank",default=1)

args = parser.parse_args()
all_files = glob.glob(args.all_files)
mem_percent = float(args.mem)
flag_flip = int(args.flipband)

# if the files don't exist exit
if not all_files:
    print('File(s) does not exist: ' + args.all_files) #, file=sys.stderr)
    sys.exit()

# check for dropped banks

row_nos=[]
for file in all_files:
    junk, header = fits.getdata(file, header=True)
    junk=None
    row_nos.append(int(header["NAXIS2"]))

# get the mode
tot_rows = max(row_nos, key=row_nos.count)

# take files with same number of rows
files=[]
for idx, file in enumerate(all_files):
    if row_nos[idx]==tot_rows:
        print('Running using files: ' + file)
        files.append(file)


# Get RAM info
meminfo = dict((i.split()[0].rstrip(':'),int(i.split()[1])) for i in open('/proc/meminfo').readlines())
mem_mib = meminfo['MemTotal']/1024


# Get nrows corresponding to a given fraction of RAM
# NOTE: Even if you give -m 1, I will scale it down to 50%
# Hence the 0.5 factor
# 100*25*7*4*20*4e-6
# 100 samples per row
# 25 chanels per bank
# 7 beams
# 4 pols
# 20 bank files
# 4e-6 sizeof(32 bit float)
if 0 < mem_percent <= 1:
     nrows=int(0.5*mem_mib*mem_percent/(100*25*7*4*20*4e-6))
else:
    raise ValueError("percent RAM usage cannot be > 1")

# Get num_runs

print "total rows : ", tot_rows
print "row step   : ", nrows
nrow_list=list(range(0,tot_rows, nrows))
nrow_list.append(tot_rows)
lrows=nrow_list[-1]-nrow_list[-2]


# This generates all the 20 bank label from A to T.
bank_labels = [chr(i) for i in range(ord('A'),ord('T')+1)]

# instead of filling zeores later, make a zeros array and fill with useful
# values later in the for loop
band_pass = np.zeros(shape=(nrows,100,500,4,7), dtype=np.float32)
# the last few rows will be left out so last_pass will contain them
#last_pass = np.zeros(shape=(lrows,100,500,4,7), dtype=np.float32)
last_pass_flag=False

# the following chunk makes up the vector to fix the frequency structure.
# as the frequencies are not contigous the following makes arrays as 
# 0-4,100-104,..,400-404 For bank A
# 5-9,105,109,..,405-409 For bank B and so on

p=0
freq_lists=np.zeros((20,25),dtype=np.int)
for _ in range(20):
    freq_lists[_,range(5)]=range(p,5+p)
    freq_lists[_,range(5,10)]=range(100+p,105+p)
    freq_lists[_,range(10,15)]=range(200+p,205+p)
    freq_lists[_,range(15,20)]=range(300+p,305+p)
    freq_lists[_,range(20,25)]=range(400+p,405+p)
    p=p+5

# dummy files to write_out_stuff
out_file_names=["BF_beam_%i.fil" %i for i in range(7)]

# a function that creates a sigproc header from the metadata in the fits files
header = getsigprocheader(files[0])


# open files to dump header
for i in range(7): #out_file_names:
    fb.create_filterbank_file(out_file_names[i],header=header[0],nbits=32)

for rows in range(len(nrow_list)-1):
    print "Manipulating rows : ", nrow_list[rows], "to" , nrow_list[rows+1]
    if rows==len(nrow_list)-2:
        last_pass = np.zeros(shape=(lrows,100,500,4,7), dtype=np.float32)
        last_pass_flag=True
        band_pass = None
    # to loop through all bank frequecies
    bank_freq_index=0
    # loops through all banks
    for abc in bank_labels:
        # open the ccorresponding file to make the frequency structure
        for file in files:
            # find which bank file to use
            if abc in file:
                # open the fits and get the data
                hdu=fits.open(file)
                data=hdu[1].data['DATA']
                # the data is this is 100 (time samples) X 25 (chans) X 4 (pol)  X 7 (beams)
                # if running for last few rows use last_pass, and raise the flag!
                if not last_pass_flag:
                    band_pass[:,:,freq_lists[bank_freq_index],:,:]=data[nrow_list[rows]:nrow_list[rows+1],:].reshape(nrows,100,25,4,7)
                else:
                    last_pass[:,:,freq_lists[bank_freq_index],:,:]=data[nrow_list[rows]:nrow_list[rows+1],:].reshape(lrows,100,25,4,7)
                hdu.close()
        bank_freq_index+=1

    print "writing rows : ", nrow_list[rows], "to" , nrow_list[rows+1]
    # write the 7 beams
   # if last_pass_flag:
   #     band_pass = last_pass
    for file_num in range(7):
        file=fb.FilterbankFile(out_file_names[file_num], mode='append')
        if last_pass_flag:
          if(flag_flip==1):
            last_pass_f = (last_pass[:,:,:,1,file_num]+last_pass[:,:,:,0,file_num]).reshape(100*lrows,500)
            last_pass_flip = np.flip(last_pass_f,axis=1)
            file.append_spectra(last_pass_flip)
          else:
            file.append_spectra((last_pass[:,:,:,1,file_num]+last_pass[:,:,:,0,file_num]).reshape(100*lrows,500))
        else:
          if (flag_flip==1):
            band_pass_f = (band_pass[:,:,:,1,file_num]+band_pass[:,:,:,0,file_num]).reshape(100*nrows,500)
            band_pass_flip = np.flip(band_pass_f,axis-1)
            file.append_spectra((band_pass[:,:,:,1,file_num]+band_pass[:,:,:,0,file_num]).reshape(100*nrows,500))
        
          else:
            file.append_spectra((band_pass[:,:,:,1,file_num]+band_pass[:,:,:,0,file_num]).reshape(100*nrows,500))
            

print "Closing Files"
for file_num in range(7):
    file=fb.FilterbankFile(out_file_names[file_num], mode='append')
    file.close()
