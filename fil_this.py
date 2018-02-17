#!/usr/bin/env python

from astropy.io import fits
import numpy as np
import sys
import glob
import argparse
import fml
import re
from itertools import chain
# change presto path here
sys.path.append("/opt/pulsar/src/PRESTOv2.7.13/lib/python/")
import filterbank as fb 
from astropy import units as u
from astropy.coordinates import SkyCoord

# written by Devansh Agarwal
# Feb 12, 2018
# mail: devanshkv@gmail.com

parser = argparse.ArgumentParser()
parser.add_argument('-f','--all_files', metavar="path", type=str,\
        help="Path to files to be merged; enclose in quotes, accepts\
        * as wildcard for directories or filenames", required=True)
parser.add_argument('-m', '--mem', type=float, help = "% RAM you \
        want to use [0-1] (def=0.6)", default = 0.6, required=False)
parser.add_argument('--noflip', nargs='?', const=False, default=True,\
        help='Flag to flip the band pass while storing to filterbank (def=True)',\
        required=False)
parser.add_argument('--ignorechan', type=str, required=False, help = "number/range\
        of channels to ignore like '0-5' or '2' or add multiple ranges like '0-2,4-6'",\
        default=None)

args = parser.parse_args()
all_files = glob.glob(args.all_files)
mem_percent = float(args.mem)
flip = args.noflip
ignorechannels=args.ignorechan


# do stuff for ignorechan
if ignorechannels:
    get_all_ranges=ignorechannels.split(',')
    def parseNumList(string):
        m = re.match(r'(\d+)(?:-(\d+))?$', string)
        start = m.group(1)
        end = m.group(2) or start
        return list(range(int(start,10), int(end,10)+1))
    ic=[]
    for ranges in get_all_ranges:
        ic.append(parseNumList(ranges))
    ignorechan=list(chain.from_iterable(ic))
    if 500 in ignorechan:
        raise ValueError("number of channels range from 0-499")
        sys.exit()
    print "ignoring channels : ", ignorechan

print "flip : ", flip
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

# Get header stuff here
def get_loc(proj_name, f_name):
    dir="/home/gbtdata/"+str(proj_name)+"/Antenna/"+str(f_name)+".fits"
    hdu=fits.open(dir)
    DEC=np.mean(hdu[2].data['DECJ2000'])
    RAJ=np.mean(hdu[2].data['RAJ2000'])
    hdu.close()
    dir="/home/gbtdata/"+str(proj_name)+"/LO1B/"+str(f_name)+".fits"
    hdu=fits.open(dir)
    CFREQ= hdu[3].data['LO1FREQ'][0]/1e6
    hdu.close()
    c = SkyCoord(ra=RAJ, dec=DEC, frame='icrs', unit='deg')
    DEC= c.dec.dms[0]*10000 + c.dec.dms[1]*100 + c.dec.dms[2]
    RAJ= c.ra.hms[0]*10000 + c.ra.hms[1]*100 + c.ra.dms[2]
    return RAJ, DEC, CFREQ 

header=fits.getheader(file)
MJD = float(header['STRTDMJD'])
proj_name=header['BCALFILE']
f_name = header['TSTAMP']
RAJ, DEC, CFREQ = get_loc(proj_name, f_name)

# get the mode
tot_rows = max(row_nos, key=row_nos.count)

# take files with same number of rows
files=[]
for idx, file in enumerate(all_files):
    if row_nos[idx]==tot_rows:
        #print('Running using files: ' + file)
        files.append(file)
print "Using ", len(files), "banks"


# Get RAM info

mem_mib = fml.FreeMemLinux(unit='MB').user_free # meminfo['MemTotal']/1024


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
     nrows=int(mem_mib*mem_percent/(100*25*7*4*20*4e-6))
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

# make header here
hdrdict={}
hdrdict["telescope_id"]=int(6)
hdrdict["machine_id"]=int(0)
hdrdict["data_type"]=int(1)
hdrdict["source_name"]=str("source")
hdrdict["barycentric"]=int(1)
hdrdict["pulsarcentric"]=int(0)
hdrdict["src_raj"]=float(RAJ)
hdrdict["src_dej"]=float(DEC)
hdrdict["nbits"]=int(32)
hdrdict["nifs"]=int(1)
hdrdict["nchans"]=int(500)
if flip:
    hdrdict["fch1"]=float(CFREQ+65.0)
    hdrdict["foff"]=float(-0.30318)
else:
    hdrdict["fch1"]=float(CFREQ-85.0)
    hdrdict["foff"]=float(0.30318)
hdrdict["tstart"]=float(MJD)
hdrdict["tsamp"]=float(0.000130)

## open files to dump header
for i in range(7): #out_file_names:
    fb.create_filterbank_file(out_file_names[i],header=hdrdict, nbits=32)
#
for rows in range(len(nrow_list)-1):
    print "Manipulating rows : ", nrow_list[rows], "to" , nrow_list[rows+1]
    if rows==len(nrow_list)-2:
        # when doing the last bit make the zeros array
        # set the flag and empty the memory taken by band_pass
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
    # band flip if needed
    # ignore chan in needed
    if ignorechannels:
        if last_pass_flag:
                last_pass[:,:,ignorechan,:,:]=0
        else:
                band_pass[:,:,ignorechan,:,:]=0

    for file_num in range(7):
        file=fb.FilterbankFile(out_file_names[file_num], mode='append')
        if last_pass_flag:
            if flip:
                file.append_spectra((last_pass[:,:,:,1,file_num]+last_pass[:,:,:,0,file_num]).reshape(100*lrows,500)[:,::-1])
            else:
                file.append_spectra((last_pass[:,:,:,1,file_num]+last_pass[:,:,:,0,file_num]).reshape(100*lrows,500))
        else:
            if flip:
                file.append_spectra((band_pass[:,:,:,1,file_num]+band_pass[:,:,:,0,file_num]).reshape(100*nrows,500)[:,::-1])
            else:
                file.append_spectra((band_pass[:,:,:,1,file_num]+band_pass[:,:,:,0,file_num]).reshape(100*nrows,500))
            

print "Closing Files"
for file_num in range(7):
    file=fb.FilterbankFile(out_file_names[file_num], mode='append')
    file.close()
