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

parser = argparse.ArgumentParser()
parser.add_argument('-f','--all_files', metavar="path", type=str,\
        help="Path to files to be merged; enclose in quotes, accepts * as wildcard for directories or filenames")
parser.add_argument('-m', '--mem', type=float, help = "% RAM you want to use [0-1], def = 0.1", default = 0.1)

args = parser.parse_args()
files = glob.glob(args.all_files)
mem_percent = float(args.mem)

if not files:
    print('File(s) does not exist: ' + args.all_files) #, file=sys.stderr)
for file in files:
    print('Running using files: ' + file)


# Get RAM info
meminfo = dict((i.split()[0].rstrip(':'),int(i.split()[1])) for i in open('/proc/meminfo').readlines())
mem_mib = meminfo['MemTotal']/1024

# Get nrows corresponding to a given fraction of RAM
# NOTE: Even if you give -m 1, I will scale it down to 75%
# Hence the 0.75 factor
# 100*25*7*4*20*4e-6
# 100 samples per row
# 25 chanels per bank
# 7 beams
# 4 pols
# 20 bank files
# 4e-6 sizeof(32 bit float)
if 0 < mem_percent <= 1:
     nrows=int(0.75*mem_mib*mem_percent/(100*25*7*4*20*4e-6))
else:
    raise ValueError("percent RAM usage cannot be > 1")

# Get num_runs

#hdu=fits.open(file)
data, header = fits.getdata(file, header=True)
tot_rows=int(header["NAXIS2"])
print "total rows : ", tot_rows

nrow_list=list(range(0,tot_rows, nrows))
nrow_list.append(tot_rows)


# This generates all the 20 bank label from A to T.
bank_labels = [chr(i) for i in range(ord('A'),ord('T')+1)]

# instead of filling zeores later, make a zeros array and fill with useful
# values later in the for loop
band_pass = np.zeros(shape=(nrows,100,500,4,7), dtype=np.float32)
# the last few rows will be left out so last_pass will contain them
last_pass = np.zeros(shape=(nrow_list[-1]-nrow_list[-2],100,500,4,7), dtype=np.float32)
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


# get dummy_header
fil_file='/lustre/projects/flag/survey_filterbank/S02/B1933+16_1/BF_pulsar_0.fil'
header=fb.read_header(fil_file)

# open files to dump header
for i in range(7): #out_file_names:
    fb.create_filterbank_file(out_file_names[i],header=header[0],nbits=32)

for rows in range(len(nrow_list)-1):
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
                try:
                    band_pass[:,:,freq_lists[bank_freq_index],:,:]=data[nrow_list[rows]:nrow_list[rows+1],:].reshape(nrows,100,25,4,7)
                except ValueError:
                    last_pass[:,:,freq_lists[bank_freq_index],:,:]=data[nrow_list[rows]:nrow_list[rows+1],:].reshape(nrows,100,25,4,7)
                    last_pass_flag=True
                hdu.close()
        bank_freq_index+=1

    print "writing rows : ", nrow_list[rows], "to" , nrow_list[rows+1]
    # write the 7 beams
    if last_pass_flag:
        band_pass = last_pass
    for file_num in range(7):
        file=fb.FilterbankFile(out_file_names[file_num], mode='append')
        file.append_spectra((band_pass[:,:,:,1,file_num]+band_pass[:,:,:,0,file_num]).reshape(100*nrows,500))


print "Closing Files"
for file_num in range(7):
    file=fb.FilterbankFile(out_file_names[file_num], mode='append')
    file.close()
