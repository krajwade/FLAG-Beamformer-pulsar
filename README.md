This is a short description of the data reduction code for FLAG-beamformer 
backend. The code runs on VEGAS engineering fits files and processes them to 
store beamformed spectra out to a filterbank format. The FLAG-Beamformer 
outputs non-contiguous frequency channels for all seven beams which the reduction code takes into account before writing out the filterbank files for all 
seven beams. The code can be found on GitHub.

The way to install the code is as follows:

1. clone the repo from the GitHub Repository () 

2. change the 


# Things to be done.

1. Flip the band to make the file consistent with the expected format for SIGPROC and PRESTO

2. Change the rw_header subroutine to access the GB ancilliary fits files for relevant header information. Only RA and DEC needs to be updated.
# KMR @ WVU, 8th June 2017 
