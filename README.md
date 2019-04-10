This is a short description of the data reduction code for FLAG-beamformer 
backend. The code runs on VEGAS engineering fits files and processes them to 
store beamformed spectra out to a filterbank format. The FLAG-Beamformer 
outputs non-contiguous frequency channels for all seven beams which the reduction code takes into account before writing out the filterbank files for all 
seven beams. The code can be found on GitHub.
Note that the following dependencies are needed:

1. cfitsio
2. cuda 8.0 toolkit

# Installation and Usage

1. clone the repo from the GitHub Repository (https://github.com/krajwade/FLAG-Beamformer-pulsar) 

2. change the bin directory if need be to where you want the binary files to be installed.i Make sure you have a bin directory already made.

3. Compile the code by running make in the folder with the Makefile.i If fitsio.h is not in standard path, change the path of the include file in include/rw_header.h. Another path you need to change in the Makefile is the LFLAGS_CUDART to the path where cuda libraries are installed.

4. Usage: path to the binary file basename of the fits file. for e.g. /users/krajwade/bf/Comm/bf_fits2fil 2017_05_24_20:07:21

5. The user has the option to: a. Qunatize to 8-bits, b) Flip the channel order and c) Run a GPU-based MAD filter (only on 8-bit data).

Current status:

a) A bit of cleanup still needed.
b) The DC component is removed in a very hackey way and only in the quantize block. Need to implement a better way to get rid of the DC component in 8 and 32-bit modes.



