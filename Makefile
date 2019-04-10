# Compiler
CC = gcc

NVCC= nvcc

# Directories
SRCDIR = src
BINDIR = bin
IDIR = include
# All the required flags

CFLAGS = -std=gnu99 -Wall -O3 -g

NVFLAGS = -g -Xptxas -arch=sm_52

IFLAGS = -I$(IDIR)

LFLAGS_CUDART = -L /usr/local/cuda/lib64  -lcudart -lcurand

LFLAGS_MATH = -lm

LFLAGS_FITSIO =  -L /usr/lib64/ -lcfitsio

OPENMP= -fopenmp

BINARY = bf_fits2fil_opt
# Final compilation

all: bf_fits2fil_opt


bf_fits2fil_opt: $(SRCDIR)/bf_fits2fil_opt.c $(SRCDIR)/rw_header.c cuda_utils.o madfilter_small.o
	$(CC) $(CFLAGS) $(OPENMP)  $^ $(LFLAGS_FITSIO) $(LFLAGS_CUDART) $(LFLAGS) $(IFLAGS) -o $(BINDIR)/$@

cuda_utils.o: $(SRCDIR)/cuda_utils.cu
	$(NVCC) -c $^ $(IFLAGS)

madfilter_small.o: $(SRCDIR)/madfilter_small.cu
	$(NVCC) -c $^ $(IFLAGS)


clean:
	rm -rf $(BINDIR)/$(BINARY)
	rm -rf *.o
