# Compiler
CC = gcc

# Directories
SRCDIR = src
BINDIR = bin
# All the required flags

CFLAGS = -std=gnu99 -Wall -O3  -g

LFLAGS_MATH = -lm

LFLAGS_FITSIO =  -L /usr/lib64/ -lcfitsio

OPENMP= -fopenmp

BINARY = bf_fits2fil
# Final compilation
all: bf_fits2fil

bf_fits2fil: $(SRCDIR)/bf_fits2fil.c $(SRCDIR)/rw_header.c
	$(CC) $(CFLAGS) $(OPENMP) $^ $(LFLAGS_FITSIO) $(LFLAGS_MATH) -o $(BINDIR)/$@

clean:
	rm -rf $(BINDIR)/$(BINARY)
