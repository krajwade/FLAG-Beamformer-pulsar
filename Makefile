# Compiler
CC = gcc

# Directories
SRCDIR = src
BINDIR = bin
# All the required flags

CFLAGS = -std=gnu99 -Wall

LFLAGS_MATH = -lm

LFLAGS_FITSIO =  -L /usr/lib64/ -lcfitsio

# Final compilation
all: bf_fits2fil

bf_fits2fil: $(SRCDIR)/bf_fits2fil.c $(SRCDIR)/rw_header.c
	$(CC) $(CFLAGS) $^ $(LFLAGS_FITSIO) $(LFLAGS_MATH) -o $(BINDIR)/$@
