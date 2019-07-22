#ifndef _CUDA_UTILS_H_
#define _CUDA_UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

void run_quant(float* h_pfbuf, int n, unsigned char* h_pcbuf);

void channel_flip_float(float* h_A, int size, int nchans);

void checkCUDAError(const char* msg);

#ifdef __cplusplus
}
#endif

#endif

