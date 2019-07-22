#ifndef _MADFILTER_SMALL_H_
#define _MADFILTER_SMALL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include<stdio.h>
#include<sys/time.h>
#include<math.h>
#include<string.h>

void run_madfilter(unsigned char* h_cdata, int size, int bsize);

#ifdef __cplusplus
}
#endif

#endif
