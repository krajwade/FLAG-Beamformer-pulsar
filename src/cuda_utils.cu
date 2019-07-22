#include "cuda_utils.h"

/* Kernel for quantization */
__global__ void Quant(float pfbuf[], int pibuf[] , int n, float fMax, float fMin) {
   /* blockDim.x = threads_per_block                            */
   /* First block gets first threads_per_block components.      */
   /* Second block gets next threads_per_block components, etc. */
   int i = blockDim.x * blockIdx.x + threadIdx.x;

   float frange;
   float fIntMax = (float) (powf(2.0 , 8) - 1.0);
   frange = fMax - fMin;
   /* block_count*threads_per_block may be >= n */
   if (i < n) {
  	 pibuf[i] = (int) roundf(((pfbuf[i]-fMin)/frange)*fIntMax);
   }
}  /* Quant */


__global__ void reverseArrayBlockFloat(float *d_b , float *d_a )
{
 int bx = blockIdx.x , tx = threadIdx.x ; 
 int old_id = blockDim.x * bx+ tx ;


// GridDim.x gives no. of block in grid in X dimension
 int new_id = (blockDim.x * gridDim.x) - 1 -  old_id ; 

 
 d_b[old_id] = d_a[new_id ]; 

}

__global__ void reverseArrayBlockNewFloat(float *d_b , float *d_a, int binsize, int bins )
{
 int bx = blockIdx.x , tx = threadIdx.x ; 
 int tid = blockDim.x * bx+ tx ;
 if (tid < bins)
 {
 int lw = tid * (binsize);
 int uw = lw + binsize;
 int i;
 for (i = 0; i < binsize; i++)
 {
   d_b[i + lw] = d_a[uw -1 - i];
 } 

 }
}


/* Qunatization GPU-based. Still debugging! */
void run_quant(float* h_pfbuf,int nchans, unsigned char* h_pcbuf) 
{
   int i;
   float fMin,fMax;
   int *d_pibuf, *h_pibuf;
   float *d_pfbuf;
   int threads_per_block;
   int block_count;
   size_t size_f,size_i;

   size_f = nchans*sizeof(float);
   size_i = nchans*sizeof(int);
   
  /* find Min and max */

   for (i = 0; i < nchans; ++i){
      if (h_pfbuf[i] > fMax)
           fMax = h_pfbuf[i];
      if (h_pfbuf[i] < fMin)
           fMin = h_pfbuf[i];

   }
   i = 0;
   /* Allocate vectors in device memory */
   h_pibuf = (int*)malloc(nchans * sizeof(int));
   cudaMalloc(&d_pfbuf, size_f);
   cudaMalloc(&d_pibuf, size_i);

   float* ptr = h_pfbuf;

   /* Copy vectors from host memory to device memory */
   cudaMemcpy(d_pfbuf, ptr, size_f, cudaMemcpyHostToDevice);

   checkCUDAError("memcpy");
   /* Define block size */
   threads_per_block = 250;

   /* Define grid size.  If we just computed n/threads_per_block */
   /* we might get fewer threads than vector components.  Using  */
   /* ceil(n/threads_per_block) guarantees at least one thread   */
   /* per vector component.  The following formula is a kludge   */
   /* since it appears that the CUDA ceil function doesn't work  */
   /* correctly.                                                 */
   block_count = (nchans + threads_per_block - 1)/threads_per_block;

   /* Invoke kernel using block_count blocks, each of which  */
   /* contains threads_per_block threads                     */
   Quant<<<block_count, threads_per_block>>>(d_pfbuf, d_pibuf, nchans, fMax,fMin);

   /* Wait for the kernel to complete */
   cudaThreadSynchronize();

  
   checkCUDAError("kernel Invocation");

   unsigned char* cptr = h_pcbuf;
   /* Copy result from device memory to host memory */
   /* h_z contains the result in host memory        */
   cudaMemcpy(h_pibuf, d_pibuf, size_i, cudaMemcpyDeviceToHost);

   checkCUDAError("memcpy");

   for (i = 0; i < nchans; ++i)
   {
	*(cptr) = (h_pibuf[i] & (255));
	++cptr;
   }
   
   /* Free device memory */
   cudaFree(d_pfbuf);
   cudaFree(d_pibuf);
   /* Free host memory */
   free(h_pibuf); 
   return;
} 

void channel_flip_float( float *h_a, int size, int nchans)
{
    
    // pointer for device memory
    float *d_b, *d_a;
    float time_kernel;
    // define grid and block size
    int numThreadsPerBlock = 10;

    int bins = (int)size/nchans;
    // Part 1 of 2: compute number of blocks needed based on array size and desired block size
    int numBlocks = (bins + numThreadsPerBlock -1)/numThreadsPerBlock ;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // allocate host and device memory
    size_t memSize = size * sizeof(float);
    cudaMalloc( (void **) &d_a, memSize );
    cudaMalloc( (void **) &d_b, memSize );

    // Initialize input array on host

    cudaEventRecord(start, 0);   // to start timing
    // Copy host array to device array
    cudaMemcpy( d_a, h_a, memSize, cudaMemcpyHostToDevice );



    // launch kernel
    dim3 dimGrid(numBlocks);
    dim3 dimBlock(numThreadsPerBlock);
 
    reverseArrayBlockNewFloat<<< dimGrid, dimBlock >>>( d_b, d_a, nchans, bins );

    // block until the device has completed
    cudaDeviceSynchronize();

    // check if kernel execution generated an error
    // Check for any CUDA errors
    checkCUDAError("kernel invocation");

    // device to host copy
    cudaMemcpy( h_a, d_b, memSize, cudaMemcpyDeviceToHost );

    // Check for any CUDA errors
    checkCUDAError("memcpy");

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);   // makes sure gpu is done, and is part of the timing module. but can synchronise in other ways instead 
    //if not insterested in timing
    cudaEventElapsedTime( &time_kernel, start, stop);
    //printf("Time for executing channel flip = %f msec\n", time_kernel);


    // free device memory
    cudaFree(d_a);
    cudaFree(d_b);

    return;
}

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }
}

