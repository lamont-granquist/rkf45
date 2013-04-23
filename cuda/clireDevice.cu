#include <stdio.h>

#include "clireDevice.hu"

void HandleError(cudaError_t error, const char *file, int line) {
  if (error != cudaSuccess) {
    printf("%s in %s at line %d\n", cudaGetErrorString(error), file, line);
    exit(EXIT_FAILURE);
  }
}

// Prepares the selected device. If device is negative, select the fastest device
void setupCuda(int device) {
  cudaDeviceProp  prop;
    
  cudaFree(0);
  
  int count;
  HANDLE_ERROR(cudaGetDeviceCount(&count));
  
  if (count == 0) {
    printf("No CUDA device found\nd");
    exit(-1);
  }
      
  if (device >= 0) {
    HANDLE_ERROR(cudaSetDevice(device));
  }
  else {
    int lastclock = 0;
        
    for (int i=0; i< count; i++) {
      HANDLE_ERROR(cudaGetDeviceProperties(&prop, i));
       
      // Choose the fastest card
      if (prop.clockRate > lastclock) {
        HANDLE_ERROR(cudaSetDevice(i));
        lastclock = prop.clockRate;
      }    
    }
  }

  int usedDevice;
  HANDLE_ERROR(cudaGetDevice(&usedDevice));
  HANDLE_ERROR(cudaGetDeviceProperties(&prop, usedDevice)); 
  HANDLE_ERROR(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
  
  printf("Device: %s\n", prop.name);
  size_t available, total;
  cudaMemGetInfo(&available, &total);
    
  printf("Device memory: %zuMB free of %zuMB total\n\n", 
         available / 1000000, total / 1000000);
}

// Given a number of jobs, try to choose sensible block and thread numbers
void distributeJobs(int jobs, int* blockResult, int* threadResult) {
  int device;
  cudaDeviceProp prop;
  HANDLE_ERROR(cudaGetDevice(&device));
  HANDLE_ERROR(cudaGetDeviceProperties(&prop, device));
    
  int multiProcessors = prop.multiProcessorCount;
  int warpSize        = prop.warpSize;
  int maxThreads      = prop.maxThreadsPerBlock;
     
  int blocks, threads;
    
  if (jobs <= multiProcessors) {
    blocks  = jobs;
    threads = 1;
  }
  else if (jobs <= multiProcessors * 2 * maxThreads) {
    blocks  = multiProcessors * 2;
    threads = ((jobs / multiProcessors / 2) / warpSize + 1) * 32;
        
    if (threads > 512) threads = 512; 
  }
  else {
    blocks  = multiProcessors * 2;
    threads = maxThreads;
  }
    
  *blockResult  = blocks;
  *threadResult = threads;
}
