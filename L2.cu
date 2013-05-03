//Main
#include <stdio.h>

int main(int argc, char const *argv[]) {
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  int device;
  for (device = 0; device < deviceCount; ++device) {
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, device);
      printf("Device %d has L2 cache size %d.\n",
             device, deviceProp.l2CacheSize);
      printf("Number of multiprocessors:     %d\n",  deviceProp.multiProcessorCount);
  }
}
