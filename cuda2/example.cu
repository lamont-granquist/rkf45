#include <stdio.h>

__device__ int get_id(void) {
  // Find the ID for this thread, based on which block it is in.
  int idx = threadIdx.x + blockIdx.x * blockDim.x; //thread x coordinate
  int idy = threadIdx.y + blockIdx.y * blockDim.y; //thread y coordinate
  int idz = threadIdx.z + blockIdx.z * blockDim.z; //thread z coordinate

  int size_1d = blockDim.x * gridDim.x;            //n.o. threads on x side
  int size_2d = size_1d * blockDim.y * gridDim.y;  //n.o. thread on x * y side

  return idx + idy * size_1d + idz * size_2d;    //unique id
}

__device__ int get_n(void) {
  return blockDim.x * blockDim.y * blockDim.z * gridDim.x * gridDim.y * gridDim.z;
}

//
__device__ int deviceFunction3(int* id,int* variable1,int* variable2,int* variable3,int* variable4) {
  *variable1 += 8;
  *variable4 += 7;
  *variable2 += 9;
  *variable3 += 2;

  return *variable1 + *variable2 + *variable3;
}

__device__ int deviceFunction2(int* id,int* variable1,int* variable2,int* variable3,int* variable4) {
  *variable3 += 8; 
  *variable1 += deviceFunction3(id,variable1,variable2,variable3,variable4);
  *variable4 += deviceFunction3(id,variable1,variable2,variable3,variable4);

  return *variable3 + *variable4;
}

__device__ int deviceFunction1(int* id,int* variable1,int* variable2,int* variable3,int* variable4) {
  *variable1 += *id;
  *variable4 += 2;
  *variable2 += deviceFunction2(id,variable1,variable2,variable3,variable4);
  *variable3 += *variable2 + *variable4;
  return *variable1 + *variable2 + *variable3 + *variable4;
}

// Kernel
__global__ void kernel(int *dev_a, int *dev_b, int *dev_c) {
  int id = get_id();
  int variable1 = 3;
  int variable2 = 5;
  int variable3 = 8;
  int variable4 = 8;

  dev_c[id] = deviceFunction1(&id,&variable1,&variable2,&variable3,&variable4);
}

// Host code
int main(int argc, char const *argv[]) {

  dim3 block_dim(2,1,1); //Number of threads per block
  dim3 grid_dim(2,3,1);  //Number of blocks per grid (cc. 1.2 only supports 2d)
  //Number of kernels:
  int nsize = grid_dim.x * grid_dim.y * grid_dim.z * block_dim.x * block_dim.y * block_dim.z; 

  // Data on the host and the device, respectively
  int a[nsize], b[nsize], c[nsize]; // host
  int *dev_a , *dev_b , *dev_c;     // device

  // Fill the arrays on the host
  for(int i = 0; i < nsize; i++) {
    a[i] = i * i;
    b[i] = 2;
  }

  // Allocate memory on the device
  cudaMalloc((void**)&dev_a, sizeof(int) * nsize);
  cudaMalloc((void**)&dev_b, sizeof(int) * nsize);
  cudaMalloc((void**)&dev_c, sizeof(int) * nsize);

  // Copy data to the device
  cudaMemcpy(dev_a, a, sizeof(int) * nsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b, b, sizeof(int) * nsize, cudaMemcpyHostToDevice);

  // Launch the kernel with 10 blocks, each with 1 thread
  kernel <<<grid_dim, block_dim>>>(dev_a, dev_b, dev_c);

  // Copy the result back from the device
  cudaMemcpy(c, dev_c, sizeof(int) * nsize, cudaMemcpyDeviceToHost);

  // Print the result: "0 5 20 45 80 125 180 245 320 405"
  for(int i = 0; i < nsize; i++) {
    printf("%d ", c[i]);
  }

  printf("\n");

  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);

  return 0;
}
