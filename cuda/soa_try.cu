#include <stdio.h>
typedef struct
{
  int *neqn;
} CUSTOMERS;

const int MAX_KERNELS = 72;

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

__device__ int id;
__device__ int gpu_array[MAX_KERNELS];

// Device code
__global__ void kernel(CUSTOMERS customers, int *dev_b, int *result) {
  int tid = get_id();

  result[tid] = customers.neqn[tid];//id;
}

const int N = 12;

// Host code
int main(int argc, char const *argv[]) {

  dim3 block_dim(2,1,1); //Number of threads per block
  dim3 grid_dim(2,3,1);  //Number of blocks per grid (cc. 1.2 only supports 2d)
  //Number of kernels:
  int nsize = grid_dim.x * grid_dim.y * grid_dim.z * block_dim.x * block_dim.y * block_dim.z; 

  CUSTOMERS customers;

  // Data on the host and the device, respectively
  int b[nsize], c[nsize]; // host
  int neqn[N];
  int *dev_b , *result;     // device
  int *dev_neqn;

  // Fill the arrays on the host
  for(int i = 0; i < nsize; i++) {
    b[i] = 2;
  }

  for(int i=0;i<N;i++) {
    neqn[i] = i;
  }
  
  // Allocate memory on the device
  cudaMalloc((void**)&dev_neqn, sizeof(int) * N);
  cudaMalloc((void**)&dev_b, sizeof(int) * nsize);
  cudaMalloc((void**)&result, sizeof(int) * nsize);

  // Copy data to the device
  cudaMemcpy(dev_neqn, neqn, sizeof(int) * nsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b, b, sizeof(int) * nsize, cudaMemcpyHostToDevice);

  //Point device pointer in host struct
  customers.neqn = dev_neqn;

  // Launch the kernel with 10 blocks, each with 1 thread
  kernel <<<grid_dim, block_dim>>>(customers, dev_b, result);

  // Copy the result back from the device
  cudaMemcpy(c, result, sizeof(int) * nsize, cudaMemcpyDeviceToHost);

  // Print the result: "0 5 20 45 80 125 180 245 320 405"
  for(int i = 0; i < nsize; i++) {
    printf("%d ", c[i]);
  }

  printf("\n");

  cudaFree(dev_neqn);
  cudaFree(dev_b);
  cudaFree(result);

  return 0;
}
