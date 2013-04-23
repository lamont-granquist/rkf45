#include <stdio.h>
typedef struct
{
  int *neqn;
  int *policy;
} CUSTOMERS;

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

// Device code
__global__ void kernel(CUSTOMERS customers, int *result) {
  int tid = get_id();

  result[tid] = customers.policy[tid]*customers.neqn[tid];
}

// Host code
int main(int argc, char const *argv[]) {

  /********** 0. SETUP **********/
  dim3 block_dim(2,1,1); //Number of threads per block
  dim3 grid_dim(2,3,1);  //Number of blocks per grid (cc. 1.2 only supports 2d)
  //Number of kernels:
  int nsize = grid_dim.x * grid_dim.y * grid_dim.z * block_dim.x * block_dim.y * block_dim.z; 

  /********** 1. MALLOC HOST  **********/
  // Data on the host and the device, respectively
  int result[nsize]; // host
  int neqn[nsize];
  int policy[nsize];

  // Fill the arrays on the host
  for(int i=0;i<nsize;i++) {
    neqn[i] = i;
    policy[i] = nsize-i;
  }
  
  /********** 2. MALLOC DEVICE  **********/
  // Allocate memory on the device
  int *dev_result;     // device
  int *dev_neqn;
  int *dev_policy;
  cudaMalloc((void**)&dev_neqn, sizeof(int) * nsize);
  cudaMalloc((void**)&dev_policy, sizeof(int) * nsize);
  cudaMalloc((void**)&dev_result, sizeof(int) * nsize);

  /********** 3. COPY HOST TO DEVICE  **********/
  // Copy data to the device
  cudaMemcpy(dev_neqn, neqn, sizeof(int) * nsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_policy, policy, sizeof(int) * nsize, cudaMemcpyHostToDevice);

  /********** 4. CUSTOMERS HOLDS POINTERS TO DEVICE **********/
  //Used to hold the pointers
  CUSTOMERS customers;
  customers.neqn = dev_neqn;
  customers.policy = dev_policy;

  /********** 5. LAUNCH WITH CUSTOMERS AND RESULT *********/
  // Launch the kernel with 10 blocks, each with 1 thread
  kernel <<<grid_dim, block_dim>>>(customers, dev_result);

  /********** 6. COPY RESULT FROM DEVICE TO HOST *********/
  // Copy the result back from the device
  cudaMemcpy(result, dev_result, sizeof(int) * nsize, cudaMemcpyDeviceToHost);

  /********** 7. PRINT HOST RESULT *********/
  // Print the result: "0 5 20 45 80 125 180 245 320 405"
  for(int i = 0; i < nsize; i++) {
    printf("%d ", result[i]);
  }

  printf("\n");

  /********** 8. FREE DEVICE MEMORY *********/
  cudaFree(dev_neqn);
  cudaFree(dev_policy);
  cudaFree(dev_result);

  return 0;
}
