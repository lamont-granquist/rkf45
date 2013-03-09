#include <stdio.h>

typedef struct
{
    int a, b;
} point;

/**************************** DEVICE ******************************/

//Calculate the id
__device__ int get_id(void) {
  // Find the ID for this thread, based on which block it is in.
  int idx = threadIdx.x + blockIdx.x * blockDim.x; //thread x coordinate
  int idy = threadIdx.y + blockIdx.y * blockDim.y; //thread y coordinate
  int idz = threadIdx.z + blockIdx.z * blockDim.z; //thread z coordinate

  int size_1d = blockDim.x * gridDim.x;            //n.o. threads on x side
  int size_2d = size_1d * blockDim.y * gridDim.y;  //n.o. thread on x * y side

  return idx + idy * size_1d + idz * size_2d;    //unique id
}

//Calculate the number of kernels
__device__ int get_n_device(void) {
  return blockDim.x * blockDim.y * blockDim.z * gridDim.x * gridDim.y * gridDim.z;
}

// Device code
__global__ void kernel(int *customer,point *myPoint, int *result) {
  int id = get_id();

  //construct(2);

  // Use the thread ID as an index
  result[id] = myPoint->a; //neqn; //get_n();//idz; //dev_a[tid] + dev_b[tid];
}

/**************************** HOST ******************************/

//Calculate the number of kernels
int get_n_host(dim3 block_dim,dim3 grid_dim) {
  return grid_dim.x * grid_dim.y * grid_dim.z * block_dim.x * block_dim.y * block_dim.z; 
}

// Host code
int main(int argc, char const *argv[]) {

  //HOST:
  //construct (decide neqn)

  //Set public vars.
  //(Must be viewable from device)

  //Start kernel
  //(estimate);

  dim3 block_dim(2,2,3); //Number of threads per block
  dim3 grid_dim(2,3,1);  //Number of blocks per grid (cc. 1.2 only supports 2d)
  int nsize = get_n_host(block_dim,grid_dim); 

  // Data on the host and the device, respectively
  int customer[nsize], result[nsize]; // host
  int *dev_customer , *dev_result;    // device pointers

  point *dev_myPoint; 

  point *myPoint = (point*) malloc(sizeof(point));
  myPoint->a = 9;
  myPoint->b = 4;
  // Fill the arrays on the host
  for(int i = 0; i < nsize; i++) {
    customer[i] = i * i;
  }
  
  // Allocate memory on the device
  cudaMalloc((void**)&dev_customer, sizeof(int) * nsize);
  cudaMalloc((void**)&dev_myPoint, sizeof(point));
  cudaMalloc((void**)&dev_result, sizeof(int) * nsize);

  // Copy data to the device
  cudaMemcpy(dev_customer, customer, sizeof(int) * nsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_myPoint, myPoint, sizeof(point), cudaMemcpyHostToDevice);

  // Launch the kernel with 10 blocks, each with 1 thread
  kernel <<<grid_dim, block_dim>>>(dev_customer,dev_myPoint,dev_result);

  // Copy the result back from the device
  cudaMemcpy(result, dev_result, sizeof(int) * nsize, cudaMemcpyDeviceToHost);

  // Print the result: "0 5 20 45 80 125 180 245 320 405"
  for(int i = 0; i < nsize; i++) {
    printf("%d ", result[i]);
  }

  printf("\n");

  cudaFree(dev_customer);
  cudaFree(dev_result);

  return 0;
}
