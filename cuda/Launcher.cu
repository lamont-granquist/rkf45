#include <stdio.h>
#include "Rkf45.hu"
#include "Customers.hu"

/**************************** HOST ******************************/

//Calculate the number of kernels
int get_n_host(dim3 block_dim,dim3 grid_dim) {
  return grid_dim.x * grid_dim.y * grid_dim.z * block_dim.x * block_dim.y * block_dim.z; 
}

// Host code
int main(int argc, char const *argv[]) {

  dim3 block_dim(2,2,3); //Number of threads per block
  dim3 grid_dim(2,3,1);  //Number of blocks per grid (cc. 1.2 only supports 2d)
  int nsize = get_n_host(block_dim,grid_dim); 

  // Data on the host and the device, respectively
  int result[nsize];
  int *dev_result;
  CUSTOMERS *dev_customers; 
  CUSTOMERS *customers;

  customers = (CUSTOMERS*) malloc(sizeof(CUSTOMERS)*nsize);

  for(int i = 0;i < nsize;i++) {
    customers[i].neqn = 1;
    customers[i].policy = 1;
    customers[i].end_year = 40;
  }

  // Allocate memory on the device
  cudaMalloc((void**)&dev_customers, sizeof(CUSTOMERS) * nsize);
  cudaMalloc((void**)&dev_result, sizeof(int) * nsize);

  // Copy data to the device
  cudaMemcpy(dev_customers, customers, sizeof(CUSTOMERS) * nsize, cudaMemcpyHostToDevice);

  // Launch the kernel with 10 blocks, each with 1 thread
  //kernel <<<grid_dim, block_dim>>>(dev_customers,dev_result);
  test_kernel <<<grid_dim, block_dim>>>(dev_customers, dev_result);

  // Copy the result back from the device
  cudaMemcpy(result, dev_result, sizeof(int) * nsize, cudaMemcpyDeviceToHost);

  // Print the result: "0 5 20 45 80 125 180 245 320 405"
  for(int i = 0; i < nsize; i++) {
    printf("%d ", result[i]);
  }

  printf("\n");

  cudaFree(dev_result);

  return 0;
}
