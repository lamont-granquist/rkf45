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

  dim3 block_dim(1,2,2); //Number of threads per block
  dim3 grid_dim(2,1,1);  //Number of blocks per grid (cc. 1.2 only supports 2d)
  int nsize = get_n_host(block_dim,grid_dim); 

  // Data on the host and the device, respectively
  float result[nsize];
  float result_cpu[51];
  float *dev_result;
  CUSTOMERS *dev_customers; 
  CUSTOMERS *customers;

  customers = (CUSTOMERS*) malloc(sizeof(CUSTOMERS)*nsize);

  for(int i = 0;i < nsize;i++) {
    customers[i].neqn = 1;
    customers[i].policy = 2;
    customers[i].end_year = 50;
    customers[i].start_year = 0;
  }

  customers[0].policy = 2;

  customers[1].policy = 2;
  customers[2].policy = 3;
  customers[3].policy = 4;
  customers[4].policy = 5;
  customers[5].policy = 6;

  // Allocate memory on the device
  cudaMalloc((void**)&dev_customers, sizeof(CUSTOMERS) * nsize);
  cudaMalloc((void**)&dev_result, sizeof(float) * nsize);

  // Copy data to the device
  cudaMemcpy(dev_customers, customers, sizeof(CUSTOMERS) * nsize, cudaMemcpyHostToDevice);

  // Launch the kernel with 10 blocks, each with 1 thread
  //test_kernel <<<grid_dim, block_dim>>>(dev_customers,dev_result); // GPU
  cpu_kernel(customers,result_cpu); //CPU

  // Copy the result back from the device
  cudaMemcpy(result, dev_result, sizeof(float) * nsize, cudaMemcpyDeviceToHost);

  // Print the result
  /*
  for(int i = 0; i < nsize; i++) {
    printf("%i: %.7f\n",i, result[i]);
  }
  */
  for(int i = 0; i < 51; i++) {
    printf("%i: %.7f\n",i, result_cpu[i]);
  }

  printf("\n");

  cudaFree(dev_result);

  return 0;
}
