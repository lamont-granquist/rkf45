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

  dim3 block_dim(8,3,6); //Number of threads per block
  dim3 grid_dim(5,4,1);  //Number of blocks per grid (cc. 1.2 only supports 2d)
  int nsize = get_n_host(block_dim,grid_dim); 

  // Data on the host and the device, respectively
  float result[nsize];
  float result_cpu[51];
  float *dev_result;
  CUSTOMERS *dev_customers; 
  CUSTOMERS *customers;

  customers = (CUSTOMERS*) malloc(sizeof(CUSTOMERS)*nsize);

  srand(19); //seed

  for(int i = 0;i < nsize;i++) {
    customers[i].policy = 1+rand()%6;
    customers[i].neqn = 1;
    if (customers[i].policy >= 5) {
      customers[i].neqn = 2;
    }
    customers[i].age = 5 + rand()%30;
    customers[i].end_year = 50;
    customers[i].start_year = 0;
  }

  // Allocate memory on the device
  cudaMalloc((void**)&dev_customers, sizeof(CUSTOMERS) * nsize);
  cudaMalloc((void**)&dev_result, sizeof(float) * nsize);

  // Copy data to the device
  cudaMemcpy(dev_customers, customers, sizeof(CUSTOMERS) * nsize, cudaMemcpyHostToDevice);

  // Launch the kernel with 10 blocks, each with 1 thread
  test_kernel <<<grid_dim, block_dim>>>(dev_customers,dev_result); // GPU
  //cpu_kernel(customers,result_cpu); //CPU

  // Copy the result back from the device
  cudaMemcpy(result, dev_result, sizeof(float) * nsize, cudaMemcpyDeviceToHost);

  // Print the result
  for(int i = 0; i < nsize; i++) {
    printf("%i: %11.7f, policy: %i, age: %i \n",i, result[i],customers[i].policy,customers[i].age);
  }
  /*
  for(int i = 0; i < 51; i++) {
    printf("%i: %.7f\n",i, result_cpu[i]);
  }
  */

  printf("\n");

  cudaFree(dev_result);

  return 0;
}
