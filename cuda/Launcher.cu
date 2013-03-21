#include <stdio.h>
#include "Rkf45.hu"
#include "Customers.hu"
#include <time.h>

/**************************** HOST ******************************/

//Calculate the number of kernels
int get_n_host(dim3 block_dim,dim3 grid_dim) {
  return grid_dim.x * grid_dim.y * grid_dim.z * block_dim.x * block_dim.y * block_dim.z; 
}

// Host code
int main(int argc, char const *argv[]) {

  /********** 0. SETUP **********/
  //dim3 block_dim(8,8,8); //Number of threads per block
  //dim3 grid_dim(64,24,1);  //Number of blocks per grid (cc. 1.2 only supports 2d)
  dim3 block_dim(2,2,1); //Number of threads per block
  dim3 grid_dim(2,1,1);  //Number of blocks per grid (cc. 1.2 only supports 2d)

  int nsize = get_n_host(block_dim,grid_dim); 

  /********** 1. MALLOC HOST  **********/
  // Data on the host and the device, respectively
  float result[nsize];
  int neqn[nsize];
  int policy[nsize];
  int age[nsize];
  int end_year[nsize];
  int start_year[nsize];

  srand(19); //seed
  for(int i = 0;i < nsize;i++) {
    policy[i] = 1+rand()%6;
    neqn[i] = 1;
    if (policy[i] >= 5) {
      neqn[i] = 2;
    }
    age[i] = 5 + rand()%30;
    end_year[i] = 50;
    start_year[i] = 0;
  }

  /********** 2. MALLOC DEVICE  **********/

  float *dev_result;
  int *dev_neqn;
  int *dev_policy;
  int *dev_age;
  int *dev_end_year;
  int *dev_start_year;
  cudaMalloc((void**)&dev_result, sizeof(float) * nsize);
  cudaMalloc((void**)&dev_neqn, sizeof(int) * nsize);
  cudaMalloc((void**)&dev_policy, sizeof(int) * nsize);
  cudaMalloc((void**)&dev_age, sizeof(int) * nsize);
  cudaMalloc((void**)&dev_end_year, sizeof(int) * nsize);
  cudaMalloc((void**)&dev_start_year, sizeof(int) * nsize);

  /********** 3. COPY HOST TO DEVICE  **********/
  // Copy data to the device
  cudaMemcpy(dev_neqn, neqn, sizeof(int) * nsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_policy, policy, sizeof(int) * nsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_age, age, sizeof(int) * nsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_end_year, end_year, sizeof(int) * nsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_start_year, start_year, sizeof(int) * nsize, cudaMemcpyHostToDevice);

  /********** 4. CUSTOMERS HOLDS POINTERS TO DEVICE **********/
  //Used to hold the pointers
  CUSTOMERS customers;
  customers.neqn = dev_neqn;
  customers.policy = dev_policy;
  customers.age = dev_age;
  customers.end_year = dev_end_year;
  customers.start_year = dev_start_year;

  //********* 5. TIMING START ************/
  //Normal timing
  clock_t start = clock();

  //Cuda timing
  cudaEvent_t cuda_start, cuda_stop;
  float cuda_time;
  cudaEventCreate(&cuda_start);
  cudaEventCreate(&cuda_stop);
  cudaEventRecord( cuda_start, 0 );

  /********** 6. LAUNCH WITH CUSTOMERS AND RESULT *********/
  gpu_kernel <<<grid_dim, block_dim>>>(customers,dev_result); // GPU
  //cpu_kernel(customers,result_cpu); //CPU

  /********** 7. TIMING ENDS *********/
  //Cuda timing
  cudaEventRecord( cuda_stop, 0 );
  cudaEventSynchronize( cuda_stop );
  cudaEventElapsedTime( &cuda_time, cuda_start, cuda_stop );
  cudaEventDestroy( cuda_start );
  cudaEventDestroy( cuda_stop );
  
  /********** 8. COPY RESULT FROM DEVICE TO HOST *********/
  // Copy the result back from the device
  cudaMemcpy(result, dev_result, sizeof(float) * nsize, cudaMemcpyDeviceToHost);

  /********** 8,5. EXTRA TIMING *********/
  //Normal timing
  clock_t end = clock();
  float time = (float) (end - start) * 1000.0f / CLOCKS_PER_SEC;

  /********** 9. PRINT HOST RESULT  *********/
  // Print the result
  for(int i = nsize-7; i < nsize; i++) {
    printf("%i: %11.7f, policy: %i, age: %i \n",i, result[i],policy[i],age[i]);
  }

  /*
  for(int i = 0; i < 51; i++) {
    printf("%i: %.7f\n",i, result_cpu[i]);
  }
  */

  printf("TIME: %f, CUDA_TIME: %f\n",time,cuda_time);

  /********** 10. FREE MEMORY   *********/
  cudaFree(dev_result);
  cudaFree(dev_policy);
  cudaFree(dev_neqn);
  cudaFree(dev_age);
  cudaFree(dev_end_year);
  cudaFree(dev_start_year);

  return 0;
}
