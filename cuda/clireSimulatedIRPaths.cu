#include <cuda.h>
#include <curand_kernel.h>

#include "clireSimulatedIRPaths.hu"
#include "clireDevice.hu"

__device__ double kappa = 1.155618, theta = 0.02384924, sigma = 0.08425505;

// Simulates an interest rate path with start rate r0.
// CIR model
__device__ void interestPath(const double r0, const int years, 
                             curandState *localState, 
                             double *dev_yieldCurves, const int n_yc,int yc) {    
  double dt = 1.0f;
  double r = r0;
  
  for (int y = 0; y <= years; y++) {
    //dev_yieldCurves[(years-y)*n_yc + yc] = r;                           //<- standard

    //if (y!=0) { dev_yieldCurves[(years-y)*2*n_yc + yc*2] = r; }          //<- formula A
    //if (y!=years) { dev_yieldCurves[(years-y-1)*2*n_yc + yc*2+1] = r;}

    dev_yieldCurves[(years-y) + yc*(years+1)] = r;                           //<- formula B

    r += MAXZERO(kappa * (theta - r) * dt + sigma * sqrtf(r * dt) * curand_normal(localState));
  }
}

// Sets up randomness states for some number of future launces of another kernel
// To put this in a seperate kernel is recommended by NVidia
__global__ void randSetupKernel(curandState *state, const int n_yc, 
                                const long long int seed) {
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  // Each thread gets same seed, a different sequence number, no offset.
  // This insures the random state for each thread is independent.
  while(tid < n_yc) {
    curand_init(seed, tid, 0, &state[tid]);
        
    tid += blockDim.x * gridDim.x;
  }
}

// Returns an array of a number of simulated interest rate paths 
// with length years*stepsPerYear+1 using the randomness from curandStates.
// The IR paths are "backwards". The value at position 0 in the array is for 
// the last year, while the last position holds the value for time 0.
__global__ void irPathsKernel(const int n_yc, const int years, 
                              curandState *curandStates, 
                              double *dev_yieldCurves) {
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  while(tid < n_yc) {
    // Copy state to local memory for efficiency. Recommended
    curandState localState = curandStates[tid];
    
    // Simulate an interest rate path
    interestPath(0.02f, years,  &localState, 
                 dev_yieldCurves, n_yc,tid);

    // Copy state back to global memory 
    curandStates[tid] = localState;
        
    tid += blockDim.x * gridDim.x;
  }
}


// Sets up the generation of irPaths and writes to the pointer dev_yieldCurves
void generateIRPaths(const int n_yc, const int years, 
                     double **dev_yieldCurves, 
                     const long long int seed) {

  curandState *dev_curandStates;
  // Allocate space for prng states on device 
  HANDLE_ERROR(cudaMalloc((void**)&dev_curandStates, n_yc * sizeof(curandState)));
  HANDLE_ERROR(cudaMemset(dev_curandStates, 0.0f,    n_yc * sizeof(curandState)));
    
  // Get an intelligent number of threads and blocks when 
  // generating a number of paths
  int blocks, threads;
  distributeJobs(n_yc, &blocks, &threads);

  // Setup prng states 
  randSetupKernel<<<blocks, threads>>>(dev_curandStates, n_yc, seed);
  HANDLE_ERROR(cudaGetLastError());
    
  HANDLE_ERROR(cudaMalloc((void**)dev_yieldCurves, n_yc * (1+years) * sizeof(double) * 2));
  HANDLE_ERROR(cudaMemset(*dev_yieldCurves, 0.0f,  n_yc * (1+years) * sizeof(double) * 2));
  
  // Call the kernel
  irPathsKernel<<<blocks, threads>>>(n_yc, years, 
                                     dev_curandStates, *dev_yieldCurves);
  HANDLE_ERROR(cudaGetLastError());

  // Dispose
  HANDLE_ERROR(cudaFree(dev_curandStates));
}
