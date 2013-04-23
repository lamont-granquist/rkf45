#include <cuda.h>
#include <curand_kernel.h>

#include "clireSimulatedIRPaths.hu"
#include "clireDevice.hu"

__device__ float kappa = 1.155618f, theta = 0.02384924f, sigma = 0.08425505f;

// Simulates an interest rate path with start rate r0.
// CIR model
__device__ void interestPath(const float r0, const int years, 
                             const int stepsPerYear, curandState *localState, 
                             float *result, const int offsetInc) {    
  float dt = 1.0f / stepsPerYear;
  float r = r0;
  
  // The last position for this interest path in the array
  int offset = years * stepsPerYear * offsetInc;
  int totalSteps = stepsPerYear * years;
  
  // Start at time 0, store it at the end of the array. As the time increases,
  // move towards the start of the array. This is beacause we calculate the
  // reserves moving from future to present.
  for (int step = 0; step < totalSteps; step++) {
    result[offset] = r;
          
    //float deltaR = kappa * (theta - r) * dt + sigma * sqrtf(r * dt) * curand_normal(localState);
    //r += (deltaR + fabsf(deltaR)) / 2.0f;         
    //r += deltaR * (deltaR > 0);
    r += MAXZERO(kappa * (theta - r) * dt + sigma * sqrtf(r * dt) * curand_normal(localState));
          
    offset -= offsetInc;
  }
  //Remember when t = 0.0f;
  result[offset] = r;
}

// Sets up randomness states for some number of future launces of another kernel
// To put this in a seperate kernel is recommended by NVidia
__global__ void randSetupKernel(curandState *state, const int numberOfPaths, 
                                const long long int seed) {
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;

  // Each thread gets same seed, a different sequence number, no offset.
  // This insures the random state for each thread is independent.
  while(tid < numberOfPaths) {
    curand_init(seed, tid, 0, &state[tid]);
        
    tid += blockDim.x * gridDim.x;
  }
}

// Returns an array of a number of simulated interest rate paths 
// with length years*stepsPerYear+1 using the randomness from curandStates.
// The IR paths are "backwards". The value at position 0 in the array is for 
// the last year, while the last position holds the value for time 0.
__global__ void irPathsKernel(const int numberOfPaths, const int years, 
                              const int stepsPerYear, curandState *curandStates, 
                              float *dev_yieldCurves) {
  unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
  
  while(tid < numberOfPaths) {
    // Copy state to local memory for efficiency. Recommended
    curandState localState = curandStates[tid];
    
    // Simulate an interest rate path
    interestPath(0.02f, years, stepsPerYear, &localState, 
                 dev_yieldCurves + tid, numberOfPaths);

    // Copy state back to global memory 
    curandStates[tid] = localState;
        
    tid += blockDim.x * gridDim.x;
  }
}


// Sets up the generation of irPaths and writes to the pointer dev_yieldCurves
void generateIRPaths(const int numberOfPaths, const int years, 
                     const int stepsPerYear, float **dev_yieldCurves, 
                     const long long int seed) {
  // Size of a single ir path
  const int singleIrPathSize = (years * stepsPerYear + 1) * sizeof(float); 
  // Size of all paths
  const int irPathsSize = numberOfPaths * singleIrPathSize; 

  curandState *dev_curandStates;
  // Allocate space for prng states on device 
  HANDLE_ERROR(cudaMalloc((void**)&dev_curandStates, numberOfPaths * sizeof(curandState)));
  HANDLE_ERROR(cudaMemset(dev_curandStates, 0,       numberOfPaths * sizeof(curandState)));
    
  // Get an intelligent number of threads and blocks when 
  // generating a number of paths
  int blocks, threads;
  distributeJobs(numberOfPaths, &blocks, &threads);

  // Setup prng states 
  randSetupKernel<<<blocks, threads>>>(dev_curandStates, numberOfPaths, seed);
  HANDLE_ERROR(cudaGetLastError());
    
  HANDLE_ERROR(cudaMalloc((void**)dev_yieldCurves, irPathsSize));
  HANDLE_ERROR(cudaMemset(*dev_yieldCurves, 0,     irPathsSize));
  
  // Call the kernel
  irPathsKernel<<<blocks, threads>>>(numberOfPaths, years, stepsPerYear, 
                                     dev_curandStates, *dev_yieldCurves);
  HANDLE_ERROR(cudaGetLastError());

  // Dispose
  HANDLE_ERROR(cudaFree(dev_curandStates));
}
