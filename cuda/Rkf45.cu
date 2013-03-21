/*
 * C implementation of the Rkf45 algoritm.
 */

const int MAX_NEQN = 2;
const float relerr = 1e-7;
const float abserr = 1e-7;
const float FloatEpsilon = 0.00000011920928955078125000f; //TODO: Calculate this constant

/********************* INIT *******************/

//Library inclusion
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //TODO: Delete

#include "Customers.hu"
#include "Rkf45.hu"

//Max,min,sign functions
/*#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a < _b ? _a : _b; })*/
#define sign(x)  ((x > 0) - ( x < 0))

__device__ __host__
void dy(int policy,int age,float t, float* V,float* result);
__device__ __host__
void bj_ii(int policy, float t, float* result);
//Declare functions
__device__ __host__
static bool local_start_to_be_reached(float t,int local_start_year,float* stepsize);
__device__ __host__
static void calculate_solutions(int policy,int age, int neqn,float t,float stepsize,float* y, float *y_diff,float* y_plus_one, float* y_plus_one_alternative);
__device__ __host__
static float calculate_solution_error(int neqn,float stepsize,float* y,float* y_plus_one, float* y_plus_one_alternative);
__device__ __host__
static void local_estimate(int policy,int age, int neqn,int local_end_year,int local_start_year,float* stepsize,float* y,float* y_diff);
__device__ __host__
static float calculate_initial_stepsize(int neqn,int start_year,int end_year,float* y, float* y_diff);
__device__ __host__
static float scale_from_error(float error,bool stepsize_decreased);
//static float FindFloatEpsilon();

/********************** Solve *********************/

/* Calculate the actual and the alternative solutions */
//y_plus_one and y_plus_one_alternative will be set
__device__ __host__
static void calculate_solutions(int policy,int age, int neqn,float t,float stepsize,float* y,float* y_diff,float* y_plus_one,float* y_plus_one_alternative) {

  float f1[MAX_NEQN];
  float f2[MAX_NEQN];
  float f3[MAX_NEQN];
  float f4[MAX_NEQN];
  float f5[MAX_NEQN];
  float f_swap[MAX_NEQN];

  float lcd_stepsize = stepsize / 4.0f; //lowest common denominator of stepsize

  //f1
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * y_diff[i];
  dy (policy,age, t + lcd_stepsize, f_swap, f1 );

  //f2
  lcd_stepsize = 3.0f * stepsize / 32.0f;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( y_diff[i] + 3.0f * f1[i] );
  dy (policy,age, t + 3.0f * stepsize / 8.0f, f_swap, f2 );


  /*printf("f_swap!   :           %.7f\n",f_swap[0]);
  printf("!     :               %.7f\n",t + 3.0f * stepsize / 8.0f);*/

  //f3
  lcd_stepsize = stepsize / 2197.0f;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( 1932.0f * y_diff[i] + ( 7296.0f * f2[i] - 7200.0f * f1[i] ) );
  dy (policy,age, t + 12.0f * stepsize / 13.0f, f_swap, f3 );

  //f4
  lcd_stepsize = stepsize / 4104.0f;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( ( 8341.0f * y_diff[i] - 845.0f * f3[i] ) + 
        ( 29440.0f * f2[i] - 32832.0f * f1[i] ) );
  dy (policy,age, t + stepsize, f_swap, f4 );

  //f5
  lcd_stepsize = stepsize / 20520.0f;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( ( -6080.0f * y_diff[i] + 
          ( 9295.0f * f3[i] - 5643.0f * f4[i] ) ) + ( 41040.0f * f1[i] - 28352.0f * f2[i] ) );
  dy (policy,age, t + stepsize / 2.0f, f_swap, f5 );

  //Calculate solution
  lcd_stepsize = stepsize / 7618050.0f;
  for (int i = 0; i < neqn; i++ )
    y_plus_one[i] = y[i] + lcd_stepsize * ( ( 902880.0f * y_diff[i] + 
          ( 3855735.0f * f3[i] - 1371249.0f * f4[i] ) ) + ( 3953664.0f * f2[i] + 277020.0f * f5[i] ) );

  //Calculate alternative solution
  for (int i = 0; i < neqn; i++ )
    y_plus_one_alternative[i] = ( -2090.0f * y_diff[i] + ( 21970.0f * f3[i] - 15048.0f * f4[i] ) ) + ( 22528.0f * f2[i] - 27360.0f * f5[i] );

  /*if(f1[0] != 0) {
    printf("y     :           %.7f\n",y[0]);
    printf("y_diff:           %.7f\n",y_diff[0]);
    printf("f_swap:           %.7f\n",f_swap[0]);
    printf("t + lcd_stepsize: %.7f\n",t + lcd_stepsize);
    printf("t:                %.7f\n",t);
    printf("f1[0]:            %.7f\n",f1[0]);
    printf("f2[0]:            %.7f\n",f2[0]);
    printf("f3[0]:            %.7f\n",f3[0]);
    printf("f4[0]:            %.7f\n",f4[0]);
    exit(0);
  }*/
}

/* Calculate the error of the solution */
//Pure
__device__ __host__
static float calculate_solution_error(int neqn,float stepsize,float* y,float* y_plus_one, float* y_plus_one_alternative) {

  //Used in calculations
  float scale = 2.0f / relerr;

  //Calculate the biggest_difference
  float biggest_difference = 0.0f;
  for (int i = 0; i < neqn; i++ )
  {
    float et = fabsf( y[i] ) + fabsf( y_plus_one[i] ) + scale * abserr;
    float ee = fabsf( y_plus_one_alternative[i] );

    biggest_difference = max ( biggest_difference, ee / et );
  }
  
  //Return the error
  return fabsf( stepsize ) * biggest_difference * scale / 752400.0f;
}

/******************* Local estimation ***********************/

/* Move from current position to local_start_year, and update all values */
// Updates y, h
__device__ __host__
static void local_estimate(int policy,int age, int neqn,int local_end_year,int local_start_year,float *stepsize,float* y,float* y_diff) {
  float t = (float)local_end_year;
  
  //Step by step integration.
  bool local_start_reached = false;
  while (!local_start_reached)
  {
    //Variables used in calculations
    bool stepsize_descresed = false;
    float hmin = 26.0f * FloatEpsilon * fabsf( t );

    local_start_reached = local_start_to_be_reached(t,local_start_year,stepsize);

    float y_plus_one[MAX_NEQN];
    float y_plus_one_alternative[MAX_NEQN];

    calculate_solutions(policy,age,neqn,t,*stepsize,y,y_diff,y_plus_one,y_plus_one_alternative);
    float error = calculate_solution_error(neqn,*stepsize,y,y_plus_one,y_plus_one_alternative);

    //Integreate 1 step
    while(error > 1.0f)
    {
      stepsize_descresed = true;
      local_start_reached = false;

      //Scale down.
      float s = max(0.1f,0.9f / powf( error, 0.2f ));
      *stepsize = s * *stepsize;  

      //Try again.
      calculate_solutions(policy,age,neqn,t,*stepsize,y,y_diff,y_plus_one,y_plus_one_alternative);
      error = calculate_solution_error(neqn,*stepsize,y,y_plus_one,y_plus_one_alternative);
    }

    //Advance in time
    t = t + *stepsize; 

    //Apply solution
    for (int i = 0; i < neqn; i++ )
      y[i] = y_plus_one[i];

    //Update y_diff
    dy (policy,age, t, y, y_diff );

    //Apply scale to stepsize
    float scale = scale_from_error(error,stepsize_descresed);
    *stepsize = sign ( *stepsize ) * max ( scale * fabsf( *stepsize ), hmin );
  }
}

/**************** Local estimation help functions ****************/


/* React if the "local start year" is about to be reached */
//Effects stepsize, returns whether the start year is reached
__device__ __host__
static bool local_start_to_be_reached(float t,int local_start_year,float* stepsize) {
    float dt = local_start_year - t;
    if ( 2.0f * fabsf( *stepsize ) > fabsf( dt ) )
    {
      if ( fabsf( dt ) <= fabsf( *stepsize ) ) //Final step?
      {
        *stepsize = dt;                   //Let stepsize hit output point
        return true;
      }
      else
      {
        *stepsize = 0.5f * dt; // If not final step, set stepsize to be second final step. (evens out)
      }
    }
    return false;
}

/* Calculate stepsize's startvalue */
__device__ __host__
static float calculate_initial_stepsize(int neqn,int start_year,int end_year,float* y,float *y_diff)
{
  //Calculate the start value of stepsize
  float s = fabsf( start_year - end_year );


  for (int k = 0; k < neqn; k++ )
  {
    float tol = relerr * fabsf( y[k] ) + abserr;
    if ( 0.0f < tol )
    {
      float ypk = fabsf( y_diff[k] );
      if ( tol < ypk * powf( s, 5.0f ) )
      {
        s = powf( ( tol / ypk ), 0.2f );
        //printf("this should not happen.\n");
      }
    }
  }

  return  max( s, 26.0f * FloatEpsilon * max( fabsf( end_year ), fabsf( start_year - end_year ) ) );
}

/* Scale from error calculations */
__device__ __host__
static float scale_from_error(float error,bool stepsize_decreased) {
  float scale = min(5.0f,0.9f / powf( error, 0.2f ));

  if (stepsize_decreased)
    scale = min( scale, 1.0f );

  return scale;
}

/*********************** Estimate **************************/

/* Estimate range */
__device__ __host__
void estimate(int policy,int age, int neqn, int end_year, int start_year,float* y,float* result0, float* result1) {

  float y_diff[MAX_NEQN];
  dy(policy,age, (float) end_year, y, y_diff);
  float stepsize = calculate_initial_stepsize(neqn,start_year,end_year,y,y_diff); 

  //Solve for one year at a time
  for (int year=end_year; year>start_year; year--) {

    //Add this years benefit to y
    bj_ii(policy,year,y);

    // Integrate over [year,year-1]
    local_estimate(policy,age,neqn,year,year-1,&stepsize,y,y_diff);

    //Copy y to results
    result0[year-start_year-1] = y[0];
    result1[year-start_year-1] = y[1];
  }
}

/*************************** Auxiliaries ****************************/

/* Find float epsilon */
/*static float FindFloatEpsilon() {
  float r = 1.0f;
  while (1.0f < (1.0f + r))
    r = r / 2.0f;
  return 2.0f * r;
}*/

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

/***** DEVICE ******/

__global__
void gpu_kernel(CUSTOMERS *customers,float *result) {

  int id = get_id();

  float y[MAX_NEQN];
  for(int i = 0;i<MAX_NEQN;i++)
    y[i] = 0.0f;

  float result0[51];
  for(int i = 0;i<51;i++)
    result0[i] = 0.0f;

  float result1[51];
  for(int i = 0;i<51;i++)
    result1[i] = 0.0f;

  estimate(
           customers[id].policy,
           customers[id].age,
           customers[id].neqn,
           customers[id].end_year,
           customers[id].start_year,
           y,
           result0,
           result1
          );

  result[id] = result0[0];

}


void cpu_kernel(CUSTOMERS *customers,float *result) {
  int cpu_id = 0;
  
  float y[MAX_NEQN];
  for(int i = 0;i<MAX_NEQN;i++)
    y[i] = 0.0f;

  float result0[51];
  for(int i = 0;i<51;i++)
    result0[i] = 0.0f;

  float result1[51];
  for(int i = 0;i<51;i++)
    result1[i] = 0.0f;

  estimate(
           customers[cpu_id].policy,
           customers[cpu_id].age,
           customers[cpu_id].neqn,
           customers[cpu_id].end_year,
           customers[cpu_id].start_year,
           y,
           result0,
           result1
          );
  
  for(int i = 0;i<51;i++)
    result[i] = result0[i];
}

/**************** RK_LIBRARY *****************/

__device__ __host__
float GM(int age, float t) {
    return 0.0005f + powf(10.0f, 5.728f - 10.0f + 0.038f*((float)age + t));
}

// Interest
__device__ __host__
float r(float t) {
  float interestrate = 0.05f;
    return interestrate;
}

__device__ __host__
float indicator(int b) {
    return b ? 1.0f : 0.0f;
}

/**************** PRODUCT, PURE ENDOWMENT ***************************/
__device__ __host__
static float b_0_PureEndowment(float t) {
    return 0.0f;
}

__device__ __host__
static float mu_01_PureEndowment(int age,float t) {
    return GM(age,t);
}

__device__ __host__
static float bj_00_PureEndowment(float t) {
    float bpension = 1.0f;
    float pensiontime = 35.0f;
    return t == pensiontime ? bpension: 0.0f;
}

__device__ __host__
static float bj_01_PureEndowment(float t) {
    return 0.0f; 
}

__device__ __host__
void bj_ii_PureEndowment(float t, float* result) {
  result[0] += bj_00_PureEndowment(t);
}

__device__ __host__
void dy_PureEndowment(int age, float t, float* V,float* result)
{
    result[0] = r(t) * V[0] - b_0_PureEndowment(t) - mu_01_PureEndowment(age,t) * (0 - V[0] + bj_01_PureEndowment(t));
}

/**************** PRODUCT, DEFFEREDLIFEANNUITY ***************************/
__device__ __host__
static float b_0_DeferredTemporaryLifeAnnuity(float t) {
    int m = 35;
    int n = 10;
    float bpension = 1.0f;
    return bpension * indicator(t > m) * indicator(t < m + n);
}

__device__ __host__
static float mu_01_DeferredTemporaryLifeAnnuity(int age, float t) {
    return GM(age,t);
}

__device__ __host__
static float bj_00_DeferredTemporaryLifeAnnuity(float t) {
    return 0.0f;
}

__device__ __host__
static float bj_01_DeferredTemporaryLifeAnnuity(float t) {
    return 0.0f; 
}

__device__ __host__
void bj_ii_DeferredTemporaryLifeAnnuity(float t, float* result) {
  result[0] += bj_00_DeferredTemporaryLifeAnnuity(t);
}

__device__ __host__
void dy_DeferredTemporaryLifeAnnuity(int age, float t, float* V,float* result)
{
    result[0] = r(t) * V[0] - b_0_DeferredTemporaryLifeAnnuity(t) - mu_01_DeferredTemporaryLifeAnnuity(age,t) * (0 - V[0] + bj_01_DeferredTemporaryLifeAnnuity(t));
}

/**************** PRODUCT, TemporaryLifeAnnuityPremium ***************************/
__device__ __host__
static float b_0_TemporaryLifeAnnuityPremium(float t) {
    int n = 35;
    int bpremium = 1;
    return -bpremium * indicator(t >= 0) * indicator(t < n);
}

__device__ __host__
static float mu_01_TemporaryLifeAnnuityPremium(int age, float t) {
    return GM(age,t);
}

__device__ __host__
static float bj_00_TemporaryLifeAnnuityPremium(float t) {
    return 0.0f;
}

__device__ __host__
static float bj_01_TemporaryLifeAnnuityPremium(float t) {
    return 0.0f; 
}

__device__ __host__
void bj_ii_TemporaryLifeAnnuityPremium(float t, float* result) {
  result[0] += bj_00_TemporaryLifeAnnuityPremium(t);
}

__device__ __host__
void dy_TemporaryLifeAnnuityPremium(int age, float t, float* V,float* result)
{
    result[0] = r(t) * V[0] - b_0_TemporaryLifeAnnuityPremium(t) - mu_01_TemporaryLifeAnnuityPremium(age,t) * (0 - V[0] + bj_01_TemporaryLifeAnnuityPremium(t));
}

/**************** PRODUCT, TermInsurance ***************************/
__device__ __host__
static float b_0_TermInsurance(float t) {
    return 0.0f;
}

__device__ __host__
static float mu_01_TermInsurance(int age, float t) {
    return GM(age,t);
}

__device__ __host__
static float bj_00_TermInsurance(float t) {
    return 0.0f;
}

__device__ __host__
static float bj_01_TermInsurance(float t) {
    int bdeath = 1;
    int n = 35;
    return bdeath * indicator(t > 0) * indicator(t < n);
}

__device__ __host__
void bj_ii_TermInsurance(float t, float* result) {
  result[0] += bj_00_TermInsurance(t);
}

__device__ __host__
void dy_TermInsurance(int age, float t, float* V,float* result)
{
    result[0] = r(t) * V[0] - b_0_TermInsurance(t) - mu_01_TermInsurance(age,t) * (0 - V[0] + bj_01_TermInsurance(t));
}

/**************** PRODUCT, DisabilityAnnuity ***************************/
__device__ __host__
static float b_0_DisabilityAnnuity(float t) {
    return 0.0f;
}

__device__ __host__
static float b_1_DisabilityAnnuity(float t) {
  int n = 35;
  int bdisabled = 1;
  //return 0.0f;
  return bdisabled * indicator(t > 0) * indicator(t < n);
}

__device__ __host__
static float GM01_DisabilityAnnuity(int age, float t) {
  return 0.0006f + powf(10.0f, 4.71609f - 10.0f + 0.06f*((float)age + t));
}

__device__ __host__
static float GM02_DisabilityAnnuity(int age, float t) {
  return GM(age,t);
}

__device__ __host__
static float GM12_DisabilityAnnuity(int age, float t) {
  return GM(age,t);
}

__device__ __host__
static float mu_01_DisabilityAnnuity(int age, float t) {
    return GM01_DisabilityAnnuity(age,t);
}

__device__ __host__
static float mu_02_DisabilityAnnuity(int age, float t) {
    return GM02_DisabilityAnnuity(age,t);
}

__device__ __host__
static float mu_12_DisabilityAnnuity(int age, float t) {
    return GM12_DisabilityAnnuity(age,t);
}

__device__ __host__
static float bj_00_DisabilityAnnuity(float t) {
    return 0.0f;
}

__device__ __host__
static float bj_01_DisabilityAnnuity(float t) {
  //int n = 35;
  //int bdisabled = 1;
  //return bdisabled * indicator(t > 0) * indicator(t < n);
  return 0.0f;
}

__device__ __host__
static float bj_02_DisabilityAnnuity(float t) {
  return 0.0f;
}

__device__ __host__
static float bj_11_DisabilityAnnuity(float t) {
  return 0.0f;
}

__device__ __host__
static float bj_12_DisabilityAnnuity(float t) {
  return 0.0f;
}

__device__ __host__
void bj_ii_DisabilityAnnuity(float t, float* result) {
  result[0] += bj_00_DisabilityAnnuity(t);
  result[1] += bj_11_DisabilityAnnuity(t);
}

__device__ __host__
void dy_DisabilityAnnuity(int age, float t, float* V,float* result)
{
  result[0] = r(t) * V[0] - b_0_DisabilityAnnuity(t) - mu_01_DisabilityAnnuity(age,t) * (V[1] - V[0] + bj_01_DisabilityAnnuity(t)) - mu_02_DisabilityAnnuity(age,t) * (0 - V[0] + bj_02_DisabilityAnnuity(t));
  result[1] = r(t) * V[1] - b_1_DisabilityAnnuity(t) - mu_12_DisabilityAnnuity(age,t) * (0 - V[1] + bj_12_DisabilityAnnuity(t)); 
}

/**************** PRODUCT, DisabilityTermInsurance ***************************/
__device__ __host__
static float b_0_DisabilityTermInsurance(float t) {
    return 0.0f;
}

__device__ __host__
static float b_1_DisabilityTermInsurance(float t) {
    return 0.0f;
}

__device__ __host__
static float GM01_DisabilityTermInsurance(int age, float t) {
  return 0.0006f + powf(10.0f, 4.71609f - 10.0f + 0.06f*((float)age + t));
}

__device__ __host__
static float GM02_DisabilityTermInsurance(int age,float t) {
  return GM(age,t);
}

__device__ __host__
static float GM12_DisabilityTermInsurance(int age, float t) {
  return GM(age,t);
}

__device__ __host__
static float mu_01_DisabilityTermInsurance(int age, float t) {
    return GM01_DisabilityTermInsurance(age,t);
}

__device__ __host__
static float mu_02_DisabilityTermInsurance(int age, float t) {
    return GM02_DisabilityTermInsurance(age,t);
}

__device__ __host__
static float mu_12_DisabilityTermInsurance(int age, float t) {
    return GM12_DisabilityTermInsurance(age,t);
}

__device__ __host__
static float bj_00_DisabilityTermInsurance(float t) {
    return 0.0f;
}

__device__ __host__
static float bj_01_DisabilityTermInsurance(float t) {
  int n = 35;
  int bdisabled = 1;
  return bdisabled * indicator(t > 0) * indicator(t < n);
}

__device__ __host__
static float bj_02_DisabilityTermInsurance(float t) {
  return 0.0f;
}

__device__ __host__
static float bj_11_DisabilityTermInsurance(float t) {
  return 0.0f;
}

__device__ __host__
static float bj_12_DisabilityTermInsurance(float t) {
  return 0.0f;
}

__device__ __host__
void bj_ii_DisabilityTermInsurance(float t, float* result) {
  result[0] += bj_00_DisabilityTermInsurance(t);
  result[1] += bj_11_DisabilityTermInsurance(t);
}

__device__ __host__
void dy_DisabilityTermInsurance(int age, float t, float* V,float* result)
{
  result[0] = r(t) * V[0] - b_0_DisabilityTermInsurance(t) - mu_01_DisabilityTermInsurance(age,t) * (V[1] - V[0] + bj_01_DisabilityTermInsurance(t)) - mu_02_DisabilityTermInsurance(age,t) * (0 - V[0] + bj_02_DisabilityTermInsurance(t));
  result[1] = r(t) * V[1] - b_1_DisabilityTermInsurance(t) - mu_12_DisabilityTermInsurance(age,t) * (0 - V[1] + bj_12_DisabilityTermInsurance(t)); 
}

/**** Policy distributor ****/

__device__ __host__
void dy(int policy,int age, float t, float* V, float* result) {
  switch(policy)
  {
    case 1:
      dy_PureEndowment(age,t,V,result);
    break;
    case 2:
      dy_DeferredTemporaryLifeAnnuity(age,t,V,result);
    break;
    case 3:
      dy_TemporaryLifeAnnuityPremium(age,t,V,result);
    break;
    case 4:
      dy_TermInsurance(age,t,V,result);
    break;
    case 5:
      dy_DisabilityAnnuity(age,t,V,result);
    break; 
    case 6:
      dy_DisabilityTermInsurance(age,t,V,result);
    break; 
  };
}

__device__ __host__
void bj_ii(int policy, float t, float* result) {
  switch(policy)
  {
    case 1:
      bj_ii_PureEndowment(t,result);
    break;
    case 2:
      bj_ii_DeferredTemporaryLifeAnnuity(t,result);
    break;
    case 3:
      bj_ii_TemporaryLifeAnnuityPremium(t,result);
    break;
    case 4:
      bj_ii_TermInsurance(t,result);
    break;
    case 5:
      bj_ii_DisabilityAnnuity(t,result);
    break;
    case 6:
      bj_ii_DisabilityTermInsurance(t,result);
    break;

  };
}
