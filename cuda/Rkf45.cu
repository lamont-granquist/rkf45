/*
 * C implementation of the Rkf45 algoritm.
 */

const int MAX_NEQN = 2;
const float relerr = 1e-11;
const float abserr = 1e-11;
const float DoubleEpsilon = 0; //TODO: Calculate this constant
/********************* INIT *******************/
//Library inclusion
#include <stdlib.h>
#include <stdio.h>

__device__ void dy(float t, float* V,float* result);
__device__ void bj_ii(float t, float* result);
/*
#include <assert.h>
#include <math.h>
#include "Matrix_Library.h"
#include <unistd.h> //only for sleep?
#include <string.h> //only for sleep?

//Max,min,sign functions /REMOVED/

//Declare functions
bool is_equal();
static void test_all();
static void time_all();
static float time_one();
static bool local_start_to_be_reached();
static void xpy();
static void calculate_solutions();
static void local_estimate();
static float calculate_initial_stepsize();
static float scale_from_error();
static float FindDoubleEpsilon();
static void allocate_equation_space();

//Declare Estimator variables

//Public variables
int neqn;
int start_year;
int end_year;
// dy       //policy
// bj_ii
float relerr;
float abserr;
float* end_year_y; 

//Private variables
static float t;
static float stepsize;
static float* f1;
static float* f2;
static float* f3;
static float* f4;
static float* f5;
static float* f_swap;
static float* y;
static float* y_diff;
static float* y_plus_one;
static float* y_plus_one_alternative;
static int local_start_year;
static int local_end_year;
static float DoubleEpsilon;

static int m; //Result length;
*/
/******************* Constructor *********************/

/* Construct */
/*
void construct(int n) {
  //DoubleEpsilon = FindDoubleEpsilon();
  allocate_equation_space();
}*/

/* Allocate equation space */
/*
static void allocate_equation_space() {
  //Global for the class
  y_plus_one             = (float) malloc(sizeof(float)*neqn);
  end_year_y             = (float) malloc(sizeof(float)*neqn);
  y                      = (float) malloc(sizeof(float)*neqn);
  y_diff                 = (float) malloc(sizeof(float)*neqn);

  //Temporary for the solve method
  f1                     = malloc(sizeof(float)*neqn);
  f2                     = malloc(sizeof(float)*neqn);
  f3                     = malloc(sizeof(float)*neqn);
  f4                     = malloc(sizeof(float)*neqn);
  f5                     = malloc(sizeof(float)*neqn);
  f_swap                 = malloc(sizeof(float)*neqn);
  y_plus_one_alternative = malloc(sizeof(float)*neqn);
}*/

/********************** Solve *********************/

/* Calculate the actual and the alternative solutions */
//y_plus_one and y_plus_one_alternative will be set
/*static void calculate_solutions() {

  float lcd_stepsize = stepsize / 4.0; //lowest common denominator of stepsize

  //f1
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * y_diff[i];
  dy ( t + lcd_stepsize, f_swap, f1 );

  //f2
  lcd_stepsize = 3.0 * stepsize / 32.0;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( y_diff[i] + 3.0 * f1[i] );
  dy ( t + 3.0 * stepsize / 8.0, f_swap, f2 );

  //f3
  lcd_stepsize = stepsize / 2197.0;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( 1932.0 * y_diff[i] + ( 7296.0 * f2[i] - 7200.0 * f1[i] ) );
  dy ( t + 12.0 * stepsize / 13.0, f_swap, f3 );

  //f4
  lcd_stepsize = stepsize / 4104.0;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( ( 8341.0 * y_diff[i] - 845.0 * f3[i] ) + 
        ( 29440.0 * f2[i] - 32832.0 * f1[i] ) );
  dy ( t + stepsize, f_swap, f4 );

  //f5
  lcd_stepsize = stepsize / 20520.0;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( ( -6080.0 * y_diff[i] + 
          ( 9295.0 * f3[i] - 5643.0 * f4[i] ) ) + ( 41040.0 * f1[i] - 28352.0 * f2[i] ) );
  dy ( t + stepsize / 2.0, f_swap, f5 );

  //Calculate solution
  lcd_stepsize = stepsize / 7618050.0;
  for (int i = 0; i < neqn; i++ )
    y_plus_one[i] = y[i] + lcd_stepsize * ( ( 902880.0 * y_diff[i] + 
          ( 3855735.0 * f3[i] - 1371249.0 * f4[i] ) ) + ( 3953664.0 * f2[i] + 277020.0 * f5[i] ) );

  //Calculate alternative solution
  for (int i = 0; i < neqn; i++ )
    y_plus_one_alternative[i] = ( -2090.0 * y_diff[i] + ( 21970.0 * f3[i] - 15048.0 * f4[i] ) ) + ( 22528.0 * f2[i] - 27360.0 * f5[i] );

}
*/
/* Calculate the error of the solution */
//Pure
/*static float calculate_solution_error() {

  //Used in calculations
  float scale = 2.0 / relerr;

  //Calculate the biggest_difference
  float biggest_difference = 0.0;
  for (int i = 0; i < neqn; i++ )
  {
    float et = fabs( y[i] ) + fabs( y_plus_one[i] ) + scale * abserr;
    float ee = fabs( y_plus_one_alternative[i] );

    biggest_difference = max ( biggest_difference, ee / et );
  }
  
  //Return the error
  return fabs( stepsize ) * biggest_difference * scale / 752400.0;
}
*/
/******************* Local estimation ***********************/

/* Move from current position to local_start_year, and update all values */
// Updates y, h
__device__ void local_estimate(float local_end_year,float local_start_year) {
  float t = local_end_year;
  
  //Step by step integration.
  bool local_start_reached = false;

  while (!local_start_reached)
  {
    //Variables used in calculations
    bool stepsize_descresed = false;
    float hmin = 26.0f * DoubleEpsilon * fabs( t );

    local_start_reached = local_start_to_be_reached();

    /*
    calculate_solutions();
    float error = calculate_solution_error();

    //Integreate 1 step
    while(error > 1.0f)
    {
      stepsize_descresed = true;
      local_start_reached = false;

      //Scale down.
      float s = max(0.1f,0.9f / pow( error, 0.2f ));
      stepsize = s * stepsize;  

      //Try again.
      calculate_solutions();
      error = calculate_solution_error();
    }

    //Advance in time
    t = t + stepsize; 

    //Apply solution
    for (int i = 0; i < neqn; i++ )
      y[i] = y_plus_one[i];

    //Update y_diff
    dy ( t, y, y_diff );

    //Apply scale to stepsize
    float scale = scale_from_error(error,stepsize_descresed);
    stepsize = sign ( stepsize ) * max ( scale * fabs( stepsize ), hmin );
    */
  }
}

/**************** Local estimation help functions ****************/


/* React if the "local start year" is about to be reached */
//Effects stepsize, returns whether the start year is reached
/*static bool local_start_to_be_reached() {
    float dt = local_end_year - t;
    if ( 2.0 * fabs( stepsize ) > fabs( dt ) )
    {
      if ( fabs( dt ) <= fabs( stepsize ) ) //Final step?
      {
        stepsize = dt;                   //Let stepsize hit output point
        return true;
      }
      else
      {
        stepsize = 0.5 * dt; // If not final step, set stepsize to be second final step. (evens out)
      }
    }
    return false;
}*/

/* Calculate stepsize's startvalue */
__device__ float calculate_initial_stepsize(int neqn,float* y, float* y_diff,float t)
{
  int start_year = 0; //NOTE: Assumption

  // s = stepsize, local var
  float s = abs((float)start_year - t);

  for (int k = 0; k < neqn; k++ )
  {
    float tol = relerr * fabs( y[k] ) + abserr;
    if ( 0.0f < tol )
    {
      float ypk = abs(y_diff[k]);
      if ( tol < ypk * __powf( s, 5 ) )
      {
        s = __powf( ( tol / ypk ), 0.2f );
        //printf("this should not happen.\n");//
      }
    }
  }

  return max( s, 26.0f * DoubleEpsilon * max( fabs( t ), fabs(start_year - t)));
}

/* Scale from error calculations */
/*
static float scale_from_error(float error,bool stepsize_decreased) {
  float scale = min(5.0,0.9 / pow( error, 0.2 ));

  if (stepsize_decreased)
    scale = min( scale, 1.0 );

  return scale;
}
*/

/*********************** Estimate **************************/

/* Estimate range */
__device__ void estimate(int neqn,int policy,int end_year,float *y) {

  float y_diff[MAX_NEQN];

  dy((float) end_year, y, y_diff);         // y_diff

  float stepsize = calculate_initial_stepsize(neqn,y,y_diff,(float) end_year); // stepsize

  //Solve for one year at a time
  for (int year=end_year; year>0; year--) { //start_year = 0;
    //Add this years benefit to y
    bj_ii(year,y);

    local_estimate((float)year,(float)year-1);
  }
}

/*************************** Auxiliaries ****************************/

/* Find float epsilon */
/*static float FindDoubleEpsilon() {
  float r = 1.0;
  while (1.0 < (1.0 + r))
    r = r / 2.0;
  return 2.0 * r;
}*/

/**************************** DEVICE ******************************/

#include "Customers.hu"
#include "Rkf45.hu"

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
__global__ void test_kernel(CUSTOMERS *customers,float *result) {
  int id = get_id();

  float y[MAX_NEQN];
  for(int i = 0;i<MAX_NEQN;i++)
    y[i] = 0.0f;

  estimate(
           customers[id].neqn,
           customers[id].policy,
           customers[id].end_year,
           y
          );

  result[id] = y[0];
}

/**************** RK_LIBRARY *****************/

__device__ float age = 30.0f;
__device__ float interestrate = 0.05f;
__device__ float bpension = 1.0f;
__device__ float pensiontime = 35.0f;

__device__ float GM(float t) {
    return 0.0005f + __powf(10.0f, 5.728f - 10.0f + 0.038f*(age + t));
}

// Interest
__device__ float r(float t) {
    return interestrate;
}

__device__ float indicator(int b) {
    return b ? 1.0f : 0.0f;
}

/**************** PRODUCT, PURE ENDOWMENT ***************************/
__device__ static float b_0(float t) {
    return 0.0f;
}

__device__ static float mu_01(float t) {
    return GM(t);
}

__device__ static float bj_00(float t) {
    return t == pensiontime ? bpension: 0.0f;
}

__device__ static float bj_01(float t) {
    return 0.0f; 
}

__device__ void bj_ii(float t, float* result) {
  result[0] += bj_00(t);
}

__device__ void dy(float t, float* V,float* result)
{
    result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
