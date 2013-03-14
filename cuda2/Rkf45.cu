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
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a < _b ? _a : _b; })
#define sign(x)  ((x > 0) - ( x < 0))

void dy(float t, float* V,float* result);
void bj_ii(float t, float* result);
//Declare functions
static bool local_start_to_be_reached(float t,int local_start_year,float* stepsize);
static void calculate_solutions(int neqn,float t,float stepsize,float* y, float *y_diff,float* y_plus_one, float* y_plus_one_alternative);
static float calculate_solution_error(int neqn,float stepsize,float* y,float* y_plus_one, float* y_plus_one_alternative);
static void local_estimate(int neqn,int local_end_year,int local_start_year,float* stepsize,float* y,float* y_diff);
static float calculate_initial_stepsize(int neqn,int start_year,int end_year,float* y, float* y_diff);
static float scale_from_error(float error,bool stepsize_decreased);
//static float FindFloatEpsilon();

/********************** Solve *********************/

/* Calculate the actual and the alternative solutions */
//y_plus_one and y_plus_one_alternative will be set
static void calculate_solutions(int neqn,float t,float stepsize,float* y,float* y_diff,float* y_plus_one,float* y_plus_one_alternative) {

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
  dy ( t + lcd_stepsize, f_swap, f1 );

  //f2
  lcd_stepsize = 3.0f * stepsize / 32.0f;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( y_diff[i] + 3.0f * f1[i] );
  dy ( t + 3.0f * stepsize / 8.0f, f_swap, f2 );


  /*printf("f_swap!   :           %.7f\n",f_swap[0]);
  printf("!     :               %.7f\n",t + 3.0f * stepsize / 8.0f);*/

  //f3
  lcd_stepsize = stepsize / 2197.0f;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( 1932.0f * y_diff[i] + ( 7296.0f * f2[i] - 7200.0f * f1[i] ) );
  dy ( t + 12.0f * stepsize / 13.0f, f_swap, f3 );

  //f4
  lcd_stepsize = stepsize / 4104.0f;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( ( 8341.0f * y_diff[i] - 845.0f * f3[i] ) + 
        ( 29440.0f * f2[i] - 32832.0f * f1[i] ) );
  dy ( t + stepsize, f_swap, f4 );

  //f5
  lcd_stepsize = stepsize / 20520.0f;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + lcd_stepsize * ( ( -6080.0f * y_diff[i] + 
          ( 9295.0f * f3[i] - 5643.0f * f4[i] ) ) + ( 41040.0f * f1[i] - 28352.0f * f2[i] ) );
  dy ( t + stepsize / 2.0f, f_swap, f5 );

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
static void local_estimate(int neqn,int local_end_year,int local_start_year,float *stepsize,float* y,float* y_diff) {
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

    calculate_solutions(neqn,t,*stepsize,y,y_diff,y_plus_one,y_plus_one_alternative);
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
      calculate_solutions(neqn,t,*stepsize,y,y_diff,y_plus_one,y_plus_one_alternative);
      error = calculate_solution_error(neqn,*stepsize,y,y_plus_one,y_plus_one_alternative);
    }

    //Advance in time
    t = t + *stepsize; 

    //Apply solution
    for (int i = 0; i < neqn; i++ )
      y[i] = y_plus_one[i];

    //Update y_diff
    dy ( t, y, y_diff );

    //Apply scale to stepsize
    float scale = scale_from_error(error,stepsize_descresed);
    *stepsize = sign ( *stepsize ) * max ( scale * fabsf( *stepsize ), hmin );
  }
}

/**************** Local estimation help functions ****************/


/* React if the "local start year" is about to be reached */
//Effects stepsize, returns whether the start year is reached
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
        printf("this should not happen.\n");
      }
    }
  }

  return  max( s, 26.0f * FloatEpsilon * max( fabsf( end_year ), fabsf( start_year - end_year ) ) );
}

/* Scale from error calculations */
static float scale_from_error(float error,bool stepsize_decreased) {
  float scale = min(5.0f,0.9f / powf( error, 0.2f ));

  if (stepsize_decreased)
    scale = min( scale, 1.0f );

  return scale;
}

/*********************** Estimate **************************/

/* Estimate range */
void estimate(int neqn, int end_year, int start_year,float* y,float* result0) { //TODO: yy

  float y_diff[MAX_NEQN];
  dy((float) end_year, y, y_diff);
  float stepsize = calculate_initial_stepsize(neqn,start_year,end_year,y,y_diff); 

  //Solve for one year at a time
  for (int year=end_year; year>start_year; year--) {

    //Add this years benefit to y
    bj_ii(year,y);

    // Integrate over [year,year-1]
    local_estimate(neqn,year,year-1,&stepsize,y,y_diff);

    //Copy y to results
    result0[year-start_year-1] = y[0];
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

__global__ void test_kernel(CUSTOMERS *customers,float *result0) {
  int id = get_id();

  float y[MAX_NEQN];
  for(int i = 0;i<MAX_NEQN;i++)
    y[i] = 0.0f;


};

void cpu_kernel(CUSTOMERS *customers,float *result0) {
  int id = 0;

  float y[MAX_NEQN];
  for(int i = 0;i<MAX_NEQN;i++)
    y[i] = 0.0f;

  estimate(
           customers[id].neqn,
           customers[id].end_year,
           customers[id].start_year,
           y,
           result0
          );
}

/**************** RK_LIBRARY *****************/

float age = 30.0f;
float interestrate = 0.05f;
float bpension = 1.0f;
float pensiontime = 35.0f;

float GM(float t) {
    return 0.0005f + powf(10.0f, 5.728f - 10.0f + 0.038f*(age + t));
}

// Interest
float r(float t) {
    return interestrate;
}

float indicator(int b) {
    return b ? 1.0f : 0.0f;
}

/**************** PRODUCT, PURE ENDOWMENT ***************************/
static float b_0(float t) {
    return 0.0f;
}

static float mu_01(float t) {
    return GM(t);
}

static float bj_00(float t) {
    return t == pensiontime ? bpension: 0.0f;
}

static float bj_01(float t) {
    return 0.0f; 
}

void bj_ii(float t, float* result) {
  result[0] += bj_00(t);
}

void dy(float t, float* V,float* result)
{
    result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
