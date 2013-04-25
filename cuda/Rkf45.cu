/*
 * C implementation of the Rkf45 algoritm.
 */

__device__
const int MAX_NEQN = 2;
__device__
const double relerr = 1e-11;
__device__
const double abserr = 1e-11;
//__device__
const double DoubleEpsilon = 0.00000011920928955078125000; 

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

__device__ 
void dy(int policy,int age,double t, double* V,double* result);
__device__ 
void bj_ii(int policy, double t, double* result);
//Declare functions
__device__ 
static bool local_start_to_be_reached(double t,int local_start_year,double* stepsize);
__device__ 
static void calculate_solutions(int policy,int age, int neqn,double t,double stepsize,double* y, double *y_diff,double* y_plus_one, double* y_plus_one_alternative);
__device__ 
static double calculate_solution_error(int neqn,double stepsize,double* y,double* y_plus_one, double* y_plus_one_alternative);
__device__ 
static void local_estimate(int policy,int age, int neqn,int local_end_year,int local_start_year,double* stepsize,double* y,double* y_diff);
__device__ 
static double calculate_initial_stepsize(int neqn,int start_year,int end_year,double* y, double* y_diff);
__device__ 
static double scale_from_error(double error,bool stepsize_decreased);
//static double FindDoubleEpsilon();

/********************** Solve *********************/

/* Calculate the actual and the alternative solutions */
//y_plus_one and y_plus_one_alternative will be set
__device__ 
static void calculate_solutions(int policy,int age, int neqn,double t,double h,double* y,double* y_diff,double* y_plus_one,double* y_plus_one_alternative) {

  double f1[MAX_NEQN];
  double f2[MAX_NEQN];
  double f3[MAX_NEQN];
  double f4[MAX_NEQN];
  double f5[MAX_NEQN];
  double f_swap[MAX_NEQN];

  //f1
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + h*(1.0/4.0)*y_diff[i];
  dy (policy,age, t + (1.0/4.0)*h, f_swap, f1 );

  //f2
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + h*(3.0/32.0)*y_diff[i]
                     + h*(9.0/32.0)*f1[i];
  dy (policy,age, t + (3.0/8.0)*h, f_swap, f2 );

  //f3
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + h*(1932.0/2197.0)* y_diff[i]
                     - h*(7200.0/2197.0)* f1[i]
                     + h*(7296.0/2197.0)* f2[i];
  dy (policy,age, t + (12.0/13.0)* h, f_swap, f3 );

  //f4
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + h*(439.0/216.0) * y_diff[i]
                     - h*(8.0/1.0)     * f1[i]
                     + h*(3680.0/513.0)* f2[i]
                     - h*(845.0/4104.0)* f3[i];
  dy (policy,age, t + h, f_swap, f4 );

  //f5
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] - h*(8.0/27.0) * y_diff[i]
                     + h*(2.0/1.0) * f1[i]
                     - h*(3544.0/2565.0) * f2[i]
                     + h*(1859.0/4104.0) * f3[i]
                     - h*(11.0/40.0) * f4[i];
  dy (policy,age, t + h*(1.0/2.0), f_swap, f5 );

  //Calculate solution
  for (int i = 0; i < neqn; i++ )
    y_plus_one[i] = y[i] + h*(16.0/135.0) * y_diff[i]
                         + h*(6656.0/12825.0) * f2[i]
                         + h*(28561.0/56430.0) * f3[i]
                         - h*(9.0/50.0) * f4[i]
                         + h*(2.0/55.0) * f5[i];

  //Calculate alternative (worse) solution
  for (int i = 0; i < neqn; i++ )
    y_plus_one_alternative[i] =  - 2090.0 * y_diff[i]
                                 + 22528.0 * f2[i]
                                 + 21970.0 * f3[i]
                                 - 15048.0 * f4[i]
                                 - 27360.0 * f5[i];

}

/* Calculate the error of the solution */
//Pure
__device__ 
static double calculate_solution_error(int neqn,double stepsize,double* y,double* y_plus_one, double* y_plus_one_alternative) {

  //Used in calculations
  double scale = 2.0 / relerr;

  //Calculate the biggest_difference
  double biggest_difference = 0.0;
  for (int i = 0; i < neqn; i++ )
  {
    double et = fabs( y[i] ) + fabs( y_plus_one[i] ) + scale * abserr;
    double ee = fabs( y_plus_one_alternative[i] );

    biggest_difference = max ( biggest_difference, ee / et );
  }
  
  //Return the error
  return fabs( stepsize ) * biggest_difference * scale / 752400.0;
}

/******************* Local estimation ***********************/

/* Move from current position to local_start_year, and update all values */
// Updates y, h
__device__ 
static void local_estimate(int policy,int age, int neqn,int local_end_year,int local_start_year,double *stepsize,double* y,double* y_diff) {
  double t = (double)local_end_year;
  
  //Step by step integration.
  bool local_start_reached = false;
  while (!local_start_reached)
  {
    //Variables used in calculations
    bool stepsize_descresed = false;
    double hmin = 26.0 * DoubleEpsilon * fabs( t );

    local_start_reached = local_start_to_be_reached(t,local_start_year,stepsize);

    double y_plus_one[MAX_NEQN];
    double y_plus_one_alternative[MAX_NEQN];

    calculate_solutions(policy,age,neqn,t,*stepsize,y,y_diff,y_plus_one,y_plus_one_alternative);
    double error = calculate_solution_error(neqn,*stepsize,y,y_plus_one,y_plus_one_alternative);

    //Integreate 1 step
    while(error > 1.0)
    {
      stepsize_descresed = true;
      local_start_reached = false;

      //Scale down.
      double s = max(0.1,0.9 / __pow( error, 0.2 ));
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
    double scale = scale_from_error(error,stepsize_descresed);
    *stepsize = sign ( *stepsize ) * max ( scale * fabs( *stepsize ), hmin );
  }
}

/**************** Local estimation help functions ****************/


/* React if the "local start year" is about to be reached */
//Effects stepsize, returns whether the start year is reached
__device__ 
static bool local_start_to_be_reached(double t,int local_start_year,double* stepsize) {
    double dt = local_start_year - t;
    if ( 2.0 * fabs( *stepsize ) > fabs( dt ) )
    {
      if ( fabs( dt ) <= fabs( *stepsize ) ) //Final step?
      {
        *stepsize = dt;                   //Let stepsize hit output point
        return true;
      }
      else
      {
        *stepsize = 0.5 * dt; // If not final step, set stepsize to be second final step. (evens out)
      }
    }
    return false;
}

/* Calculate stepsize's startvalue */
__device__ 
static double calculate_initial_stepsize(int neqn,int start_year,int end_year,double* y,double *y_diff)
{
  //Calculate the start value of stepsize
  double s = fabs( start_year - end_year );


  for (int k = 0; k < neqn; k++ )
  {
    double tol = relerr * fabs( y[k] ) + abserr;
    if ( 0.0 < tol )
    {
      double ypk = fabs( y_diff[k] );
      if ( tol < ypk * __pow( s, 5.0 ) )
      {
        s = __pow( ( tol / ypk ), 0.2 );
        //printf("this should not happen.\n");
      }
    }
  }

  return  max( s, 26.0 * DoubleEpsilon * max( fabs( end_year ), fabs( start_year - end_year ) ) );
}

/* Scale from error calculations */
__device__ 
static double scale_from_error(double error,bool stepsize_decreased) {
  double scale = min(5.0,0.9 / __pow( error, 0.2 ));

  if (stepsize_decreased)
    scale = min( scale, 1.0 );

  return scale;
}

/*********************** Estimate **************************/

/* Estimate range */
__device__ 
void estimate(int policy,int age, int neqn, int end_year, int start_year,double* y,double* result0, double* result1) {

  double y_diff[MAX_NEQN];
  dy(policy,age, (double) end_year, y, y_diff);
  double stepsize = calculate_initial_stepsize(neqn,start_year,end_year,y,y_diff); 

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

/* Find double epsilon */
/*static double FindDoubleEpsilon() {
  double r = 1.0;
  while (1.0 < (1.0 + r))
    r = r / 2.0;
  return 2.0 * r;
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
void gpu_kernel(int offset, CUSTOMERS customers,double *result,double *dev_yieldCurves) {

  int id = get_id()+offset;

  double y[MAX_NEQN];
  for(int i = 0;i<MAX_NEQN;i++)
    y[i] = 0.0;

  double result0[51];
  for(int i = 0;i<51;i++)
    result0[i] = 0.0;

  double result1[51];
  for(int i = 0;i<51;i++)
    result1[i] = 0.0;

  //To make sure of loading time of variables.
  int c_policy = customers.policy[id];
  int c_age = customers.age[id];
  int c_neqn = customers.neqn[id];
  int c_end_year = customers.end_year[id];
  int c_start_year = customers.start_year[id];

  estimate(
           c_policy,
           c_age,
           c_neqn,
           c_end_year,
           c_start_year,
           y,
           result0,
           result1
          );

  result[id] = result0[0];

}

// The Danish FSA yield curve (Finanstilsynets rentekurve).
// Data from 2011-11-16 

__device__
const float ts[] = { 
    0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 
    15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0
};

__device__
const float rs[] = { 
    1.146677033, 1.146677033, 1.146677033, 1.340669678, 1.571952911, 1.803236144, 
    2.034519377, 2.265802610, 2.497085843, 2.584085843, 2.710085843, 2.805085843, 
    2.871485843, 2.937885843, 3.004285843, 3.070685843, 3.137085843, 3.136485843, 
    3.135885843, 3.135285843, 3.134685843, 3.134085843, 3.113185843, 3.092285843, 
    3.071385843, 3.050485843, 3.029585843, 3.008685843, 2.987785843, 2.966885843, 
    2.945985843, 2.925085843
};

__device__
double interpolate(double t,int m) {
  double tm = ts[m], tm1 = ts[m+1];
  double rm = rs[m] / 100.0, rm1 = rs[m+1] / 100.0;
  double Rt = (rm * (tm1 - t) + rm1 * (t - tm)) / (tm1 - tm);
  return __logf(1.0 + Rt) + t / (tm1 - tm) * (rm1 - rm) / (1.0 + Rt);
}

__device__
double rFsa(double t) { 
    // Extrapolation:
    if (t <= 0.25) {
        //return __logf(1.0 + rs[0]/100.0);
    }
    if (t >= 30.0) {
        //return __logf(1.0 + rs[31]/100.0);
    }
    int m = 1 + floor(t);
    // Interpolation:
    if (t >= 0.25 && t < 0.5) {
        m = 0;
    }
    if (t >= 0.5 && t < 1.0) {
        m = 1;
    }
    //General case:
    return 0.05;//interpolate(t,m);
}

/**************** RK_LIBRARY *****************/

__device__ 
double GM(int age, double t) {
    return 0.0005 + __pow(10.0, 5.728 - 10.0 + 0.038*((double)age + t));
}

// Interest
__device__ 
double r(double t) {
    return rFsa(t);
}

__device__ 
double indicator(int b) {
    return b ? 1.0 : 0.0;
}

/**************** PRODUCT, PURE ENDOWMENT ***************************/
__device__ 
static double b_0_PureEndowment(double t) {
    return 0.0;
}

__device__ 
static double mu_01_PureEndowment(int age,double t) {
    return GM(age,t);
}

__device__ 
static double bj_00_PureEndowment(double t) {
    double bpension = 1.0;
    double pensiontime = 35.0;
    return t == pensiontime ? bpension: 0.0;
}

__device__ 
static double bj_01_PureEndowment(double t) {
    return 0.0; 
}

__device__ 
void bj_ii_PureEndowment(double t, double* result) {
    result[0] += bj_00_PureEndowment(t);
}

    __device__ 
void dy_PureEndowment(int age, double t, double* V,double* result)
{
    result[0] = r(t) * V[0] - b_0_PureEndowment(t) - mu_01_PureEndowment(age,t) * (0 - V[0] + bj_01_PureEndowment(t));
}

/**************** PRODUCT, DEFFEREDLIFEANNUITY ***************************/
__device__ 
static double b_0_DeferredTemporaryLifeAnnuity(double t) {
    int m = 35;
    int n = 10;
    double bpension = 1.0;
    return bpension * indicator(t > m) * indicator(t < m + n);
}

__device__ 
static double mu_01_DeferredTemporaryLifeAnnuity(int age, double t) {
    return GM(age,t);
}

__device__ 
static double bj_00_DeferredTemporaryLifeAnnuity(double t) {
    return 0.0;
}

__device__ 
static double bj_01_DeferredTemporaryLifeAnnuity(double t) {
    return 0.0; 
}

__device__ 
void bj_ii_DeferredTemporaryLifeAnnuity(double t, double* result) {
    result[0] += bj_00_DeferredTemporaryLifeAnnuity(t);
}

    __device__ 
void dy_DeferredTemporaryLifeAnnuity(int age, double t, double* V,double* result)
{
    result[0] = r(t) * V[0] - b_0_DeferredTemporaryLifeAnnuity(t) - mu_01_DeferredTemporaryLifeAnnuity(age,t) * (0 - V[0] + bj_01_DeferredTemporaryLifeAnnuity(t));
}

/**************** PRODUCT, TemporaryLifeAnnuityPremium ***************************/
__device__ 
static double b_0_TemporaryLifeAnnuityPremium(double t) {
    int n = 35;
    int bpremium = 1;
    return -bpremium * indicator(t >= 0) * indicator(t < n);
}

__device__ 
static double mu_01_TemporaryLifeAnnuityPremium(int age, double t) {
    return GM(age,t);
}

__device__ 
static double bj_00_TemporaryLifeAnnuityPremium(double t) {
    return 0.0;
}

__device__ 
static double bj_01_TemporaryLifeAnnuityPremium(double t) {
    return 0.0; 
}

__device__ 
void bj_ii_TemporaryLifeAnnuityPremium(double t, double* result) {
    result[0] += bj_00_TemporaryLifeAnnuityPremium(t);
}

    __device__ 
void dy_TemporaryLifeAnnuityPremium(int age, double t, double* V,double* result)
{
    result[0] = r(t) * V[0] - b_0_TemporaryLifeAnnuityPremium(t) - mu_01_TemporaryLifeAnnuityPremium(age,t) * (0 - V[0] + bj_01_TemporaryLifeAnnuityPremium(t));
}

/**************** PRODUCT, TermInsurance ***************************/
__device__ 
static double b_0_TermInsurance(double t) {
    return 0.0;
}

__device__ 
static double mu_01_TermInsurance(int age, double t) {
    return GM(age,t);
}

__device__ 
static double bj_00_TermInsurance(double t) {
    return 0.0;
}

__device__ 
static double bj_01_TermInsurance(double t) {
    int bdeath = 1;
    int n = 35;
    return bdeath * indicator(t > 0) * indicator(t < n);
}

__device__ 
void bj_ii_TermInsurance(double t, double* result) {
    result[0] += bj_00_TermInsurance(t);
}

    __device__ 
void dy_TermInsurance(int age, double t, double* V,double* result)
{
    result[0] = r(t) * V[0] - b_0_TermInsurance(t) - mu_01_TermInsurance(age,t) * (0 - V[0] + bj_01_TermInsurance(t));
}

/**************** PRODUCT, DisabilityAnnuity ***************************/
__device__ 
static double b_0_DisabilityAnnuity(double t) {
    return 0.0;
}

__device__ 
static double b_1_DisabilityAnnuity(double t) {
    int n = 35;
    int bdisabled = 1;
    //return 0.0;
    return bdisabled * indicator(t > 0) * indicator(t < n);
}

__device__ 
static double GM01_DisabilityAnnuity(int age, double t) {
    return 0.0006 + __pow(10.0, 4.71609 - 10.0 + 0.06*((double)age + t));
}

__device__ 
static double GM02_DisabilityAnnuity(int age, double t) {
    return GM(age,t);
}

__device__ 
static double GM12_DisabilityAnnuity(int age, double t) {
    return GM(age,t);
}

__device__ 
static double mu_01_DisabilityAnnuity(int age, double t) {
    return GM01_DisabilityAnnuity(age,t);
}

__device__ 
static double mu_02_DisabilityAnnuity(int age, double t) {
    return GM02_DisabilityAnnuity(age,t);
}

__device__ 
static double mu_12_DisabilityAnnuity(int age, double t) {
    return GM12_DisabilityAnnuity(age,t);
}

__device__ 
static double bj_00_DisabilityAnnuity(double t) {
    return 0.0;
}

__device__ 
static double bj_01_DisabilityAnnuity(double t) {
    //int n = 35;
    //int bdisabled = 1;
    //return bdisabled * indicator(t > 0) * indicator(t < n);
    return 0.0;
}

__device__ 
static double bj_02_DisabilityAnnuity(double t) {
    return 0.0;
}

__device__ 
static double bj_11_DisabilityAnnuity(double t) {
    return 0.0;
}

__device__ 
static double bj_12_DisabilityAnnuity(double t) {
    return 0.0;
}

__device__ 
void bj_ii_DisabilityAnnuity(double t, double* result) {
    result[0] += bj_00_DisabilityAnnuity(t);
    result[1] += bj_11_DisabilityAnnuity(t);
}

    __device__ 
void dy_DisabilityAnnuity(int age, double t, double* V,double* result)
{
    result[0] = r(t) * V[0] - b_0_DisabilityAnnuity(t) - mu_01_DisabilityAnnuity(age,t) * (V[1] - V[0] + bj_01_DisabilityAnnuity(t)) - mu_02_DisabilityAnnuity(age,t) * (0 - V[0] + bj_02_DisabilityAnnuity(t));
    result[1] = r(t) * V[1] - b_1_DisabilityAnnuity(t) - mu_12_DisabilityAnnuity(age,t) * (0 - V[1] + bj_12_DisabilityAnnuity(t)); 
}

/**************** PRODUCT, DisabilityTermInsurance ***************************/
__device__ 
static double b_0_DisabilityTermInsurance(double t) {
    return 0.0;
}

__device__ 
static double b_1_DisabilityTermInsurance(double t) {
    return 0.0;
}

__device__ 
static double GM01_DisabilityTermInsurance(int age, double t) {
    return 0.0006 + __pow(10.0, 4.71609 - 10.0 + 0.06*((double)age + t));
}

__device__ 
static double GM02_DisabilityTermInsurance(int age,double t) {
    return GM(age,t);
}

__device__ 
static double GM12_DisabilityTermInsurance(int age, double t) {
    return GM(age,t);
}

__device__ 
static double mu_01_DisabilityTermInsurance(int age, double t) {
    return GM01_DisabilityTermInsurance(age,t);
}

__device__ 
static double mu_02_DisabilityTermInsurance(int age, double t) {
    return GM02_DisabilityTermInsurance(age,t);
}

__device__ 
static double mu_12_DisabilityTermInsurance(int age, double t) {
    return GM12_DisabilityTermInsurance(age,t);
}

__device__ 
static double bj_00_DisabilityTermInsurance(double t) {
    return 0.0;
}

__device__ 
static double bj_01_DisabilityTermInsurance(double t) {
    int n = 35;
    int bdisabled = 1;
    return bdisabled * indicator(t > 0) * indicator(t < n);
}

__device__ 
static double bj_02_DisabilityTermInsurance(double t) {
    return 0.0;
}

__device__ 
static double bj_11_DisabilityTermInsurance(double t) {
    return 0.0;
}

__device__ 
static double bj_12_DisabilityTermInsurance(double t) {
    return 0.0;
}

__device__ 
void bj_ii_DisabilityTermInsurance(double t, double* result) {
    result[0] += bj_00_DisabilityTermInsurance(t);
    result[1] += bj_11_DisabilityTermInsurance(t);
}

    __device__ 
void dy_DisabilityTermInsurance(int age, double t, double* V,double* result)
{
    result[0] = r(t) * V[0] - b_0_DisabilityTermInsurance(t) - mu_01_DisabilityTermInsurance(age,t) * (V[1] - V[0] + bj_01_DisabilityTermInsurance(t)) - mu_02_DisabilityTermInsurance(age,t) * (0 - V[0] + bj_02_DisabilityTermInsurance(t));
    result[1] = r(t) * V[1] - b_1_DisabilityTermInsurance(t) - mu_12_DisabilityTermInsurance(age,t) * (0 - V[1] + bj_12_DisabilityTermInsurance(t)); 
}

/**** Policy distributor ****/

__device__ 
void dy(int policy,int age, double t, double* V, double* result) {
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

__device__ 
void bj_ii(int policy, double t, double* result) {
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
