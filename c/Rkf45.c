/*
 * C implementation of the Rkf45 algoritm.
 */

/********************* INIT *******************/

//Library inclusion
#include <stdlib.h>
#include <stdio.h>
#include "Policy_Distributor.h" 
#include <assert.h>
#include <math.h>
#include "Matrix_Library.h"
#include <unistd.h> //only for sleep?
#include <string.h> //only for sleep?
#include <time.h>

//Boolean values
typedef int bool;
#define false 0
#define true 1

//Max,min,sign functions
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a < _b ? _a : _b; })
#define sign(x)  ((x > 0) - ( x < 0))

//Declare functions
bool is_equal();
static void test_all();
static void time_all();
static double time_one();
static bool local_start_to_be_reached();
static void construct();
static void allocate_equation_space();
static double** estimate();
static double** compute();
static void xpy();
static double calculate_solutions();
static void local_estimate();
static double calculate_initial_stepsize();
static double scale_from_error();
static double FindDoubleEpsilon();

//Declare Estimator variables
static double const err = 1e-11;

//Public variables
static int neqn;
static int start_year;
static int end_year;
// dy
// bj_ii
static double relerr;
static double abserr;
static double* end_year_y; 

//Private variables
static double t;
static double stepsize;
static double* f1;
static double* f2;
static double* f3;
static double* f4;
static double* f5;
static double* f_swap;
static double* y;
static double* y_diff;
static double* y_plus_one;
static double* y_plus_one_alternative;
static int local_start_year;
static int local_end_year;
static double DoubleEpsilon;

static int m; //Result length;

/******************* Constructor *********************/

/* Construct */
static void construct(int n) {
  neqn = n;
  DoubleEpsilon = FindDoubleEpsilon();
  allocate_equation_space();
}

/* Allocate equation space */
static void allocate_equation_space() {
  //Global for the class
  y_plus_one             = malloc(sizeof(double)*neqn);
  end_year_y             = malloc(sizeof(double)*neqn);
  y                      = malloc(sizeof(double)*neqn);
  y_diff                     = malloc(sizeof(double)*neqn);

  //Temporary for the solve method
  f1                     = malloc(sizeof(double)*neqn);
  f2                     = malloc(sizeof(double)*neqn);
  f3                     = malloc(sizeof(double)*neqn);
  f4                     = malloc(sizeof(double)*neqn);
  f5                     = malloc(sizeof(double)*neqn);
  f_swap                 = malloc(sizeof(double)*neqn);
  y_plus_one_alternative = malloc(sizeof(double)*neqn);
}

/********************** Solve *********************/

/* Calculate the actual and the alternative solutions */
//y_plus_one and y_plus_one_alternative will be set
static double calculate_solutions() {

  double lcd_stepsize = stepsize / 4.0; //lowest common denominator of stepsize

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

/* Calculate the error of the solution */
//Pure
static double calculate_solution_error() {

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
static void local_estimate() {
  
  //Step by step integration.
  bool local_start_reached = false;
  while (!local_start_reached)
  {
    //Variables used in calculations
    bool stepsize_descresed = false;
    double hmin = 26.0 * DoubleEpsilon * fabs( t );

    local_start_reached = local_start_to_be_reached();

    calculate_solutions();
    double error = calculate_solution_error();

    //Integreate 1 step
    while(error > 1.0)
    {
      stepsize_descresed = true;
      local_start_reached = false;

      //Scale down.
      double s = max(0.1,0.9 / pow( error, 0.2 ));
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
    double scale = scale_from_error(error,stepsize_descresed);
    stepsize = sign ( stepsize ) * max ( scale * fabs( stepsize ), hmin );
  }
}

/**************** Local estimation help functions ****************/


/* React if the "local start year" is about to be reached */
//Effects stepsize, returns whether the start year is reached
static bool local_start_to_be_reached() {
    double dt = local_end_year - t;
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
}

/* Calculate stepsize's startvalue */
static double calculate_initial_stepsize()
{
  //Calculate the start value of stepsize
  double stepsize = fabs( start_year - t );


  for (int k = 0; k < neqn; k++ )
  {
    double tol = relerr * fabs( y[k] ) + abserr;
    if ( 0.0 < tol )
    {
      double ypk = fabs( y_diff[k] );
      if ( tol < ypk * pow( stepsize, 5 ) )
      {
        stepsize = pow( ( tol / ypk ), 0.2 );
        printf("this should not happen.\n");
      }
    }
  }

  return  max( stepsize, 26.0 * DoubleEpsilon * max( fabs( t ), fabs( start_year - t ) ) );
}

/* Scale from error calculations */
static double scale_from_error(double error,bool stepsize_decreased) {
  double scale = min(5.0,0.9 / pow( error, 0.2 ));

  if (stepsize_decreased)
    scale = min( scale, 1.0 );

  return scale;
}

/*********************** Estimate **************************/

/* Estimate range */
static double** estimate() {

  for(int i = 0;i<neqn;i++)                // y
    y[i] = end_year_y[i];
  t = (double) end_year;                   // t
  dy( t, y, y_diff);                       // y_diff
  stepsize = calculate_initial_stepsize(); // stepsize

  //Allocate result matrix, calculate m (length of result)
  m = end_year-start_year+1;
  double** result = allocate_double_matrix(m,neqn);

  //Solve for one year at a time
  for (int year=end_year; year>start_year; year--) {

    //Add this years benefit to y
    bj_ii(year,y);

    // Integrate over [year,year-1]
    local_start_year = year;
    local_end_year = year-1;
    local_estimate();

    //Copy y to results
    for(int i=0;i<neqn;i++)
      result[year-start_year-1][i] = y[i];
  }

  return result;
}

/*************************** Auxiliaries ****************************/

/* Find double epsilon */
static double FindDoubleEpsilon() {
  double r = 1.0;
  while (1.0 < (1.0 + r))
    r = r / 2.0;
  return 2.0 * r;
}

/* Main test function */
int main(int argc, char const *argv[]) {

  //Construct the estimator
  construct(1);

  test_all();
  time_all(12288);

  return 0;
}

void test_all() {
  policy = 1;
  assert(is_equal(test_values(),compute(),41,1));
  policy = 2;
  assert(is_equal(test_values(),compute(),51,1));
  policy = 3;
  assert(is_equal(test_values(),compute(),51,1));
  policy = 4;
  assert(is_equal(test_values(),compute(),51,1));
  printf("Tests passed\n");
}

void time_all(int customers) {
  policy = 1;
  printf("PureEndowment:                %f\n",time_one(customers));
  policy = 2;
  printf("DeferredTemporaryLifeAnnuity: %f\n",time_one(customers));
  policy = 3;
  printf("TemporaryLifeAnnuityPremium:  %f\n",time_one(customers));
  policy = 4;
  printf("TermInsurance:                %f\n",time_one(customers));
  //printf("DisabilityAnnuity:            %f\n",time_one(customers));
  //printf("DisabilityTermInsurance:      %f\n",time_one(customers));
}

double time_one(int customers) {
  clock_t start = clock();
  for (int i = 0;i<customers;i++)
    compute();
  clock_t end = clock();
  return (double) (end - start) * 1000 / CLOCKS_PER_SEC;
}

static double** compute() {
  //Set estimator variables (Term insurrance)
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  return estimate();
}

/************************** To be removed?? *********************/

/* Does two matrixes have the same values */
bool is_equal(double** a,double** b,int m,int n) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      if (fabs(a[i][j] - b[i][j]) > err)
        printf("is_equal failed: %f != %f\n",a[i][j],b[i][j]);
        //return false;
    }
  }
  return true;
}

