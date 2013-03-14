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
#include "Policy_Distributor.h" 
#include <assert.h>
#include <math.h>
#include "Matrix_Library.h"
#include <unistd.h> //only for sleep?
#include <string.h> //only for sleep?
#include "Boolean.h"

//Max,min,sign functions
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a < _b ? _a : _b; })
#define sign(x)  ((x > 0) - ( x < 0))

//Declare functions
bool is_equal();
static void test_all();
static void time_all();
static float time_one();
static bool local_start_to_be_reached(float t,int local_start_year,float* stepsize);
static void allocate_equation_space();
static void xpy();
static void calculate_solutions(int neqn,float t,float stepsize,float *y_diff);
static void local_estimate(int neqn,int local_end_year,int local_start_year,float* stepsize,float* y_diff);
static float calculate_initial_stepsize(int neqn,int start_year,float t,float* y_diff);
static float scale_from_error(float error,bool stepsize_decreased);
static float FindDoubleEpsilon();

//Declare Estimator variables

//Private variables
static float* f_swap;
static float* y;
static float* y_plus_one;
static float* y_plus_one_alternative;
static float DoubleEpsilon;

/******************* Constructor *********************/

/* Construct */
void construct(int n) {
  DoubleEpsilon = FindDoubleEpsilon();
  allocate_equation_space();
}

/* Allocate equation space */
static void allocate_equation_space() {

  //Global for the class
  y_plus_one             = malloc(sizeof(float)*MAX_NEQN);
  y_plus_one_alternative = malloc(sizeof(float)*MAX_NEQN);
  y                      = malloc(sizeof(float)*MAX_NEQN);

}

/********************** Solve *********************/

/* Calculate the actual and the alternative solutions */
//y_plus_one and y_plus_one_alternative will be set
static void calculate_solutions(int neqn,float t,float stepsize,float* y_diff) {

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
static float calculate_solution_error(int neqn,float stepsize) {

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
static void local_estimate(int neqn,int local_end_year,int local_start_year,float *stepsize,float* y_diff) {
  float t = (float)local_end_year;
  
  //Step by step integration.
  bool local_start_reached = false;
  while (!local_start_reached)
  {
    //Variables used in calculations
    bool stepsize_descresed = false;
    float hmin = 26.0f * DoubleEpsilon * fabsf( t );

    local_start_reached = local_start_to_be_reached(t,local_start_year,stepsize);

    calculate_solutions(neqn,t,*stepsize,y_diff);
    float error = calculate_solution_error(neqn,*stepsize);

    //Integreate 1 step
    while(error > 1.0f)
    {
      stepsize_descresed = true;
      local_start_reached = false;

      //Scale down.
      float s = max(0.1f,0.9f / powf( error, 0.2f ));
      *stepsize = s * *stepsize;  

      //Try again.
      calculate_solutions(neqn,t,*stepsize,y_diff);
      error = calculate_solution_error(neqn,*stepsize);
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
static float calculate_initial_stepsize(int neqn,int start_year,float t,float *y_diff)
{
  //Calculate the start value of stepsize
  float s = fabsf( start_year - t );


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

  return  max( s, 26.0f * DoubleEpsilon * max( fabsf( t ), fabsf( start_year - t ) ) );
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
float** estimate(int neqn, int end_year, int start_year,float* yy) { //TODO: yy

  float t = (float) end_year;                   // t
  float y_diff[MAX_NEQN];
  dy( t, y, y_diff);                       // y_diff
  float stepsize = calculate_initial_stepsize(neqn,start_year,t,y_diff); // stepsize

  //Allocate result matrix, calculate (length of result)
  float** result = allocate_float_matrix(end_year-start_year+1,neqn);

  //Solve for one year at a time
  for (int year=end_year; year>start_year; year--) {

    //Add this years benefit to y
    bj_ii(year,y);

    // Integrate over [year,year-1]
    local_estimate(neqn,year,year-1,&stepsize,y_diff);

    //Copy y to results
    for(int i=0;i<neqn;i++)
      result[year-start_year-1][i] = y[i];
  }

  return result;
}

/*************************** Auxiliaries ****************************/

/* Find float epsilon */
static float FindDoubleEpsilon() {
  float r = 1.0f;
  while (1.0f < (1.0f + r))
    r = r / 2.0f;
  return 2.0f * r;
}


