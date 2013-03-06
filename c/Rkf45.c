#include <stdlib.h>
#include <stdio.h>
#include "Policy_Distributor.h" 
#include <assert.h>
#include <math.h>

#include <unistd.h> //only for sleep?
#include <string.h> //only for sleep?

//Boolean values
typedef int bool;
#define false 0
#define true 1

//Max,min,sign functions
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a < _b ? _a : _b; })
#define sign(x)  ((x > 0) - ( x < 0))

//Declare functions
static bool local_start_to_be_reached();
static bool is_equal();
static void construct();
static void allocate_equation_space();
static double** estimate();
static double** allocate_double_matrix();
static void print_matrix();
static void xpy();
static double calculate_solutions();
static void local_estimate();
static double calculate_initial_stepsize();
static double scale_from_error();
static double FindDoubleEpsilon();

//Declare Estimator variables
static double const err = 1e-11;

//Public
static int neqn;
static int start_year;
static int end_year;
// dy
// bj_ii
static double relerr;
static double abserr;
static double* end_year_y; 
//Private
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

/* Main test function */
int main(int argc, char const *argv[]) {

  //Construct the estimator
  construct(1);

  //Set estimator variables (Term insurrance)
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;

  assert(is_equal(test_values(),estimate(),51,1));
  printf("test succesful\n");

  return 0;
}


/* Initiate estimator */
static void construct(int n) {
  neqn = n;
  DoubleEpsilon = FindDoubleEpsilon();
  allocate_equation_space();
}

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

/* Solve */
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

     /*
     if (test_values == 16) {
     printf("%.16lf\n",f1[0]);
     printf("%.16lf\n",f2[0]);
     printf("%.16lf\n",f3[0]);
     printf("%.16lf\n",f4[0]);
     printf("%.16lf\n",f5[0]);
     }
     test_add();
     */

  //Calculate solution
  lcd_stepsize = stepsize / 7618050.0;
  for (int i = 0; i < neqn; i++ )
    y_plus_one[i] = y[i] + lcd_stepsize * ( ( 902880.0 * y_diff[i] + 
          ( 3855735.0 * f3[i] - 1371249.0 * f4[i] ) ) + ( 3953664.0 * f2[i] + 277020.0 * f5[i] ) );

  //Calculate alternative solution
  for (int i = 0; i < neqn; i++ )
    y_plus_one_alternative[i] = ( -2090.0 * y_diff[i] + ( 21970.0 * f3[i] - 15048.0 * f4[i] ) ) + ( 22528.0 * f2[i] - 27360.0 * f5[i] );

}

//Calculate the error.
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

/* Move */
static void local_estimate() {
  
  //Step by step integration.
  bool local_start_reached = false;

  while (!local_start_reached)
  {
    //Variables used in calculations
    bool hfaild = false;
    double hmin = 26.0 * DoubleEpsilon * fabs( t );

    local_start_reached = local_start_to_be_reached();

    calculate_solutions();
    double error = calculate_solution_error();

    //Integreate 1 step
    while(error > 1.0)
    {
      hfaild = true;
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
    double scale = scale_from_error(error,hfaild);
    stepsize = sign ( stepsize ) * max ( scale * fabs( stepsize ), hmin );
  }
}

/**************** Move help functions ****************/

static bool local_start_to_be_reached() {
    double dt = local_end_year - t;
    //Reaction if stepsize is going to the endpoint.
    //Look 2.0 steps ahead, so stepsize is not 'suddenly' decreased.
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
      /* test startvalues
         printf("abserr: %.40lf\n",abserr);
         printf("tol: %.40lf\n",tol);
         printf("ypk: %.40lf\n",ypk);
         printf("y[k]: %.40lf\n",y[k]);
         printf("y_diff[k]: %.40lf\n",y_diff[k]);
         printf("stepsize: %.40lf\n",stepsize);
         */
    }
  }

  return  max( stepsize, 26.0 * DoubleEpsilon * max( fabs( t ), fabs( start_year - t ) ) );
}

/* Scale from error calculations */
static double scale_from_error(double error,bool hfailed) {
  double scale = min(5.0,0.9 / pow( error, 0.2 ));

  if (hfailed)
    scale = min( scale, 1.0 );

  return scale;
}

/* Find double epsilon */
static double FindDoubleEpsilon() {
  double r = 1.0;
  while (1.0 < (1.0 + r))
    r = r / 2.0;
  return 2.0 * r;
}

/* Estimate range */
static double** estimate() {

  //Set the initial values
  memcpy(y,end_year_y,neqn);      // y
  t = (double) end_year;          // t
  dy( t, y, y_diff);
  stepsize = calculate_initial_stepsize();

  //Allocate result matrix, calculate m (length of result)
  m = end_year-start_year+1;
  double** result = allocate_double_matrix(m,neqn);

  //Allocate benefit array
  double* benefit = malloc(sizeof(double)*neqn);

  //Solve for one year at a time
  for (int year=end_year; year>start_year; year--) {
    //calcuate this years benefit
    bj_ii(year,benefit);

    //add benefit
    xpy(y,benefit);

    //Integate
    local_start_year = year;
    local_end_year = year-1;
    local_estimate();

    //Copy y to results
    for(int i=0;i<neqn;i++)
      result[year-start_year-1][i] = y[i];
  }

  return result;
}

/*******************      Matrix functions      *******************/
//TODO: Put in a file for themself.

/* Addtion of all elements in array b to array a */
static void xpy(double* x, double* y,int length) {
  for (int i=0; i<length; i++)
    x[i] = x[i]+y[i];
}

/*Allocate Matrix */
static double** allocate_double_matrix(int m, int n) {
  /* Allocate memory for the elements */
  double *mem = malloc(m * n * sizeof(double));
  /* Allocate memory for the matrix array */
  double **mat = malloc(m * sizeof(double *));
  /* Setup array */
  if (mem != NULL && mat != NULL) {
    int i;
    for (i = 0; i < m; ++i) {
      mat[i] = &mem[i * n];
    }
  } else {
    printf("Out of memory!\n"); exit(-1);
  }
  return mat;
}

/* Print Matrix */
static void print_matrix(int m,int n,double** mat) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      printf("%.16lf, ",mat[i][j]);
    }
    printf("\n");
  }
}

/* Does two matrixes have the same values */
static bool is_equal(double** a,double** b,int m,int n) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      //printf("%.16lf, ",a[i][j]);
      if (fabs(a[i][j] - b[i][j]) > err)
        return false;
    }
  }
  return true;
}

/* Free Matrix */
static void free_double_matrix(double **mat) {
  /* Free elements */
  free(*mat);
  /* Free matrix */
  free(mat);
}

/* Testing TODO: Delete */
/*static void test() {
  assert(err < 0.00000001);

//Test xpy
int a[5] = {1,2,3,4,5};
int b[5] = {5,4,3,2,1};
xpy(a,b,5);

for(int i=0;i<5;i++)
assert(a[i]==6);

//Test fabs
assert(fabs(-100.12434588557878543)==100.12434588557878543);

//Test max
assert(max(10,5)==10);
assert(max(-10,5)==5);

//Test max
assert(min(10,5)==5);
assert(min(-10,5)==-10);
}*/
