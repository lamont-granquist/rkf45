#include <stdlib.h>
#include <stdio.h>
#include "Policy_Distributor.h" 
#include <assert.h>
#include <math.h>

#include <unistd.h> //only for sleep?

//Boolean values
typedef int bool;
#define false 0
#define true 1

//Max,min,sign functions
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b);_a < _b ? _a : _b; })
#define sign(x)  ((x > 0) - ( x < 0))

//Declare functions - TODO: make static
static bool is_equal();
void construct();
double** estimate();
void test();
double** allocate_double_matrix();
void print_matrix();
void xpy();
double solve();
void move();
double h_startvalue();
double scale_from_error();
bool test_bool();
double FindDoubleEpsilon();

//Declare Estimator variables
double const err = 1e-11;

//Private
static bool first_move = true;
static int neqn;
static double t;
static double h;
static double* yp;
static double* f1;
static double* f2;
static double* f3;
static double* f4;
static double* f5;
static double* f_swap;
static double* y_plus_one;
static double* y_plus_one_alternative;
static double DoubleEpsilon;

//Public
static double* y;
static double relerr;
static double abserr;

static int m; //Result length;

/* Main test function */
int main(int argc, char const *argv[]) {

  //Construct the estimator
  construct(1);

  //Set estimator variables (Term insurrance)
  relerr = err;
  abserr = err;
  y[0] = 0.0;

  //Some testing, might be deleted
  test();

  //Start estimator
  //estimate(0,50);
  //Print result:
  //print_matrix(51,1,estimate(0,50));

  //print_matrix(51,1,test_values());

  assert(is_equal(test_values(),estimate(0,50),51,1));
  printf("test succesful\n");

  return 0;
}


/* Initiate estimator */
void construct(int n) {
  neqn = n;
  DoubleEpsilon = FindDoubleEpsilon();
  y                      = malloc(sizeof(double)*neqn);
  yp                     = malloc(sizeof(double)*neqn);
  f1                     = malloc(sizeof(double)*neqn);
  f2                     = malloc(sizeof(double)*neqn);
  f3                     = malloc(sizeof(double)*neqn);
  f4                     = malloc(sizeof(double)*neqn);
  f5                     = malloc(sizeof(double)*neqn);
  f_swap                 = malloc(sizeof(double)*neqn);
  y_plus_one             = malloc(sizeof(double)*neqn);
  y_plus_one_alternative = malloc(sizeof(double)*neqn);
}

/* test function */
int TEST_MAX = 50;
static int NO_tests = 1;

void test_add() {
  NO_tests++;
  if (NO_tests> TEST_MAX)
    exit(0);
}

/* Solve */
double solve() {

  double ch = h / 4.0;


  //f1
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + ch * yp[i];
  dy ( t + ch, f_swap, f1 );

  //f2
  ch = 3.0 * h / 32.0;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + ch * ( yp[i] + 3.0 * f1[i] );
  dy ( t + 3.0 * h / 8.0, f_swap, f2 );

  //f3
  ch = h / 2197.0;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + ch * ( 1932.0 * yp[i] + ( 7296.0 * f2[i] - 7200.0 * f1[i] ) );
  dy ( t + 12.0 * h / 13.0, f_swap, f3 );

  //f4
  ch = h / 4104.0;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + ch * ( ( 8341.0 * yp[i] - 845.0 * f3[i] ) + 
        ( 29440.0 * f2[i] - 32832.0 * f1[i] ) );
  dy ( t + h, f_swap, f4 );

  //f5
  ch = h / 20520.0;
  for (int i = 0; i < neqn; i++ )
    f_swap[i] = y[i] + ch * ( ( -6080.0 * yp[i] + 
          ( 9295.0 * f3[i] - 5643.0 * f4[i] ) ) + ( 41040.0 * f1[i] - 28352.0 * f2[i] ) );
  dy ( t + h / 2.0, f_swap, f5 );

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
  ch = h / 7618050.0;
  for (int i = 0; i < neqn; i++ )
    y_plus_one[i] = y[i] + ch * ( ( 902880.0 * yp[i] + 
          ( 3855735.0 * f3[i] - 1371249.0 * f4[i] ) ) + ( 3953664.0 * f2[i] + 277020.0 * f5[i] ) );

  //Calculate alternative solution
  for (int i = 0; i < neqn; i++ )
    y_plus_one_alternative[i] = ( -2090.0 * yp[i] + ( 21970.0 * f3[i] - 15048.0 * f4[i] ) ) + ( 22528.0 * f2[i] - 27360.0 * f5[i] );

  /*
     printf("yp:   %.16lf\n",yp[0]);
     printf("ypo:  %.16lf\n",y_plus_one[0]);
     printf("ypoa: %.16lf\n",y_plus_one_alternative[0]);
     */

  //Calculate the error.
  double biggest_difference = 0.0;

  double scale = 2.0 / relerr; //scale
  double ae = scale * abserr;  //absolute error

  for (int i = 0; i < neqn; i++ )
  {
    double et = fabs( y[i] ) + fabs( y_plus_one[i] ) + ae;
    double ee = fabs( y_plus_one_alternative[i] );

    /*
       printf("et:   %.14lf\n",et);
       printf("ee:   %.14lf\n",ee);
       */

    biggest_difference = max ( biggest_difference, ee / et );
  }

  /*printf("bg:   %.16lf\n",biggest_difference);
    printf("yp1:  %.16lf\n",y_plus_one[0]);
    printf("yp1a: %.16lf\n",y_plus_one_alternative[0]);
    printf("err:  %.16lf\n",fabs( h ) * biggest_difference * scale / 752400.0);*/

  //Return the err
  return fabs( h ) * biggest_difference * scale / 752400.0;
}

/* Move */
void move(double t_end) {

  //Init
  if (first_move) {
    first_move = false;

    dy( t, y, yp);

    h = h_startvalue(t_end);
  }

  //Step by step integration.
  bool end_reached = false;

  while (!end_reached)
  {
    //Variables used in calculations
    bool hfaild = false;
    double dt = t_end - t;
    double hmin = 26.0 * DoubleEpsilon * fabs( t );

    //Reaction if h is going to the endpoint.
    //Look 2.0 steps ahead, so stepsize is not 'suddenly' decreased.
    if ( 2.0 * fabs( h ) > fabs( dt ) )
    {
      if ( fabs( dt ) <= fabs( h ) ) //Final step?
      {
        end_reached = true; //Return output
        h = dt;                   //Let h hit output point
      }
      else
      {
        h = 0.5 * dt; // If not final step, set h to be second final step. (evens out)
      }
    }

    double error = solve();

    //Integreate 1 step
    while(error > 1.0)
    {
      hfaild = true;
      end_reached = false;

      //Scale down.
      double s = max(0.1,0.9 / pow( error, 0.2 ));
      h = s * h;  

      //Try again.
      error = solve();
    }

    //Advance in time
    t = t + h; 

    //Apply solution
    for (int i = 0; i < neqn; i++ )
      y[i] = y_plus_one[i];

    /*Print test data
      if (y[0] != 0) {
      printf("t: %.20lf\n",t);
      printf("h: %.20lf\n",h);
      printf("y[0]: %.20lf\n",y[0]);
      sleep(1);
      }*/

    //Update yp
    dy ( t, y, yp );

    //Apply scale to stepsize
    double scale = scale_from_error(error,hfaild);
    h = sign ( h ) * max ( scale * fabs( h ), hmin );
  }
}

/**************** Move help functions ****************/

/* Calculate h's startvalue */
double h_startvalue(double t_end)
{
  //Calculate the start value of h
  double h = fabs( t_end - t );

  for (int k = 0; k < neqn; k++ )
  {
    double tol = relerr * fabs( y[k] ) + abserr;
    if ( 0.0 < tol )
    {
      double ypk = fabs( yp[k] );
      if ( tol < ypk * pow( h, 5 ) )
      {
        h = pow( ( tol / ypk ), 0.2 );
        printf("this should not happen.\n");
      }
      /* test startvalues
         printf("abserr: %.40lf\n",abserr);
         printf("tol: %.40lf\n",tol);
         printf("ypk: %.40lf\n",ypk);
         printf("y[k]: %.40lf\n",y[k]);
         printf("yp[k]: %.40lf\n",yp[k]);
         printf("h: %.40lf\n",h);
         */
    }
  }

  return  max( h, 26.0 * DoubleEpsilon * max( fabs( t ), fabs( t_end - t ) ) );
}

/* Scale from error calculations */
double scale_from_error(double error,bool hfailed) {
  double scale = min(5.0,0.9 / pow( error, 0.2 ));

  if (hfailed)
    scale = min( scale, 1.0 );

  return scale;
}

/* Find double epsilon */
double FindDoubleEpsilon() {
  double r = 1.0;
  while (1.0 < (1.0 + r))
    r = r / 2.0;
  return 2.0 * r;
}

/* Estimate range */
double** estimate(int start_year,int end_year) {

  //Start at the end year
  t = (double) end_year;

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
    move(year-1);

    //Copy y to results
    for(int i=0;i<neqn;i++)
      result[year-start_year-1][i] = y[i];
  }

  return result;
}


/* Testing */
void test() {
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
}

/* Addtion of all elements in array b to array a */
void xpy(double* x, double* y,int length) {
  for (int i=0; i<length; i++)
    x[i] = x[i]+y[i];
}

/*Allocate Matrix */
double** allocate_double_matrix(int m, int n) {
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
void print_matrix(int m,int n,double** mat) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      printf("%.16lf, ",mat[i][j]);
    }
    printf("\n");
  }
}

/* IsEqual? (2 Matrixes) */
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
void free_double_matrix(double **mat) {
  /* Free elements */
  free(*mat);
  /* Free matrix */
  free(mat);
}
