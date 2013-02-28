#include <stdlib.h>
#include <stdio.h>
#include "Policy_Distributor.h" 
#include <assert.h>
#include <stdbool.h>

//Declare functions
void construct();
double** estimate();
void test();
double** allocate_double_matrix();
void print_matrix();
void xpy();

//Declare Estimator variables
double const err = 1e-11;

//Private
static bool first_move = true;
static int neqn;
static double t;
static double* yp,f1,f2,f3,f4,f5,f_swap,y_plus_one,y_plus_one_alternative;

//Public
static double* y;
static double relerr;
static double abserr;

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
  print_matrix(10,1,estimate(0,10));

  printf("test succesful\n");
  return 0;
}

/* Initiate estimator */
void construct(int n) {
  neqn = n;

  y = malloc(sizeof(double)*neqn);
  yp = malloc(sizeof(double)*neqn);
}

/* Estimate range */
double** estimate(int start_year,int end_year) {

  //Start at the end year
  t = (double) end_year;

  //Allocate result matrix
  int result_length = end_year-start_year+1;
  double** result = allocate_double_matrix(result_length,neqn);

  //Allocate benefit array
  double* benefit = malloc(sizeof(double)*neqn);

  //Solve for one year at a time
  for (int year=end_year; year>start_year; year--) {
    //calcuate this years benefit
    bj_ii(year,benefit);

    //add benefit
    xpy(y,benefit);

    //Integate
    //TODO

    //Copy results
  }

  //dy(0.0,y,result_length,neqn,result);

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
      printf("%.2lf, ",mat[i][j]);
    }
    printf("\n");
  }
}

/* Free Matrix */
void free_double_matrix(double **mat) {
  /* Free elements */
  free(*mat);
  /* Free matrix */
  free(mat);
}
