#include <stdlib.h>
#include <stdio.h>
#include "Policy_Distributor.h" 
#include <assert.h>
#include <stdbool.h>

//Declare functions
void construct();
void estimate();
void test();
double** allocate_double_matrix();

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
  construct(neqn);

  //Set estimator variables (Term insurrance)
  relerr = err;
  abserr = err;
  y[0] = 0.0;

  //Some testing, might be deleted
  test();

  //Start estimator
  estimate(0,10);

  return 0;
}

/* Initiate estimator */
void construct(int n) {
  neqn = n;

  y = malloc(sizeof(double)*neqn);
  yp = malloc(sizeof(double)*neqn);
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

/* Estimate range */
void estimate(int start_year,int end_year) {

  t = (double) end_year;
  int result_length = end_year-start_year+1;

  construct(1);

  double** result = allocate_double_matrix(result_length,neqn);

  double V[10];

  dy(0.0,V,result_length,neqn,result);

  print_matrix(result_length,neqn,result);
  printf("test succesful\n");
}

/* Testing */
void test() {
  assert(err < 0.00000001);
}

/*Help functions */
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
