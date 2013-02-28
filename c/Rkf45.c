#include <stdlib.h>
#include <stdio.h>
#include "Policy_Distributor.h" 
#include <assert.h>

//Declare functions
void estimate();
void test();

//Declare Estimator variables
double const err = 1e-11;

static int neqn;
static double relerr;
static double abserr;
static double t;

/* Main test function */
int main(int argc, char const *argv[]) {
  
  //Set estimator variables (Term insurrance)
  neqn = 1;
  relerr = err;
  abserr = err;
  //y = { 0.0 };
  
  //Some testing, might be refactored
  test();

  //Start estimator
  estimate(0,40);

  return 0;
}

// Device code
void estimate(int start_year,int end_year) {

  t = (double) end_year;
  int result_length = end_year-start_year+1;

  double result[result_length][neqn];

  double V[10];

  dy(0.0,V,result_length,neqn,result);

  for (int i = 0;i < result_length;i++)
    for (int j = 0;j < neqn;j++)
      printf("%d: %.2lf, ",i, result[i][j]);
}

/* Testing */
void test() {
  assert(err < 0.00000001);
}
