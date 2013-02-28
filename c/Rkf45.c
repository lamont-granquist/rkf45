#include <stdlib.h>
#include <stdio.h>
#include "Policy_Distributor.h" 

//Declare functions
void estimate();

//Declare Estimator variables
static int neqn;
static double t;

/* Main test function */
int main(int argc, char const *argv[]) {
  
  //Set estimator variables
  neqn = 1;

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

  dV(0.0,V,result_length,neqn,result);

  for (int i = 0;i < result_length;i++)
    for (int j = 0;j < neqn;j++)
      printf("%d: %.2lf, ",i, result[i][j]);
}

