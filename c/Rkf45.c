#include <stdlib.h>
#include <stdio.h>
#include "Policy_Distributor.h" 

//Declare functions
void kernel();

/* Main test function */
int main(int argc, char const *argv[]) {

  //Launch kernel
  kernel();

  return 0;
}

// Device code
void kernel(int startyear,int endyear,double* start_values,int neqn) {

  double result[10];
  double V[10];

  double t = 0.0;

  dV(t,V,result);

  printf("Result: %lf\n", result[0]);
}
