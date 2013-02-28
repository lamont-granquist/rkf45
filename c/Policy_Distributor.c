#include <stdlib.h>
#include <stdio.h>
#include "TermInsurance.h"
#include "Policy_Distributor.h"

int policy = 1;

void dV(double t,double* V, int m, int n, double result[m][n]) {
  switch(policy)
  {
    case 1:
      dV_TermInsurance(t,V,m,n,result);
    break;
  };
}
