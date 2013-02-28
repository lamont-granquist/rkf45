#include <stdlib.h>
#include <stdio.h>
#include "TermInsurance.h"
#include "Policy_Distributor.h"

int policy = 1;

void dV(double t,double* V, double* result) {
  switch(policy)
  {
    case 1:
      //Test();
      dV_TermInsurance(t,V,result);
    break;
  };
}
