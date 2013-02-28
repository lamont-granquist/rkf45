#include <stdlib.h>
#include <stdio.h>
#include "TermInsurance.h"
#include "Policy_Distributor.h"

int policy = 1;

void dy(double t, double* V, double* result) {
  switch(policy)
  {
    case 1:
      dy_TermInsurance(t,V,result);
    break;
  };
}

void bj_ii(double t, double* result) {
  switch(policy)
  {
    case 1:
      bj_ii_TermInsurance(t,result);
    break;
  };
}
