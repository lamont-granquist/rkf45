#include <stdlib.h>
#include <stdio.h>
#include "TermInsurance.h"
#include "PureEndowment.h"
#include "Policy_Distributor.h"

int policy = 0;

void dy(double t, double* V, double* result) {
  switch(policy)
  {
    case 1:
      dy_PureEndowment(t,V,result);
    break;
    case 4:
      dy_TermInsurance(t,V,result);
    break;
  };
}

void bj_ii(double t, double* result) {
  switch(policy)
  {
    case 1:
      bj_ii_PureEndowment(t,result);
    break;
    case 4:
      bj_ii_TermInsurance(t,result);
    break;
  };
}

double** test_values() {
  switch(policy)
  {
    case 1:
      return tv_PureEndowment();
    break;
    case 4:
      return tv_TermInsurance();
    break;
  }
  return tv_PureEndowment();
}
