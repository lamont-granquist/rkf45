#include <stdlib.h>
#include <stdio.h>
#include "Policy_Distributor.hu"
//#include "Policies/PureEndowment.h"
/*#include "Policies/DeferredTemporaryLifeAnnuity.h"
#include "Policies/TemporaryLifeAnnuityPremium.h"
#include "Policies/TermInsurance.h"
#include "Policies/DisabilityAnnuity.h"
#include "Policies/DisabilityTermInsurance.h"*/

__device__ void dy(int policy,  t, float* V, float* result) {
  switch(policy)
  {
    case 1:
      //dy_PureEndowment(t,V,result);
    break;
    case 2:
      //dy_DeferredTemporaryLifeAnnuity(t,V,result);
    break;
    case 3:
      //dy_TemporaryLifeAnnuityPremium(t,V,result);
    break;
    case 4:
      //dy_TermInsurance(t,V,result);
    break;
    case 5:
      //dy_DisabilityAnnuity(t,V,result);
    break; 
    case 6:
      //dy_DisabilityTermInsurance(t,V,result);
    break; 
  };
}

__device__ void bj_ii(int policy,  t, float* result) {
  switch(policy)
  {
    case 1:
      //bj_ii_PureEndowment(t,result);
    break;
    case 2:
      //bj_ii_DeferredTemporaryLifeAnnuity(t,result);
    break;
    case 3:
      //bj_ii_TemporaryLifeAnnuityPremium(t,result);
    break;
    case 4:
      //bj_ii_TermInsurance(t,result);
    break;
    case 5:
      //bj_ii_DisabilityAnnuity(t,result);
    break;
    case 6:
      //bj_ii_DisabilityTermInsurance(t,result);
    break;

  };
}

float ** test_values(int policy) {
  switch(policy)
  {
    case 1:
      return //tv_PureEndowment();
    break;
    case 2:
      return //tv_DeferredTemporaryLifeAnnuity();
    break;
    case 3:
      return //tv_TemporaryLifeAnnuityPremium();
    break;
    case 4:
      return //tv_TermInsurance();
    break;
    case 5:
      return //tv_DisabilityAnnuity();
    break;
    case 6:
      return //tv_DisabilityTermInsurance();
    break;
  }
  return //tv_PureEndowment();
}
