#include <stdio.h>
#include "minunit.h"
#include "Rkf45.h"
#include "Policy_Distributor.h"
#include <math.h>
#include "Boolean.h"

//Accepted error
static double const err = 1e-11;
//Number of tests run
int tests_run = 0;

/*************** Help functions *************/

bool doubles_ewt_print(double a,double b) { //Equal within tolerance
  printf("%.16f = %.16f\n",a,b);
  return fabs(a - b) < err;
}

bool doubles_ewt(double a,double b) { //Equal within tolerance
  return fabs(a - b) < err;
}

/* Does two matrixes have the same values */
bool matrix_ewt(double** a,double** b,int m,int n) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      if (!doubles_ewt(a[i][j],b[i][j])) {
        return false;
      }
    }
  }
  return true;
}

bool matrix_ewt_print(double** a,double** b,int m,int n) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      printf("%.16f = %.16f, ",a[i][j],b[i][j]);
      if (!doubles_ewt(a[i][j],b[i][j])) {
        //return false;
      }
    }
    printf("\n");
  }
  return true;
}

/*************** Estimator tests ***************/

static char* test_PureEndowment() {
  policy = 1;
  neqn = 1;
  start_year = 0;
  end_year = 40;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  mu_assert("PureEndowment failed",matrix_ewt(test_values(),estimate(),41,1));
  return 0;
}

static char* test_DeferredTemporaryLifeAnnuity() {
  neqn = 1;
  policy = 2;
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  mu_assert("DeferredTemporaryLifeAnnuity failed",matrix_ewt(test_values(),estimate(),51,1));
  return 0;
}

static char* test_TemporaryLifeAnnuityPremium() {
  policy = 3;
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  mu_assert("TemporaryLifeAnnuityPremium failed",matrix_ewt(test_values(),estimate(),51,1));
  return 0;
}

static char* test_TermInsurance() {
  neqn = 1;
  policy = 4;
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  mu_assert("TermInsurance failed",matrix_ewt(test_values(),estimate(),51,1));
  return 0;
}

static char* test_DisabilityAnnuity() {
  neqn = 2;
  policy = 5;
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  end_year_y[1] = 0.0;
  mu_assert("DisabilityAnnuity failed",matrix_ewt(test_values(),estimate(),51,2));
  return 0;
}

static char* test_DisabilityTermInsurance() {
  neqn = 2;
  policy = 6;
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  end_year_y[1] = 0.0;
  mu_assert("DisabilityTermInsurance failed",matrix_ewt(test_values(),estimate(),51,2));
  return 0;
}

/********************* Policy tests ************/

static char* test_PureEndowment_dy() {
  double V[1];
  double result[1];
  policy = 1;
  neqn = 1;

  V[0] = 0.0;
  dy(0,V,result); 
  mu_assert("PureEndowment_dy failed",doubles_ewt(result[0],0.0));

  V[0] = 1.0;
  dy(0,V,result); 
  mu_assert("PureEndowment_dy failed",doubles_ewt(result[0],0.0512379042301291));

  V[0] = 1.0;
  dy(1,V,result); 
  mu_assert("PureEndowment_dy failed",doubles_ewt(result[0],0.0513053784411991));
  
  V[0] = -2.0;
  dy(1,V,result); 
  mu_assert("PureEndowment_dy failed",doubles_ewt(result[0],-0.102610756882398));
  return 0;
}

static char* test_DisabilityAnnuity_dy() {
  double V[2];
  double result[2];
  policy = 5;
  neqn = 2;
  V[1] = 0;

  dy(0,V,result); 
  mu_assert("DisabilityAnnuity_dy1 failed",doubles_ewt(result[0],0.0));

  V[0] = 1.0;
  dy(0,V,result); 
  mu_assert("DisabilityAnnuity_dy2 failed",doubles_ewt(result[0],0.0521660675223476));

  V[0] = 1.0;
  dy(1,V,result); 
  mu_assert("DisabilityAnnuity_dy3 failed",doubles_ewt(result[0],0.0522821603136021));
  
  V[0] = -2.0;
  dy(1,V,result); 
  mu_assert("DisabilityAnnuity_dy4 failed",doubles_ewt(result[0],-0.104564320627204));

  V[0] = 0;
  dy(34.625,V,result); 
  mu_assert("DisabilityAnnuity_dy5 failed",doubles_ewt(result[1],-1.0));
  return 0;
}

/********************* Main ********************/

static char* all_tests() {
  mu_run_test(test_PureEndowment_dy);
  mu_run_test(test_DisabilityAnnuity_dy);
  mu_run_test(test_PureEndowment);
  mu_run_test(test_DeferredTemporaryLifeAnnuity);
  mu_run_test(test_TemporaryLifeAnnuityPremium);
  mu_run_test(test_TermInsurance);
  mu_run_test(test_DisabilityAnnuity);
  mu_run_test(test_DisabilityTermInsurance);
  return 0;
}

static void init_tests() {
  construct(2);
}

int main(int argc, char **argv) {
  init_tests();
  char *result = all_tests();
  if (result != 0) {
    printf("%s\n", result);
  }
  else {
    printf("ALL TESTS PASSED\n");
  }
  printf("Tests run: %d\n", tests_run);

  return result != 0;
}
