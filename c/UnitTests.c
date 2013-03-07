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

/*************** Estimator tests ***************/

/* Does two matrixes have the same values */
bool is_equal(double** a,double** b,int m,int n) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      //printf(": %f = %f\n",a[i][j],b[i][j]);
      if (fabs(a[i][j] - b[i][j]) > err) {
        //printf("is_equal failed: %f != %f\n",a[i][j],b[i][j]);
        return false;
      }
    }
  }
  return true;
}

bool is_equal_print(double** a,double** b,int m,int n) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      printf(": %f = %f\n",a[i][j],b[i][j]);
      if (fabs(a[i][j] - b[i][j]) > err) {
        //printf("is_equal failed: %f != %f\n",a[i][j],b[i][j]);
        return false;
      }
    }
  }
  return true;
}

static char* test_PureEndowment() {
  policy = 1;
  start_year = 0;
  end_year = 40;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  mu_assert("PureEndowment failed",is_equal(test_values(),estimate(),41,1));
  return 0;
}

static char* test_DeferredTemporaryLifeAnnuity() {
  policy = 2;
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  mu_assert("DeferredTemporaryLifeAnnuity failed",is_equal(test_values(),estimate(),51,1));
  return 0;
}

static char* test_TemporaryLifeAnnuityPremium() {
  policy = 3;
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  mu_assert("TemporaryLifeAnnuityPremium failed",is_equal(test_values(),estimate(),51,1));
  return 0;
}

static char* test_TermInsurance() {
  policy = 4;
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  mu_assert("TermInsurance failed",is_equal(test_values(),estimate(),51,1));
  return 0;
}

static char* test_DisabilityAnnuity() {
  policy = 5;
  start_year = 0;
  end_year = 50;
  relerr = err;
  abserr = err;
  end_year_y[0] = 0.0;
  mu_assert("DisabilityAnnuity failed",is_equal_print(test_values(),estimate(),51,1));
  return 0;
}

/********************* Policy tests ************/

static char* test_PureEndowment_dy() {
  double V[1];
  double result[1];
  V[0] = 1.0;

  policy = 1;
  neqn = 1;
  dy(0,V,result); 
  printf("something: %.16lf\n",result[0]);
  //mu_assert("PureEndowment dy function failed",result[0] == 0.051238);
  return 0;
}

/********************* Main ********************/

static char* all_tests() {
  mu_run_test(test_PureEndowment_dy);
  mu_run_test(test_PureEndowment);
  mu_run_test(test_DeferredTemporaryLifeAnnuity);
  mu_run_test(test_TemporaryLifeAnnuityPremium);
  mu_run_test(test_TermInsurance);
  //mu_run_test(test_DisabilityAnnuity);
  return 0;
}

static void init_tests() {
  construct(1);
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
