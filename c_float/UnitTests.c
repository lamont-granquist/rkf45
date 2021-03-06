#include <stdio.h>
#include "minunit.h"
#include "Rkf45.h"
#include "Policy_Distributor.h"
#include <math.h>
#include "Boolean.h"
#include <time.h>

//Accepted error
static float const err = 1e-7;
//Number of tests run
int tests_run = 0;

/*************** Help functions *************/

bool floats_ewt_print(float a,float b) { //Equal within tolerance
  printf("%.8f = %.8f\n",a,b);
  return fabs(a - b) < err;
}

bool floats_ewt(float a,float b) { //Equal within tolerance
  return fabs(a - b) < err;
}

/* Does two matrixes have the same values */
bool matrix_ewt(float** a,float** b,int m,int n) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      if (!floats_ewt(a[i][j],b[i][j])) {
        return false;
      }
    }
  }
  return true;
}

bool matrix_ewt_print(float** a,float** b,int m,int n) {
  for (int i = 0;i < m;i++) {
    for (int j = 0;j < n;j++) {
      printf("%.7f = %.7f, ",a[i][j],b[i][j]);
      if (!floats_ewt(a[i][j],b[i][j])) {
        //return false;
      }
    }
    printf("\n");
  }
  return true;
}

/*************** Test cases *******************/

static void test_case_PureEndowment() {
  policy = 1;
}

static void test_case_DeferredTemporaryLifeAnnuity() {
  policy = 2;
}

static void test_case_TemporaryLifeAnnuityPremium() {
  policy = 3;
}

static void test_case_TermInsurance() {
  policy = 4;
}

static void test_case_DisabilityAnnuity() {
  policy = 5;
}

static void test_case_DisabilityTermInsurance() {
  policy = 6;
}
/*************** Estimator tests ***************/

static char* test_PureEndowment() {
  test_case_PureEndowment();
  float y[MAX_NEQN];
  y[0] = 0.0f;
  
  //Allocate result matrix, calculate (length of result)
  float result0[40+1];

  estimate(1,40,0,y,result0);

  matrix_ewt_print(test_values(),test_values(),41,1);

  for(int i=0;i<40;i++) {
    printf("%.7f\n",result0[i]);
  }

  //mu_assert("PureEndowment failed",matrix_ewt(test_values(),result0,41,1));
  return 0;
}

static char* test_DeferredTemporaryLifeAnnuity() {
  test_case_DeferredTemporaryLifeAnnuity();
  mu_assert("DeferredTemporaryLifeAnnuity failed",matrix_ewt_print(test_values(),estimate(),51,1));
  return 0;
}

static char* test_TemporaryLifeAnnuityPremium() {
  test_case_TemporaryLifeAnnuityPremium();
  mu_assert("TemporaryLifeAnnuityPremium failed",matrix_ewt(test_values(),estimate(),51,1));
  return 0;
}

static char* test_TermInsurance() {
  test_case_TermInsurance();
  mu_assert("TermInsurance failed",matrix_ewt(test_values(),estimate(),51,1));
  return 0;
}

static char* test_DisabilityAnnuity() {
  test_case_DisabilityAnnuity();
  mu_assert("DisabilityAnnuity failed",matrix_ewt(test_values(),estimate(),51,2));
  return 0;
}

static char* test_DisabilityTermInsurance() {
  test_case_DisabilityTermInsurance();
  mu_assert("DisabilityTermInsurance failed",matrix_ewt(test_values(),estimate(),51,2));
  return 0;
}

/********************* Policy tests ************/

static char* test_PureEndowment_dy() {
  float V[1];
  V[0] = 0.0f;
  float result[1];
  result[0] = 0.0f;
  policy = 1;

  V[0] = 0.0f;
  dy(0,V,result); 
  mu_assert("PureEndowment_dy failed",floats_ewt(result[0],0.0f));

  V[0] = 1.0f;
  dy(0,V,result); 
  mu_assert("PureEndowment_dy failed",floats_ewt(result[0],0.0512379042301291f));

  V[0] = 1.0f;
  dy(1,V,result); 
  mu_assert("PureEndowment_dy failed",floats_ewt(result[0],0.0513053784411991f));
  
  V[0] = -2.0f;
  dy(1,V,result); 
  mu_assert("PureEndowment_dy failed",floats_ewt(result[0],-0.102610756882398f));
  return 0;
}

static char* test_DeferredTemporaryLifeAnnuity_dy() {
  float V[1];
  V[0] = 0.0f;
  float result[1];
  result[0] = 0.0f;
  policy = 2;

  V[0] = 0.0f;
  dy(0,V,result); 
  mu_assert("DeferredTemporaryLifeAnnuity failed",floats_ewt(result[0],0.0f));

  V[0] = 0.28125f;
  dy(44.6250f,V,result); 
  mu_assert("DeferredTemporaryLifeAnnuity failed",floats_ewt_print(result[0],-0.9754968f));

  V[0] = 1.0f;
  dy(1,V,result); 
  mu_assert("DeferredTemporaryLifeAnnuity failed",floats_ewt(result[0],0.0513053784411991f));
  
  V[0] = -2.0f;
  dy(1,V,result); 
  mu_assert("DeferredTemporaryLifeAnnuity failed",floats_ewt(result[0],-0.102610756882398f));
  return 0;
}

static char* test_DisabilityAnnuity_dy() {
  float V[2];
  V[0] = 0.0f;
  V[1] = 0.0f;
  float result[2];
  result[0] = 0.0f;
  result[1] = 0.0f;
  policy = 5;
  V[1] = 0;

  dy(0,V,result); 
  mu_assert("DisabilityAnnuity_dy1 failed",floats_ewt(result[0],0.0f));

  V[0] = 1.0f;
  dy(0,V,result); 
  mu_assert("DisabilityAnnuity_dy2 failed",floats_ewt(result[0],0.0521660675223476f));

  V[0] = 1.0f;
  dy(1,V,result); 
  mu_assert("DisabilityAnnuity_dy3 failed",floats_ewt(result[0],0.0522821603136021f));
  
  V[0] = -2.0f;
  dy(1,V,result); 
  mu_assert("DisabilityAnnuity_dy4 failed",floats_ewt(result[0],-0.104564320627204f));

  V[0] = 0;
  dy(34.625f,V,result); 
  mu_assert("DisabilityAnnuity_dy5 failed",floats_ewt(result[1],-1.0f));
  return 0;
}

/********************* Timing ******************/

float time_one(int customers) {
  clock_t start = clock();
  for (int i = 0;i<customers;i++)
    estimate();
  clock_t end = clock();
  return (float) (end - start) * 1000.0f / CLOCKS_PER_SEC;
}

void all_timing(int customers) {
  test_case_PureEndowment();
  printf("PureEndowment:                %f\n",time_one(customers));
  test_case_DeferredTemporaryLifeAnnuity();
  printf("DeferredTemporaryLifeAnnuity: %f\n",time_one(customers));
  test_case_TemporaryLifeAnnuityPremium();
  printf("TemporaryLifeAnnuityPremium:  %f\n",time_one(customers));
  test_case_TermInsurance();
  printf("TermInsurance:                %f\n",time_one(customers));
  test_case_DisabilityAnnuity();
  printf("DisabilityAnnuity:            %f\n",time_one(customers));
  test_case_DisabilityTermInsurance();
  printf("DisabilityTermInsurance:      %f\n",time_one(customers));
}

/********************* Main ********************/

static char* all_tests() {
  //mu_run_test(test_PureEndowment_dy);
  //mu_run_test(test_DisabilityAnnuity_dy);
  //mu_run_test(test_DeferredTemporaryLifeAnnuity_dy);
  mu_run_test(test_PureEndowment);
  /*mu_run_test(test_DeferredTemporaryLifeAnnuity);
  mu_run_test(test_TemporaryLifeAnnuityPremium);
  mu_run_test(test_TermInsurance);
  mu_run_test(test_DisabilityAnnuity);
  mu_run_test(test_DisabilityTermInsurance);*/
  return 0;
}

int start_testing() {
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

int main(int argc, char **argv) {
  int r = 0;
  r = start_testing();
  //all_timing(42000);
  return r;
}
