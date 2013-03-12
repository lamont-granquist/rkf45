#ifndef RK_C_DisabilityTermInsurance_h
#define RK_C_DisabilityTermInsurance_h

double** tv_DisabilityTermInsurance();
void dy_DisabilityTermInsurance(double t,double* V,double* result);
void bj_ii_DisabilityTermInsurance(double t, double* result);

#endif
