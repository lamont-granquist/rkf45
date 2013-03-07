#ifndef RK_C_TermInsurance_h
#define RK_C_TermInsurance_h

double** tv_TermInsurance();
void dy_TermInsurance(double t,double* V,double* result);
void bj_ii_TermInsurance(double t, double* result);

#endif
