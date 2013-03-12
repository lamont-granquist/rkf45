#ifndef RK_C_DisabilityAnnuity_h
#define RK_C_DisabilityAnnuity_h

double** tv_DisabilityAnnuity();
void dy_DisabilityAnnuity(double t,double* V,double* result);
void bj_ii_DisabilityAnnuity(double t, double* result);

#endif
