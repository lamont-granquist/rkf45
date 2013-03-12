#ifndef RK_C_DeferredTemporaryLifeAnnuity_h
#define RK_C_DeferredTemporaryLifeAnnuity_h

double** tv_DeferredTemporaryLifeAnnuity();
void dy_DeferredTemporaryLifeAnnuity(double t,double* V,double* result);
void bj_ii_DeferredTemporaryLifeAnnuity(double t, double* result);

#endif
