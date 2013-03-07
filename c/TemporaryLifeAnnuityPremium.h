#ifndef RK_C_TemporaryLifeAnnuityPremium_h
#define RK_C_TemporaryLifeAnnuityPremium_h

double** tv_TemporaryLifeAnnuityPremium();
void dy_TemporaryLifeAnnuityPremium(double t,double* V,double* result);
void bj_ii_TemporaryLifeAnnuityPremium(double t, double* result);

#endif
