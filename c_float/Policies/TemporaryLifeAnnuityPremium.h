#ifndef RK_C_TemporaryLifeAnnuityPremium_h
#define RK_C_TemporaryLifeAnnuityPremium_h

float** tv_TemporaryLifeAnnuityPremium();
void dy_TemporaryLifeAnnuityPremium(float t,float* V,float* result);
void bj_ii_TemporaryLifeAnnuityPremium(float t, float* result);

#endif
