#ifndef RK_C_DisabilityTermInsurance_h
#define RK_C_DisabilityTermInsurance_h

float** tv_DisabilityTermInsurance();
void dy_DisabilityTermInsurance(float t,float* V,float* result);
void bj_ii_DisabilityTermInsurance(float t, float* result);

#endif
