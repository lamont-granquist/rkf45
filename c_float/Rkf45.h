#ifndef Rkf45_h
#define Rkf45_h

extern const int MAX_NEQN;
// dy
// bj_ii
extern float* end_year_y; 

float** estimate();
void construct(int n);

#endif
