#ifndef Rkf45_h
#define Rkf45_h

extern int neqn;
extern int start_year;
extern int end_year;
// dy
// bj_ii
extern float relerr;
extern float abserr;
extern float* end_year_y; 

float** estimate();
void construct(int n);

#endif
