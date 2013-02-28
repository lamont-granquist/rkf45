//
//  RK_Library.h
//  RK_C
//
//  Created by Nicolai Dahl on 31/01/12.
//  Copyright (c) 2012 ITU. All rights reserved.
//



#ifndef RK_C_RK_Library_h
#define RK_C_RK_Library_h


extern double age, interestrate, bpension, pensiontime;
extern int states;


double GM(double t);
double r(double t);
double indicator(int b);


void sax(double a, double* x, double* result);
void saxpy(double a, double* x, double* y, double* result);

#endif
