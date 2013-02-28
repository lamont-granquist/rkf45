//
//  RK_Library.h
//  RK_C
//
//  Created by Nicolai Dahl on 31/01/12.
//  Copyright (c) 2012 ITU. All rights reserved.
//



#ifndef RK_C_RK_Library_h
#define RK_C_RK_Library_h


extern float age, interestrate, bpension, pensiontime;
extern int states;


float GM(float t);
float r(float t);
float indicator(int b);


void sax(float a, float* x, float* result);
void saxpy(float a, float* x, float* y, float* result);

#endif
