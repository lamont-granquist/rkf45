//
//  TermInsurance.c
//  RK_C
//
//  Created by Nicolai Dahl on 05/02/12.
//  Copyright (c) 2012 ITU. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "RK_Library.h"
#include "TermInsurance.h"


static int n = 35;
static double bdeath = 1;

static double b_0(double t) {
    return 0.0;
}

static double mu_01(double t) {
    return GM(t);
}

static double bj_00(double t) {
    return 0.0;
}

static double bj_01(double t) {
    return bdeath * indicator(t > 0) * indicator(t < n);
}

void bj_ii_TermInsurance(double t, double* result) {
  result[0] = bj_00(t);
}

void dy_TermInsurance(double t, double* V,double* result)
{
    result[0] = 5.0; //r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
