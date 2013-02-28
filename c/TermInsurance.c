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
static float bdeath = 1;

static float b_0(float t) {
    return 0.0;
}

static float mu_01(float t) {
    return GM(t);
}

static void bj_00(float t, float *result) {
    result[0] = 0.0;
}

static float bj_01(float t) {
    return bdeath * indicator(t > 0) * indicator(t < n);
}

void dy_TermInsurance(double t, double* V,int m, int n, double result[m][n])
{
    result[0][0] = 5.0; //r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
