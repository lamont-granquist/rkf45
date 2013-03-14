//
//  RK_Library.c
//  RK_C
//
//  Created by Nicolai Dahl on 31/01/12.
//  Copyright (c) 2012 ITU. All rights reserved.
//

#include <stdio.h>
#include "RK_Library.h"
#include <math.h>

double age = 30;
double interestrate = 0.05;
double bpension = 1;
double pensiontime = 35;

double GM(double t) {
    return 0.0005 + pow(10.0, 5.728 - 10.0 + 0.038*(age + t));
}

// Interest
double r(double t) {
    return interestrate;
}

double indicator(int b) {
    return b ? 1.0 : 0.0;
}
