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


float age = 30.0f;
float interestrate = 0.05f;
float bpension = 1.0f;
float pensiontime = 35.0f;
int states = 0;


float GM(float t) {
    return 0.0005f + pow(10.0f, 5.728f - 10.0f + 0.038f*(age + t));
}

static float ts[] = { 
    0.25, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
    15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 };

static float rs[] = { 
    1.146677033, 1.146677033, 1.146677033, 1.340669678, 1.571952911, 1.803236144, 
    2.034519377, 2.26580261, 2.497085843, 2.584085843, 2.710085843, 2.805085843, 
    2.871485843, 2.937885843, 3.004285843, 3.070685843, 3.137085843, 3.136485843, 
    3.135885843, 3.135285843, 3.134685843, 3.134085843, 3.113185843, 3.092285843, 
    3.071385843, 3.050485843, 3.029585843, 3.008685843, 2.987785843, 2.966885843, 
    2.945985843, 2.925085843
};


static float rFsa(float t) { 
    int tslen = sizeof(ts)/sizeof(float);
    
    // Requires ts non-empty and elements strictly increasing.
    int last = tslen-1;
    if (t <= ts[0])
        return log(1 + rs[0]/100);
    else if (t >= ts[last])
        return log(1 + rs[last]/100);
    else {
        int a = 0, b = last;
        // Now a < b (bcs. ts must have more than 1 element) and ts[a] < t < ts[b]
        while (a+1 < b) {
            // Now a < b and ts[a] <= t < ts[b]
            int i = (a+b)/2;
            if (ts[i] <= t)
                a = i;
            else // t < ts[i]
                b = i;
        }
        // Now a+1>=b and ts[a] <= t < ts[b]; so a!=b and hence a+1 == b <= last
        int m = a;
        float tm = ts[m], tm1 = ts[m+1];
        float rm = rs[m] / 100, rm1 = rs[m+1] / 100;
        float Rt = (rm * (tm1 - t) + rm1 * (t - tm)) / (tm1 - tm);
        return log(1 + Rt) + t / (tm1 - tm) * (rm1 - rm) / (1 + Rt);
    }
}

// Interest
float r(float t) {
    //return interestrate;
    return rFsa(t);
}

float indicator(int b) {
    return b ? 1.0 : 0.0;
}


// sax = scalar a times x array
void sax(float a, float* x, float* result) {
    for (int i=0; i<states; i++, x++)
		result[i] = a * *x;
}

// saxpy = scalar a times x array plus y array
void saxpy(float a, float* x, float* y, float* result) {
    for (int i=0; i<states; i++, x++, y++)
		result[i] = a * *x + *y;
    
}
