//
//  DeferredTemporaryLifeAnnuity.c
//  RK_C
//
//  Created by Nicolai Dahl on 05/02/12.
//  Copyright (c) 2012 ITU. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../RK_Library.h"
#include "DeferredTemporaryLifeAnnuity.h"

static float *matrix[51];

static float row0[]  = {1.026561f};
static float row1[]  = {1.080566f};
static float row2[]  = {1.137493f};
static float row3[]  = {1.197511f};
static float row4[]  = {1.260802f};
static float row5[]  = {1.327560f};
static float row6[]  = {1.397992f};
static float row7[]  = {1.472322f};
static float row8[]  = {1.550788f};
static float row9[]  = {1.633649f};
static float row10[] = {1.721182f};
static float row11[] = {1.813685f};
static float row12[] = {1.911484f};
static float row13[] = {2.014928f};
static float row14[] = {2.124398f};
static float row15[] = {2.240309f};
static float row16[] = {2.363111f};
static float row17[] = {2.493297f};
static float row18[] = {2.631407f};
static float row19[] = {2.778035f};
static float row20[] = {2.933833f};
static float row21[] = {3.099520f};
static float row22[] = {3.275893f};
static float row23[] = {3.463836f};
static float row24[] = {3.664332f};
static float row25[] = {3.878480f};
static float row26[] = {4.107506f};
static float row27[] = {4.352792f};
static float row28[] = {4.615890f};
static float row29[] = {4.898557f};
static float row30[] = {5.202783f};
static float row31[] = {5.530833f};
static float row32[] = {5.885293f};
static float row33[] = {6.269127f};
static float row34[] = {6.685742f};
static float row35[] = {7.139072f};
static float row36[] = {6.599296f};
static float row37[] = {6.031972f};
static float row38[] = {5.434192f};
static float row39[] = {4.802500f};
static float row40[] = {4.132774f};
static float row41[] = {3.420072f};
static float row42[] = {2.658446f};
static float row43[] = {1.840702f};
static float row44[] = {0.958105f};
static float row45[] = {0.000000f};
static float row46[] = {0.000000f};
static float row47[] = {0.000000f};
static float row48[] = {0.000000f};
static float row49[] = {0.000000f};
static float row50[] = {0.000000f};

/* tv */
float** tv_DeferredTemporaryLifeAnnuity() {

  matrix[0] = row0;
  matrix[1] = row1;
  matrix[2] = row2;
  matrix[3] = row3;
  matrix[4] = row4;
  matrix[5] = row5;
  matrix[6] = row6;
  matrix[7] = row7;
  matrix[8] = row8;
  matrix[9] = row9;
  matrix[10] = row10;
  matrix[11] = row11;
  matrix[12] = row12;
  matrix[13] = row13;
  matrix[14] = row14;
  matrix[15] = row15;
  matrix[16] = row16;
  matrix[17] = row17;
  matrix[18] = row18;
  matrix[19] = row19;
  matrix[20] = row20;
  matrix[21] = row21;
  matrix[22] = row22;
  matrix[23] = row23;
  matrix[24] = row24;
  matrix[25] = row25;
  matrix[26] = row26;
  matrix[27] = row27;
  matrix[28] = row28;
  matrix[29] = row29;
  matrix[30] = row30;
  matrix[31] = row31;
  matrix[32] = row32;
  matrix[33] = row33;
  matrix[34] = row34;
  matrix[35] = row35;
  matrix[36] = row36;
  matrix[37] = row37;
  matrix[38] = row38;
  matrix[39] = row39;
  matrix[40] = row40;
  matrix[41] = row41;
  matrix[42] = row42;
  matrix[43] = row43;
  matrix[44] = row44;
  matrix[45] = row45;
  matrix[46] = row46;
  matrix[47] = row47;
  matrix[48] = row48;
  matrix[49] = row49;
  matrix[50] = row50;
  
  return matrix;
}

static int n = 10;
static int m = 35;
static float bdeath = 1.0f;

static float b_0(float t) {
    return bpension * indicator(t > m) * indicator(t < m + n);
}

static float mu_01(float t) {
    return GM(t);
}

static float bj_00(float t) {
    return 0.0f;
}

static float bj_01(float t) {
    return 0.0f;
}

void bj_ii_DeferredTemporaryLifeAnnuity(float t, float* result) {
  result[0] += bj_00(t);
}

void dy_DeferredTemporaryLifeAnnuity(float t, float* V,float* result)
{
    result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
