//
//  PureEndowment.c
//  RK_C
//
//  Created by Nicolai Dahl on 05/02/12.
//  Copyright (c) 2012 ITU. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../RK_Library.h"
#include "PureEndowment.h"

static float *matrix[41];

static float row0[]  = {0.1437959f};
static float row1[]  = {0.1513607f};
static float row2[]  = {0.1593348f};
static float row3[]  = {0.1677418f};
static float row4[]  = {0.1766073f};
static float row5[]  = {0.1859584f};
static float row6[]  = {0.1958242f};
static float row7[]  = {0.2062360f};
static float row8[]  = {0.2172272f};
static float row9[]  = {0.2288340f};
static float row10[] = {0.2410951f};
static float row11[] = {0.2540526f};
static float row12[] = {0.2677518f};
static float row13[] = {0.2822418f};
static float row14[] = {0.2975759f};
static float row15[] = {0.3138121f};
static float row16[] = {0.3310136f};
static float row17[] = {0.3492495f};
static float row18[] = {0.3685954f};
static float row19[] = {0.3891343f};
static float row20[] = {0.4109577f};
static float row21[] = {0.4341664f};
static float row22[] = {0.4588719f};
static float row23[] = {0.4851981f};
static float row24[] = {0.5132827f};
static float row25[] = {0.5432795f};
static float row26[] = {0.5753605f};
static float row27[] = {0.6097189f};
static float row28[] = {0.6465725f};
static float row29[] = {0.6861671f};
static float row30[] = {0.7287817f};
static float row31[] = {0.7747335f};
static float row32[] = {0.8243847f};
static float row33[] = {0.8781503f};
static float row34[] = {0.9365078f};
static float row35[] = {0.0000000f};
static float row36[] = {0.0000000f};
static float row37[] = {0.0000000f};
static float row38[] = {0.0000000f};
static float row39[] = {0.0000000f};
static float row40[] = {0.0000000f};

/* tv */
float** tv_PureEndowment() {

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
  
  return matrix;
}

static int n = 35;
static float bdeath = 1.0f;

static float b_0(float t) {
    return 0.0f;
}

static float mu_01(float t) {
    return GM(t);
}

static float bj_00(float t) {
    return t == pensiontime ? bpension: 0.0f;
}

static float bj_01(float t) {
    return 0.0f; 
}

void bj_ii_PureEndowment(float t, float* result) {
  result[0] += bj_00(t);
}

void dy_PureEndowment(float t, float* V,float* result)
{
    result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
