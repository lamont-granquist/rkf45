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
#include "../RK_Library.h"
#include "TermInsurance.h"

static float *matrix[51];

static float row0[] =  {0.05761727f};
static float row1[] =  {0.05934345f};
static float row2[] =  {0.06109342f};
static float row3[] =  {0.06286154f};
static float row4[] =  {0.06464227f};
static float row5[] =  {0.06642924f};
static float row6[] =  {0.06821509f};
static float row7[] =  {0.06999138f};
static float row8[] =  {0.07174849f};
static float row9[] =  {0.07347544f};
static float row10[] = {0.07515980f};
static float row11[] = {0.07678745f};
static float row12[] = {0.07834245f};
static float row13[] = {0.07980680f};
static float row14[] = {0.08116022f};
static float row15[] = {0.08237988f};
static float row16[] = {0.08344008f};
static float row17[] = {0.08431196f};
static float row18[] = {0.08496310f};
static float row19[] = {0.08535707f};
static float row20[] = {0.08545298f};
static float row21[] = {0.08520492f};
static float row22[] = {0.08456129f};
static float row23[] = {0.08346410f};
static float row24[] = {0.08184814f};
static float row25[] = {0.07963999f};
static float row26[] = {0.07675687f};
static float row27[] = {0.07310540f};
static float row28[] = {0.06857999f};
static float row29[] = {0.06306107f};
static float row30[] = {0.05641296f};
static float row31[] = {0.04848135f};
static float row32[] = {0.03909031f};
static float row33[] = {0.02803866f};
static float row34[] = {0.01509573f};
static float row35[] = {0.00000000f};
static float row36[] = {0.00000000f};
static float row37[] = {0.00000000f};
static float row38[] = {0.00000000f};
static float row39[] = {0.00000000f};
static float row40[] = {0.00000000f};
static float row41[] = {0.00000000f};
static float row42[] = {0.00000000f};
static float row43[] = {0.00000000f};
static float row44[] = {0.00000000f};
static float row45[] = {0.00000000f};
static float row46[] = {0.00000000f};
static float row47[] = {0.00000000f};
static float row48[] = {0.00000000f};
static float row49[] = {0.00000000f};
static float row50[] = {0.00000000f};

/* tv */
float** tv_TermInsurance() {

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

static int n = 35;
static float bdeath = 1.0f;

static float b_0(float t) {
    return 0.0f;
}

static float mu_01(float t) {
    return GM(t);
}

static float bj_00(float t) {
    return 0.0f;
}

static float bj_01(float t) {
    return bdeath * indicator(t > 0) * indicator(t < n);
}

void bj_ii_TermInsurance(float t, float* result) {
  result[0] += bj_00(t);
}

void dy_TermInsurance(float t, float* V,float* result)
{
    result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
