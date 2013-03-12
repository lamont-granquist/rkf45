//
//  TemporaryLifeAnnuityPremium.c
//  RK_C
//
//  Created by Nicolai Dahl on 05/02/12.
//  Copyright (c) 2012 ITU. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../RK_Library.h"
#include "TemporaryLifeAnnuityPremium.h"

static float *matrix[51];

static float row0[] =  {-15.97177f};
static float row1[] =  {-15.78593f};
static float row2[] =  {-15.59145f};
static float row3[] =  {-15.38794f};
static float row4[] =  {-15.17502f};
static float row5[] =  {-14.95226f};
static float row6[] =  {-14.71923f};
static float row7[] =  {-14.47547f};
static float row8[] =  {-14.22050f};
static float row9[] =  {-13.95383f};
static float row10[] = {-13.67492f};
static float row11[] = {-13.38322f};
static float row12[] = {-13.07813f};
static float row13[] = {-12.75905f};
static float row14[] = {-12.42530f};
static float row15[] = {-12.07619f};
static float row16[] = {-11.71095f};
static float row17[] = {-11.32880f};
static float row18[] = {-10.92886f};
static float row19[] = {-10.51020f};
static float row20[] = {-10.07182f};
static float row21[] = {-9.612608f};
static float row22[] = {-9.131372f};
static float row23[] = {-8.626795f};
static float row24[] = {-8.097425f};
static float row25[] = {-7.541655f};
static float row26[] = {-6.957700f};
static float row27[] = {-6.343563f};
static float row28[] = {-5.697003f};
static float row29[] = {-5.015492f};
static float row30[] = {-4.296166f};
static float row31[] = {-3.535767f};
static float row32[] = {-2.730567f};
static float row33[] = {-1.876291f};
static float row34[] = {-0.968004f};
static float row35[] = {0.0000000f};
static float row36[] = {0.0000000f};
static float row37[] = {0.0000000f};
static float row38[] = {0.0000000f};
static float row39[] = {0.0000000f};
static float row40[] = {0.0000000f};
static float row41[] = {0.0000000f};
static float row42[] = {0.0000000f};
static float row43[] = {0.0000000f};
static float row44[] = {0.0000000f};
static float row45[] = {0.0000000f};
static float row46[] = {0.0000000f};
static float row47[] = {0.0000000f};
static float row48[] = {0.0000000f};
static float row49[] = {0.0000000f};
static float row50[] = {0.0000000f};

/* tv */
float** tv_TemporaryLifeAnnuityPremium() {

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
static float bpremium = 1.0f;

static float b_0(float t) {
    return -bpremium * indicator(t >= 0) * indicator(t < n);
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

void bj_ii_TemporaryLifeAnnuityPremium(float t, float* result) {
  result[0] += bj_00(t);
}

void dy_TemporaryLifeAnnuityPremium(float t, float* V,float* result)
{
    result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
