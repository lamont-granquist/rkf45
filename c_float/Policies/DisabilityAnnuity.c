//
//  DisabilityAnnuity.c
//  RK_C
//
//  Created by Nicolai Dahl on 05/02/12.
//  Copyright (c) 2012 ITU. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../RK_Library.h"
#include "DisabilityAnnuity.h"
#include <math.h>

static float *matrix[51];

static float row0[]  = {0.55552600f, 15.97177000f};
static float row1[]  = {0.56979390f, 15.78593000f};
static float row2[]  = {0.58424580f, 15.59145000f};
static float row3[]  = {0.59881110f, 15.38794000f};
static float row4[]  = {0.61340610f, 15.17502000f};
static float row5[]  = {0.62793190f, 14.95226000f};
static float row6[]  = {0.64227370f, 14.71923000f};
static float row7[]  = {0.65629850f, 14.47547000f};
static float row8[]  = {0.66985320f, 14.22050000f};
static float row9[]  = {0.68276310f, 13.95383000f};
static float row10[] = {0.69483000f, 13.67492000f};
static float row11[] = {0.70583020f, 13.38322000f};
static float row12[] = {0.71551300f, 13.07813000f};
static float row13[] = {0.72359900f, 12.75905000f};
static float row14[] = {0.72977920f, 12.42530000f};
static float row15[] = {0.73371460f, 12.07619000f};
static float row16[] = {0.73503580f, 11.71095000f};
static float row17[] = {0.73334430f, 11.32880000f};
static float row18[] = {0.72821520f, 10.92886000f};
static float row19[] = {0.71920140f, 10.51020000f};
static float row20[] = {0.70584070f, 10.07182000f};
static float row21[] = {0.68766540f,  9.61260800f};
static float row22[] = {0.66421740f,  9.13137200f};
static float row23[] = {0.63506820f,  8.62679500f};
static float row24[] = {0.59984780f,  8.09742500f};
static float row25[] = {0.55828260f,  7.54165500f};
static float row26[] = {0.51024840f,  6.95770000f};
static float row27[] = {0.45584220f,  6.34356300f};
static float row28[] = {0.39547940f,  5.69700300f};
static float row29[] = {0.33002640f,  5.01549200f};
static float row30[] = {0.26098230f,  4.29616600f};
static float row31[] = {0.19072930f,  3.53576700f};
static float row32[] = {0.12287920f,  2.73056700f};
static float row33[] = {0.06275875f,  1.87629100f};
static float row34[] = {0.01809598f,  0.96800480f};
static float row35[] = {0.00000000f,  0.00000000f};
static float row36[] = {0.00000000f,  0.00000000f};
static float row37[] = {0.00000000f,  0.00000000f};
static float row38[] = {0.00000000f,  0.00000000f};
static float row39[] = {0.00000000f,  0.00000000f};
static float row40[] = {0.00000000f,  0.00000000f};
static float row41[] = {0.00000000f,  0.00000000f};
static float row42[] = {0.00000000f,  0.00000000f};
static float row43[] = {0.00000000f,  0.00000000f};
static float row44[] = {0.00000000f,  0.00000000f};
static float row45[] = {0.00000000f,  0.00000000f};
static float row46[] = {0.00000000f,  0.00000000f};
static float row47[] = {0.00000000f,  0.00000000f};
static float row48[] = {0.00000000f,  0.00000000f};
static float row49[] = {0.00000000f,  0.00000000f};
static float row50[] = {0.00000000f,  0.00000000f};

/* tv */
float** tv_DisabilityAnnuity() {

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
static int bdisabled = 1;

static float b_0(float t) {
  return 0.0f;
}

static float b_1(float t) {
  //printf("b_1: %.16f t: %.16f\n",bdisabled * indicator(t > 0) * indicator(t < n),t);
  return bdisabled * indicator(t > 0) * indicator(t < n);
}

static float GM01(float t) {
  return 0.0006f + powf(10.0f, 4.71609f - 10.0f + 0.06f*(age + t));
}

static float GM02(float t) {
  return GM(t);
}

static float GM12(float t) {
  return GM(t);
}

static float mu_01(float t) {
  return GM01(t);
}

static float mu_02(float t) {
  return GM02(t);
}

static float mu_12(float t) {
  return GM12(t);
}

static float bj_00(float t) {
  return 0.0f;
}

static float bj_01(float t) {
  return 0.0f;
}

static float bj_02(float t) {
  return 0.0f;
}

static float bj_11(float t) {
  return 0.0f;
}

static float bj_12(float t) {
  return 0.0f;
}

void bj_ii_DisabilityAnnuity(float t, float* result) {
  result[0] += bj_00(t); // 0.0
  result[1] += bj_11(t); // 0.0
}

void dy_DisabilityAnnuity(float t, float* V,float* result)
{
  result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (V[1] - V[0] + bj_01(t)) - mu_02(t) * (0 - V[0] + bj_02(t));
  //result[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t)); 
}
