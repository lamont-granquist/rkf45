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

static double *matrix[51];

static double row0[]  = {0.5555261079604120,   15.9717676673750000};
static double row1[]  = {0.5697939362470290,   15.7859295725873000};
static double row2[]  = {0.5842458860490700,   15.5914495774393000};
static double row3[]  = {0.5988112559939020,   15.3879467041578000};
static double row4[]  = {0.6134061668653150,   15.1750230434743000};
static double row5[]  = {0.6279320187147260,   14.9522626687161000};
static double row6[]  = {0.6422738476905750,   14.7192304105506000};
static double row7[]  = {0.6562985942478250,   14.4754704664814000};
static double row8[]  = {0.6698533011953220,   14.2205048162601000};
static double row9[]  = {0.6827632688462940,   13.9538314093175000};
static double row10[] = {0.6948302058887720,   13.6749220843703000};
static double row11[] = {0.7058304291697920,   13.3832201743564000};
static double row12[] = {0.7155131842695250,   13.0781377416023000};
static double row13[] = {0.7235991826741960,   12.7590523783855000};
static double row14[] = {0.7297794820426480,   12.4253034965300000};
static double row15[] = {0.7337148754958810,   12.0761880160438000};
static double row16[] = {0.7350360067043140,   11.7109553465694000};
static double row17[] = {0.7333444934112320,   11.3288015361409000};
static double row18[] = {0.7282154278216300,   10.9288624386901000};
static double row19[] = {0.7192017347676550,   10.5102057241642000};
static double row20[] = {0.7058410171284570,   10.0718215219888000};
static double row21[] = {0.6876657158018430,    9.6126114486939200};
static double row22[] = {0.6642176772198330,    9.1313757222565700};
static double row23[] = {0.6350685815395070,    8.6267980071455300};
static double row24[] = {0.5998481774680660,    8.0974275627005400};
static double row25[] = {0.5582829507417500,    7.5416581802122700};
static double row26[] = {0.5102488039940220,    6.9577032868915200};
static double row27[] = {0.4558426666237880,    6.3435664627186100};
static double row28[] = {0.3954798644641650,    5.6970064523866500};
static double row29[] = {0.3300268327704200,    5.0154955507223200};
static double row30[] = {0.2609827682846260,    4.2961699850940100};
static double row31[] = {0.1907297303473690,    3.5357705981010400};
static double row32[] = {0.1228795253423570,    2.7305717294870500};
static double row33[] = {0.0627590447157113,    1.8762956830912200};
static double row34[] = {0.0180961562710709,    0.9680095100190450};
static double row35[] = {0.0000000000000000,    0.0000000000000000};
static double row36[] = {0.0000000000000000,    0.0000000000000000};
static double row37[] = {0.0000000000000000,    0.0000000000000000};
static double row38[] = {0.0000000000000000,    0.0000000000000000};
static double row39[] = {0.0000000000000000,    0.0000000000000000};
static double row40[] = {0.0000000000000000,    0.0000000000000000};
static double row41[] = {0.0000000000000000,    0.0000000000000000};
static double row42[] = {0.0000000000000000,    0.0000000000000000};
static double row43[] = {0.0000000000000000,    0.0000000000000000};
static double row44[] = {0.0000000000000000,    0.0000000000000000};
static double row45[] = {0.0000000000000000,    0.0000000000000000};
static double row46[] = {0.0000000000000000,    0.0000000000000000};
static double row47[] = {0.0000000000000000,    0.0000000000000000};
static double row48[] = {0.0000000000000000,    0.0000000000000000};
static double row49[] = {0.0000000000000000,    0.0000000000000000};
static double row50[] = {0.0000000000000000,    0.0000000000000000};

/* tv */
double** tv_DisabilityAnnuity() {

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

static double b_0(double t) {
  return 0.0;
}

static double b_1(double t) {
  return bdisabled * indicator(t > 0) * indicator(t < n);
}

static double GM01(double t) {
  return 0.0006 + pow(10, 4.71609 - 10 + 0.06*(age + t));
}

static double GM02(double t) {
  return GM(t);
}

static double GM12(double t) {
  return GM(t);
}

static double mu_01(double t) {
  return GM01(t);
}

static double mu_02(double t) {
  return GM02(t);
}

static double mu_12(double t) {
  return GM12(t);
}

static double bj_00(double t) {
  return 0.0;
}

static double bj_01(double t) {
  return 0.0;
}

static double bj_02(double t) {
  return 0.0;
}

static double bj_11(double t) {
  return 0.0;
}

static double bj_12(double t) {
  return 0.0;
}

void bj_ii_DisabilityAnnuity(double t, double* result) {
  result[0] += bj_00(t);
  result[1] += bj_11(t);
}

void dy_DisabilityAnnuity(double t, double* V,double* result)
{
  result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
  result[1] = r(t) * V[1] - b_1(t) - mu_12(t) * (0 - V[1] + bj_12(t));
}
