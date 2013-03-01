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

static double *matrix[51];

static double row0[] = {0.0576169193132673};
static double row1[] = {0.0593440338503396};
static double row2[] = {0.0610940381466785};
static double row3[] = {0.0628621872268886};
static double row4[] = {0.0646429589229975};
static double row5[] = {0.0664299642301075};
static double row6[] = {0.0682158480098718};
static double row7[] = {0.0699921788564634};
static double row8[] = {0.0717493268311646};
static double row9[] = {0.0734763275934875};
static double row10[] = {0.0751607312303756};
static double row11[] = {0.0767884338351089};
static double row12[] = {0.0783434895820515};
static double row13[] = {0.0798079006843388};
static double row14[] = {0.0811613821939618};
static double row15[] = {0.0823810980933246};
static double row16[] = {0.0834413645163828};
static double row17[] = {0.0843133152038705};
static double row18[] = {0.0849645234136384};
static double row19[] = {0.0853585734399381};
static double row20[] = {0.0854545736024947};
static double row21[] = {0.0852066009948634};
static double row22[] = {0.0845630663660149};
static double row23[] = {0.0834659851665421};
static double row24[] = {0.0818501379168519};
static double row25[] = {0.0796420995167974};
static double row26[] = {0.0767591127460391};
static double row27[] = {0.0731077757869581};
static double row28[] = {0.0685825068600615};
static double row29[] = {0.0630637406431617};
static double row30[] = {0.0564158005824136};
static double row31[] = {0.0484843779035996};
static double row32[] = {0.0390935313045591};
static double row33[] = {0.0280420999246517};
static double row34[] = {0.0150993948779494};
static double row35[] = {0.0000000000000000};
static double row36[] = {0.0000000000000000};
static double row37[] = {0.0000000000000000};
static double row38[] = {0.0000000000000000};
static double row39[] = {0.0000000000000000};
static double row40[] = {0.0000000000000000};
static double row41[] = {0.0000000000000000};
static double row42[] = {0.0000000000000000};
static double row43[] = {0.0000000000000000};
static double row44[] = {0.0000000000000000};
static double row45[] = {0.0000000000000000};
static double row46[] = {0.0000000000000000};
static double row47[] = {0.0000000000000000};
static double row48[] = {0.0000000000000000};
static double row49[] = {0.0000000000000000};
static double row50[] = {0.0000000000000000};

/* tv */
double** tv_TermInsurance() {

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
    result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
