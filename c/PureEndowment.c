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
#include "RK_Library.h"
#include "PureEndowment.h"

static double *matrix[41];

static double row0[]  = {0.1437946974886250};
static double row1[]  = {0.1513594875720590};
static double row2[]  = {0.1593334830357050};
static double row3[]  = {0.1677404776222450};
static double row4[]  = {0.1766058889630630};
static double row5[]  = {0.1859569023966960};
static double row6[]  = {0.1958226315281160};
static double row7[]  = {0.2062342978884670};
static double row8[]  = {0.2172254324285080};
static double row9[]  = {0.2288321020171980};
static double row10[] = {0.2410931646317720};
static double row11[] = {0.2540505575320670};
static double row12[] = {0.2677496234274150};
static double row13[] = {0.2822394804908960};
static double row14[] = {0.2975734430792770};
static double row15[] = {0.3138095012097790};
static double row16[] = {0.3310108682663330};
static double row17[] = {0.3492466081065300};
static double row18[] = {0.3685923547759590};
static double row19[] = {0.3891311404830350};
static double row20[] = {0.4109543504366150};
static double row21[] = {0.4341628267168120};
static double row22[] = {0.4588681476758340};
static double row23[] = {0.4851941146396860};
static double row24[] = {0.5132784841210240};
static double row25[] = {0.5432749916555080};
static double row26[] = {0.5753557231029890};
static double row27[] = {0.6097139012822690};
static double row28[] = {0.6465671707381610};
static double row29[] = {0.6861614820516400};
static double row30[] = {0.7287757004081980};
static double row31[] = {0.7747270924521970};
static double row32[] = {0.8243778824987340};
static double row33[] = {0.8781431162166330};
static double row34[] = {0.9365001299367070};
static double row35[] = {0.0000000000000000};
static double row36[] = {0.0000000000000000};
static double row37[] = {0.0000000000000000};
static double row38[] = {0.0000000000000000};
static double row39[] = {0.0000000000000000};
static double row40[] = {0.0000000000000000};

/* tv */
double** tv_PureEndowment() {

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
static double bdeath = 1;

static double b_0(double t) {
    return 0.0;
}

static double mu_01(double t) {
    return GM(t);
}

static double bj_00(double t) {
    return t == pensiontime ? bpension: 0.0;
}

static double bj_01(double t) {
    return 0.0; 
}

void bj_ii_PureEndowment(double t, double* result) {
  result[0] += bj_00(t);
}

void dy_PureEndowment(double t, double* V,double* result)
{
    result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
