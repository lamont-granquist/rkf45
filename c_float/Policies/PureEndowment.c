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

static float row0[]  = {0.1437946974886250f};
static float row1[]  = {0.1513594875720590f};
static float row2[]  = {0.1593334830357050f};
static float row3[]  = {0.1677404776222450f};
static float row4[]  = {0.1766058889630630f};
static float row5[]  = {0.1859569023966960f};
static float row6[]  = {0.1958226315281160f};
static float row7[]  = {0.2062342978884670f};
static float row8[]  = {0.2172254324285080f};
static float row9[]  = {0.2288321020171980f};
static float row10[] = {0.2410931646317720f};
static float row11[] = {0.2540505575320670f};
static float row12[] = {0.2677496234274150f};
static float row13[] = {0.2822394804908960f};
static float row14[] = {0.2975734430792770f};
static float row15[] = {0.3138095012097790f};
static float row16[] = {0.3310108682663330f};
static float row17[] = {0.3492466081065300f};
static float row18[] = {0.3685923547759590f};
static float row19[] = {0.3891311404830350f};
static float row20[] = {0.4109543504366150f};
static float row21[] = {0.4341628267168120f};
static float row22[] = {0.4588681476758340f};
static float row23[] = {0.4851941146396860f};
static float row24[] = {0.5132784841210240f};
static float row25[] = {0.5432749916555080f};
static float row26[] = {0.5753557231029890f};
static float row27[] = {0.6097139012822690f};
static float row28[] = {0.6465671707381610f};
static float row29[] = {0.6861614820516400f};
static float row30[] = {0.7287757004081980f};
static float row31[] = {0.7747270924521970f};
static float row32[] = {0.8243778824987340f};
static float row33[] = {0.8781431162166330f};
static float row34[] = {0.9365001299367070f};
static float row35[] = {0.0000000000000000f};
static float row36[] = {0.0000000000000000f};
static float row37[] = {0.0000000000000000f};
static float row38[] = {0.0000000000000000f};
static float row39[] = {0.0000000000000000f};
static float row40[] = {0.0000000000000000f};

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
static float bdeath = 1;

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
