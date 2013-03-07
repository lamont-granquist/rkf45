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
#include "RK_Library.h"
#include "TemporaryLifeAnnuityPremium.h"

static double *matrix[51];

static double row0[] =  {-15.9717676660001000};
static double row1[] =  {-15.7859295725898000};
static double row2[] =  {-15.5914495774420000};
static double row3[] =  {-15.3879467041606000};
static double row4[] =  {-15.1750230434772000};
static double row5[] =  {-14.9522626687192000};
static double row6[] =  {-14.7192304105538000};
static double row7[] =  {-14.4754704664848000};
static double row8[] =  {-14.2205048162637000};
static double row9[] =  {-13.9538314093213000};
static double row10[] = {-13.6749220843743000};
static double row11[] = {-13.3832201743607000};
static double row12[] = {-13.0781377416068000};
static double row13[] = {-12.7590523783886000};
static double row14[] = {-12.4253034965330000};
static double row15[] = {-12.0761880160465000};
static double row16[] = {-11.7109553465719000};
static double row17[] = {-11.3288015361432000};
static double row18[] = {-10.9288624386922000};
static double row19[] = {-10.5102057241661000};
static double row20[] = {-10.0718215219905000};
static double row21[] = { -9.6126114486954800};
static double row22[] = { -9.1313757222581100};
static double row23[] = { -8.6267980071471600};
static double row24[] = { -8.0974275627022600};
static double row25[] = { -7.5416581802140900};
static double row26[] = { -6.9577032868934500};
static double row27[] = { -6.3435664627206500};
static double row28[] = { -5.6970064523888100};
static double row29[] = { -5.0154955507238000};
static double row30[] = { -4.2961699850952200};
static double row31[] = { -3.5357705981019300};
static double row32[] = { -2.7305717294876500};
static double row33[] = { -1.8762956830915000};
static double row34[] = { -0.9680095100191570};
static double row35[] = {  0.0000000000000000};
static double row36[] = {  0.0000000000000000};
static double row37[] = {  0.0000000000000000};
static double row38[] = {  0.0000000000000000};
static double row39[] = {  0.0000000000000000};
static double row40[] = {  0.0000000000000000};
static double row41[] = {  0.0000000000000000};
static double row42[] = {  0.0000000000000000};
static double row43[] = {  0.0000000000000000};
static double row44[] = {  0.0000000000000000};
static double row45[] = {  0.0000000000000000};
static double row46[] = {  0.0000000000000000};
static double row47[] = {  0.0000000000000000};
static double row48[] = {  0.0000000000000000};
static double row49[] = {  0.0000000000000000};
static double row50[] = {  0.0000000000000000};

/* tv */
double** tv_TemporaryLifeAnnuityPremium() {

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
static double bpremium = 1;

static double b_0(double t) {
    return -bpremium * indicator(t >= 0) * indicator(t < n);
}

static double mu_01(double t) {
    return GM(t);
}

static double bj_00(double t) {
    return 0.0;
}

static double bj_01(double t) {
    return 0.0; 
}

void bj_ii_TemporaryLifeAnnuityPremium(double t, double* result) {
  result[0] += bj_00(t);
}

void dy_TemporaryLifeAnnuityPremium(double t, double* V,double* result)
{
    result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
