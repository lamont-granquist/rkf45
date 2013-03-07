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
#include "RK_Library.h"
#include "DeferredTemporaryLifeAnnuity.h"

static double *matrix[51];

static double row0[]  = {1.0265607676014400};
static double row1[]  = {1.0805663523022000};
static double row2[]  = {1.1374932838714300};
static double row3[]  = {1.1975114275626800};
static double row4[]  = {1.2608022416886800};
static double row5[]  = {1.3275598043516300};
static double row6[]  = {1.3979919596875600};
static double row7[]  = {1.4723216004702400};
static double row8[]  = {1.5507881065874600};
static double row9[]  = {1.6336489620308100};
static double row10[] = {1.7211815767162100};
static double row11[] = {1.8136853437812100};
static double row12[] = {1.9114839681141600};
static double row13[] = {2.0149281079127500};
static double row14[] = {2.1243983782341300};
static double row15[] = {2.2403087740143500};
static double row16[] = {2.3631105801829600};
static double row17[] = {2.4932968486250200};
static double row18[] = {2.6314075362739100};
static double row19[] = {2.7780354160836500};
static double row20[] = {2.9338328936850100};
static double row21[] = {3.0995198879959800};
static double row22[] = {3.2758929649601700};
static double row23[] = {3.4638359512170500};
static double row24[] = {3.6643323004950200};
static double row25[] = {3.8784795419259000};
static double row26[] = {4.1075062089359300};
static double row27[] = {4.3527917332333000};
static double row28[] = {4.6158898949982800};
static double row29[] = {4.8985565532549200};
static double row30[] = {5.2027825467745500};
static double row31[] = {5.5308328651267400};
static double row32[] = {5.8852934539509200};
static double row33[] = {6.2691273543593000};
static double row34[] = {6.6857423050148000};
static double row35[] = {7.1390724853580000};
static double row36[] = {6.5992994744162900};
static double row37[] = {6.0319762072565800};
static double row38[] = {5.4341959778729600};
static double row39[] = {4.8025044157278400};
static double row40[] = {4.1327786917263300};
static double row41[] = {3.4200771644266400};
static double row42[] = {2.6584512106845400};
static double row43[] = {1.8407083534208200};
static double row44[] = {0.9581122236639260};
static double row45[] = {0.0000000000000000};
static double row46[] = {0.0000000000000000};
static double row47[] = {0.0000000000000000};
static double row48[] = {0.0000000000000000};
static double row49[] = {0.0000000000000000};
static double row50[] = {0.0000000000000000};

/* tv */
double** tv_DeferredTemporaryLifeAnnuity() {

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
static double bdeath = 1;

static double b_0(double t) {
    return bpension * indicator(t > m) * indicator(t < m + n);
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

void bj_ii_DeferredTemporaryLifeAnnuity(double t, double* result) {
  result[0] += bj_00(t);
}

void dy_DeferredTemporaryLifeAnnuity(double t, double* V,double* result)
{
    result[0] = r(t) * V[0] - b_0(t) - mu_01(t) * (0 - V[0] + bj_01(t));
}
