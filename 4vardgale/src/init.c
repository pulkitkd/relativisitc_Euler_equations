#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "dg.h"
#include "dg1d.h"

void SetTestCaseData();

// For static mesh, domain is (xmin0, xmax0)
// For moving mesh, domain may be extended to (xmin,xmax)
// (xmin0,xmax0) is used to find dx, actual computational
// domain is (xmin,xmax)
// If velocity near boundaries is zero, as in many test cases,
// then these two domains are same.
void SetTestCaseData() {
  if (test_case == SOD) {
    printf("Setting SOD test case\n");
    xmin0 = 0.0;
    xmin = xmin0;
    xmax0 = 1.0;
    xmax = xmax0;
    XS = 0.5;
    finaltime = 0.2;

    d_left = 1.0;
    u1_left = 0.0;
    u2_left = 0.0;
    p_left = 1.0;
    d_right = 0.125;
    u1_right = 0.0;
    u2_right = 0.0;
    p_right = 0.1;
    bc_left = FREE;
    bc_right = FREE;
  }
  else if (test_case == MSOD) {  // Entropy problem from Toro

    printf("Setting MSOD test case\n");
    xmin0 = 0.0;
    xmin = xmin0;
    xmax0 = 1.0;
    xmax = xmax0;
    XS = 0.3;
    finaltime = 0.5;

    d_left = 1.0;
    u1_left = 0.75;
    u2_left = 0.0;
    p_left = 1.0;
    d_right = 0.125;
    u1_right = 0.0;
    u2_right = 0.0;
    p_right = 0.1;
    bc_left = FREE;
    bc_right = FREE;
  }
  else if (test_case == SSOD) {   // See Springel paper

    printf("Setting SSOD test case\n");
    xmin0 = -10.0;
    xmin = xmin0;
    xmax0 = 10.0;
    xmax = xmax0;
    XS = 0.0;
    finaltime = 1.0;

    d_left = 1.0;
    u1_left = 0.0;
    u2_left = 0.0;
    p_left = 1.0;
    d_right = 0.25;
    u1_right = 0.0;
    u2_right = 0.0;
    p_right = 0.1795;
    bc_left = FREE;
    bc_right = FREE;
  }
  else if (test_case == BLAST) { //only BLAST is configured for relativistic case as of now
    printf("Setting BLAST test case\n");
    xmin0 = 0.0;
    xmin = xmin0;
    xmax0 = 1.0;
    xmax = xmax0;
    finaltime = 0.2;
    bc_left = PERIODIC;
    bc_right = PERIODIC;
  }
  else if (test_case == LAX) {
    printf("Setting LAX test case\n");
    xmin0 = -0.0;
    xmin = xmin0;
    xmax0 = +10.0;
    xmax = xmax0;
    XS = 5.0;
    finaltime = 1.3;
    d_left = 0.445;
    d_right = 0.5;

    u1_left = 0.698;
    u2_left = 1.0;
    u1_right = 0.0;
    u2_right = 0.0;

    p_left = 3.528;
    p_right = 0.571;
    bc_left = FREE;
    bc_right = FREE;
  }
  else if (test_case == LOWD) {
    printf("Setting LOWD test case\n");
    xmin0 = 0.0;
    xmin = xmin0;
    xmax0 = 1.0;
    xmax = xmax0;
    XS = 0.5;
    finaltime = 0.15;
    d_left = 1.0;
    d_right = 1.0;

    u1_left = -2.0;
    u2_left = 1.0;
    u1_right = 2.0;
    u2_right = 1.0;

    p_left = 0.4;
    p_right = 0.4;
    bc_left = FREE;
    bc_right = FREE;
  }
  else if (test_case == SHUOSHER) {
    printf("Setting SHUOSHER test case\n");
    xmin0 = -5.0;
    xmin = xmin0;
    xmax0 = 5.0;
    xmax = xmax0;
    finaltime = 1.8;
    bc_left = FREE;
    bc_right = FREE;
    if (ALE) {
      xmin = -10.0;
    }
  }
  else if (test_case == PULSE) {
    printf("Setting PULSE test case \n");
    xmin0 = -10.0;
    xmin = xmin0;
    xmax0 = 10.0;
    xmax = xmax0;
    XS = 0.5;
    finaltime = 1;

    d_left = 1.0;
    u1_left = 0.0;
    u2_left = 0.0;
    p_left = 1.0;
    d_right = 0.125;
    u1_right = 0.0;
    u2_right = 0.0;
    p_right = 0.1;
    bc_left = FREE;
    bc_right = FREE;
  }
  else if (test_case == POLY) {
    xmin0 = -10.0;
    xmin = xmin0;
    xmax0 = 10.0;
    xmax = xmax0;
    XS = 0.5;
    finaltime = 0.009;

    d_left = 1.0;
    u1_left = 0.0;
    u2_left = 0.0;
    p_left = 1.0;
    d_right = 0.125;
    u1_right = 0.0;
    u2_right = 0.0;
    p_right = 0.1;
    bc_left = FREE;
    bc_right = FREE;
  }
  else if (test_case == CONTACT) {
    xmin0 = 0.0;
    xmax0 = 1.0;
    xmin = -1.0;
    xmax = xmax0;
    XS = 0.5;
    finaltime = 0.5;

    d_left = 2.0;
    u1_left = 1.0;
    u2_left = 0.0;
    p_left = 1.0;
    d_right = 1.0;
    u1_right = 1.0;
    u2_right = 0.0;
    p_right = 1.0;
    bc_left = FREE;
    bc_right = FREE;
  }
  else if (test_case == NOH) {
    xmin0 = 0.0;
    xmax0 = 1.0;
    xmin = -1.0;
    xmax = 2.0;
    XS = 0.5;
    finaltime = 1.0;

    d_left = 1.0;
    u1_left = 1.0;
    u2_left = 0.0;
    p_left = 1.0e-5;
    d_right = 1.0;
    u1_right = -1.0;
    u2_right = 0.0;
    p_right = 1.0e-5;
    bc_left = FREE;
    bc_right = FREE;
  }
  else if (test_case == RELBLAST) {
    printf("Setting Relativistic Blast test case\n");
    xmin0 = 0.0;
    xmin = xmin0;
    xmax0 = 1.0;
    xmax = xmax0;
    finaltime = 0.5;
    bc_left = FIXED;
    bc_right = FIXED;

  }
  else {
    printf("Error in SetTestCaseData\n");
    exit(0);
  }

  assert(xmin0 >= xmin);
  assert(xmax0 <= xmax);
}

/* Allocate memory for main structure and set initial conditions */
CELL *Init() {
  void ReadInput();
  void InitCondEuler(REAL, REAL *);
  UINT i, j, k, l, ncl, ncr;
  REAL dx, U[NVAR], v;
  CELL *cell;

  ReadInput();
  SetTestCaseData();

  // Note that PORD = degree + 1
  // Number of Gauss quadrature points = degree + 1
  NG = PORD;

  // Gauss-Lobatto points used in positivity limiter
  NGLL = PORD;
  if (NGLL > 5) {
    printf("Update GLL points in gauss.c\n");
    exit(0);
  }

  printf("Allocating memory and setting initial condition ...\n");

  dx = (xmax0 - xmin0) / NC;

  // For some test cases, we extend domain size,
  // so add more cells
  ncl = (xmin0 - xmin) / dx;
  ncr = (xmax - xmax0) / dx;
  NC += ncl + ncr;
  NF = NC + 1;

  printf("No of cells = %d\n", NC);
  printf("No of faces = %d\n", NF);
  printf("         dx = %f\n", dx);

  n_extra_cells = 4 * NC;

  NCMAX = NC + n_extra_cells;
  NFMAX = NCMAX + 1;

  free_cell_index = NC;
  free_face_index = NF;

  printf("Maximum no of cells = %d\n", NCMAX);
  printf("Maximum no of faces = %d\n", NFMAX);

  cell = (CELL *)calloc(NCMAX, sizeof(CELL));
  if (cell == NULL) {
    printf("Init: Could not allocate cell\n");
    exit(0);
  }

  // Initialize cells
  for (i = 0; i < NCMAX; i++) {
    if (i < NC) {
      cell[i].active = true;
    } else {
      cell[i].active = false;
    }

    cell[i].level = 0;
    cell[i].xl = xmin + i * dx;
    cell[i].xr = cell[i].xl + dx;
    cell[i].x = 0.5 * (cell[i].xl + cell[i].xr);
    cell[i].h = cell[i].xr - cell[i].xl;

    cell[i].wl = 0.0;
    cell[i].wr = 0.0;

    cell[i].p = PORD;

    cell[i].ng = NG;
    cell[i].ngll = NGLL;
    cell[i].xg = (REAL *)calloc(cell[i].ng, sizeof(REAL));
    GaussPoints(&cell[i]);

    cell[i].U = (REAL **)calloc(NVAR, sizeof(REAL *));
    cell[i].Re = (REAL **)calloc(NVAR, sizeof(REAL *));
    for (j = 0; j < NVAR; j++) {
      cell[i].U[j] = (REAL *)calloc(cell[i].p, sizeof(REAL));
      cell[i].Re[j] = (REAL *)calloc(cell[i].p, sizeof(REAL));
    }

    // Use this for Periodic BC
/*    if (i == 0) {
      cell[i].lcell = &cell[NC - 1];
      cell[i].rcell = &cell[i + 1];
    } else if (i == NC - 1) {
      cell[i].lcell = &cell[i - 1];
      cell[i].rcell = &cell[0];
    } else {
      cell[i].lcell = &cell[i - 1];
      cell[i].rcell = &cell[i + 1];
    }*/

    // Use this for Fixed BC
    if (i == 0) {
      cell[i].lcell = &cell[i];
      cell[i].rcell = &cell[i + 1];
    } else if (i == NC - 1) {
      cell[i].lcell = &cell[i - 1];
      cell[i].rcell = &cell[i];
    } else {
      cell[i].lcell = &cell[i - 1];
      cell[i].rcell = &cell[i + 1];
    }
  }

  /* Set initial condition by L2 projection */
  for (i = 0; i < NC; i++) {

    for (j = 0; j < NVAR; j++)
      for (k = 0; k < cell[i].p; k++)
        cell[i].U[j][k] = 0.0;

    for (j = 0; j < cell[i].p; j++)
      for (k = 0; k < cell[i].ng; k++) {
        InitCondEuler(cell[i].xg[k], U);
        v = ShapeFun(cell[i].xg[k], &cell[i], j);
        for (l = 0; l < NVAR; l++){
          cell[i].U[l][j] += 0.5 * U[l] * v * wg[cell[i].ng - 1][k];
          //printf(" Debug: conserved variables in init() %lf \n", cell[i].U[l][j]);
          }
      }
//      printf(" Debug: conserved variables in init() %lf %lf %lf %lf \n",
//       cell[i].U[0][0], cell[i].U[1][0], cell[i].U[2][0], cell[i].U[3][0]);
  }

  return cell;
}

FACE *InitFaces(CELL *cell) {
  UINT i;
  REAL dx;
  FACE *face;

  face = (FACE *)calloc(NFMAX, sizeof(FACE));
  if (face == NULL) {
    printf("Init: Could not allocate face\n");
    exit(0);
  }

  dx = (xmax - xmin) / NC;
  dx_init = dx;

  // Initialize faces
  face[0].x = xmin;
  face[0].w = 0.0;
//  face[0].lcell = &cell[NC - 1]; // for Periodic BC without adaptation
  face[0].lcell = &cell[0];      // for Fixed BC
  face[0].rcell = &cell[0];
  cell[0].lface = &face[0];
  face[0].active = true;
  for (i = 1; i < NFMAX; i++)
    if (i != NF - 1) {
      {
        if (i < NF - 1)
          face[i].active = true;
        else
          face[i].active = false;
        face[i].x = xmin + i * dx;
        face[i].w = 0.0;
        face[i].lcell = &cell[i - 1];
        face[i].rcell = &cell[i];
        cell[i - 1].rface = &face[i];
        cell[i].lface = &face[i];
      }
    }
  face[NF - 1].x = xmax;
  face[NF - 1].w = 0.0;
  face[NF - 1].lcell = &cell[NC - 1]; // change this for periodic bc
//  face[NF - 1].rcell = &cell[0];  // Use this for Periodic BC without adaptation
  face[NF - 1].rcell = &cell[NC - 1];
  cell[NC - 1].rface = &face[NF - 1];
  face[NF - 1].active = true;

  return face;
}

void ReadInput() {
  FILE *fp;
  char dummy[100];
  fp = fopen("inp.dat", "r");
  if (fp == NULL) {
    printf("Error: Could not open inp.dat\n");
    exit(0);
  }

  fscanf(fp, "%s%lf", dummy, &cfl);
  fscanf(fp, "%s%d", dummy, &NC);
  fscanf(fp, "%s%d", dummy, &PORD);
  fscanf(fp, "%s%d", dummy, &NPLT);
  fscanf(fp, "%s%d", dummy, &FLUX);
  fscanf(fp, "%s%lf", dummy, &Mfact);
  fscanf(fp, "%s%d", dummy, &ALE);
  fscanf(fp, "%s%d", dummy, &test_case);
  fscanf(fp, "%s%d", dummy, &pos_lim);
  fscanf(fp, "%s%d", dummy, &tvd_lim);
  fscanf(fp, "%s%d", dummy, &adaptation);
  fscanf(fp, "%s%lf", dummy, &hmin);
  fscanf(fp, "%s%lf", dummy, &hmax);
  fscanf(fp, "%s%d", dummy, &siter);
  fscanf(fp, "%s%lf", dummy, &chi);
  fscanf(fp, "%s%lf", dummy, &eta);

  fclose(fp);

  switch (PORD) {
  case 1:
  case 2:
    predictor_method = taylor_1p;
    break;
  case 3:
    predictor_method = cerk_nodal_2p;
    break;
  case 4:
    predictor_method = cerk_nodal_3p;
    break;
  default:
    printf("Invalid PORD value. \n");
    exit(1);
  }
}
// Initial condition for Euler equation
// Assigns values to Primitives which are
// converted to Conserved quantities later on
void InitCondEuler(REAL x, REAL *U) {
  REAL V[NVAR];

  if (test_case == SOD || test_case == LAX || test_case == MSOD ||
      test_case == SSOD || test_case == LOWD || test_case == CONTACT ||
      test_case == NOH)
    InitCondShocktube(x, V);
  else if (test_case == BLAST)
    InitCondBlast(x, V);
  else if (test_case == SHUOSHER)
    InitCondShuOsher(x, V);
  else if (test_case == PULSE) {
    InitCondPulse(x, V);
  } else if (test_case == POLY) {
    InitCondPoly(x, V);
  } //else if (test_case == RELBLAST)
    //InitCondRelBlast(x, V);
    else {
    printf("Error: Unknown test case\n");
    exit(0);
  }
  // Conserved quantities
  REAL rho = V[0];
  REAL u1 = V[1];
  REAL u2 = V[2];
  REAL p = V[3];
  REAL alpha = GAMMA / (GAMMA -1);
  REAL h = 1 + p * alpha / rho;
  REAL BETAsq = 1 / (1 - u1 * u1 - u2 * u2);

  U[0] = rho * sqrt(BETAsq); //D
  U[1] = BETAsq * rho * h * u1; //m1
  U[2] = BETAsq * rho * h * u2; //m2
  U[3] = BETAsq * rho * h - p; //E
  //for(int q=0; q<NVAR; q++)
//printf("Debug: conserved variables in initcondeuler() %lf \n", U[q]);

}

void InitCondShocktube(REAL x, REAL *V) {
  if (x < XS) {
    V[0] = d_left;
    V[1] = u1_left;
    V[2] = u2_left;
    V[3] = p_left;
    }
    else {
    V[0] = d_right;
    V[1] = u1_right;
    V[2] = u2_right;
    V[3] = p_right;
    }
}

void InitCondBlast(REAL x, REAL *V) {
  if(x < 0.5){
    V[0] = 1.0;
    V[1] = 0.0;
    V[2] = 0.0;
    V[3] = 1000.0;
    }
//  else if (x > 0.66){
//    V[0] = 1.0;
//    V[1] = 0.0;
//    V[2] = 0.0;
//    V[3] = 1.0;
//    }
  else {
    V[0] = 1.0;
    V[1] = 0.0;
    V[2] = 0.0;
    V[3] = 0.01;
  }
//  if (x < 0.1) {
//    V[0] = 1.0;
//    V[1] = 0.0;
//    V[2] = 0.0;
//    V[3] = 1000.0;
//  } else if (x > 0.9) {
//    V[0] = 1.0;
//    V[1] = 0.0;
//    V[2] = 0.0;
//    V[3] = 100.0;
//  } else {
//    V[0] = 1.0;
//    V[1] = 0.0;
//    V[2] = 0.0;
//    V[3] = 0.01;
//  }
/*  if(x<0.5){
  V[0] = 1 + exp(-1000 * (x-0.25) * (x-0.25));
  V[1] = 0.9;
  V[2] = 0;
  V[3] = 1;
  }
  else{
  V[0] = 1 + exp(-1000 * (x-0.75) * (x-0.75));
  V[1] = -0.9;
  V[2] = 0;
  V[3] = 1;
  }*/
/*  V[0] = sin(x);
  V[1] = 0.5;
  V[2] = 0;
  V[3] = 1;*/
}

void InitCondShuOsher(REAL x, REAL *V) {
  if (x < -4.0) {
    V[0] = 3.857143;
    V[1] = 2.629369;
    V[2] = 0.0;
    V[3] = 10.333333;
  } else {
    V[0] = 1 + 0.5 * exp(-5.0 * x * x);
    V[1] = 0;
    V[2] = 0;
    V[3] = 1;
  }
}

void InitCondPulse(REAL x, REAL *V) {
  V[0] = 1 + exp(-10 * x * x);
  V[1] = 1;
  V[2] = 0;
  V[3] = 1;
}

void InitCondPoly(REAL x, REAL *V) {
  V[1] = 1;
  V[2] = 0;
  V[3] = 1;
  int order = 3;
  if (x < 0.0)
    V[0] = 1;
  else if (x >= 0.0 && x <= 1.0)
    V[0] = 1 + pow(x, order);
  else if (x >= 1.0 && x <= 2.0)
    V[0] = 1 + 2 - pow(2 - x, order);
  else if (x >= 2.0 && x <= 3.0)
    V[0] = 1 + 2 - pow(x - 2, order);
  else if (x >= 3.0 && x <= 4.0)
    V[0] = 1 + pow(4 - x, order);
  else
    V[0] = 1;
}

//Initial conditions for Moderate Blast Wave - A relativistic test case
//------------------------------------------------------------------------------
//void InitCondRelBlast(REAL x, REAL *V){
//  if (x < 0.5){
//    V[0] = 10.0;
//    V[1] = 0.0;
//    V[2] = 0.0;
//    V[3] = 13.33;
//    }
//  else{
//    V[0] = 1.0;
//    V[1] = 0.0;
//    V[2] = 0.0;
//    V[3] = 0.0;
//    }
//}
//------------------------------------------------------------------------------

