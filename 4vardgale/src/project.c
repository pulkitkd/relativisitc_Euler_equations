#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dg.h"
#include "dg1d.h"

REAL minmod(REAL, REAL, REAL);
void ApplyPositivityLimiter(CELL *cell);
void ApplyTVDLimiter(CELL *cell);
void ApplyMomentLimiter(CELL *cell);

void ApplyLimiter(CELL *cell) {

  if (tvd_lim == 0)
    return;
  else if (tvd_lim == 1)
    ApplyTVDLimiter(cell);
  else if (tvd_lim == 2)
    ApplyMomentLimiter(cell);
  else {
    printf("Unknown value of tvd limiter = %d\n", tvd_lim);
    exit(1);
  }
}

// TVD/TVB limiter
// If Mfact=0, then it becomes TVD limiter.
// For Mfact>0, we have TVB limiter.
void ApplyTVDLimiter(CELL *cell) {
  UINT i, j, k;
  REAL u[NVAR], ux[NVAR], uxb[NVAR], dul[NVAR], dur[NVAR], R[NVAR][NVAR],
    Ri[NVAR][NVAR], fact;
  REAL dxl, dxr;

  fact = sqrt(3.0);

  for (i = 0; i < NCMAX; i++)
    if (cell[i].active) {
      dxl = 0.5 * (cell[i].h + cell[i].lcell->h);
      dxr = 0.5 * (cell[i].h + cell[i].rcell->h);
      for (j = 0; j < NVAR; j++) {
        dul[j] = (cell[i].U[j][0] - cell[i].lcell->U[j][0]) / dxl;
        dur[j] = (cell[i].rcell->U[j][0] - cell[i].U[j][0]) / dxr;
        u[j] = cell[i].U[j][0];
        ux[j] = fact * cell[i].U[j][1] / cell[i].h;
      }

      EigMat(u, R, Ri);
      Multi(Ri, ux);
      Multi(Ri, dul);
      Multi(Ri, dur);
      for (j = 0; j < NVAR; j++) {
        if (fabs(ux[j]) <= Mfact * cell[i].h)
          uxb[j] = ux[j];
        else
          uxb[j] = minmod(ux[j], dul[j], dur[j]);
      }
      Multi(R, uxb);

      for (j = 0; j < NVAR; j++) {
        uxb[j] = cell[i].h * uxb[j] / fact;
        if (fabs(cell[i].U[j][1] - uxb[j]) > 1.0e-6) {
          cell[i].U[j][1] = uxb[j];
          for (k = 2; k < cell[i].p; k++)
            cell[i].U[j][k] = 0.0;
        }
      }
      if (pos_lim == TRUE)
        ApplyPositivityLimiter(&cell[i]);
    }
}

// Moment limiter of Biswas, Devine, Flaherty
void ApplyMomentLimiter(CELL *cell) {
  UINT i, j, k;
  REAL u[NVAR], ux[NVAR], uxb[NVAR], dul[NVAR], dur[NVAR], R[NVAR][NVAR],
    Ri[NVAR][NVAR], f1, f2;
  REAL dxl, dxr;
  REAL diff;
  bool tolimit[NCMAX];

  for (i = 0; i < NCMAX; i++)
    tolimit[i] = true;

  // start limiting from the last dof
  // compare U[][k] with finite differences of U[][k-1]
  for (k = PORD - 1; k > 0; k--) {
    f1 = (2.0 * k - 1.0) * sqrt(2.0 * k + 1.0);
    f2 = sqrt(2.0 * k - 1.0);
    for (i = 0; i < NCMAX; i++)
      if (cell[i].active && tolimit[i]) {
        dxl = 0.5 * (cell[i].h + cell[i].lcell->h);
        dxr = 0.5 * (cell[i].h + cell[i].rcell->h);
        // copy cell average value and compute eigenvector matrices
        for (j = 0; j < NVAR; j++)
          u[j] = cell[i].U[j][0];
        EigMat(u, R, Ri);

        for (j = 0; j < NVAR; j++) {
          dul[j] = (cell[i].U[j][k - 1] - cell[i].lcell->U[j][k - 1]) / dxl;
          dur[j] = (cell[i].rcell->U[j][k - 1] - cell[i].U[j][k - 1]) / dxr;
          ux[j] = cell[i].U[j][k] / cell[i].h;
        }

        Multi(Ri, ux);
        Multi(Ri, dul);
        Multi(Ri, dur);

        diff = 0.0;
        for (j = 0; j < NVAR; j++) {
          uxb[j] = minmod(f1 * ux[j], f2 * dul[j], f2 * dur[j]) / f1;
          diff += (fabs(ux[j] - uxb[j]) > 1.0e-2*fabs(ux[j]));
        }

        if (diff > 0) {
          Multi(R, uxb);

          for (j = 0; j < NVAR; j++)
            cell[i].U[j][k] = cell[i].h * uxb[j];
          tolimit[i] = true;
        } else {
          tolimit[i] = false;
        }
      }
   }

   if (pos_lim == TRUE)
     for (i = 0; i < NCMAX; i++)
       if (cell[i].active)
         ApplyPositivityLimiter(&cell[i]);
}

/* minmod limiter function */
REAL minmod(REAL a, REAL b, REAL c) {
  REAL sgn, m;

  if (a * b <= 0.0 || b * c <= 0.0)
    return 0.0;

  sgn = (a > 0.0) ? 1.0 : -1.0;
  a = fabs(a);
  b = fabs(b);
  c = fabs(c);
  m = (a < b) ? a : b;
  m = (c < m) ? c : m;
  return sgn * m;
}

/* Eigenvector matrix */
void EigMat(REAL *U, REAL R[][NVAR], REAL Ri[][NVAR]) {

    double V[NVAR], UV[NVAR][NVAR], R1[NVAR][NVAR], R2[NVAR][NVAR];
    double UVi[NVAR][NVAR], temp[NVAR][NVAR];
    int i, j;
    con2prim(U, V);

    double r, u1, u2, p;  // rest mass density, pressure and velocities
    r = V[0]; //defined as rho usually
    u1 = V[1];
    u2 = V[2];
    p = V[3];
    REAL a = GAMMA / (GAMMA - 1); // defined as alpha usually
    REAL h = 1 + p * a / r;
    REAL u1sq = u1 * u1;
    REAL u2sq = u2 * u2;
    REAL usq = u1sq + u2sq;
    REAL BETA = 1 / sqrt(1 - usq);
    REAL BETAsq = BETA * BETA;
    REAL H = a * a * p - r - a * (-r + p * (1 + usq));

    /* The dU/dV matrix */
    UV[0][0] = BETA;
    UV[1][0] = u1 * BETAsq;
    UV[2][0] = u2 * BETAsq;
    UV[3][0] = BETAsq;

    UV[0][1] = r * u1 * BETAsq * BETA;
    UV[1][1] = r * h * (1 + u1 * u1 - u2 * u2) * BETAsq * BETAsq;
    UV[2][1] = 2 * r * h * u1 * u2 * BETAsq * BETAsq;
    UV[3][1] = 2 * r * h * u1 * BETAsq * BETAsq;

    UV[0][2] = r * u2 * BETAsq * BETA;
    UV[1][2] = 2 * r * h * u1 * u2 * BETAsq * BETAsq;;
    UV[2][2] = r * h * (1 - u1 * u1 + u2 * u2) * BETAsq * BETAsq;
    UV[3][2] = 2 * r * h * u2 * BETAsq * BETAsq;

    UV[0][3] = 0;
    UV[1][3] = a * u1 * BETAsq;
    UV[2][3] = a * u2 * BETAsq;
    UV[3][3] = -1 + a * BETAsq;

    /*The inverse of dU/dV matrix*/
    UVi[0][0] = (h * (1 - a + usq)) / (BETA * H);
    UVi[0][1] = -u1 / (H * BETAsq * BETA);
    UVi[0][2] = -u2 / (H * BETAsq * BETA);
    UVi[0][3] = h / (BETA * H);

    UVi[1][0] = (r * u1 * (-1 + a + usq)) / H;
    UVi[1][1] = -(a * a * p + r * (-1 + u1 * u1) +
                a * (r + p * (-1 + u1 * u1 - u2 * u2))) /
                (BETAsq * H * h);
    UVi[1][2] = -(2 * a * p + r) * u1 * u2 / (BETAsq * H * h);
    UVi[1][3] = ((2 * a * p + r) * u1) / H;

    UVi[2][1] = (r * u2 * (-1 + a + u1sq)) / H;
    UVi[2][2] = -((2 * a * p + r) * u1 * u2) / (BETAsq * H * h);
    UVi[2][3] = -(a * a * p + r - r * u2 * u2 +
                a * (-r + p * (1 + u1 * u1 - u2 * u2))) / (BETAsq * H * h);
    UVi[2][4] = ((2 * a * p + r) * u2) / H;

    UVi[3][1] = -((a * r * usq)) / H;
    UVi[3][2] = a * u1 / (BETAsq * H);
    UVi[3][3] = a * u2 / (BETAsq * H);
    UVi[3][4] = (-h * (1 + usq)) / H;


  /* R1 stores primitive right eigenvectors */
    R1[0][0] = 0;
    R1[1][0] = 1;
    R1[2][0] = r / ( p * GAMMA);
    R1[3][0] = r / ( p * GAMMA);

    R1[0][1] = 0;
    R1[1][1] = 0;
    R1[2][1] = sqrt(-a * p * ((-1 + a) * h * (-1 + u1sq) + a * p * u2sq)) / (a * p * h * BETA);
    R1[3][1] = -sqrt(-a * p * ((-1 + a) * h * (-1 + u1sq) + a * p * u2sq)) / (a * p * h * BETA);

    R1[0][2] = 1;
    R1[1][2] = 0;
    R1[2][2] = -a * p * h * (-1 + u1sq) * u2 - a * p * h * u2sq * u2 + u1 * u2 *
		sqrt(-a * p * ((-1 + a) * h * (-1 + u1sq) + a * p * u2sq)) / (a * p * h * BETA *
		(-1 + u1sq));

    R1[3][2] =  -a * p * h * (-1 + u1sq) * u2 - a * p * h * u2sq * u2 - u1 * u2 *
		sqrt(-a * p * ((-1 + a) * h * (-1 + u1sq) + a * p * u2sq)) / (a * p * h * BETA *
		(-1 + u1sq));

    R1[0][3] = 0;
    R1[1][3] = 0;
    R1[2][3] = 1;
    R1[3][3] = 1;

    /* R2 stores inverse of R1 */
    /* Mathematically, it should represent the Matrix of Left Eigenvectors*/
    R2[0][0] = -GAMMA * p / r;
    R2[0][1] = -a * p * u2 / ((a - 1) * BETAsq * r * h * (-1 + u1sq));
    R2[0][2] = 0;
    R2[0][3] = 0;

    R2[1][0] = 0;
    R2[1][1] = u1 * u2 / (1 - u1sq);
    R2[1][2] = (a * p * h * BETA) / sqrt(-a * p * (r - r * u1 * u1 + a * a * p * (-1 + u1sq) +
                 a * (r * (-1 + u1sq) + p * (1 - u1sq + u2sq))));
    R2[1][3] = -(a * p * h * BETA) / sqrt(-a * p * (r - r * u1 * u1 + a * a * p * (-1 + u1sq) +
                 a * (r * (-1 + u1sq) + p * (1 - u1sq + u2sq))));

    R2[2][0] = 0;
    R2[2][1] = 1;
    R2[2][2] = 0;
    R2[2][3] = 0;

    R2[3][0] = 1;
    R2[3][1] = 0;
    R2[3][2] = 1;
    R2[3][3] = 1;

    /*We now multiply them with du/dv to get the conserved eigenvectors*/

   MultiplyMatrix(UV, R1, R);
   MultiplyMatrix(R2, UVi, Ri);



//    MultiplyMatrix(R1, R2, temp);
//    for(i = 0; i < NVAR; i++){
//      for(j = 0; j < NVAR; j++){
//        printf("%lf \t", temp[i][j]);
//        }
//        printf("\n");
//     }
//     printf("\n");
//  g1 = GAMMA - 1.0;
//  g2 = g1 / 2.0;
//
//  d = U[0];
//  //u1 = U[1] / d;
//  v = U[1] / d;
//  u2 = U[2] / d;
//  p = (GAMMA - 1.0) * (U[3] - 0.5 * d * u1 * u1);
//  c = sqrt(GAMMA * p / d);
//  h = c * c / g1 + 0.5 * u1 * u1;
//  f = d / c / 2.0;
//
//  /* Inverse eigenvector-matrix */
//  Ri[0][0] = 1.0 - g2 * v * v / c / c;
//  Ri[1][0] = (g2 * v * v - v * c) / d / c;
//  Ri[2][0] = -(g2 * v * v + v * c) / d / c;
//
//  Ri[0][1] = g1 * v / c / c;
//  Ri[1][1] = (c - g1 * v) / d / c;
//  Ri[2][1] = (c + g1 * v) / d / c;
//
//  Ri[0][2] = -g1 / c / c;
//  Ri[1][2] = g1 / d / c;
//  Ri[2][2] = -g1 / d / c;
//
//  /* Eigenvector matrix */
//  R[0][0] = 1.0;
//  R[1][0] = v;
//  R[2][0] = v * v / 2.0;
//
//  R[0][1] = f;
//  R[1][1] = (v + c) * f;
//  R[2][1] = (h + v * c) * f;
//
//  R[0][2] = -f;
//  R[1][2] = -(v - c) * f;
//  R[2][2] = -(h - v * c) * f;
}

/* Multiply matrix R and vector U */
void Multi(REAL R[][NVAR], REAL *U) {
  UINT i, j;
  REAL Ut[NVAR];

  for (i = 0; i < NVAR; i++)
    Ut[i] = U[i];

  for (i = 0; i < NVAR; i++) {
    U[i] = 0.0;
    for (j = 0; j < NVAR; j++)
      U[i] += R[i][j] * Ut[j];
  }
}

/* Multiply matrix R1 and R2 and store the result in A */
void MultiplyMatrix(REAL R1[][NVAR], REAL R2[][NVAR], REAL A[NVAR][NVAR]) {
  UINT i, j, k;

  for(k = 0; k < NVAR; k++)
    for (i = 0; i < NVAR; i++)
    {
      A[i][k] = 0.0;
      for (j = 0; j < NVAR; j++)
        A[i][k] += R1[i][j] * R2[j][k];
    }
}

// Apply positivity limiter in cell
void ApplyPositivityLimiter(CELL *cell) {
  UINT i, j;
  const REAL eps = 1.0e-12;
  REAL theta1, theta2, rho_min, rat, pre, t1, t2, t, a1, b1, c1, D;
  REAL drho, dm, dE;
  REAL **U;

  // First order scheme, nothing to do
  if (cell->p == 1)
    return;

  U = (REAL **)calloc(cell->ngll, sizeof(REAL *));
  for (i = 0; i < cell->ngll; i++)
    U[i] = (REAL *)calloc(NVAR, sizeof(REAL));

  // First, limit density

  // Compute solution at GLL nodes
  UatGLL(cell, U);

  // Find minimum value of density
  rho_min = 1.0e20;
  for (i = 0; i < cell->ngll; ++i)
    rho_min = MIN(rho_min, U[i][0]);

  rat = fabs(cell->U[0][0] - eps) / (fabs(cell->U[0][0] - rho_min) + 1.0e-14);
  theta1 = MIN(1.0, rat);

  // Apply limiter, dont change mean value
  for (i = 0; i < NVAR; ++i)
    for (j = 1; j < cell->p; ++j)
      cell->U[i][j] *= theta1;

  // Now limit the pressure

  // Compute solution at GLL nodes
  UatGLL(cell, U);

  theta2 = 1.0;
  for (i = 0; i < cell->ngll; ++i) {
    pre = (GAMMA - 1.0) * (U[i][2] - 0.5 * U[i][1] * U[i][1] / U[i][0]);
    if (pre < eps) {
      drho = U[i][0] - cell->U[0][0];
      dm = U[i][1] - cell->U[1][0];
      dE = U[i][2] - cell->U[2][0];
      a1 = 2.0 * drho * dE - dm * dm;
      b1 = 2.0 * drho * (cell->U[2][0] - eps / (GAMMA - 1.0)) +
           2.0 * cell->U[0][0] * dE - 2.0 * cell->U[1][0] * dm;
      c1 = 2.0 * cell->U[0][0] * cell->U[2][0] - pow(cell->U[1][0], 2) -
           2.0 * eps * cell->U[0][0] / (GAMMA - 1.0);
      // Divide by a1 to avoid round-off error
      b1 /= a1;
      c1 /= a1;
      D = sqrt(fabs(b1 * b1 - 4.0 * c1));
      t1 = 0.5 * (-b1 - D);
      t2 = 0.5 * (-b1 + D);
      if (t1 > -1.0e-12 && t1 < 1.0 + 1.0e-12)
        t = t1;
      else if (t2 > -1.0e-12 && t2 < 1.0 + 1.0e-12)
        t = t2;
      else {
        printf("Fatal error\n");
        printf("Mean rho = %e", cell->U[0][0]);
        printf("pre at gll point = %e\n", pre);
        printf("t1 = %e, t2 = %e\n", t1, t2);
        exit(0);
      }
      t = MIN(1.0, t);
      t = MAX(0.0, t);
      // Need t < 1.0. If t==1 upto machine precision
      // then we are suffering from round off error.
      // In this case we take the cell average value, t=0.
      if (fabs(1.0 - t) < 1.0e-14)
        t = 0.0;
      theta2 = MIN(theta2, t);
    }
  }

  // Apply limiter, dont change mean value
  for (i = 0; i < NVAR; ++i)
    for (j = 1; j < cell->p; ++j)
      cell->U[i][j] *= theta2;

  // Release memory
  for (i = 0; i < cell->ngll; i++)
    free(U[i]);
  free(U);
}
