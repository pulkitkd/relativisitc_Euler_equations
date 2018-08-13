#include <stdio.h>
#include <math.h>
#include "dg.h"
#include "dg1d.h"

/* Gauss integration points and weights */
void GaussInit() {
  UINT n = 10, i, j;

  printf("Calculating Gauss integration points and weights ...\n");

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      xg[i][j] = 0.0;
      wg[i][j] = 0.0;
    }

  /* Gauss integration points */
  xg[0][0] = 0.0;

  xg[1][0] = -1.0 / sqrt(3.0);
  xg[1][1] = 1.0 / sqrt(3.0);

  xg[2][0] = -sqrt(15.0) / 5.0;
  xg[2][1] = 0.0;
  xg[2][2] = sqrt(15.0) / 5.0;

  xg[3][0] = -sqrt(525.0 + 70.0 * sqrt(30.0)) / 35.0;
  xg[3][1] = -sqrt(525.0 - 70.0 * sqrt(30.0)) / 35.0;
  xg[3][2] = sqrt(525.0 - 70.0 * sqrt(30.0)) / 35.0;
  xg[3][3] = sqrt(525.0 + 70.0 * sqrt(30.0)) / 35.0;

  xg[4][0] = -sqrt(245.0 + 14.0 * sqrt(70.0)) / 21.0;
  xg[4][1] = -sqrt(245.0 - 14.0 * sqrt(70.0)) / 21.0;
  xg[4][2] = 0.0;
  xg[4][3] = sqrt(245.0 - 14.0 * sqrt(70.0)) / 21.0;
  xg[4][4] = sqrt(245.0 + 14.0 * sqrt(70.0)) / 21.0;

  xg[5][0] = -0.9324695142031521;
  xg[5][1] = -0.6612093864662645;
  xg[5][2] = -0.2386191860831969;
  xg[5][3] = 0.2386191860831969;
  xg[5][4] = 0.6612093864662645;
  xg[5][5] = 0.9324695142031521;

  xg[6][0] = -0.9491079123427585;
  xg[6][1] = -0.7415311855993945;
  xg[6][2] = -0.4058451513773972;
  xg[6][3] = 0.0000000000000000;
  xg[6][4] = 0.4058451513773972;
  xg[6][5] = 0.7415311855993945;
  xg[6][6] = 0.9491079123427585;

  xg[7][0] = -0.9602898564975363;
  xg[7][1] = -0.7966664774136267;
  xg[7][2] = -0.5255324099163290;
  xg[7][3] = -0.1834346424956498;
  xg[7][4] = 0.1834346424956498;
  xg[7][5] = 0.5255324099163290;
  xg[7][6] = 0.7966664774136267;
  xg[7][7] = 0.9602898564975363;

  /* Gaussian weights */
  wg[0][0] = 2.0;

  wg[1][0] = 1.0;
  wg[1][1] = 1.0;

  wg[2][0] = 5.0 / 9.0;
  wg[2][1] = 8.0 / 9.0;
  wg[2][2] = 5.0 / 9.0;

  wg[3][0] = (18.0 - sqrt(30.0)) / 36.0;
  wg[3][1] = (18.0 + sqrt(30.0)) / 36.0;
  wg[3][2] = (18.0 + sqrt(30.0)) / 36.0;
  wg[3][3] = (18.0 - sqrt(30.0)) / 36.0;

  wg[4][0] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
  wg[4][1] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
  wg[4][2] = 128.0 / 225.0;
  wg[4][3] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
  wg[4][4] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;

  wg[5][0] = 0.1713244923791704;
  wg[5][1] = 0.3607615730481386;
  wg[5][2] = 0.4679139345726910;
  wg[5][3] = 0.4679139345726910;
  wg[5][4] = 0.3607615730481386;
  wg[5][5] = 0.1713244923791704;

  wg[6][0] = 0.1294849661688697;
  wg[6][1] = 0.2797053914892766;
  wg[6][3] = 0.3818300505051189;
  wg[6][3] = 0.4179591836734694;
  wg[6][4] = 0.3818300505051189;
  wg[6][5] = 0.2797053914892766;
  wg[6][6] = 0.1294849661688697;

  wg[7][0] = 0.1012285362903763;
  wg[7][1] = 0.2223810344533745;
  wg[7][2] = 0.3137066458778873;
  wg[7][3] = 0.3626837833783620;
  wg[7][4] = 0.3626837833783620;
  wg[7][5] = 0.3137066458778873;
  wg[7][6] = 0.2223810344533745;
  wg[7][7] = 0.1012285362903763;

  /* GLL points */
  xgll[0][0] = 0.0;

  xgll[1][0] = -1.0;
  xgll[1][1] = +1.0;

  xgll[2][0] = -1.0;
  xgll[2][1] = 0.0;
  xgll[2][2] = +1.0;

  xgll[2][0] = -1.0;
  xgll[2][1] = 0.0;
  xgll[2][2] = +1.0;

  xgll[3][0] = -1.0;
  xgll[3][1] = -1.0 / sqrt(5.0);
  xgll[3][2] = 1.0 / sqrt(5.0);
  xgll[3][3] = 1.0;

  xgll[4][0] = -1.0;
  xgll[4][1] = -sqrt(3.0 / 7.0);
  xgll[4][2] = 0.00;
  xgll[4][3] = sqrt(3.0 / 7.0);
  xgll[4][4] = 1.0;

  /* GLL weights */
  wgll[0][0] = 2.0;

  wgll[1][0] = 1.0;
  wgll[1][1] = 1.0;

  wgll[2][0] = 1.0 / 3.0;
  wgll[2][1] = 4.0 / 3.0;
  wgll[2][2] = 1.0 / 3.0;

  wgll[2][0] = 1.0 / 3.0;
  wgll[2][1] = 4.0 / 3.0;
  wgll[2][2] = 1.0 / 3.0;

  wgll[3][0] = 1.0 / 6.0;
  wgll[3][1] = 5.0 / 6.0;
  wgll[3][2] = 5.0 / 6.0;
  wgll[3][3] = 1.0 / 6.0;

  wgll[4][0] = 1.0 / 10.0;
  wgll[4][1] = 49.0 / 90.0;
  wgll[4][2] = 32.00 / 45.00;
  wgll[4][3] = 49.0 / 90.0;
  wgll[4][4] = 1.0;
}

/* Find physical coordinates of Gauss points in each cell */
void GaussPoints(CELL *cell) {
  REAL xl, xr;
  UINT i;

  xl = cell->xl;
  xr = cell->xr;

  for (i = 0; i < cell->ng; i++)
    cell->xg[i] = 0.5 * (xl * (1.0 - xg[cell->ng - 1][i]) +
                         xr * (1.0 + xg[cell->ng - 1][i]));
}

void GetGaussPoints(double xl, double xr, int ng, double xgauss[ng]) {
  for (int i = 0; i < ng; i++)
    xgauss[i] = 0.5 * (xl * (1.0 - xg[ng - 1][i]) + xr * (1.0 + xg[ng - 1][i]));
}

/* Perform Gauss quadrature */
REAL GaussQuadrature(REAL *f, UINT ng) {
  REAL integral = 0.0;
  UINT i;

  for (i = 0; i < ng; i++)
    integral += f[i] * wg[ng - 1][i];

  return integral;
}

// Reinitialize Cell Values after change in dimensions
// This function only considers change in dimension of the cell.
// The following quantities do not change
// p, ng, ngll, active, lface, rface, lcell, rcell
// size of *xg
void initialize_cell_points(
  CELL *cell, double xl, double xr, double wl, double wr) {
  // Initializing cell co-ordinates
  cell->xl = xl;
  cell->xr = xr;
  cell->x = (xl + xr) / 2.0;

  // Cell Length
  cell->h = xr - xl;

  // Cell Boundary Velocity
  cell->wl = wl;
  cell->wr = wr;

  // Cell Gauss Points
  GaussPoints(cell);
}
