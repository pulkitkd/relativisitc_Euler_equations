#include "dg.h"
#include "dg1d.h"
#include <stdio.h>
#include <stdlib.h>

void UatGauss(CELL *cell, REAL **U) {
  UINT iv, ig, ip;

  for (ig = 0; ig < cell->ng; ig++)
    for (iv = 0; iv < NVAR; iv++) {
      U[ig][iv] = 0.0;
      for (ip = 0; ip < cell->p; ip++)
        U[ig][iv] += cell->U[iv][ip] * ShapeFun(cell->xg[ig], cell, ip);
    }
}

// Compute solution at GLL nodes
void UatGLL(CELL *cell, REAL **U) {
  UINT iv, ig, ip;
  REAL xi, x;

  for (ig = 0; ig < cell->ngll; ig++) {
    xi = xgll[cell->ngll - 1][ig];
    x = 0.5 * (1.0 - xi) * cell->xl + 0.5 * (1.0 + xi) * cell->xr;
    for (iv = 0; iv < NVAR; iv++) {
      U[ig][iv] = 0.0;
      for (ip = 0; ip < cell->p; ip++)
        U[ig][iv] += cell->U[iv][ip] * ShapeFun(x, cell, ip);
    }
  }
}

// Compute solution at given value of x
void Uvect(CELL *cell, REAL x, REAL *U) {
  UINT iv, ip;

  for (iv = 0; iv < NVAR; iv++) {
    U[iv] = 0.0;
    for (ip = 0; ip < cell->p; ip++)
      U[iv] += cell->U[iv][ip] * ShapeFun(x, cell, ip);
  }
}

// Matrices for Nodal to Modal Conversion
static const double invlegvander3[3][3] = {
  {0.277777777777777741, 0.444444444444444518, 0.27777777777777774},
  {-0.37267799624996494363, 0.0, 0.37267799624996494363},
  {0.24845199749997658797, -0.49690399499995317593, 0.24845199749997658797}};

static const double invlegvander4[4][4] = {{0.173927422568726937796443442174,
                                            0.326072577431273062203556557826,
                                            0.326072577431273062203556557826,
                                            0.173927422568726937796443442174},
                                           {-0.259418289292771246565670191469,
                                            -0.192012546066861755434737210425,
                                            0.192012546066861755434737210425,
                                            0.259418289292771246565670191469},
                                           {0.238144836103920094175365781383,
                                            -0.238144836103920094175365781383,
                                            -0.238144836103920094175365781383,
                                            0.238144836103920094175365781383},
                                           {-0.140235025813009784988733810551,
                                            0.35520060651490708071591257681,
                                            -0.35520060651490708071591257681,
                                            0.140235025813009784988733810551}};

// Convert nodal values to modal values in relation to the dgale1d problem.
void nodal_to_modal(int modes,
                    double nodal[modes][NVAR],
                    double modal[modes][NVAR]) {
  switch (modes) {
  case 3:
    for (int i = 0; i < NVAR; ++i) {
      for (int j = 0; j < modes; ++j) {
        modal[j][i] = 0;
        for (int k = 0; k < modes; ++k) {
          modal[j][i] += invlegvander3[j][k] * nodal[k][i];
        }
      }
    }
    break;
  case 4:
    for (int i = 0; i < NVAR; ++i) {
      for (int j = 0; j < modes; ++j) {
        modal[j][i] = 0;
        for (int k = 0; k < modes; ++k) {
          modal[j][i] += invlegvander4[j][k] * nodal[k][i];
        }
      }
    }
    break;
  default:
    printf("No Nodal to Modal Converter for %d modes. \n", modes);
    exit(1);
  }
}
