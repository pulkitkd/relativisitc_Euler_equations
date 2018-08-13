#include "dg.h"
#include "dg1d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

void exact_solution(double x, double t, double U[NVAR]) {
  double V[NVAR];
  if (test_case == SOD || test_case == LAX || test_case == MSOD ||
      test_case == SSOD || test_case == LOWD) {
    InitCondShocktube(x, V);
  } else if (test_case == BLAST) {
    InitCondBlast(x, V);
  } else if (test_case == SHUOSHER) {
    InitCondShuOsher(x, V);
  } else if (test_case == PULSE) {
    double effx = x - t;
    InitCondPulse(effx, V);
  } else if (test_case == POLY) {
    double effx = x - t;
    InitCondPoly(effx, V);
  } else {
    printf("Error: Unknown test case\n");
    exit(0);
  }
  U[0] = V[0];
  U[1] = V[0] * V[1];
  U[2] = V[0] * V[2];
  U[3] = V[3] / (GAMMA - 1.0) + 0.5 * V[0] * (V[1] * V[1]);
}

void l2error(CELL *cell, double time, double error[NVAR], int order) {
  for (int i = 0; i < NVAR; ++i)
    error[i] = 0;

  double max_error = -20;
  for (int i = 0; i < NCMAX; ++i)
    if (cell[i].active) {
      double error_xg[order];
      double error_wg[order];
      double xl = cell[i].xl;
      double xr = cell[i].xr;

      for (int j = 0; j < order; ++j) {
        error_xg[j] =
          0.5 * (xl * (1.0 - xg[order - 1][j]) + xr * (1.0 + xg[order - 1][j]));
        error_wg[j] = wg[order - 1][j];
      }
      double cell_error[NVAR];
      for (int k = 0; k < NVAR; ++k)
        cell_error[k] = 0;

      for (int j = 0; j < order; ++j) {
        double exact[NVAR];
        double solution[NVAR];
        exact_solution(error_xg[j], time, exact);
        Uvect(&cell[i], error_xg[j], solution);
        for (int k = 0; k < NVAR; ++k) {
          max_error = fmax(max_error, fabs(exact[k] - solution[k]));
          cell_error[k] +=
            pow(exact[k] - solution[k], 2) * error_wg[j] * cell[i].h * 0.5;
        }
      }
      for (int k = 0; k < NVAR; ++k) {
        error[k] += cell_error[k];
      }
    }
  for (int i = 0; i < NVAR; ++i)
    error[i] = sqrt(error[i]);
  //   printf("Max Error = %.16e \n", max_error);
}
