#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "dg.h"
#include "dg1d.h"

void Result(CELL *cell,
             double time, int iter, int final) {
  FILE *fp1, *fp2, *fp3, *fp4;
  UINT i, j;
  REAL dx, x, U[NVAR], d, u1, u2, p, c, m1, m2, w, h,
        V[NVAR], Uavg[NVAR], Vavg[NVAR];

  char solbuf[20];
  char avgbuf[20];
  char hbuf[20];

  // We see if a directory exists, and if it doesn't we create it.
  struct stat st = {0};

  if (stat("video", &st) == -1) {
    mkdir("video", 0700);
  }

  snprintf(solbuf, sizeof(solbuf), "video/%d.sol", iter);
  snprintf(avgbuf, sizeof(avgbuf), "video/%d.avg", iter);
  snprintf(hbuf, sizeof(hbuf), "video/%d.h", iter);

  if (final == 1) {
    snprintf(solbuf, sizeof(solbuf), "sol");
    snprintf(avgbuf, sizeof(avgbuf), "avg");
    snprintf(hbuf, sizeof(hbuf), "h");
  }

  fp1 = fopen(avgbuf, "w");
  fp2 = fopen(solbuf, "w");
  fp3 = fopen(hbuf, "w");// TODO (pulkit#1#): compute average solution in terms of primitives
// TODO (pulkit#2#): Plot velocity as sqrt(u1 * u1 + u2 * u2).
//                   Presently velocity is only being plotted as u1.


  for (i = 0; i < NCMAX; i++)
    if (cell[i].active) {
      // average solution
      Uavg[0] = cell[i].U[0][0];
      Uavg[1] = cell[i].U[1][0];
      Uavg[2] = cell[i].U[2][0];
      Uavg[3] = cell[i].U[3][0];

      con2prim(Uavg, Vavg);

      fprintf(fp1, "%e %e %e %e %e \n",
                 cell[i].x, Vavg[0], Vavg[1], Vavg[3], Vavg[2]);

      // more detailed solution evaluated at NPLT points inside cell
      if (NPLT == 1)
        dx = 0.5 * cell[i].h;
      else
        dx = cell[i].h / (NPLT - 1);
      for (j = 0; j < NPLT; j++) {
        if (NPLT == 1)
          x = cell[i].xl + dx;
        else
          x = cell[i].xl + dx * j;
        Uvect(&cell[i], x, U);
//        d = U[0];
//        u1 = U[1] / d;
//        u2 = U[2] / d;
//        p = (GAMMA - 1.0) * (U[3] - 0.5 * d * u1 * u1);
        con2prim(U, V);
        d = V[0];
        u1 = V[1];
        u2 = V[2];
        p = V[3];
        h = 1 + p * GAMMA / (d * (GAMMA -1));
        c = sqrt(GAMMA * p / (h * d));
        m1 = u1 / c;
        m2 = u2 / c;
        w = GetMeshVel(&cell[i], x);
        fprintf(fp2, "%e %e %e %e %e %e %e %e \n", x, d, u1, p, m1, m2, w, u2);
        //printf("in result: rho=%lf u1=%e u2=%e p=%lf \n",d, u1, u2, p);
      }
      fprintf(fp2, "\n");

      // mesh size, mesh velocity
      fprintf(fp3,
              "%d %e %e %e\n",
              i,
              cell[i].x,
              cell[i].h,
              0.5 * (cell[i].wl + cell[i].wr));
    }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);

  // Error computation is done only for problems with exact solution
  if (test_case == PULSE) {

    double error[NVAR];
    for (i = 0; i < NVAR; ++i) {
      error[i] = 0;
    }

    fp4 = fopen("error", "a");
    l2error(cell, time, error, error_order);
    fprintf(fp4,
            "%d, %d, %.16e, %.16e, %.16e\n",
            NC,
            predictor_method,
            error[0],
            error[1],
            error[2]);
    fclose(fp4);
  }
}
