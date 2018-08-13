#include <math.h>
#include "dg.h"
#include "dg1d.h"

void TimeStep(CELL *cell) {
  UINT i;
  REAL dtc, w, ucbyh, ubyh, MaxEigen;
  REAL U[NVAR], V[NVAR], Eigenvalues[NVAR];
  REAL alpha = GAMMA / (GAMMA -1);
  REAL beta = 0.1;
  dt = 1.0e20;
  ucbyh = 0.0;
  ubyh = 0.0;

  for (i = 0; i < NCMAX; i++)
    if (cell[i].active) {
      REAL x = (cell[i].xl + cell[i].xr) / 2;
      Uvect(cell, x , U);
      con2prim(U, V);
      REAL rho = V[0];
      REAL u1 = V[1];
      REAL u2 = V[2];
      REAL p = V[3];
      REAL h = 1 + p * alpha / rho;
      REAL c = sqrt(GAMMA * p / (rho * h));
      REAL u1sq = u1 * u1;
      REAL u2sq = u2 * u2;
      REAL usq = u1sq + u2sq;
      Eigenvalues[0] = fabs(u1);
      Eigenvalues[1] = fabs(u1);
      Eigenvalues[2] = fabs((u1 * (1 - c * c) +
                 c * sqrt((1 - usq) * (1 - u1sq - u2sq * c * c))) /
                 (1 - usq * c * c));
      Eigenvalues[2] = fabs((u1 * (1 - c * c) -
                 c * sqrt((1 - usq) * (1 - u1sq - u2sq * c * c))) /
                 (1 - usq * c * c));



      w = 0.5 * (cell[i].wl + cell[i].wr);
      // Check for the largest eigenvalue and store it in MaxEigen
      MaxEigen = fabs(Eigenvalues[0] - w);
  //    printf("Debug 4 \n");
      for(int j=1; j < NVAR; j++)
        if(fabs(Eigenvalues[j] - w) > MaxEigen)
          MaxEigen = fabs(Eigenvalues[j] - w);

//printf("Debug: MaxEigen %lf   dt %lf \n", MaxEigen, dt);

      dtc = cfl * cell[i].h / MaxEigen;
      dt = MIN(dtc, dt);
//printf("Debug: dtc %lf   dt %lf \n", dtc, dt);
      dtc = beta * cell[i].h / (fabs(cell[i].wr - cell[i].wl) + 1.0e-14);
      dt = MIN(dtc, dt);
      //printf("Debug: MaxEigen %lf dtc %lf   dt %lf \n", MaxEigen, dtc, dt);
      // to compute usual cfl number
      //ucbyh = MAX((fabs(u1) + c) / cell[i].h, ucbyh);
      //ubyh = MAX(fabs(u1) / cell[i].h, ubyh);
    }
//printf("Debug 6 \n");
  // This is the usual cfl number
  //cfl_f = dt * ucbyh * (2 * (PORD - 1) + 1); // full cfl
  //cfl_u = dt * ubyh * (2 * (PORD - 1) + 1);  // advective cfl
}
