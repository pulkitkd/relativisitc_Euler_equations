#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "dg.h"
#include "dg1d.h"

void Update(CELL *cell, FACE *face) {
  UINT i, j, k;
  REAL U[NVAR], V[NVAR];
  for (i = 0; i < NCMAX; i++)
    if (cell[i].active) {
      for (j = 0; j < NVAR; j++)
        for (k = 0; k < cell[i].p; k++)
          cell[i].U[j][k] = cell[i].h * cell[i].U[j][k] - dt * cell[i].Re[j][k];
    }

  MoveGrid(cell, face);

  for (i = 0; i < NCMAX; i++)
    if (cell[i].active) {
      for (j = 0; j < NVAR; j++)
        for (k = 0; k < cell[i].p; k++)
          cell[i].U[j][k] /= cell[i].h;

       //Check positivity
// TODO (pulkit#1#): Check Positivity here and in Taylor Predictor

//       REAL x = (cell[i].xl + cell[i].xr) / 2;
//       Uvect(cell, x , U);
//       con2prim(U, V);
//       REAL rho = V[0];
//       REAL p = V[3];
//       if (rho < 0) {
//         printf("Density(prim) is negative in cell %d \n", i);
//         exit(0); }
//
//       if(p < 0) {
//         printf("Pressure(prim) is negative in cell %d \n", i);
//         exit(0); }

//      if (cell[i].U[0][0] < 0.0) {
//        printf("Density is negative in cell %d !!!\n", i);
//        exit(0);
//      }
//      p = (GAMMA - 1.0) *
//          (cell[i].U[3][0] - 0.5 * pow(cell[i].U[1][0], 2) / cell[i].U[0][0]);
//      if (p < 0.0) {
//        printf("Pressure is negative in cell %d !!!\n", i);
//        exit(0);
    //}


    }
}
