#include <math.h>
#include "dg.h"
#include "dg1d.h"
#include "utils.h"

// Mid point in time, so n_time_points = 1
// Taylor expansion is around (cell->x, t_n)
void taylor_1p_predictor(double tstep,
                         int n_time_points,
                         double t[n_time_points],
                         CELL *cell,
                         int nodes,
                         double x_nodes[nodes],
                         double u_p[nodes][n_time_points][NVAR]) {

  REAL Vl[NVAR], Vr[NVAR], U[NVAR], V[NVAR];
//printf("Debug: 2.1.1  \n");

  Uvect(cell, x_nodes[0], U);
  con2prim(U, Vl);
//printf("Vl: rho = %lf u1 = %lf u2 = %lf p = %lf \n",
//         Vl[0], Vl[1], Vl[2], Vl[3]);

  Uvect(cell, x_nodes[nodes - 1], U);
  con2prim(U, Vr);
//printf("Vr: rho = %lf u1 = %lf u2 = %lf p = %lf \n",
//         Vr[0], Vr[1], Vr[2], Vr[3]);

  for(int i = 0; i < NVAR; ++i)
    V[i] = 0.5 * (Vl[i] + Vr[i]);
  //nodes = ng + 2 = ng + faces

  double A[NVAR][NVAR];
  Jacobian(V, A);


//printf("rho = %lf u1 = %lf u2 = %lf p = %lf \n",
//           vb[0], vb[1], vb[2], vb[3]);
//printf("d = %lf m1 = %lf m2 = %lf e = %lf \n",
//           ub[0], ub[1], ub[2], ub[3]);
  // compute gradient of solution. it is constant in cell since
  // solution is linear.
  double v_x[NVAR];
  for (int i = 0; i < NVAR; ++i)
    v_x[i] = (Vr[i] - Vl[i]) / cell->h;
 //   printf("Gradients: %lf, %lf, %lf, %lf \n",
//         v_x[0], v_x[1], v_x[2], v_x[3]);

  double dtp = tstep * t[0];

  for (int k = 0; k < nodes; ++k) {

    // x_nodes[k] are at time t_n, so we have to move them
    double w = GetMeshVel(cell, x_nodes[k]);
    double x = x_nodes[k] + w * dtp;
    double dx = x - cell->x;

    // Finally, we compute the predictor
    for (int i = 0; i < NVAR; ++i) {
      u_p[k][0][i] = V[i] + dx * v_x[i];
      for (int j = 0; j < NVAR; ++j) {
        u_p[k][0][i] -= dtp * A[i][j] * v_x[j];

      }
    }
//printf("predicted Sol : rho = %e u1 = %e u2 = %e p = %e \n",
//           u_p[k][0][0], u_p[k][0][1], u_p[k][0][2], u_p[k][0][3]);
/* TODO (pulkit#1#): Correct  checkPositivity arguments and implementation*/


//    if (checkPositivity(u_p[k][0]) == FALSE) {
//      for (int i = 0; i < NVAR; ++i)
//        u_p[k][0][i] = V[i];
//        printf("rho = %e u1 = %e u2 = %e p = %e \n",
//           u_p[k][0][0], u_p[k][0][1], u_p[k][0][2], u_p[k][0][3]);
//    }
  }
}
