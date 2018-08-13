#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dg.h"
#include "dg1d.h"

void get_predictor(int predictor,
                   double tstep,
                   int n_time_points,
                   double t[n_time_points],
                   CELL *cell_p,
                   int nodes,
                   double x_nodes[nodes],
                   double u_p[nodes][n_time_points][NVAR]) {
  switch (predictor) {
  case taylor_1p:
    taylor_1p_predictor(tstep, n_time_points, t, cell_p, nodes, x_nodes, u_p);
    break;
  case cerk_nodal_2p:
    cerk_nodal_2p_predictor(
      tstep, n_time_points, t, cell_p, nodes, x_nodes, u_p);
    break;
  case cerk_nodal_3p:
    cerk_nodal_3p_predictor(
      tstep, n_time_points, t, cell_p, nodes, x_nodes, u_p);
    break;
  default:
    printf("Unknown Predictor Method \n");
    exit(1);
  }
}
