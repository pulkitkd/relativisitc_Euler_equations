#include "dg.h"
#include "dg1d.h"
#include <math.h>
#include "utils.h"

// The Second Order Method CERK Method
static const double A2[2][2] = {{0, 0}, {1, 0}};
static const double B2[2][3] = {{0, 1, -0.5}, {0, 0, 0.5}};
static const double C2[2] = {0, 1.0};

// The Third Order Method CERK Method
static const double A3[4][4] = {
  {0, 0, 0, 0},
  {12.00 / 23.00, 0, 0, 0},
  {-68.00 / 375.00, 368.00 / 375.00, 0, 0},
  {31.00 / 144.00, 529.00 / 1152.00, 125.00 / 384.00, 0}};

static const double B3[4][4] = {{0.00, 1.00, -65.00 / 48.00, 41.00 / 72.00},
                                {0.00, 0.00, 529.00 / 384.00, -529.00 / 576.00},
                                {0.00, 0.00, 125.00 / 128.00, -125.00 / 192.00},
                                {0.00, 0.00, -1.00, 1.00}};

static const double C3[4] = {0, 12.00 / 23.00, 300.00 / 375.00, 1.0};

void nodal_flux(CELL *cell_p,
                double x,
                int modes,
                double U_m[modes][NVAR],
                double U[NVAR],
                double K[NVAR],
                double h) {
  double A[NVAR][NVAR];
  Jacobian(U, A);

  // mesh velocity at x
  double w = GetMeshVel(cell_p, x);

  double u_x[NVAR];
  for (int i = 0; i < NVAR; ++i) {
    A[i][i] -= w;
    u_x[i] = 0.0;
  }

  // compute gradient of solution at x
  for (int i = 0; i < cell_p->p; i++) {
    double vx = ShapeFunDeriv(x, cell_p, i) / h * cell_p->h;
    for (int j = 0; j < NVAR; j++)
      u_x[j] += vx * U_m[i][j];
  }

  for (int i = 0; i < NVAR; ++i) {
    K[i] = 0;
    for (int j = 0; j < NVAR; ++j) {
      K[i] -= A[i][j] * u_x[j];
    }
  }
}

void cerk_nodal_2p_predictor(double tstep,
                             int n_time_points,
                             double t[n_time_points],
                             CELL *cell_p,
                             int nodes,
                             double x_nodes[nodes],
                             double u_p[nodes][n_time_points][NVAR]) {

  // We need to know the cell length at the given times in order to calculate
  // the derivative properly.
  const int stages = 2;
  double h[stages];
  for (int i = 0; i < stages; ++i)
    h[i] = cell_p->h + (cell_p->wr - cell_p->wl) * C2[i] * tstep;

  int modes = nodes - 2;
  double um[modes][NVAR];

  double u0[modes][NVAR];
  for (int k = 0; k < modes; ++k) {
    Uvect(cell_p, x_nodes[k + 1], u0[k]);
  }

  nodal_to_modal(modes, u0, um);

  double K0[modes][NVAR];
  for (int k = 0; k < modes; ++k) {
    nodal_flux(cell_p, x_nodes[k + 1], modes, um, u0[k], K0[k], h[0]);
  }

  double u1[modes][NVAR];
  for (int i = 0; i < modes; ++i) {
    for (int q = 0; q < NVAR; ++q) {
      u1[i][q] = u0[i][q] + tstep * A2[1][0] * K0[i][q];
    }
  }

  double K1[modes][NVAR];
  nodal_to_modal(modes, u1, um);
  for (int k = 0; k < modes; ++k) {
    nodal_flux(cell_p, x_nodes[k + 1], modes, um, u1[k], K1[k], h[1]);
  }

  double b[n_time_points][stages];

  for (int s = 0; s < n_time_points; ++s) {
    b[s][0] = tstep * (B2[0][0] + B2[0][1] * t[s] + B2[0][2] * pow(t[s], 2));
    b[s][1] = tstep * (B2[1][0] + B2[1][1] * t[s] + B2[1][2] * pow(t[s], 2));
  }

  for (int i = 0; i < modes; ++i) {
    for (int s = 0; s < n_time_points; ++s) {
      for (int q = 0; q < NVAR; ++q) {
        u_p[i + 1][s][q] = u0[i][q] + b[s][0] * K0[i][q] + b[s][1] * K1[i][q];
      }
    }
  }
  // We next want to evaluate the solution at the end points of the cell from
  // the nodal solution, this requires us to calculate the solution at the nodes
  // first and then transform it into modal solutions, which we can then use
  // to calculate the solutions at the cell boundaries.
  double un[modes][NVAR];
  for (int s = 0; s < n_time_points; ++s) {
    for (int i = 0; i < modes; ++i) {
      for (int j = 0; j < NVAR; ++j) {
        un[i][j] = u_p[i + 1][s][j];
      }
    }
    nodal_to_modal(modes, un, um);
    // We now have stored the modal solution in um, from which, we now
    // extract solution at the end points
    for (int j = 0; j < NVAR; ++j) {
      u_p[0][s][j] = 0;
      u_p[nodes - 1][s][j] = 0;
      for (int i = 0; i < modes; ++i) {
        u_p[0][s][j] += um[i][j] * ShapeFun(cell_p->xl, cell_p, i);
        u_p[nodes - 1][s][j] += um[i][j] * ShapeFun(cell_p->xr, cell_p, i);
      }
    }
  }

  // check positivity of all states
  int pos = TRUE;
  for (int s = 0; s < n_time_points; ++s)
    for (int i = 0; i < nodes; ++i) {
      if (checkPositivity(u_p[i][s]) == FALSE)
        pos = FALSE;
    }

  // if positivity lost at any point, make first order, copy cell average
  if (pos == FALSE)
    for (int s = 0; s < n_time_points; ++s)
      for (int i = 0; i < nodes; ++i)
        for (int j = 0; j < NVAR; ++j)
          u_p[i][s][j] = cell_p->U[j][0];
}

void cerk_nodal_3p_predictor(double tstep,
                             int n_time_points,
                             double t[n_time_points],
                             CELL *cell_p,
                             int nodes,
                             double x_nodes[nodes],
                             double u_p[nodes][n_time_points][NVAR]) {

  // We need to know the cell length at the given times in order to calculate
  // the derivative properly.
  const int stages = 4;
  double h[stages];
  for (int i = 0; i < stages; ++i)
    h[i] = cell_p->h + (cell_p->wr - cell_p->wl) * C3[i] * tstep;

  int modes = nodes - 2;
  double um[modes][NVAR];

  double u0[modes][NVAR];
  for (int k = 0; k < modes; ++k) {
    Uvect(cell_p, x_nodes[k + 1], u0[k]);
  }

  nodal_to_modal(modes, u0, um);

  double K0[modes][NVAR];
  for (int k = 0; k < modes; ++k) {
    nodal_flux(cell_p, x_nodes[k + 1], modes, um, u0[k], K0[k], h[0]);
  }

  double u1[modes][NVAR];
  for (int i = 0; i < modes; ++i) {
    for (int q = 0; q < NVAR; ++q) {
      u1[i][q] = u0[i][q] + tstep * A3[1][0] * K0[i][q];
    }
  }

  double K1[modes][NVAR];
  nodal_to_modal(modes, u1, um);
  for (int k = 0; k < modes; ++k) {
    nodal_flux(cell_p, x_nodes[k + 1], modes, um, u1[k], K1[k], h[1]);
  }

  double u2[modes][NVAR];
  for (int i = 0; i < modes; ++i) {
    for (int q = 0; q < NVAR; ++q) {
      u2[i][q] =
        u0[i][q] + tstep * A3[2][0] * K0[i][q] + tstep * A3[2][1] * K1[i][q];
    }
  }

  double K2[modes][NVAR];
  nodal_to_modal(modes, u2, um);
  for (int k = 0; k < modes; ++k) {
    nodal_flux(cell_p, x_nodes[k + 1], modes, um, u2[k], K2[k], h[2]);
  }

  double u3[modes][NVAR];
  for (int i = 0; i < modes; ++i) {
    for (int q = 0; q < NVAR; ++q) {
      u3[i][q] = u0[i][q] + tstep * A3[3][0] * K0[i][q] +
                 tstep * A3[3][1] * K1[i][q] + tstep * A3[3][2] * K2[i][q];
    }
  }

  double K3[modes][NVAR];
  nodal_to_modal(modes, u3, um);
  for (int k = 0; k < modes; ++k) {
    nodal_flux(cell_p, x_nodes[k + 1], modes, um, u3[k], K3[k], h[3]);
  }

  double b[n_time_points][stages];

  for (int s = 0; s < n_time_points; ++s) {
    b[s][0] = tstep * (B3[0][0] + B3[0][1] * t[s] + B3[0][2] * pow(t[s], 2) +
                       B3[0][3] * pow(t[s], 3));
    b[s][1] = tstep * (B3[1][0] + B3[1][1] * t[s] + B3[1][2] * pow(t[s], 2) +
                       B3[1][3] * pow(t[s], 3));
    b[s][2] = tstep * (B3[2][0] + B3[2][1] * t[s] + B3[2][2] * pow(t[s], 2) +
                       B3[2][3] * pow(t[s], 3));
    b[s][3] = tstep * (B3[3][0] + B3[3][1] * t[s] + B3[3][2] * pow(t[s], 2) +
                       B3[3][3] * pow(t[s], 3));
  }

  for (int i = 0; i < modes; ++i) {
    for (int s = 0; s < n_time_points; ++s) {
      for (int q = 0; q < NVAR; ++q) {
        u_p[i + 1][s][q] = u0[i][q] + b[s][0] * K0[i][q] + b[s][1] * K1[i][q] +
                           b[s][2] * K2[i][q] + b[s][3] * K3[i][q];
      }
    }
  }
  // We next want to evaluate the solution at the end points of the cell from
  // the nodal solution, this requires us to calculate the solution at the nodes
  // first and then transform it into modal solutions, which we can then use
  // to calculate the solutions at the cell boundaries.
  double un[modes][NVAR];
  for (int s = 0; s < n_time_points; ++s) {
    for (int i = 0; i < modes; ++i) {
      for (int j = 0; j < NVAR; ++j) {
        un[i][j] = u_p[i + 1][s][j];
      }
    }
    nodal_to_modal(modes, un, um);
    // We now have stored the modal solution in um, from which, we now
    // extract solution at the end points
    for (int j = 0; j < NVAR; ++j) {
      u_p[0][s][j] = 0;
      u_p[nodes - 1][s][j] = 0;
      for (int i = 0; i < modes; ++i) {
        u_p[0][s][j] += um[i][j] * ShapeFun(cell_p->xl, cell_p, i);
        u_p[nodes - 1][s][j] += um[i][j] * ShapeFun(cell_p->xr, cell_p, i);
      }
    }
  }

  // check positivity of all states
  int pos = TRUE;
  for (int s = 0; s < n_time_points; ++s)
    for (int i = 0; i < nodes; ++i) {
      if (checkPositivity(u_p[i][s]) == FALSE)
        pos = FALSE;
    }

  // if positivity lost at any point, make first order, copy cell average
  if (pos == FALSE)
    for (int s = 0; s < n_time_points; ++s)
      for (int i = 0; i < nodes; ++i)
        for (int j = 0; j < NVAR; ++j)
          u_p[i][s][j] = cell_p->U[j][0];
}
