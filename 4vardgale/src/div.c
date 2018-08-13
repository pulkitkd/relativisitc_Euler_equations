#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include "dg.h"
#include "dg1d.h"
#include "math.h"

// TODO Fix the number of points when other methods are implemented
// Number of quadrature points for the specified predictor method
int n_time_quadrature_points(const int predictor) {
  switch (predictor) {
  case taylor_1p:
    return 1;
    break;
  case cerk_nodal_2p:
    return 2;
    break;
  case cerk_nodal_3p:
    return 2;
    break;
  default:
    printf("Predictor Method Not Supported \n");
    exit(1);
  }
}

// The values of the time quadrature points for the specified predictor method
void quadrature_points(const int n_points, double points[n_points]) {
  switch (n_points) {
  case 1:
    points[0] = 0.5;
    break;
  case 2:
    points[0] = (1 - 1.0 / sqrt(3.0)) / 2.0;
    points[1] = (1 + 1.0 / sqrt(3.0)) / 2.0;
    break;
  case 3:
    points[0] = (1 - 3.0 / sqrt(15.0)) / 2.0;
    points[1] = 1.00 / 2.0;
    points[2] = (1 + 3.0 / sqrt(15.0)) / 2.0;
    break;
  default:
    printf("%d Time Quadrature Points not Supported in Predictor \n", n_points);
    exit(1);
  }
}

// The values of qudrature weights for specified predictor method
void quadrature_weights(const int n_points, double w[n_points]) {
  switch (n_points) {
  case 1:
    w[0] = 2;
    break;
  case 2:
    w[0] = 1;
    w[1] = 1;
    break;
  case 3:
    w[0] = 5.00 / 9.00;
    w[1] = 8.00 / 9.00;
    w[2] = 5.00 / 9.00;
    break;
  default:
    printf("Quadrature Weights not defined for %d time points. \n", n_points);
    exit(1);
  }
}


void Flux(CELL *cell, FACE *face, const int predictor) {

  for (int i = 0; i < NCMAX; i++)
    if (cell[i].active)
      for (int j = 0; j < NVAR; j++)
        for (int k = 0; k < cell[i].p; k++)
          cell[i].Re[j][k] = 0.0;

  // We determine parameters for the predictor
  int n_time_points = n_time_quadrature_points(predictor);
  double t[n_time_points];
  quadrature_points(n_time_points, t);

  double w[n_time_points];
  quadrature_weights(n_time_points, w);

  // Variables to store the predicted values at the faces of cell
  double u_l[NCMAX][n_time_points][NVAR]; // Left Hand Flux for the cell
  double u_r[NCMAX][n_time_points][NVAR]; // Right Hand Flux for the cell

  // Generates predicted solution and flux function using the prediction
  for (int i = 0; i < NCMAX; ++i) {
    if (cell[i].active) {
      CELL *cell_p = &cell[i];

      // Number of Gauss Quadrature Points
      int ng = cell_p->ng;

      // Number of points at which we will be calculating our predictor.
      // We are adding two more points to the Gauss points , namely the
      // boundaries of the cell.
      int nodes = ng + 2;

      // The actual Gauss Quadrature points, our cells store this, we reference
      // them here
      double *xgauss = cell_p->xg;

      // However, we need to create a new array of points, which also contains
      // the endpoints.
      double x_nodes[nodes];
      x_nodes[0] = cell_p->xl;
      x_nodes[nodes - 1] = cell_p->xr;
      for (int j = 1; j < nodes - 1; ++j) {
        x_nodes[j] = xgauss[j - 1];
      }

      double u_p[nodes][n_time_points][NVAR];
      // We set the predictor value to zero (Initialization)
      for (int j = 0; j < nodes; ++j)
        for (int k = 0; k < n_time_points; ++k)
          for (int l = 0; l < NVAR; ++l)
            u_p[j][k][l] = 0;
//printf("Debug: 2.1 - %d\n", i);
      get_predictor(
        predictor, dt, n_time_points, t, cell_p, nodes, x_nodes, u_p);
//printf("Debug: 2.2 - %d\n", i);
      // We assign the predicted value on faces to the corresponding variables.
      for (int k = 0; k < n_time_points; ++k) {
        for (int j = 0; j < NVAR; ++j) {
          u_l[i][k][j] = u_p[0][k][j];
          u_r[i][k][j] = u_p[nodes - 1][k][j];
        }
      }

      for (int k = 0; k < n_time_points; ++k) {
        for (int j = 0; j < cell[i].ng; j++) {
          double vel = GetMeshVel(&cell[i], cell[i].xg[j]);
          double flg[NVAR];
          EulerFlux(u_p[j + 1][k], vel, flg);
          for (int m = 0; m < cell[i].p; ++m) {
            double vx = ShapeFunDeriv(cell[i].xg[j], &cell[i], m);
            for (int l = 0; l < NVAR; l++)
              cell[i].Re[l][m] -= 0.5 * cell[i].h * flg[l] * vx *
                                  wg[cell[i].ng - 1][j] * w[k] * 0.5;
          }
        }
      }
    }
  }

  /* Loop over cell faces and find flux, periodic bc */
  for (int i = 0; i < NFMAX; i++) {
    if (face[i].active) {

      // We get the pointers to the face and the cells on the left and right
      // side of the face.
      FACE *face_p = &face[i];
      ptrdiff_t lcell_index = face_p->lcell - cell;
      ptrdiff_t rcell_index = face_p->rcell - cell;
      double n_flux[NVAR];

      // We now iterate over the time points at which we have calculated the
      // predicted value.
      for (int k = 0; k < n_time_points; ++k) {
        // Mirror or reflecting bc
        if (i == 0 && bc_left == FIXED) {
          u_r[lcell_index][k][1] = -u_l[rcell_index][k][1];
        }
        if (i == NF - 1 && bc_right == FIXED) {
          u_l[rcell_index][k][1] = -u_r[lcell_index][k][1];
        }


        // TODO Comment here properly to explain the apparent reversal of
        // u_l and u_r, and why this choice is made instead of the other one.
        switch (FLUX) {
        case LF:
          LFFlux(u_r[lcell_index][k], u_l[rcell_index][k], face[i].w, n_flux);
          break;
        case ROE:
          RoeFlux(u_r[lcell_index][k], u_l[rcell_index][k], face[i].w, n_flux);
          break;
        case HLLC:
          HLLCFlux(u_r[lcell_index][k], u_l[rcell_index][k], face[i].w, n_flux);
          break;
        case ALEROE:
          ALERoeFlux(
            u_r[lcell_index][k], u_l[rcell_index][k], face[i].w, n_flux);
          break;
        case HLLCPS:
          HLLCPSFlux(
            u_r[lcell_index][k], u_l[rcell_index][k], face[i].w, n_flux);
          break;
        case HLL:
          HLLFlux(u_r[lcell_index][k], u_l[rcell_index][k], face[i].w, n_flux);
          break;
        default:
          printf("Error: Flux number %d not defined\n", FLUX);
          exit(0);
        }

        /* Add interface flux to the cells */
        if (i != 0)
          for (int j = 0; j < NVAR; j++)
            for (int m = 0; m < face[i].lcell->p; ++m) {
              double v = ShapeFun(face[i].x, face[i].lcell, m);
              face[i].lcell->Re[j][m] += n_flux[j] * v * w[k] * 0.5;
            }

        if (i != NF - 1)
          for (int j = 0; j < NVAR; j++)
            for (int m = 0; m < face[i].rcell->p; ++m) {
              double v = ShapeFun(face[i].x, face[i].rcell, m);
              face[i].rcell->Re[j][m] -= n_flux[j] * v * w[k] * 0.5;
            }
      }
    }
  }
}
