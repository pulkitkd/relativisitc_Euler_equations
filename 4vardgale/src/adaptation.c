#include "dg.h"
#include "dg1d.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

enum {
  no_adapt,
  merge_to_left,
  merge_to_right,
  refine_to_left,
  refine_to_right
};

// Adaptation Enumeration
enum { no_change, merge, refine };

// Function to merge cells. It takes two cells, lcell and rcell as inputs.
// If dest_cell is equal to merge_right, it merges into rcell, else if the
// dest_cell is equal to merge_left, it merges into left cell.
void merge_cell(CELL *lcell, CELL *rcell, int dest_cell);

void refine_cell(CELL *dlcell, FACE *nface, CELL *drcell, int refine_status);

void adapt(int nc, CELL cell[nc], int nf, FACE face[nf]) {
  double dx_min = hmin;
  double dx_max = 2 * dx_init;

  // sweep 1: based on dxmin, dxmax, level
  for (int i = 0; i < nc; ++i) {
    if (cell[i].active) {
      if (cell[i].h < dx_min) {
        cell[i].do_adapt = merge;
      } else if (cell[i].h > dx_max) {
        cell[i].do_adapt = refine;
        // more aggressive: marks both neighbours to be refined
        // cell[i].lcell->do_adapt = refine;
        // cell[i].rcell->do_adapt = refine;
      } else if (cell[i].level < cell[i].lcell->level &&
                 cell[i].level < cell[i].rcell->level) {
        cell[i].do_adapt = refine;
      } else
        cell[i].do_adapt = no_change;
    }
  }

  // sweep 2
  for (int i = 0; i < nc; ++i)
    if (cell[i].active) {
      if (cell[i].lcell->do_adapt == refine &&
          cell[i].rcell->do_adapt == refine)
        cell[i].do_adapt = refine;

      if (cell[i].do_adapt != refine) {
        double hl = cell[i].lcell->h;
        if (cell[i].lcell->do_adapt == refine)
          hl *= 0.5;
        if (cell[i].h > 2.0 * hl && cell[i].h > 2.0 * dx_min)
          cell[i].do_adapt = refine;

        double hr = cell[i].rcell->h;
        if (cell[i].rcell->do_adapt == refine)
          hr *= 0.5;
        if (cell[i].h > 2.0 * hr && cell[i].h > 2.0 * dx_min)
          cell[i].do_adapt = refine;

        // what if lcell or rcell is flagged for coarsening ?
      }
    }

  // sweep 3:
  for (int i = 0; i < nc; ++i)
    if (cell[i].active) {
      if (cell[i].lcell->do_adapt == refine &&
          cell[i].rcell->do_adapt == refine)
        cell[i].do_adapt = refine;
      if (cell[i].lcell->do_adapt == merge || cell[i].rcell->do_adapt == merge)
        if (cell[i].do_adapt == refine)
          cell[i].do_adapt = no_change;
    }

  int adapt_action[nc];

  for (int i = 0; i < nc; ++i)
    adapt_action[i] = no_adapt;

  for (int i = 0; i < nc; ++i) {
    if (cell[i].active) {
      switch (cell[i].do_adapt) {
      case no_change:
        adapt_action[i] = no_adapt;
        break;
      case merge:
        if (i == 0)
          adapt_action[i] = merge_to_left;
        else if (i == NC - 1)
          adapt_action[i] = merge_to_right;
        else {
          if (cell[i].lcell->h < cell[i].rcell->h)
            adapt_action[i] = merge_to_left;
          else
            adapt_action[i] = merge_to_right;
        }
        break;
      case refine:
        if (i == 0)
          adapt_action[i] = refine_to_right;
        else
          adapt_action[i] = refine_to_left;
        break;
      default:
        printf("Invalid Adapt Status.\n");
        exit(1);
      }
    }
  }

  for (int i = 0; i < nc; ++i) {
    if (cell[i].active) {
      switch (adapt_action[i]) {
      case no_adapt:
        break;
      case merge_to_left:
        merge_cell(cell[i].lcell, &cell[i], merge_to_left);
        break;
      case merge_to_right:
        merge_cell(&cell[i], cell[i].rcell, merge_to_right);
        break;
      case refine_to_left:
        if (free_cell_index < NCMAX && free_face_index < NFMAX) {
          CELL *dcell = &cell[free_cell_index];
          ++free_cell_index;
          FACE *nface = &face[free_face_index];
          ++free_face_index;
          refine_cell(dcell, nface, &cell[i], refine_to_left);
        } else
          printf("WARNING: Out of extra refinement cells. Ignoring refine.\n");
        break;
      case refine_to_right:
        if (free_cell_index < NCMAX && free_face_index < NFMAX) {
          CELL *dcell = &cell[free_cell_index];
          ++free_cell_index;
          FACE *nface = &face[free_face_index];
          ++free_face_index;
          refine_cell(&cell[i], nface, dcell, refine_to_right);
        } else
          printf("WARNING: Out of extra refinement cells. \n");
        break;
      default:
        printf("Invalid Adapt Action.\n");
        exit(1);
      }
    }
  }
}

void refine_cell(CELL *dlcell, FACE *nface, CELL *drcell, int refine_status) {
  CELL *pcell;
  switch (refine_status) {
  case refine_to_left:
    pcell = drcell;
    break;
  case refine_to_right:
    pcell = dlcell;
    break;
  default:
    printf("Wrong Refine Type \n");
    exit(1);
  }

  // Now, we get the information about the cells
  int p = pcell->p;
  int ng = pcell->ng;
  int ngll = pcell->ngll;
  double xl = pcell->xl;
  double xr = pcell->xr;
  double x = pcell->x;
  double wl = pcell->wl;
  double wr = pcell->wr;
  double h = pcell->h;

  // Now, we define the co-ordinates for the left and right cell
  double lxl = xl;
  double lxr = x;
  double rxl = x;
  double rxr = xr;

  double lwl = wl;
  double lwr = 1 / h * ((x - xl) * wr + (xr - x) * wl);
  double rwl = lwr;
  double rwr = wr;

  // We now calculate the gauss points for the left and right
  double lxg[ng];
  double rxg[ng];

  GetGaussPoints(lxl, lxr, ng, lxg);
  GetGaussPoints(rxl, rxr, ng, rxg);

  //   printVector(ng, rxg);
  // Next, we calculate the solution at the gauss points
  double ul[ng][NVAR];
  double ur[ng][NVAR];

  for (int i = 0; i < ng; ++i) {
    Uvect(pcell, lxg[i], ul[i]);
    Uvect(pcell, rxg[i], ur[i]);
  }

  // Having got the solutions, we now initialize the cells
  initialize_cell_points(dlcell, lxl, lxr, lwl, lwr);
  initialize_cell_points(drcell, rxl, rxr, rwl, rwr);
  dlcell->ngll = ngll;
  drcell->ngll = ngll;

  nface->w = lwr;
  nface->x = x;

  // Next, we set parameters for the cells
  CELL *lcell = pcell->lcell;
  CELL *rcell = pcell->rcell;
  FACE *lface = pcell->lface;
  FACE *rface = pcell->rface;

  lcell->rcell = dlcell;
  lface->rcell = dlcell;
  rcell->lcell = drcell;
  rface->lcell = drcell;

  dlcell->lcell = lcell;
  dlcell->lface = lface;
  dlcell->rcell = drcell;
  dlcell->rface = nface;

  drcell->lcell = dlcell;
  drcell->lface = nface;
  drcell->rface = rface;
  drcell->rcell = rcell;

  nface->lcell = dlcell;
  nface->rcell = drcell;

  // Now, we transfer the solution
  for (int i = 0; i < NVAR; ++i) {
    for (int j = 0; j < p; ++j) {
      dlcell->U[i][j] = 0;
      drcell->U[i][j] = 0;
      for (int k = 0; k < ng; ++k) {
        dlcell->U[i][j] +=
          ul[k][i] * ShapeFun(lxg[k], dlcell, j) * wg[ng - 1][k] * 0.5;
        drcell->U[i][j] +=
          ur[k][i] * ShapeFun(rxg[k], drcell, j) * wg[ng - 1][k] * 0.5;
      }
    }
  }

  dlcell->level = pcell->level + 1;
  drcell->level = pcell->level + 1;
  pcell->active = false;
  nface->active = true;
  drcell->active = true;
  dlcell->active = true;
}

void merge_cell(CELL *lcell, CELL *rcell, int dest_cell) {
  // Set which cell is to change
  CELL *dcell;
  switch (dest_cell) {
  case merge_to_right:
    dcell = rcell;
    break;
  case merge_to_left:
    dcell = lcell;
    break;
  default:
    printf("Wrong Cell Merge Type.\n");
    exit(1);
  }

  // Now, we calculate the solution for the merged cell
  int p = lcell->p;
  int ng = lcell->ng;
  double hl = lcell->h;
  double hr = rcell->h;
  double xl = lcell->xl;
  double xr = rcell->xr;
  double wl = lcell->wl;
  double wr = rcell->wr;

  // Gauss Point in Global Co-Ordinate System
  double lxg[ng];
  double rxg[ng];

  for (int j = 0; j < ng; ++j) {
    lxg[j] = lcell->xg[j];
    rxg[j] = rcell->xg[j];
  }

  // Solution values from the cells
  double lbeta[NVAR][ng];
  double rbeta[NVAR][ng];
  for (int j = 0; j < ng; ++j) {
    double Ul[NVAR];
    double Ur[NVAR];
    Uvect(lcell, lcell->xg[j], Ul);
    Uvect(rcell, rcell->xg[j], Ur);
    for (int i = 0; i < NVAR; ++i) {
      lbeta[i][j] = hl * Ul[i];
      rbeta[i][j] = hr * Ur[i];
    }
  }

  // Having got the required values from the system, we now overwrite the
  // dell with the new values
  initialize_cell_points(dcell, xl, xr, wl, wr);

  for (int i = 0; i < NVAR; ++i) {
    for (int j = 0; j < p; ++j) {
      dcell->U[i][j] = 0;
    }
  }

  for (int i = 0; i < NVAR; ++i) {
    for (int j = 0; j < p; ++j) {
      dcell->U[i][j] = 0;
      for (int k = 0; k < p; ++k) {
        dcell->U[i][j] +=
          wg[ng - 1][k] * (lbeta[i][k] * ShapeFun(lxg[k], dcell, j) +
                           rbeta[i][k] * ShapeFun(rxg[k], dcell, j));
      }
    }
  }

  for (int i = 0; i < NVAR; ++i) {
    for (int j = 0; j < p; ++j) {
      dcell->U[i][j] *= 1 / (2 * (hl + hr));
    }
  }
  switch (dest_cell) {
  case merge_to_right:
    lcell->active = 0;
    lcell->rface->active = 0;
    rcell->lcell = lcell->lcell;
    rcell->lface = lcell->lface;
    rcell->lface->rcell = rcell;
    lcell->lcell->rcell = rcell;
    break;
  case merge_to_left:
    rcell->active = 0;
    rcell->lface->active = 0;
    lcell->rcell = rcell->rcell;
    lcell->rface = rcell->rface;
    lcell->rface->lcell = lcell;
    rcell->rcell->lcell = lcell;
    break;
  default:
    printf("Wrong Cell Merge Type.\n");
    exit(1);
  }
}
