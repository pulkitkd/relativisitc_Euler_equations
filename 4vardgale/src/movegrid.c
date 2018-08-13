#include <stdio.h>
#include <stdlib.h>
#include "dg.h"
#include "dg1d.h"

void MoveGrid(CELL *cell, FACE *face) {
  UINT i;

  // First face
  i = 0;
  face[i].x += face[i].w * dt;
  face[i].rcell->xl = face[i].x;

  // Interior faces
  for (i = 1; i < NFMAX; i++)
    if (face[i].active && i != NF - 1) {
      face[i].x += face[i].w * dt;

      face[i].lcell->xr = face[i].x;
      face[i].rcell->xl = face[i].x;
    }

  // Last face
  i = NF - 1;
  face[i].x += face[i].w * dt;
  face[i].lcell->xr = face[i].x;

  dxmin = 1.0e20;
  dxmax = -1.0e20;
  for (i = 0; i < NCMAX; ++i)
    if (cell[i].active) {
      cell[i].h = cell[i].xr - cell[i].xl;
      cell[i].x = 0.5 * (cell[i].xl + cell[i].xr);
      GaussPoints(&cell[i]);
      dxmin = MIN(dxmin, cell[i].h);
      dxmax = MAX(dxmax, cell[i].h);
    }
}
