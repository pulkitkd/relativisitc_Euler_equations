#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dg.h"
#include "dg1d.h"

void MeshVel(FACE *face) {
  UINT i, it;
  REAL UL[NVAR], UR[NVAR], vl, vr;

  if (ALE == 0)
    return;

  // First face
  i = 0;
  if (bc_left == FIXED) {
    face[i].w = 0.0;
  } else {
    Uvect(face[i].rcell, face[i].x, UR);
    face[i].w = UR[1] / UR[0];
  }
  face[i].rcell->wl = face[i].w;

  // Last face
  i = NF - 1;
  if (bc_right == FIXED) {
    face[i].w = 0.0;
  } else {
    Uvect(face[i].lcell, face[i].x, UL);
    face[i].w = UL[1] / UL[0];
  }
  face[i].lcell->wr = face[i].w;

  // Interior faces
  for (i = 1; i < NFMAX; i++)
    if (face[i].active && i != NF - 1) {
      Uvect(face[i].lcell, face[i].x, UL);
      Uvect(face[i].rcell, face[i].x, UR);

      vl = UL[1] / UL[0];
      vr = UR[1] / UR[0];
      face[i].w = 0.5 * (vl + vr);
      //         RoeAverage(UL, UR, UA);
      //         face[i].w = UA[1]/UA[0];

      face[i].lcell->wr = face[i].w;
      face[i].rcell->wl = face[i].w;
    }

  // apply springels velocity correction
  for (i = 1; i < NFMAX; i++)
    if (face[i].active && i != NF - 1) {
      REAL R = face[i].rcell->x - face[i].lcell->x;
      REAL s = 0.5 * (face[i].rcell->x + face[i].lcell->x);
      REAL r = face[i].x;
      REAL d = fabs(s - r);

      Uvect(face[i].lcell, face[i].x, UL);
      Uvect(face[i].rcell, face[i].x, UR);
      REAL cl = sound_speed(UL);
      REAL cr = sound_speed(UR);
      REAL c = 0.5 * (cl + cr);

      REAL rat = d / (eta * R);
      if (0.9 <= rat && rat <= 1.1)
        face[i].w +=
          chi * c * (s - r) / d * (d - 0.9 * eta * R) / (0.2 * eta * R);
      else if (rat > 1.1)
        face[i].w += chi * c * (s - r) / d;
    }

  // smooth mesh velocity
  for (it = 0; it < siter; ++it) {
    for (i = 1; i < NFMAX; i++)
      if (face[i].active && i != NF - 1) {
        face[i].w = (face[i].w + face[i].lcell->wl + face[i].rcell->wr) / 3.0;
      }

    // copy to cells
    for (i = 1; i < NFMAX; i++)
      if (face[i].active && i != NF - 1) {
        face[i].lcell->wr = face[i].w;
        face[i].rcell->wl = face[i].w;
      }
  }
}

// Do linear interpolation to compute mesh velocity at x
double GetMeshVel(CELL *cell, REAL x) {
  return ((cell->xr - x) * cell->wl + (x - cell->xl) * cell->wr) / cell->h;
}
