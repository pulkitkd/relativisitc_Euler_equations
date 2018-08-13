#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "dg.h"
#include "dg1d.h"

REAL RoeEigVal(REAL *, REAL *);
REAL MaxEigVal(REAL *, REAL *, REAL);

//-----------------------------------------------------------------------------
// compute speed of sound
//-----------------------------------------------------------------------------
REAL sound_speed(REAL *U) {
  REAL V[NVAR];
  con2prim(U, V);
  REAL rho = V[0];
  REAL p = V[3];
  REAL alpha = GAMMA / (GAMMA -1);
  REAL h = 1 + p * alpha / rho;
  return sqrt(GAMMA * p / (rho * h));
}

//-----------------------------------------------------------------------------
/* Compute flux for 1-d euler equation given conserved vector */
//-----------------------------------------------------------------------------// TODO (pulkit#1#): Modify the fluxes to use mesh velocity w



void EulerFlux(REAL *V, REAL w, REAL *flux) {
  REAL rho = V[0];
  REAL u1 = V[1];
  REAL u2 = V[2];
  REAL p = V[3];
  REAL alpha = GAMMA / (GAMMA -1);
  REAL h = 1 + p * alpha / rho;
  REAL BETAsq = 1 / (1 - u1 * u1 - u2 * u2);

  flux[0] = sqrt(BETAsq) * rho * u1; //Du1
  flux[1] = BETAsq * rho * h * u1 * u1 + p; //m1u1 + p
  flux[2] = BETAsq * rho * h * u2 * u1; //m2u1
  flux[3] = BETAsq * rho * h * u1; //(E+p)u1

//  D = U[0] = rho * BETA
//  m1 = U[1] = BETAsq * rho * h * u1
//  m2 = U[2] = BETAsq * rho * h * u2
//  E = U[3] = BETAsq * rho * h - p
}
//-----------------------------------------------------------------------------
// Compute flux Jacobian
//-----------------------------------------------------------------------------
void Jacobian(REAL *V, REAL A[][NVAR]) {
  REAL rho = V[0];
  REAL u1 = V[1];
  REAL u2 = V[2];
  REAL p = V[3];
  REAL alpha = GAMMA / (GAMMA -1);
  REAL h = 1 + p * alpha / rho;
  REAL u1sq = u1 * u1;
  REAL u2sq = u2 * u2;
  REAL u = sqrt(u1sq + u2sq);
  REAL usq = u * u;
  REAL BETA = 1 / sqrt(1 - usq);
  REAL BETAsq = BETA * BETA;
  REAL asqp = alpha * alpha * p;
  REAL H = alpha * alpha * p - rho - alpha * (-rho + p * (1 + usq));

  A[0][0] = u1;
  A[0][1] = (-1 + alpha) * rho * rho * h;
  A[0][2] = 0;
  A[0][3] = (1 - alpha) * rho * u1 / (BETAsq * H);
  A[1][0] = 0;
  A[1][1] = (asqp - rho + alpha * (-2 * p + rho)) * u1 / H;
  A[1][2] = 0;
  A[1][3] = (rho - rho * u1sq + asqp * (-1 + u1sq) +  alpha * (rho * (-1 + u1sq)
            + p * (1 - u1sq + u2sq))) / (BETAsq * H * rho * h);
  A[2][0] = 0;
  A[2][1] = -alpha * p * u2 / (BETAsq * H);
  A[2][2] = u1;
  A[2][3] = -(asqp - rho + alpha * (-2 * p + rho)) * u1 * u2
             / (BETAsq * rho * h * H);
  A[3][0] = 0;
  A[3][1] = alpha * p * rho * h / H;
  A[3][2] = 0;
  A[3][3] = (asqp - rho + alpha * (-2 * p + rho)) * u1 / H;
 }
//-----------------------------------------------------------------------------
/* Roe flux for 1-d euler equations */
//-----------------------------------------------------------------------------
void RoeFlux(REAL *Ul, REAL *Ur, REAL w, REAL *flux) {
  REAL left[NVAR], right[NVAR];
  UINT i, j;

  left[0] = Ul[0];
  left[1] = Ul[1] / Ul[0];
  left[2] = Ul[2] / Ul[0];
  left[3] = (GAMMA - 1) * (Ul[2] - 0.5 * left[0] * left[1] * left[1]);

  right[0] = Ur[0];
  right[1] = Ur[1] / Ur[0];
  right[2] = Ur[2] / Ur[0];
  right[2] = (GAMMA - 1) * (Ur[2] - 0.5 * right[0] * right[1] * right[1]);

  REAL fl = sqrt(left[0]);
  REAL fr = sqrt(right[0]);
  REAL u = (fl * left[1] + fr * right[1]) / (fl + fr);

  REAL Hl = GAMMA * left[2] / left[0] / (GAMMA - 1.0) + 0.5 * pow(left[1], 2);
  REAL Hr =
    GAMMA * right[2] / right[0] / (GAMMA - 1.0) + 0.5 * pow(right[1], 2);

  // average of fluxes
  flux[0] =
    0.5 * (left[0] * left[1] - w * Ul[0] + right[0] * right[1] - w * Ur[0]);
  flux[1] = 0.5 * (left[2] + left[0] * pow(left[1], 2) - w * Ul[1] + right[2] +
                   right[0] * pow(right[1], 2) - w * Ur[1]);
  flux[2] = 0.5 * (Hl * left[0] * left[1] - w * Ul[2] +
                   Hr * right[0] * right[1] - w * Ur[2]);

  // Add conservative dissipation
  REAL H = (fl * Hl + fr * Hr) / (fl + fr);
  REAL a = sqrt((GAMMA - 1.0) * (H - 0.5 * u * u));
  REAL R[3][3];
  R[0][0] = R[0][1] = R[0][2] = 1.0;
  R[1][0] = u - a;
  R[1][1] = u;
  R[1][2] = u + a;
  R[2][0] = H - u * a;
  R[2][1] = 0.5 * u * u;
  R[2][2] = H + u * a;

  REAL Lambda[] = {fabs(u - w - a), fabs(u - w), fabs(u - w + a)};

  REAL dU[] = {Ur[0] - Ul[0], Ur[1] - Ul[1], Ur[2] - Ul[2]};

  REAL aa[3];
  aa[1] = (GAMMA - 1.0) / (a * a) * (dU[0] * (H - u * u) + u * dU[1] - dU[3]);
  aa[0] = 0.5 / a * (dU[0] * (u + a) - dU[1] - a * aa[1]);
  aa[2] = dU[0] - aa[0] - aa[1];

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      flux[i] -= 0.5 * aa[j] * Lambda[j] * R[i][j];
}

//-----------------------------------------------------------------------------
/* ALE Roe flux for 1-d euler equations */
//-----------------------------------------------------------------------------
void ALERoeFlux(REAL *Ul, REAL *Ur, REAL w, REAL *flux) {
  REAL left[NVAR], right[NVAR];
  UINT i, j;

  left[0] = Ul[0];
  left[1] = Ul[1] / Ul[0];
  left[2] = (GAMMA - 1) * (Ul[2] - 0.5 * left[0] * left[1] * left[1]);

  right[0] = Ur[0];
  right[1] = Ur[1] / Ur[0];
  right[2] = (GAMMA - 1) * (Ur[2] - 0.5 * right[0] * right[1] * right[1]);

  REAL fl = sqrt(left[0]);
  REAL fr = sqrt(right[0]);
  REAL u = (fl * left[1] + fr * right[1]) / (fl + fr);

  REAL Hl = GAMMA * left[2] / left[0] / (GAMMA - 1.0) + 0.5 * pow(left[1], 2);
  REAL Hr =
    GAMMA * right[2] / right[0] / (GAMMA - 1.0) + 0.5 * pow(right[1], 2);

  // average of fluxes
  flux[0] =
    0.5 * (left[0] * left[1] - w * Ul[0] + right[0] * right[1] - w * Ur[0]);
  flux[1] = 0.5 * (left[2] + left[0] * pow(left[1], 2) - w * Ul[1] + right[2] +
                   right[0] * pow(right[1], 2) - w * Ur[1]);
  flux[2] = 0.5 * (Hl * left[0] * left[1] - w * Ul[2] +
                   Hr * right[0] * right[1] - w * Ur[2]);

  // Add conservative dissipation
  REAL H = (fl * Hl + fr * Hr) / (fl + fr);
  REAL a = sqrt((GAMMA - 1.0) * (H - 0.5 * u * u));
  REAL R[3][3];
  R[0][0] = R[0][1] = R[0][2] = 1.0;
  R[1][0] = u - a;
  R[1][1] = u;
  R[1][2] = u + a;
  R[2][0] = H - u * a;
  R[2][1] = 0.5 * u * u;
  R[2][2] = H + u * a;

  REAL Lambda[] = {fabs(u - w - a), fabs(u - w), fabs(u - w + a)};

  // Prevent Lambda[1] from becoming zero
  double delta = 0.1 * a;
  if (Lambda[1] < delta)
    Lambda[1] = 0.5 * (delta + Lambda[1] * Lambda[1] / delta);

  REAL dU[] = {Ur[0] - Ul[0], Ur[1] - Ul[1], Ur[2] - Ul[2]};

  REAL aa[3];
  aa[1] = (GAMMA - 1.0) / (a * a) * (dU[0] * (H - u * u) + u * dU[1] - dU[3]);
  aa[0] = 0.5 / a * (dU[0] * (u + a) - dU[1] - a * aa[1]);
  aa[2] = dU[0] - aa[0] - aa[1];

  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      flux[i] -= 0.5 * aa[j] * Lambda[j] * R[i][j];
}

//-----------------------------------------------------------------------------
/* Lax-Friedrichs flux for 1-d euler equations */
//-----------------------------------------------------------------------------
/*// TODO (pulkit#1#): Check if the Flux is implemented correctly!!!*/

void LFFlux(REAL *Vl, REAL *Vr, REAL w, REAL *flux) {
  UINT i;
  REAL Fl[NVAR], Fr[NVAR], Ul[NVAR], Ur[NVAR], lam;

  EulerFlux(Vl, w, Fl);
  EulerFlux(Vr, w, Fr);

  // lam = RoeEigVal (Ul, Ur);
  lam = MaxEigVal(Vl, Vr, w);// TODO (pulkit#1#): Check if the formula should use Conserved instead of primitives
  prim2con(Vl, Ul);
  prim2con(Vr, Ur);

  for (i = 0; i < NVAR; i++)
    flux[i] = 0.5 * (Fl[i] + Fr[i] - lam * (Ur[i] - Ul[i]));

}

//-----------------------------------------------------------------------------
/* Maximum eigenvalue of left and right states */
//-----------------------------------------------------------------------------
REAL MaxEigVal(REAL *Vl, REAL *Vr, REAL w) {
  REAL LeftEigen[NVAR], RightEigen[NVAR]; //Arrays of eigenvalues
  REAL LeftEigenMax, RightEigenMax; //Variables to store the largest left and right eigenvalues

  REAL alpha = GAMMA / (GAMMA -1);
// TODO (pulkit#1#): Modify eigenvalues to include mesh velocity w

  REAL rhol = Vl[0];
  REAL u1l = Vl[1];
  REAL u2l = Vl[2];
  REAL pl = Vl[3];
  REAL hl = 1 + pl * alpha / rhol;
  REAL cl = sqrt(GAMMA * pl / (rhol * hl));
  REAL u1lsq = u1l * u1l;
  REAL u2lsq = u2l * u2l;
  REAL ulsq = u1lsq + u2lsq;
  LeftEigen[0] = fabs(u1l);
  LeftEigen[1] = fabs(u1l);
  LeftEigen[2] = fabs((u1l * (1 - cl * cl) +
                 cl * sqrt((1 - ulsq) * (1 - u1lsq - u2lsq * cl * cl))) /
                 (1 - ulsq * cl * cl));
  LeftEigen[3] = fabs((u1l * (1 - cl * cl) -
             cl * sqrt((1 - ulsq) * (1 - u1lsq - u2lsq * cl * cl))) /
             (1 - ulsq * cl * cl));

  LeftEigenMax = fabs(LeftEigen[0] - w);
  for(int i=1; i < NVAR; i++)
    if(fabs(LeftEigen[i] - w) > LeftEigenMax)
      LeftEigenMax = fabs(LeftEigen[i] - w);

  REAL rhor = Vr[0];
  REAL u1r = Vr[1];
  REAL u2r = Vr[2];
  REAL pr = Vr[3];
  REAL hr = 1 + pr * alpha / rhor;
  REAL cr = sqrt(GAMMA * pr / (rhor * hr));
  //REAL lr = fabs(u1r - w) + cr;
  REAL u1rsq = u1r * u1r;
  REAL u2rsq = u2r * u2r;
  REAL ursq = u1rsq + u2rsq;
  RightEigen[0] = fabs(u1r);
  RightEigen[1] = fabs(u1r);
  RightEigen[2] = fabs((u1r * (1 - cr * cr) +
             cr * sqrt((1 - ursq) * (1 - u1rsq - u2rsq * cr * cr))) /
             (1 - ursq * cr * cr));
  RightEigen[3] = fabs((u1r * (1 - cr * cr) -
             cr * sqrt((1 - ursq) * (1 - u1rsq - u2rsq * cr * cr))) /
             (1 - ursq * cr * cr));
  RightEigenMax = fabs(RightEigen[0] - w);
  for(int i=1; i < NVAR; i++)
    if(fabs(RightEigen[i] - w) > RightEigenMax)
      RightEigenMax = fabs(RightEigen[i] - w);

//  dr = Ur[0];
//  u1r = Ur[1] / dr;
//  u2r = Ur[2] / dr;
//  ur = sqrt(u1r * u1r + u2r * u2r);
//  pr = (GAMMA - 1.0) * (Ur[3] - 0.5 * dr * ur * ur);
//  ar = sqrt(GAMMA * pr / dr);
//  lr = fabs(ur - w) + ar;

  if (LeftEigenMax > RightEigenMax)
    return LeftEigenMax;
  else
    return RightEigenMax;
}

//-----------------------------------------------------------------------------
/* Returns maximum eigenvalue abs(u) + a based on Roe-average */
//-----------------------------------------------------------------------------
REAL RoeEigVal(REAL *Ul, REAL *Ur) {
  REAL dl, ul, pl, hl, dr, ur, pr, hr, d1, d2, d3, u, h, a;

  dl = Ul[0];
  ul = Ul[1] / dl;
  pl = (GAMMA - 1.0) * (Ul[2] - 0.5 * dl * ul * ul);
  hl = GAMMA * pl / dl / (GAMMA - 1.0) + 0.5 * ul * ul;

  dr = Ur[0];
  ur = Ur[1] / dr;
  pr = (GAMMA - 1.0) * (Ur[2] - 0.5 * dr * ur * ur);
  hr = GAMMA * pr / dr / (GAMMA - 1.0) + 0.5 * ur * ur;

  d1 = sqrt(dl);
  d2 = sqrt(dr);
  d3 = 1.0 / (d1 + d2);

  u = (d1 * ul + d2 * ur) * d3;
  h = (d1 * hl + d2 * hr) * d3;
  a = sqrt((GAMMA - 1.0) * (h - 0.5 * u * u));

  return fabs(u) + a;
}

//-----------------------------------------------------------------------------
/* Returns maximum eigenvalue abs(u) + a based on Roe-average */
//-----------------------------------------------------------------------------
void RoeAverage(REAL *Ul, REAL *Ur, REAL *U) {
  REAL dl, ul, pl, hl, dr, ur, pr, hr, d1, d2, d3, d, u, h, p;

  dl = Ul[0];
  ul = Ul[1] / dl;
  pl = (GAMMA - 1.0) * (Ul[2] - 0.5 * dl * ul * ul);
  hl = GAMMA * pl / dl / (GAMMA - 1.0) + 0.5 * ul * ul;

  dr = Ur[0];
  ur = Ur[1] / dr;
  pr = (GAMMA - 1.0) * (Ur[2] - 0.5 * dr * ur * ur);
  hr = GAMMA * pr / dr / (GAMMA - 1.0) + 0.5 * ur * ur;

  d1 = sqrt(dl);
  d2 = sqrt(dr);
  d3 = 1.0 / (d1 + d2);

  u = (d1 * ul + d2 * ur) * d3;
  h = (d1 * hl + d2 * hr) * d3;
  d = sqrt(dl * dr);
  p = (GAMMA - 1.0) * (h - 0.5 * u * u) * d / GAMMA;

  U[0] = d;
  U[1] = d * u;
  U[3] = p / (GAMMA - 1.0) + 0.5 * d * u * u;
}

//-----------------------------------------------------------------------------
/* HLLC flux */
//-----------------------------------------------------------------------------
void HLLCFlux(REAL *Ul, REAL *Ur, REAL w, REAL *flux) {

  // We assign variables for easy reading of the flux.
  // Firstly, density
  double rho_l = Ul[0];
  double rho_r = Ur[0];

  // Next, velocity
  double v_l = Ul[1] / Ul[0];
  double v_r = Ur[1] / Ur[0];

  // We first calculate the pressure on the left and right hand side
  double p_l = (GAMMA - 1) * (Ul[2] - 0.5 * rho_l * v_l * v_l);
  double p_r = (GAMMA - 1) * (Ur[2] - 0.5 * rho_r * v_r * v_r);

  // Calculate the enthalpy
  double H_l = GAMMA / (GAMMA - 1) * p_l / rho_l + 0.5 * v_l * v_l;
  double H_r = GAMMA / (GAMMA - 1) * p_r / rho_r + 0.5 * v_r * v_r;

  // Next we calculate the speed of sound on both the side of the interface
  double c_l = sqrt(GAMMA * p_l / rho_l);
  double c_r = sqrt(GAMMA * p_r / rho_r);

  // Relative velocity of the fluid with the moving mesh
  double u_l = v_l - w;
  double u_r = v_r - w;

  // Next, we have to calculate roe average of velocity and speed of sound
  double v_roe =
    (sqrt(rho_l) * v_l + sqrt(rho_r) * v_r) / (sqrt(rho_l) + sqrt(rho_r));

  // Roe Average of Enthalpy
  double H_roe =
    (sqrt(rho_l) * H_l + sqrt(rho_r) * H_r) / (sqrt(rho_l) + sqrt(rho_r));

  // Roe Average Speed of Sound
  double c_roe = sqrt((GAMMA - 1) * (H_roe - 0.5 * v_roe * v_roe));

  //   Next, we need to calculate the signal velocities
  double s_l = fmin(u_l - c_l, v_roe - w - c_roe);
  double s_r = fmax(u_r + c_r, v_roe - w + c_roe);

  // Next, we calculate the intermedia signal velocity
  double s_m =
    (p_l - p_r + rho_r * u_r * (s_r - u_r) - rho_l * u_l * (s_l - u_l)) /
    (rho_r * (s_r - u_r) - rho_l * (s_l - u_l));

  double p_star = p_r + rho_r * (u_r - s_r) * (u_r - s_m);

  if (0 < s_l) {
    EulerFlux(Ul, w, flux);

  } else if (s_l <= 0 && 0 <= s_m) {
    double U_star[NVAR];
    double s_diff = s_l - s_m;
    U_star[0] = (s_l - u_l) * rho_l / s_diff;
    U_star[1] = ((s_l - u_l) * rho_l * v_l + p_star - p_l) / s_diff;
    U_star[2] = ((s_l - u_l) * Ul[2] - p_l * u_l + p_star * s_m) / s_diff;

    flux[0] = s_m * U_star[0];
    flux[1] = s_m * U_star[1] + p_star;
    flux[2] = s_m * U_star[2] + (s_m + w) * p_star;

  } else if (s_m <= 0 && 0 <= s_r) {
    double U_star[NVAR];
    double s_diff = s_r - s_m;
    U_star[0] = (s_r - u_r) * rho_r / s_diff;
    U_star[1] = ((s_r - u_r) * rho_r * v_r + p_star - p_r) / s_diff;
    U_star[2] = ((s_r - u_r) * Ur[2] - p_r * u_r + p_star * s_m) / s_diff;

    flux[0] = s_m * U_star[0];
    flux[1] = s_m * U_star[1] + p_star;
    flux[2] = s_m * U_star[2] + (s_m + w) * p_star;

  } else
    EulerFlux(Ur, w, flux);
}

//-----------------------------------------------------------------------------
/* HLL-CPS flux */
//-----------------------------------------------------------------------------
void HLLCPSFlux(REAL *Ul, REAL *Ur, REAL w, REAL *flux) {

   UINT i;

   // We assign variables for easy reading of the flux.
   // Firstly, density
   double rho_l = Ul[0];
   double rho_r = Ur[0];

   // Next, velocity
   double v_l = Ul[1] / Ul[0];
   double v_r = Ur[1] / Ur[0];

   // We first calculate the pressure on the left and right hand side
   double p_l = (GAMMA - 1) * (Ul[2] - 0.5 * rho_l * v_l * v_l);
   double p_r = (GAMMA - 1) * (Ur[2] - 0.5 * rho_r * v_r * v_r);

   // Calculate the enthalpy
   double H_l = GAMMA / (GAMMA - 1) * p_l / rho_l + 0.5 * v_l * v_l;
   double H_r = GAMMA / (GAMMA - 1) * p_r / rho_r + 0.5 * v_r * v_r;

   // Next we calculate the speed of sound on both the side of the interface
   double c_l = sqrt(GAMMA * p_l / rho_l);
   double c_r = sqrt(GAMMA * p_r / rho_r);
   double c = 0.5 * (c_l + c_r);

   // Relative velocity of the fluid with the moving mesh
   double u_l = v_l - w;
   double u_r = v_r - w;
   double u = 0.5 * (u_l + u_r);

   // Next, we have to calculate roe average of velocity and speed of sound
   double v_roe =
     (sqrt(rho_l) * v_l + sqrt(rho_r) * v_r) / (sqrt(rho_l) + sqrt(rho_r));

   // Roe Average of Enthalpy
   double H_roe =
     (sqrt(rho_l) * H_l + sqrt(rho_r) * H_r) / (sqrt(rho_l) + sqrt(rho_r));

   // Roe Average Speed of Sound
   double c_roe = sqrt((GAMMA - 1) * (H_roe - 0.5 * v_roe * v_roe));

   //   Next, we need to calculate the signal velocities
   double s_l = fmin(u_l - c_l, v_roe - w - c_roe);
   s_l = fmin(0.0, s_l);
   double s_r = fmax(u_r + c_r, v_roe - w + c_roe);
   s_r = fmax(0.0, s_r);
   if (fabs(u) < 1.0e-13) {
     s_l = -c;
     s_r = c;
   }

   // convective flux
   double f1[NVAR];
   if (u >= 0.0) {
     double ff = u * (u_l - s_l) / (u - s_l);
     for (i = 0; i < NVAR; ++i)
       f1[i] = ff * Ul[i];
   } else {
     double ff = u * (u_r - s_r) / (u - s_r);
     for (i = 0; i < NVAR; ++i)
       f1[i] = ff * Ur[i];
   }

   // pressure flux
   double f2[NVAR];
   f2[0] = 0.0;
   f2[1] = 0.5 * (p_l + p_r);
   f2[2] = 0.5 * (p_l * v_l + p_r * v_r);

   double a1 = 0.5 * (s_r + s_l) / (s_r - s_l);
   double a2 = s_l * s_r / (c * c * (s_r - s_l));
   double du[NVAR];

   du[0] = -a2 * (p_l - p_r);
   du[1] = a1 * (p_l - p_r) - a2 * (p_l * v_l - p_r * v_r);
   du[3] = a1 * (p_l * v_l - p_r * v_r) -
           a2 * (c * c * (p_l - p_r) / (GAMMA - 1) +
                 0.5 * (p_l * v_l * v_l - p_r * v_r * v_r));

   // total flux
   for (i = 0; i < NVAR; ++i)
     flux[i] = f1[i] + f2[i] + du[i];
}

//-----------------------------------------------------------------------------
/* HLL flux */
//-----------------------------------------------------------------------------
void HLLFlux(REAL *Ul, REAL *Ur, REAL w, REAL *flux) {
   UINT i;

   // We assign variables for easy reading of the flux.
   // Firstly, density
   double rho_l = Ul[0];
   double rho_r = Ur[0];

   // Next, velocity
   double v_l = Ul[1] / Ul[0];
   double v_r = Ur[1] / Ur[0];

   // We first calculate the pressure on the left and right hand side
   double p_l = (GAMMA - 1) * (Ul[2] - 0.5 * rho_l * v_l * v_l);
   double p_r = (GAMMA - 1) * (Ur[2] - 0.5 * rho_r * v_r * v_r);

   // Calculate the enthalpy
   double H_l = GAMMA / (GAMMA - 1) * p_l / rho_l + 0.5 * v_l * v_l;
   double H_r = GAMMA / (GAMMA - 1) * p_r / rho_r + 0.5 * v_r * v_r;

   // Next we calculate the speed of sound on both the side of the interface
   double c_l = sqrt(GAMMA * p_l / rho_l);
   double c_r = sqrt(GAMMA * p_r / rho_r);

   // Relative velocity of the fluid with the moving mesh
   double u_l = v_l - w;
   double u_r = v_r - w;

   // Next, we have to calculate roe average of velocity and speed of sound
   double v_roe =
     (sqrt(rho_l) * v_l + sqrt(rho_r) * v_r) / (sqrt(rho_l) + sqrt(rho_r));

   // Roe Average of Enthalpy
   double H_roe =
     (sqrt(rho_l) * H_l + sqrt(rho_r) * H_r) / (sqrt(rho_l) + sqrt(rho_r));

   // Roe Average Speed of Sound
   double c_roe = sqrt((GAMMA - 1) * (H_roe - 0.5 * v_roe * v_roe));

   //   Next, we need to calculate the signal velocities
   double s_l = fmin(u_l - c_l, v_roe - w - c_roe);
   double s_r = fmax(u_r + c_r, v_roe - w + c_roe);

   if (0 <= s_l)
     EulerFlux(Ul, w, flux);
   else if (s_r <= 0)
     EulerFlux(Ur, w, flux);
   else {
     double a = s_l * s_r / (s_r - s_l);
     double fl[NVAR], fr[NVAR];
     EulerFlux(Ul, w, fl);
     EulerFlux(Ur, w, fr);
     for (i = 0; i < NVAR; ++i)
       flux[i] =
         (s_r * fl[i] - s_l * fr[i]) / (s_r - s_l) + a * (Ur[i] - Ul[i]);
   }
}
