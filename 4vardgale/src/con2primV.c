#include <stdio.h>
#include <math.h>
#include "dg.h"
#include "dg1d.h"


//The function to calculate the f and f' for the function
//f(D,m,E,u). The function has coefficients in terms of m
//and not m1, m2. The output is u and not u1 or u2.
void func(double conserved[], double f[2], double u)
{
  double D, m, E;
  D = conserved[0];
  m = sqrt(conserved[1] * conserved[1] + conserved[2] * conserved[2]);
  E = conserved[3];

  double usq = u * u;
  double GE = GAMMA * E;
  double Gm = GAMMA * m;
  double c1 = GE * u - m - Gm * usq + m * usq;
  double c2 = D * (GAMMA - 1) * u;
  double c2sq = c2 * c2;
  f[0] = c1 * c1 - c2sq + c2sq * usq;
  f[1] = 2 * c1 * (GE - 2 * u * (Gm - m)) - (c2sq * (2 / u - 4 * u));


//  f[0] = (GAMMA*u*E - m - m*GAMMA*u*u + m*u*u)*
//         (GAMMA*u*E - m - m*GAMMA*u*u + m*u*u) -
//         (D*D*(GAMMA - 1)*(GAMMA - 1)*u*u*(1 - u*u));

//  f[1] = 2*(GAMMA*u*E - m - m*GAMMA*u*u + m*u*u)*
//         (GAMMA*E - 2*m*GAMMA*u + 2*m*u) -
//         (D*D*(GAMMA - 1)*(GAMMA - 1)*2*u*(1 - 2*u*u));
}

//Function to predict an initial guess for Newton-Raphson
double InitialGuess(double conserved[])
{
  double z, ul, uu, m;
  double delta = 1e-8;

  m = sqrt(conserved[1] * conserved[1] + conserved[2] * conserved[2]);
  z = 2 * m / (GAMMA * conserved[3]);
  ul = (1 - sqrt(1 - z * z * (GAMMA - 1))) / (z * (GAMMA - 1));

  if (m > 0)
    uu = fmin(1.0, m / conserved[3] - delta);
  else
    uu = fmax(-1.0, m / conserved[3] + delta);

  return (ul + uu) / 2;
}

//Newton-raphson procedure
double newton(double conserved[])
{
  double x=0, f[2], s=1, m1, m2;
  int itrtn = 0;
  m1=conserved[1];
  m2=conserved[2];
// if |m1| and |m2| both is < 1e-12, then set m = 0
  if(fabs(m1) > 1e-10 || fabs(m2) > 1e-10)
  {
    x = InitialGuess(conserved);
    func(conserved, f, x);
    while(fabs(s) > 1e-12 && fabs(f[0]) > 1e-12)
    {
      s = -f[0] / f[1];
      x = x + s;
      while(x + s > 1)
        s = s / 2;
      func(conserved,f,x);
      itrtn++;
      if(itrtn >= 30)
      {
        printf("Newton method did not converge \n");
        return 0;
      }
    }

    return x;
  }
  else
  return x;
}

void con2prim(double conserved[], double primitive[])
{
  double rho,D,m1,m2,E,m,u;

  D = conserved[0], m1 = conserved[1];
  m2 = conserved[2], E = conserved[3];
  rho = primitive[0];

  m = sqrt(m1 * m1 + m2 * m2);

  u = newton(conserved);
  if(fabs(m) < 1e-10) {
    primitive[1] = 0;
    primitive[2] = 0;
    }
  else {
    primitive[1] = m1 * u / m;
    primitive[2] = m2 * u / m;
    }
    primitive[0] = D * sqrt(1 - u * u);
    primitive[3] = (GAMMA - 1) * (E - m * u - rho);

    printf("Prim: %lf %.16lf %.16lf %lf \n" ,primitive[0], primitive[1], primitive[2], primitive[3]);
  printf("Cons: %lf %.16lf %.16lf %lf \n" ,conserved[0], conserved[1], conserved[2], conserved[3]);

 }

void prim2con(double primitive[], double conserved[])
{
  double D, rho, u1, u2, p, b, h, usq;
  int i;
  rho = primitive[0];
  u1 = primitive[1];
  u2 = primitive[2];
  p = primitive[3];

  usq = u1 * u1 + u2 * u2;
  if(usq > 1)
    {
      printf("Error! Magnitude of velocity exceeded 1."
             " The process will return Zero. \n");
      for(i=0; i<3; i++)
        conserved[i] = 0;
        rho = 0;
    }

  else {
    b = 1 / sqrt(1 - usq);

    if(fabs(rho) > 1e-10) {
      h = 1 + p * GAMMA / (rho * (GAMMA - 1));
      D = rho * b;
      conserved[0] = D; // D
      conserved[1] = D * b * h * u1; // m1
      conserved[2] = D * b * h * u2; // m2
      conserved[3] = D * b * h - p; // E
      }
    else {
      conserved[0] = 0;
      conserved[1] = 0;
      conserved[2] = 0;
      conserved[3] = 0;
      }
    }
}

//void main()
//{
//double conserved[4],primitive[4];
//
//primitive[0]=10;
//primitive[1]=0.0;
//primitive[2]=0.0;
//primitive[3]=13.33;
//
//prim2con(primitive, conserved);
//con2prim(conserved, primitive);
//printf("%lf %.16lf %.16lf %lf \n" ,primitive[0], primitive[1], primitive[2], primitive[3]);
//printf("%lf %.16lf %.16lf %lf \n" ,conserved[0], conserved[1], conserved[2], conserved[3]);
//
//}


