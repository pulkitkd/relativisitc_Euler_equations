#include <stdio.h>
#include <math.h>
#include "dg.h"
#include "dg1d.h"

//The function to calculate the f and f' for the function
//f(D,m,E,p). The function has coefficients in terms of m
//and not m1, m2.

void func(double *conserved, double *f, double p)
{
  double D, m, E, a;
  D = conserved[0];
  m = sqrt(conserved[1]*conserved[1]+conserved[2]*conserved[2]);
  E = conserved[3];
  a = sqrt(1 - m * m / ((E + p) * (E + p)));
  f[0] = E - p / (GAMMA - 1) - D * a - m * m / (E + p);
  f[1] = -1/(GAMMA-1) - D * m * m/(pow((p+E),3)*a) + m*m/((E+p)*(E+p));
}

//Function to predict an initial guess for Newton-Raphson
double InitialGuess(double *conserved)
{
  double  D, E;
  D = conserved[0];
  E = conserved[3];
  //m = sqrt(conserved[1]*conserved[1] + conserved[2]*conserved[2]);
  return (E-D)*(GAMMA-1);
}

//Newton-Raphson procedure
double newton(double *conserved)
{
  double x=0, f[2], s=1, m1, m2, m;
  int iter=0;
  m1=conserved[1];
  m2=conserved[2];
  m=sqrt(m1 * m1 + m2 * m2);

  x = InitialGuess(conserved);
  if(fabs(m) > 1e-10)
  {
    func(conserved, f, x);
    while(fabs(s) > x * 1e-12)
    {
      s = -f[0] / f[1];

      while(x+s < 0)
        s/=2;

      x = x + s;
      func(conserved, f, x);
      iter++;
      if(iter > 50)
      {
	    printf("Error! Newton Raphson did not converge! \n");
	    printf("Conserved quantities received are:"
	             "%lf \t %e \t %lf \t %lf \n", conserved[0],
	             conserved[1], conserved[2], conserved[3]);
	    return 0;
      }
    }
    return x;
  }
  else
  return x;
}

//accepts an array of conserved variables and returns an array of primitives
void con2prim(double *conserved, double *primitive)
{
  double D,m1,m2,E,m,u;
  double p;

  D = conserved[0], m1 = conserved[1];
  m2 = conserved[2], E = conserved[3];

  m = sqrt(m1 * m1 + m2 * m2);
  p = newton(conserved);
  u = m / (E + p);

  primitive[3] = p;
  primitive[0] = D * sqrt(1 - u*u);

  if(fabs(m1) > 1e-9)
    primitive[1] = m1*u/m;
  else
    primitive[1] = 0;

  if(fabs(m2) > 1e-9)
    primitive[2] = m2*u/m;
  else
    primitive[2] = 0;

//printf("Cons: %lf %lf %lf %lf \n", conserved[0], conserved[1], conserved[2], conserved[3]);
//printf("Prim: %lf %lf %lf %lf \n", primitive[0], primitive[1], primitive[2], primitive[3]);
 }

//accepts an array of primitives and returns conserved variables
void prim2con(double *primitive, double *conserved)
{
  double D, rho, u1, u2, p, BETA, h, usq;
  int i;
  rho = primitive[0];
  u1 = primitive[1];
  u2 = primitive[2];
  p = primitive[3];

  usq=u1*u1 + u2*u2;
  if(usq>1)
    {
      printf("Error! Magnitude of velocity exceeded 1."
             " The process will return Zero. \n");
      for(i=0; i<3; i++)
        conserved[i]=0;
      primitive[0]=0;
    }


  BETA = 1/sqrt(1 - usq);
  h = 1 + p*GAMMA/(rho*(GAMMA - 1));
  D = rho * BETA;

  conserved[0] = D; // D
  conserved[1] = D * BETA * h * u1; // m1
  conserved[2] = D * BETA * h * u2; // m2
  conserved[3] = D * h * BETA - p; // E

}

//void main()
//{
//double conserved[4],primitive[4];
//
//primitive[0] = 10.0;
//primitive[1] = 0.0;
//primitive[2] = 0.0;
//primitive[3] = 13.33;
//
//prim2con(primitive, conserved);
//printf("D= %lf m1= %lf m2= %lf E= %lf \n", conserved[0], conserved[1],conserved[2],conserved[3]);
//
//con2prim(conserved, primitive);
//printf("rho=%lf u1=%e u2=%e p=%lf \n",primitive[0],primitive[1],primitive[2], primitive[3]);
//}


