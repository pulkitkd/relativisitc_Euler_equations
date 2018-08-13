#include <stdio.h>
#include <math.h>
#define PI 3.1415926535897932

void main()
{
double rho, u, h;
FILE *fp;
double Tf, dt, x;
int NC, i;
x = 0;
NC = 1000;
h = 1.0 / NC;
u = 0.9;
Tf = 1;
dt = 0.005;
fp = fopen("sin-density-exact.dat","w");

//for(t = 0; t < Tf; t = t + dt) 
    while(x < 1) {
      rho = 1 + 0.5 * sin(2 * PI * (x - u * Tf));
      printf("%lf \t %lf \n",x, rho);
      fprintf(fp, "%lf \t %lf \n",x,rho);
      x = x + h;
    }
    fclose(fp);
}

