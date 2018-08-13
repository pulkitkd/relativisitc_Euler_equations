#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NVAR 4
#define GAMMA 1.4
void MultiplyMatrix(double R1[][NVAR], double R2[][NVAR], double A[NVAR][NVAR]) {
  int i, j, k;

  for(k = 0; k < NVAR; k++)
    for (i = 0; i < NVAR; i++)
    {
      A[i][k] = 0.0;
      for (j = 0; j < NVAR; j++)
        A[i][k] += R1[i][j] * R2[j][k];
    }
}

void main()
{
double R1[NVAR][NVAR];
double R2[NVAR][NVAR];
double UV[NVAR][NVAR];
double UVi[NVAR][NVAR];
double temp[NVAR][NVAR];
double R[NVAR][NVAR];
double Ri[NVAR][NVAR];

int i, j, k;
double r, u1, u2, p;
    r = 2; //defined as rho usually
    u1 = 0.5;
    u2 = 0;
    p = 1.5;
    double a = GAMMA / (GAMMA - 1); // defined as alpha usually
    double h = 1 + p * a / r;
    double u1sq = u1 * u1;
    double u2sq = u2 * u2;
    double usq = u1sq + u2sq;
    double BETA = 1 / sqrt(1 - usq);
    double BETAsq = BETA * BETA;
    double H = a * a * p - r - a * (-r + p * (1 + usq));

/* R1 stores primitive right eigenvectors */
    R1[0][0] = 0;
    R1[1][0] = 1;
    R1[2][0] = r / ( p * GAMMA);
    R1[3][0] = r / ( p * GAMMA);

    R1[0][1] = 0;
    R1[1][1] = 0;
    R1[2][1] = sqrt(-a * p * ((-1 + a) * h * (-1 + u1sq) + a * p * u2sq)) / (a * p * h * BETA);
    R1[3][1] = -sqrt(-a * p * ((-1 + a) * h * (-1 + u1sq) + a * p * u2sq)) / (a * p * h * BETA);

    R1[0][2] = 1;
    R1[1][2] = 0;
    R1[2][2] = -a * p * h * (-1 + u1sq) * u2 - a * p * h * u2sq * u2 + u1 * u2 *
		sqrt(-a * p * ((-1 + a) * h * (-1 + u1sq) + a * p * u2sq)) / (a * p * h * BETA *
		(-1 + u1sq));

    R1[3][2] =  -a * p * h * (-1 + u1sq) * u2 - a * p * h * u2sq * u2 - u1 * u2 *
		sqrt(-a * p * ((-1 + a) * h * (-1 + u1sq) + a * p * u2sq)) / (a * p * h * BETA *
		(-1 + u1sq));

    R1[0][3] = 0;
    R1[1][3] = 0; 
    R1[2][3] = 1;
    R1[3][3] = 1;

    /* R2 stores inverse of R1 */
    /* Mathematically, it should represent the Matrix of Left Eigenvectors*/
    R2[0][0] = -GAMMA * p / r;
    R2[0][1] = -a * p * u2 / ((a - 1) * BETAsq * r * h * (-1 + u1sq));
    R2[0][2] = 0;
    R2[0][3] = 0;

    R2[1][0] = 0;
    R2[1][1] = u1 * u2 / (1 - u1sq);
    R2[1][2] = (a * p * h * BETA) / sqrt(-a * p * (r - r * u1 * u1 + a * a * p * (-1 + u1sq) +
                 a * (r * (-1 + u1sq) + p * (1 - u1sq + u2sq))));
    R2[1][3] = -(a * p * h * BETA) / sqrt(-a * p * (r - r * u1 * u1 + a * a * p * (-1 + u1sq) +
                 a * (r * (-1 + u1sq) + p * (1 - u1sq + u2sq))));

    R2[2][0] = 0;
    R2[2][1] = 1;
    R2[2][2] = 0;
    R2[2][3] = 0;

    R2[3][0] = 1;
    R2[3][1] = 0;
    R2[3][2] = 1;
    R2[3][3] = 1;
    
    /*We now multiply them with du/dv to get the conserved eigenvectors*/

//   MultiplyMatrix(UV, R1, R);
//   MultiplyMatrix(R2, UVi, Ri);

    MultiplyMatrix(R2, R1, temp);
    for(i = 0; i < NVAR; i++){
      for(j = 0; j < NVAR; j++){
        printf("%lf \t", temp[i][j]);
        }
        printf("\n");
     }
     printf("\n");
}
