#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define NVAR 2
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

void main() {

double A[2][2], I[2][2], temp[2][2];
int i, j = 0;

A[0][0] = 1;
A[0][1] = 6;
A[1][0] = 3;
A[1][1] = 4;

I[0][0] = 1;
I[0][1] = 3;
I[1][0] = 5;
I[1][1] = 1;

MultiplyMatrix(A, I, temp);

 for(i = 0; i < NVAR; i++){
      for(j = 0; j < NVAR; j++){
        printf("%lf \t", temp[i][j]);
        }
        printf("\n");
     }
     printf("\n");
}


