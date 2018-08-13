#include <math.h>
#include <stdio.h>
#include "dg.h"
#include "dg1d.h"
#include "utils.h"
// This file contains functions which allow us to perform quick operations
// on data structures used in dgale1d. This is in spirit of providing common
// operators for classes in c++.

// Copy matrix A to matrix B
void copyArray(int rows, int cols, double A[rows][cols], double B[rows][cols]) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      B[i][j] = A[i][j];
    }
  }
}

// Copy matric cell A to matrix B
void copyCellArray(int rows, int cols, double **A, double B[rows][cols]) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      B[i][j] = A[i][j];
    }
  }
}

void printArray(int rows, int cols, double A[rows][cols]) {
  for (int i = 0; i < rows; ++i) {
    printf("[");
    for (int j = 0; j < cols; ++j) {
      printf("%.16e ", A[i][j]);
    }
    printf("]\n");
  }
}

void printVector(int cols, double A[cols]) {
  printf("[");
  for (int j = 0; j < cols; ++j) {
    printf("%.16e ", A[j]);
  }
  printf("]\n");
}
// Subtract array B from array A and store the result in array C
void subtractArray(int rows,
                   int cols,
                   double A[rows][cols],
                   double B[rows][cols],
                   double C[rows][cols]) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      C[i][j] = A[i][j] - B[i][j];
    }
  }
}

// Calculate the L^\infty norm of input array (returns a double)
double supArray(int rows, int cols, double A[rows][cols]) {
  double sup = fabs(A[0][0]);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      sup = fmax(fabs(sup), fabs(A[i][j]));
    }
  }
  return sup;
}

double l1Array(int rows, int cols, double A[rows][cols]) {
  double sup = 0;
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      sup += fabs(A[i][j]);
    }
  }
  return sup;
}

// Transpose of an Array
// Transpose Array A of rows rows and cols columns into array B
void transposeArray(int rows,
                    int cols,
                    double A[rows][cols],
                    double B[cols][rows]) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      B[j][i] = A[i][j];
    }
  }
}

void test_nodal_to_modal(int nc, CELL cell[nc]) {
  // Original Values
  double original[nc][NG][NVAR];
  double nodal[nc][NG][NVAR];
  double original_cell[nc][NG][NVAR];
  for (int i = 0; i < nc; ++i) {
    copyCellArray(NVAR, NG, cell[i].U, original_cell[i]);
    transposeArray(NVAR, NG, original_cell[i], original[i]);
  }

  // Generate Nodal Values

  for (int i = 0; i < nc; ++i) {
    for (int j = 0; j < NG; ++j) {
      Uvect(&cell[i], cell[i].xg[j], nodal[i][j]);
    }
  }

  // Convert nodal to modal
  double modal[nc][NG][NVAR];
  for (int i = 0; i < nc; ++i) {
    nodal_to_modal(NG, nodal[i], modal[i]);
  }

  // Compare Arrays and Print the Result
  double error[nc][NG][NVAR];
  double scalar_error = 0;
  for (int i = 0; i < nc; ++i) {
    subtractArray(NG, NVAR, original[i], modal[i], error[i]);
    scalar_error += supArray(NG, NVAR, error[i]);
  }
  scalar_error /= nc;

  printf("Scalar Error in Conversion = %.16e \n", scalar_error);
}

void printlegVander(CELL *cell) {
  double legVander[PORD][PORD];
  for (int i = 0; i < NCMAX; ++i) {
    for (int j = 0; j < NG; ++j)
      for (int k = 0; k < NG; ++k) {
        legVander[j][k] = ShapeFun(cell[i].xg[j], &cell[i], k);
      }
    printf("LegVander Matrix for Cell %d \n", i);
    printArray(NG, NG, legVander);
  }
}

//int checkPositivity(REAL *U) {
//  REAL ZERO = 1.0e-13;
//  REAL p = (GAMMA - 1) * (U[3] - 0.5 * U[1] * U[1] / U[0]);
//  if (U[0] <= ZERO || p <= ZERO)
//    return FALSE;
//  else
//    return TRUE;
//}

int checkPositivity(REAL *U) {
  REAL V[NVAR];
  con2prim(U, V);
  REAL ZERO = 1.0e-13;
  if (V[0] <= ZERO || V[3] <= ZERO)
    return FALSE;
  else
    return TRUE;
}
