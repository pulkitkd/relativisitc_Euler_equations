#ifndef UTILS_H
#define UTILS_H

// Copy Arrays
void copyArray(int rows, int cols, double A[rows][cols], double B[rows][cols]);

// Special Function to Copy Cell Arrays because of the double pointer nature
void copyCellArray(int rows, int cols, double **A, double B[rows][cols]);

// Print Array
void printArray(int rows, int cols, double A[rows][cols]);
void printVector(int cols, double A[cols]);

// Subtract array B from array A and store the result in array C
void subtractArray(int rows,
                   int cols,
                   double A[rows][cols],
                   double B[rows][cols],
                   double C[rows][cols]);

// Calculate the L^\infty norm of input array (returns a double)
double supArray(int rows, int cols, double A[rows][cols]);
double l1Array(int rows, int cols, double A[rows][cols]);

// Transpose of an Array
// Transpose Array A of rows rows and cols columns into array B
void transposeArray(int rows,
                    int cols,
                    double A[rows][cols],
                    double B[cols][rows]);

// Provides real time testing to above function
void test_nodal_to_modal(int nc, CELL cell[nc]);
void printlegVander(CELL *cell);
int checkPositivity(REAL *);
#endif // UTILS_H
