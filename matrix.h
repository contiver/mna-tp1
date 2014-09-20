#ifndef MATRIX_H
#define MATRIX_H

typedef struct Matrix *Matrix;

/* Allocate matrix of rows x cols, and initialize it with all 0 */
Matrix nullMatrix(int rows, int cols);

/* free memory used by Matrix ADT */
void freeMatrix(Matrix mat);

/* Returns a deep copy of mat */
Matrix copyMatrix(Matrix mat);

/* Returns the matrix multiplication between m1 and m2 */
Matrix matrixMult(Matrix m1, Matrix m2);

/* Returns the L matrix of the Ising model
 * Expects size of the matrix, i.e. 2*m
 */
Matrix build_L(int size);

/* Returns the K matrix of the Ising model
 * Expects size of the matrix, i.e. 2*m
 */
Matrix build_K(int size);

/* Returns a new Matrix by transposing the elements of the input */
Matrix transpose(Matrix mat);

/* Returns a Matrix (vector) containing the eigenvector for the
 * dominant eigenvalue */
Matrix powerIteration(Matrix A);

void printMatrix(Matrix mat);

void printEigenvector(Matrix mat);

/* Return column size of mat */
int cols(Matrix mat);

/* Return row size of mat */
int rows(Matrix mat);

/* ========================================================================= */

typedef struct CCSMatrix *CCSMatrix;

CCSMatrix identityCCSMatrix(int size);

void freeCCSMatrix(CCSMatrix mat);

Matrix ccsToMatrix(CCSMatrix ccs);

void printCCS(CCSMatrix ccs);
    
#endif
