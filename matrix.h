#ifndef MATRIX_H
#define MATRIX_H
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct MatrixCDT *Matrix;

typedef struct MatrixCDT{
    double  **elem;
    int       rows;
    int       cols;
} MatrixCDT;

typedef struct CCSMatrixCDT *CCSMatrix;

typedef struct CCSMatrixCDT {
    int     *row_index;
    int     *col_ptr;
    double  *val;
    int      nnz;
    int      rows;
    int      cols;
} CCSMatrixCDT;

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

Matrix identityMatrix(int size);

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
/* Compressed sparse column matrix, a.k.a Compressed column storage. Provides
 * efficient storage for sparse matrices.*/

CCSMatrix identityCCSMatrix(int size);

void freeCCSMatrix(CCSMatrix mat);

Matrix ccsToMatrix(CCSMatrix ccs);

CCSMatrix matrixToCCS(Matrix mat);

void printCCS(CCSMatrix ccs);
    
#endif
