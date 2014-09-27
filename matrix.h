#ifndef MATRIX_H
#define MATRIX_H
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846264338327

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

typedef struct CRSMatrixCDT *CRSMatrix;

typedef struct CRSMatrixCDT {
    int     *col_index;
    int     *row_ptr;
    double  *val;
    int      nnz;
    int      rows;
    int      cols;
} CRSMatrixCDT;

/* Allocate matrix of rows x cols, and initialize it with all 0 */
Matrix nullMatrix(int rows, int cols);

/* free memory used by mat */
void freeMatrix(Matrix mat);

/* Returns a deep copy of mat */
Matrix copyMatrix(Matrix mat);

/* Returns the matrix multiplication between m1 and m2 */
Matrix matrixMult(Matrix m1, Matrix m2);

/* Returns the A matrix of the Ising model
 * Expects size of the matrix, i.e. 2*m
 */
Matrix build_A(int size);

/* Returns the L matrix of the Ising model
 * Expects size of the matrix, i.e. 2*m
 */
Matrix build_L(int size);

/* Returns the K matrix of the Ising model
 * Expects size of the matrix, i.e. 2*m
 */
Matrix build_K(int size);

/* Returns a new Matrix by transposing the elements of the input */
Matrix transposeMatrix(Matrix mat);

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

CCSMatrix build_CCS_A(int size);

CCSMatrix build_CCS_K(int size);

CCSMatrix build_CCS_L(int size);

void freeCCSMatrix(CCSMatrix mat);

/* Returns a new dense matrix with the elements of the input CCSMatrix */
Matrix ccsToMatrix(CCSMatrix ccs);

/* Returns a new CCSMatrix with the elements of the input dense matrix */
CCSMatrix matrixToCCS(Matrix mat);

/* Returns the multiplication of ccs1 with ccs2 into a new CCSMatrix */
CCSMatrix ccsMult(CCSMatrix ccs1, CCSMatrix ccs2);

/* Returns the value at position (row,col) in ccs */
double ccsValueAt(int row, int col, CCSMatrix ccs);

/* Prints memory representation, mainly for debugging */
void printCCS(CCSMatrix ccs);

void print_CCSMatrix(CCSMatrix ccs);

#endif
