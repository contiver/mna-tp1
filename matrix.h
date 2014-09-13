#ifndef MATRIX_H
#define MATRIX_H

typedef struct Matrix *Matrix;

/* Returns matrix multiplication between m1 and m2 */
Matrix matrixMult(Matrix m1, Matrix m2);

/* Returns the L matrix of the Ising model
 * Expects size of the matrix, i.e. 2*m
 */
Matrix build_L(int size);

/* Returns the K matrix of the Ising model
 * Expects size of the matrix, i.e. 2*m
 */
Matrix build_K(int size);

/* Allocate matrix of rows x cols, and initialize it with all 0 */
Matrix nullMatrix(int rows, int cols);

/* Allocates a Matrix, but doesn't initialize it */
Matrix newMatrix(void);

/* free memory used by Matrix ADT */
void freeMatrix(Matrix mat);

void printMatrix(Matrix mat);

void printEigenvector(Matrix mat);


double norm2(Matrix vec);

Matrix powerIteration(Matrix A);

void* mallocTest(size_t size);

void* callocTest(size_t nmemb, size_t size);

#endif
