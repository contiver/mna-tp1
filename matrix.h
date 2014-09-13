#ifndef MATRIX_H
#define MATRIX_H

typedef struct Matrix *Matrix;

Matrix matrixMult(Matrix m1, Matrix m2);

Matrix build_L(int size);

Matrix build_K(int size);

/* Allocate matrix of size x size, and initialize it with all 0 */
Matrix nullMatrix(int rows, int cols);

void printMatrix(Matrix mat);
void printEigenvector(Matrix mat);
void freeMatrix(Matrix mat);
double norm2(Matrix vec);
Matrix powerIteration(Matrix A);
void* mallocTest(size_t size);
void* callocTest(size_t nmemb, size_t size);

#endif
