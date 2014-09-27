#include <stdarg.h>
#include "matrix.h"
#include "utils.h"

#define MIN(a, b)  ((a) > (b) ? (a) : (b))

/* Allocates a Matrix or rows x cols size, but doesn't initialize it */
Matrix newMatrix(int rows, int cols);

void   scalarProduct(double scalar, Matrix mat);
Matrix addMatrices(Matrix mat1, Matrix mat2);
Matrix substractMatrices(Matrix mat1, Matrix mat2);

static Matrix matrixFunction(double (*fn)(double,double), Matrix mat1, Matrix mat2);
static double norm2(Matrix vec);
static void   copyColumn(int col, Matrix src, double *dest);

Matrix
nullMatrix(int rows, int cols){
    Matrix ret   = xmalloc(sizeof(MatrixCDT));
    ret->elem    = xmalloc(rows * sizeof(*(ret->elem)));
    ret->elem[0] = xcalloc(sizeof(*(ret->elem[0])), rows * cols);
    ret->rows    = rows;
    ret->cols    = cols;

    for(int i = 0; i < rows; i++){
        ret->elem[i] = ret->elem[0] + cols * i;
    }
    return ret;
}

Matrix
newMatrix(int rows, int cols){
    Matrix ret   = xmalloc(sizeof(MatrixCDT));
    ret->elem    = xmalloc(rows * sizeof(*(ret->elem)));
    ret->elem[0] = xmalloc(rows * cols * sizeof(*(ret->elem[0])));
    ret->rows    = rows;
    ret->cols    = cols;

    for(int i = 0; i < rows; i++){
        ret->elem[i] = ret->elem[0] + cols * i;
    }
    return ret;
}

Matrix
build_A(int m){
    Matrix K = build_K(2*m);
    Matrix L = build_L(2*m);
    Matrix A = matrixMult(K, L);
    freeMatrix(K);
    freeMatrix(L);
    return A;
}

Matrix
build_K(int size){
    Matrix K    = nullMatrix(size, size);
    double alfa = PI / 4;
    double a    = cos(alfa);
    double b    = sin(alfa);
    double c    = -b;

    size--;
    while(size > 0){
        K->elem[size][size]     = a;
        K->elem[size-1][size-1] = a;
        K->elem[size][size-1]   = c;
        K->elem[size-1][size]   = b;
        size -= 2;
    }
    return K;
}

Matrix
build_L(int size){
    Matrix L    = nullMatrix(size, size);
    double beta = PI / 4;
    double a    = cos(beta);
    double b    = sin(beta);
    double c    = -b;

    size--;
    L->elem[0][0]       = a;
    L->elem[0][size]    = c;
    L->elem[size][0]    = b;
    L->elem[size][size] = a;

    while(size > 2){
        L->elem[size-1][size-1] = a;
        L->elem[size-1][size-2] = c;
        L->elem[size-2][size-1] = b;
        L->elem[size-2][size-2] = a;
        size -= 2;
    }
    return L;
}

/* naive algorithm: O(n^3)
 *
 * Given a matrices m1 of size n x m, and m2 of size m x p, return the
 * multiplication of both
 */
Matrix
matrixMult(Matrix mat1, Matrix mat2){
    if(mat1->cols != mat2->rows){
        die("Cannot multiply a %d x %d matrix by a %d x %d matrix",
            mat1->rows, mat1->cols, mat2->rows, mat1->cols);
    }

    Matrix ret = newMatrix(mat1->rows, mat2->cols);

    for(int i = 0; i < mat1->rows; i++){
        for(int j = 0; j < mat2->cols; j++){
            double value = 0;
            for(int k = 0; k < mat1->cols; k++){
                value += mat1->elem[i][k] * mat2->elem[k][j];
            }
            ret->elem[i][j] = value;
        }
    }

    return ret;
}

static inline double
plus(double a, double b){
    return a + b;
}

static inline double
minus(double a, double b){
    return a - b;
}

static inline Matrix
matrixFunction(double (*fn)(double,double), Matrix mat1, Matrix mat2){
    Matrix ret = newMatrix(mat1->rows, mat1->cols);

    for(int i = 0; i < mat1->rows; i++){
        for(int j = 0; j < mat1->cols; j++){
            ret->elem[i][j] = fn(mat1->elem[i][j], mat2->elem[i][j]);
        }
    }

    return ret;
}

Matrix
addMatrices(Matrix mat1, Matrix mat2){
    if(mat1->rows != mat2->rows || mat1->cols != mat2->cols){
        die("Can't add a %d x %d matrix with a %d x %d one",
             mat1->rows, mat1->cols, mat2->rows, mat2->cols);
    }
    return matrixFunction(plus, mat1, mat2);
}

Matrix
substractMatrices(Matrix mat1, Matrix mat2){
    if(mat1->rows != mat2->rows || mat1->cols != mat2->cols){
        die("Can't substract a %d x %d matrix with a %d x %d one",
             mat1->rows, mat1->cols, mat2->rows, mat2->cols);
    }
    return matrixFunction(minus, mat1, mat2);
}

void
printMatrix(Matrix mat){
    for(int i = 0; i < mat->rows; i++){
        for(int j = 0; j < mat->cols; j++){
            printf("%7.3f ", mat->elem[i][j]);
        }
        putchar('\n');
    }
}

void
freeMatrix(Matrix mat){
    free(mat->elem[0]);
    free(mat->elem);
    free(mat);
}

Matrix
copyMatrix(Matrix mat){
    Matrix ret = newMatrix(mat->rows, mat->cols);

    for(int i = 0; i < mat->rows; i++){
        for(int j = 0; j < mat->cols; j++){
            ret->elem[i][j] = mat->elem[i][j];
        }
    }
    return ret;
}

static double
norm2(Matrix vec){
    if(vec->cols != 1) die("Bad usage of norm2 function. Aborting...");

    double val = 0;
    for(int i = 0; i < vec->rows; i++){
        val += vec->elem[i][0] * vec->elem[i][0];
    }
    return sqrt(val);
}

Matrix
transposeMatrix(Matrix mat){
    Matrix ret = newMatrix(mat->cols, mat->rows);

    for(int i = 0; i < mat->cols; i++){
        for(int j = 0; j < mat->rows; j++){
            ret->elem[i][j] = mat->elem[j][i];
        }
    }
    return ret;
}

Matrix
identityMatrix(int size){
    Matrix ret = nullMatrix(size, size);
    for(int i = 0; i < size; i++){
        ret->elem[i][i] = 1.0;
    }
    return ret;
}

int
cols(Matrix mat){
    return mat->cols;
}

int
rows(Matrix mat){
    return mat->rows;
}

void
scalarProduct(double scalar, Matrix mat){
    for(int i = 0; i < mat->rows; i++){
        for(int j = 0; j < mat->cols; j++){
            mat->elem[i][j] *= scalar;
        }
    }
}

static int
sign(double x){
    return (x > 0) - (x < 0);
}

void
copyColumn(int col, Matrix src, double *dest){
    for(int row = 0; row < src->rows; row++){
        dest[row] = src->elem[row][col];
    }
}

/* TODO terminar código
 * Ejemplo de código C y explicación en:
 * http://rosettacode.org/wiki/QR_decomposition
 */
/*
void
householder(Matrix mat, Matrix Q, Matrix R){
    Matrix q[mat->rows];
    Matrix z = mat, z1;

    //double u[] = substractVector()a - norm(a)*e;

    int minDim = MIN(mat->rows-1, mat->cols);
    for(int k = 0; k < minDim; k++){
        double basisVector[mat->rows];
    }
}
*/
