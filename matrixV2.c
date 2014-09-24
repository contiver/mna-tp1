#include <stdlib.h>

typedef struct MatrixCDT *Matrix;

typedef struct MatrixCDT{
    int       rows;
    int       cols;
    double   *elem;
} MatrixCDT;

inline double dot(double *a, double *b, int len, int step);

Matrix
new_Matrix(int rows, int cols){
    Matrix ret = malloc(sizeof(MatrixCDT));
    ret->elem  = malloc(sizeof(ret->elem[0]) * rows * cols);
    ret->rows  = rows;
    ret->cols  = cols;
    return ret;
}

Matrix
null_Matrix(int rows, int cols){
    Matrix ret = malloc(sizeof(MatrixCDT));
    ret->elem  = calloc(rows * cols, sizeof(ret->elem[0]));
    ret->rows  = rows;
    ret->cols  = cols;
    return ret;
}

Matrix
identity(int size){
    Matrix ret = null_Matrix(size, size);
    for(int i = 0; i < size; i++){
        ret->elem[i * size + i] = 1;
    }
    return ret;
}

Matrix
matrix_Mult(Matrix mat1, Matrix mat2){
    if (mat1->cols != mat2->rows) return NULL;

    Matrix ret = new_Matrix(mat1->rows, mat2->cols);
    double *currElem = ret->elem;
    double *mat1Row  = mat1->elem;

    for (int row = 0; row < mat1->rows; row++, mat1Row += mat1->cols){
        for (int col = 0; col < mat2->cols; col++){
            *currElem++ = dot(mat1Row, mat2->elem + col, mat1->cols, mat2->cols);
        }
    }

    return ret;
}

inline double
dot(double *row1, double *col2, int len, int step){
    double ret = 0;
    while (len--) {
        ret  += *row1++ * *col2;
        col2 += step;
    }
    return ret;
}

double
elem_At(Matrix mat, int row, int col){
    return mat->elem[row * mat->cols + col];
}
