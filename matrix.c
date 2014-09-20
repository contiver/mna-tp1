#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "matrix.h"

#define PI 3.14159265358979323846264338327

typedef struct MatrixCDT{
    double  **elem;
    int       rows;
    int       cols;
} MatrixCDT;

static double norm2(Matrix vec);

/* Allocates a Matrix or rows x cols size, but doesn't initialize it */
static Matrix newMatrix(int rows, int cols);


Matrix
nullMatrix(int rows, int cols){
    Matrix ret   = malloc(sizeof(MatrixCDT));
    ret->elem    = malloc(rows * sizeof(*(ret->elem)));
    ret->elem[0] = calloc(sizeof(*(ret->elem[0])), rows * cols);

    for(int i = 0; i < rows; i++){
        ret->elem[i] = ret->elem[0] + cols * i;
    }

    ret->rows = rows;
    ret->cols = cols;

    return ret;
}

static Matrix
newMatrix(int rows, int cols){
    Matrix ret   = malloc(sizeof(MatrixCDT));
    ret->elem    = malloc(rows * sizeof(*(ret->elem)));
    ret->elem[0] = malloc(rows * cols * sizeof(*(ret->elem[0])));

    for(int i = 0; i < rows; i++){
        ret->elem[i] = ret->elem[0] + cols * i;
    }

    ret->rows = rows;
    ret->cols = cols;

    return ret;
}

Matrix
build_K(int size){
    Matrix K = nullMatrix(size, size);

    double alfa = PI / 4;

    double a = cos(alfa);
    double b = sin(alfa);
    double c = -b;

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
    Matrix L = nullMatrix(size, size);

    double beta = PI / 4;

    size--;
    L->elem[0][0]       = cos(beta);
    L->elem[0][size]    = -sin(beta);
    L->elem[size][0]    = sin(beta);
    L->elem[size][size] = cos(beta);

    double a = cos(beta);
    double b = sin(beta);
    double c = -b;

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
        char s[70];
        sprintf(s,"Cannot multiply %d x %d matrix by a %d x %d matrix",
                mat1->rows, mat1->cols, mat2->rows, mat1->cols);
        perror(s);
        exit(EXIT_FAILURE);
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

Matrix
matrixAddition(Matrix mat1, Matrix mat2){
    if(mat1->rows != mat2->rows || mat1->cols != mat2->cols){
        perror("Can't add matrices of different size");
        return NULL;
    }
    Matrix ret = newMatrix(mat1->rows, mat1->cols);

    for(int i = 0; i < mat1->rows; i++){
        for(int j = 0; j < mat1->cols; j++){
            ret->elem[i][j] = mat1->elem[i][j] + mat2->elem[i][j];
        }
    }

    return ret;
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
printEigenvector(Matrix mat){
    printf("Eigenvector v1:\n");
    printf("[%f", mat->elem[0][0]);

    for(int i = 1; i < mat->rows; i++){
        printf(" %f", mat->elem[i][0]);
    }
    printf("]\n");
}

void
freeMatrix(Matrix mat){
    if(mat->elem != NULL){
        free(mat->elem[0]);
        free(mat->elem);
    }
    free(mat);
}

/* Requires that:
 * A has an eigenvalue that is strictly greater in magnitude than its other
 * eigenvalues
 * The starting vector p has a nonzero component in the direction of an
 * eigenvector associated with the dominant eigenvalue.
 */
Matrix
powerIteration(Matrix A){
    Matrix p = nullMatrix(A->rows, 1);
    double norm;

    for(int i = 0; i < A->rows; i++){
        p->elem[i] = malloc(sizeof(p->elem[i]));
        p->elem[i][0] = 1.0;
    }

    /* cableadas 100 iteraciones, modificar */
    for(int i = 0; i < 100; i++){
        p = matrixMult(A, p);
        norm = norm2(p);
        for(int j = 0; j < A->rows; j++){
            p->elem[j][0] /= norm;
        }
    }
    return p;
}

Matrix
copyMatrix(Matrix mat){
    Matrix ret = newMatrix(mat->rows, mat->cols);

    for(int i = 0; i < mat->rows; i++)
        for(int j = 0; j < mat->cols; j++)
            ret->elem[i][j] = mat->elem[i][j];

    return ret;
}

static double
norm2(Matrix vec){
    if(vec->cols != 1){
        perror("Bad usage of norm2 function. Aborting...");
        exit(EXIT_FAILURE);
    }
    double val = 0;

    for(int i = 0; i < vec->rows; i++){
        val += vec->elem[i][0] * vec->elem[i][0];
    }

    return sqrt(val);
}

Matrix
transpose(Matrix mat){
    Matrix ret = newMatrix(mat->cols, mat->rows);

    for(int i = 0; i < mat->cols; i++)
        for(int j = 0; j < mat->rows; j++)
            ret->elem[i][j] = mat->elem[j][i];

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
cols(Matrix M){
    return M->cols;
}

int
rows(Matrix M){
    return M->rows;
}

/* ========================================================================= */
/* Compressed sparse column matrix, a.k.a
 * Compressed column storage */

typedef struct CCSMatrixCDT {
    int     *row_index;
    int     *col_ptr;
    double  *val;
    int      nnz;
    int      rows;
    int      cols;
} CCSMatrixCDT;

CCSMatrix
identityCCSMatrix(int size){
    CCSMatrix ret  = malloc(sizeof(CCSMatrixCDT));
    ret->val       = malloc(size * sizeof (ret->val[0]));
    ret->row_index = malloc(size * sizeof (ret->row_index[0]));
    ret->col_ptr   = malloc((size+1) * sizeof (ret->col_ptr[0]));
    ret->nnz       = size;
    ret->rows      = size;
    ret->cols      = size;

    for(int i = 0; i < size; i++){
        ret->row_index[i] = i;
        ret->col_ptr[i]   = i;
        ret->val[i]       = 1.0;
    }
    ret->col_ptr[size] = size;

    return ret;
}

void
freeCCSMatrix(CCSMatrix mat){
    free(mat->val);
    free(mat->row_index);
    free(mat->col_ptr);
    free(mat);
}

void
printCCS(CCSMatrix ccs){
    printf("val vector:\n");
    for(int i = 0; i < ccs->nnz; i++){
        printf("%.2f ", ccs->val[i]);
    }
    printf("\n");

    printf("row_index vector:\n");
    for(int i = 0; i < ccs->nnz; i++){
        printf("%d ", ccs->row_index[i]);
    }
    printf("\n");

    printf("col_ptr vector:\n");
    for(int i = 0; i <= ccs->nnz; i++){
        printf("%d ", ccs->col_ptr[i]);
    }
    printf("\n");

    printf("ccs->nnz: %d\n", ccs->nnz);
    printf("ccs->rows: %d\n", ccs->rows);
    printf("ccs->cols: %d\n", ccs->cols);
}

Matrix
ccsToMatrix(CCSMatrix ccs){
    Matrix ret = nullMatrix(ccs->rows, ccs->rows);

    int j = 0;
    int elements = 0;

    for(int col = 0; col < ccs->nnz; col++){
        elements += ccs->col_ptr[col+1] - ccs->col_ptr[col];

        while(j < elements){
            ret->elem[ccs->row_index[j]][col] = ccs->val[j];
            j++;
        }
    }
    return ret;
}

CCSMatrix
matrixToCCS(Matrix mat){
    CCSMatrix ret  = malloc(sizeof(CCSMatrixCDT));
    int memSize    = (mat->cols * mat->rows) * 0.1 + 2;
    ret->val       = malloc(memSize * sizeof(ret->val[0]));
    ret->row_index = malloc(memSize * sizeof(ret->row_index[0]));
    ret->col_ptr   = malloc((mat->cols + 1) * sizeof(ret->col_ptr[0]));
    ret->cols      = mat->cols;
    ret->rows      = mat->rows;

    for(int col = 0; col < mat->cols; col++){
        bool columnStarterNotFound = true;
        for(int row = 0; row < mat->rows; row++){
            if(mat->elem[row][col] != 0.0){
                if(ret->nnz == memSize){
                    memSize *= 1.2;     // Arbitrary incremental value
                    ret->val = realloc(ret->val, memSize * sizeof(ret->val[0]));
                    ret->row_index = realloc(ret->row_index, memSize * sizeof(ret->row_index[0]));
                }
                ret->val[ret->nnz] = mat->elem[row][col];
                ret->row_index[ret->nnz] = row;
                if(columnStarterNotFound){
                    ret->col_ptr[col] = ret->nnz;
                }
                columnStarterNotFound = false;
                ret->nnz++;
            }
        }
    }
    ret->col_ptr[ret->cols] = ret->nnz;

    return ret;
}

/*
   void
   print_CCSMatrix(CCSMatrix ccs){
   int val_index = 0;

   int i = 0, j = 0;
   while(i < ccs->rows){
   while(j < ccs->cols){
   if(ccs->row_index[ccs->col_ptr[j]] == i)
   if(ccs->row_index[val_index] == i && ccs->col_ptr[j] == ){
   printf("");
   }else{
   printf("0 ");
   }

   }
   }
   if( i == ccs->row_index[i] )

   }
   }
   }
   */
