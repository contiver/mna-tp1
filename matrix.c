#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

struct Matrix{
    double **matrix;
    int rows;
    int cols;
};

/* local functions */
static double norm2(Matrix vec);
static void* mallocTest(size_t size);
static void* callocTest(size_t nmemb, size_t size);


Matrix
newMatrix(){
    Matrix ret = mallocTest(sizeof(Matrix));
    ret->matrix = NULL;
    ret->rows = 0;
    ret->cols = 0;
    return ret;
}

Matrix
build_K(int size){
    Matrix K = nullMatrix(size, size);
    double **mat = K->matrix;

    double alfa = M_PI / 4;

    double a = cos(alfa);
    double b = sin(alfa);
    double c = -b;

    size--;
    while(size > 0){
        mat[size][size]     = a;
        mat[size-1][size-1] = a;
        mat[size][size-1]   = c;
        mat[size-1][size]   = b;
        size -= 2;
    }
    return K;
}

Matrix
build_L(int size){
    Matrix L = nullMatrix(size, size);
    double **mat = L->matrix;

    double beta = M_PI / 4;

    size--;
    mat[0][0]       = cos(beta);
    mat[0][size]    = -sin(beta);
    mat[size][0]    = sin(beta);
    mat[size][size] = cos(beta);

    double a = cos(beta);
    double b = sin(beta);
    double c = -b;

    while(size > 2){
        mat[size-1][size-1] = a;
        mat[size-1][size-2] = c;
        mat[size-2][size-1] = b;
        mat[size-2][size-2] = a;
        size -= 2;
    }
    return L;
}

Matrix
nullMatrix(int rows, int cols){
    double **matrix = callocTest(rows, sizeof(*matrix));
    int i;

    for(i = 0; i < rows; i++){
        matrix[i] = callocTest(cols, sizeof(*matrix[i]));
    }

    Matrix answ = newMatrix();
    answ->rows = rows;
    answ->cols = cols;
    answ->matrix = matrix;

    return answ;
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
    int n = mat1->rows;
    int m = mat1->cols;
    int p = mat2->cols;
    double **m1 = mat1->matrix;
    double **m2 = mat2->matrix;
    double **answ = mallocTest(n * sizeof(*answ));
    double value;
    int i, j, k;

    for(i = 0; i < n; i++){
        answ[i] = mallocTest(p * sizeof(*answ[i]));
    }

    for(i = 0; i < n; i++){
        for(j = 0; j < p; j++){
            value = 0;
            for(k = 0; k < m; k++){
                value += m1[i][k] * m2[k][j];
            }
            answ[i][j] = value;
        }
    }

    Matrix ret = newMatrix();
    ret->matrix = answ;
    ret->rows = n;
    ret->cols = p;
    return ret;
}

void
printMatrix(Matrix mat){
    int i, j;
    char s[10];

    for(i = 0; i < mat->rows; i++){
        for(j = 0; j < mat->cols; j++){
            sprintf(s, "%5f", mat->matrix[i][j]);
            printf("%.8s ", s);
        }
        putchar('\n');
    }
}

void
printEigenvector(Matrix mat){
    printf("Eigenvector v1:\n");
    printf("[%f", mat->matrix[0][0]);
    int i;
    for(i = 1; i < mat->rows; i++){
        printf(" %f", mat->matrix[i][0]);
    }
    printf("]\n");
}

void
freeMatrix(Matrix mat){
    if(mat->matrix != NULL){
        int i;
        for(i = 0; i < mat->rows; i++){
            free(mat->matrix[i]);
        }
        free(mat->matrix);
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
    int i, j;

    for(i = 0; i < A->rows; i++){
        p->matrix[i] = mallocTest(sizeof(p->matrix[i]));
        p->matrix[i][0] = 1.0;
    }

    /* cableadas 100 iteraciones, modificar */
    for(i = 0; i < 200; i++){
        p = matrixMult(A, p);
        norm = norm2(p);
        for(j = 0; j < A->rows; j++){
            p->matrix[j][0] /= norm;
        }
    }
    return p;
}

Matrix
copyMatrix(Matrix mat){
    Matrix ret = newMatrix();
    ret->cols = mat->cols;
    ret->rows = mat->rows;
    ret->matrix = malloc(mat->rows * sizeof(*(ret->matrix)));

    int i, j;
    for(i = 0; i < mat->rows; i++){
        ret->matrix[i] = malloc(mat->cols * sizeof(*(ret->matrix[i])));
    }

    for(i = 0; i < mat->rows; i++){
        for(j = 0; j < mat->cols; j++){
            ret->matrix[i][j] = mat->matrix[i][j];
        }
    }

    return ret;
}

static double
norm2(Matrix vec){
    if(vec->cols != 1){
        perror("Bad usage of norm2 function. Aborting...");
        exit(EXIT_FAILURE);
    }
    double val = 0;
    int i;

    for(i = 0; i < vec->rows; i++){
        val += vec->matrix[i][0] * vec->matrix[i][0];
    }

    return sqrt(val);
}

/* Wrapper for checking if malloc returned NULL */
static void*
mallocTest(size_t size){
    void* mem = malloc(size);
    if(mem == NULL){
        perror("Couldn't allocate memory upon malloc call. Aborting...");
        exit(EXIT_FAILURE);
    }
    return mem;
}

/* Wrapper for checking if malloc returned NULL */
static void*
callocTest(size_t nmemb, size_t size){
    void* mem = calloc(nmemb, size);
    if(mem == NULL){
        perror("Couldn't allocate memory upon calloc call. Aborting...");
        exit(EXIT_FAILURE);
    }
    return mem;
}

Matrix
transpose(Matrix M){
    double **answ = mallocTest(M->cols * sizeof(*answ));
    int i, j;

    for(i = 0; i < M->cols; i++){
        answ[i] = mallocTest(M->rows * sizeof(*answ[i]));
    }

    for(i = 0; i < M->cols; i++){
        for(j = 0; j < M->rows; j++){
            answ[i][j] = M->matrix[j][i];
        }
    }

    Matrix ret = newMatrix();
    ret->matrix = answ;
    ret->rows = M->cols;
    ret->cols = M->rows;
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
