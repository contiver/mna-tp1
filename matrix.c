#include <stdarg.h>
#include "matrix.h"


/* Allocates a Matrix or rows x cols size, but doesn't initialize it */
Matrix newMatrix(int rows, int cols);

static double norm2(Matrix vec);

void * xmalloc(size_t len);
void * xcalloc(size_t nmemb, size_t size);
void * xrealloc(void *p, size_t len);
void die(const char *errstr, ...);
Matrix matrixAddition(Matrix mat1, Matrix mat2);
void scalarProduct(double scalar, Matrix mat);

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

    for(int i = 0; i < rows; i++){
        ret->elem[i] = ret->elem[0] + cols * i;
    }

    ret->rows = rows;
    ret->cols = cols;

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
    double a = cos(beta);
    double b = sin(beta);
    double c = -b;

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
    free(mat->elem[0]);
    free(mat->elem);
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

    /* FIXME cableadas 100 iteraciones, modificar */
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

    for(int i = 0; i < mat->rows; i++){
        for(int j = 0; j < mat->cols; j++){
            ret->elem[i][j] = mat->elem[i][j];
        }
    }

    return ret;
}

static double
norm2(Matrix vec){
    if(vec->cols != 1){
        die("Bad usage of norm2 function. Aborting...");
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

void *
xrealloc(void *p, size_t len) {
    if( !(p = realloc(p, len)) )
        die("Out of memory\n");

    return p;
}

void *
xmalloc(size_t len) {
    void *p = malloc(len);

    if(!p) die("Out of memory.\n");

    return p;
}

void *
xcalloc(size_t nmemb, size_t size) {
    void *p = calloc(nmemb, size);

    if(!p) die("Out of memory.\n");

    return p;
}

void
die(const char *errstr, ...) {
    va_list ap;

    va_start(ap, errstr);
    vfprintf(stderr, errstr, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
}
