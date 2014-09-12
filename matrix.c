#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double **matrixProduct(int m, double **m1, double **m2);
double **build_L(int m);
double **build_K(int m);
void printMatrix(int size, double **matrix);
void freeMatrix(int size, double **matrix);

int main(int argc, char *argv[]){
    if (argc < 2){
        printf("Please supply the m (integer) value as a "
               "parameter to the program\n");
        return EXIT_SUCCESS;
    }
    int m;
    sscanf(argv[1], "%d", &m);
    int size = 2 * m;

    double **A;
    double **K = build_K(size);
    double **L = build_L(size);
    A = matrixProduct(size, L, K);
    printMatrix(size, A);

    freeMatrix(size, K);
    freeMatrix(size, L);
    freeMatrix(size, A);

    return EXIT_SUCCESS;
}

double **
build_K(int size){
    double **K;
    nullMatrix(size, K);

    double alfa = M_PI / 4;

    double a = cos(alfa);
    double b = sin(alfa);
    double c = -b;

    size--;
    while(size > 0){
        K[size][size]     = a;
        K[size-1][size-1] = a;
        K[size][size-1]   = c;
        K[size-1][size]   = b;
        size -= 2;
    }
    return K;
}

double **
build_L(int size){
    double **L;
    nullMatrix(size, L);

    double beta = M_PI / 4;

    size--;
    L[0][0]       = cos(beta);
    L[0][size]    = -sin(beta);
    L[size][0]    = sin(beta);
    L[size][size] = cos(beta);

    double a = cos(beta);
    double b = sin(beta);
    double c = -b;

    while(size > 2){
        L[size-1][size-1] = a;
        L[size-1][size-2] = c;
        L[size-2][size-1] = b;
        L[size-2][size-2] = a;
        size -= 2;
    }
    return L;
}

/* Allocate matrix of size x size, and initialize it with all 0 */
void
nullMatrix(int size, double **matrix){
    matrix = malloc(size * sizeof(*L));
    int i, j;

    for(i = 0; i < size; i++){
        matrix[i] = malloc(size * sizeof(*matrix[i]));
    }
    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            matrix[i][j] = 0;
        }
    }
}

/* naive algorithm: O(n^3) */
double **
matrixProduct(int size, double **m1, double **m2){
    double **answ = malloc(size * sizeof(*answ));
    int i, j, k;

    for(i = 0; i < size; i++){
        answ[i] = malloc(size * sizeof(*answ[i]));
    }

    double value;

    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            value = 0;
            for(k = 0; k < size; k++){
                value += m1[i][k] * m2[k][j];
            }
            answ[i][j] = value;
        }
    }
    return answ;
}

void
printMatrix(int size, double **matrix){
    int i, j;
    char s[10];

    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++){
            sprintf(s, "%5f", matrix[i][j]);
            printf("%.8s ", s);
        }
        putchar('\n');
    }
}

void
freeMatrix(int size, double **matrix){
    int i;
    for(i = 0; i < size; i++){
        free(matrix[i]);
    }
    free(matrix);
}
