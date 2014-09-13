#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

int
main(int argc, char *argv[]){
    if (argc < 2){
        printf("Please supply the m (integer) value as a "
               "parameter to the program\n");
        return EXIT_SUCCESS;
    }
    int m;
    sscanf(argv[1], "%d", &m);
    int size = 2 * m;

    Matrix K = build_K(size);
    Matrix L = build_L(size);
    Matrix A = matrixMult(L, K);
    printf("Matrix A for m = %d:\n", m);
    printMatrix(A);

    Matrix p = powerIteration(A);
    printEigenvector(p);

    /*
    double pt = transpose(p);
    matrixMult(size, si)
    */

    freeMatrix(K);
    freeMatrix(L);
    freeMatrix(A);

    return EXIT_SUCCESS;
}
