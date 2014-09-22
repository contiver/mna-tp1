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

    Matrix A = build_A(m);

    printf("\nMatrix A for m = %d:\n", m);
    printMatrix(A);

    //Matrix p = powerIteration(A);
    //printEigenvector(p);

    freeMatrix(A);

    return EXIT_SUCCESS;
}
