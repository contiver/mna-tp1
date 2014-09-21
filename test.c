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

    CCSMatrix K = build_CCS_K(2*m);
    CCSMatrix K2 = build_CCS_K(2*m);
    CCSMatrix M = ccsMult(K, K2);

    printCCS(M);
    Matrix R = ccsToMatrix(M);
    printMatrix(R);

    printf("The right value of K * K is\n");
    printMatrix(matrixMult(build_K(2*m), build_K(2*m)));

    return EXIT_SUCCESS;
}
