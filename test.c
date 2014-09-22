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
    CCSMatrix L = build_CCS_L(2*m);
    CCSMatrix A = ccsMult(K, L);

    printf("Matrix K is: \n");
    print_CCSMatrix(K);
    printf("Matrix L is: \n");
    print_CCSMatrix(L);
    printCCS(A);
    print_CCSMatrix(A);

    printf("The real Matrix A is: \n");
    printMatrix(build_A(m));

    return EXIT_SUCCESS;
}
