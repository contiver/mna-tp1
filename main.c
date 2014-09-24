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

    CCSMatrix A = build_CCS_A(m);

    printf("\nMatrix A for m = %d:\n", m);
    printCCS(A);
    print_CCSMatrix(A);

    freeCCSMatrix(A);

    return EXIT_SUCCESS;
}
