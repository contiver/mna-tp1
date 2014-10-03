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

    printf("\nMatrix A for m = %d:\n", m);

    CCSMatrix A = build_Fast_CCS_A(m);
    /*
    CCSMatrix A2 = build_CCS_A(m);
    print_CCSMatrix(A);
    printf("esta es A2\n");
    print_CCSMatrix(A2);
    //print_CCSMatrix(A);

//    printMatrix(A);

*/
    return EXIT_SUCCESS;
}
