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

    CCSMatrix id = identityCCSMatrix(m);

    printCCS(id);

    Matrix idMat = ccsToMatrix(id);
    printMatrix(idMat);

    //Matrix p = powerIteration(A);
    //printEigenvector(p);

    freeMatrix(idMat);


    // TODO: REVISAR ESTE FREE DEL CCS
    freeCCSMatrix(id);

    return EXIT_SUCCESS;
}
