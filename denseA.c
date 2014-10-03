#include "matrix.h"

void usage(void);

void
usage(void){
    printf("usage: densaA m [p]\n"
            "      m    number to use when creating the 2m x 2m matrix\n"
            "      p    optional. Whether to print matrix to stdout or not\n"
            "\n"
            "Example: denseA 10 p\n\n");
}

int
main(int argc, char *argv[]){
    if (argc < 2){
        usage();
        return EXIT_SUCCESS;
    }
    int m;
    sscanf(argv[1], "%d", &m);

    Matrix A = build_A(m);

    if(argc > 2){
        char p;
        sscanf(argv[2], "%c", &p);
        if(p == 'p'){
            printf("\nMatrix A for m = %d:\n", m);
            printMatrix(A);
        }
    }
    return EXIT_SUCCESS;
}
