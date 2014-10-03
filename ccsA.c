#include "matrix.h"

void usage(void);

void
usage(void){
    printf("usage: ccsA m [p]\n"
            "      m    number to use when creating the 2m x 2m matrix\n"
            "      p    optional. Whether to print matrix to stdout or not\n"
            "\n"
            "Example: ccsA 10 p\n\n");
}

int
main(int argc, char *argv[]){
    if (argc < 2){
        usage();
        return EXIT_SUCCESS;
    }
    int m;
    sscanf(argv[1], "%d", &m);

    CCSMatrix A = build_CCS_A(m);

    if(argc > 2){
        char p;
        sscanf(argv[2], "%c", &p);
        if(p == 'p'){
            printf("\nCCS Matrix A for m = %d:\n", m);
            print_CCSMatrix(A);
        }
    }
    return EXIT_SUCCESS;
}

