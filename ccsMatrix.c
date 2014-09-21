#include "matrix.h"

/* ========================================================================= */
/* Compressed sparse column matrix, a.k.a
 * Compressed column storage */

CCSMatrix
identityCCSMatrix(int size){
    CCSMatrix ret  = malloc(sizeof(CCSMatrixCDT));
    ret->val       = malloc(size * sizeof (ret->val[0]));
    ret->row_index = malloc(size * sizeof (ret->row_index[0]));
    ret->col_ptr   = malloc((size+1) * sizeof (ret->col_ptr[0]));
    ret->nnz       = size;
    ret->rows      = size;
    ret->cols      = size;

    for(int i = 0; i < size; i++){
        ret->row_index[i] = i;
        ret->col_ptr[i]   = i;
        ret->val[i]       = 1.0;
    }
    ret->col_ptr[size] = size;

    return ret;
}

void
freeCCSMatrix(CCSMatrix mat){
    free(mat->val);
    free(mat->row_index);
    free(mat->col_ptr);
    free(mat);
}

void
printCCS(CCSMatrix ccs){
    printf("val vector:\n");
    for(int i = 0; i < ccs->nnz; i++){
        printf("%.2f ", ccs->val[i]);
    }
    printf("\n");

    printf("row_index vector:\n");
    for(int i = 0; i < ccs->nnz; i++){
        printf("%d ", ccs->row_index[i]);
    }
    printf("\n");

    printf("col_ptr vector:\n");
    for(int i = 0; i <= ccs->nnz; i++){
        printf("%d ", ccs->col_ptr[i]);
    }
    printf("\n");

    printf("ccs->nnz: %d\n", ccs->nnz);
    printf("ccs->rows: %d\n", ccs->rows);
    printf("ccs->cols: %d\n", ccs->cols);
}

Matrix
ccsToMatrix(CCSMatrix ccs){
    Matrix ret = nullMatrix(ccs->rows, ccs->rows);

    int j = 0;
    int elements = 0;

    for(int col = 0; col < ccs->nnz; col++){
        elements += ccs->col_ptr[col+1] - ccs->col_ptr[col];

        while(j < elements){
            ret->elem[ccs->row_index[j]][col] = ccs->val[j];
            j++;
        }
    }
    return ret;
}

CCSMatrix
matrixToCCS(Matrix mat){
    CCSMatrix ret  = malloc(sizeof(CCSMatrixCDT));
    int memSize    = (mat->cols * mat->rows) * 0.1 + 2;
    ret->val       = malloc(memSize * sizeof(ret->val[0]));
    ret->row_index = malloc(memSize * sizeof(ret->row_index[0]));
    ret->col_ptr   = malloc((mat->cols + 1) * sizeof(ret->col_ptr[0]));
    ret->cols      = mat->cols;
    ret->rows      = mat->rows;

    double * dp;
    int * ip;

    for(int col = 0; col < mat->cols; col++){
        bool columnStarterNotFound = true;
        for(int row = 0; row < mat->rows; row++){
            if(mat->elem[row][col] != 0.0){
                if(ret->nnz == memSize){
                    memSize *= 1.2;     // Arbitrary incremental value
                    dp       = realloc(ret->val, memSize * sizeof(ret->val[0]));
                    ip       = realloc(ret->row_index, memSize * sizeof(ret->row_index[0]));
                    if(dp == NULL || ip == NULL){ 
                        perror("Not enough memory"); 
                        exit(EXIT_FAILURE);
                    }
                    ret->val = dp;
                    ret->row_index = ip;
                }
                ret->val[ret->nnz] = mat->elem[row][col];
                ret->row_index[ret->nnz] = row;
                if(columnStarterNotFound){
                    ret->col_ptr[col] = ret->nnz;
                }
                columnStarterNotFound = false;
                ret->nnz++;
            }
        }
    }
    ret->col_ptr[ret->cols] = ret->nnz;

    return ret;
}

/*
   void
   print_CCSMatrix(CCSMatrix ccs){
   int val_index = 0;

   int i = 0, j = 0;
   while(i < ccs->rows){
   while(j < ccs->cols){
   if(ccs->row_index[ccs->col_ptr[j]] == i)
   if(ccs->row_index[val_index] == i && ccs->col_ptr[j] == ){
   printf("");
   }else{
   printf("0 ");
   }

   }
   }
   if( i == ccs->row_index[i] )

   }
   }
   }
   */
