#include "matrix.h"

/* ========================================================================= */
/* Compressed sparse column matrix, a.k.a
 * Compressed column storage */

static CCSMatrix newCCSMatrix(int nnz, int rows, int cols);

CCSMatrix
identityCCSMatrix(int size){
    CCSMatrix ret = newCCSMatrix(size, size, size);

    for(int i = 0; i < size; i++){
        ret->row_index[i] = i;
        ret->col_ptr[i]   = i;
        ret->val[i]       = 1.0;
    }
    ret->col_ptr[size] = size;

    return ret;
}

static CCSMatrix
newCCSMatrix(int nnz, int rows, int cols){
    CCSMatrix ret  = malloc(sizeof(CCSMatrixCDT));
    ret->val       = malloc(nnz * sizeof (ret->val[0]));
    ret->row_index = malloc(nnz * sizeof (ret->row_index[0]));
    ret->col_ptr   = malloc((cols +1) * sizeof (ret->col_ptr[0]));
    ret->nnz       = nnz;
    ret->rows      = rows;
    ret->cols      = cols;

    ret->col_ptr[cols] = nnz;

    return ret;
}

CCSMatrix
build_CCS_K(int size){
    CCSMatrix K  = newCCSMatrix(4*(size/2), size, size);

    double alfa = PI / 4;

    double a = cos(alfa);
    double b = sin(alfa);
    double c = -b;

    int i = 0;
    int j = 0;
    int counter = 0;
    while(i < K->nnz){
        K->col_ptr[j] = counter;

        K->val[i] = a;
        K->row_index[i++] = j;
        K->val[i] = c;
        K->row_index[i++] = j + 1;

        K->col_ptr[j+1] = counter + 2;

        K->val[i] = b;
        K->row_index[i++] = j;
        K->val[i] = a;
        K->row_index[i++] = j + 1;
        j += 2;
        counter += 4;
    }

    return K;
}

/*
CCSMatrix
build_CCS_L(int size){
    CCSMatrix L  = newCCSMatrix(4 + 4*((size/2)-1), size, size);

    double beta = PI / 4;

    double a = cos(beta);
    double b = sin(beta);
    double c = -b;

    L->val[0]        = cos(beta);
    L->val[1]        = -sin(beta);
    L->val[L->nnz-2] = sin(beta);
    L->val[L->nnz-1] = cos(beta);

    int i = 0;
    int j = 0;
    int counter = 0;
    while(i < L->nnz){
        L->col_ptr[j] = counter;

        L->val[i] = a;
        L->row_index[i++] = j;
        L->val[i] = c;
        L->row_index[i++] = j + 1;

        L->col_ptr[j+1] = counter + 2;

        L->val[i] = b;
        L->row_index[i++] = j;
        L->val[i] = a;
        L->row_index[i++] = j + 1;
        j += 2;
        counter += 4;
    }

    return L;
}

Matrix
build_L(int size){
    Matrix L = nullMatrix(size, size);

    double beta = PI / 4;

    size--;
    L->elem[0][0]       = cos(beta);
    L->elem[0][size]    = -sin(beta);
    L->elem[size][0]    = sin(beta);
    L->elem[size][size] = cos(beta);

    while(size > 2){
        L->elem[size-1][size-1] = a;
        L->elem[size-1][size-2] = c;
        L->elem[size-2][size-1] = b;
        L->elem[size-2][size-2] = a;
        size -= 2;
    }
    return L;
}
*/

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
    for(int i = 0; i <= ccs->cols; i++){
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

    for(int col = 0; col < ccs->cols; col++){
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
    int elemEstimate = (mat->cols * mat->rows) * 0.1 + 6;

    CCSMatrix ret    = malloc(sizeof(CCSMatrixCDT));
    ret->val         = malloc(elemEstimate * sizeof(ret->val[0]));
    ret->row_index   = malloc(elemEstimate * sizeof(ret->row_index[0]));
    ret->col_ptr     = malloc((mat->cols + 1) * sizeof(ret->col_ptr[0]));
    ret->cols        = mat->cols;
    ret->rows        = mat->rows;

    double * dp;
    int * ip;

    for(int col = 0; col < mat->cols; col++){
        bool columnStarterNotFound = true;
        for(int row = 0; row < mat->rows; row++){
            if(mat->elem[row][col] != 0.0){
                if(ret->nnz == elemEstimate){
                    elemEstimate *= 1.3;     // Arbitrary incremental value
                    dp = realloc(ret->val, elemEstimate * sizeof(ret->val[0]));
                    ip = realloc(ret->row_index, elemEstimate * sizeof(ret->row_index[0]));
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
