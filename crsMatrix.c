#include "matrix.h"

/* ========================================================================= */
/* Compressed sparse row matrix, a.k.a Compressed row storage                */

double crsValueAt(int row, int col, CRSMatrix crs);

static CRSMatrix
newCRSMatrix(int nnz, int rows, int cols){
    CRSMatrix ret      = malloc(sizeof(CRSMatrixCDT));
    ret->val           = malloc(nnz * sizeof (ret->val[0]));
    ret->col_index     = malloc(nnz * sizeof (ret->col_index[0]));
    ret->row_ptr       = malloc((rows +1) * sizeof (ret->row_ptr[0]));
    ret->nnz           = nnz;
    ret->rows          = rows;
    ret->cols          = cols;
    ret->row_ptr[rows] = nnz;

    return ret;
}

CRSMatrix
transposeCCS(CCSMatrix ccs){
    CRSMatrix ret = newCRSMatrix(ccs->nnz, ccs->cols, ccs->rows);
    memcpy(ret->val, ccs->val, ccs->nnz * sizeof(ret->val[0]));
    memcpy(ret->col_index, ccs->row_index, ccs->nnz * sizeof(ret->col_index[0]));
    memcpy(ret->row_ptr, ccs->col_ptr, ccs->nnz * sizeof(ret->row_ptr[0]));
    
    return ret;
}

void
print_CRSMatrix(CRSMatrix crs){
    for(int i = 0; i < crs->rows; i++){
        for(int j = 0; j < crs->cols; j++){
            printf("%7.3f ", crsValueAt(i, j, crs));
        }
        putchar('\n');
    }
}

double
crsValueAt(int row, int col, CRSMatrix crs){
    for(int i = crs->row_ptr[row]; i < crs->row_ptr[row+1]; i++){
        if(crs->col_index[i] == col) return crs->val[i];
    }
    return 0.0;
}

/* Maybe traversing just once, creating a (value,row,col) tuple 
 * and inserting it in some self sorting structure like a tree
 * could make this faster */
CCSMatrix
crsToCcs(CRSMatrix crs){
    CCSMatrix ret = newCCSMatrix(crs->nnz, crs->rows, crs->cols);

    for(int i = 0, pos = 0; i < crs->cols; i++){
        bool columnStarterNotFound = true;
        for(int j = 0; j < crs->rows; j++){
            int elemsInRow = crs->row_ptr[j+1] - crs->row_ptr[j];
            if(elemsInRow == 0){
                continue;
            }
            for(int k = crs->row_ptr[j]; k < crs->row_ptr[j+1]; k++){
                if(crs->col_index[k] == j){
                    ret->val[pos] = crs->val[k];
                    ret->row_index[pos] = j;
                    if( columnStarterNotFound ){
                        ret->col_ptr[i] = pos;
                        columnStarterNotFound = false;
                    }
                    pos++;
                }
            }
        }
    }

    return ret;
}
