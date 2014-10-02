#include "matrix.h"
#include "utils.h"

#define EPSILON 0.0000000000001

/* ========================================================================= */
/* Compressed sparse column matrix, a.k.a  Compressed column storage         */

CCSMatrix newCCSMatrix(int nnz, int rows, int cols);
CCSMatrix transposeCRS(CRSMatrix crs);
double ccsValueAt(int row, int col, CCSMatrix ccs);
bool reallocCCS(CCSMatrix ccs);
void ccsScalarMult(double scalar, CCSMatrix ccs);
static inline int reallocret(int estimate, CCSMatrix ret);

CCSMatrix
newCCSMatrix(int nnz, int rows, int cols){
    CCSMatrix ret      = malloc(sizeof(CCSMatrixCDT));
    ret->val           = malloc(nnz * sizeof (ret->val[0]));
    ret->row_index     = malloc(nnz * sizeof (ret->row_index[0]));
    ret->col_ptr       = malloc((cols +1) * sizeof (ret->col_ptr[0]));
    ret->nnz           = nnz;
    ret->rows          = rows;
    ret->cols          = cols;
    ret->col_ptr[cols] = nnz;

    return ret;
}

void
freeCCSMatrix(CCSMatrix mat){
    free(mat->val);
    free(mat->row_index);
    free(mat->col_ptr);
    free(mat);
}

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

CCSMatrix
build_CCS_A(int m){
    CCSMatrix K = build_CCS_K(2*m);
    CCSMatrix L = build_CCS_L(2*m);
    CCSMatrix A = ccsMult(K,L);
    freeCCSMatrix(K);
    freeCCSMatrix(L);
    return A;
}

CCSMatrix
build_CCS_K(int size){
    CCSMatrix K  = newCCSMatrix(4*(size/2), size, size);
    double alfa  = PI / 4;
    double a     = cos(alfa);
    double b     = sin(alfa);
    double c     = -b;

    for(int i = 0, j = 0, counter = 0; i < K->nnz; j += 2, counter += 4){
        K->col_ptr[j]     = counter;
        K->val[i]         = a;
        K->row_index[i++] = j;

        K->val[i]         = c;
        K->row_index[i++] = j + 1;

        K->col_ptr[j+1]   = counter + 2;
        K->val[i]         = b;
        K->row_index[i++] = j;

        K->val[i]         = a;
        K->row_index[i++] = j + 1;
    }

    return K;
}

CCSMatrix
build_CCS_L(int size){
    CCSMatrix L  = newCCSMatrix(4 + 4*((size/2)-1), size, size);
    double beta  = PI / 4;
    double a     = cos(beta);
    double b     = sin(beta);
    double c     = -b;

    L->val[0]              = a;
    L->row_index[0]        = 0;
    L->col_ptr[0]          = 0;

    L->val[1]              = b;
    L->row_index[1]        = L->rows - 1;

    L->val[L->nnz-2]       = c;
    L->row_index[L->nnz-2] = 0;
    L->col_ptr[L->cols-1]  = L->nnz-2;

    L->val[L->nnz-1]       = a;
    L->row_index[L->nnz-1] = L->rows - 1;

    for(int i = 2, j = 1, counter = 2; i < L->nnz - 2; j += 2, counter += 4){
        L->col_ptr[j] = counter;

        L->val[i] = a;
        L->row_index[i++] = j;

        L->val[i] = c;
        L->row_index[i++] = j+1;

        L->col_ptr[j+1] = counter + 2;

        L->val[i] = b;
        L->row_index[i++] = j;
        L->val[i] = a;
        L->row_index[i++] = j+1;
    }

    return L;
}

CCSMatrix
ccsMult(CCSMatrix ccs1, CCSMatrix ccs2){
    int elemEstimate = 6 + (ccs1->rows * ccs2->cols) * 0.1;
    CCSMatrix ret    = newCCSMatrix(elemEstimate, ccs1->rows, ccs2->cols);
    int nnz = 0;
    int stepsToBacktrack = 0;
    double currentVal = 0.0;
    bool columnStarterNotFound;

    // TODO ver si esto todavía es necesario !
    //for(i = 0; i < ret->cols; i++){
        //ret->col_ptr[i] = -1;
    //}

    for(int j = 0; j < ccs2->cols; j++){
        columnStarterNotFound = true;
        for(int i = 0; i < ccs1->rows; i++){
            currentVal = 0.0;

            // If no nonzero elements are present in the column...
            if(ccs2->col_ptr[j+1] - ccs2->col_ptr[j] == 0){
                // The element is zero, so just continue with the next one.
                stepsToBacktrack++;
                break;
            }

            for(int k = ccs2->col_ptr[j]; k < ccs2->col_ptr[j+1]; k++){
                double ccs1Val = ccsValueAt(i ,ccs2->row_index[k], ccs1);
                currentVal += ccs1Val * ccs2->val[k];
                //currentVal += ccsValueAt(i ,ccs2->row_index[k], ccs1) * ccs2->val[k];
            }

            if(fabs(currentVal) < EPSILON){
                continue;
            }

            if(columnStarterNotFound){
                ret->col_ptr[j] = nnz;
                while(stepsToBacktrack > 0){
                    ret->col_ptr[j-stepsToBacktrack] = nnz;
                    stepsToBacktrack--;
                }
                columnStarterNotFound = false;
            }

            /* Make sure enough memory has been allocated */
            if(nnz == elemEstimate){
                elemEstimate = reallocret(elemEstimate, ret);
            }

            ret->val[nnz] = currentVal;
            ret->row_index[nnz] = i;
            nnz++;
        }
        // Needed for the case in which even though there are elements in the
        // column the total sum is zero.
        //if(ret->col_ptr[j] == -1) stepsToBacktrack++;
    }
    ret->nnz = nnz;
    ret->col_ptr[ret->cols] = nnz;
    return ret;
}

double
ccsValueAt(int row, int col, CCSMatrix ccs){
    for(int i = ccs->col_ptr[col]; i < ccs->col_ptr[col+1]; i++){
        if(ccs->row_index[i] == row) return ccs->val[i];
    }
    return 0.0;
}

/* reallocates memory to use just the necessary amount */
bool
reallocCCS(CCSMatrix ccs){
    double *dp;
    int    *ip;
    int    *ip2;

    dp  = realloc(ccs->val, ccs->nnz * sizeof(ccs->val[0]));
    ip  = realloc(ccs->row_index, ccs->nnz * sizeof(ccs->row_index[0]));
    ip2 = realloc(ccs->col_ptr, (ccs->cols + 1) * sizeof(ccs->col_ptr[0]));

    if(!dp || !ip || !ip2) return false;

    ccs->val       = dp;
    ccs->row_index = ip;
    ccs->col_ptr   = ip2;
    return true;
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

static inline int
reallocret(int estimate, CCSMatrix ret){
    estimate *= 1.5;
    double *dp = realloc(ret->val, estimate * sizeof(ret->val[0]));
    int *ip = realloc(ret->row_index, estimate * sizeof(ret->row_index[0]));
    if(!dp || !ip) die("Out of memory while creating a ccs matrix from a dense one");
    ret->val = dp;
    ret->row_index = ip;
    return estimate;
}

CCSMatrix
matrixToCCS(Matrix mat){
    int elemEstimate = (mat->cols * mat->rows) * 0.1 + 6;
    CCSMatrix ret = newCCSMatrix(elemEstimate, mat->rows, mat->cols);

    for(int col = 0; col < mat->cols; col++){
        bool columnStarterNotFound = true;
        for(int row = 0; row < mat->rows; row++){
            if(mat->elem[row][col] != 0.0){
                if(ret->nnz == elemEstimate){
                    elemEstimate = reallocret(elemEstimate, ret);
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

void
print_CCSMatrix(CCSMatrix ccs){
    for(int i = 0; i < ccs->rows; i++){
        for(int j = 0; j < ccs->cols; j++){
            printf("%7.3f ", ccsValueAt(i, j, ccs));
        }
        putchar('\n');
    }
}

CCSMatrix
transposeCRS(CRSMatrix crs){
    CCSMatrix ret = newCCSMatrix(crs->nnz, crs->cols, crs->rows);
    memcpy(ret->val, crs->val, crs->nnz * sizeof(ret->val[0]));
    memcpy(ret->row_index, crs->col_index, crs->nnz * sizeof(ret->row_index[0]));
    memcpy(ret->col_ptr, crs->row_ptr, crs->nnz * sizeof(ret->col_ptr[0]));

    return ret;
}

/* FIXME Naive approach! if scalar is zero or small enough that the
 * product is rounded to zero, it breaks the data structure. */
void
ccsScalarMult(double scalar, CCSMatrix ccs){
    for(int i = 0; i < ccs->nnz; i++){
        ccs->val[i] *= scalar;
    }
}
