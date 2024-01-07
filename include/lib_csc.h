#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "atlas_headers.h"

typedef struct {
    int *values;
    int *row_indices;
    int *column_pointers;
    int nnz;  // Number of non-zero elements
    int n;    // Size of the matrix (assumed square)
} CSCMatrix;

CSCMatrix convertToCSC(int *poissonMatrix, int n);
void dcscmv(const CSCMatrix *matrix, const int *x, int *y) ;