#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "atlas_headers.h"

typedef struct {
    int *values;
    int *column_indices;
    int *row_pointers;
    int nnz;  // Number of non-zero elements
    int n;    // Size of the matrix (assumed square)
} CSRMatrix;

CSRMatrix convertToCSR(int *poissonMatrix, int n);
void dcsrmv(const CSRMatrix *matrix, const int *x, int *y) ;