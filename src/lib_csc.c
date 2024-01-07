#include "lib_csc.h"

CSCMatrix convertToCSC(int *poissonMatrix, int n) {
    CSCMatrix cscMatrix;

    int nnz = 3 * n - 2;

    cscMatrix.values = (int *)malloc(nnz * sizeof(int));
    cscMatrix.row_indices = (int *)malloc(nnz * sizeof(int));
    cscMatrix.column_pointers = (int *)malloc((n + 1) * sizeof(int));
    cscMatrix.nnz = nnz;
    cscMatrix.n = n;

    int index = 0;  

    for (int i = 0; i < n; i++) {
        cscMatrix.values[index] = 4;
        cscMatrix.row_indices[index] = i;
        index++;

        if (i > 0) {
            cscMatrix.values[index] = -1;
            cscMatrix.row_indices[index] = i - 1;
            index++;
        }

        if (i < n - 1) {
            cscMatrix.values[index] = -1;
            cscMatrix.row_indices[index] = i + 1;
            index++;
        }
    }

    cscMatrix.column_pointers[0] = 0;
    for (int i = 1; i <= n; i++) {
        cscMatrix.column_pointers[i] = cscMatrix.column_pointers[i - 1] + 3;
    }

    return cscMatrix;
}

void dcscmv(const CSCMatrix *matrix, const int *x, int *y) {
    for (int j = 0; j < matrix->n; j++) {
        y[j] = 0;
        for (int i = matrix->column_pointers[j]; i < matrix->column_pointers[j + 1]; i++) {
            y[matrix->row_indices[i]] += matrix->values[i] * x[j];
        }
    }
}




