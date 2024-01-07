#include "lib_csr.h"


CSRMatrix convertToCSR(int *poissonMatrix, int n) {
    CSRMatrix csrMatrix;

    int nnz = 3 * n - 2;

    csrMatrix.values = (int *)malloc(nnz * sizeof(int));
    csrMatrix.column_indices = (int *)malloc(nnz * sizeof(int));
    csrMatrix.row_pointers = (int *)malloc((n + 1) * sizeof(int));
    csrMatrix.nnz = nnz;
    csrMatrix.n = n;

    int index = 0;  

    for (int i = 0; i < n; i++) {
        csrMatrix.values[index] = 4;
        csrMatrix.column_indices[index] = i;
        index++;

        if (i > 0) {
            csrMatrix.values[index] = -1;
            csrMatrix.column_indices[index] = i - 1;
            index++;
        }

        if (i < n - 1) {
            csrMatrix.values[index] = -1;
            csrMatrix.column_indices[index] = i + 1;
            index++;
        }
    }

    int row_start = 0;
    for (int i = 0; i < n; i++) {
        csrMatrix.row_pointers[i] = row_start;

        int row_nnz = 3;
        if (i == 0 || i == n - 1) {
            row_nnz = 2;
        }

        row_start += row_nnz;
    }
    csrMatrix.row_pointers[n] = row_start;  
    return csrMatrix;
}

void dcsrmv(const CSRMatrix *matrix, const int *x, int *y) {
    for (int i = 0; i < matrix->n; i++) {
        y[i] = 0;
        for (int j = matrix->row_pointers[i]; j < matrix->row_pointers[i + 1]; j++) {
            y[i] += matrix->values[j] * x[matrix->column_indices[j]];
        }
    }
}