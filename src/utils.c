#include <string.h>
void sparse_to_full(int n, int k, int p, const double *U, const int *iU, double *out){
    memset(out, 0, n * k * sizeof(double));

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < p; ++j){
            int col = iU[i * p + j];
            if (col == -1) break;
            out[i * k + col] = U[i * p + j];
        }
    }
}

void sparse_to_full_c(int n, int k, int p, const double *U, const int *iU, double *out){
    memset(out, 0, n * k * sizeof(double));

    for (int i = 0; i < n; ++i){
        for (int j = 0; j < p; ++j){
            int col = iU[i * p + j];
            if (col == -1) break;
            out[col * n + i] = U[i * p + j];
        }
    }
}

