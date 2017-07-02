#ifndef __SPMAT_H_
#define __SPMAT_H_

typedef struct{
    int *indx;
    int *pntr;
    double *d;
    int nnz;
    int n;
    char label[9];
    char desc[73];
} adj_matrix;

#ifdef __cplusplus
extern "C" {
#endif


int read_adj_matrix_csr(const char *, adj_matrix *, int);
void adj_matrix_destroy(adj_matrix *);

#ifdef __cplusplus
}
#endif

#endif

