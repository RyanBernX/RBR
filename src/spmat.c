#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "spmat.h"

/** @fn int read_adj_matrixcsr(const char *filename, adj_matrix *A)
 *  @brief Read the sparse matrix from file. The sparse matrix should have
 *  Rutherford Boeing format.
 *  @param filename The matrix filename.
 *  @param A Pointer to the sparse matrix object.
 *  @param one_based set 1 to generate a one-based matrix.
 *  @return (int) The function returns 0 if it ends with no errors. Otherwise,
 *  either it fails to open the file or the file is illegal.
 */

int read_adj_matrix_csr(const char *filename, adj_matrix* A, int one_based){
    FILE *fp = fopen(filename, "r");
    int nnz, row, col;
    char headerl1[82], headerl2[82], headerl3[82], headerl4[82], which[4];
    int *pntr, *indx;

    if (fp == NULL){
        fprintf(stderr, "Cannot open file %s. \n", filename);
        return -1;
    }

    fgets(headerl1, 82, fp);
    fgets(headerl2, 82, fp);
    fgets(headerl3, 82, fp);
    fgets(headerl4, 82, fp);
    sscanf(headerl3, "%3s %d %d %d", which, &row, &col, &nnz);

    /* check the shape */
    if (row != col){
        fprintf(stderr, "The adj matrix must be square.(got %d rows and %d cols)\n", row, col);
        fclose(fp);
        return -1;
    }

    if(strcmp(which, "rsa") == 0 || strcmp(which, "rua") == 0){
        strncpy(A->desc, headerl1, 72);
        strncpy(A->label, headerl1 + 72, 8);
        A->n = col; A->nnz = nnz;

        /* allocating memory */
        A->indx = (int*)malloc(nnz * sizeof(int));
        A->pntr = (int*)malloc((col + 1) * sizeof(int));
        A->d    = (double*)malloc(col * sizeof(double));
        pntr = A->pntr;
        indx = A->indx;

        /* input PNTR */
        for(int i = 0; i <= col; ++i){
            fscanf(fp, "%d", pntr + i);
            if (i == 0){
                one_based ^= pntr[0];
            }
            pntr[i] -= one_based;
            /* update degree one by one */
            if (i > 0) A->d[i-1] = pntr[i] - pntr[i-1];
            
        }
        /* input INDX */
        for(int i = 0; i < nnz; ++i){
            fscanf(fp, "%d", indx + i);
            indx[i] -= one_based;
        }
    }
    else{
        fprintf(stderr, "Format error. Rutherford Boeing CSR format matrix needed.\n");
        fclose(fp);
        return -1;
    }
    fclose(fp);
    return 0;
}

void adj_matrix_destroy(adj_matrix *A){
    free(A->pntr);
    free(A->indx);
    free(A->d);
}
