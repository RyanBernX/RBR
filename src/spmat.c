/*
 * RBR
 *
 * Copyright (C) 2019  Haoyang Liu (liuhaoyang@pku.edu.cn)
 *                     Zaiwen Wen  (wenzw@pku.edu.cn)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.\
 */

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
 *  @param weighted set 1 to indicate that this is a weighted adjacency matrix
 *  @return (int) The function returns 0 if it ends with no errors. Otherwise,
 *  either it fails to open the file or the file is illegal.
 */

int read_adj_matrix_csr(const char *filename, adj_matrix* A, int one_based, int weighted){
    FILE *fp = fopen(filename, "r");
    int nnz, row, col;
    char headerl1[82], headerl2[82], headerl3[82], headerl4[82], which[4];
    int *pntr, *indx;
    double *val;

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
            if (i > 0 && weighted == 0) A->d[i-1] = pntr[i] - pntr[i-1];
            
        }
        /* input INDX */
        for(int i = 0; i < nnz; ++i){
            fscanf(fp, "%d", indx + i);
            indx[i] -= one_based;
        }
        /* input VAL for weighted graph and A->d */
        if (weighted){
            A->val = (double*)malloc(nnz * sizeof(double));
            val = A->val;
            for (int i = 0; i < nnz; ++i){
                fscanf(fp, "%lg", val + i);
            }
            
            /* compute A->d = A * ones */
            for (int i = 0; i < col; i++){
                A->d[i] = 0;
                for (int j = pntr[i]; j < pntr[i+1]; ++j){
                    A->d[i] += val[j];
                }
            }
        }
        else {
            A->val = NULL;
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
    if (A->val != NULL) free(A->val);
}
