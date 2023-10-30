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

#include "adj_matrix.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

adj_matrix * adj_matrix_create(RBR_INT n, RBR_INT nnz){
    adj_matrix * A = (adj_matrix*)calloc(1, sizeof(adj_matrix));

    if (A == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        return NULL;
    }

    A->n = n;
    A->nnz = nnz;
    A->pntr = (RBR_INT*)malloc((n + 1) * sizeof(RBR_INT));
    A->indx = (RBR_INT*)malloc(nnz * sizeof(RBR_INT));
    A->val = (double*)malloc(nnz * sizeof(double));
    A->d = (double*)malloc(n * sizeof(double));

    if (A->pntr == NULL || A->indx == NULL || A->val == NULL || A->d == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        adj_matrix_destroy(A);
        return NULL;
    }

    return A;
}

/**
 * @fn int read_adj_matrixcsr(const char *filename, adj_matrix *A)
 * @brief Read the sparse matrix from file. The sparse matrix should have
 * Rutherford Boeing format.
 * @param filename The matrix filename.
 * @param one_based set 1 to generate a one-based matrix.
 * @param weighted set 1 to indicate that this is a weighted adjacency matrix
 * @return (pointer to adj_matrix) The function returns the pointer to adj_matrix
 * when successfully executed, otherwise returns NULL, in which case
 * either it fails to open the file or the file is illegal.
 */

adj_matrix * read_adj_matrix_csr(const char *filename, int one_based, int weighted){
    FILE *fp = fopen(filename, "r");
    RBR_INT nnz, row, col;
    char headerl1[82], headerl2[82], headerl3[82], headerl4[82], which[4];
    RBR_INT *pntr, *indx, cnt;
    double *val;
    adj_matrix * A = NULL;

    if (fp == NULL){
        fprintf(stderr, "[ERROR] %s: cannot open file %s. \n", __func__, filename);
        return NULL;
    }

    if (fgets(headerl1, 82, fp) == NULL){
        fprintf(stderr, "[ERROR] %s: Line 1 is invalid or missing.\n", __func__);
        goto cleanup;
    }
    if (fgets(headerl2, 82, fp) == NULL){
        fprintf(stderr, "[ERROR] %s: Line 2 is invalid or missing.\n", __func__);
        goto cleanup;
    }
    if (fgets(headerl3, 82, fp) == NULL){
        fprintf(stderr, "[ERROR] %s: Line 3 is invalid or missing.\n", __func__);
        goto cleanup;
    }
    if (fgets(headerl4, 82, fp) == NULL){
        fprintf(stderr, "[ERROR] %s: Line 4 is invalid or missing.\n", __func__);
        goto cleanup;
    }
    sscanf(headerl3, "%3s " _fmt_int " " _fmt_int " " _fmt_int, which, &row, &col, &nnz);

    /* check the shape */
    if (row != col){
        fprintf(stderr, "[ERROR] %s: The adj matrix must be square.(got " _fmt_int " rows and " _fmt_int " cols)\n", __func__, row, col);
        goto cleanup;
    }

    A = calloc(1, sizeof(adj_matrix));

    if (A == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        goto cleanup;
    }

    if (strcmp(which, "iua") == 0 || strcmp(which, "rua") == 0){
        strncpy(A->desc, headerl1, 72);
        strncpy(A->label, headerl1 + 72, 8);
        A->n = col;
        A->nnz = nnz;

        /* allocating memory */
        A->indx = (RBR_INT*)malloc(nnz * sizeof(RBR_INT));
        A->pntr = (RBR_INT*)malloc((col + 1) * sizeof(RBR_INT));
        A->d    = (double*)malloc(col * sizeof(double));
        if (A->indx == NULL || A->pntr == NULL || A->d == NULL){
            fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
            goto cleanup;
        }
        pntr = A->pntr;
        indx = A->indx;

        /* input PNTR */
        for(RBR_INT i = 0; i <= col; ++i){
            cnt = fscanf(fp, _fmt_int, pntr + i);

            if (cnt < 1){
                fprintf(stderr, "[ERROR] %s: invalid inputs or EOF reached\n", __func__);
                goto cleanup;
            }

            if (i == 0){
                one_based ^= pntr[0];
            }
            pntr[i] -= one_based;
            /* update degree one by one */
            if (i > 0 && weighted == 0){
                A->d[i-1] = pntr[i] - pntr[i-1];
            }
        }
        /* input INDX */
        for(RBR_INT i = 0; i < nnz; ++i){
            cnt = fscanf(fp, _fmt_int, indx + i);
            if (cnt < 1){
                fprintf(stderr, "[ERROR] %s: invalid inputs or EOF reached\n", __func__);
                goto cleanup;
            }
            indx[i] -= one_based;
        }

        /* input VAL for weighted graph and A->d */
        if (weighted){
            A->val = (double*)malloc(nnz * sizeof(double));

            if (A->val == NULL){
                fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
                goto cleanup;
            }

            val = A->val;
            for (RBR_INT i = 0; i < nnz; ++i){
                cnt = fscanf(fp, "%lg", val + i);
                if (cnt < 1){
                    fprintf(stderr, "[ERROR] %s: invalid inputs or EOF reached\n", __func__);
                    goto cleanup;
                }
            }

            /* compute A->d = A * ones */
            for (RBR_INT i = 0; i < col; i++){
                A->d[i] = 0;
                for (RBR_INT j = pntr[i]; j < pntr[i+1]; ++j){
                    A->d[i] += val[j];
                }
            }
        } else {
            A->val = NULL;
        }

    }

    // format not supported
    else{
        fprintf(stderr, "[ERROR] %s: format error. Rutherford Boeing CSR format matrix needed.\n", __func__);
        goto cleanup;
    }

    if (0){
cleanup:
        adj_matrix_destroy(A);
        fclose(fp);
        return NULL;
    }

    fclose(fp);
    return A;
}

void adj_matrix_destroy(adj_matrix *A){
    if (A == NULL) return;
    if (A->pntr != NULL) free(A->pntr);
    if (A->indx != NULL) free(A->indx);
    if (A->d != NULL) free(A->d);
    if (A->val != NULL) free(A->val);

    free(A);
}

