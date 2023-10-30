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

/*
 * ==========================================================================
 *
 *       Filename:  kmeans.c
 *
 *    Description:  perform K-means algorithm on standard data stuctures
 *
 *        Version:  1.0
 *        Created:  01/13/2019 07:13:51 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BICMR, PKU
 *
 * ==========================================================================
 */

#include "rbr_subroutines.h"
#include "rbr_blas.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// functions used only in this source file
int check_diff(RBR_INT n, const RBR_INT *labels1, const RBR_INT *labels2);
RBR_INT my_idamin(RBR_INT n, const double *x);

int kmeans(RBR_INT n, RBR_INT p, const double *X, RBR_INT K, RBR_INT *labels, double *H,
        kmeans_param opts){
    int exit_code = 0;
    
    /* memory allocation */
    RBR_INT *cluster_num = (RBR_INT*)malloc(K * sizeof(RBR_INT));
    RBR_INT *old_labels = (RBR_INT*)malloc(n * sizeof(RBR_INT));
    double *XH = (double*)malloc(n * K * sizeof(double));
    double *X_norm = (double*)malloc(n * sizeof(double));
    double *H_norm = (double*)malloc(K * sizeof(double));

    if (cluster_num == NULL || old_labels == NULL || XH == NULL || X_norm == NULL || H_norm == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        exit_code = 1;
        goto cleanup;
    }

    /* initialization */
    switch (opts.init){
        /* 
        case KMEANS_CENTER:
            mu = (double*)calloc(p, sizeof(double));
            ones = (double*)malloc(n * sizeof(double));
            for (int i = 0; i < n; ++i) ones[i] = 1;

            cblas_dgemv(CblasRowMajor, CblasTrans, n, p, 1.0, X, p, ones, 1, 0.0, mu, 1);
            cblas_dscal(p, 1.0 / n, mu, 1);
            for (int i = 0; i < K; ++i){
                LAPACKE_dlarnv(3, seed, p, H + i * p);
                cblas_daxpy(p, 1.0, mu, 1, H + i * p, 1);
            }
            free(mu); free(ones);
            break;
            */
        case KMEANS_SAMPLE:
            for (RBR_INT i = 0; i < K; ++i){
                memcpy(H + i * p, X + (rand() % n) * p, p * sizeof(double));
            }
            break;
        case KMEANS_LABELS:
            memset(cluster_num, 0, K * sizeof(RBR_INT));
            for (RBR_INT i = 0; i < n; ++i){
                if (cluster_num[labels[i]] == 0){
                    memset(H + labels[i] * p, 0, p * sizeof(double));
                }
                cblas_daxpy(p, 1.0, X + i * p, 1, H + labels[i] * p, 1);
                ++cluster_num[labels[i]];
            }
            for (RBR_INT i = 0; i < K; ++i){
                if (cluster_num[i] != 0){
                    cblas_dscal(p, 1.0 / cluster_num[i], H + i * p, 1);
                } else {
                    fprintf(stderr, "[WARNING] %s: cluster %d is empty!\n", __func__, (int)i);
                }
            }
            break;
        case KMEANS_USER:
            /* user supplied initial centroid, do nothing */
            break;
        default:
            break;
    }


    /* main loop */
    for (int iter = 0; iter < opts.maxit; ++iter){
        /* re-assign clusters */
        memcpy(old_labels, labels, n * sizeof(RBR_INT));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, K, p, -2.0,
                X, p, H, p, 0.0, XH, K);
        for (RBR_INT i = 0; i < n; ++i){
            X_norm[i] = cblas_ddot(p, X + i * p, 1, X + i * p, 1);
        }
        for (RBR_INT i = 0; i < K; ++i){
            H_norm[i] = cblas_ddot(p, H + i * p, 1, H + i * p, 1);
        }
        for (RBR_INT i = 0; i < n; ++i){
            cblas_daxpy(K, 1.0, H_norm, 1, XH + i * K, 1);
        }
        for (RBR_INT i = 0; i < K; ++i){
            cblas_daxpy(n, 1.0, X_norm, 1, XH + i, K);
        }

        for (RBR_INT i = 0; i < n; ++i){
            //labels[i] = cblas_idamin(K, XH + i * K, 1);
            labels[i] = my_idamin(K, XH + i * K);
        }

        /* compute centroids */
        memset(cluster_num, 0, K * sizeof(RBR_INT));
        for (RBR_INT i = 0; i < n; ++i){
            if (cluster_num[labels[i]] == 0){
                memset(H + labels[i] * p, 0, p * sizeof(double));
            }
            cblas_daxpy(p, 1.0, X + i * p, 1, H + labels[i] * p, 1);
            ++cluster_num[labels[i]];
        }
        for (RBR_INT i = 0; i < K; ++i){
            if (cluster_num[i] != 0){
                cblas_dscal(p, 1.0 / cluster_num[i], H + i * p, 1);
            } else {
                fprintf(stderr, "[WARNING] %s: iter %d cluster %d is empty!\n", __func__, iter, (int)i);
            }
        }


        /* check stopping rule */
        int moved = check_diff(n, old_labels, labels);
        if (!moved){
            break;
        }
        if (opts.verbose){
            printf("[INFO] %s: iter %d: moved %d\n", __func__, iter, moved);
        }
    }

cleanup:
    free(cluster_num);
    free(old_labels);
    free(XH);
    free(X_norm);
    free(H_norm);
    return exit_code;
}

int check_diff(RBR_INT n, const RBR_INT *labels1, const RBR_INT *labels2){
    for (RBR_INT i = 0; i < n; ++i){
        if (labels1[i] != labels2[i]){
	    return 1;
        }
    }
    return 0;
}

RBR_INT my_idamin(RBR_INT n, const double *x){
    if (n == 0){
        return -1;
    }

    RBR_INT ipos = 0, v = x[0];
    for (RBR_INT i = 1; i < n; ++i){
        if (fabs(x[i]) < v){
            v = fabs(x[i]);
            ipos = i;
        }
    }
    return ipos;
}
