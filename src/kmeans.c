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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_MKL
#include "mkl.h"
#else
#include "cblas.h"
#endif
#include "DCRBR.h"

void kmeans(int n, int p, const double *X, int K, int *labels, double *H,
        kmeans_param opts);
int check_diff(int n, int *labels1, int *labels2);
int my_idamin(int n, const double *x);

void kmeans(int n, int p, const double *X, int K, int *labels, double *H,
        kmeans_param opts){
    /* MKL_INT seed[] = {23, (unsigned)time(NULL) % 4096, 23, 1};
    double *mu, *ones;
     **/
    
    /* memory allocation */
    int *cluster_num = (int*)malloc(K * sizeof(int));
    int *old_labels = (int*)malloc(n * sizeof(int));
    double *XH = (double*)malloc(n * K * sizeof(double));
    double *X_norm = (double*)malloc(n * sizeof(double));
    double *H_norm = (double*)malloc(K * sizeof(double));

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
            for (int i = 0; i < K; ++i){
                memcpy(H + i * p, X + (rand() % n) * p, p * sizeof(double));
            }
            break;
        case KMEANS_LABELS:
            memset(cluster_num, 0, K * sizeof(int));
            for (int i = 0; i < n; ++i){
                if (cluster_num[labels[i]] == 0){
                    memset(H + labels[i] * p, 0, p * sizeof(double));
                }
                cblas_daxpy(p, 1.0, X + i * p, 1, H + labels[i] * p, 1);
                ++cluster_num[labels[i]];
            }
            for (int i = 0; i < K; ++i){
                if (cluster_num[i] != 0){
                    cblas_dscal(p, 1.0 / cluster_num[i], H + i * p, 1);
                } else {
                    printf("cluster %d is empty!\n", i);
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
        memcpy(old_labels, labels, n * sizeof(int));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, K, p, -2.0,
                X, p, H, p, 0.0, XH, K);
        for (int i = 0; i < n; ++i){
            X_norm[i] = cblas_ddot(p, X + i * p, 1, X + i * p, 1);
        }
        for (int i = 0; i < K; ++i){
            H_norm[i] = cblas_ddot(p, H + i * p, 1, H + i * p, 1);
        }
        for (int i = 0; i < n; ++i){
            cblas_daxpy(K, 1.0, H_norm, 1, XH + i * K, 1);
        }
        for (int i = 0; i < K; ++i){
            cblas_daxpy(n, 1.0, X_norm, 1, XH + i, K);
        }

        for (int i = 0; i < n; ++i){
            //labels[i] = cblas_idamin(K, XH + i * K, 1);
            labels[i] = my_idamin(K, XH + i * K);
        }

        /* compute centroids */
        memset(cluster_num, 0, K * sizeof(int));
        for (int i = 0; i < n; ++i){
            if (cluster_num[labels[i]] == 0){
                memset(H + labels[i] * p, 0, p * sizeof(double));
            }
            cblas_daxpy(p, 1.0, X + i * p, 1, H + labels[i] * p, 1);
            ++cluster_num[labels[i]];
        }
        for (int i = 0; i < K; ++i){
            if (cluster_num[i] != 0){
                cblas_dscal(p, 1.0 / cluster_num[i], H + i * p, 1);
            } else {
                printf("iter %d cluster %d is empty!\n", iter, i);
            }
        }


        /* check stopping rule */
	int moved = check_diff(n, old_labels, labels);
        if (!moved){
            break;
        }
	if (opts.verbose){
	    printf("iter %d: moved %d\n", iter, moved);
	}
    }

    free(cluster_num);
    free(old_labels);
    free(XH);
    free(X_norm);
    free(H_norm);
    return;
}

int check_diff(int n, int *labels1, int *labels2){
    for (int i = 0; i < n; ++i){
        if (labels1[i] != labels2[i]){
	    return 1;
        }
    }
    return 0;
}

int my_idamin(int n, const double *x){
    int ipos = 0, v = x[0];
    for (int i = 1; i < n; ++i){
        if (fabs(x[i]) < v){
            v = fabs(x[i]);
            ipos = i;
        }
    }
    return ipos;
}
