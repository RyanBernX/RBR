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

#include "RBR.h"
#include "rbr_blas.h"
#include "rbr_subroutines.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

rbr_result * rbr(const adj_matrix *A, RBR_INT k, const rbr_param * param) {
    /* k is the number of communities */
    int maxIter, roundIter;
    double lambda;
    double *U = NULL;
    RBR_INT *index = NULL, *labels = NULL, *iU = NULL;

    rbr_result * out = rbr_result_create();

    if (out == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory. Failed to create output.\n", __func__);
        return NULL;
    }

    /* Import local variables */
    RBR_INT nnzero = A->nnz;
    RBR_INT nnode = A->n;
    double *d = A->d;
    const RBR_INT *ir = A->pntr;
    const RBR_INT *jc = A->indx;
    maxIter = param->maxIter;
    roundIter = param->roundIter;
    RBR_INT p = param->p;

    if (param->verbose){
        printf("---------------- Input summary ----------------\n");
        printf("  Nodes: " _fmt_int "           Edges: " _fmt_int "\n", A->n, A->nnz / 2);
        printf("maxIter: %d               p: " _fmt_int "\n", param->maxIter, param->p);
        printf("-----------------------------------------------\n");
        if (param->extract == RBR_ROUNDING){
            printf("Extraction method: rounding\n");
        } else if (param->extract == RBR_KMEANS){
            printf("Extraction method: kmeans\n");
        } else {
            printf("Extraction method: none\n");
        }
        printf("Use full storage of U: %s\n", param->full == 1 ? "yes" : "no");
        printf("-----------------------------------------------\n");
    }
    labels = (RBR_INT *)malloc(nnode * sizeof(RBR_INT));

    /* sum(d) = nnzero, so lambda = 1/nnezero. */
    lambda = 1.0 / nnzero;

    /* timing */
    out->elapsed = omp_get_wtime();

    /* Initialize U 
     *
     * full = 0
     * U is sparse row major, each row with p elements
     * iU indicates which column each element belongs to
     *
     * full = 1
     * U is dense row major, each row with k elements */

    if (param->full){
        U = (double*)calloc(nnode * k, sizeof(double));
    } else {
        U = (double*)calloc(nnode * p, sizeof(double));
        iU = (RBR_INT*)malloc(nnode * p * sizeof(RBR_INT));
    }

    index = (RBR_INT*)malloc(nnode * sizeof(RBR_INT));

    // initialize UTB
    double UTB[k];

    // validate memory
    if (U == NULL || labels == NULL || index == NULL || (param->full == 0 && iU == NULL)){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        goto cleanup_on_error;
    }

    param->shuffle ? shuffling(nnode, index) : shuffling_false(nnode, index);

    /* Initialize obj and UTB; UTB = U'*d; */
    memset(UTB, 0, k * sizeof(double));
    for (RBR_INT i = 0; i < nnode; ++i){
        if (param->full){
            RBR_INT init_class = index[i] % k;
            U[i + nnode * init_class] = 1;
            UTB[init_class] += d[i];
        } else {
            U[i * p] = 1.0;
            iU[i * p] = index[i] % k;
            if(p > 1) iU[i * p + 1] = -1;
            UTB[iU[i * p]] += d[i];
        }
    }


    /* TODO: obj = sum(sum(full((A*U)).*U))-lambda*UTB'*UTB; */

    /* randomly choose idx to avoid false sharing. */

    for (int iter = 0; iter < maxIter; ++iter){
        if (param->verbose){
            printf("Iter: %d (%f s) \n", iter, omp_get_wtime() - out->elapsed);
        }

        if (param->shuffle) shuffling(nnode, index);
#pragma omp parallel
        {
#pragma omp for nowait
            /* j is the # of row that is being optimized. */
            for (RBR_INT j = 0; j < nnode; ++j){
                RBR_INT rid = param->shuffle ? index[j] : j;

                double midb1[k], b[k], xmid[k];
                RBR_INT ixmid[k];

                memset(b, 0, k * sizeof(double));

                if (param->full){
                    /* Sparse m-by-v product: b = -U'*A(:,idxj)
                     * Note: A is symmetric ! csr <=> csc */
                    for (RBR_INT ii = 0; ii < k; ++ii){
                        for (RBR_INT jj = ir[rid]; jj < ir[rid + 1]; ++jj){
                            b[ii] -= U[jc[jj] + ii * nnode];
                        }
                    }

                    /* b += (UTB-U(idx(j),:)'*d(idx(j)))*lambda*d(idx(j)) */
                    cblas_dcopy(k, UTB, 1, midb1, 1);
                    cblas_daxpy(k, -d[rid], U + rid, nnode, midb1, 1);
                    cblas_daxpy(k, lambda * d[rid], midb1, 1, b, 1);

                    /* b is derived, now we solve the subproblem:
                     * minimize b'x, s.t x'x=1, x>=0. */

                    /* [minVal, minPos] = min(b); */
                    double minVal = b[0]; 
                    RBR_INT minPos = 0;    
                    for (RBR_INT ii = 1; ii < k; ++ii) {
                        if (b[ii] < minVal){
                            minVal = b[ii];
                            minPos = ii;
                        }
                    }

                    if (minVal < 0){
                        for (RBR_INT iii = 0; iii < k; iii++)   {
                            if (b[iii] < 0) xmid[iii] = -b[iii];
                            else xmid[iii] = 0;
                        }
                        double norm = 1.0 / cblas_dnrm2(k, xmid, 1);
                        cblas_dscal(k, norm, xmid, 1);
                    }

                    else{
                        memset(xmid, 0, k * sizeof(double));
                        xmid[minPos] = 1;
                    }

                    /* Update UTB and U ; */
                    cblas_daxpy(k, d[rid], xmid, 1, UTB, 1);  /* UTB = UTB + (xmid - U(idx(j), :)')*d(idx(j)); */
                    cblas_daxpy(k, -d[rid], U + rid, nnode, UTB, 1); 

                    cblas_dcopy(k, xmid, 1, U + rid, nnode);  /* U(idx(j), :) = xmid; */

                } else { /* param.full = 0 */

                    /* Sparse - Sparse product: b = -U'*A(:,idxj)
                     * Note: A is symmetric ! csr <=> csc */
                    for (RBR_INT jj = ir[rid]; jj < ir[rid + 1]; ++jj){
                        /* consider row #jc[jj] in U */
                        double *rU = U + jc[jj] * p;
                        RBR_INT *irU = iU + jc[jj] * p;
                        for (RBR_INT ii = 0; ii < p; ++ii){
                            if (irU[ii] == -1) break;
                            b[irU[ii]] -= rU[ii];
                        }
                    }

                    /* b += (UTB-U(idx(j),:)'*d(idx(j)))*lambda*d(idx(j)) */
                    cblas_dcopy(k, UTB, 1, midb1, 1);
                    for (RBR_INT i = 0; i < p; ++i){
                        double *rU = U + rid * p;
                        RBR_INT *irU = iU + rid * p;
                        if (irU[i] == -1) break;
                        midb1[irU[i]] -= d[rid] * rU[i];
                    }
                    cblas_daxpy(k, lambda * d[rid], midb1, 1, b, 1);

                    /* b is derived, now we solve the subproblem
                       minimize b'x
                       s.t x'x=1, x>=0, |x|_0<=p*/
                    RBR_INT iPos;
                    cblas_dcopy(k, b, 1, xmid, 1);
                    solve_sub_U(k, p, xmid, ixmid, &iPos);

                    if (iPos > 0){
                        double norm = 1.0 / cblas_dnrm2(iPos, xmid, 1);
                        cblas_dscal(iPos, -norm, xmid, 1);
                    }

                    /* UTB = UTB + (xmid - U(idx(j), :)')*d(idx(j)); */
                    for (RBR_INT i = 0; i < iPos; ++i)
                        UTB[ixmid[i]] += d[rid] * xmid[i];
                    for (RBR_INT i = 0; i < p; ++i){
                        double *rU = U + rid * p;
                        RBR_INT *irU = iU + rid * p;
                        if (irU[i] == -1) break;
                        UTB[irU[i]] -= d[rid] * rU[i];
                    }

                    /* U(idx(j), :) = xmid; */
                    cblas_dcopy(iPos, xmid, 1, U + rid * p, 1);
                    memcpy(iU + rid * p, ixmid, iPos * sizeof(RBR_INT));
                    if (iPos < p) iU[rid * p + iPos] = -1;
                }

            }
            /* TODO: If necessary, update the objective. */

            if (param->extract == RBR_ROUNDING && iter > roundIter){
                /* Rounding Process */
#pragma omp for nowait
                for (RBR_INT j = 0; j < nnode; ++j){
                    if (param->full){
                        RBR_INT maxPos = imax(k, U + j, nnode);
                        labels[j] = maxPos;
                        for (RBR_INT jj = 0; jj < k; ++jj) U[j + jj * nnode] = 0;
                        U[j + nnode * maxPos] = 1;
                    } else {
                        RBR_INT rid = param->shuffle ? index[j] : j;
                        RBR_INT maxPos = imax_p(p, U + rid * p, iU + rid * p);
                        labels[rid] = iU[rid * p + maxPos];
                        U[rid * p] = 1.0;
                        iU[rid * p] = iU[rid * p + maxPos];
                        if (p > 1) iU[rid * p + 1] = -1;
                    }
                }

            }

        }

        /* due to the rounding process, we have to recalculate the UTB each time. */
        if (param->extract == RBR_ROUNDING && iter > roundIter){
            //cblas_dgemv(CblasColMajor, CblasTrans, nnode, k, 1, U, nnode, d, 1, 0, UTB, 1);
            memset(UTB, 0, k * sizeof(double));

            //#pragma omp parallel for
            for (RBR_INT j = 0; j < nnode; ++j)
                UTB[labels[j]] += d[j];
        }

    }

    out->elapsed = omp_get_wtime() - out->elapsed;
    /* kmeans */
    if (param->extract == RBR_KMEANS){
        double *full_U;
        /* perform one-step rounding to find the initial value */
        if (param->full){
#pragma omp for
            for (RBR_INT j = 0; j < nnode; ++j){
                RBR_INT maxPos = imax(k, U + j, nnode);
                labels[j] = maxPos;
            }

            full_U = U;
        } else {
#pragma omp for
            for (RBR_INT jj = 0; jj < nnode; ++jj){
                RBR_INT maxPos = imax_p(p, U + jj * p, iU + jj * p);
                labels[jj] = iU[jj * p + maxPos];
            }

            full_U = (double*)malloc(nnode * k * sizeof(double));
            
            if (full_U == NULL){
                fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
                goto cleanup_on_error;
            }
            sparse_to_full(nnode, k, p, U, iU, full_U);
        }

        kmeans_param kpar = param->k_param;
        kpar.init = KMEANS_LABELS;
        double *H = (double*)malloc(k * k * sizeof(double));

        if (H == NULL){
            fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
            if (!param->full){
                free(full_U);
            }
            goto cleanup_on_error;
        }
        int kret = kmeans(nnode, k, full_U, k, labels, H, kpar);
        free(H); free(full_U);

        if (kret){
            fprintf(stderr, "[ERROR] %s: kmeans returns nonzero code\n", __func__);
            goto cleanup_on_error;
        }
    }
    if (param->verbose)
        printf("Elapsed: %f sec.\n", out->elapsed);

    //for (int i = 0; i < nnode; ++i) labels[i] = iU[i * p];
    out->labels = labels;
    out->U = U;
    out->iU = iU;

    if (param->extract != RBR_NONE){
        /* re-compute UTB */
        memset(UTB, 0, k * sizeof(double));
        for (RBR_INT j = 0; j < nnode; ++j)
            UTB[labels[j]] += d[j];
        out->Q = cal_Q(A, nnode, k, labels, lambda, UTB);
        out->CC = cal_CC(A, nnode, k, labels);
        out->S  = cal_Strength(A, nnode, k, labels);
    }

    if (0){
cleanup_on_error:
        if (U) free(U);
        if (iU) free(iU);
        if (labels) free(labels);
        out->exit_code = 100;
    }
    free(index);
    return out;

}

void shuffling(RBR_INT n, RBR_INT *card){
    for (RBR_INT i = 0; i < n; ++i) card[i] = i;

    /* Shuffle elements by randomly exchanging each with one other. */
    for (RBR_INT i = 0; i < n - 1; ++i) {
        RBR_INT r = i + rand() % (n - i);
        RBR_INT temp = card[i]; card[i] = card[r]; card[r] = temp;
    }
}

void shuffling_false(RBR_INT n, RBR_INT *card){
    for (RBR_INT i = 0; i < n; ++i) card[i] = i;
}

RBR_INT imax(RBR_INT n, const double *x, RBR_INT incx){
    if (n == 0){
        return -1;
    }

    RBR_INT id, ret = 0;
    double val = x[0];

    for (id = 1; id < n; ++id){
        if (x[id * incx] > val){
            val = x[id * incx];
            ret = id;
        }
    }

    return ret;
}

RBR_INT imax_p(RBR_INT n, const double *x, const RBR_INT *ix){
    RBR_INT id, ret = 0;
    double val = x[0];

    for (id = 1; ix[id] != -1 && id < n; ++id){
        if (x[id] > val){
            val = x[id];
            ret = id;
        }
    }

    return ret;
}

