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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#ifdef HAVE_MKL
#include "mkl.h"
#else
#include _CBLAS_HEADER
#endif
#include "DCRBR.h"


DCRBR_out rbr(adj_matrix *A, int k, DCRBR_param param) {
    /* k is the number of communities */
    int nnode, nnzero, maxIter, roundIter, p;
    double lambda;
    double *U, *d;
    int *index, *labels, *ir, *jc, *iU = NULL;
    DCRBR_out out;

    /* Import local variables */
    nnzero = A->nnz; nnode = A->n; d = A->d;
    ir = A->pntr; jc = A->indx;
    maxIter = param.maxIter; roundIter = param.roundIter, p = param.p;
    if (param.verbose){
        fprintf(stderr, "Input summary:\n");
        fprintf(stderr, "Nodes: %d    Edges: %d\n", A->n, A->nnz / 2);
        fprintf(stderr, "maxIter: %d    p: %d\n", param.maxIter, param.p);
        if (param.extract == RBR_ROUNDING){
            fprintf(stderr, "extraction method: rounding\n");
        } else if (param.extract == RBR_KMEANS){
            fprintf(stderr, "extraction method: kmeans\n");
        } else {
            fprintf(stderr, "extraction method: none\n");
        }
        fprintf(stderr, "Use full storage of U: %s\n", param.full == 1 ? "yes" : "no");
    }
    labels = (int*)malloc(nnode * sizeof(int));

    /* sum(d) = nnzero, so lambda = 1/nnezero. */
    lambda = 1.0 / nnzero;

    /* timing */
    out.elapsed = omp_get_wtime();

    /* Initialize U 
     *
     * full = 0
     * U is sparse row major, each row with p elements
     * iU indicates which column each element belongs to
     *
     * full = 1
     * U is dense row major, each row with k elements */

    if (param.full){
        U = (double*)calloc(nnode * k, sizeof(double));
    } else {
        U = (double*)calloc(nnode * p, sizeof(double));
        iU = (int*)malloc(nnode * p * sizeof(int));
    }

    index = (int*)malloc(nnode * sizeof(int));
    param.shuffle ? shuffling(nnode, index) : shuffling_false(nnode, index);

    /* Initialize obj and UTB; UTB = U'*d; */
    double UTB[k];
    memset(UTB, 0, k * sizeof(double));
    for (int i = 0; i < nnode; ++i){
        if (param.full){
            int init_class = index[i] % k;
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
        if (param.verbose)
            fprintf(stderr, "Iter: %d\n", iter);

        if (param.shuffle) shuffling(nnode, index);
#pragma omp parallel
        {
#pragma omp for nowait
            /* j is the # of row that is being optimized. */
            for (int j = 0; j < nnode; ++j){
                int rid = param.shuffle ? index[j] : j;

                double midb1[k], b[k], xmid[k];
                int ixmid[k];

                memset(b, 0, k * sizeof(double));

                if (param.full){
                    /* Sparse m-by-v product: b = -U'*A(:,idxj)
                     * Note: A is symmetric ! csr <=> csc */
                    for (int ii = 0; ii < k; ++ii){
                        for (int jj = ir[rid]; jj < ir[rid + 1]; ++jj){
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
                    double minVal = b[0];  int minPos = 0;    
                    for (int ii = 1; ii < k; ++ii) {
                        if (b[ii] < minVal){
                            minVal = b[ii];
                            minPos = ii;
                        }
                    }

                    if (minVal < 0){
                        for (int iii = 0; iii < k; iii++)   {
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
                    for (int jj = ir[rid]; jj < ir[rid + 1]; ++jj){
                        /* consider row #jc[jj] in U */
                        double *rU = U + jc[jj] * p;
                        int *irU = iU + jc[jj] * p;
                        for (int ii = 0; ii < p; ++ii){
                            if (irU[ii] == -1) break;
                            b[irU[ii]] -= rU[ii];
                        }
                    }

                    /* b += (UTB-U(idx(j),:)'*d(idx(j)))*lambda*d(idx(j)) */
                    cblas_dcopy(k, UTB, 1, midb1, 1);
                    for (int i = 0; i < p; ++i){
                        double *rU = U + rid * p;
                        int *irU = iU + rid * p;
                        if (irU[i] == -1) break;
                        midb1[irU[i]] -= d[rid] * rU[i];
                    }
                    cblas_daxpy(k, lambda * d[rid], midb1, 1, b, 1);

                    /* b is derived, now we solve the subproblem
                       minimize b'x
                       s.t x'x=1, x>=0, |x|_0<=p*/
                    int iPos;
                    cblas_dcopy(k, b, 1, xmid, 1);
                    solve_sub_U(k, p, xmid, ixmid, &iPos);

                    if (iPos > 0){
                        double norm = 1.0 / cblas_dnrm2(iPos, xmid, 1);
                        cblas_dscal(iPos, -norm, xmid, 1);
                    }

                    /* UTB = UTB + (xmid - U(idx(j), :)')*d(idx(j)); */
                    for (int i = 0; i < iPos; ++i)
                        UTB[ixmid[i]] += d[rid] * xmid[i];
                    for (int i = 0; i < p; ++i){
                        double *rU = U + rid * p;
                        int *irU = iU + rid * p;
                        if (irU[i] == -1) break;
                        UTB[irU[i]] -= d[rid] * rU[i];
                    }

                    /* U(idx(j), :) = xmid; */
                    cblas_dcopy(iPos, xmid, 1, U + rid * p, 1);
                    memcpy(iU + rid * p, ixmid, iPos * sizeof(int));
                    if (iPos < p) iU[rid * p + iPos] = -1;
                }

            }
            /* TODO: If necessary, update the objective. */

            if (param.extract == RBR_ROUNDING && iter > roundIter){
                /* Rounding Process */
#pragma omp for nowait
                for (int j = 0; j < nnode; ++j){
                    if (param.full){
                        int maxPos = imax(k, U + j, nnode);
                        labels[j] = maxPos;
                        for (int jj = 0; jj < k; ++jj) U[j + jj * nnode] = 0;
                        U[j + nnode * maxPos] = 1;
                    } else {
                        int rid = param.shuffle ? index[j] : j;
                        int maxPos = imax_p(p, U + rid * p, iU + rid * p);
                        labels[rid] = iU[rid * p + maxPos];
                        U[rid * p] = 1.0;
                        iU[rid * p] = iU[rid * p + maxPos];
                        if (p > 1) iU[rid * p + 1] = -1;
                    }
                }

            }

        }

        /* due to the rounding process, we have to recalculate the UTB each time. */
        if (param.extract == RBR_ROUNDING && iter > roundIter){
            //cblas_dgemv(CblasColMajor, CblasTrans, nnode, k, 1, U, nnode, d, 1, 0, UTB, 1);
            memset(UTB, 0, k * sizeof(double));

            //#pragma omp parallel for
            for (int j = 0; j < nnode; ++j)
                UTB[labels[j]] += d[j];
        }

    }
    out.elapsed = omp_get_wtime() - out.elapsed;
    /* kmeans */
    if (param.extract == RBR_KMEANS){
        double *full_U;
        /* perform one-step rounding to find the initial value */
        if (param.full){
#pragma omp for
            for (int j = 0; j < nnode; ++j){
                int maxPos = imax(k, U + j, nnode);
                labels[j] = maxPos;
            }

            full_U = U;
        } else {
#pragma omp for
            for (int jj = 0; jj < nnode; ++jj){
                int maxPos = imax_p(p, U + jj * p, iU + jj * p);
                labels[jj] = iU[jj * p + maxPos];
            }

            full_U = (double*)malloc(nnode * k * sizeof(double));
            sparse_to_full(nnode, k, p, U, iU, full_U);
        }

        param.k_param.init = KMEANS_LABELS;
        double *H = (double*)malloc(k * k * sizeof(double));
        kmeans(nnode, k, full_U, k, labels, H, param.k_param);
        free(H); free(full_U);
    }
    if (param.verbose)
        fprintf(stderr, "Elapsed: %f sec.\n", out.elapsed);

    //for (int i = 0; i < nnode; ++i) labels[i] = iU[i * p];
    out.labels = labels;
    out.U = U; out.iU = iU;
    if (param.extract != RBR_NONE){
        /* re-compute UTB */
        memset(UTB, 0, k * sizeof(double));
        for (int j = 0; j < nnode; ++j)
            UTB[labels[j]] += d[j];
        out.Q = cal_Q(A, nnode, k, labels, lambda, UTB);
        out.CC = cal_CC(A, nnode, k, labels);
        out.S  = cal_Strength(A, nnode, k, labels);
    }

    free(index);
    return out;

}

void shuffling(int n, int *card){
    for (int i = 0; i < n; ++i) card[i] = i;

    /* Shuffle elements by randomly exchanging each with one other. */
    for (int i = 0; i < n - 1; ++i) {
        int r = i + rand() % (n - i);
        int temp = card[i]; card[i] = card[r]; card[r] = temp;
    }
}

void shuffling_false(int n, int *card){
    for (int i = 0; i < n; ++i) card[i] = i;
}

int imax(int n, const double *x, int incx){
    int id, ret = 0;
    double val = x[0];

    for (id = 1; id < n; ++id)
        if (x[id * incx] > val){
            val = x[id * incx];
            ret = id;
        }

    return ret;
}

int imax_p(int n, const double *x, const int *ix){
    int id, ret = 0;
    double val = x[0];

    for (id = 1; ix[id] != -1 && id < n; ++id)
        if (x[id] > val){
            val = x[id];
            ret = id;
        }

    return ret;
}

int* generate_labels(int nnode, int k, const double *U){
    int *ret = (int*)malloc(nnode * sizeof(int));
    int idx;

    for (int i = 0; i < nnode; ++i){
        /* find the first 1 in row i */
        for (idx = 0; idx < k; ++idx)
            if ( (int)U[idx * nnode + i] == 1) break;

        ret[i] = idx;
    }

    return ret;
}

DCRBR_out rbr_maxcut(adj_matrix *A, int k, DCRBR_param param) {
    int nnode, nnzero, maxIter, roundIter, p;
    double *U, *d, *val;
    int *index, *labels, *ir, *jc, *iU = NULL;
    DCRBR_out out;

    /* Import local variables */
    nnzero = A->nnz; nnode = A->n; d = A->d;
    ir = A->pntr; jc = A->indx; val = A->val;
    maxIter = param.maxIter; roundIter = param.roundIter, p = param.p;
    if (param.verbose){
        fprintf(stderr, "Input summary:\n");
        fprintf(stderr, "Nodes: %d    Edges: %d\n", A->n, A->nnz / 2);
        fprintf(stderr, "maxIter: %d    p: %d\n", param.maxIter, param.p);
        if (param.extract == RBR_ROUNDING){
            fprintf(stderr, "extraction method: rounding\n");
        } else if (param.extract == RBR_KMEANS){
            fprintf(stderr, "extraction method: kmeans\n");
        } else {
            fprintf(stderr, "extraction method: none\n");
        }
        fprintf(stderr, "Use full storage of U: %s\n", param.full == 1 ? "yes" : "no");
        fprintf(stderr, "Shuffle: %s\n", param.shuffle == 1 ? "yes" : "no");
    }
    labels = (int*)malloc(nnode * sizeof(int));

    /* timing */
    out.elapsed = omp_get_wtime();

    /* Initialize U 
     *
     * full = 0
     * U is sparse row major, each row with p elements
     * iU indicates which column each element belongs to
     *
     * full = 1
     * U is dense row major, each row with k elements */

    if (param.full){
        U = (double*)calloc(nnode * k, sizeof(double));
    } else {
        U = (double*)calloc(nnode * p, sizeof(double));
        iU = (int*)malloc(nnode * p * sizeof(int));
    }

    index = (int*)malloc(nnode * sizeof(int));
    param.shuffle ? shuffling(nnode, index) : shuffling_false(nnode, index);

    /* Initialize U */
    for (int i = 0; i < nnode; ++i){
        if (param.full){
            int init_class = index[i] % k;
            U[i + nnode * init_class] = 1;
        } else {
            U[i * p] = 1.0;
            iU[i * p] = index[i] % k;
            if(p > 1) iU[i * p + 1] = -1;
        }
    }


    /* randomly choose idx to avoid false sharing. */

    for (int iter = 0; iter < maxIter; ++iter){
        if (param.verbose)
            fprintf(stderr, "Iter: %d\n", iter);

        if (param.shuffle) shuffling(nnode, index);
#pragma omp parallel
        {
#pragma omp for nowait
            /* j is the # of row that is being optimized. */
            for (int j = 0; j < nnode; ++j){
                int rid = param.shuffle ? index[j] : j;

                double b[k], xmid[k];
                int ixmid[k];

                memset(b, 0, k * sizeof(double));

                if (param.full){
                    /* Sparse m-by-v product: b = -U'*A(:,idxj)
                     * Note: A is symmetric ! csr <=> csc */
                    for (int ii = 0; ii < k; ++ii){
                        for (int jj = ir[rid]; jj < ir[rid + 1]; ++jj){
                            int col = jc[jj];
                            if (val == NULL){
                                b[ii] -= U[col + ii * nnode];
                            } else {
                                b[ii] -= val[jj] * U[col + ii * nnode];
                            }
                        }
                    }

                    /* b is derived, now we solve the subproblem:
                     * minimize b'x, s.t x'x=1. */
                    double nrm = cblas_dnrm2(k, b, 1);
                    if (fabs(nrm > 1e-14)){
                        cblas_dscal(k, 1.0 / nrm, b, 1);
                        cblas_dcopy(k, b, 1, U + rid, nnode);
                    }


                } else { /* param.full = 0 */

                    /* Sparse - Sparse product: b = -U'*A(:,idxj)
                     * Note: A is symmetric ! csr <=> csc */
                    for (int jj = ir[rid]; jj < ir[rid + 1]; ++jj){
                        /* consider row #jc[jj] in U */
                        double *rU = U + jc[jj] * p;
                        int *irU = iU + jc[jj] * p;
                        for (int ii = 0; ii < p; ++ii){
                            if (irU[ii] == -1) break;
                            if (val == NULL){
                                b[irU[ii]] -= rU[ii];
                            } else {
                                b[irU[ii]] -= val[jj] * rU[ii];
                            }
                        }
                    }

                    /* solve the subproblem */
                    cblas_dcopy(k, b, 1, xmid, 1);
                    solve_sub_maxcut(k, p, xmid, ixmid);
                    double nrm = cblas_dnrm2(p, xmid, 1);
                    if (fabs(nrm > 1e-14)){
                        cblas_dscal(p, 1.0 / nrm, xmid, 1);
                        cblas_dcopy(p, xmid, 1, U + rid * p, 1);
                        memcpy(iU + rid * p, ixmid, p * sizeof(int));
                    }
                }

            }
            /* TODO: If necessary, update the objective. */

        }

    }
    out.elapsed = omp_get_wtime() - out.elapsed;


    //for (int i = 0; i < nnode; ++i) labels[i] = iU[i * p];
    out.labels = labels;
    out.U = U; out.iU = iU;
    if (param.full){
        out.funct_V = cal_maxcut_value(A, nnode, k, U);
    } else {
        double *U_tmp = (double*)malloc(nnode * k * sizeof(double));
        sparse_to_full_c(nnode, k, p, U, iU, U_tmp);
        out.funct_V = cal_maxcut_value(A, nnode, k, U_tmp);
        free (U_tmp);
    }


    if (param.verbose){
        fprintf(stderr, "Cut value: %e.\n", out.funct_V);
        fprintf(stderr, "Elapsed: %f sec.\n", out.elapsed);
    }

    free(index);
    return out;

}


void DCRBR_out_destroy(DCRBR_out *out){
    free(out->labels);
    free(out->U);
    if (out->iU != NULL) free(out->iU);
}

/** @fn double cal_obj_value()
 * @brief compute the objective function value. Requires: U column major, ldu=n
 */
double cal_obj_value(adj_matrix *A, int n, int k, const double *U, double lambda, const double *UTB){
    double *AU, obj;
    int *ir, *jc;
    AU = (double*)calloc(n * k, sizeof(double));

    ir = A->pntr; jc = A->indx;

    /* compute A * U */
    for (int j = 0; j < k; ++j){
#pragma omp parallel for
        for (int r = 0; r < n; ++r){
            for (int i = ir[r]; i < ir[r + 1]; ++i)
                AU[j * n + r] += U[jc[i] + j * n];
        }
    }
    /* compute <AU, U> */
    obj = cblas_ddot(n * k, AU, 1, U, 1);
    obj -= lambda * cblas_ddot(k, UTB, 1, UTB, 1);

    return obj;
}

double cal_maxcut_value(adj_matrix *A, int n, int k, const double *U){
    double *AU, *val, cut;
    int *ir, *jc;
    AU = (double*)calloc(n * k, sizeof(double));

    ir = A->pntr; jc = A->indx; val = A->val;

    /* compute A * U */
    for (int j = 0; j < k; ++j){
        for (int r = 0; r < n; ++r){
            for (int i = ir[r]; i < ir[r + 1]; ++i){
                if (val == NULL){
                    AU[j * n + r] += U[jc[i] + j * n];
                } else {
                    AU[j * n + r] += val[i] * U[jc[i] + j * n];
                }
            }
        }
    }

    /* compute cut value */
    double one = 1;
    cut = cblas_ddot(n, A->d, 1, &one, 0);
    cut -= cblas_ddot(n * k, AU, 1, U, 1);

    return cut / 4;
}

