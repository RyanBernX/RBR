/*
 * ===========================================================================
 *
 *       Filename:  rbr_mc.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/29/2023 10:13:10 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BDA, PKU
 *
 * ===========================================================================
 */

#include "RBR.h"
#include "rbr_blas.h"
#include "rbr_subroutines.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

rbr_result * rbr_maxcut(const adj_matrix *A, RBR_INT k, const rbr_param *param) {
    int maxIter, roundIter;
    double *U;
    RBR_INT *index, *labels, *iU = NULL;
    double fun = 0, fun_p = 0, df = 0;

    /* Import local variables */
    RBR_INT nnode = A->n;
    const RBR_INT *ir = A->pntr;
    const RBR_INT *jc = A->indx;
    const double *val = A->val;

    maxIter = param->maxIter;
    roundIter = param->roundIter;
    RBR_INT p = param->p;

    rbr_result *out = rbr_result_create();
    if (out == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory. Failed to create output.\n", __func__);
        return NULL;
    }

    if (param->verbose){
        printf("---------------- Input summary ----------------\n");
        printf("  Nodes: " _fmt_int "           Edges: " _fmt_int "\n", A->n, A->nnz / 2);
        printf("maxIter: %d               p: " _fmt_int "\n", param->maxIter, param->p);
        printf("-----------------------------------------------\n");

        if (param->extract == RBR_ROUNDING){
            printf("Extraction method: rounding\n");
        } else if (param->extract == RBR_NONE){
            printf("Extraction method: none\n");
        }

        printf("Use full storage of U: %s\n", param->full == 1 ? "yes" : "no");
        printf("Shuffle: %s\n", param->shuffle == 1 ? "yes" : "no");
        printf("-----------------------------------------------\n");
    }
    labels = (RBR_INT*)malloc(nnode * sizeof(RBR_INT));

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

    // validate memory
    if (U == NULL || labels == NULL || index == NULL || (param->full == 0 && iU == NULL)){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        goto cleanup_on_error;
    }

    param->shuffle ? shuffling(nnode, index) : shuffling_false(nnode, index);

    /* Initialize U */
    for (RBR_INT i = 0; i < nnode; ++i){
        if (param->full){
            RBR_INT init_class = index[i] % k;
            U[i * k + init_class] = 1;
        } else {
            U[i * p] = 1.0;
            iU[i * p] = index[i] % k;
            if(p > 1) iU[i * p + 1] = -1;
        }
    }

    /* initialize function */
    if (param->full){
        fun = cal_maxcut_value(A, nnode, k, U);
    } else {
        fun = cal_maxcut_value_p(A, nnode, k, p, U, iU);
    }

    /* randomly choose idx to avoid false sharing. */
    int iter;
    for (iter = 0; iter < maxIter; ++iter){
        if (param->verbose && iter % param->funct_chk_freq == 0){
            if (iter == 0){
                printf("Iter: %7d  fval: %13.7e  df: -------------  (%f s)\n", iter, fun, omp_get_wtime() - out->elapsed);
            } else {
                printf("Iter: %7d  fval: %13.7e  df: %13.7e  (%f s)\n", iter, fun, df, omp_get_wtime() - out->elapsed);
            }
        }

        if (iter > 0 && fabs(df) < param->tol) break;

        if (param->shuffle) shuffling(nnode, index);
#pragma omp parallel
        {
#pragma omp for
            /* j is the # of row that is being optimized. */
            for (RBR_INT j = 0; j < nnode; ++j){
                RBR_INT rid = param->shuffle ? index[j] : j;

                double b[k], xmid[k];
                RBR_INT ixmid[k];

                memset(b, 0, k * sizeof(double));

                if (param->full){
                    /* Sparse m-by-v product: b = -U'*A(:,idxj)
                     * Note: A is symmetric ! csr <=> csc */
                    for (RBR_INT jj = ir[rid]; jj < ir[rid + 1]; ++jj){
                        RBR_INT col = jc[jj];
                        if (val == NULL){
                            cblas_daxpy(k, 1, U + col * k, 1, b, 1);
                        } else {
                            cblas_daxpy(k, val[jj], U + col * k, 1, b, 1);
                        }
                    }

                    /* b is derived, now we solve the subproblem:
                     * minimize b'x, s.t x'x=1. */
                    double nrm = cblas_dnrm2(k, b, 1);
                    if (fabs(nrm) > 1e-14){
                        cblas_dscal(k, -1.0 / nrm, b, 1);
                        cblas_dcopy(k, b, 1, U + rid * k, 1);
                    }

                } else { /* param.full = 0 */

                    /* Sparse - Sparse product: b = -U'*A(:,idxj)
                     * Note: A is symmetric ! csr <=> csc */
                    for (RBR_INT jj = ir[rid]; jj < ir[rid + 1]; ++jj){
                        /* consider row #jc[jj] in U */
                        double *rU = U + jc[jj] * p;
                        RBR_INT *irU = iU + jc[jj] * p;
                        for (RBR_INT ii = 0; ii < p; ++ii){
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
                    if (fabs(nrm) > 1e-14){
                        cblas_dscal(p, 1.0 / nrm, xmid, 1);
                        cblas_dcopy(p, xmid, 1, U + rid * p, 1);
                        memcpy(iU + rid * p, ixmid, p * sizeof(RBR_INT));
                    }
                }

            }

        }
        /* TODO: If necessary, update the objective. */
        if (iter % param->funct_chk_freq == 0){
            fun_p = fun;
            if (param->full){
                fun = cal_maxcut_value(A, nnode, k, U);
            } else {
                fun = cal_maxcut_value_p(A, nnode, k, p, U, iU);
            }
            df = fun - fun_p;
        }

    }
    out->elapsed = omp_get_wtime() - out->elapsed;

    /* compute the best cut */
    if (param->extract == RBR_ROUNDING){
        if (param->verbose){
            printf("Perform rounding for %d times...\n", maxIter - roundIter);
        }
        out->Q = rounding_maxcut(A, nnode, k, p, maxIter - roundIter, U, iU, labels);
        out->labels = labels;
    }

    out->U = U;
    out->iU = iU;
    out->iter = iter;
    if (param->full){
        out->funct_V = cal_maxcut_value(A, nnode, k, U);
    } else {
        out->funct_V = cal_maxcut_value_p(A, nnode, k, p, U, iU);
    }

    if (param->verbose){
        printf("Funct value: %13.6e.\n", out->funct_V);
        if (param->extract == RBR_ROUNDING){
            printf("Cut value: %15.8g\n", out->Q);
        }
        printf("Elapsed: %f sec.\n", out->elapsed);
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


/** @fn double cal_obj_value()
 * @brief compute the objective function value. Requires: U column major, ldu=n
 */
double cal_obj_value(const adj_matrix *A, RBR_INT n, RBR_INT k, const double *U, double lambda, const double *UTB){
    double *AU, obj;
    AU = (double*)calloc(n * k, sizeof(double));

    if (AU == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        return NAN;
    }

    const RBR_INT *ir = A->pntr;
    const RBR_INT *jc = A->indx;

    /* compute A * U */
    for (RBR_INT j = 0; j < k; ++j){
#pragma omp parallel for
        for (RBR_INT r = 0; r < n; ++r){
            for (RBR_INT i = ir[r]; i < ir[r + 1]; ++i)
                AU[j * n + r] += U[jc[i] + j * n];
        }
    }
    /* compute <AU, U> */
    obj = cblas_ddot(n * k, AU, 1, U, 1);
    obj -= lambda * cblas_ddot(k, UTB, 1, UTB, 1);

    return obj;
}

double cal_maxcut_value(const adj_matrix *A, RBR_INT n, RBR_INT k, const double *U){
    double *AU, cut;
    AU = (double*)calloc(n * k, sizeof(double));

    const RBR_INT *ir = A->pntr;
    const RBR_INT *jc = A->indx;
    const double *val = A->val;

    if (AU == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        return NAN;
    }

    /* compute A * U */
#pragma omp parallel for
    for (RBR_INT r = 0; r < n; ++r){
        for (RBR_INT i = ir[r]; i < ir[r + 1]; ++i){
            if (val == NULL){
                cblas_daxpy(k, 1.0, U + jc[i] * k, 1, AU + r * k, 1);
            } else {
                cblas_daxpy(k, val[i], U + jc[i] * k, 1, AU + r * k, 1);
            }
        }
    }

    /* compute cut value */
    double one = 1;
    cut = cblas_ddot(n, A->d, 1, &one, 0);
    cut -= cblas_ddot(n * k, AU, 1, U, 1);

    free(AU);
    return cut / 4;
}

/** @fn double cal_maxcut_value_p()
 * @brief Compute the maxcut value. Both A and U are sparse.
 */

double cal_maxcut_value_p(const adj_matrix *A, RBR_INT n, RBR_INT k, RBR_INT p, const double *U, const RBR_INT *iU){
    double *AU, cut;
    /* maybe waste of storage 
       Note: row major index */
    AU = (double*)calloc(n * k, sizeof(double));

    const RBR_INT *ir = A->pntr;
    const RBR_INT *jc = A->indx;
    const double *val = A->val;

    if (AU == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        return NAN;
    }

    /* compute A * U */
#pragma omp parallel for
    for (RBR_INT r = 0; r < n; ++r){ /* for every row of A and U*/
        RBR_INT i_start_AU = r * k;
        for (RBR_INT i = ir[r]; i < ir[r + 1]; ++i){
            RBR_INT i_start = jc[i] * p;
            for (RBR_INT c = 0; c < p; ++c){
                RBR_INT col = iU[i_start + c];
                if (col == -1) break;
                if (val == NULL){
                    AU[i_start_AU + col] += U[i_start + c];
                } else {
                    AU[i_start_AU + col] += val[i] * U[i_start + c];
                }
            }
        }
    }

    /* compute cut value */
    double one = 1;
    cut = cblas_ddot(n, A->d, 1, &one, 0);
    for (RBR_INT r = 0; r < n; ++r){
        RBR_INT i_start = r * p;
        for (RBR_INT c = 0; c < p; ++c){
            RBR_INT col = iU[i_start + c];
            if (col == -1) break;
            cut -= U[i_start + c] * AU[r * k + col];
        }
    }

    free(AU);
    return cut / 4;
}

double cal_maxcut_value_i(const adj_matrix *A, RBR_INT n, RBR_INT *labels){
    double cut = 0;
#pragma omp parallel for shared(A,n,labels) reduction(+:cut)
    for (RBR_INT r = 0; r < n; ++r){
        for (RBR_INT i = A->pntr[r]; i < A->pntr[r+1]; ++i){
            /* r, jc[i] */
            if (labels[r] != labels[A->indx[i]]){
                if (A->val == NULL){
                    ++cut;
                } else {
                    cut += A->val[i];
                }
            }
        }
    }
    return cut / 2;
}

double rounding_maxcut(const adj_matrix *A, RBR_INT n, RBR_INT k, RBR_INT p, int ntries, const double *U, const RBR_INT *iU, RBR_INT *labels){
    double cut = 0, cut_tmp;

    RBR_INT *labels_tmp = (RBR_INT*)malloc(n * sizeof(RBR_INT));
    double *Ur = (double*)malloc(n * sizeof(double));
    double *r = (double*)malloc(k * sizeof(double));

    if (labels_tmp == NULL || Ur == NULL || r == NULL){
        fprintf(stderr, "[ERROR] %s: insufficient memory\n", __func__);
        if (labels_tmp) free(labels_tmp);
        if (Ur) free(Ur);
        if (r) free(r);

        return NAN;
    }

    for (int i = 0; i < ntries; ++i){
        /* random gaussian */
        random_gaussian_vector(k, r);
        double nrm = cblas_dnrm2(k, r, 1);
        /* normalize */
        cblas_dscal(k, 1.0 / nrm, r, 1);
        /* Ur = U * r */
        if (iU == NULL){
            cblas_dgemv(CblasRowMajor, CblasNoTrans, n, k, 1, U, k, r, 1, 0, Ur, 1);
        } else {
            memset(Ur, 0, n * sizeof(double));
            /* sparse U * r */
#pragma omp parallel for
            for (RBR_INT row = 0; row < n; ++row){
                RBR_INT i_start = row * p;
                for (RBR_INT i = 0; i < p; ++i){
                    RBR_INT col = iU[i_start + i];
                    if (col == -1) break;
                    Ur[row] += U[i_start + i] * r[col];
                }
            }
        }
        /* rounding */
        for (RBR_INT j = 0; j < n; ++j){
            labels_tmp[j] = Ur[j] > 0 ? 1 : -1;
        }
        cut_tmp = cal_maxcut_value_i(A, n, labels_tmp);
        /* record the best cut */
        if (i == 0 || cut_tmp > cut){
            cut = cut_tmp;
            memcpy(labels, labels_tmp, n * sizeof(RBR_INT));
        }
    }
    free(labels_tmp); free(Ur); free(r);
    return cut;
}

