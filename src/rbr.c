#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "mkl.h"
#include "DCRBR.h"


DCRBR_out rbr(adj_matrix *A, int k, DCRBR_param param) {
    /* k is the number of communities;   This part is Correct */
    int nnode, nnzero, maxIter, roundIter, p;
    double lambda;
    double *U, *d;
    int *index, *labels, *ir, *jc, *iU;
    DCRBR_out out;

    /* Import local variables */
    nnzero = A->nnz; nnode = A->n; d = A->d;
    ir = A->pntr; jc = A->indx;
    maxIter = param.maxIter; roundIter = param.roundIter, p = param.p;
    if (param.verbose){
        fprintf(stderr, "Input summary:\n");
        fprintf(stderr, "Nodes: %d    Edges: %d\n", A->n, A->nnz / 2);
        fprintf(stderr, "maxIter: %d    p: %d\n", param.maxIter, param.p);
        fprintf(stderr, "extraction method: %d\n", (int)param.extract);
    }
    labels = (int*)malloc(nnode * sizeof(int));

    /* sum(d) = nnzero, so lambda = 1/nnezero. */
    lambda = 1.0 / nnzero;

    /* timing */
    out.elapsed = omp_get_wtime();

    /* Initialize U 
       Note: U is sparse row major, each row with p elements
       iU indicates which column each element belongs to*/
    U = (double*)calloc(nnode * p, sizeof(double));
    iU = (int*)malloc(nnode * p * sizeof(int));
    index = (int*)malloc(nnode * sizeof(int));
    param.shuffle ? shuffling(nnode, index) : shuffling_false(nnode, index);

    /* Initialize obj and UTB; UTB = U'*d; */
    double UTB[k];
    memset(UTB, 0, k * sizeof(double));
    for (int i = 0; i < nnode; ++i){
        U[i * p] = 1.0;
        iU[i * p] = index[i] % k;
        if(p > 1) iU[i * p + 1] = -1;
        UTB[iU[i * p]] += d[i];
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
            /* TODO: If necessary, update the objective. */

            if (param.extract == RBR_ROUNDING && iter > roundIter){
#pragma omp for nowait
                /* Rounding Process */
                for (int jj = 0; jj < nnode; ++jj){
                    int rid = param.shuffle ? index[jj] : jj;
                    int maxPos = imax_p(p, U + rid * p, iU + rid * p);
                    labels[rid] = iU[rid * p + maxPos];
                    U[rid * p] = 1.0;
                    iU[rid * p] = iU[rid * p + maxPos];
                    if (p > 1) iU[rid * p + 1] = -1;
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
        /* perform one-step rounding to find the initial value */
#pragma omp for
        /* Rounding Process */
        for (int jj = 0; jj < nnode; ++jj){
            int maxPos = imax_p(p, U + jj * p, iU + jj * p);
            labels[jj] = iU[jj * p + maxPos];
        }
        param.k_param.init = KMEANS_LABELS;
        double *full_U = (double*)malloc(nnode * k * sizeof(double));
        sparse_to_full(nnode, k, p, U, iU, full_U);
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

void DCRBR_out_destroy(DCRBR_out *out){
    free(out->labels);
    free(out->U);
    free(out->iU);
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
