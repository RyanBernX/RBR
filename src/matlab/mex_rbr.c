#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "RBR.h"
#include "rbr_subroutines.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    adj_matrix *A;
    rbr_param * param;

    if (nlhs > 2){
        mexErrMsgTxt("Too many output parameters.\n");
        return;
    }

    if (nrhs != 8){
        mexErrMsgTxt("Usage: [labels, U] = mex_rbr(A, d, k, p, maxit, extract, full, prob)\n");
        return;
    }

    mwIndex *ir, *jc;
    double *pr;

    mwSize m = mxGetM(prhs[0]);
    mwSize n = mxGetN(prhs[0]);

    if (m != n){
        mexErrMsgTxt("Matrix must be square.");
        return;
    }

    if (!mxIsSparse(prhs[0])){
        mexErrMsgTxt("Matrix must be sparse.");
        return;
    }

    param = rbr_param_create();

    ir = mxGetIr(prhs[0]);
    jc = mxGetJc(prhs[0]);
    pr = mxGetPr(prhs[0]);

    mwSize nnz = (mwSize)(jc[n] - jc[0]);
    A = adj_matrix_create(n, nnz);

    if (param == NULL || A == NULL){
        if (param) rbr_param_destroy(param);
        if (A) adj_matrix_destroy(A);
        mexErrMsgTxt("RBR internal error. Failed to create RBR workspace.");
        return;
    }


    RBR_INT k = (RBR_INT)(*mxGetPr(prhs[2]));
    param->p = (RBR_INT)(*mxGetPr(prhs[3]));
    param->maxIter = (int)(*mxGetPr(prhs[4]));
    param->roundIter = param->maxIter - 5;
    switch ((int)(*mxGetPr(prhs[5]))){
        case 0:
            param->extract = RBR_ROUNDING;
            break;
        case 1:
            param->extract = RBR_KMEANS;
            break;
        case 2:
            param->extract = RBR_NONE;
            break;
        default:
            param->extract = RBR_ROUNDING;
            break;
    }
    param->full = (int)(*mxGetPr(prhs[6]));
    int prob = (int)(*mxGetPr(prhs[7]));

    /* import A */
    // copy d to A->d
    memcpy(A->d, mxGetPr(prhs[1]), n);

    for (mwIndex i = 0; i <= n; ++i){
        A->pntr[i] = jc[i];
    }

    for (mwIndex i = 0; i < nnz; ++i){
        A->indx[i] = ir[i];
    }

    rbr_result * out = NULL;
    if (prob == 0){ /* community */
        out = rbr(A, k, param);
    } else if (prob == 1){ /* maxcut */
        memcpy(A->val, pr, nnz);
        out = rbr_maxcut(A, k, param);
    }

    if (out == NULL || out->exit_code == 100){
        adj_matrix_destroy(A);
        rbr_param_destroy(param);
        rbr_result_destroy(out);
        mexErrMsgTxt("RBR internal error. Solver returned nonzero code.");
        return;
    }

    /* extract output */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    pr = mxGetPr(plhs[0]);
    for (mwIndex i = 0; i < n; ++i){
        pr[i] = out->labels[i];
    }

    if (nlhs == 2){
        plhs[1] = mxCreateDoubleMatrix(n, k, mxREAL);
        pr = mxGetPr(plhs[1]);
        if (param->full == 0){
            sparse_to_full_c(n, k, param->p, out->U, out->iU, pr);
        } else {
            memcpy(pr, out->U, n * k * sizeof(double));
        }
    }

    adj_matrix_destroy(A);
    rbr_param_destroy(param);
    rbr_result_destroy(out);

    return;
}
