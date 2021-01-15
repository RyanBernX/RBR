#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "mex.h"
#include "DCRBR.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    adj_matrix A;
    DCRBR_param param;
    DCRBR_out out;

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

    ir = mxGetIr(prhs[0]);
    jc = mxGetJc(prhs[0]);
    pr = mxGetPr(prhs[0]);
    A.d = mxGetPr(prhs[1]);

    int k = (int)(*mxGetPr(prhs[2]));
    param.p = (int)(*mxGetPr(prhs[3]));
    param.maxIter = (int)(*mxGetPr(prhs[4]));
    param.roundIter = param.maxIter - 5;
    switch ((int)(*mxGetPr(prhs[5]))){
        case 0:
            param.extract = RBR_ROUNDING;
            break;
        case 1:
            param.extract = RBR_KMEANS;
            break;
        case 2:
            param.extract = RBR_NONE;
            break;
        default:
            param.extract = RBR_ROUNDING;
            break;
    }
    param.full = (int)(*mxGetPr(prhs[6]));
    int prob = (int)(*mxGetPr(prhs[7]));

    /* import A */
    int n, nnz;
    n = mxGetM(prhs[0]);
    nnz = (int)(jc[n] - jc[0]);

    A.n = n; A.nnz = nnz;
    A.pntr = (int*)malloc((n + 1) * sizeof(int));
    A.indx = (int*)malloc(nnz * sizeof(int));
    A.val = NULL;

    for (int i = 0; i <= n; ++i) A.pntr[i] = jc[i];
    for (int i = 0; i < nnz; ++i) A.indx[i] = ir[i];

    if (prob == 0){ /* community */
        out = rbr(&A, k, param);
    } else if (prob == 1){ /* maxcut */
        A.val = pr;
        out = rbr_maxcut(&A, k, param);
    }


    /* extract output */
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    pr = mxGetPr(plhs[0]);
    for (int i = 0; i < n; ++i) pr[i] = out.labels[i];

    if (nlhs == 2){
        plhs[1] = mxCreateDoubleMatrix(n, k, mxREAL);
        pr = mxGetPr(plhs[1]);
        if (param.full == 0){
            sparse_to_full_c(n, k, param.p, out.U, out.iU, pr);
        } else {
            memcpy(pr, out.U, n * k * sizeof(double));
        }
    }

    DCRBR_out_destroy(&out);

    free(A.pntr); free(A.indx);

    return;
}
