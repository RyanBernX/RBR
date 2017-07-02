#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int m, n;
  double *U, *id;

  if (nrhs != 1){
    mexErrMsgTxt("Usage: id = mexFunction(U)");
    return;
  }

  m = mxGetM(prhs[0]); n = mxGetN(prhs[0]);
  U = mxGetPr(prhs[0]);

  plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
  id = mxGetPr(plhs[0]);

  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      if (U[j * m + i] > 0){
        id[i] = (double)j;
        break;
      }

  return;
}

