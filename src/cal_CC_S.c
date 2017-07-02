#include "mex.h"
#include "DCRBR.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  int n, *ilabels;
  double *labels, *ret, k;
  mwIndex *ir, *jc;
  adj_matrix A;

  if (nrhs != 3){
    mexErrMsgTxt("Usage: ret = cal_CC_S(A, k, labels)");
    return;
  }

  n = mxGetM(prhs[0]);
  if (n != mxGetM(prhs[2])){
    mexErrMsgTxt("Dimension must agree.");
    return;
  }

  ir = mxGetIr(prhs[0]);
  jc = mxGetJc(prhs[0]);

  /* import options */
  A.n = n; A.nnz = jc[n] - jc[0];
  A.d = (double*)malloc(n * sizeof(double));
  A.pntr = (int*)malloc((n + 1) * sizeof(int));
  A.indx = (int*)malloc(A.nnz * sizeof(int));
  ilabels = (int*)malloc(n * sizeof(int));
  for (int i = 0; i < n; ++i){
    A.d[i] = jc[i+1] - jc[i];
    A.pntr[i] = jc[i];
  }
  A.pntr[n] = jc[n];
  for (int i = 0; i < A.nnz; ++i)
    A.indx[i] = ir[i];

  k = *mxGetPr(prhs[1]);
  labels = mxGetPr(prhs[2]);
  for (int i = 0; i < n; ++i)
    ilabels[i] = (int)labels[i];

  //printf("%d  %d  %ld  %ld\n", A.n, A.nnz, jc[n], jc[0]);
  //printf("%d  %d\n", A.pntr[0], A.indx[0]);

  plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
  ret = mxGetPr(plhs[0]);

  /* call C function */
  ret[0] = cal_Strength(&A, n, (int)k, ilabels);
  ret[1] = cal_CC(&A, n, (int)k, ilabels);

  free(A.pntr); free(A.indx); free(A.d);
  free(ilabels);
}



double cal_Strength(adj_matrix *A, int n, int k, const int *labels){
  double *scores, ret = 0;
  int *deg, *pntr = A->pntr, *indx = A->indx;

  scores = (double*)malloc(k * sizeof(double));
  deg = (int*)calloc(k, sizeof(int));

  for (int i = 0; i < k; ++i) scores[i] = 1;

  for (int i = 0; i < n; ++i){ /* for each node */
    int deg_in = 0, deg_out = 0;
    for (int j = pntr[i]; j < pntr[i+1]; ++j){
      int t = indx[j]; /* destnation */
      labels[i] == labels[t] ? ++deg_in : ++deg_out;
    }
    if (deg_in <= deg_out) scores[labels[i]] = 0.5;

    deg[labels[i]] += deg_in - deg_out;
  }

  for (int i = 0; i < k; ++i){
    if (deg[i] <= 0) scores[i] = 0;
    ret += scores[i];
  }

  free(scores); free(deg);
  return ret / k;
}

double cal_CC(adj_matrix *A, int n, int k, const int *labels){
  double *scores, ret = 0;
  int *comm_cnt;

  scores = (double*)calloc(k, sizeof(double));
  comm_cnt = (int*)calloc(k, sizeof(int));

  for (int v = 0; v < n; ++v){ /* for each node */
    int label = labels[v], deg = A->d[v], CC_cnt = 0;
    int *s = A->indx + A->pntr[v];
    ++comm_cnt[label];

    if (deg <= 1){
      scores[label] += 1;
      continue;
    }

    int dom = deg * (deg - 1) / 2;
    /* enumerate all neighbours */
    for (int i = 0; i < deg; ++i)
      for (int  j = i + 1; j < deg; ++j)
        if (labels[s[i]] == label && labels[s[j]] == label && have_edge(A, s[i], s[j]))
          ++CC_cnt;

    scores[label] += 1.0 * CC_cnt / dom;
  }

  for (int i = 0; i < k; ++i){
    scores[i] = comm_cnt[i] > 0 ? scores[i] / comm_cnt[i] : 1;
    ret += scores[i];
  }

  free(scores); free(comm_cnt);
  return ret / k;
}

int have_edge(adj_matrix *A, int s, int t){
  int *indx = A->indx, *pntr = A->pntr;

  for (int j = pntr[s]; j < pntr[s+1]; ++j){
    /* TODO: use binary search if it is too slow */
    if (indx[j] == t) return 1;
  }

  return 0;
}

