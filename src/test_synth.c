#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "DCRBR.h"

int main(int argc, char **argv){
  adj_matrix A;
  DCRBR_out out;
  DCRBR_param param;
  int n, k, mis_best = -1;
  int *info;
  double tt = 0, obj_best = 0;

  if (argc < 3){
    fprintf(stderr, "Usage ./main <csr_matrix> <k>\n");
    return 1;
  }

  if (read_adj_matrix_csr(argv[1], &A, 0))
    return 1;

  k = strtol(argv[2], NULL, 10); n = A.n;
  info = (int*)malloc(n * sizeof(int));
  for (int i = 0; i < n; ++i) info[i] = i / (n / k);
  omp_set_num_threads(1);

  srand((unsigned)time(NULL));
  param.verbose = 0; param.maxIter = 15; param.shuffle = 1;
  for (int i = 0; i < 10; ++i){
    out = DCRBR_rounding(&A, k, param);
    int mis = cal_mis(n, k, out.labels, info);
    if (i == 0 || out.funct_V > obj_best){
      mis_best = mis;
      obj_best = out.funct_V;
    }
    tt += out.elapsed;
    DCRBR_out_destroy(&out);
  }

  printf("Time: %f    Mis%%: %f    Obj: %e\n", tt, 1.0 * mis_best / n, obj_best);

  free(info);
  adj_matrix_destroy(&A);
  return 0;
}