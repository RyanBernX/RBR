#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "DCRBR.h"

int main(int argc, char **argv){
  adj_matrix A;
  DCRBR_out out;
  DCRBR_param param;
  int offset = 5;
  int k;
  int threads;
  double tt = 0, obj_best = 0;

  if (argc < 4){
    fprintf(stderr, "Usage ./main <csr_matrix> <k> <threads>\n");
    return 1;
  }

  if (read_adj_matrix_csr(argv[1], &A, 0))
    return 1;


  k = strtol(argv[2], NULL, 10);
  threads = strtol(argv[3], NULL, 10);

  omp_set_num_threads(threads);

  srand((unsigned)time(NULL));
  param.verbose = 1; param.maxIter = 75; param.shuffle = 1;
  param.roundIter = param.maxIter - offset;
  for (int i = 0; i < 1; ++i){
    out = DCRBR_rounding(&A, k, param);
    if (i == 0 || out.Q > obj_best){
      obj_best = out.Q;
    }
    tt += out.elapsed;
    DCRBR_out_destroy(&out);
  }

  printf("Time: %f    CC: %f    S: %f    Q: %f\n", tt, out.CC, out.S, obj_best);

  adj_matrix_destroy(&A);
  return 0;
}
