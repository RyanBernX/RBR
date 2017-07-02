#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "DCRBR.h"

int read_info(const char *filename, int n, int *info);
void write_labels(const char *filename, int n, const int *labels);

int main(int argc, char **argv){
  adj_matrix A;
  DCRBR_out out;
  DCRBR_param param;
  int n, k, mis_best = -1;
  int *info;
  double tt = 0, obj_best = 0, mis_avg = 0;
  FILE *fp;
  int *label_best;

  if (argc < 5){
    fprintf(stderr, "Usage ./main <csr_matrix> <k> <info> <logfile> <label>\n");
    return 1;
  }

  if (read_adj_matrix_csr(argv[1], &A, 0))
    return 1;


  k = strtol(argv[2], NULL, 10); n = A.n;

  info = (int*)malloc(n * sizeof(int));
  if (read_info(argv[3], n, info)){
    free(info);
    return 1;
  }

  srand((unsigned)time(NULL));
  param.verbose = 0; param.maxIter = 10; param.shuffle = 1;
  param.roundIter = 4; param.p = k;

  fprintf(stderr, "log file is %s\n", argv[4]);

  fp = fopen(argv[4], "w");
  fprintf(fp, "mis\n");
  label_best = (int*)malloc(n * sizeof(int));
  for (int i = 0; i < 300; ++i){
    out = DCRBR_rounding(&A, k, param);
    int mis = cal_mis(n, k, out.labels, info);
    fprintf(fp, "%f\n", 1. * mis / n);
    mis_avg += 1. * mis / n;
    if (mis < mis_best || i == 0){
      mis_best = mis;
      memcpy(label_best, out.labels, n * sizeof(int));
    }
    tt += out.elapsed;
    DCRBR_out_destroy(&out);
  }

  printf("Time: %f    Mis%%: %f\n", tt, mis_avg / 300);
  write_labels(argv[5], n, label_best);
  
  fclose(fp);
  free(info); free(label_best);
  adj_matrix_destroy(&A);
  return 0;
}

int read_info(const char *filename, int n, int *info){
  FILE *fp;

  fp = fopen(filename, "r");

  if (fp == NULL){
    fprintf(stderr, "Cannot open info file %s\n", filename);
    return 1;
  }

  for (int i = 0; i < n; ++i)
    fscanf(fp, "%d", info + i);

  fclose(fp);
  return 0;
}

void write_labels(const char *filename, int n, const int *labels){
  FILE *fp;

  fp = fopen(filename, "w");

  if (fp == NULL){
    fprintf(stderr, "Cannot open file %s\n", filename);
    return;
  }

  for (int i = 0; i < n; ++i)
    fprintf(fp, "%d %d\n", i, labels[i]);

  fclose(fp);
}

