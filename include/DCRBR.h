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

#ifndef __DCRBR_H_
#define __DCRBR_H_

#include "spmat.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  KMEANS_SAMPLE,
  KMEANS_CENTER,
  KMEANS_LABELS,
  KMEANS_USER
} KMEANS_INIT;

typedef enum {
  RBR_ROUNDING,
  RBR_KMEANS,
  RBR_NONE
} rbr_extract;

typedef struct {
  int maxit;
  KMEANS_INIT init;
  int verbose;
} kmeans_param;

typedef struct {
  int *labels;
  double *U;
  int *iU;
  int iter;
  double elapsed;
  double funct_V;
  double CC;
  double S;
  double Q;
} DCRBR_out;

typedef struct {
  int verbose;
  int maxIter;
  int shuffle;
  int p;
  int roundIter;
  int full;
  int funct_chk_freq;
  double tol;
  rbr_extract extract;
  kmeans_param k_param;
} DCRBR_param;


/* from contrib */
void dsort(int n, double *data, int *index);

/* c++ plugin (nth_element) */
void solve_sub_U(int n, int p, double *xmid, int *ixmid, int *iPos);
void solve_sub_maxcut(int n, int p, double *xmid, int *ixmid);

/* c++ plugin (random gaussian) */
void random_gaussian_vector(int n, double *x);

void shuffling(int n, int *card);
void shuffling_false(int n, int *card);
int imax(int n, const double *x, int incx);
int imax_p(int n, const double *x, const int *ix);
DCRBR_out DCRBR_rounding(adj_matrix *A, int k, DCRBR_param param);
DCRBR_out DCRBR_kmeans(adj_matrix *A, int k, DCRBR_param param);
DCRBR_out DCRBR_rounding_p(adj_matrix *A, int k, DCRBR_param param);
DCRBR_out DCRBR_kmeans_p(adj_matrix *A, int k, DCRBR_param param);
DCRBR_out rbr(adj_matrix *A, int k, DCRBR_param param);
DCRBR_out rbr_maxcut(adj_matrix *A, int k, DCRBR_param param);
void DCRBR_out_destroy(DCRBR_out *out);
int *generate_labels(int n, int k, const double *U);
double cal_obj_value(adj_matrix *A, int n, int k, const double *U, double lambda, const double *UTB);
double cal_maxcut_value(adj_matrix *A, int n, int k, const double *U);
double cal_maxcut_value_p(adj_matrix *A, int n, int k, int p, const double *U, const int *iU);
double cal_maxcut_value_i(adj_matrix *A, int n, int *labels);
double rounding_maxcut(adj_matrix *A, int n, int k, int p, int ntries, const double *U, const int *iU, int *labels);

/* metric.c */
int cal_mis(int, int, const int*, const int*);
void cal_confusion_matrix(int n, int k, const int *labels, const int *info, int *mat);
double cal_Fsame(int n, int k, const int *labels, const int *info);
double cal_NMI(int n, int k, const int *labels, const int *info);
double cal_Q(adj_matrix *A, int n, int k, const int *labels, double lambda, const double *UTB);
double cal_Strength(adj_matrix *A, int n, int k, const int *labels);
double cal_CC(adj_matrix *A, int n, int k, const int *labels);
int have_edge(adj_matrix *A, int s, int t);

/* kmeans */
void kmeans(int n, int p, const double *X, int K, int *labels, double *H,
        kmeans_param opts);

/* utils.c */
void sparse_to_full(int n, int k, int p, const double *X, const int *iU,
        double *out);
void sparse_to_full_c(int n, int k, int p, const double *X, const int *iU,
        double *out);

#ifdef __cplusplus
}
#endif


#endif
