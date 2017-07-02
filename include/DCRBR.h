#ifndef __DCRBR_H_
#define __DCRBR_H_

#include "spmat.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  int *labels;
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
} DCRBR_param;

/* from contrib */
void dsort(int n, double *data, int *index);

void shuffling(int n, int *card);
void shuffling_false(int n, int *card);
int imax(int n, const double *x, int incx);
int imax_p(int n, const double *x, const int *ix);
DCRBR_out DCRBR_rounding(adj_matrix *A, int k, DCRBR_param param);
DCRBR_out DCRBR_rounding_p(adj_matrix *A, int k, DCRBR_param param);
void DCRBR_out_destroy(DCRBR_out *out);
int *generate_labels(int n, int k, const double *U);
int cal_mis(int, int, const int*, const int*);
double cal_obj_value(adj_matrix *A, int n, int k, const double *U, double lambda, const double *UTB);
double cal_Q(adj_matrix *A, int n, int k, const int *labels, double lambda, const double *UTB);
double cal_Strength(adj_matrix *A, int n, int k, const int *labels);
double cal_CC(adj_matrix *A, int n, int k, const int *labels);
int have_edge(adj_matrix *A, int s, int t);

#ifdef __cplusplus
}
#endif


#endif
