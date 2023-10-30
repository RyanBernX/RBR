/*
 * ===========================================================================
 *
 *       Filename:  rbr_subroutines.h
 *
 *    Description:  decls for subroutines
 *
 *        Version:  1.0
 *        Created:  10/27/2023 11:06:38 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BDA, PKU
 *
 * ===========================================================================
 */

#ifndef RBR_RBR_SUBROUTINES_H
#define RBR_RBR_SUBROUTINES_H

#include "adj_matrix.h"
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* c++ plugin (nth_element) */
int solve_sub_U(RBR_INT n, RBR_INT p, double *xmid, RBR_INT *ixmid, RBR_INT *iPos);
int solve_sub_maxcut(RBR_INT n, RBR_INT p, double *xmid, RBR_INT *ixmid);

/* c++ plugin (random gaussian) */
void random_gaussian_vector(RBR_INT n, double *x);

#ifdef __cplusplus
}
#endif

// defined in rbr.c
void shuffling(RBR_INT n, RBR_INT *card);
void shuffling_false(RBR_INT n, RBR_INT *card);
RBR_INT imax(RBR_INT n, const double *x, RBR_INT incx);
RBR_INT imax_p(RBR_INT n, const double *x, const RBR_INT *ix);

// defined in rbr_mc.c
double cal_obj_value(const adj_matrix *A, RBR_INT n, RBR_INT k, const double *U, double lambda, const double *UTB);
double cal_maxcut_value(const adj_matrix *A, RBR_INT n, RBR_INT k, const double *U);
double cal_maxcut_value_p(const adj_matrix *A, RBR_INT n, RBR_INT k, RBR_INT p, const double *U, const RBR_INT *iU);
double cal_maxcut_value_i(const adj_matrix *A, RBR_INT n, RBR_INT *labels);
double rounding_maxcut(const adj_matrix *A, RBR_INT n, RBR_INT k, RBR_INT p, int ntries, const double *U, const RBR_INT *iU, RBR_INT *labels);

/* metric.c */
RBR_INT cal_mis(RBR_INT n, RBR_INT k, const RBR_INT *labels, const RBR_INT *info);
void cal_confusion_matrix(RBR_INT n, RBR_INT k, const RBR_INT *labels, const RBR_INT *info, RBR_INT *mat);
double cal_Fsame(RBR_INT n, RBR_INT k, const RBR_INT *labels, const RBR_INT *info);
double cal_NMI(RBR_INT n, RBR_INT k, const RBR_INT *labels, const RBR_INT *info);
double cal_Q(const adj_matrix *A, RBR_INT n, RBR_INT k, const RBR_INT *labels, double lambda, const double *UTB);
double cal_Strength(const adj_matrix *A, RBR_INT n, RBR_INT k, const RBR_INT *labels);
double cal_CC(const adj_matrix *A, RBR_INT n, RBR_INT k, const RBR_INT *labels);
int have_edge(const adj_matrix *A, RBR_INT s, RBR_INT t);

/* kmeans */
int kmeans(RBR_INT n, RBR_INT p, const double *X, RBR_INT K, RBR_INT *labels, double *H,
        kmeans_param opts);

/* utils.c */
void sparse_to_full(RBR_INT n, RBR_INT k, RBR_INT p, const double *X, const RBR_INT *iU,
        double *out);
void sparse_to_full_c(RBR_INT n, RBR_INT k, RBR_INT p, const double *X, const RBR_INT *iU,
        double *out);

#endif

