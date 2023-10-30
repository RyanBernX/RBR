/*
 * ===========================================================================
 *
 *       Filename:  types.h
 *
 *    Description:  struct and typdefs for RBR
 *
 *        Version:  1.0
 *        Created:  10/27/2023 11:02:56 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BDA, PKU
 *
 * ===========================================================================
 */

#ifndef RBR_TYPES_H
#define RBR_TYPES_H

#include <stdint.h>
#include <inttypes.h>

#ifdef RBR_ILP64
#define RBR_INT int64_t
#define _fmt_int "%" PRId64
#else
#define RBR_INT int32_t
#define _fmt_int "%" PRId32
#endif

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
    int exit_code;
    RBR_INT *labels;
    double *U;
    RBR_INT *iU;
    int iter;
    double elapsed;
    double funct_V;
    double CC;
    double S;
    double Q;
} rbr_result;

typedef struct {
    int verbose;
    int maxIter;
    int shuffle;
    RBR_INT p;
    int roundIter;
    int full;
    int funct_chk_freq;
    double tol;
    rbr_extract extract;
    kmeans_param k_param;
} rbr_param;

#ifdef __cplusplus
}
#endif


#endif

