/*
 * ===========================================================================
 *
 *       Filename:  arrabit_blas_lapack.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  10/20/2023 03:00:17 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Haoyang Liu (@liuhy), liuhaoyang@pku.edu.cn
 *   Organization:  BDA, PKU
 *
 * ===========================================================================
 */

#ifndef RBR_BLAS_LAPACK_H
#define RBR_BLAS_LAPACK_H

#include "types.h"

// CBLAS enums
typedef enum CBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102} CBLAS_ORDER;
typedef enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, CblasConjNoTrans=114} CBLAS_TRANSPOSE;
typedef enum CBLAS_UPLO      {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
typedef enum CBLAS_DIAG      {CblasNonUnit=131, CblasUnit=132} CBLAS_DIAG;
typedef enum CBLAS_SIDE      {CblasLeft=141, CblasRight=142} CBLAS_SIDE;
typedef CBLAS_ORDER CBLAS_LAYOUT;

// CBLAS functions
double cblas_dasum (const RBR_INT n, const double *x, const RBR_INT incx);
void cblas_daxpy(const RBR_INT n, const double alpha, const double *x, const RBR_INT incx, double *y, const RBR_INT incy);
void cblas_dcopy(const RBR_INT n, const double *x, const RBR_INT incx, double *y, const RBR_INT incy);
double cblas_ddot(const RBR_INT n, const double *x, const RBR_INT incx, const double *y, const RBR_INT incy);
void cblas_dscal(const RBR_INT N, const double alpha, double *X, const RBR_INT incX);
void cblas_dgemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE trans, const RBR_INT m, const RBR_INT n, const double alpha, const double *a, const RBR_INT lda,  const double  *x, const RBR_INT incx,  const double beta,  double  *y, const RBR_INT incy);
void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_TRANSPOSE TransB, const RBR_INT M, const RBR_INT N, const RBR_INT K,
        const double alpha, const double *A, const RBR_INT lda, const double *B, const RBR_INT ldb, const double beta, double *C, const RBR_INT ldc);
double cblas_dnrm2 (const RBR_INT N, const double *X, const RBR_INT incX);
RBR_INT cblas_idamax(const RBR_INT n, const double *x, const RBR_INT incx);
void cblas_dswap(const RBR_INT n, double *x, const RBR_INT incx, double *y, const RBR_INT incy);

#endif

