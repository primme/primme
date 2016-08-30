/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2015 College of William & Mary,
 *   James R. McCombs, Eloy Romero Alcalde, Andreas Stathopoulos, Lingfei Wu
 *
 *   This file is part of PRIMME.
 *
 *   PRIMME is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   PRIMME is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *******************************************************************************
 * File: numerical.h
 *
 * Purpose - Contains prototypes for fundamental numerical functions.
 *
 ******************************************************************************/

#ifndef NUMERICAL_H
#define NUMERICAL_H

#include <limits.h>    
#include "primme.h"

// TEMP: Interface @(...) and macros

#ifndef USE_DOUBLECOMPLEX
#define USE_DOUBLECOMPLEX
#endif

#ifdef USE_DOUBLECOMPLEX
#  ifndef __cplusplus
#     include <complex.h> /* definition of creal, cabs, conj */
#     define SCALAR complex double
#     define REAL double
#     define REAL_PART(x) (creal(x))
#     define ABS(x) (cabs(x))
#     define CONJ(x) (conj(x))
#  else
#     include <complex> /* definition of real, abs, conj */
#     define SCALAR std::complex<double>
#     define REAL double
#     define REAL_PART(x) (std::real(x))
#     define ABS(x) (std::abs(x))
#     define CONJ(x) (std::conj(x))
#  endif
#else
#  define SCALAR double
#  define REAL double
#  define REAL_PART(x) (x)
#  define ABS(x) (fabs(x))
#  define CONJ(x) (x)
#endif

/* complex.h may be defined in primme.h or here; so undefine I */
#ifdef I
#   undef I
#endif

#ifndef __cplusplus
#include <tgmath.h>   /* select proper function abs from fabs, cabs... */
#endif

#define MACHINE_EPSILON 1.11e-16

#ifdef F77UNDERSCORE
#define FORTRAN_FUNCTION(X) X ## _
#else
#define FORTRAN_FUNCTION(X) X
#endif

#ifndef max
#  define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#  define min(a, b) ((a) < (b) ? (a) : (b))
#endif


#ifdef __cplusplus
extern "C" {
#endif

int Num_zhpev_zprimme(int iopt, __PRIMME_COMPLEX_DOUBLE__ *ap, double *w, __PRIMME_COMPLEX_DOUBLE__ *z, int ldz, 
   int n, __PRIMME_COMPLEX_DOUBLE__ *aux, double *rwork, int naux);
void Num_zheev_zprimme(const char *jobz, const char *uplo, int n, __PRIMME_COMPLEX_DOUBLE__ *a, int lda, 
   double *w, __PRIMME_COMPLEX_DOUBLE__ *work, int ldwork, double *rwork, int *info);
void Num_zgesvd_zprimme(const char *jobu, const char *jobvt, int m, int n, __PRIMME_COMPLEX_DOUBLE__ *a, int lda,
    double *s, __PRIMME_COMPLEX_DOUBLE__ *u, int ldu, __PRIMME_COMPLEX_DOUBLE__ *vt, int ldvt, __PRIMME_COMPLEX_DOUBLE__ *work,
    int ldwork, double *rwork, int *info);
void Num_zhetrf_zprimme(const char *uplo, int n, __PRIMME_COMPLEX_DOUBLE__ *a, int lda, int *ipivot,
   __PRIMME_COMPLEX_DOUBLE__ *work, int ldwork, int *info);
void Num_zhetrs_zprimme(const char *uplo, int n, int nrhs, __PRIMME_COMPLEX_DOUBLE__ *a, int lda, 
   int *ipivot, __PRIMME_COMPLEX_DOUBLE__ *b, int ldb, int *info);


void Num_copy_dprimme(int n, double *x, int incx, double *y, int incy);
void Num_copy_zprimme(int n, __PRIMME_COMPLEX_DOUBLE__ *x, int incx, __PRIMME_COMPLEX_DOUBLE__ *y, int incy);
__PRIMME_COMPLEX_DOUBLE__ Num_dot_zprimme(int n, __PRIMME_COMPLEX_DOUBLE__ *x, int incx, __PRIMME_COMPLEX_DOUBLE__ *y, int incy);
void Num_orgqr_zprimme(int m, int n, int k, __PRIMME_COMPLEX_DOUBLE__ *a, int lda, __PRIMME_COMPLEX_DOUBLE__ *tau,
      __PRIMME_COMPLEX_DOUBLE__ *rwork, int lrwork, int *info);
void Num_gemm_zprimme(const char *transa, const char *transb, int m, int n, int k, 
   __PRIMME_COMPLEX_DOUBLE__ alpha, __PRIMME_COMPLEX_DOUBLE__ *a, int lda, __PRIMME_COMPLEX_DOUBLE__ *b, int ldb, 
   __PRIMME_COMPLEX_DOUBLE__ beta, __PRIMME_COMPLEX_DOUBLE__ *c, int ldc);
void Num_symm_zprimme(const char *side, const char *uplo, int m, int n, __PRIMME_COMPLEX_DOUBLE__ alpha, 
   __PRIMME_COMPLEX_DOUBLE__ *a, int lda, __PRIMME_COMPLEX_DOUBLE__ *b, int ldb, __PRIMME_COMPLEX_DOUBLE__ beta, 
   __PRIMME_COMPLEX_DOUBLE__ *c, int ldc);
void Num_symv_zprimme(const char *uplo, int n, __PRIMME_COMPLEX_DOUBLE__ alpha, 
   __PRIMME_COMPLEX_DOUBLE__ *a, int lda, __PRIMME_COMPLEX_DOUBLE__ *x, int lncx, __PRIMME_COMPLEX_DOUBLE__ beta, 
   __PRIMME_COMPLEX_DOUBLE__ *y, int lncy); 
void Num_axpy_zprimme(int n, __PRIMME_COMPLEX_DOUBLE__ alpha, __PRIMME_COMPLEX_DOUBLE__ *x, int incx, 
   __PRIMME_COMPLEX_DOUBLE__ *y, int incy);
void Num_gemv_zprimme(const char *transa, int m, int n, __PRIMME_COMPLEX_DOUBLE__ alpha, __PRIMME_COMPLEX_DOUBLE__ *a,
   int lda, __PRIMME_COMPLEX_DOUBLE__ *x, int incx, __PRIMME_COMPLEX_DOUBLE__ beta, __PRIMME_COMPLEX_DOUBLE__ *y, int incy);
void Num_larnv_zprimme(int idist, int *iseed, int length, __PRIMME_COMPLEX_DOUBLE__ *x);
void Num_scal_zprimme(int n, __PRIMME_COMPLEX_DOUBLE__ alpha, __PRIMME_COMPLEX_DOUBLE__ *x, int incx);
void Num_swap_zprimme(int n, __PRIMME_COMPLEX_DOUBLE__ *x, int incx, __PRIMME_COMPLEX_DOUBLE__ *y, int incy);
void Num_copy_matrix_zprimme(__PRIMME_COMPLEX_DOUBLE__ *x, int m, int n, int ldx, __PRIMME_COMPLEX_DOUBLE__ *y, int ldy);
void Num_zero_matrix_zprimme(__PRIMME_COMPLEX_DOUBLE__ *x, int m, int n, int ldx);
void Num_copy_trimatrix_zprimme(__PRIMME_COMPLEX_DOUBLE__ *x, int m, int n, int ldx, int ul, int i0, __PRIMME_COMPLEX_DOUBLE__ *y, int ldy, int zero);
void Num_geqrf_zprimme(int m, int n, __PRIMME_COMPLEX_DOUBLE__ *a, int lda, __PRIMME_COMPLEX_DOUBLE__ *tau, __PRIMME_COMPLEX_DOUBLE__ *rwork, int lrwork, int *info);
int Num_update_VWXR_zprimme(__PRIMME_COMPLEX_DOUBLE__ *V, __PRIMME_COMPLEX_DOUBLE__ *W, int mV, int nV, int ldV,
   __PRIMME_COMPLEX_DOUBLE__ *h, int nh, int ldh, double *hVals,
   __PRIMME_COMPLEX_DOUBLE__ *X0, int nX0b, int nX0e, int ldX0,
   __PRIMME_COMPLEX_DOUBLE__ *X1, int nX1b, int nX1e, int ldX1,
   __PRIMME_COMPLEX_DOUBLE__ *X2, int nX2b, int nX2e, int ldX2,
   __PRIMME_COMPLEX_DOUBLE__ *Wo, int nWob, int nWoe, int ldWo,
   __PRIMME_COMPLEX_DOUBLE__ *R, int nRb, int nRe, int ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   __PRIMME_COMPLEX_DOUBLE__ *rwork, int lrwork, primme_params *primme);
void Num_compute_residual_zprimme(int n, double eval, __PRIMME_COMPLEX_DOUBLE__ *x, __PRIMME_COMPLEX_DOUBLE__ *Ax, __PRIMME_COMPLEX_DOUBLE__ *r);
void permute_vecs_iprimme(int *vecs, int n, int *perm_, int *iwork);
void permute_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm_,
      double *rwork, int *iwork);
void permute_vecs_zprimme(__PRIMME_COMPLEX_DOUBLE__ *vecs, int m, int n, int ld, int *perm_,
      __PRIMME_COMPLEX_DOUBLE__ *rwork, int *iwork);
double* Num_compact_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm,
      double *work, int ldwork, int avoidCopy);
__PRIMME_COMPLEX_DOUBLE__* Num_compact_vecs_zprimme(__PRIMME_COMPLEX_DOUBLE__ *vecs, int m, int n, int ld, int *perm,
      __PRIMME_COMPLEX_DOUBLE__ *work, int ldwork, int avoidCopy);
void Num_copy_compact_trimatrix_zprimme(__PRIMME_COMPLEX_DOUBLE__ *x, int m, int n, int i0, __PRIMME_COMPLEX_DOUBLE__ *y, int ldy);
void Num_copy_trimatrix_compact_zprimme(__PRIMME_COMPLEX_DOUBLE__ *x, int m, int n, int ldx, int i0, __PRIMME_COMPLEX_DOUBLE__ *y, int *ly);
void Num_copy_matrix_columns_dprimme(double *x, int m, int *xin, int n, int ldx, double *y,
      int *yin, int ldy);
void Num_copy_matrix_columns_zprimme(__PRIMME_COMPLEX_DOUBLE__ *x, int m, int *xin, int n, int ldx, __PRIMME_COMPLEX_DOUBLE__ *y,
      int *yin, int ldy);
int Num_compute_residual_columns_zprimme(int m, double *evals, __PRIMME_COMPLEX_DOUBLE__ *x, int n, int *p,
   int ldx, __PRIMME_COMPLEX_DOUBLE__ *Ax, int ldAx,
   __PRIMME_COMPLEX_DOUBLE__ *xo, int no, int ldxo, int io0, __PRIMME_COMPLEX_DOUBLE__ *ro, int ldro,
   __PRIMME_COMPLEX_DOUBLE__ *xd, int nd, int *pd, int ldxd, __PRIMME_COMPLEX_DOUBLE__ *rd, int ldrd,
   __PRIMME_COMPLEX_DOUBLE__ *rwork, int lrwork);
void Num_trmm_zprimme(const char *side, const char *uplo, const char *transa,
   const char *diag, int m, int n, __PRIMME_COMPLEX_DOUBLE__ alpha, __PRIMME_COMPLEX_DOUBLE__ *a, int lda, __PRIMME_COMPLEX_DOUBLE__ *b,
   int ldb);
void Num_trsm_zprimme(const char *side, const char *uplo, const char *transa, const char *diag,
      int m, int n, __PRIMME_COMPLEX_DOUBLE__ alpha, __PRIMME_COMPLEX_DOUBLE__ *a, int lda, __PRIMME_COMPLEX_DOUBLE__ *b, int ldb);
int compute_submatrix_zprimme(__PRIMME_COMPLEX_DOUBLE__ *X, int nX, int ldX, 
   __PRIMME_COMPLEX_DOUBLE__ *H, int nH, int ldH, __PRIMME_COMPLEX_DOUBLE__ *R, int ldR,
   __PRIMME_COMPLEX_DOUBLE__ *rwork, int lrwork);
double Num_lamch_zprimme(const char *cmach);

#define PRIMME_BLOCK_SIZE 512

#ifdef __cplusplus
}
#endif

#endif
