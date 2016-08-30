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
#  include <complex.h>
#  ifdef I
#     undef I
#  endif
#  define IMAGINARY _Complex_I
#  define SCALAR complex double
#  define REAL double
#  define REAL_PART(x) (creal(x))
#  define ABS(x) (cabs(x))
#  define CONJ(x) (conj(x))
#else
#  define IMAGINARY 0.0
#  define SCALAR double
#  define REAL double
#  define REAL_PART(x) (x)
#  define ABS(x) (fabs(x))
#  define CONJ(x) (x)
#endif

#include <tgmath.h>   /* select proper function abs from fabs, cabs... */

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

int Num_zhpev_zprimme(int iopt, complex double *ap, double *w, complex double *z, int ldz, 
   int n, complex double *aux, double *rwork, int naux);
void Num_zheev_zprimme(const char *jobz, const char *uplo, int n, complex double *a, int lda, 
   double *w, complex double *work, int ldwork, double *rwork, int *info);
void Num_zgesvd_zprimme(const char *jobu, const char *jobvt, int m, int n, complex double *a, int lda,
    double *s, complex double *u, int ldu, complex double *vt, int ldvt, complex double *work,
    int ldwork, double *rwork, int *info);
void Num_zhetrf_zprimme(const char *uplo, int n, complex double *a, int lda, int *ipivot,
   complex double *work, int ldwork, int *info);
void Num_zhetrs_zprimme(const char *uplo, int n, int nrhs, complex double *a, int lda, 
   int *ipivot, complex double *b, int ldb, int *info);


void Num_copy_dprimme(int n, double *x, int incx, double *y, int incy);
void Num_copy_zprimme(int n, complex double *x, int incx, complex double *y, int incy);
complex double Num_dot_zprimme(int n, complex double *x, int incx, complex double *y, int incy);
void Num_orgqr_zprimme(int m, int n, int k, complex double *a, int lda, complex double *tau,
      complex double *rwork, int lrwork, int *info);
void Num_gemm_zprimme(const char *transa, const char *transb, int m, int n, int k, 
   complex double alpha, complex double *a, int lda, complex double *b, int ldb, 
   complex double beta, complex double *c, int ldc);
void Num_symm_zprimme(const char *side, const char *uplo, int m, int n, complex double alpha, 
   complex double *a, int lda, complex double *b, int ldb, complex double beta, 
   complex double *c, int ldc);
void Num_symv_zprimme(const char *uplo, int n, complex double alpha, 
   complex double *a, int lda, complex double *x, int lncx, complex double beta, 
   complex double *y, int lncy); 
void Num_axpy_zprimme(int n, complex double alpha, complex double *x, int incx, 
   complex double *y, int incy);
void Num_gemv_zprimme(const char *transa, int m, int n, complex double alpha, complex double *a,
   int lda, complex double *x, int incx, complex double beta, complex double *y, int incy);
void Num_larnv_zprimme(int idist, int *iseed, int length, complex double *x);
void Num_scal_zprimme(int n, complex double alpha, complex double *x, int incx);
void Num_swap_zprimme(int n, complex double *x, int incx, complex double *y, int incy);
void Num_copy_matrix_zprimme(complex double *x, int m, int n, int ldx, complex double *y, int ldy);
void Num_zero_matrix_zprimme(complex double *x, int m, int n, int ldx);
void Num_copy_trimatrix_zprimme(complex double *x, int m, int n, int ldx, int ul, int i0, complex double *y, int ldy, int zero);
void Num_geqrf_zprimme(int m, int n, complex double *a, int lda, complex double *tau, complex double *rwork, int lrwork, int *info);
int Num_update_VWXR_zprimme(complex double *V, complex double *W, int mV, int nV, int ldV,
   complex double *h, int nh, int ldh, double *hVals,
   complex double *X0, int nX0b, int nX0e, int ldX0,
   complex double *X1, int nX1b, int nX1e, int ldX1,
   complex double *X2, int nX2b, int nX2e, int ldX2,
   complex double *Wo, int nWob, int nWoe, int ldWo,
   complex double *R, int nRb, int nRe, int ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   complex double *rwork, int lrwork, primme_params *primme);
void Num_compute_residual_zprimme(int n, double eval, complex double *x, complex double *Ax, complex double *r);
void permute_vecs_iprimme(int *vecs, int n, int *perm_, int *iwork);
void permute_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm_,
      double *rwork, int *iwork);
void permute_vecs_zprimme(complex double *vecs, int m, int n, int ld, int *perm_,
      complex double *rwork, int *iwork);
double* Num_compact_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm,
      double *work, int ldwork, int avoidCopy);
complex double* Num_compact_vecs_zprimme(complex double *vecs, int m, int n, int ld, int *perm,
      complex double *work, int ldwork, int avoidCopy);
void Num_copy_compact_trimatrix_zprimme(complex double *x, int m, int n, int i0, complex double *y, int ldy);
void Num_copy_trimatrix_compact_zprimme(complex double *x, int m, int n, int ldx, int i0, complex double *y, int *ly);
void Num_copy_matrix_columns_dprimme(double *x, int m, int *xin, int n, int ldx, double *y,
      int *yin, int ldy);
void Num_copy_matrix_columns_zprimme(complex double *x, int m, int *xin, int n, int ldx, complex double *y,
      int *yin, int ldy);
int Num_compute_residual_columns_zprimme(int m, double *evals, complex double *x, int n, int *p,
   int ldx, complex double *Ax, int ldAx,
   complex double *xo, int no, int ldxo, int io0, complex double *ro, int ldro,
   complex double *xd, int nd, int *pd, int ldxd, complex double *rd, int ldrd,
   complex double *rwork, int lrwork);
void Num_trmm_zprimme(const char *side, const char *uplo, const char *transa,
   const char *diag, int m, int n, complex double alpha, complex double *a, int lda, complex double *b,
   int ldb);
void Num_trsm_zprimme(const char *side, const char *uplo, const char *transa, const char *diag,
      int m, int n, complex double alpha, complex double *a, int lda, complex double *b, int ldb);
int compute_submatrix_zprimme(complex double *X, int nX, int ldX, 
   complex double *H, int nH, int ldH, complex double *R, int ldR,
   complex double *rwork, int lrwork);
double Num_lamch_zprimme(const char *cmach);

#define PRIMME_BLOCK_SIZE 512

#ifdef __cplusplus
}
#endif

#endif
