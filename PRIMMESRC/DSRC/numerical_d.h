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

#ifndef USE_DOUBLE
#define USE_DOUBLE
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

// TEMP: Error management system

#define CHKERR(ERRN, RETURN) { \
   int err = (ERRN);\
   if (err) {\
      fprintf(stderr, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", err, __LINE__, #ERRN );\
      return (RETURN);\
   }\
}
#define CHKERRM(ERRN, RETURN, ...) { \
   int err = (ERRN);\
   if (err) {\
      fprintf(stderr, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", err, __LINE__, #ERRN );\
      fprintf(stderr, "PRIMME: " __VA_ARGS__);\
      return (RETURN);\
   }\
}

#define TO_INT(X) ((X) < INT_MAX ? (X) : INT_MAX)

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

#ifdef USE_DOUBLECOMPLEX
int Num_hpev_zprimme(int iopt, SCALAR *ap, double *w, SCALAR *z, int ldz, 
   int n, SCALAR *aux, double *rwork, int naux);
void Num_heev_zprimme(const char *jobz, const char *uplo, int n, SCALAR *a, int lda, 
   double *w, SCALAR *work, int ldwork, double *rwork, int *info);
void Num_gesvd_zprimme(const char *jobu, const char *jobvt, int m, int n, SCALAR *a, int lda,
    double *s, SCALAR *u, int ldu, SCALAR *vt, int ldvt, SCALAR *work,
    int ldwork, double *rwork, int *info);
#endif

#ifdef USE_DOUBLE
int Num_hpev_dprimme(int iopt, double *ap, double *w, double *z, int ldz, 
   int n, double *aux, int naux);
void Num_heev_dprimme(const char *jobz, const char *uplo, int n, double *a, int lda, 
   double *w, double *work, int ldwork, int *info);
void Num_gesvd_dprimme(const char *jobu, const char *jobvt, int m, int n, double *a, int lda,
    double *s, double *u, int ldu, double *vt, int ldvt, double *work,
    int ldwork, int *info);
#endif

void Num_hetrf_dprimme(const char *uplo, int n, SCALAR *a, int lda, int *ipivot,
   SCALAR *work, int ldwork, int *info);
void Num_hetrs_dprimme(const char *uplo, int n, int nrhs, SCALAR *a, int lda, 
   int *ipivot, SCALAR *b, int ldb, int *info);
void Num_copy_dprimme(int n, double *x, int incx, double *y, int incy);
#ifdef USE_DOUBLECOMPLEX
void Num_copy_zprimme(int n, SCALAR *x, int incx, SCALAR *y, int incy);
#endif
SCALAR Num_dot_dprimme(int n, SCALAR *x, int incx, SCALAR *y, int incy);
void Num_gemm_dprimme(const char *transa, const char *transb, int m, int n, int k, 
   SCALAR alpha, SCALAR *a, int lda, SCALAR *b, int ldb, 
   SCALAR beta, SCALAR *c, int ldc);
void Num_hemm_dprimme(const char *side, const char *uplo, int m, int n,
SCALAR alpha, 
   SCALAR *a, int lda, SCALAR *b, int ldb, SCALAR beta, 
   SCALAR *c, int ldc);
void Num_hemv_dprimme(const char *uplo, int n, SCALAR alpha, 
   SCALAR *a, int lda, SCALAR *x, int lncx, SCALAR beta, 
   SCALAR *y, int lncy); 
void Num_axpy_dprimme(int n, SCALAR alpha, SCALAR *x, int incx, 
   SCALAR *y, int incy);
void Num_gemv_dprimme(const char *transa, int m, int n, SCALAR alpha, SCALAR *a,
   int lda, SCALAR *x, int incx, SCALAR beta, SCALAR *y, int incy);
void Num_larnv_dprimme(int idist, PRIMME_INT *iseed, PRIMME_INT length, SCALAR *x);
void Num_scal_dprimme(int n, SCALAR alpha, SCALAR *x, int incx);
void Num_swap_dprimme(int n, SCALAR *x, int incx, SCALAR *y, int incy);
void Num_copy_matrix_dprimme(SCALAR *x, int m, int n, int ldx, SCALAR *y, int ldy);
void Num_zero_matrix_dprimme(SCALAR *x, int m, int n, int ldx);
void Num_copy_trimatrix_dprimme(SCALAR *x, int m, int n, int ldx, int ul, int i0, SCALAR *y, int ldy, int zero);
int Num_update_VWXR_dprimme(SCALAR *V, SCALAR *W, int mV, int nV, int ldV,
   SCALAR *h, int nh, int ldh, double *hVals,
   SCALAR *X0, int nX0b, int nX0e, int ldX0,
   SCALAR *X1, int nX1b, int nX1e, int ldX1,
   SCALAR *X2, int nX2b, int nX2e, int ldX2,
   SCALAR *Wo, int nWob, int nWoe, int ldWo,
   SCALAR *R, int nRb, int nRe, int ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   SCALAR *rwork, int lrwork, primme_params *primme);
void Num_compute_residual_dprimme(int n, double eval, SCALAR *x, SCALAR *Ax, SCALAR *r);
void permute_vecs_iprimme(int *vecs, int n, int *perm_, int *iwork);
void permute_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm_,
      double *rwork, int *iwork);
#ifdef USE_DOUBLECOMPLEX
void permute_vecs_zprimme(SCALAR *vecs, int m, int n, int ld, int *perm_,
      SCALAR *rwork, int *iwork);
#endif
double* Num_compact_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm,
      double *work, int ldwork, int avoidCopy);
#ifdef USE_DOUBLECOMPLEX
SCALAR* Num_compact_vecs_dprimme(SCALAR *vecs, int m, int n, int ld, int *perm,
      SCALAR *work, int ldwork, int avoidCopy);
#endif
void Num_copy_compact_trimatrix_dprimme(SCALAR *x, int m, int n, int i0, SCALAR *y, int ldy);
void Num_copy_trimatrix_compact_dprimme(SCALAR *x, int m, int n, int ldx, int i0, SCALAR *y, int *ly);
void Num_copy_matrix_columns_dprimme(double *x, int m, int *xin, int n, int ldx, double *y,
      int *yin, int ldy);
#ifdef USE_DOUBLECOMPLEX
void Num_copy_matrix_columns_dprimme(SCALAR *x, int m, int *xin, int n, int
ldx, SCALAR *y,
      int *yin, int ldy);
#endif
int Num_compute_residual_columns_dprimme(int m, double *evals, SCALAR *x, int n, int *p,
   int ldx, SCALAR *Ax, int ldAx,
   SCALAR *xo, int no, int ldxo, int io0, SCALAR *ro, int ldro,
   SCALAR *xd, int nd, int *pd, int ldxd, SCALAR *rd, int ldrd,
   SCALAR *rwork, int lrwork);
void Num_trmm_dprimme(const char *side, const char *uplo, const char *transa,
   const char *diag, int m, int n, SCALAR alpha, SCALAR *a, int lda, SCALAR *b,
   int ldb);
void Num_trsm_dprimme(const char *side, const char *uplo, const char *transa, const char *diag,
      int m, int n, SCALAR alpha, SCALAR *a, int lda, SCALAR *b, int ldb);
int compute_submatrix_dprimme(SCALAR *X, int nX, int ldX, 
   SCALAR *H, int nH, int ldH, SCALAR *R, int ldR,
   SCALAR *rwork, size_t *lrwork);
double Num_lamch_dprimme(const char *cmach);

#define PRIMME_BLOCK_SIZE 512

#ifdef __cplusplus
}
#endif

#endif
