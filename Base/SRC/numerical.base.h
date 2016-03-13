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

#include "primme.h"
#include "common_numerical.h"
#ifdefarithm L_DEFCPLX
#include "Complexz.h"
#endifarithm

#ifdef __cplusplus
extern "C" {
#endif

#ifdefarithm L_DEFCPLX
int Num_zhpev_zprimme(int iopt, @(type) *ap, double *w, @(type) *z, int ldz, 
   int n, @(type) *aux, double *rwork, int naux);
void Num_zheev_zprimme(const char *jobz, const char *uplo, int n, @(type) *a, int lda, 
   double *w, @(type) *work, int ldwork, double *rwork, int *info);
void Num_zgesvd_zprimme(const char *jobu, const char *jobvt, int m, int n, @(type) *a, int lda,
    double *s, @(type) *u, int ldu, @(type) *vt, int ldvt, @(type) *work,
    int ldwork, double *rwork, int *info);
void Num_zhetrf_zprimme(const char *uplo, int n, @(type) *a, int lda, int *ipivot,
   @(type) *work, int ldwork, int *info);
void Num_zhetrs_zprimme(const char *uplo, int n, int nrhs, @(type) *a, int lda, 
   int *ipivot, @(type) *b, int ldb, int *info);
#endifarithm

#ifdefarithm L_DEFREAL
int Num_dspev_dprimme(int iopt, double *ap, double *w, double *z, int ldz, 
   int n, double *aux, int naux);
void Num_dsyev_dprimme(const char *jobz, const char *uplo, int n, double *a, int lda, 
   double *w, double *work, int ldwork, int *info);
void Num_dgesvd_dprimme(const char *jobu, const char *jobvt, int m, int n, double *a, int lda,
    double *s, double *u, int ldu, double *vt, int ldvt, double *work,
    int ldwork, int *info);
void Num_dsytrf_dprimme(const char *uplo, int n, double *a, int lda, int *ipivot, 
   double *work, int ldwork, int *info);
void Num_dsytrs_dprimme(const char *uplo, int n, int nrhs, double *a, int lda, 
   int *ipivot, double *b, int ldb, int *info);
#endifarithm

void Num_dcopy_dprimme(int n, double *x, int incx, double *y, int incy);
#ifdefarithm L_DEFCPLX
void Num_zcopy_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy);
#endifarithm
@(type) Num_dot_@(pre)primme(int n, @(type) *x, int incx, @(type) *y, int incy);
void Num_orgqr_@(pre)primme(int m, int n, int k, @(type) *a, int lda, @(type) *tau,
      @(type) *rwork, int lrwork, int *info);
void Num_gemm_@(pre)primme(const char *transa, const char *transb, int m, int n, int k, 
   @(type) alpha, @(type) *a, int lda, @(type) *b, int ldb, 
   @(type) beta, @(type) *c, int ldc);
void Num_symm_@(pre)primme(const char *side, const char *uplo, int m, int n, @(type) alpha, 
   @(type) *a, int lda, @(type) *b, int ldb, @(type) beta, 
   @(type) *c, int ldc);
void Num_symv_@(pre)primme(const char *uplo, int n, @(type) alpha, 
   @(type) *a, int lda, @(type) *x, int lncx, @(type) beta, 
   @(type) *y, int lncy); 
void Num_axpy_@(pre)primme(int n, @(type) alpha, @(type) *x, int incx, 
   @(type) *y, int incy);
void Num_gemv_@(pre)primme(const char *transa, int m, int n, @(type) alpha, @(type) *a,
   int lda, @(type) *x, int incx, @(type) beta, @(type) *y, int incy);
void Num_larnv_@(pre)primme(int idist, int *iseed, int length, @(type) *x);
void Num_scal_@(pre)primme(int n, @(type) alpha, @(type) *x, int incx);
void Num_swap_@(pre)primme(int n, @(type) *x, int incx, @(type) *y, int incy);
void Num_copy_matrix_@(pre)primme(@(type) *x, int m, int n, int ldx, @(type) *y, int ldy);
void Num_copy_trimatrix_@(pre)primme(@(type) *x, int m, int n, int ldx, int ul, int i0, @(type) *y, int ldy, int zero);
void Num_geqrf_@(pre)primme(int m, int n, @(type) *a, int lda, @(type) *tau, @(type) *rwork, int lrwork, int *info);
int Num_update_VWXR_@(pre)primme(@(type) *V, @(type) *W, int mV, int nV, int ldV,
   @(type) *h, int nh, int ldh, double *hVals,
   @(type) *X0, int nX0b, int nX0e, int ldX0,
   @(type) *X1, int nX1b, int nX1e, int ldX1,
   @(type) *X2, int nX2b, int nX2e, int ldX2,
   @(type) *Wo, int nWob, int nWoe, int ldWo,
   @(type) *R, int nRb, int nRe, int ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   @(type) *rwork, int lrwork, primme_params *primme);
void Num_compute_residual_@(pre)primme(int n, double eval, @(type) *x, @(type) *Ax, @(type) *r);
void permute_vecs_iprimme(int *vecs, int n, int *perm_, int *iwork);
void permute_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm_,
      double *rwork, int *iwork);
#ifdefarithm L_DEFCPLX
void permute_vecs_zprimme(Complex_Z *vecs, int m, int n, int ld, int *perm_,
      Complex_Z *rwork, int *iwork);
#endifarithm
double* Num_compact_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm,
      double *work, int ldwork, int avoidCopy);
#ifdefarithm L_DEFCPLX
@(type)* Num_compact_vecs_@(pre)primme(@(type) *vecs, int m, int n, int ld, int *perm,
      @(type) *work, int ldwork, int avoidCopy);
#endifarithm
void Num_copy_compact_trimatrix_@(pre)primme(@(type) *x, int m, int n, int i0, @(type) *y, int ldy);
void Num_copy_trimatrix_compact_@(pre)primme(@(type) *x, int m, int n, int ldx, int i0, @(type) *y, int *ly);
void Num_copy_matrix_i_@(pre)primme(@(type) *x, int m, int *xin, int n, int ldx, @(type) *y,
      int *yin, int ldy);
int Num_compute_residual_i_@(pre)primme(int m, double *evals, @(type) *x, int n, int *p,
   int ldx, @(type) *Ax, int ldAx,
   @(type) *xo, int no, int ldxo, @(type) *ro, int ldro,
   @(type) *xd, int nd, int *pd, int ldxd, @(type) *rd, int ldrd,
   @(type) *rwork, int lrwork);
void Num_trmm_@(pre)primme(const char *side, const char *uplo, const char *transa,
   const char *diag, int m, int n, @(type) alpha, @(type) *a, int lda, @(type) *b,
   int ldb);
void Num_trsm_@(pre)primme(const char *side, const char *uplo, const char *transa, const char *diag,
      int m, int n, @(type) alpha, @(type) *a, int lda, @(type) *b, int ldb);

#define PRIMME_BLOCK_SIZE 512

#ifdef __cplusplus
}
#endif

#endif
