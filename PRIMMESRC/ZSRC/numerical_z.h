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
#include "Complexz.h"

#ifdef __cplusplus
extern "C" {
#endif

int Num_zhpev_zprimme(int iopt, Complex_Z *ap, double *w, Complex_Z *z, int ldz, 
   int n, Complex_Z *aux, double *rwork, int naux);
void Num_zheev_zprimme(const char *jobz, const char *uplo, int n, Complex_Z *a, int lda, 
   double *w, Complex_Z *work, int ldwork, double *rwork, int *info);
void Num_zgesvd_zprimme(const char *jobu, const char *jobvt, int m, int n, Complex_Z *a, int lda,
    double *s, Complex_Z *u, int ldu, Complex_Z *vt, int ldvt, Complex_Z *work,
    int ldwork, double *rwork, int *info);
void Num_zhetrf_zprimme(const char *uplo, int n, Complex_Z *a, int lda, int *ipivot,
   Complex_Z *work, int ldwork, int *info);
void Num_zhetrs_zprimme(const char *uplo, int n, int nrhs, Complex_Z *a, int lda, 
   int *ipivot, Complex_Z *b, int ldb, int *info);


void Num_dcopy_dprimme(int n, double *x, int incx, double *y, int incy);
void Num_zcopy_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy);
Complex_Z Num_dot_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy);
void Num_orgqr_zprimme(int m, int n, int k, Complex_Z *a, int lda, Complex_Z *tau,
      Complex_Z *rwork, int lrwork, int *info);
void Num_gemm_zprimme(const char *transa, const char *transb, int m, int n, int k, 
   Complex_Z alpha, Complex_Z *a, int lda, Complex_Z *b, int ldb, 
   Complex_Z beta, Complex_Z *c, int ldc);
void Num_symm_zprimme(const char *side, const char *uplo, int m, int n, Complex_Z alpha, 
   Complex_Z *a, int lda, Complex_Z *b, int ldb, Complex_Z beta, 
   Complex_Z *c, int ldc);
void Num_symv_zprimme(const char *uplo, int n, Complex_Z alpha, 
   Complex_Z *a, int lda, Complex_Z *x, int lncx, Complex_Z beta, 
   Complex_Z *y, int lncy); 
void Num_axpy_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx, 
   Complex_Z *y, int incy);
void Num_gemv_zprimme(const char *transa, int m, int n, Complex_Z alpha, Complex_Z *a,
   int lda, Complex_Z *x, int incx, Complex_Z beta, Complex_Z *y, int incy);
void Num_larnv_zprimme(int idist, int *iseed, int length, Complex_Z *x);
void Num_scal_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx);
void Num_swap_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy);
void Num_copy_matrix_zprimme(Complex_Z *x, int m, int n, int ldx, Complex_Z *y, int ldy);
void Num_copy_trimatrix_zprimme(Complex_Z *x, int m, int n, int ldx, int ul, int i0, Complex_Z *y, int ldy, int zero);
void Num_geqrf_zprimme(int m, int n, Complex_Z *a, int lda, Complex_Z *tau, Complex_Z *rwork, int lrwork, int *info);
int Num_update_VWXR_zprimme(Complex_Z *V, Complex_Z *W, int mV, int nV, int ldV,
   Complex_Z *h, int nh, int ldh, double *hVals,
   Complex_Z *X0, int nX0b, int nX0e, int ldX0,
   Complex_Z *X1, int nX1b, int nX1e, int ldX1,
   Complex_Z *X2, int nX2b, int nX2e, int ldX2,
   Complex_Z *Wo, int nWob, int nWoe, int ldWo,
   Complex_Z *R, int nRb, int nRe, int ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   Complex_Z *rwork, int lrwork, primme_params *primme);
void Num_compute_residual_zprimme(int n, double eval, Complex_Z *x, Complex_Z *Ax, Complex_Z *r);
void permute_vecs_iprimme(int *vecs, int n, int *perm_, int *iwork);
void permute_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm_,
      double *rwork, int *iwork);
void permute_vecs_zprimme(Complex_Z *vecs, int m, int n, int ld, int *perm_,
      Complex_Z *rwork, int *iwork);
double* Num_compact_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm,
      double *work, int ldwork, int avoidCopy);
Complex_Z* Num_compact_vecs_zprimme(Complex_Z *vecs, int m, int n, int ld, int *perm,
      Complex_Z *work, int ldwork, int avoidCopy);
void Num_copy_compact_trimatrix_zprimme(Complex_Z *x, int m, int n, int i0, Complex_Z *y, int ldy);
void Num_copy_trimatrix_compact_zprimme(Complex_Z *x, int m, int n, int ldx, int i0, Complex_Z *y, int *ly);
void Num_copy_matrix_i_zprimme(Complex_Z *x, int m, int *xin, int n, int ldx, Complex_Z *y,
      int *yin, int ldy);
int Num_compute_residual_i_zprimme(int m, double *evals, Complex_Z *x, int n, int *p,
   int ldx, Complex_Z *Ax, int ldAx,
   Complex_Z *xo, int no, int ldxo, Complex_Z *ro, int ldro,
   Complex_Z *xd, int nd, int *pd, int ldxd, Complex_Z *rd, int ldrd,
   Complex_Z *rwork, int lrwork);
void Num_trmm_zprimme(const char *side, const char *uplo, const char *transa,
   const char *diag, int m, int n, Complex_Z alpha, Complex_Z *a, int lda, Complex_Z *b,
   int ldb);
void Num_trsm_zprimme(const char *side, const char *uplo, const char *transa, const char *diag,
      int m, int n, Complex_Z alpha, Complex_Z *a, int lda, Complex_Z *b, int ldb);

#define PRIMME_BLOCK_SIZE 512

#ifdef __cplusplus
}
#endif

#endif
