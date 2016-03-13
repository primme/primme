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

#ifdef __cplusplus
extern "C" {
#endif


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

void Num_dcopy_dprimme(int n, double *x, int incx, double *y, int incy);
double Num_dot_dprimme(int n, double *x, int incx, double *y, int incy);
void Num_orgqr_dprimme(int m, int n, int k, double *a, int lda, double *tau,
      double *rwork, int lrwork, int *info);
void Num_gemm_dprimme(const char *transa, const char *transb, int m, int n, int k, 
   double alpha, double *a, int lda, double *b, int ldb, 
   double beta, double *c, int ldc);
void Num_symm_dprimme(const char *side, const char *uplo, int m, int n, double alpha, 
   double *a, int lda, double *b, int ldb, double beta, 
   double *c, int ldc);
void Num_symv_dprimme(const char *uplo, int n, double alpha, 
   double *a, int lda, double *x, int lncx, double beta, 
   double *y, int lncy); 
void Num_axpy_dprimme(int n, double alpha, double *x, int incx, 
   double *y, int incy);
void Num_gemv_dprimme(const char *transa, int m, int n, double alpha, double *a,
   int lda, double *x, int incx, double beta, double *y, int incy);
void Num_larnv_dprimme(int idist, int *iseed, int length, double *x);
void Num_scal_dprimme(int n, double alpha, double *x, int incx);
void Num_swap_dprimme(int n, double *x, int incx, double *y, int incy);
void Num_copy_matrix_dprimme(double *x, int m, int n, int ldx, double *y, int ldy);
void Num_copy_trimatrix_dprimme(double *x, int m, int n, int ldx, int ul, int i0, double *y, int ldy, int zero);
void Num_geqrf_dprimme(int m, int n, double *a, int lda, double *tau, double *rwork, int lrwork, int *info);
int Num_update_VWXR_dprimme(double *V, double *W, int mV, int nV, int ldV,
   double *h, int nh, int ldh, double *hVals,
   double *X0, int nX0b, int nX0e, int ldX0,
   double *X1, int nX1b, int nX1e, int ldX1,
   double *X2, int nX2b, int nX2e, int ldX2,
   double *Wo, int nWob, int nWoe, int ldWo,
   double *R, int nRb, int nRe, int ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   double *rwork, int lrwork, primme_params *primme);
void Num_compute_residual_dprimme(int n, double eval, double *x, double *Ax, double *r);
void permute_vecs_iprimme(int *vecs, int n, int *perm_, int *iwork);
void permute_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm_,
      double *rwork, int *iwork);
double* Num_compact_vecs_dprimme(double *vecs, int m, int n, int ld, int *perm,
      double *work, int ldwork, int avoidCopy);
void Num_copy_compact_trimatrix_dprimme(double *x, int m, int n, int i0, double *y, int ldy);
void Num_copy_trimatrix_compact_dprimme(double *x, int m, int n, int ldx, int i0, double *y, int *ly);
void Num_copy_matrix_i_dprimme(double *x, int m, int *xin, int n, int ldx, double *y,
      int *yin, int ldy);
int Num_compute_residual_i_dprimme(int m, double *evals, double *x, int n, int *p,
   int ldx, double *Ax, int ldAx,
   double *xo, int no, int ldxo, double *ro, int ldro,
   double *xd, int nd, int *pd, int ldxd, double *rd, int ldrd,
   double *rwork, int lrwork);
void Num_trmm_dprimme(const char *side, const char *uplo, const char *transa,
   const char *diag, int m, int n, double alpha, double *a, int lda, double *b,
   int ldb);
void Num_trsm_dprimme(const char *side, const char *uplo, const char *transa, const char *diag,
      int m, int n, double alpha, double *a, int lda, double *b, int ldb);

#define PRIMME_BLOCK_SIZE 512

#ifdef __cplusplus
}
#endif

#endif
