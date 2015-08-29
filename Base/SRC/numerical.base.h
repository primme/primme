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
void Num_dsytrf_dprimme(const char *uplo, int n, double *a, int lda, int *ipivot, 
   double *work, int ldwork, int *info);
void Num_dsytrs_dprimme(const char *uplo, int n, int nrhs, double *a, int lda, 
   int *ipivot, double *b, int ldb, int *info);
#endifarithm

void Num_@(pre)copy_@(pre)primme(int n, @(type) *x, int incx, @(type) *y, int incy);
@(type) Num_dot_@(pre)primme(int n, @(type) *x, int incx, @(type) *y, int incy);
void Num_gemm_@(pre)primme(const char *transa, const char *transb, int m, int n, int k, 
   @(type) alpha, @(type) *a, int lda, @(type) *b, int ldb, 
   @(type) beta, @(type) *c, int ldc);
void Num_symm_@(pre)primme(const char *side, const char *uplo, int m, int n, @(type) alpha, 
   @(type) *a, int lda, @(type) *b, int ldb, @(type) beta, 
   @(type) *c, int ldc);
void Num_axpy_@(pre)primme(int n, @(type) alpha, @(type) *x, int incx, 
   @(type) *y, int incy);
void Num_gemv_@(pre)primme(const char *transa, int m, int n, @(type) alpha, @(type) *a,
   int lda, @(type) *x, int incx, @(type) beta, @(type) *y, int incy);
void Num_larnv_@(pre)primme(int idist, int *iseed, int length, @(type) *x);
void Num_scal_@(pre)primme(int n, @(type) alpha, @(type) *x, int incx);
void Num_swap_@(pre)primme(int n, @(type) *x, int incx, @(type) *y, int incy);

#ifdef __cplusplus
}
#endif

#endif


