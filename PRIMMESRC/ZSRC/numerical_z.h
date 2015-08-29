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
#include "Complexz.h"

#ifdef __cplusplus
extern "C" {
#endif

int Num_zhpev_zprimme(int iopt, Complex_Z *ap, double *w, Complex_Z *z, int ldz, 
   int n, Complex_Z *aux, double *rwork, int naux);
void Num_zheev_zprimme(const char *jobz, const char *uplo, int n, Complex_Z *a, int lda, 
   double *w, Complex_Z *work, int ldwork, double *rwork, int *info);
void Num_zhetrf_zprimme(const char *uplo, int n, Complex_Z *a, int lda, int *ipivot,
   Complex_Z *work, int ldwork, int *info);
void Num_zhetrs_zprimme(const char *uplo, int n, int nrhs, Complex_Z *a, int lda, 
   int *ipivot, Complex_Z *b, int ldb, int *info);


void Num_zcopy_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy);
Complex_Z Num_dot_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy);
void Num_gemm_zprimme(const char *transa, const char *transb, int m, int n, int k, 
   Complex_Z alpha, Complex_Z *a, int lda, Complex_Z *b, int ldb, 
   Complex_Z beta, Complex_Z *c, int ldc);
void Num_symm_zprimme(const char *side, const char *uplo, int m, int n, Complex_Z alpha, 
   Complex_Z *a, int lda, Complex_Z *b, int ldb, Complex_Z beta, 
   Complex_Z *c, int ldc);
void Num_axpy_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx, 
   Complex_Z *y, int incy);
void Num_gemv_zprimme(const char *transa, int m, int n, Complex_Z alpha, Complex_Z *a,
   int lda, Complex_Z *x, int incx, Complex_Z beta, Complex_Z *y, int incy);
void Num_larnv_zprimme(int idist, int *iseed, int length, Complex_Z *x);
void Num_scal_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx);
void Num_swap_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy);

#ifdef __cplusplus
}
#endif

#endif


