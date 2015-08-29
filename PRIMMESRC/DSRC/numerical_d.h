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

#ifdef __cplusplus
extern "C" {
#endif


int Num_dspev_dprimme(int iopt, double *ap, double *w, double *z, int ldz, 
   int n, double *aux, int naux);
void Num_dsyev_dprimme(const char *jobz, const char *uplo, int n, double *a, int lda, 
   double *w, double *work, int ldwork, int *info);
void Num_dsytrf_dprimme(const char *uplo, int n, double *a, int lda, int *ipivot, 
   double *work, int ldwork, int *info);
void Num_dsytrs_dprimme(const char *uplo, int n, int nrhs, double *a, int lda, 
   int *ipivot, double *b, int ldb, int *info);

void Num_dcopy_dprimme(int n, double *x, int incx, double *y, int incy);
double Num_dot_dprimme(int n, double *x, int incx, double *y, int incy);
void Num_gemm_dprimme(const char *transa, const char *transb, int m, int n, int k, 
   double alpha, double *a, int lda, double *b, int ldb, 
   double beta, double *c, int ldc);
void Num_symm_dprimme(const char *side, const char *uplo, int m, int n, double alpha, 
   double *a, int lda, double *b, int ldb, double beta, 
   double *c, int ldc);
void Num_axpy_dprimme(int n, double alpha, double *x, int incx, 
   double *y, int incy);
void Num_gemv_dprimme(const char *transa, int m, int n, double alpha, double *a,
   int lda, double *x, int incx, double beta, double *y, int incy);
void Num_larnv_dprimme(int idist, int *iseed, int length, double *x);
void Num_scal_dprimme(int n, double alpha, double *x, int incx);
void Num_swap_dprimme(int n, double *x, int incx, double *y, int incy);

#ifdef __cplusplus
}
#endif

#endif


