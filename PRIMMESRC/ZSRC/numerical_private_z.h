/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2005  James R. McCombs,  Andreas Stathopoulos
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
 * File: numerical_private.h
 *
 * Purpose - Contains definitions and prototypes for exclusive use with 
 *           numerical.c.  There are various definitions for use with Sun,
 *           IBM, and Cray.
 *
 * Module name      : %M%
 * SID              : %I%
 * Date             : %G%
 ******************************************************************************/

#ifndef NUMERICAL_H
#define NUMERICAL_H

#if !defined(NUM_SUN) && !defined(NUM_IBM) && !defined(NUM_CRAY)
#define NUM_SUN
#endif

#ifdef NUM_SUN

#define ZDOTCSUB  zdotcsub_
#define ZCOPY  zcopy_
#define ZSWAP  zswap_
#define ZGEMM  zgemm_
#define ZHEMM  zhemm_
#define ZAXPY  zaxpy_
#define ZGEMV  zgemv_
#define ZSCAL  zscal_
#define ZLARNV zlarnv_
#define ZHEEV  zheev_
#define ZHETRF zhetrf_
#define ZHETRS zhetrs_

#define DCOPY  dcopy_
#define DSWAP  dswap_
#define DGEMM  dgemm_
#define DSYMM  dsymm_
#define DAXPY  daxpy_
#define DGEMV  dgemv_
#define DDOT   ddot_
#define DSCAL  dscal_
#define DLARNV dlarnv_
#define DSYEV  dsyev_
#define DSYTRF dsytrf_
#define DSYTRS dsytrs_

#elif defined(NUM_IBM)

#define ZDOTCSUB  zdotcsub
#define ZCOPY  zcopy
#define ZSWAP  zswap
#define ZGEMM  zgemm
#define ZHEMM  zhemm
#define ZAXPY  zaxpy
#define ZGEMV  zgemv
#define ZSCAL  zscal
#define ZLARNV zlarnv
#define ZHEEV  zheev
#define ZHETRF zhetrf
#define ZHETRS zhetrs

#define DCOPY  dcopy
#define DSWAP  dswap
#define DGEMM  dgemm
#define DSYMM  dsymm
#define DAXPY  daxpy
#define DGEMV  dgemv
#define DDOT   ddot
#define DSCAL  dscal
#define DLARNV dlarnv
#define DSYEV  dsyev
#define DSYTRF dsytrf
#define DSYTRS dsytrs

#ifdef NUM_ESSL
#include <essl.h>
#endif

#elif defined(NUM_CRAY)
#include <fortran.h>
#include <string.h>

#define ZDOTCSUB  zdotcsub
#define ZCOPY  zcopy
#define ZSWAP  zswap
#define ZGEMM  zgemm
#define ZHEMM  zhemm
#define ZAXPY  zaxpy
#define ZGEMV  zgemv
#define ZSCAL  zscal
#define ZLARNV zlarnv
#define ZHEEV  zheev
#define ZHETRF zhetrf
#define ZHETRS zhetrs

#define DCOPY  SCOPY
#define DSWAP  SSWAP
#define DGEMM  SGEMM
#define DSYMM  DSYMM
#define DAXPY  SAXPY
#define DGEMV  SGEMV
#define DDOT   SDOT
#define DLAMCH SLAMCH
#define DSCAL  SSCAL
#define DLARNV SLARNV
#define DSYEV  SSYEV
#define DSYTRF DSYTRF
#define DSYTRS DSYTRS

#endif
#ifdef Cplusplus
extern "C" {
#endif /* Cplusplus */

#ifndef NUM_CRAY

void DCOPY(int *n, double *x, int *incx, double *y, int *incy);
void DSWAP(int *n, double *x, int *incx, double *y, int *incy);
void DGEMM(char *transa, char *transb, int *m, int *n, int *k, double *alpha, 
   double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void DSYMM(char *side, char *uplo, int *m, int *n, double *alpha, double *a,
   int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void DAXPY(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
void DGEMV(char *transa, int *m, int *n, double *alpha, double *a, int *lda, 
   double *x, int *incx, double *beta, double *y, int *incy);
double DDOT(int *n, double *x, int *incx, double *y, int *incy);
double DLAMCH(char *cmach);
void DSCAL(int *n, double *alpha, double *x, int *incx);
void DLARNV(int *idist, int *iseed, int *n, double *x);
void DSYEV(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
   double *work, int *ldwork, int *info);
void DSYTRF(char *uplo, int *n, double *a, int *lda, int *ipivot, double *work,
   int *ldwork, int *info);
void DSYTRS(char *uplo, int *n, int *nrhs, double *a, int *lda, int *ipivot,
   double *b, int *ldb, int *info);

void   ZCOPY(int *n, void *x, int *incx, void *y, int *incy);
void   ZSWAP(int *n, void *x, int *incx, void *y, int *incy);
void   ZGEMM(char *transa, char *transb, int *m, int *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
void   ZHEMM(char *side, char *uplo, int *m, int *n, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
void   ZAXPY(int *n, void *alpha, void *x, int *incx, void *y, int *incy);
void   ZGEMV(char *transa, int *m, int *n, void *alpha, void *a, int *lda, void *x, int *incx, void *beta, void *y, int *incy);
void   ZSCAL(int *n, void *alpha, void *x, int *incx);
void   ZLARNV(int *idist, int *iseed, int *n, void *x);
void   ZHEEV(char *jobz, char *uplo, int *n, void *a, int *lda, double *w, void *work, int *ldwork, double *rwork, int *info);
void   ZHETRF(char *uplo, int *n, void *a, int *lda, int *ipivot, void *work, int *ldwork, int *info);
void   ZHETRS(char *uplo, int *n, int *nrhs, void *a, int *lda, int *ipivot, void *b, int *ldb, int *info);
void   ZDOTCSUB(void *dot, int *n, void *x, int *incx, void *y, int *incy);

#ifdef NUM_ESSL
int dspev(int iopt, double *ap, double *w, double *z, int ldz, int n, double *aux, int naux);
int zhpev(int iopt, void *ap, double *w, void *z, int ldz, int n, void *aux, int naux);
#endif

#else

void DCOPY(int *n, double *x, int *incx, double *y, int *incy);
void DSWAP(int *n, double *x, int *incx, double *y, int *incy);
void DGEMM(_fcd transa_fcd, _fcd transb_fcd, int *m, int *n, int *k, 
   double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, 
   double *c, int *ldc);
void DSYMM(_fcd side_fcd, _fcd uplo_fcd, int *m, int *n, double *alpha, 
   double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void DAXPY(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
void DGEMV(_fcd transa_fcd, int *m, int *n, double *alpha, double *a, int *lda, 
   double *x, int *incx, double *beta, double *y, int *incy);
double DDOT(int *n, double *x, int *incx, double *y, int *incy);
double DLAMCH(_fcd cmach_fcd);
void DSCAL(int *n, double *alpha, double *x, int *incx);
void DLARNV(int *idist, int *iseed, int *n, double *x);
void DSYEV(_fcd jobz_fcd, _fcd uplo_fcd, int *n, double *a, int *lda, double *w,
   double *work, int *ldwork, int *info);

void DSYTRF(_fcd uplo, int *n, double *a, int *lda, int *ipivot, double *work,
   int *ldwork, int *info);
void DSYTRS(_fcd uplo, int *n, int *nrhs, double *a, int *lda, int *ipivot,
   double *b, int *ldb, int *info);

void   ZCOPY(int *n, void *x, int *incx, void *y, int *incy);
void   ZSWAP(int *n, void *x, int *incx, void *y, int *incy);
void   ZGEMM(_fcd transa, _fcd transb, int *m, int *n, int *k, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
void   ZHEMM(_fcd side, _fcd uplo, int *m, int *n, void *alpha, void *a, int *lda, void *b, int *ldb, void *beta, void *c, int *ldc);
void   ZAXPY(int *n, void *alpha, void *x, int *incx, void *y, int *incy);
void   ZGEMV(_fcd transa, int *m, int *n, void *alpha, void *a, int *lda, void *x, int *incx, void *beta, void *y, int *incy);
void   ZSCAL(int *n, void *alpha, void *x, int *incx);
void   ZLARNV(int *idist, int *iseed, int *n, void *x);
void   ZHEEV(_fcd jobz, _fcd uplo, int *n, void *a, int *lda, double *w, void *work, int *ldwork, double *rwork, int *info);
void   ZDOTCSUB(void *dot, int *n, void *x, int *incx, void *y, int *incy);

void   ZHETRF(_fcd uplo, int *n, void *a, int *lda, int *ipivot, void *work, int *ldwork, int *info);
void   ZHETRS(_fcd uplo, int *n, int *nrhs, void *a, int *lda, int *ipivot, void *b, int *ldb, int *info);

#endif /* else (if not cray)*/

#ifdef Cplusplus
}
#endif /* Cplusplus */

#endif /* NUMERICAL_H */
