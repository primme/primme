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
 * File: numerical_private.h
 *
 * Purpose - Contains definitions and prototypes for exclusive use with 
 *           numerical.c.  There are various definitions for use with Sun,
 *           IBM, and Cray.
 *
 ******************************************************************************/

#ifndef NUMERICAL_PRIVATE_H
#define NUMERICAL_PRIVATE_H

#include "common_numerical.h"

#ifdef F77UNDERSCORE
#define FORTRAN_FUNCTION(X) X ## _
#else
#define FORTRAN_FUNCTION(X) X
#endif

#ifndef NUM_CRAY

#define ZDOTCSUB  FORTRAN_FUNCTION(zdotcsub)
#define ZCOPY     FORTRAN_FUNCTION(zcopy)
#define ZSWAP     FORTRAN_FUNCTION(zswap)
#define ZGEMM     FORTRAN_FUNCTION(zgemm)
#define ZHEMM     FORTRAN_FUNCTION(zhemm)
#define ZAXPY     FORTRAN_FUNCTION(zaxpy)
#define ZGEMV     FORTRAN_FUNCTION(zgemv)
#define ZSCAL     FORTRAN_FUNCTION(zscal)
#define ZLARNV    FORTRAN_FUNCTION(zlarnv)
#define ZHEEV     FORTRAN_FUNCTION(zheev)
#define ZHETRF    FORTRAN_FUNCTION(zhetrf)
#define ZHETRS    FORTRAN_FUNCTION(zhetrs)

#define DCOPY     FORTRAN_FUNCTION(dcopy)
#define DSWAP     FORTRAN_FUNCTION(dswap)
#define DGEMM     FORTRAN_FUNCTION(dgemm)
#define DSYMM     FORTRAN_FUNCTION(dsymm)
#define DAXPY     FORTRAN_FUNCTION(daxpy)
#define DGEMV     FORTRAN_FUNCTION(dgemv)
#define DDOT      FORTRAN_FUNCTION(ddot)
#define DSCAL     FORTRAN_FUNCTION(dscal)
#define DLARNV    FORTRAN_FUNCTION(dlarnv)
#define DSYEV     FORTRAN_FUNCTION(dsyev)
#define DSYTRF    FORTRAN_FUNCTION(dsytrf)
#define DSYTRS    FORTRAN_FUNCTION(dsytrs)

#ifdef NUM_ESSL
#include <essl.h>
#endif

#else /* NUM_CRAY */

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

#endif /* NUM_CRAY */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef NUM_CRAY

void DCOPY(PRIMME_BLASINT *n, double *x, PRIMME_BLASINT *incx, double *y, PRIMME_BLASINT *incy);
void DSWAP(PRIMME_BLASINT *n, double *x, PRIMME_BLASINT *incx, double *y, PRIMME_BLASINT *incy);
void DGEMM(const char *transa, const char *transb, PRIMME_BLASINT *m, PRIMME_BLASINT *n, PRIMME_BLASINT *k, double *alpha, 
   double *a, PRIMME_BLASINT *lda, double *b, PRIMME_BLASINT *ldb, double *beta, double *c, PRIMME_BLASINT *ldc);
void DSYMM(const char *side, const char *uplo, PRIMME_BLASINT *m, PRIMME_BLASINT *n, double *alpha, double *a,
   PRIMME_BLASINT *lda, double *b, PRIMME_BLASINT *ldb, double *beta, double *c, PRIMME_BLASINT *ldc);
void DAXPY(PRIMME_BLASINT *n, double *alpha, double *x, PRIMME_BLASINT *incx, double *y, PRIMME_BLASINT *incy);
void DGEMV(const char *transa, PRIMME_BLASINT *m, PRIMME_BLASINT *n, double *alpha, double *a, PRIMME_BLASINT *lda, 
   double *x, PRIMME_BLASINT *incx, double *beta, double *y, PRIMME_BLASINT *incy);
double DDOT(PRIMME_BLASINT *n, double *x, PRIMME_BLASINT *incx, double *y, PRIMME_BLASINT *incy);
double DLAMCH(const char *cmach);
void DSCAL(PRIMME_BLASINT *n, double *alpha, double *x, PRIMME_BLASINT *incx);
void DLARNV(PRIMME_BLASINT *idist, PRIMME_BLASINT *iseed, PRIMME_BLASINT *n, double *x);
void DSYEV(const char *jobz, const char *uplo, PRIMME_BLASINT *n, double *a, PRIMME_BLASINT *lda, double *w,
   double *work, PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info);
void DSYTRF(const char *uplo, PRIMME_BLASINT *n, double *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot, double *work,
   PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info);
void DSYTRS(const char *uplo, PRIMME_BLASINT *n, PRIMME_BLASINT *nrhs, double *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot,
   double *b, PRIMME_BLASINT *ldb, PRIMME_BLASINT *info);

void   ZCOPY(PRIMME_BLASINT *n, void *x, PRIMME_BLASINT *incx, void *y, PRIMME_BLASINT *incy);
void   ZSWAP(PRIMME_BLASINT *n, void *x, PRIMME_BLASINT *incx, void *y, PRIMME_BLASINT *incy);
void   ZGEMM(const char *transa, const char *transb, PRIMME_BLASINT *m, PRIMME_BLASINT *n, PRIMME_BLASINT *k, void *alpha, void *a, PRIMME_BLASINT *lda, void *b, PRIMME_BLASINT *ldb, void *beta, void *c, PRIMME_BLASINT *ldc);
void   ZHEMM(const char *side, const char *uplo, PRIMME_BLASINT *m, PRIMME_BLASINT *n, void *alpha, void *a, PRIMME_BLASINT *lda, void *b, PRIMME_BLASINT *ldb, void *beta, void *c, PRIMME_BLASINT *ldc);
void   ZAXPY(PRIMME_BLASINT *n, void *alpha, void *x, PRIMME_BLASINT *incx, void *y, PRIMME_BLASINT *incy);
void   ZGEMV(const char *transa, PRIMME_BLASINT *m, PRIMME_BLASINT *n, void *alpha, void *a, PRIMME_BLASINT *lda, void *x, PRIMME_BLASINT *incx, void *beta, void *y, PRIMME_BLASINT *incy);
void   ZSCAL(PRIMME_BLASINT *n, void *alpha, void *x, PRIMME_BLASINT *incx);
void   ZLARNV(PRIMME_BLASINT *idist, PRIMME_BLASINT *iseed, PRIMME_BLASINT *n, void *x);
void   ZHEEV(const char *jobz, const char *uplo, PRIMME_BLASINT *n, void *a, PRIMME_BLASINT *lda, double *w, void *work, PRIMME_BLASINT *ldwork, double *rwork, PRIMME_BLASINT *info);
void   ZHETRF(const char *uplo, PRIMME_BLASINT *n, void *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot, void *work, PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info);
void   ZHETRS(const char *uplo, PRIMME_BLASINT *n, PRIMME_BLASINT *nrhs, void *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot, void *b, PRIMME_BLASINT *ldb, PRIMME_BLASINT *info);
void   ZDOTCSUB(void *dot, PRIMME_BLASINT *n, void *x, PRIMME_BLASINT *incx, void *y, PRIMME_BLASINT *incy);

#ifdef NUM_ESSL
PRIMME_BLASINT dspev(PRIMME_BLASINT iopt, double *ap, double *w, double *z, PRIMME_BLASINT ldz, PRIMME_BLASINT n, double *aux, PRIMME_BLASINT naux);
PRIMME_BLASINT zhpev(PRIMME_BLASINT iopt, void *ap, double *w, void *z, PRIMME_BLASINT ldz, PRIMME_BLASINT n, void *aux, PRIMME_BLASINT naux);
#endif

#else /* NUM_CRAY */

void DCOPY(PRIMME_BLASINT *n, double *x, PRIMME_BLASINT *incx, double *y, PRIMME_BLASINT *incy);
void DSWAP(PRIMME_BLASINT *n, double *x, PRIMME_BLASINT *incx, double *y, PRIMME_BLASINT *incy);
void DGEMM(_fcd transa_fcd, _fcd transb_fcd, PRIMME_BLASINT *m, PRIMME_BLASINT *n, PRIMME_BLASINT *k, 
   double *alpha, double *a, PRIMME_BLASINT *lda, double *b, PRIMME_BLASINT *ldb, double *beta, 
   double *c, PRIMME_BLASINT *ldc);
void DSYMM(_fcd side_fcd, _fcd uplo_fcd, PRIMME_BLASINT *m, PRIMME_BLASINT *n, double *alpha, 
   double *a, PRIMME_BLASINT *lda, double *b, PRIMME_BLASINT *ldb, double *beta, double *c, PRIMME_BLASINT *ldc);
void DAXPY(PRIMME_BLASINT *n, double *alpha, double *x, PRIMME_BLASINT *incx, double *y, PRIMME_BLASINT *incy);
void DGEMV(_fcd transa_fcd, PRIMME_BLASINT *m, PRIMME_BLASINT *n, double *alpha, double *a, PRIMME_BLASINT *lda, 
   double *x, PRIMME_BLASINT *incx, double *beta, double *y, PRIMME_BLASINT *incy);
double DDOT(PRIMME_BLASINT *n, double *x, PRIMME_BLASINT *incx, double *y, PRIMME_BLASINT *incy);
double DLAMCH(_fcd cmach_fcd);
void DSCAL(PRIMME_BLASINT *n, double *alpha, double *x, PRIMME_BLASINT *incx);
void DLARNV(PRIMME_BLASINT *idist, PRIMME_BLASINT *iseed, PRIMME_BLASINT *n, double *x);
void DSYEV(_fcd jobz_fcd, _fcd uplo_fcd, PRIMME_BLASINT *n, double *a, PRIMME_BLASINT *lda, double *w,
   double *work, PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info);

void DSYTRF(_fcd uplo, PRIMME_BLASINT *n, double *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot, double *work,
   PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info);
void DSYTRS(_fcd uplo, PRIMME_BLASINT *n, PRIMME_BLASINT *nrhs, double *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot,
   double *b, PRIMME_BLASINT *ldb, PRIMME_BLASINT *info);

void   ZCOPY(PRIMME_BLASINT *n, void *x, PRIMME_BLASINT *incx, void *y, PRIMME_BLASINT *incy);
void   ZSWAP(PRIMME_BLASINT *n, void *x, PRIMME_BLASINT *incx, void *y, PRIMME_BLASINT *incy);
void   ZGEMM(_fcd transa, _fcd transb, PRIMME_BLASINT *m, PRIMME_BLASINT *n, PRIMME_BLASINT *k, void *alpha, void *a, PRIMME_BLASINT *lda, void *b, PRIMME_BLASINT *ldb, void *beta, void *c, PRIMME_BLASINT *ldc);
void   ZHEMM(_fcd side, _fcd uplo, PRIMME_BLASINT *m, PRIMME_BLASINT *n, void *alpha, void *a, PRIMME_BLASINT *lda, void *b, PRIMME_BLASINT *ldb, void *beta, void *c, PRIMME_BLASINT *ldc);
void   ZAXPY(PRIMME_BLASINT *n, void *alpha, void *x, PRIMME_BLASINT *incx, void *y, PRIMME_BLASINT *incy);
void   ZGEMV(_fcd transa, PRIMME_BLASINT *m, PRIMME_BLASINT *n, void *alpha, void *a, PRIMME_BLASINT *lda, void *x, PRIMME_BLASINT *incx, void *beta, void *y, PRIMME_BLASINT *incy);
void   ZSCAL(PRIMME_BLASINT *n, void *alpha, void *x, PRIMME_BLASINT *incx);
void   ZLARNV(PRIMME_BLASINT *idist, PRIMME_BLASINT *iseed, PRIMME_BLASINT *n, void *x);
void   ZHEEV(_fcd jobz, _fcd uplo, PRIMME_BLASINT *n, void *a, PRIMME_BLASINT *lda, double *w, void *work, PRIMME_BLASINT *ldwork, double *rwork, PRIMME_BLASINT *info);
void   ZDOTCSUB(void *dot, PRIMME_BLASINT *n, void *x, PRIMME_BLASINT *incx, void *y, PRIMME_BLASINT *incy);

void   ZHETRF(_fcd uplo, PRIMME_BLASINT *n, void *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot, void *work, PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info);
void   ZHETRS(_fcd uplo, PRIMME_BLASINT *n, PRIMME_BLASINT *nrhs, void *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot, void *b, PRIMME_BLASINT *ldb, PRIMME_BLASINT *info);

#endif /* NUM_CRAY */

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUMERICAL_PRIVATE_H */
