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

#include "numerical_@(pre).h"

#if !defined(PRIMME_BLASINT_SIZE) || PRIMME_BLASINT_SIZE == 64
#  include <stdint.h>
#  include <inttypes.h>
#  define PRIMME_BLASINT int64_t
#  define PRIMME_BLASINT_P PRId64
#  define PRIMME_BLASINT_MAX INT64_MAX
#elif PRIMME_BLASINT_SIZE == 0
#  include <limits.h>
#  define PRIMME_BLASINT int
#  define PRIMME_BLASINT_P "d"
#  define PRIMME_BLASINT_MAX INT_MAX
#elif PRIMME_INT == 32
#  include <stdint.h>
#  include <inttypes.h>
#  define PRIMME_BLASINT int32_t
#  define PRIMME_BLASINT_P PRId32
#  define PRIMME_BLASINT_MAX INT32_MAX
#else
#  define PRIMME_BLASINT PRIMME_INT_SIZE
#  define PRIMME_BLASINT_P "d"
#  define PRIMME_BLASINT_MAX INT_MAX
#endif

#ifndef NUM_CRAY

#ifdef USE_DOUBLE
#  define LAPACK_FUNCTION(S,C,D,Z) FORTRAN_FUNCTION(D)
#elif defined(USE_DOUBLECOMPLEX)
#  define LAPACK_FUNCTION(S,C,D,Z) FORTRAN_FUNCTION(Z)
#endif

#define XCOPY     LAPACK_FUNCTION(scopy , ccopy , dcopy , zcopy )   
#define XSWAP     LAPACK_FUNCTION(sswap , cswap , dswap , zswap )
#define XGEMM     LAPACK_FUNCTION(sgemm , cgemm , dgemm , zgemm )
#define XTRMM     LAPACK_FUNCTION(strmm , ctrmm , dtrmm , ztrmm )
#define XTRSM     LAPACK_FUNCTION(strsm , ctrsm , dtrsm , ztrsm )
#define XHEMM     LAPACK_FUNCTION(ssymm , chemm , dsymm , zhemm )
#define XHEMV     LAPACK_FUNCTION(ssymv , chemv , dsymv , zhemv )
#define XAXPY     LAPACK_FUNCTION(saxpy , caxpy , daxpy , zaxpy )
#define XGEMV     LAPACK_FUNCTION(sgemv , cgemv , dgemv , zgemv )
#define XDOT      LAPACK_FUNCTION(sdot  ,       , ddot  ,       )
#define XSCAL     LAPACK_FUNCTION(sscal , cscal , dscal , zscal )
#define XLARNV    LAPACK_FUNCTION(slarnv, clarnv, dlarnv, zlarnv)
#define XHEEV     LAPACK_FUNCTION(ssyev , cheev , dsyev , zheev )
#define XGESVD    LAPACK_FUNCTION(sgesvd, cgesvd, dgesvd, zgesvd)
#define XHETRF    LAPACK_FUNCTION(ssytrf, chetrf, dsytrf, zhetrf)
#define XHETRS    LAPACK_FUNCTION(ssytrs, chetrs, dsytrs, zhetrs)
#define XLAMCH    LAPACK_FUNCTION(slamch, clamch, dlamch, zlamch)

#ifdef NUM_ESSL
#include <essl.h>
#endif

#else /* NUM_CRAY */

#include <fortran.h>
#include <string.h>

#ifdef USE_DOUBLE
#  define LAPACK_FUNCTION(D,Z) D
#elif defined(USE_DOUBLECOMPLEX)
#  define LAPACK_FUNCTION(D,Z) Z
#endif

#define XCOPY  LAPACK_FUNCTION(SCOPY  , zcopy )
#define XSWAP  LAPACK_FUNCTION(SSWAP  , zswap )
#define XGEMM  LAPACK_FUNCTION(SGEMM  , zgemm )
#define XTRMM  LAPACK_FUNCTION(STRMM  , ztrmm )
#define XTRSM  LAPACK_FUNCTION(STRSM  , ztrsm )
#define XSYMM  LAPACK_FUNCTION(DSYMM  , zhemm )
#define XSYMV  LAPACK_FUNCTION(DSYMV  , zhemv )
#define XAXPY  LAPACK_FUNCTION(SAXPY  , zaxpy )
#define XGEMV  LAPACK_FUNCTION(SGEMV  , zgemv )
#define XDOT   LAPACK_FUNCTION(SDOT   ,       )
#define XLAMCH LAPACK_FUNCTION(SLAMCH , zlarnv)
#define XSCAL  LAPACK_FUNCTION(SSCAL  , zscal )
#define XLARNV LAPACK_FUNCTION(SLARNV ,       )
#define XSYEV  LAPACK_FUNCTION(SSYEV  , zheev )
#define XGESVD LAPACK_FUNCTION(SGESVD , zhetrf)
#define XSYTRF LAPACK_FUNCTION(SSYTRF , zgesvd)
#define XSYTRS LAPACK_FUNCTION(SSYTRS , zhetrs)

#endif /* NUM_CRAY */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef NUM_CRAY
#  define STRING const char * 
#else
#  define STRING _fcd
#endif

void XCOPY(PRIMME_BLASINT *n, SCALAR *x, PRIMME_BLASINT *incx, SCALAR *y, PRIMME_BLASINT *incy);
void XSWAP(PRIMME_BLASINT *n, SCALAR *x, PRIMME_BLASINT *incx, SCALAR *y, PRIMME_BLASINT *incy);
void XGEMM(STRING transa, STRING transb, PRIMME_BLASINT *m, PRIMME_BLASINT *n, PRIMME_BLASINT *k, SCALAR *alpha, SCALAR *a, PRIMME_BLASINT *lda, SCALAR *b, PRIMME_BLASINT *ldb, SCALAR *beta, SCALAR *c, PRIMME_BLASINT *ldc);
void XGEMV(STRING transa, PRIMME_BLASINT *m, PRIMME_BLASINT *n, SCALAR *alpha, SCALAR *a, PRIMME_BLASINT *lda, SCALAR *x, PRIMME_BLASINT *incx, SCALAR *beta, SCALAR *y, PRIMME_BLASINT *incy);
void XTRMM(STRING side, STRING uplo, STRING transa, STRING diag, PRIMME_BLASINT *m, PRIMME_BLASINT *n, SCALAR *alpha, SCALAR *a, PRIMME_BLASINT *lda, SCALAR *b, PRIMME_BLASINT *ldb);
void XTRSM(STRING side, STRING uplo, STRING transa, STRING diag, PRIMME_BLASINT *m, PRIMME_BLASINT *n, SCALAR *alpha, SCALAR *a, PRIMME_BLASINT *lda, SCALAR *b, PRIMME_BLASINT *ldb);
void XHEMM(STRING side, STRING uplo, PRIMME_BLASINT *m, PRIMME_BLASINT *n, SCALAR *alpha, SCALAR *a, PRIMME_BLASINT *lda, SCALAR *b, PRIMME_BLASINT *ldb, SCALAR *beta, SCALAR *c, PRIMME_BLASINT *ldc);
void XHEMV(STRING uplo, PRIMME_BLASINT *n, SCALAR *alpha, SCALAR *a, PRIMME_BLASINT *lda, SCALAR *x, PRIMME_BLASINT *lncx, SCALAR *beta, SCALAR *y, PRIMME_BLASINT *lncy);
void XAXPY(PRIMME_BLASINT *n, SCALAR *alpha, SCALAR *x, PRIMME_BLASINT *incx, SCALAR *y, PRIMME_BLASINT *incy);
#ifdef USE_DOUBLE
SCALAR XDOT(PRIMME_BLASINT *n, SCALAR *x, PRIMME_BLASINT *incx, SCALAR *y, PRIMME_BLASINT *incy);
void XHEEV(STRING jobz, STRING uplo, PRIMME_BLASINT *n, SCALAR *a, PRIMME_BLASINT *lda, SCALAR *w, SCALAR *work, PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info);
void XGESVD(STRING jobu, STRING jobvt, PRIMME_BLASINT *m, PRIMME_BLASINT *n, SCALAR *a, PRIMME_BLASINT *lda, SCALAR *s, SCALAR *u, PRIMME_BLASINT *ldu, SCALAR *vt, PRIMME_BLASINT *ldvt, SCALAR *work, PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info); 
#elif defined(USE_DOUBLECOMPLEX)
void XHEEV(STRING jobz, STRING uplo, PRIMME_BLASINT *n, SCALAR *a, PRIMME_BLASINT *lda, REAL *w, SCALAR *work, PRIMME_BLASINT *ldwork, REAL *rwork, PRIMME_BLASINT *info);
void XGESVD(STRING jobu, STRING jobvt, PRIMME_BLASINT *m, PRIMME_BLASINT *n, SCALAR *a, PRIMME_BLASINT *lda, REAL *s, SCALAR *u, PRIMME_BLASINT *ldu, SCALAR *vt, PRIMME_BLASINT *ldvt, SCALAR *work, PRIMME_BLASINT *ldwork, REAL *rwork, PRIMME_BLASINT *info);
#endif
SCALAR XLAMCH(STRING cmach);
void XSCAL(PRIMME_BLASINT *n, SCALAR *alpha, SCALAR *x, PRIMME_BLASINT *incx);
void XLARNV(PRIMME_BLASINT *idist, PRIMME_BLASINT *iseed, PRIMME_BLASINT *n, SCALAR *x);
void XHETRF(STRING uplo, PRIMME_BLASINT *n, SCALAR *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot, SCALAR *work, PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info);
void XHETRS(STRING uplo, PRIMME_BLASINT *n, PRIMME_BLASINT *nrhs, SCALAR *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot, SCALAR *b, PRIMME_BLASINT *ldb, PRIMME_BLASINT *info);

#ifdef NUM_ESSL
#  ifdef USE_DOUBLE
PRIMME_BLASINT dspev(PRIMME_BLASINT iopt, SCALAR *ap, SCALAR *w, SCALAR *z, PRIMME_BLASINT ldz, PRIMME_BLASINT n, SCALAR *aux, PRIMME_BLASINT naux);
#  elif defined(USE_DOUBLECOMPLEX)
PRIMME_BLASINT zhpev(PRIMME_BLASINT iopt, void *ap, SCALAR *w, void *z, PRIMME_BLASINT ldz, PRIMME_BLASINT n, void *aux, PRIMME_BLASINT naux);
#  endif
#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUMERICAL_PRIVATE_H */
