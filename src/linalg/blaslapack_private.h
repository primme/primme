/*******************************************************************************
 * Copyright (c) 2017, College of William & Mary
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the College of William & Mary nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * PRIMME: https://github.com/primme/primme
 * Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
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

#include "numerical.h"

#if !defined(PRIMME_BLASINT_SIZE) || PRIMME_BLASINT_SIZE == 32
#  include <stdint.h>
#  include <inttypes.h>
#  define PRIMME_BLASINT int32_t
#  define PRIMME_BLASINT_P PRId32
#  define PRIMME_BLASINT_MAX INT32_MAX
#elif PRIMME_BLASINT_SIZE == 0
#  include <limits.h>
#  define PRIMME_BLASINT int
#  define PRIMME_BLASINT_P "d"
#  define PRIMME_BLASINT_MAX INT_MAX
#elif PRIMME_BLASINT_SIZE == 64
#  include <stdint.h>
#  include <inttypes.h>
#  define PRIMME_BLASINT int64_t
#  define PRIMME_BLASINT_P PRId64
#  define PRIMME_BLASINT_MAX INT64_MAX
#else
#  define PRIMME_BLASINT PRIMME_BLASINT_SIZE
#  define PRIMME_BLASINT_P "d"
#  define PRIMME_BLASINT_MAX ((PRIMME_BLASINT_SIZE)INT_MAX)*INT_MAX
#endif

#ifndef NUM_CRAY

#ifdef USE_DOUBLE
#  define LAPACK_FUNCTION(S,C,D,Z) FORTRAN_FUNCTION(D)
#elif defined(USE_DOUBLECOMPLEX)
#  define LAPACK_FUNCTION(S,C,D,Z) FORTRAN_FUNCTION(Z)
#elif defined(USE_FLOAT)
#  define LAPACK_FUNCTION(S,C,D,Z) FORTRAN_FUNCTION(S)
#elif defined(USE_FLOATCOMPLEX)
#  define LAPACK_FUNCTION(S,C,D,Z) FORTRAN_FUNCTION(C)
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
#define XHEEVX    LAPACK_FUNCTION(ssyevx, cheevx, dsyevx, zheevx)
#define XGESVD    LAPACK_FUNCTION(sgesvd, cgesvd, dgesvd, zgesvd)
#define XHETRF    LAPACK_FUNCTION(ssytrf, chetrf, dsytrf, zhetrf)
#define XHETRS    LAPACK_FUNCTION(ssytrs, chetrs, dsytrs, zhetrs)

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
#define XSCAL  LAPACK_FUNCTION(SSCAL  , zscal )
#define XLARNV LAPACK_FUNCTION(SLARNV ,       )
#define XHEEV  LAPACK_FUNCTION(SSYEV  , zheev )
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
#ifndef USE_COMPLEX
SCALAR XDOT(PRIMME_BLASINT *n, SCALAR *x, PRIMME_BLASINT *incx, SCALAR *y, PRIMME_BLASINT *incy);
void XHEEVX(STRING jobz, STRING range, STRING uplo, PRIMME_BLASINT *n, SCALAR *a, PRIMME_BLASINT *lda, SCALAR *vl, SCALAR *vu, PRIMME_BLASINT *il, PRIMME_BLASINT *iu,  SCALAR *abstol, PRIMME_BLASINT *m,  SCALAR *w, SCALAR *z, PRIMME_BLASINT *ldz, SCALAR *work, PRIMME_BLASINT *ldwork, PRIMME_BLASINT *iwork, PRIMME_BLASINT *ifail, PRIMME_BLASINT *info);
void XGESVD(STRING jobu, STRING jobvt, PRIMME_BLASINT *m, PRIMME_BLASINT *n, SCALAR *a, PRIMME_BLASINT *lda, SCALAR *s, SCALAR *u, PRIMME_BLASINT *ldu, SCALAR *vt, PRIMME_BLASINT *ldvt, SCALAR *work, PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info); 
#else
void XHEEVX(STRING jobz, STRING range, STRING uplo, PRIMME_BLASINT *n, SCALAR *a, PRIMME_BLASINT *lda, REAL *vl, REAL *vu, PRIMME_BLASINT *il, PRIMME_BLASINT *iu, REAL *abstol, PRIMME_BLASINT *m,  REAL *w, SCALAR *z, PRIMME_BLASINT *ldz, SCALAR *work, PRIMME_BLASINT *ldwork, REAL *rwork, PRIMME_BLASINT *iwork, PRIMME_BLASINT *ifail, PRIMME_BLASINT *info);
void XGESVD(STRING jobu, STRING jobvt, PRIMME_BLASINT *m, PRIMME_BLASINT *n, SCALAR *a, PRIMME_BLASINT *lda, REAL *s, SCALAR *u, PRIMME_BLASINT *ldu, SCALAR *vt, PRIMME_BLASINT *ldvt, SCALAR *work, PRIMME_BLASINT *ldwork, REAL *rwork, PRIMME_BLASINT *info);
#endif
void XSCAL(PRIMME_BLASINT *n, SCALAR *alpha, SCALAR *x, PRIMME_BLASINT *incx);
void XLARNV(PRIMME_BLASINT *idist, PRIMME_BLASINT *iseed, PRIMME_BLASINT *n, SCALAR *x);
void XHETRF(STRING uplo, PRIMME_BLASINT *n, SCALAR *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot, SCALAR *work, PRIMME_BLASINT *ldwork, PRIMME_BLASINT *info);
void XHETRS(STRING uplo, PRIMME_BLASINT *n, PRIMME_BLASINT *nrhs, SCALAR *a, PRIMME_BLASINT *lda, PRIMME_BLASINT *ipivot, SCALAR *b, PRIMME_BLASINT *ldb, PRIMME_BLASINT *info);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NUMERICAL_PRIVATE_H */
