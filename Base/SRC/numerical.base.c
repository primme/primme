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
 * File: numerical.c
 *
 * Purpose - This file contains for the most part C wrapper routines for
 *    calling various BLAS and LAPACK FORTRAN routines.
 *
 ******************************************************************************/


#include <stdarg.h>
#ifdefarithm L_DEFCPLX
#include "Complex.h"
#endifarithm
#include "numerical_private_@(pre).h"
#include "numerical_@(pre).h"

/******************************************************************************/
void Num_@(pre)copy_@(pre)primme(int n, @(type) *x, int incx, @(type) *y, int incy) {

#ifdefarithm L_DEFCPLX
   ZCOPY(&n, x, &incx, y, &incy);
#endifarithm
#ifdefarithm L_DEFREAL
   DCOPY(&n, x, &incx, y, &incy);
#endifarithm
}
/******************************************************************************/

void Num_gemm_@(pre)primme(char *transa, char *transb, int m, int n, int k, 
   @(type) alpha, @(type) *a, int lda, @(type) *b, int ldb, 
   @(type) beta, @(type) *c, int ldc) {

#ifdef NUM_CRAY
   _fcd transa_fcd, transb_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
   transb_fcd = _cptofcd(transb, strlen(transb));
#ifdefarithm L_DEFCPLX
   ZGEMM(transa_fcd, transb_fcd, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, 
         c, &ldc);
#endifarithm
#ifdefarithm L_DEFREAL
   DGEMM(transa_fcd, transb_fcd, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, 
         c, &ldc);
#endifarithm
#else
#ifdefarithm L_DEFCPLX
   ZGEMM(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endifarithm
#ifdefarithm L_DEFREAL
   DGEMM(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endifarithm
#endif

}

/******************************************************************************/
void Num_symm_@(pre)primme(char *side, char *uplo, int m, int n, @(type) alpha, 
   @(type) *a, int lda, @(type) *b, int ldb, @(type) beta, 
   @(type) *c, int ldc) {

#ifdef NUM_CRAY
   _fcd side_fcd, uplo_fcd;

   side_fcd = _cptofcd(side, strlen(side));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));
#ifdefarithm L_DEFCPLX
   ZHEMM(side_fcd, uplo_fcd, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endifarithm
#ifdefarithm L_DEFREAL
   DSYMM(side_fcd, uplo_fcd, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endifarithm
#else
#ifdefarithm L_DEFCPLX
   ZHEMM(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endifarithm
#ifdefarithm L_DEFREAL
   DSYMM(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endifarithm
#endif   

}

/******************************************************************************/
void Num_axpy_@(pre)primme(int n, @(type) alpha, @(type) *x, int incx, 
   @(type) *y, int incy) {

#ifdefarithm L_DEFCPLX
   ZAXPY(&n, &alpha, x, &incx, y, &incy);
#endifarithm
#ifdefarithm L_DEFREAL
   DAXPY(&n, &alpha, x, &incx, y, &incy);
#endifarithm

}

/******************************************************************************/
void Num_gemv_@(pre)primme(char *transa, int m, int n, @(type) alpha, @(type) *a,
   int lda, @(type) *x, int incx, @(type) beta, @(type) *y, int incy) {

#ifdef NUM_CRAY
   _fcd transa_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
#ifdefarithm L_DEFCPLX
   ZGEMV(transa_fcd, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endifarithm
#ifdefarithm L_DEFREAL
   DGEMV(transa_fcd, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endifarithm
#else
#ifdefarithm L_DEFCPLX
   ZGEMV(transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endifarithm
#ifdefarithm L_DEFREAL
   DGEMV(transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endifarithm

#endif

}

/******************************************************************************/
@(type) Num_dot_@(pre)primme(int n, @(type) *x, int incx, @(type) *y, int incy) {

#ifdefarithm L_DEFCPLX
/* ---- Explicit implementation of the zdotc() --- */
   int i;
   Complex_Z zdotc = {+0.0e+00,+0.0e00};
   if (n <= 0) return(zdotc);
   if (incx==1 && incy==1) {
      for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
        zdotc.r += x[i].r*y[i].r + x[i].i*y[i].i;
        zdotc.i += x[i].r*y[i].i - x[i].i*y[i].r;
      }
   }
   else {
      int ix,iy;
      ix = 0;
      iy = 0;
      if(incx <= 0) ix = (-n+1)*incx;
      if(incy <= 0) iy = (-n+1)*incy;
      for (i=0;i<n;i++) {
        zdotc.r += x[ix].r*y[iy].r + x[ix].i*y[iy].i;
        zdotc.i += x[ix].r*y[iy].i - x[ix].i*y[iy].r;
        ix += incx;
        iy += incy;
      }
   }
   return(zdotc);
/* -- end of explicit implementation of the zdotc() - */
#endifarithm
#ifdefarithm L_DEFREAL
   return(DDOT(&n, x, &incx, y, &incy));
#endifarithm

}

/******************************************************************************/
void Num_larnv_@(pre)primme(int idist, int *iseed, int length, @(type) *x) {
#ifdefarithm L_DEFCPLX
   ZLARNV(&idist, iseed, &length, x);
#endifarithm
#ifdefarithm L_DEFREAL
   DLARNV(&idist, iseed, &length, x);
#endifarithm

}

/******************************************************************************/
void Num_scal_@(pre)primme(int n, @(type) alpha, @(type) *x, int incx) {

#ifdefarithm L_DEFCPLX
   ZSCAL(&n, &alpha, x, &incx);
#endifarithm
#ifdefarithm L_DEFREAL
   DSCAL(&n, &alpha, x, &incx);
#endifarithm

}

/******************************************************************************/
void Num_swap_@(pre)primme(int n, @(type) *x, int incx, @(type) *y, int incy) {

#ifdefarithm L_DEFCPLX
   ZSWAP(&n, x, &incx, y, &incy);
#endifarithm
#ifdefarithm L_DEFREAL
   DSWAP(&n, x, &incx, y, &incy);
#endifarithm

}

/******************************************************************************/
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdef NUM_ESSL
#ifdefarithm L_DEFREAL
int Num_dspev_dprimme(int iopt, double *ap, double *w, double *z, int ldz, 
   int n, double *aux, int naux) {

   int ret;

   ret = dspev(iopt, ap, w, z, ldz, n, aux, naux);
   return (ret);
}
#endifarithm
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdefarithm L_DEFCPLX
int Num_zhpev_zprimme(int iopt, Complex_Z *ap, double *w, Complex_Z *z, int ldz,
   int n, Complex_Z *aux, int naux) {

   int ret;

   ret = zhpev(iopt, ap, w, z, ldz, n, aux, naux);
   return (ret);
}
#endifarithm
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#else
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdefarithm L_DEFREAL
void Num_dsyev_dprimme(char *jobz, char *uplo, int n, double *a, int lda, 
   double *w, double *work, int ldwork, int *info) {

#ifdef NUM_CRAY
   _fcd jobz_fcd, uplo_fcd;

   jobz_fcd = _cptofcd(jobz, strlen(jobz));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));

   DSYEV(jobz_fcd, uplo_fcd, &n, a, &lda, w, work, &ldwork, info); 

#else

   DSYEV(jobz, uplo, &n, a, &lda, w, work, &ldwork, info); 

#endif
}
#endifarithm

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdefarithm L_DEFCPLX
void Num_zheev_zprimme(char *jobz, char *uplo, int n, Complex_Z *a, int lda,
   double *w, Complex_Z *work, int ldwork, double *rwork, int *info) {

#ifdef NUM_CRAY
   _fcd jobz_fcd, uplo_fcd;

   jobz_fcd = _cptofcd(jobz, strlen(jobz));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));

   ZHEEV(jobz_fcd, uplo_fcd, &n, a, &lda, w, work, &ldwork, rwork, info); 

#else

   ZHEEV(jobz, uplo, &n, a, &lda, w, work, &ldwork, rwork, info); 

#endif
}
#endifarithm

#endif

#ifdefarithm L_DEFREAL
/******************************************************************************/
void Num_dsytrf_dprimme(char *uplo, int n, double *a, int lda, int *ipivot, 
   double *work, int ldwork, int *info) {

#ifdef NUM_CRAY
   _fcd uplo_fcd;

   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   DSYTRF(uplo_fcd, &n, a, &lda, ipivot, work, &ldwork, info);
#else

   DSYTRF(uplo, &n, a, &lda, ipivot, work, &ldwork, info);

#endif

}
#endifarithm

#ifdefarithm L_DEFCPLX
/******************************************************************************/
void Num_zhetrf_zprimme(char *uplo, int n, @(type) *a, int lda, int *ipivot,
   @(type) *work, int ldwork, int *info) {

#ifdef NUM_CRAY
        _fcd uplo_fcd;

        uplo_fcd = _cptofcd(uplo, strlen(uplo));
        ZHETRF(uplo_fcd, &n, a, &lda, ipivot, work, &ldwork, info);
#else

        ZHETRF(uplo, &n, a, &lda, ipivot, work, &ldwork, info);

#endif

}
#endifarithm

#ifdefarithm L_DEFREAL
/******************************************************************************/
void Num_dsytrs_dprimme(char *uplo, int n, int nrhs, double *a, int lda, 
   int *ipivot, double *b, int ldb, int *info) {

#ifdef NUM_CRAY
   _fcd uplo_fcd;

   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   DSYTRS(uplo_fcd, &n, &nrhs, a, &lda, ipivot, b, &ldb, info);
#else

   DSYTRS(uplo, &n, &nrhs, a, &lda, ipivot, b, &ldb, info);
#endif

}
#endifarithm

#ifdefarithm L_DEFCPLX
/******************************************************************************/
void Num_zhetrs_zprimme(char *uplo, int n, int nrhs, @(type) *a, int lda, 
   int *ipivot, @(type) *b, int ldb, int *info) {

#ifdef NUM_CRAY
        _fcd uplo_fcd;

        uplo_fcd = _cptofcd(uplo, strlen(uplo));
        ZHETRS(uplo_fcd, &n, &nrhs, a, &lda, ipivot, b, &ldb, info);
#else

        ZHETRS(uplo, &n, &nrhs, a, &lda, ipivot, b, &ldb, info);
#endif

}
  
#endifarithm
