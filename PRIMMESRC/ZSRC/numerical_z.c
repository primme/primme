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
#include "Complex.h"
#include "numerical_private_z.h"
#include "numerical_z.h"

/******************************************************************************/
void Num_zcopy_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy) {

   ZCOPY(&n, x, &incx, y, &incy);
}
/******************************************************************************/

void Num_gemm_zprimme(char *transa, char *transb, int m, int n, int k, 
   Complex_Z alpha, Complex_Z *a, int lda, Complex_Z *b, int ldb, 
   Complex_Z beta, Complex_Z *c, int ldc) {

#ifdef NUM_CRAY
   _fcd transa_fcd, transb_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
   transb_fcd = _cptofcd(transb, strlen(transb));
   ZGEMM(transa_fcd, transb_fcd, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, 
         c, &ldc);
#else
   ZGEMM(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif

}

/******************************************************************************/
void Num_symm_zprimme(char *side, char *uplo, int m, int n, Complex_Z alpha, 
   Complex_Z *a, int lda, Complex_Z *b, int ldb, Complex_Z beta, 
   Complex_Z *c, int ldc) {

#ifdef NUM_CRAY
   _fcd side_fcd, uplo_fcd;

   side_fcd = _cptofcd(side, strlen(side));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   ZHEMM(side_fcd, uplo_fcd, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#else
   ZHEMM(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif   

}

/******************************************************************************/
void Num_axpy_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx, 
   Complex_Z *y, int incy) {

   ZAXPY(&n, &alpha, x, &incx, y, &incy);

}

/******************************************************************************/
void Num_gemv_zprimme(char *transa, int m, int n, Complex_Z alpha, Complex_Z *a,
   int lda, Complex_Z *x, int incx, Complex_Z beta, Complex_Z *y, int incy) {

#ifdef NUM_CRAY
   _fcd transa_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
   ZGEMV(transa_fcd, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#else
   ZGEMV(transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

#endif

}

/******************************************************************************/
Complex_Z Num_dot_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy) {

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

}

/******************************************************************************/
void Num_larnv_zprimme(int idist, int *iseed, int length, Complex_Z *x) {
   ZLARNV(&idist, iseed, &length, x);

}

/******************************************************************************/
void Num_scal_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx) {

   ZSCAL(&n, &alpha, x, &incx);

}

/******************************************************************************/
void Num_swap_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy) {

   ZSWAP(&n, x, &incx, y, &incy);

}

/******************************************************************************/
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdef NUM_ESSL
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
int Num_zhpev_zprimme(int iopt, Complex_Z *ap, double *w, Complex_Z *z, int ldz,
   int n, Complex_Z *aux, int naux) {

   int ret;

   ret = zhpev(iopt, ap, w, z, ldz, n, aux, naux);
   return (ret);
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#else
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
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

#endif


/******************************************************************************/
void Num_zhetrf_zprimme(char *uplo, int n, Complex_Z *a, int lda, int *ipivot,
   Complex_Z *work, int ldwork, int *info) {

#ifdef NUM_CRAY
        _fcd uplo_fcd;

        uplo_fcd = _cptofcd(uplo, strlen(uplo));
        ZHETRF(uplo_fcd, &n, a, &lda, ipivot, work, &ldwork, info);
#else

        ZHETRF(uplo, &n, a, &lda, ipivot, work, &ldwork, info);

#endif

}


/******************************************************************************/
void Num_zhetrs_zprimme(char *uplo, int n, int nrhs, Complex_Z *a, int lda, 
   int *ipivot, Complex_Z *b, int ldb, int *info) {

#ifdef NUM_CRAY
        _fcd uplo_fcd;

        uplo_fcd = _cptofcd(uplo, strlen(uplo));
        ZHETRS(uplo_fcd, &n, &nrhs, a, &lda, ipivot, b, &ldb, info);
#else

        ZHETRS(uplo, &n, &nrhs, a, &lda, ipivot, b, &ldb, info);
#endif

}
  
