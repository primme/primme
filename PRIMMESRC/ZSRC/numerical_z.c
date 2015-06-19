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

#include "primme.h"
#include <stdarg.h>
/*#include "Complex.h"*/
#include "numerical_private_z.h"
#include "numerical_z.h"

/******************************************************************************/
void Num_zcopy_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy) {
   long long int ln = n;
   long long int lincx = incx;
   long long int lincy = incy;

   ZCOPY(&ln, x, &lincx, y, &lincy);
}
/******************************************************************************/

void Num_gemm_zprimme(char *transa, char *transb, int m, int n, int k, 
   Complex_Z alpha, Complex_Z *a, int lda, Complex_Z *b, int ldb, 
   Complex_Z beta, Complex_Z *c, int ldc) {
   long long int lm = m;
   long long int ln = n;
   long long int lk = k;
   long long int llda = lda;
   long long int lldb = ldb;
   long long int lldc = ldc;

#ifdef NUM_CRAY
   _fcd transa_fcd, transb_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
   transb_fcd = _cptofcd(transb, strlen(transb));
   ZGEMM(transa_fcd, transb_fcd, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, 
	 c, &lldc);
#else
   ZGEMM(transa, transb, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endif

}

/******************************************************************************/
void Num_symm_zprimme(char *side, char *uplo, int m, int n, Complex_Z alpha, 
   Complex_Z *a, int lda, Complex_Z *b, int ldb, Complex_Z beta, 
   Complex_Z *c, int ldc) {
   long long int lm = m;
   long long int ln = n;
   long long int llda = lda;
   long long int lldb = ldb;
   long long int lldc = ldc;

#ifdef NUM_CRAY
   _fcd side_fcd, uplo_fcd;

   side_fcd = _cptofcd(side, strlen(side));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   ZHEMM(side_fcd, uplo_fcd, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#else
   ZHEMM(side, uplo, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endif   

}

/******************************************************************************/
void Num_axpy_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx, 
   Complex_Z *y, int incy) {
   long long int ln = n;
   long long int lincx = incx;
   long long int lincy = incy;

   ZAXPY(&ln, &alpha, x, &lincx, y, &lincy);

}

/******************************************************************************/
void Num_gemv_zprimme(char *transa, int m, int n, Complex_Z alpha, Complex_Z *a,
   int lda, Complex_Z *x, int incx, Complex_Z beta, Complex_Z *y, int incy) {
   long long int lm = m;
   long long int ln = n;
   long long int llda = lda;
   long long int lincx = incx;
   long long int lincy = incy;

#ifdef NUM_CRAY
   _fcd transa_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
   ZGEMV(transa_fcd, &lm, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);
#else
   ZGEMV(transa, &lm, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);

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
   long long int lidist = idist;
   long long int llength = length;
   long long int *liseed;
   long long int temp[4];
   int i;
   for(i=0;i<4;i++)
   {
      temp[i] = iseed[i];
   }
   liseed = temp;

   ZLARNV(&lidist, liseed, &llength, x);

   for(i=0;i<4;i++)
   {
      iseed[i] = liseed[i];
   }

}

/******************************************************************************/
void Num_scal_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx) {
   long long int ln = n;
   long long int lincx = incx;

   ZSCAL(&ln, &alpha, x, &lincx);

}

/******************************************************************************/
void Num_swap_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy) {
   long long int ln = n;
   long long int lincx = incx;
   long long int lincy = incy;

   ZSWAP(&ln, x, &lincx, y, &lincy);

}

/******************************************************************************/
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef NUM_ESSL
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int Num_zhpev_zprimme(int iopt, Complex_Z *ap, double *w, Complex_Z *z, int ldz,
   int n, Complex_Z *aux, int naux) {
   long long int liopt = iopt;
   long long int lldz = ldz;
   long long int ln = n;
   long long int lnaux = naux;

   int ret;

   ret = zhpev(liopt, ap, w, z, lldz, ln, aux, lnaux);
   return (ret);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Num_zheev_zprimme(char *jobz, char *uplo, int n, Complex_Z *a, int lda,
   double *w, Complex_Z *work, int ldwork, double *rwork, int *info) {
   long long int ln = n;
   long long int llda = lda;
   long long int lldwork = ldwork;
   long long int linfo = (long long int)*info;

#ifdef NUM_CRAY
   _fcd jobz_fcd, uplo_fcd;

   jobz_fcd = _cptofcd(jobz, strlen(jobz));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));

   ZHEEV(jobz_fcd, uplo_fcd, &ln, a, &llda, w, work, &lldwork, rwork, &linfo); 

#else

   ZHEEV(jobz, uplo, &ln, a, &llda, w, work, &lldwork, rwork, &linfo); 

#endif
    *info = (int) linfo;
}

#endif


/******************************************************************************/
void Num_zhetrf_zprimme(char *uplo, int n, Complex_Z *a, int lda, int *ipivot,
   Complex_Z *work, int ldwork, int *info) {
   long long int ln = n;
   long long int llda = lda;

   long long int *lipivot = (long long int *)primme_calloc(n, sizeof(long long int), "lipivot array");
   int i;
   for(i=0;i<n;i++)
   {
      lipivot[i] = ipivot[i];
   }
   
   long long int lldwork = ldwork;
   long long int linfo = (long long int)*info;

#ifdef NUM_CRAY
	_fcd uplo_fcd;

	uplo_fcd = _cptofcd(uplo, strlen(uplo));
	ZHETRF(uplo_fcd, &ln, a, &llda, lipivot, work, &lldwork, &linfo);
#else

	ZHETRF(uplo, &ln, a, &llda, lipivot, work, &lldwork, &linfo);

#endif

   for(i=0;i<n;i++)
   {
      ipivot[i] = lipivot[i];
   }
    *info = (int) linfo;

}


/******************************************************************************/
void Num_zhetrs_zprimme(char *uplo, int n, int nrhs, Complex_Z *a, int lda, 
   int *ipivot, Complex_Z *b, int ldb, int *info) {
   long long int ln = n;
   long long int lnrhs = nrhs;
   long long int llda = lda;

   long long int *lipivot = (long long int *)primme_calloc(n, sizeof(long long int), "lipivot array");
   int i;
   for(i=0;i<n;i++)
   {
      lipivot[i] = ipivot[i];
   }
   
   long long int lldb = ldb;
   long long int linfo = (long long int)*info;

#ifdef NUM_CRAY
	_fcd uplo_fcd;

	uplo_fcd = _cptofcd(uplo, strlen(uplo));
	ZHETRS(uplo_fcd, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
#else

	ZHETRS(uplo, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
#endif

    for(i=0;i<n;i++)
   {
      ipivot[i] = lipivot[i];
   }
    *info = (int) linfo;

}
  
