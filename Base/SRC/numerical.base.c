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
 * File: numerical.c
 *
 * Purpose - This file contains mostly C wrapper routines for
 *           calling various BLAS and LAPACK FORTRAN routines.
 *
 ******************************************************************************/


#include <stdarg.h>
#ifdefarithm L_DEFCPLX
#include "Complexz.h"
#endifarithm
#include "numerical_private_@(pre).h"
#include "numerical_@(pre).h"
#include "primme.h"
#include <stdlib.h>   /* free */

/******************************************************************************/
void Num_@(pre)copy_@(pre)primme(int n, @(type) *x, int incx, @(type) *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

#ifdefarithm L_DEFCPLX
   ZCOPY(&ln, x, &lincx, y, &lincy);
#endifarithm
#ifdefarithm L_DEFREAL
   DCOPY(&ln, x, &lincx, y, &lincy);
#endifarithm
}
/******************************************************************************/

void Num_gemm_@(pre)primme(const char *transa, const char *transb, int m, int n, int k, 
   @(type) alpha, @(type) *a, int lda, @(type) *b, int ldb, 
   @(type) beta, @(type) *c, int ldc) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lk = k;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldb = ldb;
   PRIMME_BLASINT lldc = ldc;


#ifdef NUM_CRAY
   _fcd transa_fcd, transb_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
   transb_fcd = _cptofcd(transb, strlen(transb));
#ifdefarithm L_DEFCPLX
   ZGEMM(transa_fcd, transb_fcd, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, 
         c, &lldc);
#endifarithm
#ifdefarithm L_DEFREAL
   DGEMM(transa_fcd, transb_fcd, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, 
         c, &lldc);
#endifarithm
#else
#ifdefarithm L_DEFCPLX
   ZGEMM(transa, transb, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endifarithm
#ifdefarithm L_DEFREAL
   DGEMM(transa, transb, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endifarithm
#endif

}

/******************************************************************************/
void Num_symm_@(pre)primme(const char *side, const char *uplo, int m, int n, @(type) alpha, 
   @(type) *a, int lda, @(type) *b, int ldb, @(type) beta, 
   @(type) *c, int ldc) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldb = ldb;
   PRIMME_BLASINT lldc = ldc;

#ifdef NUM_CRAY
   _fcd side_fcd, uplo_fcd;

   side_fcd = _cptofcd(side, strlen(side));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));
#ifdefarithm L_DEFCPLX
   ZHEMM(side_fcd, uplo_fcd, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endifarithm
#ifdefarithm L_DEFREAL
   DSYMM(side_fcd, uplo_fcd, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endifarithm
#else
#ifdefarithm L_DEFCPLX
   ZHEMM(side, uplo, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endifarithm
#ifdefarithm L_DEFREAL
   DSYMM(side, uplo, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endifarithm
#endif 

}

/******************************************************************************/
void Num_axpy_@(pre)primme(int n, @(type) alpha, @(type) *x, int incx, 
   @(type) *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

#ifdefarithm L_DEFCPLX
   ZAXPY(&ln, &alpha, x, &lincx, y, &lincy);
#endifarithm
#ifdefarithm L_DEFREAL
   DAXPY(&ln, &alpha, x, &lincx, y, &lincy);
#endifarithm

}

/******************************************************************************/
void Num_gemv_@(pre)primme(const char *transa, int m, int n, @(type) alpha, @(type) *a,
   int lda, @(type) *x, int incx, @(type) beta, @(type) *y, int incy) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

#ifdef NUM_CRAY
   _fcd transa_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
#ifdefarithm L_DEFCPLX
   ZGEMV(transa_fcd, &lm, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);
#endifarithm
#ifdefarithm L_DEFREAL
   DGEMV(transa_fcd, &lm, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);
#endifarithm
#else
#ifdefarithm L_DEFCPLX
   ZGEMV(transa, &lm, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);
#endifarithm
#ifdefarithm L_DEFREAL
   DGEMV(transa, &lm, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);
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
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

   return(DDOT(&ln, x, &lincx, y, &lincy));
#endifarithm

}

/******************************************************************************/
void Num_larnv_@(pre)primme(int idist, int *iseed, int length, @(type) *x) {

   PRIMME_BLASINT lidist = idist;
   PRIMME_BLASINT llength = length;
   PRIMME_BLASINT temp[4];
   PRIMME_BLASINT *liseed = temp;
   int i;

   if (sizeof(int) == sizeof(PRIMME_BLASINT)) {
      liseed = (PRIMME_BLASINT*)iseed; /* cast avoid compiler warning */
   } else {
      liseed = temp;
      for(i=0; i<4; i++)
         liseed[i] = (PRIMME_BLASINT)iseed[i];
   }

#ifdefarithm L_DEFCPLX
   ZLARNV(&lidist, liseed, &llength, x);
#endifarithm
#ifdefarithm L_DEFREAL
   DLARNV(&lidist, liseed, &llength, x);
#endifarithm

   if (sizeof(int) != sizeof(PRIMME_BLASINT))
      for(i=0; i<4; i++)
         iseed[i] = (int)liseed[i];

}

/******************************************************************************/
void Num_scal_@(pre)primme(int n, @(type) alpha, @(type) *x, int incx) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;

#ifdefarithm L_DEFCPLX
   ZSCAL(&ln, &alpha, x, &lincx);
#endifarithm
#ifdefarithm L_DEFREAL
   DSCAL(&ln, &alpha, x, &lincx);
#endifarithm

}

/******************************************************************************/
void Num_swap_@(pre)primme(int n, @(type) *x, int incx, @(type) *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

#ifdefarithm L_DEFCPLX
   ZSWAP(&ln, x, &lincx, y, &lincy);
#endifarithm
#ifdefarithm L_DEFREAL
   DSWAP(&ln, x, &lincx, y, &lincy);
#endifarithm

}

/******************************************************************************/
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdef NUM_ESSL
#ifdefarithm L_DEFREAL
int Num_dspev_dprimme(int iopt, double *ap, double *w, double *z, int ldz, 
   int n, double *aux, int naux) {

   PRIMME_BLASINT liopt = iopt;
   PRIMME_BLASINT lldz = ldz;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lnaux = naux;

   return dspev(liopt, ap, w, z, lldz, ln, aux, lnaux);
}
#endifarithm
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdefarithm L_DEFCPLX
int Num_zhpev_zprimme(int iopt, Complex_Z *ap, double *w, Complex_Z *z, int ldz,
   int n, Complex_Z *aux, int naux) {

   PRIMME_BLASINT liopt = iopt;
   PRIMME_BLASINT lldz = ldz;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lnaux = naux;

   return zhpev(liopt, ap, w, z, lldz, ln, aux, lnaux);
}
#endifarithm
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#else
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdefarithm L_DEFREAL
void Num_dsyev_dprimme(const char *jobz, const char *uplo, int n, double *a, int lda, 
   double *w, double *work, int ldwork, int *info) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldwork = ldwork;
   PRIMME_BLASINT linfo = 0;
 
#ifdef NUM_CRAY
   _fcd jobz_fcd, uplo_fcd;

   jobz_fcd = _cptofcd(jobz, strlen(jobz));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));

   DSYEV(jobz_fcd, uplo_fcd, &ln, a, &llda, w, work, &lldwork, &linfo);

#else

   DSYEV(jobz, uplo, &ln, a, &llda, w, work, &lldwork, &linfo);

#endif
   *info = linfo;
}
#endifarithm

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdefarithm L_DEFCPLX
void Num_zheev_zprimme(const char *jobz, const char *uplo, int n, Complex_Z *a, int lda,
   double *w, Complex_Z *work, int ldwork, double *rwork, int *info) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldwork = ldwork;
   PRIMME_BLASINT linfo = 0;

#ifdef NUM_CRAY
   _fcd jobz_fcd, uplo_fcd;

   jobz_fcd = _cptofcd(jobz, strlen(jobz));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));

   ZHEEV(jobz_fcd, uplo_fcd, &ln, a, &llda, w, work, &lldwork, rwork, &linfo); 

#else

   ZHEEV(jobz, uplo, &ln, a, &llda, w, work, &lldwork, rwork, &linfo);

#endif
   *info = linfo;
}
#endifarithm

#endif

#ifdefarithm L_DEFREAL
/******************************************************************************/
void Num_dsytrf_dprimme(const char *uplo, int n, double *a, int lda, int *ipivot, 
   double *work, int ldwork, int *info) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT *lipivot;
   PRIMME_BLASINT lldwork = ldwork;
   PRIMME_BLASINT linfo = 0; 
   int i;

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      lipivot = (PRIMME_BLASINT *)primme_calloc(n, sizeof(PRIMME_BLASINT), "lipivot array");
   } else {
      lipivot = (PRIMME_BLASINT *)ipivot; /* cast avoid compiler warning */
   }

#ifdef NUM_CRAY
   _fcd uplo_fcd;

   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   DSYTRF(uplo_fcd, &ln, a, &llda, lipivot, work, &lldwork, &linfo);
#else

   DSYTRF(uplo, &ln, a, &llda, lipivot, work, &lldwork, &linfo);

#endif

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      for(i=0; i<n; i++)
         ipivot[i] = (int)lipivot[i];
      free(lipivot);
   }
   *info = (int)linfo;

}
#endifarithm

#ifdefarithm L_DEFCPLX
/******************************************************************************/
void Num_zhetrf_zprimme(const char *uplo, int n, @(type) *a, int lda, int *ipivot,
   @(type) *work, int ldwork, int *info) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT *lipivot;
   PRIMME_BLASINT lldwork = ldwork;
   PRIMME_BLASINT linfo = 0; 
   int i;

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      lipivot = (PRIMME_BLASINT *)primme_calloc(n, sizeof(PRIMME_BLASINT), "lipivot array");
   } else {
      lipivot = (PRIMME_BLASINT *)ipivot; /* cast avoid compiler warning */
   }

#ifdef NUM_CRAY
        _fcd uplo_fcd;

        uplo_fcd = _cptofcd(uplo, strlen(uplo));
        ZHETRF(uplo_fcd, &ln, a, &llda, lipivot, work, &lldwork, &linfo);
#else

        ZHETRF(uplo, &ln, a, &llda, lipivot, work, &lldwork, &linfo);

#endif

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      for(i=0; i<n; i++)
         ipivot[i] = (int)lipivot[i];
      free(lipivot);
   }
   *info = (int)linfo;

}
#endifarithm

#ifdefarithm L_DEFREAL
/******************************************************************************/
void Num_dsytrs_dprimme(const char *uplo, int n, int nrhs, double *a, int lda, 
   int *ipivot, double *b, int ldb, int *info) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lnrhs = nrhs;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT *lipivot;
   PRIMME_BLASINT lldb = ldb;
   PRIMME_BLASINT linfo = 0; 
   int i;

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      lipivot = (PRIMME_BLASINT *)primme_calloc(n, sizeof(PRIMME_BLASINT), "lipivot array");
   } else {
      lipivot = (PRIMME_BLASINT *)ipivot; /* cast avoid compiler warning */
   }

#ifdef NUM_CRAY
   _fcd uplo_fcd;

   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   DSYTRS(uplo_fcd, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
#else

   DSYTRS(uplo, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
#endif

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      for(i=0; i<n; i++)
         ipivot[i] = (int)lipivot[i];
      free(lipivot);
   }
   *info = (int)linfo;

}
#endifarithm

#ifdefarithm L_DEFCPLX
/******************************************************************************/
void Num_zhetrs_zprimme(const char *uplo, int n, int nrhs, @(type) *a, int lda, 
   int *ipivot, @(type) *b, int ldb, int *info) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lnrhs = nrhs;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT *lipivot;
   PRIMME_BLASINT lldb = ldb;
   PRIMME_BLASINT linfo = 0; 
   int i;

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      lipivot = (PRIMME_BLASINT *)primme_calloc(n, sizeof(PRIMME_BLASINT), "lipivot array");
   } else {
      lipivot = (PRIMME_BLASINT *)ipivot; /* cast avoid compiler warning */
   }

#ifdef NUM_CRAY
        _fcd uplo_fcd;

        uplo_fcd = _cptofcd(uplo, strlen(uplo));
        ZHETRS(uplo_fcd, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
#else

        ZHETRS(uplo, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
#endif

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      for(i=0; i<n; i++)
         ipivot[i] = (int)lipivot[i];
      free(lipivot);
   }
   *info = (int)linfo;

}
 
#endifarithm
