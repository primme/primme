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
 * File: blaslapack.c
 *
 * Purpose - This file contains mostly C wrapper routines for
 *           calling various BLAS and LAPACK FORTRAN routines.
 *
 ******************************************************************************/

#include <stdlib.h>   /* free */
#include <string.h>   /* memmove */
#include <assert.h>
#include <math.h>
#include "template.h"
#include "blaslapack_private.h"
#include "blaslapack.h"

/*******************************************************************************
 * Subroutine Num_copy_Sprimme - y(0:n*incy-1:incy) = x(0:n*incx-1:incx)
 ******************************************************************************/
 
TEMPLATE_PLEASE
void Num_copy_Sprimme(PRIMME_INT n, SCALAR *x, int incx, SCALAR *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

   while(n > 0) {
      ln = (PRIMME_BLASINT)min(n, PRIMME_BLASINT_MAX-1);
      XCOPY(&ln, x, &lincx, y, &lincy);
      n -= (PRIMME_INT)ln;
      x += ln;
      y += ln;
   }
}

/*******************************************************************************
 * Subroutine Num_gemm_Sprimme - C = op(A)*op(B), with C size m x n
 ******************************************************************************/

TEMPLATE_PLEASE
void Num_gemm_Sprimme(const char *transa, const char *transb, int m, int n,
      int k, SCALAR alpha, SCALAR *a, int lda, SCALAR *b, int ldb, SCALAR beta,
      SCALAR *c, int ldc) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lk = k;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldb = ldb;
   PRIMME_BLASINT lldc = ldc;

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return;

   /* Quick exit */
   if (k == 0) {
      Num_zero_matrix_Sprimme(c, m, n, ldc);
      return;
   }

#ifdef NUM_CRAY
   _fcd transa_fcd, transb_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
   transb_fcd = _cptofcd(transb, strlen(transb));
   XGEMM(transa_fcd, transb_fcd, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, 
         c, &lldc);
#else
   XGEMM(transa, transb, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endif

}

/*******************************************************************************
 * Subroutine Num_gemm_Sprimme - C = A*B or B*A where A is Hermitian,
 *    where C size m x n.
 ******************************************************************************/

TEMPLATE_PLEASE
void Num_hemm_Sprimme(const char *side, const char *uplo, int m, int n,
      SCALAR alpha, SCALAR *a, int lda, SCALAR *b, int ldb, SCALAR beta, 
      SCALAR *c, int ldc) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldb = ldb;
   PRIMME_BLASINT lldc = ldc;

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return;

#ifdef NUM_CRAY
   _fcd side_fcd, uplo_fcd;

   side_fcd = _cptofcd(side, strlen(side));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   XHEMM(side_fcd, uplo_fcd, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#else
   XHEMM(side, uplo, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endif 

}

/*******************************************************************************
 * Subroutine Num_trmm_Sprimme - C = A*B or B*A where A is triangular,
 *    with C size m x n.
 ******************************************************************************/

TEMPLATE_PLEASE
void Num_trmm_Sprimme(const char *side, const char *uplo,
      const char *transa, const char *diag, int m, int n, SCALAR alpha,
      SCALAR *a, int lda, SCALAR *b, int ldb) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldb = ldb;

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return;

#ifdef NUM_CRAY
   _fcd side_fcd, uplo_fcd, transa_fcd, diag_fcd;

   side_fcd = _cptofcd(side, strlen(side));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   transa_fcd = _cptofcd(transa, strlen(transa));
   diag_fcd = _cptofcd(diag, strlen(diag));
   XTRMM(side_fcd, uplo_fcd, transa_fcd, diag_fcd, &lm, &ln, &alpha, a, &llda, b, &lldb);
#else
   XTRMM(side, uplo, transa, diag, &lm, &ln, &alpha, a, &llda, b, &lldb);
#endif

}

/*******************************************************************************
 * Subroutine Num_gemv_Sprimme - y = alpha*A*x + beta*y, with A size m x n
 ******************************************************************************/

TEMPLATE_PLEASE
void Num_gemv_Sprimme(const char *transa, PRIMME_INT m, int n, SCALAR alpha,
      SCALAR *a, int lda, SCALAR *x, int incx, SCALAR beta, SCALAR *y,
      int incy) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

   /* Zero dimension matrix may cause problems */
   if (n == 0) return;

   while(m > 0) {
      lm = (PRIMME_BLASINT)min(m, PRIMME_BLASINT_MAX-1);
#ifdef NUM_CRAY
      _fcd transa_fcd;

      transa_fcd = _cptofcd(transa, strlen(transa));
      XGEMV(transa_fcd, &lm, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);
#else
      XGEMV(transa, &lm, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);
#endif
      m -= (PRIMME_INT)lm;
      a += lm;
      if (transa[0] == 'n' || transa[0] == 'N') {
         y += lm;
      }
      else {
         x += lm;
         beta = 1.0;
      }
   }
}

/*******************************************************************************
 * Subroutine Num_hemv_Sprimme - y = alpha*A*x + beta*y where A is Hermitian
 ******************************************************************************/

TEMPLATE_PLEASE
void Num_hemv_Sprimme(const char *uplo, int n, SCALAR alpha, 
   SCALAR *a, int lda, SCALAR *x, int incx, SCALAR beta, 
   SCALAR *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

   /* Zero dimension matrix may cause problems */
   if (n == 0) return;

#ifdef NUM_CRAY
   _fcd uplo_fcd;

   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   XHEMV(uplo_fcd, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);
#else
   XHEMV(uplo, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);
#endif

}

/*******************************************************************************
 * Subroutine Num_axpy_Sprimme - y += alpha*x
 ******************************************************************************/

TEMPLATE_PLEASE
void Num_axpy_Sprimme(PRIMME_INT n, SCALAR alpha, SCALAR *x, int incx, 
   SCALAR *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

   while(n > 0) {
      ln = (PRIMME_BLASINT)min(n, PRIMME_BLASINT_MAX-1);
      XAXPY(&ln, &alpha, x, &lincx, y, &lincy);
      n -= (PRIMME_INT)ln;
      x += ln;
      y += ln;
   }
}

/*******************************************************************************
 * Subroutine Num_dot_Sprimme - y'*x
 ******************************************************************************/

TEMPLATE_PLEASE
SCALAR Num_dot_Sprimme(PRIMME_INT n, SCALAR *x, int incx, SCALAR *y, int incy) {

/* NOTE: vecLib doesn't follow BLAS reference for sdot */
#if defined(USE_COMPLEX) || (defined(USE_FLOAT) && (defined(__APPLE__) || defined(__MACH__)))
/* ---- Explicit implementation of the zdotc() --- */
   PRIMME_INT i;
   SCALAR zdotc = 0.0;
   if (n <= 0) return(zdotc);
   if (incx == 1 && incy == 1) {
      for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
         zdotc += CONJ(x[i]) * y[i];
      }
   }
   else {
      for (i=0;i<n;i++) { /* zdotc = zdotc + dconjg(x(i))* y(i) */
         zdotc += CONJ(x[i*incx]) * y[i*incy];
      }
   }
   return zdotc;
/* -- end of explicit implementation of the zdotc() - */
#else
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;
   SCALAR r = 0.0;

   while(n > 0) {
      ln = (PRIMME_BLASINT)min(n, PRIMME_BLASINT_MAX-1);
      r += XDOT(&ln, x, &lincx, y, &lincy);
      n -= (PRIMME_INT)ln;
      x += ln;
      y += ln;
   }

   return r;
#endif
}

/*******************************************************************************
 * Subroutine Num_larnv_Sprimme - x(0:n*incy-1:incy) = rand(0:n-1)
 ******************************************************************************/
 
TEMPLATE_PLEASE
void Num_larnv_Sprimme(int idist, PRIMME_INT *iseed, PRIMME_INT length,
      SCALAR *x) {

   PRIMME_BLASINT lidist = idist;
   PRIMME_BLASINT llength;
   PRIMME_BLASINT temp[4];
   PRIMME_BLASINT *liseed;
   int i;

   if (sizeof(PRIMME_INT) == sizeof(PRIMME_BLASINT)) {
      liseed = (PRIMME_BLASINT*)iseed; /* cast avoid compiler warning */
   } else {
      liseed = temp;
      for(i=0; i<4; i++)
         liseed[i] = (PRIMME_BLASINT)iseed[i];
   }

   while(length > 0) {
      llength = (PRIMME_BLASINT)min(length, PRIMME_BLASINT_MAX-1);
      XLARNV(&lidist, liseed, &llength, x);
      length -= (PRIMME_INT)llength;
      x += llength;
   }

   if (sizeof(PRIMME_INT) != sizeof(PRIMME_BLASINT))
      for(i=0; i<4; i++)
         iseed[i] = (int)liseed[i];
}

/*******************************************************************************
 * Subroutine Num_scal_Sprimme - x(0:n*incx-1:incx) *= alpha
 ******************************************************************************/
 
TEMPLATE_PLEASE
void Num_scal_Sprimme(PRIMME_INT n, SCALAR alpha, SCALAR *x, int incx) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;

   while(n > 0) {
      ln = (PRIMME_BLASINT)min(n, PRIMME_BLASINT_MAX-1);
      XSCAL(&ln, &alpha, x, &lincx);
      n -= (PRIMME_INT)ln;
      x += ln;
   }
}

/*******************************************************************************
 * Subroutine Num_swap_Sprimme - swap x(0:n*incx-1:incx) and y(0:n*incy-1:incy)
 ******************************************************************************/
 
TEMPLATE_PLEASE
void Num_swap_Sprimme(PRIMME_INT n, SCALAR *x, int incx, SCALAR *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

   while(n > 0) {
      ln = (PRIMME_BLASINT)min(n, PRIMME_BLASINT_MAX-1);
      XSWAP(&ln, x, &lincx, y, &lincy);
      n -= (PRIMME_INT)ln;
      x += ln;
      y += ln;
   }
}

/*******************************************************************************
 * Subroutines for dense eigenvalue decomposition
 * NOTE: xheevx is used instead of xheev because xheev is not in ESSL
 ******************************************************************************/
 
TEMPLATE_PLEASE
void Num_heev_Sprimme(const char *jobz, const char *uplo, int n, SCALAR *a,
      int lda, REAL *w, SCALAR *work, int ldwork, int *info) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldwork = ldwork;
   PRIMME_BLASINT linfo = 0;
   SCALAR *z;
   REAL abstol=0.0;
#ifdef USE_COMPLEX
   REAL *rwork;
#endif
   PRIMME_BLASINT *iwork, *ifail;
   SCALAR dummys=0;
   REAL   dummyr=0;
   PRIMME_BLASINT dummyi=0;

   /* Zero dimension matrix may cause problems */
   if (n == 0) return;

   /* NULL matrices and zero leading dimension may cause problems */
   if (a == NULL) a = &dummys;
   if (llda < 1) llda = 1;
   if (w == NULL) w = &dummyr;

   /* Borrow space from work for z, rwork and iwork or set dummy values */
   if (ldwork != -1) {
      if (
               WRKSP_MALLOC_PRIMME(n*n, &z, &work, &lldwork) 
#ifdef USE_COMPLEX
            || WRKSP_MALLOC_PRIMME(7*n, &rwork, &work, &lldwork)
#endif
            || WRKSP_MALLOC_PRIMME(5*n, &iwork, &work, &lldwork)
            || WRKSP_MALLOC_PRIMME(n, &ifail, &work, &lldwork)
         ) {
         *info = -1;
         return;
      }
   }
   else {
      z = &dummys;
#ifdef USE_COMPLEX
      rwork = &dummyr;
#endif
      iwork = &dummyi;
      ifail = &dummyi;
   }

#ifdef NUM_CRAY
   _fcd jobz_fcd, range_fcd, uplo_fcd;

   jobz_fcd = _cptofcd(jobz, strlen(jobz));
   range_fcd = _cptofcd("A", 1);
   uplo_fcd = _cptofcd(uplo, strlen(uplo));

   XHEEVX(jobz_fcd, range_fcd, uplo_fcd, &ln, a, &llda, &dummyr, &dummyr,
         &dummyi, &dummyi, &abstol, &dummyi, w, z, &ln, work, &lldwork,
#  ifdef USE_COMPLEX
         rwork,
#  endif
         iwork, ifail, &linfo);
#else
   XHEEVX(jobz, "A", uplo, &ln, a, &llda, &dummyr, &dummyr,
         &dummyi, &dummyi, &abstol, &dummyi, w, z, &ln, work, &lldwork,
#  ifdef USE_COMPLEX
         rwork,
#  endif
         iwork, ifail, &linfo);
#endif

   /* Copy z to a or add the extra space for z, iwork and ifail */
   if (ldwork != -1) {
      Num_copy_matrix_Sprimme(z, n, n, n, a, lda);
   }
   else {
      work[0] += (REAL)n*n + sizeof(PRIMME_BLASINT)*6*n/sizeof(SCALAR) + 6.0;
#ifdef USE_COMPLEX
      work[0] += (REAL)sizeof(REAL)*7*n/sizeof(SCALAR) + 2.0;
#endif
   }
   *info = (int)linfo;
}

/*******************************************************************************
 * Subroutines for dense singular value decomposition
 ******************************************************************************/
 
TEMPLATE_PLEASE
#ifndef USE_COMPLEX
void Num_gesvd_Sprimme(const char *jobu, const char *jobvt, int m, int n,
      SCALAR *a, int lda, REAL *s, SCALAR *u, int ldu, SCALAR *vt, int ldvt,
      SCALAR *work, int ldwork, int *info){

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldu = ldu;
   PRIMME_BLASINT lldvt = ldvt;
   PRIMME_BLASINT lldwork = ldwork;
   PRIMME_BLASINT linfo = 0;
   SCALAR dummys=0;
   REAL   dummyr=0;

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return;

   /* NULL matrices and zero leading dimension may cause problems */
   if (a == NULL) a = &dummys;
   if (llda < 1) llda = 1;
   if (s == NULL) s = &dummyr;
   if (u == NULL) u = &dummys;
   if (lldu < 1) lldu = 1;
   if (vt == NULL) vt = &dummys;
   if (lldvt < 1) lldvt = 1;

#ifdef NUM_CRAY
   _fcd jobu_fcd, jobvt_fcd;

   jobu_fcd = _cptofcd(jobu, strlen(jobu));
   jobvt_fcd = _cptofcd(jobvt, strlen(jobvt));
   XGESVD(jobu_fcd, jobvt_fcd, &lm, &ln, a, &llda, s, u, &lldu, vt, &lldvt, work,
#else
   XGESVD(jobu, jobvt, &lm, &ln, a, &llda, s, u, &lldu, vt, &lldvt, work,
          &lldwork, &linfo);
#endif
   *info = (int)linfo;
}

#else
void Num_gesvd_Sprimme(const char *jobu, const char *jobvt, int m, int n,
   SCALAR *a, int lda, REAL *s, SCALAR *u, int ldu, SCALAR *vt, int ldvt,
   SCALAR *work, int ldwork, REAL *rwork, int *info){

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldu = ldu;
   PRIMME_BLASINT lldvt = ldvt;
   PRIMME_BLASINT lldwork = ldwork;
   PRIMME_BLASINT linfo = 0;
   SCALAR dummys=0;
   REAL   dummyr=0;

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return;

   /* NULL matrices and zero leading dimension may cause problems */
   if (a == NULL) a = &dummys;
   if (llda < 1) llda = 1;
   if (s == NULL) s = &dummyr;
   if (u == NULL) u = &dummys;
   if (lldu < 1) lldu = 1;
   if (vt == NULL) vt = &dummys;
   if (lldvt < 1) lldvt = 1;

#ifdef NUM_CRAY
   _fcd jobu_fcd, jobvt_fcd;

   jobu_fcd = _cptofcd(jobu, strlen(jobu));
   jobvt_fcd = _cptofcd(jobvt, strlen(jobvt));
   XGESVD(jobu_fcd, jobvt_fcd, &lm, &ln, a, &llda, s, u, &lldu, vt, &lldvt, work,
          &lldwork, rwork, &linfo);
#else
   XGESVD(jobu, jobvt, &lm, &ln, a, &llda, s, u, &lldu, vt, &lldvt, work,
          &lldwork, rwork, &linfo);
#endif
   *info = (int)linfo;
}
#endif

/*******************************************************************************
 * Subroutine Num_hetrf_Sprimme - LL^H factorization with pivoting
 ******************************************************************************/

TEMPLATE_PLEASE
void Num_hetrf_Sprimme(const char *uplo, int n, SCALAR *a, int lda, int *ipivot,
   SCALAR *work, int ldwork, int *info) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT *lipivot;
   PRIMME_BLASINT lldwork = ldwork;
   PRIMME_BLASINT linfo = 0; 
   int i;
   SCALAR dummys=0;
   PRIMME_BLASINT dummyi=0;

   /* Zero dimension matrix may cause problems */
   if (n == 0) return;

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      if (MALLOC_PRIMME(n, &lipivot) != 0) {
         *info = -1;
         return;
      }
   } else {
      lipivot = (PRIMME_BLASINT *)ipivot; /* cast avoid compiler warning */
   }

   /* NULL matrices and zero leading dimension may cause problems */
   if (a == NULL) a = &dummys;
   if (llda < 1) llda = 1;
   if (lipivot == NULL) lipivot = &dummyi;

#ifdef NUM_CRAY
   _fcd uplo_fcd;

   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   XHETRF(uplo_fcd, &ln, a, &llda, lipivot, work, &lldwork, &linfo);
#else
   XHETRF(uplo, &ln, a, &llda, lipivot, work, &lldwork, &linfo);
#endif

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      if (ipivot) for(i=0; i<n; i++)
         ipivot[i] = (int)lipivot[i];
      free(lipivot);
   }
   *info = (int)linfo;

}

/*******************************************************************************
 * Subroutine Num_hetrs_Sprimme - b = A\b where A stores a LL^H factorization
 ******************************************************************************/
 
TEMPLATE_PLEASE
void Num_hetrs_Sprimme(const char *uplo, int n, int nrhs, SCALAR *a,
      int lda, int *ipivot, SCALAR *b, int ldb, int *info) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lnrhs = nrhs;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT *lipivot;
   PRIMME_BLASINT lldb = ldb;
   PRIMME_BLASINT linfo = 0; 
   int i;

   /* Zero dimension matrix may cause problems */
   if (n == 0 || nrhs == 0) return;

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      if (MALLOC_PRIMME(n, &lipivot) != 0) {
         *info = -1;
         return;
      }
   } else {
      lipivot = (PRIMME_BLASINT *)ipivot; /* cast avoid compiler warning */
   }

#ifdef NUM_CRAY
   _fcd uplo_fcd;

   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   XHETRS(uplo_fcd, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
#else
   XHETRS(uplo, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
#endif

   if (sizeof(int) != sizeof(PRIMME_BLASINT)) {
      if (ipivot) for(i=0; i<n; i++)
         ipivot[i] = (int)lipivot[i];
      free(lipivot);
   }
   *info = (int)linfo;
}

/*******************************************************************************
 * Subroutine Num_trsm_Sprimme - b = op(A)\b
 ******************************************************************************/
 
TEMPLATE_PLEASE
void Num_trsm_Sprimme(const char *side, const char *uplo, const char *transa,
      const char *diag, int m, int n, SCALAR alpha, SCALAR *a, int lda,
      SCALAR *b, int ldb) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldb = ldb;

   /* Zero dimension matrix may cause problems */
   if (m == 0 || n == 0) return;

#ifdef NUM_CRAY
   _fcd side_fcd, uplo_fcd, transa_fcd, diag_fcd;

   side_fcd = _cptofcd(side, strlen(side));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   transa_fcd = _cptofcd(transa, strlen(transa));
   diag_fcd = _cptofcd(diag, strlen(diag));
   XTRSM(side_fcd, uplo_fcd, transa_fcd, diag_fcd, &lm, &ln, &alpha, a, &llda, b, &lldb);
#else
   XTRSM(side, uplo, transa, diag, &lm, &ln, &alpha, a, &llda, b, &lldb);
#endif
}
