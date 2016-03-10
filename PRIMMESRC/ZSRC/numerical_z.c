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
#include "Complexz.h"
#include "numerical_private_z.h"
#include "numerical_z.h"
#include "primme.h"
#include <stdlib.h>   /* free */
#include <string.h>   /* memmove */
#include <assert.h>
#include <math.h>

/******************************************************************************/
void Num_zcopy_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

   ZCOPY(&ln, x, &lincx, y, &lincy);
}
/******************************************************************************/

void Num_gemm_zprimme(const char *transa, const char *transb, int m, int n, int k, 
   Complex_Z alpha, Complex_Z *a, int lda, Complex_Z *b, int ldb, 
   Complex_Z beta, Complex_Z *c, int ldc) {

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
   ZGEMM(transa_fcd, transb_fcd, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, 
         c, &lldc);
#else
   ZGEMM(transa, transb, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endif

}

/******************************************************************************/
void Num_symm_zprimme(const char *side, const char *uplo, int m, int n, Complex_Z alpha, 
   Complex_Z *a, int lda, Complex_Z *b, int ldb, Complex_Z beta, 
   Complex_Z *c, int ldc) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldb = ldb;
   PRIMME_BLASINT lldc = ldc;

#ifdef NUM_CRAY
   _fcd side_fcd, uplo_fcd;

   side_fcd = _cptofcd(side, strlen(side));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   ZHEMM(side_fcd, uplo_fcd, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#else
   ZHEMM(side, uplo, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endif 

}

void Num_trmm_zprimme(const char *side, const char *uplo, const char *transa,
   const char *diag, int m, int n, Complex_Z alpha, Complex_Z *a, int lda, Complex_Z *b,
   int ldb) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldb = ldb;

   ZTRMM(side, uplo, transa, diag, &lm, &ln, &alpha, a, &llda, b, &lldb);

}


/******************************************************************************/
void Num_symv_zprimme(const char *uplo, int n, Complex_Z alpha, 
   Complex_Z *a, int lda, Complex_Z *x, int incx, Complex_Z beta, 
   Complex_Z *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;


   ZHEMV(uplo, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);

}

/******************************************************************************/
void Num_axpy_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx, 
   Complex_Z *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

   ZAXPY(&ln, &alpha, x, &lincx, y, &lincy);

}

/******************************************************************************
 * Function Num_compute_residual - This subroutine performs the next operation
 *    in a cache-friendly way:
 *
 *    r = Ax - eval*x
 *
 * PARAMETERS
 * ---------------------------
 * n           The number of rows of x, Ax and r
 * eval        The value to compute the residual vector r
 * x           The vector x
 * Ax          The vector Ax
 * r           On output r = Ax - eval*x
 *
 ******************************************************************************/

void Num_compute_residual_zprimme(int n, double eval, Complex_Z *x, 
   Complex_Z *Ax, Complex_Z *r) {

   Complex_Z ztmp = {+0.0e+00,+0.0e00};
   int k, M=min(n,PRIMME_BLOCK_SIZE);
   *(double*)&ztmp = -eval;

   for (k=0; k<n; k+=M, M=min(M,n-k)) {
      Num_zcopy_zprimme(M, &Ax[k], 1, &r[k], 1);
      Num_axpy_zprimme(M, ztmp, &x[k], 1, &r[k], 1);
   }

}

/******************************************************************************
 * Function Num_compute_residual_i - This subroutine performs the next operations
 *    in a cache-friendly way:
 *
 *    X = X(p); Ax = Ax(p)
 *    j = k = 0; XD = RD = []
 *    for i=0:nd-1
 *       if pd(i) == j
 *          XD = [XD XO(j)]; RD = [RD RO(j)]; j++
 *       else
 *          XD = [XD X(p(k)]; RD = [RD AX(p(k)) - evals(p(k))*X(p(k))]; k++
 *       end if
 *    end for
 *
 * NOTE: X and XD *can* overlap, but X(0:n-1) and XD *cannot* overlap (same for R and RD)
 *       XO and XD *can* overlap (same for RO and RD)
 *       p should be a list of increasing indices
 *       pd should be a merge of two increasing lists
 *
 * PARAMETERS
 * ---------------------------
 * m           The number of rows of matrices x, Ax, xo, ro, xd and rd
 * evals       The values to compute the residual
 * x           The matrix that does x = x(p)
 * n           The number of columns of the output x
 * p           The columns to copy back to x and Ax
 * ldx         The leading dimension of x
 * Ax          The matrix that does Ax = Ax(p)
 * ldAx        The leading dimension of Ax
 * xo          Alternative source of columns for xd
 * no          The number of columns in xo
 * ldxo        The leading dimension of xo
 * ro          Alternative source of columns for rd
 * ldro        The leading dimension of ro
 * xd          The matrix that will have columns from x and xo
 * nd          The maximum size of xd
 * pd          The indices of the columns to generate xd and rd
 * ldxd        The leading dimension of xd
 * rd          The matrix that will have columns from r and ro
 * ldrd        The leading dimension of rd
 * rwork       Workspace
 * lrwork      The size of rwork
 *
 ******************************************************************************/

int Num_compute_residual_i_zprimme(int m, double *evals, Complex_Z *x, int n, int *p, 
   int ldx, Complex_Z *Ax, int ldAx,
   Complex_Z *xo, int no, int ldxo, Complex_Z *ro, int ldro,
   Complex_Z *xd, int nd, int *pd, int ldxd, Complex_Z *rd, int ldrd,
   Complex_Z *rwork, int lrwork) {

   int i, id, k, io, M=min(m,PRIMME_BLOCK_SIZE);
   Complex_Z *X0, *R0;

   /* Return memory requirement */

   if (evals == NULL) {
      return nd*M*2;
   }

   /* Quick exit */

   if (n == 0) {
      Num_copy_matrix_zprimme(xo, m, min(no,nd), ldxo, xd, ldxd);
      Num_copy_matrix_zprimme(ro, m, min(no,nd), ldro, rd, ldrd);
      return 0;
   }

   X0 = rwork;
   R0 = X0+nd*M;
   assert(nd*M*2 <= lrwork);

   for (k=0; k<m; k+=M, M=min(M,m-k)) {
      for (i=id=io=0; i < n || id < nd; id++) {
         if (id < nd && io < no && pd[id] == io) {
            Num_copy_matrix_zprimme(&xo[io*ldxo+k], M, 1, ldxo, &X0[id*M], M);
            Num_copy_matrix_zprimme(&ro[io*ldro+k], M, 1, ldro, &R0[id*M], M);
            io++;
         }
         else {
            assert(id >= nd || i < n);
            Num_copy_matrix_zprimme(&x[p[i]*ldx+k],   M, 1, ldx,  &x[i*ldx +k],  ldx);
            Num_copy_matrix_zprimme(&Ax[p[i]*ldAx+k], M, 1, ldAx, &Ax[i*ldAx+k], ldAx);
            if (id < nd) {
               Num_copy_matrix_zprimme(&x[p[i]*ldx+k], M, 1, ldx, &X0[id*M], M);
               Num_compute_residual_zprimme(M, evals[p[i]], &x[p[i]*ldx+k], &Ax[p[i]*ldAx+k], &R0[id*M]);
            }
            i++;
         }
      }
      assert(id >= nd);
      Num_copy_matrix_zprimme(X0, M, nd, M, &xd[k], ldxd);
      Num_copy_matrix_zprimme(R0, M, nd, M, &rd[k], ldrd);
   }

   return 0;
}

/******************************************************************************/
void Num_gemv_zprimme(const char *transa, int m, int n, Complex_Z alpha, Complex_Z *a,
   int lda, Complex_Z *x, int incx, Complex_Z beta, Complex_Z *y, int incy) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

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

   ZLARNV(&lidist, liseed, &llength, x);

   if (sizeof(int) != sizeof(PRIMME_BLASINT))
      for(i=0; i<4; i++)
         iseed[i] = (int)liseed[i];

}

/******************************************************************************/
void Num_scal_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;

   ZSCAL(&ln, &alpha, x, &lincx);

}

/******************************************************************************/
void Num_swap_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy) {

   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

   ZSWAP(&ln, x, &lincx, y, &lincy);

}

/******************************************************************************/
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#ifdef NUM_ESSL
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
int Num_zhpev_zprimme(int iopt, Complex_Z *ap, double *w, Complex_Z *z, int ldz,
   int n, Complex_Z *aux, int naux) {

   PRIMME_BLASINT liopt = iopt;
   PRIMME_BLASINT lldz = ldz;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lnaux = naux;

   return zhpev(liopt, ap, w, z, lldz, ln, aux, lnaux);
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
#else
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
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
   *info = (int)linfo;
}

#endif



/******************************************************************************/
void Num_zgesvd_zprimme(const char *jobu, const char *jobvt, int m, int n, Complex_Z *a,
   int lda, double *s, Complex_Z *u, int ldu, Complex_Z *vt, int ldvt,
   Complex_Z *work, int ldwork, double *rwork, int *info){

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT lldu = ldu;
   PRIMME_BLASINT lldvt = ldvt;
   PRIMME_BLASINT lldwork = ldwork;
   PRIMME_BLASINT linfo = 0;

   ZGESVD(jobu, jobvt, &lm, &ln, a, &llda, s, u, &lldu, vt, &lldvt, work,
          &lldwork, rwork, &linfo);

   *info = (int)linfo;
}

/******************************************************************************/

void Num_geqrf_zprimme(int m, int n, Complex_Z *a, int lda, Complex_Z *tau,
      Complex_Z *rwork, int lrwork, int *info) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT llrwork = lrwork;
   PRIMME_BLASINT linfo = 0;

   ZGEQRF
      (&lm, &ln, a, &llda, tau, rwork, &llrwork, &linfo);

   *info = (int)linfo;
}

void Num_orgqr_zprimme(int m, int n, int k, Complex_Z *a, int lda, Complex_Z *tau,
      Complex_Z *rwork, int lrwork, int *info) {

   PRIMME_BLASINT lm = m;
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lk = k;
   PRIMME_BLASINT llda = lda;
   PRIMME_BLASINT llrwork = lrwork;
   PRIMME_BLASINT linfo = 0;

   ZUNGQR(&lm, &ln, &lk, a, &llda, tau, rwork, &llrwork, &linfo);

   *info = (int)linfo;
}



/******************************************************************************/
void Num_zhetrf_zprimme(const char *uplo, int n, Complex_Z *a, int lda, int *ipivot,
   Complex_Z *work, int ldwork, int *info) {

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


/******************************************************************************/
void Num_zhetrs_zprimme(const char *uplo, int n, int nrhs, Complex_Z *a, int lda, 
   int *ipivot, Complex_Z *b, int ldb, int *info) {

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
 

/******************************************************************************
 * Function Num_copy_matrix - Copy the matrix x into y
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * y           On output y = x
 * ldy         The leading dimension of y
 *
 * NOTE: x and y *can* overlap
 *
 ******************************************************************************/

void Num_copy_matrix_zprimme(Complex_Z *x, int m, int n, int ldx, Complex_Z *y, int ldy) {
   int i,j;

   assert(ldx >= m && ldy >= m);

   /* Do nothing if x and y are the same matrix */
   if (x == y && ldx == ldy) return;

   /* Copy a contiguous memory region */
   if (ldx == ldy && ldx == m) {
      memmove(y, x, sizeof(Complex_Z)*m*n);
   }

   /* Copy matrix some rows down or up */
   else if (ldx == ldy && (y > x ? y-x : x-y) < ldx) {
      for (i=0; i<n; i++)
         memmove(&y[i*ldy], &x[i*ldx], sizeof(Complex_Z)*m);
   }

   /* Copy matrix some columns forward */
   else if (ldx == ldy && y > x && y-x > ldx) {
      for (i=n-1; i>=0; i--)
         for (j=0; j<m; j++)
            y[i*ldy+j] = x[i*ldx+j];
   }

   /* Copy matrix some columns backward and the general case */
   else {
      /* TODO: assert x and y don't overlap */
      for (i=0; i<n; i++)
         for (j=0; j<m; j++)
            y[i*ldy+j] = x[i*ldx+j];
   }

}

/******************************************************************************
 * Function Num_copy_matrix_i - Copy the matrix x(xin) into y(yin)
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * m           The number of rows of x
 * xin         The column indices to copy
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * y           On output y(yin) = x(xin)
 * yin         The column indices of y to be modified
 * ldy         The leading dimension of y
 *
 * NOTE: x(xin) and y(yin) *cannot* overlap
 *
 ******************************************************************************/

void Num_copy_matrix_i_zprimme(Complex_Z *x, int m, int *xin, int n, int ldx, Complex_Z *y,
      int *yin, int ldy) {

   int i,j;

   /* TODO: assert x and y don't overlap */
   for (i=0; i<n; i++)
      for (j=0; j<m; j++)
         y[(yin?yin[i]:i)*ldy+j] = x[(xin?xin[i]:i)*ldx+j];
} 


/******************************************************************************
 * Function Num_copy_trimatrix - Copy the upper/lower triangular part of the
 *    matrix x into y
 *
 * PARAMETERS
 * ---------------------------
 * x           The source matrix
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * ul          if 0, copy the upper part; otherwise copy the lower part 
 * i0          The row index that diagonal starts from
 * y           On output y = x
 * ldy         The leading dimension of y
 * zero        If nonzero, zero the triangular part not copied
 *
 * NOTE: the distance between x and y can be less than ldx, or
 *       x and y *cannot* overlap at all
 *
 ******************************************************************************/

void Num_copy_trimatrix_zprimme(Complex_Z *x, int m, int n, int ldx, int ul,
      int i0, Complex_Z *y, int ldy, int zero) {

   int i, j, jm;
   Complex_Z tzero = {+0.0e+00,+0.0e00};             /*constants*/

   assert(m == 0 || n == 0 || (ldx >= m && ldy >= m));
   if (x == y) return;
   if (ul == 0) {
      /* Copy upper part */

      if (ldx == ldy && (y > x ? y-x : x-y) < ldx) {
         /* x and y overlap */
         for (i=0; i<n; i++) {
            memmove(&y[i*ldy], &x[i*ldx], sizeof(Complex_Z)*min(i0+i+1, m));
            /* zero lower part*/
            if (zero) for (j=min(i0+i+1, m); j<m; j++) y[i*ldy+j] = tzero;
         }
      }
      else {
         /* x and y don't overlap */
         for (i=0; i<n; i++) {
            for (j=0, jm=min(i0+i+1, m); j<jm; j++)
               y[i*ldy+j] = x[i*ldx+j];
            /* zero lower part*/
            if (zero) for (j=min(i0+i+1, m); j<m; j++) y[i*ldy+j] = tzero;
         }
      }
   }
   else {
      /* Copy lower part */

      if (ldx == ldy && (y > x ? y-x : x-y) < ldx) {
         /* x and y overlap */
         for (i=0; i<n; i++) {
            memmove(&y[i*ldy+i0+i], &x[i*ldx+i0+i], sizeof(Complex_Z)*(m-min(i0+i, m)));
            /* zero upper part*/
            if (zero) for (j=0, jm=min(i0+i, m); j<jm; j++) y[i*ldy+j] = tzero;
         }
      }
      else {
         /* x and y don't overlap */
         for (i=0; i<n; i++) {
            for (j=i+i0; j<m; j++)
               y[i*ldy+j] = x[i*ldx+j];
            /* zero upper part*/
            if (zero) for (j=0, jm=min(i0+i, m); j<jm; j++) y[i*ldy+j] = tzero;
         }
      }
   }
}


/******************************************************************************
 * Function Num_copy_trimatrix - Copy the upper triangular part of the matrix x
 *    into y contiguously, i.e., y has all columns of x row-stacked
 *
 * PARAMETERS
 * ---------------------------
 * x           The source upper triangular matrix
 * m           The number of rows of x
 * n           The number of columns of x
 * ldx         The leading dimension of x
 * i0          The row index that diagonal starts from
 * y           On output y = x and nonzero elements of y are contiguous
 * ly          Output the final length of y
 *
 * NOTE: x and y *cannot* overlap
 *
 ******************************************************************************/

void Num_copy_trimatrix_compact_zprimme(Complex_Z *x, int m, int n, int ldx, int i0, Complex_Z *y, int *ly) {
   int i, j, k;

   assert(m == 0 || n == 0 || ldx >= m);

   for (i=0, k=0; i<n; i++)
      for (j=0; j<=i+i0; j++)
         y[k++] = x[i*ldx+j];
   if (ly) *ly = k;
}

/******************************************************************************
 * Function Num_copy_trimatrix - Copy y into the upper triangular part of the
 *    matrix x
 *
 * PARAMETERS
 * ---------------------------
 * x           The source vector
 * m           The number of rows of y
 * n           The number of columns of y
 * i0          The row index that diagonal starts from
 * y           On output the upper triangular part of y has x
 * ldy         The leading dimension of y
 *
 * NOTE: x and y *cannot* overlap
 *
 ******************************************************************************/

void Num_copy_compact_trimatrix_zprimme(Complex_Z *x, int m, int n, int i0, Complex_Z *y, int ldy) {

   int i, j, k;

   assert(m == 0 || n == 0 || (ldy >= m && m >= n));

   for (i=n-1, k=(n+1)*n/2+i0*n-1; i>=0; i--)
      for (j=i+i0; j>=0; j--)
         y[i*ldy+j] = x[k--];
}


/******************************************************************************
 * Function Num_update_VWXR - This subroutine performs the next operations:
 *
 *    X0 = V*h(nX0b+1:nX0e), X1 = V*h(nX1b+1:nX1e), X2 = V*h(nX2b+1:nX2e)
 *    Wo = W*h(nWob+1:nWoe),
 *    R = W*h(nRb+1:nRe) - W*h(nRb+1:nRe)*diag(hVals(nRb+1:nRe)),
 *    Rnorms = norms(R),
 *    rnorms = norms(W*h(nrb+1:nre) - W*h(nrb+1:nre)*diag(hVals(nrb+1:nre)))
 *
 * NOTE: if Rnorms and rnorms are requested, nRb-nRe+nrb-nre < mV
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V, W        input basis
 * mV,nV,ldV   number of rows and columns and leading dimension of V and W
 * h           input rotation matrix
 * nh          Number of columns of h
 * ldh         The leading dimension of h
 * hVals       Array of values
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * X0          Output matrix V*h(nX0b:nX0e-1) (optional)
 * nX0b, nX0e  Range of columns of h
 * X1          Output matrix V*h(nX1b:nX1e-1) (optional)
 * nX1b, nX1e  Range of columns of h
 * X2          Output matrix V*h(nX2b:nX2e-1) (optional)
 * nX2b, nX2e  Range of columns of h
 * Wo          Output matrix W*h(nWob:nWoe-1) (optional)
 * nWob, nWoe  Range of columns of h
 * R           Output matrix (optional)
 * nRb, nRe    Range of columns of h and hVals
 * Rnorms      Output array with the norms of R (optional)
 * rnorms      Output array with the extra residual vector norms (optional)
 * nrb, nre    Columns of residual vector to compute the norm
 * 
 * NOTE: n*e, n*b are zero-base indices of ranges where the first value is
 *       included and the last isn't.
 *
 ******************************************************************************/

int Num_update_VWXR_zprimme(Complex_Z *V, Complex_Z *W, int mV, int nV, int ldV,
   Complex_Z *h, int nh, int ldh, double *hVals,
   Complex_Z *X0, int nX0b, int nX0e, int ldX0,
   Complex_Z *X1, int nX1b, int nX1e, int ldX1,
   Complex_Z *X2, int nX2b, int nX2e, int ldX2,
   Complex_Z *Wo, int nWob, int nWoe, int ldWo,
   Complex_Z *R, int nRb, int nRe, int ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   Complex_Z *rwork, int lrwork, primme_params *primme) {

   int i, j;         /* Loop variables */
   int m=min(PRIMME_BLOCK_SIZE, mV);   /* Number of rows in the cache */
   int nXb, nXe, nYb, nYe, ldX, ldY;
   Complex_Z *X, *Y;
   double *tmp, *tmp0;
   Complex_Z tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};

   /* Return memory requirements */
   if (V == NULL) {
      return 2*m*nV;
   }

   /* R or Rnorms or rnorms imply W */
   assert(!(R || Rnorms || rnorms) || W);

   nXb = min(min(min(min(X0?nX0b:INT_MAX, X1?nX1b:INT_MAX), X2?nX2b:INT_MAX),
         R?nRb:INT_MAX), rnorms?nrb:INT_MAX);
   nXe = max(max(max(X0?nX0e:0, X1?nX1e:0), R?nRe:0), rnorms?nre:0);
   nYb = min(min(Wo?nWob:INT_MAX, R?nRb:INT_MAX), rnorms?nrb:INT_MAX);
   nYe = max(max(Wo?nWoe:0, R?nRe:0), rnorms?nre:0);

   assert((nXe-nXb+nYe-nYb)*nV <= lrwork); /* Check workspace for X and Y */
   assert(2*(nRe-nRb+nre-nrb) <= lrwork); /* Check workspace for tmp and tmp0 */

   X = rwork;
   Y = rwork + m*(nXe-nXb);
   ldX = ldY = m;

   if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = 0.0;
   if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = 0.0;

   for (i=0; i < mV; i+=m, m=min(m,mV-i)) {
      /* X = V*h(nXb:nXe-1) */
      Num_gemm_zprimme("N", "N", m, nXe-nXb, nV, tpone,
         &V[i], ldV, &h[nXb*ldh], ldh, tzero, X, ldX);

      /* X0 = X(nX0b-nXb:nX0e-nXb-1) */
      if (X0) Num_copy_matrix_zprimme(&X[ldX*(nX0b-nXb)], m, nX0e-nX0b,
            ldX, &X0[i], ldX0);

      /* X1 = X(nX1b-nXb:nX1e-nXb-1) */
      if (X1) Num_copy_matrix_zprimme(&X[ldX*(nX1b-nXb)], m, nX1e-nX1b,
            ldX, &X1[i], ldX1);

      /* X2 = X(nX2b-nXb:nX2e-nXb-1) */
      if (X2) Num_copy_matrix_zprimme(&X[ldX*(nX2b-nXb)], m, nX2e-nX2b,
            ldX, &X2[i], ldX2);

      /* Y = W*h(nYb:nYe-1) */
      if (nYb < nYe) Num_gemm_zprimme("N", "N", m, nYe-nYb, nV,
            tpone, &W[i], ldV, &h[nYb*ldh], ldh, tzero, Y, ldY);

      /* Wo = Y(nWob-nYb:nWoe-nYb-1) */
      if (Wo) Num_copy_matrix_zprimme(&Y[ldY*(nWob-nYb)], m, nWoe-nWob,
            ldY, &Wo[i], ldWo);

      /* R = Y(nRb-nYb:nRe-nYb-1) - X(nRb-nYb:nRe-nYb-1)*diag(nRb:nRe-1) */
      for (j=nRb; j<nRe; j++) {
         Num_compute_residual_zprimme(m, hVals[j], &X[ldX*(j-nXb)], &Y[ldY*(j-nYb)],
               &R[i+ldR*(j-nRb)]);
         if (Rnorms) {
            Complex_Z ztmp;
            ztmp = Num_dot_zprimme(m, &R[i+ldR*(j-nRb)], 1, &R[i+ldR*(j-nRb)], 1);
            Rnorms[j-nRb] += *(double*)&ztmp;
         }
      }

      /* rnorms = Y(nrb-nYb:nre-nYb-1) - X(nrb-nYb:nre-nYb-1)*diag(nrb:nre-1) */
      if (rnorms) for (j=nrb; j<nre; j++) {
         Complex_Z ztmp;
         Num_compute_residual_zprimme(m, hVals[j], &X[ldX*(j-nXb)], &Y[ldY*(j-nYb)],
               &Y[ldY*(j-nYb)]);
         ztmp = Num_dot_zprimme(m, &Y[ldY*(j-nYb)], 1, &Y[ldY*(j-nYb)], 1);
         rnorms[j-nrb] += *(double*)&ztmp;
      }
   }

   /* Reduce Rnorms and rnorms and sqrt the results */

   if (primme->globalSumDouble) {
      tmp = (double*)rwork;
      j = 0;
      if (Rnorms) for (i=nRb; i<nRe; i++) tmp[j++] = Rnorms[i-nRb];
      if (rnorms) for (i=nrb; i<nre; i++) tmp[j++] = rnorms[i-nrb];
      tmp0 = tmp+j;
      if (j) primme->globalSumDouble(tmp, tmp0, &j, primme);
      j = 0;
      if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = sqrt(tmp0[j++]);
      if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = sqrt(tmp0[j++]);
   }
   else {
      if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = sqrt(Rnorms[i-nRb]);
      if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = sqrt(rnorms[i-nrb]);
   }

   return 0; 
}

/******************************************************************************
 * Subroutine permute_vecs - This routine permutes a set of vectors according
 *            to a permutation array perm.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * m, n, ld    The number of rows and columns and the leading dimension of vecs
 * perm        The permutation of the columns
 * rwork       Temporary space of size the number of rows
 * iwork       Temporary space of size the number of columns
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * vecs        The matrix whose columns will be reordered
 *
 ******************************************************************************/

void permute_vecs_zprimme(Complex_Z *vecs, int m, int n, int ld, int *perm_,
      Complex_Z *rwork, int *iwork) {

   int currentIndex;     /* Index of vector in sorted order                   */
   int sourceIndex;      /* Position of out-of-order vector in original order */
   int destinationIndex; /* Position of out-of-order vector in sorted order   */
   int tempIndex;        /* Used to swap                                      */
   int *perm=iwork;      /* A copy of perm_                                   */

   assert((perm_>iwork?perm_-iwork:iwork-perm_) >= n);

   /* Copy of perm_ into perm, to avoid to modify the input permutation */

   for (tempIndex=0; tempIndex<n; tempIndex++)
      perm[tempIndex] = perm_[tempIndex];

   /* Continue until all vectors are in the sorted order */

   currentIndex = 0;
   while (1) {

      /* Find a vector that does not belong in its original position */
      while ((currentIndex < n) && (perm[currentIndex] == currentIndex)) {
         currentIndex++;
      }

      /* Return if they are in the sorted order */
      if (currentIndex >= n) {
         return;
      }

      /* Copy the vector to a buffer for swapping */
      Num_zcopy_zprimme(m, &vecs[currentIndex*ld], 1, rwork, 1);

      destinationIndex = currentIndex;
      /* Copy vector perm[destinationIndex] into position destinationIndex */

      while (perm[destinationIndex] != currentIndex) {

         sourceIndex = perm[destinationIndex];
         Num_zcopy_zprimme(m, &vecs[sourceIndex*ld], 1, 
            &vecs[destinationIndex*ld], 1);
         tempIndex = perm[destinationIndex];
         perm[destinationIndex] = destinationIndex;
         destinationIndex = tempIndex;
      }

      /* Copy the vector from the buffer to where it belongs */
      Num_zcopy_zprimme(m, rwork, 1, &vecs[destinationIndex*ld], 1);
      perm[destinationIndex] = destinationIndex;

      currentIndex++;
   }

   /* Check permutation */
   for (currentIndex=0; currentIndex < n; currentIndex++)
      assert(perm[currentIndex] == currentIndex);

}



/******************************************************************************
 * Subroutine Num_compact_vecs - copy certain columns of matrix into another
 *       matrix, i.e., work = vecs(perm). If avoidCopy and perm indices are
 *       consecutive the routine returns a reference in vecs and doesn't copy.
 *            
 *
 * PARAMETERS
 * ---------------------------
 * 
 * vecs        The matrix
 * m           The number of rows of vecs
 * n           The number of columns of to copy
 * ld          The leading dimension of vecs
 * perm        The indices of columns to copy
 * work        The columns are copied to this matrix
 * ldwork      The leading dimension of work
 * avoidCopy   If nonzero, the copy is avoid
 *
 * return      Reference of a matrix with the columns ordered as perm
 *
 ******************************************************************************/

Complex_Z* Num_compact_vecs_zprimme(Complex_Z *vecs, int m, int n, int ld, int *perm,
      Complex_Z *work, int ldwork, int avoidCopy) {

   int i;

   if (avoidCopy) {
      for (i=0; i<n-1 && perm[i]+1 == perm[i+1]; i++);
      if (i >= n-1) return &vecs[ld*perm[0]];
   }

   for (i=0; i < n; i++) {
      Num_copy_matrix_zprimme(&vecs[perm[i]*ld], m, 1, ld, &work[i*ldwork], ld);
   }
   return work;
}
