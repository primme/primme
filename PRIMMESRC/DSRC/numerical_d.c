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
#include "numerical_private_d.h"
#include "numerical_d.h"

/******************************************************************************/
void Num_dcopy_dprimme(int n, double *x, int incx, double *y, int incy) {
   long long int ln = n;
   long long int lincx = incx;
   long long int lincy = incy;

   DCOPY(&ln, x, &lincx, y, &lincy);
}
/******************************************************************************/

void Num_gemm_dprimme(char *transa, char *transb, int m, int n, int k, 
   double alpha, double *a, int lda, double *b, int ldb, 
   double beta, double *c, int ldc) {
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
   DGEMM(transa_fcd, transb_fcd, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, 
	 c, &lldc);
#else
   DGEMM(transa, transb, &lm, &ln, &lk, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endif

}

/******************************************************************************/
void Num_symm_dprimme(char *side, char *uplo, int m, int n, double alpha, 
   double *a, int lda, double *b, int ldb, double beta, 
   double *c, int ldc) {
   long long int lm = m;
   long long int ln = n;
   long long int llda = lda;
   long long int lldb = ldb;
   long long int lldc = ldc;

#ifdef NUM_CRAY
   _fcd side_fcd, uplo_fcd;

   side_fcd = _cptofcd(side, strlen(side));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   DSYMM(side_fcd, uplo_fcd, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#else
   DSYMM(side, uplo, &lm, &ln, &alpha, a, &llda, b, &lldb, &beta, c, &lldc);
#endif   

}

/******************************************************************************/
void Num_axpy_dprimme(int n, double alpha, double *x, int incx, 
   double *y, int incy) {
   long long int ln = n;
   long long int lincx = incx;
   long long int lincy = incy;

   DAXPY(&ln, &alpha, x, &lincx, y, &lincy);

}

/******************************************************************************/
void Num_gemv_dprimme(char *transa, int m, int n, double alpha, double *a,
   int lda, double *x, int incx, double beta, double *y, int incy) {
   long long int lm = m;
   long long int ln = n;
   long long int llda = lda;
   long long int lincx = incx;
   long long int lincy = incy;

#ifdef NUM_CRAY
   _fcd transa_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
   DGEMV(transa_fcd, &lm, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);
#else
   DGEMV(transa, &lm, &ln, &alpha, a, &llda, x, &lincx, &beta, y, &lincy);

#endif

}

/******************************************************************************/
double Num_dot_dprimme(int n, double *x, int incx, double *y, int incy) {
   long long int ln = n;
   long long int lincx = incx;
   long long int lincy = incy;

   return(DDOT(&ln, x, &lincx, y, &lincy));

}

/******************************************************************************/
void Num_larnv_dprimme(int idist, int *iseed, int length, double *x) {
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

   DLARNV(&lidist, liseed, &llength, x);
   
   for(i=0;i<4;i++)
   {
      iseed[i] = liseed[i];
   }

}

/******************************************************************************/
void Num_scal_dprimme(int n, double alpha, double *x, int incx) {
   long long int ln = n;
   long long int lincx = incx;

   DSCAL(&ln, &alpha, x, &lincx);

}

/******************************************************************************/
void Num_swap_dprimme(int n, double *x, int incx, double *y, int incy) {
   long long int ln = n;
   long long int lincx = incx;
   long long int lincy = incy;

   DSWAP(&ln, x, &lincx, y, &lincy);

}

/******************************************************************************/
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef NUM_ESSL
int Num_dspev_dprimme(int iopt, double *ap, double *w, double *z, int ldz, 
   int n, double *aux, int naux) {
   long long int liopt = iopt;
   long long int lldz = ldz;
   long long int ln = n;
   long long int lnaux = naux;

   int ret;

   ret = dspev(liopt, ap, w, z, lldz, ln, aux, lnaux);
   return (ret);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Num_dsyev_dprimme(char *jobz, char *uplo, int n, double *a, int lda, 
   double *w, double *work, int ldwork, int *info) {
   long long int ln = n;
   long long int llda = lda;
   long long int lldwork = ldwork;
   long long int linfo = (long long int) *info;
   

#ifdef NUM_CRAY
   _fcd jobz_fcd, uplo_fcd;

   jobz_fcd = _cptofcd(jobz, strlen(jobz));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));

   DSYEV(jobz_fcd, uplo_fcd, &ln, a, &llda, w, work, &lldwork, &linfo); 

#else

   DSYEV(jobz, uplo, &ln, a, &llda, w, work, &lldwork, &linfo); 

#endif
    *info = (int)linfo;

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

/******************************************************************************/
void Num_dsytrf_dprimme(char *uplo, int n, double *a, int lda, int *ipivot, 
   double *work, int ldwork, int *info) {
   long long int ln = n;
   long long int llda = lda;

   long long int *lipivot = (long long int *)primme_calloc(n, sizeof(long long int), "lipivot array");
   int i;
   for(i=0;i<n;i++)
   {
      lipivot[i] = ipivot[i];
   }
   
   long long int lldwork = ldwork;
   long long int linfo = (long long int) *info;   

#ifdef NUM_CRAY
   _fcd uplo_fcd;

   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   DSYTRF(uplo_fcd, &ln, a, &llda, lipivot, work, &lldwork, &linfo);
#else

   DSYTRF(uplo, &ln, a, &llda, lipivot, work, &lldwork, &linfo);

#endif

   for(i=0;i<n;i++)
   {
      ipivot[i] = lipivot[i];
   }
   *info = (int) linfo;
}


/******************************************************************************/
void Num_dsytrs_dprimme(char *uplo, int n, int nrhs, double *a, int lda, 
   int *ipivot, double *b, int ldb, int *info) {
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
   long long int linfo = (long long int) *info;

#ifdef NUM_CRAY
   _fcd uplo_fcd;

   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   DSYTRS(uplo_fcd, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
#else

   DSYTRS(uplo, &ln, &lnrhs, a, &llda, lipivot, b, &lldb, &linfo);
#endif

    for(i=0;i<n;i++)
   {
      ipivot[i] = lipivot[i];
   }
    *info = (int) linfo;
}

