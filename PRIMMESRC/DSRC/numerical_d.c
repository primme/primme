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
#include "numerical_private_d.h"
#include "numerical_d.h"

/******************************************************************************/
void Num_dcopy_dprimme(int n, double *x, int incx, double *y, int incy) {

   DCOPY(&n, x, &incx, y, &incy);
}
/******************************************************************************/

void Num_gemm_dprimme(char *transa, char *transb, int m, int n, int k, 
   double alpha, double *a, int lda, double *b, int ldb, 
   double beta, double *c, int ldc) {

#ifdef NUM_CRAY
   _fcd transa_fcd, transb_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
   transb_fcd = _cptofcd(transb, strlen(transb));
   DGEMM(transa_fcd, transb_fcd, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, 
	 c, &ldc);
#else
   DGEMM(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif

}

/******************************************************************************/
void Num_symm_dprimme(char *side, char *uplo, int m, int n, double alpha, 
   double *a, int lda, double *b, int ldb, double beta, 
   double *c, int ldc) {

#ifdef NUM_CRAY
   _fcd side_fcd, uplo_fcd;

   side_fcd = _cptofcd(side, strlen(side));
   uplo_fcd = _cptofcd(uplo, strlen(uplo));
   DSYMM(side_fcd, uplo_fcd, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#else
   DSYMM(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif   

}

/******************************************************************************/
void Num_axpy_dprimme(int n, double alpha, double *x, int incx, 
   double *y, int incy) {

   DAXPY(&n, &alpha, x, &incx, y, &incy);

}

/******************************************************************************/
void Num_gemv_dprimme(char *transa, int m, int n, double alpha, double *a,
   int lda, double *x, int incx, double beta, double *y, int incy) {

#ifdef NUM_CRAY
   _fcd transa_fcd;

   transa_fcd = _cptofcd(transa, strlen(transa));
   DGEMV(transa_fcd, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#else
   DGEMV(transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);

#endif

}

/******************************************************************************/
double Num_dot_dprimme(int n, double *x, int incx, double *y, int incy) {

   return(DDOT(&n, x, &incx, y, &incy));

}

/******************************************************************************/
void Num_larnv_dprimme(int idist, int *iseed, int length, double *x) {
   DLARNV(&idist, iseed, &length, x);

}

/******************************************************************************/
void Num_scal_dprimme(int n, double alpha, double *x, int incx) {

   DSCAL(&n, &alpha, x, &incx);

}

/******************************************************************************/
void Num_swap_dprimme(int n, double *x, int incx, double *y, int incy) {

   DSWAP(&n, x, &incx, y, &incy);

}

/******************************************************************************/
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#ifdef NUM_ESSL
int Num_dspev_dprimme(int iopt, double *ap, double *w, double *z, int ldz, 
   int n, double *aux, int naux) {

   int ret;

   ret = dspev(iopt, ap, w, z, ldz, n, aux, naux);
   return (ret);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#else
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif

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

