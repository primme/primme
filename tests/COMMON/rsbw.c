/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2017 College of William & Mary,
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
 * File: rsbw.c
 * 
 * Purpose - librsb wrapper.
 * 
 ******************************************************************************/

#include <rsb.h>        /* for rsb_lib_init */
#include <blas_sparse.h>
#include <stdio.h>
#include <assert.h>
#include "rsbw.h"
#include "num.h"
#include <math.h>

int readMatrixRSB(const char* matrixFileName, blas_sparse_matrix *matrix, double *fnorm) {
#if !defined(RSB_NUMERICAL_TYPE_DOUBLE) || !defined(RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX)
   fprintf(stderr, "Needed librsb with support for 'double' and 'double complex'.\n");
   return -1;
#else
#  ifdef USE_DOUBLECOMPLEX
   rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX;
#  else
   rsb_type_t typecode = RSB_NUMERICAL_TYPE_DOUBLE;
#  endif
   if (rsb_lib_init(RSB_NULL_INIT_OPTIONS) != RSB_ERR_NO_ERROR) assert(0);
   *matrix = blas_invalid_handle;
   if ((rsb_perror(NULL, rsb_lib_init(RSB_NULL_INIT_OPTIONS)))!=RSB_ERR_NO_ERROR) {
     fprintf(stderr, "Error while initializing librsb.\n");
     return -1;
   }

   *matrix = rsb_load_spblas_matrix_file_as_matrix_market(matrixFileName, typecode);
   if ( *matrix == blas_invalid_handle) {
      fprintf(stderr, "ERROR: Could not read matrix file\n");
      return -1;
   }

   if (BLAS_ussp(*matrix, blas_rsb_autotune_next_operation) != 0) assert(0);
   if (BLAS_dusget_infinity_norm(*matrix, fnorm, blas_no_trans) != 0) assert(0);

   return 0;
#endif
}

void RSBMatvec(void *x, void *y, int *blockSize, primme_params *primme) {
   int i;
   SCALAR *xvec, *yvec;
   blas_sparse_matrix *matrix;
#ifdef USE_DOUBLECOMPLEX
   const SCALAR one=(SCALAR)1;
#endif
  
   if (*blockSize <= 0) return; 

   matrix = (blas_sparse_matrix *)primme->matrix;
   xvec = (SCALAR *)x;
   yvec = (SCALAR *)y;

   for (i=0; i<*blockSize*primme->nLocal; i++)
      yvec[i] = 0;
#ifndef USE_DOUBLECOMPLEX
   if(BLAS_dusmm(blas_colmajor, blas_no_trans, *blockSize, 1.0, *matrix, xvec, primme->nLocal, yvec, primme->nLocal) != 0) assert(0);
#else
   if(BLAS_zusmm(blas_colmajor, blas_no_trans, *blockSize, &one, *matrix, xvec, primme->nLocal, yvec, primme->nLocal) != 0) assert(0);
#endif
}

void RSBMatvecSVD(void *x, int *ldx, void *y, int *ldy, int *blockSize,
                  int *trans, primme_svds_params *primme_svds) {
   
   int i, j;
   SCALAR *xvec, *yvec;
   blas_sparse_matrix *matrix;
#ifdef USE_DOUBLECOMPLEX
   const SCALAR one=(SCALAR)1;
#endif
   
   matrix = (blas_sparse_matrix *)primme_svds->matrix;
   xvec = (SCALAR *)x;
   yvec = (SCALAR *)y;

   for (i=0; i<(*blockSize); i++) {
     if (*trans == 0){
      for (j=0; j<primme_svds->mLocal; j++) yvec[i*(*ldy)+j] = 0;
     } else {
      for (j=0; j<primme_svds->nLocal; j++) yvec[i*(*ldy)+j] = 0;
     }
   }

#ifndef USE_DOUBLECOMPLEX
   if (BLAS_dusmm(blas_colmajor, *trans == 0 ? blas_no_trans : blas_trans, *blockSize, 1.0, *matrix, xvec, *ldx, yvec, *ldy) != 0) assert(0);
#else
   if (BLAS_zusmm(blas_colmajor, *trans == 0 ? blas_no_trans : blas_trans, *blockSize, &one, *matrix, xvec, *ldx, yvec, *ldy) != 0) assert(0);
#endif
}


/******************************************************************************
 * Generates the diagonal of A.
 *
 *    P = Diag(A)
 *
 * This will be used with solver provided shifts as (P-shift_i)^(-1) 
******************************************************************************/
static void getDiagonal(blas_sparse_matrix matrix, double *diag) {
#ifndef USE_DOUBLECOMPLEX
   if (BLAS_dusget_diag(matrix, diag) != 0) assert(0);
#else
   int n = BLAS_usgp(matrix, blas_num_rows), i;
   SCALAR *d = (SCALAR *)primme_calloc(n, sizeof(PRIMM_NUM), "aux");
   if (BLAS_zusget_diag(matrix, d) != 0) assert(0);
   for (i=0; i<n; i++) diag[i] = REAL_PART(d[i]);
#endif
}

int createInvDiagPrecRSB(blas_sparse_matrix matrix, double shift, double **prec) {
   int i, n = BLAS_usgp(matrix, blas_num_rows);
   double *diag;

   diag = (double*)primme_calloc(n, sizeof(double), "diag");
   getDiagonal(matrix, diag);
   for (i=0; i<n; i++)
      diag[i] -= shift;
   *prec = diag;
   return 1;
}


/******************************************************************************
 * Generates sum of square values per rows and then per columns 
 *
******************************************************************************/
static void getSumSquares(blas_sparse_matrix matrix, double *diag) {
   int nnz = BLAS_usgp(matrix, blas_num_nonzeros);
   int m = BLAS_usgp(matrix, blas_num_rows);
   int n = BLAS_usgp(matrix, blas_num_cols);
   int i, *AI, *AJ;
   double *sumr = diag, *sumc = &diag[m], v;
   SCALAR *A;

   for (i=0; i < m + n; i++) {
      diag[i] = 0.0;
   }

   A = (SCALAR *)primme_calloc(nnz, sizeof(SCALAR), "A");
   AI = (int *)primme_calloc(nnz*2, sizeof(int), "AI AJ");
   AJ = AI + nnz;
   for (i=0; i<nnz; i++) {
      v = A[i]*CONJ(A[i]);
      sumr[AI[i]] += v;
      sumc[AJ[i]] += v;
   }
}

int createInvNormalPrecRSB(blas_sparse_matrix matrix, double shift, double **prec) {
   int i;
   double *diag, minDenominator=1e-14;
   int n = BLAS_usgp(matrix, blas_num_rows)+BLAS_usgp(matrix, blas_num_cols);

   diag = (double*)primme_calloc(n, sizeof(double), "diag");
   getSumSquares(matrix, diag);
   for (i=0; i<n; i++) {
      diag[i] -= shift*shift;
      if (fabs(diag[i]) < minDenominator)
         diag[i] = copysign(minDenominator, diag[i]);
   }
   *prec = diag;
   return 1;
}
