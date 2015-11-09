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
 * File: mat.c
 * 
 * Purpose - Wrapper implemented with routines that handle matrices in
 *           CSR format.
 * 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "native.h"

static void getDiagonal(const CSRMatrix *matrix, double *diag);

#ifdef __cplusplus
extern "C" {
#endif

#ifndef USE_DOUBLECOMPLEX
void FORTRAN_FUNCTION(amux)(int*, double*, double*, double*, int*, int*);
void FORTRAN_FUNCTION(ilut)(int*, double*, int*, int*, int*, double*, double*, int*, int*, int*,
                            double*, double*, int*, int*, int*, int*);
void FORTRAN_FUNCTION(lusol0)(int*, double*, double*, double*, int*, int*);
#else
void FORTRAN_FUNCTION(zamux)(int*, PRIMME_NUM*, PRIMME_NUM*, PRIMME_NUM*, int*, int*);
void FORTRAN_FUNCTION(zilut)(int*, PRIMME_NUM*, int*, int*, int*, double*, PRIMME_NUM*, int*, int*, int*,
                             PRIMME_NUM*, int*, int*);
void FORTRAN_FUNCTION(zlusol)(int*, PRIMME_NUM*, PRIMME_NUM*, PRIMME_NUM*, int*, int*);
#endif

#ifdef __cplusplus
}
#endif

/******************************************************************************
 * Applies the matrix vector multiplication on a block of vectors.
 * Because a block function is not available, we call blockSize times
 * the SPARSKIT function amux(). Note the (void *) parameters x, y that must 
 * be cast as doubles for use in amux()
 *
******************************************************************************/
void CSRMatrixMatvec(void *x, void *y, int *blockSize, primme_params *primme) {
   
   int i;
   PRIMME_NUM *xvec, *yvec;
   CSRMatrix *matrix;
   
   matrix = (CSRMatrix *)primme->matrix;
   xvec = (PRIMME_NUM *)x;
   yvec = (PRIMME_NUM *)y;

   for (i=0;i<*blockSize;i++) {
#ifndef USE_DOUBLECOMPLEX
      FORTRAN_FUNCTION(amux)
#else
      FORTRAN_FUNCTION(zamux)
#endif
            (&primme->n, &xvec[primme->nLocal*i], &yvec[primme->nLocal*i], 
             matrix->AElts, matrix->JA, matrix->IA);
   }
}

/******************************************************************************
 * Applies the (already inverted) diagonal preconditioner
 *
 *    y(i) = P*x(i), i=1:blockSize, 
 *    with P = (Diag(A)-shift)^(-1)
 *
******************************************************************************/

int createInvDiagPrecNative(const CSRMatrix *matrix, double shift, double **prec) {
   int i;
   double *diag, minDenominator=1e-14;

   diag = (double*)primme_calloc(matrix->n, sizeof(double), "diag");
   getDiagonal(matrix, diag);
   for (i=0; i<matrix->n; i++)
      diag[i] -= shift;
   for (i=0; i<matrix->n; i++)
      if (fabs(diag[i]) < minDenominator)
         diag[i] = copysign(minDenominator, diag[i]);
   *prec = diag;
   return 1;
}

void ApplyInvDiagPrecNative(void *x, void *y, int *blockSize, 
                                        primme_params *primme) {
   int i, j;
   double *diag;
   PRIMME_NUM *xvec, *yvec;
   const int nLocal = primme->nLocal, bs = *blockSize;
   
   diag = (double *)primme->preconditioner;
   xvec = (PRIMME_NUM *)x;
   yvec = (PRIMME_NUM *)y;

   for (i=0; i<bs; i++)
      for (j=0; j<nLocal; j++)
         yvec[primme->n*i+j] = xvec[primme->n*i+j]/diag[j];
}

/******************************************************************************
 * Applies a Davidson type preconditioner
 *
 *    x(i) = (Diag(A) - primme.Shifts(i) I)^(-1) * y(i),   i=1:blockSize
 *    
 * NOTE that each block vector may have its own shift provided by dprimme
 * in the array primme->ShiftsForPreconditioner
 *
 * To avoid division with too small numbers we limit how small relatively 
 * to ||A|| the denominators can be. In the absense of ||A|| we use 1e-14.
 *
******************************************************************************/

int createInvDavidsonDiagPrecNative(const CSRMatrix *matrix, double **prec) {
   double *diag;

   diag = (double*)primme_calloc(matrix->n, sizeof(double), "diag");
   getDiagonal(matrix, diag);
   *prec = diag;
   return 1;
}

void ApplyInvDavidsonDiagPrecNative(void *x, void *y, int *blockSize, 
                                        primme_params *primme) {
   int i, j;
   double *diag, shift, d, minDenominator;
   PRIMME_NUM *xvec, *yvec;
   const int nLocal = primme->nLocal, bs = *blockSize;
   
   diag = (double *)primme->preconditioner;
   xvec = (PRIMME_NUM *)x;
   yvec = (PRIMME_NUM *)y;
   minDenominator = 1e-14*(primme->aNorm >= 0.0L ? primme->aNorm : 1.);

   for (i=0; i<bs; i++) {
      shift = primme->ShiftsForPreconditioner[i];
      for (j=0; j<nLocal; j++) {
         d = diag[j] - shift;
         d = (fabs(d) > minDenominator) ? d : copysign(minDenominator, d);
         yvec[primme->n*i+j] = xvec[primme->n*i+j]/d;
      }
   }
}

/******************************************************************************
 * Applies the ILUT preconditioner 
 *
 *    y(i) = U^(-1)*( L^(-1)*x(i)), i=1:blockSize, 
 *    with L,U = ilut(A-shift) 
 * 
 * It calls the SPARSKIT lusol0 function for each block vector.
 *
******************************************************************************/

int createILUTPrecNative(const CSRMatrix *matrix, double shift, int level,
                         double threshold, double filter, CSRMatrix **prec) {
#ifdef USE_DOUBLECOMPLEX
   int ierr;
   int lenFactors;
   PRIMME_NUM *W;
   int *iW;
   CSRMatrix *factors;

   if (shift != 0.0) {
      shiftCSRMatrix(-shift, (CSRMatrix*)matrix);
   }

   /* Work arrays */
   W = (PRIMME_NUM *)primme_calloc(matrix->n+1, sizeof(PRIMME_NUM), "W");
   iW = (int *)primme_calloc(matrix->n*2, sizeof(int), "iW");

   /* Max size of factorization */
   lenFactors = 9*matrix->nnz;

   factors = (CSRMatrix *)primme_calloc(1,  sizeof(CSRMatrix), "factors");
   factors->AElts = (PRIMME_NUM *)primme_calloc(lenFactors,
                                sizeof(PRIMME_NUM), "iluElts");
   factors->JA = (int *)primme_calloc(lenFactors, sizeof(int), "Jilu");
   factors->IA = (int *)primme_calloc(matrix->n+1, sizeof(int), "Iilu");
   factors->n = matrix->n;
   factors->nnz = lenFactors;
   
   FORTRAN_FUNCTION(zilut)
         ((int*)&matrix->n, (PRIMME_NUM*)matrix->AElts, (int*)matrix->JA,
          (int*)matrix->IA, &level, &threshold,
          factors->AElts, factors->JA, factors->IA, &lenFactors, W, iW, &ierr);
   
   if (ierr != 0)  {
      fprintf(stderr, "ZILUT factorization could not be completed\n");
      return(-1);
   }

   if (shift != 0.0L) {
      shiftCSRMatrix(shift, (CSRMatrix*)matrix);
   }

   /* free workspace */
   free(W); free(iW);

   *prec = factors;
   return 0;
#else
   int ierr;
   int lenFactors;
   double *W1, *W2;
   int *iW1, *iW2, *iW3;
   CSRMatrix *factors;

   if (shift != 0.0) {
      shiftCSRMatrix(-shift, (CSRMatrix*)matrix);
   }

   /* Work arrays */
   W1 = (double *)primme_calloc( matrix->n+1,  sizeof(double), "W1");
   W2 = (double *)primme_calloc( matrix->n,  sizeof(double), "W2");
   iW1 = (int *)primme_calloc( matrix->n,  sizeof(int), "iW1");
   iW2 = (int *)primme_calloc( matrix->n,  sizeof(int), "iW2");
   iW3 = (int *)primme_calloc( matrix->n,  sizeof(int), "iW2");
   /* Max size of factorization */
   lenFactors = 9*matrix->nnz;
   factors = (CSRMatrix *)primme_calloc(1,  sizeof(CSRMatrix), "factors");
   factors->AElts = (double *)primme_calloc(lenFactors,
                                sizeof(double), "iluElts");
   factors->JA = (int *)primme_calloc(lenFactors, sizeof(int), "Jilu");
   factors->IA = (int *)primme_calloc(matrix->n+1, sizeof(int), "Iilu");
   factors->n = matrix->n;
   factors->nnz = lenFactors;
   
   FORTRAN_FUNCTION(ilut)
        ((int*)&matrix->n, (double*)matrix->AElts, (int*)matrix->JA,
         (int*)matrix->IA, &level, &threshold,
         factors->AElts, factors->JA, factors->IA, &lenFactors, 
         W1, W2, iW1, iW2, iW3, &ierr);
   
   if (ierr != 0)  {
      fprintf(stderr, "ILUT factorization could not be completed\n");
      return(-1);
   }

   if (shift != 0.0L) {
      shiftCSRMatrix(shift, (CSRMatrix*)matrix);
   }

   /* free workspace */
   free(W1); free(W2); free(iW1); free(iW2); free(iW3);

   *prec = factors;
   return 0;
#endif
}

void ApplyILUTPrecNative(void *x, void *y, int *blockSize, primme_params *primme) {
   int i;
   PRIMME_NUM *xvec, *yvec;
   CSRMatrix *prec;
   
   prec = (CSRMatrix *)primme->preconditioner;
   xvec = (PRIMME_NUM *)x;
   yvec = (PRIMME_NUM *)y;

   for (i=0; i<*blockSize; i++) {
#ifdef USE_DOUBLECOMPLEX
      FORTRAN_FUNCTION(zlusol)
#else
      FORTRAN_FUNCTION(lusol0)
#endif
             (&primme->n, &xvec[primme->n*i], &yvec[primme->n*i],
              prec->AElts, prec->JA, prec->IA);
   }
}


/******************************************************************************
 * Generates the diagonal of A.
 *
 *    P = Diag(A)
 *
 * This will be used with solver provided shifts as (P-shift_i)^(-1) 
******************************************************************************/
static void getDiagonal(const CSRMatrix *matrix, double *diag) {
   int i, j;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   for (i=0; i < matrix->n; i++) {
      diag[i] = 0.;
      for (j=matrix->IA[i]; j <= matrix->IA[i+1]-1; j++) {
         if (matrix->JA[j-1]-1 == i) {
            diag[i] = REAL_PART(matrix->AElts[j-1]);
         }
      }
   }
}
