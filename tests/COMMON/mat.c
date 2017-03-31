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
void FORTRAN_FUNCTION(atmuxr)(int*, int*, double*, double*, double*, int*, int*);
void FORTRAN_FUNCTION(ilut)(int*, double*, int*, int*, int*, double*, double*, int*, int*, int*,
                            double*, double*, int*, int*, int*, int*);
void FORTRAN_FUNCTION(lusol0)(int*, double*, double*, double*, int*, int*);
#else
void FORTRAN_FUNCTION(zamux)(int*, SCALAR*, SCALAR*, SCALAR*, int*, int*);
void FORTRAN_FUNCTION(zatmuxr)(int*, int*, SCALAR*, SCALAR*, SCALAR*, int*, int*);
void FORTRAN_FUNCTION(zilut)(int*, SCALAR*, int*, int*, int*, double*, SCALAR*, int*, int*, int*,
                             SCALAR*, int*, int*);
void FORTRAN_FUNCTION(zlusol)(int*, SCALAR*, SCALAR*, SCALAR*, int*, int*);
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
void CSRMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {
   
   int i;
   int n = (int)primme->n;
   SCALAR *xvec, *yvec;
   CSRMatrix *matrix;
   
   matrix = (CSRMatrix *)primme->matrix;
   xvec = (SCALAR *)x;
   yvec = (SCALAR *)y;

   for (i=0;i<*blockSize;i++) {
#ifndef USE_DOUBLECOMPLEX
      FORTRAN_FUNCTION(amux)
#else
      FORTRAN_FUNCTION(zamux)
#endif
            (&n, &xvec[*ldx*i], &yvec[*ldy*i], 
             matrix->AElts, matrix->JA, matrix->IA);
   }
   *ierr = 0;
}

void CSRMatrixMatvecSVD(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, int *trans, primme_svds_params *primme_svds, int *ierr) {
   
   int i;
   int m = (int)primme_svds->m;
   int n = (int)primme_svds->n;
   SCALAR *xvec, *yvec;
   CSRMatrix *matrix;
   
   matrix = (CSRMatrix *)primme_svds->matrix;
   xvec = (SCALAR *)x;
   yvec = (SCALAR *)y;

   if (*trans == 0) {
      for (i=0;i<*blockSize;i++) {
#ifndef USE_DOUBLECOMPLEX
         FORTRAN_FUNCTION(amux)
#else
         FORTRAN_FUNCTION(zamux)
#endif
              (&m, &xvec[(*ldx)*i], &yvec[(*ldy)*i], 
               matrix->AElts, matrix->JA, matrix->IA);
      }
   } else {
      for (i=0;i<*blockSize;i++) {
#ifndef USE_DOUBLECOMPLEX
         FORTRAN_FUNCTION(atmuxr)
#else
         FORTRAN_FUNCTION(zatmuxr)
#endif
           (&n, &m, &xvec[(*ldx)*i], 
            &yvec[(*ldy)*i], matrix->AElts, matrix->JA, matrix->IA);
      }
   }
   *ierr = 0;
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
   double *diag;

   diag = (double*)primme_calloc(matrix->n, sizeof(double), "diag");
   getDiagonal(matrix, diag);
   for (i=0; i<matrix->n; i++)
      diag[i] -= shift;
   *prec = diag;
   return 1;
}

static void ApplyInvDiagPrecNativeGen(SCALAR *xvec, int ldx, SCALAR *yvec,
      int ldy, int nLocal, int bs, double *diag, double *shifts, double aNorm) {
   int i, j;
   const double minDenominator = 1e-14*(aNorm >= 0.0L ? aNorm : 1.);

   for (i=0; i<bs; i++) {
      double shift = shifts ? shifts[i] : 0.0;
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (j=0; j<nLocal; j++) {
         double d = diag[j] - shift;
         d = (fabs(d) > minDenominator) ? d : copysign(minDenominator, d);
         yvec[ldy*i+j] = xvec[ldx*i+j]/d;
      }
   }
}


void ApplyInvDiagPrecNative(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {
   ApplyInvDiagPrecNativeGen((SCALAR*)x, *ldx, (SCALAR*)y, *ldy,
      primme->nLocal, *blockSize, (double*)primme->preconditioner, NULL, primme->aNorm);
   *ierr = 0;
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

void ApplyInvDavidsonDiagPrecNative(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, 
                                        primme_params *primme, int *ierr) {
   ApplyInvDiagPrecNativeGen((SCALAR*)x, *ldx, (SCALAR*)y, *ldy,
      primme->nLocal, *blockSize, (double*)primme->preconditioner,
      primme->ShiftsForPreconditioner, primme->aNorm);
   *ierr = 0;
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
   SCALAR *W;
   int *iW;
   CSRMatrix *factors;

   if (shift != 0.0) {
      shiftCSRMatrix(-shift, (CSRMatrix*)matrix);
   }

   /* Work arrays */
   W = (SCALAR *)primme_calloc(matrix->n+1, sizeof(SCALAR), "W");
   iW = (int *)primme_calloc(matrix->n*2, sizeof(int), "iW");

   /* Max size of factorization */
   lenFactors = 9*matrix->nnz;

   factors = (CSRMatrix *)primme_calloc(1,  sizeof(CSRMatrix), "factors");
   factors->AElts = (SCALAR *)primme_calloc(lenFactors,
                                sizeof(SCALAR), "iluElts");
   factors->JA = (int *)primme_calloc(lenFactors, sizeof(int), "Jilu");
   factors->IA = (int *)primme_calloc(matrix->n+1, sizeof(int), "Iilu");
   factors->n = matrix->n;
   factors->nnz = lenFactors;
   
   FORTRAN_FUNCTION(zilut)
         ((int*)&matrix->n, (SCALAR*)matrix->AElts, (int*)matrix->JA,
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

void ApplyILUTPrecNative(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {
   int i;
   int n = (int)primme->n;
   SCALAR *xvec, *yvec;
   CSRMatrix *prec;
   
   prec = (CSRMatrix *)primme->preconditioner;
   xvec = (SCALAR *)x;
   yvec = (SCALAR *)y;

   for (i=0; i<*blockSize; i++) {
#ifdef USE_DOUBLECOMPLEX
      FORTRAN_FUNCTION(zlusol)
#else
      FORTRAN_FUNCTION(lusol0)
#endif
             (&n, &xvec[*ldx*i], &yvec[*ldy*i],
              prec->AElts, prec->JA, prec->IA);
   }
   *ierr = 0;
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

/******************************************************************************
 * Generates sum of square values per rows and then per columns 
 *
******************************************************************************/
static void getSumSquares(const CSRMatrix *matrix, double *diag) {
   int i, j;
   double *sumr = diag, *sumc = &diag[matrix->m], v;

   for (i=0; i < matrix->m + matrix->n; i++) {
      diag[i] = 0.0;
   }

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   for (i=0; i < matrix->m; i++) {
      for (j=matrix->IA[i]; j <= matrix->IA[i+1]-1; j++) {
         v = REAL_PART(matrix->AElts[j-1]*CONJ(matrix->AElts[j-1]));
         sumr[i]   += v;
         sumc[matrix->JA[j-1]-1] += v;
      }
   }
}


/******************************************************************************
 * Applies the diagonal preconditioner over A'A or AA'
 *
******************************************************************************/

int createInvNormalPrecNative(const CSRMatrix *matrix, double shift, double **prec) {
   int i;
   double *diag, minDenominator=1e-14;

   diag = (double*)primme_calloc(matrix->m+matrix->n, sizeof(double), "diag");
   getSumSquares(matrix, diag);
   for (i=0; i<matrix->m+matrix->n; i++) {
      diag[i] -= shift*shift;
      if (fabs(diag[i]) < minDenominator)
         diag[i] = copysign(minDenominator, diag[i]);
   }
   *prec = diag;
   return 1;
}

static void ApplyInvNormalPrecNativeSvdsGen(void *x, PRIMME_INT *ldx, void *y,
      PRIMME_INT *ldy, int *blockSize, int *mode,
      primme_svds_params *primme_svds, double *shifts) {

   double *diag = (double *)primme_svds->preconditioner;
   SCALAR *xvec = (SCALAR *)x, *yvec = (SCALAR *)y;
   const int bs = *blockSize;
   double *sumr = diag, *sumc = &diag[primme_svds->mLocal];
   
   if (*mode == primme_svds_op_AtA) {
      ApplyInvDiagPrecNativeGen(xvec, *ldx, yvec, *ldy, primme_svds->nLocal, bs,
         sumc, shifts, primme_svds->aNorm);
   }
   else if (*mode == primme_svds_op_AAt) {
      ApplyInvDiagPrecNativeGen(xvec, *ldx, yvec, *ldy, primme_svds->mLocal, bs,
         sumr, shifts, primme_svds->aNorm);
   }
   else if (*mode == primme_svds_op_augmented) {
      ApplyInvDiagPrecNativeGen(xvec, *ldx, yvec, *ldy, primme_svds->nLocal, bs,
         sumc, shifts, primme_svds->aNorm);
      ApplyInvDiagPrecNativeGen(xvec+primme_svds->nLocal, *ldx,
         yvec+primme_svds->nLocal, *ldy, primme_svds->mLocal, bs,
         sumr, shifts, primme_svds->aNorm);
   } 
}

void ApplyInvNormalPrecNative(void *x, PRIMME_INT *ldx, void *y,
      PRIMME_INT *ldy, int *blockSize, int *mode,
      primme_svds_params *primme_svds, int *ierr) {

   ApplyInvNormalPrecNativeSvdsGen(x, ldx, y, ldy, blockSize, mode, primme_svds, NULL);
   *ierr = 0;
}

void ApplyInvDavidsonNormalPrecNative(void *x, PRIMME_INT *ldx, void *y,
      PRIMME_INT *ldy, int *blockSize, int *mode,
      primme_svds_params *primme_svds, int *ierr) {

   primme_params *primme =
      primme_svds->method == *mode ? &primme_svds->primme : &primme_svds->primmeStage2;
   ApplyInvNormalPrecNativeSvdsGen(x, ldx, y, ldy, blockSize, mode, primme_svds,
      primme->ShiftsForPreconditioner);
   *ierr = 0;
}
