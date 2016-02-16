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
 * File: update_projection.c
 *
 * Purpose - Adds blockSize new columns and rows to H.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <assert.h>
#include "primme.h"
#include "const.h"
#include "update_projection_@(pre).h"
#include "numerical_@(pre).h"

/*******************************************************************************
 * Subroutine update_projection - Z = X'*Y. It assumes Z is a hermitian matrix 
 *    whose columns will be updated with blockSize vectors.  Even though space 
 *    for the entire Z is allocated, only the upper triangular portion is 
 *    stored. 
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * X           Matrix with size nLocal x numCols+blockSize
 * ldX         The leading dimension of X
 * Y           Matrix with size nLocal x numCols+blockSize
 * ldY         The leading dimension of Y
 * Z           Matrix with size nLocal x numCols+blockSize
 * numCols     The number of columns that hasn't changed
 * blockSize   The number of columns that has changed
 * rwork       Workspace
 * lrwork      Size of rwork
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * Z           Matrix with size numCols+blockSize with value X'*Y, only
 *             Z(1:numCols+blockSize,numCols:numCols+blockSize) will be updated
 * ldZ         The leading dimension of Z
 *
 ******************************************************************************/

int update_projection_@(pre)primme(@(type) *X, int ldX, @(type) *Y, int ldY,
   @(type) *Z, int ldZ, int nLocal, int numCols, int blockSize, @(type) *rwork,
   int lrwork, primme_params *primme) {

   int count, count_doubles, m;
   @(type) tpone = @(tpone), tzero = @(tzero);

   /* -------------------------- */
   /* Return memory requirements */
   /* -------------------------- */

   if (X == NULL) {
      return (numCols+blockSize)*numCols*2;
   }

   assert(ldX >= nLocal && ldY >= nLocal && ldZ >= numCols+blockSize);

   /* ------------ */
   /* Quick return */
   /* ------------ */

   if (blockSize <= 0) return 0;

   /* --------------------------------------------------------------------- */
   /* Grow Z by blockSize number of rows and columns all at once            */
   /* --------------------------------------------------------------------- */

   m = numCols+blockSize;
   Num_gemm_@(pre)primme("C", "N", m, blockSize, nLocal, tpone, 
      X, ldX, &Y[ldY*numCols], ldY, tzero, &Z[ldZ*numCols], ldZ);

   /* -------------------------------------------------------------- */
   /* Alternative to the previous call:                              */
   /*    Compute next the additional rows of each new column vector. */
   /*    Only the upper triangular portion is computed and stored.   */
   /* -------------------------------------------------------------- */

   /*
   for (j = numCols; j < numCols+blockSize; j++) {
      Num_gemv_@(pre)primme("C", primme->nLocal, j-numCols+1, tpone,
         &X[primme->nLocal*numCols], primme->nLocal, &Y[primme->nLocal*j], 1, 
         tzero, &rwork[maxCols*(j-numCols)+numCols], 1);  
   }
   */

   if (primme->n > 1) {
      /* --------------------------------------------------------------------- */
      /* Reduce the upper triangular part of the new columns in Z.             */
      /* --------------------------------------------------------------------- */

      Num_copy_trimatrix_compact_@(pre)primme(&Z[ldZ*numCols], m, blockSize, ldZ,
            numCols, rwork, &count);
      assert(count*2 <= lrwork);

      count_doubles = count;
#ifdefarithm L_DEFCPLX
      count_doubles *= 2;
#endifarithm
      primme->globalSumDouble(rwork, (double*)&rwork[count], &count_doubles, primme);

      Num_copy_compact_trimatrix_@(pre)primme(&rwork[count], m, blockSize, numCols,
            &Z[ldZ*numCols], ldZ);
   }

   return 0;
}
