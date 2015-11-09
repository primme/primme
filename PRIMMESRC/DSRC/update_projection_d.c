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
#include "primme.h"
#include "const.h"
#include "update_projection_d.h"
#include "numerical_d.h"

/*******************************************************************************
 * Subroutine update_projection - Z = X'*Y. It assumes Z is a hermitian matrix 
 *    whose columns will be updated with blockSize vectors.  Even though space 
 *    for the entire Z is allocated, only the upper triangular portion is 
 *    stored. 
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * X             Some nLocal x numCols matrix
 * Y             Some nLocal x numCols matrix
 * numCols       Number of rows and columns in Z
 * maxCols       Maximum (leading) dimension of Z
 * blockSize     Number of rows and columns to be added to Z
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * Z      X'*Y where Y is some nLocal x numCols matrix.
 * rwork  Must be at least maxCols*blockSize in length
 ******************************************************************************/

void update_projection_dprimme(double *X, double *Y, double *Z, 
   int numCols, int maxCols, int blockSize, double *rwork, 
   primme_params *primme) {

   int j;    /* Loop variable  */ 
   int count;
   double tpone = +1.0e+00, tzero = +0.0e+00;

   /* --------------------------------------------------------------------- */
   /* Zero the work array to prevent floating point traps during all-reduce */
   /* --------------------------------------------------------------------- */

   for (j = 0; j < maxCols*blockSize; j++) {
      rwork[j] = tzero;
   }

   /* --------------------------------------------------------------------- */
   /* Grow Z by blockSize number of rows and columns all at once            */
   /* --------------------------------------------------------------------- */

   Num_gemm_dprimme("C", "N", numCols+blockSize, blockSize, primme->nLocal, tpone, 
      X, primme->nLocal, &Y[primme->nLocal*numCols], primme->nLocal, 
      tzero, rwork, maxCols);

   /* -------------------------------------------------------------- */
   /* Alternative to the previous call:                              */
   /*    Compute next the additional rows of each new column vector. */
   /*    Only the upper triangular portion is computed and stored.   */
   /* -------------------------------------------------------------- */

   /*
   for (j = numCols; j < numCols+blockSize; j++) {
      Num_gemv_dprimme("C", primme->nLocal, j-numCols+1, tpone,
         &X[primme->nLocal*numCols], primme->nLocal, &Y[primme->nLocal*j], 1, 
         tzero, &rwork[maxCols*(j-numCols)+numCols], 1);  
   }
   */
   
   count = maxCols*blockSize;
   (*primme->globalSumDouble)(rwork, &Z[maxCols*numCols], &count, primme);
}
