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
 * File: factorize.c
 *
 * Purpose - Functions to factorize and back-solve a hermitian matrix M.
 *  
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "primme.h"
#include "factorize_d.h"
#include "numerical_d.h"

/******************************************************************************
 * Function UDUDecompose - This function computes an UDU decomposition of the
 *   matrix M.  See LAPACK routine dsytrf for more information on how the
 *   decomposition is performed.
 *
 *
 * Input Parameters
 * ----------------
 * M  A (numOrthoConst+numEvals) x (numOrthoConst+numEvals) array that contains 
 *    the upper triagular portion of a dimM x dimM hermitian matrix.  
 *    The leading dimension of the array is numOrthoConst+numEvals.
 *
 * dimM  The dimension of the matrix M
 *
 * rwork Real work array of dimension at least dimM.  Optimal size is dimM*NB 
 *       where NB is the block size returned by LAPACK routine ilaenv.
 * 
 *
 * Output Parameters
 * -----------------
 * UDU  Array of dimension dimM x dimM containing the UDU decomposition of M.
 *
 * ipivot  Integer array of length dimM containing pivot mapping
 *
 *
 * Return Value
 * ------------
 * int error code: 0 upon success
 *                 dsytrf error code
 ******************************************************************************/
 
int UDUDecompose_dprimme(double *M, double *UDU, int *ipivot, int dimM, 
   double *rwork, int rworkSize, primme_params *primme) {

   int i, j;
   int info;

   /* Quick return for M of dimension 1 */

   if (dimM <= 1) {
      *UDU = *M;
      info = 0;
   }
   else {

      /* Copy the upper triangular portion of M into UDU */

      for (j = 0; j < dimM; j++) {
         for (i = 0; i <= j; i++) {
            UDU[dimM*j+i] = M[(primme->numOrthoConst+primme->numEvals)*j+i];
         }
      }

      /* Perform the decomposition */
      Num_dsytrf_dprimme("U", dimM, UDU, dimM, ipivot, rwork, rworkSize, &info);
   }

   return info;
}

/******************************************************************************
 * Function UDUSolve - This function solves a dense hermitian linear system
 *   given a right hand side (rhs) and a UDU factorization.
 *
 *
 * Input Parameters
 * ----------------
 * UDU     Two-dimensional of dimension dim and leading dimension dim.
 *         Contains block diagonal and multipliers necessary to construct
 *         the upper triangular matrix U.  See LAPACK routine dsytrf for more
 *         details.
 *
 * ipivot  Permutation array that determines how rows and columns of the
 *         factorization were permuted for stability.
 *
 * dim     The dimension of the linear system
 *
 * rhs     The right hand side of the linear system
 *
 * primme  Structure containing various solver parameters
 *
 *
 * Output Parameters
 * -----------------
 * sol     The solution of the linear system 
 *
 ******************************************************************************/

int UDUSolve_dprimme(double *UDU, int *ipivot, int dim, double *rhs, 
   double *sol) {

   int info;

   if (dim == 1) {
      *sol = *rhs/(*UDU); 
      info = 0;
   }
   else {
      Num_dcopy_dprimme(dim, rhs, 1, sol, 1);
      Num_dsytrs_dprimme("U", dim, 1, UDU, dim, ipivot, sol, dim, &info);
   }

   return info;

}
