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
#include "factorize_@(pre).h"
#include "numerical_@(pre).h"

/******************************************************************************
 * Function UDUDecompose - This function computes an UDU decomposition of the
 *   matrix M.  See LAPACK routine dsytrf for more information on how the
 *   decomposition is performed.
 *
 *
 * Input Parameters
 * ----------------
 * M  A (numOrthoConst+numEvals) x (numOrthoConst+numEvals) array that contains 
 *    the upper triangular portion of a dimM x dimM hermitian matrix.  
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
 
int UDUDecompose_@(pre)primme(@(type) *M, int ldM, @(type) *UDU, int ldUDU,
   int *ipivot, int dimM, @(type) *rwork, int rworkSize, primme_params *primme) {

   int info;

   /* Quick return for M with dimension 0 */

   if (dimM == 0) return 0;

   /* Return memory requirement */

   if (M == NULL) {
      @(type) w;
#ifdefarithm L_DEFCPLX
      Num_zhetrf_zprimme("U", dimM, UDU, ldUDU, ipivot, &w, -1, &info);
#endifarithm
#ifdefarithm L_DEFREAL
      Num_dsytrf_dprimme("U", dimM, UDU, ldUDU, ipivot, &w, -1, &info);
#endifarithm
      return (int)*(double*)&w;
    }

   /* if ld is zero, change by the matrix size */
   if (ldUDU == 0) ldUDU = dimM;

   /* Quick return for M with dimension 1 */

   if (dimM <= 1) {
      *UDU = *M;
      info = 0;
   }
   else {

      /* Copy the upper triangular portion of M into UDU */

      Num_copy_trimatrix_@(pre)primme(M, dimM, dimM, ldM, 0 /* up */, 0, UDU, ldUDU, 0);

      /* Perform the decomposition */
#ifdefarithm L_DEFCPLX
      Num_zhetrf_zprimme("U", dimM, UDU, ldUDU, ipivot, rwork, rworkSize, &info);
#endifarithm
#ifdefarithm L_DEFREAL
      Num_dsytrf_dprimme("U", dimM, UDU, ldUDU, ipivot, rwork, rworkSize, &info);
#endifarithm
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

int UDUSolve_@(pre)primme(@(type) *UDU, int *ipivot, int dim, @(type) *rhs, 
   @(type) *sol) {

   int info;

   if (dim == 1) {
#ifdefarithm L_DEFCPLX
      z_div_primme(sol, rhs, UDU); 
#endifarithm
#ifdefarithm L_DEFREAL
      *sol = *rhs/(*UDU); 
#endifarithm
      info = 0;
   }
   else {
      Num_@(pre)copy_@(pre)primme(dim, rhs, 1, sol, 1);
#ifdefarithm L_DEFCPLX
      Num_zhetrs_zprimme("U", dim, 1, UDU, dim, ipivot, sol, dim, &info);
#endifarithm
#ifdefarithm L_DEFREAL
      Num_dsytrs_dprimme("U", dim, 1, UDU, dim, ipivot, sol, dim, &info);
#endifarithm
   }

   return info;

}
