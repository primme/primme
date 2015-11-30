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
 * File: init.c
 *
 * Purpose - Generate the basis that will be used during the first
 *           iteration of the method.
 *  
 ******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "primme.h"
#include "const.h"
#include "init_d.h"
#include "init_private_d.h"
#include "update_projection_d.h"
#include "update_W_d.h"
#include "ortho_d.h"
#include "factorize_d.h"
#include "numerical_d.h"
#include "wtime.h"                       /* Needed for CostModel */


/*******************************************************************************
 * subroutine init_basis - This subroutine is used to 
 *    initialize the basis V.
 *
 * Here we initialize the basis V as well as the array of locked
 * vectors if locking has been enabled.  The cases handled are:
 * 
 *    I. No initial vectors are provided.  BlockSize orthonormal initial
 *       vectors are created and used to form an orthonormal block Krylov
 *       basis.  The size of the basis will be minRestartSize.
 * 
 *   II. Initial vectors are provided.
 *
 *       1.a if locking is disabled, minRestartSize or greater number of initial
 *          vectors are provided and they are orthonormalized by calling ortho.
 * 
 *       1.b if locking is enable, up to minRestartSize of initial vectors are
 *          copied to V and orthonormalized. 
 *       2. There are fewer than minRestartSize initial vectors provided.
 *          A Krylov subspace of dimension restartSize - initSize vectors
 *          is created so that restartSize initial vectors will be available.
 * 
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * evecs      Array containing the orthogonalization constraints and initial 
 *            guesses. Holds as many as primme->numOrthoConst + primme->numEvals
 *            vectors.
 *
 * machEps    double machine precision
 * 
 * rwork      Double precision work array needed by other subroutines called
 *            by initialize_basis
 *
 * rworkSize  At most the maximum of (maximum size required by the 
 *            orthogonalization routine, maximum worksize by UDUDecompose)
 *
 * primme       Structure containing various solver parameters
 * 
 * 
 * 
 * OUTPUT PARAMETERS
 * -----------------
 * *basisSize   The size of the resulting basis V
 *
 * *nextGuess   The index of the next initial guess stored in the evecs array
 *
 * *numGuesses  When locking is enabled, the number of remaining initial guesses
 * 
 * *timeForMV   Time estimate for applying the matvec operator.
 *              Measured only if primme.dynamicMethodSwitch is on.
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * V            The orthonormal basis
 *
 * W            A*V
 *
 * evecsHat     K^{-1}*evecs, given a preconditioner K
 *
 * M            evecs'*evecsHat.  Its dimension is as large as 
 *              (primme->numOrthoConst + primme->numEvals).
 *
 * UDU          The factorization of M
 *
 * ipivot       The pivots of the UDU factorization
 *
 * Return value
 * ------------
 *  0 - Successful return
 * -1 - Orthogonalization failure
 * -2 - Failure to initialize block Krylov subspace
 * -3 - Failure to initialize Krylov subspace
 * -4 - Failure to UDU decompose M
 *
 ******************************************************************************/

int init_basis_dprimme(double *V, double *W, double *evecs, 
   double *evecsHat, double *M, double *UDU, int *ipivot, 
   double machEps,  double *rwork, int rworkSize, int *basisSize, 
   int *nextGuess, int *numGuesses, double *timeForMV,
   primme_params *primme) {

   int ret;          /* Return value                              */
   int currentSize;
   int initSize;
   int random;

   /*-----------------------------------------------------------------------*/
   /* Orthogonalize the orthogonalization constraints provided by the user. */
   /* If a preconditioner is given and inner iterations are to be           */
   /* performed, then initialize M.                                         */
   /*-----------------------------------------------------------------------*/

   if (primme->numOrthoConst > 0) {
   /* lingfei: primme_svds. change ortho function for returning Q and R */
      ret = ortho_dprimme(evecs, NULL, primme->nLocal, 0, 
        primme->numOrthoConst - 1, NULL, 0, 0, primme->nLocal, 
        primme->iseed, machEps, rwork, rworkSize, primme);

      /* Push an error message onto the stack trace if an error occured */
      if (ret < 0) {
         primme_PushErrorMessage(Primme_init_basis, Primme_ortho, ret, 
                         __FILE__, __LINE__, primme);
         return ORTHO_FAILURE;
      }

      /* Initialize evecsHat, M, and its factorization UDU,ipivot. This   */
      /* allows the orthogonalization constraints to be included in the   */
      /* projector (I-QQ'). Only needed if there is preconditioning, and  */
      /* JDqmr inner iterations with a right, skew projector. Only in     */
      /* that case, is UDU not NULL                                       */

      if (UDU != NULL) {

         (*primme->applyPreconditioner)
            (evecs, evecsHat, &primme->numOrthoConst, primme); 
         primme->stats.numPreconds += primme->numOrthoConst;

         update_projection_dprimme(evecs, evecsHat, M, 0, 
            primme->numOrthoConst+primme->numEvals, primme->numOrthoConst, 
            rwork, primme);

         ret = UDUDecompose_dprimme(M, UDU, ipivot, primme->numOrthoConst, 
            rwork, rworkSize, primme);

         if (ret != 0) {
            primme_PushErrorMessage(Primme_init_basis, Primme_ududecompose, ret,
               __FILE__, __LINE__, primme);
            return UDUDECOMPOSE_FAILURE;
         }

      }  /* if evecsHat and M=evecs'evecsHat, UDU are needed */

   }  /* if numOrthoCont >0 */


   /* Handle case when some or all initial guesses are provided by */ 
   /* the user                                                     */
   if (!primme->locking) {
      initSize = primme->initSize;
   }
   else {
      initSize = min(primme->minRestartSize, primme->initSize);
   }
   *numGuesses = primme->initSize - initSize;
   *nextGuess = primme->numOrthoConst + initSize;

   /* Copy over the initial guesses provided by the user */
   Num_dcopy_dprimme(primme->nLocal*initSize, 
         &evecs[primme->numOrthoConst*primme->nLocal], 1, V, 1);

   switch(primme->InitBasisMode) {
   case primme_init_krylov:
      random = 0;
      break;
   case primme_init_random:
      random = max(0,primme->minRestartSize-initSize);
      break;
   case primme_init_user:
      random = max(primme->maxBlockSize-initSize, 0);
      break;
   default:
      assert(0);
   }
   Num_larnv_dprimme(2, primme->iseed, primme->nLocal*random,
      &V[primme->nLocal*initSize]);
   *basisSize = initSize + random;

   /* Orthonormalize the guesses provided by the user */ 
   ret = ortho_dprimme(V, NULL, primme->nLocal, 0, *basisSize-1, 
         evecs, primme->nLocal, primme->numOrthoConst, primme->nLocal, 
         primme->iseed, machEps, rwork, rworkSize, primme);

   /* Push an error message onto the stack trace if an error occurred */
   if (ret < 0) {
      primme_PushErrorMessage(Primme_init_basis, Primme_ortho, ret, 
            __FILE__, __LINE__, primme);
      return ORTHO_FAILURE;
   }

   update_W_dprimme(V, W, 0, *basisSize, primme);

   if (primme->InitBasisMode == primme_init_krylov) {
      ret = init_block_krylov(V, W, *basisSize, primme->minRestartSize - 1, evecs, 
            primme->numOrthoConst, machEps, rwork, rworkSize, primme); 

      /* Push an error message onto the stack trace if an error occurred */
      if (ret < 0) {
         primme_PushErrorMessage(Primme_init_basis, Primme_init_block_krylov,
               ret, __FILE__, __LINE__, primme);
         return INIT_BLOCK_KRYLOV_FAILURE;
      }

      *basisSize = primme->minRestartSize;
   }

   /* ----------------------------------------------------------- */
   /* If time measurements are needed, waste one MV + one Precond */
   /* Put dummy results in the first open space of W (currentSize)*/
   /* ----------------------------------------------------------- */
   if (primme->dynamicMethodSwitch) {
      currentSize = primme->nLocal*(*basisSize);
      ret = 1;
      *timeForMV = primme_wTimer(0);
       (*primme->matrixMatvec)(V, &W[currentSize], &ret, primme);
      *timeForMV = primme_wTimer(0) - *timeForMV;
      primme->stats.numMatvecs += 1;
   }
      
   return 0;
}


/*******************************************************************************
 * Subroutine init_block_krylov - Initializes the basis as an orthonormal 
 *    block Krylov subspace.  
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * dv1, dv2    Range of indices over which the basis will be generated
 * 
 * locked      The array of locked Ritz vectors
 * 
 * numLocked   The number of vectors in the locked array
 *
 * machEps     machine precision needed in ortho()
 *
 * rwork       Real work array
 *
 * rworkSize   Size of rwork array
 *
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * V  The orthonormal basis
 * 
 * W  A*V
 *
 * Return value
 * ------------
 * int -  0 upon success
 *       -1 if orthogonalization failed
 * 
 ******************************************************************************/

static int init_block_krylov(double *V, double *W, int dv1, int dv2, 
   double *locked, int numLocked, double machEps, double *rwork, 
   int rworkSize, primme_params *primme) {

   int i;               /* Loop variables */
   int numNewVectors;   /* Number of vectors to be generated */
   int ret;             /* Return code.                      */  
   int ONE = 1;         /* Used for passing it by reference in matrixmatvec */
   int blockSize;       /* blockSize used in practice */
   
   numNewVectors = dv2 - dv1 + 1;

   /*----------------------------------------------------------------------*/
   /* Generate a single Krylov space if there are only a few vectors to be */
   /* generated, else generate a block Krylov space with                   */
   /* primme->maxBlockSize as the block Size.                              */ 
   /*----------------------------------------------------------------------*/

   blockSize = numNewVectors <= primme->maxBlockSize ? 1 : primme->maxBlockSize;

   /*----------------------------------------------------------------------*/
   /* Generate the initial vectors.                                        */
   /*----------------------------------------------------------------------*/

   if (dv1+blockSize-1 <= dv2) {
      Num_larnv_dprimme(2, primme->iseed, primme->nLocal*blockSize,
         &V[primme->nLocal*dv1]);
   }
   ret = ortho_dprimme(V, NULL, primme->nLocal, dv1, 
      dv1+blockSize-1, locked, primme->nLocal, numLocked, 
      primme->nLocal, primme->iseed, machEps, rwork, rworkSize, primme);

   /* Generate the remaining vectors in the sequence */

   for (i = dv1+blockSize; i <= dv2; i++) {
      (*primme->matrixMatvec)(&V[primme->nLocal*(i-blockSize)], 
         &V[primme->nLocal*i], &ONE, primme);
      Num_dcopy_dprimme(primme->nLocal, &V[primme->nLocal*i], 1,
         &W[primme->nLocal*(i-blockSize)], 1);

       /* lingfei: primme_svds. change ortho function for returning Q and R */
       ret = ortho_dprimme(V, NULL, primme->nLocal, i, i, locked, 
         primme->nLocal, numLocked, primme->nLocal, primme->iseed, machEps,
         rwork, rworkSize, primme);

      if (ret < 0) {
         primme_PushErrorMessage(Primme_init_block_krylov, Primme_ortho, 
                         ret, __FILE__, __LINE__, primme);
         return ORTHO_FAILURE;
      }
   }

   primme->stats.numMatvecs += dv2-(dv1+blockSize)+1;
   update_W_dprimme(V, W, dv2-blockSize+1, blockSize, primme);

   return 0;
}
