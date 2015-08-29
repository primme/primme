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
#include "primme.h"
#include "const.h"
#include "init_z.h"
#include "init_private_z.h"
#include "update_projection_z.h"
#include "update_W_z.h"
#include "ortho_z.h"
#include "factorize_z.h"
#include "numerical_z.h"
#include "wtime.h"                       /* Needed for CostModel */


/*******************************************************************************
 * subroutine init_basis - This subroutine is used to 
 *    initialize the basis V.
 *
 * Here we initialize the basis V as well as the array of locked
 * vectors if locking has been enabled.  The cases handled are:
 * 
 * I. Locking is disabled
 *    A. No initial vectors are provided.  BlockSize orthonormal initial
 *       vectors are created and used to form an orthonormal block Krylov
 *       basis.  The size of the basis will be minRestartSize.
 * 
 *    B. Initial vectors are provided.
 *      Subcases:
 *       1. minRestartSize or greater number of initial vectors are provided
 *          and they are orthonormalized by calling ortho.
 * 
 *       2. Fewer than minRestartSize initial vectors were provided.
 *          Additional vectors are placed in the basis by forming
 *          a block Krylov subspace.  The dimension of the block Krylov space
 *          will be minRestartSize - initSize.
 * 
 * II. Locking is enabled
 * 
 *     A. There are initSize initial vectors provided.  They are
 *        then orthonormalized and copied to V.
 * 
 *     B. There are fewer than minRestartSize initial vectors provided.
 *        A Krylov subspace of dimension restartSize - initSize vectors
 *        is created so that restartSize initial vectors will be available.
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

int init_basis_zprimme(Complex_Z *V, Complex_Z *W, Complex_Z *evecs, 
   Complex_Z *evecsHat, Complex_Z *M, Complex_Z *UDU, int *ipivot, 
   double machEps,  Complex_Z *rwork, int rworkSize, int *basisSize, 
   int *nextGuess, int *numGuesses, double *timeForMV,
   primme_params *primme) {

   int ret;          /* Return value                              */
   int currentSize;

   /*-----------------------------------------------------------------------*/
   /* Orthogonalize the orthogonalization constraints provided by the user. */
   /* If a preconditioner is given and inner iterations are to be           */
   /* performed, then initialize M.                                         */
   /*-----------------------------------------------------------------------*/

   if (primme->numOrthoConst > 0) {
      ret = ortho_zprimme(evecs, primme->nLocal, 0, 
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

         update_projection_zprimme(evecs, evecsHat, M, 0, 
            primme->numOrthoConst+primme->numEvals, primme->numOrthoConst, 
            rwork, primme);

         ret = UDUDecompose_zprimme(M, UDU, ipivot, primme->numOrthoConst, 
            rwork, rworkSize, primme);

         if (ret != 0) {
            primme_PushErrorMessage(Primme_init_basis, Primme_ududecompose, ret,
               __FILE__, __LINE__, primme);
            return UDUDECOMPOSE_FAILURE;
         }

      }  /* if evecsHat and M=evecs'evecsHat, UDU are needed */

   }  /* if numOrthoCont >0 */


   /*-----------------------------------------------------------------------*/
   /* No locking                                                            */
   /*-----------------------------------------------------------------------*/
   if (!primme->locking) {

      /* Handle case when no initial guesses are provided by the user */
      if (primme->initSize == 0) {

         ret = init_block_krylov(V, W, 0, primme->minRestartSize - 1, evecs, 
            primme->numOrthoConst, machEps, rwork, rworkSize, primme); 

         /* Push an error message onto the stack trace if an error occured */
         if (ret < 0) {
            primme_PushErrorMessage(Primme_init_basis, Primme_init_block_krylov,
                            ret, __FILE__, __LINE__, primme);
            return INIT_BLOCK_KRYLOV_FAILURE;
         }

         *basisSize = primme->minRestartSize;

      }
      else {
      /* Handle case when some or all initial guesses are provided by */ 
      /* the user                                                     */

         /* Copy over the initial guesses provided by the user */
         Num_zcopy_zprimme(primme->nLocal*primme->initSize, 
            &evecs[primme->numOrthoConst*primme->nLocal], 1, V, 1);

         /* Orthonormalize the guesses provided by the user */ 

         ret = ortho_zprimme(V, primme->nLocal, 0, primme->initSize-1, 
            evecs, primme->nLocal, primme->numOrthoConst, primme->nLocal, 
            primme->iseed, machEps, rwork, rworkSize, primme);

         /* Push an error message onto the stack trace if an error occured */
         if (ret < 0) {
            primme_PushErrorMessage(Primme_init_basis, Primme_ortho, ret, 
                            __FILE__, __LINE__, primme);
            return ORTHO_FAILURE;
         }

         update_W_zprimme(V, W, 0, primme->initSize, primme);

         /* An insufficient number of initial guesses were provided by */
         /* the user.  Generate a block Krylov space to fill the       */
         /* remaining vacancies.                                       */

         if (primme->initSize < primme->minRestartSize) {

            ret = init_block_krylov(V, W, primme->initSize, 
               primme->minRestartSize - 1, evecs, primme->numOrthoConst, 
               machEps, rwork, rworkSize, primme);

            /* Push an error message onto the stack trace if an error occured */
            if (ret < 0) {
               primme_PushErrorMessage(Primme_init_basis, 
                  Primme_init_block_krylov, ret, __FILE__, __LINE__, primme);
               return INIT_KRYLOV_FAILURE;
            }

            *basisSize = primme->minRestartSize;
         }
         else {
            *basisSize = primme->initSize;
         }

      }

      *numGuesses = 0;
      *nextGuess = 0;

   }
   else {
   /*-----------------------------------------------------------------------*/
   /* Locking                                                               */
   /*-----------------------------------------------------------------------*/

      *numGuesses = primme->initSize;
      *nextGuess = primme->numOrthoConst;

      /* If some initial guesses are available, copy them to the basis       */
      /* and orthogonalize them against themselves and the orthogonalization */
      /* constraints.                                                        */

      if (primme->initSize > 0) {
         currentSize = min(primme->initSize, primme->minRestartSize);
         Num_zcopy_zprimme(primme->nLocal*currentSize, 
            &evecs[primme->numOrthoConst*primme->nLocal], 1, V, 1);

         ret = ortho_zprimme(V, primme->nLocal, 0, currentSize-1, evecs,
            primme->nLocal, primme->numOrthoConst, primme->nLocal,
            primme->iseed, machEps, rwork, rworkSize, primme);

         if (ret < 0) {
            primme_PushErrorMessage(Primme_init_basis, Primme_ortho, ret,
                   __FILE__, __LINE__, primme);
            return ORTHO_FAILURE;
         }
      
         update_W_zprimme(V, W, 0, currentSize, primme);
         *numGuesses = *numGuesses - currentSize;
         *nextGuess = *nextGuess + currentSize;
         
      }
      else {
         currentSize = 0;
      }

      /* If an insufficient number of guesses was provided, then fill */
      /* the remaining vacancies with a block Krylov space.           */

      if (currentSize < primme->minRestartSize) {
         
         ret = init_block_krylov(V, W, currentSize, primme->minRestartSize - 1,
            evecs, primme->numOrthoConst, machEps, rwork, rworkSize, primme);

         if (ret < 0) {
            primme_PushErrorMessage(Primme_init_basis, Primme_init_block_krylov,
                            ret, __FILE__, __LINE__, primme);
            return INIT_BLOCK_KRYLOV_FAILURE;
         }

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
 * dv1, dv2    Range of indicies over which the basis will be generated
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

static int init_block_krylov(Complex_Z *V, Complex_Z *W, int dv1, int dv2, 
   Complex_Z *locked, int numLocked, double machEps, Complex_Z *rwork, 
   int rworkSize, primme_params *primme) {

   int i;               /* Loop variables */
   int numNewVectors;   /* Number of vectors to be generated */
   int ret;             /* Return code.                      */  
   int ONE = 1;         /* Used for passing it by reference in matrixmatvec */
   
   numNewVectors = dv2 - dv1 + 1;

   /*----------------------------------------------------------------------*/
   /* Generate a single Krylov space if there are only a few vectors to be */
   /* generated, else generate a block Krylov space with                   */
   /* primme->maxBlockSize as the block Size.                              */ 
   /*----------------------------------------------------------------------*/

   if (numNewVectors <= primme->maxBlockSize) {

      /* Create and orthogonalize the inital vectors */

      Num_larnv_zprimme(2, primme->iseed,primme->nLocal,&V[primme->nLocal*dv1]);
      ret = ortho_zprimme(V, primme->nLocal, dv1, dv1, locked, 
         primme->nLocal, numLocked, primme->nLocal, primme->iseed, machEps, 
         rwork, rworkSize, primme);

      if (ret < 0) {
         primme_PushErrorMessage(Primme_init_block_krylov, Primme_ortho, ret, 
            __FILE__, __LINE__, primme);
         return ORTHO_FAILURE;
      }

      /* Generate the remainder of the Krylov space. */

      for (i = dv1; i < dv2; i++) {
         (*primme->matrixMatvec)
           (&V[primme->nLocal*i], &V[primme->nLocal*(i+1)], &ONE, primme);
         Num_zcopy_zprimme(primme->nLocal, &V[primme->nLocal*(i+1)], 1,
            &W[primme->nLocal*i], 1);
         ret = ortho_zprimme(V, primme->nLocal, i+1, i+1, locked, 
            primme->nLocal, numLocked, primme->nLocal, primme->iseed, machEps,
            rwork, rworkSize, primme);
      
         if (ret < 0) {
            primme_PushErrorMessage(Primme_init_block_krylov, Primme_ortho, 
                            ret, __FILE__, __LINE__, primme);
            return ORTHO_FAILURE;
         }
      }

      primme->stats.numMatvecs += dv2-dv1;
      update_W_zprimme(V, W, dv2, 1, primme);

   }
   else {
   /*----------------------------------------------------------------------*/
   /* Generate the initial vectors.                                        */
   /*----------------------------------------------------------------------*/

      Num_larnv_zprimme(2, primme->iseed, primme->nLocal*primme->maxBlockSize,
         &V[primme->nLocal*dv1]);
      ret = ortho_zprimme(V, primme->nLocal, dv1, 
         dv1+primme->maxBlockSize-1, locked, primme->nLocal, numLocked, 
         primme->nLocal, primme->iseed, machEps, rwork, rworkSize, primme);

      /* Generate the remaining vectors in the sequence */

      for (i = dv1+primme->maxBlockSize; i <= dv2; i++) {
         (*primme->matrixMatvec)(&V[primme->nLocal*(i-primme->maxBlockSize)], 
            &V[primme->nLocal*i], &ONE, primme);
         Num_zcopy_zprimme(primme->nLocal, &V[primme->nLocal*i], 1,
            &W[primme->nLocal*(i-primme->maxBlockSize)], 1);

         ret = ortho_zprimme(V, primme->nLocal, i, i, locked, 
            primme->nLocal, numLocked, primme->nLocal, primme->iseed, machEps,
            rwork, rworkSize, primme);

         if (ret < 0) {
            primme_PushErrorMessage(Primme_init_block_krylov, Primme_ortho, 
                            ret, __FILE__, __LINE__, primme);
            return ORTHO_FAILURE;
         }

      }

      primme->stats.numMatvecs += dv2-(dv1+primme->maxBlockSize)+1;
      update_W_zprimme(V, W, dv2-primme->maxBlockSize+1, primme->maxBlockSize,
         primme);

   }
         
   return 0;
}
