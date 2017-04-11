/*******************************************************************************
 * Copyright (c) 2017, College of William & Mary
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the College of William & Mary nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * PRIMME: https://github.com/primme/primme
 * Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
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
#include "numerical.h"
#include "init.h"
#include "update_projection.h"
#include "update_W.h"
#include "ortho.h"
#include "factorize.h"
#include "auxiliary_eigs.h"
#include "wtime.h"                       /* Needed for CostModel */

static int init_block_krylov(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *W, PRIMME_INT ldW, int dv1, int dv2, SCALAR *locked,
      PRIMME_INT ldlocked, int numLocked, double machEps, SCALAR *rwork,
      size_t *rworkSize, primme_params *primme);

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
 *  error code
 ******************************************************************************/

TEMPLATE_PLEASE
int init_basis_Sprimme(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *W, PRIMME_INT ldW, SCALAR *evecs, PRIMME_INT ldevecs,
      SCALAR *evecsHat, PRIMME_INT ldevecsHat, SCALAR *M, int ldM, SCALAR *UDU,
      int ldUDU, int *ipivot, double machEps, SCALAR *rwork, size_t *rworkSize,
      int *basisSize, int *nextGuess, int *numGuesses, primme_params *primme) {

   int i;
   int initSize;
   int random=0;

   /* Return memory requirement */

   if (V == NULL) {
      update_projection_Sprimme(NULL, 0, NULL, 0, NULL, 0, nLocal,
            0, primme->numOrthoConst, NULL, rworkSize, 1/*symmetric*/, primme);
      UDUDecompose_Sprimme(NULL, 0, NULL, 0, NULL,
            primme->numOrthoConst, NULL, rworkSize, primme);
      ortho_Sprimme(NULL, 0, NULL, 0, 0, 
            primme->numOrthoConst-1, NULL, 0, 0, nLocal, 
            NULL, 0.0, NULL, rworkSize, primme);
      ortho_Sprimme(NULL, 0, NULL, 0, 0, *basisSize-1, 
            NULL, 0, primme->numOrthoConst, nLocal, 
            NULL, 0.0, NULL, rworkSize, primme);
      return 0;
   }

   /*-----------------------------------------------------------------------*/
   /* Orthogonalize the orthogonalization constraints provided by the user. */
   /* If a preconditioner is given and inner iterations are to be           */
   /* performed, then initialize M.                                         */
   /*-----------------------------------------------------------------------*/

   if (primme->numOrthoConst > 0) {
      CHKERR(ortho_Sprimme(evecs, ldevecs, NULL, 0, 0, 
        primme->numOrthoConst - 1, NULL, 0, 0, nLocal, 
        primme->iseed, machEps, rwork, rworkSize, primme), -1);

      /* Initialize evecsHat, M, and its factorization UDU,ipivot. This   */
      /* allows the orthogonalization constraints to be included in the   */
      /* projector (I-QQ'). Only needed if there is preconditioning, and  */
      /* JDqmr inner iterations with a right, skew projector. Only in     */
      /* that case, is UDU not NULL                                       */

      if (UDU != NULL) {

         CHKERR(applyPreconditioner_Sprimme(evecs, primme->nLocal, ldevecs,
                  evecsHat, ldevecsHat, primme->numOrthoConst, primme), -1);

         CHKERR(update_projection_Sprimme(evecs, ldevecs, evecsHat,
                  ldevecsHat, M, ldM, nLocal, 0, primme->numOrthoConst, rwork,
                  rworkSize, 1/*symmetric*/, primme), -1);

         CHKERR(UDUDecompose_Sprimme(M, ldM, UDU, ldUDU, ipivot,
                  primme->numOrthoConst, rwork, rworkSize, primme), -1);

      }  /* if evecsHat and M=evecs'evecsHat, UDU are needed */

   }  /* if numOrthoCont >0 */


   /* Handle case when some or all initial guesses are provided by */ 
   /* the user                                                     */
   if (!primme->locking) {
      initSize = min(primme->maxBasisSize, primme->initSize);
   }
   else {
      initSize = min(primme->minRestartSize, primme->initSize);
   }
   *numGuesses = primme->initSize - initSize;
   *nextGuess = primme->numOrthoConst + initSize;

   /* Copy over the initial guesses provided by the user */
   Num_copy_matrix_Sprimme(&evecs[primme->numOrthoConst*ldevecs],
         nLocal, initSize, ldevecs, V, ldV);

   switch(primme->initBasisMode) {
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
   for (i=0; i<random; i++) {
      Num_larnv_Sprimme(2, primme->iseed, nLocal,
            &V[ldV*(initSize+i)]);
   }
   *basisSize = initSize + random;

   /* Orthonormalize the guesses provided by the user */ 
   CHKERR(ortho_Sprimme(V, ldV, NULL, 0, 0, *basisSize-1, 
         evecs, ldevecs, primme->numOrthoConst, nLocal, 
         primme->iseed, machEps, rwork, rworkSize, primme), -1)

   CHKERR(matrixMatvec_Sprimme(V, nLocal, ldV, W, ldW, 0, *basisSize,
            primme), -1);

   if (primme->initBasisMode == primme_init_krylov) {
      CHKERR(init_block_krylov(V, nLocal, ldV, W, ldW, *basisSize,
            primme->minRestartSize-1, evecs, ldevecs, primme->numOrthoConst,
            machEps, rwork, rworkSize, primme), -1); 

      *basisSize = primme->minRestartSize;
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

static int init_block_krylov(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *W, PRIMME_INT ldW, int dv1, int dv2, SCALAR *locked,
      PRIMME_INT ldlocked, int numLocked, double machEps, SCALAR *rwork,
      size_t *rworkSize, primme_params *primme) {

   int i;               /* Loop variables */
   int numNewVectors;   /* Number of vectors to be generated */
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
      for (i=dv1; i<dv1+blockSize; i++) {
         Num_larnv_Sprimme(2, primme->iseed, nLocal, &V[ldV*i]);
      }
   }
   CHKERR(ortho_Sprimme(V, ldV, NULL, 0, dv1, 
            dv1+blockSize-1, locked, ldlocked, numLocked, 
            nLocal, primme->iseed, machEps, rwork, rworkSize, primme), -1);

   /* Generate the remaining vectors in the sequence */

   for (i = dv1+blockSize; i <= dv2; i++) {
      CHKERR(matrixMatvec_Sprimme(&V[ldV*(i-blockSize)], nLocal, ldV,
               &V[ldV*i], ldV, 0, 1, primme), -1);

      Num_copy_Sprimme(nLocal, &V[ldV*i], 1,
         &W[ldW*(i-blockSize)], 1);

      CHKERR(ortho_Sprimme(V, ldV, NULL, 0, i, i, locked, 
               ldlocked, numLocked, nLocal, primme->iseed, machEps,
               rwork, rworkSize, primme), -1);
   }

   CHKERR(matrixMatvec_Sprimme(V, nLocal, ldV, W, ldW, dv2-blockSize+1,
            blockSize, primme), -1);

   return 0;
}
