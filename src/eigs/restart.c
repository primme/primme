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
 * File: restart.c
 *
 * Purpose - Restart V and related matrices (eg. W, H, Q, R, QtV...).
 *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "const.h"
#include "numerical.h"
#include "auxiliary_eigs.h"
#include "restart.h"
#include "ortho.h"
#include "solve_projection.h"
#include "factorize.h"
#include "update_projection.h"
#include "update_W.h"
#include "convergence.h"
#include "globalsum.h"
#include "wtime.h"

static int restart_soft_locking_Sprimme(int *restartSize, SCALAR *V,
       SCALAR *W, PRIMME_INT nLocal, int basisSize, PRIMME_INT ldV, SCALAR **X,
       SCALAR **R, SCALAR *hVecs, int ldhVecs, int *restartPerm,
       REAL *hVals, int *flags, int *iev, int *ievSize, REAL *blockNorms,
       SCALAR *evecs, REAL *evals, REAL *resNorms, SCALAR *evecsHat,
       PRIMME_INT ldevecsHat, SCALAR *M, int ldM, int *numConverged,
       int *numConvergedStored, int numPrevRetained, int *indexOfPreviousVecs,
       int *hVecsPerm, int reset, double machEps, SCALAR *rwork,
       size_t *rworkSize, int *iwork, int iworkSize, primme_params *primme);

static int restart_locking_Sprimme(int *restartSize, SCALAR *V, SCALAR *W, 
      PRIMME_INT nLocal, int basisSize, PRIMME_INT ldV, SCALAR **X, SCALAR **R,
      SCALAR *hVecs, int ldhVecs, int *restartPerm, REAL *hVals, int *flags,
      int *iev, int *ievSize, REAL *blockNorms, SCALAR *evecs,
      PRIMME_INT ldevecs, REAL *evals, int *numConverged, int *numLocked,
      REAL *resNorms, int *lockedFlags, int *evecsperm, int numPrevRetained,
      int *indexOfPreviousVecs, int *hVecsPerm, int reset, double machEps,
      SCALAR *rwork, size_t *rworkSize, int *iwork, int iworkSize,
      primme_params *primme);

static int restart_projection_Sprimme(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *H, int ldH, SCALAR *VtBV, int ldVtBV, SCALAR *Q,
      PRIMME_INT ldQ, PRIMME_INT nLocal, SCALAR *R, int ldR, SCALAR *QtV,
      int ldQtV, SCALAR *hU, int ldhU, int newldhU,
      int indexOfPreviousVecsBeforeRestart, SCALAR *hVecs, int ldhVecs,
      int newldhVecs, REAL *hVals, REAL *hSVals, int *restartPerm,
      int *hVecsPerm, int restartSize, int basisSize, int numPrevRetained,
      int indexOfPreviousVecs, SCALAR *evecs, int *evecsSize,
      PRIMME_INT ldevecs, SCALAR *evecsHat, PRIMME_INT ldevecsHat, SCALAR *M,
      int ldM, SCALAR *UDU, int ldUDU, int *ipivot, int *targetShiftIndex,
      int numConverged, int *numArbitraryVecs, SCALAR *hVecsRot, int ldhVecsRot,
      size_t *rworkSize, SCALAR *rwork, int iworkSize, int *iwork,
      double machEps, primme_params *primme);

static int restart_RR(SCALAR *H, int ldH, SCALAR *hVecs, int ldhVecs,
      int newldhVecs, REAL *hVals, int restartSize,
      int basisSize, int numLocked, int numPrevRetained,
      int indexOfPreviousVecs, int *hVecsPerm, int *targetShiftIndex,
      double machEps, size_t *rworkSize, SCALAR *rwork, int iworkSize,
      int *iwork, primme_params *primme);

static int restart_RR_explicit(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, PRIMME_INT nLocal, SCALAR *H, int ldH, SCALAR *VtBV,
      int ldVtBV, SCALAR *hVecs, int ldhVecs, REAL *hVals, int restartSize,
      int numLocked, int *targetShiftIndex, double machEps, size_t *rworkSize,
      SCALAR *rwork, int iworkSize, int *iwork, primme_params *primme);

static int restart_refined(SCALAR *V, PRIMME_INT ldV, SCALAR *W, PRIMME_INT ldW,
      SCALAR *H, int ldH, SCALAR *Q, PRIMME_INT ldQ, PRIMME_INT nLocal,
      SCALAR *R, int ldR, SCALAR *hU, int ldhU, int newldhU,
      int indexOfPreviousVecsBeforeRestart, SCALAR *hVecs, int ldhVecs,
      int newldhVecs, REAL *hVals, REAL *hSVals, int *restartPerm,
      int *hVecsPerm, int restartSize, int basisSize, int numPrevRetained,
      int indexOfPreviousVecs, int *targetShiftIndex, int numConverged,
      int *numArbitraryVecs, SCALAR *hVecsRot, int ldhVecsRot,
      size_t *rworkSize, SCALAR *rwork, int iworkSize, int *iwork,
      double machEps, primme_params *primme);

static int restart_harmonic(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *H, int ldH, SCALAR *Q, PRIMME_INT nLocal,
      PRIMME_INT ldQ, SCALAR *R, int ldR, SCALAR *QtV, int ldQtV, SCALAR *hU,
      int ldhU, int newldhU, SCALAR *hVecs, int ldhVecs, int newldhVecs,
      REAL *hVals, REAL *hSVals, int *restartPerm, int *hVecsPerm,
      int restartSize, int basisSize, int numPrevRetained, 
      int indexOfPreviousVecs, int *targetShiftIndex, int numConverged,
      int *numArbitraryVecs, SCALAR *hVecsRot, int ldhVecsRot,
      size_t *rworkSize, SCALAR *rwork, int iworkSize, int *iwork,
      double machEps, primme_params *primme);

static int dtr_Sprimme(int numLocked, SCALAR *hVecs, REAL *hVals, int *flags, 
  int basisSize, int numFree, int *iev, SCALAR *rwork, primme_params *primme);

static int ortho_coefficient_vectors_Sprimme(SCALAR *hVecs, int basisSize,
      int ldhVecs, int indexOfPreviousVecs, SCALAR *hU, int ldhU, SCALAR *R,
      int ldR, SCALAR *VtBV, int ldVtBV, int *numPrevRetained, double machEps,
      SCALAR *rwork, size_t *rworkSize, primme_params *primme);

static void insertionSort(REAL newVal, REAL *evals, REAL newNorm,
   REAL *resNorms, int newFlag, int *flags, int *perm, int numLocked,
   primme_params *primme);

static int compute_residual_columns(PRIMME_INT m, REAL *evals, SCALAR *x,
      int n, int *p, PRIMME_INT ldx, SCALAR *Ax, PRIMME_INT ldAx,
      SCALAR *xo, int no, PRIMME_INT ldxo, int io0, SCALAR *ro, PRIMME_INT ldro,
      SCALAR *xd, int nd, int *pd, PRIMME_INT ldxd, SCALAR *rd, PRIMME_INT ldrd,
      SCALAR *rwork, PRIMME_INT lrwork);


/*******************************************************************************
 * Subroutine: restart - This routine replaces V with V*c, some subset
 *             of the Ritz vectors of the current and the previous iteration.
 *             Related bases and matrices are updated accordingly.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 *
 * nLocal           Number of rows of V, W, Q, evecs and evecsHat assigned to the node
 *
 * ldV              The leading dimension of V, W, Q, evecs and evecsHat
 *
 * basisSize        The number of columns in V, W and Q
 *
 * resNorms         The residual norms of the converged eigenpairs
 *
 * previousHVecs    Coefficient vectors retained from the previous iteration
 *
 * ldpreviousHVecs  The leading dimension of previousHVecs
 *
 * numGuesses       Number of remaining initial guesses
 *
 * rwork            Real work array
 *
 * rworkSize        Must be of size 
 *
 * iwork            Integer work array
 *                  
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSizeOutput Output the number of columns of the restarted V.
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * W                A*V
 *
 * hU               The left singular vectors of R or the eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of the input hU
 *
 * newldhU          The leading dimension of the output hU
 *
 * hVecs            The coefficient vectors
 *
 * ldhVecs          The leading dimension of the input hVecs
 *
 * newldhVecs       The leading dimension of the output hVecs
 *
 * hVals            The Rayleigh quotient of candidates
 *
 * hSVals           The singular values of R
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * iev              Array of size blockSize indicating which Ritz vectors are
 *                  targeted in the block
 *
 * ievSize          The length of iev
 *
 * blockNorms       The residual norms of the eigenpairs in the block
 *
 * evecs            The converged Ritz vectors. Without locking, all converged
 *                  eigenvectors are copied from V to evecs if skew projections
 *                  are required
 *
 * evals            The converged Ritz values
 *
 * resNorms         The residual norms of the converged eigenpairs
 *
 * evecsHat         K^{-1}evecs
 *
 * prevRitzVals     Projected values retained from the previous iteration
 *
 * numPrevRitzVals  Length of the vector prevRitzVals
 *
 * ldevecsHat       The leading dimension of evecsHat
 *
 * M, ldM           evecs'*evecsHat and the leading dimension of M
 *
 * UDU              The factorization of M
 *
 * ldUDU            The leading dimension of UDU
 *
 * ipivot           The pivot array of the UDU factorization
 *
 * H                The projection V'*A*V
 *
 * ldH              The leading dimension of H
 *
 * VtBV             The projection V'*B*V
 *
 * ldVtBV           The leading dimension of VtBV
 *
 * Q, R             The factors of the QR decomposition of (A - targetShift*B)*V
 *
 * ldQ, ldR         The leading dimension of Q and R
 *
 * numConverged     The number of converged eigenpairs
 *
 * numLocked        The number of locked eigenpairs
 *
 * lockedFlags      The flags of the locked pairs
 *
 * numConvergedStored The # of converged vectors copied to evecs
 *
 * numPrevRetained  As input the number of columns of previousHVecs. As output the
 *                  number of columns added to V
 *
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 *
 * numArbitraryVecs The number of columns of hVecsRot
 *
 * hVecsRot         hVecs = hV*hVecsRot, where hV are the original coefficient
 *                  vectors returned by solve_H
 *
 * ldhVecsRot       The leading dimension of hVecsRot
 *
 * restartsSinceReset Number of restarts since last reset of V and W
 *
 * reset            flag to reset V and W at this restart
 *
 *
 * Return value
 * ------------
 * int   > 0 the restart size   
 *        -2 restart_H failed
 *        -4 factorization of M failed
 *        -5 flags do not correspond to converged pairs in pseudolocking
 *       
 ******************************************************************************/
 
TEMPLATE_PLEASE
int restart_Sprimme(SCALAR *V, SCALAR *W, PRIMME_INT nLocal, int basisSize,
       PRIMME_INT ldV, REAL *hVals, REAL *hSVals, int *flags, int *iev,
       int *ievSize, REAL *blockNorms, SCALAR *evecs, PRIMME_INT ldevecs,
       int *evecsPerm, REAL *evals, REAL *resNorms, SCALAR *evecsHat,
       PRIMME_INT ldevecsHat, SCALAR *M, int ldM, SCALAR *UDU, int ldUDU,
       int *ipivot, int *numConverged, int *numLocked, int *lockedFlags,
       int *numConvergedStored, SCALAR *previousHVecs, int *numPrevRetained,
       int ldpreviousHVecs, int numGuesses, REAL *prevRitzVals,
       int *numPrevRitzVals, SCALAR *H, int ldH, SCALAR *VtBV, int ldVtBV,
       SCALAR *Q, PRIMME_INT ldQ, SCALAR *R, int ldR, SCALAR* QtV, int ldQtV,
       SCALAR *hU, int ldhU, int newldhU, SCALAR *hVecs, int ldhVecs,
       int newldhVecs, int *restartSizeOutput, int *targetShiftIndex,
       int *numArbitraryVecs, SCALAR *hVecsRot, int ldhVecsRot,
       int *restartsSinceReset, int *reset, double machEps, SCALAR *rwork,
       size_t *rworkSize, int *iwork, int iworkSize, primme_params *primme) {

   int i;                   /* Loop indices */
   int restartSize;         /* Basis size after restarting                   */
   int indexOfPreviousVecs; /* Column index in hVecs with previous vecs      */
   int *iwork0;             /* Temporal integer workspace pointer            */
   int iworkSize0;          /* Size of iwork0                                */
   int *restartPerm;        /* Permutation of hVecs used to restart V        */
   int *hVecsPerm;          /* Permutation of hVecs to sort as primme.target */
   int indexOfPreviousVecsBeforeRestart=0;/* descriptive enough name, isn't? */
   double aNorm = primme?max(primme->aNorm, primme->stats.estimateLargestSVal):0.0;

   /* Return memory requirement */

   if (V == NULL) {
      iworkSize0 = 0;
      CHKERR(ortho_coefficient_vectors_Sprimme(NULL, basisSize, 0,
            basisSize, NULL, 0, NULL, 0, VtBV, 0, numPrevRetained, 0.0, NULL,
            rworkSize, primme), -1);

      if (primme->locking) {
         CHKERR(restart_locking_Sprimme(&basisSize, NULL, NULL, nLocal,
                  basisSize, 0, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL,
                  ievSize, NULL, NULL, 0, NULL, numConverged, numConverged,
                  NULL, NULL, NULL, *numPrevRetained, NULL, NULL, 0, 0.0, NULL,
                  rworkSize, &iworkSize0, 0, primme), -1);
      }
      else {
         CHKERR(restart_soft_locking_Sprimme(&basisSize, NULL, NULL,
               nLocal, basisSize, 0, NULL, NULL, NULL, 0, NULL, NULL, NULL,
               NULL, ievSize, NULL, NULL, NULL, NULL, evecsHat, 0, NULL, 0,
               numConverged, numConverged, *numPrevRetained, NULL, NULL, 0, 0.0,
               NULL, rworkSize, &iworkSize0, 0, primme), -1);
      }

      CHKERR(restart_projection_Sprimme(NULL, 0, NULL, 0, NULL, 0, VtBV, 0,
               NULL, 0, 0, NULL, 0, NULL, 0, NULL, 0, 0, 0,
               NULL, 0, 0, NULL, NULL, NULL, NULL, basisSize, basisSize,
               *numPrevRetained, basisSize, NULL,
               NULL, 0, NULL, 0, NULL, 0, NULL, 0, NULL, 0, 0, NULL, NULL, 0,
               rworkSize, NULL, 0, &iworkSize0, 0.0, primme), -1);

      iworkSize0 += 2*basisSize; /* for restartPerm and hVecsPerm */
      *iwork = max(*iwork, iworkSize0);

      return 0;
   }

   /* ----------------------------------------------------------- */
   /* Remove the SKIP_UNTIL_RESTART flags.                        */
   /* ----------------------------------------------------------- */

   for (i=0, *numConverged=*numLocked; i<basisSize; i++) {
      if (flags[i] == SKIP_UNTIL_RESTART) {
         flags[i] = UNCONVERGED;
      }
      else if (flags[i] != UNCONVERGED &&
            /* Don't check more than numEvals */
               *numConverged < primme->numEvals &&
            /* Check only the first pairs, except if finding closest_leq/geq  */
            /* because with refined extraction the pairs may not be ordered   */
            /* by this criterion.                                             */
               (i < primme->numEvals-*numLocked
                || primme->target == primme_closest_geq
                || primme->target == primme_closest_leq)) {
         (*numConverged)++;
      }
   }

   /* ----------------------------------------------------------- */
   /* Special case: If (basisSize+numLocked) is the entire space, */
   /* then everything should be converged. Do not test, just flag */
   /* everything as converged. But only do that when multiple     */
   /* shifts are not involved.                                    */
   /* ----------------------------------------------------------- */

   if (basisSize + *numLocked + primme->numOrthoConst >= primme->n) {
      if (primme->target == primme_largest || primme->target == primme_smallest
            || (primme->target != primme_closest_leq
               && primme->target != primme_closest_geq
               && *numLocked >= primme->numTargetShifts-1)) {
         for (i = 0; i < basisSize; i++)
            flags[i] = CONVERGED;
         *numConverged = basisSize + *numLocked;
      }
      restartSize = basisSize;
      *numPrevRetained = 0;
   }
   /* --------------------------------------------------------------------- */
   /* If basis isn't full, restart with the current basis size.             */
   /* --------------------------------------------------------------------- */
   else if (basisSize <= primme->maxBasisSize - primme->maxBlockSize) {
      restartSize = basisSize;
      *numPrevRetained = 0;
   }
   /* --------------------------------------------------------------------- */
   /* If dynamic thick restarting is to be used, then determine the minimum */
   /* number of free spaces to be maintained and call the DTR routine.      */
   /* The DTR routine will determine how many coefficient vectors from the  */
   /* left and right of H-spectrum to retain at restart. If DTR is not used */
   /* then set the restart size to the minimum restart size.                */
   /* --------------------------------------------------------------------- */
   else if (primme->restartingParams.scheme == primme_dtr) {
      int numFree = *numPrevRetained+max(3, primme->maxBlockSize);
      restartSize = dtr_Sprimme(*numLocked, hVecs, hVals, flags, basisSize, numFree, 
            iev, rwork, primme);
   }
   else {
      restartSize = min(basisSize, primme->minRestartSize);
   }

   /* --------------------------------------------------------------------- */
   /* Experimentally observed that ||V'*V - I|| grows with the number of    */
   /* restarts as much as sqrt(restarts)*machEps, and ||A*V - W|| as        */
   /* sqrt(restarts)*machEps*aNorm. For now the last norm is tracked here   */
   /* and used by check_convergence to consider the error computing the     */
   /* residual norm. V and W are asked to be reset when the error is as     */
   /* much as the current residual norm. If using refined, only W is        */
   /* reset: we haven't seen any benefit by resetting V also.               */
   /* --------------------------------------------------------------------- */

   if (!*reset) {
      ++*restartsSinceReset;
   }
   else {
      *restartsSinceReset = 0;
      if (!Q) *reset = 2; /* if Q, only reset W, not V */
      if (primme->printLevel >= 5 && primme->procID == 0) {
         fprintf(primme->outputFile, 
               "Resetting V, W and QR.\n");
         fflush(primme->outputFile);
      }
   }

   /* ----------------------------------------------------------------------- */
   /* Insert as many initial guesses as eigenpairs have converged.            */
   /* Leave sufficient restarting room in the restarted basis so that to      */
   /* insert (in main_iter) as many initial guesses as the number of          */
   /* eigenpairs that converged.                                              */
   /* ----------------------------------------------------------------------- */

   restartSize -= min(min(numGuesses, *numConverged-*numLocked), restartSize);

   /* ----------------------------------------------------------------------- */
   /* Limit restartSize so that it plus 'to be locked' plus previous Ritz     */
   /* vectors do not exceed basisSize.                                        */
   /* ----------------------------------------------------------------------- */

   if (primme->locking)
      restartSize = min(restartSize, basisSize-(*numConverged-*numLocked));

   /* ----------------------------------------------------------------------- */
   /* Limit the number of previous retained vectors such that the final basis */
   /* size isn't larger than the current basis size.                          */
   /* ----------------------------------------------------------------------- */

   *numPrevRetained = max(0, min(min(
         *numPrevRetained,
         restartSize - (*numConverged-*numLocked)), 
         basisSize - (restartSize+*numConverged-*numLocked)));

   /* ----------------------------------------------------------------------- */
   /* Restarting with a small number of coefficient vectors from the previous */
   /* iteration can accelerate convergence.  The previous                     */
   /* coefficient vectors must be combined with the current coefficient       */
   /* vectors by first orthogonalizing the previous ones versus the current   */
   /* restartSize ones.  The orthogonalized previous vectors are then         */
   /* inserted into the hVecs array at hVecs(:,indexOfPreviousVecs).          */
   /* When using refined, avoid to replace linear dependent previous coeffi-  */
   /* cient vectors by random vectors, it may increase the norm of R after    */
   /* restarting.                                                             */
   /* ----------------------------------------------------------------------- */

   if (primme->locking)
      indexOfPreviousVecs = restartSize+*numConverged-*numLocked;
   else
      indexOfPreviousVecs = restartSize;
   indexOfPreviousVecsBeforeRestart = indexOfPreviousVecs;

   Num_copy_matrix_Sprimme(previousHVecs, basisSize, *numPrevRetained,
         ldpreviousHVecs, &hVecs[ldhVecs*indexOfPreviousVecs], ldhVecs);

   CHKERR(ortho_coefficient_vectors_Sprimme(hVecs, basisSize, ldhVecs,
         indexOfPreviousVecs, hU, ldhU, R, ldR, VtBV, ldVtBV, numPrevRetained,
         machEps, rwork, rworkSize, primme), -1);

   /* ----------------------------------------------------------------------- */
   /* Restart V and W, and compute X and residual vectors for next candidates */
   /* ----------------------------------------------------------------------- */

   restartPerm = iwork;
   hVecsPerm = &restartPerm[basisSize];
   iwork0 = &hVecsPerm[basisSize];
   iworkSize0 = iworkSize - 2*basisSize;

   if (!primme->locking) {
      SCALAR *X, *Res;
      CHKERR(restart_soft_locking_Sprimme(&restartSize, V, W, nLocal,
               basisSize, ldV, &X, &Res, hVecs, ldhVecs, restartPerm, hVals,
               flags, iev, ievSize, blockNorms, evecs, evals, resNorms,
               evecsHat, ldevecsHat, M, ldM, numConverged, numConvergedStored,
               *numPrevRetained, &indexOfPreviousVecs, hVecsPerm, *reset,
               machEps, rwork, rworkSize, iwork0, iworkSize0, primme), -1);
   }
   else {
      SCALAR *X, *Res;
      CHKERR(restart_locking_Sprimme(&restartSize, V, W, nLocal, basisSize,
               ldV, &X, &Res, hVecs, ldhVecs, restartPerm, hVals, flags, iev,
               ievSize, blockNorms, evecs, ldevecs, evals, numConverged,
               numLocked, resNorms, lockedFlags, evecsPerm, *numPrevRetained,
               &indexOfPreviousVecs, hVecsPerm, *reset, machEps, rwork,
               rworkSize, iwork0, iworkSize0, primme), -1);
   }

   *reset = 0;

   /* Rearrange prevRitzVals according to restartPerm */

   if (primme->target != primme_smallest && primme->target != primme_largest) {
      permute_vecs_Rprimme(prevRitzVals, 1, basisSize, 1, restartPerm,
            (REAL*)rwork, iwork0);
      for (i=0; i<restartSize; i++) {
         if (restartPerm[i] >= *numPrevRitzVals) {
            prevRitzVals[i] = hVals[i];
         }
      }
      permute_vecs_Rprimme(prevRitzVals, 1, restartSize, 1, hVecsPerm,
            (REAL*)rwork, iwork0);
      *numPrevRitzVals = restartSize;
   }

   if (newldhVecs == 0) newldhVecs = restartSize;
   if (newldhU == 0) newldhU = restartSize;
   CHKERR(restart_projection_Sprimme(V, ldV, W, ldV, H, ldH, VtBV, ldVtBV, Q,
            ldQ, nLocal, R, ldR, QtV, ldQtV, hU, ldhU, newldhU,
            indexOfPreviousVecsBeforeRestart, hVecs, ldhVecs, newldhVecs, hVals,
            hSVals, restartPerm, hVecsPerm, restartSize, basisSize,
            *numPrevRetained, indexOfPreviousVecs, evecs, numConvergedStored,
            primme->nLocal, evecsHat, ldevecsHat, M, ldM, UDU, ldUDU, ipivot,
            targetShiftIndex, *numConverged, numArbitraryVecs, hVecsRot,
            ldhVecsRot, rworkSize, rwork, iworkSize0, iwork0, machEps, primme),
         -1);

   /* If all request eigenpairs converged, force the converged vectors at the */
   /* beginning of V                                                          */

   if (*numConverged >= primme->numEvals && !primme->locking) {
      permute_vecs_Sprimme(V, nLocal, restartSize, ldV, hVecsPerm, rwork,
            iwork0);
   }

   *restartSizeOutput = restartSize; 

   /* If VtBV is computed, estimate the orthogonality error as                */
   /* ||I - V'*B*V||_F * ||A||. Use this value to estimate the residual error */
   /* and to bound the minimum residual error norm the solver can converge    */

   if (VtBV) {
      REAL n = 0.0;
      int i,j;
      for (i=0; i<restartSize; i++) {
         for (j=0; j<i; j++) n += 2*REAL_PART(CONJ(VtBV[i*ldVtBV+j])*VtBV[i*ldVtBV+j]);
         n += REAL_PART(CONJ(VtBV[i*ldVtBV+i] - (SCALAR)1.0)*(VtBV[i*ldVtBV+i] - (SCALAR)1.0));
      }
      if (*restartsSinceReset <= 1) {
         primme->stats.maxConvTol = max(primme->stats.maxConvTol,
               sqrt(n)*primme->stats.estimateLargestSVal);
      }
      primme->stats.estimateResidualError =
         sqrt((double)*restartsSinceReset) * sqrt(n) * aNorm;
   }
   else {
      primme->stats.estimateResidualError =
         2 * sqrt((double)*restartsSinceReset) * machEps * aNorm;
   }


   return 0;
}

/*******************************************************************************
 * Subroutine: restart_soft_locking - This routine replaces V with V*c, some
 *             subset of hVecs. Also it may include components from vectors 
 *             from the previous iteration.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 *
 * nLocal           Number of rows of V assigned to the node
 *
 * ldV              The leading dimension of V and W
 *
 * hR               The factors of the QR decomposition of (A - targetShift*B)*V
 *
 * ldhR             The leading dimension of Q and R
 *
 * hU               The eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of hU
 *
 * basisSize        Size of the basis V
 *
 * resNorms         The residual norms of the converged eigenpairs
 *
 * previousHVecs    Coefficient vectors retained from the previous iteration
 *
 * rwork            Real work array
 *
 * rworkSize        Must be of size 
 *
 * iwork            Integer work array
 *                  
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSize      The number of vectors to restart with.
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * W                A*V
 *
 * X                Reference to the Ritz vectors of the eigenpairs in the block
 *
 * R                Reference to the residual vectors of the eigenpairs in the block
 *
 * hVecs, ldhVecs   The eigenvectors of H and the leading dimension of H
 *
 * restartPerm      The permutation applied to the columns of hVecs before restarting
 *
 * hVals            The eigenvalues of H
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * iev              Array of size blockSize indicating which Ritz vectors are
 *                  targeted in the block
 *
 * ievSize          The length of iev
 *
 * blockNorms       The residual norms of the eigenpairs in the block
 *
 * evecs            The converged Ritz vectors. Without locking, all converged
 *                  eigenvectors are copied from V to evecs if skew projections
 *                  are required
 *
 * evals            The converged Ritz values
 *
 * resNorms         The residual norms of the converged eigenpairs
 *
 * evecsHat         K^{-1}evecs
 *
 * ldevecsHat       The leading dimension of evecsHat
 *
 * M, ldM           evecs'*evecsHat and the leading dimension of M
 *
 * numConverged     The number of converged eigenpairs
 *
 * numConvergedStored The # of converged vectors copied to evecs
 *
 * numPrevRetained  As input the number of columns of previousHVecs. As output the
 *                  number of columns added to V
 *
 * indexOfPreviousVecs The first column in the output V that has a vector from previousHVecs
 *
 * hVecsPerm        The permutation that orders the output hVals and hVecs as primme.target
 *
 * 
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * reset            flag to reset V and W at this restart
 * 
 *
 * Return value
 * ------------
 * int   > 0 the restart size   
 *        -2 restart_H failed
 *        -4 factorization of M failed
 *        -5 flags do not correspond to converged pairs in pseudolocking
 *       
 ******************************************************************************/
 
static int restart_soft_locking_Sprimme(int *restartSize, SCALAR *V,
       SCALAR *W, PRIMME_INT nLocal, int basisSize, PRIMME_INT ldV, SCALAR **X,
       SCALAR **R, SCALAR *hVecs, int ldhVecs, int *restartPerm,
       REAL *hVals, int *flags, int *iev, int *ievSize, REAL *blockNorms,
       SCALAR *evecs, REAL *evals, REAL *resNorms, SCALAR *evecsHat,
       PRIMME_INT ldevecsHat, SCALAR *M, int ldM, int *numConverged,
       int *numConvergedStored, int numPrevRetained, int *indexOfPreviousVecs,
       int *hVecsPerm, int reset, double machEps, SCALAR *rwork,
       size_t *rworkSize, int *iwork, int iworkSize, primme_params *primme) {

   int i, j, k;               /* loop indices */
   int wholeSpace=0;          /* if all pairs in V are marked as converged */
   double aNorm;

   /* Return memory requirement */

   if (V == NULL) {
      SCALAR t;
      REAL d;
      *rworkSize = max(*rworkSize, (size_t)basisSize); /* permute_vecs for hVecs */
      CHKERR(Num_reset_update_VWXR_Sprimme(NULL, NULL, nLocal, basisSize, 0, &t,
            *restartSize, 0, NULL,
            &t, 0, *restartSize, 0,
            &t, *numConverged, *numConverged+*ievSize, 0,
            NULL, 0, 0, 0, 0,
            &t, 0, *restartSize, 0,
            &t, *numConverged, *numConverged+*ievSize, 0, &d,
            NULL, 0, 0,
            0, 0.0, NULL, rworkSize, primme), -1);
      /* if evecsHat, permutation matrix & compute_submatrix workspace */
      if (evecsHat) {
         *rworkSize = max(*rworkSize, 
               (size_t)(primme->numOrthoConst+*numConverged)*
               (size_t)(primme->numOrthoConst+*numConverged)*2);
      }
      /* int workspace for permute_vecs and checking_convergence */
      *iwork = max(*iwork, basisSize);
      return 0;
   }

   aNorm = max(primme->stats.estimateLargestSVal, primme->aNorm);

   /* -------------------------------------------------------------------------- */ 
   /* Check if any of the previous flagged converged eigenvalues seems           */
   /* to have become unconverged by checking hVals[i]-evals[i] < tol.            */
   /* If it fails, flag it UNCONVERGED and let it be targeted again. This avoids */  
   /* early converged but unwanted evs preventing wanted from being targeted.    */
   /* Update maxConvTol for the remaining converged pairs.                       */
   /* Update the number of converged values also.                                */
   /* -------------------------------------------------------------------------- */

   if (basisSize + primme->numOrthoConst < primme->n) {
      primme->stats.maxConvTol = 0.0;
      *numConverged = 0;
      for (i=0; i<primme->numEvals; i++) {
         if (flags[i] != UNCONVERGED && fabs(hVals[i]-evals[i]) > resNorms[i]) {
            flags[i] = UNCONVERGED;
         }
         else if (flags[i] != UNCONVERGED) {
            primme->stats.maxConvTol = max(primme->stats.maxConvTol, resNorms[i]);
            if (i < primme->numEvals) (*numConverged)++;
         }
      }
   }

   /* -------------------------------------------------------------- */
   /* Indicate that all pairs are marked as converged                */
   /* -------------------------------------------------------------- */

   else {
      wholeSpace = 1;
   }

   /* -------------------------------------------------------------- */
   /* Restart V and W by replacing it with the current Ritz vectors. */
   /* Result of restartPerm (modified part indicated with ---)       */
   /*                                                                */
   /*      converged | non-conv | prevRitzVecs | X & R               */
   /* V: [-----------|----------|--------------|- X ---|    )        */
   /* W: [-----------|----------|--------------|- R ---|    )        */
   /*                           ^ indexOfPreviousVecs                */
   /*    [-----------) numConverged            ^ restartSize         */
   /*                           [--------------) numPrevRetained     */
   /*                                  ievSize [-------)             */
   /*                                                                */
   /* X & R have the eigenvectors and residual vectors of the        */
   /* first ievSize candidates pairs to be targeted after restart.   */
   /* Their computation is performed more efficiently here together  */
   /* with the V, W                                                  */
   /* -------------------------------------------------------------- */

   *indexOfPreviousVecs = *restartSize;

   *restartSize += numPrevRetained;

   *ievSize = max(0, min(min(min(
                  primme->maxBlockSize,
                  primme->numEvals-*numConverged+1),
                  primme->maxBasisSize-*restartSize-numPrevRetained),
                  basisSize-*numConverged));
   *ievSize = max(0, min(*ievSize, primme->minRestartSize - *numConverged));

   /* Generate restartPerm */

   for (i=j=k=0; i<basisSize; i++) {
      if (k >= *numConverged || flags[i] == UNCONVERGED) {
         restartPerm[*numConverged + j++] = i;
      }
      else {
         restartPerm[k++] = i;
      }
   }
   assert(j+k == basisSize);
 
   /* Permute hVals and hVecs */

   permute_vecs_Rprimme(hVals, 1, basisSize, 1, restartPerm, (REAL*)rwork, iwork);
   permute_vecs_Sprimme(hVecs, basisSize, basisSize, ldhVecs, restartPerm, rwork,
         iwork);

   /* -------------------------------------------------------------- */
   /* Restart V and W by replacing it with the current Ritz vectors. */
   /* Compute X, R, blockNorms for the next values in the block.     */
   /* -------------------------------------------------------------- */

   *X = &V[*restartSize*ldV];
   *R = &W[*restartSize*ldV];

   CHKERR(Num_reset_update_VWXR_Sprimme(V, W, nLocal, basisSize, ldV,
            hVecs, *restartSize, ldhVecs, hVals,
            V, 0, *restartSize, ldV,
            *X, *numConverged, *numConverged+*ievSize, ldV,
            NULL, 0, 0, 0, 0,
            W, 0, *restartSize, ldV,
            *R, *numConverged, *numConverged+*ievSize, ldV, blockNorms,
            NULL, 0, 0,
            reset, machEps, rwork, rworkSize, primme), -1);

   if (!wholeSpace) {
      /* ----------------------------------------------------------------- */
      /* Generate the permutation hVecsPerm that undoes restartPerm        */
      /* ----------------------------------------------------------------- */

      for (i=0; i<basisSize; i++)
         hVecsPerm[restartPerm[i]] = i;
   }

   else {
      REAL *fakeResNorms = (REAL*)rwork;
      size_t rworkSize0 = *rworkSize;

      for (i=0; i<*restartSize; i++)
         fakeResNorms[i] = aNorm*machEps;
      primme->stats.estimateResidualError = 0;
      assert(rworkSize0 >= (size_t)*restartSize);
      rworkSize0 -= (size_t)*restartSize;
      CHKERR(check_convergence_Sprimme(V, nLocal, ldV, NULL, 0, NULL, 0, 0,
               0, *restartSize, flags, fakeResNorms, hVals, NULL, machEps,
               rwork+*restartSize, &rworkSize0, iwork, iworkSize, primme), -1);

      *numConverged = 0;
      for (i=0; i<*restartSize && *numConverged<primme->numEvals; i++) {
         if (flags[i] != UNCONVERGED) {
            (*numConverged)++;
         }
      }
      for (i=j=k=0; i<*restartSize; i++) {
         if (flags[i] != UNCONVERGED && k < *numConverged) {
            hVecsPerm[k++] = i;
         }
         else {
            hVecsPerm[(*numConverged) + j++] = i;
         }
      }

      assert(iworkSize >= *restartSize);
      permute_vecs_Rprimme(hVals, 1, *restartSize, 1, hVecsPerm, (REAL*)rwork,
            iwork);
      permute_vecs_Sprimme(hVecs, basisSize, *restartSize, ldhVecs,
            hVecsPerm, rwork, iwork);
      permute_vecs_iprimme(restartPerm, *restartSize, hVecsPerm, iwork);
      for (i=0; i<*restartSize; i++)
         hVecsPerm[i] = i;
   }

   /* Compute iev */

   for (i=0; i<*ievSize; i++)
      for (j=0; j<*restartSize; j++)
         if (hVecsPerm[j] == *numConverged+i)
            iev[i] = j;

   /* --------------------------------------------------------------------- */
   /* If the user requires (I-QQ') projectors in JDQMR without locking,     */
   /* the converged eigenvectors are copied temporarily to evecs. There     */
   /* they stay locked  for use in (I-QQ') and (I-K^{-1}Q () Q') projectors.*/
   /* The Ritz vectors remain in the basis, and they will overwrite evecs   */
   /* the end.                                                              */
   /* We recommend against this type of usage. It's better to use locking.  */
   /* --------------------------------------------------------------------- */

   /* Andreas NOTE: is done inefficiently for the moment. We should only */
   /* add the recently converged. But we need to differentiate them      */
   /* from flags...                                                      */
   /* Eloy NOTE: that will not be enough, the Ritz vectors change from   */
   /* one restart to the next even if they weren't targeted. Probably    */
   /* the best we can do is to always copy all converged vectors.        */

   if (evecsHat) {
      int newNumConvergedStored=0, oldSizeM, newSizeM;
      size_t rworkSize0;

      /* Pack evecs and evecsHat for the converged pairs restartPerm[0:numConverged] */

      for (i=0; i < *numConverged && restartPerm[i] < *numConvergedStored; i++) {
         Num_copy_matrix_Sprimme(&evecs[(restartPerm[i]+primme->numOrthoConst)*nLocal],
               nLocal, 1, nLocal,
               &evecs[(newNumConvergedStored+primme->numOrthoConst)*nLocal],
               nLocal);
         Num_copy_matrix_Sprimme(&evecsHat[(restartPerm[i]+primme->numOrthoConst)*ldevecsHat],
               nLocal, 1, ldevecsHat,
               &evecsHat[(newNumConvergedStored+primme->numOrthoConst)*ldevecsHat],
               ldevecsHat);
         newNumConvergedStored++;
      }

      /* Apply restartPerm to rows and columns of M */

      oldSizeM = *numConvergedStored + primme->numOrthoConst;
      newSizeM = newNumConvergedStored + primme->numOrthoConst;
      for (i=0; i<oldSizeM*newSizeM; i++)
         rwork[i] = 0.0;
      for (i=0; i < primme->numOrthoConst; i++)
         rwork[oldSizeM*i + i] = 1.0;
      for (; i < newSizeM; i++)
         rwork[oldSizeM*i + restartPerm[i]+primme->numOrthoConst] = 1.0;
      assert(*rworkSize >= (size_t)oldSizeM*newSizeM);
      rworkSize0 = *rworkSize - (size_t)oldSizeM*newSizeM;
      compute_submatrix_Sprimme(rwork, newSizeM, oldSizeM, M, oldSizeM, ldM,
         M, ldM, rwork+oldSizeM*newSizeM, &rworkSize0);

      *numConvergedStored = newNumConvergedStored;
   }

   return 0;
}

/*******************************************************************************
 * Subroutine: restart_locking - This routine is only called when locking and
 *    replaces V (W) with V*c (W*c), where c is some subset of hVecs 
 *    and prevRitzVecs.
 *
 *    This subroutine locks converged Ritz pairs. The converged Ritz vectors
 *    are copied to evecs, and the converged Ritz values are copied to evals.
 *    The evals array is maintained in sorted order; however the evecs array
 *    is not.  Instead, a permutation array is maintained so that the locked
 *    Ritz vectors may be sorted upon return to the user.  
 *
 *    Vectors that have remained converged are locked to the evecs array and
 *    may be replaced afterward by initial guesses if there are any remaining.
 *
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal           Number of rows of V assigned to the node
 *
 * ldV              The leading dimension of V and W
 *
 * ldhR             The leading dimension of Q and R
 *
 * hU               The eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of hU
 *
 * basisSize        Size of the basis V
 *
 * numGuesses       Number of remaining initial guesses
 *
 * previousHVecs    Coefficient vectors retained from the previous iteration
 *
 * ldpreviousHVecs  The leading dimension of previousHVecs
 *
 * rwork            Real work array
 *
 * rworkSize        Must be of size 
 *                  
 * iwork            Integer work array
 *                  
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSize      The number of vectors to restart with
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * W                A*V
 *
 * X                Reference to the Ritz vectors of the eigenpairs in the block
 *
 * R                Reference to the residual vectors of the eigenpairs in the block
 *
 * hVecs, ldhVecs   The eigenvectors of H and the leading dimension of H
 *
 * restartPerm      The permutation applied to the columns of hVecs before restarting
 *
 * hVals            The eigenvalues of H
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * iev              Array of size blockSize indicating which Ritz vectors are
 *                  targeted in the block
 *
 * ievSize          The length of iev
 *
 * blockNorms       The residual norms of the eigenpairs in the block
 *
 * evecs            The converged Ritz vectors
 *
 * numConverged     The number of converged eigenpairs
 *
 * numLocked        The number of eigenpairs to be locked
 *
 * resNorms         The residual norms of the converged eigenpairs
 *
 * evecsperm        The permutation that orders the converged pairs as primme.target
 *
 * lockedFlags      The flags of the locked pairs
 *
 * numPrevRetained  As input the number of columns of previousHVecs. As output the
 *                  number of columns added to V
 *
 * indexOfPreviousVecs The first column in the output V that has a vector from previousHVecs
 *
 * hVecsPerm        The permutation that orders the output hVals and hVecs as primme.target
 *
 * numArbitraryVecs On input, the number of leading coefficient vectors in
 *                  hVecs that do not come from solving the projected problem.
 *                  On output, the number of such vectors that are in the
 *                  restarted basis
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * reset            flag to reset V and W in the next restart
 *
 *
 * Return value
 * ------------
 * int   > 0 the restart size   
 *        -2 restart_H failed
 *        -4 factorization of M failed
 *        -5 flags do not correspond to converged pairs in pseudolocking
 *       
 ******************************************************************************/
 
static int restart_locking_Sprimme(int *restartSize, SCALAR *V, SCALAR *W, 
      PRIMME_INT nLocal, int basisSize, PRIMME_INT ldV, SCALAR **X, SCALAR **R,
      SCALAR *hVecs, int ldhVecs, int *restartPerm, REAL *hVals, int *flags,
      int *iev, int *ievSize, REAL *blockNorms, SCALAR *evecs,
      PRIMME_INT ldevecs, REAL *evals, int *numConverged, int *numLocked,
      REAL *resNorms, int *lockedFlags, int *evecsperm, int numPrevRetained,
      int *indexOfPreviousVecs, int *hVecsPerm, int reset, double machEps,
      SCALAR *rwork, size_t *rworkSize, int *iwork, int iworkSize,
      primme_params *primme) {

   int i, j, k;             /* Loop variables                                 */
   int numPacked;           /* The number of coefficient vectors moved to the */
                            /* end of the hVecs array.                        */
   int maxBlockSize;        /* max block size for the next iteration          */
   int sizeBlockNorms;      /* Size of the block                              */
   int failed;              /* Number of vectors to be locked that didn't pass the conv. test */
   int *ifailed;            /* Indices of vectors to be locked vectors that failed */
   int left;                /* Index of the first column with a vector to be locked */
   int overbooking;         /* If basis has more converged pairs than numEvals */
   REAL *lockedResNorms;  /* residual norms for locked vectors              */
   int numLocked0 = *numLocked; /* aux variables                         */
   REAL *blockNorms0;
   size_t rworkSize0 = *rworkSize;

   /* Return memory requirement */
   if (V == NULL) {
      SCALAR t;
      REAL d;
      /* for permute_vecs and fakeResNorms */
      *rworkSize = max(*rworkSize, (size_t)basisSize*2);
      *rworkSize = max(*rworkSize, (size_t)
            compute_residual_columns(nLocal, NULL, NULL,
               basisSize, NULL, 0, NULL, 0, NULL, primme->maxBlockSize, 0, 0,
               NULL, 0, NULL, primme->maxBlockSize, NULL, 0, NULL, 0, NULL, 0));
      CHKERR(Num_reset_update_VWXR_Sprimme(NULL, NULL, nLocal, basisSize,
               0, NULL, *restartSize, 0, NULL,
               &t, 0, *restartSize+*numLocked, 0,
               &t, 0, *ievSize, 0,
               &t, 0, *restartSize, *restartSize+*numLocked, 0,
               &t, 0, *restartSize, 0,
               &t, 0, *ievSize, 0, &d,
               &d, *restartSize, *numLocked,
               0, 0.0, NULL, rworkSize, primme), -1);
      CHKERR(check_convergence_Sprimme(NULL, nLocal, 0, NULL, 0, NULL,
               *numLocked, 0, *restartSize, *restartSize+*numLocked, NULL, NULL,
               NULL, NULL, 0.0, NULL, rworkSize, iwork, 0, primme), -1);
      /* for permute_vecs and ifailed */
      *iwork = max(*iwork, 2*basisSize);
      return 0;
   }

   overbooking = (*numConverged > primme->numEvals);

   /* -------------------------------------------------------------- */
   /* Rearrange hVals and hVecs so that V*hVecs and W*hVecs will     */
   /* look like this (modified part indicated with ---):             */
   /*                                                                */
   /*      restarted  | prevRitzVecs | to be locked | X & R          */
   /* V: [------------|--------------|--------------|- X -|    )     */
   /* W: [------------|--------------|--------------|- R -|    )     */
   /*                                ^ left         ^ restartSize    */
   /*                 ^ indexOfPreviousVecs                          */
   /*          [--------------) numPrevRetained                      */
   /*                                [--------------) numPacked      */
   /*                         [----) sizeBlockNorms [-----)          */
   /*                                maxBlockSize   [-------)        */
   /*                                                                */
   /*                           to be locked                         */
   /* evecs: [----------------|--------------)                       */
   /*                         ^ numLocked                            */
   /*                                                                */
   /* X & R have the eigenvectors and residual vectors of the        */
   /* first sizeBlockNorms candidates pairs to be targeted after     */
   /* restart. Their computation is performed more efficiently here  */
   /* together with the V, W                                         */
   /* -------------------------------------------------------------- */

   maxBlockSize = max(0, min(min(
         *restartSize,
         primme->maxBlockSize),
         primme->numEvals-*numConverged+1));

   sizeBlockNorms = max(0, min(
         maxBlockSize,
         primme->maxBasisSize-*restartSize-numPrevRetained
               -*numConverged+*numLocked));

   *indexOfPreviousVecs = *restartSize;

   left = *restartSize + numPrevRetained;

   *restartSize += numPrevRetained + *numConverged - *numLocked;

   for (i=k=numPacked=0; i<basisSize; i++) {
      if (flags[i] != UNCONVERGED
            /* Try to converge all pairs in case of overbooking               */
            && (overbooking || (
            /* Otherwise don't check more than numEvals-numLocked             */
                  numPacked < primme->numEvals-*numLocked
            /* Check only the first pairs, except if finding closest_leq/geq  */
            /* because with refined extraction the pairs may not be ordered   */
            /* by this criterion.                                             */
                  && (i < primme->numEvals-*numLocked
                     || primme->target == primme_closest_geq
                     || primme->target == primme_closest_leq)))) {
         restartPerm[left + numPacked++] = i;
      }
      else if (k < left) {
         restartPerm[k++] = i;
      }
      else if (!overbooking) {
         restartPerm[min(*numConverged, primme->numEvals)-*numLocked + k++] = i;
      }
      else {
         restartPerm[*numConverged-*numLocked + k++] = i;
      }
   }

   assert(iworkSize >= basisSize);
   permute_vecs_Rprimme(hVals, 1, basisSize, 1, restartPerm, (REAL*)rwork, iwork);
   permute_vecs_Sprimme(hVecs, basisSize, basisSize, ldhVecs, restartPerm, rwork, iwork);

   /* -------------------------------------------------------------- */
   /* Restart V and W by replacing it V*hVecs and W*hVecs, copy      */
   /* vectors to be locked into evecs and compute X and R.           */
   /* -------------------------------------------------------------- */

   if (!overbooking) {
      lockedResNorms = &resNorms[*numLocked];
   }
   else {
      lockedResNorms = (REAL*)rwork;
      rwork += basisSize;
      assert(rworkSize0 >= (size_t)basisSize);
      rworkSize0 -= basisSize;
   }
   *X = &V[*restartSize*ldV];
   *R = &W[*restartSize*ldV];
   CHKERR(Num_reset_update_VWXR_Sprimme(V, W, nLocal, basisSize, ldV,
            hVecs, *restartSize, ldhVecs, hVals,
            V, 0, *restartSize, ldV,
            *X, 0, sizeBlockNorms, ldV,
            evecs, *numLocked+primme->numOrthoConst, left,
                !overbooking?left+numPacked:left, ldevecs,
            W, 0, *restartSize, ldV,
            *R, 0, sizeBlockNorms, ldV, blockNorms,
            lockedResNorms, left, *restartSize,
            reset, machEps, rwork, &rworkSize0, primme), -1);
 
   /* -------------------------------------------------------------- */
   /* Recompute flags for the vectors to be locked.                  */
   /* NOTE: the eigenvals flagged as practically converged will keep */
   /*       it if their residual norm is still less than             */
   /*       sqrt(numLocked)*tol                                      */
   /* -------------------------------------------------------------- */

   permute_vecs_iprimme(flags, basisSize, restartPerm, iwork);
   CHKERR(check_convergence_Sprimme(&V[ldV*left],
            nLocal, ldV, NULL, 0, NULL, *numLocked, 0, left,
            left+numPacked, flags, lockedResNorms, hVals, NULL,
            machEps, rwork, &rworkSize0, iwork, iworkSize, primme), -1);

   /* -------------------------------------------------------------- */
   /* Copy the values for the converged values into evals, and in    */
   /* case of overbooking also the converged vectors into evecs.     */
   /* -------------------------------------------------------------- */

   for (i=left, j=0; i < left+numPacked; i++) {
      if (flags[i] != UNCONVERGED && *numLocked+j < primme->numEvals) {
         if (overbooking) {
            Num_copy_matrix_Sprimme(&V[i*ldV], nLocal, 1, ldV,
                  &evecs[(*numLocked+primme->numOrthoConst+j)*ldevecs],
                  ldevecs);
         }
         evals[*numLocked+j++] = hVals[i];
      }
      else {
         flags[i] = UNCONVERGED;
      }
   }

   /* -------------------------------------------------------------- */
   /* Merge back the pairs that failed to be locked into the         */
   /* restarted basis.                                               */
   /* -------------------------------------------------------------- */

   /* Generate the permutation that separates the pairs that failed  */
   /* to be locked to the ones that passed.                          */

   assert(iworkSize >= numPacked);
   ifailed = iwork;
   for (i=left, failed=0; i < left+numPacked; i++)
      if (flags[i] == UNCONVERGED) ifailed[failed++] = i-left;
   for (i=left, j=0; i < left+numPacked; i++)
      if (flags[i] != UNCONVERGED) ifailed[failed+j++] = i-left;

   if (1 /* Put zero to disable the new feature */) {
      /* New feature: the pairs that failed to be locked are         */
      /* rearranged with the rest of the restarted vectors and they  */
      /* are considered to be included in the block.                 */

      /* Generate hVecsPerm merging candidate vectors and vectors to */
      /* be locked that failed to pass the convergence test with the */
      /* restarted pairs. The order is determined by the original    */
      /* position of the pairs.                                      */

      maxBlockSize = max(0, min(
            maxBlockSize,
            primme->maxBasisSize-*restartSize-numPrevRetained
                  -*numConverged+*numLocked));

      blockNorms0 = (REAL*)rwork;
      for (i=0; i<sizeBlockNorms; i++) blockNorms0[i] = blockNorms[i];
      for (i=j=k=0; i<*indexOfPreviousVecs || j<failed; k++) {

         /* If there is no more failed vectors or the candidate      */
         /* vector was before the failed, take the candidate         */

         if (i < *indexOfPreviousVecs && (j >= failed 
                  || restartPerm[i] < restartPerm[left+ifailed[j]])) {

            if (k < maxBlockSize && i < sizeBlockNorms)
               blockNorms[k] = blockNorms0[i];
            hVecsPerm[k] = i++;

         }

         /* Otherwise take the failed */

         else {
            if (k < maxBlockSize)
               blockNorms[k] = resNorms[numLocked0+ifailed[j]];
            hVecsPerm[k] = left + j++;
         }
      }

      /* Generate the rest of the permutation of hVecsPerm           */

      for (i=0; i<numPrevRetained; i++) hVecsPerm[k++] = i+*indexOfPreviousVecs;
      assert(k == *indexOfPreviousVecs + numPrevRetained + failed);
      for (; k < basisSize; k++) hVecsPerm[k] = -1;

      /* -------------------------------------------------------------- */
      /* Update X and R with the failed vectors and place the failed    */
      /* vectors after the candidate vectors (modified part indicated   */
      /* with ---):                                                     */
      /*                                                                */
      /*      restarted  | prevRitzVecs | failed locked | X & R         */
      /* V: [            |              |---------------|- X ---|    )  */
      /* W: [            |              |---------------|- R ---|    )  */
      /*                                ^ left                          */
      /*                         failed [---------------)               */
      /*                                 maxBlockSize   [---------)     */
      /*                                                                */
      /* -------------------------------------------------------------- */

      compute_residual_columns(nLocal, &hVals[left],
            &V[left*ldV], failed, ifailed, ldV, &W[left*ldV], ldV, *X,
            sizeBlockNorms, ldV, 0, *R, ldV,
            &V[(left+failed)*ldV], maxBlockSize, hVecsPerm, ldV,
            &W[(left+failed)*ldV], ldV, rwork, TO_INT(rworkSize0));
   }
   else {
      /* The failed pairs are not rearranged with the rest of           */
      /* non-converged pairs and they may not be in the block in the    */
      /* next iteration. This was the behaviour of locking in previous  */
      /* versions than 2.0.                                             */

      for (i=0; i<basisSize; i++) hVecsPerm[i] = i;

      /* Pack together the failed vectors and place them after          */
      /* candidates.                                                    */

      Num_compact_vecs_Sprimme(&V[left*ldV], nLocal, failed, ldV,
            ifailed, &V[left*ldV], ldV, 0);
      Num_compact_vecs_Sprimme(&W[left*ldV], nLocal, failed, ldV,
            ifailed, &W[left*ldV], ldV, 0);

      /* Copy X and R after the failed vectors */

      Num_copy_matrix_Sprimme(*X, nLocal, sizeBlockNorms, ldV,
            &V[(left+failed)*ldV], ldV);
      Num_copy_matrix_Sprimme(*R, nLocal, sizeBlockNorms, ldV,
            &W[(left+failed)*ldV], ldV);
   }

   /* Modify hVals, hVecs and restartPerm to add the pairs failed to */
   /* be locked just after the candidates                            */

   Num_compact_vecs_Sprimme(&hVecs[left*ldhVecs], basisSize, failed,
         ldhVecs, ifailed, &hVecs[left*ldhVecs], ldhVecs, 0);
   Num_compact_vecs_Rprimme(&hVals[left], 1, failed, 1, ifailed, &hVals[left],
         1, 0);
   assert(iworkSize >= numPacked*2);
   permute_vecs_iprimme(&restartPerm[left], numPacked, ifailed,
         ifailed+numPacked);

   *X = &V[(left+failed)*ldV];
   *R = &W[(left+failed)*ldV];

   /* Pack those Ritz vectors that are actually converged from */
   /* those that have been copied into evecs. Also insertSort  */
   /* the converged Ritz value within the evals array.         */ 

   for (i=left; i < left+numPacked; i++) {
       if (flags[i] != UNCONVERGED && *numLocked < primme->numEvals) {
         REAL resNorm = resNorms[*numLocked] = lockedResNorms[i-left];
         REAL eval = evals[*numLocked];
         if (!overbooking) {
            Num_copy_matrix_Sprimme(
                  &evecs[(numLocked0+i-left+primme->numOrthoConst)*ldevecs],
                  nLocal, 1, ldevecs,
                  &evecs[(*numLocked+primme->numOrthoConst)*ldevecs], ldevecs);
         }

         (*numLocked)++;

         /* Report a pair was hard locked */
         /* NOTE: do this before sorting evals */
         if (primme->monitorFun) {
            primme_event EVENT_LOCKED = primme_event_locked;
            int err;
            lockedFlags[*numLocked-1] = flags[i];
            primme->stats.elapsedTime = primme_wTimer(0);
            CHKERRM((primme->monitorFun(NULL, NULL, NULL, NULL, NULL, NULL,
                        NULL, evals, numLocked, lockedFlags, resNorms, NULL, NULL,
                        &EVENT_LOCKED, primme, &err), err), -1,
                  "Error returned by monitorFun: %d", err);
         }

         insertionSort(eval, evals, resNorm, resNorms, flags[i], lockedFlags,
               evecsperm, *numLocked-1, primme);

         /* Update maxConvTol if it wasn't practically converged */
         if (flags[i] == CONVERGED) {
            primme->stats.maxConvTol = max(primme->stats.maxConvTol, resNorm);
         }
      }
   }

   *restartSize = left + failed;
   *ievSize = min(maxBlockSize, sizeBlockNorms+failed);
   *numConverged = *numLocked;
   for (i=0; i<*ievSize; i++) iev[i] = i;

   /* ----------------------------------------------------------------- */
   /* There is no converged pair in the restarted basis                 */
   /* ----------------------------------------------------------------- */
   
   for (i=0; i<basisSize; i++) flags[i] = UNCONVERGED;

   return 0;
}


/******************************************************************************
 * Function Num_reset_update_VWXR - This subroutine performs the next operations:
 *
 *    X0 = V*h(nX0b+1:nX0e), X1 = V*h(nX1b+1:nX1e), evecs(evecsSize:) = V*h(nX2b+1:nX2e)
 *    if reset, then
 *       reortho evecs(evecsSize:)
 *       reortho X0
 *    end
 *    Wo = A*X0(nWob+1-nX0b:nWoe-nX0b),
 *    R = Wo(nRb+1-nWob:nRe-nWob) - X0(nRb+1-nX0b:nRe-nX0b)*diag(hVals(nRb+1:nRe)),
 *    Rnorms = norms(R),
 *    rnorms = norms(Wo(nrb+1-nWob:nre-nWob) - X0(nrb+1-nX0b:nre-nX0b)*diag(hVals(nrb+1:nre))),
 *
 * NOTE: if Rnorms and rnorms are requested, nRb-nRe+nrb-nre < mV
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V, W        input basis
 * mV,nV,ldV   number of rows and columns and leading dimension of V and W
 * h           input rotation matrix
 * nh          Number of columns of h
 * ldh         The leading dimension of h
 * hVals       Array of values
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * X0          Output matrix V*h(nX0b:nX0e-1) (optional)
 * nX0b, nX0e  Range of columns of h
 * X1          Output matrix V*h(nX1b:nX1e-1) (optional)
 * nX1b, nX1e  Range of columns of h
 * evecs       Output matrix V*h(nX2b:nX2e-1) (optional)
 * evecsSize   First column where to copy V*h(nX2b:nX2e-1) (optional)
 * nX2b, nX2e  Range of columns of h
 * Wo          Output matrix W*h(nWob:nWoe-1) (optional)
 * nWob, nWoe  Range of columns of h
 * R           Output matrix (optional)
 * nRb, nRe    Range of columns of h and hVals
 * Rnorms      Output array with the norms of R (optional)
 * rnorms      Output array with the extra residual vector norms (optional)
 * nrb, nre    Columns of residual vector to compute the norm
 * reset       if reset>1, reothogonalize Xi; if reset>0, recompute Wo=A*X0
 * 
 * NOTE: n*e, n*b are zero-base indices of ranges where the first value is
 *       included and the last isn't.
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_reset_update_VWXR_Sprimme(SCALAR *V, SCALAR *W, PRIMME_INT mV,
   int nV, PRIMME_INT ldV,
   SCALAR *h, int nh, int ldh, REAL *hVals,
   SCALAR *X0, int nX0b, int nX0e, PRIMME_INT ldX0,
   SCALAR *X1, int nX1b, int nX1e, PRIMME_INT ldX1,
   SCALAR *evecs, int evecsSize, int nX2b, int nX2e, PRIMME_INT ldevecs,
   SCALAR *Wo, int nWob, int nWoe, PRIMME_INT ldWo,
   SCALAR *R, int nRb, int nRe, PRIMME_INT ldR, REAL *Rnorms,
   REAL *rnorms, int nrb, int nre,
   int reset, double machEps, SCALAR *rwork, size_t *lrwork,
   primme_params *primme) {

   int i, j;         /* Loop variables */
   REAL *tmp, *tmp0;

   /* Return memory requirements */
   if (V == NULL) {
      *lrwork = max(*lrwork,
            (size_t)Num_update_VWXR_Sprimme(
               V, W, mV, nV, ldV, h, nh, ldh, hVals,
               X0, nX0b, nX0e, ldX0,
               X1, nX1b, nX1e, ldX1,
               evecs, nX2b, nX2e, ldevecs,
               Wo, nWob, nWoe, ldWo,
               R, nRb, nRe, ldR, Rnorms,
               rnorms, nrb, nre,
               NULL, 0, primme));
      return 0;
   }

   /* Quick exit */
   if (reset == 0) {
      CHKERR(Num_update_VWXR_Sprimme(
               V, W, mV, nV, ldV, h, nh, ldh, hVals,
               X0, nX0b, nX0e, ldX0,
               X1, nX1b, nX1e, ldX1,
               evecs?&evecs[ldevecs*evecsSize]:NULL, nX2b, nX2e, ldevecs,
               Wo, nWob, nWoe, ldWo,
               R, nRb, nRe, ldR, Rnorms,
               rnorms, nrb, nre,
               rwork, TO_INT(*lrwork), primme), -1);
      return 0;
   }


   /* R or Rnorms or rnorms imply W */
   assert(!(R || Rnorms || rnorms) || W);

   assert((size_t)(nre-nrb)*2 <= *lrwork); /* Check workspace for tmp and tmp0 */

   /* X_i = V*h(nX_ib:nX_ie-1) */

   assert(!reset || !evecs || (nX0b <= nX2b && nX2e <= nX0e));
   Num_update_VWXR_Sprimme(V, NULL, mV, nV, ldV, h, nh, ldh, NULL,
         X0, nX0b, nX0e, ldX0,
         X1, nX1b, nX1e, ldX1,
         evecs?&evecs[ldevecs*evecsSize]:NULL, nX2b, nX2e, ldevecs,
         NULL, 0, 0, 0,
         NULL, 0, 0, 0, NULL,
         NULL, 0, 0,
         rwork, TO_INT(*lrwork), primme);

   /* Reortho [evecs(evecSize:) X0] against evecs if asked */

   if (reset > 1) {
      CHKERR(ortho_Sprimme(evecs, ldevecs, NULL, 0, evecsSize, 
               evecsSize+nX2e-nX2b-1, NULL, 0, 0, mV, primme->iseed, 
               machEps, rwork, lrwork, primme), -1);
      CHKERR(ortho_Sprimme(X0, ldX0, NULL, 0, 0, nX2b-nX0b-1, evecs,
               ldevecs, evecsSize+nX2e-nX2b, mV, primme->iseed, 
               machEps, rwork, lrwork, primme), -1);
      Num_copy_matrix_Sprimme(&evecs[ldevecs*evecsSize], mV, nX2e-nX2b,
            ldevecs, &X0[ldX0*(nX2b-nX0b)], ldX0);
      CHKERR(ortho_Sprimme(X0, ldX0, NULL, 0, nX2e-nX0b, nX0e-nX0b-1,
               evecs, ldevecs, evecsSize, mV, primme->iseed, machEps, rwork,
               lrwork, primme), -1);
      assert(!X1 || (nX0b <= nX1b && nX1e <= nX0e));
      if (X1) Num_copy_matrix_Sprimme(&X0[ldX0*(nX1b-nX0b)], mV, nX1e-nX1b,
            ldX0, X1, ldX1);
   }

   /* Compute W = A*V for the orthogonalized corrections */

   assert(nWob == nX0b && nWoe == nX0e);
   CHKERR(matrixMatvec_Sprimme(X0, mV, ldX0, Wo, ldWo, 0, nWoe-nWob,
            primme), -1);
 
   /* R = Y(nRb-nYb:nRe-nYb-1) - X(nRb-nYb:nRe-nYb-1)*diag(nRb:nRe-1) */
   for (j=nRb; j<nRe; j++) {
      Num_compute_residual_Sprimme(mV, hVals[j], &X0[ldX0*(j-nX0b)],
            &Wo[ldWo*(j-nWob)], &R[ldR*(j-nRb)]);
      if (Rnorms) {
         Rnorms[j-nRb] = REAL_PART(Num_dot_Sprimme(mV, &R[ldR*(j-nRb)], 1,
                  &R[ldR*(j-nRb)], 1));
      }
   }

   /* rnorms = Y(nrb-nYb:nre-nYb-1) - X(nrb-nYb:nre-nYb-1)*diag(nrb:nre-1) */
   if (rnorms) for (j=nrb; j<nre; j++) {
      int m=min(PRIMME_BLOCK_SIZE, mV);   /* Number of rows in the cache */
      rnorms[j-nrb] = 0.0;
      for (i=0; i < mV; i+=m, m=min(m,mV-i)) {
         Num_compute_residual_Sprimme(m, hVals[j], &X0[ldX0*(j-nX0b)],
               &Wo[i+ldWo*(j-nWob)], rwork);
         rnorms[j-nrb] +=
            REAL_PART(Num_dot_Sprimme(m, rwork, 1, rwork, 1));
      }
   }

   /* Reduce Rnorms and rnorms and sqrt the results */

   if (primme->globalSumReal) {
      tmp = (REAL*)rwork;
      j = 0;
      if (Rnorms) for (i=nRb; i<nRe; i++) tmp[j++] = Rnorms[i-nRb];
      if (rnorms) for (i=nrb; i<nre; i++) tmp[j++] = rnorms[i-nrb];
      tmp0 = tmp+j;
      if (j) CHKERR(globalSum_Rprimme(tmp, tmp0, j, primme), -1);
      j = 0;
      if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = sqrt(tmp0[j++]);
      if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = sqrt(tmp0[j++]);
   }
   else {
      if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = sqrt(Rnorms[i-nRb]);
      if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = sqrt(rnorms[i-nrb]);
   }

   return 0; 
}


/*******************************************************************************
 * Subroutine: restart_projection - This routine updates Q, R, H and QtV to match the
 *    changes in V.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal           Number of rows of V assigned to the node
 *
 * basisSize        Size of the basis V
 *
 * numConverged     The number of converged eigenpairs without locking
 *
 * numLocked        The number of Ritz vectors that have been locked 
 *
 * numPrevRetained  The number of coefficient vectors from previousHVecs in hVecs
 *  
 * indexOfPreviousVecs The first column in hVecs that has a vector from previousHVecs
 *
 * indexOfPreviousVecsBeforeRestart The number of columns in hVecs that previousHVecs
 *                  where orthogonalized against
 *
 * rwork            Real work array
 *
 * rworkSize        Size of rwork
 *                  
 * iwork            Integer work array
 *
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSize      The number of vectors to restart with. If negative,
 *                  use dynamic thick restart.
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * ldV              The leading dimension of V
 *
 * W                A*V
 *
 * ldW              The leading dimension of W
 *
 * H                The projection V'*A*V
 *
 * ldH              The leading dimension of H
 *
 * VtBV             The projection V'*B*V
 *
 * ldVtBV           The leading dimension of VtBV
 *
 * Q, R             The factors of the QR decomposition of (A - targetShift*B)*V
 *
 * ldQ, ldR         The leading dimension of Q and R
 *
 * QtV              = Q'*V
 *
 * ldQtV            The leading dimension of QtV
 *
 * hU               The left coefficient vectors or the eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of the input hU
 *
 * newldhU          The leading dimension of the output hU
 *
 * hVecs            The coefficient vectors
 *
 * ldhVecs          The leading dimension of the input hVecs
 *
 * newldhVecs       The leading dimension of the output hVecs
 *
 * hVals            The eigenvalues of H
 *
 * hSVals           The singular values of R
 *
 * hVecsPerm        The permutation that orders hVals and hVecs as primme.target
 *
 * restartPerm      The permutation applied to the columns of hVecs before restarting
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * evecs            The converged Ritz vectors. Without locking, all converged
 *                  eigenvectors are copied from V to evecs if skew projections
 *                  are required.
 *
 * ldevecs          The leading dimension of evecs
 *
 * evecsHat         K^{-1}evecs
 *
 * ldevecsHat       The leading dimension of evecsHat
 *
 * M                evecs'*evecsHat
 *
 * ldM              The leading dimension of M
 *
 * UDU              The factorization of M
 *
 * ldUDU            The leading dimension of UDU
 *
 * ipivot           The pivot array of the UDU factorization
 *
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 *
 * numArbitraryVecs On input, the number of coefficients vectors that do
 *                  not correspond to solutions of the projected problem.
 *                  They appear in the first numArbitraryVecs positions of hVecs.
 *                  On output, the number of such vectors that are in the
 *                  restarted basis.
 *
 * Return value
 * ------------
 * int   > 0 the restart size   
 *        -2 restart_H failed
 *        -4 factorization of M failed
 *        -5 flags do not correspond to converged pairs in pseudolocking
 *       
 ******************************************************************************/
 
static int restart_projection_Sprimme(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *H, int ldH, SCALAR *VtBV, int ldVtBV, SCALAR *Q,
      PRIMME_INT ldQ, PRIMME_INT nLocal, SCALAR *R, int ldR, SCALAR *QtV,
      int ldQtV, SCALAR *hU, int ldhU, int newldhU,
      int indexOfPreviousVecsBeforeRestart, SCALAR *hVecs, int ldhVecs,
      int newldhVecs, REAL *hVals, REAL *hSVals, int *restartPerm,
      int *hVecsPerm, int restartSize, int basisSize, int numPrevRetained,
      int indexOfPreviousVecs, SCALAR *evecs, int *evecsSize,
      PRIMME_INT ldevecs, SCALAR *evecsHat, PRIMME_INT ldevecsHat, SCALAR *M,
      int ldM, SCALAR *UDU, int ldUDU, int *ipivot, int *targetShiftIndex,
      int numConverged, int *numArbitraryVecs, SCALAR *hVecsRot, int ldhVecsRot,
      size_t *rworkSize, SCALAR *rwork, int iworkSize, int *iwork,
      double machEps, primme_params *primme) {

   /* -------------------------------------------------------- */
   /* Restart projected problem matrices H and R               */
   /* -------------------------------------------------------- */

   switch (primme->projectionParams.projection) {
   case primme_proj_RR:
      if (!VtBV) {
         CHKERR(restart_RR(H, ldH, hVecs, ldhVecs, newldhVecs, hVals,
                  restartSize, basisSize, numConverged, numPrevRetained,
                  indexOfPreviousVecs, hVecsPerm, targetShiftIndex, machEps,
                  rworkSize, rwork, iworkSize, iwork, primme), -1);
      }
      else {
         CHKERR(restart_RR_explicit(V, ldV, W, ldW, nLocal, H, ldH, VtBV,
                  ldVtBV, hVecs, newldhVecs, hVals, restartSize,
                  numConverged, targetShiftIndex, machEps, rworkSize, rwork,
                  iworkSize, iwork, primme), -1);
      }
      break;

   case primme_proj_harmonic:
      CHKERR(restart_harmonic(V, ldV, W, ldW, H, ldH, Q, ldQ, nLocal, R, ldR,
            QtV, ldQtV, hU, ldhU, newldhU, hVecs, ldhVecs, newldhVecs, hVals,
            hSVals, restartPerm, hVecsPerm, restartSize, basisSize,
            numPrevRetained, indexOfPreviousVecs, targetShiftIndex,
            numConverged, numArbitraryVecs, hVecsRot, ldhVecsRot, rworkSize,
            rwork, iworkSize, iwork, machEps, primme), -1);
      break;

   case primme_proj_refined:
      CHKERR(restart_refined(V, ldV, W, ldW, H, ldH, Q, ldQ, nLocal, R, ldR, hU,
            ldhU, newldhU, indexOfPreviousVecsBeforeRestart,
            hVecs, ldhVecs, newldhVecs, hVals, hSVals, restartPerm, hVecsPerm,
            restartSize, basisSize, numPrevRetained, indexOfPreviousVecs,
            targetShiftIndex, numConverged, numArbitraryVecs, hVecsRot,
            ldhVecsRot, rworkSize, rwork, iworkSize, iwork, machEps, primme),
            -1);
      break;

   default:
      assert(0);
   }

   if (evecsHat) {
      int numRecentlyConverged = numConverged - *evecsSize;

      /* Return memory requirement */
      if (H == NULL) {
         CHKERR(update_projection_Sprimme(NULL, 0, NULL, 0, NULL, 0, nLocal,
                  *evecsSize, basisSize, NULL, rworkSize, 1/*symmetric*/,
                  primme), -1);
         CHKERR(UDUDecompose_Sprimme(NULL, 0, NULL, 0, NULL, *evecsSize,
                  NULL, rworkSize, primme), -1);
         return 0;
      }

      /* Compute K^{-1}x for all newly locked eigenvectors */

      /* TODO: primme.shiftsForPreconditioner is undefined at that point;
         maybe it makes sense to always set NULL shiftsForPreconditioner
         when SkewQ is enabled to force the same preconditioner. */
      CHKERR(applyPreconditioner_Sprimme(
               &evecs[ldevecs*(*evecsSize+primme->numOrthoConst)], nLocal,
               ldevecs,
               &evecsHat[ldevecsHat*(*evecsSize+primme->numOrthoConst)],
               ldevecsHat, numRecentlyConverged, primme), -1);

      /* Update the projection evecs'*evecsHat now that evecs and evecsHat   */
      /* have been expanded by numRecentlyConverged columns.  Required       */
      /* workspace is numLocked*numEvals.  The most ever needed would be     */
      /* maxBasisSize*numEvals.                                              */

      CHKERR(update_projection_Sprimme(evecs, primme->nLocal, evecsHat,
               primme->nLocal, M, ldM, nLocal, *evecsSize+primme->numOrthoConst,
               numRecentlyConverged, rwork, rworkSize, 1/*symmetric*/, primme),
               -1);
      *evecsSize = numConverged;

      CHKERR(UDUDecompose_Sprimme(M, ldM, UDU, ldUDU, ipivot,
               *evecsSize+primme->numOrthoConst, rwork, rworkSize, primme), -1);
   }

   return 0;
}

/*******************************************************************************
 * Function restart_RR - This routine is used to recompute H = V'*A*V once V 
 *   has been restarted.  If no coefficient vectors from the previous iteration
 *   have been retained, then the restarted H will be diagonal and the 
 *   new eigenvectors (coefficient vectors) of H will be the standard basis
 *   vectors.  If previous coefficient vectors have been retained, then H will 
 *   contain the numPrevRetained x numPrevRetained submatrix 
 *   previousHVecs'*H*previousHvecs and the rest of the elements of H will be
 *   along the diagonal.
 *   
 *
 * INPUT PARAMETERS
 * ----------------
 * restartSize   Number of vectors the basis was restarted with
 * 
 * basisSize     Maximum size of the basis V
 *
 * numLocked     The number of Ritz vectors that have been locked 
 *
 * numPrevRetained The number of vectors retained from the previous iteration
 *
 * indexOfPreviousVecs  The index within hVecs where the previous vectors were
 *                      inserted.  Its also where the overlap matrix
 *                      previousHVecs'*H*previousHvecs will be inserted within
 *                      the restarted H.
 *
 * hVecsPerm     The permutation that orders hVals and hVecs as primme.target
 *
 * iwork         Integer work array
 *
 * rwork         Work array.  Necessary only when coefficient vectors from the
 *               previous iteration are retained.
 *
 * rworkSize     Size of rwork
 *
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * H      Will contain H = V'*A*V given the restarted V.  H will be
 *        diagonal since V will contain only Ritz vectors, unless previous
 *        vectors are retained.  In that case, it will be diagonal except for
 *        a numPrevRetained x numPrevRetained submatrix inserted at
 *        H(numPrevRetained, numPrevRetained)
 *
 * ldH    The leading dimension of H
 *
 * VtBV   Will contain VtBV = V'*B*V given the restarted V. VtBV will be
 *        diagonal since V will contain only Ritz vectors, unless previous
 *        vectors are retained.  In that case, it will be diagonal except for
 *        a numPrevRetained x numPrevRetained submatrix inserted at
 *        VtBV(numPrevRetained, numPrevRetained)
 *
 * ldVtBV The leading dimension of VtBV
 *
 * hVecs  If the new H is diagonal, then it will contain the standard basis
 *        vectors.  If previous coefficient vectors are retained, then 
 *        restartSize - numPrevRetained of the vectors will be standard basis
 *        vectors.  The remaining numPrevRetained vectors will contain
 *        numPrevRetained non-zero elements corresponding to the 
 *        numPrevRetined x numPrevRetained submatrix.
 * 
 * ldhVecs  The leading dimension of the input hVecs
 *
 * newldhVecs  The leading dimension of the output hVecs
 *
 * hVals  The eigenvalues of the restarted H
 * 
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 * 
 * Return value
 * ------------
 * Error code: 0 upon success
 *            -1 eigenvalues of submatrix could not be computed
 *
 ******************************************************************************/

static int restart_RR(SCALAR *H, int ldH, SCALAR *hVecs, int ldhVecs,
      int newldhVecs, REAL *hVals, int restartSize,
      int basisSize, int numLocked, int numPrevRetained,
      int indexOfPreviousVecs, int *hVecsPerm, int *targetShiftIndex,
      double machEps, size_t *rworkSize, SCALAR *rwork, int iworkSize,
      int *iwork, primme_params *primme) {

   int i, j;          /* Loop variables                                       */
   int orderedIndexOfPreviousVecs;  /* index of prev. vecs after applying hVecsPerm */
   double aNorm = primme?max(primme->aNorm, primme->stats.estimateLargestSVal):0.0;

   /* Return memory requirement */

   if (H == NULL) {
      CHKERR(compute_submatrix_Sprimme(NULL, numPrevRetained, 0, NULL,
               basisSize, 0, NULL, 0, NULL, rworkSize), -1);
      CHKERR(solve_H_Sprimme(NULL, numPrevRetained, 0, NULL, 0, NULL, 0, NULL,
               0, NULL, 0, NULL, 0, NULL, NULL, numLocked, 0.0, rworkSize, NULL,
               0, iwork, primme), -1);
      return 0;
   }

   /* ---------------------------------------------------------------------- */
   /* If coefficient vectors from the previous iteration were retained, then */
   /* insert the computed overlap matrix into the restarted H                */
   /* ---------------------------------------------------------------------- */

   CHKERR(compute_submatrix_Sprimme(&hVecs[ldhVecs*indexOfPreviousVecs],
            numPrevRetained, ldhVecs, H,
            basisSize, ldH, &H[ldH*indexOfPreviousVecs+indexOfPreviousVecs],
            ldH, rwork, rworkSize), -1);

   /* ----------------------------------------------------------------- */
   /* Y = V*hVecs([0:indexOfPreviousVecs-1 \                            */
   /*                 indexOfPreviousVecs+numPrevRetained:restartSize]) */
   /* are the new Ritz vectors so the new H and VtBV will be diag(hVals)*/
   /* and I respectively. Note it is only stored the upper triangular   */
   /* part of H and VtBV.                                               */
   /* ----------------------------------------------------------------- */

   for (j=0; j < indexOfPreviousVecs; j++) {
      for (i=0; i < j; i++) {
         H[ldH*j+i] = 0.0;
      }
      H[ldH*j+j] = hVals[j];
   }
   for (j=indexOfPreviousVecs; j<indexOfPreviousVecs+numPrevRetained; j++) {
      for (i=0; i < indexOfPreviousVecs; i++) {
         H[ldH*j+i] = 0.0;
      }
   }
   for (j=indexOfPreviousVecs+numPrevRetained; j < restartSize; j++) {
      for (i=0; i < j; i++) {
         H[ldH*j+i] = 0.0;
      }
      H[ldH*j+j] = hVals[j];
   }

   /* ---------------------------------------------------------------------- */
   /* Solve the whole matrix when the targetShift has to change              */ 
   /* ---------------------------------------------------------------------- */

   if (targetShiftIndex && primme->targetShifts && (*targetShiftIndex < 0
            || fabs(primme->targetShifts[*targetShiftIndex]
               - primme->targetShifts[min(primme->numTargetShifts-1,
                  numLocked)])) > machEps*aNorm) {

      *targetShiftIndex = min(primme->numTargetShifts-1, numLocked);

      CHKERR(solve_H_Sprimme(H, restartSize, ldH, NULL, 0, NULL, 0, NULL,
               0, NULL, 0, hVecs, newldhVecs, hVals, NULL, numLocked, machEps,
               rworkSize, rwork, iworkSize, iwork, primme), -1);

      return 0;
   }


   /* --------------------------------------------------------------------- */
   /* Find the index of the first column with a retained vector. All        */
   /* retained vectors should be together.                                  */
   /* --------------------------------------------------------------------- */

   for (i=0, orderedIndexOfPreviousVecs=restartSize; i<restartSize; i++) {
      if (hVecsPerm[i] == indexOfPreviousVecs) {
         orderedIndexOfPreviousVecs = i;
         break;
      }
   }
   assert(orderedIndexOfPreviousVecs != restartSize
         || indexOfPreviousVecs >= restartSize);

   for (i=orderedIndexOfPreviousVecs+1;
         i<orderedIndexOfPreviousVecs+numPrevRetained; i++)
      assert(hVecsPerm[i-1]+1 == hVecsPerm[i]);

   /* --------------------------------------------------------------------- */
   /* Given the above H, we know the eigenvectors of H will be the          */
   /* canonical basis except for the retained vectors.                      */
   /* --------------------------------------------------------------------- */

   for (j=0; j < restartSize; j++) {
      for (i=0; i < restartSize; i++) {
          hVecs[newldhVecs*j+i] = 0.0;
      }
      hVecs[newldhVecs*j+hVecsPerm[j]] = 1.0;
   }

   /* Apply permutation hVecsPerm to hVals */
   permute_vecs_Rprimme(hVals, 1, restartSize, 1, hVecsPerm, (REAL*)rwork,
         iwork);

   /* ---------------------------------------------------------------------- */
   /* Solve the overlap matrix corresponding for the retained vectors to     */ 
   /* compute the coefficient vectors.                                       */
   /* ---------------------------------------------------------------------- */

   CHKERR(solve_H_Sprimme(
            &H[ldH*indexOfPreviousVecs+indexOfPreviousVecs], numPrevRetained,
            ldH, NULL, 0, NULL, 0, NULL, 0, NULL, 0,
            &hVecs[newldhVecs*orderedIndexOfPreviousVecs+indexOfPreviousVecs],
            newldhVecs, &hVals[orderedIndexOfPreviousVecs], NULL, numLocked,
            machEps, rworkSize, rwork, iworkSize, iwork, primme), -1);

   return 0;
}

/*******************************************************************************
 * Function restart_RR_explicit - This routine is used to recompute H = V'*A*V
 *   VtBV = V'*B*V explcitly.
 *
 * INPUT PARAMETERS
 * ----------------
 * restartSize   Number of vectors the basis was restarted with
 * 
 * basisSize     Maximum size of the basis V
 *
 * numLocked     The number of Ritz vectors that have been locked 
 *
 * numPrevRetained The number of vectors retained from the previous iteration
 *
 * indexOfPreviousVecs  The index within hVecs where the previous vectors were
 *                      inserted.  Its also where the overlap matrix
 *                      previousHVecs'*H*previousHvecs will be inserted within
 *                      the restarted H.
 *
 * hVecsPerm     The permutation that orders hVals and hVecs as primme.target
 *
 * iwork         Integer work array
 *
 * rwork         Work array.  Necessary only when coefficient vectors from the
 *               previous iteration are retained.
 *
 * rworkSize     Size of rwork
 *
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * H      Will contain H = V'*A*V given the restarted V.  H will be
 *        diagonal since V will contain only Ritz vectors, unless previous
 *        vectors are retained.  In that case, it will be diagonal except for
 *        a numPrevRetained x numPrevRetained submatrix inserted at
 *        H(numPrevRetained, numPrevRetained)
 *
 * ldH    The leading dimension of H
 *
 * VtBV   Will contain VtBV = V'*B*V given the restarted V. VtBV will be
 *        diagonal since V will contain only Ritz vectors, unless previous
 *        vectors are retained.  In that case, it will be diagonal except for
 *        a numPrevRetained x numPrevRetained submatrix inserted at
 *        VtBV(numPrevRetained, numPrevRetained)
 *
 * ldVtBV The leading dimension of VtBV
 *
 * hVecs  If the new H is diagonal, then it will contain the standard basis
 *        vectors.  If previous coefficient vectors are retained, then 
 *        restartSize - numPrevRetained of the vectors will be standard basis
 *        vectors.  The remaining numPrevRetained vectors will contain
 *        numPrevRetained non-zero elements corresponding to the 
 *        numPrevRetined x numPrevRetained submatrix.
 * 
 * ldhVecs  The leading dimension of the input hVecs
 *
 * newldhVecs  The leading dimension of the output hVecs
 *
 * hVals  The eigenvalues of the restarted H
 * 
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 * 
 * Return value
 * ------------
 * Error code: 0 upon success
 *            -1 eigenvalues of submatrix could not be computed
 *
 ******************************************************************************/

static int restart_RR_explicit(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, PRIMME_INT nLocal, SCALAR *H, int ldH, SCALAR *VtBV,
      int ldVtBV, SCALAR *hVecs, int ldhVecs, REAL *hVals, int restartSize,
      int numLocked, int *targetShiftIndex, double machEps, size_t *rworkSize,
      SCALAR *rwork, int iworkSize, int *iwork, primme_params *primme) {

   CHKERR(update_projection_Sprimme(V, ldV, W, ldW, H, ldH, nLocal, 0,
            restartSize, rwork, rworkSize, 1/*symmetric*/, primme), -1);

   if (VtBV) CHKERR(update_projection_Sprimme(V, ldV, V, ldV, VtBV, ldVtBV,
            nLocal, 0, restartSize, rwork, rworkSize, 1/*symmetric*/, primme),
         -1);

   if (primme->numTargetShifts > 0 && targetShiftIndex)
      *targetShiftIndex = min(primme->numTargetShifts-1, numLocked);

   CHKERR(solve_H_Sprimme(H, restartSize, ldH, VtBV, ldVtBV, NULL, 0, NULL, 0,
            NULL, 0, hVecs, ldhVecs, hVals, NULL,
            targetShiftIndex?*targetShiftIndex:0, machEps,
            rworkSize, rwork, iworkSize, iwork, primme), -1);

   return 0;
}

/*******************************************************************************
 * Function restart_refined - This routine is used to recompute the QR decomposition
 *    of W (=A*V), after V being replaced by V*hVecs. The coefficient vectors
 *    except hVecs(indexOfPrevRetained:indexOfPrevRetained+numPrevRetained) are
 *    right singular vectors of R. The output R for these columns will be
 *    a diagonal matrix with the singular values as diagonal elements. The
 *    output Q for these columns will be the update with the left singular
 *    vectors. The other columns of R and Q are recomputed properly without
 *    requiring explicitly recompute the QR factorization.
 *
 *    Also H = V'*A*V and QtV = Q'*V are recomputed properly.
 *   
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal           Number of rows of V assigned to the node
 *
 * basisSize        Size of the basis V
 *
 * numConverged     The number of converged eigenpairs without locking
 *
 * numPrevRetained  The number of coefficient vectors from previousHVecs in hVecs
 *  
 * indexOfPreviousVecs The first column in hVecs that has a vector from previousHVecs
 *
 * indexOfPreviousVecsBeforeRestart The number of columns in hVecs that previousHVecs
 *                  where orthogonalized against
 *
 * rwork            Real work array
 *
 * rworkSize        Size of rwork
 *                  
 * iwork            Integer work array
 *
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSize      The number of vectors to restart with. If negative,
 *                  use dynamic thick restart.
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * ldV              The leading dimension of V
 *
 * W                A*V
 *
 * ldW              The leading dimension of W
 *
 * H                The projection V'*A*V
 *
 * ldH              The leading dimension of H
 *
 * Q, R             The factors of the QR decomposition of (A - targetShift*B)*V
 *
 * ldQ, ldR         The leading dimension of Q and R
 *
 * hU               The left singular vectors of R
 *
 * ldhU             The leading dimension of the input hU
 *
 * newldhU          The leading dimension of the output hU
 *
 * hVecs            The eigenvectors of H
 *
 * ldhVecs          The leading dimension of the input hVecs
 *
 * newldhVecs       The leading dimension of the output hVecs
 *
 * hVals            The eigenvalues of H
 *
 * hSVals           The singular values of R
 *
 * hVecsPerm        The permutation that orders hVals and hVecs as primme.target
 *
 * restartPerm      The permutation applied to the columns of hVecs before restarting
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 *
 * numArbitraryVecs On input, the number of coefficients vectors that do
 *                  not correspond to solutions of the projected problem.
 *                  They appear in the first numArbitraryVecs positions of hVecs.
 *                  On output, the number of such vectors that are in the
 *                  restarted basis.
 * 
 * Return value
 * ------------
 * Error code: 0 upon success
 *            -1 eigenvalues of submatrix could not be computed
 *
 ******************************************************************************/

static int restart_refined(SCALAR *V, PRIMME_INT ldV, SCALAR *W, PRIMME_INT ldW,
      SCALAR *H, int ldH, SCALAR *Q, PRIMME_INT ldQ, PRIMME_INT nLocal,
      SCALAR *R, int ldR, SCALAR *hU, int ldhU, int newldhU,
      int indexOfPreviousVecsBeforeRestart, SCALAR *hVecs, int ldhVecs,
      int newldhVecs, REAL *hVals, REAL *hSVals, int *restartPerm,
      int *hVecsPerm, int restartSize, int basisSize, int numPrevRetained,
      int indexOfPreviousVecs, int *targetShiftIndex, int numConverged,
      int *numArbitraryVecs, SCALAR *hVecsRot, int ldhVecsRot,
      size_t *rworkSize, SCALAR *rwork, int iworkSize, int *iwork,
      double machEps, primme_params *primme) {

   int i, j;          /* Loop variables                                       */
   int mhVecsRot0;    /* Number of rows of hVecsRot0                          */
   int newNumArbitraryVecs; /* arbitrary vectors after restarting             */
   int nRegular;      /* restarted vectors that weren't retained              */
   SCALAR *hVecsRot0;  /* Part of hVecsRot                                   */
   SCALAR *RPrevhVecs;
   SCALAR *rwork0;
   size_t rworkSize0;
   int *restartPerm0;
   int *invhVecsPerm;
   int *iwork0;
   double aNorm = primme?max(primme->aNorm, primme->stats.estimateLargestSVal):0.0;

   /* Return memory requirement */
 
   if (V == NULL) {
      CHKERR(compute_submatrix_Sprimme(NULL, basisSize, 0, NULL, basisSize,
               0, NULL, 0, NULL, rworkSize), -1);
      CHKERR(update_Q_Sprimme(NULL, nLocal, 0, NULL, 0, NULL, 0, NULL, 0,
               0.0, 0, basisSize, NULL, rworkSize, 0.0, primme), -1);
      CHKERR(update_projection_Sprimme(NULL, 0, NULL, 0, NULL, 0, nLocal, 0, basisSize,
            NULL, rworkSize, 0/*unsymmetric*/, primme), -1);
      /* Workspace for  R(indexOfPrevVecs:) = R * hVecs(indexOfPrevVecs:) */
      /* The workspace for permute_vecs(hU) is basisSize */
      *rworkSize = max(*rworkSize, (size_t)basisSize*(size_t)basisSize);
      *rworkSize = max(*rworkSize,
            (size_t)Num_update_VWXR_Sprimme(NULL, NULL, nLocal, basisSize,
               0, NULL, basisSize, 0, NULL,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0, NULL,
               NULL, 0, 0,
               NULL, 0, primme));
      CHKERR(solve_H_Sprimme(NULL, basisSize, 0, NULL, 0, NULL, 0, NULL, 0,
               NULL, 0, NULL, 0, NULL, NULL, numConverged, 0.0, rworkSize,
               NULL, 0, iwork, primme), -1);
      return 0;
   }

   /* ------------------------------- */
   /* Replace H by hVecs' * H * hVecs */
   /* ------------------------------- */

   if (H) compute_submatrix_Sprimme(hVecs, restartSize, ldhVecs, H,
         basisSize, ldH, H, ldH, rwork, rworkSize);

   /* -------------------------------------- */
   /* Quick exit if the target has changed   */
   /* -------------------------------------- */

   /* NOTE: Force to pass the next condition if you want to rebuild the QR    */
   /* factorization at every restart.                                         */

   /* NOTE: keep the same condition here as in main_iter */

   if (*targetShiftIndex < 0 || fabs(primme->targetShifts[*targetShiftIndex]
            - primme->targetShifts[min(primme->numTargetShifts-1, numConverged)])
         > machEps*aNorm) {

      *targetShiftIndex = min(primme->numTargetShifts-1, numConverged);

      CHKERR(update_Q_Sprimme(V, nLocal, ldV, W, ldW, Q, ldQ, R, ldR,
               primme->targetShifts[*targetShiftIndex], 0,
               restartSize, rwork, rworkSize, machEps, primme), -1);

      CHKERR(solve_H_Sprimme(H, restartSize, ldH, NULL, 0, R, ldR, NULL, 0, hU,
               newldhU, hVecs, newldhVecs, hVals, hSVals, numConverged, machEps,
               rworkSize, rwork, iworkSize, iwork, primme), -1);

      *numArbitraryVecs = 0;

      return 0;
   }

   /* restartPerm0 = hVecsPerm(restartPerm)                             */

   restartPerm0 = iwork;
   for (i=0; i<restartSize; i++) {
      restartPerm0[i] = restartPerm[hVecsPerm[i]];
      assert(restartPerm0[i] >= 0 && restartPerm0[i] < basisSize);
   }

   /* ----------------------------------------------------------------- */
   /* Update numArbitraryVecs as the number of arbitrary vectors in     */
   /* the restarted basis.                                              */
   /* ----------------------------------------------------------------- */

   for (i=newNumArbitraryVecs=0; i < restartSize-numPrevRetained; i++) {
      if (restartPerm0[i] < *numArbitraryVecs) {
         newNumArbitraryVecs++;
      }
   }

   /* Check all arbitrary vectors will be together and at the beginning after restart */

   for (i=0; i<newNumArbitraryVecs; i++)
      assert(restartPerm0[i] < *numArbitraryVecs);

   /* -------------------------------------------------------------------- */
   /* Note that hVecs=[hV*hVecsRot prevhVecs], where hV are the right      */
   /* singular vectors of R and prevhVecs are the retained coefficient     */
   /* vectors after orthogonalizing them against the first                 */
   /* indexOfPreviousVecsBeforeRestart coefficient vectors in hVecs.       */
   /*                                                                      */
   /* In restart V is replaced by by V*Y and W by W*Y.                     */
   /* To restart the QR decomposition of (A-\tau I)*V, Q is replaced by    */
   /* Q*Z and R by T, where Z and T are the QR decomposition of R*Y.       */
   /* We avoid computing R*Y because the errors go into the direction of   */
   /* of the largest singular values. Note that                            */
   /*                                                                      */
   /* R*Y = R*[ hV*hVecsRot prevhVecs ]                                    */
   /*     = [ hU*diag(hSVals)*hVecsRot R*prevhVecs ].                      */
   /*                                                                      */
   /* Note that the code works with R(:,restartPerm0) instead of R,        */
   /* because the new arbitrary columns are all together at the leading    */
   /* columns of R(:,restartPerm0).                                        */
   /* -------------------------------------------------------------------- */


   /* RPrevhVecs = R * hVecs(indexOfPreviousVecs:indexOfPreviousVecs+numPrevRetained) */

   assert(*rworkSize >= (size_t)(numPrevRetained*basisSize));
   RPrevhVecs = rwork; rwork0 = rwork + numPrevRetained*basisSize;
   rworkSize0 = *rworkSize - (size_t)numPrevRetained*basisSize;
   Num_gemm_Sprimme("N", "N", basisSize, numPrevRetained, basisSize, 1.0,
         R, ldR, &hVecs[ldhVecs*indexOfPreviousVecs], ldhVecs, 0.0,
         RPrevhVecs, basisSize);
 
   /* Do hVecsRot0 = diag(hSvals)*hVecsRot(hVecsPerm(restartPerm)) limited */
   /* to the columns 0:newNumArbitraryVecs-1 (in 3 steps)                  */
   
   nRegular = restartSize - numPrevRetained;
   for (i=0, mhVecsRot0=*numArbitraryVecs; i < nRegular; i++)
      mhVecsRot0 = max(mhVecsRot0, restartPerm0[i]+1);
   assert(rworkSize0 >= (size_t)mhVecsRot0*nRegular);
   hVecsRot0 = rwork0; rwork0 += mhVecsRot0*nRegular;
   rworkSize0 -= (size_t)mhVecsRot0*nRegular;

   /* 1) hVecsRot0 = hVecsRot(hVecsPerm(restartPerm(0:newNumArbitraryVecs-1))) */
   Num_zero_matrix_Sprimme(hVecsRot0, mhVecsRot0, nRegular, mhVecsRot0);
   Num_copy_matrix_columns_Sprimme(hVecsRot, *numArbitraryVecs,
         restartPerm0, newNumArbitraryVecs, ldhVecsRot, hVecsRot0, NULL,
         mhVecsRot0);

   /* 2) hVecsRot0 = diag(hSVals)*hVecsRot0 */
   for (i=0; i<newNumArbitraryVecs; i++) {
      for (j=0; j<*numArbitraryVecs; j++) {
         hVecsRot0[mhVecsRot0*i+j] *= hSVals[j];
      }
   }

   /* 3) hVecsRot0(:,c) = diag(hSVals)*I(:,restartPerm0(c))  */
   /*  for c = newNumArbVecs:nRegular-1                   */
   for (i=newNumArbitraryVecs; i<nRegular; i++) {
      hVecsRot0[mhVecsRot0*i+restartPerm0[i]] = hSVals[restartPerm0[i]];
   }

   /* R = zeros() */
   Num_zero_matrix_Sprimme(R, primme->maxBasisSize, primme->maxBasisSize,
         ldR);
   
   /* [hVecsRot0, R] = ortho(hVecsRot0) */
   CHKERR(ortho_Sprimme(hVecsRot0, mhVecsRot0, R, ldR, 0, nRegular-1, NULL,
         0, 0, mhVecsRot0, primme->iseed, machEps, rwork0, &rworkSize0, NULL),
         -1);

   /* hU = hU * hVecsRot0 */
   Num_gemm_Sprimme("N", "N", basisSize, nRegular, mhVecsRot0, 1.0, hU,
         ldhU, hVecsRot0, mhVecsRot0, 0.0, rwork0, basisSize);
   Num_copy_matrix_Sprimme(rwork0, basisSize, nRegular, basisSize, hU,
         basisSize);

   /* hU = [hU RPrevhVecs] */
   Num_copy_matrix_Sprimme(RPrevhVecs, basisSize, numPrevRetained,
         basisSize, &hU[basisSize*nRegular], basisSize);

   /* [hU, R] = ortho(hU, nRegular:nRegular+numPrevRetained-1) */
   CHKERR(ortho_Sprimme(hU, basisSize, R, ldR, nRegular,
            nRegular+numPrevRetained-1, NULL, 0, 0, basisSize, primme->iseed,
            machEps, rwork, rworkSize, NULL), -1);

   /* Zero upper triangular part of R(:,newNumArbVecs:nRegular-1) for the     */
   /* columns that restartPerm0[i] >= numArbVecs                              */
   for (i=newNumArbitraryVecs; i<nRegular; i++) {
      if (restartPerm0[i] >= *numArbitraryVecs) {
         for (j=0; j<=i; j++) {
            R[ldR*i+j] = 0.0;
         }
         R[ldR*i+i] = hSVals[restartPerm0[i]];
      }
   }

   /* When the retained coefficient vectors were orthogonalized against all   */
   /* arbitrary vectors then R*prevhVecs is orthogonal to                     */
   /* hU(0:indexOfPreviousVecsBeforeRestart) and that rows of R corresponding */
   /* to the retained vectors should be zero.                                 */

    if (*numArbitraryVecs <= indexOfPreviousVecsBeforeRestart) {
      /* Zero R(0:nRegular-1,nRegular:restartSize) */
      Num_zero_matrix_Sprimme(&R[ldR*nRegular], nRegular, numPrevRetained,
            ldR);
   }

   /* ----------------------------------- */
   /* Restart Q by replacing it with Q*hU */
   /* ----------------------------------- */

   CHKERR(Num_update_VWXR_Sprimme(Q, NULL, nLocal, basisSize, ldQ, hU,
            restartSize,
            basisSize, NULL,
            Q, 0, restartSize, ldQ,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0, NULL,
            NULL, 0, 0,
            rwork, TO_INT(*rworkSize), primme), -1);

   /* ---------------------------------------------------------------------- */
   /* R may lost the triangular structure after the previous permutation, so */
   /* the easiest way to recompute the projected matrices hVecs, hU and      */
   /* and hSVals is to call solve_H. For LOBPCG style, it means to           */
   /* solve twice the projected problem, and it may be a performance issue.  */
   /*                                                                        */
   /* TODO: solve only the columns of R corresponding to the new             */
   /*       arbitrary vectors and the retained vectors.                      */
   /* ---------------------------------------------------------------------- */

   /* Compute the singular value decomposition of R.                         */
   /* Notice that hVals aren't changed. They are updated as a permutation    */
   /* over the original hVals.                                               */
   /* NOTE: the eigenvalues of the retained pairs aren't correctly computed. */

   assert(*rworkSize >= (size_t)restartSize);
   rworkSize0 = *rworkSize - (size_t)restartSize;
   CHKERR(solve_H_Sprimme(H, restartSize, ldH, NULL, 0, R, ldR, NULL, 0, hU,
            newldhU, hVecs, newldhVecs, (REAL*)rwork, hSVals, numConverged,
            machEps, &rworkSize0, rwork+restartSize, iworkSize, iwork, primme),
         -1);

   permute_vecs_Rprimme(hVals, 1, restartSize, 1, hVecsPerm, (REAL*)rwork, iwork);

   /* ---------------------------------------------------------------------- */
   /* Permute back the columns of R                                          */
   /* ---------------------------------------------------------------------- */

   /* hVecsPerm(invhVecsPerm) = 1:n */
   invhVecsPerm = iwork; iwork0 = iwork + restartSize;
   for (i=0; i<restartSize; i++) {
      invhVecsPerm[hVecsPerm[i]] = i;
   }

   /* R = R(:, invhVecsPerm) */
   permute_vecs_Sprimme(R, restartSize, restartSize, ldR, invhVecsPerm,
         rwork, iwork0);

   /* ----------------------------------------------------------------- */
   /* After all the changes in hVecs and R new arbitrary vectors may    */
   /* have been introduced. When the retained coefficient vectors are   */
   /* orthogonalized against all arbitrary vectors then the new number  */
   /* of arbitrary vectors is at most the largest index out of order    */
   /* in hVecsPerm, plus one. Otherwise the easiest way is to consider  */
   /* all restarted vectors as arbitrary vectors.                       */
   /* ----------------------------------------------------------------- */

   if (*numArbitraryVecs <= indexOfPreviousVecsBeforeRestart) {
      for (i=*numArbitraryVecs=newNumArbitraryVecs; i<restartSize; i++)
         if (hVecsPerm[i] != i) *numArbitraryVecs=i+1;
   }
   else {
      *numArbitraryVecs = restartSize;
   }

   /* ----------------------------------------------------------------------- */
   /* We want to compute hVecsRot so that it corresponds to the rotation of   */
   /* the right singular vectors of R that produces I(:,hVecsPerm), because   */
   /* the desired coefficient vectors were already decided in the last        */
   /* iteration. The current hVecs is the right singular vectors of           */
   /* R*I(:,hVecsPerm). Notice R*I(:,hVecsPerm) = hU * S * hV'*I(:,hVecsPerm).*/
   /* Setting hVecsRot as current hVecs'=hV'*I(:,hVecsPerm) will satisfy that */
   /* hVecs = hV * hVecsRot, with hVecs=I(:,hVecsPerm).                       */
   /* ----------------------------------------------------------------------- */

   /* hVecsRot <- hVecs' for arbitrary vectors */

   Num_zero_matrix_Sprimme(hVecsRot, primme->maxBasisSize, primme->maxBasisSize,
         ldhVecsRot);
   for (j=0; j < *numArbitraryVecs; j++) {
      for (i=0; i<restartSize; i++) {
         hVecsRot[ldhVecsRot*j+i] = CONJ(hVecs[newldhVecs*i+j]);
      }
   }
 
   /* hVecs <- I for arbitrary vectors */

   Num_zero_matrix_Sprimme(hVecs, restartSize, *numArbitraryVecs, newldhVecs);
   for (j=0; j < *numArbitraryVecs; j++) {
      hVecs[newldhVecs*j+hVecsPerm[j]] = 1.0;
   }

   return 0;
}

/*******************************************************************************
 * Function restart_harmonic - This routine is used to recompute the QR decomposition
 *    of W (=A*V), after V being replaced by V*hVecs. 
 *    Also H = V'*A*V and QtV = Q'*V are recomputed properly.
 *   
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal           Number of rows of V assigned to the node
 *
 * basisSize        Size of the basis V
 *
 * numConverged     The number of converged eigenpairs without locking
 *
 * numPrevRetained  The number of coefficient vectors from previousHVecs in hVecs
 *  
 * indexOfPreviousVecs The first column in hVecs that has a vector from previousHVecs
 *
 * rwork            Real work array
 *
 * rworkSize        Size of rwork
 *                  
 * iwork            Integer work array
 *
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * restartSize      The number of vectors to restart with. If negative,
 *                  use dynamic thick restart.
 *
 * V                The orthonormal basis. After restart, contains Ritz vectors
 *                  plus the orthogonal components from numPrevRetained Ritz 
 *                  vectors from the penultimate step.
 *
 * ldV              The leading dimension of V
 *
 * W                A*V
 *
 * ldW              The leading dimension of W
 *
 * H                The projection V'*A*V
 *
 * ldH              The leading dimension of H
 *
 * Q, R             The factors of the QR decomposition of (A - targetShift*B)*V
 *
 * ldQ, ldR         The leading dimension of Q and R
 *
 * QtV              = Q'*V
 *
 * ldQtV            The leading dimension of QtV
 *
 * hU               The eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of the input hU
 *
 * newldhU          The leading dimension of the output hU
 *
 * hVecs            The eigenvectors of H
 *
 * ldhVecs          The leading dimension of the input hVecs
 *
 * newldhVecs       The leading dimension of the output hVecs
 *
 * hVals            The eigenvalues of H
 *
 * hSVals           The singular values of R
 *
 * hVecsPerm        The permutation that orders hVals and hVecs as primme.target
 *
 * restartPerm      The permutation applied to the columns of hVecs before restarting
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 *
 * numArbitraryVecs On input, the number of coefficients vectors that do
 *                  not correspond to solutions of the projected problem.
 *                  They appear in the first numArbitraryVecs positions of hVecs.
 *                  On output, the number of such vectors that are in the
 *                  restarted basis.
 * 
 * Return value
 * ------------
 * Error code: 0 upon success
 *            -1 eigenvalues of submatrix could not be computed
 *
 ******************************************************************************/

static int restart_harmonic(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *H, int ldH, SCALAR *Q, PRIMME_INT ldQ,
      PRIMME_INT nLocal, SCALAR *R, int ldR, SCALAR *QtV, int ldQtV, SCALAR *hU,
      int ldhU, int newldhU, SCALAR *hVecs, int ldhVecs, int newldhVecs,
      REAL *hVals, REAL *hSVals, int *restartPerm, int *hVecsPerm,
      int restartSize, int basisSize, int numPrevRetained, 
      int indexOfPreviousVecs, int *targetShiftIndex, int numConverged,
      int *numArbitraryVecs, SCALAR *hVecsRot, int ldhVecsRot,
      size_t *rworkSize, SCALAR *rwork, int iworkSize, int *iwork,
      double machEps, primme_params *primme) {

   (void)ldhU; /* unused parameter */
   (void)restartPerm; /* unused parameter */
   (void)hVecsPerm; /* unused parameter */
   (void)numPrevRetained; /* unused parameter */
   (void)indexOfPreviousVecs; /* unused parameter */
   (void)hVecsRot; /* unused parameter */
   (void)ldhVecsRot; /* unused parameter */

   /* Return memory requirement */
 
   if (V == NULL) {
      CHKERR(compute_submatrix_Sprimme(NULL, basisSize, 0, NULL, basisSize,
               0, NULL, 0, NULL, rworkSize), -1);
      CHKERR(update_Q_Sprimme(NULL, nLocal, 0, NULL, 0, NULL, 0, NULL, 0,
               0.0, 0, basisSize, NULL, rworkSize, 0.0, primme), -1);
      CHKERR(update_projection_Sprimme(NULL, 0, NULL, 0, NULL, 0, nLocal,
               0, basisSize, NULL, rworkSize, 0/*unsymmetric*/, primme), -1);
      CHKERR(solve_H_Sprimme(NULL, basisSize, 0, NULL, 0, NULL, 0, NULL, 0,
               NULL, 0, NULL, 0, NULL, NULL, numConverged, 0.0, rworkSize,
               NULL, 0, iwork, primme), -1);
      return 0;
   }

   /* ------------------------------- */
   /* Replace H by hVecs' * H * hVecs */
   /* ------------------------------- */

   CHKERR(compute_submatrix_Sprimme(hVecs, restartSize, ldhVecs, H,
            basisSize, ldH, H, ldH, rwork, rworkSize), -1);

   /* ------------------------------- */
   /* Update targetShiftIndex         */
   /* ------------------------------- */

   *targetShiftIndex = min(primme->numTargetShifts-1, numConverged);

   /* ------------------------------- */
   /* Compute QR                      */
   /* ------------------------------- */

   CHKERR(update_Q_Sprimme(V, nLocal, ldV, W, ldW, Q, ldQ, R, ldR,
         primme->targetShifts[*targetShiftIndex], 0,
         restartSize, rwork, rworkSize, machEps, primme), -1);

   /* ------------------------------- */
   /* Update QtV                      */
   /* ------------------------------- */

   CHKERR(update_projection_Sprimme(Q, ldQ, V, ldV, QtV, ldQtV, nLocal, 0,
            restartSize, rwork, rworkSize, 0/*unsymmetric*/, primme), -1);

   /* ------------------------------- */
   /* Solve the projected problem     */
   /* ------------------------------- */

   CHKERR(solve_H_Sprimme(H, restartSize, ldH, NULL, 0, R, ldR, QtV, ldQtV, hU,
            newldhU, hVecs, newldhVecs, hVals, hSVals, numConverged, machEps,
            rworkSize, rwork, iworkSize, iwork, primme), -1);

   *numArbitraryVecs = 0;

   return 0;
}


/*******************************************************************************
 * Function dtr - This function determines the number of coefficient vectors
 *    to retain from both the left and right side of the spectrum.  The vectors
 *    are then copied so that they are contiguous in memory.
 *
 * Input parameters
 * ----------------
 * numLocked  The number of Ritz vectors that have been locked
 * 
 * basisSize  The current size of the basis
 *
 * numFree    Number of vacancies to be left in the basis
 *
 * iev        Array of size blockSize that determines index Ritz value index 
 *            each block vector is associated with.       
 *
 * rwork      REAL work array of size maxBasisSize^2. Use for hVec swapping.
 *
 * primme       Structure containing various solver parameters       
 *
 *
 * Input/output parameters
 * -----------------------
 * hVecs      The eigenvectors (coefficient vectors) of the projection V'*H*V
 *
 * hVals      The eigenvalues (Ritz values) of H
 *
 * flags      Array of size basisSize indicating the convergence of the Ritz
 *            vectors
 * 
 * Return value
 * ------------
 * int  The new restart size
 *
 ******************************************************************************/


static int dtr_Sprimme(int numLocked, SCALAR *hVecs, REAL *hVals, int *flags, 
  int basisSize, int numFree, int *iev, SCALAR *rwork, primme_params *primme)
{

   int i;                 /* Loop variable */
   int l, lOpt, lMin;     /* Determine how many left side vectors to retain   */
   int r, rOpt;           /* Determine how many right side vectors to retain  */
   int maxIndex;          /* basisSize - 1                                    */
   int restartSize;       /* The new restart size                             */
   REAL currentRitzVal; /* The current Ritz value the solver is computing   */
   REAL newVal, optVal; /* Used to find the optimum gap ratio               */

   /* ---------------------------------------------------------------- */
   /* Compute lOpt and rOpt with respect to the first Ritz value being */
   /* targeted by the block.                                           */
   /* ---------------------------------------------------------------- */

   currentRitzVal = hVals[iev[0]];
   maxIndex = basisSize-1;

   /* If locking is engaged, then lMin must be large enough to retain */
   /* the coefficient vector associated with a converged target.      */
   /* lMin should be no smaller than primme->minRestartSize.            */

   if (primme->locking) {

      lMin = 0;

      /* Determine the largest index of any converged but unlocked target */
      /* Ritz vector.                                                     */

      for (l = 0; l < basisSize; l++) {
         if (flags[l] != UNCONVERGED && numLocked + l < primme->numEvals) {
            lMin = l;
         }
      }

      lMin = max(lMin, min(basisSize, primme->minRestartSize));

   }
   else {
      lMin = min(basisSize, primme->minRestartSize);
   }

   
   lOpt = lMin;
   rOpt = 0;   
   optVal = 0.0L;

   if (primme->printLevel >= 5 && primme->procID == 0) {
      fprintf(primme->outputFile,"DTR basisSize: %d\n", basisSize);
   }

   /* ---------------------------------------------------------------------- */
   /* Compute lOpt and rOpt that maximize the function.                      */
   /* maximize the function (basisSize-numFree-lMin-rMin)*                   */
   /*                       sqrt((currentRitzVal - hVals[l+1])/              */
   /*                            (hVals[l+1]-hVals[basisSize-1-r]))          */
   /* ---------------------------------------------------------------------- */

   for (l = lMin; l < basisSize - numFree; l++) {
      for (r = 0; r < basisSize - l - numFree; r++)       {
         if ((basisSize - l - r) % primme->maxBlockSize == 0) {
            newVal = (basisSize - l - r)
                     * sqrt((currentRitzVal - hVals[l+1])/
                            (hVals[l+1]-hVals[maxIndex-r]));

            if (newVal > optVal) {
               optVal = newVal;
               lOpt = l;
               rOpt = r;
            }

         }
      }
   }


   restartSize = lOpt + rOpt;

   /* --------------------------------------------------------------- */
   /* Swap the rOpt vectors from the right hand side so that they are */
   /* contiguous with the vectors from the left hand side.            */
   /* --------------------------------------------------------------- */

   i = basisSize - restartSize; 

   Num_copy_Sprimme(i*basisSize, &hVecs[basisSize*lOpt], 1, rwork, 1);
   Num_copy_Sprimme(rOpt*basisSize, &hVecs[basisSize*(basisSize-rOpt)], 1,
      &hVecs[basisSize*lOpt], 1);
   Num_copy_Sprimme(i*basisSize, rwork, 1, &hVecs[basisSize*restartSize], 1);

   /* Do the same with the eigenvalues of H */

   Num_copy_Rprimme(i, &hVals[lOpt], 1, (REAL *) rwork, 1);
   Num_copy_Rprimme(rOpt, &hVals[(basisSize-rOpt)], 1, &hVals[lOpt], 1);
   Num_copy_Rprimme(i, (REAL *) rwork, 1, &hVals[restartSize], 1);

   /* Set only those flags lower than restartSize. The rest will be reset */
   for (i = 0; i < rOpt; i++) {
      flags[lOpt + i] = flags[basisSize-rOpt + i];
   }

   if (primme->printLevel >= 5 && primme->procID == 0) {
      fprintf(primme->outputFile,"DTR restart size: %d L: %d R: %d\n", 
         restartSize, lOpt, rOpt);
   }

   for (i=restartSize; i<primme->maxBasisSize; i++) {
      flags[i] = UNCONVERGED;
   }

   return restartSize;

}

/*******************************************************************************
 * Subroutine ortho_coefficient_vectors - Orthogonalize properly the columns of
 *    hVecs from indexOfPreviousVecs plus numPrevRetained.
 *    
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * hVecs            The eigenvectors of H
 * ldhVecs          The leading dimension of the input hVecs
 * basisSize        Rows in hVecs and size of R and hU
 * hU               The eigenvectors of QtV/R
 * ldhU             The leading dimension of hU
 * R                The factors of the QR decomposition of (A - targetShift*B)*V
 * ldR              The leading dimension of R
 * outR             Rotations of orthogonalizing hVecs
 * iwork            Integer work array
 * rwork            Work array
 * rworkSize        Length of the work array
 *
 ******************************************************************************/

static int ortho_coefficient_vectors_Sprimme(SCALAR *hVecs, int basisSize,
      int ldhVecs, int indexOfPreviousVecs, SCALAR *hU, int ldhU, SCALAR *R,
      int ldR, SCALAR *VtBV, int ldVtBV, int *numPrevRetained, double machEps,
      SCALAR *rwork, size_t *rworkSize, primme_params *primme) {

   int i;
   SCALAR *rwork0;
   int newNumPrevRetained=0;

   /* Return memory requirement */

   if (!hVecs) {
      size_t rworkSize0 = 0;
      CHKERR(ortho_Sprimme(NULL, ldhVecs, NULL, 0, 0, basisSize, NULL, 0, 0,
               basisSize, NULL, 0.0, rwork, &rworkSize0, NULL), -1);
      CHKERR(Bortho_local_Sprimme(&hVecs[ldhVecs*indexOfPreviousVecs], ldhVecs,
               NULL, *numPrevRetained, 0, *numPrevRetained-1,
               hVecs, ldhVecs, indexOfPreviousVecs, basisSize, VtBV, ldVtBV,
               primme->iseed, machEps, rwork0, &rworkSize0, NULL), -1);
      /* for dummyR and broadcasting previous retained vectors in hVecs */
      rworkSize0 += (size_t)basisSize*(*numPrevRetained)*2+2;
      *rworkSize = max(*rworkSize, rworkSize0);
      return 0;
   }

   if (primme->projectionParams.projection == primme_proj_harmonic) {

      /* TODO: pending to explain this, see solve_H_Harm for some guidance */

      CHKERR(ortho_Sprimme(hVecs?&hVecs[ldhVecs*indexOfPreviousVecs]:NULL,
               ldhVecs, NULL, 0, 0, *numPrevRetained-1,
               &hU[ldhU*indexOfPreviousVecs], ldhU, indexOfPreviousVecs,
               basisSize, primme->iseed, machEps, rwork, rworkSize, NULL), -1);
      if (hVecs) Num_trsm_Sprimme("L", "U", "N", "N", basisSize,
            *numPrevRetained, 1.0, R, ldR, &hVecs[ldhVecs*indexOfPreviousVecs],
            ldhVecs);

   }

   if (primme->procID == 0) {
      /* Avoid that ortho replaces linear dependent vectors by random vectors.*/
      /* The random vectors make the restarting W with large singular value.  */
      /* This change has shown benefit finding the smallest singular values.  */

      SCALAR *dummyR = rwork;
      rwork0 = rwork + (*numPrevRetained)*(*numPrevRetained);
      assert(*rworkSize >= (size_t)(*numPrevRetained)*(*numPrevRetained));
      size_t rworkSize0 = *rworkSize - (size_t)(*numPrevRetained)*(*numPrevRetained);

      CHKERR(Bortho_local_Sprimme(&hVecs[ldhVecs*indexOfPreviousVecs], ldhVecs,
               dummyR, *numPrevRetained, 0, *numPrevRetained-1,
               hVecs, ldhVecs, indexOfPreviousVecs, basisSize, VtBV, ldVtBV,
               primme->iseed, machEps, rwork0, &rworkSize0, NULL), -1);

      for (i=0; i<*numPrevRetained; i++) {
         if (REAL_PART(dummyR[(*numPrevRetained)*i+i]) == 0.0) {
            break;
         }
      }
      newNumPrevRetained = i;
   }

   /* Broadcast hVecs(indexOfPreviousVecs:indexOfPreviousVecs+numPrevRetained) */

   if (primme->procID == 0) {
      rwork[0] = (SCALAR)newNumPrevRetained;
      Num_copy_matrix_Sprimme(&hVecs[ldhVecs*indexOfPreviousVecs], basisSize,
            *numPrevRetained, ldhVecs, rwork+1, basisSize);
   }
   else {
      rwork[0] = 0.0;
      Num_zero_matrix_Sprimme(rwork+1, basisSize, *numPrevRetained, basisSize);
   }

   rwork0 = rwork + basisSize*(*numPrevRetained) + 1;
   assert(*rworkSize >= 2u*basisSize*(*numPrevRetained) + 2);
   CHKERR(globalSum_Sprimme(rwork, rwork0, basisSize*(*numPrevRetained)+1,
            primme), -1);
   *numPrevRetained = (int)REAL_PART(rwork0[0]);
   Num_copy_matrix_Sprimme(rwork0+1, basisSize, *numPrevRetained, basisSize,
         &hVecs[ldhVecs*indexOfPreviousVecs], ldhVecs);

  return 0;
}


/*******************************************************************************
 * Function retain_previous_coefficients - This function is part of the
 *    recurrence-based restarting method for accelerating convergence.
 *    This function is called one iteration before restart so that the
 *    coefficients (eigenvectors of the projection H) corresponding to 
 *    a few of the target Ritz vectors may be retained at restart. 
 *    The desired coefficients are copied to a separate storage space so
 *    that they may be preserved until the restarting routine is called.
 *
 *
 * Input parameters
 * ----------------
 * hVecs            The coefficients (eigenvectors of the projection H).
 *
 * ldhVecs          The leading dimension of the input hVecs
 *
 * hV               The right singular vectors or R
 *
 * ldhV             The leading dimension of the input hV
 *
 * hU               The left singular vectors of R or the eigenvectors of QtV/R
 *
 * ldhU             The leading dimension of the input hU
 *
 * hSVals           The singular values of R
 *
 * mprevious        The number of rows of previousHVecs
 *
 * basisSize        The current size of the basis
 *
 * iev              Array of size block size.  It maps the block index to the Ritz
 *                  value index each block vector corresponds to.
 *
 * blockSize        The number of block vectors generated during the current iteration
 *
 * primme           Structure containing various solver parameters
 *
 *
 * Output parameters
 * -----------------
 * previousHVecs    Coefficient vectors retained
 *
 * ldpreviousHVecs  The leading dimension of previousHVecs
 *
 * numPrevRetained  The number of vectors retained
 *
 * prevhSVals       The retained singular values
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int retain_previous_coefficients_Sprimme(SCALAR *hVecs, int ldhVecs,
   SCALAR *hU, int ldhU, SCALAR *previousHVecs, int ldpreviousHVecs,
   int mprevious, int basisSize, int *iev, int blockSize, int *flags,
   int *numPrevRetained, int *iwork, int iworkSize, primme_params *primme) {

   int i;               /* Loop indices                                  */
   int index;           /* The index of some coefficient vector in hVecs */ 
   int *cols;           /* Indices of the column retained                */

   /* Return memory requirement */

   if (hVecs == NULL) {
      *iwork = max(*iwork, primme->restartingParams.maxPrevRetain);
      return 0;
   }

   assert(iworkSize >= primme->restartingParams.maxPrevRetain);
   cols = iwork;

   /* First, retain coefficient vectors corresponding to current block */
   /* vectors.  If all of those have been retained, then retain the    */ 
   /* the next unconverged coefficient vectors beyond iev[blockSize-1].*/

   *numPrevRetained = 0;
   for (i=0, index=0; i < primme->restartingParams.maxPrevRetain
         && index < basisSize && i <= blockSize; index++) {

      if (i < blockSize) {
         index = cols[i] = iev[i];
         i++;
      }
      else if (flags[index] == UNCONVERGED) {
         cols[i] = index;
         i++;
      }
   }
   *numPrevRetained = i;

   if (primme->projectionParams.projection == primme_proj_RR
         || primme->projectionParams.projection == primme_proj_refined) {
      Num_copy_matrix_columns_Sprimme(hVecs, basisSize, cols,
            *numPrevRetained, ldhVecs, previousHVecs, NULL, ldpreviousHVecs);

      /* Zero the maxBasisSize-basisSize last elements of the matrix */

      Num_zero_matrix_Sprimme(previousHVecs+basisSize, mprevious-basisSize,
            *numPrevRetained, ldpreviousHVecs);
   }
   else if (primme->projectionParams.projection == primme_proj_harmonic) {
      /* If harmonic, the coefficient vectors (i.e., the eigenvectors of the  */
      /* projected problem) are in hU; so retain them.                        */

      Num_copy_matrix_columns_Sprimme(hU, basisSize, cols,
            *numPrevRetained, ldhU, previousHVecs, NULL, ldpreviousHVecs);

      /* Zero the maxBasisSize-basisSize last elements of the matrix */

      Num_zero_matrix_Sprimme(previousHVecs+basisSize, mprevious-basisSize,
            *numPrevRetained, ldpreviousHVecs);
   }

   if (primme->printLevel >= 5 && primme->procID == 0) {
      fprintf(primme->outputFile, "retain_previous: numPrevRetained: %d\n",
            *numPrevRetained);
   }

   return 0;
}

/******************************************************************************
 * Subroutine insertionSort -- This subroutine locks a converged Ritz value
 *   by insertion sorting it into the evals array.  A permutation array, perm,
 *   is maintained to keep track of the position the value would have been
 *   placed in had sorting not been performed.  This allows the locked Ritz
 *   vectors to be sorted at a later time upon return to the user.
 *   The order is ascending or descending for smallest/largest respectively.
 *   For interior, it is the order of convergence except for the same shifts.
 *   In that case, Ritz values that satisfy the criterion closer come first.
 *
 *
 * Input parameters
 * ----------------
 * newVal   The Ritz value to be locked
 *
 * newNorm  The residual norm of the Ritz value to be locked
 *
 * newFlag     The current flag
 *
 * numLocked  The current number of locked Ritz vectors
 * 
 * primme  Structure containing various solver parameters
 *
 *
 * Input/Output parameters
 * -----------------------
 * evals    The sorted list of locked Ritz values
 *
 * resNorms The residual norms corresponding to the locked Ritz values
 *
 * flags    The flags of the locked pairs
 *
 * perm     The permutation array indicating each Ritz values original
 *          unsorted position.
 *
 ******************************************************************************/

static void insertionSort(REAL newVal, REAL *evals, REAL newNorm,
   REAL *resNorms, int newFlag, int *flags, int *perm, int numLocked,
   primme_params *primme) {

   int i, current; /* Indices used for sorting */
   REAL ithShift, currentShift;

   /* ------------------------------------------------------------------ */
   /* Find smallest index to insert the Ritz value. The eigenvalue order */
   /* depends on how we target eigenvalues.                              */
   /* ------------------------------------------------------------------ */

   if ( primme->target == primme_smallest ) {

      for (i = numLocked; i > 0; i--) {
         if (newVal >= evals[i-1]) break; 
      }
   }
   else if ( primme->target == primme_largest ) {

      for (i = numLocked; i > 0; i--) {
         if (newVal <= evals[i-1]) break; 
      }
   }
   else {
   /* For interior cases maintain convergence order except for the same shift *
    * Example: s1 s2 s2 s3 s4. Only eigenvalues converged for s2 may switch.  *
    * Therefore, we only need to look back as long as the shift is the same.  */

      currentShift =
        primme->targetShifts[min(primme->numTargetShifts-1, numLocked)];

      if ( primme->target == primme_closest_geq ) {
         for (i = numLocked; i > 0; i--) {
            ithShift =primme->targetShifts[min(primme->numTargetShifts-1, i-1)];
            if ( ithShift != currentShift || 
            newVal-currentShift >= evals[i-1]-currentShift ) break; 
         }
      }
      else if ( primme->target == primme_closest_leq ) {
         for (i = numLocked; i > 0; i--) {
            ithShift =primme->targetShifts[min(primme->numTargetShifts-1, i-1)];
            if ( ithShift != currentShift || 
            currentShift-newVal >= currentShift-evals[i-1] ) break; 
         }
      }
      else if ( primme->target == primme_closest_abs ) {
         for (i = numLocked; i > 0; i--) {
            ithShift =primme->targetShifts[min(primme->numTargetShifts-1, i-1)];
            if ( ithShift != currentShift || 
            fabs(newVal-currentShift) >= fabs(evals[i-1]-currentShift) ) break; 
         }
      }
      else if ( primme->target == primme_largest_abs ) {
         for (i = numLocked; i > 0; i--) {
            ithShift =primme->targetShifts[min(primme->numTargetShifts-1, i-1)];
            if ( ithShift != currentShift || 
            fabs(newVal-currentShift) <= fabs(evals[i-1]-currentShift) ) break; 
         }
      }
      else {
         /* This should never happen */
         assert(0);
         i = 0; /* Avoid warning */
      }
   }

   /* Shift the array to make room for the new Ritz value */

   for (current = numLocked-1; current >= i; current--) {
      evals[current+1] = evals[current];
      resNorms[current+1] = resNorms[current];
      perm[current+1] = perm[current];
      flags[current+1] = flags[current];
   }

   /* Insert the new value */

   evals[i] = newVal;
   resNorms[i] = newNorm;
   perm[i] = numLocked;
   flags[i] = newFlag;
   return;

}

/******************************************************************************
 * Function compute_residual_columns - This subroutine performs the next
 *    operations in a cache-friendly way:
 *
 *    X = X(p); Ax = Ax(p)
 *    j = k = 0; XD = RD = []
 *    for i=0:nd-1
 *       if pd(i) == j+io0
 *          XD = [XD XO(j)]; RD = [RD RO(j)]; j++
 *       else
 *          XD = [XD X(p(k)]; RD = [RD AX(p(k)) - evals(p(k))*X(p(k))]; k++
 *       end if
 *    end for
 *
 *           n        nd             no
 * X:  [-----|-------------]           (input/output)
 * XO:                      [--------| (input)
 * XD:       [--------|                (output)
 *
 * NOTE: X and XD *can* overlap, but X(0:n-1) and XD *cannot* overlap (same for R and RD)
 *       XO and XD *can* overlap (same for RO and RD)
 *       p should be a list of increasing indices
 *       pd should be a merge of two increasing lists
 *
 * PARAMETERS
 * ---------------------------
 * m           The number of rows of matrices x, Ax, xo, ro, xd and rd
 * evals       The values to compute the residual
 * x           The matrix that does x = x(p)
 * n           The number of columns of the output x
 * p           The columns to copy back to x and Ax
 * ldx         The leading dimension of x
 * Ax          The matrix that does Ax = Ax(p)
 * ldAx        The leading dimension of Ax
 * xo          Alternative source of columns for xd
 * no          The number of columns in xo
 * ldxo        The leading dimension of xo
 * io0         The index of the first column in xo
 * ro          Alternative source of columns for rd
 * ldro        The leading dimension of ro
 * xd          The matrix that will have columns from x and xo
 * nd          The maximum size of xd
 * pd          The indices of the columns to generate xd and rd
 * ldxd        The leading dimension of xd
 * rd          The matrix that will have columns from r and ro
 * ldrd        The leading dimension of rd
 * rwork       Workspace
 * lrwork      The size of rwork
 *
 ******************************************************************************/

static int compute_residual_columns(PRIMME_INT m, REAL *evals, SCALAR *x,
      int n, int *p, PRIMME_INT ldx, SCALAR *Ax, PRIMME_INT ldAx,
      SCALAR *xo, int no, PRIMME_INT ldxo, int io0, SCALAR *ro, PRIMME_INT ldro,
      SCALAR *xd, int nd, int *pd, PRIMME_INT ldxd, SCALAR *rd, PRIMME_INT ldrd,
      SCALAR *rwork, PRIMME_INT lrwork) {

   int i, id, k, io, M=min(m,PRIMME_BLOCK_SIZE);
   SCALAR *X0, *R0;

   /* Return memory requirement */

   if (evals == NULL) {
      return nd*M*2;
   }

   /* Quick exit */

   if (n == 0) {
      Num_copy_matrix_Sprimme(xo, m, min(no,nd), ldxo, xd, ldxd);
      Num_copy_matrix_Sprimme(ro, m, min(no,nd), ldro, rd, ldrd);
      return 0;
   }

   X0 = rwork;
   R0 = X0+nd*M;
   assert(nd*M*2 <= lrwork);

   for (k=0; k<m; k+=M, M=min(M,m-k)) {
      for (i=id=io=0; i < n || id < nd; id++) {
         if (id < nd && io < no && pd[id] == io+io0) {
            Num_copy_matrix_Sprimme(&xo[io*ldxo+k], M, 1, ldxo, &X0[id*M], M);
            Num_copy_matrix_Sprimme(&ro[io*ldro+k], M, 1, ldro, &R0[id*M], M);
            io++;
         }
         else {
            assert(id >= nd || i < n);
            Num_copy_matrix_Sprimme(&x[p[i]*ldx+k],   M, 1, ldx,  &x[i*ldx +k],  ldx);
            Num_copy_matrix_Sprimme(&Ax[p[i]*ldAx+k], M, 1, ldAx, &Ax[i*ldAx+k], ldAx);
            if (id < nd) {
               Num_copy_matrix_Sprimme(&x[p[i]*ldx+k], M, 1, ldx, &X0[id*M], M);
               Num_compute_residual_Sprimme(M, evals[p[i]], &x[p[i]*ldx+k], &Ax[p[i]*ldAx+k], &R0[id*M]);
            }
            i++;
         }
      }
      assert(id >= nd);
      Num_copy_matrix_Sprimme(X0, M, nd, M, &xd[k], ldxd);
      Num_copy_matrix_Sprimme(R0, M, nd, M, &rd[k], ldrd);
   }

   return 0;
}
