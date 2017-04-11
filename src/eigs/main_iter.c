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
 * File: main_iter.c
 *
 * Purpose - This is the main Davidson-type iteration 
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "const.h"
#include "wtime.h"
#include "numerical.h"
#include "main_iter.h"
#include "main_iter_private.h"
#include "convergence.h"
#include "correction.h"
#include "init.h"
#include "ortho.h"
#include "restart.h"
#include "solve_projection.h"
#include "update_projection.h"
#include "update_W.h"
#include "globalsum.h"
#include "auxiliary_eigs.h"

static int verify_norms(SCALAR *V, PRIMME_INT ldV, SCALAR *W, PRIMME_INT ldW,
      REAL *hVals, int basisSize, REAL *resNorms, int *flags, int *converged,
      double machEps, SCALAR *rwork, size_t *rworkSize, int *iwork,
      int iworkSize, primme_params *primme);

/******************************************************************************
 * Subroutine main_iter - This routine implements a more general, parallel, 
 *    block (Jacobi)-Davidson outer iteration with a variety of options.
 *
 * A coarse outline of the algorithm performed is as follows:
 *
 *  1. Initialize basis V
 *  2. while (converged Ritz vectors have become unconverged) do
 *  3.    Update W=A*V, and H=V'*A*V, and solve the H eigenproblem
 *  4.    while (not all Ritz pairs have converged) do
 *  5.       while (maxBasisSize has not been reached) 
 *  6.          Adjust the block size if necessary
 *  7.          Compute the residual norms of each block vector and 
 *                check for convergence.  If a Ritz vector has converged,
 *                target an unconverged Ritz vector to replace it.
 *  8.          Solve the correction equation for each residual vector
 *  9.          Insert corrections in V, orthogonalize, update W and H 
 * 10.          Solve the new H eigenproblem to obtain the next Ritz pairs    
 * 11.       endwhile
 * 12.       Restart V with appropriate vectors 
 * 13.       If (locking) Lock out Ritz pairs that converged since last restart
 * 14.    endwhile
 * 15.    if Ritz vectors have become unconverged, reset convergence flags
 * 16. endwhile
 *    
 *
 * INPUT only parameters
 * ---------------------
 * machEps machine precision 
 *
 * OUTPUT arrays and parameters
 * ----------------------------
 * evals    The converged Ritz values.  It is maintained in sorted order.
 *          If locking is not engaged, values are copied to this array only 
 *          upon return.  If locking is engaged, then values are copied to this 
 *          array as they converge.
 *
 * perm     A permutation array that maps the converged Ritz vectors to
 *          their corresponding Ritz values in evals.  
 *
 * resNorms The residual norms corresponding to the converged Ritz vectors
 *
 * intWork  Integer work array
 *
 * realWork SCALAR work array
 *
 * INPUT/OUTPUT arrays and parameters
 * ----------------------------------
 * evecs    Stores initial guesses. Upon return, it contains the converged Ritz
 *          vectors.  If locking is not engaged, then converged Ritz vectors 
 *          are copied to this array just before return.  
 *
 * primme.initSize: On output, it stores the number of converged eigenvectors. 
 *           If smaller than numEvals and locking is used, there are
 *              only primme.initSize vectors in evecs.
 *           Without locking all numEvals approximations are in evecs
 *              but only the initSize ones are converged.
 *           During the execution, access to primme.initSize gives 
 *              the number of converged pairs up to that point. The pairs
 *              are available in evals and evecs, but only when locking is used
 *
 * Return Value
 * ------------
 * int -  0 the requested number of Ritz values converged
 *       -1 if solver did not converge within required number of iterations
 *       -2 if initialization failed
 *       -3 if orthogonalization failed
 *       -4 if solve_H failed
 *       -5 if solve_correction failed
 *       -6 if restart failed
 *       -7 if lock_vectors failed
 *       
 ******************************************************************************/

TEMPLATE_PLEASE
int main_iter_Sprimme(REAL *evals, int *perm, SCALAR *evecs, PRIMME_INT ldevecs,
   REAL *resNorms, double machEps, int *intWork, void *realWork, 
   primme_params *primme) {
         
   int i;                   /* Loop variable                                 */
   int blockSize;           /* Current block size                            */
   int availableBlockSize;  /* There is enough space in basis for this block */
   int basisSize;           /* Current size of the basis V                   */
   int numLocked;           /* Number of locked Ritz vectors                 */
   int numGuesses;          /* Current number of initial guesses in evecs    */
   int nextGuess;           /* Index of the next initial guess in evecs      */ 
   int numConverged;        /* Number of converged Ritz pairs                */
   int targetShiftIndex;    /* Target shift used in the QR factorization     */
   int recentlyConverged;   /* Number of target Ritz pairs that have         */
                            /*    converged during the current iteration.    */
   int maxRecentlyConverged;/* Maximum value allowed for recentlyConverged   */
   int numConvergedStored;  /* Numb of Ritzvecs temporarily stored in evecs  */
                            /*    to allow for skew projectors w/o locking   */
   int converged;           /* True when all required Ritz vals. converged   */
   int LockingProblem;      /* Flag==1 if practically converged pairs locked */
   int restartLimitReached; /* True when maximum restarts performed          */
   int numPrevRetained;     /* Number of vectors retained using recurrence-  */
                            /* based restarting.                             */
   int numArbitraryVecs;    /* Columns in hVecs computed with RR instead of  */
                            /* the current extraction method.                */
   int maxEvecsSize;        /* Maximum capacity of evecs array               */
   size_t rworkSize;        /* Size of rwork array                           */
   int iworkSize;           /* Size of iwork array                           */
   int numPrevRitzVals = 0; /* Size of the prevRitzVals updated in correction*/
   int ret;                 /* Return value                                  */
   int touch=0;             /* param used in inner solver stopping criteria  */

   int *iwork;              /* Integer workspace pointer                     */
   int *flags;              /* Indicates which Ritz values have converged    */
   int *lockedFlags=NULL;   /* Flags for the locked pairs                    */
   int *ipivot;             /* The pivot for the UDU factorization of M      */
   int *iev;                /* Evalue index each block vector corresponds to */

   SCALAR *V;               /* Basis vectors                                 */
   PRIMME_INT ldV;          /* The leading dimension of V                    */
   SCALAR *W;               /* Work space storing A*V                        */
   PRIMME_INT ldW;          /* The leading dimension of W                    */
   SCALAR *H;               /* Upper triangular portion of V'*A*V            */
   SCALAR *M = NULL;        /* The projection Q'*K*Q, where Q = [evecs, x]   */
                            /* x is the current Ritz vector and K is a       */
                            /* hermitian preconditioner.                     */
   SCALAR *UDU = NULL;      /* The factorization of M=Q'KQ                   */
   SCALAR *evecsHat = NULL; /* K^{-1}evecs                                   */
   PRIMME_INT ldevecsHat=0; /* The leading dimension of evecsHat             */
   SCALAR *rwork;           /* Real work space.                              */
   SCALAR *hVecs;           /* Eigenvectors of H                             */
   SCALAR *hU=NULL;         /* Left singular vectors of R                    */
   SCALAR *previousHVecs;   /* Coefficient vectors retained by               */
                            /* recurrence-based restarting                   */

   int numQR;               /* Maximum number of QR factorizations           */
   SCALAR *Q = NULL;        /* QR decompositions for harmonic or refined     */
   PRIMME_INT ldQ;          /* The leading dimension of Q                    */
   SCALAR *R = NULL;        /* projection: (A-target[i])*V = QR              */
   SCALAR *QtV = NULL;      /* Q'*V                                          */
   SCALAR *hVecsRot=NULL;   /* transformation of hVecs in arbitrary vectors  */

   REAL *hVals;             /* Eigenvalues of H                              */
   REAL *hSVals=NULL;       /* Singular values of R                          */
   REAL *prevRitzVals;      /* Eigenvalues of H at previous outer iteration  */
                            /* by robust shifting algorithm in correction.c  */
   REAL *basisNorms;        /* Residual norms of basis at pairs              */
   REAL *blockNorms;        /* Residual norms corresponding to current block */
                            /* vectors.                                      */
   double smallestResNorm;  /* the smallest residual norm in the block       */
   int reset=0;             /* Flag to reset V and W                         */
   int restartsSinceReset=0;/* Restart since last reset of V and W           */
   int wholeSpace=0;        /* search subspace reach max size                */

   /* Runtime measurement variables for dynamic method switching             */
   primme_CostModel CostModel; /* Structure holding the runtime estimates of */
                            /* the parameters of the model.Only visible here */
   double tstart=0.0;       /* Timing variable for accumulative time spent   */

   /* -------------------------------------------------------------- */
   /* Subdivide the workspace                                        */
   /* -------------------------------------------------------------- */

   maxEvecsSize = primme->numOrthoConst + primme->numEvals;
   if (primme->projectionParams.projection != primme_proj_RR) {
      numQR = 1;
   }
   else {
      numQR = 0;
   }

   /* Use leading dimension ldOPs for the large dimension mats: V, W and Q */

   ldV = ldW = ldQ = primme->ldOPs;
   rwork         = (SCALAR *) realWork;
   V             = rwork; rwork += primme->ldOPs*primme->maxBasisSize;
   W             = rwork; rwork += primme->ldOPs*primme->maxBasisSize;
   if (numQR > 0) {
      Q          = rwork; rwork += primme->ldOPs*primme->maxBasisSize*numQR;
      R          = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize*numQR;
      hU         = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize*numQR;
   }
   if (primme->projectionParams.projection == primme_proj_harmonic) {
      QtV        = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize*numQR;
   }
   H             = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize;
   hVecs         = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize;
   previousHVecs = rwork; rwork += primme->maxBasisSize*primme->restartingParams.maxPrevRetain;
   if (primme->projectionParams.projection == primme_proj_refined
       || primme->projectionParams.projection == primme_proj_harmonic) {
      hVecsRot   = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize*numQR;
   }

   if (primme->correctionParams.precondition && 
         primme->correctionParams.maxInnerIterations != 0 &&
         primme->correctionParams.projectors.RightQ &&
         primme->correctionParams.projectors.SkewQ           ) {
      ldevecsHat = primme->nLocal;
      evecsHat   = rwork; rwork += ldevecsHat*maxEvecsSize;
      M          = rwork; rwork += maxEvecsSize*maxEvecsSize;
      UDU        = rwork; rwork += maxEvecsSize*maxEvecsSize;
   }

   #ifndef USE_COMPLEX
   #  define TO_REAL(X) X
   #else
   #  define TO_REAL(X) (X+1)/2
   #endif
   hVals         = (REAL *)rwork; rwork += TO_REAL(primme->maxBasisSize);
   if (numQR > 0) {
      hSVals     = (REAL *)rwork; rwork += TO_REAL(primme->maxBasisSize);
   }
   prevRitzVals  = (REAL *)rwork; rwork += TO_REAL(primme->maxBasisSize+primme->numEvals);
   blockNorms    = (REAL *)rwork; rwork += TO_REAL(primme->maxBlockSize);
   basisNorms    = (REAL *)rwork; rwork += TO_REAL(primme->maxBasisSize);
   #undef TO_REAL

   rworkSize     = primme->realWorkSize/sizeof(SCALAR) - (rwork - (SCALAR*)realWork);

   /* Integer workspace */

   iwork = intWork; iworkSize = (int)(primme->intWorkSize/sizeof(int));
   if (primme->locking) {
      lockedFlags = iwork; iwork += primme->numEvals; iworkSize -= primme->numEvals;
   }
   flags = iwork; iwork += primme->maxBasisSize; iworkSize -= primme->maxBasisSize;
   iev = iwork; iwork += primme->maxBlockSize; iworkSize -= primme->maxBlockSize;
   ipivot = iwork; iwork += maxEvecsSize; iworkSize -= maxEvecsSize;

   /* -------------------------------------------------------------- */
   /* Initialize counters and flags                                  */
   /* -------------------------------------------------------------- */

   primme->stats.numOuterIterations            = 0; 
   primme->stats.numRestarts                   = 0;
   primme->stats.numMatvecs                    = 0;
   primme->stats.numPreconds                   = 0;
   primme->stats.numGlobalSum                  = 0;
   primme->stats.volumeGlobalSum               = 0;
   primme->stats.numOrthoInnerProds            = 0.0;
   primme->stats.elapsedTime                   = 0.0;
   primme->stats.timeMatvec                    = 0.0;
   primme->stats.timePrecond                   = 0.0;
   primme->stats.timeOrtho                     = 0.0;
   primme->stats.timeGlobalSum                 = 0.0;
   primme->stats.estimateMinEVal               = HUGE_VAL;
   primme->stats.estimateMaxEVal               = -HUGE_VAL;
   primme->stats.estimateLargestSVal           = -HUGE_VAL;
   primme->stats.maxConvTol                    = 0.0;
   primme->stats.estimateResidualError         = 0.0;

   numLocked = 0;
   converged = FALSE;
   LockingProblem = 0;

   numPrevRetained = 0;
   blockSize = 0; 

   for (i=0; i<primme->numEvals; i++) perm[i] = i;

   /* -------------------------------------- */
   /* Quick return for matrix of dimension 1 */
   /* -------------------------------------- */

   if (primme->numEvals == 0) {
      primme->initSize = 0;
      return 0;
   }

   /* -------------------------------------- */
   /* Quick return for matrix of dimension 1 */
   /* -------------------------------------- */

   if (primme->n == 1) {
      evecs[0] = 1.0;
      CHKERR(matrixMatvec_Sprimme(&evecs[0], primme->nLocal, ldevecs,
            W, ldW, 0, 1, primme), -1);
      evals[0] = REAL_PART(W[0]);
      V[0] = 1.0;

      resNorms[0] = 0.0L;
      primme->stats.numMatvecs++;
      primme->initSize = 1;
      return 0;
   }

   /* ------------------------------------------------ */
   /* Especial configuration for matrix of dimension 2 */
   /* ------------------------------------------------ */

   if (primme->n == 2) {
      primme->minRestartSize = 2;
      primme->restartingParams.maxPrevRetain = 0;
   }

   /* -------------------- */
   /* Initialize the basis */
   /* -------------------- */

   CHKERR(init_basis_Sprimme(V, primme->nLocal, ldV, W, ldW, evecs, ldevecs,
            evecsHat, primme->nLocal, M, maxEvecsSize, UDU, 0, ipivot, machEps,
            rwork, &rworkSize, &basisSize, &nextGuess, &numGuesses, primme),
         -1);

   /* Now initSize will store the number of converged pairs */
   primme->initSize = 0;

   /* ----------------------------------------------------------- */
   /* Dynamic method switch means we need to decide whether to    */
   /* allow inner iterations based on runtime timing measurements */
   /* ----------------------------------------------------------- */
   if (primme->dynamicMethodSwitch > 0) {
      initializeModel(&CostModel, primme);
      CostModel.MV = primme->stats.timeMatvec/primme->stats.numMatvecs;
      if (primme->numEvals < 5)
         primme->dynamicMethodSwitch = 1;   /* Start tentatively GD+k */
      else
         primme->dynamicMethodSwitch = 3;   /* Start GD+k for 1st pair */
      primme->correctionParams.maxInnerIterations = 0; 
   }

   /* ---------------------------------------------------------------------- */
   /* Outer most loop                                                        */
   /* Without locking, restarting can cause converged Ritz values to become  */
   /* unconverged. Keep performing JD iterations until they remain converged */
   /* ---------------------------------------------------------------------- */
   while (!converged &&
          ( primme->maxMatvecs == 0 || 
            primme->stats.numMatvecs < primme->maxMatvecs ) &&
          ( primme->maxOuterIterations == 0 ||
            primme->stats.numOuterIterations < primme->maxOuterIterations) ) {

      /* Reset convergence flags. This may only reoccur without locking */

      primme->initSize = numConverged = numConvergedStored = 0;
      for (i=0; i<primme->maxBasisSize; i++)
         flags[i] = UNCONVERGED;

      /* Compute the initial H and solve for its eigenpairs */

      targetShiftIndex = 0;
      if (Q) CHKERR(update_Q_Sprimme(V, primme->nLocal, ldV, W, ldW, Q, ldQ, R,
               primme->maxBasisSize, primme->targetShifts[targetShiftIndex], 0,
               basisSize, rwork, &rworkSize, machEps, primme), -1);

      if (H) CHKERR(update_projection_Sprimme(V, ldV, W, ldW, H,
               primme->maxBasisSize, primme->nLocal, 0, basisSize, rwork,
               &rworkSize, 1/*symmetric*/, primme), -1);

      if (QtV) CHKERR(update_projection_Sprimme(Q, ldQ, V, ldV, QtV,
               primme->maxBasisSize, primme->nLocal, 0, basisSize, rwork,
               &rworkSize, 0/*unsymmetric*/, primme), -1);

      CHKERR(solve_H_Sprimme(H, basisSize, primme->maxBasisSize, R,
               primme->maxBasisSize, QtV, primme->maxBasisSize, hU, basisSize,
               hVecs, basisSize, hVals, hSVals, numConverged, machEps,
               &rworkSize, rwork, iworkSize, iwork, primme), -1);
      
      numArbitraryVecs = 0;
      maxRecentlyConverged = availableBlockSize = blockSize = 0;
      smallestResNorm = HUGE_VAL;

      /* -------------------------------------------------------------- */
      /* Begin the iterative process.  Keep restarting until all of the */
      /* required eigenpairs have been found (no verification)          */
      /* -------------------------------------------------------------- */
      while (numConverged < primme->numEvals &&
             ( primme->maxMatvecs == 0 || 
               primme->stats.numMatvecs < primme->maxMatvecs ) &&
             ( primme->maxOuterIterations == 0 ||
               primme->stats.numOuterIterations < primme->maxOuterIterations) ) {
 
            numPrevRetained = 0;
         /* ----------------------------------------------------------------- */
         /* Main block Davidson loop.                                         */
         /* Keep adding vectors to the basis V until the basis has reached    */
         /* maximum size or the basis plus the locked vectors span the entire */
         /* space. Once this happens, restart with a smaller basis.           */
         /* ----------------------------------------------------------------- */
         while (basisSize < primme->maxBasisSize &&
                basisSize < primme->n - primme->numOrthoConst - numLocked &&
                ( primme->maxMatvecs == 0 || 
                  primme->stats.numMatvecs < primme->maxMatvecs) &&
                ( primme->maxOuterIterations == 0 ||
                  primme->stats.numOuterIterations < primme->maxOuterIterations) ) {

            primme->stats.numOuterIterations++;

            /* When QR are computed and there are more than one target shift, */
            /* limit blockSize and the converged values to one.               */

            if (primme->numTargetShifts > numConverged+1 && Q) {
               availableBlockSize = 1;
               maxRecentlyConverged = numConverged-numLocked+1;
            }
            else {
               availableBlockSize = primme->maxBlockSize;
               maxRecentlyConverged = max(0, primme->numEvals-numConverged);
            }

            /* Limit blockSize to vacant vectors in the basis */

            availableBlockSize = min(availableBlockSize, primme->maxBasisSize-basisSize);

            /* Limit blockSize to remaining values to converge plus one */

            availableBlockSize = min(availableBlockSize, maxRecentlyConverged+1);

            /* Limit basisSize to the matrix dimension */

            availableBlockSize = max(0, min(availableBlockSize, 
                  primme->n - basisSize - numLocked - primme->numOrthoConst));

            /* Set the block with the first unconverged pairs */
            if (availableBlockSize > 0) {
               prepare_candidates_Sprimme(V, ldV, W, ldW, primme->nLocal, H,
                  primme->maxBasisSize, basisSize,
                  &V[basisSize*ldV], &W[basisSize*ldW],
                  hVecs, basisSize, hVals, hSVals, flags,
                  maxRecentlyConverged, blockNorms, blockSize,
                  availableBlockSize, evecs, numLocked, ldevecs, evals,
                  resNorms, targetShiftIndex, machEps, iev, &blockSize,
                  &recentlyConverged, &numArbitraryVecs, &smallestResNorm,
                  hVecsRot, primme->maxBasisSize, numConverged, basisNorms,
                  &reset, rwork, &rworkSize, iwork, iworkSize, primme);
               assert(recentlyConverged >= 0);
            }
            else {
               blockSize = recentlyConverged = 0;
            }

            /* If the total number of converged pairs, including the     */
            /* recentlyConverged ones, are greater than or equal to the  */
            /* target number of eigenvalues, attempt to restart, verify  */
            /* their convergence, lock them if necessary, and return.    */
            /* For locking interior, restart and lock now any converged. */
            /* If Q, restart after an eigenpair converged to recompute   */
            /* QR with a different shift.                                */
            /* Also if it has been converged as many pairs as initial    */
            /* guesses has introduced in V, then restart and introduce   */
            /* new guesses.                                              */

            numConverged += recentlyConverged;

            /* Report iteration */

            if (primme->monitorFun) {
               primme_event EVENT_OUTER_ITERATION = primme_event_outer_iteration;
               primme->stats.elapsedTime = primme_wTimer(0);
               int err;
               CHKERRM((primme->monitorFun(hVals, &basisSize, flags, iev,
                           &blockSize, basisNorms, &numConverged, evals,
                           &numLocked, lockedFlags, resNorms, NULL, NULL,
                           &EVENT_OUTER_ITERATION, primme, &err),
                        err), -1, "Error returned by monitorFun: %d", err);
            }

            /* Reset touch every time an eigenpair converges */

            if (recentlyConverged > 0) touch = 0;

            if (numConverged >= primme->numEvals ||
                (primme->locking && recentlyConverged > 0
                  && primme->target != primme_smallest
                  && primme->target != primme_largest
                  && primme->projectionParams.projection == primme_proj_RR) ||
                targetShiftIndex < 0 ||
                  (blockSize == 0 && recentlyConverged > 0) ||
                /* NOTE: use the same condition as in restart_refined */
                (Q && fabs(primme->targetShifts[targetShiftIndex] -
                  primme->targetShifts[
                     min(primme->numTargetShifts-1, numConverged)]) >= 
                        max(primme->aNorm, primme->stats.estimateLargestSVal))
               || (numConverged >= nextGuess-primme->numOrthoConst
                  && numGuesses > 0)) {

               break;

            }

            /* If the block size is zero, the whole basis spans an exact     */
            /* (converged) eigenspace. Then, since not all needed evecs have */
            /* been found, we must generate a new set of vectors to proceed. */
            /* This set should be of size AvailableBlockSize, and random     */
            /* as there is currently no locking to bring in new guesses.     */
            /* We zero out the V(AvailableBlockSize), avoid any correction   */
            /* and let ortho create the random vectors.                      */

            if (blockSize == 0) {
               blockSize = availableBlockSize;
               Num_scal_Sprimme(blockSize*primme->nLocal, 0.0,
                  &V[ldV*basisSize], 1);
            }
            else {

               /* Solve the correction equations with the new blockSize Ritz */
               /* vectors and residuals.                                     */

               /* - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
               /* If dynamic method switching, time the inner method     */
               if (primme->dynamicMethodSwitch > 0) {
                  tstart = primme_wTimer(0); /* accumulate correction time */

                  if (CostModel.resid_0 == -1.0L)       /* remember the very */
                     CostModel.resid_0 = blockNorms[0]; /* first residual */

                  /*if some pairs converged OR we evaluate jdqmr at every step*/
                  /* update convergence statistics and consider switching */
                  if (recentlyConverged > 0 || primme->dynamicMethodSwitch == 2)
                  {
                     CostModel.MV =
                        primme->stats.timeMatvec/primme->stats.numMatvecs;
                     ret = update_statistics(&CostModel, primme, tstart, 
                        recentlyConverged, 0, numConverged, blockNorms[0], 
                        primme->stats.estimateLargestSVal); 

                     if (ret) switch (primme->dynamicMethodSwitch) {
                        /* for few evals (dyn=1) evaluate GD+k only at restart*/
                        case 3:
                           CHKERR(switch_from_GDpk(&CostModel,primme), -1);
                           break;
                        case 2: case 4:
                           CHKERR(switch_from_JDQMR(&CostModel,primme), -1);
                     } /* of if-switch */
                  } /* of recentlyConv > 0 || dyn==2 */
               } /* dynamic switching */
               /* - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

               CHKERR(solve_correction_Sprimme(V, ldV, W, ldW, evecs, ldevecs,
                        evecsHat, ldevecsHat, UDU, ipivot, evals, numLocked,
                        numConvergedStored, hVals, prevRitzVals,
                        &numPrevRitzVals, flags, basisSize, blockNorms, iev,
                        blockSize, &touch, machEps, rwork, &rworkSize, iwork, iworkSize,
                        primme), -1);

               /* ------------------------------------------------------ */
               /* If dynamic method switch, accumulate inner method time */
               /* ------------------------------------------------------ */
               if (primme->dynamicMethodSwitch > 0) 
                  CostModel.time_in_inner += primme_wTimer(0) - tstart;

              
            } /* end of else blocksize=0 */

            /* Orthogonalize the corrections with respect to each other */
            /* and the current basis.                                   */
            CHKERR(ortho_Sprimme(V, ldV, NULL, 0, basisSize, 
               basisSize+blockSize-1, evecs, ldevecs, 
               primme->numOrthoConst+numLocked, primme->nLocal, primme->iseed, 
               machEps, rwork, &rworkSize, primme), -1);

            /* Compute W = A*V for the orthogonalized corrections */

            CHKERR(matrixMatvec_Sprimme(V, primme->nLocal, ldV, W, ldW,
                     basisSize, blockSize, primme), -1);

            if (Q) CHKERR(update_Q_Sprimme(V, primme->nLocal, ldV, W, ldW, Q,
                     ldQ, R, primme->maxBasisSize,
                     primme->targetShifts[targetShiftIndex], basisSize,
                     blockSize, rwork, &rworkSize, machEps, primme), -1);

            /* Extend H by blockSize columns and rows and solve the */
            /* eigenproblem for the new H.                          */

            if (H) CHKERR(update_projection_Sprimme(V, ldV, W, ldW, H,
                     primme->maxBasisSize, primme->nLocal, basisSize, blockSize,
                     rwork, &rworkSize, 1/*symmetric*/, primme), -1);

            if (QtV) CHKERR(update_projection_Sprimme(Q, ldQ, V, ldV, QtV,
                     primme->maxBasisSize, primme->nLocal, basisSize, blockSize,
                     rwork, &rworkSize, 0/*unsymmetric*/, primme), -1);

            if (basisSize+blockSize >= primme->maxBasisSize) {
               CHKERR(retain_previous_coefficients_Sprimme(hVecs,
                        basisSize, hU, basisSize, previousHVecs,
                        primme->maxBasisSize, primme->maxBasisSize, basisSize,
                        iev, blockSize, flags, &numPrevRetained, iwork,
                        iworkSize, primme), -1);
            }


            basisSize += blockSize;
            blockSize = 0;

            CHKERR(solve_H_Sprimme(H, basisSize, primme->maxBasisSize, R,
                     primme->maxBasisSize, QtV, primme->maxBasisSize, hU,
                     basisSize, hVecs, basisSize, hVals, hSVals, numConverged,
                     machEps, &rworkSize, rwork, iworkSize, iwork, primme), -1);

            numArbitraryVecs = 0;

            /* If QR decomposition accumulates so much error, force it to     */
            /* reset by setting targetShiftIndex to -1. We use the next       */
            /* heuristic. Note that if (s_0, u_0, y_0) is the smallest triplet*/
            /* of R, (A-tau*I)*V = Q*R, and l_0 is the Rayleigh quotient of   */
            /* V*y_0, then                                                    */
            /*    |l_0-tau| = |y_0'*V'*(A-tau*I)*V*y_0| = |y_0'*V'*Q*R*y_0| = */
            /*       = |y_0'*V'*Q*u_0*s_0| <= s_0.                            */
            /* So when |l_0-tau|-machEps*|A| > s_0, we consider to reset the  */
            /* QR factorization. machEps*|A| is the error computing l_0.      */
            /* NOTE: rarely observed |l_0-tau|-machEps*|A| slightly greater   */
            /* than s_0 after resetting. The condition restartsSinceReset > 0 */
            /* avoids infinite loop in those cases.                           */

            if (primme->projectionParams.projection == primme_proj_refined &&
                  basisSize > 0 && restartsSinceReset > 0 &&
                  fabs(primme->targetShifts[targetShiftIndex]-hVals[0])
                    -max(primme->aNorm, primme->stats.estimateLargestSVal)
                      *machEps > hSVals[0]) {

               availableBlockSize = 0;
               targetShiftIndex = -1;
               reset = 2;
               if (primme->printLevel >= 5 && primme->procID == 0) {
                  fprintf(primme->outputFile, 
                        "Resetting V, W and QR: Some errors in QR detected.\n");
                  fflush(primme->outputFile);
               }

               break;
            }

           /* --------------------------------------------------------------- */
         } /* while (basisSize<maxBasisSize && basisSize<n-orthoConst-numLocked)
            * --------------------------------------------------------------- */

         wholeSpace = basisSize >= primme->n - primme->numOrthoConst - numLocked;

         /* ----------------------------------------------------------------- */
         /* Restart basis will need the final coefficient vectors in hVecs    */
         /* to lock out converged vectors and to compute X and R for the next */
         /* iteration. The function prepare_vecs will make sure that hVecs    */
         /* has proper coefficient vectors.                                   */
         /* ----------------------------------------------------------------- */

         if (targetShiftIndex >= 0) {
            /* -------------------------------------------------------------- */
            /* TODO: merge all this logic in the regular main loop. After all */
            /* restarting is a regular step that also shrinks the basis.      */
            /* -------------------------------------------------------------- */

            int j,k,l,m;
            double *dummySmallestResNorm, dummyZero = 0.0;

            if (blockSize > 0) {
               availableBlockSize = blockSize;
               maxRecentlyConverged = 0;
            }

            /* When there are more than one target shift,                     */
            /* limit blockSize and the converged values to one.               */

            else if (primme->numTargetShifts > numConverged+1) {
               if (primme->locking) {
                  maxRecentlyConverged =
                     max(min(primme->numEvals, numLocked+1) - numConverged, 0);
               }
               else {
                  maxRecentlyConverged =
                     max(min(primme->numEvals, numConverged+1) - numConverged, 0);
               }
               availableBlockSize = maxRecentlyConverged;
            }

            else {
               maxRecentlyConverged = max(0, primme->numEvals-numConverged);

               /* Limit blockSize to vacant vectors in the basis */

               availableBlockSize = max(0, min(primme->maxBlockSize, primme->maxBasisSize-(numConverged-numLocked)));

               /* Limit blockSize to remaining values to converge plus one */

               availableBlockSize = min(availableBlockSize, maxRecentlyConverged+1);
            }

            /* Limit basisSize to the matrix dimension */

            availableBlockSize = max(0, min(availableBlockSize, 
                  primme->n - numLocked - primme->numOrthoConst));

            /* -------------------------------------------------------------- */
            /* NOTE: setting smallestResNorm to zero may overpass the inner   */
            /* product condition (ip) and result in a smaller                 */
            /* numArbitraryVecs before restarting. This helps or is needed to */
            /* pass testi-100-LOBPCG_OrthoBasis_Window-100-primme_closest_abs */
            /* -primme_proj_refined.F                                         */
            /* -------------------------------------------------------------- */

            if (primme->target == primme_closest_abs ||
                  primme->target == primme_largest_abs) {
               dummySmallestResNorm = &dummyZero;
            }
            else {
               dummySmallestResNorm = &smallestResNorm;
            }

            prepare_candidates_Sprimme(V, ldV, W, ldW, primme->nLocal, H,
                  primme->maxBasisSize, basisSize,
                  NULL, NULL,
                  hVecs, basisSize, hVals, hSVals, flags,
                  maxRecentlyConverged, blockNorms, blockSize,
                  availableBlockSize, evecs, numLocked, ldevecs, evals,
                  resNorms, targetShiftIndex, machEps, iev, &blockSize,
                  &recentlyConverged, &numArbitraryVecs, dummySmallestResNorm,
                  hVecsRot, primme->maxBasisSize, numConverged, basisNorms,
                  &reset, rwork, &rworkSize, iwork, iworkSize, primme);
            assert(recentlyConverged >= 0);

            /* When QR is computed and there are more than one target shift   */
            /* and some eigenpair has converged values to one, don't provide  */
            /* candidate for the next iteration, because it may be the closest*/
            /* to a different target.                                         */

            if (Q && numConverged+recentlyConverged > numLocked
                  && primme->numTargetShifts > numLocked+1) {
               blockSize = 0;
            }

            /* Updated the number of converged pairs */
            /* Intentionally we include the pairs flagged SKIP_UNTIL_RESTART */

            for (i=0, numConverged=numLocked; i<basisSize; i++) {
               if (flags[i] != UNCONVERGED && numConverged < primme->numEvals
                     && (i < primme->numEvals-numLocked
                        /* Refined and prepare_vecs may not completely    */
                        /* order pairs considering closest_leq/geq; so we */
                        /* find converged pairs beyond the first remaining*/
                        /* pairs to converge.                             */
                        || primme->target == primme_closest_geq
                        || primme->target == primme_closest_leq)) {

                  numConverged++;

               }
            }

            /* Move the converged pairs and the ones in iev at the beginning */

            for (i=k=l=m=0; i<basisSize; i++) {
               int iIsInIev = 0;
               for (j=0; j<blockSize; j++)
                  if (iev[j] == i) iIsInIev = 1;
               if ((flags[i] != UNCONVERGED && m++ < numConverged-numLocked)
                     || iIsInIev) {
                  iwork[k++] = i;
               }
               else {
                  iwork[numConverged-numLocked+blockSize+l++] = i;
               }
            }
            permute_vecs_Rprimme(hVals, 1, basisSize, 1, iwork, (REAL*)rwork, iwork+basisSize);
            permute_vecs_Sprimme(hVecs, basisSize, basisSize, basisSize, iwork, rwork, iwork+basisSize);
            permute_vecs_iprimme(flags, basisSize, iwork, iwork+basisSize);
            if (hVecsRot) {
               Num_zero_matrix_Sprimme(&hVecsRot[numArbitraryVecs*primme->maxBasisSize], primme->maxBasisSize,
                     basisSize-numArbitraryVecs, primme->maxBasisSize);
               for (i=numArbitraryVecs; i<basisSize; i++)
                  hVecsRot[primme->maxBasisSize*i+i] = 1.0;
               permute_vecs_Sprimme(hVecsRot, basisSize, basisSize, primme->maxBasisSize, iwork, rwork, iwork+basisSize);
               for (i=j=0; i<basisSize; i++)
                  if (iwork[i] != i) j=i+1;
               numArbitraryVecs = max(numArbitraryVecs, j);
            }
         }

         /* ------------------ */
         /* Restart the basis  */
         /* ------------------ */

         int oldNumLocked = numLocked;
         assert(ldV == ldW); /* this function assumes ldV == ldW */
         restart_Sprimme(V, W, primme->nLocal, basisSize, ldV, hVals, hSVals,
               flags, iev, &blockSize, blockNorms, evecs, ldevecs, perm,
               evals, resNorms, evecsHat, primme->nLocal, M, maxEvecsSize, UDU,
               0, ipivot, &numConverged, &numLocked, lockedFlags,
               &numConvergedStored, previousHVecs, &numPrevRetained,
               primme->maxBasisSize, numGuesses, prevRitzVals, &numPrevRitzVals,
               H, primme->maxBasisSize, Q, ldQ, R, primme->maxBasisSize,
               QtV, primme->maxBasisSize, hU, basisSize, 0, hVecs, basisSize, 0,
               &basisSize, &targetShiftIndex, &numArbitraryVecs, hVecsRot,
               primme->maxBasisSize, &restartsSinceReset, &reset, machEps,
               rwork, &rworkSize, iwork, iworkSize, primme);

         /* If there are any initial guesses remaining, then copy it */
         /* into the basis.                                          */

         if (numGuesses > 0) {
            /* Try to keep minRestartSize guesses in the search subspace */

            int numNew = max(0, min(
                     primme->minRestartSize + numConverged 
                     - (nextGuess - primme->numOrthoConst), numGuesses));

            /* Don't make the resulting basis size larger than maxBasisSize */

            numNew = max(0,
                  min(basisSize+numNew, primme->maxBasisSize) - basisSize);

            /* Don't increase basis size beyond matrix dimension */

            numNew = max(0,
                  min(basisSize+numNew+primme->numOrthoConst+numLocked,
                     primme->n)
                  - primme->numOrthoConst - numLocked - basisSize);

            Num_copy_matrix_Sprimme(&evecs[nextGuess*ldevecs], primme->nLocal,
                  numNew, ldevecs, &V[basisSize*ldV], ldV);

            nextGuess += numNew;
            numGuesses -= numNew;

            CHKERR(ortho_Sprimme(V, ldV, NULL, 0, basisSize, basisSize+numNew-1,
                     evecs, ldevecs, numLocked+primme->numOrthoConst,
                     primme->nLocal, primme->iseed, machEps, rwork, &rworkSize,
                     primme), -1);

            /* Compute W = A*V for the orthogonalized corrections */

            CHKERR(matrixMatvec_Sprimme(V, primme->nLocal, ldV, W, ldW,
                     basisSize, numNew, primme), -1);

            if (Q) CHKERR(update_Q_Sprimme(V, primme->nLocal, ldV, W, ldW, Q,
                     ldQ, R, primme->maxBasisSize,
                     primme->targetShifts[targetShiftIndex], basisSize, numNew,
                     rwork, &rworkSize, machEps, primme), -1);

            /* Extend H by numNew columns and rows and solve the */
            /* eigenproblem for the new H.                       */

            if (H) CHKERR(update_projection_Sprimme(V, ldV, W, ldW, H,
                     primme->maxBasisSize, primme->nLocal, basisSize, numNew,
                     rwork, &rworkSize, 1/*symmetric*/, primme), -1);
            if (QtV) CHKERR(update_projection_Sprimme(Q, ldQ, V, ldV, QtV,
                     primme->maxBasisSize, primme->nLocal, basisSize, numNew,
                     rwork, &rworkSize, 0/*unsymmetric*/, primme), -1);
            basisSize += numNew;
            CHKERR(solve_H_Sprimme(H, basisSize, primme->maxBasisSize, R,
                  primme->maxBasisSize, QtV, primme->maxBasisSize, hU,
                  basisSize, hVecs, basisSize, hVals, hSVals, numConverged,
                  machEps, &rworkSize, rwork, iworkSize, iwork, primme), -1);

         }
 
         primme->stats.numRestarts++;

         primme->initSize = numConverged;

         /* ------------------------------------------------------------- */
         /* If dynamic method switching == 1, update model parameters and */
         /* evaluate whether to switch from GD+k to JDQMR. This is after  */
         /* restart. GD+k is also evaluated if a pair converges.          */
         /* ------------------------------------------------------------- */
         if (primme->dynamicMethodSwitch == 1 ) {
            tstart = primme_wTimer(0);
            CostModel.MV = primme->stats.timeMatvec/primme->stats.numMatvecs;
            ret = update_statistics(&CostModel, primme, tstart, 0, 1,
               numConverged, blockNorms[0], primme->stats.estimateMaxEVal); 
            CHKERR(switch_from_GDpk(&CostModel, primme), -1);
         } /* ---------------------------------------------------------- */

         /* ----------------------------------------------------------- */
         /* After having exhausted the search subspace, everything      */
         /* should be converged. However the code can only mark one     */
         /* pair as converged per iteration when multiple shifts still  */
         /* remain active. In this situation, exit when no pair         */
         /* converges at some iteration.                                */
         /* ----------------------------------------------------------- */

         if (wholeSpace && (
                  primme->target == primme_largest
                  || primme->target == primme_smallest
                  || (primme->target != primme_closest_leq
                     && primme->target != primme_closest_geq
                     && oldNumLocked >= primme->numTargetShifts-1)
                  || oldNumLocked == numLocked))
               break;

        /* ----------------------------------------------------------- */
      } /* while ((numConverged < primme->numEvals)  (restarting loop)
         * ----------------------------------------------------------- */

      /* ------------------------------------------------------------ */
      /* If locking is enabled, check to make sure the required       */
      /* number of eigenvalues have been computed, else make sure the */
      /* residual norms of the converged Ritz vectors have remained   */
      /* converged by calling verify_norms.                           */
      /* ------------------------------------------------------------ */

      if (primme->locking) {

         /* if dynamic method, give method recommendation for future runs */
         if (primme->dynamicMethodSwitch > 0 ) {
            if (CostModel.accum_jdq_gdk < 0.96) 
               primme->dynamicMethodSwitch = -2;  /* Use JDQMR_ETol */
            else if (CostModel.accum_jdq_gdk > 1.04) 
               primme->dynamicMethodSwitch = -1;  /* Use GD+k */
            else
               primme->dynamicMethodSwitch = -3;  /* Close call. Use dynamic */
         }

         /* Return flag showing if there has been a locking problem */
         intWork[0] = LockingProblem;

         /* If all of the target eigenvalues have been computed, */
         /* then return success, else return with a failure.     */
 
         if (numConverged == primme->numEvals || wholeSpace) {
            if (primme->aNorm <= 0.0L) primme->aNorm = primme->stats.estimateLargestSVal;
            return 0;
         }
         else {
            CHKERRNOABORTM(-1, -1, "Maximum iterations or matvecs reached");
         }

      }
      else {      /* no locking. Verify that everything is converged  */

         /* Determine if the maximum number of matvecs has been reached */

         restartLimitReached = primme->maxMatvecs > 0 && 
                               primme->stats.numMatvecs >= primme->maxMatvecs;

         /* ---------------------------------------------------------- */
         /* The norms of the converged Ritz vectors must be recomputed */
         /* before return.  This is because the restarting may cause   */
         /* some converged Ritz vectors to become slightly unconverged.*/
         /* If some have become unconverged, then further iterations   */
         /* must be performed to force these approximations back to a  */
         /* converged state.                                           */
         /* ---------------------------------------------------------- */

         CHKERR(verify_norms(V, ldV, W, ldW, hVals, numConverged, resNorms,
                  flags, &converged, machEps, rwork, &rworkSize, iwork,
                  iworkSize, primme), -1);

         /* ---------------------------------------------------------- */
         /* If the convergence limit is reached or the target vectors  */
         /* have remained converged, then copy the current Ritz values */
         /* and vectors to the output arrays and return, else continue */
         /* iterating.                                                 */
         /* ---------------------------------------------------------- */

         if (restartLimitReached || converged || wholeSpace) {
            for (i=0; i < primme->numEvals; i++) {
               evals[i] = hVals[i];
               perm[i] = i;
            }

            Num_copy_matrix_Sprimme(V, primme->nLocal, primme->numEvals, ldV,
               &evecs[ldevecs*primme->numOrthoConst], ldevecs);

            /* The target values all remained converged, then return */
            /* successfully, else return with a failure code.        */
            /* Return also the number of actually converged pairs    */
 
            primme->initSize = numConverged;

            /* if dynamic method, give method recommendation for future runs */
            if (primme->dynamicMethodSwitch > 0 ) {
               if (CostModel.accum_jdq_gdk < 0.96) 
                  primme->dynamicMethodSwitch = -2;  /* Use JDQMR_ETol */
               else if (CostModel.accum_jdq_gdk > 1.04) 
                  primme->dynamicMethodSwitch = -1;  /* Use GD+k */
               else
                  primme->dynamicMethodSwitch = -3;  /* Close call.Use dynamic*/
            }

            if (converged) {
               if (primme->aNorm <= 0.0L) primme->aNorm = primme->stats.estimateLargestSVal;
               return 0;
            }
            else {
	       CHKERRNOABORTM(-1, -1, "Maximum iterations or matvecs reached");
            }

         }
         else if (!converged) {
            /* ------------------------------------------------------------ */
            /* Reorthogonalize the basis, recompute W=AV, and continue the  */
            /* outer while loop, resolving the epairs. Slow, but robust!    */
            /* ------------------------------------------------------------ */
            CHKERR(ortho_Sprimme(V, ldV, NULL, 0, 0,
                     basisSize-1, evecs, ldevecs,
                     primme->numOrthoConst+numLocked, primme->nLocal,
                     primme->iseed, machEps, rwork, &rworkSize, primme), -1);
            CHKERR(matrixMatvec_Sprimme(V, primme->nLocal, ldV, W, ldW, 0,
                     basisSize, primme), -1);

            if (primme->printLevel >= 2 && primme->procID == 0) {
               fprintf(primme->outputFile, 
                 "Verifying before return: Some vectors are unconverged.\n");
               fflush(primme->outputFile);
            }

            restartsSinceReset = 0;
            reset = 0;
            primme->stats.estimateResidualError = 0.0;

           /* ------------------------------------------------------------ */
         } /* End of elseif(!converged). Restart and recompute all epairs
            * ------------------------------------------------------------ */

        /* ------------------------------------------------------------ */
      } /* End of non locking
         * ------------------------------------------------------------ */

     /* -------------------------------------------------------------- */
   } /* while (!converged)  Outer verification loop
      * -------------------------------------------------------------- */

   if (primme->aNorm <= 0.0L) primme->aNorm = primme->stats.estimateLargestSVal;

   return 0;

}

/******************************************************************************
          Some basic functions within the scope of main_iter
*******************************************************************************/

/*******************************************************************************
 * Subroutine prepare_candidates - This subroutine puts into the block the first
 *    unconverged Ritz pairs, up to maxBlockSize. If needed, compute residuals
 *    and rearrange the coefficient vectors in hVecs.
 * 
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V              The orthonormal basis
 * W              A*V
 * nLocal         Local length of vectors in the basis
 * basisSize      Size of the basis V and W
 * ldV            The leading dimension of V and X
 * ldW            The leading dimension of W and R
 * hVecs          The projected vectors
 * ldhVecs        The leading dimension of hVecs
 * hVals          The Ritz values
 * maxBasisSize   maximum allowed size of the basis
 * numSoftLocked  Number of vectors that have converged (not updated here)
 * remainedEvals  Remained number of eigenpairs that the user wants computed
 * blockNormsSize Number of already computed residuals
 * maxBlockSize   maximum allowed size of the block
 * evecs          Converged eigenvectors
 * evecsSize      The size of evecs
 * numLocked      The number of vectors currently locked (if locking)
 * numConverged   Number of converged pairs (soft+hard locked)
 * rwork          Real work array, used by check_convergence and Num_update_VWXR
 * primme         Structure containing various solver parameters
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * X             The eigenvectors put in the block
 * R             The residual vectors put in the block
 * flags         Array indicating which eigenvectors have converged     
 * iev           indicates which eigenvalue each block vector corresponds to
 * blockNorms    Residual norms of the Ritz vectors being computed during the
 *               current iteration
 * blockSize     Dimension of the block
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * recentlyConverged  Number of pairs converged
 * reset         flag to reset V and W in the next restart
 * 
 ******************************************************************************/

TEMPLATE_PLEASE
int prepare_candidates_Sprimme(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, PRIMME_INT nLocal, SCALAR *H, int ldH, int basisSize,
      SCALAR *X, SCALAR *R, SCALAR *hVecs, int ldhVecs, REAL *hVals,
      REAL *hSVals, int *flags, int remainedEvals, REAL *blockNorms,
      int blockNormsSize, int maxBlockSize, SCALAR *evecs, int numLocked,
      PRIMME_INT ldevecs, REAL *evals, REAL *resNorms, int targetShiftIndex,
      double machEps, int *iev, int *blockSize, int *recentlyConverged,
      int *numArbitraryVecs, double *smallestResNorm, SCALAR *hVecsRot,
      int ldhVecsRot, int numConverged, REAL *basisNorms, int *reset,
      SCALAR *rwork, size_t *rworkSize, int *iwork, int iworkSize,
      primme_params *primme) {

   int i, blki;         /* loop variables */
   REAL *hValsBlock;    /* contiguous copy of the hVals to be tested */
   SCALAR *hVecsBlock;  /* contiguous copy of the hVecs columns to be tested */     
   int *flagsBlock;     /* contiguous copy of the flags to be tested */
   REAL *hValsBlock0;   /* workspace for hValsBlock */
   SCALAR *hVecsBlock0; /* workspace for hVecsBlock */
   double targetShift;  /* current target shift */
   size_t rworkSize0;   /* current size of rwork */
   int lasti;           /* last tested pair */

   /* -------------------------- */
   /* Return memory requirements */
   /* -------------------------- */

   if (V == NULL) {
      SCALAR t;
      REAL d;
      size_t lrw=0;
      int liw=0;

      CHKERR(check_convergence_Sprimme(NULL, nLocal, 0, NULL, 0, NULL,
               numLocked, 0, 0, basisSize, NULL, NULL,
               NULL, NULL, 0.0, NULL, &lrw, &liw, 0, primme), -1);
      lrw = max(lrw,
            (size_t)Num_update_VWXR_Sprimme(NULL, NULL, nLocal, basisSize,
               0, NULL, 0, 0, NULL, &t, basisSize-maxBlockSize, basisSize, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               &t, basisSize-maxBlockSize, basisSize, 0, &d,
               NULL, 0, 0,
               NULL, 0, primme));
      CHKERR(prepare_vecs_Sprimme(basisSize, 0, maxBlockSize, NULL, 0,
               NULL, NULL, NULL, 0, 0, NULL, 0.0, NULL, 0, NULL, 0, 0.0, &lrw,
               NULL, 0, &liw, primme), -1);
      *rworkSize = max(*rworkSize,
            (size_t)maxBlockSize+(size_t)maxBlockSize*(size_t)basisSize+lrw);
      *iwork = max(*iwork, liw + basisSize);
      return 0;
   }

   *blockSize = 0;
   hValsBlock0 = (REAL*)rwork;
   hVecsBlock0 = &rwork[maxBlockSize];
   rwork += maxBlockSize + ldhVecs*maxBlockSize;
   assert(*rworkSize >= (size_t)(maxBlockSize + ldhVecs*maxBlockSize));
   rworkSize0 = *rworkSize - maxBlockSize - ldhVecs*maxBlockSize;
   flagsBlock = iwork;
   iwork += maxBlockSize;
   iworkSize -= maxBlockSize;
   assert(iworkSize >= 0);
   targetShift = primme->targetShifts ? primme->targetShifts[targetShiftIndex] : 0.0;
   lasti = -1;

   /* Pack hVals for already computed residual pairs */

   hValsBlock = Num_compact_vecs_Rprimme(hVals, 1, blockNormsSize, 1, &iev[*blockSize],
         hValsBlock0, 1, 1 /* avoid copy */);

   /* If some residual norms have already been computed, set the minimum   */
   /* of them as the smallest residual norm. If not, use the value from    */
   /* previous iteration.                                                  */

   if (blockNormsSize > 0) {
      for (*smallestResNorm = HUGE_VAL, i=0; i < blockNormsSize; i++) {
         *smallestResNorm = min(*smallestResNorm, blockNorms[i]);
      }
   }

   *recentlyConverged = 0;
   while (1) {
      /* Recompute flags in iev(*blockSize:*blockSize+blockNormsize) */
      for (i=*blockSize; i<blockNormsSize; i++)
         flagsBlock[i-*blockSize] = flags[iev[i]];
      CHKERR(check_convergence_Sprimme(X?&X[(*blockSize)*ldV]:NULL, nLocal,
            ldV, R?&R[(*blockSize)*ldW]:NULL, ldW, evecs, numLocked,
            ldevecs, 0, blockNormsSize, flagsBlock,
            &blockNorms[*blockSize], hValsBlock, reset, machEps, rwork,
            &rworkSize0, iwork, iworkSize, primme), -1);

      /* Compact blockNorms, X and R for the unconverged pairs in    */
      /* iev(*blockSize:*blockSize+blockNormsize). Do the proper     */
      /* actions for converged pairs.                                */

      for (blki=*blockSize, i=0; i < blockNormsSize && *blockSize < maxBlockSize; i++, blki++) {
         /* Write back flags and residual norms */
         flags[iev[blki]] = flagsBlock[i];
         basisNorms[iev[blki]] = blockNorms[blki];

         /* Ignore some cases */
         if ((primme->target == primme_closest_leq
                  && hVals[iev[blki]]-blockNorms[blki] > targetShift) ||
               (primme->target == primme_closest_geq
                && hVals[iev[blki]]+blockNorms[blki] < targetShift)) {
         }
         else if (flagsBlock[i] != UNCONVERGED
                         && *recentlyConverged < remainedEvals
                         && (iev[blki] < primme->numEvals-numLocked
                            /* Refined and prepare_vecs may not completely    */
                            /* order pairs considering closest_leq/geq; so we */
                            /* find converged pairs beyond the first remaining*/
                            /* pairs to converge.                             */
                            || primme->target == primme_closest_geq
                            || primme->target == primme_closest_leq)) {

            /* Write the current Ritz value in evals and the residual in resNorms;  */
            /* it will be checked by restart routine later.                         */
            /* Also print the converged eigenvalue.                                 */

            if (!primme->locking) {
               evals[iev[blki]] = hVals[iev[blki]];
               resNorms[iev[blki]] = blockNorms[blki];
               primme->stats.maxConvTol = max(primme->stats.maxConvTol, blockNorms[blki]);
            }

            /* Count the new solution */
            (*recentlyConverged)++;

            /* Reset smallestResNorm if some pair converged */
            if (*blockSize == 0) {
               *smallestResNorm = HUGE_VAL;
            }

            /* Report a pair was soft converged */
            if (primme->monitorFun) {
               int ONE = 1, numConverged0 = numConverged+*recentlyConverged;
               primme_event EVENT_CONVERGED = primme_event_converged;
               int err;
               CHKERRM((primme->monitorFun(hVals, &basisSize, flags, &iev[blki],
                           &ONE, basisNorms, &numConverged0, NULL, NULL, NULL,
                           NULL, NULL, NULL, &EVENT_CONVERGED, primme, &err),
                        err), -1, "Error returned by monitorFun: %d", err);
            }
         }
         else if (flagsBlock[i] == UNCONVERGED) {
            /* Update the smallest residual norm */
            if (*blockSize == 0) {
               *smallestResNorm = HUGE_VAL;
            }
            *smallestResNorm = min(*smallestResNorm, blockNorms[blki]);

            blockNorms[*blockSize] = blockNorms[blki];
            iev[*blockSize] = iev[blki];
            if (X) Num_copy_matrix_Sprimme(&X[blki*ldV], nLocal, 1, ldV,
                  &X[(*blockSize)*ldV], ldV);
            if (R) Num_copy_matrix_Sprimme(&R[blki*ldW], nLocal, 1, ldW,
                  &R[(*blockSize)*ldW], ldW);
            (*blockSize)++;
         }

         lasti = iev[blki];
      }

      /* Generate well conditioned coefficient vectors; start from the last   */
      /* position visited (variable i)                                        */

      blki = *blockSize;
      prepare_vecs_Sprimme(basisSize, lasti+1, maxBlockSize-blki, H, ldH, hVals,
            hSVals, hVecs, ldhVecs, targetShiftIndex, numArbitraryVecs,
            *smallestResNorm, flags, 1, hVecsRot, ldhVecsRot, machEps,
            &rworkSize0, rwork, iworkSize, iwork, primme);

      /* Find next candidates, starting from iev(*blockSize)+1 */

      for (i=lasti+1; i<basisSize && blki < maxBlockSize; i++)
         if (flags[i] == UNCONVERGED) iev[blki++] = i;

      /* If no new candidates or all required solutions converged yet, go out */

      if (blki == *blockSize || *recentlyConverged >= remainedEvals) break;
      blockNormsSize = blki - *blockSize;

      /* Pack hVals & hVecs */

      hValsBlock = Num_compact_vecs_Rprimme(hVals, 1, blockNormsSize, 1, &iev[*blockSize],
         hValsBlock0, 1, 1 /* avoid copy */);
      hVecsBlock = Num_compact_vecs_Sprimme(hVecs, basisSize, blockNormsSize, ldhVecs, &iev[*blockSize],
         hVecsBlock0, ldhVecs, 1 /* avoid copy */);

      /* Compute X, R and residual norms for the next candidates                                   */
      /* X(basisSize:) = V*hVecs(*blockSize:*blockSize+blockNormsize)                              */
      /* R(basisSize:) = W*hVecs(*blockSize:*blockSize+blockNormsize) - X(basisSize:)*diag(hVals)  */
      /* blockNorms(basisSize:) = norms(R(basisSize:))                                             */

      assert(ldV == ldW); /* This functions only works in this way */
      CHKERR(Num_update_VWXR_Sprimme(V, W, nLocal, basisSize, ldV,
               hVecsBlock, basisSize, ldhVecs, hValsBlock,
               X?&X[(*blockSize)*ldV]:NULL, 0, blockNormsSize, ldV,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               R?&R[(*blockSize)*ldV]:NULL, 0, blockNormsSize, ldV, &blockNorms[*blockSize],
               !R?&blockNorms[*blockSize]:NULL, 0, blockNormsSize,
               rwork, rworkSize0, primme), -1);
   }

   return 0;
}

/*******************************************************************************
 * Function verify_norms - This subroutine computes the residual norms of the 
 *    target eigenvectors before the Davidson-type main iteration terminates. 
 *    This is done to insure that the eigenvectors have remained converged. 
 *    If any have become unconverged, this routine performs restarting before 
 *    the eigensolver iteration is once again performed.  
 *
 * Note: This routine assumes it is called immediately after a call to the 
 *       restart subroutine.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V            The orthonormal basis
 *
 * W            A*V
 *
 * hVals        The eigenvalues of V'*A*V
 *
 * basisSize    Size of the basis V
 *
 * rworkSize    Length of rwork
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * rwork   Must be at least 2*primme->numEvals in size
 *
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * resNorms      Residual norms for each eigenpair
 *
 * flag          Indicates which Ritz pairs have converged
 *
 * converged     Return 1 if the basisSize values are converged, and 0 otherwise
 *
 ******************************************************************************/
   
static int verify_norms(SCALAR *V, PRIMME_INT ldV, SCALAR *W, PRIMME_INT ldW,
      REAL *hVals, int basisSize, REAL *resNorms, int *flags, int *converged,
      double machEps, SCALAR *rwork, size_t *rworkSize, int *iwork,
      int iworkSize, primme_params *primme) {

   int i;         /* Loop variable                                     */
   REAL *dwork = (REAL *) rwork; /* pointer to cast rwork to REAL*/

   /* Compute the residual vectors */

   for (i=0; i < basisSize; i++) {
      Num_axpy_Sprimme(primme->nLocal, -hVals[i], &V[ldV*i], 1, &W[ldW*i], 1);
      dwork[i] = REAL_PART(Num_dot_Sprimme(primme->nLocal, &W[ldW*i],
               1, &W[ldW*i], 1));
   }
      
   CHKERR(globalSum_Rprimme(dwork, resNorms, basisSize, primme), -1);
   for (i=0; i < basisSize; i++)
      resNorms[i] = sqrt(resNorms[i]);

   /* Check for convergence of the residual norms. */

   CHKERR(check_convergence_Sprimme(V, primme->nLocal, ldV, W, ldW, NULL, 0,
            0, 0, basisSize, flags, resNorms, hVals, NULL, machEps, rwork,
            rworkSize, iwork, iworkSize, primme), -1);

   /* Set converged to 1 if the first basisSize pairs are converged */

   *converged = 1;
   for (i=0; i<basisSize; i++) {
      if (flags[i] == UNCONVERGED) {
         *converged = 0;
      }
   }
    
   return 0;
}

/******************************************************************************
           Dynamic Method Switching uses the following functions 
    ---------------------------------------------------------------------
     If primme->dynamicMethodSwitch > 0 find which of GD+k,JDQMR is best
     NOTE: JDQMR requires additional memory for inner iteration, so we 
       assume either the user has called primme_set_method(DYNAMIC) 
       or has allocated the appropriate space.
*******************************************************************************/

/******************************************************************************
 * Function switch_from_JDQMR - 
 *    If primme->dynamicMethodSwitch=2,4 try to switch from JDQMR_ETol to GD+k.
 *    Assumes that the CostModel has been updated through runtime measurements.
 *    Based on this CostModel, the switch occurs only if
 *
 *                          expected_JDQMR_ETol_time
 *                ratio =  --------------------------  > 1.05
 *                             expected_GD+k_time
 *
 *    There are two cases determining when this function is called.
 *
 *    primme->dynamicMethodSwitch = 2 
 *        numEvals < 5 (few eigenvalues). We must dynamically decide the 
 *        best method before an eigenvalue converges. Because a very slow
 *        and expensive inner iteration may take too many inner steps, we 
 *        must re-evaluate the ratio at every outer step, before the call 
 *        to solve the correction equation. 
 *        After a switch, dynamicMethodSwitch is 1.
 *        
 *    primme->dynamicMethodSwitch = 4 
 *        numEvals > 4 (many eigenvalues). We can afford to check how both 
 *        methods converge to an eigenvalue, and then measure the statistics.
 *        In this case we re-evaluate the ratio every time one or more 
 *        eigenvalues converge, just before the next correction equation.
 *        After a switch, dynamicMethodSwitch is 3.
 *
 * INPUT/OUTPUT
 * ------------
 * model        The CostModel that contains all the model relevant parameters
 *              (accum_jdq_gdk, accum_jdq, accum_gdk updated)
 *
 * primme       The solver parameters (dynamicMethodSwitch, maxInnerIterations
 *              changed if there is a switch)
 *
 ******************************************************************************/

static int switch_from_JDQMR(primme_CostModel *model, primme_params *primme) {

   int switchto=0;
   double est_slowdown, est_ratio_MV_outer, ratio, globalRatio; 

   /* ----------------------------------------------------------------- */
   /* Asymptotic evaluation of the JDQMR versus GD+k for small numEvals */
   /* ----------------------------------------------------------------- */
   if (primme->dynamicMethodSwitch == 2) {

     /* For numEvals<4, (dyn=2), after we first estimate timings, decide if 
      * must always use GD+k (eg., because the operator is very expensive)
      * Use a best case scenario for JDQMR (small slowdown and many inner its)*/
      est_slowdown       = 1.1;     /* a small enough slowdown */
      est_ratio_MV_outer = 1000;    /* practically all time is in inner its */
      ratio = ratio_JDQMR_GDpk(model, 0, est_slowdown, est_ratio_MV_outer);

      /* If more many procs, make sure that all have the same ratio */
      if (primme->numProcs > 1) {
         CHKERR(globalSum_dprimme(&ratio, &globalRatio, 1, primme), -1); 
         ratio = globalRatio/primme->numProcs;
      }

      if (ratio > 1.05) { 
         /* Always use GD+k. No further model updates */
         primme->dynamicMethodSwitch = -1;
         primme->correctionParams.maxInnerIterations = 0;
         if (primme->printLevel >= 3 && primme->procID == 0) 
            fprintf(primme->outputFile, 
            "Ratio: %e Switching permanently to GD+k\n", ratio);
         return 0;
      }
   }

   /* Select method to switch to if needed: 2->1 and 4->3 */
   switch (primme->dynamicMethodSwitch) {
      case 2: switchto = 1; break;
      case 4: switchto = 3; 
   }

   /* ----------------------------------------------------------------- *
    * Compute the ratio of expected times JDQMR/GD+k. To switch to GD+k, the 
    * ratio must be > 1.05. Update accum_jdq_gdk for recommendation to user
    * ----------------------------------------------------------------- */

   ratio = ratio_JDQMR_GDpk(model, 0, model->JDQMR_slowdown, 
                                  model->ratio_MV_outer);

   /* If more many procs, make sure that all have the same ratio */
   if (primme->numProcs > 1) {
      CHKERR(globalSum_dprimme(&ratio, &globalRatio, 1, primme), -1);
      ratio = globalRatio/primme->numProcs;
   }
   
   if (ratio > 1.05) {
      primme->dynamicMethodSwitch = switchto; 
      primme->correctionParams.maxInnerIterations = 0;
   }

   model->accum_jdq += model->gdk_plus_MV_PR*ratio;
   model->accum_gdk += model->gdk_plus_MV_PR;
   model->accum_jdq_gdk = model->accum_jdq/model->accum_gdk;

   if (primme->printLevel >= 3 && primme->procID == 0) 
      switch (primme->correctionParams.maxInnerIterations) {
         case 0: fprintf(primme->outputFile, 
            "Ratio: %e JDQMR switched to GD+k\n", ratio); break;
         case -1: fprintf(primme->outputFile, 
            "Ratio: %e Continue with JDQMR\n", ratio);
      }

   return 0;

} /* end of switch_fromJDQMR() */

/******************************************************************************
 * Function switch_from_GDpk - 
 *    If primme->dynamicMethodSwitch=1,3 try to switch from GD+k to JDQMR_ETol.
 *    Assumes that the CostModel has been updated through runtime measurements.
 *    If no JDQMR measurements exist (1st time), switch to it unconditionally.
 *    Otherwise, based on the CostModel, the switch occurs only if
 *
 *                          expected_JDQMR_ETol_time
 *                ratio =  --------------------------  < 0.95
 *                             expected_GD+k_time
 *
 *    There are two cases determining when this function is called.
 *
 *    primme->dynamicMethodSwitch = 1
 *        numEvals < 5 (few eigenvalues). GD+k does not have an inner iteration
 *        so it is not statistically meaningful to check every outer step. 
 *        For few eigenvalues, we re-evaluate the ratio immediately after
 *        restart of the method. This allows also the cost of restart to 
 *        be included in the measurements.
 *        After a switch, dynamicMethodSwitch is 2.
 *        
 *    primme->dynamicMethodSwitch = 3 
 *        numEvals > 4 (many eigenvalues). As with JDQMR, we can check how both 
 *        methods converge to an eigenvalue, and then measure the statistics.
 *        In this case, we re-evaluate the ratio every time one or more 
 *        eigenvalues converge, just before the next preconditioner application
 *        After a switch, dynamicMethodSwitch is 4.
 *
 * INPUT/OUTPUT
 * ------------
 * model        The CostModel that contains all the model relevant parameters
 *              (accum_jdq_gdk, accum_jdq, accum_gdk updated)
 *
 * primme       The solver parameters (dynamicMethodSwitch, maxInnerIterations
 *              changed if there is a switch)
 *
 ******************************************************************************/
static int switch_from_GDpk(primme_CostModel *model, primme_params *primme) {

   int switchto=0;
   double ratio, globalRatio;

   /* if no restart has occurred (only possible under dyn=3) current timings */
   /* do not include restart costs. Remain with GD+k until a restart occurs */
   if (primme->stats.numRestarts == 0) return 0;

   /* Select method to switch to if needed: 1->2 and 3->4 */
   switch (primme->dynamicMethodSwitch) {
      case 1: switchto = 2; break;
      case 3: switchto = 4; 
   }

   /* If JDQMR never run before, switch to it to get first measurements */
   if (model->qmr_only == 0.0) {
      primme->dynamicMethodSwitch = switchto;
      primme->correctionParams.maxInnerIterations = -1;
      if (primme->printLevel >= 3 && primme->procID == 0) 
         fprintf(primme->outputFile, 
         "Ratio: N/A  GD+k switched to JDQMR (first time)\n");
      return 0;
   }

   /* ------------------------------------------------------------------- *
    * Compute the ratio of expected times JDQMR/GD+k. To switch to JDQMR, the
    * ratio must be < 0.95. Update accum_jdq_gdk for recommendation to user
    * ------------------------------------------------------------------- */

   ratio = ratio_JDQMR_GDpk(model, 0, model->JDQMR_slowdown, 
                                  model->ratio_MV_outer);

   /* If more many procs, make sure that all have the same ratio */
   if (primme->numProcs > 1) {
      CHKERR(globalSum_dprimme(&ratio, &globalRatio, 1, primme), -1); 
      ratio = globalRatio/primme->numProcs;
   }

   if (ratio < 0.95) {
      primme->dynamicMethodSwitch = switchto;
      primme->correctionParams.maxInnerIterations = -1;
   } 

   model->accum_jdq += model->gdk_plus_MV_PR*ratio;
   model->accum_gdk += model->gdk_plus_MV_PR;
   model->accum_jdq_gdk = model->accum_jdq/model->accum_gdk;

   if (primme->printLevel >= 3 && primme->procID == 0) 
      switch (primme->correctionParams.maxInnerIterations) {
         case 0: fprintf(primme->outputFile, 
            "Ratio: %e Continue with GD+k\n", ratio); break;
         case -1: fprintf(primme->outputFile, 
            "Ratio: %e GD+k switched to JDQMR\n", ratio);
      }

   return 0;
}

/******************************************************************************
 * Function update_statistics - 
 *
 *    Performs runtime measurements and updates the cost model that describes
 *    the average cost of Matrix-vector and Preconditioning operations, 
 *    the average cost of running 1 full iteration of GD+k and of JDQMR, 
 *    the number of inner/outer iterations since last update,
 *    the current convergence rate measured for each of the two methods,
 *    and based on these rates, the expected slowdown of JDQMR over GD+k
 *    in terms of matrix-vector operations.
 *    Times are averaged with one previous measurement, and convergence
 *    rates are averaged over a window which is reset over 10 converged pairs.
 *
 *    The function is called right before switch_from_JDQMR/switch_from_GDpk.
 *    Depending on the current method running, this is at two points:
 *       If some eigenvalues just converged, called before solve_correction()
 *       If primme->dynamicMethodSwitch = 2, called before solve_correction()
 *       If primme->dynamicMethodSwitch = 1, called after restart()
 *
 *    The Algorithm
 *    -------------
 *    The dynamic switching algorithm starts with GD+k. After the 1st restart
 *    (for dyn=1) or after an eigenvalue converges (dyn=3), we collect the 
 *    first measurements for GD+k and switch to JDQMR_ETol. We collect the 
 *    first measurements for JDQMR after 1 outer step (dyn=2) or after an
 *    eigenvalue converges (dyn=4). From that time on, the method is chosen
 *    dynamically by switch_from_JDQMR/switch_from_GDpk using the ratio.
 *
 *    Cost/iteration breakdown
 *    ------------------------
 *    kout: # of outer iters. kinn: # of inner QMR iters. nMV = # of Matvecs
 *    All since last call to update_statistics.
 *    1 outer step: costJDQMR = costGD_outer_method + costQMR_Iter*kinn + mv+pr
 *    
 *        <---------1 step JDQMR----------->
 *       (GDout)(---------QMR--------------)
 *        gd mv  pr q+mv+pr .... q+mv+pr mv  
 *                  <-------kinn------->      kinn = nMV/kout - 2
 *        (-----)(------time_in_inner------)  
 *
 *    The model
 *    ---------
 *    time_in_inner = kout (pr+kinn*(q+mv+pr)+mv) 
 *                  = kout (pr+mv) + nMV*(q+mv+pr) - 2kout(q+mv+pr)
 *    time_in_outer = elapsed_time-time_in_inner = kout*(gd+mv)
 *    JDQMR_time = time_in_inner+time_in_outer = 
 *               = kout(gd+mv) + kout(pr+mv) + nMV(q+mv+pr) -2kout(q+mv+pr)
 *               = kout (gd -2q -pr) + nMV(q+mv+pr)
 *    GDpk_time  = gdOuterIters*(gd+mv+pr)
 *
 *    Letting slowdown = (nMV for JDQMR)/(gdOuterIters) we have the ratio:
 *
 *          JDQMR_time     q+mv+pr + kout/nMV (gd-2q-pr)
 *          ----------- = ------------------------------- slowdown
 *           GDpk_time               gd+mv+pr
 *
 *    Because the QMR of JDQMR "wastes" one MV per outer step, and because
 *    its number of outer steps cannot be more than those of GD+k
 *           (kinn+2)/(kinn+1) < slowdown < kinn+2
 *    However, in practice 1.1 < slowdown < 2.5. 
 *
 *    For more details see papers [1] and [2] (listed in primme_S.c)
 *
 *
 * INPUT
 * -----
 * primme           Structure containing the solver parameters
 * current_time     Time stamp taken before update_statistics is called
 * recentConv       Number of converged pairs since last update_statistics
 * calledAtRestart  True if update_statistics is called at restart by dyn=1
 * numConverged     Total number of converged pairs
 * currentResNorm   Residual norm of the next unconverged epair
 * aNormEst         Estimate of ||A||_2. Conv Tolerance = aNormEst*primme.eps
 *
 * INPUT/OUTPUT
 * ------------
 * model            The model parameters updated
 *
 * Return value
 * ------------
 *     0    Not enough iterations to update model. Continue with current method
 *     1    Model updated. Proceed with relative evaluation of the methods
 *
 ******************************************************************************/
static int update_statistics(primme_CostModel *model, primme_params *primme,
   double current_time, int recentConv, int calledAtRestart, int numConverged, 
   double currentResNorm, double aNormEst) {

   double low_res, elapsed_time, time_in_outer, kinn;
   int kout, nMV;

   /* ------------------------------------------------------- */
   /* Time in outer and inner iteration since last update     */
   /* ------------------------------------------------------- */
   elapsed_time = current_time - model->timer_0;
   time_in_outer = elapsed_time - model->time_in_inner;

   /* ------------------------------------------------------- */
   /* Find # of outer, MV, inner iterations since last update */
   /* ------------------------------------------------------- */
   kout = primme->stats.numOuterIterations - model->numIt_0;
   nMV  = primme->stats.numMatvecs - model->numMV_0;
   if (calledAtRestart) kout++; /* Current outer iteration is complete,  */
                                /*   but outerIterations not incremented yet */
   if (kout == 0) return 0;     /* No outer iterations. No update or evaluation
                                 *   Continue with current method */
   kinn = ((double) nMV)/kout - 2;

   if (primme->correctionParams.maxInnerIterations == -1 && 
       (kinn < 1.0 && model->qmr_only == 0.0L) )  /* No inner iters yet, and */
      return 0;              /* no previous QMR timings. Continue with JDQMR */

   /* --------------------------------------------------------------- */
   /* After one or more pairs converged, currentResNorm corresponds   */
   /* to the next unconverged pair. To measure the residual reduction */
   /* during the previous step, we must use the convergence tolerance */
   /* Also, update how many evals each method found since last reset  */
   /* --------------------------------------------------------------- */
   if (recentConv > 0) {
      /* Use tolerance as the lowest residual norm to estimate conv rate */
      if (primme->aNorm > 0.0L) 
         low_res = primme->eps*primme->aNorm;
      else 
         low_res = primme->eps*aNormEst;
      /* Update num of evals found */
      if (primme->correctionParams.maxInnerIterations == -1)
          model->nevals_by_jdq += recentConv;
      else
          model->nevals_by_gdk += recentConv;
   }
   else 
      /* For dyn=1 at restart, and dyn=2 at every step. Use current residual */
      low_res = currentResNorm;

   /* ------------------------------------------------------- */
   /* Update model timings and parameters                     */
   /* ------------------------------------------------------- */

   /* update outer iteration time for both GD+k,JDQMR.Average last two updates*/
   if (model->gdk_plus_MV == 0.0L) 
      model->gdk_plus_MV = time_in_outer/kout;
   else 
      model->gdk_plus_MV = (model->gdk_plus_MV + time_in_outer/kout)/2.0L;

   /* ---------------------------------------------------------------- *
    * Reset the conv rate averaging window every 10 converged pairs. 
    * See also Notes 2 and 3 below.
    * Note 1: For large numEvals, we should average the convergence 
    * rate over only the last few converged pairs. To avoid a more 
    * expensive moving window approach, we reset the window when >= 10 
    * additional pairs converge. To avoid a complete reset, we consider 
    * the current average rate to be the new rate for the "last" pair.
    * By scaling down the sums: sum_logResReductions and sum_MV, this 
    * rate does not dominate the averaging of subsequent measurements.
    * ---------------------------------------------------------------- */

   if (numConverged/10 >= model->nextReset) {
      model->gdk_sum_logResReductions /= model->nevals_by_gdk;
      model->gdk_sum_MV /= model->nevals_by_gdk;
      model->jdq_sum_logResReductions /= model->nevals_by_jdq;
      model->jdq_sum_MV /= model->nevals_by_jdq;
      model->nextReset = numConverged/10+1;
      model->nevals_by_gdk = 1;
      model->nevals_by_jdq = 1;
   }

   switch (primme->dynamicMethodSwitch) {

      case 1: case 3: /* Currently running GD+k */
        /* Update Precondition times */
        if (model->PR == 0.0L) 
           model->PR          = model->time_in_inner/kout;
        else 
           model->PR          = (model->PR + model->time_in_inner/kout)/2.0L;
        model->gdk_plus_MV_PR = model->gdk_plus_MV + model->PR;
        model->MV_PR          = model->MV + model->PR;

        /* update convergence rate.
         * ---------------------------------------------------------------- *
         * Note 2: This is NOT a geometric average of piecemeal rates. 
         * This is the actual geometric average of all rates per MV, or 
         * equivalently the total rate as TotalResidualReductions over
         * the corresponding nMVs. If the measurement intervals (in nMVs)
         * are identical, the two are equivalent. 
         * Note 3: In dyn=1,2, we do not record residual norm increases.
         * This overestimates convergence rates a little, but otherwise, a 
         * switch would leave the current method with a bad rate estimate.
         * ---------------------------------------------------------------- */

        if (low_res <= model->resid_0) 
           model->gdk_sum_logResReductions += log(low_res/model->resid_0);
        model->gdk_sum_MV += nMV;
        model->gdk_conv_rate = exp( model->gdk_sum_logResReductions
                                   /model->gdk_sum_MV );
        break;

      case 2: case 4: /* Currently running JDQMR */
        /* Basic timings for QMR iteration (average of last two updates) */
        if (model->qmr_plus_MV_PR == 0.0L) {
           model->qmr_plus_MV_PR=(model->time_in_inner/kout- model->MV_PR)/kinn;
           model->ratio_MV_outer = ((double) nMV)/kout;
        }
        else {
           if (kinn != 0.0) model->qmr_plus_MV_PR = (model->qmr_plus_MV_PR  + 
              (model->time_in_inner/kout - model->MV_PR)/kinn )/2.0L;
           model->ratio_MV_outer =(model->ratio_MV_outer+((double) nMV)/kout)/2;
        }
        model->qmr_only = model->qmr_plus_MV_PR - model->MV_PR;

        /* Update the cost of a hypothetical GD+k, as measured outer + PR */
        model->gdk_plus_MV_PR = model->gdk_plus_MV + model->PR;

        /* update convergence rate */
        if (low_res <= model->resid_0) 
           model->jdq_sum_logResReductions += log(low_res/model->resid_0);
        model->jdq_sum_MV += nMV;
        model->jdq_conv_rate = exp( model->jdq_sum_logResReductions
                                   /model->jdq_sum_MV);
        break;
   }
   update_slowdown(model);

   /* ------------------------------------------------------- */
   /* Reset counters to measure statistics at the next update */
   /* ------------------------------------------------------- */
   model->numIt_0 = primme->stats.numOuterIterations;
   if (calledAtRestart) model->numIt_0++; 
   model->numMV_0 = primme->stats.numMatvecs;
   model->timer_0 = current_time;      
   model->time_in_inner = 0.0;
   model->resid_0 = currentResNorm;

   return 1;
}

/******************************************************************************
 * Function ratio_JDQMR_GDpk -
 *    Using model parameters, computes the ratio of expected times:
 *
 *          JDQMR_time     q+mv+pr + kout/nMV (gd-2q-pr)
 *          ----------- = ------------------------------- slowdown
 *           GDpk_time               gd+mv+pr
 *
 ******************************************************************************/
static double ratio_JDQMR_GDpk(primme_CostModel *model, int numLocked,
   double estimate_slowdown, double estimate_ratio_MV_outer) {
   
   return estimate_slowdown* 
     ( model->qmr_plus_MV_PR + model->project_locked*numLocked + 
       (model->gdk_plus_MV - model->qmr_only - model->qmr_plus_MV_PR  
        + (model->reortho_locked - model->project_locked)*numLocked  
       )/estimate_ratio_MV_outer
     ) / (model->gdk_plus_MV_PR + model->reortho_locked*numLocked);
}

/******************************************************************************
 * Function update_slowdown -
 *    Given the model measurements for convergence rates, computes the slowdown
 *           log(GDpk_conv_rate)/log(JDQMR_conv_rate)
 *    subject to the bounds 
 *    max(1.1, (kinn+2)/(kinn+1)) < slowdown < min(2.5, kinn+2)
 *
 ******************************************************************************/
static void update_slowdown(primme_CostModel *model) {
  double slowdown;

  if (model->gdk_conv_rate < 1.0) {
     if (model->jdq_conv_rate < 1.0) 
        slowdown = log(model->gdk_conv_rate)/log(model->jdq_conv_rate);
     else if (model->jdq_conv_rate == 1.0)
        slowdown = 2.5;
     else
        slowdown = -log(model->gdk_conv_rate)/log(model->jdq_conv_rate);
  }
  else if (model->gdk_conv_rate == 1.0) 
        slowdown = 1.1;
  else { /* gdk > 1 */
     if (model->jdq_conv_rate < 1.0) 
        slowdown = log(model->gdk_conv_rate)/log(model->jdq_conv_rate);
     else if (model->jdq_conv_rate == 1.0)
        slowdown = 1.1;
     else  /* both gdk, jdq > 1 */
        slowdown = log(model->jdq_conv_rate)/log(model->gdk_conv_rate);
  }

  /* Slowdown cannot be more than the matvecs per outer iteration */
  /* nor less than MV per outer iteration/(MV per outer iteration-1) */
  slowdown = max(model->ratio_MV_outer/(model->ratio_MV_outer-1.0),
                                    min(slowdown, model->ratio_MV_outer));
  /* Slowdown almost always in [1.1, 2.5] */
  model->JDQMR_slowdown = max(1.1, min(slowdown, 2.5));
}

/******************************************************************************
 * Function initializeModel - Initializes model members
 ******************************************************************************/
static void initializeModel(primme_CostModel *model, primme_params *primme) {
   model->MV_PR          = 0.0L;
   model->MV             = 0.0L;
   model->PR             = 0.0L;
   model->qmr_only       = 0.0L;
   model->qmr_plus_MV_PR = 0.0L;
   model->gdk_plus_MV_PR = 0.0L;
   model->gdk_plus_MV    = 0.0L;
   model->project_locked = 0.0L;
   model->reortho_locked = 0.0L;

   model->gdk_conv_rate  = 0.0001L;
   model->jdq_conv_rate  = 0.0001L;
   model->JDQMR_slowdown = 1.5L;
   model->ratio_MV_outer = 0.0L;

   model->nextReset      = 1;
   model->gdk_sum_logResReductions = 0.0;
   model->jdq_sum_logResReductions = 0.0;
   model->gdk_sum_MV     = 0.0;
   model->jdq_sum_MV     = 0.0;
   model->nevals_by_gdk  = 0;
   model->nevals_by_jdq  = 0;

   model->numMV_0 = primme->stats.numMatvecs;
   model->numIt_0 = primme->stats.numOuterIterations+1;
   model->timer_0 = primme_wTimer(0);
   model->time_in_inner  = 0.0L;
   model->resid_0        = -1.0L;

   model->accum_jdq      = 0.0L;
   model->accum_gdk      = 0.0L;
   model->accum_jdq_gdk  = 1.0L;
}

#if 0
/******************************************************************************
 *
 * Function to display the model parameters --- For debugging purposes only
 *
 ******************************************************************************/
static void displayModel(primme_CostModel *model){
   fprintf(stdout," MV %e\n", model->MV);
   fprintf(stdout," PR %e\n", model->PR);
   fprintf(stdout," MVPR %e\n", model->MV_PR);
   fprintf(stdout," QMR %e\n", model->qmr_only);
   fprintf(stdout," QMR+MVPR %e\n", model->qmr_plus_MV_PR);
   fprintf(stdout," GD+MVPR %e\n", model->gdk_plus_MV_PR);
   fprintf(stdout," GD+MV %e\n\n", model->gdk_plus_MV);
   fprintf(stdout," project %e\n", model->project_locked);
   fprintf(stdout," reortho %e\n", model->reortho_locked);
   fprintf(stdout," GD convrate %e\n", model->gdk_conv_rate);
   fprintf(stdout," JD convrate %e\n", model->jdq_conv_rate);
   fprintf(stdout," slowdown %e\n", model->JDQMR_slowdown);
   fprintf(stdout," outer/Totalmv %e\n", model->ratio_MV_outer);
   
   fprintf(stdout," gd sum_logResRed %e\n",model->gdk_sum_logResReductions);
   fprintf(stdout," jd sum_logResRed %e\n",model->jdq_sum_logResReductions);
   fprintf(stdout," gd sum_MV %e\n", model->gdk_sum_MV);
   fprintf(stdout," jd sum_MV %e\n", model->jdq_sum_MV);
   fprintf(stdout," gd nevals %d\n", model->nevals_by_gdk);
   fprintf(stdout," jd nevals %d\n", model->nevals_by_jdq);

   fprintf(stdout," Accumul jd %e\n", model->accum_jdq);
   fprintf(stdout," Accumul gd %e\n", model->accum_gdk);
   fprintf(stdout," Accumul ratio %e\n", model->accum_jdq_gdk);

   fprintf(stdout,"0: MV %d IT %d timer %e timInn %e res %e\n", model->numMV_0,
      model->numIt_0,model->timer_0,model->time_in_inner, model->resid_0);
   fprintf(stdout," ------------------------------\n");
}
#endif
