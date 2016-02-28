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
 * File: main_iter.c
 *
 * Purpose - This is the main Davidson-type iteration 
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "primme.h"
#include "const.h"
#include "wtime.h"
#include "main_iter_@(pre).h"
#include "main_iter_private_@(pre).h"
#include "convergence_@(pre).h"
#include "correction_@(pre).h"
#include "init_@(pre).h"
#include "ortho_@(pre).h"
#include "restart_@(pre).h"
#include "locking_@(pre).h"
#include "solve_H_@(pre).h"
#include "update_projection_@(pre).h"
#include "update_W_@(pre).h"
#include "numerical_@(pre).h"

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
 * realWork @(type) work array
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

int main_iter_@(pre)primme(double *evals, int *perm, @(type) *evecs, 
   double *resNorms, double machEps, int *intWork, void *realWork, 
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
   int maxRecentlyConverged;/* Maximum converged values per iteration        */
   int numConvergedStored;  /* Numb of Ritzvecs temporarily stored in evecs  */
                            /*    to allow for skew projectors w/o locking   */
   int converged;           /* True when all required Ritz vals. converged   */
   int LockingProblem;      /* Flag==1 if practically converged pairs locked */
   int restartLimitReached; /* True when maximum restarts performed          */
   int numPrevRetained;     /* Number of vectors retained using recurrence-  */
                            /* based restarting.                             */
   int maxEvecsSize;        /* Maximum capacity of evecs array               */
   int rworkSize;           /* Size of rwork array                           */
   int numPrevRitzVals = 0; /* Size of the prevRitzVals updated in correction*/
   int restartSize;         /* Basis size after restarting                   */
   int indexOfPreviousVecs; /* Column index in hVecs with previous vecs      */
   int ret;                 /* Return value                                  */

   int *iwork;              /* Integer workspace pointer                     */
   int *iwork0;             /* Temporal integer workspace pointer            */
   int *flags;              /* Indicates which Ritz values have converged    */
   int *ipivot;             /* The pivot for the UDU factorization of M      */
   int *iev;                /* Evalue index each block vector corresponds to */
   int *restartPerm;        /* Permutation of hVecs used to restart V        */
   int *hVecsPerm;          /* Permutation of hVecs to sort as primme.target */

   double tol;              /* Required tolerance for residual norms         */
   @(type) *V;              /* Basis vectors                                 */
   @(type) *W;              /* Work space storing A*V                        */
   @(type) *H;              /* Upper triangular portion of V'*A*V            */
   @(type) *M = NULL;       /* The projection Q'*K*Q, where Q = [evecs, x]   */
                            /* x is the current Ritz vector and K is a       */
                            /* hermitian preconditioner.                     */
   @(type) *UDU = NULL;     /* The factorization of M=Q'KQ                   */
   @(type) *evecsHat = NULL;/* K^{-1}evecs                                   */
   @(type) *rwork;          /* Real work space.                              */
   @(type) *hVecs;          /* Eigenvectors of H                             */
   @(type) *hU=NULL;        /* Left singular vectors of R                    */
   @(type) *previousHVecs;  /* Coefficient vectors retained by               */
                            /* recurrence-based restarting                   */

   int numQR;               /* Maximum number of QR factorizations           */
   @(type) *Q = NULL;       /* QR decompositions for harmonic or refined     */
   @(type) *R = NULL;       /* projection: (A-target[i])*V = QR              */
   @(type) *QV = NULL;      /* Q'*V                                          */

   double *hVals;           /* Eigenvalues of H                              */
   double *hSVals=NULL;     /* Singular values of R                          */
   double *prevRitzVals;    /* Eigenvalues of H at previous outer iteration  */
                            /* by robust shifting algorithm in correction.c  */
   double *blockNorms;      /* Residual norms corresponding to current block */
                            /* vectors.                                      */
   @(type) tpone = @(tpone);/* constant 1.0 of type @(type) */
   @(type) tzero = @(tzero);/* constant 0.0 of type @(type) */

   /* Runtime measurement variables for dynamic method switching             */
   primme_CostModel CostModel; /* Structure holding the runtime estimates of */
                            /* the parameters of the model.Only visible here */
   double timeForMV;        /* Measures time for 1 matvec operation          */
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

   rwork         = (@(type) *) realWork;
   V             = rwork; rwork += primme->nLocal*primme->maxBasisSize;
   W             = rwork; rwork += primme->nLocal*primme->maxBasisSize;
   if (numQR > 0) {
      Q          = rwork; rwork += primme->nLocal*primme->maxBasisSize*numQR;
      R          = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize*numQR;
      hU         = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize*numQR;
   }
   if (primme->projectionParams.projection == primme_proj_Harm) {
      QV         = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize*numQR;
   }
   H             = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize;
   hVecs         = rwork; rwork += primme->maxBasisSize*primme->maxBasisSize;
   previousHVecs = rwork; rwork += primme->maxBasisSize*primme->restartingParams.maxPrevRetain;

   if (primme->correctionParams.precondition && 
         primme->correctionParams.maxInnerIterations != 0 &&
         primme->correctionParams.projectors.RightQ &&
         primme->correctionParams.projectors.SkewQ           ) {
      evecsHat   = rwork; rwork += primme->nLocal*maxEvecsSize;
      M          = rwork; rwork += maxEvecsSize*maxEvecsSize;
      UDU        = rwork; rwork += maxEvecsSize*maxEvecsSize;
   }

   hVals         = (double *)rwork; rwork += primme->maxBasisSize*sizeof(double)/sizeof(@(type)) + 1;
   if (numQR > 0) {
      hSVals     = (double *)rwork; rwork += primme->maxBasisSize*sizeof(double)/sizeof(@(type)) + 1;
   }
   prevRitzVals  = (double *)rwork; rwork += (primme->maxBasisSize+primme->numEvals)*sizeof(double)/sizeof(@(type)) + 1;
   blockNorms    = (double *)rwork; rwork += primme->maxBlockSize*sizeof(double)/sizeof(@(type)) + 1;

   rworkSize     = primme->realWorkSize/sizeof(@(type)) - (rwork - (@(type)*)realWork);

   /* Integer workspace */

   flags = intWork;
   iev = flags + primme->maxBasisSize;
   ipivot = iev + primme->maxBlockSize;
   iwork = ipivot + maxEvecsSize;

   /* -------------------------------------------------------------- */
   /* Initialize counters and flags                                  */
   /* -------------------------------------------------------------- */

   primme->stats.numOuterIterations = 0;
   primme->stats.numRestarts = 0;
   primme->stats.numMatvecs = 0;
   numLocked = 0;
   converged = FALSE;
   LockingProblem = 0;

   numPrevRetained = 0;
   blockSize = 0; 

   /* ---------------------------------------- */
   /* Set the tolerance for the residual norms */
   /* ---------------------------------------- */

   primme->stats.estimateMaxEVal   = -HUGE_VAL;
   primme->stats.estimateMinEVal   = HUGE_VAL;
   primme->stats.estimateLargestSVal = -HUGE_VAL;
   primme->stats.maxConvTol        = 0.0L;
   if (primme->aNorm > 0.0L) {
      tol = primme->eps*primme->aNorm;
   }
   else {
      tol = primme->eps; /* tol*estimateLargestSVal will be checked */
   }

   /* -------------------------------------- */
   /* Quick return for matrix of dimension 1 */
   /* -------------------------------------- */

   if (primme->n == 1) {
      evecs[0] = tpone;
      matrixMatvec_@(pre)primme(&evecs[0], primme->nLocal, primme->nLocal,
            W, primme->nLocal, 0, 1, primme);
#ifdefarithm L_DEFREAL
      evals[0] = W[0];
#endifarithm
#ifdefarithm L_DEFCPLX
      evals[0] = W[0].r;
#endifarithm
      V[0] = tpone;

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

   ret = init_basis_@(pre)primme(V, primme->nLocal, primme->nLocal, W,
         primme->nLocal, evecs, primme->nLocal, evecsHat, primme->nLocal,
         M, maxEvecsSize, UDU, 0, ipivot, machEps, rwork, rworkSize, &basisSize,
         &nextGuess, &numGuesses, &timeForMV, primme);

   if (ret < 0) {
      primme_PushErrorMessage(Primme_main_iter, Primme_init_basis, ret, 
                      __FILE__, __LINE__, primme);
      return INIT_FAILURE;
   }

   /* ----------------------------------------------------------- */
   /* Dynamic method switch means we need to decide whether to    */
   /* allow inner iterations based on runtime timing measurements */
   /* ----------------------------------------------------------- */
   if (primme->dynamicMethodSwitch > 0) {
      initializeModel(&CostModel, primme);
      CostModel.MV = timeForMV;
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
      reset_flags_@(pre)primme(flags, 0, primme->maxBasisSize-1);

      /* Compute the initial H and solve for its eigenpairs */

      targetShiftIndex = 0;
      if (Q) update_Q_@(pre)primme(V, primme->nLocal, primme->nLocal, W, primme->nLocal, Q,
            primme->nLocal, R, primme->maxBasisSize, primme->targetShifts[targetShiftIndex], 0,
            basisSize, rwork, rworkSize, machEps, primme);

      if (H) update_projection_@(pre)primme(V, primme->nLocal, W, primme->nLocal, H,
            primme->maxBasisSize, primme->nLocal, 0, basisSize, rwork,
            rworkSize, 1/*symmetric*/, primme);

      if (QV) update_projection_@(pre)primme(Q, primme->nLocal, V, primme->nLocal, QV,
            primme->maxBasisSize, primme->nLocal, 0, basisSize, rwork,
            rworkSize, 0/*unsymmetric*/, primme);

      ret = solve_H_@(pre)primme(H, basisSize, primme->maxBasisSize, R,
            primme->maxBasisSize, QV, primme->maxBasisSize, hU, basisSize, hVecs,
            basisSize, hVals, hSVals, numConverged, machEps, rworkSize, rwork,
            iwork, primme);
      
      if (ret != 0) {
         primme_PushErrorMessage(Primme_main_iter, Primme_solve_h, ret, 
                         __FILE__, __LINE__, primme);
         return SOLVE_H_FAILURE;
      }

      /* -------------------------------------------------------------- */
      /* Begin the iterative process.  Keep restarting until all of the */
      /* required eigenpairs have been found (no verification)          */
      /* -------------------------------------------------------------- */
      while (numConverged < primme->numEvals &&
             ( primme->maxMatvecs == 0 || 
               primme->stats.numMatvecs < primme->maxMatvecs ) &&
             ( primme->maxOuterIterations == 0 ||
               primme->stats.numOuterIterations < primme->maxOuterIterations) ) {
 
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
            numPrevRetained = 0;

            /* When QR are computed and there are more than one target shift,	*/
            /* limit blockSize and the converged values to one.               */

            if (primme->numTargetShifts > numConverged+1 && R) {
               availableBlockSize = 1;
               maxRecentlyConverged = numConverged-numLocked+1;
            }
            else {
               availableBlockSize = primme->maxBlockSize;
               maxRecentlyConverged = primme->numEvals-numLocked;
            }

            /* Limit blockSize to vacant vectors in the basis */

            availableBlockSize = min(availableBlockSize, primme->maxBasisSize-basisSize);

            assert(blockSize <= availableBlockSize);

            /* Set the block with the first unconverged pairs */
            prepare_candidates_@(pre)(V, W, primme->nLocal, basisSize, primme->nLocal,
               &V[basisSize*primme->nLocal], &W[basisSize*primme->nLocal], hVecs, basisSize,
               hVals, flags, numConverged-numLocked, maxRecentlyConverged, blockNorms,
               blockSize, availableBlockSize, evecs, numLocked, evals, resNorms, machEps,
               iev, &blockSize, &recentlyConverged, rwork, rworkSize, iwork, primme);

            /* print residuals */
            print_residuals(hVals, blockNorms, numConverged, numLocked, iev, blockSize,
               primme);

            /* If the total number of converged pairs, including the     */
            /* recentlyConverged ones, are greater than or equal to the  */
            /* target number of eigenvalues, attempt to restart, verify  */
            /* their convergence, lock them if necessary, and return.    */
            /* For locking interior, restart and lock now any converged. */
            /* Restart if converging a new value implies recompute QR.   */

            numConverged += recentlyConverged;

            if (numConverged >= primme->numEvals ||
                (primme->locking && recentlyConverged > 0
                    && primme->target != primme_smallest
                    && primme->target != primme_largest) ||
                (Q && min(primme->numTargetShifts, numConverged) !=
                        min(primme->numTargetShifts, numConverged-recentlyConverged))
                              ) {
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
               Num_scal_@(pre)primme(blockSize*primme->nLocal, tzero,
                  &V[primme->nLocal*basisSize], 1);
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
                     ret = update_statistics(&CostModel, primme, tstart, 
                        recentlyConverged, 0, numConverged, blockNorms[0], 
                        primme->stats.estimateLargestSVal); 

                     if (ret) switch (primme->dynamicMethodSwitch) {
                        /* for few evals (dyn=1) evaluate GD+k only at restart*/
                        case 3: switch_from_GDpk(&CostModel,primme); break;
                        case 2: case 4: switch_from_JDQMR(&CostModel,primme);
                     } /* of if-switch */
                  } /* of recentlyConv > 0 || dyn==2 */
               } /* dynamic switching */
               /* - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

               ret = solve_correction_@(pre)primme(V, W, evecs, evecsHat, UDU, 
                 ipivot, evals, numLocked, numConvergedStored, hVals, 
                 prevRitzVals, &numPrevRitzVals, flags, basisSize, blockNorms, 
                 iev, blockSize, tol, machEps, primme->stats.estimateLargestSVal,
                 rwork, iwork, rworkSize, primme);

               if (ret != 0) {
                  primme_PushErrorMessage(Primme_main_iter, 
                     Primme_solve_correction, ret, __FILE__, __LINE__, primme);
                  return SOLVE_CORRECTION_FAILURE;
               }

               /* ------------------------------------------------------ */
               /* If dynamic method switch, accumulate inner method time */
               /* ------------------------------------------------------ */
               if (primme->dynamicMethodSwitch > 0) 
                  CostModel.time_in_inner += primme_wTimer(0) - tstart;

              
            } /* end of else blocksize=0 */

            /* Orthogonalize the corrections with respect to each other */
            /* and the current basis.                                   */
            ret = ortho_@(pre)primme(V, primme->nLocal, NULL, 0, basisSize, 
               basisSize+blockSize-1, evecs, primme->nLocal, 
               primme->numOrthoConst+numLocked, primme->nLocal, primme->iseed, 
               machEps, rwork, rworkSize, primme);

            if (ret < 0) {
               primme_PushErrorMessage(Primme_main_iter, Primme_ortho, ret,
                               __FILE__, __LINE__, primme);
               return ORTHO_FAILURE;
            }
           
            /* Compute W = A*V for the orthogonalized corrections */

            matrixMatvec_@(pre)primme(V, primme->nLocal, primme->nLocal, W, primme->nLocal,
                  basisSize, blockSize, primme);

            if (Q) update_Q_@(pre)primme(V, primme->nLocal, primme->nLocal, W, primme->nLocal, Q,
                  primme->nLocal, R, primme->maxBasisSize,
                  primme->targetShifts[targetShiftIndex], basisSize,
                  blockSize, rwork, rworkSize, machEps, primme);

            /* If harmonic, the coefficient vectors (i.e., the eigenvectors of the  */
            /* projected problem) are in hU; so retain them.                        */

            numPrevRetained = retain_previous_coefficients(QV ? hU : hVecs, 
               previousHVecs, basisSize, iev, blockSize, primme);

            /* Extend H by blockSize columns and rows and solve the */
            /* eigenproblem for the new H.                          */

            if (H) update_projection_@(pre)primme(V, primme->nLocal, W, primme->nLocal, H,
                  primme->maxBasisSize, primme->nLocal, basisSize, blockSize, rwork,
                  rworkSize, 1/*symmetric*/, primme);

            if (QV) update_projection_@(pre)primme(Q, primme->nLocal, V, primme->nLocal, QV,
                  primme->maxBasisSize, primme->nLocal, basisSize, blockSize, rwork,
                  rworkSize, 0/*unsymmetric*/, primme);

            basisSize += blockSize;
            blockSize = 0;

            ret = solve_H_@(pre)primme(H, basisSize, primme->maxBasisSize, R,
                  primme->maxBasisSize, QV, primme->maxBasisSize, hU, basisSize, hVecs,
                  basisSize, hVals, hSVals, numConverged, machEps, rworkSize, rwork,
                  iwork, primme);

            if (ret != 0) {
               primme_PushErrorMessage(Primme_main_iter, Primme_solve_h, ret,
                               __FILE__, __LINE__, primme);
               return SOLVE_H_FAILURE;
            }
            
           /* --------------------------------------------------------------- */
         } /* while (basisSize<maxBasisSize && basisSize<n-orthoConst-numLocked)
            * --------------------------------------------------------------- */

         /* ------------------ */
         /* Restart the basis  */
         /* ------------------ */

         /* ----------------------------------------------------------- */
         /* Special case: If (basisSize+numLocked) is the entire space, */
         /* then everything should be converged. Do not test, just flag */
         /* everything as converged                                     */
         /* ----------------------------------------------------------- */

         if (basisSize + numLocked + primme->numOrthoConst >= primme->n) {
            for (i = 0; i < basisSize && numConverged < primme->numEvals; i++)
               if (flags[i] == UNCONVERGED) { flags[i] = CONVERGED; numConverged++; }
            restartSize = basisSize;
            numPrevRetained = 0;
         }

         /* --------------------------------------------------------------------- */
         /* If basis isn't full, restart with the current basis size.             */
         /* If dynamic thick restarting is to be used, then determine the minimum */
         /* number of free spaces to be maintained and call the DTR routine.      */
         /* The DTR routine will determine how many coefficient vectors from the  */
         /* left and right of H-spectrum to retain at restart. If DTR is not used */
         /* then set the restart size to the minimum restart size.                */
         /* --------------------------------------------------------------------- */
      
         else if (basisSize <= primme->maxBasisSize - primme->maxBlockSize) {
             restartSize = basisSize;
         }
         else if (primme->restartingParams.scheme == primme_dtr) {
            int numFree = numPrevRetained+max(3, primme->maxBlockSize);
            restartSize = dtr_@(pre)(numLocked, hVecs, hVals, flags, basisSize, numFree, 
                  iev, rwork, primme);
         }
         else {
            restartSize = min(basisSize, primme->minRestartSize);
         }

         restartPerm = iwork;
         hVecsPerm = &restartPerm[basisSize];
         iwork0 = &hVecsPerm[basisSize];

         if (!primme->locking) {
            @(type) *X, *Res;
            ret = restart_@(pre)primme(&restartSize, V, W, primme->nLocal, R, primme->maxBasisSize,
               hU, basisSize, basisSize, primme->nLocal, &X, &Res, hVecs, basisSize, restartPerm,
               hVals, flags, iev, &blockSize, blockNorms, evecs, evals, resNorms, evecsHat,
               primme->nLocal, M, maxEvecsSize, &numConverged, &numConvergedStored, previousHVecs,
               &numPrevRetained, primme->maxBasisSize, &indexOfPreviousVecs, hVecsPerm, machEps,
               rwork, rworkSize, iwork0, primme);
         }
         else {
            @(type) *X, *Res;
            ret = restart_locking_@(pre)primme(&restartSize, V, W, primme->nLocal, R,
               primme->maxBasisSize, hU, basisSize, basisSize, primme->nLocal, &X, &Res, hVecs, basisSize,
               restartPerm, hVals, flags, iev, &blockSize, blockNorms, evecs, evals, &numConverged,
               &numLocked, resNorms, perm, numGuesses, previousHVecs, &numPrevRetained,
               primme->maxBasisSize, &indexOfPreviousVecs, hVecsPerm, machEps, rwork, rworkSize,
               iwork0, primme);
         }

         if (ret != 0) {
            primme_PushErrorMessage(Primme_main_iter, Primme_restart, 
                            basisSize, __FILE__, __LINE__, primme);
            return RESTART_FAILURE;
         }

         /* Rearrange prevRitzVals according to restartPerm */

         if (primme->target != primme_smallest && primme->target != primme_largest) {
            permute_vecs_d(prevRitzVals, 1, basisSize, 1, restartPerm, (double*)rwork, iwork0);
            permute_vecs_d(prevRitzVals, 1, restartSize, 1, hVecsPerm, (double*)rwork, iwork0);
            numPrevRitzVals = restartSize;
         }

         after_restart_@(pre)primme(V, primme->nLocal, W, primme->nLocal, H,
               primme->maxBasisSize, Q, primme->nLocal, primme->nLocal, R,
               primme->maxBasisSize, QV, primme->maxBasisSize, hU, basisSize, restartSize,
               hVecs, basisSize, restartSize, hVals , hSVals, restartPerm, hVecsPerm,
               restartSize, basisSize, numPrevRetained, indexOfPreviousVecs, evecs,
               &numConvergedStored, primme->nLocal, evecsHat, primme->nLocal, M,
               maxEvecsSize, UDU, 0, ipivot, targetShiftIndex, numConverged,
               rworkSize, rwork, iwork0, machEps, primme);

         targetShiftIndex = min(primme->numTargetShifts-1, numConverged);
         basisSize = restartSize;


         /* If there are any initial guesses remaining, then copy it */
         /* into the basis, else flag the vector as locked so it may */
         /* be discarded later.                                      */

         if (numGuesses > 0) {
            int numNew = min(primme->minRestartSize-basisSize, numGuesses);

            Num_copy_matrix_@(pre)primme(&evecs[nextGuess*primme->nLocal], primme->nLocal,
                  numNew, primme->nLocal, &V[basisSize*primme->nLocal], primme->nLocal);

            nextGuess += numNew;
            numGuesses -= numNew;

            ret = ortho_@(pre)primme(V, primme->nLocal, NULL, 0, basisSize, 
                  basisSize+numNew-1, evecs, primme->nLocal, numLocked, 
                  primme->nLocal, primme->iseed, machEps, rwork, rworkSize, primme);

            if (ret < 0) {
               primme_PushErrorMessage(Primme_lock_vectors, Primme_ortho, ret, 
                     __FILE__, __LINE__, primme);
               return ORTHO_FAILURE;
            }   

            /* Compute W = A*V for the orthogonalized corrections */

            matrixMatvec_@(pre)primme(V, primme->nLocal, primme->nLocal, W, primme->nLocal,
                  basisSize, numNew, primme);

            if (Q) update_Q_@(pre)primme(V, primme->nLocal, primme->nLocal, W, primme->nLocal, Q,
                  primme->nLocal, R, primme->maxBasisSize,
                  primme->targetShifts[targetShiftIndex], basisSize,
                  numNew, rwork, rworkSize, machEps, primme);

            /* If harmonic, the coefficient vectors (i.e., the eigenvectors of the  */
            /* projected problem) are in hU; so retain them.                        */

            numPrevRetained = retain_previous_coefficients(QV ? hU : hVecs, 
               previousHVecs, basisSize, iev, numNew, primme);

            /* Extend H by numNew columns and rows and solve the */
            /* eigenproblem for the new H.                       */

            if (H) update_projection_@(pre)primme(V, primme->nLocal, W, primme->nLocal, H,
                  primme->maxBasisSize, primme->nLocal, basisSize, numNew, rwork,
                  rworkSize, 1/*symmetric*/, primme);
            if (QV) update_projection_@(pre)primme(Q, primme->nLocal, V, primme->nLocal, QV,
                  primme->maxBasisSize, primme->nLocal, basisSize, numNew, rwork,
                  rworkSize, 0/*unsymmetric*/, primme);
            basisSize += numNew;
            ret = solve_H_@(pre)primme(H, basisSize, primme->maxBasisSize, R,
                  primme->maxBasisSize, QV, primme->maxBasisSize, hU, basisSize, hVecs,
                  basisSize, hVals, hSVals, numConverged, machEps, rworkSize, rwork, iwork,
                  primme);

            if (ret != 0) {
               primme_PushErrorMessage(Primme_main_iter, Primme_solve_h, ret,
                               __FILE__, __LINE__, primme);
               return SOLVE_H_FAILURE;
            }
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
            ret = update_statistics(&CostModel, primme, tstart, 0, 1,
               numConverged, blockNorms[0], primme->stats.estimateMaxEVal); 
            switch_from_GDpk(&CostModel, primme);
         } /* ---------------------------------------------------------- */


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
 
         if (numConverged == primme->numEvals) {
            if (primme->aNorm <= 0.0L) primme->aNorm = primme->stats.estimateMaxEVal;
            return 0;
         }
         else {
            return MAX_ITERATIONS_REACHED;
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

         converged = verify_norms(V, W, hVecs, hVals, basisSize, resNorms, 
            flags, tol, primme->stats.estimateMaxEVal, rwork, &numConverged, primme);

         /* ---------------------------------------------------------- */
         /* If the convergence limit is reached or the target vectors  */
         /* have remained converged, then copy the current Ritz values */
         /* and vectors to the output arrays and return, else continue */
         /* iterating.                                                 */
         /* ---------------------------------------------------------- */

         if (restartLimitReached || converged) {
            for (i=0; i < primme->numEvals; i++) {
               evals[i] = hVals[i];
               perm[i] = i;
            }

            Num_@(pre)copy_@(pre)primme(primme->nLocal*primme->numEvals, V, 1, 
               &evecs[primme->nLocal*primme->numOrthoConst], 1);

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
               if (primme->aNorm <= 0.0L) primme->aNorm = primme->stats.estimateMaxEVal;
               return 0;
            }
            else {
               return MAX_ITERATIONS_REACHED;
            }

         }
         else if (!converged) {
            /* ------------------------------------------------------------ */
            /* Reorthogonalize the basis, recompute W=AV, and continue the  */
            /* outer while loop, resolving the epairs. Slow, but robust!    */
            /* ------------------------------------------------------------ */
            ret = ortho_@(pre)primme(V, primme->nLocal, NULL, 0, 0, basisSize-1, evecs, 
               primme->nLocal, primme->numOrthoConst+numLocked, primme->nLocal,
               primme->iseed, machEps, rwork, rworkSize, primme);
            if (ret < 0) {
               primme_PushErrorMessage(Primme_main_iter, Primme_ortho, ret,
                               __FILE__, __LINE__, primme);
               return ORTHO_FAILURE;
            }
            matrixMatvec_@(pre)primme(V, primme->nLocal, primme->nLocal, W,
                  primme->nLocal, 0, basisSize, primme);

            if (primme->printLevel >= 2 && primme->procID == 0) {
               fprintf(primme->outputFile, 
                 "Verifying before return: Some vectors are unconverged. ");
               fprintf(primme->outputFile, "Restarting at #MV %d\n",
                 primme->stats.numMatvecs);
               fflush(primme->outputFile);
            }

           /* ------------------------------------------------------------ */
         } /* End of elseif(!converged). Restart and recompute all epairs
            * ------------------------------------------------------------ */

        /* ------------------------------------------------------------ */
      } /* End of non locking
         * ------------------------------------------------------------ */

     /* -------------------------------------------------------------- */
   } /* while (!converged)  Outer verification loop
      * -------------------------------------------------------------- */

   if (primme->aNorm <= 0.0L) primme->aNorm = primme->stats.estimateMaxEVal;

   return 0;

}

/******************************************************************************
          Some basic functions within the scope of main_iter
*******************************************************************************/

/*******************************************************************************
 * Subroutine prepare_candidates - This subroutine put into the block the first
 *    unconverged Ritz pairs, up to maxBlockSize.
 * 
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V              The orthonormal basis
 * W              A*V
 * nLocal         Local length of vectors in the basis
 * basisSize      Size of the basis V and W
 * ldV            The leading dimension of V, W, X and R
 * hVecs          The projected vectors
 * ldhVecs        The leading dimension of hVecs
 * hVals          The Ritz values
 * maxBasisSize   maximum allowed size of the basis
 * numSoftLocked  Number of vectors that have converged (not updated here)
 * numEvals       Remained number of eigenpairs that the user wants computed
 * blockNormsSize Number of already computed residuals
 * maxBlockSize   maximum allowed size of the block
 * evecs          Converged eigenvectors
 * evecsSize      The size of evecs
 * numLocked      The number of vectors currently locked (if locking)
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
 ******************************************************************************/

int prepare_candidates_@(pre)(@(type) *V, @(type) *W, int nLocal, int basisSize,
   int ldV, @(type) *X, @(type) *R, @(type) *hVecs, int ldhVecs, double *hVals,
   int *flags, int numSoftLocked, int numEvals, double *blockNorms,
   int blockNormsSize, int maxBlockSize, @(type) *evecs, int numLocked,
   double *evals, double *resNorms, double machEps, int *iev, int *blockSize,
   int *recentlyConverged, @(type) *rwork, int rworkSize, int *iwork,
   primme_params *primme) {

   int i, blki, ret;
   double *hValsBlock, *hValsBlock0;
   @(type) *hVecsBlock, *hVecsBlock0;
   int *flagsBlock;

   /* -------------------------- */
   /* Return memory requirements */
   /* -------------------------- */

   if (V == NULL) {
      @(type) t;

      return maxBlockSize+maxBlockSize*basisSize+max(
         check_convergence_@(pre)primme(NULL, nLocal, 0, NULL, 0, NULL, numLocked, 0,
               basisSize-maxBlockSize, basisSize, NULL, NULL, NULL, 0.0, NULL, 0, NULL, primme),
         Num_update_VWXR_@(pre)(NULL, NULL, nLocal, basisSize, 0, NULL, 0, 0, NULL,
               &t, basisSize-maxBlockSize, basisSize, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               NULL, 0, 0, 0,
               &t, basisSize-maxBlockSize, basisSize, 0, &machEps,
               NULL, 0, 0,
               NULL, 0, primme));
   }

   *blockSize = 0;
   hValsBlock0 = (double*)rwork;
   hVecsBlock0 = &rwork[maxBlockSize];
   rwork += maxBlockSize + ldhVecs*maxBlockSize;
   rworkSize -= maxBlockSize + ldhVecs*maxBlockSize;
   flagsBlock = iwork;
   iwork += maxBlockSize;

   /* Pack hVals */

   hValsBlock = Num_compact_vecs_d(hVals, 1, blockNormsSize, 1, &iev[*blockSize],
         hValsBlock0, 1, 1 /* avoid copy */);

   *recentlyConverged = 0;
   while (1) {

      /* Recompute flags in iev */
      ret = check_convergence_@(pre)primme(&X[(*blockSize)*ldV], nLocal, ldV,
         &R[(*blockSize)*ldV], ldV, evecs, numLocked, primme->nLocal, 0, blockNormsSize, flagsBlock,
         &blockNorms[*blockSize], hValsBlock, machEps, rwork, rworkSize, iwork, primme);
      if (ret != 0) return ret;

      /* Compact blockNorms, X and R for the unconverged pairs between left      */
      /* and right                                                               */

      for (blki=*blockSize, i=0; i < blockNormsSize && *blockSize < maxBlockSize; i++, blki++) {
         if (flagsBlock[i] == UNCONVERGED) {
            blockNorms[*blockSize] = blockNorms[blki];
            iev[*blockSize] = iev[blki];
            Num_copy_matrix_@(pre)primme(&X[blki*ldV], nLocal, 1, ldV,
                  &X[(*blockSize)*ldV], ldV);
            Num_copy_matrix_@(pre)primme(&R[blki*ldV], nLocal, 1, ldV,
                  &R[(*blockSize)*ldV], ldV);
            (*blockSize)++;
         }

         else {
            /* Write the current Ritz value in evals and the residual in resNorms;  */
            /* it will be checked by restart routine later.                         */
            /* Also print the converged eigenvalue.                                 */

            if (!primme->locking) {
               if (primme->procID == 0 && primme->printLevel >= 2)
                  fprintf(primme->outputFile, 
                        "#Converged %d eval[ %d ]= %e norm %e Mvecs %d Time %g\n",
                        iev[blki]-*blockSize, iev[blki], hVals[iev[blki]], blockNorms[blki],
                        primme->stats.numMatvecs, primme_wTimer(0));
               evals[iev[blki]] = hVals[iev[blki]];
               resNorms[iev[blki]] = blockNorms[blki];
            }

            /* Write back flags and count the new solution */

            flags[iev[blki]] = flagsBlock[i];
            (*recentlyConverged)++;
         }
      }

      /* Find next candidates */

      for (i=0, blki=*blockSize; i<basisSize && i<numEvals && blki < maxBlockSize; i++)
         if (flags[i] == UNCONVERGED) iev[blki++] = i;

      /* If no new candidates, go out */

      if (blki == *blockSize) break;
      blockNormsSize = blki - *blockSize;

      /* Pack hVals & hVecs */

      hValsBlock = Num_compact_vecs_d(hVals, 1, blockNormsSize, 1, &iev[*blockSize],
         hValsBlock0, 1, 1 /* avoid copy */);
      hVecsBlock = Num_compact_vecs_@(pre)(hVecs, basisSize, blockNormsSize, ldhVecs, &iev[*blockSize],
         hVecsBlock0, ldhVecs, 1 /* avoid copy */);

      /* Compute X, R and residual norms for the next candidates */
      /* X(basisSize:) = V*hVecs(left:right-1)                                   */
      /* R(basisSize:) = W*hVecs(left:right-1) - X(basisSize:)*diag(hVals)       */
      /* blockNorms(basisSize:) = norms(R(basisSize:))                           */

      ret = Num_update_VWXR_@(pre)(V, W, nLocal, basisSize, ldV, hVecsBlock, basisSize,
         ldhVecs, hValsBlock,
         &X[(*blockSize)*ldV], 0, blockNormsSize, ldV,
         NULL, 0, 0, 0,
         NULL, 0, 0, 0,
         NULL, 0, 0, 0,
         &R[(*blockSize)*ldV], 0, blockNormsSize, ldV, &blockNorms[*blockSize],
         NULL, 0, 0,
         rwork, rworkSize, primme);
      if (ret != 0) return ret;

   }

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
 * hVecs      The coefficients (eigenvectors of the projection H).
 *
 * basisSize  The current size of the basis
 *
 * iev        Array of size block size.  It maps the block index to the Ritz
 *            value index each block vector corresponds to.
 *
 * blockSize  The number of block vectors generated during the current iteration
 *
 * primme       Structure containing various solver parameters
 *
 *
 * Output parameters
 * -----------------
 * previousHVecs  The coefficients to be retained
 *
 *
 * Return value
 * ------------
 * The number of vectors retained
 *
 ******************************************************************************/

static int retain_previous_coefficients(@(type) *hVecs, @(type) *previousHVecs, 
   int basisSize, int *iev, int blockSize, primme_params *primme) {

   int i, j;            /* Loop indices                                  */
   int index;           /* The index of some coefficient vector in hVecs */ 
   int numPrevRetained; /* The number of coefficent vectors retained     */
   @(type) tzero = @(tzero);

   numPrevRetained = 0;

   /* If coefficient vectors are to be retained and its the iteration  */
   /* before restart occurs, then retain at most maxPrevRetain vectors */
 
   if (primme->restartingParams.maxPrevRetain > 0 && 
       basisSize+blockSize >= primme->maxBasisSize)
   {
      index = -1;

      /* ------------------------------------------------------------- */
      /* Retain as many coefficients corresponding to unconverged Ritz */
      /* vectors as possible.                                          */
      /* ------------------------------------------------------------- */

      for (i = 0; i < primme->restartingParams.maxPrevRetain; i++) {

         /* First, retain coefficient vectors corresponding to current block */
         /* vectors.  If all of those have been retained, then retain the    */ 
         /* the next coefficient beyond iev[blockSize-1].                    */

         if (i < blockSize) {
            index = iev[i];
         }
         else {
            index++;
         }

         /* If there is a coefficient vector at index index, then retain it */

         if (index < basisSize) {
            Num_@(pre)copy_@(pre)primme(basisSize, &hVecs[basisSize*index], 1, 
               &previousHVecs[primme->maxBasisSize*numPrevRetained], 1);

            /* Zero the maxBasisSize-basisSize last elements of the buffer */

            for (j = basisSize; j < primme->maxBasisSize; j++) {
               previousHVecs[primme->maxBasisSize*numPrevRetained+j] = tzero;
            } 

            numPrevRetained++;
         }
         else {
            break;
         }
      }

      if (primme->printLevel >= 5 && primme->procID == 0) {
         fprintf(primme->outputFile, "retain_previous: numPrevRetained: %d\n",
                 numPrevRetained);
      }

   }

   return numPrevRetained;
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
 * hVecs        The eigenvectors of V'*A*V
 *
 * hVals        The eigenvalues of V'*A*V
 *
 * basisSize    Size of the basis V
 *
 * tol          Required tolerance for the residual norms
 *
 * aNormEstimate if primme->aNorm<=0, use tol*aNormEstimate
 * 
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
 * numConverged  Number of eigenpairs that remained converged
 *
 *
 * RETURN VALUE
 * ------------
 * TRUE if the numReqEvals Ritz values have remained converged after restarting,
 * FALSE otherwise.
 ******************************************************************************/
   
static int verify_norms(@(type) *V, @(type) *W, @(type) *hVecs, 
   double *hVals, int basisSize, double *resNorms, int *flag, double tol, 
   double aNormEstimate, void *rwork, int *numConverged, primme_params *primme){

   int i;         /* Loop varible                                      */
   int converged; /* True when all requested Ritz values are converged */
   int nev, n;    /* convenience integers for numEvals and nLocal      */
   double *dwork = (double *) rwork; /* pointer to cast rwork to double*/
#ifdefarithm L_DEFCPLX
   @(type) ztmp;  /* temp complex var */
#endifarithm

   nev = primme->numEvals;
   n   = primme->nLocal;

   /* Set up the tolerance if necessary */

   if (primme->aNorm <= 0.0L) {
      tol = tol * aNormEstimate;
   }

   /* Compute the residual vectors */

   for (i=0; i < nev; i++) {
#ifdefarithm L_DEFCPLX
      {ztmp.r = -hVals[i]; ztmp.i = 0.0L;}
      Num_axpy_zprimme(n, ztmp, &V[n*i], 1, &W[n*i], 1);
      ztmp = Num_dot_zprimme(n, &W[n*i], 1, &W[n*i], 1);
      dwork[nev+i] = ztmp.r;
#endifarithm
#ifdefarithm L_DEFREAL
      Num_axpy_@(pre)primme(n, -hVals[i], &V[n*i], 1, &W[n*i], 1);
      dwork[nev+i] = Num_dot_dprimme(n, &W[n*i], 1, &W[n*i], 1);
#endifarithm

   }
      
   (*primme->globalSumDouble)(&dwork[nev], &dwork[0], &nev, primme); 
   converged = 1;

   /* Check for convergence of the residual norms. */

   for (i=0; i < nev; i++) {
      dwork[i] = sqrt(dwork[i]);
      
      if (dwork[i] < tol) {
         flag[i] = CONVERGED;
      }
      else {
         converged = 0;
         flag[i] = UNCONVERGED;
         *numConverged = *numConverged - 1;
      }

      resNorms[i] = dwork[i];
   }

   return converged;
}

/*******************************************************************************
 * Subroutine print_residuals - This function displays the residual norms of 
 *    each Ritz vector computed at this iteration.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * ritzValues     The Ritz values
 * blockNorms     Residual norms of the corresponding Ritz vectors
 * numConverged   Number of already converged epairs
 * numLocked      Number of locked (converged) epairs
 * iev            indicates which eigenvalue each block vector corresponds to
 * blockSize      Number of pairs in the block
 * primme         Structure containing various solver parameters
 ******************************************************************************/

static void print_residuals(double *ritzValues, double *blockNorms,
   int numConverged, int numLocked, int *iev, int blockSize, 
   primme_params *primme) {

   int i;  /* Loop variable */
   int found;  /* Loop variable */

   if (primme->printLevel >= 3 && primme->procID == 0) {

      if (primme->locking) 
         found = numLocked;
      else 
         found = numConverged;

      for (i=0; i < blockSize; i++) {
         fprintf(primme->outputFile, 
            "OUT %d conv %d blk %d MV %d Sec %E EV %13E |r| %.3E\n",
         primme->stats.numOuterIterations, found, i, primme->stats.numMatvecs,
         primme_wTimer(0), ritzValues[iev[i]], blockNorms[i]);
      }

      fflush(primme->outputFile);
   }

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

static void switch_from_JDQMR(primme_CostModel *model, primme_params *primme) {

   int switchto=0, one=1;
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
         (*primme->globalSumDouble)(&ratio, &globalRatio, &one, primme); 
         ratio = globalRatio/primme->numProcs;
      }

      if (ratio > 1.05) { 
         /* Always use GD+k. No further model updates */
         primme->dynamicMethodSwitch = -1;
         primme->correctionParams.maxInnerIterations = 0;
         if (primme->printLevel >= 3 && primme->procID == 0) 
            fprintf(primme->outputFile, 
            "Ratio: %e Switching permanently to GD+k\n", ratio);
         return;
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
      (*primme->globalSumDouble)(&ratio, &globalRatio, &one, primme); 
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

   return;

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
static void switch_from_GDpk(primme_CostModel *model, primme_params *primme) {

   int switchto=0, one = 1;
   double ratio, globalRatio;

   /* if no restart has occurred (only possible under dyn=3) current timings */
   /* do not include restart costs. Remain with GD+k until a restart occurs */
   if (primme->stats.numRestarts == 0) return;

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
      return;
   }

   /* ------------------------------------------------------------------- *
    * Compute the ratio of expected times JDQMR/GD+k. To switch to JDQMR, the
    * ratio must be < 0.95. Update accum_jdq_gdk for recommendation to user
    * ------------------------------------------------------------------- */

   ratio = ratio_JDQMR_GDpk(model, 0, model->JDQMR_slowdown, 
                                  model->ratio_MV_outer);

   /* If more many procs, make sure that all have the same ratio */
   if (primme->numProcs > 1) {
      (*primme->globalSumDouble)(&ratio, &globalRatio, &one, primme); 
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

   return;
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
 *    For more details see papers [1] and [2] (listed in primme_@(pre).c)
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
