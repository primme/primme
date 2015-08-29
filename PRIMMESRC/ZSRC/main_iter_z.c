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
#include "primme.h"
#include "const.h"
#include "wtime.h"
#include "main_iter_z.h"
#include "main_iter_private_z.h"
#include "convergence_z.h"
#include "correction_z.h"
#include "init_z.h"
#include "ortho_z.h"
#include "restart_z.h"
#include "locking_z.h"
#include "solve_H_z.h"
#include "update_projection_z.h"
#include "update_W_z.h"
#include "numerical_z.h"

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
 * realWork Complex_Z work array
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

int main_iter_zprimme(double *evals, int *perm, Complex_Z *evecs, 
   double *resNorms, double machEps, int *intWork, void *realWork, 
   primme_params *primme) {
         
   int i;                   /* Loop variable                                 */
   int blockSize;           /* Current block size                            */
   int AvailableBlockSize;  /* There is enough space in basis for this block */
   int ievMax;              /* Index of next eigenvalue to be approximated   */
   int basisSize;           /* Current size of the basis V                   */
   int numLocked;           /* Number of locked Ritz vectors                 */
   int numGuesses;          /* Current number of initial guesses in evecs    */
   int nextGuess;           /* Index of the next initial guess in evecs      */ 
   int numConverged;        /* Number of converged Ritz pairs                */
   int recentlyConverged;   /* Number of target Ritz pairs that have         */
                            /*    converged during the current iteration.    */
   int numConvergedStored;  /* Numb of Ritzvecs temporarily stored in evecs  */
                            /*    to allow for skew projectors w/o locking   */
   int converged;           /* True when all required Ritz vals. converged   */
   int LockingProblem;      /* Flag==1 if practically converged pairs locked */
   int restartLimitReached; /* True when maximum restarts performed          */
   int numPrevRetained;     /* Number of vectors retained using recurrence-  */
                            /* based restarting.                             */
   int maxEvecsSize;        /* Maximum capacity of evecs array               */
   int doubleSize;          /* sizeof the three double arrays hVals,         */
                            /*                      prevRitzVals, blockNorms */
   int rworkSize;           /* Size of rwork array                           */
   int numPrevRitzVals = 0; /* Size of the prevRitzVals updated in correction*/
   int ret;                 /* Return value                                  */

   int *iwork;              /* Integer workspace pointer                     */
   int *flag;               /* Inidicates which Ritz values have converged   */
   int *ipivot;             /* The pivot for the UDU factorization of M      */
   int *iev;                /* Evalue index each block vector corresponds to */
   int ONE = 1;             /* To be passed by reference in matrixMatvec     */

   double largestRitzValue; /* The largest modulus of any Ritz value computed*/
   double tol;              /* Required tolerance for residual norms         */
   double maxConvTol;       /* Max locked residual norm (see convergence.c)  */
   Complex_Z *V;               /* Basis vectors                               */
   Complex_Z *W;               /* Work space storing A*V                      */
   Complex_Z *H;               /* Upper triangular portion of V'*A*V          */
   Complex_Z *M;               /* The projection Q'*K*Q, where Q = [evecs, x] */
                            /* x is the current Ritz vector and K is a       */
                            /* hermitian preconditioner.                     */
   Complex_Z *UDU;             /* The factorization of M=Q'KQ                 */
   Complex_Z *evecsHat;       /* K^{-1}evecs                                   */
   Complex_Z *rwork;          /* Real work space.                              */
   Complex_Z *hVecs;          /* Eigenvectors of H                             */
   Complex_Z *previousHVecs;   /* Coefficient vectors retained by            */
                            /* recurrence-based restarting                   */
   double *hVals;           /* Eigenvalues of H                              */
   double *prevRitzVals;    /* Eigenvalues of H at previous outer iteration  */
                            /* by robust shifting algorithm in correction.c  */
   double *blockNorms;      /* Residual norms corresponding to current block */
                            /* vectors.                                      */
   Complex_Z tpone = {+1.0e+00,+0.0e00};/* constant 1.0 of type Complex_Z */
   Complex_Z tzero = {+0.0e+00,+0.0e00};/* constant 0.0 of type Complex_Z */

   /* Runtime measurement variables for dynamic method switching             */
   primme_CostModel CostModel; /* Structure holding the runtime estimates of */
                            /* the parameters of the model.Only visible here */
   double timeForMV;        /* Measures time for 1 matvec operation          */
   double tstart;           /* Timing variable for accumulative time spent   */

   /* -------------------------------------------------------------- */
   /* Subdivide the workspace                                        */
   /* -------------------------------------------------------------- */

   maxEvecsSize = primme->numOrthoConst + primme->numEvals;

   V             = (Complex_Z *) realWork;
   W             = V + primme->nLocal*primme->maxBasisSize;
   H             = W + primme->nLocal*primme->maxBasisSize;
   hVecs         = H + primme->maxBasisSize*primme->maxBasisSize;
   previousHVecs = hVecs + primme->maxBasisSize*primme->maxBasisSize;
   if (! (primme->correctionParams.precondition && 
          primme->correctionParams.maxInnerIterations != 0 &&
          primme->correctionParams.projectors.RightQ &&
          primme->correctionParams.projectors.SkewQ           ) ) {
      evecsHat   = NULL;
      M          = NULL;
      UDU        = NULL;
      rwork      = previousHVecs + primme->restartingParams.maxPrevRetain*
                           primme->maxBasisSize;
   }
   else {
      evecsHat   = previousHVecs + primme->restartingParams.maxPrevRetain*
                           primme->maxBasisSize;
      M          = evecsHat + primme->nLocal*maxEvecsSize;
      UDU        = M + maxEvecsSize*maxEvecsSize; 
      rwork      = UDU + maxEvecsSize*maxEvecsSize; 
   }
   /* Size of three double arrays that go at the end */
   doubleSize    = (2*primme->maxBasisSize + primme->numEvals + 
                                primme->maxBlockSize)*sizeof(double);
   /* rwork size is (Total - endarrays) - frontarrays */
   rworkSize     = (primme->realWorkSize - doubleSize)/sizeof(Complex_Z);
   rworkSize     = rworkSize - (rwork-V);

   hVals         = (double *)(rwork + rworkSize);
   prevRitzVals  = hVals + primme->maxBasisSize;
   blockNorms    = prevRitzVals + (primme->numEvals+primme->maxBasisSize);

   /* Integer workspace */

   flag = intWork;
   iev = flag + primme->maxBasisSize;
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
   blockSize = primme->maxBlockSize; 
   ievMax = primme->maxBlockSize;      

   for (i=0; i < primme->maxBlockSize; i++) {
      iev[i] = i;
   }

   /* ---------------------------------------- */
   /* Set the tolerance for the residual norms */
   /* ---------------------------------------- */

   largestRitzValue = 0.0L;
   if (primme->aNorm > 0.0L) {
      tol = primme->eps*primme->aNorm;
   }
   else {
      tol = primme->eps; /* tol*largestRitzValue will be checked */
   }
   maxConvTol = tol;

   /* -------------------------------------- */
   /* Quick return for matrix of dimension 1 */
   /* -------------------------------------- */

   if (primme->n == 1) {
      evecs[0] = tpone;
      (*primme->matrixMatvec)(&evecs[0], W, &ONE, primme);
      evals[0] = W[0].r;
      V[0] = tpone;

      resNorms[0] = 0.0L;
      primme->stats.numMatvecs++;
      return 0;
   }

   /* -------------------- */
   /* Initialize the basis */
   /* -------------------- */

   ret = init_basis_zprimme(V, W, evecs, evecsHat, M, UDU, ipivot, machEps,
           rwork, rworkSize, &basisSize, &nextGuess, &numGuesses, &timeForMV,
           primme);

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
      reset_flags_zprimme(flag, 0, primme->maxBasisSize-1);

      /* Compute the initial H and solve for its eigenpairs */
   
      update_projection_zprimme(V, W, H, 0,primme->maxBasisSize,basisSize,
         hVecs,primme);
      ret = solve_H_zprimme(H, hVecs, hVals, basisSize, primme->maxBasisSize,
         &largestRitzValue, numLocked, rworkSize, rwork, iwork, primme);

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

            /* Adjust the block size if necessary. Remember the available for */
            /* expansion slots in the basis, as blockSize may be reduced later*/

            adjust_blockSize(iev, flag, &blockSize, primme->maxBlockSize, 
               &ievMax, basisSize, primme->maxBasisSize, numLocked, 
               numConverged, primme->numEvals, primme->n);
            AvailableBlockSize = blockSize;

            /* Check the convergence of the blockSize Ritz vectors computed */

            recentlyConverged = check_convergence_zprimme(V, W, hVecs, 
               hVals, flag, basisSize, iev, &ievMax, blockNorms, &blockSize, 
               numConverged, numLocked, evecs, tol, maxConvTol, 
               largestRitzValue, rwork, primme);

            /* If the total number of converged pairs, including the     */
            /* recentlyConverged ones, are greater than or equal to the  */
            /* target number of eigenvalues, attempt to restart, verify  */
            /* their convergence, lock them if necessary, and return.    */
            /* For locking interior, restart and lock now any converged. */

            numConverged += recentlyConverged;

            if (numConverged >= primme->numEvals ||
                (primme->locking && recentlyConverged > 0
                    && primme->target != primme_smallest
                    && primme->target != primme_largest)) {
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
               blockSize = AvailableBlockSize;
               Num_scal_zprimme(blockSize*primme->nLocal, tzero,
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
                        largestRitzValue); 

                     if (ret) switch (primme->dynamicMethodSwitch) {
                        /* for few evals (dyn=1) evaluate GD+k only at restart*/
                        case 3: switch_from_GDpk(&CostModel,primme); break;
                        case 2: case 4: switch_from_JDQMR(&CostModel,primme);
                     } /* of if-switch */
                  } /* of recentlyConv > 0 || dyn==2 */
               } /* dynamic switching */
               /* - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

               ret = solve_correction_zprimme(V, W, evecs, evecsHat, UDU, 
                 ipivot, evals, numLocked, numConvergedStored, hVals, 
                 prevRitzVals, &numPrevRitzVals, flag, basisSize, blockNorms, 
                 iev, blockSize, tol, machEps, largestRitzValue, rwork, iwork, 
                 rworkSize, primme);

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

            ret = ortho_zprimme(V, primme->nLocal, basisSize, 
               basisSize+blockSize-1, evecs, primme->nLocal, 
               primme->numOrthoConst+numLocked, primme->nLocal, primme->iseed, 
               machEps, rwork, rworkSize,primme);

            if (ret < 0) {
               primme_PushErrorMessage(Primme_main_iter, Primme_ortho, ret,
                               __FILE__, __LINE__, primme);
               return ORTHO_FAILURE;
            }
           
            /* Compute W = A*V for the orthogonalized corrections */

            update_W_zprimme(V, W, basisSize, blockSize, primme);

            /* Extend H by blockSize columns and rows and solve the */
            /* eigenproblem for the new H.                          */

            numPrevRetained = retain_previous_coefficients(hVecs, 
               previousHVecs, basisSize, iev, blockSize, primme);
            update_projection_zprimme(V, W, H, basisSize, 
               primme->maxBasisSize, blockSize, hVecs, primme);
            basisSize = basisSize + blockSize;
            ret = solve_H_zprimme(H, hVecs, hVals, basisSize, 
               primme->maxBasisSize, &largestRitzValue, numLocked, rworkSize, 
               rwork, iwork, primme);

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

         basisSize = restart_zprimme(V, W, H, hVecs, hVals, flag, iev, 
            evecs, evecsHat, M, UDU, ipivot, basisSize, numConverged, 
            &numConvergedStored, numLocked, numGuesses, previousHVecs, 
            numPrevRetained, machEps, rwork, rworkSize, primme);

         if (basisSize <= 0) {
            primme_PushErrorMessage(Primme_main_iter, Primme_restart, 
                            basisSize, __FILE__, __LINE__, primme);
            return RESTART_FAILURE;
         }

         /* ----------------------------------------------------------- */
         /* If locking is engaged, then call the lock vectors routine,  */
         /* else mark all non target Ritz values as unconverged.        */
         /* ----------------------------------------------------------- */

         if (primme->locking) {
            ret = lock_vectors_zprimme(tol, &largestRitzValue, &maxConvTol,
               &basisSize, &numLocked, &numGuesses, &nextGuess, V, W, H, 
               evecsHat, M, UDU, ipivot, hVals, hVecs, evecs, evals, perm, 
               machEps, resNorms, &numPrevRitzVals, prevRitzVals, flag, 
               rwork, rworkSize, iwork, &LockingProblem, primme);
            numConverged = primme->initSize = numLocked;

            if (ret < 0) {
               primme_PushErrorMessage(Primme_main_iter, Primme_lock_vectors,
                               ret, __FILE__, __LINE__, primme);
               return LOCK_VECTORS_FAILURE;
            }
            
         }
         else {
            /* reset_flags_zprimme(flag, primme->numEvals, primme->maxBasisSize-1);*/
            check_reset_flags_zprimme(flag, &numConverged, hVals, 
               prevRitzVals, numPrevRitzVals, tol, largestRitzValue, primme);
         }

         primme->stats.numRestarts++;

         /* ------------------------------------------------------------- */
         /* If dynamic method switching == 1, update model parameters and */
         /* evaluate whether to switch from GD+k to JDQMR. This is after  */
         /* restart. GD+k is also evaluated if a pair converges.          */
         /* ------------------------------------------------------------- */
         if (primme->dynamicMethodSwitch == 1 ) {
            tstart = primme_wTimer(0);
            ret = update_statistics(&CostModel, primme, tstart, 0, 1,
               numConverged, blockNorms[0], largestRitzValue); 
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
            if (primme->aNorm <= 0.0L) primme->aNorm = largestRitzValue;
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
            flag, tol, largestRitzValue, rwork, &numConverged, primme);

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

            Num_zcopy_zprimme(primme->nLocal*primme->numEvals, V, 1, 
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
               if (primme->aNorm <= 0.0L) primme->aNorm = largestRitzValue;
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

            ret = ortho_zprimme(V, primme->nLocal,0, basisSize-1, evecs, 
               primme->nLocal, primme->numOrthoConst+numLocked, primme->nLocal,
               primme->iseed, machEps, rwork, rworkSize, primme);
            if (ret < 0) {
               primme_PushErrorMessage(Primme_main_iter, Primme_ortho, ret,
                               __FILE__, __LINE__, primme);
               return ORTHO_FAILURE;
            }
            update_W_zprimme(V, W, 0, basisSize, primme);

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

   if (primme->aNorm <= 0.0L) primme->aNorm = largestRitzValue;

   return 0;

}

/******************************************************************************
          Some basic functions within the scope of main_iter
*******************************************************************************/

/*******************************************************************************
 * Subroutine adjust_blockSize - This subroutine increases or decreases the 
 *    block size depending on how many vectors are available to 
 *    put into the block.  
 * 
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * iev             indicates which eigenvalue each block vector corresponds to
 * flag            array of flags indicating which eigenvalues have converged
 * basisSize       current size of the basis
 * maxBasisSize    maximum allowed size of the basis
 * maxBlockSize    maximum allowed size of the block
 * numLocked       number of locked evals (if locking)
 * numConverged    number of converged evals (if locking=numLocked)
 * numWantedEvs    number of eigenvalues saught
 * matrixDimension the dimension of the large eigenproblem matrix 
 * 
 * INPUT/OUTPUT PARAMETERS
 * -----------------------
 * blockSize  The size of the block.
 * ievMax     index of next vector to be targeted by the solver
 ******************************************************************************/

static void adjust_blockSize(int *iev, int *flag, int *blockSize, 
   int maxBlockSize, int *ievMax, int basisSize, int maxBasisSize, 
   int numLocked, int numConverged, int numWantedEvs, int matrixDimension) {

   /* If the block size is larger than the number of vacancies in V, */
   /* reduce the block size else, if the blockSize is smaller than   */
   /* maxBlockSize, try to find unconverged vectors to increase the  */
   /* block size. Do not exceed matrix dim, or target more than Evs  */

   *blockSize = 0;
   *ievMax = 0;

   while ((*blockSize < maxBlockSize) && (*blockSize < basisSize) && 
      (*blockSize < (maxBasisSize - basisSize)) &&
      (*blockSize+numConverged < numWantedEvs + 1) &&   
      (*ievMax < basisSize) &&
      (*blockSize < (matrixDimension - numLocked - basisSize)) ){
   
      if (flag[*ievMax] == UNCONVERGED) {
         *blockSize = *blockSize + 1;
         iev[*blockSize-1] = *ievMax;
      }

      *ievMax = *ievMax + 1;

   }
}

/*******************************************************************************
 * Function retain_previous_coefficients - This function is part of the
 *    recurrence-based restarting method for accelerating convergence.
 *    This function is called one iteration before restart so that the
 *    coefficients (eigenvectors of the projection H) corresponding to 
 *    a few of the target Ritz vectors may be retained at restart. 
 *    The desired coefficients are copied to a seperate storage space so
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

static int retain_previous_coefficients(Complex_Z *hVecs, Complex_Z *previousHVecs, 
   int basisSize, int *iev, int blockSize, primme_params *primme) {

   int i, j;            /* Loop indices                                  */
   int index;           /* The index of some coefficient vector in hVecs */ 
   int numPrevRetained; /* The number of coefficent vectors retained     */
   Complex_Z tzero = {+0.0e+00,+0.0e00};

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
            Num_zcopy_zprimme(basisSize, &hVecs[basisSize*index], 1, 
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
 * Function check_reset_flags - This subroutine is called only in soft-locking
 *    after restart. It mainly calls reset_flags for the unwanted evals. But
 *    it also checks if any of the previous flagged converged eigenvalues seems
 *    to have become unconverged by checking hVals[i]-prevRitzVals[i] < tol
 *    If not, it flags it UNCONVERGED and lets it be targeted again. This avoids 
 *    early converged but unwanted evs preventing wanted from being targeted.
 ******************************************************************************/

void check_reset_flags_zprimme(int *flag, int *numConverged, 
   double *hVals, double *prevRitzVals, int numPrevRitzVals,
   double tol, double aNormEstimate, primme_params *primme) {

   int i;

   reset_flags_zprimme(flag, primme->numEvals, primme->maxBasisSize-1);

   if (primme->aNorm <= 0.0L) {
      tol = tol * aNormEstimate;
   }
   for (i=0;i<primme->numEvals;i++) {
      if (i >= numPrevRitzVals) break;
      if ((flag[i] == CONVERGED) && (fabs(hVals[i]-prevRitzVals[i]) > tol)) {
         (*numConverged)--;
         flag[i] = UNCONVERGED;
      }
   }

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
   
static int verify_norms(Complex_Z *V, Complex_Z *W, Complex_Z *hVecs, 
   double *hVals, int basisSize, double *resNorms, int *flag, double tol, 
   double aNormEstimate, void *rwork, int *numConverged, primme_params *primme){

   int i;         /* Loop varible                                      */
   int converged; /* True when all requested Ritz values are converged */
   int nev, n;    /* convenience integers for numEvals and nLocal      */
   double *dwork = (double *) rwork; /* pointer to cast rwork to double*/
   Complex_Z ztmp;  /* temp complex var */

   nev = primme->numEvals;
   n   = primme->nLocal;

   /* Set up the tolerance if necessary */

   if (primme->aNorm <= 0.0L) {
      tol = tol * aNormEstimate;
   }

   /* Compute the residual vectors */

   for (i=0; i < nev; i++) {
      {ztmp.r = -hVals[i]; ztmp.i = 0.0L;}
      Num_axpy_zprimme(n, ztmp, &V[n*i], 1, &W[n*i], 1);
      ztmp = Num_dot_zprimme(n, &W[n*i], 1, &W[n*i], 1);
      dwork[nev+i] = ztmp.r;

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

   int switchto, one=1;
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

   int switchto, one = 1;
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
 *    For more details see papers [1] and [2] (listed in primme_z.c)
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
