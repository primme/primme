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
 * File: locking.c
 *
 * Purpose - This file contains routines for locking converged Ritz pairs.
 *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include "primme.h"
#include "wtime.h"
#include "const.h"
#include "locking_d.h"
#include "locking_private_d.h"
#include "ortho_d.h"
#include "update_projection_d.h"
#include "update_W_d.h"
#include "solve_H_d.h"
#include "update_W_d.h"
#include "restart_d.h"
#include "factorize_d.h"
#include "numerical_d.h"
#include <assert.h>

/******************************************************************************
 * Function lock_vectors - This subroutine locks converged Ritz pairs.  The
 *    converged Ritz vectors are copied to evecs, and the converged Ritz values
 *    are copied to evals.  The evals array is maintained in sorted order; 
 *    however the evecs array is not.  Instead, a permutation array is 
 *    maintained so that the locked Ritz vectors may be sorted upon return to
 *    the user.  
 *
 *    This routine assumes that it is called after return from the restart
 *    routine and that the converged Ritz vectors are at the end of the basis
 *    V.  This way they can easily be replaced by new initial guesses.  Before
 *    this can be done though, the vectors must be checked for convergence (or 
 *    practical convergence) a second time. Because Ritz vectors are not locked
 *    when they first converge but only after restart, we need to make sure 
 *    that they have not become unconverged (or not practically converged).
 *
 *    Vectors that have remain converged are locked to the evecs array and
 *    are replaced by initial guesses if there are any remaining.  The initial
 *    guesses are orthogonalized.  If any vectors were locked, the eigenproblem
 *    for the projected matrix H is then solved for the modified basis.
 *
 *    
 * Input arrays and parameters
 * ---------------------------
 * tol        Convergence tolerance for residual norms
 *
 * aNormEsimate if (primme->aNorm<=0), use tol*aNormEstimate (=largestRitzValue)
 *
 * maxConvTol The maximum residual norm > tol, of any locked Ritz vectors.
 *            Vectors with ||r|| > tol can be locked if an accuracy problem 
 *            has been detected in check_convergence(). Accuracy is not lost
 *            however, because the missing components can be recovered from 
 *            the rest of the locked vectors with a Rayleigh Ritz at the end.
 *
 * machEps    Double machine precision
 *
 * rwork      Real work array
 *
 * rworkSize  Size of rwork.  Must be at least: 
 *              2*maxBasisSize + MAX( orthoSize, 3*maxBasisSize, 
 *                (primme->numOrthoConst+primme->numEvals)*maxBasisSize)
 *
 * iwork      Integer work array of size maxBasisSize used as perm in solve_H
 *
 * primme       Structure containing various solver parameters
 *
 * 
 * Input/Output parameters
 * -----------------------
 * basisSize  The size of the basis V
 *
 * numLocked  The number of vectors that have been locked this far
 *
 * numGuesses The number of initial guesses remaining in evecs
 *
 * nextGuess  The index of the next initial guess in the evecs array
 * 
 * V          The basis (w/o converged and possibly with new guesses)
 *
 * W          A*V
 * 
 * H          The projection V'*A*V
 *
 * evecsHat   The K*{-1}*evecs updated for the newly locked vectors (if needed)
 *
 * M          the matrix evecs'*evecsHat
 *
 * UDU        the factorization of the matrix M 
 *
 * ipivot     the pivot array for the UDU factorization
 *
 * hVals      The eigenvalues of H (the Ritz values)
 *
 * hVecs      The eigenvectors of H (the Ritz vectors)
 *
 * evecs      Contains the initial guesses and stores the locked vectors
 *
 * evals      Contains the locked Ritz values
 *
 * perm       A permutation array that maps locked Ritz vectors to the
 *            sorted order of the evals array
 *
 * resNorms   Residual norms of the locked Ritz vectors
 *
 * numPrevRitzVals  Number of Ritz values from iteration (current-1) 
 *
 * prevRitzVals The Ritz values of iteration (current-1) closest to locked ones 
 *            will be removed.
 *
 * LockingProblem Set to 1 if some evector "practically" converged, ie, some of 
 *            its components are in evecs. A Rayleigh-Ritz may be needed.
 *
 * flag       Flag array indicating the status of the restart vectors:
 *               UNCONVERGED  : A Ritz vector that is unconverged
 *               LOCK_IT      : A Ritz vector that is to be locked
 *               UNCONDITIONAL_LOCK_IT : A Ritz vector that is to be locked
 *                              without checking its norm.
 *               INITIAL_GUESS: A new initial guess that must be orthogonalized
 *
 * Return value
 * ------------
 * int  Error code: 0 if locking occured successfully 
 *                 -1 if orthogonalization of initial guesses failed
 *                 -2 if a block Krylov basis could not be generated
 *                 -3 if the call to solve_H failed
 *                 -4 if the call to UDUDecompose failed
 *
 ******************************************************************************/

int lock_vectors_dprimme(double tol, double *aNormEstimate, double *maxConvTol, 
   int *basisSize, int *numLocked, int *numGuesses, int *nextGuess,
   double *V, double *W, double *H, double *evecsHat, double *M, 
   double *UDU, int *ipivot, double *hVals, double *hVecs, 
   double *evecs, double *evals, int *perm, double machEps, 
   double *resNorms, int *numPrevRitzVals, double *prevRitzVals, 
   int *flag, double *rwork, int rworkSize, int *iwork, 
   int *LockingProblem, primme_params *primme) {

   int i;             /* Loop counter                                       */
   int numCandidates; /* Number of targeted Ritz vectors converged before   */
                      /* restart.                                           */
   int newStart;      /* Index in evecs where the locked vectors were added */
   int numNewVectors; /* Number of vectors added to the basis to replace    */
                      /* locked vectors.                                    */
   int candidate;     /* Index of Ritz vector to be checked for convergence */
   int numDeflated;   /* The number of vectors actually locked              */
   int numReplaced;   /* The number of locked vectors that were replaced by */
                      /* initial guesses.                                   */
   int numRecentlyLocked; /* Number of vectors locked.                      */
   int evecsSize;     /* The number of orthogonalization constraints plus   */
                      /* the number of locked vectors.                      */
   int ret;           /* Used to store return values.                       */
   int workinW;       /* Flag whether an active W vector is used as tempwork*/
   int entireSpace = (*basisSize+*numLocked >= primme->n); /* bool if entire*/
                      /* space is built, so current ritzvecs are accurate.  */

   double *norms, *tnorms; /* Array of residual norms, and temp array       */
   double attainableTol;   /* Used to verify a practical convergence problem*/
   double *residual;  /* Stores residual vector                             */
   double ztmp;       /* temp variable */

   /* ----------------------------------------*/
   /* Assign temporary work space for residual*/
   /* ----------------------------------------*/

   if (*basisSize < primme->maxBasisSize) {
      /* compute residuals in the next open slot of W */
      residual = &W[*basisSize*primme->nLocal];
      workinW = 0;
   }
   else {
      /* This basiSize==maxBasisSize, immediately after restart, can only occur
       * if the basisSize + numLocked = n (at which case we lock everything)
       * OR if (numConverged + restartSize + numPrevRetain > basisSize), ie.
       * too many converged. Since we do not know which evec will be locked
       * we use the W[LAST] as temporary space, but only after W[LAST] has 
       * been used to compute residual(LAST) -the while loop starts from LAST.
       * After all lockings, if the LAST evec was not locked, we must  
       * recompute W[LAST]=Av. This matvec event is extremely infrequent */
      residual = &W[(*basisSize-1)*primme->nLocal];
      workinW = 1;
   }

   /* -------------------------------------*/
   /* Set the tolerance, and attainableTol */
   /* -------------------------------------*/
   
   if (primme->aNorm <= 0.0L) {
      tol = tol * (*aNormEstimate);
   }
   attainableTol=max(tol,sqrt(primme->numOrthoConst+*numLocked)*(*maxConvTol));

   /* -------------------------------------------------------- */
   /* Determine how many Ritz vectors converged before restart */
   /* -------------------------------------------------------- */

   i = *basisSize - 1;
   while ((flag[i] == LOCK_IT ||flag[i] == UNCONDITIONAL_LOCK_IT) && i >= 0) {
      i--;
   }
      
   numCandidates = *basisSize - i - 1;

   if (numCandidates == 0) {
      return 0;
   }

   /* --------------------------------- */
   /* Compute residuals and their norms */
   /* --------------------------------- */

   tnorms = (double *) rwork;
   norms  = tnorms + numCandidates;

   for (i = *basisSize-1, candidate = numCandidates-1;  
      i >= *basisSize-numCandidates; i--, candidate--) {
      Num_dcopy_dprimme(primme->nLocal, &W[primme->nLocal*i], 1, residual, 1);
      ztmp = -hVals[i];
      Num_axpy_dprimme(primme->nLocal, ztmp, &V[primme->nLocal*i],1,residual,1);
      tnorms[candidate] = Num_dot_dprimme(primme->nLocal,residual,1,residual,1);
   }

   /* Global sum the dot products */
   (*primme->globalSumDouble)(tnorms, norms, &numCandidates, primme); 

   numRecentlyLocked = 0;

   /* ------------------------------------------------------------------- */
   /* Check the convergence of each residual norm.  If the Ritz vector is */
   /* converged, then lock it.                                            */
   /* ------------------------------------------------------------------- */

   for (i = *basisSize - numCandidates, candidate = 0; i < *basisSize; i++, 
      candidate++) {

      norms[candidate] = sqrt(norms[candidate]);

      /* If the vector has become (regularly or practically) unconverged, */
      /* then flag it, else lock it and replace it with an initial guess, */
      /* if one is available. Exception: If the entire space is spanned,  */
      /* we can't do better, so lock it.                                  */


      if ((flag[i]!=UNCONDITIONAL_LOCK_IT && norms[candidate] >= tol 
                                   && !entireSpace ) ||
          (flag[i]==UNCONDITIONAL_LOCK_IT && norms[candidate] >= attainableTol
                                   && !entireSpace )) {
         flag[i] = UNCONVERGED;
      }
      else {
        /* If an unconditional lock has become converged, show it and */
        /* record the max converged tolerance accordingly             */
        if (norms[candidate]<tol) {
           flag[i]=LOCK_IT;
           *maxConvTol = max(*maxConvTol, tol);
        }
        else {
           *maxConvTol = max(*maxConvTol, norms[candidate]);
           *LockingProblem = 1;
        }

         if (primme->printLevel >= 2 && primme->procID == 0) { 
            fprintf(primme->outputFile, 
            "Lock epair[ %d ]= %e norm %.4e Mvecs %d Time %.4e Flag %d\n",
                  *numLocked+1, hVals[i], norms[candidate], 
                   primme->stats.numMatvecs,primme_wTimer(0),flag[i]);
            fflush(primme->outputFile);
         }

         /* Copy the converged Ritz vector to the evecs array and  */
         /* insert the converged Ritz value in sorted order within */
         /* the evals array.                                       */

         Num_dcopy_dprimme(primme->nLocal, &V[primme->nLocal*i], 1, 
            &evecs[primme->nLocal*(primme->numOrthoConst + *numLocked)], 1);
         insertionSort(hVals[i], evals, norms[candidate], resNorms, perm, 
            *numLocked, primme);

         /* If there are any initial guesses remaining, then copy it */
         /* into the basis, else flag the vector as locked so it may */
         /* be discarded later.                                      */

         if (*numGuesses > 0) {
            Num_dcopy_dprimme(primme->nLocal, 
               &evecs[primme->nLocal*(*nextGuess)], 1, &V[primme->nLocal*i], 1);
            flag[i] = INITIAL_GUESS;
            *numGuesses = *numGuesses - 1;
            *nextGuess = *nextGuess + 1;
         }
         else {
            flag[i] = LOCKED;
         }

         *numLocked = *numLocked + 1;
         numRecentlyLocked++;
            
      }
   }
      
   evecsSize = primme->numOrthoConst + *numLocked;

   /* -------------------------------------------------------------------- */
   /* If a W vector was used as workspace for residual AND its evec has    */
   /* not been locked out,  recompute it, W = A*v. This is rare.           */
   /* -------------------------------------------------------------------- */
      if (workinW && flag[*basisSize-1] != LOCKED) {
         update_W_dprimme(V, W, *basisSize-1, 1, primme);
      }
   
   /* -------------------------------------------------------------------- */
   /* Return IF all target Ritz vectors have been locked, ELSE update the  */
   /* evecsHat array by applying the preconditioner (if preconditioning is */
   /* needed, and JDQMR with right, skew Q projector is applied            */
   /* -------------------------------------------------------------------- */

   if (*numLocked >= primme->numEvals) {
      return 0;
   }
   else if (UDU != NULL) {

      /* Compute K^{-1}x for all newly locked eigenvectors */

      newStart = primme->nLocal*(evecsSize - numRecentlyLocked);
      (*primme->applyPreconditioner)( &evecs[newStart], &evecsHat[newStart], 
                                    &numRecentlyLocked, primme);
      primme->stats.numPreconds += numRecentlyLocked;

      /* Update the projection evecs'*evecsHat now that evecs and evecsHat   */
      /* have been expanded by numRecentlyLocked columns.  Required          */
      /* workspace is numLocked*numEvals.  The most ever needed would be     */
      /* maxBasisSize*numEvals.                                              */

      update_projection_dprimme(evecs, evecsHat, M, 
         evecsSize-numRecentlyLocked, primme->numOrthoConst+primme->numEvals, 
         numRecentlyLocked, rwork, primme);

      ret = UDUDecompose_dprimme(M, UDU, ipivot, evecsSize, rwork, 
         rworkSize, primme);

      if (ret != 0) {
         primme_PushErrorMessage(Primme_lock_vectors, Primme_ududecompose, ret,
            __FILE__, __LINE__, primme);
         return UDUDECOMPOSE_FAILURE;
      }

   }


   /* --------------------------------------------------------------------- */
   /* Swap, towards the end of the basis, vectors that were locked but not  */
   /* replaced by new initial guesses.                                      */
   /* --------------------------------------------------------------------- */

   numDeflated = swap_flagVecs_toEnd(*basisSize, LOCKED, V, W, H, hVals, flag, 
      primme);

   /* --------------------------------------------------------------------- */
   /* Reduce the basis size by numDeflated and swap the new initial guesses */
   /* towards the end of the basis.                                         */
   /* --------------------------------------------------------------------- */
  
   numReplaced = swap_flagVecs_toEnd(*basisSize-numDeflated, INITIAL_GUESS, 
      V, W, H, hVals, flag, primme);

   *basisSize = *basisSize - (numDeflated + numReplaced);

   if (primme->printLevel >= 5 && primme->procID == 0) {
      fprintf(primme->outputFile, "numDeflated: %d numReplaced: %d \
              basisSize: %d\n", numDeflated, numReplaced, *basisSize);
   }

   /* ---------------------------------------------------------------------- */
   /* If there are new initial guesses, then orthogonalize them and update W */ 
   /* ---------------------------------------------------------------------- */

   if (numReplaced > 0) {
      ret = ortho_dprimme(V, primme->nLocal, *basisSize, 
         *basisSize+numReplaced-1, evecs, primme->nLocal, evecsSize, 
         primme->nLocal, primme->iseed, machEps, rwork, rworkSize, primme);

      if (ret < 0) {
         primme_PushErrorMessage(Primme_lock_vectors, Primme_ortho, ret, 
                         __FILE__, __LINE__, primme);
         return ORTHO_FAILURE;
      }   

      update_W_dprimme(V, W, *basisSize, numReplaced, primme);
   }


   numNewVectors = numReplaced;

   if (primme->printLevel >= 5 && primme->procID == 0) {
      fprintf(primme->outputFile, "Number of new vectors: %d\n", numNewVectors);
   }

   /* ---------------------------------------------------------------- */
   /* If new vectors were added to the basis, then extend the rows and */
   /* columns of H by numNewVectors.                                   */
   /* ---------------------------------------------------------------- */

   if (numNewVectors > 0) {
      update_projection_dprimme(V, W, H, *basisSize, primme->maxBasisSize, 
         numNewVectors, hVecs, primme);
      *basisSize = *basisSize + numNewVectors;
   }

   /* ----------------------------------------------------------------- */
   /* Because vectors have been removed from the basis and possibly new */
   /* ones have been added, we must solve the eigenproblem for H.       */
   /* ----------------------------------------------------------------- */

   ret = solve_H_dprimme(H, hVecs, hVals, *basisSize, primme->maxBasisSize,
      aNormEstimate, *numLocked, rworkSize, rwork, iwork, primme);
   reset_flags_dprimme(flag, 0, primme->maxBasisSize - 1);

   if (ret < 0) {
      primme_PushErrorMessage(Primme_lock_vectors, Primme_solve_h, ret, 
                      __FILE__, __LINE__, primme);
      return SOLVE_H_FAILURE;
   } 

   /* ----------------------------------------------------------------- */
   /* Remove locked evals from prevRitzVals. This is difficult since    */
   /* some may have converged out of order. Also, if numReplaced>0 they */
   /* don't match hVals. Their role is not critical (Olsen eps), so     */
   /* easier to assume the first ones (1:numRecentlyLocked) are locked. */
   /* ----------------------------------------------------------------- */

   if (numRecentlyLocked>0) {
      *numPrevRitzVals = *numPrevRitzVals-numRecentlyLocked ;
      for(i=0;i<*numPrevRitzVals;i++)
         prevRitzVals[i] = prevRitzVals[i+numRecentlyLocked];
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
 * perm     The permutation array indicating each Ritz values original
 *          unsorted position.
 *
 ******************************************************************************/

static void insertionSort(double newVal, double *evals, double newNorm,
   double *resNorms, int *perm, int numLocked, primme_params *primme) {

   int i, current; /* Indices used for sorting */
   double ithShift, currentShift;

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
   }

   /* Insert the new value */

   evals[i] = newVal;
   resNorms[i] = newNorm;
   perm[i] = numLocked;
   return;

}

/******************************************************************************
 * Function swap_flagVecs_toEnd - This function swaps Ritz vectors with flag 
 *    value flagValue toward the end of the basis V.  This allows locked vectors
 *    to be discarded from V and initial guesses to be orthogonalized more
 *    efficiently by keeping them contiguous in memory.
 *
 * 
 * Input parameters
 * ----------------
 * basisSize   The current size of the basis V
 *
 * flagValue   Ritz vectors with this flag value will be swapped
 *
 * 
 * Input/Output parameters
 * -----------------------
 * V           The basis
 *
 * W           A*V
 *
 * H           The projection matrix V'*A*V
 *
 * hVals       The eigenvalues of H
 *
 * flag        Values indicating the state of each Ritz vector
 *
 *
 * Return value
 * ------------
 * int   The number of vectors whose flag value is flagValue
 ******************************************************************************/

static int swap_flagVecs_toEnd(int basisSize, int flagValue, double *V, 
  double *W, double *H, double *hVals, int *flag, primme_params *primme) {

   int left, right; /* Search indices                                   */
   int numFlagged;  /* Number of Ritz vectors with flag value flagValue */
   int itemp;       /* Temporary value used for swapping                */
   double dtemp;    /* Temporary value used for swapping                */
   double ztmp;    /* Temporary value used for swapping                */
  
   right = basisSize - 1;
   numFlagged = 0;

   /* Search for values that have flag value flagValue and swap */
   /* them towards the end of the basis.                        */
 
   while (right > 0) {

      /* Find a Ritz vector that doesn't have the necessary flag value */
      while (right >= 0 && flag[right] == flagValue) {
         numFlagged++;
         right--;
      }

      /* There are no vectors to swap with, so swapping is complete */
      if (right <= 0) {
         return numFlagged;
      }

      left = right - 1;

      /* Find a Ritz vector whose flag value is flagValue */

      while (left >= 0 && flag[left] != flagValue) { 
         left--;
      }

      if (left < 0) {
         return numFlagged;
      }

      /* Swap the two columns of V and W */

      Num_swap_dprimme(primme->nLocal, &V[primme->nLocal*left], 1, 
                                       &V[primme->nLocal*right], 1);
      Num_swap_dprimme(primme->nLocal, &W[primme->nLocal*left], 1, 
                                       &W[primme->nLocal*right], 1);

      /* Swap Ritz values */

      dtemp = hVals[left];
      hVals[left] = hVals[right];
      hVals[right] = dtemp;

      /* After restarting, the eigenvectors of H are the standard */
      /* basis vectors (H is diagonal).  Thus, they don't need to */
      /* be swapped.  Just swap the diagonal elements of H.       */
      
      ztmp = H[primme->maxBasisSize*left+left];
      H[primme->maxBasisSize*left+left] = H[primme->maxBasisSize*right+right];
      H[primme->maxBasisSize*right+right] = ztmp;

      itemp = flag[left];
      flag[left] = flag[right];
      flag[right] = itemp;

   }

   return numFlagged;
}
