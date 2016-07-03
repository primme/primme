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
 * Purpose - This file contains the routines to restart with locking.
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
#include "convergence_d.h"
#include "update_projection_d.h"
#include "update_W_d.h"
#include "solve_H_d.h"
#include "restart_d.h"
#include "factorize_d.h"
#include "numerical_d.h"
#include <assert.h>

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
 
int restart_locking_dprimme(int *restartSize, double *V, double *W, 
      int nLocal, int basisSize, int ldV, double **X, double **R,
      double *hVecs, int ldhVecs, int *restartPerm, double *hVals, int *flags,
      int *iev, int *ievSize, double *blockNorms, double *evecs, double *evals,
      int *numConverged, int *numLocked, double *resNorms, int *evecsperm,
      int numPrevRetained, int *indexOfPreviousVecs, int *hVecsPerm,
      int reset, double machEps, double *rwork, int rworkSize, int *iwork, 
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
   double *lockedResNorms;  /* residual norms for locked vectors              */
   int ret, numLocked0 = *numLocked; /* aux variables                         */
   double *blockNorms0;

   /* Return memory requirement */
   if (V == NULL) {
      double t;
      double d;
      return max(max(max(
            /* for permute_vecs and fakeResNorms */
            2*basisSize,
            Num_compute_residual_columns_dprimme(nLocal, NULL, NULL, basisSize, NULL, 0,
               NULL, 0, NULL, primme->maxBlockSize, 0, 0, NULL, 0, NULL,
               primme->maxBlockSize, NULL, 0, NULL, 0, NULL, 0)),
            Num_reset_update_VWXR_dprimme(NULL, NULL, 0, basisSize, 0, NULL,
               *restartSize, 0, NULL,
               &t, 0, *restartSize+*numLocked, 0,
               &t, 0, *ievSize, 0,
               &t, 0, *restartSize, *restartSize+*numLocked, 0,
               &t, 0, *restartSize, 0,
               &t, 0, *ievSize, 0, &d,
               &d, *restartSize, *numLocked,
               0, 0.0, NULL, 0, primme)),
            check_convergence_dprimme(NULL,
               nLocal, 0, NULL, 0, NULL, *numLocked, 0, *restartSize,
               *restartSize+*numLocked, NULL, NULL, NULL, NULL, 0.0, NULL, 0,
               NULL, primme));
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
   /* X & R will have the eigenvectors and residual vectors of the   */
   /* first sizeBlockNorms candidates pairs.                         */
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

   for (i=j=k=numPacked=0; i<basisSize; i++) {
      if (flags[i] != UNCONVERGED
            && (i < primme->numEvals-*numLocked
            /* Refined and prepare_vecs may not completely order pairs        */
            /* considering closest_leq/geq; so we find converged pairs beyond */
            /* the first remaining pairs to converge.                         */
               || primme->target == primme_closest_geq
               || primme->target == primme_closest_leq
            /* Try to converge all pairs in case of overbooking               */
                || overbooking)) {
         restartPerm[left + numPacked++] = i;
      }
      else if (k < left) {
         restartPerm[k++] = i;
      }
      else {
         restartPerm[*numConverged-*numLocked + k++] = i;
      }
   }

   assert(numPacked == *numConverged-*numLocked || overbooking);
 
   permute_vecs_dprimme(hVals, 1, basisSize, 1, restartPerm, (double*)rwork, iwork);
   permute_vecs_dprimme(hVecs, basisSize, basisSize, ldhVecs, restartPerm, rwork, iwork);

   /* -------------------------------------------------------------- */
   /* Restart V and W by replacing it V*hVecs and W*hVecs, copy      */
   /* vectors to be locked into evecs and compute X and R.           */
   /* -------------------------------------------------------------- */

   if (!overbooking) {
      lockedResNorms = &resNorms[*numLocked];
   }
   else {
      lockedResNorms = (double*)rwork;
      rwork += basisSize;
      rworkSize -= basisSize;
   }
   *X = &V[*restartSize*ldV];
   *R = &W[*restartSize*ldV];
   ret = Num_reset_update_VWXR_dprimme(V, W, nLocal, basisSize, ldV, hVecs,
         *restartSize, ldhVecs, hVals,
         V, 0, *restartSize, ldV,
         *X, 0, sizeBlockNorms, ldV,
         evecs, *numLocked+primme->numOrthoConst, left,
               !overbooking?*restartSize:left, primme->nLocal,
         W, 0, *restartSize, ldV,
         *R, 0, sizeBlockNorms, ldV, blockNorms,
         lockedResNorms, left, *restartSize,
         reset, machEps, rwork, rworkSize, primme);
   if (ret != 0) return ret;
 
   /* -------------------------------------------------------------- */
   /* Recompute flags for the vectors to be locked.                  */
   /* NOTE: the eigenvals flagged as practically converged will keep */
   /*       it if their residual norm is still less than             */
   /*       sqrt(numLocked)*tol                                      */
   /* -------------------------------------------------------------- */

   permute_vecs_iprimme(flags, basisSize, restartPerm, iwork);
   ret = check_convergence_dprimme(&V[ldV*left],
         nLocal, ldV, NULL, 0, NULL, *numLocked, 0, left,
         *restartSize, flags, lockedResNorms, hVals, &reset /*dummy*/, machEps,
         rwork, rworkSize, iwork, primme);
   if (ret != 0) return ret;

   /* -------------------------------------------------------------- */
   /* Copy the values for the converged pairs into evals       */
   /* -------------------------------------------------------------- */

   for (i=left, j=0; i < *restartSize; i++) {
      if (flags[i] != UNCONVERGED && *numLocked+j < primme->numEvals) {
         evals[*numLocked+j++] = hVals[i];
      }
   }

   /* -------------------------------------------------------------- */
   /* Merge back the pairs that failed to be locked into the         */
   /* restarted basis.                                               */
   /* -------------------------------------------------------------- */

   /* Generate the permutation that separates the pairs that failed  */
   /* to be locked to the ones that passed.                          */

   ifailed = iwork;
   for (i=left, failed=0; i < *restartSize; i++)
      if (flags[i] == UNCONVERGED) ifailed[failed++] = i-left;
   for (i=left, j=0; i < *restartSize; i++)
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

      blockNorms0 = (double*)rwork;
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

      Num_compute_residual_columns_dprimme(nLocal, &hVals[left],
            &V[left*ldV], failed, ifailed, ldV, &W[left*ldV], ldV, *X,
            sizeBlockNorms, ldV, 0, *R, ldV,
            &V[(left+failed)*ldV], maxBlockSize, hVecsPerm, ldV,
            &W[(left+failed)*ldV], ldV, rwork, rworkSize);
   }
   else {
      /* The failed pairs are not rearranged with the rest of           */
      /* non-converged pairs and they may not be in the block in the    */
      /* next iteration. This was the behaviour of locking in previous  */
      /* versions than 2.0.                                             */

      for (i=0; i<basisSize; i++) hVecsPerm[i] = i;

      /* Pack together the failed vectors and place them after          */
      /* candidates.                                                    */

      Num_compact_vecs_dprimme(&V[left*ldV], nLocal, failed, ldV,
            ifailed, &V[left*ldV], ldV, 0);
      Num_compact_vecs_dprimme(&W[left*ldV], nLocal, failed, ldV,
            ifailed, &W[left*ldV], ldV, 0);

      /* Copy X and R after the failed vectors */

      Num_copy_matrix_dprimme(*X, nLocal, sizeBlockNorms, ldV,
            &V[(left+failed)*ldV], ldV);
      Num_copy_matrix_dprimme(*R, nLocal, sizeBlockNorms, ldV,
            &W[(left+failed)*ldV], ldV);
   }

   /* Modify hVals, hVecs and restartPerm to add the pairs failed to */
   /* be locked just after the candidates                            */

   Num_compact_vecs_dprimme(&hVecs[left*ldhVecs], basisSize, failed,
         ldhVecs, ifailed, &hVecs[left*ldhVecs], ldhVecs, 0);
   Num_compact_vecs_dprimme(&hVals[left], 1, failed, 1, ifailed, &hVals[left],
         1, 0);
   permute_vecs_iprimme(&restartPerm[left], numPacked, ifailed,
         ifailed+numPacked);

   *X = &V[(left+failed)*ldV];
   *R = &W[(left+failed)*ldV];

   /* Pack those Ritz vectors that are actually converged from */
   /* those that have been copied into evecs. Also insertSort  */
   /* the converged Ritz value within the evals array.         */ 

   for (i=left; i < *restartSize; i++) {
       if (flags[i] != UNCONVERGED && *numLocked < primme->numEvals) {
         double resNorm = lockedResNorms[i-left];
         double eval = evals[*numLocked];
         if (!overbooking)
            Num_copy_matrix_dprimme(
                  &evecs[(numLocked0+i-left+primme->numOrthoConst)*primme->nLocal],
                  nLocal, 1, primme->nLocal,
                  &evecs[(*numLocked+primme->numOrthoConst)*primme->nLocal],
                  primme->nLocal);
         else
            Num_copy_matrix_dprimme(&V[i*nLocal], nLocal, 1, primme->nLocal,
                  &evecs[(*numLocked+primme->numOrthoConst)*primme->nLocal],
                  primme->nLocal);
         insertionSort(eval, evals, resNorm, resNorms, evecsperm,
            *numLocked, primme);
         (*numLocked)++;

         /* Update maxConvTol */
         primme->stats.maxConvTol = max(primme->stats.maxConvTol, resNorm);

         if (primme->printLevel >= 2 && primme->procID == 0) { 
            fprintf(primme->outputFile, 
                  "Lock epair[ %d ]= %e norm %.4e Mvecs %d Time %.4e Flag %d\n",
                  *numLocked, eval, resNorm, 
                  primme->stats.numMatvecs, primme_wTimer(0), flags[i]);
            fflush(primme->outputFile);
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
