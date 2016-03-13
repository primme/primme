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
#include "locking_z.h"
#include "locking_private_z.h"
#include "ortho_z.h"
#include "convergence_z.h"
#include "update_projection_z.h"
#include "update_W_z.h"
#include "solve_H_z.h"
#include "restart_z.h"
#include "factorize_z.h"
#include "numerical_z.h"
#include <assert.h>

/*******************************************************************************
 * Subroutine: restart_locking - This routine replaces V with V*c, some subset
 *             of the Ritz vectors, corresponding to the restartSize chosen
 *             eigenvalues of V'*A*V. It may include components from the 
 *             Ritz vectors from the (maxBasisSize-1) step (i.e., recurrence
 *             restarting).
 *
 *    This subroutine locks converged Ritz pairs.  The
 *    converged Ritz vectors are copied to evecs, and the converged Ritz values
 *    are copied to evals.  The evals array is maintained in sorted order; 
 *    however the evecs array is not.  Instead, a permutation array is 
 *    maintained so that the locked Ritz vectors may be sorted upon return to
 *    the user.  
 *
 *    Vectors that have remain converged are locked to the evecs array and
 *    are replaced by initial guesses if there are any remaining.  The initial
 *    guesses are orthogonalized.  If any vectors were locked, the eigenproblem
 *    for the projected matrix H is then solved for the modified basis.
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
 * hU               The eigenvectors of QV/R
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
 * numLocked        The number of locked eigenpairs
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
 *
 * Return value
 * ------------
 * int   > 0 the restart size   
 *        -2 restart_H failed
 *        -4 factorization of M failed
 *        -5 flags do not correspond to converged pairs in pseudolocking
 *       
 ******************************************************************************/
 
int restart_locking_zprimme(int *restartSize, Complex_Z *V, Complex_Z *W,
   int nLocal, Complex_Z *hR, int ldhR, Complex_Z *hU, int ldhU,
   int basisSize, int ldV, Complex_Z **X, Complex_Z **R, Complex_Z *hVecs,
   int ldhVecs, int *restartPerm, double *hVals, int *flags, int *iev, int *ievSize,
   double *blockNorms, Complex_Z *evecs, double *evals, int *numConverged, int *numLocked,
   double *resNorms, int *evecsperm, int numGuesses, Complex_Z *previousHVecs,
   int *numPrevRetained, int ldpreviousHVecs, int *indexOfPreviousVecs, int *hVecsPerm,
   double machEps, Complex_Z *rwork, int rworkSize, int *iwork, primme_params *primme) {

   int i, j, k;             /* Loop variables                       */
   int numPacked;           /* The number of coefficient vectors moved to the */
                            /* end of the hVecs array.                        */
   int maxBlockSize;          /* max block size for the next iteration        */
   int sizeBlockNorms;        /* size of blockNorms                           */
   int failed;                /* Number of locked vectors that didn't pass the conv. test */
   int *ifailed;              /* Indices of the locked vectors that failed    */
   int left, ret, numLocked0 = *numLocked; /* aux variables                */
   double *blockNorms0;

   /* Return memory requirement */
   if (V == NULL) {
      Complex_Z t;
      double d;
      return max(max(max(max(
            /* for permute_vecs */
            basisSize,
            Num_compute_residual_i_zprimme(nLocal, NULL, NULL, basisSize, NULL, 0,
               NULL, 0, NULL, primme->maxBlockSize, 0, NULL, 0, NULL,
               primme->maxBlockSize, NULL, 0, NULL, 0, NULL, 0)),
            ortho_coefficient_vectors_zprimme(NULL, basisSize, 0, 0, *restartSize, NULL, NULL,
               0, NULL, 0, *numPrevRetained, 0.0, NULL, NULL, 0, primme)),
            Num_update_VWXR_zprimme(NULL, NULL, 0, basisSize, 0, NULL,
               *restartSize, 0, NULL,
               &t, 0, *restartSize+*numLocked, 0,
               &t, 0, *ievSize, 0,
               &t, *restartSize, *restartSize+*numLocked, 0,
               &t, 0, *restartSize, 0,
               &t, 0, *ievSize, 0, &d,
               &d, *restartSize, *numLocked,
               NULL, 0, primme)),
            check_convergence_zprimme(NULL,
               nLocal, 0, NULL, 0, NULL, *numLocked, 0, *restartSize,
               *restartSize+*numLocked, NULL, NULL, NULL, 0.0, NULL, 0,
               NULL, primme));
   }

   /* ----------------------------------------------------------------------- */
   /* Limit restartSize so that it plus 'to be locked' plus previous Ritz     */
   /* vectors do not exceed basisSize.                                        */
   /* ----------------------------------------------------------------------- */

   *restartSize = min(*restartSize, basisSize-(*numConverged-*numLocked));

   /* ----------------------------------------------------------------------- */
   /* Insert as many initial guesses as eigenpairs have converged.            */
   /* ----------------------------------------------------------------------- */

   *restartSize -= min(min(numGuesses, *numConverged-*numLocked), *restartSize);

   /* ----------------------------------------------------------------------- */
   /* Insert vectors from previous iteration.                                 */
   /* ----------------------------------------------------------------------- */

   *numPrevRetained = min(basisSize, *restartSize+*numConverged-*numLocked + *numPrevRetained)
      - (*restartSize+*numConverged-*numLocked);
   *indexOfPreviousVecs = *restartSize;

   /* ----------------------------------------------------------------------- */
   /* Swap coefficient vectors and hVals corresponding to                     */
   /* converged Ritz vectors to the end of the hVecs(:, restartSize) subarray.*/
   /* This allows the converged Ritz vectors to be stored contiguously in     */
   /* memory after restart.  This significantly reduces the amount of data    */
   /* movement the locking routine would have to perform otherwise.           */
   /* ----------------------------------------------------------------------- */

   for (i=j=0; i<basisSize; i++) {
      if (flags[i] == UNCONVERGED) {
         restartPerm[j < *restartSize+*numPrevRetained ? j : *numConverged-*numLocked+j] = i;
         j++;
      }
   }
   assert(j >= *restartSize+*numPrevRetained);
   for (i=numPacked=0, j=*restartSize+*numPrevRetained; i<basisSize; i++) {
      if (flags[i] != UNCONVERGED) {
         restartPerm[j++] = i;
         numPacked++;
      }
   }

   assert(numPacked == *numConverged-*numLocked);
 
   permute_vecs_dprimme(hVals, 1, basisSize, 1, restartPerm, (double*)rwork, iwork);
   permute_vecs_zprimme(hVecs, basisSize, basisSize, ldhVecs, restartPerm, rwork, iwork);

   *restartSize += numPacked + *numPrevRetained;

   /* ----------------------------------------------------------------------- */
   /* Restarting with a small number of coefficient vectors from the previous */
   /* iteration can be retained to accelerate convergence.  The previous      */
   /* coefficient vectors must be combined with the current coefficient       */
   /* vectors by first orthogonalizing the previous ones versus the current   */
   /* restartSize ones.  The orthogonalized previous vectors are then         */
   /* inserted into the hVecs array at hVecs(:,indexOfPreviousVecs).          */
   /* ----------------------------------------------------------------------- */

   Num_copy_matrix_zprimme(previousHVecs, basisSize, *numPrevRetained,
         ldpreviousHVecs, &hVecs[ldhVecs*(*indexOfPreviousVecs)], ldhVecs);

   ret = ortho_coefficient_vectors_zprimme(hVecs, basisSize, ldhVecs, *indexOfPreviousVecs,
         *restartSize, restartPerm, hU, ldhU, hR, ldhR, *numPrevRetained, machEps,
         iwork, rwork, rworkSize, primme);
   if (ret != 0) return ret;

   /* -------------------------------------------------------------- */
   /* Restart V and W by replacing it with the current Ritz vectors. */
   /* -------------------------------------------------------------- */

   maxBlockSize = min(min(primme->maxBlockSize, primme->numEvals-*numConverged),
         primme->maxBasisSize-*restartSize);
   sizeBlockNorms = X ? min(maxBlockSize, *indexOfPreviousVecs) : 0;
   left = *restartSize-*numConverged+*numLocked;
   if (X) {
      *X = &V[*restartSize*ldV];
      *R = &W[*restartSize*ldV];
   }
   ret = Num_update_VWXR_zprimme(V, W, nLocal, basisSize, ldV, hVecs,
         *restartSize, ldhVecs, hVals,
         V, 0, *restartSize, ldV,
         X?*X:NULL, 0, sizeBlockNorms, ldV,
         &evecs[primme->nLocal*(*numLocked+primme->numOrthoConst)], left, *restartSize, primme->nLocal,
         W, 0, *restartSize, ldV,
         X?*R:NULL, 0, sizeBlockNorms, ldV, blockNorms,
         &resNorms[*numLocked], left, *restartSize,
         rwork, rworkSize, primme);
   if (ret != 0) return ret;
 
   /* ----------------------------------------------------------------------------- */
   /* Recompute flags for the new locked Ritz vectors                               */
   /* NOTE: the eigenvals flagged as practically converged will keep it if their    */
   /*       residual norm is still less than sqrt(numLocked)*tol                    */
   /* ----------------------------------------------------------------------------- */

   permute_vecs_iprimme(flags, basisSize, restartPerm, iwork);
   ret = check_convergence_zprimme(&V[ldV*left],
         nLocal, ldV, NULL, 0, NULL, *numLocked, 0, left,
         *restartSize, flags, &resNorms[*numLocked], hVals, machEps, rwork, rworkSize,
         iwork, primme);
   if (ret != 0) return ret;

   for (i=left; i < *restartSize && i-left+*numLocked < primme->numEvals; i++)
      evals[*numLocked+i-left] = hVals[i];

   /* ----------------------------------------------------------------------------- */
   /* Pack hVals, hVecs, blockNorms and restartPerm for the unconverged pairs.      */
   /* When there are more unconverged than vacancies in the block overwrite the     */
   /* residual vectors computed previously in Num_update_VWXR.                      */ 
   /* ----------------------------------------------------------------------------- */

   ifailed = iwork;
   for (i=left, failed=0; i < *restartSize; i++)
      if (flags[i] == UNCONVERGED) ifailed[failed++] = i-left;
   for (i=left, j=0; i < *restartSize; i++)
      if (flags[i] != UNCONVERGED) ifailed[failed+j++] = i-left;

   if (1 /* Put zero to disable the new feature */) {
      /* New feature: the candidate pairs to be locked that failed the              */
      /* convergence test are rearrange with the rest of non-converged pairs        */
      /* and are included in the block.                                             */

      /* Generate hVecsPerm merging back the locked pairs that failed to pass the   */
      /* convergence test with the restarted pairs.                                 */
      blockNorms0 = (double*)rwork;
      for (i=0; i<sizeBlockNorms; i++) blockNorms0[i] = blockNorms[i];
      for (i=j=k=0; i<*indexOfPreviousVecs || j<failed; k++) {
         if (i < *indexOfPreviousVecs && (j >= failed || restartPerm[i] < restartPerm[left+ifailed[j]])) {
            if (k < maxBlockSize && i < sizeBlockNorms) blockNorms[k] = blockNorms0[i];
            hVecsPerm[k] = i++;
         }
         else {
            if (k < maxBlockSize) blockNorms[k] = resNorms[numLocked0+ifailed[j]];
            hVecsPerm[k] = *indexOfPreviousVecs + *numPrevRetained + j++;
         }
      }

      /* Generate the rest of the permutation of hVecsPerm                          */
      for (i=0; i<*numPrevRetained; i++) hVecsPerm[k++] = *indexOfPreviousVecs+i;
      assert(k == *indexOfPreviousVecs + *numPrevRetained + failed);
      for (; k < basisSize; k++) hVecsPerm[k] = -1;

      /* Pack X and R for the unconverged pairs.                                     */
      Num_compute_residual_i_zprimme(nLocal, &hVals[left], &V[left*ldV], failed, ifailed, 
            ldV, &W[left*ldV], ldV,
            X?*X:NULL, sizeBlockNorms, ldV, X?*R:NULL, ldV,
            &V[(left+failed)*ldV], maxBlockSize, hVecsPerm, ldV, &W[(left+failed)*ldV], ldV,
            rwork, rworkSize);
   }
   else {
      /* The failed candidates are not rearrange with the rest of non-converged pairs  */
      /* and they may not be in block in the next iteration. This was the behaviour of */
      /* locking in previous versions than 2.0.                                        */

      for (i=0; i<basisSize; i++) hVecsPerm[i] = i;

      /* Pack V and W for the unconverged pairs.                                       */
      Num_compact_vecs_zprimme(&V[left*ldV], nLocal, failed, ldV, ifailed, &V[left*ldV],
            ldV, 0);
      Num_compact_vecs_zprimme(&W[left*ldV], nLocal, failed, ldV, ifailed, &W[left*ldV],
            ldV, 0);

      Num_copy_matrix_zprimme(*X, nLocal, sizeBlockNorms, ldV, &V[(left+failed)*ldV], ldV);
      Num_copy_matrix_zprimme(*R, nLocal, sizeBlockNorms, ldV, &W[(left+failed)*ldV], ldV);
   }

   /* Pack hVals, hVecs and restartPerm for the failed pairs  */
   Num_compact_vecs_zprimme(&hVecs[left*ldhVecs], basisSize, failed, ldhVecs, ifailed,
         &hVecs[left*ldhVecs], ldhVecs, 0);
   Num_compact_vecs_dprimme(&hVals[left], 1, failed, 1, ifailed, &hVals[left], 1, 0);
   permute_vecs_iprimme(&restartPerm[left], numPacked, ifailed, ifailed+numPacked);

   if (X) {
      *X = &V[(left+failed)*ldV];
      *R = &W[(left+failed)*ldV];
   }

   /* Pack the converged Ritz vector into the evecs array and */
   /* insert the converged Ritz value in sorted order within  */
   /* the evals array.                                        */

   for (i=left; i < *restartSize; i++) {
       if (flags[i] != UNCONVERGED) {
         double resNorm = resNorms[i-left+numLocked0];
         double eval = evals[i-left+numLocked0];
         Num_copy_matrix_zprimme(&evecs[(numLocked0+i-left+primme->numOrthoConst)*primme->nLocal], nLocal, 1,
               primme->nLocal, &evecs[(*numLocked+primme->numOrthoConst)*primme->nLocal], primme->nLocal);
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
