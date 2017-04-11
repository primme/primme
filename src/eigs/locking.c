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
 * File: locking.c
 *
 * Purpose - This file contains the routines to restart with locking.
 *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "wtime.h"
#include "const.h"
#include "numerical.h"
#include "locking.h"
#include "convergence.h"
#include "auxiliary_eigs.h"
#include "restart.h"

static void insertionSort(REAL newVal, REAL *evals, REAL newNorm,
   REAL *resNorms, int newFlag, int *flags, int *perm, int numLocked,
   primme_params *primme);

static int compute_residual_columns(PRIMME_INT m, REAL *evals, SCALAR *x,
      int n, int *p, PRIMME_INT ldx, SCALAR *Ax, PRIMME_INT ldAx,
      SCALAR *xo, int no, PRIMME_INT ldxo, int io0, SCALAR *ro, PRIMME_INT ldro,
      SCALAR *xd, int nd, int *pd, PRIMME_INT ldxd, SCALAR *rd, PRIMME_INT ldrd,
      SCALAR *rwork, PRIMME_INT lrwork);

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
 
TEMPLATE_PLEASE
int restart_locking_Sprimme(int *restartSize, SCALAR *V, SCALAR *W, 
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


