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
 * File: restart.c
 *
 * Purpose - Compute the Ritz vectors corresponding to the restartSize
 *           smallest eigenvalues.
 *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "primme.h"
#include "const.h"
#include "restart_@(pre).h"
#include "restart_private_@(pre).h"
#include "ortho_@(pre).h"
#include "solve_H_@(pre).h"
#include "factorize_@(pre).h"
#include "update_projection_@(pre).h"
#include "update_W_@(pre).h"
#include "convergence_@(pre).h"
#include "numerical_@(pre).h"


/*******************************************************************************
 * Subroutine: restart - This routine replaces V with V*c, some subset
 *             of the Ritz vectors, corresponding to the restartSize chosen
 *             eigenvalues of V'*A*V. It may include components from the 
 *             Ritz vectors from the (maxBasisSize-1) step (i.e., recurrence
 *             restarting).
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 *
 * nLocal           Number of rows of V assigned to the node
 *
 * ldV              The leading dimension of V and W
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
 * hVecsperm        The permutation applied to the columns of hVecs before restarting
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
 * Vperm            The permutation that orders the output hVals and hVecs as primme.target
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
 
int restart_@(pre)primme(int *restartSize, @(type) *V, @(type) *W, int nLocal,
   int basisSize, int ldV, @(type) **X, @(type) **R, @(type) *hVecs, int ldhVecs,
   int *hVecsperm, double *hVals, int *flags, int *iev, int *ievSize,
   double *blockNorms, @(type) *evecs, double *evals, double *resNorms,
   @(type) *evecsHat, int ldevecsHat, @(type) *M, int ldM, int *numConverged,
   int *numConvergedStored, @(type) *previousHVecs, int *numPrevRetained,
   int ldpreviousHVecs, int *indexOfPreviousVecs, int *Vperm, double machEps,
   @(type) *rwork, int rworkSize, int *iwork, primme_params *primme) {
 
   int i, j, ret;
   @(type) tpone = @(tpone), tzero = @(tzero);             /*constants*/

   /* Return memory requirement */

   if (V == NULL) {
      @(type) t;
      double d;
      return max(max(max(
                  nLocal,      /* permute_vecs for hVecs */
                  Num_update_VWXR_@(pre)(NULL, NULL, nLocal, basisSize, 0, &t,
                     *restartSize, 0, NULL,
                     &t, 0, *restartSize, 0,
                     &t, *numConverged, *numConverged+*ievSize, 0,
                     NULL, 0, 0, 0,
                     &t, 0, *restartSize, 0,
                     &t, *numConverged, *numConverged+*ievSize, 0, &d,
                     NULL, 0, 0,
                     NULL, 0, primme)),
                  /* if evecsHat, permutation matrix & compute_submatrix workspace */
                  evecsHat ? (primme->numOrthoConst+*numConverged)*
                     (primme->numOrthoConst+*numConverged)*2 : 0),
                  ortho_@(pre)primme(NULL, 0, NULL, 0, 0,
                     *numPrevRetained-1, NULL, 0, *restartSize, 0,
                     NULL, 0, NULL, 0, NULL));
   }

   /* -------------------------------------------------------------------------- */ 
   /* Check if any of the previous flagged converged eigenvalues seems           */
   /* to have become unconverged by checking hVals[i]-evals[i] < tol.            */
   /* If not, flag it UNCONVERGED and let it be targeted again. This avoids      */  
   /* early converged but unwanted evs preventing wanted from being targeted.    */
   /* -------------------------------------------------------------------------- */

   if (basisSize != primme->n) {
      for (i=0; i<primme->numEvals; i++) {
         if (flags[i] != UNCONVERGED && fabs(hVals[i]-evals[i]) > resNorms[i]) {
            flags[i] = UNCONVERGED;
         }
      }
   }

   /* ----------------------------------------------------------------- */
   /* Generate the permutation that put converged Ritz vectors first.   */
   /* Also update the number of converged values.                       */
   /* ----------------------------------------------------------------- */

   for (i=j=0, *numConverged=0; i<basisSize; i++)
      if (flags[i] != UNCONVERGED) hVecsperm[j++] = i, (*numConverged)++;
   *numConverged = min(*numConverged, primme->numEvals);
   for (i=0; i<basisSize; i++)
      if (flags[i] == UNCONVERGED) hVecsperm[j++] = i;

   permute_vecs_@(pre)(hVecs, basisSize, basisSize, ldhVecs, hVecsperm, rwork,
         iwork);

   /* Permute hVals */
   permute_vecs_d(hVals, 1, basisSize, 1, hVecsperm, (double*)rwork, iwork);

   /* ----------------------------------------------------------------------- */
   /* Restarting with a small number of coefficient vectors from the previous */
   /* iteration can be retained to accelerate convergence.  The previous      */
   /* coefficient vectors must be combined with the current coefficient       */
   /* vectors by first orthogonalizing the previous ones versus the current   */
   /* restartSize ones.  The orthogonalized previous vectors are then         */
   /* inserted into the hVecs array at hVecs(:,indexOfPreviousVecs).          */
   /* ----------------------------------------------------------------------- */

   *numPrevRetained = min(basisSize, *restartSize + *numPrevRetained)
      - *restartSize;
   *indexOfPreviousVecs = *restartSize;

   Num_copy_matrix_@(pre)primme(previousHVecs, basisSize, *numPrevRetained,
         ldpreviousHVecs, &hVecs[ldhVecs*(*indexOfPreviousVecs)], ldhVecs);

   ret = ortho_@(pre)primme(hVecs, ldhVecs, NULL, 0, *indexOfPreviousVecs,
         *indexOfPreviousVecs+*numPrevRetained-1, NULL, 0, 0, basisSize,
         primme->iseed, machEps, rwork, rworkSize, NULL);
   if (ret != 0) return ret;

   *restartSize += *numPrevRetained;

   /* -------------------------------------------------------------- */
   /* Restart V and W by replacing it with the current Ritz vectors. */
   /* Compute X, R, blockNorms for the next values in the block.     */
   /* -------------------------------------------------------------- */

   if (X) {
      *ievSize = min(min(min(primme->maxBlockSize, primme->numEvals-*numConverged),
               *indexOfPreviousVecs-*numConverged), primme->maxBasisSize-*restartSize);
      *X = &V[*restartSize*ldV];
      *R = &W[*restartSize*ldV];
   }
   else {
      *ievSize = 0;
   }

   ret = Num_update_VWXR_@(pre)(V, W, nLocal, basisSize, ldV, hVecs,
         *restartSize, ldhVecs, hVals,
         V, 0, *restartSize, ldV,
         X?*X:NULL, *numConverged, *numConverged+*ievSize, ldV,
         NULL, 0, 0, 0,
         W, 0, *restartSize, ldV,
         X?*R:NULL, *numConverged, *numConverged+*ievSize, ldV, blockNorms,
         NULL, 0, 0,
         rwork, rworkSize, primme);
   if (ret != 0) return ret;

   if (X)
      for (i=0; i<*ievSize; i++)
         iev[i] = hVecsperm[i+*numConverged];

   /* Undo the reordering of hVecsperm */

   for (i=0; i<basisSize; i++)
      Vperm[i] = hVecsperm[i];
   permute_vecs_i(Vperm, basisSize, hVecsperm, iwork);
   

   /* --------------------------------------------------------------------- */
   /* If the user requires (I-QQ') projectors in JDQMR without locking,     */
   /* the converged eigenvectors are copied temporarily to evecs. There     */
   /* they stay locked  for use in (I-QQ') and (I-K^{-1}Q () Q') projectors.*/
   /* NOTE THIS IS NOT LOCKING! The Ritz vectors remain in the basis, and   */
   /* they will overwrite evecs at the end.                                 */
   /* We recommend against this type of usage. It's better to use locking.  */
   /* --------------------------------------------------------------------- */

   /* Andreas NOTE: is done inefficiently for the moment. We should only */
   /* add the recently converged. But we need to differentiate them      */
   /* from flags...                                                      */

   if (evecsHat) {
      int newNumConvergedStored=0, oldSizeM, newSizeM;

      /* Pack evecs and evecsHat for the converged pairs hVecsperm[0:numConverged] */

      for (i=0; i < *numConverged && hVecsperm[i] < *numConvergedStored; i++) {
         Num_copy_matrix_@(pre)primme(&evecs[(hVecsperm[i]+primme->numOrthoConst)*nLocal],
               nLocal, 1, nLocal,
               &evecs[(newNumConvergedStored+primme->numOrthoConst)*nLocal],
               nLocal);
         Num_copy_matrix_@(pre)primme(&evecsHat[(hVecsperm[i]+primme->numOrthoConst)*ldevecsHat],
               nLocal, 1, ldevecsHat,
               &evecsHat[(newNumConvergedStored+primme->numOrthoConst)*ldevecsHat],
               ldevecsHat);
         newNumConvergedStored++;
      }

      /* Apply hVecsperm to rows and columns of M */

      oldSizeM = *numConvergedStored + primme->numOrthoConst;
      newSizeM = newNumConvergedStored + primme->numOrthoConst;
      for (i=0; i<oldSizeM*newSizeM; i++)
         rwork[i] = tzero;
      for (i=0; i < primme->numOrthoConst; i++)
         rwork[oldSizeM*i + i] = tpone;
      for (; i < newSizeM; i++)
         rwork[oldSizeM*i + hVecsperm[i]+primme->numOrthoConst] = tpone;
      compute_submatrix(rwork, newSizeM, oldSizeM, M, oldSizeM, ldM,
         M, ldM, rwork+oldSizeM*newSizeM, rworkSize-oldSizeM*newSizeM);

      *numConvergedStored = newNumConvergedStored;
   }

   return 0;
}


/*******************************************************************************
 * Subroutine: restart_projected_problem - This routine replaces V with V*c, some subset
 *             of the Ritz vectors, corresponding to the restartSize chosen
 *             eigenvalues of V'*A*V. It may include components from the 
 *             Ritz vectors from the (maxBasisSize-1) step (i.e., recurrence
 *             restarting).
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal           Number of rows of V assigned to the node
 *
 * iev              Array of size blockSize indicating which Ritz vectors are
 *                  being targeted by the block.
 *
 * basisSize        Size of the basis V
 *
 * numConverged     The number of converged eigenpairs without locking
 *
 * numLocked        The number of Ritz vectors that have been locked 
 *
 * numGuesses       Number of remaining initial guesses
 *
 * previousHVecs    Coefficient vectors retained from the previous iteration
 *
 * numPrevRetained  The number of coefficient vectors in previousHVecs
 *  
 * rwork            Real work array
 *
 * rworkSize        Must be of size 
 *                  (primme->restartingParams.maxPrevRetain)^2 + 
 *                  MAX(
 *                  primme->maxBasisSize*primme->restartingParams.maxPrevRetain,
#ifdefarithm L_DEFCPLX
 *                  5*primme->maxPrevRetain, primme->maxBasisSize^2)
#endifarithm
#ifdefarithm L_DEFREAL
 *                  3*primme->maxPrevRetain, primme->maxBasisSize^2)
#endifarithm
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
 * W                A*V
 *
 * H                The projection V'*A*V
 *
 * hVecs            The eigenvectors of H
 *
 * hVals            The eigenvalues of H
 *
 * flags            Array indicating the convergence of the Ritz vectors
 *
 * evecs            The converged Ritz vectors. Without locking, all converged
 *                  eigenvectors are copied from V to evecs if skew projections
 *                  are required.
 *
 * evecsHat         K^{-1}evecs
 *
 * M                evecs'*evecsHat
 *
 * UDU              The factorization of M
 *
 * ipivot           The pivot array of the UDU factorization
 *
 * numConvergedStored The # of converged vectors copied to evecs w/o locking
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
 
int after_restart_@(pre)primme(@(type) *V, int ldV, @(type) *W, int ldW,
   @(type) *H, int ldH, @(type) *Q, int nLocal, int ldQ, @(type) *R, int ldR,
   @(type) *hU, int ldhU, int newldhU, @(type) *hVecs, int ldhVecs, int newldhVecs,
   double *hVals, double *hSVals, int *hVecsperm, int *Vperm,
   int restartSize, int basisSize, int numPrevRetained,
   int indexOfPreviousVecs, @(type) *evecs, int *evecsSize,
   int ldevecs, @(type) *evecsHat, int ldevecsHat, @(type) *M, int ldM, @(type) *UDU,
   int ldUDU, int *ipivot, int numConvergedBeforeRestart, int numConverged,
   int rworkSize, @(type) *rwork, int *iwork, double machEps, primme_params *primme) {

   int ret;

   /* -------------------------------------------------------- */
   /* Restart projected problem matrices H and R               */
   /* -------------------------------------------------------- */

   switch (primme->projectionParams.projection) {
   case primme_proj_RR:
      ret = restart_RR(H, ldH, hVecs, ldhVecs, newldhVecs, hVals, restartSize,
            basisSize, numConverged, numPrevRetained, indexOfPreviousVecs,
            Vperm, rworkSize, rwork, iwork, primme);
      break;

   case primme_proj_ref:
      ret = restart_ref(V, ldV, W, ldW, H, ldH, Q, nLocal, ldQ, R, ldR, hU, ldhU, newldhU, hVecs,
            ldhVecs, newldhVecs, hVals, hSVals, hVecsperm, Vperm, restartSize, basisSize,
            numPrevRetained, indexOfPreviousVecs, numConvergedBeforeRestart, numConverged,
            rworkSize, rwork, iwork, machEps, primme);
      break;

   default:
      assert(0);
   }

   if (H && ret != 0) {
      primme_PushErrorMessage(Primme_restart, Primme_restart_h, ret, __FILE__, 
         __LINE__, primme);
      return RESTART_H_FAILURE;
   }

   if (evecsHat) {
      int numRecentlyConverged = numConverged - *evecsSize;

      /* Return memory requirement */
      if (H == NULL) {
         return max(max(
               /* Workspace for restart_RR or restart_ref */
               ret,
               update_projection_@(pre)primme(NULL, 0, NULL, 0, NULL, 0, nLocal,
                  *evecsSize, basisSize, NULL, 0, primme)),
               UDUDecompose_@(pre)primme(NULL, 0, NULL, 0, NULL, *evecsSize, NULL, 
                  0, primme));
      }

      /* Compute K^{-1}x for all newly locked eigenvectors */

      /* TODO: primme.shiftsForPreconditioner is undefined at that point;
         maybe it makes sense to always set NULL shiftsForPreconditioner
         when SkewQ is enabled to force the same preconditioner. */
      (*primme->applyPreconditioner)(&evecs[primme->nLocal*(*evecsSize+primme->numOrthoConst)],
            &evecsHat[primme->nLocal*(*evecsSize+primme->numOrthoConst)], &numRecentlyConverged,
            primme);
      primme->stats.numPreconds += numRecentlyConverged;

      /* Update the projection evecs'*evecsHat now that evecs and evecsHat   */
      /* have been expanded by numRecentlyConverged columns.  Required       */
      /* workspace is numLocked*numEvals.  The most ever needed would be     */
      /* maxBasisSize*numEvals.                                              */

      update_projection_@(pre)primme(evecs, primme->nLocal, evecsHat, primme->nLocal,
            M, ldM, nLocal, *evecsSize+primme->numOrthoConst, numRecentlyConverged, rwork,
            rworkSize, primme);
      *evecsSize = numConverged;

      ret = UDUDecompose_@(pre)primme(M, ldM, UDU, ldUDU, ipivot, *evecsSize+primme->numOrthoConst,
            rwork, rworkSize, primme);

      if (ret != 0) {
         primme_PushErrorMessage(Primme_lock_vectors, Primme_ududecompose, ret,
            __FILE__, __LINE__, primme);
         return UDUDECOMPOSE_FAILURE;
      }

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
 * previousHVecs   Coefficient vectors from the previous iteration.  They are
 *                 orthonormal to the current coefficient vectors.
 *
 * numPrevRetained The number of vectors retained from the previous iteration
 *
 * indexOfPreviousVecs  The index within hVecs where the previous vectors were
 *                      inserted.  Its also where the overlap matrix
 *                      previousHVecs'*H*previousHvecs will be inserted within
 *                      the restarted H.
 *
 * rwork         Work array.  Necessary only when coefficient vectors from the
 *               previous iteration are retained.
 *
 * rworkSize     Can be zero if no previous vectors are retained.  Otherwise,
 *               it must be at least 
 *               numPrevRetained*numPrevRetained + 
#ifdefarithm L_DEFCPLX
 *               max(basisSize*numPrevRetained, 5*numPrevRetained) 
#endifarithm
#ifdefarithm L_DEFREAL
 *               max(basisSize*numPrevRetained, 3*numPrevRetained) 
#endifarithm
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
 * hVecs  If the new H is diagonal, then it will contain the standard basis
 *        vectors.  If previous coefficient vectors are retained, then 
 *        restartSize - numPrevRetained of the vectors will be standard basis
 *        vectors.  The remaining numPrevRetained vectors will contain
 *        numPrevRetained non-zero elements corresponding to the 
 *        numPrevRetined x numPrevRetained submatrix.
 * 
 * hVals  The eigenvalues of the restarted H
 * 
 * 
 * Return value
 * ------------
 * Error code: 0 upon success
 *            -1 eigenvalues of submatrix could not be computed
 *
 ******************************************************************************/

static int restart_RR(@(type) *H, int ldH, @(type) *hVecs, int ldhVecs,
   int newldhVecs, double *hVals, int restartSize, int basisSize, int numLocked,
   int numPrevRetained, int indexOfPreviousVecs, int *Vperm,
   int rworkSize, @(type) *rwork, int *iwork, primme_params *primme) {

   int i, j;          /* Loop variables                                       */
   int ret;           /* Return value                                         */
   int orderedIndexOfPreviousVecs;  /* index of prev. vecs after applying Vperm */
   @(type) tpone = @(tpone), tzero = @(tzero);             /*constants*/

   /* Return memory requirement */

   if (H == NULL) {
      return max(
            compute_submatrix(NULL, numPrevRetained, 0, NULL, basisSize, 0, NULL,
               0, NULL, 0),
            solve_H_RR_@(pre)primme(NULL, 0, NULL, 0, NULL, numPrevRetained,
               numLocked, 0, NULL, NULL, primme));
   }

   /* ---------------------------------------------------------------------- */
   /* If coefficient vectors from the previous iteration were retained, then */
   /* set up work space for computing the numPrevRetained x numPrevRetained  */
   /* submatrix, then compute the submatrix.                                 */
   /* ---------------------------------------------------------------------- */

   compute_submatrix(&hVecs[ldhVecs*indexOfPreviousVecs], numPrevRetained, ldhVecs, H,
         basisSize, ldH, &H[ldH*indexOfPreviousVecs+indexOfPreviousVecs],
         ldH, rwork, rworkSize);

   /* ----------------------------------------------------------------- */
   /* V*hVecs yields a diagonal matrix composed of the Ritz values of A */
   /* with respect to V.  Set H to a diagonal matrix with the Ritz      */
   /* values on the diagonal.                                           */
   /* ----------------------------------------------------------------- */

   for (j=0; j < indexOfPreviousVecs; j++) {
      for (i=0; i <= j; i++) {
         H[ldH*j+i] = tzero;
      }
      *(double*)&H[ldH*j+j] = hVals[j];
   }
   for (j=indexOfPreviousVecs; j<indexOfPreviousVecs+numPrevRetained; j++) {
      for (i=0; i < indexOfPreviousVecs; i++) {
         H[ldH*j+i] = tzero;
      }
   }
   for (j=indexOfPreviousVecs+numPrevRetained; j < restartSize; j++) {
      for (i=0; i <= j; i++) {
         H[ldH*j+i] = tzero;
      }
      *(double*)&H[ldH*j+j] = hVals[j];
   }

   /* --------------------------------------------------------------------- */
   /* Given the above H, we know the eigenvectors of H will be the standard */
   /* basis vectors if no previous coefficient vectors are retained         */
   /* --------------------------------------------------------------------- */

   for (j=0; j < restartSize; j++) {
      for (i=0; i < restartSize; i++) {
          hVecs[newldhVecs*j+i] = tzero;
      }
      hVecs[newldhVecs*j+Vperm[j]] = tpone;
   }

   /* Apply permutation Vperm to hVals */
   permute_vecs_d(hVals, 1, restartSize, 1, Vperm, (double*)rwork, iwork);

   /* ---------------------------------------------------------------------- */
   /* If coefficient vectors from the previous iteration have been retained, */
   /* then insert the computed overlap matrix into the restarted H and solve */ 
   /* the resulting eigenproblem for the resulting H.                        */
   /* ---------------------------------------------------------------------- */

   for (i=0, orderedIndexOfPreviousVecs=restartSize; i<restartSize; i++) {
      if (Vperm[i] == indexOfPreviousVecs) {
         orderedIndexOfPreviousVecs = i;
         break;
      }
   }
   assert(orderedIndexOfPreviousVecs != restartSize || indexOfPreviousVecs >= restartSize);

   ret = solve_H_RR_@(pre)primme(
         &H[ldH*indexOfPreviousVecs+indexOfPreviousVecs], ldH,
         &hVecs[newldhVecs*orderedIndexOfPreviousVecs+orderedIndexOfPreviousVecs],
         newldhVecs, &hVals[orderedIndexOfPreviousVecs], numPrevRetained,
         numLocked, rworkSize, rwork, iwork, primme);

   /* TODO: reorder all hVals */

   if (ret != 0) {
      primme_PushErrorMessage(Primme_restart_h, Primme_insert_submatrix, 
            ret, __FILE__, __LINE__, primme);
      return INSERT_SUBMATRIX_FAILURE;
   }

   return 0;
}

/*******************************************************************************
 * Function restart_ref - This routine is used to recompute H = V'*A*V once V 
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
 * previousHVecs   Coefficient vectors from the previous iteration.  They are
 *                 orthonormal to the current coefficient vectors.
 *
 * numPrevRetained The number of vectors retained from the previous iteration
 *
 * indexOfPreviousVecs  The index within hVecs where the previous vectors were
 *                      inserted.  Its also where the overlap matrix
 *                      previousHVecs'*H*previousHvecs will be inserted within
 *                      the restarted H.
 *
 * rwork         Work array.  Necessary only when coefficient vectors from the
 *               previous iteration are retained.
 *
 * rworkSize     Can be zero if no previous vectors are retained.  Otherwise,
 *               it must be at least 
 *               numPrevRetained*numPrevRetained + 
#ifdefarithm L_DEFCPLX
 *               max(basisSize*numPrevRetained, 5*numPrevRetained) 
#endifarithm
#ifdefarithm L_DEFREAL
 *               max(basisSize*numPrevRetained, 3*numPrevRetained) 
#endifarithm
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
 * hVecs  If the new H is diagonal, then it will contain the standard basis
 *        vectors.  If previous coefficient vectors are retained, then 
 *        restartSize - numPrevRetained of the vectors will be standard basis
 *        vectors.  The remaining numPrevRetained vectors will contain
 *        numPrevRetained non-zero elements corresponding to the 
 *        numPrevRetined x numPrevRetained submatrix.
 * 
 * hVals  The eigenvalues of the restarted H
 * 
 * 
 * Return value
 * ------------
 * Error code: 0 upon success
 *            -1 eigenvalues of submatrix could not be computed
 *
 ******************************************************************************/

static int restart_ref(@(type) *V, int ldV, @(type) *W, int ldW, @(type) *H,
   int ldH, @(type) *Q, int nLocal, int ldQ, @(type) *R, int ldR, @(type) *hU,
   int ldhU, int newldhU, @(type) *hVecs, int ldhVecs, int newldhVecs,
   double *hVals, double *hSVals, int *perm, int *Vperm, int restartSize,
   int basisSize, int numPrevRetained, int indexOfPreviousVecs,
   int numConvergedBeforeRestart, int numConverged, int rworkSize,
   @(type) *rwork, int *iwork, double machEps, primme_params *primme) {

   int i, j;          /* Loop variables                                       */
   int ret;           /* Return value                                         */
   @(type) tpone = @(tpone), tzero = @(tzero);             /*constants*/
   int orderedIndexOfPreviousVecs;  /* index of prev. vecs after applying Vperm */

   /* Return memory requirement */
 
   if (H == NULL) {
      @(type) t;
      int geqrfSize;    /* Workspace used by Num_geqrf_@(pre)primme */
      int orgqrSize;    /* Workspace used by Num_orgqr_@(pre)primme */

      Num_geqrf_@(pre)primme(basisSize, numPrevRetained, NULL, basisSize, NULL, &t, -1, &ret);
      geqrfSize = *(double*)&t;
      Num_orgqr_@(pre)primme(basisSize, numPrevRetained, numPrevRetained, NULL, basisSize, NULL,
         &t, -1, &ret);
      orgqrSize = *(double*)&t;

      return max(max(max(max(max(max(
         compute_submatrix(NULL, basisSize, 0, NULL, basisSize, 0, NULL, 0,
            NULL, 0),
         update_Q_@(pre)primme(NULL, nLocal, 0, NULL, 0, NULL, 0, NULL, 0, 0.0, 0,
            basisSize, NULL, 0, 0.0, primme)),
         /* Workspace for  R(indexOfPrevVecs:) = R * hVecs(indexOfPrevVecs:) */
         basisSize*basisSize),
         /* Workspace for permute_vecs(hU) */
         basisSize),
         basisSize+max(geqrfSize, orgqrSize)),
         Num_update_VWXR_@(pre)(NULL, NULL, nLocal, basisSize, 0, NULL, basisSize, 0, NULL,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0,
            NULL, 0, 0, 0, NULL,
            NULL, 0, 0,
            NULL, 0, primme)),
         solve_H_Ref_@(pre)primme(NULL, 0, NULL, 0, NULL, 0, NULL, NULL, 0, NULL,
            numPrevRetained, 0, NULL, primme));
   }

   /* ------------------------------- */
   /* Replace H by hVecs' * H * hVecs */
   /* ------------------------------- */

   compute_submatrix(hVecs, restartSize, ldhVecs, H, basisSize, ldH, H, ldH,
      rwork, rworkSize);

   /* -------------------------------------- */
   /* Quick exit if the target has changed   */
   /* -------------------------------------- */

   if (primme->targetShifts[min(primme->numTargetShifts-1, numConvergedBeforeRestart)]   
         != primme->targetShifts[min(primme->numTargetShifts-1, numConverged)]) {

      ret = update_Q_@(pre)primme(V, nLocal, ldV, W, ldW, Q, ldQ, R, ldR,
            primme->targetShifts[min(primme->numTargetShifts-1, numConverged)], 0,
            restartSize, rwork, rworkSize, machEps, primme);
      if (ret != 0) return ret;

      ret = solve_H_Ref_@(pre)primme(H, ldH, hVecs, newldhVecs, hU, newldhU, hSVals, 
            R, ldR, hVals, restartSize, rworkSize, rwork, primme);
      if (ret != 0) return ret;

      return 0;

   }
   
   /* -------------------------------------------------------------------- */
   /* During the restart V and W has replaced by V*hVecs and W*hVecs.      */
   /* Currently the QR decomposition correspond to W before restarting.    */
   /* To update the QR decomposition to the new W, replace Q by Q*Qn and   */
   /* R by Rn where Qn and Rn are the QR decomposition of R*hVecs = Qn*Rn. */
   /* -------------------------------------------------------------------- */

   /* -------------------------------------------------------------------- */
   /* R(indexOfPrevVecs:) = R * hVecs(indexOfPrevVecs:)                    */
   /* -------------------------------------------------------------------- */

   Num_copy_matrix_@(pre)primme(&hVecs[indexOfPreviousVecs*ldhVecs], basisSize,
      numPrevRetained, ldhVecs, rwork, basisSize);
   Num_trmm_@(pre)primme("L", "U", "N", "N", basisSize, numPrevRetained, tpone,
      R, ldR, rwork, basisSize);
   Num_copy_matrix_@(pre)primme(rwork, basisSize, numPrevRetained, basisSize,
      &R[ldR*indexOfPreviousVecs], ldR);

   /* -------------------------------------------------------------------- */
   /* hVecs(0:indexOfPrevVecs) are the perm right singular vectors of R.   */
   /* So R*hVecs(0:indexOfPrevVecs) = U(perm)*diag(hSVals(perm)).          */
   /* -------------------------------------------------------------------- */

   for (j=0; j < indexOfPreviousVecs; j++) {
      for (i=0; i < primme->maxBasisSize; i++) {
         R[ldR*j+i] = tzero;
      }
      *(double*)&R[ldR*j+j] = hSVals[perm[j]];
   }

   permute_vecs_@(pre)(hU, basisSize, basisSize, ldhU, perm, rwork, iwork);

   /* -------------------------------------------------------------------- */
   /* Compute the QR decomposition of R(indexOfPrevVecs:restartSize-1)     */
   /* -------------------------------------------------------------------- */

   Num_geqrf_@(pre)primme(basisSize, numPrevRetained, &R[ldR*indexOfPreviousVecs],
      ldR, rwork, rwork+basisSize, rworkSize-basisSize, &ret);
   assert(ret == 0);

   /* -------------------------------------------------------------------- */
   /* Place the Q factor besides hU and apply permutation Vperm            */
   /* -------------------------------------------------------------------- */

   Num_copy_matrix_@(pre)primme(&R[ldR*indexOfPreviousVecs], basisSize,
      numPrevRetained, ldR, &hU[ldhU*indexOfPreviousVecs], ldhU);
   Num_orgqr_@(pre)primme(basisSize, numPrevRetained, numPrevRetained,
      &hU[ldhU*indexOfPreviousVecs], ldhU, rwork, rwork+basisSize,
      rworkSize-basisSize, &ret);
   assert(ret == 0);

   /* -------------------------------------------------------------------- */
   /* Move the R factor to diagonal part of                                */
   /* R(indexOfPrevVecs:restartSize-1) and zero the lower triangular part  */
   /* of R(indexOfPrevVecs:restartSize-1, indexOfPrevVecs:restartSize-1).  */
   /* -------------------------------------------------------------------- */

   Num_copy_trimatrix_@(pre)primme(&R[ldR*indexOfPreviousVecs], numPrevRetained,
      numPrevRetained, ldR, 0 /* copy upper part */, 0,
      &R[ldR*indexOfPreviousVecs+indexOfPreviousVecs],
      primme->maxBasisSize, 1 /* zero lower part */);

   /* -------------------------------------------------------------------- */
   /* Zero R(0:indexOfPrevVecs-1, indexOfPrevVecs:restartSize-1)           */
   /* -------------------------------------------------------------------- */

   for (j=indexOfPreviousVecs; j<indexOfPreviousVecs+numPrevRetained; j++) {
      for (i=0; i < indexOfPreviousVecs; i++) {
         R[ldR*j+i] = tzero;
      }
   }
   for (j=indexOfPreviousVecs+numPrevRetained; j < restartSize; j++) {
      for (i=0; i <= j; i++) {
         R[ldR*j+i] = tzero;
      }
      *(double*)&R[ldR*j+j] = hSVals[perm[j]];
   }

   /* ----------------------------------- */
   /* Restart Q by replacing it with Q*hU */
   /* ----------------------------------- */

   Num_update_VWXR_@(pre)(Q, NULL, nLocal, basisSize, ldQ, hU, restartSize,
      ldhU, NULL,
      Q, 0, restartSize, ldQ,
      NULL, 0, 0, 0,
      NULL, 0, 0, 0,
      NULL, 0, 0, 0,
      NULL, 0, 0, 0, NULL,
      NULL, 0, 0,
      rwork, rworkSize, primme);

   /* ----------------------------------------------------------------------- */
   /* Given the above R, we know the right vectors will be the standard       */
   /* basis vectors if no previous coefficient vectors are retained           */
   /* ----------------------------------------------------------------------- */

   for (i=0, orderedIndexOfPreviousVecs=restartSize; i<restartSize; i++) {
      if (Vperm[i] == indexOfPreviousVecs) {
         orderedIndexOfPreviousVecs = i;
         break;
      }
   }
   assert(orderedIndexOfPreviousVecs != restartSize || indexOfPreviousVecs >= restartSize);

   for (j=0; j < restartSize; j++) {
      for (i=0; i < restartSize; i++) {
          hVecs[newldhVecs*j+i] = tzero;
          hU[newldhU*j+i] = tzero;
      }
      if (j < orderedIndexOfPreviousVecs || j >= orderedIndexOfPreviousVecs+numPrevRetained) {
         hVecs[newldhVecs*j+Vperm[j]] = tpone;
         hU[newldhU*j+Vperm[j]] = tpone;
      }
   }

   permute_vecs_d(hSVals, 1, basisSize, 1, perm, (double*)rwork, iwork);
   permute_vecs_d(hSVals, 1, restartSize, 1, Vperm, (double*)rwork, iwork);
   permute_vecs_d(hVals, 1, restartSize, 1, Vperm, (double*)rwork, iwork);

   /* ---------------------------------------------------------------------- */
   /* If coefficient vectors from the previous iteration have been retained, */
   /* then insert the computed overlap matrix into the restarted R and solve */ 
   /* the resulting projected problem.                                       */
   /* ---------------------------------------------------------------------- */

   if (numPrevRetained > 0) {
      ret = solve_H_Ref_@(pre)primme(
         &H[ldH*indexOfPreviousVecs+indexOfPreviousVecs], ldH,
         &hVecs[newldhVecs*orderedIndexOfPreviousVecs+orderedIndexOfPreviousVecs], newldhVecs,
         &hU[newldhU*orderedIndexOfPreviousVecs+orderedIndexOfPreviousVecs], newldhU,
         &hSVals[orderedIndexOfPreviousVecs],
         &R[ldR*indexOfPreviousVecs+indexOfPreviousVecs], ldR,
         &hVals[orderedIndexOfPreviousVecs], numPrevRetained,
         rworkSize, rwork, primme);

      /* TODO: reorder all hVals when numPacked != 0 */

      if (ret != 0) {
         primme_PushErrorMessage(Primme_restart_h, Primme_insert_submatrix, 
            ret, __FILE__, __LINE__, primme);
         return INSERT_SUBMATRIX_FAILURE;
      }
   }

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
 * rwork      double work array of size maxBasisSize^2. Use for hVec swapping.
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


int dtr_@(pre)(int numLocked, @(type) *hVecs, double *hVals, int *flags, 
  int basisSize, int numFree, int *iev, @(type) *rwork, primme_params *primme)
{

   int i;                 /* Loop variable */
   int l, lOpt, lMin;     /* Determine how many left side vectors to retain   */
   int r, rOpt;           /* Determine how many right side vectors to retain  */
   int maxIndex;          /* basisSize - 1                                    */
   int restartSize;       /* The new restart size                             */
   double currentRitzVal; /* The current Ritz value the solver is computing   */
   double newVal, optVal; /* Used to find the optimum gap ratio               */

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

   Num_@(pre)copy_@(pre)primme(i*basisSize, &hVecs[basisSize*lOpt], 1, rwork, 1);
   Num_@(pre)copy_@(pre)primme(rOpt*basisSize, &hVecs[basisSize*(basisSize-rOpt)], 1,
      &hVecs[basisSize*lOpt], 1);
   Num_@(pre)copy_@(pre)primme(i*basisSize, rwork, 1, &hVecs[basisSize*restartSize], 1);

   /* Do the same with the eigenvalues of H */

   Num_dcopy_primme(i, &hVals[lOpt], 1, (double *) rwork, 1);
   Num_dcopy_primme(rOpt, &hVals[(basisSize-rOpt)], 1, &hVals[lOpt], 1);
   Num_dcopy_primme(i, (double *) rwork, 1, &hVals[restartSize], 1);

   /* Set only those flags lower than restartSize. The rest will be reset */
   for (i = 0; i < rOpt; i++) {
      flags[lOpt + i] = flags[basisSize-rOpt + i];
   }

   if (primme->printLevel >= 5 && primme->procID == 0) {
      fprintf(primme->outputFile,"DTR restart size: %d L: %d R: %d\n", 
         restartSize, lOpt, rOpt);
   }

   reset_flags_@(pre)primme(flags, restartSize, primme->maxBasisSize);
   return restartSize;

}

 
/*******************************************************************************
 * Subroutine compute_submatrix - This subroutine computes the nX x nX submatrix
 *    R = X'*H*X, where H stores the upper triangular part of a symmetric matrix.
 *    
 * Input parameters
 * ----------------
 * X        The coefficient vectors retained from the previous iteration
 *
 * nX       Number of columns of X
 *
 * H        Matrix
 *
 * nH       Dimension of H
 *
 * ldH      Leading dimension of H
 *
 * rwork    Work array.  Must be of size nH x nX
 *
 * lrwork   Length of the work array
 *
 * ldR      Leading dimension of R
 *
 * Output parameters
 * -----------------
 * R - nX x nX matrix computed 
 *
 ******************************************************************************/

static int compute_submatrix(@(type) *X, int nX, int ldX, 
   @(type) *H, int nH, int ldH, @(type) *R, int ldR,
   @(type) *rwork, int lrwork) {

   @(type) tpone = @(tpone), tzero = @(tzero);

   /* Return memory requirement */
   if (X == NULL) {
      return nH*nX;
   }

   if (nH == 0 || nX == 0) return 0;

   assert(lrwork >= nH*nX);

   Num_symm_@(pre)primme("L", "U", nH, nX, tpone, H, ldH, X, ldX, tzero, rwork, nH);
   
   Num_gemm_@(pre)primme("C", "N", nX, nX, nH, tpone, X, ldX, rwork, nH, tzero, R, 
      ldR);

   return 0;
}


/*******************************************************************************
 * Subroutine reset_flags - Marks a series of Ritz values as unconverged
 *
 * INPUT PARAMETERS
 * ----------------
 * first Index of first Ritz value to be marked
 * last  Index of last Ritz value to be marked 
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * flags The flags indicating which Ritz values have converged
 ******************************************************************************/
     
void reset_flags_@(pre)primme(int *flags, int first, int last) {

   int i;  /* Loop variable */

   for (i = 0; i <= last; i++) {
      flags[i] = UNCONVERGED;
   }

}
