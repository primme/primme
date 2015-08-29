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
#include "primme.h"
#include "const.h"
#include "restart_d.h"
#include "restart_private_d.h"
#include "ortho_d.h"
#include "factorize_d.h"
#include "update_projection_d.h"
#include "numerical_d.h"


/*******************************************************************************
 * Subroutine: restart - This routine replaces V with V*c, some subset
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
 *                  3*primme->maxPrevRetain, primme->maxBasisSize^2)
 *                  
 * primme           Structure containing various solver parameters
 *
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
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
 
int restart_dprimme(double *V, double *W, double *H, double *hVecs,
   double *hVals, int *flags, int *iev, double *evecs, double *evecsHat, 
   double *M, double *UDU, int *ipivot, int basisSize, int numConverged, 
   int *numConvergedStored, int numLocked, int numGuesses, 
   double *previousHVecs, int numPrevRetained, double machEps, 
   double *rwork, int rworkSize, primme_params *primme) {
  
   int numFree;             /* The number of basis vectors to be left free    */
   int numPacked;           /* The number of coefficient vectors moved to the */
                            /* end of the hVecs array.                        */
   int restartSize;         /* The number of vectors to restart with          */
   int indexOfPreviousVecs=0; /* Position within hVecs array the previous       */
                            /* coefficient vectors will be stored             */
   int i, n, eStart;        /* various variables                              */
   int ret;                 /* Return value                                   */

   numPacked = 0;

   /* --------------------------------------------------------------------- */
   /* If dynamic thick restarting is to be used, then determine the minimum */
   /* number of free spaces to be maintained and call the DTR routine.      */
   /* The DTR routine will determine how many coefficient vectors from the  */
   /* left and right of H-spectrum to retain at restart. If DTR is not used */
   /* then set the restart size to the minimum restart size.                */
   /* --------------------------------------------------------------------- */

   if (primme->restartingParams.scheme == primme_dtr) {
      numFree = numPrevRetained+max(3, primme->maxBlockSize);
      restartSize = dtr(numLocked, hVecs, hVals, flags, basisSize, numFree, 
                        iev, rwork, primme);
   }
   else {
      restartSize = min(basisSize, primme->minRestartSize);
   }

   /* ----------------------------------------------------------------------- */
   /* If locking is engaged, then swap coefficient vectors corresponding to   */
   /* converged Ritz vectors to the end of the hVecs(:, restartSize) subarray.*/
   /* This allows the converged Ritz vectors to be stored contiguously in     */
   /* memory after restart.  This significantly reduces the amount of data    */
   /* movement the locking routine would have to perform otherwise.           */
   /* The following function also covers some limit cases where restartSize   */
   /* plus 'to be locked' and previous Ritz vectors may exceed the basisSize  */
   /* ----------------------------------------------------------------------- */

   if (primme->locking) {
      numPacked = pack_converged_coefficients(&restartSize, basisSize, 
         &numPrevRetained, numLocked, numGuesses, hVecs, hVals, flags, primme);
   }

   /* ----------------------------------------------------------------------- */
   /* Restarting with a small number of coefficient vectors from the previous */
   /* iteration can be retained to accelerate convergence.  The previous      */
   /* coefficient vectors must be combined with the current coefficient       */
   /* vectors by first orthogonalizing the previous ones versus the current   */
   /* restartSize ones.  The orthogonalized previous vectors are then         */
   /* inserted into the hVecs array at hVecs(:,indexOfPreviousVecs).          */
   /* ----------------------------------------------------------------------- */

   if (numPrevRetained > 0) {
      indexOfPreviousVecs = combine_retained_vectors(hVals, flags, hVecs,
         basisSize, &restartSize, numPacked, previousHVecs, 
         &numPrevRetained, machEps, rwork, primme);
   }

   /* -------------------------------------------------------- */
   /* Restart V by replacing it with the current Ritz vectors. */
   /* -------------------------------------------------------- */

   restart_X(V, hVecs, primme->nLocal, basisSize, restartSize, rwork,rworkSize);
   
   /* ------------------------------------------------------------ */
   /* Restart W by replacing it with W times the eigenvectors of H */
   /* ------------------------------------------------------------ */

   restart_X(W, hVecs, primme->nLocal, basisSize, restartSize, rwork,rworkSize);

   /* ---------------------------------------------------------------- */
   /* Because we have replaced V by the Ritz vectors, V'*A*V should be */
   /* diagonal with the Ritz values on the diagonal.  The eigenvectors */
   /* of the new matrix V'*A*V become the standard basis vectors.      */
   /* ---------------------------------------------------------------- */

   ret = restart_H(H, hVecs, hVals, restartSize, basisSize, previousHVecs, 
      numPrevRetained, indexOfPreviousVecs, rworkSize, rwork, primme);

   if (ret != 0) {
      primme_PushErrorMessage(Primme_restart, Primme_restart_h, ret, __FILE__, 
         __LINE__, primme);
      return RESTART_H_FAILURE;
   }

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

   if (!primme->locking && primme->correctionParams.maxInnerIterations != 0 && 
        numConverged > 0 &&
        (primme->correctionParams.projectors.LeftQ ||
         primme->correctionParams.projectors.RightQ )  ) {

       n = primme->nLocal;
       *numConvergedStored = 0;
       eStart = primme->numOrthoConst;

       for (i=0;i<primme->numEvals;i++) {
           if (flags[i] == CONVERGED) {
              if (*numConvergedStored < numConverged) {
                 Num_dcopy_dprimme(n, &V[i*n], 1, 
                              &evecs[(eStart+*numConvergedStored)*n], 1);
                 (*numConvergedStored)++;
              }
           } /* if converged */
       } /* for */
       if (*numConvergedStored != numConverged) {
          if (primme->printLevel >= 1 && primme->procID == 0) {
             fprintf(primme->outputFile, 
             "Flags and converged eigenpairs do not correspond %d %d\n",
                numConverged, *numConvergedStored);
          }
          return PSEUDOLOCK_FAILURE;
       }

      /* Update also the M = K^{-1}evecs and its udu factorization if needed */
      if (UDU != NULL) {

         apply_preconditioner_block(&evecs[eStart*n], &evecsHat[eStart*n], 
                                    numConverged, primme );
         /* rwork must be maxEvecsSize*numEvals! */
         update_projection_dprimme(evecs, evecsHat, M, eStart*n,
           primme->numOrthoConst+primme->numEvals, numConverged, rwork, primme);

         ret = UDUDecompose_dprimme(M, UDU, ipivot, eStart+numConverged, 
                         rwork, rworkSize, primme);
         if (ret != 0) {
            primme_PushErrorMessage(Primme_lock_vectors,Primme_ududecompose,ret,
               __FILE__, __LINE__, primme);
            return UDUDECOMPOSE_FAILURE;
         }
      } /* if UDU factorization is needed */
   } /* if this pseudo locking should take place */

   return restartSize;
}


/*******************************************************************************
 * Subroutine restart_X - This subroutine computes X*hVecs and places 
 *    the result in X.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal       Number of rows of V assigned to the node
 *
 * basisSize    Current size of the basis V
 *
 * restartSize  Number of Ritz vectors V/W will be restarted with 
 *
 * rwork        Work array that must be at least of size restartSize
 *
 * rworkSize    The size availble in rwork. Matrix multiply blocks of X with 
 *              hVecs, producing blocks of the new X of size
 *              (AvailRows * restartSize) = rworkSize
 *              Therefore rworkSize must be at least restartSize.
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * X      Holds either V or W before and after restarting
 *
 * hVecs  The eigenvectors of V'*A*V before and after restarting
 *
 ******************************************************************************/
  
static void restart_X(double *X, double *hVecs, int nLocal, 
   int basisSize, int restartSize, double *rwork, int rworkSize) {

   int i, k;  /* Loop variables */
   int AvailRows = min(rworkSize/restartSize, nLocal);
   double tpone = +1.0e+00, tzero = +0.0e+00;
   i = 0;

   while (i < nLocal) {
      /* Block matrix multiply */
      Num_gemm_dprimme("N", "N", AvailRows, restartSize, basisSize, tpone,
         &X[i], nLocal, hVecs, basisSize, tzero, rwork, AvailRows );

      /* Copy the result in the desired location of X */
      for (k=0; k < restartSize; k++) {
         Num_dcopy_dprimme(AvailRows, &rwork[AvailRows*k],1, &X[i+nLocal*k], 1);
      }
      i = i+AvailRows;
      AvailRows = min(AvailRows, nLocal-i);
   }
}


/*******************************************************************************
 * Function restart_H - This routine is used to recompute H = V'*A*V once V 
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
 *               max(basisSize*numPrevRetained, 3*numPrevRetained) 
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

static int restart_H(double *H, double *hVecs, double *hVals, 
   int restartSize, int basisSize, double *previousHVecs, 
   int numPrevRetained, int indexOfPreviousVecs, int rworkSize, 
   double *rwork, primme_params *primme) {

   int i, j;          /* Loop variables                                       */
   int workSpaceSize; /* Workspace size needed by insert_submatrix            */
   int ret;           /* Return value                                         */
   double *subMatrix;/* Contains the submatrix previousHVecs'*H*previousHvecs*/
   double *workSpace;/* Workspace size needed                              */
   double tpone = +1.0e+00, tzero = +0.0e+00;             /*constants*/
   
   /* ---------------------------------------------------------------------- */
   /* If coefficient vectors from the previous iteration were retained, then */
   /* set up work space for computing the numPrevRetained x numPrevRetained  */
   /* submatrix, then compute the submatrix.                                 */
   /* ---------------------------------------------------------------------- */

   if (numPrevRetained > 0) {
      subMatrix = rwork;
      workSpace = &rwork[numPrevRetained*numPrevRetained];
      workSpaceSize = rworkSize - numPrevRetained*numPrevRetained;
      compute_submatrix(previousHVecs, numPrevRetained, H, basisSize, 
         primme->maxBasisSize, subMatrix, workSpace);
   }
      
   /* ----------------------------------------------------------------- */
   /* V*hVecs yields a diagonal matrix composed of the Ritz values of A */
   /* with respect to V.  Set H to a diagonal matrix with the Ritz      */
   /* values on the diagonal.                                           */
   /* ----------------------------------------------------------------- */

   for (j=0; j < restartSize; j++) {
      for (i=0; i < j; i++) {
         H[primme->maxBasisSize*j+i] = tzero;
      }
      H[primme->maxBasisSize*j+j] = hVals[j];
   }

   /* --------------------------------------------------------------------- */
   /* Given the above H, we know the eigenvectors of H will be the standard */
   /* basis vectors if no previous coefficient vectors are retained         */
   /* --------------------------------------------------------------------- */

   for (j=0; j < restartSize; j++) {
      for (i=0; i < j; i++) {
          hVecs[restartSize*j+i] = tzero;
          hVecs[restartSize*i+j] = tzero;
      }
      hVecs[restartSize*j+j] = tpone;
   }      

   /* ---------------------------------------------------------------------- */
   /* If coefficient vectors from the previous iteration have been retained, */
   /* then insert the computed overlap matrix into the restarted H and solve */ 
   /* the resulting eigenproblem for the resulting H.                        */
   /* ---------------------------------------------------------------------- */

   if (numPrevRetained > 0) {
      ret = insert_submatrix(H, hVals, hVecs, restartSize, subMatrix, 
        numPrevRetained, indexOfPreviousVecs, workSpaceSize, workSpace, primme);

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


static int dtr(int numLocked, double *hVecs, double *hVals, int *flags, 
  int basisSize, int numFree, int *iev, double *rwork, primme_params *primme)
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
         if (  (flags[l] == CONVERGED || flags[l] == PRACTICALLY_CONVERGED) 
             && (numLocked + l < primme->numEvals)) {
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

   Num_dcopy_dprimme(i*basisSize, &hVecs[basisSize*lOpt], 1, rwork, 1);
   Num_dcopy_dprimme(rOpt*basisSize, &hVecs[basisSize*(basisSize-rOpt)], 1,
      &hVecs[basisSize*lOpt], 1);
   Num_dcopy_dprimme(i*basisSize, rwork, 1, &hVecs[basisSize*restartSize], 1);

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

   reset_flags_dprimme(flags, restartSize, primme->maxBasisSize);
   return restartSize;

}

/******************************************************************************
 * Function pack_converged_coefficients - This function is called before
 *    restart so that the coefficient vectors (eigenvectors of H) are packed
 *    towards the end of the hVecs[0..restartSize-1] subarray.  This ensures 
 *    that the Ritz vectors to be locked will be contiguous in memory.  This 
 *    allows the basis to be updated more easily with new initial guesses while
 *    the converged Ritz vectors are locked. 
 *    The routine is ONLY called if locking is used.
 *
 *    Note: to guarantee restartSize vectors will be in the basis after lacking
 *    and restart have occurred, we should restart with 
 *    restartSize+numFlagged-anyInitialGuessesStillAvailable. 
 *    All converged vectors are retained beyond that.
 *
 *    Note: if too many vectors have converged, and adding numPrevRetained
 *    ones may increase the basisSize, then we simply do not add any previous
 *    ones, but compute all the Ritz vectors throwing nothing away this time
 *
 * Input parameters
 * ----------------
 * basisSize     The current size of the basis
 *
 * numLocked     The number of vectors that have been locked
 *
 * numGuesses    Number of initial guesses remaining
 *
 * primme        Structure containing various solver parameters
 *
 * 
 * Input/Output parameters
 * -----------------------
 * restartSize   The number of (nonconverged) vectors to restart the basis with
 *
 * numPrevRetained The number of vectors from the previous iteration to include
 *
 * hVecs         The eigenvectors of the projection H
 *
 * hVals         The eigenvalues of the projection H
 *
 * flag          Array indicating the convergence status of each the
 *               basisSize current Ritz vectors.  The array is of size 
 *               basisSize.
 *
 *
 * Return value
 * ------------
 * int  The number of vectors to be locked
 ******************************************************************************/

static int pack_converged_coefficients(int *restartSize, int basisSize, 
   int *numPrevRetained, int numLocked, int numGuesses, double *hVecs, 
   double *hVals, int *flag, primme_params *primme) {

   int i;            /* Loop variable                        */
   int left, right;  /* Search indices                       */
   int itemp;        /* Temporary variable used for swapping */
   int numFlagged;   /* Number of target converged Ritz vectors that have */
                     /* converged since the last time this function was   */
                     /* called.                                           */
   double dtemp;     /* Temporary variable used for swapping */


   /* ---------------------------------------------- */
   /* Only converged target vectors should be locked */
   /* ---------------------------------------------- */

   for (i = 0, numFlagged = 0; i < basisSize; i++) {
      /* Make sure the vector is converged and it's a target Ritz vector */

      if (flag[i] == CONVERGED && (numLocked + i < primme->numEvals)) {
         flag[i] = LOCK_IT;
         numFlagged++;
      }
      else if (flag[i] == PRACTICALLY_CONVERGED && 
                                  (numLocked + i < primme->numEvals)) {
         flag[i] = UNCONDITIONAL_LOCK_IT;
         numFlagged++;
      }
   }

   /* ----------------------------------------------------------- */
   /* Special case: If (basisSize+numLocked) is the entire space, */
   /* then everything should be converged. Do not test, just flag */
   /* everything as converged to be locked.                       */
   /* ----------------------------------------------------------- */

   if (basisSize + numLocked + primme->numOrthoConst == primme->n) {
      for (numFlagged = 0; 
           numFlagged < min(primme->numEvals-numLocked, basisSize); 
           numFlagged++) {
      flag[numFlagged] = LOCK_IT;
      }
   }

   /* ------------------------------------------------------------------- */
   /* Redefine restartSize so that we have at least restartSize in the    */
   /* basis when we restart and after initial guesses are substituted in  */ 
   /* If more than primme.restartSize have converged we need to keep them */
   /* If necessary do not throw any vectors (ie., if the above>basisSize) */
   /* In that case, there is no need to keep any previous Ritz vectors.   */
   /* ------------------------------------------------------------------- */

   itemp = min(numFlagged, numGuesses);
   if (itemp >= *restartSize)
      *restartSize = numFlagged;
   else
      *restartSize = *restartSize + numFlagged - itemp;

   if (*restartSize + *numPrevRetained >= basisSize) {
      *restartSize = basisSize;
      *numPrevRetained = 0;
   }

   /* ------------------------------------------------------------------ */
   /* The right index starts at the end of the flags[0..restartSize-1]   */
   /* subarray and stops decreasing when it finds a vector that is not   */
   /* to be locked.  The left index is then used to find a replacement   */
   /* to swap with.  This replacement must be a vector targeted for      */
   /* locking.  If no replacement can be found, the packing is finished. */
   /* ------------------------------------------------------------------ */

   right = *restartSize - 1;

   while (right >= 0) {

      /* Find a vector that is not to be locked */
 
      while (right >= 0 && (flag[right] == LOCK_IT || 
                            flag[right] == UNCONDITIONAL_LOCK_IT) ) {
         right--;
      }

      left = right - 1;

      /* Find a vector that is to be locked */

      while (left >= 0 && flag[left] != LOCK_IT && 
                          flag[left] != UNCONDITIONAL_LOCK_IT) { 
         left--;
      }

      /* If no such vector could be found, packing is complete */

      if (left < 0) {
         return numFlagged;
      }

      /* Swap the coefficient vectors corresponding to left and right. */
      Num_swap_dprimme(basisSize, &hVecs[basisSize*left], 1, 
         &hVecs[basisSize*right], 1);

      /* Swap the Ritz values */
      dtemp = hVals[left];
      hVals[left] = hVals[right];
      hVals[right] = dtemp;

      /* Swap the flag values */
      itemp = flag[left];
      flag[left] = flag[right];
      flag[right] = itemp;
 
   }

   return numFlagged;

}

/*******************************************************************************
 * Function combine_retained_vectors -- This function combines the current
 *   coefficient vectors with the ones retained from the previous iteration.
 *   The retained coefficient vectors are first orthogonalized versus themselves
 *   and the current coefficient vectors.  The previous coefficients are then
 *   copied to hVecs(:,restartSize-numPacked).  This is because the locking 
 *   routine (if locking is engaged) requires coefficient vectors coresponding
 *   to converged Ritz vectors be stored at the tail of the hVecs array.
 *
 * Input parameters
 * ----------------
 * basisSize    The current basis size
 *
 * numPacked    The number of coefficient vectors corresponding to converged
 *              Ritz vectors that are to be locked.
 * 
 * previousHVecs  The coefficients retained from the previous iteration
 *
 * rwork        Real workspace needed by the orthogonalization routine.
 *              basisSize is a sufficient size.
 *
 *
 * Input/output parameters
 * -----------------------
 * hVals  The eigenvalues (Ritz values) of the projection matrix H
 *
 * flags  Array indicating the convergence of the Ritz values
 *
 * hVecs  The eigenvectors (coefficient vectors) of the projection matrix H
 *
 * restartSize  The current restart size
 *
 * numPrevRetained  The actual number of previous coefficient vectors to be
 *                  retained at restart.  May be reduced due to 
 *                  orthogonalization difficulties, or because they do not 
 *                  all fit within the basisSize.
 *
 * Return value
 * ------------
 * The index within hVecs where the retained vectors were inserted.
 *
 ******************************************************************************/

static int combine_retained_vectors(double *hVals, int *flags, double *hVecs,
   int basisSize, int *restartSize, int numPacked, double *previousHVecs, 
   int *numPrevRetained, double machEps, double *rwork, 
   primme_params *primme) {

   int i;                    /* Loop variable */
   int indexOfPreviousVecs;  /* The index within hVecs where the previous */
                             /* coefficients are inserted.                */

   if (primme->printLevel >= 5 && primme->procID == 0) {
      fprintf(primme->outputFile,
  "combine: basisSize: %d restartSize: %d numPrevRetained: %d numPacked: %d \n",
         basisSize, *restartSize, *numPrevRetained, numPacked);
   }

   /* ------------------------------------------------------------------ */
   /* Orthogonalize the coefficents from the previous iteration with the */
   /* current coefficients.  This may annihilate some of the previous    */
   /* vectors, if the have much overlap with the current vectors.  Thus, */
   /* numPrevRetained may be reduced.                                    */
   /* ------------------------------------------------------------------ */

   *numPrevRetained = ortho_retained_vectors_dprimme(hVecs, basisSize, 
      *restartSize, previousHVecs, *numPrevRetained, machEps, rwork);

   /* --------------------------------------------------------------------- */
   /* If locking is engaged and there exist retained previous coefficent    */
   /* vectors, then move the last numPacked vectors forward numPrevRetained */
   /* spaces.  This will make room for the orthogonalized previous vectors. */
   /* --------------------------------------------------------------------- */

   if (primme->locking && *numPrevRetained > 0) {
      indexOfPreviousVecs = *restartSize - numPacked;
      /* WARNING: dcopy's -1 step is not implemented appropriately in Mac's 
       * Veclib (and other libs).  Perform vector by vector instead.
      Num_dcopy_dprimme(basisSize*numPacked, 
         &hVecs[basisSize*(*restartSize-numPacked)], -1, 
         &hVecs[basisSize*(*restartSize-numPacked+*numPrevRetained)], -1);
      */
      for (i=1;i<=numPacked;i++) Num_dcopy_dprimme(basisSize,
                     &hVecs[basisSize*(*restartSize-i)],1,
                     &hVecs[basisSize*(*restartSize+*numPrevRetained-i)],1);

      /* Move the Ritz values and flag values forward as well */

      for (i = numPacked-1; i >= 0; i--) {
         hVals[*restartSize-numPacked+*numPrevRetained+i] = 
            hVals[*restartSize-numPacked+i];      
         flags[*restartSize-numPacked+*numPrevRetained+i] = 
            flags[*restartSize-numPacked+i];      
      }

   }
   else {
      indexOfPreviousVecs = *restartSize;
   }
   
   /* ----------------------------------------------------------------*/
   /* Copy the orthogonalized previous coefficents to the hVecs array */
   /* ----------------------------------------------------------------*/

   Num_dcopy_dprimme(basisSize*(*numPrevRetained), previousHVecs, 1, 
      &hVecs[basisSize*indexOfPreviousVecs], 1);

   /* Initialize the Ritz and flag values */
   for (i = 0; i < *numPrevRetained; i++) {
      hVals[indexOfPreviousVecs+i] = 0.0L;
      flags[indexOfPreviousVecs+i] = UNCONVERGED;
   }

   
   if (primme->printLevel >= 5 && primme->procID == 0) {
      fprintf(primme->outputFile, "numPrevRetained: %d restartSize: %d\n", 
         *numPrevRetained, *restartSize);
   }

   /* Increase the restart size with the previous vectors added. */

   *restartSize = *restartSize + *numPrevRetained;

   return indexOfPreviousVecs;
}
   
      
/*******************************************************************************
 * Subroutine compute_submatrix - This subroutine computes the 
 *    numPrevRetained x numPrevRetained submatrix 
 *    previousHVecs'*H*previousHVecs.
 *    
 * Input parameters
 * ----------------
 * previousHVecs   The coefficient vectors retained from the previous iteration
 *
 * numPrevRetained  Number of previous vectors retained
 *
 * H               The projection matrix V'*A*V
 *
 * basisSize       The current size of the basis and dimension of H
 *
 * maxBasisSize    The maximum basis size and leading dimension of H
 *
 * rwork           Work array.  Must be of size basisSize x numPrevRetained
 *
 * 
 * Output parameters
 * -----------------
 * subMatrix - numPrevRetained x numPrevRetained submatrix to be computed    
 *
 ******************************************************************************/

static void compute_submatrix(double *previousHVecs, int numPrevRetained, 
   double *H, int basisSize, int maxBasisSize, double *subMatrix, 
   double *rwork) {

   double tpone = +1.0e+00, tzero = +0.0e+00;

   Num_symm_dprimme("L", "U", basisSize, numPrevRetained, tpone, H, 
      maxBasisSize, previousHVecs, maxBasisSize, tzero, rwork, basisSize);
   
   Num_gemm_dprimme("C", "N", numPrevRetained, numPrevRetained, basisSize,
      tpone, previousHVecs, basisSize, rwork, basisSize, tzero, subMatrix, 
      numPrevRetained);
}


/*******************************************************************************
 * Function insert_submatrix -- This function inserts the submatrix
 *    previousHVecs'*H*previousHVecs into the restarted H at 
 *    H(indexOfRetainedVecs, indexOfRetainedVecs).  It then solves
 *    the eigenproblem for the submatrix and constructs the numPrevRetained
 *    eigenvectors of H corresponding to the submatrix.
 *
 * Input parameters
 * ----------------
 * restartSize  The basis size after restart
 *
 * numPrevRetained  Number of coefficient vectors retained from the previous
 *                  iteration
 * 
 * indexOfPreviousVecs  The position within hVecs where the coefficient vectors
 *                      corresponding to the submatrix will be inserted.  It
 *                      also indicates where the submatrix will be inserted
 *                      within H.
 *
 * rworkSize        Size of the work array.  Must be at least 3*maxPrevRetain
 *
 * rwork            Workspace needed to solve the eigenproblem for the submatrix
 *
 * primme             Structure containing various solve parameters
 *
 *
 * Input/Output parameters
 * -----------------------
 * H      The projection matrix V'*A*V
 *
 * hVals  The eigenvalues (Ritz values) of H
 *
 * hVecs  The eigenvectors of H
 *
 * submatrix        The matrix previousHVecs'*H*previousHVecs
 *
 ******************************************************************************/

static int insert_submatrix(double *H, double *hVals, double *hVecs, 
   int restartSize, double *subMatrix, int numPrevRetained, 
   int indexOfPreviousVecs, int rworkSize, double *rwork, 
   primme_params *primme) {

   int info;
   int i, j;

   /* ---------------------------------------------------------------------- */
   /* Copy the submatrix into H with the upper right corner of the submatrix */
   /* at H(indexOfPreviousVecs, indexOfPreviousVecs).                        */
   /* ---------------------------------------------------------------------- */

   for (j = indexOfPreviousVecs; j < indexOfPreviousVecs+numPrevRetained; j++) {
      for (i = indexOfPreviousVecs; i <= j; i++) {
         H[primme->maxBasisSize*j+i] = 
         subMatrix[numPrevRetained*(j-indexOfPreviousVecs)
                              +(i-indexOfPreviousVecs)];
      }
   }

   /* ----------------------------------------- */
   /* Solve the eigenproblrm for the submatrix. */
   /* ----------------------------------------- */
   Num_dsyev_dprimme("V", "U", numPrevRetained, subMatrix, numPrevRetained, 
      &hVals[indexOfPreviousVecs], rwork, rworkSize, &info);

   if (info != 0) {
      primme_PushErrorMessage(Primme_insert_submatrix, Primme_num_dsyev, info, 
         __FILE__, __LINE__, primme);
      return NUM_DSYEV_FAILURE;
   }

   /* ---------------------------------------------------------------------- */
   /* The eigenvectors of the submatrix are of dimension numPrevRetained  */
   /* and reside in the submatrix array after dsyev is called.  The       */
   /* eigenvectors of the new H corresponding to the submatrix are easily */
   /* constructed using the eigenvectors of the submatrix.                */
   /* ---------------------------------------------------------------------- */

   for (j = indexOfPreviousVecs; j < indexOfPreviousVecs+numPrevRetained; j++) {
      for (i = indexOfPreviousVecs; 
           i < indexOfPreviousVecs + numPrevRetained; i++) 
      {
         hVecs[restartSize*j+i] = 
         subMatrix[numPrevRetained*(j-indexOfPreviousVecs)+
                   (i-indexOfPreviousVecs)];
      }
   }

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
 * flag  The flags indicating which Ritz values have converged
 ******************************************************************************/
     
void reset_flags_dprimme(int *flag, int first, int last) {

   int i;  /* Loop variable */

   for (i = 0; i <= first-1; i++) {
      if (flag[i] != CONVERGED && flag[i] != PRACTICALLY_CONVERGED ) {
         flag[i] = UNCONVERGED;
      }
   }

   for (i=first; i <= last; i++) {
      flag[i] = UNCONVERGED;
   }

}

/*******************************************************************************
 * Subroutine apply_preconditioner_block - This subroutine applies the 
 *    preconditioner to a block of vectors v by computing: K^{-1}v
 *    (duplicated here as with correection.c to allow for static use)
 *
 * Input Parameters
 * ----------------
 * v         The vectors the preconditioner will be applied to.
 *
 * blockSize The number of vectors in the blocks v, result
 *
 * primme      Structure containing various solver parameters
 * 
 * Output parameters
 * -----------------
 * result    The result of the application of K^{-1}
 *
 ******************************************************************************/

static void apply_preconditioner_block(double *v, double *result, 
                int blockSize, primme_params *primme) {
         
   if (primme->correctionParams.precondition) {

      (*primme->applyPreconditioner)(v, result, &blockSize, primme);
      primme->stats.numPreconds += blockSize;
   }
   else {
      Num_dcopy_dprimme(primme->nLocal*blockSize, v, 1, result, 1);
   }

}

