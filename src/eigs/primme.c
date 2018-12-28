/*******************************************************************************
 * Copyright (c) 2018, College of William & Mary
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
 * File: primme.c
 *
 * Purpose - Real, SCALAR precision front end to the multimethod eigensolver
 *
 * For the moment please cite the following two papers: 
 *
 *  A. Stathopoulos, Nearly optimal preconditioned methods for hermitian
 *    eigenproblems under limited memory. Part I: Seeking one eigenvalue,
 *    Tech Report: WM-CS-2005-03, July, 2005. To appear in SISC.
 *  A. Stathopoulos and J. R. McCombs, Nearly optimal preconditioned methods
 *    for hermitian eigenproblems under limited memory. Part II: Seeking many
 *    eigenvalues, Tech Report: WM-CS-2006-02, June, 2006.
 *
 * Additional information on the algorithms appears in:
 *
 *  J. R. McCombs and A. Stathopoulos, Iterative Validation of Eigensolvers:
 *    A Scheme for Improving the Reliability of Hermitian Eigenvalue Solvers
 *    Tech Report: WM-CS-2005-02, July, 2005, to appear in SISC.
 *  A. Stathopoulos, Locking issues for finding a large number of eigenvectors 
 *    of hermitian matrices, Tech Report: WM-CS-2005-03, July, 2005.
 *
 *   Some papers to be listed here. For the moment contact the author:
 *
 *                       andreas@cs.wm.edu
 *
 ******************************************************************************/

#include "numerical.h"
#include "const.h"
#include "primme_interface.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "main_iter.h"
#include "auxiliary_eigs.h"
#endif

#ifdef SUPPORTED_TYPE

static int Sprimme_for_real(XREAL *evals, XSCALAR *evecs_, XREAL *resNorms,
            primme_params *primme, primme_context ctx);
static int check_input(
      REAL *evals, SCALAR *evecs, REAL *resNorms, primme_params *primme);
static void convTestFunAbsolute(double *eval, void *evec, double *rNorm, int *isConv,
   primme_params *primme, int *ierr);
static void default_monitor(void *basisEvals, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms, int *numConverged,
      void *lockedEvals, int *numLocked, int *lockedFlags, void *lockedNorms,
      int *inner_its, void *LSRes, const char *msg, double *time,
      primme_event *event, primme_params *primme, int *err);

#endif /* SUPPORTED_TYPE */

/*******************************************************************************
 * Subroutine Sprimme - This routine is a front end used to perform 
 *    error checking on the input parameters, perform validation, 
 *    and make the call to main_iter. 
 *
 *    Calling Sprimme with all evals, evecs, resNorms set to NULL
 *    returns the int and real memory required in the following primme fields:
 *            int primme->intWorkSize : bytes of int workspace needed
 *       long int primme->realWorkSize: bytes of real workspace needed
 * 
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * evals  Contains the converged Ritz values upon return.   Should be of size 
 *        primme->numEvals.
 * 
 * evecs  The local portions of the converged Ritz vectors.  The dimension of
 *        the array is at least primme->nLocal*primme->numEvals
 *
 * resNorms  The residual norms of the converged Ritz vectors.  Should be of 
 *           size primme->numEvals
 *  
 * primme  Structure containing various solver parameters and statistics
 *         See readme.txt for INPUT/OUTPUT variables in primme
 *
 * Return Value
 * ------------
 *  0 - Success
 *  1 - Reporting only memory space required
 * -1 - Failure to allocate workspace
 * -2 - Malloc failure in allocating a permutation integer array
 * -3 - main_iter encountered a problem
 * -4 ...-32 - Invalid input (parameters or primme struct) returned 
 *             by check_input()
 *
 ******************************************************************************/

int Sprimme(XREAL *evals, XSCALAR *evecs, XREAL *resNorms,
            primme_params *primme) {

#ifdef SUPPORTED_TYPE
   /* Generate context */

   primme_context ctx = primme_get_context(primme);

   /* Main call */

   int ret = Sprimme_for_real(evals, evecs, resNorms, primme, ctx);

   /* Free context */

   primme_free_context(ctx);
#else
   int ret = PRIMME_FUNCTION_UNAVAILABLE;
#endif

   return ret;
}


#ifdef SUPPORTED_TYPE

static int Sprimme_for_real(XREAL *evals, XSCALAR *evecs_, XREAL *resNorms,
            primme_params *primme, primme_context ctx) {

   int *perm;
   SCALAR *evecs = (SCALAR *)evecs_; /* Change type of evecs */

   /* zero out the timer */
   double t0 = primme_wTimer();

   /* Set some defaults for sequential programs */
   if (primme->numProcs <= 1 && evals != NULL && evecs != NULL &&
         resNorms != NULL) {
      primme->nLocal = primme->n;
      primme->procID = 0;
   }

   /* Set some defaults  */
   primme_set_defaults(primme);

   if (primme->orth == primme_orth_default) {
#ifdef USE_HOST
      /* The current code for block orthogonalization does not produce     */
      /* a machine precision orthonormal basis. So block orthogonalization */
      /* is used only when V'BV is computed explicitly.                    */

      if (primme->maxBlockSize > 1) {
         primme->orth = primme_orth_explicit_I;
      } else {
         primme->orth = primme_orth_implicit_I;
      }
#else

      /* Observed orthogonality issues finding the largest/smallest values in  */
      /* single precision. Computing V'*B*V and solving the projected problem  */
      /* V'AVx = V'BVxl mitigates the problem.                                 */

#  ifdef USE_FLOAT
      if (primme->projectionParams.projection == primme_proj_RR &&
                  (primme->target == primme_largest ||
                        primme->target == primme_smallest ||
                        primme->target == primme_largest_abs) ||
            primme->maxBlockSize > 1) {
         primme->orth = primme_orth_explicit_I;
      }
      else
#  endif
         primme->orth = primme_orth_implicit_I;
#endif
   }

   /* If we are free to choose the leading dimension of V and W, use    */
   /* a multiple of PRIMME_BLOCK_SIZE. This may improve the performance */
   /* of Num_update_VWXR_Sprimme and Bortho_block_Sprimme.              */

   if (primme->ldOPs == -1) {
      if (PRIMME_BLOCK_SIZE < INT_MAX) {
         primme->ldOPs = min(((primme->nLocal + PRIMME_BLOCK_SIZE - 1)
                  /PRIMME_BLOCK_SIZE)*PRIMME_BLOCK_SIZE, primme->nLocal);
      } else {
         primme->ldOPs = primme->nLocal;
      }
   }

   /* Deprecated input:                                              */
   if (evals == NULL && evecs == NULL && resNorms == NULL)
      return 0;

   /* Reset random number seed if inappropriate for DLARENV */
   /* Yields unique quadruples per proc if procID < 4096^3  */

   if (primme->iseed[0]<0 || primme->iseed[0]>4095) primme->iseed[0] = 
      primme->procID % 4096;
   if (primme->iseed[1]<0 || primme->iseed[1]>4095) primme->iseed[1] = 
      (int)(primme->procID/4096+1) % 4096;
   if (primme->iseed[2]<0 || primme->iseed[2]>4095) primme->iseed[2] = 
      (int)((primme->procID/4096)/4096+2) % 4096;
   if (primme->iseed[3]<0 || primme->iseed[3]>4095) primme->iseed[3] = 
      (2*(int)(((primme->procID/4096)/4096)/4096)+1) % 4096;

   /* Set default convTetFun  */

   if (!primme->convTestFun) {
      primme->convTestFun = convTestFunAbsolute;
      if (primme->eps == 0.0) {
         primme->eps = MACHINE_EPSILON*1e4;
      }
   }

   /* Set default monitor     */

   if (!primme->monitorFun) {
      primme->monitorFun = default_monitor;
   }

   /* Check primme input data for bounds, correct values etc. */

   int ret = check_input(evals, evecs, resNorms, primme);
   if (ret != 0) return ret;

   /* Cast evals and resNorms to HREAL */

   HREAL *evals0, *resNorms0;
   CHKERR(Num_matrix_astype_RHprimme(evals, 1, primme->numEvals, 1,
         PRIMME_OP_REAL, (void **)&evals0, NULL, PRIMME_OP_HREAL, 1 /* alloc */,
         0 /* not copy */, ctx));
   CHKERR(Num_matrix_astype_RHprimme(resNorms, 1, primme->numEvals, 1,
         PRIMME_OP_REAL, (void **)&resNorms0, NULL, PRIMME_OP_HREAL,
         1 /* alloc */, 0 /* not copy */, ctx));

   /* --------------------------------------------------------- */
   /* Allocate workspace that will be needed locally by Sprimme */
   /* --------------------------------------------------------- */
   CHKERR(Num_malloc_iprimme(primme->numEvals, &perm, ctx));

   /*----------------------------------------------------------------------*/
   /* Call the solver                                                      */
   /*----------------------------------------------------------------------*/

   CHKERR(main_iter_Sprimme(
         evals0, perm, evecs, primme->ldevecs, resNorms0, t0, &ret, ctx));

   /*----------------------------------------------------------------------*/
   /* If locking is engaged, the converged Ritz vectors are stored in the  */
   /* order they converged.  They must then be permuted so that they       */
   /* correspond to the sorted Ritz values in evals.                       */
   /*----------------------------------------------------------------------*/

   CHKERR(permute_vecs_Sprimme(&evecs[primme->numOrthoConst * primme->ldevecs],
            primme->nLocal, primme->initSize, primme->ldevecs, perm,
            ctx));

   CHKERR(Num_free_iprimme(perm, ctx));

   /* Copy back evals and resNorms */

   CHKERR(Num_matrix_astype_RHprimme(evals0, 1, primme->numEvals, 1,
         PRIMME_OP_HREAL, (void **)&evals, NULL, PRIMME_OP_REAL,
         0 /* no alloc */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_RHprimme(resNorms0, 1, primme->numEvals, 1,
         PRIMME_OP_HREAL, (void **)&resNorms, NULL, PRIMME_OP_REAL,
         0 /* no alloc */, 1 /* copy */, ctx));

   if (evals != (XREAL *)evals0) CHKERR(Num_free_RHprimme(evals0, ctx));
   if (resNorms != (XREAL *)resNorms0)
      CHKERR(Num_free_RHprimme(resNorms0, ctx));

   primme->stats.elapsedTime = primme_wTimer() - t0;
   return ret;
}


/******************************************************************************
 *
 * static int check_input(double *evals, SCALAR *evecs, double *resNorms, 
 *                        primme_params *primme) 
 *
 * INPUT
 * -----
 *  evals, evecs, resNorms   Output arrays for primme
 *  primme                   the main structure of parameters 
 *
 * return value -   0    If input parameters in primme are appropriate
 *              -4..-32  Inappropriate input parameters were found
 *
 ******************************************************************************/
static int check_input(
      REAL *evals, SCALAR *evecs, REAL *resNorms, primme_params *primme) {
   int ret;
   ret = 0;

   if (primme == NULL)
      ret = -4;
   else if (primme->n < 0 || primme->nLocal < 0 || primme->nLocal > primme->n) 
      ret = -5;
   else if (primme->numProcs < 1)
      ret = -6;
   else if (primme->matrixMatvec == NULL) 
      ret = -7;
   else if (primme->applyPreconditioner == NULL && 
            primme->correctionParams.precondition > 0 ) 
      ret = -8;
   else if (primme->numEvals > primme->n)
      ret = -10;
   else if (primme->numEvals < 0)
      ret = -11;
   else if (fabs(primme->eps) != 0.0L && primme->eps < MACHINE_EPSILON )
      ret = -12;
   else if ( primme->target != primme_smallest  &&
             primme->target != primme_largest  &&
             primme->target != primme_largest_abs  &&
             primme->target != primme_closest_geq  &&
             primme->target != primme_closest_leq  &&
             primme->target != primme_closest_abs    )
      ret = -13;
   else if (primme->numOrthoConst < 0 || primme->numOrthoConst > primme->n)
      ret = -16;
   else if (primme->maxBasisSize < 2 && primme->n > 2) 
      ret = -17;
   else if (primme->minRestartSize < 0 || (primme->minRestartSize == 0
                                    && primme->n > 2 && primme->numEvals > 0))
      ret = -18;
   else if (primme->maxBlockSize < 0
             || (primme->maxBlockSize == 0 && primme->numEvals > 0)) 
      ret = -19;
   else if (primme->restartingParams.maxPrevRetain < 0)
      ret = -20;
   else if (primme->initSize < 0) 
      ret = -22;
   else if (primme->locking == 0 && primme->initSize > primme->maxBasisSize)
      ret = -23;
   else if (primme->locking > 0 && primme->initSize > primme->numEvals)
      ret = -24;
   else if (primme->minRestartSize + primme->restartingParams.maxPrevRetain 
                   >= primme->maxBasisSize && primme->n > primme->maxBasisSize)
      ret = -25;
   else if (primme->minRestartSize > primme->n && primme->n > 2)
      ret = -26;
   else if (primme->printLevel < 0 || primme->printLevel > 5)
      ret = -27; 
   else if (primme->correctionParams.convTest != primme_full_LTolerance &&
            primme->correctionParams.convTest != primme_decreasing_LTolerance &&
            primme->correctionParams.convTest != primme_adaptive_ETolerance &&
            primme->correctionParams.convTest != primme_adaptive )
      ret = -28;
   else if (primme->correctionParams.convTest == primme_decreasing_LTolerance &&
            primme->correctionParams.relTolBase <= 1.0L ) 
      ret = -29;
   else if (evals == NULL)
      ret = -30;
   else if (evecs == NULL)
      ret = -31;
   else if (resNorms == NULL)
      ret = -32;
   else if (primme->locking == 0 && primme->minRestartSize < primme->numEvals &&
            primme->n > 2)
      ret = -33;
   else if (primme->ldevecs < primme->nLocal)
      ret = -34;
   else if (primme->ldOPs != 0 && primme->ldOPs < primme->nLocal)
      ret = -35;
   /* Booked -36 and -37 */
   else if (primme->locking == 0
         && (primme->target == primme_closest_leq
            || primme->target == primme_closest_geq))
      ret = -38;
   else if (primme->massMatrixMatvec &&
            primme->projectionParams.projection != primme_proj_RR)
      ret = -39;
   /* Please keep this if instruction at the end */
   else if ( primme->target == primme_largest_abs ||
             primme->target == primme_closest_geq ||
             primme->target == primme_closest_leq ||
             primme->target == primme_closest_abs   ) {
      if (primme->numTargetShifts <= 0) {
         ret = -14;
      }
      else if (primme->targetShifts == NULL ) {
         ret = -15;
      }
   }

   return ret;
  /***************************************************************************/
} /* end of check_input
   ***************************************************************************/

/*******************************************************************************
 * Subroutine convTestFunAbsolute - This routine implements primme_params.
 *    convTestFun and return an approximate eigenpair converged when           
 *    resNorm < eps*|A| for standard problems, and
 *    resNorm < (|A| + max(|\lambda|)*|B|)*eps for generalized problems
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * evec         The approximate eigenvector
 * eval         The approximate eigenvalue 
 * rNorm        The norm of the residual vector
 * primme       Structure containing various solver parameters
 *
 * OUTPUT PARAMETERS
 * ----------------------------------
 * isConv      if it isn't zero the approximate pair is marked as converged
 ******************************************************************************/

static void convTestFunAbsolute(double *eval, void *evec, double *rNorm,
      int *isConv, primme_params *primme, int *ierr) {

   (void)eval; /* unused parameter */
   (void)evec; /* unused parameter */

   if (primme->massMatrixMatvec == NULL) {
      *isConv = *rNorm < max(primme->eps, MACHINE_EPSILON * 2) *
                               problemNorm_Sprimme(0, primme);
   }
   else {
      *isConv = *rNorm < max(primme->eps, MACHINE_EPSILON) *
                               problemNorm_Sprimme(0, primme);
   }
   *ierr = 0;
}

/*******************************************************************************
 * Subroutine default_monitor - report iterations, #MV, residual norm,
 *    eigenvalues, etc. at every inner/outer iteration and when some pair
 *    converges.       
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * basisEvals   The approximate eigenvalues of the basis
 * basisSize    The size of the basis
 * basisFlags   The state of every approximate pair of the basis (see conv_flags)
 * iblock       Indices of the approximate pairs in the block
 * blockSize    The size of the block
 * basisNorms   The approximate residual norms of the pairs of the basis
 * numConverged The number of pairs converged in the basis and the locked pairs
 *              (this value isn't monotonic!)
 * lockedEvals  The locked eigenvalues
 * numLocked    The number of pairs locked
 * lockedFlags  The state of each locked eigenpair (see conv_flags)
 * lockedNorms  The residual norms of the locked pairs
 * inner_its    The number of performed QMR iterations in the current correction equation
 * LSRes        The residual norm of the linear system at the current QMR iteration
 * event        The event reported
 * primme       Structure containing various solver parameters and statistics
 *
 * OUTPUT
 * ------
 * err          Error code
 * 
 ******************************************************************************/

static void default_monitor(void *basisEvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedEvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, const char *msg, double *time,
      primme_event *event, primme_params *primme, int *err) {

   XREAL *basisEvals = (XREAL*)basisEvals_, *basisNorms = (XREAL*)basisNorms_,
        *lockedEvals = (XREAL*)lockedEvals_, *lockedNorms = (XREAL*)lockedNorms_,
        *LSRes = (XREAL*)LSRes_;

   assert(event != NULL && primme != NULL);

   /* Only print report if this is proc zero */
   if (primme->procID == 0 && primme->outputFile) {
      switch(*event) {
      case primme_event_outer_iteration:
         assert(basisSize && (!*basisSize || (basisEvals && basisFlags)) &&
                blockSize && (!*blockSize || (iblock && basisNorms)) &&
                numConverged);
         assert(!primme->locking ||
                (numLocked && (!*numLocked || (lockedEvals && lockedFlags &&
                                                    lockedNorms))));
         if (primme->printLevel >= 3) {
            int i;  /* Loop variable */
            int found;  /* Reported eigenpairs found */

            if (primme->locking) 
               found = *numLocked;
            else 
               found = *numConverged;

            for (i=0; i < *blockSize; i++) {
               fprintf(primme->outputFile, 
                     "OUT %" PRIMME_INT_P " conv %d blk %d MV %" PRIMME_INT_P " Sec %E EV %13E |r| %.3E\n",
                     primme->stats.numOuterIterations, found, i,
                     primme->stats.numMatvecs, primme->stats.elapsedTime,
                     (double)basisEvals[iblock[i]], (double)basisNorms[iblock[i]]);
            }
         }
         break;
      case primme_event_inner_iteration:
         assert(basisSize && iblock && basisNorms && inner_its && LSRes);
         (void)inner_its;
         if (primme->printLevel >= 4) {
            fprintf(primme->outputFile,
                  "INN MV %" PRIMME_INT_P " Sec %e Eval %e Lin|r| %.3e EV|r| %.3e\n",
                  primme->stats.numMatvecs, primme->stats.elapsedTime,
                  (double)basisEvals[iblock[0]], (double)*LSRes,
                  (double)basisNorms[iblock[0]]);
         }
        break;
      case primme_event_converged:
         assert(numConverged && iblock && basisEvals && basisNorms);
         if ((!primme->locking && primme->printLevel >= 2)
               || (primme->locking && primme->printLevel >= 5))
            fprintf(primme->outputFile, 
                  "#Converged %d eval[ %d ]= %e norm %e Mvecs %" PRIMME_INT_P " Time %g\n",
                  *numConverged, iblock[0], (double)basisEvals[iblock[0]],
                  (double)basisNorms[iblock[0]], primme->stats.numMatvecs,
                  primme->stats.elapsedTime);
         break;
      case primme_event_locked:
         assert(numLocked && lockedEvals && lockedNorms && lockedFlags);
         if (primme->printLevel >= 2) { 
            fprintf(primme->outputFile, 
                  "Lock epair[ %d ]= %e norm %.4e Mvecs %" PRIMME_INT_P " Time %.4e Flag %d\n",
                  *numLocked-1, (double)lockedEvals[*numLocked-1], (double)lockedNorms[*numLocked-1], 
                  primme->stats.numMatvecs, primme->stats.elapsedTime, lockedFlags[*numLocked-1]);
         }
         break;
      case primme_event_message:
         assert(msg != NULL);
         if (primme->printLevel >= 2) { 
            fprintf(primme->outputFile, 
                  "%s\n", msg);
         }
         break;
      case primme_event_profile:
         assert(msg != NULL && time != NULL);
         if (primme->printLevel >= 2) { 
            fprintf(primme->outputFile, 
                  "time for %s : %g\n", msg, *time);
         }
         break;
      default:
         break;
      }
      fflush(primme->outputFile);
   }
   *err = 0;
}

#endif /* SUPPORTED_TYPE */
