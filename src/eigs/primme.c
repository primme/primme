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

#include <stdlib.h>   /* mallocs, free */
#include <stdio.h>    
#include "const.h"
#include "wtime.h"
#include "numerical.h"
#include "main_iter.h"
#include "init.h"
#include "ortho.h"
#include "solve_projection.h"
#include "restart.h"
#include "correction.h"
#include "update_projection.h"
#include "primme_interface.h"

#define ALLOCATE_WORKSPACE_FAILURE -1
#define MALLOC_FAILURE             -2
#define MAIN_ITER_FAILURE          -3

static int allocate_workspace(primme_params *primme, int allocate);
static int check_input(REAL *evals, SCALAR *evecs, REAL *resNorms,
                       primme_params *primme);
static void convTestFunAbsolute(double *eval, void *evec, double *rNorm, int *isConv,
   primme_params *primme, int *ierr);
static void default_monitor(void *basisEvals, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms, int *numConverged,
      void *lockedEvals, int *numLocked, int *lockedFlags, void *lockedNorms,
      int *inner_its, void *LSRes, primme_event *event, primme_params *primme,
      int *err);


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
 
int Sprimme(REAL *evals, SCALAR *evecs, REAL *resNorms, 
            primme_params *primme) {
      
   int ret;
   int *perm;
   double machEps;

   /* ------------------ */
   /* zero out the timer */
   /* ------------------ */
   primme_wTimer(1);

   /* ----------------------- */
   /*  Find machine precision */
   /* ----------------------- */
   machEps = MACHINE_EPSILON;

   /* ------------------ */
   /* Set some defaults  */
   /* ------------------ */
   primme_set_defaults(primme);

   /* -------------------------------------------------------------- */
   /* If needed, we are ready to estimate required memory and return */
   /* -------------------------------------------------------------- */
   if (evals == NULL && evecs == NULL && resNorms == NULL)
       return allocate_workspace(primme, FALSE);

   /* ----------------------------------------------------- */
   /* Reset random number seed if inappropriate for DLARENV */
   /* Yields unique quadruples per proc if procID < 4096^3  */
   /* ----------------------------------------------------- */

   if (primme->iseed[0]<0 || primme->iseed[0]>4095) primme->iseed[0] = 
      primme->procID % 4096;
   if (primme->iseed[1]<0 || primme->iseed[1]>4095) primme->iseed[1] = 
      (int)(primme->procID/4096+1) % 4096;
   if (primme->iseed[2]<0 || primme->iseed[2]>4095) primme->iseed[2] = 
      (int)((primme->procID/4096)/4096+2) % 4096;
   if (primme->iseed[3]<0 || primme->iseed[3]>4095) primme->iseed[3] = 
      (2*(int)(((primme->procID/4096)/4096)/4096)+1) % 4096;

   /* ----------------------- */
   /* Set default convTetFun  */
   /* ----------------------- */

   if (!primme->convTestFun) {
      primme->convTestFun = convTestFunAbsolute;
      if (primme->eps == 0.0) {
         primme->eps = machEps*1e4;
      }
   }

   /* ----------------------- */
   /* Set default monitor     */
   /* ----------------------- */

   if (!primme->monitorFun) {
      primme->monitorFun = default_monitor;
   }

   /* ------------------------------------------------------- */
   /* Check primme input data for bounds, correct values etc. */
   /* ------------------------------------------------------- */

   CHKERRNOABORT(ret = check_input(evals, evecs, resNorms, primme), ret);

   /* ----------------------------------------------------------------------- */
   /* Compute AND allocate memory requirements for main_iter and subordinates */
   /* ----------------------------------------------------------------------- */

   CHKERRNOABORT(allocate_workspace(primme, TRUE), ALLOCATE_WORKSPACE_FAILURE);

   /* --------------------------------------------------------- */
   /* Allocate workspace that will be needed locally by Sprimme */
   /* --------------------------------------------------------- */
   CHKERRNOABORT(MALLOC_PRIMME(primme->numEvals, &perm), MALLOC_FAILURE);

   /*----------------------------------------------------------------------*/
   /* Call the solver                                                      */
   /*----------------------------------------------------------------------*/

   CHKERRNOABORT(main_iter_Sprimme(evals, perm, evecs, primme->ldevecs,
            resNorms, machEps, primme->intWork, primme->realWork, primme),
         MAIN_ITER_FAILURE);

   /*----------------------------------------------------------------------*/
   /* If locking is engaged, the converged Ritz vectors are stored in the  */
   /* order they converged.  They must then be permuted so that they       */
   /* correspond to the sorted Ritz values in evals.                       */
   /*----------------------------------------------------------------------*/

   assert(primme->realWorkSize >= sizeof(SCALAR)*primme->nLocal
         && primme->intWorkSize >= (int)sizeof(int)*primme->initSize);
   permute_vecs_Sprimme(&evecs[primme->numOrthoConst*primme->ldevecs],
         primme->nLocal, primme->initSize, primme->ldevecs, perm,
         (SCALAR*)primme->realWork, (int*)primme->intWork);

   free(perm);

   primme->stats.elapsedTime = primme_wTimer(0);
   return(0);
}


/******************************************************************************
 * Function allocate_workspace - This function computes the amount of integer 
 *    and real workspace needed by the solver and possibly allocates the space 
 *
 * Input: 
 *   allocate  If false, no allocation occurs, but the amounts of int and real 
 *                       workspaces in BYTES are returned in the primme fields 
 *                       primme.intWorkSize, and primme.realWorkSize 
 *             If  true, and if the user-provided space is not sufficient,
 *                       allocation is also performed.
 *
 * Output
 *  primme.intWorkSize   Size of integer space allocated in bytes
 *  primme.realWorkSize  Size of real space allocated in bytes (LONG INT)
 *  *(primme.intWork)    Pointer to the integer space allocated
 *  *(primme.realWork)   Pointer to the real space allocated
 *   
 * 
 * Return value
 * ------------
 * int -  0 if (allocate == true) and the given workspaces are large enough or
 *             have been allocated successfully
 *       -1 if (allocate == true) and memory allocation has failed
 *        1 if (allocate==false) 
 *
 ******************************************************************************/

static int allocate_workspace(primme_params *primme, int allocate) {

   size_t realWorkSize=0;  /* Size of real work space.                  */
   size_t rworkByteSize=0; /* Size of all real data in bytes            */
   int intWorkSize=0;/* Size of integer work space in bytes             */

   int dataSize;     /* Number of SCALAR positions allocated, excluding */
                     /* REAL (see doubleSize below) and work space.  */
   int doubleSize=0; /* Number of doubles allocated exclusively to the  */
                     /* double arrays: hVals, prevRitzVals, blockNorms  */
   int maxEvecsSize; /* Maximum number of vectors in evecs and evecsHat */
   SCALAR *evecsHat=NULL;/* not NULL when evecsHat will be used        */
   SCALAR t;        /* dummy variable */

   maxEvecsSize = primme->numOrthoConst + primme->numEvals;
 
   /* first determine real workspace */

   /*----------------------------------------------------------------------*/
   /* Compute the memory required by the main iteration data structures    */
   /*----------------------------------------------------------------------*/

   dataSize = primme->ldOPs*primme->maxBasisSize   /* Size of V            */
      + primme->ldOPs*primme->maxBasisSize         /* Size of W            */
      + primme->maxBasisSize*primme->maxBasisSize  /* Size of H            */
      + primme->maxBasisSize*primme->maxBasisSize  /* Size of hVecs        */
      + primme->restartingParams.maxPrevRetain*primme->maxBasisSize;
                                                   /* size of prevHVecs    */

   /*----------------------------------------------------------------------*/
   /* Add memory for Harmonic or Refined projection                        */
   /*----------------------------------------------------------------------*/
   if (primme->projectionParams.projection == primme_proj_harmonic ||
         primme->projectionParams.projection == primme_proj_refined) {

      dataSize += primme->ldOPs*primme->maxBasisSize     /* Size of Q      */
         + primme->maxBasisSize*primme->maxBasisSize     /* Size of R      */
         + primme->maxBasisSize*primme->maxBasisSize     /* Size of hU     */
         + primme->maxBasisSize*primme->maxBasisSize;    /* Size of hVecsRot */
      doubleSize += primme->maxBasisSize;                /* Size of hSVals */
   }
   if (primme->projectionParams.projection == primme_proj_harmonic) {
      /* Stored QtV = Q'*V */
      dataSize +=
            primme->maxBasisSize*primme->maxBasisSize;      /* Size of QtV */
   }


   /*----------------------------------------------------------------------*/
   /* Add also memory needed for JD skew projectors                        */
   /*----------------------------------------------------------------------*/
   if ( (primme->correctionParams.precondition && 
         primme->correctionParams.maxInnerIterations != 0 &&
         primme->correctionParams.projectors.RightQ &&
         primme->correctionParams.projectors.SkewQ          ) ) {

      dataSize = dataSize + 
         + primme->nLocal*maxEvecsSize             /* Size of evecsHat     */ 
         + maxEvecsSize*maxEvecsSize               /* Size of M            */
         + maxEvecsSize*maxEvecsSize;              /* Size of UDU          */
      evecsHat = &t; /* set not NULL */
   }

   /*----------------------------------------------------------------------*/
   /* Determine workspace required by init and its children                */
   /*----------------------------------------------------------------------*/

   CHKERR(init_basis_Sprimme(NULL, primme->nLocal, 0, NULL, 0, NULL, 0,
            NULL, 0, NULL, 0, NULL, 0, NULL, 0, NULL, &realWorkSize,
            &primme->maxBasisSize, NULL, NULL, primme), -1);

   /*----------------------------------------------------------------------*/
   /* Determine orthogalization workspace with and without locking.        */
   /*----------------------------------------------------------------------*/

   CHKERR(ortho_Sprimme(NULL, 0, NULL, 0, primme->maxBasisSize,
            primme->maxBasisSize+primme->maxBlockSize-1, NULL, primme->nLocal, 
            primme->locking?maxEvecsSize:primme->numOrthoConst+1, primme->nLocal,
            NULL, 0.0, NULL, &realWorkSize, primme), -1);

   /*----------------------------------------------------------------------*/
   /* Determine workspace required by solve_H and its children             */
   /*----------------------------------------------------------------------*/

   CHKERR(solve_H_Sprimme(NULL, primme->maxBasisSize, 0, NULL, 0, NULL, 0,
            NULL, 0, NULL, 0, NULL, NULL, 0, 0.0, &realWorkSize, NULL, 0,
            &intWorkSize, primme), -1);

   /*----------------------------------------------------------------------*/
   /* Determine workspace required by solve_correction and its children    */
   /*----------------------------------------------------------------------*/

   CHKERR(solve_correction_Sprimme(NULL, 0, NULL, 0, NULL, 0, NULL, 0, NULL, 
            NULL, NULL, maxEvecsSize, 0, NULL, NULL, NULL, NULL, 
            primme->maxBasisSize, NULL, NULL, primme->maxBlockSize, NULL,
            0.0, NULL, &realWorkSize, &intWorkSize, 0, primme), -1);

   /*----------------------------------------------------------------------*/
   /* Determine workspace required by restarting and its children          */
   /*----------------------------------------------------------------------*/

   CHKERR(restart_Sprimme(NULL, NULL, primme->nLocal, primme->maxBasisSize,
            0, NULL, NULL, NULL, NULL, &primme->maxBlockSize, NULL, NULL, 0,
            NULL, NULL, NULL, evecsHat, 0, NULL, 0, NULL, 0, NULL,
            &primme->numEvals, &primme->numEvals, &primme->numEvals, NULL, NULL,
            &primme->restartingParams.maxPrevRetain, primme->maxBasisSize,
            primme->initSize, NULL, &primme->maxBasisSize, NULL,
            primme->maxBasisSize, NULL, 0, NULL, 0, NULL, 0, NULL, 0, 0, NULL,
            0, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, 0.0, NULL,
            &realWorkSize, &intWorkSize, 0, primme), -1);

   /*----------------------------------------------------------------------*/
   /* Determine workspace required by main_iter and its children           */
   /*----------------------------------------------------------------------*/

   CHKERR(update_projection_Sprimme(NULL, 0, NULL, 0, NULL, 0, 0, 0,
            primme->maxBasisSize, NULL, &realWorkSize, 0, primme), -1);

   CHKERR(prepare_candidates_Sprimme(NULL, 0, NULL, 0, primme->nLocal, NULL, 0,
            primme->maxBasisSize, NULL, NULL, NULL, 0, NULL, NULL, NULL,
            primme->numEvals, NULL, 0, primme->maxBlockSize,
            NULL, primme->numEvals, 0, NULL, NULL, 0, 0.0, NULL,
            &primme->maxBlockSize, NULL, NULL, NULL, NULL, 0, 0, NULL, NULL,
            NULL, &realWorkSize, &intWorkSize, 0, primme), -1);

   CHKERR(retain_previous_coefficients_Sprimme(NULL, 0, NULL, 0, NULL, 0,
            0, 0, NULL, primme->maxBlockSize, NULL,
            &primme->restartingParams.maxPrevRetain, &intWorkSize, 0, primme),
            -1);
 
   /*----------------------------------------------------------------------*/
   /* Workspace needed by function verify_norms                            */
   /*----------------------------------------------------------------------*/
   realWorkSize = max(realWorkSize, (size_t)2*primme->numEvals);

   /*----------------------------------------------------------------------*/
   /* The following size is always allocated as REAL                       */
   /*----------------------------------------------------------------------*/

   doubleSize += 4     /* padding cause by TO_REAL aligning them to SCALAR */
      + primme->maxBasisSize                       /* Size of hVals        */
      + primme->numEvals+primme->maxBasisSize      /* Size of prevRitzVals */
      + primme->maxBlockSize                       /* Size of blockNorms   */
      + primme->maxBasisSize;                      /* Size of basisNorms   */

   /*----------------------------------------------------------------------*/
   /* Determine the integer workspace needed                               */
   /*----------------------------------------------------------------------*/

   intWorkSize += primme->maxBasisSize /* Size of flag               */
      + 2*primme->maxBlockSize         /* Size of iev and ilev       */
      + maxEvecsSize;                  /* Size of ipivot             */
   if (primme->locking) {
      intWorkSize += primme->numEvals; /* Size of lockedFlags        */
   }

   /*----------------------------------------------------------------------*/
   /* byte sizes:                                                          */
   /*----------------------------------------------------------------------*/
   
   rworkByteSize = (dataSize + realWorkSize)*sizeof(SCALAR)
                                + doubleSize*sizeof(REAL); 

   /*----------------------------------------------------------------------*/
   /* If only the amount of required workspace is needed return it in bytes*/
   /*----------------------------------------------------------------------*/

   if (!allocate) {
      primme->intWorkSize  = intWorkSize*sizeof(int);
      primme->realWorkSize = rworkByteSize;
      return 1;
   }

   /*----------------------------------------------------------------------*/
   /* Allocate the required workspace, if the user did not provide enough  */
   /*----------------------------------------------------------------------*/
   if (primme->realWork != NULL && primme->realWorkSize < rworkByteSize) {
      return -35;
   }
   else if (primme->realWork == NULL) {
      primme->realWorkSize = rworkByteSize;
      if (primme->printLevel >= 5) fprintf(primme->outputFile, 
         "Allocating real workspace: %zd bytes\n", primme->realWorkSize);
      CHKERRM(MALLOC_PRIMME(rworkByteSize, (char**)&primme->realWork), MALLOC_FAILURE,
            "Failed to allocate %zd bytes\n", rworkByteSize);
   }

   if (primme->intWork != NULL
         && primme->intWorkSize < intWorkSize*(int)sizeof(int)) {
      return -36;
   }
   else if (primme->intWork == NULL) {
      primme->intWorkSize = intWorkSize*sizeof(int);
      if (primme->printLevel >= 5) fprintf(primme->outputFile, 
         "Allocating integer workspace: %d bytes\n", primme->intWorkSize);
      CHKERRM(MALLOC_PRIMME(intWorkSize, &primme->intWork), MALLOC_FAILURE,
            "Failed to allocate %d bytes\n", primme->intWorkSize);
   }

   return 0;

  /***************************************************************************/
} /* end of allocate workspace
  ****************************************************************************/

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
static int check_input(REAL *evals, SCALAR *evecs, REAL *resNorms, 
                       primme_params *primme) {
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
   /* ret = -9 is free */
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
   else if (primme->maxBasisSize < 2 && primme->maxBasisSize != primme->n) 
      ret = -17;
   else if (primme->minRestartSize < 0 || (primme->minRestartSize == 0
                                    && primme->n > 2 && primme->numEvals > 0))
      ret = -18;
   else if (primme->maxBlockSize < 0
             || (primme->maxBlockSize == 0 && primme->numEvals > 0)) 
      ret = -19;
   else if (primme->restartingParams.maxPrevRetain < 0)
      ret = -20;
   else if (primme->restartingParams.scheme != primme_thick &&
            primme->restartingParams.scheme != primme_dtr)
      ret = -21;
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
 *    resNorm < eps*(aNorm != 0 ? aNorm : aNormEstimate) or
 *    resNorm is close to machineEpsilon * aNorm.          
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

   const double machEps = MACHINE_EPSILON;
   const double aNorm = (primme->aNorm > 0.0) ?
      primme->aNorm : primme->stats.estimateLargestSVal;
   (void)eval; /* unused parameter */
   (void)evec; /* unused parameter */
   *isConv = *rNorm < max(
               primme->eps * aNorm,
               machEps * 3.16 * primme->stats.estimateLargestSVal);
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
      int *inner_its, void *LSRes_, primme_event *event, primme_params *primme,
      int *err)
{
   REAL *basisEvals = (REAL*)basisEvals_, *basisNorms = (REAL*)basisNorms_,
        *lockedEvals = (REAL*)lockedEvals_, *lockedNorms = (REAL*)lockedNorms_,
        *LSRes = (REAL*)LSRes_;
   assert(event != NULL && primme != NULL);

   /* Only print report if this is proc zero */
   if (primme->procID == 0 && primme->outputFile) {
      switch(*event) {
      case primme_event_outer_iteration:
         assert(basisEvals && basisSize && basisFlags && iblock && blockSize
                && basisNorms && numConverged);
         assert(!primme->locking || (lockedEvals && numLocked && lockedFlags
                 && lockedNorms));
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
                     primme->stats.numMatvecs, primme_wTimer(0),
                     basisEvals[iblock[i]], (double)basisNorms[iblock[i]]);
            }
         }
         break;
      case primme_event_inner_iteration:
         assert(basisSize && iblock && basisNorms && inner_its && LSRes);
         (void)inner_its;
         if (primme->printLevel >= 4) {
            fprintf(primme->outputFile,
                  "INN MV %" PRIMME_INT_P " Sec %e Eval %e Lin|r| %.3e EV|r| %.3e\n",
                  primme->stats.numMatvecs, primme_wTimer(0),
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
                  *numConverged, iblock[0], basisEvals[iblock[0]],
                  basisNorms[iblock[0]], primme->stats.numMatvecs,
                  primme_wTimer(0));
         break;
      case primme_event_locked:
         assert(numLocked && lockedEvals && lockedNorms && lockedFlags);
         if (primme->printLevel >= 2) { 
            fprintf(primme->outputFile, 
                  "Lock epair[ %d ]= %e norm %.4e Mvecs %" PRIMME_INT_P " Time %.4e Flag %d\n",
                  *numLocked-1, lockedEvals[*numLocked-1], lockedNorms[*numLocked-1], 
                  primme->stats.numMatvecs, primme_wTimer(0), lockedFlags[*numLocked-1]);
         }
         break;
      default:
         break;
      }
      fflush(primme->outputFile);
   }
   *err = 0;
}
