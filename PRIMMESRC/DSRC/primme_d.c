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
 * File: primme.c
 *
 * Purpose - Real, double precision front end to the multimethod eigensolver
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
#include "primme.h"
#include "const.h"
#include "wtime.h"
#include "main_iter_d.h"
#include "ortho_d.h"
#include "solve_H_d.h"
#include "correction_d.h"
#include "primme_private_d.h"
#include "numerical_d.h"

/*******************************************************************************
 * Subroutine dprimme - This routine is a front end used to perform 
 *    error checking on the input parameters, perform validation, 
 *    and make the call to main_iter. 
 *
 *    Calling dprimme with all evals, evecs, resNorms set to NULL
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
 
int dprimme(double *evals, double *evecs, double *resNorms, 
            primme_params *primme) {
      
   int ret;
   int *perm;
   double machEps;

   /* ------------------ */
   /* zero out the timer */
   /* ------------------ */
   primme_wTimer(1);

   /* ---------------------------- */
   /* Clear previous error reports */
   /* ---------------------------- */
   primme_DeleteStackTrace(primme);

   /* ----------------------- */
   /*  Find machine precision */
   /* ----------------------- */
   machEps = Num_dlamch_primme("E");

   /* ----------------------------------------- */
   /* Set some defaults for sequential programs */
   /* ----------------------------------------- */
   if (primme->numProcs == 1) {
      primme->nLocal = primme->n;
      primme->procID = 0;
      if (primme->globalSumDouble == NULL) 
         primme->globalSumDouble = primme_seq_globalSumDouble;
   }

   /* --------------------------------------------------------------------- */
   /* Decide on whether to use locking (hard locking), or not (soft locking)*/
   /* --------------------------------------------------------------------- */
   if (primme->target != primme_smallest &&
       primme->target != primme_largest ) {
      /* Locking is necessary as interior Ritz values can cross shifts */
      primme->locking = 1;
   }
   else {
      if (primme->locking == 0) {
         /* use locking when not enough vectors to restart with */
         primme->locking = (primme->numEvals > primme->minRestartSize);   
      }
   }

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

   /* ------------------------------------------------------- */
   /* Check primme input data for bounds, correct values etc. */
   /* ------------------------------------------------------- */

   ret = check_input(evals, evecs, resNorms, primme);

   if (ret != 0) {
      primme_PushErrorMessage(Primme_dprimme, Primme_check_input, ret,
                      __FILE__, __LINE__, primme);
      primme->stats.elapsedTime = primme_wTimer(0);
      return ret;
   }
   
   /* ----------------------------------------------------------------------- */
   /* Compute AND allocate memory requirements for main_iter and subordinates */
   /* ----------------------------------------------------------------------- */

   ret = allocate_workspace(primme, TRUE);

   if (ret != 0) {
      primme_PushErrorMessage(Primme_dprimme, Primme_allocate_workspace, ret,
                      __FILE__, __LINE__, primme);
      primme->stats.elapsedTime = primme_wTimer(0);
      return ALLOCATE_WORKSPACE_FAILURE;
   }

   /* --------------------------------------------------------- */
   /* Allocate workspace that will be needed locally by dprimme */
   /* --------------------------------------------------------- */
   perm = (int *)primme_calloc((primme->numEvals), sizeof(int), "Perm array");

   if (perm == NULL) {
      primme_PushErrorMessage(Primme_dprimme, Primme_malloc, 0, 
                      __FILE__, __LINE__, primme);
      primme->stats.elapsedTime = primme_wTimer(0);
      return MALLOC_FAILURE;
   }

   /*----------------------------------------------------------------------*/
   /* Call the solver                                                      */
   /*----------------------------------------------------------------------*/

   ret = main_iter_dprimme(evals, perm, evecs, resNorms, machEps, 
                   primme->intWork, primme->realWork, primme);

   if (ret < 0) {
      primme_PushErrorMessage(Primme_dprimme, Primme_main_iter, 
                      ret, __FILE__, __LINE__, primme);
      primme->stats.elapsedTime = primme_wTimer(0);
      return MAIN_ITER_FAILURE;
   }
   /*----------------------------------------------------------------------*/

   /*----------------------------------------------------------------------*/
   /* If locking is engaged, the converged Ritz vectors are stored in the  */
   /* order they converged.  They must then be permuted so that they       */
   /* correspond to the sorted Ritz values in evals.                       */
   /*----------------------------------------------------------------------*/

   permute_evecs_dprimme(&evecs[primme->numOrthoConst], perm, 
        (double *) primme->realWork, primme->numEvals, primme->nLocal);

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

   long int realWorkSize;  /* Size of real work space.                  */
   long int rworkByteSize; /* Size of all real data in bytes            */

   int dataSize;     /* Number of double positions allocated, excluding */
                     /* doubles (see doubleSize below) and work space.  */
   int doubleSize;   /* Number of doubles allocated exclusively to the  */
                     /* double arrays: hVals, prevRitzVals, blockNorms  */
   int maxEvecsSize; /* Maximum number of vectors in evecs and evecsHat */
   int intWorkSize;  /* Size of integer work space in bytes             */
   int orthoSize;    /* Amount of work space required by ortho routine  */
   int solveCorSize; /* work space for solve_correction and inner_solve */

   maxEvecsSize = primme->numOrthoConst + primme->numEvals;

   /* first determine real workspace */

   /*----------------------------------------------------------------------*/
   /* Compute the memory required by the main iteration data structures    */
   /*----------------------------------------------------------------------*/

   dataSize = primme->nLocal*primme->maxBasisSize  /* Size of V            */
      + primme->nLocal*primme->maxBasisSize        /* Size of W            */
      + primme->maxBasisSize*primme->maxBasisSize  /* Size of H            */
      + primme->maxBasisSize*primme->maxBasisSize  /* Size of hVecs        */
      + primme->restartingParams.maxPrevRetain*primme->maxBasisSize;
                                                   /* size of prevHVecs    */

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
   }

   /*----------------------------------------------------------------------*/
   /* Determine orthogalization workspace with and without locking.        */
   /*----------------------------------------------------------------------*/

   if (primme->locking) {
      orthoSize = ortho_dprimme(NULL, primme->nLocal, primme->maxBasisSize,
         primme->maxBasisSize+primme->maxBlockSize-1, NULL, primme->nLocal, 
         maxEvecsSize, primme->nLocal, NULL, 0.0, NULL, 0, primme);
   }
   else {
      orthoSize = ortho_dprimme(NULL, primme->nLocal, primme->maxBasisSize,
         primme->maxBasisSize+primme->maxBlockSize-1, NULL, primme->nLocal, 
         primme->numOrthoConst+1, primme->nLocal, NULL, 0.0, NULL, 0, primme);
   }

   /*----------------------------------------------------------------------*/
   /* Determine workspace required by solve_correction and its children    */
   /*----------------------------------------------------------------------*/

   solveCorSize = solve_correction_dprimme(NULL, NULL, NULL, NULL, NULL, 
                  NULL, NULL, maxEvecsSize, 0, NULL, NULL, NULL, NULL, 
                  primme->maxBasisSize, NULL, NULL, primme->maxBlockSize, 
                  1.0, 0.0, 1.0, NULL, NULL, 0, primme);

   /*----------------------------------------------------------------------*/
   /* Workspace is reused in many functions. Allocate the max needed by any*/
   /*----------------------------------------------------------------------*/
   realWorkSize = Num_imax_primme(8,

      /* Workspace needed by init_basis */
      Num_imax_primme(3, 
         maxEvecsSize*primme->numOrthoConst, maxEvecsSize, orthoSize),

      /* Workspace needed by solve_correction and its child inner_solve */
      solveCorSize, 

      /* Workspace needed by function solve_H */
#ifdef ESSL
      2*primme->maxBasisSize +
         primme->maxBasisSize*(primme->maxBasisSize + 1)/2,
#else
      3*primme->maxBasisSize,
#endif
   
      /* Workspace needed by function check_convergence */ 
      max(primme->maxBasisSize*primme->maxBlockSize + primme->maxBlockSize,
                2*maxEvecsSize*primme->maxBlockSize),

      /* Workspace needed by function restart*/
      primme->restartingParams.maxPrevRetain*
      primme->restartingParams.maxPrevRetain  /* for submatrix of prev hvecs */
      + Num_imax_primme(4, primme->maxBasisSize, 
           3*primme->restartingParams.maxPrevRetain,
           primme->maxBasisSize*primme->restartingParams.maxPrevRetain,
           primme->maxBasisSize*primme->maxBasisSize,     /* for DTR copying */
           maxEvecsSize*primme->numEvals), /*this one is for UDU w/o locking */

      /* Workspace needed by function verify_norms */
      2*primme->numEvals,

      /* space needed by lock vectors (no need w/o lock but doesn't add any) */
      (2*primme->maxBasisSize) + Num_imax_primme(3, 
          maxEvecsSize*primme->maxBasisSize, orthoSize, 3*primme->maxBasisSize),

      /* maximum workspace needed by ortho */ 
      orthoSize);

   /*----------------------------------------------------------------------*/
   /* The following size is always alloced as double                       */
   /*----------------------------------------------------------------------*/

   doubleSize = primme->maxBasisSize               /* Size of hVals        */
      + primme->numEvals+primme->maxBasisSize      /* Size of prevRitzVals */
      + primme->maxBlockSize;                      /* Size of blockNorms   */

   /*----------------------------------------------------------------------*/
   /* Determine the integer workspace needed                               */
   /*----------------------------------------------------------------------*/

   intWorkSize = primme->maxBasisSize /* Size of flag               */
      + 2*primme->maxBlockSize        /* Size of iev and ilev       */
      + maxEvecsSize                  /* Size of ipivot             */
      + 2*primme->maxBasisSize;       /* Size of 2 perms in solve_H */

   /*----------------------------------------------------------------------*/
   /* byte sizes:                                                          */
   /*----------------------------------------------------------------------*/
   
   rworkByteSize = (dataSize + realWorkSize)*sizeof(double)
                                + doubleSize*sizeof(double); 

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
   if (primme->realWorkSize < rworkByteSize || primme->realWork == NULL) {
      if (primme->realWork != NULL) {
         free(primme->realWork);
      }
      primme->realWorkSize = rworkByteSize;
      primme->realWork = (void *) primme_valloc(rworkByteSize,"Real Alloc");
      if (primme->printLevel >= 5) fprintf(primme->outputFile, 
         "Allocating real workspace: %ld bytes\n", primme->realWorkSize);
   }

   if (primme->intWorkSize < intWorkSize*(int)sizeof(int) || primme->intWork==NULL) {
      if (primme->intWork != NULL) {
         free(primme->intWork);
      }
      primme->intWorkSize = intWorkSize*sizeof(int);
      primme->intWork= (int *)primme_valloc(primme->intWorkSize ,"Int Alloc");
      if (primme->printLevel >= 5) fprintf(primme->outputFile, 
         "Allocating integer workspace: %d bytes\n", primme->intWorkSize);
   }

   if (primme->intWork == NULL || primme->realWork == NULL) {

      primme_PushErrorMessage(Primme_allocate_workspace, Primme_malloc, 0, 
         __FILE__, __LINE__, primme);
      return MALLOC_FAILURE;
   }
      
   return 0;

  /***************************************************************************/
} /* end of allocate workspace
  ****************************************************************************/

/******************************************************************************
 *
 * static int check_input(double *evals, double *evecs, double *resNorms, 
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
static int check_input(double *evals, double *evecs, double *resNorms, 
                       primme_params *primme) {
   int ret;
   ret = 0;

   if (primme == NULL)
      ret = -4;
   else if (primme->n <= 0 || primme->nLocal <= 0) 
      ret = -5;
   else if (primme->numProcs < 1)
      ret = -6;
   else if (primme->matrixMatvec == NULL) 
      ret = -7;
   else if (primme->applyPreconditioner == NULL && 
            primme->correctionParams.precondition ) 
      ret = -8;
   else if (primme->globalSumDouble == NULL)
      ret = -9;
   else if (primme->numEvals > primme->n)
      ret = -10;
   else if (primme->numEvals < 0)
      ret = -11;
   else if (primme->eps > 0.0L && primme->eps < Num_dlamch_primme("E") )
      ret = -12;
   else if ( primme->target != primme_smallest  &&
             primme->target != primme_largest  &&
             primme->target != primme_closest_geq  &&
             primme->target != primme_closest_leq  &&
             primme->target != primme_closest_abs    )
      ret = -13;
   else if ( primme->target == primme_closest_geq ||
             primme->target == primme_closest_leq ||
             primme->target == primme_closest_abs   ) {
      if (primme->numTargetShifts <= 0) {
         ret = -14;
      }
      else if (primme->targetShifts == NULL ) {
         ret = -15;
      }
   }
   else if (primme->numOrthoConst < 0 || primme->numOrthoConst >=primme->n)
      ret = -16;
   else if (primme->maxBasisSize < 2) 
      ret = -17;
   else if (primme->minRestartSize <= 0) 
      ret = -18;
   else if (primme->maxBlockSize <= 0) 
      ret = -19;
   else if (primme->restartingParams.maxPrevRetain < 0)
      ret = -20;
   else if (primme->restartingParams.scheme != primme_thick &&
            primme->restartingParams.scheme != primme_dtr)
      ret = -21;
   else if (primme->initSize < 0) 
      ret = -22;
   else if (!primme->locking && primme->initSize > primme->maxBasisSize)
      ret = -23;
   else if (primme->locking && primme->initSize > primme->numEvals)
      ret = -24;
   else if (primme->minRestartSize + primme->restartingParams.maxPrevRetain 
                   >= primme->maxBasisSize)
      ret = -25;
   else if (primme->minRestartSize >= primme->n)
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

   return ret;
  /***************************************************************************/
} /* end of check_input
   ***************************************************************************/
