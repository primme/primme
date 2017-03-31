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
 * File: correction.c
 *
 * Purpose - Computes the correction corresponding to each block Ritz vectors.
 *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "const.h"
#include "numerical.h"
#include "correction.h"
#include "inner_solve.h"
#include "globalsum.h"
#include "auxiliary_eigs.h"

static REAL computeRobustShift(int blockIndex, double resNorm, 
   REAL *prevRitzVals, int numPrevRitzVals, REAL *sortedRitzVals, 
   REAL *approxOlsenShift, int numSorted, int *ilev, primme_params *primme);

static void mergeSort(REAL *lockedEvals, int numLocked, REAL *ritzVals, 
   int *flags, int basisSize, REAL *sortedRitzVals, int *ilev, int blockSize,
   primme_params *primme);
 
static int Olsen_preconditioner_block(SCALAR *r, PRIMME_INT ldr, SCALAR *x,
      PRIMME_INT ldx, int blockSize, SCALAR *rwork, primme_params *primme);

static int setup_JD_projectors(SCALAR *x, SCALAR *evecs, PRIMME_INT ldevecs,
      SCALAR *evecsHat, PRIMME_INT ldevecsHat, SCALAR *Kinvx, SCALAR *xKinvx, 
      SCALAR **Lprojector, PRIMME_INT *ldLprojector, SCALAR **RprojectorQ,
      PRIMME_INT *ldRprojectorQ, SCALAR **RprojectorX,
      PRIMME_INT *ldRprojectorX,  int *sizeLprojector, int *sizeRprojectorQ,
      int *sizeRprojectorX, int numLocked, int numConverged,
      primme_params *primme);


/*******************************************************************************
 * Subroutine solve_correction - This routine solves the correction equation
 *    for each Ritz vector and residual in the block.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * W              The last blockSize vectors of W contain the residuals
 *
 * evecs          The converged Ritz vectors.  Array is of dimension numLocked.
 *
 * evecsHat       K^{-1}evecs given a preconditioner K. 
 *                (accessed only if skew projector is requested)
 *
 * UDU            The factorization of the matrix evecs'*evecsHat
 * ipivot         The pivots for the UDU factorization
 *                (UDU, ipivot accessed only if skew projector is requested)
 *
 * lockedEvals    The locked eigenvalues
 *
 * numLocked      The number of locked eigenvalues (zero if not locking)
 *
 * numConvergedStored  Number of converged eigenvectors in V that have been
 *                copied in evecs when no locking is employed, to accommodate 
 *                for a requested skew projection with evecs. 
 *                   
 * ritzVals       Array of size basisSize. The Ritz values corresponding to 
 *                the current basis V
 *
 * flags          Array indicating which eigenvectors have converged
 *
 * basisSize      The size of the basis V
 *
 * iev            Array of size blockSize.  The index of each block Ritz value
 *                in ritzVals (and each Ritz vector in V)
 *
 * blockSize      The current block size
 *
 * machEps        machine precision 
 *
 * rwork          Real workspace of size          
 *                3*maxEvecsSize + 2*primme->maxBlockSize 
 *                + (primme->numEvals+primme->maxBasisSize)
 *                        *----------------------------------------------------*
 *                        | The following are optional and mutually exclusive: |
 *                        *------------------------------+                     |
 *                + 4*primme->nLocal + primme->nLocal    | For QMR work and sol|
 *                + primme->nLocal*primme->maxBlockSize  | OLSEN for Kinvx     |
 *                                                       *---------------------*
 *
 * rworkSize      the size of rwork. If less than needed, func returns needed.
 *
 * iwork          Integer workspace of size maxBlockSize 
 * 
 * primme         Structure containing various solver parameters 
 * 
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * V              The orthonormal basis.  The last blockSize vectors of V
 *                contain the Ritz vectors
 *
 * blockNorms     On input, residual norms of the Ritz vectors computed 
 *                during the current outer iteration.
 *                On output, the approximate residual norms of the corrected
 *                Ritz vectors (through JDQMR only).
 *
 *
 * prevRitzVals   Array of size numPrevRitzVals.  The Ritz values from the
 *                previous iteration and converged Ritz values
 *
 * numPrevRitzVals  The of size prevRitzVals updated every outer step
 *
 * touch            Parameter used in inner solve stopping criteria
 *
 * Return Value
 * ------------
 * int  Error code: 0 upon success, nonzero otherwise
 *                 -1 innner solver failure
 *                 >0 the needed rworkSize, if the given was not enough
 *
 ******************************************************************************/
 
TEMPLATE_PLEASE
int solve_correction_Sprimme(SCALAR *V, PRIMME_INT ldV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *evecs, PRIMME_INT ldevecs, SCALAR *evecsHat,
      PRIMME_INT ldevecsHat, SCALAR *UDU, int *ipivot, REAL *lockedEvals, 
      int numLocked, int numConvergedStored, REAL *ritzVals, 
      REAL *prevRitzVals, int *numPrevRitzVals, int *flags, int basisSize, 
      REAL *blockNorms, int *iev, int blockSize, int *touch, double machEps,
      SCALAR *rwork, size_t *rworkSize, int *iwork, int iworkSize,
      primme_params *primme) {

   int blockIndex;         /* Loop index.  Ranges from 0..blockSize-1.       */
   int ritzIndex;          /* Ritz value index blockIndex corresponds to.    */
                           /* Possible values range from 0..basisSize-1.     */
   int sortedIndex;        /* Ritz value index in sortedRitzVals, blockIndex */
                           /* corresponds to. Range 0..numLocked+basisSize-1 */
   size_t neededRsize;     /* Needed size for rwork. If not enough return    */
   size_t linSolverRWorkSize;/* Size of the linSolverRWork array.            */
   int *ilev;              /* Array of size blockSize.  Maps the target Ritz */
                           /* values to their positions in the sortedEvals   */
                           /* array.                                         */
   int sizeLprojector;     /* Sizes of the various left/right projectors     */
   int sizeRprojectorQ;    /* These will be 0/1/or numOrthConstr+numLocked   */
   int sizeRprojectorX;    /* or numOrthConstr+numConvergedStored w/o locking*/

   SCALAR *r, *x, *sol;  /* Residual, Ritz vector, and correction.         */
   SCALAR *linSolverRWork;/* Workspace needed by linear solver.            */
   REAL *sortedRitzVals; /* Sorted array of current and converged Ritz     */
                           /* values.  Size of array is numLocked+basisSize. */
   double *blockOfShifts;  /* Shifts for (A-shiftI) or (if needed) (K-shiftI)*/
   REAL *approxOlsenEps; /* Shifts for approximate Olsen implementation    */
   SCALAR *Kinvx;         /* Workspace to store K^{-1}x                     */
   SCALAR *Lprojector;   /* Q pointer for (I-Q*Q'). Usually points to evecs*/
   SCALAR *RprojectorQ;  /* May point to evecs/evecsHat depending on skewQ */
   SCALAR *RprojectorX;  /* May point to x/Kinvx depending on skewX        */
   PRIMME_INT ldLprojector;  /* The leading dimension of Lprojector     */
   PRIMME_INT ldRprojectorQ; /* The leading dimension of RprojectorQ    */
   PRIMME_INT ldRprojectorX; /* The leading dimension of RprojectorL    */


   SCALAR xKinvx;                        /* Stores x'*K^{-1}x if needed    */
   REAL eval, shift, robustShift;       /* robust shift values.           */

   /*------------------------------------------------------------*/
   /* Subdivide the workspace with pointers, and figure out      */
   /* the total amount of needed real workspace (neededRsize)    */
   /*------------------------------------------------------------*/

   /* needed worksize */
   neededRsize = 0;
   Kinvx       = rwork;
   /* Kinvx will have nonzero size if precond and both RightX and SkewX */
   if (primme->correctionParams.projectors.RightX &&  
       primme->correctionParams.projectors.SkewX ) { 

      /* OLSEN's method requires a block, but JDQMR is vector by vector */
      if (primme->correctionParams.maxInnerIterations == 0) {    
         sol = Kinvx + primme->ldOPs*blockSize;
         neededRsize = neededRsize + primme->ldOPs*blockSize;
      }
      else {
         sol = Kinvx + primme->nLocal;
         neededRsize = neededRsize + primme->nLocal;
      }
   }
   else {
      sol = Kinvx + 0;
   }
   if (primme->correctionParams.maxInnerIterations == 0) {    
      linSolverRWork = sol + 0;                   /* sol not needed for GD */
      linSolverRWorkSize = 0;                     /* No inner solver used  */
   }
   else {
      linSolverRWork = sol + primme->nLocal;      /* sol needed in innerJD */
      neededRsize = neededRsize + primme->nLocal;
      linSolverRWorkSize =                        /* Inner solver worksize */
              4*primme->nLocal + 2*(primme->numOrthoConst+primme->numEvals);
      neededRsize = neededRsize + linSolverRWorkSize;
   }
   sortedRitzVals = (REAL *)(linSolverRWork + linSolverRWorkSize);
   blockOfShifts  = ALIGN(sortedRitzVals + (numLocked+basisSize), double);
   approxOlsenEps = ALIGN(blockOfShifts  + blockSize, REAL);
   neededRsize = neededRsize + numLocked+basisSize
      + blockSize*(1+sizeof(double)/sizeof(REAL)) + 2;

   /* Return memory requirements */
   if (V == NULL) {
      *rworkSize = max(*rworkSize, neededRsize);
      *iwork = max(*iwork, primme->numEvals+primme->maxBasisSize);
      return 0;
   }
   assert(neededRsize <= *rworkSize);

   /* Subdivide also the integer work space */
   ilev = iwork;       /* of size blockSize */

   /*------------------------------------------------------------*/
   /*  Figuring out preconditioning shifts  (robust, Olsen, etc) */
   /*------------------------------------------------------------*/
   /* blockOfShifts will contain the preconditioning shifts:               */
   /* either Ritz values or robustShifts computed below. These shifts      */
   /* will be used in the correction equations or in inverting (K-sigma I) */
   /* approxOlsenEps will contain error approximations for eigenavalues    */
   /* to be used for Olsen's method (when innerIterations =0).             */
    
   if (primme->locking && 
      (primme->target == primme_smallest || primme->target == primme_largest)) {
      /* Combine the sorted list of locked Ritz values with the sorted  */
      /* list of current Ritz values, ritzVals.  The merging of the two */
      /* lists lockedEvals and ritzVals is stored in sortedRitzVals.    */

      assert(iworkSize >= numLocked+blockSize);
      mergeSort(lockedEvals, numLocked, ritzVals, flags, basisSize, 
                   sortedRitzVals, ilev, blockSize, primme);
   }
   else {
      /* In the case of soft-locking or when we look for interior ones  */
      /* the sorted evals are simply the ritzVals, targeted as iev      */

      sortedRitzVals = ritzVals;
      ilev = iev;
   }

   /*-----------------------------------------------------------------*/
   /* For interior pairs use not the robust, but user provided shifts */
   /*-----------------------------------------------------------------*/

   if (primme->target != primme_smallest && primme->target != primme_largest) {

      for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {

         /* Considering |Ritz value - exact eigenvalue| <= residual norm, then      */
         /* we take the closest point in the interval Ritz value +- residual norm   */
         /* to the user shift as the proper shift.                                  */

         sortedIndex = ilev[blockIndex];
         if (sortedRitzVals[sortedIndex] - blockNorms[blockIndex]
               <  primme->targetShifts[min(primme->numTargetShifts-1,numLocked)]
             &&   primme->targetShifts[min(primme->numTargetShifts-1,numLocked)]
               < sortedRitzVals[sortedIndex] + blockNorms[blockIndex])
            blockOfShifts[blockIndex] = 
               primme->targetShifts[min(primme->numTargetShifts-1, numLocked)];
         else
            blockOfShifts[blockIndex] = sortedRitzVals[sortedIndex]
               + blockNorms[blockIndex] * (
                  primme->targetShifts[min(primme->numTargetShifts-1,numLocked)]
                     < sortedRitzVals[sortedIndex] ? -1 : 1);
         
         if (sortedIndex < *numPrevRitzVals) {
            approxOlsenEps[blockIndex] = 
            fabs(prevRitzVals[sortedIndex] - sortedRitzVals[sortedIndex]);
         }  
         else {
            approxOlsenEps[blockIndex] = blockNorms[blockIndex];
         }  
      } /* for loop */

      /* Remember the previous ritz values*/
      *numPrevRitzVals = basisSize;
      Num_copy_Rprimme(*numPrevRitzVals, sortedRitzVals, 1, prevRitzVals, 1);

   } /* user provided shifts */
   else {    
   /*-----------------------------------------------------------------*/
   /* else it is primme_smallest or primme_largest                    */
   /*-----------------------------------------------------------------*/

      if (primme->correctionParams.robustShifts) { 
         /*----------------------------------------------------*/
         /* Subtract/add a robust shift from/to the Ritz value */
         /*----------------------------------------------------*/

         /* Find the robust shift for each block vector */
         for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
   
            sortedIndex = ilev[blockIndex];
            eval = sortedRitzVals[sortedIndex];
   
            robustShift = computeRobustShift(blockIndex, 
              blockNorms[blockIndex], prevRitzVals, *numPrevRitzVals, 
              sortedRitzVals, &approxOlsenEps[blockIndex], 
              numLocked+basisSize, ilev, primme);
   
            /* Subtract/add the shift if looking for the smallest/largest  */
            /* eigenvalues, Do not go beyond the previous computed eigval  */
       
            if (primme->target == primme_smallest) {
               blockOfShifts[blockIndex] = eval - robustShift;
               if (sortedIndex > 0) blockOfShifts[blockIndex] = 
                  max(blockOfShifts[blockIndex], sortedRitzVals[sortedIndex-1]);
            }
            else {
               blockOfShifts[blockIndex] = eval + robustShift;
               if (sortedIndex > 0) blockOfShifts[blockIndex] = 
                  min(blockOfShifts[blockIndex], sortedRitzVals[sortedIndex-1]);
            } /* robust shifting */
   
         }  /* for loop */
   
      }  /* endif robust shifts */
      else {
         /*--------------------------------------------------------------*/
         /* Otherwise, the shifts for both preconditioner and correction */
         /* equation should be just the Ritz values. For Olsen's method, */
         /* the shifts for r-eps*x, are chosen as the difference in Ritz */
         /* value between successive iterations.                         */
         /*--------------------------------------------------------------*/
   
         for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
            ritzIndex   =  iev[blockIndex];
            sortedIndex = ilev[blockIndex];
            blockOfShifts[blockIndex] = ritzVals[ritzIndex];
            if (sortedIndex < *numPrevRitzVals) {
               approxOlsenEps[blockIndex] = 
               fabs(prevRitzVals[sortedIndex] - sortedRitzVals[sortedIndex]);
            }
            else {
               approxOlsenEps[blockIndex] = blockNorms[blockIndex]; 
            }
         } /* for loop */
      } /* else no robust shifts */

      /* Remember the previous ritz values*/
      *numPrevRitzVals = numLocked+basisSize;
      Num_copy_Rprimme(*numPrevRitzVals, sortedRitzVals, 1, prevRitzVals, 1);

   } /* else primme_smallest or primme_largest */

   /* Equip the primme struct with the blockOfShifts, in case the user */
   /* wants to precondition (K-sigma_i I)^{-1} with a different shift  */
   /* for each vector                                                  */

   primme->ShiftsForPreconditioner = blockOfShifts;

   /*------------------------------------------------------------ */
   /*  Generalized Davidson variants -- No inner iterations       */
   /*------------------------------------------------------------ */
   if (primme->correctionParams.maxInnerIterations == 0) {
      /* This is Generalized Davidson or approximate Olsen's method. */
      /* Perform block preconditioning (with or without projections) */
      
      r = &W[ldW*basisSize];    /* All the block residuals    */
      x = &V[ldV*basisSize];    /* All the block Ritz vectors */
      
      if ( primme->correctionParams.projectors.RightX &&
           primme->correctionParams.projectors.SkewX    ) {    
           /* Compute exact Olsen's projected preconditioner. This is */
          /* expensive and rarely improves anything! Included for completeness*/
          
          Olsen_preconditioner_block(r, ldW, x, ldV, blockSize, Kinvx, primme);
      }
      else {
         if ( primme->correctionParams.projectors.RightX ) {   
            /*Compute a cheap approximation to OLSENS, where (x'Kinvr)/xKinvx */
            /*is approximated by e: Kinvr-e*Kinvx=Kinv(r-e*x)=Kinv(I-ct*x*x')r*/

            for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
               /* Compute r_i = r_i - err_i * x_i */
               Num_axpy_Sprimme(primme->nLocal, -approxOlsenEps[blockIndex],
               &x[ldV*blockIndex],1,&r[ldW*blockIndex],1);
            } /* for */
         }

         /* GD: compute K^{-1}r , or approx.Olsen: K^{-1}(r-ex) */

         CHKERR(applyPreconditioner_Sprimme(r, primme->nLocal, ldW,
                  x, ldV, blockSize, primme), -1);
      }
   }
   /* ------------------------------------------------------------ */
   /*  JDQMR --- JD inner-outer variants                           */
   /* ------------------------------------------------------------ */
   else {  /* maxInnerIterations > 0  We perform inner-outer JDQMR */
      int touch0 = *touch;

      /* Solve the correction for each block vector. */

      for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {

         r = &W[ldW*(basisSize+blockIndex)];
         x = &V[ldV*(basisSize+blockIndex)];

         /* Set up the left/right/skew projectors for JDQMR.        */
         /* The pointers Lprojector, Rprojector(Q/X) point to the   */
         /* appropriate arrays for use in the projection step       */

         CHKERR(setup_JD_projectors(x, evecs, ldevecs, evecsHat, ldevecsHat,
                  Kinvx, &xKinvx, &Lprojector, &ldLprojector, &RprojectorQ,
                  &ldRprojectorQ, &RprojectorX, &ldRprojectorX, &sizeLprojector,
                  &sizeRprojectorQ, &sizeRprojectorX, numLocked,
                  numConvergedStored, primme), -1);

         /* Map the index of the block vector to its corresponding eigenvalue */
         /* index, and the shift for the correction equation. Also make the   */
         /* shift available to primme, in case (K-shift I)^-1 is needed       */

         ritzIndex = iev[blockIndex];
         shift = blockOfShifts[blockIndex];
         primme->ShiftsForPreconditioner = &blockOfShifts[blockIndex];

         /* Pass the original value of touch and update touch as the maximum  */
         /* value that takes for all inner_solve calls                        */
         int touch1 = touch0;

         CHKERR(inner_solve_Sprimme(x, r, &blockNorms[blockIndex], evecs,
                  ldevecs, UDU, ipivot, &xKinvx,
                  Lprojector, ldLprojector, RprojectorQ, ldRprojectorQ,
                  RprojectorX, ldRprojectorX, sizeLprojector, sizeRprojectorQ,
                  sizeRprojectorX, sol, ritzVals[ritzIndex], shift, &touch1,
                  machEps, linSolverRWork, linSolverRWorkSize,
                  primme), -1);
         *touch = max(*touch, touch1);

         Num_copy_Sprimme(primme->nLocal, sol, 1, 
            &V[ldV*(basisSize+blockIndex)], 1);

      } /* end for each block vector */
   } /* JDqmr variants */

   return 0;

}
      

/*******************************************************************************
 * Subroutine computeRobustShift - This function computes the robust shift
 *    to be used in the correction equation.  The standard shift is the current
 *    Ritz value.  The calling routine eventually adds the robust shift, 
 *    epsilon, to the standard shift to accelerate the convergence of the outer 
 *    method.  
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * blockIndex    Index of the block vector for which the correction will be 
 *               solved.
 *
 * resNorm       The residual norm of the current Ritz vector
 *
 * prevRitzVals  Sorted array of locked Ritz values and Ritz values from the 
 *               previous outer iteration.  Array size is numPrevRitzVals.
 *
 * numPrevRitzVals The number of values in prevRitzVals
 * 
 * sortedRitzVals  Sorted array of locked Ritz values and current Ritz values.
 *                 The array is of size numLocked+basisSize
 *
 * approxOlsenShift The shift to be used to modify r-shift x, for approx Olsen's
 *
 * numSorted     The size of sortedRitzVals.
 *
 * ilev          Array of size blockSize that maps targeted Ritz values to
 *               their position with the array sortedRitzVals.
 *
 * Return Value
 * ------------
 * double  The value to be added to the current shift.  The calling routine
 *         will compute shift := shift +/- epsilon.
 ******************************************************************************/

static REAL computeRobustShift(int blockIndex, double resNorm, 
   REAL *prevRitzVals, int numPrevRitzVals, REAL *sortedRitzVals, 
   REAL *approxOlsenShift, int numSorted, int *ilev, primme_params *primme) {

   int sortedIndex;                 /* Index of the current Ritz value in */
                                    /* sortedRitzVals.                    */
   REAL gap, lowerGap, upperGap;  /* Gaps between the current and       */
                                    /* neighboring Ritz Values.           */ 
   REAL delta;                    /* The difference between the current */
                                    /* Ritz value and its value in the    */
                                    /* previous iteration.                */
   REAL epsilon;                  /* The return value.  The calling     */
                                    /* routine will modify the shift by   */
                                    /* this amount.                       */    

   /* If this is the 1st outer iteration, there are no previous */
   /* Ritz values. Return the residual norm                     */

   if (primme->stats.numOuterIterations <= 1) {
      *approxOlsenShift = resNorm;
      return resNorm;
   }

   sortedIndex = ilev[blockIndex];

   /* Compute the gap when the first eigenvalue with respect to the */
   /* current basis is to be computed.                              */

   if (sortedIndex == 0 && numSorted >= 2) {
      lowerGap = HUGE_VAL;
      gap = fabs(sortedRitzVals[1] - sortedRitzVals[0]);
   }
   else if (sortedIndex > 0 && numSorted >= 2 && sortedIndex+1 < numSorted) {

      /* Take the smaller of the two gaps if an interior eigenvalue is */
      /* targeted.                                                     */

      lowerGap = fabs(sortedRitzVals[sortedIndex] - 
                      sortedRitzVals[sortedIndex-1]);
      upperGap = fabs(sortedRitzVals[sortedIndex+1] - 
                      sortedRitzVals[sortedIndex]);
      gap = min(lowerGap, upperGap);
   }
   else {
      lowerGap = fabs(sortedRitzVals[sortedIndex] -
                      sortedRitzVals[sortedIndex-1]);
      gap = lowerGap;
   }
   
   /* Compute the change in a Ritz value between successive iterations */

   if (sortedIndex < numPrevRitzVals) {
      delta = fabs(prevRitzVals[sortedIndex] - sortedRitzVals[sortedIndex]);
   }
   else {
      delta = HUGE_VAL;
   }

   /* Compute epsilon according to the theory of Davis and Kahan described */
   /* in The Symmetric Eigenvalue Problem by B.N. Parlett.                 */

   if (gap > resNorm) {
      epsilon = min(delta, min(resNorm*resNorm/gap, lowerGap));
   }
   else {
      epsilon = min(resNorm, lowerGap);
   }

   *approxOlsenShift = min(delta, epsilon);

   /* If the above is too large a shift set it to a milder shift */
   /* epsilon = min(delta, epsilon); */

   return epsilon;

}


/*******************************************************************************
 * Subroutine mergeSort -- This routine merges the sorted lockedEvals array and 
 *   the sorted ritzVals array into a single sorted array, sortedRitzVals.
 *   It is only called for extreme (largest, smallest) eigenvalues, not for 
 *   interior ones. Thus it does not cover the interior case. 
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * lockedEvals  Array of size numLocked.  Ritz values that have been locked.
 *
 * numLocked    The size of array lockedEvals.
 *
 * ritzVals     The current Ritz values.  They are the eigenvalues of the
 *              projection V'*A*V.
 *
 * flags        Array of flags indicating which Ritz values in ritzVals are
 *              converged.
 *
 * basisSize    The current size of the basis V.
 *
 * ilev         Maps the blockSize targeted Ritz values to their position
 *              within sortedRitzVals.
 * 
 * blockSize    The number of Ritz values targeted during the current outer
 *              iteration.
 *
 * primme       A structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * sortedRitzVals  The array lockedEvals and ritzVals are merged into.
 *
 ******************************************************************************/


static void mergeSort(REAL *lockedEvals, int numLocked, REAL *ritzVals, 
   int *flags, int basisSize, REAL *sortedRitzVals, int *ilev, int blockSize,
   primme_params *primme) {
   
   int count;         /* The number of Ritz values merged        */
   int eval, ritzVal; /* The index of the next locked Ritz value */
                      /*   and current Ritz values to be merged  */
   int blockIndex;    /* Counter used to index ilev              */

   count = 0;
   eval = 0;
   ritzVal = 0;
   blockIndex = 0;
   
   /* Continue merging until all elements have been merged */

   while (count < numLocked+basisSize) {

      /* I. The list of locked eigenvalues has been completely merged, or  */
      /* if ritzVals[ritzVal] is greater/less than lockedEvals[eval], then */
      /* merge ritzVals[ritzVal] into the list of sorted Ritz values.      */
      /*                                                                   */
      /* II. If the list of current Ritz values has been completely merged,*/
      /* or if lockedEvals[eval] greater/less than ritzVals[ritzVal], then */
      /* merge lockedEvals[eval] into the list of sorted Ritz values.      */
  
      if (eval >= numLocked || (primme->target == primme_largest && 
                               ritzVals[ritzVal] >= lockedEvals[eval])
                           || (primme->target == primme_smallest &&
                               ritzVals[ritzVal] <= lockedEvals[eval])
         ) 
      {
         sortedRitzVals[count] = ritzVals[ritzVal];

         /* If the Ritz value just merged is unconverged and there is room */
         /* enough in the block to target it, then keep track of its index */
         /* in the sorted list.                                            */

         if (blockIndex < blockSize && flags[ritzVal] == UNCONVERGED) {
            ilev[blockIndex] = count;
            blockIndex++;
         }

         ritzVal++;
      }
      else if (ritzVal >= basisSize || (primme->target == primme_largest &&
                                       lockedEvals[eval] >= ritzVals[ritzVal])
                                   || (primme->target == primme_smallest &&
                                       lockedEvals[eval] <= ritzVals[ritzVal])
              )
      {
         sortedRitzVals[count] = lockedEvals[eval];
         eval++;
      } 

      count++;
   }

}

/*******************************************************************************
 * Subroutine Olsen_preconditioner_block - This subroutine applies the projected
 *    preconditioner to a block of blockSize vectors r by computing:
 *       (I - (K^{-1}x_i)x_i^T / (x_i^T K^{-1}x_i) ) K^{-1}r_i
 *
 * Input Parameters
 * ----------------
 * r          The vectors the preconditioner and projection will be applied to.
 *
 * blockSize  The number of vectors in r, x
 *
 * rwork      SCALAR work array of size (primme.nLocal + 4*blockSize)
 *
 * primme       Structure containing various solver parameters
 *
 * Output parameters
 * -----------------
 * x The current Ritz vector used in the projection. It also carries the result.
 *
 ******************************************************************************/

static int Olsen_preconditioner_block(SCALAR *r, PRIMME_INT ldr, SCALAR *x,
      PRIMME_INT ldx, int blockSize, SCALAR *rwork, primme_params *primme) {

   int blockIndex;
   SCALAR alpha;
   SCALAR *Kinvx, *xKinvx, *xKinvr, *xKinvx_local, *xKinvr_local;

   /*------------------------------------------------------------------*/
   /* Subdivide workspace                                              */
   /*------------------------------------------------------------------*/
   Kinvx = rwork;
   xKinvx_local = Kinvx + ldx*blockSize;
   xKinvr_local = xKinvx_local + blockSize;
   xKinvx       = xKinvr_local + blockSize;
   xKinvr       = xKinvx + blockSize;
          
   /*------------------------------------------------------------------ */
   /* Compute K^{-1}x for block x. Kinvx memory requirement (blockSize*nLocal)*/
   /*------------------------------------------------------------------ */

   CHKERR(applyPreconditioner_Sprimme(x, primme->nLocal, ldx, Kinvx, ldx,
            blockSize, primme), -1);

   /*------------------------------------------------------------------ */
   /* Compute local x^TK^{-1}x and x^TK^{-1}r = (K^{-1}x)^Tr for each vector */
   /*------------------------------------------------------------------ */

   for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
      xKinvx_local[blockIndex] =
        Num_dot_Sprimme(primme->nLocal, &x[ldx*blockIndex],1, 
                           &Kinvx[ldx*blockIndex],1);
      xKinvr_local[blockIndex] =
        Num_dot_Sprimme(primme->nLocal, &Kinvx[ldx*blockIndex],1,
                                   &r[ldr*blockIndex],1);
   }      
   CHKERR(globalSum_Sprimme(xKinvx_local, xKinvx, 2*blockSize, primme), -1);

   /*------------------------------------------------------------------*/
   /* Compute K^{-1}r                                                  */
   /*------------------------------------------------------------------*/

   CHKERR(applyPreconditioner_Sprimme(r, primme->nLocal, ldr, x, ldx,
            blockSize, primme), -1);

   /*------------------------------------------------------------------*/
   /* Compute K^(-1)r  - ( xKinvr/xKinvx ) K^(-1)r for each vector     */
   /*------------------------------------------------------------------*/

   for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
      if (ABS(xKinvx[blockIndex]) > 0.0L) 
         alpha = - xKinvr[blockIndex]/xKinvx[blockIndex];
      else   
         alpha = 0.0;

      Num_axpy_Sprimme(primme->nLocal,alpha,&Kinvx[ldx*blockIndex],
                                       1, &x[ldx*blockIndex],1);
   } /*for*/

   return 0;

} /* of Olsen_preconditiner_block */

/*******************************************************************************
 *   subroutine setup_JD_projectors()
 *
 *   Sets up whether to use any of the following projections in the 
 *   correction equation:
 *
 *  INPUT
 *  -----
 *   x                The Ritz vector
 *   evecs            Converged locked eigenvectors (denoted as Q herein)
 *   evecsHat         K^{-1}*evecs
 *   numLocked        Number of locked eigenvectors (if locking)
 *   numConverged     Number of converged e-vectors copied in evecs (no locking)
 *   primme           The main data structures that contains the choices for
 *
 *       primme->LeftQ  : evecs in the left projector
 *       primme->LeftX  : x in the left projector
 *       primme->RightQ : evecs in the right projector
 *       primme->RightX : x in the right projector
 *       primme->SkewQ  : evecs in the right one, should be skewed
 *       primme->SkewX  : x in the right one, should be skewed
 *                   
 *  OUTPUT
 *  ------
 *  *Kinvx            The result of K^{-1}x (if needed, otherwise NULL)
 * **Lprojector       Pointer to the left projector for [Q x] (could be NULL)
 * **RprojectorQ      Pointer to the right projector for Q (could be NULL)
 * **RprojectorX      Pointer to the right projector for X (could be NULL)
 *   sizeLprojector   Size of the left projector(numConverged/numLocked)/+1 or 0
 *   sizeRprojectorQ  Size of the Q right projectr (numConverged/numLocked or 0)
 *   sizeRprojectorX  Size of the X right projectr (1 or 0)
 *
 * ============================================================================
 * Functionality:
 *
 *   Q = locked vectors
 *   x = Ritz vector
 *   K = preconditioner
 *   r = a vector to apply the 
 *
 *   (I-QQ')(I-xx')(A-shift I)(I-Kx(x'Kx)^(-1)x')(I-KQ(Q'K Q)^(-1)Q') K 
 *    -----  -----             -----------------  ------------------
 *     Lq     Lx                   Rx  Sx              Rq Sq
 *
 * Control parameters from primme (shortened here, see above for full names):
 *
 *     Lq = 1/0 whether to apply the left projector I-QQ'
 *     Lx = 1/0 whether to apply the left projector I-xx'
 *     Rx = 1/0 whether to apply the right projector with x (skew or orthogonal)
 *     Rq = 1/0 whether to apply the right projector with Q (skew or orthogonal)
 *     if (Rx == 1)
 *        Sx = 1 applies the right skew projector (I-Kx(x'Kx)^(-1)x')
 *        Sx = 0 applies the right orthogonal projector (I-xx')
 *     if (Rq == 1)
 *        Sq = 1 applies the right skew projector (I-KQ(Q'K Q)^(-1)Q') 
 *        Sq = 0 applies the right orthogonal projector (I-QQ')
 *
 * Examples
 * Lq,Lx,Rx,Sx,Rq,Sq
 *  1  1  1  1  1  1  (the above equation)                   full, expensive JD
 *  0  0  1  1  1  1  (I-Kx(x'Kx)^(-1)x')(I-KQ(Q'K Q)^(-1)Q')        classic JD
 *  0  1  1  1  0  0  (I-xx')(A-shift I)(I-Kx(x'Kx)^(-1)x')     JD w/o deflaton
 *  0  0  0  0  0  0  (A-shift I) K                                  classic GD
 *  1  1  1  1  0  0  (I-QQ')(I-xx')(A-shift I)(I-Kx(x'Kx)^(-1)x')  RECOMMENDED
 *  1  1  0  0  0  0  (I-QQ')(I-xx')(A-shift I)                         similar
 *      Others...
 *                    Researchers can experiment with other projection schemes,
 *                    although our experience says they are rarely beneficial
 *
 * The left orthogonal projector for x and Q can be performed as one block
 * containing [Q x]. However, the right projections (if either is skew)
 * are performed separately for Q and x. There are memory reasons for 
 * doing so, but also we do not have to factor (Q'KQ) at every outer step;
 * only when an eval converges. 
 *
 ******************************************************************************/

static int setup_JD_projectors(SCALAR *x, SCALAR *evecs, PRIMME_INT ldevecs,
      SCALAR *evecsHat, PRIMME_INT ldevecsHat, SCALAR *Kinvx, SCALAR *xKinvx, 
      SCALAR **Lprojector, PRIMME_INT *ldLprojector, SCALAR **RprojectorQ,
      PRIMME_INT *ldRprojectorQ, SCALAR **RprojectorX,
      PRIMME_INT *ldRprojectorX,  int *sizeLprojector, int *sizeRprojectorQ,
      int *sizeRprojectorX, int numLocked, int numConverged,
      primme_params *primme) {

   int n, sizeEvecs;
   SCALAR xKinvx_local;

   *sizeLprojector  = 0;
   *sizeRprojectorQ = 0;
   *sizeRprojectorX = 0;
   *ldLprojector  = 0;
   *ldRprojectorQ = 0;
   *ldRprojectorX = 0;
   *Lprojector  = NULL;
   *RprojectorQ = NULL;
   *RprojectorX = NULL;

   n = primme->nLocal;
   if (primme->locking) 
      sizeEvecs = primme->numOrthoConst+numLocked;
   else
      sizeEvecs = primme->numOrthoConst+numConverged;
   
   /* --------------------------------------------------------*/
   /* Set up the left projector arrays. Make x adjacent to Q. */
   /* --------------------------------------------------------*/
   
   if (primme->correctionParams.projectors.LeftQ) {
   
         *sizeLprojector = sizeEvecs;
         *Lprojector = evecs;
         *ldLprojector = ldevecs;
         if (primme->correctionParams.projectors.LeftX) {
            Num_copy_Sprimme(n, x, 1, &evecs[sizeEvecs*ldevecs], 1);
            *sizeLprojector = *sizeLprojector + 1;
         }
   }
   else {
      if (primme->correctionParams.projectors.LeftX) {
         *Lprojector = x;
         *ldLprojector = n;
         *sizeLprojector = 1;
      }
   }
      
   /* --------------------------------------------------------*/
   /* Set up the right projector arrays. Q and x separately   */
   /* --------------------------------------------------------*/
   
   /* ------------*/
   /* First for Q */
   /* ------------*/
   if (primme->correctionParams.projectors.RightQ) {
   
      if (primme->correctionParams.precondition    &&
          primme->correctionParams.projectors.SkewQ) {
         *RprojectorQ = evecsHat;       /* Use the K^(-1)evecs array */
         *ldRprojectorQ = ldevecsHat;
      }                               
      else {  /* Right Q but not SkewQ */
         *RprojectorQ = evecs;          /* Use just the evecs array. */
         *ldRprojectorQ = ldevecs;
      }
      *sizeRprojectorQ = sizeEvecs;
   
   } 
   else { /* if no RightQ projector */
   
      *RprojectorQ = NULL;       
      *sizeRprojectorQ = 0;
   }
   
   /* ------------*/
   /* Then for x  */
   /* ------------*/
   if (primme->correctionParams.projectors.RightX) {
   
      if (primme->correctionParams.precondition   &&
          primme->correctionParams.projectors.SkewX) {
         CHKERR(applyPreconditioner_Sprimme(x, primme->nLocal, primme->nLocal,
                  Kinvx, primme->nLocal, 1, primme), -1);
         primme->stats.numPreconds += 1;
         *RprojectorX  = Kinvx;
         xKinvx_local = Num_dot_Sprimme(primme->nLocal, x, 1, Kinvx, 1);
         CHKERR(globalSum_Sprimme(&xKinvx_local, xKinvx, 1, primme), -1);
      }      
      else {
         *RprojectorX = x;
         *xKinvx = 1.0;
      }
      *sizeRprojectorX = 1;
      *ldRprojectorX  = n;
   }
   else { 
         *RprojectorX = NULL;
         *sizeRprojectorX = 0;
         *xKinvx = 1.0;
   }

   return 0;

} /* setup_JD_projectors */
