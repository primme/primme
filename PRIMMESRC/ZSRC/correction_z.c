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
 * File: correction.c
 *
 * Purpose - Computes the correction corresponding to each block Ritz vectors.
 *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include "primme.h"
#include "const.h"
#include "correction_z.h"
#include "correction_private_z.h"
#include "inner_solve_z.h"
#include "numerical_z.h"

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
 * eresTol        the eigenvalue residual tolerance
 *
 * machEps        machine precision 
 *
 * aNormEstimate if primme->aNorm<=0, eresTol*aNormEstimate (=largestRitzValue)
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
 *
 * Return Value
 * ------------
 * int  Error code: 0 upon success, nonzero otherwise
 *                 -1 innner solver failure
 *                 >0 the needed rworkSize, if the given was not enough
 *
 ******************************************************************************/
 

int solve_correction_zprimme(Complex_Z *V, Complex_Z *W, Complex_Z *evecs, 
   Complex_Z *evecsHat, Complex_Z *UDU, int *ipivot, double *lockedEvals, 
   int numLocked, int numConvergedStored, double *ritzVals, 
   double *prevRitzVals, int *numPrevRitzVals, int *flags, int basisSize, 
   double *blockNorms, int *iev, int blockSize, double eresTol, 
   double machEps, double aNormEstimate, Complex_Z *rwork, int *iwork, 
   int rworkSize, primme_params *primme) {

   int blockIndex;         /* Loop index.  Ranges from 0..blockSize-1.       */
   int ritzIndex;          /* Ritz value index blockIndex corresponds to.    */
                           /* Possible values range from 0..basisSize-1.     */
   int sortedIndex;        /* Ritz value index in sortedRitzVals, blockIndex */
                           /* corresponds to. Range 0..numLocked+basisSize-1 */
   int neededRsize;        /* Needed size for rwork. If not enough return    */
   int linSolverRWorkSize; /* Size of the linSolverRWork array.              */
   int *ilev;              /* Array of size blockSize.  Maps the target Ritz */
                           /* values to their positions in the sortedEvals   */
                           /* array.                                         */
   int sizeLprojector;     /* Sizes of the various left/right projectors     */
   int sizeRprojectorQ;    /* These will be 0/1/or numOrthConstr+numLocked   */
   int sizeRprojectorX;    /* or numOrthConstr+numConvergedStored w/o locking*/

   int ret;                /* Return code.                                   */
   Complex_Z *r, *x, *sol;  /* Residual, Ritz vector, and correction.         */
   Complex_Z *linSolverRWork;/* Workspace needed by linear solver.            */
   double *sortedRitzVals; /* Sorted array of current and converged Ritz     */
                           /* values.  Size of array is numLocked+basisSize. */
   double *blockOfShifts;  /* Shifts for (A-shiftI) or (if needed) (K-shiftI)*/
   double *approxOlsenEps; /* Shifts for approximate Olsen implementation    */
   Complex_Z *Kinvx;         /* Workspace to store K^{-1}x                     */
   Complex_Z *Lprojector;   /* Q pointer for (I-Q*Q'). Usually points to evecs*/
   Complex_Z *RprojectorQ;  /* May point to evecs/evecsHat depending on skewQ */
   Complex_Z *RprojectorX;  /* May point to x/Kinvx depending on skewX        */

   Complex_Z xKinvx;                        /* Stores x'*K^{-1}x if needed    */
   double eval, shift, robustShift;       /* robust shift values.           */
   Complex_Z tmpShift;                      /* Temp shift for daxpy           */

   /*------------------------------------------------------------*/
   /* Subdivide the workspace with pointers, and figure out      */
   /* the total amount of needed real workspace (neededRsize)    */
   /*------------------------------------------------------------*/

   /* needed worksize */
   neededRsize = 0;
   Kinvx       = rwork;
   /* Kinvx will have nonzero size if precond and both RightX and SkewX */
   if (primme->correctionParams.precondition &&        
       primme->correctionParams.projectors.RightX &&  
       primme->correctionParams.projectors.SkewX ) { 

      /* OLSEN's method requires a block, but JDQMR is vector by vector */
      if (primme->correctionParams.maxInnerIterations == 0) {    
         sol = Kinvx + primme->nLocal*blockSize;
         neededRsize = neededRsize + primme->nLocal*blockSize;
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
   sortedRitzVals = (double *)(linSolverRWork + linSolverRWorkSize);
   blockOfShifts  = sortedRitzVals + (numLocked+basisSize);
   approxOlsenEps = blockOfShifts  + blockSize;
   neededRsize = neededRsize + numLocked+basisSize + 2*blockSize;

   if (neededRsize > rworkSize) {
      return(neededRsize);
   }

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
         sortedIndex = ilev[blockIndex];
         blockOfShifts[blockIndex] = 
            primme->targetShifts[ min(primme->numTargetShifts-1, numLocked) ];
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
      Num_dcopy_primme(*numPrevRitzVals, sortedRitzVals, 1, prevRitzVals, 1);

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
      Num_dcopy_primme(*numPrevRitzVals, sortedRitzVals, 1, prevRitzVals, 1);

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
      
      r = &W[primme->nLocal*basisSize];    /* All the block residuals    */
      x = &V[primme->nLocal*basisSize];    /* All the block Ritz vectors */
      
      if ( primme->correctionParams.projectors.RightX &&
           primme->correctionParams.projectors.SkewX    ) {    
           /* Compute exact Olsen's projected preconditioner. This is */
          /* expensive and rarely improves anything! Included for completeness*/
          
          Olsen_preconditioner_block(r, x, blockSize, Kinvx, primme);
      }
      else {
         if ( primme->correctionParams.projectors.RightX ) {   
            /*Compute a cheap approximation to OLSENS, where (x'Kinvr)/xKinvx */
            /*is approximated by e: Kinvr-e*Kinvx=Kinv(r-e*x)=Kinv(I-ct*x*x')r*/

            for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
               /* Compute r_i = r_i - err_i * x_i */
               {tmpShift.r = -approxOlsenEps[blockIndex]; tmpShift.i = 0.0L;}
               Num_axpy_zprimme(primme->nLocal, tmpShift,
               &x[primme->nLocal*blockIndex],1,&r[primme->nLocal*blockIndex],1);
            } /* for */
         }

         /* GD: compute K^{-1}r , or approx.Olsen: K^{-1}(r-ex) */

         apply_preconditioner_block(r, x, blockSize, primme );
      }
   }
   /* ------------------------------------------------------------ */
   /*  JDQMR --- JD inner-outer variants                           */
   /* ------------------------------------------------------------ */
   else {  /* maxInnerIterations > 0  We perform inner-outer JDQMR */

      /* Solve the correction for each block vector. */

      for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {

         r = &W[primme->nLocal*(basisSize+blockIndex)];
         x = &V[primme->nLocal*(basisSize+blockIndex)];

         /* Set up the left/right/skew projectors for JDQMR.        */
         /* The pointers Lprojector, Rprojector(Q/X) point to the   */
         /* appropriate arrays for use in the projection step       */

         setup_JD_projectors(x, r, evecs, evecsHat, Kinvx, &xKinvx, 
            &Lprojector, &RprojectorQ, &RprojectorX, 
            &sizeLprojector, &sizeRprojectorQ, &sizeRprojectorX,
            numLocked, numConvergedStored, primme);

         /* Map the index of the block vector to its corresponding eigenvalue */
         /* index, and the shift for the correction equation. Also make the   */
         /* shift available to primme, in case (K-shift I)^-1 is needed       */

         ritzIndex = iev[blockIndex];
         shift = blockOfShifts[blockIndex];
         primme->ShiftsForPreconditioner = &blockOfShifts[blockIndex];

         ret = inner_solve_zprimme(x, r, &blockNorms[blockIndex], evecs, 
            evecsHat, UDU, ipivot, &xKinvx, Lprojector, RprojectorQ, 
            RprojectorX, sizeLprojector, sizeRprojectorQ, sizeRprojectorX,
            sol, ritzVals[ritzIndex], shift, eresTol, aNormEstimate, 
            machEps, linSolverRWork, linSolverRWorkSize, primme);

         if (ret != 0) {
            primme_PushErrorMessage(Primme_solve_correction, Primme_inner_solve,
                            ret, __FILE__, __LINE__, primme);
            return (INNER_SOLVE_FAILURE);
         }

         Num_zcopy_zprimme(primme->nLocal, sol, 1, 
            &V[primme->nLocal*(basisSize+blockIndex)], 1);

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

static double computeRobustShift(int blockIndex, double resNorm, 
   double *prevRitzVals, int numPrevRitzVals, double *sortedRitzVals, 
   double *approxOlsenShift, int numSorted, int *ilev, primme_params *primme) {

   int sortedIndex;                 /* Index of the current Ritz value in */
                                    /* sortedRitzVals.                    */
   double gap, lowerGap, upperGap;  /* Gaps between the current and       */
                                    /* neighboring Ritz Values.           */ 
   double delta;                    /* The difference between the current */
                                    /* Ritz value and its value in the    */
                                    /* previous iteration.                */
   double epsilon;                  /* The return value.  The calling     */
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
      gap = Num_fmin_primme(2, lowerGap, upperGap);
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
      epsilon = Num_fmin_primme(3, delta, resNorm*resNorm/gap, lowerGap);
   }
   else {
      epsilon = Num_fmin_primme(2, resNorm, lowerGap);
   }

   *approxOlsenShift = Num_fmin_primme(2, delta, epsilon);

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


static void mergeSort(double *lockedEvals, int numLocked, double *ritzVals, 
   int *flags, int basisSize, double *sortedRitzVals, int *ilev, int blockSize,
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
 * Subroutine apply_preconditioner_block - This subroutine applies the 
 *    preconditioner to a block of vectors v by computing: K^{-1}v
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

static void apply_preconditioner_block(Complex_Z *v, Complex_Z *result, 
                int blockSize, primme_params *primme) {
         
   if (primme->correctionParams.precondition) {

      (*primme->applyPreconditioner)(v, result, &blockSize, primme);
      primme->stats.numPreconds += blockSize;
   }
   else {
      Num_zcopy_zprimme(primme->nLocal*blockSize, v, 1, result, 1);
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
 * rwork      Complex_Z work array of size (primme.nLocal + 4*blockSize)
 *
 * primme       Structure containing various solver parameters
 *
 * Output parameters
 * -----------------
 * x The current Ritz vector used in the projection. It also carries the result.
 *
 ******************************************************************************/

static void Olsen_preconditioner_block(Complex_Z *r, Complex_Z *x,
                int blockSize, Complex_Z *rwork, primme_params *primme) {

   int blockIndex, count;
   Complex_Z alpha;
   Complex_Z *Kinvx, *xKinvx, *xKinvr, *xKinvx_local, *xKinvr_local;
   Complex_Z ztmp;
   Complex_Z tzero = {+0.0e+00,+0.0e00};

   /*------------------------------------------------------------------*/
   /* Subdivide workspace                                              */
   /*------------------------------------------------------------------*/
   Kinvx = rwork;
   xKinvx_local = Kinvx + primme->nLocal*blockSize;
   xKinvr_local = xKinvx_local + blockSize;
   xKinvx       = xKinvr_local + blockSize;
   xKinvr       = xKinvx + blockSize;
          
   /*------------------------------------------------------------------ */
   /* Compute K^{-1}x for block x. Kinvx memory requirement (blockSize*nLocal)*/
   /*------------------------------------------------------------------ */

   apply_preconditioner_block(x, Kinvx, blockSize, primme );

   /*------------------------------------------------------------------ */
   /* Compute local x^TK^{-1}x and x^TK^{-1}r = (K^{-1}x)^Tr for each vector */
   /*------------------------------------------------------------------ */

   for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
      xKinvx_local[blockIndex] =
        Num_dot_zprimme(primme->nLocal, &x[primme->nLocal*blockIndex],1, 
                           &Kinvx[primme->nLocal*blockIndex],1);
      xKinvr_local[blockIndex] =
        Num_dot_zprimme(primme->nLocal, &Kinvx[primme->nLocal*blockIndex],1,
                                   &r[primme->nLocal*blockIndex],1);
   }      
   count = 4*blockSize;
   (*primme->globalSumDouble)(xKinvx_local, xKinvx, &count, primme);

   /*------------------------------------------------------------------*/
   /* Compute K^{-1}r                                                  */
   /*------------------------------------------------------------------*/

   apply_preconditioner_block(r, x, blockSize, primme );

   /*------------------------------------------------------------------*/
   /* Compute K^(-1)r  - ( xKinvr/xKinvx ) K^(-1)r for each vector     */
   /*------------------------------------------------------------------*/

   for (blockIndex = 0; blockIndex < blockSize; blockIndex++) {
      if (z_abs_primme(xKinvx[blockIndex]) > 0.0L) 
      {
         ztmp.r = -xKinvr[blockIndex].r;
         ztmp.i = -xKinvr[blockIndex].i;
         z_div_primme(&alpha, &ztmp, &xKinvx[blockIndex]);
      }
      else   
         alpha = tzero;

      Num_axpy_zprimme(primme->nLocal,alpha,&Kinvx[primme->nLocal*blockIndex],
                                       1, &x[primme->nLocal*blockIndex],1);
   } /*for*/

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
 *   r                The residual vector for x
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

static void setup_JD_projectors(Complex_Z *x, Complex_Z *r, Complex_Z *evecs, 
   Complex_Z *evecsHat, Complex_Z *Kinvx, Complex_Z *xKinvx, 
   Complex_Z **Lprojector, Complex_Z **RprojectorQ, Complex_Z **RprojectorX, 
   int *sizeLprojector, int *sizeRprojectorQ, int *sizeRprojectorX, 
   int numLocked, int numConverged, primme_params *primme) {

   int n, sizeEvecs;
   int ONE = 1;
   /* In Complex, the size of the array to globalSum is twice as large */
   int count = 2;
   Complex_Z xKinvx_local;
   Complex_Z tpone = {+1.0e+00,+0.0e00};

   *sizeLprojector  = 0;
   *sizeRprojectorQ = 0;
   *sizeRprojectorX = 0;
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
         if (primme->correctionParams.projectors.LeftX) {
            Num_zcopy_zprimme(n, x, 1, &evecs[sizeEvecs*n], 1);
            *sizeLprojector = *sizeLprojector + 1;
         }
   }
   else {
      if (primme->correctionParams.projectors.LeftX) {
         *Lprojector = x;
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
      }                               
      else {  /* Right Q but not SkewQ */
         *RprojectorQ = evecs;          /* Use just the evecs array. */
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
         (*primme->applyPreconditioner)(x, Kinvx, &ONE, primme);
         primme->stats.numPreconds += 1;
         *RprojectorX  = Kinvx;
         xKinvx_local = Num_dot_zprimme(primme->nLocal, x, 1, Kinvx, 1);
         (*primme->globalSumDouble)(&xKinvx_local, xKinvx, &count, primme);
      }      
      else {
         *RprojectorX = x;
         *xKinvx = tpone;
      }
      *sizeRprojectorX = 1;
   }
   else { 
         *RprojectorX = NULL;
         *sizeRprojectorX = 0;
         *xKinvx = tpone;
   }
         
} /* setup_JD_projectors */
