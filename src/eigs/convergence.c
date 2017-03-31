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
 * File: convergence.c
 *
 * Purpose - Checks for convergence of the block vectors.
 *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "const.h"
#include "numerical.h"
#include "convergence.h"
#include "ortho.h"
#include "auxiliary_eigs.h"


static int check_practical_convergence(SCALAR *R, PRIMME_INT nLocal,
      PRIMME_INT ldR, SCALAR *evecs, int evecsSize, PRIMME_INT ldevecs,
      int left, int *iev, int numToProject, int *flags, REAL *blockNorms,
      double tol, SCALAR *rwork, size_t *rworkSize, primme_params *primme);

/*******************************************************************************
 * Subroutine check_convergence - This procedure checks the block vectors for  
 *    convergence.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * X              The Ritz vectors in the block
 * nLocal, ldX    The number of rows and the leading dimension of X
 * ldR            The leading dimension of R
 * evecs          The locked vectors
 * numLocked      The number of columns in evecs besides numOrthoConst
 * ldevecs        The leading dimension of evecs
 * left, right    Range of vectors to be checked for convergence
 * blockNorms     Residual norms of the Ritz vectors starting from left
 * hVals          The Ritz values
 * rwork          Real work array that must be of size 
 * rworkSize      The size of rwork
 * primme         Structure containing various solver parameters
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * R              The residual vectors of the Ritz vectors in the block
 *                (this routine may remove the evecs directions in some vectors)
 * flags          Array indicating which eigenvectors have converged     
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * reset          flag to reset V and W in the next restart
 ******************************************************************************/

TEMPLATE_PLEASE
int check_convergence_Sprimme(SCALAR *X, PRIMME_INT nLocal, PRIMME_INT ldX,
      SCALAR *R, PRIMME_INT ldR, SCALAR *evecs, int numLocked,
      PRIMME_INT ldevecs, int left, int right, int *flags, REAL *blockNorms,
      REAL *hVals, int *reset, double machEps, SCALAR *rwork,
      size_t *rworkSize, int *iwork, int iworkSize, primme_params *primme) {

   int i;                  /* Loop variable                                      */
   int numToProject;       /* Number of vectors with potential accuracy problem  */
   int *toProject = iwork; /* Indices from left with potential accuracy problem  */
   double tol;             /* Residual tolerance                                 */
   double attainableTol=0; /* Used in locking to check near convergence problem  */
   int isConv;             /* return of convTestFun                              */
   double targetShift;     /* target shift */

   /* -------------------------- */
   /* Return memory requirements */
   /* -------------------------- */

   if (flags == NULL) {
      CHKERR(check_practical_convergence(NULL, 0, 0, NULL, numLocked, 0, left,
            NULL, right-left, NULL, NULL, 0, NULL, rworkSize, primme), -1);
      *iwork = max(*iwork, right-left); /* for toProject */
      return 0;
   }
 
   /* Check enough space for toProject */
   assert(iworkSize >= right-left);

   targetShift = primme->numTargetShifts > 0 ?
      primme->targetShifts[min(primme->initSize, primme->numTargetShifts-1)] : 0.0;
 
   /* -------------------------------------------- */
   /* Tolerance based on our dynamic norm estimate */
   /* -------------------------------------------- */

   tol = max(machEps * max(primme->stats.estimateLargestSVal, primme->aNorm),
               primme->stats.maxConvTol);

   /* ---------------------------------------------------------------------- */
   /* If locking, set tol beyond which we need to check for accuracy problem */
   /* ---------------------------------------------------------------------- */
   if (primme->locking) {
      attainableTol = sqrt((double)(primme->numOrthoConst+numLocked))*tol;
   }

   /* ----------------------------------------------------------------- */
   /* Determine which Ritz vectors have converged < tol and flag them.  */
   /* ----------------------------------------------------------------- */

   numToProject = 0;
   for (i=left; i < right; i++) {
       
      /* Don't trust any residual norm below estimateResidualError */
      blockNorms[i-left] = max(blockNorms[i-left], primme->stats.estimateResidualError);

      /* Refine doesn't order the pairs considering closest_leq/gep. */
      /* Then ignore values so that value +-residual is completely   */
      /* outside of the desired region.                              */

      if ((primme->target == primme_closest_leq
               && hVals[i]-blockNorms[i-left] > targetShift) ||
            (primme->target == primme_closest_geq
             && hVals[i]+blockNorms[i-left] < targetShift)) {
         flags[i] = UNCONVERGED;
         continue;
      }

      CHKERR(convTestFun_Sprimme(hVals[i], X?&X[ldX*(i-left)]:NULL,
               blockNorms[i-left], &isConv, primme), -1);

      if (isConv) {
         flags[i] = CONVERGED;
      }

      /* ----------------------------------------------------------------- */
      /* If residual norm is around the bound of the error in the          */
      /* residual norm, then stop converging this value and force reset    */
      /* of V and W in the next restart.                                   */
      /* ----------------------------------------------------------------- */

      else if (blockNorms[i-left] <= primme->stats.estimateResidualError && reset) {
         flags[i] = SKIP_UNTIL_RESTART;
         *reset = 1;
      }

      /* ----------------------------------------------------------------- */
      /* If locking there may be an accuracy problem close to convergence. */
      /* Check if there is danger if R is provided. If the Ritz vector was */
      /* flagged practically converged before and R is not provided then   */
      /* consider converged still.                                         */
      /* ----------------------------------------------------------------- */

      else if (primme->locking && numLocked > 0 && blockNorms[i-left] < attainableTol ) {
         if (R) {
            toProject[numToProject++] = i-left;
         }
         else if (flags[i] != PRACTICALLY_CONVERGED) {
            flags[i] = UNCONVERGED;
         }
      }

      else {
         flags[i] = UNCONVERGED;
      }

   }

   /* --------------------------------------------------------------- */
   /* Project the TO_BE_PROJECTED residuals and check for practical   */
   /* convergence among them.                                         */
   /* --------------------------------------------------------------- */

   if (numToProject > 0) {
      CHKERR(check_practical_convergence(R, nLocal, ldR, evecs,
               primme->numOrthoConst+numLocked, ldevecs, left, toProject,
               numToProject, flags, blockNorms, tol, rwork, rworkSize, primme),
               -1);
   }

   return 0;

}

/*******************************************************************************
 * Subroutine check_practical_convergence(): for pairs whose residual norm is
 *    less than tol*sqrt(numConverged) but greater than tol, they will be
 *    flagged as converged if || R(i) - (I-QQ')R(i) || = || Q'R(i) || > tol and
 *    || (I-QQ')R(i) || < tol*tol/||R(i)||/2, where Q are the locked vectors.
 *
 *    NOTE: the routine removes evecs directions from the residual vectors,
 *          but blockNorms isn't changed.
 *
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal          The local length of the vectors in the basis
 * evecs           The locked eigenvectors
 * evecsSize       The number of locked eigenvectors
 * ldevecs         The leading dimension of evecs
 * left            Base index indicating which flags are to be recomputed
 * iev             Indices of flags to recompute
 * numToProject    Size of iev
 * blockNorms      The norms of the residual vectors starting by index 'left'
 * tol             The required convergence tolerance
 * rwork           real work array of size: 2*maxEvecsSize*primme->maxBlockSize
 * rworkSize       The size of rwork
 * primme          Structure containing various solver parameters
 *
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * R               The residual vectors
 * ldR             The leading dimension of R
 * flags           Indicates which Ritz pairs have converged
 ******************************************************************************/
static int check_practical_convergence(SCALAR *R, PRIMME_INT nLocal,
      PRIMME_INT ldR, SCALAR *evecs, int evecsSize, PRIMME_INT ldevecs,
      int left, int *iev, int numToProject, int *flags, REAL *blockNorms,
      double tol, SCALAR *rwork, size_t *rworkSize, primme_params *primme) {

   int i;
   REAL *overlaps;
   size_t rworkSize0;

   /* -------------------------- */
   /* Return memory requirements */
   /* -------------------------- */

   if (R == NULL) {
      size_t lrw=0;
      CHKERR(ortho_single_iteration_Sprimme(NULL, nLocal, evecsSize, 0,
               NULL, NULL, numToProject, 0, NULL, NULL, NULL, &lrw, primme),
               -1);
      /* for overlaps */
      lrw += numToProject;
      *rworkSize = max(*rworkSize, lrw);
      return 0;
   }

   /* ------------------------------------------------------------------ */
   /* Compute norms of the projected res and the differences from res    */ 
   /* note: ||residual - (I-QQ')residual||=||Q'*r||=||overlaps||         */
   /* ------------------------------------------------------------------ */

   /* R = (I-evecs*evecs)*R */
   /* overlaps(i) = || evecs'*R(i) || */
   /* newResiduals(i) = || (I-evecs*evecs')*R(i) || */

   overlaps = (REAL*)rwork;

   assert(*rworkSize >= (size_t)numToProject);
   rworkSize0 = *rworkSize - numToProject;
   CHKERR(ortho_single_iteration_Sprimme(evecs, nLocal, evecsSize, ldevecs,
            R, iev, numToProject, ldR, overlaps, NULL, rwork+numToProject,
            &rworkSize0, primme), -1);

   /* ------------------------------------------------------------------ */
   /* For each projected residual check whether there is an accuracy     */
   /* problem and, if so, declare it CONVERGED to lock later.            */
   /* normDiff is a lower bound to the attainable accuracy for this pair */
   /* so problems exist only if normDiff > tol. Then, we stop if Tol is  */
   /* the geometric mean of normPr and r.                                */
   /* ------------------------------------------------------------------ */

   for (i=0; i < numToProject; i++) {

      double normPr   = sqrt(max(0.0, blockNorms[iev[i]]*blockNorms[iev[i]]
                               - overlaps[i]*overlaps[i]));  /* || (I-QQ')res || */
      double normDiff = overlaps[i];                         /* || res - (I-QQ')res || */
      double blockNorm = blockNorms[iev[i]];

      /* ------------------------------------------------------------------ */
      /* NOTE: previous versions than 2.0 used the next criterion instead:  */
      /*                                                                    */
      /* if (normDiff >= tol && normPr < tol*tol/blockNorm/2)               */
      /*                                                                    */
      /* Due to rounding-off errors in computing normDiff, the test         */
      /* normDiff >= tol may fail even it is satisfied in exact arithmetic. */
      /* ------------------------------------------------------------------ */
      
      if (normPr < normDiff*normDiff/blockNorm/2) {
         if (primme->printLevel >= 5 && primme->procID == 0) {
            fprintf(primme->outputFile,
               " PRACTICALLY_CONVERGED %d norm(I-QQt)r %e bound %e\n",
                left+iev[i],normPr,tol*tol/normDiff);
                  fflush(primme->outputFile);
         }
         flags[left+iev[i]] = PRACTICALLY_CONVERGED;
      }
      else {
         flags[left+iev[i]] = UNCONVERGED;
      }
      blockNorms[iev[i]] = normPr;

   }

   return 0;
}
