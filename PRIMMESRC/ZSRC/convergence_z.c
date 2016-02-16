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
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boion, MA  02110-1301  USA
 *
 *******************************************************************************
 * File: convergence.c
 *
 * Purpose - Checks for convergence of the block vectors.
 *
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "primme.h"
#include "const.h"
#include "wtime.h"
#include "convergence_z.h"
#include "convergence_private_z.h"
#include "ortho_z.h"
#include "numerical_z.h"

/* Extra estates for flags */
#define PRACTICALLY_CONVERGED  2


/*******************************************************************************
 * Subroutine check_convergence - This procedure checks the block vectors for  
 *    convergence.  For each of the Ritz vectors that has converged, an
 *    unconverged one will be chosen to replace it. 
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
 ******************************************************************************/

int check_convergence_zprimme(Complex_Z *X, int nLocal, int ldX, Complex_Z *R,
   int ldR, Complex_Z *evecs, int numLocked, int ldevecs, int left, int right,
   int *flags, double *blockNorms, double *hVals, double machEps, Complex_Z *rwork,
   int rworkSize, int *iwork, primme_params *primme) {

   int i;                  /* Loop variable                                      */
   int numToProject;       /* Number of vectors with potential accuracy problem  */
   int *toProject = iwork; /* Indices from left with potential accuracy problem  */
   double tol;             /* Residual tolerance                                 */
   double attainableTol=0; /* Used in locking to check near convergence problem  */
   int isConv;             /* return of convTestFun                              */
   int ret=0;

   /* -------------------------- */
   /* Return memory requirements */
   /* -------------------------- */

   if (X == NULL) {
      return R ? check_practical_convergence(NULL, 0, 0, NULL, numLocked, 0, left,
         NULL, right-left, NULL, NULL, 0, NULL, 0, primme) : 0;
   }
   
   /* -------------------------------------------- */
   /* Tolerance based on our dynamic norm estimate */
   /* -------------------------------------------- */

   tol = max(machEps * max(primme->stats.estimateLargestSVal, primme->aNorm),
               primme->stats.maxConvTol);

   /* ---------------------------------------------------------------------- */
   /* If locking, set tol beyond which we need to check for accuracy problem */
   /* ---------------------------------------------------------------------- */
   if (primme->locking) {
      attainableTol = sqrt(primme->numOrthoConst+numLocked)*tol;
   }

   /* ----------------------------------------------------------------- */
   /* Determine which Ritz vectors have converged < tol and flag them.  */
   /* ----------------------------------------------------------------- */

   numToProject = 0;
   for (i=left; i < right; i++) {
       
      primme->convTestFun(&hVals[i], &X[ldX*(i-left)], &blockNorms[i-left],
            &isConv, primme);

      if (isConv) {
         flags[i] = CONVERGED;
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
   /* convergence among them. Those practically converged evecs are   */
   /* swapped just before the converged ones at the end of the block. */
   /* numVacancies and recentlyConverged are also updated             */
   /* --------------------------------------------------------------- */
   if (numToProject > 0) {
      ret = check_practical_convergence(R, nLocal, ldR, evecs,
         primme->numOrthoConst+numLocked, ldevecs, left, toProject,
         numToProject, flags, blockNorms, tol, rwork, rworkSize, primme);
   }

   return ret;

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
static int check_practical_convergence(Complex_Z *R, int nLocal, int ldR,
   Complex_Z *evecs, int evecsSize, int ldevecs, int left, int *iev,
   int numToProject, int *flags, double *blockNorms, double tol, Complex_Z *rwork,
   int rworkSize, primme_params *primme) {

   int i, ret;
   double *overlaps;

   /* -------------------------- */
   /* Return memory requirements */
   /* -------------------------- */

   if (R == NULL) {
      return numToProject + ortho_single_iteration_zprimme(NULL, nLocal,
         evecsSize, 0, NULL, NULL, numToProject, 0, NULL, NULL, NULL, 0, primme);
   }

   /* ------------------------------------------------------------------ */
   /* Compute norms of the projected res and the differences from res    */ 
   /* note: ||residual - (I-QQ')residual||=||Q'*r||=||overlaps||         */
   /* ------------------------------------------------------------------ */

   /* R = (I-evecs*evecs)*R */
   /* overlaps(i) = || evecs'*R(i) || */
   /* newResiduals(i) = || (I-evecs*evecs')*R(i) || */

   overlaps = (double*)rwork;

   ret = ortho_single_iteration_zprimme(evecs, nLocal, evecsSize, ldevecs,
      R, iev, numToProject, ldR, overlaps, NULL, rwork+numToProject,
      rworkSize-numToProject, primme);
   if (ret != 0) return ret;

   /* ------------------------------------------------------------------ */
   /* For each projected residual check whether there is an accuracy     */
   /* problem and, if so, declare it CONVERGED to lock later.            */
   /* normDiff is a lower bound to the attainable accuracy for this pair */
   /* so problems exist only if normDiff > tol. Then, we stop if Tol is  */
   /* the geometric mean of normPr and r.                                */
   /* ------------------------------------------------------------------ */

   for (i=0; i < numToProject; i++) {

      double normPr   = sqrt(blockNorms[iev[i]]*blockNorms[iev[i]]
                               - overlaps[i]*overlaps[i]);   /* || (I-QQ')res || */
      double normDiff = overlaps[i];                         /* || res - (I-QQ')res || */
      double blockNorm = blockNorms[iev[i]];

      assert(blockNorms[iev[i]] >= overlaps[i]);

      if (/*normDiff >= tol &&*/ normPr < normDiff*normDiff/blockNorm/2) {
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

   }

   return 0;
}
