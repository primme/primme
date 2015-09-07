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
#include "primme.h"
#include "const.h"
#include "wtime.h"
#include "convergence_z.h"
#include "convergence_private_z.h"
#include "numerical_z.h"

/*******************************************************************************
 * Subroutine check_convergence - This procedure checks the block vectors for  
 *    convergence.  For each of the Ritz vectors that has converged, an
 *    unconverged one will be chosen to replace it. 
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V              The orthonormal basis
 * W              A*V
 * hVecs          The eigenvectors of V'*A*V
 * hVals          The Ritz values
 * basisSize      Size of the basis V
 * numReqEvals    Total number of eigenpairs the user wants computed
 * numConverged   Number of vectors that have converged (not updated here)
 * numLocked      The number of vectors currently locked (if locking)
 * evecs          
 * tol            Tolerance used to determine convergence of residual norms
 * maxConvTol     The max residual norm > tol for any locked eigenpair 
 *                that has been determined to have an accuracy problem 
 * aNormEstimate  If primme->aNorm<=0, use tol*aNormEstimate (=largestRitzValue)
 * rwork          Real work array that must be of size 
 *                MAX(2*maxEvecsSize*primme->maxBlockSize, primme->maxBlockSize+
 *                    primme->maxBasisSize*primme->maxBlockSize);
 * primme           Structure containing various solver parameters
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * flags         Array indicating which eigenvectors have converged     
 * iev           indicates which eigenvalue each block vector corresponds to
 * ievMax        index of next eigenvalue to be targeted by the block 
 * blockNorms    Residual norms of the Ritz vectors being computed during the
 *               current iteration
 * blockSize     Dimension of the block
 ******************************************************************************/

int check_convergence_zprimme(Complex_Z *V, Complex_Z *W, Complex_Z *hVecs, 
   double *hVals, int *flags, int basisSize, int *iev, int *ievMax, 
   double *blockNorms, int *blockSize, int numConverged, int numLocked, 
   Complex_Z *evecs, double tol, double maxConvTol, double aNormEstimate, 
   Complex_Z *rwork, primme_params *primme) {

   int i;             /* Loop variable                                        */
   int left, right;   /* Range of block vectors to be checked for convergence */
   int start;         /* starting index in block of converged/tobeProject vecs*/
   int numVacancies;  /* Number of vacant positions between left and right    */
   int recentlyConverged; /* The number of Ritz values declared converged     */
                          /* since the last iteration                         */
   int numToProject;      /* Number of vectors with potential accuracy problem*/
   double attainableTol;  /* Used in locking to check near convergence problem*/

   /* -------------------------------------------- */
   /* Tolerance based on our dynamic norm estimate */
   /* -------------------------------------------- */

   if (primme->aNorm <= 0.0L) {
      tol = tol * aNormEstimate;
   }

   /* ---------------------------------------------------------------------- */
   /* If locking, set tol beyond which we need to check for accuracy problem */
   /* ---------------------------------------------------------------------- */
   if (primme->locking) {
      attainableTol = sqrt(primme->numOrthoConst+numLocked)*maxConvTol;
   }   

      
   /* --------------------------------------------------------------- */
   /* Compute each Ritz vector and its corresponding residual vector. */
   /* The Ritz vector and residual are stored temporarily in V and W  */
   /* respectively.  For each Ritz vector, determine if it has        */
   /* converged.  If it has, try to replace it with one that hasn't.  */
   /* --------------------------------------------------------------- */

   recentlyConverged = 0;
   left = 0;
   right = *blockSize - 1;
   numVacancies = 1;

   while (numVacancies > 0 && 
          (numConverged + recentlyConverged) < primme->numEvals) {

      /* Consider the newly added vectors in the block and reset counters */
      numVacancies = 0;
      numToProject = 0;

      /* Copy needed hvecs into the front of the work array. */

      for (i=left; i <= right; i++) {
         Num_zcopy_zprimme(basisSize, &hVecs[basisSize*iev[i]], 1, 
            &rwork[basisSize*(i-left)], 1);
      }
            
      /* ----------------------------------------------------------------- */
      /* Compute the Ritz vectors, residuals, and norms for the next       */
      /* blockSize unconverged Ritz vectors.  The Ritz vectors will be     */
      /* placed from V(0,lft) to V(0,rgt) and the residual vectors from    */
      /* W(0,lft) to W(0,rgt).                                             */
      /* ----------------------------------------------------------------- */
      /* rwork must be maxBasisSize*maxBlockSize + maxBlockSize in size,   */
      /* maxBasisSize*maxBlockSize holds selected hVecs to facilitate      */
      /* blocking, and maxBlockSize to hold the residual norms             */
      /* ----------------------------------------------------------------- */

      compute_resnorms(V, W, rwork, hVals, basisSize, blockNorms,
         iev, left, right, &rwork[basisSize*(right-left+1)], primme);

      print_residuals(hVals, blockNorms, numConverged, numLocked, iev, 
         left, right, primme);

      /* ----------------------------------------------------------------- */
      /* Determine which Ritz vectors have converged < tol and flag them.  */
      /* ----------------------------------------------------------------- */

      for (i=left; i <= right; i++) {
       
         /* ------------------------------------*/
         /* If the vector is converged, flag it */
         /* ------------------------------------*/
         if (blockNorms[i] < tol) {
            flags[iev[i]] = CONVERGED;
            numVacancies++;

            if ((!primme->locking && iev[i] < primme->numEvals) || 
               (primme->locking && ((numLocked + iev[i]) < primme->numEvals))) {

               recentlyConverged++;

               if (!primme->locking && primme->procID == 0 && 
                   primme->printLevel >= 2) { fprintf(primme->outputFile, 
                  "#Converged %d eval[ %d ]= %e norm %e Mvecs %d Time %g\n",
                  numConverged+recentlyConverged, iev[i], hVals[iev[i]], 
                  blockNorms[i], primme->stats.numMatvecs,primme_wTimer(0));
                  fflush(primme->outputFile);
               } /* printf */
            } /*if */
         } /*if converged */
         /* ---------------------------------------------------------------- */
         /* If locking there may be an accuracy problem close to convergence */
         /* Check if there is danger and set these Ritz vecs for projection  */
         /* ---------------------------------------------------------------- */
         else if (primme->locking && numLocked > 0 &&
                  blockNorms[i] < attainableTol ) {

            flags[iev[i]] = TO_BE_PROJECTED;
            numToProject++;
         }

      } /* for */

      /* ---------------------------------------------------------------- */
      /* If some of the Ritz vectors in the block have converged, or need */
      /* to be projected against evecs, move those flagged Ritz vectors   */
      /* and residuals towards the end of the block [left,right]. Also    */
      /* swap iev, and blockNorms for the targeted block.                 */
      /* ---------------------------------------------------------------- */
      if (numVacancies > 0 || numToProject > 0) {

         swap_UnconvVecs(V, W, primme->nLocal, basisSize, iev, flags, 
            blockNorms, primme->numOrthoConst + numLocked, *blockSize, left);
      }
      /* --------------------------------------------------------------- */
      /* Project the TO_BE_PROJECTED residuals and check for practical   */
      /* convergence among them. Those practically converged evecs are   */
      /* swapped just before the converged ones at the end of the block. */
      /* numVacancies and recentlyConverged are also updated             */
      /* --------------------------------------------------------------- */
      if (numToProject > 0) {

         start = *blockSize - numVacancies - numToProject;

         check_practical_convergence(V, W, evecs, numLocked, basisSize, 
            *blockSize, start, numToProject, iev, flags, blockNorms, tol, 
            &recentlyConverged, &numVacancies, rwork, primme);
      }

      /* ---------------------------------------------------------------- */
      /* Replace the vacancies, with as many unconverged vectors beyond   */
      /* ievMax as possible. If not enough are available reduce blockSize */
      /* ---------------------------------------------------------------- */

      if (numVacancies > 0) {
         replace_vectors(iev, flags, *blockSize, basisSize, numVacancies, 
                         &left, &right, ievMax); 
         numVacancies = right - left + 1;
         *blockSize = left + numVacancies;
      }

   } /* while there are vacancies */

   return recentlyConverged;
}


/*******************************************************************************
 * Subroutine compute_resnorms - This routine computes the Ritz vectors, the
 *    corresponding residual vectors, and the residual norms. The Ritz vectors 
 *    are stored in V(0,nv+left) through V(0,nv+right), the residual vectors 
 *    are stored in W(0,nv+left) through W(0,nv+right), and the residual norms 
 *    are stored in the blockNorms array.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * hVecs        The eigenvectors of V'*A*V
 * hVals        The eigenvalues of V'*A*V
 * nLocal       Number of rows of V assigned to the node
 * basisSize    Number of vectors in the basis V
 * iev          indicates which eigenvalue each block vector correspondso to
 * left, right  Ritz vectors left through right will be computed
 * primme       Structure containing various solver parameters
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * V           The orthonormal basis.
 * W           A*V
 * blockNorms  Norms of the residual vectors 
 * rwork       Must be at least maxBlockSize in length
 ******************************************************************************/

static void compute_resnorms(Complex_Z *V, Complex_Z *W, Complex_Z *hVecs, 
   double *hVals, int basisSize, double *blockNorms, int *iev, int left, 
   int right, void *rwork, primme_params *primme) {

   int i;            /* Loop variable                             */
   int numResiduals; /* Number of residual vectors to be computed */
   double *dwork = (double *) rwork;  /* pointer casting rwork to double */
   Complex_Z ztmp;     /* temp var holding shift                    */
   Complex_Z tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00};       /* constants */

   numResiduals = right - left + 1;

   /* We want to compute residuals r = Ax-hVal*x for the Ritz vectors x */
   /* Eqivalently, r = A*V*hVec - hval*V*hVec = W*hVec - hVal*V*hVec.   */

   /* Compute the Ritz vectors */

   Num_gemm_zprimme("N", "N", primme->nLocal, numResiduals, basisSize, 
      tpone, V, primme->nLocal, hVecs, basisSize, tzero,
      &V[primme->nLocal*(basisSize+left)], primme->nLocal);

   /* Compute W*hVecs */

   Num_gemm_zprimme("N", "N", primme->nLocal, numResiduals, basisSize, 
      tpone, W, primme->nLocal, hVecs, basisSize, tzero,
      &W[primme->nLocal*(basisSize+left)], primme->nLocal);

   /* Compute the residuals */

   for (i=left; i <= right; i++) {
      {ztmp.r = -hVals[iev[i]]; ztmp.i = 0.0L;}
      Num_axpy_zprimme(primme->nLocal, ztmp, &V[primme->nLocal*(basisSize+i)],
       1, &W[primme->nLocal*(basisSize+i)], 1);
   }

   /* Compute the residual norms */

   for (i=left; i <= right; i++) {
      ztmp = Num_dot_zprimme(primme->nLocal, &W[primme->nLocal*(basisSize+i)],
         1, &W[primme->nLocal*(basisSize+i)] , 1);
      dwork[i] = ztmp.r;
   }
   
   (*primme->globalSumDouble)(&dwork[left], &blockNorms[left], &numResiduals,
                              primme);

   for (i=left; i <= right; i++) {
      blockNorms[i] = sqrt(blockNorms[i]);
   }

}


/*******************************************************************************
 * Subroutine print_residuals - This function displays the residual norms of 
 *    each Ritz vector computed at this iteration.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * ritzValues     The Ritz values
 * blockNorms     Residual norms of the corresponding Ritz vectors
 * numConverged   Number of already converged epairs
 * numLocked      Number of locked (converged) epairs
 * iev            indicates which eigenvalue each block vector corresponds to
 * left, right    Ritz values left through right will be displayed
 * numOuterIterations  Number of JD iterations performed
 * numMatvecs     Number of matrix-vector multiplications performed
 * primme         Structure containing various solver parameters
 ******************************************************************************/

static void print_residuals(double *ritzValues, double *blockNorms,
   int numConverged, int numLocked, int *iev, int left, int right, 
   primme_params *primme) {

   int i;  /* Loop variable */
   int found;  /* Loop variable */

   if (primme->printLevel >= 3 && primme->procID == 0) {

      if (primme->locking) 
         found = numLocked;
      else 
         found = numConverged;

      for (i=left; i <= right; i++) {
         fprintf(primme->outputFile, 
            "OUT %d conv %d blk %d MV %d Sec %E EV %13E |r| %.3E\n",
         primme->stats.numOuterIterations, found, i, primme->stats.numMatvecs,
         primme_wTimer(0), ritzValues[iev[i]], blockNorms[i]);
      }

      fflush(primme->outputFile);
   }


}



/*******************************************************************************
 * Subroutine swap_UnconvVecs - This procedure copies unconverged Ritz vectors 
 *    towards the beginning of the block starting at V(:,basisSize+left).
 *    Flagged vectors are copied to the right end of the block, ending at 
 *    basisSize + blockSize-1. The corresponding residuals in W are also 
 *    copied towards the front and end correspondingly of the block 
 *    W(:,basisSize+left:basisSize+blockSize-1). This ensures that the 
 *    unconverged Ritz vectors as well as the residuals are contiguous.  Also, 
 *    new vectors selected to replace the converged ones will be contiguous too.
 * 
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * nLocal        Number of rows of V assigned to the node
 * basisSize     Number of vectors in the basis V
 * flags         Indicates which of the Ritz vectors have converged/flagged
 * blockSize     Number of block vectors
 * left, right   Indicies of the block vectors to be swapped  
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * V             The orthonormal basis
 * W             A*V
 * iev           the index of each block vector
 * blockNorms    the norms for each vector in the block
 ******************************************************************************/
     
static void swap_UnconvVecs(Complex_Z *V, Complex_Z *W, int nLocal, 
   int basisSize, int *iev, int *flags, double *blockNorms, int dimEvecs, 
   int blockSize, int left) {

   int right; /* holds the right swapping position */
   int temp;  /* used to swap integers */
   double dtemp; /* used to swap doubles */

   /* Search from left to right within the block looking for flagged     */
   /* Ritz vectors.  If a flagged one is found, find also an unconverged */
   /* Ritz vector from the right side of the block. If the flagged was   */
   /* converged, replace it with the unconverged, otherwise if it was    */
   /* was flagged TO_BE_PROJECTED, swap it with the unconverged vector.  */

   while (left < blockSize) {

      /* If block vector left is unconverged, move left one vector. */ 

      if (flags[iev[left]] == UNCONVERGED) {
         left++;
      }
      else {

         /* If block vector left is converged, find an unconverged */
         /* vector somewhere between left+1 and blockSize.         */

         right = left + 1;
   
         /* If the end of the block has been reached, there are */
         /* no more unconverged vectors left.  If not, keep     */
         /* searching right for another unconverged vector to   */
         /* perform the replacement.                            */

         while (right < blockSize && flags[iev[right]] != UNCONVERGED) {
            right++;
         }

         /* No unconverged block vector could be found */

         if (right == blockSize) {
            return;
         }

         /* An unconverged Ritz vector was found and should */
         /* replace or be swapped with block vector left.   */

         if (flags[iev[left]] != TO_BE_PROJECTED) { 
                /* replace */
            Num_zcopy_zprimme(nLocal, &V[nLocal*(basisSize+right)], 1,
               &V[nLocal*(basisSize+left)], 1);
            Num_zcopy_zprimme(nLocal, &W[nLocal*(basisSize+right)], 1,
               &W[nLocal*(basisSize+left)], 1);
            temp = iev[left];
            iev[left] = iev[right];
            iev[right] = temp;
            blockNorms[left] = blockNorms[right];
         }
         else { /* swap */
            Num_swap_zprimme(nLocal, &V[nLocal*(basisSize+left)], 1, 
                              &V[nLocal*(basisSize+right)], 1);
            Num_swap_zprimme(nLocal, &W[nLocal*(basisSize+left)], 1, 
                              &W[nLocal*(basisSize+right)], 1);
            temp = iev[left];
            iev[left] = iev[right];
            iev[right] = temp;
            dtemp = blockNorms[left];
            blockNorms[left] = blockNorms[right];
            blockNorms[right] = dtemp;
         } /* end of swaps */

         left++;

      }  /* looking for replacement */

   } /* while left < blockSize */
}


/*******************************************************************************
 * Subroutine replace_vectors - This routine determines which eigenvectors
 *            are to be targeted in the next iteration.  If one is marked
 *            as converged, another is targeted in its place if one is
 *            available.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * flags           Indicates which Ritz pairs have converged
 * blockSize       The number of block vectors
 * basisSize       Number of vectors in the basis
 * numVacancies    Number of Ritz values between left and right that were
 *                 declared converged  
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * iev             Indicates which Ritz value each block vector corresponds to
 * left, right     Indices indicating which vectors are to be replaced
 * ievMax          Index of the next Ritz value to be targeted by the block
 ******************************************************************************/

static void replace_vectors(int *iev, int *flags, int blockSize, int basisSize,
   int numVacancies, int *left, int *right, int *ievMax) {

   int i;  /* Counter variable */

   /* The unconverged Ritz vectors were copied towards the front of the    */
   /* block by the swap_UnconvVecs routine, so there are now vacant spaces */
   /* towards the end of the block that must be replaced.                  */

   *left = blockSize - numVacancies;
   
   /* Target unconverged Ritz vectors.  Replace the convergedCount      */
   /* converged Ritz vectors with unconverged vectors.  If there are    */
   /* fewer than convergedCount unconverged vectors available, then     */
   /* replace as many as possible.                                      */

   i = *left;

   while ((i < blockSize) && (*ievMax < basisSize)) {

      /* If the Ritz vector is unconverged, target it */

      if (flags[*ievMax] == UNCONVERGED) {
         iev[i] = *ievMax;
         i++;
      }

      *ievMax = *ievMax + 1;

   }

   *right = i - 1;

}
/*******************************************************************************
 * Subroutine check_practical_convergence()
 *
 *       This function is called after swaping has pushed converged (C)
 *       and to be projected (P) evecs at the end of blockSize.
 *       Makes C, P contiguous, projects P, checks and finds practically 
 *       converged evecs (p), and makes p contiguous to C, updating
 *       flags[], blockNorms[], and numVacancies. Eg:
 *
 *          blockSize 
 *       |------CPCCPPCP|  -> |------PPPPCCCC|  ->  |------PpPpCCCC| ->
 *       |------PPppCCCC|, numVacancies = 6
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                    
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * evecs           The locked eigenvectors
 * numLocked       The number of locked eigenvectors
 * basisSize       Number of vectors in the basis
 * blockSize       The number of block vectors
 * start           Starting index in V,W of vectors converged or to be projected
 * numToProject    The number of vectors to project. 
 * tol             The required convergence tolerance
 * rwork           real work array of size: 2*maxEvecsSize*primme->maxBlockSize
 * primme          Structure containing various solver parameters
 *
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * V               The basis vectors
 * W               A*V
 * iev             Indicates which Ritz value each block vector corresponds to
 * flags           Indicates which Ritz pairs have converged
 * blockNorms      The norms of the block vectors to be targeted
 * recentlyConverged Number of converged vectors in the whole basis V
 *                   = converged+practicallyConverged 
 * numVacancies    Number of Ritz values between left and right that were
 *                 declared converged or practicallyConverged. [left, right]
 *                 can be smaller than the whole V. See while-loop in 
 *                 check_convergence.
 * left, right     Indices indicating which vectors are to be replaced
 * ievMax          Index of the next Ritz value to be targeted by the block
 ******************************************************************************/
static void check_practical_convergence(Complex_Z *V, Complex_Z *W, 
   Complex_Z *evecs, int numLocked, int basisSize, int blockSize, int start, 
   int numToProject, int *iev, int *flags, double *blockNorms, double tol, 
   int *recentlyConverged, int *numVacancies, Complex_Z *rwork, 
   primme_params *primme) {

   int i, n, dimEvecs;
   int count; 
   double normPr; 
   double normDiff;
   Complex_Z *overlaps;  /* Pointer in rwork to keep Q'*r */
   Complex_Z tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00}, tmone = {-1.0e+00,+0.0e00}; /* constants */

   /* convenience variables */
   n        = primme->nLocal;
   dimEvecs = primme->numOrthoConst + numLocked;

   /* Subdivide rwork */
   overlaps = rwork + dimEvecs*primme->maxBlockSize;

   /* --------------------------------------------------------------- */
   /* Reset any TO_BE_PROJECTED flags back to UNCONVERGED, and  swap  */
   /* again CONVERGED flags toward the end of the block. Now swapping */
   /* occurs in the block basisSize+[start:blockSize]:                */
   /*        [ . . . . . . . P P P P C C C C ]     blockSize          */
   /*                   start^       ^start+numToProject              */
   /* --------------------------------------------------------------- */

   for (i=start; i < blockSize; i++)
      if (flags[iev[i]] == TO_BE_PROJECTED)
         flags[iev[i]] = UNCONVERGED;

   if (*numVacancies > 0)
      swap_UnconvVecs(V, W, primme->nLocal, basisSize, iev, flags, 
                   blockNorms, dimEvecs, blockSize, start);

   /* ------------------------------------------------------------------ */
   /* Project the numToProject residuals agaist (I-evecs*evecs')         */
   /* ------------------------------------------------------------------ */

   /* overlaps = evecs'*residuals */

   Num_gemm_zprimme("C", "N", dimEvecs, numToProject, n, tpone, evecs, n, 
                  &W[(basisSize+start)*n], n, tzero, rwork, dimEvecs);

   /* In Complex, the size of the array to globalSum is twice as large */
   count = 2*(dimEvecs*numToProject);
   (*primme->globalSumDouble)(rwork, overlaps, &count, primme);

   /* residuals = residuals - evecs*overlaps */

   Num_gemm_zprimme("N", "N", n, numToProject, dimEvecs, tmone, evecs, n, 
                  overlaps, dimEvecs, tpone, &W[(basisSize+start)*n], n);

   /* ------------------------------------------------------------------ */
   /* Compute norms^2 of the projected res and the differences from res  */ 
   /* note: ||residual - (I-QQ')residual||=||Q'*r||=||overlaps||         */
   /* ------------------------------------------------------------------ */

   for (i=0; i < numToProject; i++) {
      /* || res - (I-QQ')res || */
      rwork[i] = Num_dot_zprimme(dimEvecs, &overlaps[dimEvecs*i], 1, 
                                &overlaps[dimEvecs*i], 1);
      /* || (I-QQ')res || */
      rwork[i+numToProject] = Num_dot_zprimme(n, &W[(basisSize+start+i)*n], 1,
                                &W[(basisSize+start+i)*n], 1);
   }
   /* global sum ||overlaps|| and ||(I-QQ')r|| */
   /* In Complex, the size of the array to globalSum is twice as large */
   count = 2*(2*numToProject);
   (*primme->globalSumDouble)(rwork, &rwork[count], &count, primme);

   /* ------------------------------------------------------------------ */
   /* For each projected residual check whether there is an accuracy     */
   /* problem and, if so, declare it PRACTICALLY_CONVERGED to lock later.*/
   /* normDiff is a lower bound to the attainable accuracy for this pair */
   /* so problems exist only if normDiff > tol. Then, we stop if Tol is  */
   /* the geometric mean of normPr and r.                                */
   /* ------------------------------------------------------------------ */

   for (i=start; i < start+numToProject; i++) {

      normDiff = sqrt(rwork[i-start].r);
      normPr   = sqrt(rwork[i-start+numToProject].r);

      /* printf(" R Pr |R-Pr| %e %e %e \n", blockNorms[i],normPr,normDiff); */

      if (normDiff >= tol && normPr < tol*tol/blockNorms[i]/2) {
         if (primme->printLevel >= 5 && primme->procID == 0) {
            fprintf(primme->outputFile, 
               " PRACTICALLY_CONVERGED %d norm(I-QQt)r %e bound %e\n",
                iev[i],normPr,tol*tol/normDiff);
                  fflush(primme->outputFile);
         }
         flags[iev[i]] = PRACTICALLY_CONVERGED;
         (*numVacancies)++;
         if (numLocked + iev[i] < primme->numEvals){
            recentlyConverged++;
         }
      } /* if practically converged */

   } /* for each projected residual */

   /* ------------------------------------------------------------------ */
   /* Finally swap all practically converged toward the end of the block */
   /* to be contiguous with the converged ones. Swap all relevant arrays */
   /* ------------------------------------------------------------------ */

   start = blockSize - *numVacancies;

   swap_UnconvVecs(V, W, primme->nLocal, basisSize, iev, flags, blockNorms, 
                         dimEvecs, blockSize, start);

}
