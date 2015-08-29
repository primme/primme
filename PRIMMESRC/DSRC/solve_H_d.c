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
 * File: solve_H.c
 * 
 * Purpose - Solves the eigenproblem for the matrix V'*A*V.
 *
 ******************************************************************************/

#include <math.h>
#include "primme.h"
#include "solve_H_d.h"
#include "solve_H_private_d.h"
#include "numerical_d.h"

/*******************************************************************************
 * Subroutine solve_H - This procedure solves the eigenproblem for the
 *            matrix H.
 *        
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * H             The matrix V'*A*V
 * basisSize     Current size of the orthonormal basis V
 * maxBasisSize  The maximum size of the basis V
 * numLocked     Number of eigenvalues locked, to determine ordering shift.
 * lrwork        Length of the work array rwork
 * primme          Strucuture containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hVecs             The eigenvectors of H
 * hVals             The eigenvalues of H
 * largestRitzValue  Maintains the largest in absolute value Ritz value seen
 * rwork             Must be of size at least 3*maxBasisSize
 * iwork             Permutation array for evecs/evals with desired targeting 
 *                   order. hVecs/hVals are permuted in the right order.
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev was unsuccsessful
 ******************************************************************************/

int solve_H_dprimme(double *H, double *hVecs, double *hVals, 
   int basisSize, int maxBasisSize, double *largestRitzValue, int numLocked, 
   int lrwork, double *rwork, int *iwork, primme_params *primme) {

   int i, j; /* Loop variables    */
   int info; /* dsyev error value */
   int index;
   int *permu, *permw;
   double targetShift;


   /* ---------------------- */
   /* Divide the iwork space */
   /* ---------------------- */
   permu  = iwork;
   permw = permu + basisSize;

#ifdef NUM_ESSL
   int apSize, idx;
#endif


   /* ------------------------------------------------------------------- */
   /* Copy the upper triangular portion of H into hvecs.  We need to do   */
   /* this since DSYEV overwrites the input matrix with the eigenvectors. */  
   /* Note that H is maxBasisSize-by-maxBasisSize and the basisSize-by-   */
   /* basisSize submatrix of H is copied into hvecs.                      */
   /* ------------------------------------------------------------------- */

#ifdef NUM_ESSL
   idx = 0;

   if (primme->target != primme_largest) { /* smallest or any of closest_XXX */
      for (j=0; j < basisSize; j++) {
         for (i=0; i <= j; i++) {
            rwork[idx] = H[maxBasisSize*j+i];
            idx++;
         }
      }
   }
   else { /* (primme->target == primme_largest)  */
      for (j=0; j < basisSize; j++) {
         for (i=0; i <= j; i++) {
            rwork[idx] = -H[maxBasisSize*j+i];
            idx++;
         }
      }
   }
   
   apSize = basisSize*(basisSize + 1)/2;
   lrwork = lrwork - apSize;

   info = Num_dspev_dprimme(21, rwork, hVals, hVecs, basisSize, basisSize, 
      &rwork[apSize], lrwork);

   if (info != 0) {
      primme_PushErrorMessage(Primme_solve_h, Primme_num_dspev, info, __FILE__, 
         __LINE__, primme);
      return NUM_DSPEV_FAILURE;
   }

#else
   if (primme->target != primme_largest) {
      for (j=0; j < basisSize; j++) {
         for (i=0; i <= j; i++) { 
            hVecs[basisSize*j+i] = H[maxBasisSize*j+i];
         }
      }      
   }
   else { /* (primme->target == primme_largest) */
      for (j=0; j < basisSize; j++) {
         for (i=0; i <= j; i++) { 
            hVecs[basisSize*j+i] = -H[maxBasisSize*j+i];
         }
      }
   }

   Num_dsyev_dprimme("V", "U", basisSize, hVecs, basisSize, hVals, rwork, 
                lrwork, &info);

   if (info != 0) {
      primme_PushErrorMessage(Primme_solve_h, Primme_num_dsyev, info, __FILE__, 
         __LINE__, primme);
      return NUM_DSYEV_FAILURE;
   }

#endif

   /* ----------------------------------------------------------------------- */
   /* Update the largest absolute Ritz value ever seen as an estimate of ||A||
    * ----------------------------------------------------------------------- */
   *largestRitzValue = Num_fmax_primme(3, 
           *largestRitzValue, fabs(hVals[0]), fabs(hVals[basisSize-1]));

   /* ---------------------------------------------------------------------- */
   /* ORDER the eigenvalues and their eigenvectors according to the desired  */
   /* target:  smallest/Largest or interior closest abs/leq/geq to a shift   */
   /* ---------------------------------------------------------------------- */

   if (primme->target == primme_smallest) 
      return 0;

   if (primme->target == primme_largest) {
      for (i = 0; i < basisSize; i++) {
         hVals[i] = -hVals[i];
      }
   }
   else { 
      /* ---------------------------------------------------------------- */
      /* Select the interior shift. Use the first unlocked shift, and not */
      /* higher ones, even if some eigenpairs in the basis are converged. */
      /* Then order the ritz values based on the closeness to the shift   */
      /* from the left, from right, or in absolute value terms            */
      /* ---------------------------------------------------------------- */

      targetShift = 
        primme->targetShifts[min(primme->numTargetShifts-1, numLocked)];

      if (primme->target == primme_closest_geq) {
   
         /* ---------------------------------------------------------------- */
         /* find hVal closest to the right of targetShift, i.e., closest_geq */
         /* ---------------------------------------------------------------- */
         for (j=0;j<basisSize;j++) 
              if (hVals[j]>=targetShift) break;
           
         /* figure out this ordering */
         index = 0;
   
         for (i=j; i<basisSize; i++) {
            permu[index++]=i;
         }
         for (i=0; i<j; i++) {
            permu[index++]=i;
         }
      }
      else if (primme->target == primme_closest_leq) {
         /* ---------------------------------------------------------------- */
         /* find hVal closest_leq to targetShift                             */
         /* ---------------------------------------------------------------- */
         for (j=basisSize-1; j>=0 ;j--) 
             if (hVals[j]<=targetShift) break;
           
         /* figure out this ordering */
         index = 0;
   
         for (i=j; i>=0; i--) {
            permu[index++]=i;
         }
         for (i=basisSize-1; i>j; i--) {
            permu[index++]=i;
         }
      }
      else if (primme->target == primme_closest_abs) {

         /* ---------------------------------------------------------------- */
         /* find hVal closest but geq than targetShift                       */
         /* ---------------------------------------------------------------- */
         for (j=0;j<basisSize;j++) 
             if (hVals[j]>=targetShift) break;

         i = j-1;
         index = 0;
         while (i>=0 && j<basisSize) {
            if (fabs(hVals[i]-targetShift) < fabs(hVals[j]-targetShift)) 
               permu[index++] = i--;
            else 
               permu[index++] = j++;
         }
         if (i<0) {
            for (i=j;i<basisSize;i++) 
                    permu[index++] = i;
         }
         else if (j>=basisSize) {
            for (j=i;j>=0;j--)
                    permu[index++] = j;
         }
      }

      /* ---------------------------------------------------------------- */
      /* Reorder hVals and hVecs according to the permutation             */
      /* ---------------------------------------------------------------- */
      for (i=0;i<basisSize;i++) 
          permw[i] = permu[i];
      permute_evecs_dprimme(hVals, permu, rwork, basisSize, 1);
      permute_evecs_dprimme(hVecs, permw, rwork, basisSize, basisSize);
   }


   return 0;   
}

/******************************************************************************
 * Subroutine permute_evecs- This routine permutes a set of vectors according
 *            to a permutation array perm. It is supposed to be called on 
 *            the eigenvectors after the eigenvalues have been sorted, so that 
 *            the vectors are in the same order as the sorted eigenvalues.
 *
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * rwork  Used as temorary space to swap the vectors
 * nev    Number of eigenvectors/eigenvalues
 * nLocal Number of rows of each vector this process stores
 *
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * perm  The array indicating what order the eigenvectors must be permuted to.
 *       perm[i] indicates that the ith eigenvector in the sorted ordering
 *       should be the perm[i]-th vector from the original ordering.
 * evecs  The eigenvectors to be put in sorted order.
 *
 ******************************************************************************/
     
void permute_evecs_dprimme(double *evecs, int *perm, double *rwork, int nev, 
   int nLocal) {

   int currentIndex;     /* Index of eigenvector in sorted order              */
   int sourceIndex;      /* Position of out-of-order vector in original order */
   int destinationIndex; /* Position of out-of-order vector in sorted order   */
   int tempIndex;        /* Used to swap                                      */
   
   currentIndex = 0;

   /* Continue until all eigenvectors are in the sorted order */

   while (1) {

      /* Find a vector that does not belong in its original position */
      while ((currentIndex < nev) && (perm[currentIndex] == currentIndex)) {
         currentIndex++;
      }

      /* Return if they are in the sorted order */
      if (currentIndex >= nev) {
         return;
      }

      /* Copy the vector to a buffer for swapping */
      Num_dcopy_primme(nLocal, &evecs[currentIndex*nLocal], 1, rwork, 1);

      destinationIndex = currentIndex;
      /* Copy vector perm[destinationIndex] into position destinationIndex */

      while (perm[destinationIndex] != currentIndex) {

         sourceIndex = perm[destinationIndex];
         Num_dcopy_primme(nLocal, &evecs[sourceIndex*nLocal], 1, 
            &evecs[destinationIndex*nLocal], 1);
         tempIndex = perm[destinationIndex];
         perm[destinationIndex] = destinationIndex;
         destinationIndex = tempIndex;
      }

      /* Copy the vector from the buffer to where it belongs */
      Num_dcopy_primme(nLocal, rwork, 1, &evecs[destinationIndex*nLocal], 1);
      perm[destinationIndex] = destinationIndex;

      currentIndex++;
   }

  /***************************************************************************/
} /* end of permute_evecs
   ***************************************************************************/

