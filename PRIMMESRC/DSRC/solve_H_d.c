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
#include <assert.h>
#include "primme.h"
#include "solve_H_d.h"
#include "solve_H_private_d.h"
#include "numerical_d.h"
#include "ortho_d.h"

/*******************************************************************************
 * Subroutine solve_H - This procedure solves the project problem and return
 *       the projected vectors (hVecs) and values (hVals) in the order according
 *       to primme.target.
 *        
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * H              The matrix V'*A*V
 * basisSize      The dimension of H, R, hU
 * ldH            The leading dimension of H
 * R              The factor R for the QR decomposition of (A - target*I)*V
 * ldR            The leading dimension of R
 * numConverged   Number of eigenvalues converged to determine ordering shift
 * lrwork         Length of the work array rwork
 * primme         Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hU                The left singular vectors of R
 * hVecs             The eigenvectors of H or the right singular vectors
 * hVals             The Ritz values
 * hSVals            The singular values of R
 * rwork             Workspace
 * iwork             Workspace in integers
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev was unsuccessful
 ******************************************************************************/

int solve_H_dprimme(double *H, int basisSize, int ldH, double *R, int ldR,
   double *hU, int ldhU, double *hVecs, int ldhVecs, double *hVals, double *hSVals,
   int numConverged, int lrwork, double *rwork, int *iwork, primme_params *primme) {

   int i, ret;

   switch (primme->projectionParams.projection) {
   case primme_proj_RR:
      ret = solve_H_RR_dprimme(H, ldH, hVecs, ldhVecs, hVals, basisSize,
            numConverged, lrwork, rwork, iwork, primme);
      break;

   case primme_proj_harmonic:
      assert(0);
      break;

   case primme_proj_refined:
      ret = solve_H_Ref_dprimme(H, ldH, hVecs, ldhVecs, hU, ldhU, hSVals, 
            R, ldR, hVals, basisSize, lrwork, rwork, primme);
      break;

   default:
      assert(0);
   }

   /* Return memory requirements */

   if (H == NULL) {
      return ret;
   }

   if (ret != 0) return ret;

   /* -------------------------------------------------------- */
   /* Update the leftmost and rightmost Ritz values ever seen  */
   /* -------------------------------------------------------- */
   for (i=0; i<basisSize; i++) {
      primme->stats.estimateMinEVal = min(primme->stats.estimateMinEVal,
            hVals[i]); 
      primme->stats.estimateMaxEVal = max(primme->stats.estimateMaxEVal,
            hVals[i]); 
   }
   primme->stats.estimateLargestSVal = max(fabs(primme->stats.estimateMinEVal),
                                           fabs(primme->stats.estimateMaxEVal));

   return 0;
}


/*******************************************************************************
 * Subroutine solve_H_RR - This procedure solves the eigenproblem for the
 *            matrix H.
 *        
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * H              The matrix V'*A*V
 * basisSize      The dimension of H, R, hU
 * ldH            The leading dimension of H
 * R              The factor R for the QR decomposition of (A - target*I)*V
 * ldR            The leading dimension of R
 * numConverged   Number of eigenvalues converged to determine ordering shift
 * lrwork         Length of the work array rwork
 * primme         Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hU                The left singular vectors of R
 * hVecs             The eigenvectors of H or the right singular vectors
 * hVals             The Ritz values
 * hSVals            The singular values of R
 * rwork             Workspace
 * iwork             Workspace in integers
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev was unsuccessful
 ******************************************************************************/

int solve_H_RR_dprimme(double *H, int ldH, double *hVecs,
   int ldhVecs, double *hVals, int basisSize, int numConverged, int lrwork,
   double *rwork, int *iwork, primme_params *primme) {

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

   /* Some LAPACK implementations don't like zero-size matrices */
   if (basisSize == 0) return 0;

   /* Return memory requirements */
   if (H == NULL) {
#ifdef NUM_ESSL
      return 2*basisSize + basisSize*(basisSize + 1)/2;
#else
      double rwork0;
      lrwork = 0;
      Num_dsyev_dprimme("V", "U", basisSize, hVecs, basisSize, hVals, &rwork0, 
            -1, &info);

      if (info != 0) {
         primme_PushErrorMessage(Primme_solve_h, Primme_num_dsyev, info, __FILE__, 
               __LINE__, primme);
         return NUM_DSYEV_FAILURE;
      }
      lrwork += (int)*(double*)&rwork0;
      return lrwork;
#endif
   }

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
            rwork[idx] = H[ldH*j+i];
            idx++;
         }
      }
   }
   else { /* (primme->target == primme_largest)  */
      for (j=0; j < basisSize; j++) {
         for (i=0; i <= j; i++) {
            rwork[idx] = -H[ldH*j+i];
            idx++;
         }
      }
   }
   
   apSize = basisSize*(basisSize + 1)/2;
   lrwork = lrwork - apSize;

   info = Num_dspev_dprimme(21, rwork, hVals, hVecs, ldhVecs, basisSize, 
      &rwork[apSize], lrwork);

   if (info != 0) {
      primme_PushErrorMessage(Primme_solve_h, Primme_num_dspev, info, __FILE__, 
         __LINE__, primme);
      return NUM_DSPEV_FAILURE;
   }

#else /* NUM_ESSL */
   if (primme->target != primme_largest) {
      for (j=0; j < basisSize; j++) {
         for (i=0; i <= j; i++) { 
            hVecs[ldhVecs*j+i] = H[ldH*j+i];
         }
      }      
   }
   else { /* (primme->target == primme_largest) */
      for (j=0; j < basisSize; j++) {
         for (i=0; i <= j; i++) { 
            hVecs[ldhVecs*j+i] = -H[ldH*j+i];
         }
      }
   }

   Num_dsyev_dprimme("V", "U", basisSize, hVecs, ldhVecs, hVals, rwork, 
                lrwork, &info);

   if (info != 0) {
      primme_PushErrorMessage(Primme_solve_h, Primme_num_dsyev, info, __FILE__, 
         __LINE__, primme);
      return NUM_DSYEV_FAILURE;
   }
#endif /* NUM_ESSL */

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

      /* TODO: order properly when numTargetShifts > 1 */

      targetShift = 
        primme->targetShifts[min(primme->numTargetShifts-1, numConverged)];

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
      permute_vecs_d(hVals, 1, basisSize, 1, permu, (double*)rwork, permw);
      permute_vecs_d(hVecs, basisSize, basisSize, ldhVecs, permu, rwork, permw);
   }

   return 0;   
}

/*******************************************************************************
 * Subroutine solve_H_Ref - This procedure solves the singular value
 *            decomposition of matrix R
 *        
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * H             The matrix V'*A*V
 * ldH           The leading dimension of H
 * R             The R factor for the QR decomposition of (A - target*I)*V
 * ldR           The leading dimension of R
 * basisSize     Current size of the orthonormal basis V
 * lrwork        Length of the work array rwork
 * primme        Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hVecs             The right singular vectors of R
 * hU                The left singular vectors of R
 * hSVals            The singular values of R
 * hVals             The Ritz values of the vectors in hVecs
 * rwork             Must be of size at least 3*maxBasisSize
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev was unsuccessful
 ******************************************************************************/

int solve_H_Ref_dprimme(double *H, int ldH, double *hVecs,
   int ldhVecs, double *hU, int ldhU, double *hSVals, double *R, int ldR,
   double *hVals, int basisSize, int lrwork, double *rwork, primme_params *primme) {

   int i, j; /* Loop variables    */
   int info; /* dsyev error value */
   double tpone = +1.0e+00, tzero = +0.0e+00, ztmp;

   /*lingfei: primme_svds. if the refined projection is used after
   the RayRitz projection.  It is designed for interior eigenvalue
   problem, especially for singular value problem with augmented
   method. In this case, we need compute:
   1) hVecs = Vm**T*A**T*A*Vm - 2*theta*Vm**T*A*Vm + theta*theta*Im
            = WTWm - 2*theta*Hm + theta*theta*Im;
   2) Another approach is for two-stages SVD problem;
      (AVm - theta*Vm) = Wm - thetaVm = QR; */

   /* Some LAPACK implementations don't like zero-size matrices */
   if (basisSize == 0) return 0;

   /* Return memory requirements */
   if (H == NULL) {
      double rwork0;
      lrwork = 0;
      Num_dgesvd_dprimme("S", "O", basisSize, basisSize, R, basisSize, 
            NULL, NULL, basisSize, hVecs, basisSize, &rwork0, -1, &info);

      if (info != 0) {
         primme_PushErrorMessage(Primme_solve_h, Primme_num_dgesvd, info, __FILE__, 
               __LINE__, primme);
         return NUM_DGESVD_FAILURE;
      }
      lrwork += (int)*(double*)&rwork0;
      return lrwork;
   }

   /* Copy upper triangular part of R into hVecs and zero lower triangular
      part of hVecs */
   Num_copy_trimatrix_dprimme(R, basisSize, basisSize, ldR, 0, 0, hVecs, ldhVecs, 1);

   /*Since Ritz vectors in hVecs is not needed, we use hVecs to hold refined
     Ritz vectors. Note gesvd returns transpose(V) rather than V and sorted in
     descending order of the singular values */

   Num_dgesvd_dprimme(hU?"S":"N", "O", basisSize, basisSize, hVecs, ldhVecs,
         hSVals?hSVals:hVals, hU, ldhU, hVecs, ldhVecs, rwork, lrwork, &info);

   if (info != 0) {
      primme_PushErrorMessage(Primme_solve_h, Primme_num_dgesvd, info, __FILE__, 
            __LINE__, primme);
      return NUM_DGESVD_FAILURE;
   }

   /* Transpose back V and rearrange V, hSVals and hU in ascending order
      of singular values */
   for (j=0; j < basisSize; j++) {
      for (i=0; i < basisSize; i++) { 
         rwork[basisSize*j+i] = hVecs[ldhVecs*i + (basisSize-1-j)];
      }
   }      
   Num_copy_matrix_dprimme(rwork, basisSize, basisSize, basisSize, hVecs, ldhVecs);

   if (hU) {
      Num_copy_matrix_dprimme(hU, basisSize, basisSize, ldhU, rwork, basisSize);
      for (j=0; j < basisSize; j++) {
         for (i=0; i < basisSize; i++) {
            hU[ldhU*j+i] = rwork[basisSize*(basisSize-j-1) + i];
         }
      }
   }

   if (hSVals) {
      Num_dcopy_dprimme(basisSize, hSVals, 1, (double*)rwork, 1);
      for (i=0; i < basisSize; i++)
         hSVals[i] = ((double*)rwork)[basisSize-i-1];
   }

   /* compute Rayleigh quotient lambda_i = x_i'*H*x_i */
   Num_symm_dprimme("L", "U", basisSize, basisSize, tpone, H,
      ldH, hVecs, ldhVecs, tzero, rwork, basisSize);

   for (i=0; i<basisSize; i++) {
      ztmp = Num_dot_dprimme(basisSize, &hVecs[ldhVecs*i], 1, &rwork[basisSize*i], 1);
      hVals[i] = ztmp;
   }

   return 0;
}
