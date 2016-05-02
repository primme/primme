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
#include "const.h"
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
 * basisSize      The dimension of H, R, QtV and hU
 * ldH            The leading dimension of H
 * R              The factor R for the QR decomposition of (A - target*I)*V
 * ldR            The leading dimension of R
 * QtV            Q'*V
 * ldQtV          The leading dimension of QtV
 * numConverged   Number of eigenvalues converged to determine ordering shift
 * lrwork         Length of the work array rwork
 * primme         Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hU             The left singular vectors of R or the eigenvectors of QtV/R
 * ldhU           The leading dimension of hU
 * hVecs          The coefficient vectors such as V*hVecs will be the Ritz vectors
 * ldhVecs        The leading dimension of hVecs
 * hVals          The Ritz values
 * hSVals         The singular values of R
 * rwork          Workspace
 * iwork          Workspace in integers
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev was unsuccessful
 ******************************************************************************/

int solve_H_dprimme(double *H, int basisSize, int ldH, double *R, int ldR,
   double *QtV, int ldQtV, double *hU, int ldhU, double *hVecs, int ldhVecs,
   double *hVals, double *hSVals, int numConverged, double machEps, int lrwork,
   double *rwork, int *iwork, primme_params *primme) {

   int i, ret;

   switch (primme->projectionParams.projection) {
   case primme_proj_RR:
      ret = solve_H_RR_dprimme(H, ldH, hVecs, ldhVecs, hVals, basisSize,
            numConverged, lrwork, rwork, iwork, primme);
      break;

   case primme_proj_harmonic:
      ret = solve_H_Harm_dprimme(H, ldH, QtV, ldQtV, R, ldR, hVecs, ldhVecs, hU,
            ldhU, hVals, basisSize, numConverged, machEps, lrwork, rwork, iwork, primme);
      break;

   case primme_proj_refined:
      ret = solve_H_Ref_dprimme(H, ldH, hVecs, ldhVecs, hU, ldhU, hSVals, 
            R, ldR, hVals, basisSize, numConverged, lrwork, rwork, iwork, primme);
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
 * numConverged   Number of eigenvalues converged to determine ordering shift
 * lrwork         Length of the work array rwork
 * primme         Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * hVecs          The eigenvectors of H or the right singular vectors
 * ldhVecs        The leading dimension of hVecs
 * hVals          The Ritz values
 * hSVals         The singular values of R
 * rwork          Workspace
 * iwork          Workspace in integers
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev was unsuccessful
 ******************************************************************************/

static int solve_H_RR_dprimme(double *H, int ldH, double *hVecs,
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
      else if (primme->target == primme_largest_abs) {

         j = 0;
         i = basisSize-1;
         index = 0;
         while (i>=j) {
            if (fabs(hVals[i]-targetShift) > fabs(hVals[j]-targetShift)) 
               permu[index++] = i--;
            else 
               permu[index++] = j++;
         }

      }

      /* ---------------------------------------------------------------- */
      /* Reorder hVals and hVecs according to the permutation             */
      /* ---------------------------------------------------------------- */
      permute_vecs_dprimme(hVals, 1, basisSize, 1, permu, (double*)rwork, permw);
      permute_vecs_dprimme(hVecs, basisSize, basisSize, ldhVecs, permu, rwork, permw);
   }

   return 0;   
}

/*******************************************************************************
 * Subroutine solve_H_Harm - This procedure implements the harmonic extraction
 *    in a novelty way. In standard harmonic the next eigenproblem is solved:
 *       V'*(A-s*I)'*(A-s*I)*V*X = V'*(A-s*I)'*V*X*L,
 *    where (L_{i,i},X_i) are the harmonic-Ritz pairs. In practice, it is
 *    computed (A-s*I)*V = Q*R and it is solved instead:
 *       R*X = Q'*V*X*L,
 *    which is a generalized non-Hermitian problem. Instead of dealing with
 *    complex solutions, which are unnatural in context of Hermitian problems,
 *    we propose the following. Note that,
 *       (A-s*I)*V = Q*R -> Q'*V*inv(R) = Q'*inv(A-s*I)*Q.
 *    And note that Q'*V*inv(R) is Hermitian if A is, and also that
 *       Q'*V*inv(R)*Y = Y*inv(L) ->  Q'*V*X*L = R*X,
 *    with Y = R*X. So this routine computes X by solving the Hermitian problem
 *    Q'*V*inv(R).
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
 * hVecs         The orthogonal basis of inv(R) * eigenvectors of QtV/R
 * ldhVecs       The leading dimension of hVecs
 * hU            The eigenvectors of QtV/R
 * ldhU          The leading dimension of hU
 * hVals         The Ritz values of the vectors in hVecs
 * rwork         Workspace
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev was unsuccessful
 ******************************************************************************/

static int solve_H_Harm_dprimme(double *H, int ldH, double *QtV, int ldQtV,
   double *R, int ldR, double *hVecs, int ldhVecs, double *hU, int ldhU,
   double *hVals, int basisSize, int numConverged, double machEps, int lrwork,
   double *rwork, int *iwork, primme_params *primme) {

   int i, ret;
   double tzero = +0.0e+00, tpone = +1.0e+00;
   double *oldTargetShifts, zero=0.0;
   primme_target oldTarget;

   /* Some LAPACK implementations don't like zero-size matrices */
   if (basisSize == 0) return 0;

   /* Return memory requirements */
   if (QtV == NULL) {
      return solve_H_RR_dprimme(QtV, ldQtV, hVecs, ldhVecs, hVals, basisSize,
         0, lrwork, rwork, iwork, primme);
   }

   /* QAQ = QtV*inv(R) */

   Num_copy_matrix_dprimme(QtV, basisSize, basisSize, ldQtV, hVecs, ldhVecs);
   Num_trsm_dprimme("R", "U", "N", "N", basisSize, basisSize, tpone, R, ldR,
         hVecs, ldhVecs);

   /* Compute eigenpairs of QAQ */

   oldTargetShifts = primme->targetShifts;
   oldTarget = primme->target;
   primme->targetShifts = &zero;
   switch(primme->target) {
      case primme_closest_geq:
         primme->target = primme_largest;
         break;
      case primme_closest_leq:
         primme->target = primme_smallest;
         break;
      case primme_closest_abs:
         primme->target = primme_largest_abs;
         break;
      default:
         assert(0);
   }
   ret = solve_H_RR_dprimme(hVecs, ldhVecs, hVecs, ldhVecs, hVals, basisSize,
         0, lrwork, rwork, iwork, primme);
   primme->targetShifts = oldTargetShifts;
   primme->target = oldTarget;
   if (ret != 0) return ret;

   Num_copy_matrix_dprimme(hVecs, basisSize, basisSize, ldhVecs, hU, ldhU);

   /* Transfer back the eigenvectors to V, hVecs = R\hVecs */

   Num_trsm_dprimme("L", "U", "N", "N", basisSize, basisSize, tpone, R, ldR,
         hVecs, ldhVecs);
   ret = ortho_dprimme(hVecs, ldhVecs, NULL, 0, 0, basisSize-1, NULL, 0, 0,
         basisSize, primme->iseed, machEps, rwork, lrwork, primme);
   if (ret != 0) return ret;
 
   /* Compute Rayleigh quotient lambda_i = x_i'*H*x_i */

   Num_symm_dprimme("L", "U", basisSize, basisSize, tpone, H,
      ldH, hVecs, ldhVecs, tzero, rwork, basisSize);

   for (i=0; i<basisSize; i++) {
      double ztmp = Num_dot_dprimme(basisSize, &hVecs[ldhVecs*i], 1, &rwork[basisSize*i], 1);
      hVals[i] = ztmp;
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
 * hVecs         The right singular vectors of R
 * ldhVecs       The leading dimension of hVecs
 * hU            The left singular vectors of R
 * ldhU          The leading dimension of hU
 * hSVals        The singular values of R
 * hVals         The Ritz values of the vectors in hVecs
 * rwork         Workspace
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 Num_dsyev was unsuccessful
 ******************************************************************************/

static int solve_H_Ref_dprimme(double *H, int ldH, double *hVecs,
   int ldhVecs, double *hU, int ldhU, double *hSVals, double *R, int ldR,
   double *hVals, int basisSize, int targetShiftIndex, int lrwork, double *rwork,
   int *iwork, primme_params *primme) {

   int i, j; /* Loop variables    */
   int info; /* dsyev error value */
   double tpone = +1.0e+00, tzero = +0.0e+00, ztmp;

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

   /* Since Ritz vectors in hVecs is not needed, we use hVecs to hold refined
      Ritz vectors. Note gesvd returns transpose(V) rather than V and sorted in
      descending order of the singular values */

   Num_dgesvd_dprimme("S", "O", basisSize, basisSize, hVecs, ldhVecs,
         hSVals, hU, ldhU, hVecs, ldhVecs, rwork, lrwork, &info);

   if (info != 0) {
      primme_PushErrorMessage(Primme_solve_h, Primme_num_dgesvd, info, __FILE__, 
            __LINE__, primme);
      return NUM_DGESVD_FAILURE;
   }

   /* Transpose back V */

   for (j=0; j < basisSize; j++) {
      for (i=0; i < basisSize; i++) { 
         rwork[basisSize*j+i] = hVecs[ldhVecs*i+j];
      }
   }      
   Num_copy_matrix_dprimme(rwork, basisSize, basisSize, basisSize, hVecs, ldhVecs);

   /* Rearrange V, hSVals and hU in ascending order of singular value   */
   /* if target is not largest abs.                                     */

   if (primme->target == primme_closest_abs 
         || primme->target == primme_closest_leq
         || primme->target == primme_closest_geq) {
      int *perm = iwork;
      int *iwork0 = iwork + basisSize;

      for (i=0; i<basisSize; i++) perm[i] = basisSize-1-i;
      permute_vecs_dprimme(hSVals, 1, basisSize, 1, perm, (double*)rwork, iwork0);
      permute_vecs_dprimme(hVecs, basisSize, basisSize, ldhVecs, perm, rwork, iwork0);
      permute_vecs_dprimme(hU, basisSize, basisSize, ldhU, perm, rwork, iwork0);
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

/*******************************************************************************
 * Function prepare_vecs - This subroutine checks that the
 *    conditioning of the coefficient vectors are good enough to converge
 *    with the requested accuracy. For now refined extraction is the only one that
 *    may present problems: two similar singular values in the projected problem
 *    may correspond to distinct eigenvalues in the original problem. If that is
 *    the case, the singular vector may have components of both eigenvectors,
 *    which prevents the residual norm be lower than some degree. Don't dealing
 *    with this may lead into stagnation.
 *
 *    It is checked that the next upper bound about the angle of the right
 *    singular vector v of A and the right singular vector vtilde of A+E,
 *
 *      sin(v,vtilde) <= sqrt(2)*||E||/sval_gap <= sqrt(2)*||A||*machEps/sval_gap,
 *
 *    is less than the upper bound about the angle of exact eigenvector u and
 *    the approximate eigenvector utilde,
 *
 *      sin(u,utilde) <= ||r||/eval_gap <= ||A||*eps/eval_gap.
 *
 *    (see pp. 211 in Matrix Algorithms vol. 2 Eigensystems, G. W. Steward).
 *
 *    If the inequality doesn't hold, do Rayleigh-Ritz onto the subspace
 *    spanned by both vectors.
 *
 *    we have found cases where this is not enough or the performance improves
 *    if Rayleigh-Ritz is also done when the candidate vector has a small
 *    angle with the last vector in V and when the residual norm is larger than
 *    the singular value.
 *
 *    When only one side of the shift is targeted (primme_closest_leq/geq), we
 *    allow to take eigenvalue of the other side but close to the shift. In the
 *    current heuristic they shouldn't be farther than the smallest residual
 *    norm in the block. This heuristic obtained good results in solving the
 *    augmented problem with shifts from solving the normal equations.
 *
 * NOTE: this function assumes hSVals are arranged in increasing order.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * basisSize    Projected problem size
 * i0           Index of the first pair to check
 * blockSize    Number of candidates wanted
 * H            The matrix V'*A*V
 * ldH          The leading dimension of H
 * hVals        The Ritz values
 * hSVals       The singular values of R
 * hVecs        The coefficient vectors
 * ldhVecs      The leading dimension of hVecs
 * targetShiftIndex The target shift used in (A - targetShift*B) = Q*R
 * arbitraryVecs The number of vectors modified (input/output)
 * smallestResNorm The smallest residual norm in the block
 * flags        Array indicating the convergence of the Ritz vectors
 * RRForAll     If false compute Rayleigh-Ritz only in clusters with
 *              candidates. If true, compute it in every cluster.
 * machEps      Machine precision
 * rworkSize    The length of rwork
 * rwork        Workspace
 * iwork        Integer workspace
 * primme       Structure containing various solver parameters
 *
 ******************************************************************************/

int prepare_vecs_dprimme(int basisSize, int i0, int blockSize,
      double *H, int ldH, double *hVals, double *hSVals, double *hVecs,
      int ldhVecs, int targetShiftIndex, int *arbitraryVecs,
      double smallestResNorm, int *flags, int RRForAll, double machEps,
      int rworkSize, double *rwork, int *iwork, primme_params *primme) {

   int i, j, k;         /* Loop indices */
   int candidates;      /* Number of eligible pairs */
   int someCandidate;   /* If there is an eligible pair in the cluster */
   double tpone = +1.0e+00, tzero = +0.0e+00;
   double targetShift, aNorm;
   int left, right, *perm, ret;

   /* Quick exit */

   if (primme->projectionParams.projection != primme_proj_refined
         || basisSize == 0) {
      if (arbitraryVecs) *arbitraryVecs = 0;
      return 0;
   }

   /* Return memory requirement */

   if (H == NULL) {
      return basisSize*basisSize*2 + /* aH, ahVecs */
         max(
               compute_submatrix_dprimme(NULL, basisSize, 0, NULL, basisSize, 0,
                  NULL, 0, NULL, 0),
               solve_H_RR_dprimme(NULL, 0, NULL, 0, NULL, basisSize, 0, 0, NULL,
                  NULL, primme));
   }

   /* Special case: If (basisSize+numLocked) is the entire space, */
   /* then everything should be converged. Just do RR with the    */
   /* entire space.                                               */
 
   if (basisSize + (primme->locking?primme->initSize:0) 
         + primme->numOrthoConst >= primme->n) {

      /* Compute and sort eigendecomposition aH*ahVecs = ahVecs*diag(hVals(j:i-1)) */
      ret = solve_H_RR_dprimme(H, ldH, hVecs, ldhVecs, hVals, basisSize,
            targetShiftIndex, rworkSize, rwork, iwork, primme);
      if (ret != 0) return ret;

      *arbitraryVecs = 0;

      return 0;
   }

   targetShift = primme->targetShifts[targetShiftIndex];
   aNorm = (primme->aNorm <= 0.0) ?
      primme->stats.estimateLargestSVal : primme->aNorm;

   for (candidates=0, i=min(*arbitraryVecs,basisSize), j=i0;
         j < basisSize && candidates < blockSize; ) {

      /* -------------------------------------------------------------------- */
      /* Count all eligible values (candidates) from j up to i.               */
      /* -------------------------------------------------------------------- */

      for ( ; j < i; j++)
         if ((!flags || flags[j] == UNCONVERGED) && (
               (primme->target == primme_closest_leq
                && hVals[j]-smallestResNorm <= targetShift) ||
               (primme->target == primme_closest_geq
                && hVals[j]+smallestResNorm >= targetShift) ||
               (primme->target != primme_closest_leq
                && primme->target != primme_closest_geq)))
            candidates++;
     
      if (candidates >= blockSize) break;
 
      /* -------------------------------------------------------------------- */
      /* Find the first i-th vector i>j with enough good conditioning, ie.,   */
      /* the singular value is separated enough from the rest (see the header */
      /* comment in this function). Also check if there is an unconverged     */
      /* value in the block.                                                  */
      /* -------------------------------------------------------------------- */

      for (i=j+1, someCandidate=0; i<basisSize; i++) {

         /* Check that this approximation:                                    */
         /* sin singular vector: max(hSVals)*machEps/(hSvals[i]-hSVals[i+1])  */
         /* is less than the next ones:                                       */
         /* sin eigenvector    : aNorm*eps/(hVals[i]-hVals[i+1])          */
         /* sin current and previous hVecs(:,i): hVecs(end,i)                 */
         /* NOTE: we don't want to check hVecs(end,i) just after restart, so  */
         /* we don't use the value when it is zero.                           */
         /* Also check that singular value is larger than the residual norm   */

         double ip0 = fabs(*(double*)&hVecs[(i-1)*ldhVecs+basisSize-1]);
         double ip = (flags && ip0 != 0.0) ? ip0 : HUGE_VAL;
         double minDiff = sqrt(2.0)*hSVals[basisSize-1]*machEps/
            min(ip, aNorm*primme->eps/fabs(hVals[i]-hVals[i-1]));

         if (!flags || flags[i-1] == UNCONVERGED) someCandidate = 1;

         if (fabs(hSVals[i]-hSVals[i-1]) >= minDiff 
               && fabs(hVals[i-1]-targetShift) < hSVals[i-1]+machEps*hSVals[basisSize-1]
               && hSVals[i-1] >= smallestResNorm)
            break;
      }
      i = min(i, basisSize);

      /* ----------------------------------------------------------------- */
      /* If the cluster j:i-1 is larger than one vector and there is some  */
      /* unconverged pair in there, compute the approximate eigenvectors   */
      /* with Rayleigh-Ritz. If RRForAll do also this when there is no     */
      /* candidate in the cluster.                                         */
      /* ----------------------------------------------------------------- */

      if (i-j > 1 && (someCandidate || RRForAll)) {
         double *rwork0 = rwork, *aH, *ahVecs;
         int rworkSize0 = rworkSize;
         int aBasisSize = i-j;
         aH = rwork0; rwork0 += aBasisSize*aBasisSize; rworkSize0 -= aBasisSize*aBasisSize;
         ahVecs = rwork0; rwork0 += aBasisSize*aBasisSize; rworkSize0 -= aBasisSize*aBasisSize;
         assert(rworkSize0 >= 0);

         /* aH = hVecs(:,j:i-1)'*H*hVecs(:,j:i-1) */
         compute_submatrix_dprimme(&hVecs[ldhVecs*j], aBasisSize,
               ldhVecs, H, basisSize, ldH, aH, aBasisSize, rwork0,
               rworkSize0);

         /* Compute and sort eigendecomposition aH*ahVecs = ahVecs*diag(hVals(j:i-1)) */
         ret = solve_H_RR_dprimme(aH, aBasisSize, ahVecs, aBasisSize,
               &hVals[j], aBasisSize, targetShiftIndex, rworkSize0, rwork0,
               iwork, primme);
         if (ret != 0) return ret;

         /* hVecs(:,j:i-1) = hVecs(:,j:i-1)*ahVecs */
         Num_gemm_dprimme("N", "N", basisSize, aBasisSize, aBasisSize,
               tpone, &hVecs[ldhVecs*j], ldhVecs, ahVecs, aBasisSize, tzero,
               rwork0, basisSize);
         Num_copy_matrix_dprimme(rwork0, basisSize, aBasisSize, basisSize,
               &hVecs[ldhVecs*j], ldhVecs);

         /* Indicate that before i may not be singular vectors */
         *arbitraryVecs = i;

         /* Remove converged flags from j upto i */
         if (flags && !RRForAll) for (k=j; k<i; k++) flags[k] = UNCONVERGED;
      }
   }

   if (primme->target != primme_closest_leq && primme->target != primme_closest_geq)
      return 0;

   /* -------------------------------------------------------------------- */
   /* Rearrange hVals and hVecs from 0 up to i-th putting the candidates   */
   /* first and then the unwanted pairs (eg. hVals[i]>targetshift if       */
   /* target is primme_closest_leq).                                       */
   /* -------------------------------------------------------------------- */

   /* Count all eligible values (candidates) from 0 up to i */

   for (j=0, candidates=0; j < i; j++)
      if (  (primme->target == primme_closest_leq && hVals[j]-smallestResNorm <= targetShift)
          ||(primme->target == primme_closest_geq && hVals[j]+smallestResNorm >= targetShift))
         candidates++;

   perm = iwork;
   iwork += basisSize;
   for (j=right=left=0; j < i; j++) {
      if (  (primme->target == primme_closest_leq && hVals[j]-smallestResNorm <= targetShift)
          ||(primme->target == primme_closest_geq && hVals[j]+smallestResNorm >= targetShift))
         perm[right++] = j;
      else
         perm[candidates+left++] = j;
   }

   permute_vecs_dprimme(hVals, 1, i, 1, perm, (double*)rwork, iwork);
   permute_vecs_dprimme(hSVals, 1, i, 1, perm, (double*)rwork, iwork);
   permute_vecs_dprimme(hVecs, basisSize, i, ldhVecs, perm, rwork, iwork);

   /* If something has changed between arbitraryVecs and i, notify */

   for (j=*arbitraryVecs; j<i; j++)
      if (perm[j] != j) *arbitraryVecs = j+1;

   return 0;
}
