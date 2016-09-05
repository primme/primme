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
 * File: update_W.c
 *
 * Purpose - Computes A*V(:,0) through A*V(:,blockSize-1) after
 *           V has been expanded by blksze correction vectors.
 *
 ******************************************************************************/

#include <assert.h>
#include "primme.h"
#include "numerical_d.h"
#include "update_W_d.h"
#include "ortho_d.h"


/*******************************************************************************
 * Subroutine matrixMatvec_ - Computes A*V(:,nv+1) through A*V(:,nv+blksze)
 *           where V(:,nv+1:nv+blksze) are the new correction vectors.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V          The orthonormal basis
 * nLocal     Number of rows of each vector stored on this node
 * ldV        The leading dimension of V
 * ldW        The leading dimension of W
 * basisSize  Number of vectors in V
 * blockSize  The current block size
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * W          A*V
 ******************************************************************************/

int matrixMatvec_dprimme(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *W, PRIMME_INT ldW, int basisSize, int blockSize,
      primme_params *primme) {

   int i, ONE=1, ierr=0;

   if (blockSize <= 0) return 0;

   /* W(:,c) = A*V(:,c) for c = basisSize:basisSize+blockSize-1 */
   if (ldV == nLocal && ldW == nLocal) {
      primme->matrixMatvec(&V[ldV*basisSize], &W[ldW*basisSize], &blockSize,
            primme);
   }
   else {
      for (i=0; i<basisSize; i++) {
         primme->matrixMatvec(&V[ldV*(basisSize+i)], &W[ldW*(basisSize+i)], &ONE,
               primme);
      }
   }

   primme->stats.numMatvecs += blockSize;

   return ierr;

}

/*******************************************************************************
 * Subroutine update_QR - Computes the QR factorization (A-targetShift*I)*V
 *    updating only the columns nv:nv+blockSize-1 of Q and R.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V          The orthonormal basis
 * nLocal     Number of rows of each vector stored on this node
 * ldV        The leading dimension of V
 * W          A*V
 * ldW        The leading dimension of W
 * basisSize  Number of vectors in V
 * blockSize  The current block size
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * Q          The Q factor
 * R          The R factor
 ******************************************************************************/

int update_Q_dprimme(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *W, PRIMME_INT ldW, SCALAR *Q, PRIMME_INT ldQ, SCALAR *R, int ldR,
      double targetShift, int basisSize, int blockSize, SCALAR *rwork,
      size_t *rworkSize, double machEps, primme_params *primme) {

   int i, j, ret;

   /* Return memory requirement */
   if (V == NULL) {
      ortho_dprimme(NULL, 0, NULL, 0, basisSize,
         basisSize+blockSize-1, NULL, 0, 0, primme->nLocal, 
         NULL, machEps, NULL, rworkSize, primme);
      return 0;
   }

   /* Quick exit */

   if (blockSize <= 0 || Q == NULL || R == NULL) return 0;

   assert(ldV >= nLocal && ldW >= nLocal && ldQ >= nLocal && ldR >= basisSize+blockSize);   

   /* Q(:,c) = W(:,c) - V(:,c)*target for c = basisSize:basisSize+blockSize-1 */
   for (i=basisSize; i<basisSize+blockSize; i++) {
      Num_compute_residual_dprimme(nLocal, targetShift, &V[ldV*i], &W[ldW*i],
            &Q[ldQ*i]);
   }

   /* Ortho Q(:,c) for c = basisSize:basisSize+blockSize-1 */
   ret = ortho_dprimme(Q, ldQ, R, ldR, basisSize, basisSize+blockSize-1, NULL,
         0, 0, nLocal, primme->iseed, machEps, rwork, rworkSize, primme);

   /* Zero the lower triangular part of R */
   for (i=basisSize; i<basisSize+blockSize; i++) {
      for (j=i+1; j<ldR; j++) {
         R[ldR*i+j] = 0.0;
      }
   }

   return ret;
}
