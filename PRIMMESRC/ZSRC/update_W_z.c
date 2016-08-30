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
#include "update_W_z.h"
#include "ortho_z.h"
#include "numerical_z.h"


/*******************************************************************************
 * Subroutine matrixMatvec_ - Computes A*V(:,nv+1) through A*V(:,nv+blksze)
 *           where V(:,nv+1:nv+blksze) are the new correction vectors.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V          The orthonormal basis
 * basisSize  Number of vectors in V
 * blockSize  The current block size
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * W  A*V
 ******************************************************************************/

void matrixMatvec_zprimme(__PRIMME_COMPLEX_DOUBLE__ *V, int nLocal, int ldV, __PRIMME_COMPLEX_DOUBLE__ *W,
   int ldW, int basisSize, int blockSize, primme_params *primme) {

   int i, ONE=1;

   if (blockSize <= 0) return;

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

}

int update_Q_zprimme(__PRIMME_COMPLEX_DOUBLE__ *V, int nLocal, int ldV, __PRIMME_COMPLEX_DOUBLE__ *W, int ldW,
      __PRIMME_COMPLEX_DOUBLE__ *Q, int ldQ, __PRIMME_COMPLEX_DOUBLE__ *R, int ldR, double targetShift, int basisSize,
      int blockSize, __PRIMME_COMPLEX_DOUBLE__ *rwork, int rworkSize, double machEps, primme_params *primme) {

   int i, j, ret;

   /* Return memory requirement */
   if (V == NULL) {
      return ortho_zprimme(NULL, 0, NULL, 0, basisSize,
         basisSize+blockSize-1, NULL, 0, 0, primme->nLocal, 
         NULL, machEps, NULL, 0, primme);
   }

   /* Quick exit */

   if (blockSize <= 0 || Q == NULL || R == NULL) return 0;

   assert(ldV >= nLocal && ldW >= nLocal && ldQ >= nLocal && ldR >= basisSize+blockSize);   

   /* Q(:,c) = W(:,c) - V(:,c)*target for c = basisSize:basisSize+blockSize-1 */
   for (i=basisSize; i<basisSize+blockSize; i++) {
      Num_compute_residual_zprimme(nLocal, targetShift, &V[ldV*i], &W[ldW*i],
            &Q[ldQ*i]);
   }

   /* Ortho Q(:,c) for c = basisSize:basisSize+blockSize-1 */
   ret = ortho_zprimme(Q, ldQ, R, ldR, basisSize, basisSize+blockSize-1, NULL,
         0, 0, nLocal, primme->iseed, machEps, rwork, rworkSize, primme);

   /* Zero the lower triangular part of R */
   for (i=basisSize; i<basisSize+blockSize; i++) {
      for (j=i+1; j<ldR; j++) {
         R[ldR*i+j] = 0.0;
      }
   }

   return ret;
}
