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
#include "update_W_d.h"
#include "ortho_d.h"
#include "numerical_d.h"


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

void matrixMatvec_dprimme(double *V, int nLocal, int ldV, double *W,
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

int update_Q_dprimme(double *V, int nLocal, int ldV, double *W, int ldW,
      double *Q, int ldQ, double *R, int ldR, double targetShift, int basisSize,
      int blockSize, double *rwork, int rworkSize, double machEps, primme_params *primme) {

   int i, ret;

   /* Return memory requirement */
   if (V == NULL) {
      return ortho_dprimme(NULL, 0, NULL, 0, basisSize,
         basisSize+blockSize-1, NULL, 0, 0, primme->nLocal, 
         NULL, machEps, NULL, 0, primme);
   }

   /* Quick exit */

   if (blockSize <= 0 || Q == NULL || R == NULL) return 0;

   assert(ldV >= nLocal && ldW >= nLocal && ldQ >= nLocal && ldR >= basisSize+blockSize);   

   /* Q(:,c) = W(:,c) - V(:,c)*target for c = basisSize:basisSize+blockSize-1 */
   for (i=basisSize; i<basisSize+blockSize; i++) {
      Num_compact_res_dprimme(nLocal, targetShift, &V[ldV*i], &W[ldW*i],
            NULL, NULL, NULL, &Q[ldQ*i]);
   }

   /* Ortho Q(:,c) for c = basisSize:basisSize+blockSize-1 */
   ret = ortho_dprimme(Q, ldQ, R, ldR, basisSize, basisSize+blockSize-1, NULL,
         0, 0, nLocal, primme->iseed, machEps, rwork, rworkSize, primme);

   return ret;
}
