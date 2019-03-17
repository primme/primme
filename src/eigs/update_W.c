/*******************************************************************************
 * Copyright (c) 2018, College of William & Mary
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
 * File: update_W.c
 *
 * Purpose - Computes A*V(:,0) through A*V(:,blockSize-1) after
 *           V has been expanded by blksze correction vectors.
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/update_W.c"
#endif


#include "numerical.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "update_W.h"
#include "auxiliary_eigs.h"
#include "ortho.h"
#endif

#ifdef SUPPORTED_TYPE

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

TEMPLATE_PLEASE
int matrixMatvec_Sprimme(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *W, PRIMME_INT ldW, int basisSize, int blockSize,
      primme_context ctx) {

   primme_params *primme = ctx.primme;

   if (blockSize <= 0)
      return 0;

   assert(ldV >= nLocal && ldW >= nLocal);
   assert(primme->ldOPs == 0 || primme->ldOPs >= nLocal);

   double t0 = primme_wTimer();

   /* Cast V and W */

   SCALAR *Vb = &V[ldV * basisSize], *Wb = &W[ldW * basisSize];
   void *V0, *W0;
   PRIMME_INT ldV0, ldW0;
   CHKERR(Num_matrix_astype_Sprimme(Vb, nLocal, blockSize, ldV,
         PRIMME_OP_SCALAR, &V0, &ldV0, primme->matrixMatvec_type, 1 /* alloc */,
         1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_Sprimme(Wb, nLocal, blockSize, ldW,
         PRIMME_OP_SCALAR, &W0, &ldW0, primme->matrixMatvec_type, 1 /* alloc */,
         0 /* no copy */, ctx));

   /* W(:,c) = A*V(:,c) for c = basisSize:basisSize+blockSize-1 */

   int ierr = 0;
   CHKERRM(
         (primme->matrixMatvec(V0, &ldV0, W0, &ldW0, &blockSize, primme, &ierr),
               ierr),
         PRIMME_USER_FAILURE, "Error returned by 'matrixMatvec' %d", ierr);

   /* Copy back W */

   CHKERR(Num_matrix_astype_Sprimme(W0, nLocal, blockSize, ldW0,
         primme->matrixMatvec_type, (void **)&Wb, &ldW,
         PRIMME_OP_SCALAR, 0 /* not alloc */, 1 /* copy */, ctx));

   if (Vb != V0) CHKERR(Num_free_Sprimme((SCALAR*)V0, ctx));
   if (Wb != W0) CHKERR(Num_free_Sprimme((SCALAR*)W0, ctx));

   primme->stats.timeMatvec += primme_wTimer() - t0;
   primme->stats.numMatvecs += blockSize;

   return 0;
}

/*******************************************************************************
 * Subroutine massMatrixMatvec - Computes B*V(:,nv+1:nv+blksze)
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
 * BV          B*V
 ******************************************************************************/

TEMPLATE_PLEASE
int massMatrixMatvec_Sprimme(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *BV, PRIMME_INT ldBV, int basisSize, int blockSize,
      primme_context ctx) {

   primme_params *primme = ctx.primme;

   if (blockSize <= 0)
      return 0;

   assert(ldV >= nLocal && ldBV >= nLocal);
   assert(primme->ldOPs == 0 || primme->ldOPs >= nLocal);

   double t0 = primme_wTimer();

   /* Cast V and BV */

   SCALAR *Vb = &V[ldV * basisSize], *BVb = &BV[ldBV * basisSize];
   void *V0, *BV0;
   PRIMME_INT ldV0, ldBV0;
   CHKERR(Num_matrix_astype_Sprimme(Vb, nLocal, blockSize, ldV,
         PRIMME_OP_SCALAR, &V0, &ldV0, primme->massMatrixMatvec_type,
         1 /* alloc */, 1 /* copy */, ctx));
   CHKERR(Num_matrix_astype_Sprimme(BVb, nLocal, blockSize, ldBV,
         PRIMME_OP_SCALAR, &BV0, &ldBV0, primme->massMatrixMatvec_type,
         1 /* alloc */, 0 /* no copy */, ctx));

   /* BV(:,c) = B*V(:,c) for c = basisSize:basisSize+blockSize-1 */

   int ierr = 0;
   CHKERRM((primme->massMatrixMatvec(
                  V0, &ldV0, BV0, &ldBV, &blockSize, primme, &ierr),
                 ierr),
         PRIMME_USER_FAILURE, "Error returned by 'massMatrixMatvec' %d", ierr);

   /* Copy back BV */

   CHKERR(Num_matrix_astype_Sprimme(BV0, nLocal, blockSize, ldBV0,
         primme->matrixMatvec_type, (void **)&BVb, &ldBV, PRIMME_OP_SCALAR,
         0 /* not alloc */, 1 /* copy */, ctx));

   if (Vb != V0) CHKERR(Num_free_Sprimme((SCALAR*)V0, ctx));
   if (BVb != BV0) CHKERR(Num_free_Sprimme((SCALAR*)BV0, ctx));

   primme->stats.timeMatvec += primme_wTimer() - t0;
   primme->stats.numMatvecs += blockSize;

   return 0;
}

/*******************************************************************************
 * Subroutine update_QR - Computes the QR factorization (A-targetShift*B)*V
 *    updating only the columns nv:nv+blockSize-1 of Q and R.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * BV         B*V
 * nLocal     Number of rows of each vector stored on this node
 * ldBV       The leading dimension of BV
 * W          A*V
 * ldW        The leading dimension of W
 * basisSize  Number of vectors in V
 * blockSize  The current block size
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * Q          The Q factor
 * R          The R factor
 * QtQ        Q'Q
 ******************************************************************************/

TEMPLATE_PLEASE
int update_Q_Sprimme(SCALAR *BV, PRIMME_INT nLocal, PRIMME_INT ldBV, SCALAR *W,
      PRIMME_INT ldW, SCALAR *Q, PRIMME_INT ldQ, HSCALAR *R, int ldR,
      HSCALAR *QtQ, int ldQtQ, double targetShift, int basisSize,
      int blockSize, int *nQ, primme_context ctx) {

   int i;

   /* Quick exit */

   if (blockSize <= 0 || R == NULL) return 0;

   assert(ldBV >= nLocal && ldW >= nLocal && ldQ >= nLocal &&
          ldR >= basisSize + blockSize);

   /* Q(:,c) = W(:,c) - BV(:,c)*target for c = basisSize:basisSize+blockSize-1 */
   for (i=basisSize; i<basisSize+blockSize; i++) {
      Num_compute_residual_Sprimme(
            nLocal, targetShift, &BV[ldBV * i], &W[ldW * i], &Q[ldQ * i], ctx);
   }

   /* Ortho Q(:,c) for c = basisSize:basisSize+blockSize-1 */

   CHKERR(ortho_block_Sprimme(Q, ldQ, QtQ, ldQtQ, R,
         ldR, *nQ, *nQ + blockSize - 1, NULL, 0, 0, NULL, 0, nLocal,
         ctx.primme->maxBasisSize, nQ, ctx));

   /* Zero the lower-left part of R */

   Num_zero_matrix_SHprimme(&R[basisSize], blockSize, basisSize, ldR, ctx);

   return 0;
}

#endif /* SUPPORTED_TYPE */
