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
 * File: auxiliary.c
 *
 * Purpose - Miscellanea functions used by PRIMME EIGS
 *
 ******************************************************************************/

#include <assert.h>
#include <math.h>
#include "const.h"
#include "numerical.h"
#include "globalsum.h"
#include "auxiliary_eigs.h"
#include "wtime.h"

/******************************************************************************
 * Function Num_compute_residual - This subroutine performs the next operation
 *    in a cache-friendly way:
 *
 *    r = Ax - eval*x
 *
 * PARAMETERS
 * ---------------------------
 * n           The number of rows of x, Ax and r
 * eval        The value to compute the residual vector r
 * x           The vector x
 * Ax          The vector Ax
 * r           On output r = Ax - eval*x
 *
 ******************************************************************************/

TEMPLATE_PLEASE
void Num_compute_residual_Sprimme(PRIMME_INT n, SCALAR eval, SCALAR *x, 
   SCALAR *Ax, SCALAR *r) {

   int k, M=min(n,PRIMME_BLOCK_SIZE);

   for (k=0; k<n; k+=M, M=min(M,n-k)) {
      Num_copy_Sprimme(M, &Ax[k], 1, &r[k], 1);
      Num_axpy_Sprimme(M, -eval, &x[k], 1, &r[k], 1);
   }

}

/******************************************************************************
 * Function Num_update_VWXR - This subroutine performs the next operations:
 *
 *    X0 = V*h(nX0b+1:nX0e), X1 = V*h(nX1b+1:nX1e), X2 = V*h(nX2b+1:nX2e)
 *    Wo = W*h(nWob+1:nWoe),
 *    R = W*h(nRb+1:nRe) - W*h(nRb+1:nRe)*diag(hVals(nRb+1:nRe)),
 *    Rnorms = norms(R),
 *    rnorms = norms(W*h(nrb+1:nre) - W*h(nrb+1:nre)*diag(hVals(nrb+1:nre)))
 *
 * NOTE: if Rnorms and rnorms are requested, nRb-nRe+nrb-nre < mV
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V, W        input basis
 * mV,nV,ldV   number of rows and columns and leading dimension of V and W
 * h           input rotation matrix
 * nh          Number of columns of h
 * ldh         The leading dimension of h
 * hVals       Array of values
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * X0          Output matrix V*h(nX0b:nX0e-1) (optional)
 * nX0b, nX0e  Range of columns of h
 * X1          Output matrix V*h(nX1b:nX1e-1) (optional)
 * nX1b, nX1e  Range of columns of h
 * X2          Output matrix V*h(nX2b:nX2e-1) (optional)
 * nX2b, nX2e  Range of columns of h
 * Wo          Output matrix W*h(nWob:nWoe-1) (optional)
 * nWob, nWoe  Range of columns of h
 * R           Output matrix (optional)
 * nRb, nRe    Range of columns of h and hVals
 * Rnorms      Output array with the norms of R (optional)
 * rnorms      Output array with the extra residual vector norms (optional)
 * nrb, nre    Columns of residual vector to compute the norm
 * 
 * NOTE: n*e, n*b are zero-base indices of ranges where the first value is
 *       included and the last isn't.
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_update_VWXR_Sprimme(SCALAR *V, SCALAR *W, PRIMME_INT mV, int nV,
      PRIMME_INT ldV, SCALAR *h, int nh, int ldh, REAL *hVals,
      SCALAR *X0, int nX0b, int nX0e, PRIMME_INT ldX0,
      SCALAR *X1, int nX1b, int nX1e, PRIMME_INT ldX1,
      SCALAR *X2, int nX2b, int nX2e, PRIMME_INT ldX2,
      SCALAR *Wo, int nWob, int nWoe, PRIMME_INT ldWo,
      SCALAR *R, int nRb, int nRe, PRIMME_INT ldR, REAL *Rnorms,
      REAL *rnorms, int nrb, int nre,
      SCALAR *rwork, int lrwork, primme_params *primme) {

   PRIMME_INT i;     /* Loop variables */
   int j;            /* Loop variables */
   int m=min(PRIMME_BLOCK_SIZE, mV);   /* Number of rows in the cache */
   int nXb, nXe, nYb, nYe, ldX, ldY;
   SCALAR *X, *Y;
   REAL *tmp, *tmp0;

   /* Return memory requirements */
   if (V == NULL) {
      return 2*m*nV;
   }

   /* R or Rnorms or rnorms imply W */
   assert(!(R || Rnorms || rnorms) || W);

   nXb = min(min(min(min(X0?nX0b:INT_MAX, X1?nX1b:INT_MAX), X2?nX2b:INT_MAX),
         R?nRb:INT_MAX), rnorms?nrb:INT_MAX);
   nXe = max(max(max(X0?nX0e:0, X1?nX1e:0), R?nRe:0), rnorms?nre:0);
   nYb = min(min(Wo?nWob:INT_MAX, R?nRb:INT_MAX), rnorms?nrb:INT_MAX);
   nYe = max(max(Wo?nWoe:0, R?nRe:0), rnorms?nre:0);

   assert(nXe <= nh || nXb >= nXe); /* Check dimension */
   assert(nYe <= nh || nYb >= nYe); /* Check dimension */
   assert((nXe-nXb+nYe-nYb)*m <= lrwork); /* Check workspace for X and Y */
   assert(2*(nRe-nRb+nre-nrb) <= lrwork); /* Check workspace for tmp and tmp0 */

   X = rwork;
   Y = rwork + m*(nXe-nXb);
   ldX = ldY = m;

   if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = 0.0;
   if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = 0.0;

   for (i=0; i < mV; i+=m, m=min(m,mV-i)) {
      /* X = V*h(nXb:nXe-1) */
      Num_gemm_Sprimme("N", "N", m, nXe-nXb, nV, 1.0,
         &V[i], ldV, &h[nXb*ldh], ldh, 0.0, X, ldX);

      /* X0 = X(nX0b-nXb:nX0e-nXb-1) */
      if (X0) Num_copy_matrix_Sprimme(&X[ldX*(nX0b-nXb)], m, nX0e-nX0b,
            ldX, &X0[i], ldX0);

      /* X1 = X(nX1b-nXb:nX1e-nXb-1) */
      if (X1) Num_copy_matrix_Sprimme(&X[ldX*(nX1b-nXb)], m, nX1e-nX1b,
            ldX, &X1[i], ldX1);

      /* X2 = X(nX2b-nXb:nX2e-nXb-1) */
      if (X2) Num_copy_matrix_Sprimme(&X[ldX*(nX2b-nXb)], m, nX2e-nX2b,
            ldX, &X2[i], ldX2);

      /* Y = W*h(nYb:nYe-1) */
      if (nYb < nYe) Num_gemm_Sprimme("N", "N", m, nYe-nYb, nV,
            1.0, &W[i], ldV, &h[nYb*ldh], ldh, 0.0, Y, ldY);

      /* Wo = Y(nWob-nYb:nWoe-nYb-1) */
      if (Wo) Num_copy_matrix_Sprimme(&Y[ldY*(nWob-nYb)], m, nWoe-nWob,
            ldY, &Wo[i], ldWo);

      /* R = Y(nRb-nYb:nRe-nYb-1) - X(nRb-nYb:nRe-nYb-1)*diag(nRb:nRe-1) */
      if (R) for (j=nRb; j<nRe; j++) {
         Num_compute_residual_Sprimme(m, hVals[j], &X[ldX*(j-nXb)], &Y[ldY*(j-nYb)],
               &R[i+ldR*(j-nRb)]);
         if (Rnorms) {
            Rnorms[j-nRb] +=
               REAL_PART(Num_dot_Sprimme(m, &R[i+ldR*(j-nRb)], 1,
                        &R[i+ldR*(j-nRb)], 1));
         }
      }

      /* rnorms = Y(nrb-nYb:nre-nYb-1) - X(nrb-nYb:nre-nYb-1)*diag(nrb:nre-1) */
      if (rnorms) for (j=nrb; j<nre; j++) {
         Num_compute_residual_Sprimme(m, hVals[j], &X[ldX*(j-nXb)], &Y[ldY*(j-nYb)],
               &Y[ldY*(j-nYb)]);
         rnorms[j-nrb] += 
            REAL_PART(Num_dot_Sprimme(m, &Y[ldY*(j-nYb)], 1,
                     &Y[ldY*(j-nYb)], 1));
      }
   }

   /* Reduce Rnorms and rnorms and sqrt the results */

   if (primme->numProcs > 1) {
      tmp = (REAL*)rwork;
      j = 0;
      if (R && Rnorms) for (i=nRb; i<nRe; i++) tmp[j++] = Rnorms[i-nRb];
      if (rnorms) for (i=nrb; i<nre; i++) tmp[j++] = rnorms[i-nrb];
      tmp0 = tmp+j;
      if (j) CHKERR(globalSum_Rprimme(tmp, tmp0, j, primme), -1);
      j = 0;
      if (R && Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = sqrt(tmp0[j++]);
      if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = sqrt(tmp0[j++]);
   }
   else {
      if (R && Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = sqrt(Rnorms[i-nRb]);
      if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = sqrt(rnorms[i-nrb]);
   }

   return 0; 
}

/*******************************************************************************
 * Subroutine applyPreconditioner - apply preconditioner to V
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V          The orthonormal basis
 * nLocal     Number of rows of each vector stored on this node
 * ldV        The leading dimension of V
 * ldW        The leading dimension of W
 * blockSize  The number of columns of V and W.
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * W          M*V
 ******************************************************************************/

TEMPLATE_PLEASE
int applyPreconditioner_Sprimme(SCALAR *V, PRIMME_INT nLocal, PRIMME_INT ldV,
      SCALAR *W, PRIMME_INT ldW, int blockSize, primme_params *primme) {

   int i, ONE=1, ierr=0;
   double t0;

   if (blockSize <= 0) return 0;
   assert(primme->nLocal == nLocal);

   t0 = primme_wTimer(0);

   if (primme->correctionParams.precondition) {
      if (primme->ldOPs == 0
            || (ldV == primme->ldOPs && ldW == primme->ldOPs)) {
         CHKERRM((primme->applyPreconditioner(V, &ldV, W, &ldW, &blockSize,
                     primme, &ierr), ierr), -1,
               "Error returned by 'applyPreconditioner' %d", ierr);
      }
      else {
         for (i=0; i<blockSize; i++) {
            CHKERRM((primme->applyPreconditioner(&V[ldV*i], &primme->ldOPs,
                        &W[ldW*i], &primme->ldOPs, &ONE, primme, &ierr), ierr),
                  -1, "Error returned by 'applyPreconditioner' %d", ierr);
         }
      }
      primme->stats.numPreconds += blockSize;
   }
   else {
      Num_copy_matrix_Sprimme(V, nLocal, blockSize, ldV, W, ldW);
   }

   primme->stats.timePrecond += primme_wTimer(0) - t0;

   return 0;
}

/*******************************************************************************
 * Subroutine convTestFun - wrapper around primme.convTestFun; evaluate if the
 *    the approximate eigenpair eval, evec with given residual norm is
 *    considered as converged.
 *
 * INPUT PARAMETERS
 * ----------------
 * eval     the eigenvalue
 * evec     the eigenvector
 * rNorm    the residual vector norm
 * 
 * OUTPUT
 * ------
 * isconv   if non-zero, the pair is considered converged.
 ******************************************************************************/

TEMPLATE_PLEASE
int convTestFun_Sprimme(REAL eval, SCALAR *evec, REAL rNorm, int *isconv, 
      struct primme_params *primme) {

   int ierr=0;
   double evald = eval, rNormd = rNorm;

   CHKERRM((primme->convTestFun(&evald, evec, &rNormd, isconv, primme, &ierr),
            ierr), -1, "Error returned by 'convTestFun' %d", ierr);

   return 0;
}
