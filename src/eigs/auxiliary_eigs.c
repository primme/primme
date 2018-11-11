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
 * File: auxiliary.c
 *
 * Purpose - Miscellanea functions used by PRIMME EIGS
 *
 ******************************************************************************/

#include <string.h> /* memset */
#include "const.h"
#include "numerical.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "auxiliary_eigs.h"
#endif

#ifdef USE_DOUBLE

/******************************************************************************
 * Function monitor_report - pass to the monitor the reports
 *
 * PARAMETERS
 * ---------------------------
 * fun      function name or message to report
 * time     time spent on the call
 * ctx      primme context
 *
 ******************************************************************************/

static int monitor_report(const char *fun, double time, primme_context ctx) {
   if (ctx.primme && ctx.primme->monitorFun) {
      int err;
      primme_event event =
            (time >= 0.0 ? primme_event_profile : primme_event_message);

#ifdef PRIMME_PROFILE
      /* Avoid profiling this function. It will turn out in a recursive call */
      ctx.path = NULL;
#endif

      CHKERRM((ctx.primme->monitorFun(NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, fun, &time, &event, ctx.primme, &err),
                    err),
            -1, "Error returned by monitorFun: %d", err);
   }
   return 0;
}

/******************************************************************************
 * Function primme_get_context - return a context from the primme_params
 *
 * PARAMETERS
 * ---------------------------
 * primme      primme_params struct
 *
 ******************************************************************************/

TEMPLATE_PLEASE
primme_context primme_get_context(primme_params *primme) {
   primme_context ctx;
   memset(&ctx, 0, sizeof(primme_context));
   if (primme) {
      ctx.primme = primme;
      ctx.printLevel = primme->printLevel;
      ctx.outputFile = primme->outputFile;
      ctx.numProcs = primme->numProcs;
      ctx.procID = primme->procID;
      ctx.mpicomm = primme->commInfo;
      ctx.queue = primme->queue;
      ctx.report = monitor_report;
#ifdef PRIMME_PROFILE
      if (primme->profile) {
         /* Compile regex. If there is no errors, set path to a nonzero       */
         /* value. Set ctx.report to the function that will channel the       */
         /* reports to the monitor. Report errors if they are.                */

         int ierr = regcomp(&ctx.profile, primme->profile, REG_NOSUB);
         if (ierr || MALLOC_PRIMME(1, &ctx.timeoff)) {
            char errmsg[100];
            regerror(ierr, &ctx.profile, errmsg, 100);
            if (ctx.report && ierr != 0) ctx.report(errmsg, -1, ctx);
            regfree(&ctx.profile);
            ctx.path = NULL;
         } else {
            *ctx.timeoff = 0.0;
            ctx.path = "";
         }
      } else {
         ctx.path = NULL;
      }
#endif
   }

   return ctx;
} 

/******************************************************************************
 * Function primme_free_context - free memory associated to the context
 *
 * PARAMETERS
 * ---------------------------
 * ctx         context
 *
 ******************************************************************************/

TEMPLATE_PLEASE
void primme_free_context(primme_context ctx) {

   /* Deregister the allocation of the current frame */

   primme_frame *curr = ctx.mm;
   Mem_deregister_alloc(curr, ctx);

   /* Pop the current frame */

   Mem_pop_frame(&ctx);

   /* Free the current frame */

   if (curr) free(curr);

   /* Free profiler */
#ifdef PRIMME_PROFILE
   if (ctx.path) regfree(&ctx.profile);
   if (ctx.timeoff) free(ctx.timeoff);
#endif
}

#endif /* USE_DOUBLE */

/******************************************************************************
 * Function Num_compute_residual - This subroutine performs the next operation
 *    in a cache-friendly way:
 *
 *    r = Ax - eval*Bx
 *
 * PARAMETERS
 * ---------------------------
 * n           The number of rows of x, Ax and r
 * eval        The value to compute the residual vector r
 * Bx          The vector Bx
 * Ax          The vector Ax
 * r           On output r = Ax - eval*Bx
 *
 ******************************************************************************/

TEMPLATE_PLEASE
void Num_compute_residual_Sprimme(PRIMME_INT n, HSCALAR eval, SCALAR *Bx,
   SCALAR *Ax, SCALAR *r, primme_context ctx) {

   int k, M=min(n,PRIMME_BLOCK_SIZE);

   for (k=0; k<n; k+=M, M=min(M,n-k)) {
      Num_copy_Sprimme(M, &Ax[k], 1, &r[k], 1, ctx);
      Num_axpy_Sprimme(M, -eval, &Bx[k], 1, &r[k], 1, ctx);
   }

}

/******************************************************************************
 * Function Num_update_VWXR - This subroutine performs the next operations:
 *
 *    X0 = V*h(nX0b+1:nX0e), X1 = V*h(nX1b+1:nX1e), X2 = V*h(nX2b+1:nX2e)
 *    Wo = W*h(nWob+1:nWoe),
 *    BX0 = BV*h(nBX0b+1:nBX0e), BX1 = BV*h(nBX1b+1:nBX1e)
 *    R = W*h(nRb+1:nRe) - BV*h(nRb+1:nRe)*diag(hVals(nRb+1:nRe)),
 *    Rnorms = norms(R),
 *    rnorms = norms(W*h(nrb+1:nre) - BV*h(nrb+1:nre)*diag(hVals(nrb+1:nre)))
 *    xnorms = norms(V*h(nxb+1:nxe))
 *    G = (V*h(1:nG))'*BV*h(1:nG)
 *    H = (V*h(1:nH))'*W*h(1:nH)
 *
 * NOTE: if Rnorms and rnorms are requested, nRb-nRe+nrb-nre < mV
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * V, W, BV    input basis
 * mV,nV,ldV   number of rows and columns and leading dimension of V, W and BV
 * h           input rotation matrix
 * nh          Number of columns of h
 * ldh         The leading dimension of h
 * hVals       Array of values
 *
 * OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------
 * X0             Output matrix V*h(nX0b:nX0e-1) (optional)
 * nX0b, nX0e     Range of columns of h
 * X1             Output matrix V*h(nX1b:nX1e-1) (optional)
 * nX1b, nX1e     Range of columns of h
 * X2             Output matrix V*h(nX2b:nX2e-1) (optional)
 * nX2b, nX2e     Range of columns of h
 * Wo             Output matrix W*h(nWob:nWoe-1) (optional)
 * nWob, nWoe     Range of columns of h
 * BX0            Output matrix BV*h(nBX0b:nBX0e-1) (optional)
 * nBX0b, nBX0e   Range of columns of h
 * BX1            Output matrix BV*h(nBX1b:nBX1e-1) (optional)
 * nBX1b, nBX1e   Range of columns of h
 * BX2            Output matrix BV*h(nBX2b:nBX2e-1) (optional)
 * nBX2b, nBX2e   Range of columns of h
 * R              Output matrix (optional)
 * nRb, nRe       Range of columns of h and hVals
 * Rnorms         Output array with the norms of R (optional)
 * rnorms         Output array with the extra residual vector norms (optional)
 * nrb, nre       Columns of residual vector to compute the norm
 * xnorms         Output array with V*h(nxb:nxe-1) vector norms (optional)
 * nxb, nxe       Columns of V*h to compute the norm
 * 
 * NOTE: n*e, n*b are zero-base indices of ranges where the first value is
 *       included and the last isn't.
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_update_VWXR_Sprimme(SCALAR *V, SCALAR *W, SCALAR *BV, PRIMME_INT mV,
      int nV, PRIMME_INT ldV, HSCALAR *h, int nh, int ldh, HREAL *hVals,
      SCALAR *X0, int nX0b, int nX0e, PRIMME_INT ldX0,
      SCALAR *X1, int nX1b, int nX1e, PRIMME_INT ldX1,
      SCALAR *X2, int nX2b, int nX2e, PRIMME_INT ldX2,
      SCALAR *Wo, int nWob, int nWoe, PRIMME_INT ldWo,
      SCALAR *R, int nRb, int nRe, PRIMME_INT ldR, HREAL *Rnorms,
      SCALAR *BX0, int nBX0b, int nBX0e, PRIMME_INT ldBX0,
      SCALAR *BX1, int nBX1b, int nBX1e, PRIMME_INT ldBX1,
      SCALAR *BX2, int nBX2b, int nBX2e, PRIMME_INT ldBX2,
      HREAL *rnorms, int nrb, int nre,
      HSCALAR *G, int nG, int ldG,
      HSCALAR *H, int nH, int ldH,
      HREAL *xnorms, int nxb, int nxe,
      primme_context ctx) {

   PRIMME_INT i;     /* Loop variables */
   int j;            /* Loop variables */
   int m=min(PRIMME_BLOCK_SIZE, mV);   /* Number of rows in the cache */
   int nXb, nXe, nYb, nYe, nBXb, nBXe, ldG0=0, ldH0=0;
   PRIMME_INT ldX, ldY, ldBX;
   SCALAR *X, *Y, *BX;
   HSCALAR *G0=NULL, *H0=NULL, *workGH = NULL;

   assert(mV <= ldV && nh <= ldh && (!G || nG <= ldG) && (!H || nH <= ldH));
   assert(Rnorms == NULL || nRb >= nRe || nrb >= nre || Rnorms != rnorms);

   /* Figure out which columns of V*h, W*h and BV*h to compute */

   nXb = nYb = nBXb = INT_MAX;
   nXe = nYe = nBXe = 0;

   if (X0 && nX0b < nX0e) nXb = min(nXb, nX0b), nXe = max(nXe, nX0e);
   if (X1 && nX1b < nX1e) nXb = min(nXb, nX1b), nXe = max(nXe, nX1e);
   if (X2 && nX2b < nX2e) nXb = min(nXb, nX2b), nXe = max(nXe, nX2e);
   if (R && nRb < nRe && !BV) nXb = min(nXb, nRb), nXe = max(nXe, nRe);
   if (G && nG > 0) nXb = min(nXb, 0), nXe = max(nXe, nG);
   if (H && nH > 0) nXb = min(nXb, 0), nXe = max(nXe, nH);
   if (rnorms && nrb < nre && !BV) nXb = min(nXb, nrb), nXe = max(nXe, nre);
   if (xnorms && nxb < nxe) nXb = min(nXb, nxb), nXe = max(nXe, nxe);
   nXe = max(nXe, nXb);

   if (Wo && nWob < nWoe) nYb = min(nYb, nWob), nYe = max(nYe, nWoe);
   if (R && nRb < nRe) nYb = min(nYb, nRb), nYe = max(nYe, nRe);
   if (rnorms && nrb < nre) nYb = min(nYb, nrb), nYe = max(nYe, nre);
   if (H && nH > 0) nYb = min(nYb, 0), nYe = max(nYe, nH);
   nYe = max(nYe, nYb);

   if (BV) {
      if (BX0 && nBX0b < nBX0e) nBXb = min(nBXb, nBX0b), nBXe = max(nBXe, nBX0e);
      if (BX1 && nBX1b < nBX1e) nBXb = min(nBXb, nBX1b), nBXe = max(nBXe, nBX1e);
      if (BX2 && nBX2b < nBX2e) nBXb = min(nBXb, nBX2b), nBXe = max(nBXe, nBX2e);
      if (R && nRb < nRe) nBXb = min(nBXb, nRb), nBXe = max(nBXe, nRe);
      if (G && nG > 0) nBXb = min(nBXb, 0), nBXe = max(nBXe, nG);
      if (rnorms && nrb < nre) nBXb = min(nBXb, nrb), nBXe = max(nBXe, nre);
   }
   nBXe = max(nBXe, nBXb);

   assert(nXe <= nh || nXb >= nXe); /* Check dimension */
   assert(nYe <= nh || nYb >= nYe); /* Check dimension */
   assert(nBXe <= nh || nBXb >= nXe); /* Check dimension */

   CHKERR(Num_malloc_Sprimme(m * (nXe - nXb), &X, ctx));
   CHKERR(Num_malloc_Sprimme(m * (nYe - nYb), &Y, ctx));
   CHKERR(Num_malloc_Sprimme(m * (nBXe - nBXb), &BX, ctx));
   ldX = ldY = ldBX = m;
   Num_zero_matrix_Sprimme(X, m, nXe - nXb, ldX, ctx);
   Num_zero_matrix_Sprimme(Y, m, nYe - nYb, ldY, ctx);
   Num_zero_matrix_Sprimme(BX, m, nBXe - nBXb, ldBX, ctx);

   int nGH = (G ? nG * nG : 0) + (H ? nH * nH : 0);
   if (ctx.numProcs > 1) {
      CHKERR(Num_malloc_SHprimme(nGH, &workGH, ctx));
      if (G) {
         G0 = workGH;
         ldG0 = nG;
      }
      if (H) {
         H0 = workGH + (G ? nG * nG : 0);
         ldH0 = nH;
      }
   }
   else {
      G0 = G; ldG0 = ldG;
      H0 = H; ldH0 = ldH;
   }

   if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = 0.0;
   if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = 0.0;
   if (G) Num_zero_matrix_SHprimme(G0, nG, nG, ldG0, ctx);
   if (H) Num_zero_matrix_SHprimme(H0, nH, nH, ldH0, ctx);
   if (xnorms) for (i=nxb; i<nxe; i++) xnorms[i-nxb] = 0.0;

   for (i=0; i < mV; i+=m, m=min(m,mV-i)) {
      /* X = V*h(nXb:nXe-1) */
      CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nXe-nXb, nV, 1.0,
         &V[i], ldV, &h[nXb*ldh], ldh, 0.0, X, ldX, ctx));

      /* X0 = X(nX0b-nXb:nX0e-nXb-1) */
      if (X0) Num_copy_matrix_Sprimme(&X[ldX*(nX0b-nXb)], m, nX0e-nX0b,
            ldX, &X0[i], ldX0, ctx);

      /* X1 = X(nX1b-nXb:nX1e-nXb-1) */
      if (X1) Num_copy_matrix_Sprimme(&X[ldX*(nX1b-nXb)], m, nX1e-nX1b,
            ldX, &X1[i], ldX1, ctx);

      /* X2 = X(nX2b-nXb:nX2e-nXb-1) */
      if (X2) Num_copy_matrix_Sprimme(&X[ldX*(nX2b-nXb)], m, nX2e-nX2b,
            ldX, &X2[i], ldX2, ctx);

      /* Y = W*h(nYb:nYe-1) */
      if (nYb < nYe) CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nYe-nYb, nV,
            1.0, &W[i], ldV, &h[nYb*ldh], ldh, 0.0, Y, ldY, ctx));

      /* Wo = Y(nWob-nYb:nWoe-nYb-1) */
      if (Wo) Num_copy_matrix_Sprimme(&Y[ldY*(nWob-nYb)], m, nWoe-nWob,
            ldY, &Wo[i], ldWo, ctx);

      /* BX = BV*h(nBXb:nBXe-1) */

      CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nBXe-nBXb, nV, 1.0,
         &BV[i], ldV, &h[nBXb*ldh], ldh, 0.0, BX, ldBX, ctx));

      /* BX0 = BX(nX0b-nXb:nX0e-nXb-1) */
      if (BX0) Num_copy_matrix_Sprimme(&BX[ldBX*(nBX0b-nBXb)], m, nBX0e-nBX0b,
            ldBX, &BX0[i], ldBX0, ctx);

      /* BX1 = BX(nBX1b-nBXb:nBX1e-nBXb-1) */
      if (BX1) Num_copy_matrix_Sprimme(&BX[ldBX*(nBX1b-nBXb)], m, nBX1e-nBX1b,
            ldBX, &BX1[i], ldBX1, ctx);

      /* BX2 = BX(nBX2b-nBXb:nBX2e-nBXb-1) */
      if (BX2) Num_copy_matrix_Sprimme(&BX[ldBX*(nBX2b-nBXb)], m, nBX2e-nBX2b,
            ldBX, &BX2[i], ldBX2, ctx);

      /* G += X(:,0:nG-1)'*X(:,0:nG-1) */

      if (G) {
         CHKERR(Num_gemm_ddh_Sprimme("C", "N", nG, nG, m, 1.0, X, ldX,
               BV ? BX : X, ldX, 1.0, G0, ldG0, ctx));
      }

      /* H = X(:,0:nH-1)'*Y(:,0:nH-1) */

      if (H) {
         CHKERR(Num_gemm_ddh_Sprimme(
               "C", "N", nH, nH, m, 1.0, X, ldX, Y, ldY, 1.0, H0, ldH0, ctx));
      }

      /* xnorms = norm(X(nxb-nXb:nxe-nXb-1)) */
      if (xnorms) for (j=nxb; j<nxe; j++) {
            xnorms[j - nxb] += REAL_PART(Num_dot_Sprimme(
                  m, &X[ldX * (j - nXb)], 1, &X[ldX * (j - nXb)], 1, ctx));
      }

      /* R = Y(nRb-nYb:nRe-nYb-1) - BX(nRb-nYb:nRe-nYb-1)*diag(nRb:nRe-1) */
      if (R) for (j=nRb; j<nRe; j++) {
            Num_compute_residual_Sprimme(m, hVals[j],
                  BV ? &BX[ldBX * (j - nBXb)] : &X[ldX * (j - nXb)],
                  &Y[ldY * (j - nYb)], &R[i + ldR * (j - nRb)], ctx);
            if (Rnorms) {
               Rnorms[j - nRb] +=
                     REAL_PART(Num_dot_Sprimme(m, &R[i + ldR * (j - nRb)], 1,
                           &R[i + ldR * (j - nRb)], 1, ctx));
         }
      }

      /* rnorms = Y(nrb-nYb:nre-nYb-1) - BX(nrb-nYb:nre-nYb-1)*diag(nrb:nre-1) */
      if (rnorms) for (j=nrb; j<nre; j++) {
            Num_compute_residual_Sprimme(m, hVals[j],
                  BV ? &BX[ldBX * (j - nBXb)] : &X[ldX * (j - nXb)],
                  &Y[ldY * (j - nYb)], &Y[ldY * (j - nYb)], ctx);
            rnorms[j - nrb] += REAL_PART(Num_dot_Sprimme(
                  m, &Y[ldY * (j - nYb)], 1, &Y[ldY * (j - nYb)], 1, ctx));
      }
   }

   /* Reduce and copy back G0 and H0 */

   if (ctx.numProcs > 1) {
      CHKERR(globalSum_SHprimme(workGH, workGH, nGH, ctx));
   }
   if (G) Num_copy_matrix_SHprimme(G0, nG, nG, ldG0, G, ldG, ctx);
   if (H) Num_copy_matrix_SHprimme(H0, nH, nH, ldH0, H, ldH, ctx);

   /* Reduce Rnorms, rnorms and xnorms and sqrt the results */

   if (ctx.numProcs > 1) {
      HREAL *tmp;
      CHKERR(Num_malloc_RHprimme(nRe - nRb + nre - nrb + nxe - nxb, &tmp, ctx));
      j = 0;
      if (Rnorms) for (i=nRb; i<nRe; i++) tmp[j++] = Rnorms[i-nRb];
      if (rnorms) for (i=nrb; i<nre; i++) tmp[j++] = rnorms[i-nrb];
      if (xnorms) for (i=nxb; i<nxe; i++) tmp[j++] = xnorms[i-nxb];
      if (j) CHKERR(globalSum_RHprimme(tmp, tmp, j, ctx));
      j = 0;
      if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = sqrt(tmp[j++]);
      if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = sqrt(tmp[j++]);
      if (xnorms) for (i=nxb; i<nxe; i++) xnorms[i-nxb] = sqrt(tmp[j++]);
      CHKERR(Num_free_RHprimme(tmp, ctx));
   }
   else {
      if (Rnorms) for (i=nRb; i<nRe; i++) Rnorms[i-nRb] = sqrt(Rnorms[i-nRb]);
      if (rnorms) for (i=nrb; i<nre; i++) rnorms[i-nrb] = sqrt(rnorms[i-nrb]);
      if (xnorms) for (i=nxb; i<nxe; i++) xnorms[i-nxb] = sqrt(xnorms[i-nxb]);
   }

   CHKERR(Num_free_Sprimme(X, ctx));
   CHKERR(Num_free_Sprimme(Y, ctx));
   CHKERR(Num_free_Sprimme(BX, ctx));
   CHKERR(Num_free_SHprimme(workGH, ctx));

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
      SCALAR *W, PRIMME_INT ldW, int blockSize, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, ONE=1, ierr=0;
   double t0;

   if (blockSize <= 0) return 0;
   assert(primme->nLocal == nLocal);

   t0 = primme_wTimer();

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
      Num_copy_matrix_Sprimme(V, nLocal, blockSize, ldV, W, ldW, ctx);
   }

   primme->stats.timePrecond += primme_wTimer() - t0;

   return 0;
}

/*******************************************************************************
 * Subroutine convTestFun - wrapper around primme.convTestFun; evaluate if the
 *    the approximate eigenpair eval, evec with given residual norm is
 *    considered as converged.
 *
 * INPUT PARAMETERS
 * ----------------
 * eval       the eigenvalue
 * evec       the eigenvector
 * givenEvec  whether eigenvector is provided
 * rNorm      the residual vector norm
 * 
 * OUTPUT
 * ------
 * isconv   if non-zero, the pair is considered converged.
 ******************************************************************************/

TEMPLATE_PLEASE
int convTestFun_Sprimme(HREAL eval, SCALAR *evec, int givenEvec, HREAL rNorm,
      int *isconv, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int ierr=0;
   double evald = eval, rNormd = rNorm;
   /* If an evec is going to be passed to convTestFun, but nLocal is 0,       */
   /* then fake the evec with a nonzero pointer in order to not be mistaken   */
   /* by not passing a vector.                                                */

   if (primme->nLocal == 0 && givenEvec) evec = (SCALAR *)0 + 1;

   CHKERRM((primme->convTestFun(&evald, givenEvec ? evec : NULL, &rNormd,
                  isconv, primme, &ierr),
                 ierr),
         -1, "Error returned by 'convTestFun' %d", ierr);

   return 0;
}

#ifdef USE_HOST

TEMPLATE_PLEASE
int globalSum_Sprimme(SCALAR *sendBuf, SCALAR *recvBuf, int count, 
   primme_context ctx) {

   primme_params *primme = ctx.primme;
   int ierr;
   double t0=0.0;

   if (primme && primme->globalSumReal) {
      t0 = primme_wTimer();

      /* If it is a complex type, count real and imaginary part */
#ifdef USE_COMPLEX
      count *= 2;
#endif
      CHKERRM((primme->globalSumReal(sendBuf, recvBuf, &count, primme, &ierr),
               ierr), PRIMME_USER_FAILURE,
            "Error returned by 'globalSumReal' %d", ierr);

      primme->stats.numGlobalSum++;
      primme->stats.timeGlobalSum += primme_wTimer() - t0;
      primme->stats.volumeGlobalSum += count;
   }
   else {
      Num_copy_Sprimme(count, sendBuf, 1, recvBuf, 1, ctx);
   }

   return 0;
}

#endif /* USE_HOST */

/*******************************************************************************
 * Subroutine problemNorm - return an estimation of |B\A|
 * 
 * INPUT PARAMETERS
 * ----------------
 * overrideUserEstimations    if nonzero, use estimations of |A| and |inv(B)| if
 *                            they are larger than primme.aNorm and primme.invBNorm 
 * 
 * OUTPUT
 * ------
 * return                     estimation of |B\A|
 ******************************************************************************/

TEMPLATE_PLEASE
HREAL problemNorm_Sprimme(
      int overrideUserEstimations, struct primme_params *primme) {

   if (!overrideUserEstimations) {
      if (!primme->massMatrixMatvec) {
         return primme->aNorm > 0.0 ? primme->aNorm
                                    : primme->stats.estimateLargestSVal;
      } else {
         return primme->aNorm > 0.0 && primme->invBNorm > 0.0
                      ? primme->aNorm * primme->invBNorm
                      : primme->stats.estimateLargestSVal;
      }
   }
   else {
      if (!primme->massMatrixMatvec) {
         return max(primme->aNorm > 0.0 ? primme->aNorm : 0.0,
               primme->stats.estimateLargestSVal);
      } else {
         return max(primme->aNorm > 0.0 && primme->invBNorm > 0.0
                          ? primme->aNorm * primme->invBNorm
                          : 0.0,
               primme->stats.estimateLargestSVal);
      }
   }
}

/*******************************************************************************
 * Subroutine deltaEig - return an estimation of the minimum distance that
 * that two distinct eigenvalues can have. We estimate that as the smallest
 * e so that the eigenpair (\lambda+e,x) has residual vector norm smaller than 
 * (|Ax| + max(|\lambda|)|Bx|)*\epsilon := |(A,B)|*\epsilon.
 *
 * If (\lambda,x) is an exact eigenpair, then |e| is constrain as
 *
 *    |Ax - (\lambda+e)Bx| = |e|*|Bx| <= |(A,B)|*\epsilon
 *
 * If x is B-normal, then x'*B*x = x'*chol(B)*chol(B)'*x = 1 and
 * |chol(B)'*x| = 1. Therefore
 *
 *    |Bx| = |chol(B)*chol(B)'x| >= minsval(chol(B)) * |chol(B)'*x|
 *                               >= sqrt(minsval(B))
 *
 * Therefore, |e| <= sqrt(|inv(B)|) * |(A,B)| * \epsilon
 * 
 * INPUT PARAMETERS
 * ----------------
 * overrideUserEstimations    if nonzero, use estimations of |A| and |B| if
 *                            they are larger than primme.aNorm and primme.BNorm 
 * 
 * OUTPUT
 * ------
 * return                     estimation of the minimum distance
 ******************************************************************************/

TEMPLATE_PLEASE
HREAL deltaEig_Sprimme(
      int overrideUserEstimations, struct primme_params *primme)
{
   HREAL BNorm;

   if (overrideUserEstimations) {
      BNorm = max(primme->BNorm, primme->stats.estimateBNorm);
   }
   else {
      BNorm =
            (primme->BNorm > 0.0 ? primme->BNorm : primme->stats.estimateBNorm);
   }

   return problemNorm_Sprimme(overrideUserEstimations, primme) /
          sqrt(BNorm) * MACHINE_EPSILON;
}

/*******************************************************************************
 * Function dist_dots - Computes several dot products in parallel
 *
 * Input Parameters
 * ----------------
 * x, y     Operands of the dot product operations
 *
 * ldx, ldy Leading dimension of x and y
 *
 * m        Length of the vectors x and y
 *
 * n        Number of columns in x and y
 *
 * ctx      Structure containing various solver parameters
 *
 * result   The inner products
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_dist_dots_Sprimme(SCALAR *x, PRIMME_INT ldx, SCALAR *y, PRIMME_INT ldy,
      PRIMME_INT m, int n, HSCALAR *result, primme_context ctx) {

   int i;
   for (i=0; i<n; i++) {
      result[i] = Num_dot_Sprimme(m, &x[ldx * i], 1, &y[ldy * i], 1, ctx);
   }
   CHKERR(globalSum_SHprimme(result, result, n, ctx));

   return 0;
}

/*******************************************************************************
 * Function dist_dots_real - Computes several dot products in parallel.
 *    Returns only the real part.
 *
 * Input Parameters
 * ----------------
 * x, y     Operands of the dot product operations
 *
 * ldx, ldy Leading dimension of x and y
 *
 * m        Length of the vectors x and y
 *
 * n        Number of columns in x and y
 *
 * ctx      Structure containing various solver parameters
 *
 * result   The inner products
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int Num_dist_dots_real_Sprimme(SCALAR *x, PRIMME_INT ldx, SCALAR *y,
      PRIMME_INT ldy, PRIMME_INT m, int n, HREAL *result, primme_context ctx) {

   int i;
   for (i=0; i<n; i++) {
      result[i] =
            REAL_PART(Num_dot_Sprimme(m, &x[ldx * i], 1, &y[ldy * i], 1, ctx));
   }
   CHKERR(globalSum_RHprimme(result, result, n, ctx));

   return 0;
}
