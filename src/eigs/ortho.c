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
 * File: ortho.c
 *
 * Purpose - Orthonormalizes a block of vectors, vector by vector, 
 *           against two bases and among themselves. Gram-Schmidt is used 
 *           with reorthogonalization based on Daniel's test. 
 *           For the purpose of the test, the norm of the resulting vector 
 *           is computed without synchronizations. Because of floating point,
 *           this norm is not accurate if reortho is needed, but it is more 
 *           than sufficient for the test. If at least one reortho has been
 *           performed, the actual norm must be computed. 
 *
 * Note on numerical properties
 *           The relative error in the implicit s1 is 
 *                    O( machineEps * (s0/s1)^2 ), 
 *           so when Daniel's test succeeds error(s1) = O(2*machineEps)).
 *
 *           For almost linearly dependent vectors a large error in s1
 *           might fail to detect a problem during the second reortho.
 *           Let s1 = s1hat*(1+err), with err = O(macheps*(s0/s1)^2).
 *              eg. s0/s1 = 1e8 => (s1-s1hat)/s1hat = O(1)
 *              eg. s0/s1 = 1e7 => (s1-s1hat)/s1hat = O(1e-2)
 *              eg. s0/s1 = 1e4 => (s1-s1hat)/s1hat = O(1e-8)
 *           When s0 is reset to s00 = s1, and perform again the Daniel's test
 *              s11/s00 = s11/s1hat(1+err) < 0.75 <=> s11/s1hat < .75*(1+err)
 *           if err is around 1e-2 it does not affect Daniels test a lot.
 *           Thus, if s0/s1 > 1.49e7, we must compute the actual s1.
 *
 * Note on linear dependency
 *           If the original vector norm is reduced more than machEps times,
 *           it is considered linear dependent. If R is not returned, the
 *           vector is replaced by a random vector. Otherwise the vector is
 *           zeroed.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "numerical.h"
#include "auxiliary_eigs.h"
#include "ortho.h"
#include "const.h"
#include "wtime.h"
 

/**********************************************************************
 * Function ortho - This routine orthonormalizes
 * a block of of vectors (from b1 to including b2 in basis)
 * against other vectors in the same array (from 0 to b1-1 in basis),
 * against a set of locked vectors (from 0 to numLocked-1 in locked),
 * and themselves.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * ldV        Leading dimension of the basis
 * b1, b2     Range of indices of vectors to be orthonormalized
 *            (b1 can be zero, but b1 must be <= b2)
 * ldR        Leading dimension in R
 * locked     Array that holds locked vectors if they are in-core
 * ldLocked   Leading dimension of locked
 * numLocked  Number of vectors in locked
 * nLocal     Number of rows of each vector stored on this node
 * machEps    Double machine precision
 *
 * rworkSize  Length of rwork array
 * primme     Primme struct. Contains globalSumDouble and Parallelism info
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * V       Basis vectors
 * R       Rotations done in the basis: V_input = V_output * R
 * iseed   Seeds used to generate random vectors
 * rwork   Contains buffers and other necessary work arrays
 *
 * Return Value
 * ------------
 * >0  - Insufficient workspace provided.  Positive number indicates
 *       appropriate amount.
 *  0  - success
 * -1  - some size or leading dimension was < 0
 * -2  - b1 > b2
 * -3  - A limit number of randomizations has been performed without
 *       yielding an orthogonal direction
 * 
 **********************************************************************/

static int Bortho_gen_Sprimme(SCALAR *V, PRIMME_INT ldV, SCALAR *R, int ldR,
      int b1, int b2, SCALAR *locked, PRIMME_INT ldLocked,
      int numLocked, PRIMME_INT nLocal, int (*B)(SCALAR*,PRIMME_INT,SCALAR*,
         PRIMME_INT,int,void*), void *Bctx,
      PRIMME_INT *iseed, primme_context ctx) {
              
   primme_params *primme = ctx.primme;
   int i, j;                /* Loop indices */
   int messages = 1;        /* messages = 1 prints the intermediate results */
   /* TODO: replace by a dynamic criterion when to stop orthogonalizing local */
   /* vectors. Observed performance improvement when maxNumOrthos increases.  */
   int maxNumOrthos = primme?3:7; /* We let 2 reorthogonalizations before randomize */
                                  /* for nLocal length vectors, and 6 orthogonalisations */
                                  /* for the rest */
   int maxNumRandoms = 10;  /* We do not allow more than 10 randomizations */
   double tol = sqrt(2.0L)/2.0L; /* We set Daniel et al. test to .707 */
   double t0;

   messages = (primme && primme->procID == 0 && primme->printLevel >= 3
         && primme->outputFile);

   /*----------------------------------*/
   /* input and workspace verification */
   /*----------------------------------*/
   assert(nLocal >= 0 && numLocked >= 0 &&
          ldV >= nLocal && (numLocked == 0 || ldLocked >= nLocal) &&
          (R == NULL || ldR > b2));

   tol = sqrt(2.0L)/2.0L;

   /* Zero the columns from b1 to b2 of R */
   if (R)
      for(i=b1; i <= b2; i++)
         for (j=0; j <= i; j++)
            R[ldR*i+j] = 0.0;

   /*---------------------------------------------------*/
   /* main loop to orthogonalize new vectors one by one */
   /*---------------------------------------------------*/

   t0 = primme_wTimer(0);

   // Allocate overlaps and Bx

   SCALAR *overlaps;
   CHKERR(Num_malloc_SHprimme(b2+1 + numLocked, &overlaps, ctx));
   SCALAR *Bx = NULL;
   CHKERR(Num_malloc_Sprimme(B?nLocal:0, &Bx, ctx));

   for(i=b1; i <= b2; i++) {
    
      int nOrth;  // number of orthogonalizations of the current vector
      int randomizations;  // times the current vector has been replaced by a
                           // random vector
      int reorth = 1;   // flag to keep iterating
      int updateR = (R ? 1 : 0); // flag to keep updating R, after the first
      int Bx_update = 0;   // flag indicating if Bx = B*V[i]
                                 // randomization it is set to zero
      REAL s0=0.0;   // B-norm of the current vector before deflating V and locked
      REAL s02=0.0;  // s0 squared
      REAL s1=0.0;   // B-norm of the current vector after deflating V and locked
      REAL s12=0.0;  // s1 squared

      for (nOrth=0, randomizations=0; reorth; ) {

         if (nOrth >= maxNumOrthos) {
            /* Stop updating R when replacing one of the columns of the V */
            /* with a random vector                                           */

            if (updateR) {
               R[ldR*i + i] = 0.0;
               updateR = 0;
            }

            if (randomizations >= maxNumRandoms) {
               return -3;
            }
            if (messages){
               fprintf(primme->outputFile, "Randomizing in ortho: %d, vector size of %" PRIMME_INT_P "\n", i, nLocal);
            }

            Num_larnv_Sprimme(2, iseed, nLocal, &V[ldV*i], ctx); 
            randomizations++;
            nOrth = 0;
            Bx_update = 0;    // V[i] has changed, so Bx != B*V[i]
         }

         nOrth++;

         // Compute B*V[i]

         if (B && !Bx_update) {
            CHKERR(B(&V[ldV*i], ldV, Bx, nLocal, 1, Bctx));
         }
         else if (!B) {
            Bx = &V[ldV*i];
         }

         // Compute the B norm of the current vector, V[i], if it wasn't computed
         // in previous iteration

         if (nOrth == 1) {
            s02 = REAL_PART(Num_dot_Sprimme(nLocal, &V[ldV*i], 1, Bx, 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += 1;
         }

         // Compute overlaps = [V[0:i-1]'*B*V[i] 
 
         if (i > 0) {
            CHKERR(Num_gemv_ddh_Sprimme(
                  "C", nLocal, i, 1.0, V, ldV, Bx, 1, 0.0, overlaps, 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += i;
         }

         if (numLocked > 0) {
            CHKERR(Num_gemv_ddh_Sprimme("C", nLocal, numLocked, 1.0, locked,
                  ldLocked, Bx, 1, 0.0, &overlaps[i], 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += numLocked;
         }

         overlaps[i+numLocked] = s02;
         CHKERR(globalSum_SHprimme(overlaps, overlaps, i + numLocked + 1,
                  ctx));

         if (updateR) {
             Num_axpy_SHprimme(i, 1.0, overlaps, 1, &R[ldR*i], 1, ctx);
         }

         if (numLocked > 0) { /* locked array most recently accessed */
            // Compute V[i] = V[i] - locked'*overlaps[i:i+numLocked-1]
            CHKERR(Num_gemv_dhd_Sprimme("N", nLocal, numLocked, -1.0, locked,
                  ldLocked, &overlaps[i], 1, 1.0, &V[ldV * i], 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += numLocked;
         }

         if (i > 0) {
            CHKERR(Num_gemv_dhd_Sprimme("N", nLocal, i, -1.0, V, ldV, overlaps,
                  1, 1.0, &V[ldV * i], 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += i;
         }

         Bx_update = 0;    // V[i] has changed, so Bx != B*V[i]
 
         if (nOrth == 1) {
            s0 = sqrt(s02 = REAL_PART(overlaps[i+numLocked]));
         }

         /* Compute the norm of the resulting vector implicitly */
         
         {
            REAL temp = REAL_PART(Num_dot_SHprimme(i+numLocked,overlaps,1,overlaps,1,ctx));
            s1 = sqrt(s12 = max(0.0L, s02-temp));
         }
         
         /* If s1 decreased too much, its implicit computation may have       */
         /* problem. Compute s1 explicitly in that cases                      */
         
         int s1_update = 0;      // flag if s1 has been computed explicitly
         if ( s1 < s0*sqrt(MACHINE_EPSILON) || nOrth > 1 || !primme) {  
            if (B) {
               CHKERR(B(&V[ldV*i], ldV, Bx, nLocal, 1, Bctx));
               Bx_update = 1;
            }
            REAL temp =
                REAL_PART(Num_dot_Sprimme(nLocal, &V[ldV * i], 1, Bx, 1, ctx));
            if (primme) primme->stats.numOrthoInnerProds += 1;
            CHKERR(globalSum_RHprimme(&temp, &s12, 1, ctx));
            s1 = sqrt(s12);
            s1_update = 1;
         }

         if (s1 <= MACHINE_EPSILON*s0) {
            if (messages) {
               fprintf(primme->outputFile, 
                 "Vector %d lost all significant digits in ortho\n", i-b1);
            }
            nOrth = maxNumOrthos;
         }
         else if (s1 <= tol*s0 || (!primme && nOrth < maxNumOrthos)) {
            /* No numerical benefit in normalizing the vector before reortho */
            s0 = s1;
            s02 = s12;
         }
         else {
            if (updateR) {
               if (!s1_update) {
                  if (B && !Bx_update) {
                     CHKERR(B(&V[ldV*i], ldV, Bx, nLocal, 1, Bctx));
                  }
                  REAL temp = REAL_PART(Num_dot_Sprimme(nLocal,
                           &V[ldV*i], 1, Bx, 1, ctx));
                  if (primme) primme->stats.numOrthoInnerProds += 1;
                  CHKERR(globalSum_RHprimme(&temp, &s1, 1, ctx));
                  s1 = sqrt(max((REAL)0, s1));
               }
               R[ldR*i + i] = s1;
            }

            if (ISFINITE((REAL)(1.0/s1))) {
               Num_scal_Sprimme(nLocal, 1.0/s1, &V[ldV*i], 1, ctx);
               break;
            }
            else {
               if (messages) {
                  fprintf(primme->outputFile, 
                        "Vector %d lost all significant digits in ortho\n", i-b1);
               }
               nOrth = maxNumOrthos;
            }
         } 
      }
   }

   if (primme) primme->stats.timeOrtho += primme_wTimer(0) - t0;

   CHKERR(Num_free_SHprimme(overlaps, ctx));
   if (B) CHKERR(Num_free_Sprimme(Bx, ctx));

   /* Check orthogonality */
   /*
   if (numLocked) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*numLocked*numLocked);
      Num_gemm_Sprimme("C", "N", numLocked, numLocked, nLocal, 1.0, locked,
            ldLocked, locked, ldLocked, 0.0, H, numLocked);
      for(i=0; i < numLocked; i++) {
         for(j=0; j < i; j++) assert(ABS(H[numLocked*i+j]) < 1e-13);
         assert(fabs(1 - ABS(H[numLocked*i+i])) < 1e-13);
      }
      free(H);
   }
   if (b2+1) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*(b2+1)*(b2+1));
      Num_gemm_Sprimme("C", "N", b2+1, b2+1, nLocal, 1.0, V,
            ldV, V, ldV, 0.0, H, b2+1);
      for(i=0; i < b2+1; i++) {
         for(j=0; j < i; j++) assert(ABS(H[(b2+1)*i+j]) < 1e-13);
         assert(H[(b2+1)*i+i] == 0.0 || fabs(1 - ABS(H[(b2+1)*i+i])) < 1e-13);
      }
      free(H);
   }
   if (numLocked) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*(b2+1)*numLocked);
      Num_gemm_Sprimme("C", "N", numLocked, b2+1, nLocal, 1.0, locked,
            ldLocked, V, ldV, 0.0, H, numLocked);
      for(i=0; i < b2+1; i++) {
         for(j=0; j < numLocked; j++) assert(ABS(H[numLocked*i+j]) < 1e-13);
      }
      free(H);
   }
   */

   return 0;
}

TEMPLATE_PLEASE
int ortho_Sprimme(SCALAR *V, PRIMME_INT ldV, SCALAR *R, int ldR, int b1, int b2,
                  SCALAR *locked, PRIMME_INT ldLocked, int numLocked,
                  PRIMME_INT nLocal, PRIMME_INT *iseed, primme_context ctx) {

  return Bortho_gen_Sprimme(V, ldV, R, ldR, b1, b2, locked, ldLocked, numLocked,
                            nLocal, NULL, NULL, iseed, ctx);
}

#ifdef USE_HOST
struct local_matvec_ctx { SCALAR *B; int n, ldB; };

static int local_matvec(SCALAR *x, PRIMME_INT ldx, SCALAR *y, PRIMME_INT ldy ,int bs,
      void *ctx_)
{
   struct local_matvec_ctx *ctx = (struct local_matvec_ctx*)ctx_;
   Num_hemm_SHprimme(
         "L", "U", ctx->n, bs, 1.0, ctx->B, ctx->ldB, x, ldx, 0.0, y, ldy);
   return 0;
}

TEMPLATE_PLEASE
int Bortho_local_Sprimme(SCALAR *V, int ldV, SCALAR *R,
      int ldR, int b1, int b2, SCALAR *locked, int ldLocked,
      int numLocked, PRIMME_INT nLocal, SCALAR *B, int ldB, PRIMME_INT *iseed,
      primme_context ctx) {

   (void)ctx; 
   struct local_matvec_ctx Bctx = {B, (int)nLocal, ldB};
   return Bortho_gen_Sprimme(V, ldV, R, ldR, b1, b2, locked, ldLocked,
                             numLocked, nLocal, B ? local_matvec : NULL, &Bctx,
                             iseed, primme_get_context(NULL));
}
#endif /* USE_HOST */

/**********************************************************************
 * Function ortho_single_iteration -- This function orthogonalizes
 *    applies ones the projector (I-QQ') on X. Optionally returns
 *    the norms ||Q'X(i)|| and ||(I-QQ')X(i)||.
 *   
 * ARRAYS AND PARAMETERS
 * ----------------
 * Q               The basis of the projector I-QQ'.
 * mQ, nQ, ldQ     Rows, columns and leading dimension of Q.
 * X               The vectors to apply the projector I-QQ'.
 * inX             Column indices to apply the projector (optional).
 * nX              Number of columns to apply the projector.
 * ldX             The leading dimension of X.
 * overlaps        The norms of Q'X(i) (optional).
 * norms           The norms of (I-QQ')X(i) (optional).
 * rwork           Auxiliary space
 * lrwork          Available rwork
 *
 ****************************************************************************/

TEMPLATE_PLEASE
int ortho_single_iteration_Sprimme(SCALAR *Q, PRIMME_INT mQ, PRIMME_INT nQ,
      PRIMME_INT ldQ, SCALAR *X, int *inX, int nX, PRIMME_INT ldX,
      REAL *overlaps, REAL *norms, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j, M=PRIMME_BLOCK_SIZE, m=min(M, mQ);
   SCALAR *y, *y0, *X0;

   double t0 = primme_wTimer(0);

   CHKERR(Num_malloc_SHprimme(nQ*nX, &y, ctx));
   CHKERR(Num_malloc_SHprimme(nQ*nX, &y0, ctx));
   CHKERR(Num_malloc_Sprimme(m*nX, &X0, ctx));

   /* Check if the indices of inX are contiguous */

   if (inX && nX > 0) {
      for (i=0, j=inX[0]; i<nX && j==inX[i]; i++, j++);
      if (i >= nX) {
         X = &X[inX[0]*ldX];
         inX = NULL;
      }
   }

   /* y = Q'*X */
   if (!inX) {
      CHKERR(Num_gemm_ddh_Sprimme(
            "C", "N", nQ, nX, mQ, 1.0, Q, ldQ, X, ldX, 0.0, y, nQ, ctx));
   }
   else {
      Num_zero_matrix_SHprimme(y, nQ, nX, nQ, ctx);
      for (i=0, m=min(M,mQ); i < mQ; i+=m, m=min(m,mQ-i)) {
        Num_copy_matrix_columns_Sprimme(&X[i], m, inX, nX, ldX, X0, NULL, m,
                                        ctx);
        CHKERR(Num_gemm_ddh_Sprimme(
              "C", "N", nQ, nX, m, 1.0, &Q[i], ldQ, X0, m, 1.0, y, nQ, ctx));
      }
   }
   primme->stats.numOrthoInnerProds += nQ*nX;

   /* Store the reduction of y in y0 */
   CHKERR(globalSum_SHprimme(y, y0, nQ*nX, ctx));
   
   /* overlaps(i) = norm(y0(:,i))^2 */
   for (i=0; i<nX; i++) {
      overlaps[i] =
         sqrt(REAL_PART(Num_dot_SHprimme(nQ, &y0[nQ*i], 1, &y0[nQ*i], 1, ctx)));
   }

   /* X = X - Q*y0; norms(i) = norms(X(i))^2 */
   if (norms) for (i=0; i<nX; i++) norms[i] = 0.0;
   for (i=0, m=min(M,mQ); i < mQ; i+=m, m=min(m,mQ-i)) {
      if (inX) {
        Num_copy_matrix_columns_Sprimme(&X[i], m, inX, nX, ldX, X0, NULL, m,
                                        ctx);
      }
      CHKERR(Num_gemm_dhd_Sprimme("N", "N", m, nX, nQ, -1.0, &Q[i], ldQ, y0, nQ,
            1.0, inX ? X0 : &X[i], inX ? m : ldX, ctx));
      if (inX) {
         Num_copy_matrix_columns_Sprimme(X0, m, NULL, nX, m, &X[i], inX, ldX, ctx);
      }
      if (norms) for (j=0; j<nX; j++) {
         SCALAR *v = inX ? &X0[j*m] : &X[j*ldX+i];
         norms[j] += REAL_PART(Num_dot_Sprimme(m, v, 1, v, 1, ctx));
      }
   }

   if (norms) {
      /* Store the reduction of norms */
      CHKERR(globalSum_RHprimme(norms, norms, nX, ctx));
 
      for (i=0; i<nX; i++) norms[i] = sqrt(norms[i]);
      primme->stats.numOrthoInnerProds += nX;
   }

   CHKERR(Num_free_SHprimme(y, ctx));
   CHKERR(Num_free_SHprimme(y0, ctx));
   CHKERR(Num_free_Sprimme(X0, ctx));

   primme->stats.timeOrtho += primme_wTimer(0) - t0;

   return 0;
}
