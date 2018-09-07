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
#include "const.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "factorize.h"
#include "auxiliary_eigs.h"
#include "ortho.h"
#include "wtime.h"
#endif
 
static int Num_ortho_kernel(SCALAR *Q, PRIMME_INT M, int nQ, PRIMME_INT ldQ,
      SCALAR *V, int nV, PRIMME_INT ldV, SCALAR *X, int nX, PRIMME_INT ldX,
      HSCALAR *A, int ldA, HREAL *D, HSCALAR *Y, int ldY,
      SCALAR *W, int nW, PRIMME_INT ldW, SCALAR *Z, int nZ, PRIMME_INT ldZ,
      HSCALAR *B, int ldB, primme_context ctx);

static int eig(HSCALAR *H, int n, int ldH, HSCALAR *Y, int ldhVecs, HREAL *svals,
      primme_context ctx);
static int rank_estimation(HSCALAR *V, int n0, int n1, int n, int ldV);


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
 * RLocked    R for locked vectors
 * ldRLocked  The leading dimension of RLocked
 * nLocal     Number of rows of each vector stored on this node
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

static int Bortho_gen_Sprimme(SCALAR *V, PRIMME_INT ldV, HSCALAR *R, int ldR,
      int b1, int b2, SCALAR *locked, PRIMME_INT ldLocked, int numLocked,
      HSCALAR *RLocked, int ldRLocked, PRIMME_INT nLocal,
      int (*B)(SCALAR *, PRIMME_INT, SCALAR *, PRIMME_INT, int, void *),
      void *Bctx, PRIMME_INT *iseed, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i;                   /* Loop indices */
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
   if (R) Num_zero_matrix_SHprimme(&R[b1 * ldR], b2 + 1, b2 - b1 + 1, ldR, ctx);
   if (RLocked) {
      Num_zero_matrix_SHprimme(RLocked, numLocked, b2 - b1 + 1, ldRLocked, ctx);
   }

   /*---------------------------------------------------*/
   /* main loop to orthogonalize new vectors one by one */
   /*---------------------------------------------------*/

   t0 = primme_wTimer(0);

   // Allocate overlaps and Bx

   HSCALAR *overlaps;
   CHKERR(Num_malloc_SHprimme(b2+1 + numLocked, &overlaps, ctx));
   SCALAR *Bx = NULL;
   CHKERR(Num_malloc_Sprimme(B?nLocal:0, &Bx, ctx));

   for(i=b1; i <= b2; i++) {
    
      int nOrth;  // number of orthogonalizations of the current vector
      int randomizations;  // times the current vector has been replaced by a
                           // random vector
      int reorth = 1;   // flag to keep iterating
      int updateR =
            (R || RLocked ? 1 : 0); // flag to keep updating R, after the first
      int Bx_update = 0;   // flag indicating if Bx = B*V[i]
                                 // randomization it is set to zero
      HREAL s0=0.0;   // B-norm of the current vector before deflating V and locked
      HREAL s02=0.0;  // s0 squared
      HREAL s1=0.0;   // B-norm of the current vector after deflating V and locked
      HREAL s12=0.0;  // s1 squared

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
             if (R) Num_axpy_SHprimme(i, 1.0, overlaps, 1, &R[ldR*i], 1, ctx);
             if (RLocked) {
                Num_axpy_SHprimme(numLocked, 1.0, overlaps + i, 1,
                      &RLocked[ldRLocked * (i - b1)], 1, ctx);
             }
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
            HREAL temp = REAL_PART(Num_dot_SHprimme(i+numLocked,overlaps,1,overlaps,1,ctx));
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
            HREAL temp =
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
            if (updateR && R) {
               if (!s1_update) {
                  if (B && !Bx_update) {
                     CHKERR(B(&V[ldV*i], ldV, Bx, nLocal, 1, Bctx));
                  }
                  HREAL temp = REAL_PART(Num_dot_Sprimme(nLocal,
                           &V[ldV*i], 1, Bx, 1, ctx));
                  if (primme) primme->stats.numOrthoInnerProds += 1;
                  CHKERR(globalSum_RHprimme(&temp, &s1, 1, ctx));
                  s1 = sqrt(max((HREAL)0, s1));
               }
               R[ldR*i + i] = s1;
            }

            if (ISFINITE((HREAL)(1.0/s1))) {
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
int ortho_Sprimme(SCALAR *V, PRIMME_INT ldV, HSCALAR *R, int ldR, int b1, int b2,
                  SCALAR *locked, PRIMME_INT ldLocked, int numLocked,
                  PRIMME_INT nLocal, PRIMME_INT *iseed, primme_context ctx) {

  return Bortho_gen_Sprimme(V, ldV, R, ldR, b1, b2, locked, ldLocked, numLocked,
                            NULL, 0, nLocal, NULL, NULL, iseed, ctx);
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

   /* Remove MPI communications from the context */

   ctx.numProcs = 1;
   ctx.procID = 0;
   ctx.mpicomm = NULL;
   ctx.primme = NULL;

   /* Call orthogonalization */

   struct local_matvec_ctx Bctx = {B, (int)nLocal, ldB};
   return Bortho_gen_Sprimme(V, ldV, R, ldR, b1, b2, locked, ldLocked,
         numLocked, NULL, 0, nLocal, B ? local_matvec : NULL, &Bctx, iseed,
         ctx);
}

#endif /* USE_HOST */

/**********************************************************************
 * Function ortho_block - This routine orthonormalizes
 * a block of of vectors (from b1 to including b2 in basis)
 * against other vectors in the same array (from 0 to b1-1 in basis),
 * against a set of locked vectors (from 0 to numLocked-1 in locked),
 * and themselves.
 *
 * The following conditions must always be met: 
 * ldBasis > 0, nLocal > 0, b1 >= 0, b2 >= 0, b2 >= b1, numLocked >= 0, 
 * rworkSize > 0
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
 * maxRank    largest rank of the basis being orthogonalized
 *
 * ctx        primme context
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * V       Basis vectors stored in core memory
 * R       Rotations done in the basis: input_basis = output_basis * R
 * b2_out  The number of linear independent columns
 *
 * Return Value
 * ------------
 *  error code
 * 
 **********************************************************************/

TEMPLATE_PLEASE
int ortho_block_Sprimme(SCALAR *V, PRIMME_INT ldV, HSCALAR *VLtVL, int ldVLtVL,
      HSCALAR *R, PRIMME_INT ldR, int b1, int b2, SCALAR *locked,
      PRIMME_INT ldLocked, int numLocked, HSCALAR *RLocked, int ldRLocked,
      PRIMME_INT nLocal, int maxRank, int *b2_out, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j;               /* loop indices */
   HSCALAR *A, *C, *Y;      /* auxiliary local matrices */
   HSCALAR *VLtVLdA;        /* auxiliary local matrices */
   HSCALAR *fVLtVL;         /* auxiliary local matrices */
   HREAL *D, *N;            /* singular values */
   int ldA;                /* leading dimension of A */
   b2++; /* NOTE: Let's use C range convention */

   /* Quick exit */

   if (b2 <= b1) {
      *b2_out = b2;
      return 0;
   }

   if (VLtVL == NULL) {
      CHKERR(Bortho_gen_Sprimme(V, ldV, R, ldR, b1, b2 - 1, locked, ldLocked,
            numLocked, RLocked, ldRLocked, nLocal, NULL, NULL, primme->iseed,
            ctx));
      *b2_out = b2;
      return 0;
   }

   // TEMP!!!
   // CHKERR(Bortho_gen_Sprimme(V, ldV, R, ldR, b1, b2 - 1, locked, ldLocked,
   //          numLocked, RLocked, ldRLocked, nLocal, NULL, NULL, primme->iseed,
   //          ctx));
   // Num_zero_matrix_SHprimme(&VLtVL[ldVLtVL*(numLocked+b1)], numLocked+b2, b2-b1, ldVLtVL, ctx);
   // for (i=numLocked+b1; i<numLocked+b2; i++) {
   //    VLtVL[ldVLtVL*i+i] = 1.0;
   // }
   // *b2_out = b2;
   // return 0;
   
 
   /* input and workspace verification */

   assert(nLocal >= 0 && numLocked >= 0 &&
          ldV >= nLocal && (numLocked == 0 || ldLocked >= nLocal) &&
          (R == NULL || ldR >= b2) && (!R || numLocked == 0));

   /* Zero the columns from b1 to b2 of R */

   HSCALAR *r = NULL; /* The diagonal block of R(b1:b2) */
   int ldr = 0;

   if (R) {
      Num_zero_matrix_SHprimme(&R[ldR * b1], b1, b2 - b1, ldR, ctx);
      r = &R[ldR * b1 + b1];
      ldr = ldR;
   }

   if (RLocked) {
      Num_zero_matrix_SHprimme(RLocked, numLocked, b2 - b1, ldRLocked, ctx);
      if (!r) {
         CHKERR(Num_malloc_SHprimme((b2 - b1) * (b2 - b1), &r, ctx));
         ldr = b2 - b1;
      }
   }
   if (r) {
      Num_zero_matrix_SHprimme(r, b2 - b1, b2 - b1, ldr, ctx);
      for (i=0; i<b2 - b1; i++) {
         r[ldr*i+i] = 1.0;
      }
   }

   /* Allocate workspace */

   CHKERR(Num_malloc_RHprimme(b2-b1, &D, ctx));
   CHKERR(Num_malloc_RHprimme(b2-b1, &N, ctx));
   A = &VLtVL[ldVLtVL*(b1+numLocked)];
   ldA = ldVLtVL;
   CHKERR(Num_malloc_SHprimme((numLocked+b2)*b2, &VLtVLdA, ctx));
   CHKERR(Num_malloc_SHprimme((b2-b1)*(b2-b1), &Y, ctx));
   int nVL = b1+numLocked;
   CHKERR(Num_malloc_SHprimme(nVL*nVL, &fVLtVL, ctx));
   int *pVLtVL;
   CHKERR(Num_malloc_iprimme(nVL, &pVLtVL, ctx));
   CHKERR(Num_malloc_SHprimme((b2-b1)*(b2-b1), &C, ctx));

   /* Factor VLtVL */

   CHKERR(UDUDecompose_SHprimme(VLtVL, ldVLtVL, fVLtVL, nVL, pVLtVL, nVL, ctx));

   /* Main loop to orthogonalize new vectors one by one. Just kidding. Although
    * seeming complicated it is just doing iterative CG-SVQB */

   double t0 = primme_wTimer(0);

   *b2_out = b2; 
   int its;
   int maxits = 5, plus2 = 5;
   for (its=0; its<maxits; its++) {
      /* Notation:                                          */
      /* Vp = [locked V(0:b1-1)]; Xp = V(b1:b2-1) */
      /* Xc = V(b1:b2-1)                                    */
      /* Xn is Xp in the next iteration                     */

      /* Do Xp <= (Xp - Vp*VLtVLdA)/N*Y/D if it isn't the first iteration     */
      /* Do A <= [locked V(0:b2)]'*Xc always                                  */

      CHKERR(Num_ortho_kernel(locked, nLocal, numLocked, ldLocked, V, b1, ldV,
            &V[b1 * ldV], b2 - b1, ldV, its == 0 ? NULL : VLtVLdA, nVL, D, Y,
            b2 - b1, V, b2, ldV, &V[ldV * b1], b2 - b1, ldV, A, ldA, ctx));
      if (primme) primme->stats.numOrthoInnerProds += (numLocked+b2)*(b2-b1);

      /* Check convergence */

      if (rank_estimation(VLtVL, numLocked + b1,
                   numLocked + b2, maxRank, ldVLtVL) == numLocked + b2) {
         if (its >= plus2 - 1) break;
         else plus2 = min(its + 2, plus2);
      }

      if (ctx.procID == 0) {
         /* C = Xn'*Xn = Xc'*(I-Vc/(Vc'*Vc)*Vc')*Xc = Xc'*Xc -
            (Xc'*Vc)*((Vc'*Vc)\(Vc'*Xc))       */
         /*   = A(b1:b2-1,:) -
            (A(0:numLocked+b1-1,:)'*(VLtVL\A(0:numLocked+b1-1,:))   */

         Num_copy_matrix_SHprimme(
               &A[numLocked + b1], b2 - b1, b2 - b1, ldA, C, b2 - b1, ctx);
         CHKERR(UDUSolve_SHprimme(
               fVLtVL, pVLtVL, nVL, A, b2 - b1, ldA, VLtVLdA, nVL, ctx));
         Num_gemm_SHprimme("C", "N", b2 - b1, b2 - b1, numLocked + b1, -1.0, A,
               ldA, VLtVLdA, nVL, 1.0, C, b2 - b1, ctx);

         /* N(i) = ||Xc(:,i)||; normXc = max ||N(i)|| */

         for (i=0; i<b2-b1; i++) {
            N[i] = sqrt(max(ABS(C[(b2-b1)*i+i]), MACHINE_EPSILON));
         }
         for (i=0; i<b2-b1; i++) {
            for (j=0; j<=i; j++) {
               C[(b2 - b1) * i + j] =
                     (C[(b2 - b1) * i + j] + CONJ(C[(b2 - b1) * j + i])) /
                     (HSCALAR)2.0 / N[i] / N[j];
            }
         }

         /* [D, Y] = eig(C) */

         CHKERR(eig(C, b2 - b1, b2 - b1, &Y[(b2 - b1) * (b1 - b1) + b1 - b1],
               b2 - b1, D, ctx));

         /* D = sqrt(D) */

         for (i=0; i<b2-b1; i++) {
            D[i] = sqrt(max(D[i], MACHINE_EPSILON * (b2 - b1)));
         }
      } else {
         Num_zero_matrix_SHprimme(VLtVLdA, nVL, b2-b1, nVL, ctx);
         for (i=0; i<b2-b1; i++) D[i] = 0.0;
         for (i=0; i<b2-b1; i++) N[i] = 0.0;
         Num_zero_matrix_SHprimme(&Y[(b2 - b1) * (b1 - b1) + b1 - b1], b2 - b1,
               b2 - b1, b2 - b1, ctx);
      }

      CHKERR(globalSum_SHprimme(VLtVLdA, VLtVLdA, nVL*(b2-b1), ctx));
      CHKERR(globalSum_RHprimme(N, N, b2-b1, ctx));
      CHKERR(globalSum_RHprimme(D, D, b2-b1, ctx));
      CHKERR(globalSum_SHprimme(Y, Y, (b2 - b1) * (b2 - b1), ctx));

      if (RLocked) {
         /* R(0:b1-1,:) += VtBV\Vc'*Xc*gi = VtBV\A(0:b1-1,:)*R(b1:b2-1,:) */

         Num_gemm_SHprimme("N", "N", numLocked, b2 - b1, b2 - b1, 1.0, VLtVLdA,
               nVL, r, ldr, 1.0, RLocked, ldRLocked, ctx);
      }
      if (R) {
         /* R(0:b1-1,:) += VtBV\Vc'*Xc*gi = VtBV\A(0:b1-1,:)*R(b1:b2-1,:) */

         Num_gemm_SHprimme("N", "N", b1, b2 - b1, b2 - b1, 1.0,
               &VLtVLdA[numLocked], nVL, &R[ldR * b1 + b1], ldR, 1.0,
               &R[ldR * b1], ldR, ctx);
      }
      if (r) {
         /* R(b1:b2-1,b1:b2-1) = N * R(b1:b2-1,b1:b2-1) */

         for (i = 0; i < b2 - b1; i++) {
            for (j = 0; j < b2 - b1; j++) { r[ldr * i + j] *= N[j]; }
         }

         /* C = Y' * R(b1:b2-1,b1:b2-1) */

         Num_gemm_SHprimme("C", "N", b2 - b1, b2 - b1, b2 - b1, 1.0, Y, b2 - b1,
               r, ldr, 0.0, C, b2 - b1, ctx);

         /* R(b1:b2-1,b1:b2-1) = D * C */

         for (i = 0; i < b2 - b1; i++) {
            for (j = 0; j < b2 - b1; j++) {
               r[ldr * i + j] = D[j] * C[(b2 - b1) * i + j];
            }
         }
      }

      /* Y = N \ Y */

      for (i=0; i<b2-b1; i++) {
         for (j=0; j<b2-b1; j++) {
            Y[(b2-b1)*i+j] /= N[j];
         }
      }
   } /* end while, I hope you enjoyed the loop */

   if (primme) primme->stats.timeOrtho += primme_wTimer(0) - t0;

   *b2_out = rank_estimation(VLtVL, numLocked + b1, numLocked + b2,
                        maxRank, ldVLtVL) - numLocked;
   
   if (!R && RLocked) CHKERR(Num_free_SHprimme(r, ctx));
   CHKERR(Num_free_RHprimme(D, ctx));
   CHKERR(Num_free_RHprimme(N, ctx));
   CHKERR(Num_free_SHprimme(VLtVLdA, ctx));
   CHKERR(Num_free_SHprimme(Y, ctx));
   CHKERR(Num_free_iprimme(pVLtVL, ctx));
   CHKERR(Num_free_SHprimme(fVLtVL, ctx));
   CHKERR(Num_free_SHprimme(C, ctx));

   /* Check orthogonality */
   /*if (numLocked) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*numLocked*numLocked);
      Num_gemm_Sprimme("C", "N", numLocked, numLocked, nLocal, 1.0, locked,
            ldLocked, locked, ldLocked, 0.0, H, numLocked);
      for(i=0; i < numLocked; i++) {
         for(j=0; j < i; j++) assert(fabs(*(double*)&H[numLocked*i+j]) < 1e-13);
         assert(fabs(1 - *(double*)&H[numLocked*i+i]) < 1e-13);
      }
      free(H);
   }
   if (b2) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*(b2)*(b2));
      Num_gemm_Sprimme("C", "N", b2, b2, nLocal, 1.0, V,
            ldV, V, ldV, 0.0, H, b2);
      for(i=0; i < b2; i++) {
         for(j=0; j < i; j++) assert(ABS(H[b2*i+j]) < 1e-13);
         assert(ABS(H[b2*i+i]) == 0.0 || ABS(1.0 - H[b2*i+i]) < 1e-13);
      }
      free(H);
   }
   if (numLocked) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*b2*numLocked);
      Num_gemm_Sprimme("C", "N", numLocked, b2, nLocal, 1.0, locked,
            ldLocked, basis, ldBasis, 0.0, H, numLocked);
      for(i=0; i < b2; i++) {
         for(j=0; j < numLocked; j++) assert(fabs(*(double*)&H[numLocked*i+j]) < 1e-13);
      }
      free(H);
   }*/

   return 0;
}

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
      HREAL *overlaps, HREAL *norms, primme_context ctx) {

   primme_params *primme = ctx.primme;
   int i, j, M=PRIMME_BLOCK_SIZE, m=min(M, mQ);
   HSCALAR *y, *y0;
   SCALAR *X0;

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

/******************************************************************************
 * Function Num_ortho_kernel - This subroutine performs the next operations:
 *
 *    X = (X - [Q V]*A)*Y/D if A, D and Y provided
 *    B = [Q W]'*Z
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * Q           input matrix
 * m           number of rows on Q, V, W, Z, X
 * nQ          number of columns of Q
 * ldQ         leading dimension of Q
 * V           input matrix of size m x nV
 * ldV         leading dimension of V
 * A           input matrix of size nQ+nV x nX
 * ldA         leading dimension of A
 * D           diagonal of a diagonal matrix (optional)
 * Y           input matrix of size nX x nX (optional)
 * ldY         leading dimension of Y (optional)
 * W           input matrix of size m x nW
 * ldW         leading dimension of W
 * Z           input matrix of size m x nZ
 * ldZ         leading dimension of Z
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * X           input/output matrix (optional)
 * nX          number of columns of X (optional)
 * ldX         leading dimension of X (optional)
 * B           output matrix of size nQ+nW x nZ
 * ldB         leading dimension of B
 * rwork       scalar workspace
 * ldrwork     size of rwork
 *
 ******************************************************************************/

static int Num_ortho_kernel(SCALAR *Q, PRIMME_INT M, int nQ, PRIMME_INT ldQ,
      SCALAR *V, int nV, PRIMME_INT ldV, SCALAR *X, int nX, PRIMME_INT ldX,
      HSCALAR *A, int ldA, HREAL *D, HSCALAR *Y, int ldY,
      SCALAR *W, int nW, PRIMME_INT ldW, SCALAR *Z, int nZ, PRIMME_INT ldZ,
      HSCALAR *B, int ldB, primme_context ctx)
{

   PRIMME_INT i;     /* Loop variables */
   // int m=min(PRIMME_BLOCK_SIZE, M);   /* Number of rows in the cache */
   int m;
   //int m=min(max((1024*1024 - 2*(nQ+nX)*(nQ+nX))/((nQ+nX)*3),1), M);   /* Number of rows in the cache */
   // int m=M;   /* Number of rows in the cache */
   SCALAR *Xo;
   HSCALAR *Bo;

   if (A && D && Y) {
      m = min(PRIMME_BLOCK_SIZE, M);   /* Number of rows in the cache */
   }
   else {
      m = M;
   }

   CHKERR(Num_malloc_SHprimme((nQ+nW)*nZ, &Bo, ctx));
   CHKERR(Num_malloc_Sprimme(m*nX, &Xo, ctx));

   /* Bo = zeros() */
   Num_zero_matrix_SHprimme(Bo, nQ + nW, nZ, nQ + nW, ctx);

   /* D[i] = 1/D[i] */
   for (i=0; i < nX; i++) D[i] = 1.0/D[i];

   for (i=0; i < M; i+=m, m=min(m,M-i)) {
      if (A && D && Y) {
         /* X = X - Q*A(0:nQ,:) */
         Num_gemm_dhd_Sprimme("N", "N", m, nX, nQ, -1.0, &Q[i], ldQ, A, ldA,
               1.0, &X[i], ldX, ctx);

         /* X = X - V*A(nQ:nQ+nV) */
         Num_gemm_dhd_Sprimme("N", "N", m, nX, nV, -1.0, &V[i], ldV, A + nQ,
               ldA, 1.0, &X[i], ldX, ctx);

         /* Xo = X*Y */
         Num_gemm_dhd_Sprimme(
               "N", "N", m, nX, nX, 1.0, &X[i], ldX, Y, ldY, 0.0, Xo, m, ctx);

         /* X = Xo*D */
         Num_scale_matrix_Sprimme(Xo, m, nX, m, D, &X[i], ldX, ctx);
      }

      /* Bo(0:nQ-1,:) += Q'*Z */
      Num_gemm_ddh_Sprimme("C", "N", nQ, nZ, m, 1.0, &Q[i], ldQ, &Z[i], ldZ,
            1.0, Bo, nQ + nW, ctx);

      /* Bo(nQ:nQ+nW-1) += W'*Z */
      Num_gemm_ddh_Sprimme("C", "N", nW, nZ, m, 1.0, &W[i], ldW, &Z[i], ldZ,
            1.0, Bo + nQ, nQ + nW, ctx);
   }

   /* B = globalSum(Bo) */
   CHKERR(globalSum_SHprimme(Bo, Bo, (nQ + nW) * nZ, ctx));
   Num_copy_matrix_SHprimme(Bo, nQ+nW, nZ, nQ+nW, B, ldB, ctx);

   CHKERR(Num_free_SHprimme(Bo, ctx));
   CHKERR(Num_free_Sprimme(Xo, ctx));

   return 0;
}

/*******************************************************************************
 * Subroutine eig - This procedure compute the eigen-decomposition of H
 *        
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * H             The input matrix of size n x n
 * n             Size of H
 * ldH           The leading dimension of H
 * primme        Structure containing various solver parameters
 * 
 * INPUT/OUTPUT ARRAYS
 * -------------------
 * Y             Right singular vectors of H
 * ldY           The leading dimension of Y
 * evals         The eigenvalues of H, non-increasing order
 * rwork         Workspace
 * lrwork        Length of the work array rwork
 *
 * Return Value
 * ------------
 * int -  0 upon successful return
 *     - -1 was unsuccessful
 ******************************************************************************/

static int eig(HSCALAR *H, int n, int ldH, HSCALAR *Y, int ldY, HREAL *evals,
      primme_context ctx) {

   int i, j; /* Loop variables    */

   /* XHEEV returns the eigenvalues in ascending order. This function will    */
   /* return them in descending order. An easy way is to compute the          */
   /* eigenpairs of -H instead. Do Y = -H.                                    */

   for (i=0; i < n; i++) {
      for (j = 0; j <= i; j++) { Y[ldY * i + j] = -H[ldH * i + j]; }
   }
 
   CHKERR(Num_heev_SHprimme("V", "U", n, Y, ldY, evals, ctx));

   /* evals = -evals */

   for (i=0; i < n; i++) { 
      evals[i] *= -1.0;
   }

   return 0;
}

static int rank_estimation(HSCALAR *V, int n0, int n1, int n, int ldV) {

   (void)n0;

   int i, j;

   for(i=0; i<n1; i++) {
      HREAL norm1 = 0.0;
      for (j = 0; j < i; j++)
         norm1 += ABS(V[i * ldV + j]) /
                  sqrt(ABS(V[i * ldV + i]) * ABS(V[j * ldV + j]));
      if (norm1 > (i+1)*1.0/n*.8) break;
      // TEMP!!!
      for (j = 0; j < i; j++)
      if (ABS(V[i * ldV + j]) /
                  sqrt(ABS(V[i * ldV + i]) * ABS(V[j * ldV + j])) > .8/n) break;
   }

   return i;
}
