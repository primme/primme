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
 * File: main_iter.c
 *
 * Purpose - This is the main Lanczos-type iteration 
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/lanczos.c"
#endif


#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "lanczos.h"
#include "convergence.h"
#include "correction.h"
#include "init.h"
#include "ortho.h"
#include "restart.h"
#include "solve_projection.h"
#include "update_projection.h"
#include "update_W.h"
#include "auxiliary_eigs.h"
#include "auxiliary_eigs_normal.h"
#endif

#ifdef SUPPORTED_TYPE

/******************************************************************************
 * Subroutine lanczos - This routine implements a more general, parallel, 
 *    block Lanczos outer iteration
 *
 * A coarse outline of the algorithm performed is as follows:
 *
 *  1. Initialize basis V_1
 *     Initial iteration step:
 *  2. W_1' = A * V_1
 *  3. a_1 = W_1' * V_1
 *  4. W_1 = W_1' - a_1 * V_1
 *  2. for j = 2, 3, ...., basisSize do
 *  3.    b_j = ||W_{j-1}||
 *  4.    V_j = W_{j-1} / b_j
 *  5.    W_j' = A * V_j
 *  6.    a_1 = W_j' * V_j
 *  7.    W_j = W_j' - a_j * V_j - b_j * V_{j-1}
 *  8. endfor
 *  9. Let T = [a_1   b_2                              0  ]
 *             [b_2   a_2   b_3                           ]
 *             [      b_3   a_3   b_4                     ]
 *             [            ...   ...   ...               ]
 *             [                  b_{m-1}   a_{m-1}   b_m ]
 *             [  0                         b_m       a_m}]
 *
 * 
 * INPUT/OUTPUT arrays and parameters
 * ----------------------------------
 * evecs    Stores initial guesses. Upon return, it contains the converged Ritz
 *          vectors.  If locking is not engaged, then converged Ritz vectors 
 *          are copied to this array just before return.  
 *
 * primme.initSize: On output, it stores the number of converged eigenvectors. 
 *           If smaller than numEvals and locking is used, there are
 *              only primme.initSize vectors in evecs.
 *           Without locking all numEvals approximations are in evecs
 *              but only the initSize ones are converged.
 *           During the execution, access to primme.initSize gives 
 *              the number of converged pairs up to that point. The pairs
 *              are available in evals and evecs, but only when locking is used
 * ret      successful state
 * numRet   number of returned pairs in evals, evecs, and resNorms
 *
 * Return Value
 * ------------
 * error code      
 ******************************************************************************/

TEMPLATE_PLEASE
int lanczos_Sprimme(HEVAL *evals, SCALAR *evecs, PRIMME_INT ldevecs,
      HREAL *resNorms, double startTime, int *ret, int *numRet,
      primme_context ctx) {

   /* Default error is something is wrong in this function */
   *ret = PRIMME_MAIN_ITER_FAILURE;

   primme_params *primme = ctx.primme;

                            /* primme parameters */
   SCALAR *V;                 /* Basis vectors                                 */
   SCALAR *AV;                /* Projected Basis vectors (AV = A*V)            */
   SCALAR *W;                 /* Work space storing A*V                        */
   SCALAR *Va;                /* Temporary work array of size nLocalxblockSize */
   SCALAR *VtBV;              /* V'*B*V */
   SCALAR *scales;            /* Scale basis vector blocks                     */

   SCALAR *H;                 /* Upper triangular portion of V'*A*V            */
   HSCALAR *hVecs;            /* Eigenvectors of H                             */

   PRIMME_INT ldV;            /* The leading dimension of V                    */
   PRIMME_INT ldAV;           /* The leading dimension of AV                   */
   PRIMME_INT ldW;            /* The leading dimension of W                    */
   PRIMME_INT ldVa;           /* The leading dimension of Va                   */
   PRIMME_INT ldVtBV;         /* The leading dimension of VtBV                 */
   PRIMME_INT ldH;            /* The leading dimension of H                    */
   PRIMME_INT ldhVecs;        /* The leading dimension of hVecs                */

   int i, j, k;                  /* Loop variable                                 */
   int blockSize;             /* Current block size                            */
   int maxEvecsSize;          /* Maximum capacity of evecs array               */
   int numConverged;          /* Number of converged Ritz pairs                */
   int *flags;                /* Indicates which Ritz values have converged    */
   int maxRank;               /* maximum size of the main space being orthonormalize */

   double smallestResNorm;    /* the smallest residual norm in the block       */
   HEVAL *hVals;              /* Eigenvalues of H                              */

   HREAL *hSVals;             /* Singular values of R                          */
   HREAL *basisNorms;         /* Residual norms of basis at pairs              */
   HREAL *blockNorms;         /* Residual norms corresponding to current block */
                              /* vectors.                                      */

   numConverged = 0;
   PRIMME_INT nLocal = primme->nLocal;
   PRIMME_INT maxBasisSize = primme->maxBasisSize;
   int maxBlockSize = min(primme->maxBlockSize, maxBasisSize);    /* In case user enters in too big a block size */
   blockSize = maxBlockSize;
   if (primme->locking)
      maxRank = primme->numOrthoConst + primme->numEvals + primme->maxBasisSize;
   else
      maxRank = primme->numOrthoConst + primme->maxBasisSize;
  
   maxEvecsSize = primme->numOrthoConst + primme->numEvals;

   /* -------------------------------------------------------------- */
   /* Allocate objects                                               */
   /* -------------------------------------------------------------- */

   /* Setting up leading dimensions for matrices */
   ldV = ldW = ldVa = ldAV = primme->ldOPs;
   ldH = ldhVecs = ldVtBV = maxBasisSize;

   CHKERR(Num_malloc_Sprimme(ldV*maxBasisSize, &V, ctx));
   CHKERR(Num_malloc_Sprimme(ldAV*maxBasisSize, &AV, ctx));
   CHKERR(Num_malloc_Sprimme(ldW*maxBlockSize, &W, ctx));
   CHKERR(Num_malloc_Sprimme(ldVa*blockSize, &Va, ctx));
   CHKERR(Num_malloc_Sprimme(maxBasisSize*maxBasisSize, &VtBV, ctx));
   CHKERR(Num_malloc_Sprimme(maxBlockSize, &scales, ctx));

   CHKERR(Num_malloc_Sprimme(maxBasisSize*maxBasisSize, &H, ctx));
   CHKERR(Num_malloc_SHprimme(maxBasisSize*maxBasisSize, &hVecs, ctx));
   CHKERR(Num_malloc_SHprimme(maxBasisSize, &hSVals, ctx));

   CHKERR(KIND(Num_malloc_RHprimme, Num_malloc_SHprimme)(maxBasisSize, &hVals, ctx)); 
   CHKERR(Num_malloc_RHprimme(maxBlockSize, &blockNorms, ctx));
   CHKERR(Num_malloc_RHprimme(maxBasisSize, &basisNorms, ctx));
   CHKERR(Num_malloc_iprimme(maxBasisSize, &flags, ctx));

   CHKERR(Num_zero_matrix_Sprimme(H, maxBasisSize, maxBasisSize, ldH, ctx));
   CHKERR(Num_zero_matrix_SHprimme(hVecs, maxBasisSize, maxBasisSize, ldhVecs, ctx));
   CHKERR(Num_zero_matrix_SHprimme(hSVals, 1, maxBasisSize, 1, ctx));
   CHKERR(Num_zero_matrix_RHprimme(basisNorms, 1, maxBasisSize, 1, ctx));

   for(i = 0; i < maxBasisSize; i++)
      flags[i] = UNCONVERGED;

  /* Initialize counters and flags ---------------------------- */

   primme->stats.numOuterIterations            = 0; 
   primme->stats.numRestarts                   = 0;
   primme->stats.numMatvecs                    = 0;
   primme->stats.numPreconds                   = 0;
   primme->stats.numGlobalSum                  = 0;
   primme->stats.numBroadcast                  = 0;
   primme->stats.volumeGlobalSum               = 0;
   primme->stats.volumeBroadcast               = 0;
   primme->stats.flopsDense                    = 0;
   primme->stats.numOrthoInnerProds            = 0.0;
   primme->stats.elapsedTime                   = 0.0;
   primme->stats.timeMatvec                    = 0.0;
   primme->stats.timePrecond                   = 0.0;
   primme->stats.timeOrtho                     = 0.0;
   primme->stats.timeGlobalSum                 = 0.0;
   primme->stats.timeBroadcast                 = 0.0;
   primme->stats.timeDense                     = 0.0;
   primme->stats.estimateMinEVal               = HUGE_VAL;
   primme->stats.estimateMaxEVal               = -HUGE_VAL;
   primme->stats.estimateLargestSVal           = -HUGE_VAL;
   primme->stats.estimateBNorm                 = primme->massMatrixMatvec ? -HUGE_VAL : 1.0;
   primme->stats.estimateInvBNorm              = primme->massMatrixMatvec ? -HUGE_VAL : 1.0;
   primme->stats.maxConvTol                    = 0.0;
   primme->stats.estimateResidualError         = 0.0;
   primme->stats.lockingIssue                  = 0;

   /* -------------------------------------------------------------- */
   /* Lanczos algorithm - Build basis (no restarting)                */
   /* -------------------------------------------------------------- */

   /* Setup first basis vector for Lanczos ------------------------------------------------ */
   //for( i=0; i < blockSize; i++)
   CHKERR(Num_larnv_Sprimme(2, primme->iseed, nLocal, &V[0], ctx));
   CHKERR(ortho_Sprimme(&V[0], ldV, NULL, 0, 0, blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));

   printf("First column after ortho: ");
   for(i = 0; i < nLocal; i++)
      printf("%f, ", V[i]);
   printf("\n");  

   /* Initial Iteration Step -------------------------------------------------------------- */
   CHKERR(matrixMatvec_Sprimme(&V[0], nLocal, ldV, &W[0], ldW, 0, blockSize, ctx));                                                                         /* W = A*V_{0}*/
   CHKERR(update_projection_Sprimme(&W[0], ldW, &V[0], ldV, &H[0], blockSize, nLocal, 0, blockSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));      /* a_0 = W'*V_{0}*/ 
   CHKERR(Num_gemm_Sprimme("N", "T", nLocal, blockSize, 0, 1.0, &V[0], ldV, &H[0], blockSize, 1.0, &Va[0], ldVa, ctx));                                     /* Va = V_{0}*a_{0}'           */
   CHKERR(Num_axpy_Sprimme(nLocal, -1.0, &Va[0], 1, &W[0], 1, ctx));                                                                                        /* W = W - (V_{0}*a_0)         */

   blockSize = min(maxBlockSize, maxBasisSize-blockSize);      /* Adjust block size if needed   */

   CHKERR(update_projection_Sprimme(&V[0], ldV, &V[0], ldV, &VtBV[0], 1, nLocal, 0, 1, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));  /* VtBV = V'*B*V */
   printf("Iteration 0 --------\n");
   printf("%f \n", VtBV[0]);

   /* Main loop of Lanczos ----------------------------------------- */
   for(i = maxBlockSize; i < maxBasisSize; i += blockSize) {
      blockSize = min(blockSize, maxBasisSize-i);                                                                                                                                /* Adjust block size if needed   */
 
      CHKERR(ortho_Sprimme(&W[0], ldW, &H[i*ldH + (i-blockSize)], blockSize, 0, blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));                                           /* [W, b_i] = QR(W)              */
      CHKERR(Num_copy_matrix_Sprimme(&W[0], nLocal, blockSize, ldW, &V[i*ldV], ldV, ctx));                                                                                       /* V_{j} = W                     */

      CHKERR(matrixMatvec_Sprimme(&V[(i-blockSize)*ldV], nLocal, ldV, &W[0], ldW, 0, blockSize, ctx));                                                                           /* W = A*V_{i-1}                 */
      CHKERR(update_projection_Sprimme(&W[0], ldW, &V[i*ldV], ldV, &H[i*ldH+i], blockSize, nLocal, 0, blockSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));              /* a_i = W'*V_{i}*/ 

      CHKERR(Num_gemm_Sprimme("N", "T", nLocal, blockSize, 0, 1.0, &V[ldV*i], ldV, &H[i*blockSize + i], blockSize, 1.0, &Va[0], ldVa, ctx));                                     /* Va = V_{i}*a_{i}'             */
      CHKERR(Num_axpy_Sprimme(nLocal, -1.0, &Va[0], 1, &W[0], 1, ctx));                                                                                                          /* W = W - Va                    */

      CHKERR(Num_gemm_Sprimme("N", "T", nLocal, blockSize, 0, 1.0, &V[ldV*(i-blockSize)], ldV, &H[i*ldH + (i-blockSize)], blockSize, 1.0, &Va[0], ldVa, ctx));                   /* Va = V_{i-1}*b_{i}'           */
      CHKERR(Num_axpy_Sprimme(nLocal, -1.0, &Va[0], 1, &W[0], 1, ctx));                                                                                                          /* W = W - Va                    */

      CHKERR(update_projection_Sprimme(&V[0], ldV, &V[0], ldV, &VtBV[0], i+blockSize, nLocal, 0, i+blockSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));  /* VtBV = V'*B*V */
      printf("Iteration %d (After) --------\n", i);
      for(j = 0; j < i+blockSize; j++){
         for(k = 0; k < i+blockSize; k++){
            printf("%f, ", VtBV[k*(i+blockSize)+j]);
         }
         printf("\n");
      }
      printf("\n\n");

   }

   /* ---------------------------------------------------------- */
   /* Moving onto eigenvalue problem                             */
   /* ---------------------------------------------------------- */

   /* Quick return for matrix of dimension 1 */
   if (primme->numEvals == 0) {
      primme->initSize = 0;
      *ret = 0;
      goto clean;
   }

   /* Now initSize will store the number of converged pairs */
   primme->initSize = 0;

   /* H = V'AV */
//   CHKERR(matrixMatvec_Sprimme(&V[0], nLocal, ldV, &AV[0], ldAV, 0, maxBasisSize, ctx));                                                            /* AV = A*V */
//   CHKERR(update_projection_Sprimme(&V[0], ldV, &AV[0], ldAV, &H[0], maxBasisSize, nLocal, 0, maxBasisSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));  /* H = V'*A*V */

   /* VtBV = V'BV */
   CHKERR(update_projection_Sprimme(&V[0], ldV, &V[0], ldV, &VtBV[0], maxBasisSize, nLocal, 0, maxBasisSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));  /* VtBV = V'*B*V */

   printf("THIS IS H ------------ \n");
   for(i = 0; i < maxBasisSize; i++){
      for(j = 0; j < maxBasisSize; j++){
         printf("%f, ", H[j*ldH+i]);
      }
      printf("\n");
   }
      printf("\n\n\n");
   printf("THIS IS VtBV ------------ \n");
   for(i = 0; i < maxBasisSize; i++){
      for(j = 0; j < maxBasisSize; j++){
         printf("%f, ", VtBV[j*maxBasisSize+i]);
      }
      printf("\n");
   }
      printf("\n\n\n");
   /* Solve the projected problem */
   CHKERR(solve_H_Sprimme(H, maxBasisSize, ldH, VtBV, ldVtBV, NULL, 0, NULL, 0, NULL, 0, NULL, 0, hVecs, ldhVecs, hVals, hSVals, numConverged, ctx));   

   /* Check the convergence of the Ritz vectors */
   //CHKERR(check_convergence_Sprimme(V, ldV, 1 /* given X */, NULL, 0, 0, evecs, 0, ldevecs, Bevecs, ldBevecs, VtBV, ldVtBV, 0, maxBasisSize, flags, resNorms, hVals, NULL, 0, ctx));

   /* ---------------------------------------------------------- */
   /* Deallocate arrays                                          */
   /* ---------------------------------------------------------- */
clean:
   CHKERR(Num_free_Sprimme(V, ctx));
   CHKERR(Num_free_Sprimme(Va, ctx));
   CHKERR(Num_free_Sprimme(W, ctx));
   CHKERR(Num_free_Sprimme(VtBV, ctx));
   CHKERR(Num_free_Sprimme(H, ctx));
   CHKERR(Num_free_SHprimme(hVecs, ctx));
   CHKERR(Num_free_RHprimme(blockNorms, ctx));
   CHKERR(Num_free_RHprimme(basisNorms, ctx));

   *ret = 0;      /* Set return value to PASSED */

return 0;
}

#endif   // SUPPORT_TYPE
