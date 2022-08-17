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
   int i;                   /* Loop variable                                 */
   int basisSize;           /* Current size of the basis V                   */

   SCALAR *V;               /* Basis vectors                                 */
   PRIMME_INT ldV;          /* The leading dimension of V                    */
   SCALAR *W;               /* Work space storing A*V                        */
   PRIMME_INT ldW;          /* The leading dimension of W                    */
   SCALAR *Va;              /* Temporary work array of size nLocalxblockSize */
   PRIMME_INT ldVa;         /* The leading dimension of Va                   */
   HSCALAR *T;              /* Holds the tridiagonal matrix from Lanczos     XXX: I don't define "T" right now, I have the alphas stored in SCALAR *a, and the betas in SCALAR *B*/
   PRIMME_INT ldT;          /* The leading dimension of T                    */
   SCALAR *a;               /* Inner product of W'V                          */
   PRIMME_INT lda;          /* The leading dimension of a                    */
   SCALAR *B;               /* R factor of W                                 */
   PRIMME_INT ldB;          /* The leading dimension of B                    */


   PRIMME_INT nLocal = primme->nLocal;
   PRIMME_INT maxBasisSize = primme->maxBasisSize;
   int blockSize = primme->maxBlockSize;

   /* -------------------------------------------------------------- */
   /* Allocate objects                                               */
   /* -------------------------------------------------------------- */

   /* Use leading dimension ldOPs for the large dimension mats: V, W, T, a, and B */

   ldV = ldVa = ldW = ldT = primme->ldOPs;
   lda = ldB = blockSize;
   CHKERR(Num_malloc_Sprimme(ldV*maxBasisSize, &V, ctx));
   CHKERR(Num_malloc_Sprimme(ldVa*blockSize, &Va, ctx));
   CHKERR(Num_malloc_Sprimme(ldW*blockSize, &W, ctx));
   CHKERR(Num_malloc_SHprimme(ldT*ldT, &T, ctx));
   CHKERR(Num_malloc_Sprimme(maxBasisSize*blockSize, &a, ctx));
   CHKERR(Num_malloc_Sprimme(maxBasisSize*blockSize, &B, ctx));

   /* -------------------------------------------------------------- */
   /* Lanczos algorithm                                              */
   /* -------------------------------------------------------------- */

   /* Write first block of orthonormal random vectors ---------- */
   for( i=0; i < blockSize; i++) { 
      CHKERR(Num_larnv_Sprimme(2, primme->iseed, nLocal, &V[i], ctx));
   }
   CHKERR(ortho_Sprimme(V, ldV, NULL, 0, 0, blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));

   /* Initial iteration step ---------------------------------- */
   CHKERR(matrixMatvec_Sprimme(&V[0], nLocal, ldV, &W[0], ldW, maxBasisSize, blockSize, ctx));                             // W = A*V_1 
   CHKERR(Num_gemm_Sprimme("T", "N", blockSize, blockSize, 0, 1.0, &W[0], ldW, &V[0], ldV, 1.0, &a[0], lda, ctx));         // a = W' * V_1;
   CHKERR(Num_gemm_Sprimme("N", "N", nLocal, blockSize, 0, 1.0, &V[0], ldV, &a[0], lda, 1.0, &Va[0], ldVa, ctx));          // Va = V_1*a
   CHKERR(Num_axpy_Sprimme(nLocal, -1.0, &Va[0], 1, &W[0], 1, ctx));                                                       // W = W - Va
   CHKERR(ortho_Sprimme(&W[0], ldW, &B[0], ldB, 0, blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));                  // [W, B] = QR(W)
   CHKERR(Num_copy_matrix_Sprimme(&W[0], nLocal, blockSize, ldW, &V[0], ldV, ctx));                                        // V_2 = W
 
   /* Main loop of Lanczos ------------------------------------ */
   // XXX: This assumes the block size divides evenly the size of the basis
   for( i=2*blockSize; i < maxBasisSize; i+=blockSize) {
      CHKERR(matrixMatvec_Sprimme(&V[ldV*(i-blockSize)], nLocal, ldV, &W[0], ldW, maxBasisSize, blockSize, ctx));                                           // W = A*V_{i-1} 
      CHKERR(Num_gemm_Sprimme("N", "T", nLocal, blockSize, 0, 1.0, &V[ldV*(i-2*blockSize)], ldV, &B[ldB*(i-2*blockSize)], ldB, 1.0, &Va[0], ldVa, ctx));    // Va = V_{i-2}*B_{i-1}'
      CHKERR(Num_axpy_Sprimme(nLocal, -1.0, &Va[0], 1, &W[0], 1, ctx));                                                                                     // W = W - Va
      CHKERR(Num_gemm_Sprimme("T", "N", blockSize, blockSize, 0, 1.0, &W[0], ldW, &V[ldV*(i-blockSize)], ldV, 1.0, &a[lda*(i-blockSize)], lda, ctx));       // a = W' * V_i-1;
      CHKERR(Num_gemm_Sprimme("N", "N", nLocal, blockSize, 0, 1.0, &V[ldV*(i-blockSize)], ldV, &a[lda*(i-blockSize)], lda, 1.0, &Va[0], ldVa, ctx));        // Va = V_{i-1}*a_{i}'
      CHKERR(Num_axpy_Sprimme(nLocal, -1.0, &Va[0], 1, &W[0], 1, ctx));                                                                                     // W = W - Va
      CHKERR(ortho_Sprimme(&W[0], ldW, &B[ldB*(i-blockSize)], ldB, 0, blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));                                // [W, B_i] = QR(W)
      CHKERR(Num_copy_matrix_Sprimme(&W[0], nLocal, blockSize, ldW, &V[i], ldV, ctx));                                                                      // V_j+1 = W
   }

   /* ---------------------------------------------------------- */
   /* Deallocate arrays                                          */
   /* ---------------------------------------------------------- */
   CHKERR(Num_free_Sprimme(V, ctx));
   CHKERR(Num_free_Sprimme(Va, ctx));
   CHKERR(Num_free_Sprimme(W, ctx));
   CHKERR(Num_free_SHprimme(T, ctx));
   CHKERR(Num_free_Sprimme(a, ctx));
   CHKERR(Num_free_Sprimme(B, ctx));


 }
