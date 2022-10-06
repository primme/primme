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
   SCALAR *VtBV;              /* V'*B*V */

   HSCALAR *H;                /* Upper triangular portion of V'*A*V            */
   HSCALAR *hVecs;            /* Eigenvectors of H                             */
   HSCALAR *RhVecs;           /* Ritz vectors                                  */

   PRIMME_INT ldV;            /* The leading dimension of V                    */
   PRIMME_INT ldH;            /* The leading dimension of H                    */
   PRIMME_INT ldVtBV;         /* The leading dimension of VtBV                 */
   PRIMME_INT ldhVecs;        /* The leading dimension of hVecs                */
   PRIMME_INT ldRhVecs;       /* The leading dimension of RhVecs               */

   int i, j, k;               /* Loop variable                                 */
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
   int maxBlockSize = blockSize = min(primme->maxBlockSize, maxBasisSize);    /* In case user enters in too big a block size */

   if (primme->locking)
      maxRank = primme->numOrthoConst + primme->numEvals + primme->maxBasisSize;
   else
      maxRank = primme->numOrthoConst + primme->maxBasisSize;
  
   maxEvecsSize = primme->numOrthoConst + primme->numEvals;

   /* -------------------------------------------------------------- */
   /* Allocate objects                                               */
   /* -------------------------------------------------------------- */

   /* Setting up leading dimensions for matrices */
   ldV = ldRhVecs = primme->ldOPs;
   ldH = ldhVecs = ldVtBV = maxBasisSize;

   CHKERR(Num_malloc_Sprimme(ldV*(maxBasisSize+maxBlockSize), &V, ctx));
   CHKERR(Num_malloc_Sprimme(maxBasisSize*maxBasisSize, &VtBV, ctx));

   CHKERR(Num_malloc_SHprimme(maxBasisSize*maxBasisSize, &H, ctx));
   CHKERR(Num_malloc_SHprimme(maxBasisSize*maxBasisSize, &hVecs, ctx));
   CHKERR(Num_malloc_SHprimme(nLocal*maxBasisSize, &RhVecs, ctx));
   CHKERR(Num_malloc_SHprimme(maxBasisSize, &hSVals, ctx));

   CHKERR(KIND(Num_malloc_RHprimme, Num_malloc_SHprimme)(maxBasisSize, &hVals, ctx)); 
   CHKERR(Num_malloc_RHprimme(maxBlockSize, &blockNorms, ctx));
   CHKERR(Num_malloc_RHprimme(maxBasisSize, &basisNorms, ctx));
   CHKERR(Num_malloc_iprimme(maxBasisSize, &flags, ctx));

   CHKERR(Num_zero_matrix_SHprimme(H, maxBasisSize, maxBasisSize, ldH, ctx));
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
   CHKERR(Num_larnv_Sprimme(2, primme->iseed, nLocal*blockSize, &V[0], ctx));
   //for(i = 0; i < blockSize*nLocal; i++)
   //   V[i] = i+1;
   CHKERR(ortho_Sprimme(&V[0], ldV, NULL, 0, 0, blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));

   /* Initial iteration before for loop --------------------------------------------------- */
   /* Adjust block size first */   
   blockSize = min(blockSize, maxBasisSize - maxBlockSize);

   /* Put next chuck of vectors in the basis */
   CHKERR(matrixMatvec_Sprimme(&V[0], nLocal, ldV, &V[maxBlockSize*ldV], ldV, 0, blockSize, ctx));

   /* Compute and subtract first alpha */
   CHKERR(update_projection_Sprimme(&V[maxBlockSize*ldV], ldV, &V[0], ldV, &H[0], ldH, nLocal, 0, blockSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));
   CHKERR(Num_gemm_Sprimme("N", "N", nLocal, blockSize, blockSize, -1.0, &V[0], ldV, &H[0], ldH, 1.0, &V[maxBlockSize*ldV], ldV, ctx)); 

   for(i = maxBlockSize; i < maxBasisSize; i += blockSize) {

      /* Adjust block size in case maxBasisSize % maxBlockSize != 0 */
      blockSize = min(blockSize, maxBasisSize - i);

      /* Ortho to prepare for next iteration */
      CHKERR(ortho_Sprimme(&V[i*ldV], ldV, &H[i*ldH + (i-blockSize)], ldH, 0, blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));

      /* Put next chunk of vectors in the basis */
      CHKERR(matrixMatvec_Sprimme(&V[i*ldV], nLocal, ldV, &V[(i+blockSize)*ldV], ldV, 0, blockSize, ctx));

      /* Subtract beta term from new chunk of vectors */
      CHKERR(Num_gemm_Sprimme("N", "T", nLocal, blockSize, blockSize, -1.0, &V[(i-blockSize)*ldV], ldV, &H[i*ldH + (i-blockSize)], ldH, 1.0, &V[(i+blockSize)*ldV], ldV, ctx)); 

      /* Find and subtract alpha term */
      CHKERR(update_projection_Sprimme(&V[(i+blockSize)*ldV], ldV, &V[i*ldV], ldV, &H[i*ldH+i], ldH, nLocal, 0, blockSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));
      CHKERR(Num_gemm_Sprimme("N", "N", nLocal, blockSize, blockSize, -1.0, &V[i*ldV], ldV, &H[i*ldH+i], ldH, 1.0, &V[(i+blockSize)*ldV], ldV, ctx)); 
      
      CHKERR(update_projection_Sprimme(&V[0], ldV, &V[0], ldV, &VtBV[0], ldVtBV, nLocal, 0, maxBasisSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));                      /* VtBV = V'*B*V */
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

   /* VtBV = V'BV */
   CHKERR(update_projection_Sprimme(&V[0], ldV, &V[0], ldV, &VtBV[0], ldVtBV, nLocal, 0, maxBasisSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx));                      /* VtBV = V'*B*V */

   /* Solve the projected problem */
   CHKERR(solve_H_Sprimme(&H[0], maxBasisSize, ldH, VtBV, ldVtBV, NULL, 0, NULL, 0, NULL, 0, NULL, 0, &hVecs[0], ldhVecs, &hVals[0], &hSVals[0], numConverged, ctx));   

   /* Find Ritz Vectors */
   CHKERR(Num_gemm_Sprimme("N", "N", nLocal, maxBasisSize, maxBasisSize, 1.0, &V[0], ldV, &hVecs[0], ldhVecs, 0.0, &RhVecs[0], ldRhVecs, ctx));                 /* Va = V_{i-1}*b_{i}'           */

   /* Check the convergence of the Ritz vectors */
   //CHKERR(check_convergence_Sprimme(&RhVecs[0], ldRhVecs, 1, NULL, 0, 0, NULL, 0, 0, SCALAR *Bevecs, PRIMME_INT ldBevecs, &VtBV[0], ldVtBV, 0, maxBasisSize-1, &flags[0], &blockNorms[0], HEVAL *hVals, int *reset, int practConvCheck, primme_context ctx));

   /* ---------------------------------------------------------- */
   /* Deallocate arrays                                          */
   /* ---------------------------------------------------------- */
clean:
   CHKERR(Num_free_Sprimme(V, ctx));
   CHKERR(Num_free_Sprimme(VtBV, ctx));

   CHKERR(Num_free_SHprimme(H, ctx));
   CHKERR(Num_free_SHprimme(hVecs, ctx));
   CHKERR(Num_free_SHprimme(RhVecs, ctx));
   CHKERR(Num_free_SHprimme(hSVals, ctx));

   CHKERR(Num_free_RHprimme(blockNorms, ctx));
   CHKERR(Num_free_RHprimme(basisNorms, ctx));
   CHKERR(KIND(Num_free_RHprimme, Num_free_SHprimme)(hVals, ctx)); 

   CHKERR(Num_free_iprimme(flags, ctx));

   *ret = 0;      /* Set return value to PASSED */

return 0;
}

#endif   // SUPPORT_TYPE
