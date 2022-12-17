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
 * File: lanczos.c
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
#include "sketch.h"
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
      HREAL *resNorms, int *ret, int *numRet, int fullOrtho,
      primme_context ctx) {

   /* Default error is something is wrong in this function */
   *ret = PRIMME_MAIN_ITER_FAILURE;

   primme_params *primme = ctx.primme;

                            /* primme parameters */
   SCALAR *V;                 /* Basis vectors                                 */
   SCALAR *AVhVecs;           /* A*V*hVecs                                     */
   SCALAR *rwork;             /* temporary work vector                         */
   SCALAR *identity;          /* Used to copy beta from the lower triangular to the upper  */

   HSCALAR *H;                /* Upper triangular portion of V'*A*V            */
   HSCALAR *hVecs;            /* Eigenvectors of H                             */
   HSCALAR *hVals;

   PRIMME_INT ldV;            /* The leading dimension of V                    */
   PRIMME_INT ldAVhVecs;      /* The leading dimension of AVhVecs              */
   PRIMME_INT ldidentity;     /* The leading dimension of identity             */
   PRIMME_INT ldH;            /* The leading dimension of H                    */
   PRIMME_INT ldhVecs;        /* The leading dimension of hVecs                */
   PRIMME_INT ldrwork;        /* The leading dimension of rwork                */

   int i, j;                  /* Loop variable                                 */
   int blockSize;             /* Current block size                            */
   int numConverged;          /* Number of converged Ritz pairs                */
   int *flags;                /* Indicates which Ritz values have converged    */
   int reset = 0;             /* reset variable for check_convergence          */

   HREAL *basisNorms;         /* Residual norms of basis at pairs              */
   HREAL *blockNorms;         /* Residual norms corresponding to current block */
                              /* vectors.                                      */

   PRIMME_INT nLocal = primme->nLocal;
   PRIMME_INT maxBasisSize = primme->maxBasisSize;
   int maxBlockSize = blockSize = min(primme->maxBlockSize, maxBasisSize);    /* In case user enters in too big a block size */


   /* -------------------------------------------------------------- */
   /* Allocate objects                                               */
   /* -------------------------------------------------------------- */

   /* Setting up leading dimensions for matrices */
   ldV = ldAVhVecs = ldrwork = primme->ldOPs;
   ldhVecs = maxBasisSize;
   ldidentity = maxBlockSize;

   if(primme->projectionParams.projection == primme_proj_sketched)
   {
      ldH = maxBlockSize + maxBasisSize;

      CHKERR(Num_malloc_SHprimme((maxBasisSize+blockSize)*maxBasisSize, &H, ctx));
      CHKERR(Num_zero_matrix_SHprimme(H, maxBasisSize+maxBlockSize, maxBasisSize, ldH, ctx));
   }else{
      ldH = maxBasisSize;

      CHKERR(Num_malloc_SHprimme(maxBasisSize*maxBasisSize, &H, ctx));
      CHKERR(Num_zero_matrix_SHprimme(H, maxBasisSize, maxBasisSize, ldH, ctx));
   }

   CHKERR(Num_malloc_Sprimme(ldV*(maxBasisSize+maxBlockSize), &V, ctx));
   CHKERR(Num_malloc_Sprimme(ldAVhVecs*primme->numEvals, &AVhVecs, ctx));
   CHKERR(Num_malloc_Sprimme(ldrwork*maxBasisSize, &rwork, ctx));
   CHKERR(Num_malloc_Sprimme(maxBlockSize*maxBlockSize, &identity, ctx));

   CHKERR(Num_malloc_SHprimme(maxBasisSize*maxBasisSize, &hVecs, ctx));
   CHKERR(Num_malloc_SHprimme(maxBasisSize, &hVals, ctx));

   CHKERR(Num_malloc_RHprimme(maxBlockSize, &blockNorms, ctx));
   CHKERR(Num_malloc_RHprimme(maxBasisSize, &basisNorms, ctx));
   CHKERR(Num_malloc_iprimme(maxBasisSize, &flags, ctx));

   CHKERR(Num_zero_matrix_SHprimme(identity, maxBlockSize, maxBlockSize, ldidentity, ctx));

   for(i = 0; i < maxBasisSize; i++)
      flags[i] = UNCONVERGED;
   for(i = 0; i < maxBlockSize; i++)
      identity[i*ldidentity+i] = 1.0;

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

   /* Let's insert random vectors into the basis ------------------------------------------ */
   CHKERR(Num_larnv_Sprimme(3, primme->iseed, nLocal*blockSize, &V[0], ctx));
   CHKERR(ortho_Sprimme(V, nLocal, NULL, 0, 0, blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));

   /* Initial iteration before for loop --------------------------------------------------- */
   blockSize = min(blockSize, maxBasisSize - maxBlockSize); /* Adjust block size first */

   CHKERR(matrixMatvec_Sprimme(V, nLocal, nLocal, &V[maxBlockSize*ldV], nLocal, 0, blockSize, ctx)); /* W = A*V_0 */

   /* Compute and subtract first alpha (W = W - V*(W'V))*/
   CHKERR(update_projection_Sprimme(&V[blockSize*ldV], nLocal, &V[0], nLocal, &H[0], ldH, nLocal, 0, blockSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx)); 
   CHKERR(Num_gemm_Sprimme("N", "N", nLocal, blockSize, blockSize, 1.0, &V[0], nLocal, &H[0], ldH, 0.0, &rwork[0], ldrwork, ctx)); 
   CHKERR(Num_axpy_Sprimme(ldV*blockSize, -1.0, &rwork[0], 1, &V[ldV*blockSize], 1, ctx)); 

   for(i = maxBlockSize; i < maxBasisSize; i += blockSize) {

      blockSize = min(blockSize, maxBasisSize - i); /* Adjust block size */

      CHKERR(ortho_Sprimme(&V[i*ldV], nLocal, &H[(i-blockSize)*ldH + i], ldH, 0, blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));   /* [V_i, b_i] = qr(V_i) */
      if(fullOrtho)
         CHKERR(ortho_Sprimme(V, nLocal, NULL, 0, i, i+blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));   /* V_i = cgs(V(:, 0:i-1), V_i) */

      CHKERR(Num_gemm_Sprimme("N", "C", blockSize, blockSize, blockSize, 1.0, &identity[0], ldidentity, &H[(i-blockSize)*ldH + i], ldH, 0.0, &H[i*ldH + (i-blockSize)], ldH, ctx));  

      CHKERR(matrixMatvec_Sprimme(&V[i*nLocal], nLocal, nLocal, &V[(i+blockSize)*nLocal], nLocal, 0, blockSize, ctx));                             /* V_{i+1} = AV_i */

      /* Subtract beta term from new chunk of vectors (W = W - V_{i-1}*b_i)*/
      CHKERR(Num_gemm_Sprimme("N", "C", nLocal, blockSize, blockSize, 1.0, &V[(i-blockSize)*nLocal], nLocal, &H[(i-blockSize)*ldH + i], ldH, 0.0, &rwork[0], ldrwork, ctx));  
      CHKERR(Num_axpy_Sprimme(nLocal*blockSize, -1.0, &rwork[0], 1, &V[nLocal*(i+blockSize)], 1, ctx)); 

      /* Find and subtract alpha term */
      CHKERR(update_projection_Sprimme(&V[(i+blockSize)*nLocal], nLocal, &V[i*nLocal], nLocal, &H[i*ldH+i], ldH, nLocal, 0, blockSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx)); /* a_i = W'*V_i */
      CHKERR(Num_gemm_Sprimme("N", "N", nLocal, blockSize, blockSize, 1.0, &V[i*nLocal], nLocal, &H[i*ldH+i], ldH, 0.0, &rwork[0], ldrwork, ctx)); /* rwork = V_i*a_i */
      CHKERR(Num_axpy_Sprimme(ldV*blockSize, -1.0, &rwork[0], 1, &V[nLocal*(i+blockSize)], 1, ctx)); /* W = W - rwork */

   }

   /* Moving on to the eigenvalue problem */
   primme->initSize = 0;

   /* Quick return for matrix of dimension 1 */
   if (primme->numEvals == 0) goto clean;

   if(primme->projectionParams.projection == primme_proj_sketched)
   {

      /* Adding a row to H */
      CHKERR(ortho_Sprimme(&V[ldV*maxBasisSize], nLocal, &H[(maxBasisSize-blockSize)*ldH + maxBasisSize], ldH, 0, blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));   /* [V_i, b_i] = qr(V_i) */
      if(fullOrtho)
         CHKERR(ortho_Sprimme(V, nLocal, NULL, 0, maxBasisSize, maxBasisSize+blockSize-1, NULL, 0, 0, nLocal, primme->iseed, ctx));   /* Orthogonalized the last block of V against the rest of the basis */

      /* Getting our sketched basis and projected sketched basis */
      CHKERR(apply_sketching_Sprimme(H, ldH, V, ldV, hVecs, ldhVecs, hVals, maxBasisSize, blockSize, ctx));

   } else {
      CHKERR(solve_H_Sprimme(&H[0], maxBasisSize, ldH, NULL, 0, NULL, 0, NULL, 0, NULL, 0, NULL, 0, hVecs, ldhVecs, hVals, NULL, 0, ctx));
   }

   for(i = 0; i < primme->numEvals; i++) evals[i] = REAL_PART(hVals[i]);

   /* Find Ritz Vectors */
   CHKERR(Num_gemm_Sprimme("N", "N", nLocal, primme->numEvals, maxBasisSize, 1.0, V, nLocal, hVecs, ldhVecs, 0.0, evecs, ldevecs, ctx));    /* evecs = V*hVecs */

   /* Find residual norms */
   CHKERR(matrixMatvec_Sprimme(evecs, nLocal, ldevecs, AVhVecs, ldAVhVecs, 0, primme->numEvals, ctx)); /* AVhVecs = A*V*hVecs */
   CHKERR(Num_compute_residuals_Sprimme(nLocal, primme->numEvals, evals, evecs, ldevecs, AVhVecs, ldAVhVecs, rwork, ldrwork, ctx));
   for(j = 0; j < primme->numEvals; j++) resNorms[j] = sqrt(Num_dot_Sprimme(ldrwork, &rwork[j*ldrwork], 1, &rwork[j*ldrwork], 1, ctx));

   /* Check the convergence of the Ritz vectors */
   CHKERR(check_convergence_Sprimme(evecs, ldevecs, 1 /* given X */, NULL, 0, 0 /* not given R */, NULL, 0, 0, NULL, 0, NULL, 0, 0, primme->numEvals, flags, resNorms, hVals, &reset, -1, ctx));

   /* Find number of converged eigenpairs */
   numConverged = 0;
   for(i = 0; i < primme->numEvals && flags[i] == CONVERGED; i++) numConverged = i+1;

   primme->initSize = numConverged;
   *numRet = numConverged;

   /* ---------------------------------------------------------- */
   /* Deallocate arrays                                          */
   /* ---------------------------------------------------------- */
   clean:

   CHKERR(Num_free_Sprimme(V, ctx));
   CHKERR(Num_free_Sprimme(AVhVecs, ctx));
   CHKERR(Num_free_Sprimme(rwork, ctx));
   CHKERR(Num_free_Sprimme(identity, ctx));

   CHKERR(Num_free_SHprimme(H, ctx));
   CHKERR(Num_free_SHprimme(hVecs, ctx));
   CHKERR(Num_free_SHprimme(hVals, ctx));

   CHKERR(Num_free_RHprimme(blockNorms, ctx));
   CHKERR(Num_free_RHprimme(basisNorms, ctx));

   CHKERR(Num_free_iprimme(flags, ctx));

   *ret = 0;

return 0;
}

#endif   // SUPPORT_TYPE
