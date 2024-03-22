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

TEMPLATE_PLEASE
int print_lanczos_timings_Sprimme(PRIMME_INT basisSize, primme_context ctx) {

   primme_params *primme = ctx.primme;

   printf("Basis Size: %-" PRIMME_INT_P "\n", basisSize); 
   printf("Iterations: %-" PRIMME_INT_P "\n", primme->stats.numOuterIterations); 
   printf("Restarts  : %-" PRIMME_INT_P "\n", primme->stats.numRestarts);
   printf("Matvecs   : %-" PRIMME_INT_P "\n", primme->stats.numMatvecs);
   printf("Preconds  : %-" PRIMME_INT_P "\n", primme->stats.numPreconds);
   printf("Elapsed Time        : %-22.10E\n", primme->stats.elapsedTime);
   printf("MatVec Time         : %-22.10E\n", primme->stats.timeMatvec);
   printf("Precond Time        : %-22.10E\n", primme->stats.timePrecond);
   printf("Ortho Time          : %-22.10E\n", primme->stats.timeOrtho);
   printf("GlobalSum Time      : %-22.10E\n", primme->stats.timeGlobalSum);
   printf("Broadcast Time      : %-22.10E\n", primme->stats.timeBroadcast);
   printf("SketchedMatvec Time : %-22.10E\n", primme->stats.timeSketchMatvec);
   printf("Rayleigh-Ritz Time  : %-22.10E\n", primme->stats.timeRR);
   printf("Stabilization Time  : %-22.10E\n", primme->stats.timeStabilization);
   printf("Residual Time       : %-22.10E\n", primme->stats.timeResiduals);

   return 0;
}


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

#ifdef USE_HERMITIAN

   /* Default error is something is wrong in this function */
   *ret = PRIMME_MAIN_ITER_FAILURE;

   primme_params *primme = ctx.primme;

   double elapsed_time = primme_wTimer();
   double solve_timer;
   double resid_timer;
   primme->stats.elapsedTime = 0.0;  // XXX: Added this for debugging purposes - Heather

                            /* primme parameters */
   SCALAR *V;                 /* Basis vectors                                 */
   SCALAR *AVhVecs;           /* A*V*hVecs                                     */
   SCALAR *rwork;             /* temporary work vector                         */

   HSCALAR *H;                /* Upper triangular portion of V'*A*V            */
   HSCALAR *hVecs;            /* Eigenvectors of H                             */
   HREAL *hVals;              /* The Ritz values */

   PRIMME_INT ldV;            /* The leading dimension of V                    */
   PRIMME_INT ldAVhVecs;      /* The leading dimension of AVhVecs              */
   PRIMME_INT ldH;            /* The leading dimension of H                    */
   PRIMME_INT ldhVecs;        /* The leading dimension of hVecs                */
   PRIMME_INT ldrwork;        /* The leading dimension of rwork                */

   PRIMME_INT i, j;           /* Loop variables                                 */
   int blockSize;             /* Current block size                            */
   int numConverged;          /* Number of converged Ritz pairs                */
   int *flags;                /* Indicates which Ritz values have converged    */
   int *eval_perm;            /* For sorting the returned Ritz pairs           */
   int *global_start;         /* For building the initial starting vectors     */
   int reset = 0;             /* reset variable for check_convergence          */

   int maxBlockSize = blockSize = min(primme->maxBlockSize, primme->maxBasisSize);    /* In case user enters in too big a block size */
   primme->maxBasisSize -= primme->maxBasisSize % maxBlockSize; /* Ensure the maximum basis size is a multiple of the maximum blocksize */

   /* Sketching Variables -----------------------------------------------------*/
   SCALAR *V_temp;            /* Copy of basis vectors for sketching           */
   SCALAR *SV;                /* The sketched basis                            */
   SCALAR *T;                 /* The "R" factor in the QR decomposition of SV  */
   SCALAR *SW;                /* The projected sketched basis                  */
   SCALAR *S_vals;            /* CSC Formatted Values                          */
   PRIMME_INT *S_rows;        /* CSC Formatted Rows                            */
   PRIMME_INT nnzPerCol, ldSV, ldSW, ldT;  /* Size and nnz of the sketching matrix */
   REAL *normalize_evecs;

   /* -------------------------------------------------------------- */
   /* Allocate objects                                               */
   /* -------------------------------------------------------------- */

   /* Setting up leading dimensions for matrices */
   ldV = ldAVhVecs = ldrwork = primme->ldOPs;
   ldhVecs = primme->maxBasisSize;

   if(primme->projectionParams.projection == primme_proj_sketched)
   {
      ldH = maxBlockSize + primme->maxBasisSize;

      CHKERR(Num_malloc_SHprimme(ldH*primme->maxBasisSize, &H, ctx));
      CHKERR(Num_zero_matrix_SHprimme(H, ldH, primme->maxBasisSize, ldH, ctx));
   }else{
      ldH = primme->maxBasisSize;

      CHKERR(Num_malloc_SHprimme(ldH*ldH, &H, ctx));
      CHKERR(Num_zero_matrix_SHprimme(H, ldH, ldH, ldH, ctx));
   }

   CHKERR(Num_malloc_Sprimme(ldV*(primme->maxBasisSize+maxBlockSize), &V, ctx));
   CHKERR(Num_malloc_Sprimme(ldAVhVecs*primme->numEvals, &AVhVecs, ctx));
   CHKERR(Num_malloc_Sprimme(ldrwork*primme->maxBasisSize, &rwork, ctx));

   CHKERR(Num_malloc_SHprimme(ldhVecs*ldhVecs, &hVecs, ctx));
   CHKERR(Num_malloc_RHprimme(ldhVecs, &hVals, ctx));

   CHKERR(Num_malloc_iprimme(ldhVecs, &eval_perm, ctx));
   CHKERR(Num_malloc_iprimme(primme->numEvals, &flags, ctx));
   CHKERR(Num_malloc_iprimme(primme->numProcs, &global_start, ctx));

   for(i = 0; i < primme->numEvals; i++) flags[i] = UNCONVERGED;
  
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
   primme->stats.timeRR                        = 0.0;
   primme->stats.timeSketchMatvec              = 0.0;
   primme->stats.timeStabilization             = 0.0;
   primme->stats.timeResiduals                 = 0.0;
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

   /* Quick return for matrix of dimension 1 */
   if (primme->numEvals == 0) goto clean;

   /* Inserting random initial vectors into the basis ------------------------------------------ */
   blockSize = min(maxBlockSize, primme->maxBasisSize - maxBlockSize); /* Adjust block size first */

   // Set initial vectors
   for(i = 0; i < ldV*blockSize; i++) V[i] = 1.0;
   CHKERR(ortho_Sprimme(V, ldV, NULL, 0, 0, blockSize-1, NULL, 0, 0, primme->nLocal, primme->iseed, ctx));

   /* Build/update the sketching matrix if needed ------------------------------ */
   if(primme->projectionParams.projection == primme_proj_sketched)
   { 
      /* Default settings for sketch size and nnz per column. Based on Yuji and Tropp's manuscript */
      ldSV = ldSW = primme->sketchingParams.sketchSize;
      ldT = primme->maxBasisSize+maxBlockSize;
      nnzPerCol = primme->sketchingParams.nnzPerCol; 

      S_rows = (PRIMME_INT*)malloc(nnzPerCol*primme->nLocal*sizeof(PRIMME_INT));

      CHKERR(Num_malloc_Sprimme(nnzPerCol*primme->nLocal, &S_vals, ctx));
      CHKERR(Num_malloc_Sprimme(ldSV*(primme->maxBasisSize+maxBlockSize), &SV, ctx));
      CHKERR(Num_malloc_Sprimme(ldT*ldT, &T, ctx));
      CHKERR(Num_malloc_Sprimme(ldSW*primme->maxBasisSize, &SW, ctx));
      CHKERR(Num_malloc_Sprimme(ldV*(primme->maxBasisSize+maxBlockSize), &V_temp, ctx));

      CHKERR(Num_malloc_Rprimme(primme->numEvals, &normalize_evecs, ctx));

      /* Build Sketch CSR Locally */
      CHKERR(build_sketch_Sprimme(S_rows, S_vals, ctx));

   } /* End sketching matrix build */
   
   /* INITIAL ITERATION --------------------------------------------------- */
   CHKERR(matrixMatvec_Sprimme(V, ldV, primme->nLocal, &V[blockSize*ldV], ldV, 0, blockSize, ctx)); /* W = A*V_0 */

   /* Compute and subtract first alpha (W = W - V*(W'V))*/
   CHKERR(update_projection_Sprimme(V, ldV, &V[blockSize*ldV], ldV, H, ldH, primme->nLocal, 0, blockSize, 0, ctx)); 
   CHKERR(Num_gemm_Sprimme("N", "N", ldrwork, blockSize, blockSize, 1.0, V, ldV, H, ldH, 0.0, rwork, ldrwork, ctx)); 
   CHKERR(Num_axpy_Sprimme(ldV*blockSize, -1.0, rwork, 1, &V[ldV*blockSize], 1, ctx)); 

   /* Update sketched basis if sketching is turned on */
   if(primme->projectionParams.projection == primme_proj_sketched)
      CHKERR(sketch_basis_Sprimme(V, ldV, SV, ldSV, T, ldT, 0, blockSize, S_rows, S_vals, ctx));

   primme->stats.numOuterIterations++;

   /* ---------------------------------------------------------------
    * Main Lanczos Loop
    * --------------------------------------------------------------- */
   i = blockSize;
   while(i < primme->maxBasisSize && (primme->maxOuterIterations == 0 || primme->stats.numOuterIterations < primme->maxOuterIterations) && primme->stats.numMatvecs < primme->maxMatvecs) {

      blockSize = min(blockSize, primme->maxBasisSize - i); /* Adjust block size if needed */

      CHKERR(ortho_Sprimme(&V[i*ldV], ldV, &H[(i-blockSize)*ldH + i], ldH, 0, blockSize-1, NULL, 0, 0, primme->nLocal, primme->iseed, ctx));   /* [V_i, b_i] = qr(V_i) */
      if(fullOrtho)
         CHKERR(ortho_Sprimme(V, ldV, NULL, 0, i, i+blockSize-1, NULL, 0, 0, primme->nLocal, primme->iseed, ctx));   /* V_i = cgs(V(:, 0:i-1), V_i) */

      /* Symmetrize H */
      CHKERR(Num_copy_matrix_conj_Sprimme(&H[(i-blockSize)*ldH + i], blockSize, blockSize, ldH, &H[i*ldH + (i-blockSize)], ldH, ctx));

      CHKERR(matrixMatvec_Sprimme(&V[i*ldV], primme->nLocal, ldV, &V[(i+blockSize)*ldV], ldV, 0, blockSize, ctx));                             /* V_{i+1} = AV_i */

      /* Subtract beta term from new chunk of vectors (W_new = W - V_{i-1}*b_i')*/
      CHKERR(Num_gemm_Sprimme("N", "N", ldrwork, blockSize, blockSize, 1.0, &V[(i-blockSize)*ldV], ldV, &H[i*ldH + (i-blockSize)], ldH, 0.0, rwork, ldrwork, ctx));  
      CHKERR(Num_axpy_Sprimme(ldV*blockSize, -1.0, rwork, 1, &V[ldV*(i+blockSize)], 1, ctx)); 

      /* Find and subtract alpha term (W_new = W_new - V_i*a_i)*/
      CHKERR(update_projection_Sprimme(&V[i*ldV], ldV, &V[(i+blockSize)*ldV], ldV, &H[i*ldH+i], ldH, primme->nLocal, 0, blockSize, KIND(1 /*symmetric*/, 0 /* unsymmetric */), ctx)); /* a_i = W'*V_i */
      CHKERR(Num_gemm_Sprimme("N", "N", ldrwork, blockSize, blockSize, 1.0, &V[i*ldV], ldV, &H[i*ldH+i], ldH, 0.0, rwork, ldrwork, ctx)); /* rwork = V_i*a_i */
      CHKERR(Num_axpy_Sprimme(ldV*blockSize, -1.0, rwork, 1, &V[ldV*(i+blockSize)], 1, ctx)); /* W = W - rwork */

      /* Update sketched basis if sketching is turned on */
      if(primme->projectionParams.projection == primme_proj_sketched)
         CHKERR(sketch_basis_Sprimme(V, ldV, SV, ldSV, T, ldT, i, blockSize, S_rows, S_vals, ctx));

      if(primme->printLevel >= 2 && (PRIMME_INT)(i+blockSize) % 100 == 0)
      {
         /* Moving on to the eigenvalue problem */
         primme->initSize = 0;

         PRIMME_INT numEvals = min(i+blockSize, primme->numEvals);

         if(primme->projectionParams.projection == primme_proj_sketched)
         {
            /* Adding a row to H */
            CHKERR(Num_copy_matrix_Sprimme(V, ldV, i+2*blockSize, ldV, V_temp, ldV, ctx));
            CHKERR(ortho_Sprimme(&V_temp[ldV*(i+blockSize)], ldV, &H[i*ldH + (i+blockSize)], ldH, 0, blockSize-1, NULL, 0, 0, primme->nLocal, primme->iseed, ctx));   /* [V_i, b_i] = qr(V_i) */

            /* SW = SV*H */
            CHKERR(sketch_basis_Sprimme(V_temp, ldV, SV, ldSV, T, ldT, i+blockSize, blockSize, S_rows, S_vals, ctx));
            CHKERR(Num_gemm_Sprimme("N", "N", ldSV, i+blockSize, i+2*blockSize, 1.0, SV, ldSV, H, ldH, 0.0, SW, ldSW, ctx));

            /* Getting our sketched basis and projected sketched basis */
            CHKERR(sketched_RR_Sprimme(SV, ldSV, T, ldT, SW, ldSW, hVecs, i+blockSize, hVals, i+blockSize, ctx));

         } else { /* End sketching */ 
            solve_timer = primme_wTimer();
            CHKERR(solve_H_Sprimme(H, i+blockSize, ldH, NULL, 0, NULL, 0, NULL, 0, NULL, 0, NULL, 0, hVecs, i+blockSize, hVals, NULL, 0, ctx));
            primme->stats.timeRR += primme_wTimer() - solve_timer;
         } /* End non-sketching */

         /* Check how many pairs have converged basis in approximate residual */
         resid_timer = primme_wTimer();
         if(primme->procID == 0){
            for(j = 0; j < primme->numEvals; j++) {
               if(j < numEvals) {
                  resNorms[j] = fabs(H[(i-blockSize)*ldH + i]*hVecs[(i+blockSize)*j+(i+blockSize)-1]);
               } else {
                  resNorms[j] = 1.0;
               }
            }
         }
         primme->stats.timeResiduals += primme_wTimer() - resid_timer;
         numConverged = 0;
         CHKERR(broadcast_Rprimme(resNorms, numEvals, ctx));
         for(j = 0; j < primme->numEvals; j++) if(resNorms[j] < primme->aNorm*primme->eps) numConverged++;

         /* FOR TESTING: Print residual information */
         if(primme->procID == 0){
            printf("BasisSize %ld: NumConverged = %d, Convergence Tolerance = %.6E (%.6E x %.6E)\n", i+blockSize, numConverged, primme->aNorm*primme->eps, primme->aNorm, primme->eps);
            for(j = 0; j < primme->numEvals; j++){
               if(j < numEvals) {
                  printf("BasisSize %ld Eval[%ld] = %.10E, ResNorm[%ld] = %.10E\n", i+blockSize, j, hVals[j], j, resNorms[j]);
               } else {
                  printf("BasisSize %ld Eval[%ld] = NaN, ResNorm[%ld] = NaN\n", i+blockSize, j, j);
               }
            }
         }
         
         /* If everything is marked as converged, find the ACTUAL residuals and make sure */
         if(numConverged == primme->numEvals) {
            CHKERR(Num_gemm_Sprimme("N", "N", primme->nLocal, numEvals, i+blockSize, 1.0, V, ldV, hVecs, i+blockSize, 0.0, evecs, ldevecs, ctx));    /* evecs = V*hVecs */
            for(j = 0; j < primme->numEvals; j++) evals[j] = hVals[j];
           
            /* If sketching, scale the eigenvectors to ensure they are normal before computing their actual residual norms */ 
            if(primme->projectionParams.projection == primme_proj_sketched) {
               for(j = 0; j < primme->numEvals; j++) normalize_evecs[j] = Num_dot_Sprimme(primme->nLocal, &evecs[j*ldevecs], 1, &evecs[j*ldevecs], 1, ctx);
               CHKERR(globalSum_Rprimme(normalize_evecs, numEvals, ctx));
               for(j = 0; j < numEvals; j++) CHKERR(Num_scal_Sprimme(primme->nLocal, 1/sqrt(normalize_evecs[j]), &evecs[j*ldevecs], 1, ctx));
            }
            
            resid_timer = primme_wTimer(); 
            /* Compute residual vectors */
            CHKERR(matrixMatvec_Sprimme(evecs, primme->nLocal, ldevecs, AVhVecs, ldAVhVecs, 0, numEvals, ctx)); /* AVhVecs = A*V*hVecs */
            CHKERR(Num_compute_residuals_Sprimme(primme->nLocal, numEvals, evals, evecs, ldevecs, AVhVecs, ldAVhVecs, rwork, ldrwork, ctx)); /* Compute residual vectors */

            /* Compute residual norms */
            for(j = 0; j < numEvals; j++) resNorms[j] = Num_dot_Sprimme(ldrwork, &rwork[j*ldrwork], 1, &rwork[j*ldrwork], 1, ctx);
            CHKERR(globalSum_Rprimme(resNorms, numEvals, ctx));  
            for(j = 0; j < numEvals; j++) resNorms[j] = sqrt(resNorms[j]);
            primme->stats.timeResiduals += primme_wTimer() - resid_timer;
            
            /* Check the convergence of the Ritz vectors */
            CHKERR(check_convergence_Sprimme(evecs, ldevecs, 1 /* given X */, NULL, 0, 0 /* not given R */, NULL, 0, 0, NULL, 0, NULL, 0, 0, numEvals, flags, resNorms, hVals, &reset, -1, ctx));
            /* Find number of converged eigenpairs */
            numConverged = 0;
            for(j = 0; j < numEvals && flags[j] == CONVERGED; j++) numConverged = j+1;

            if(numConverged == primme->numEvals)
            {
               primme->initSize = numConverged;
               *numRet = numConverged;
               goto clean;
            } 

         } /* End if numConverged == primme->numEvals */

      } /* End the checking of residuals */

      /* Update loop variables */
      i += blockSize;
      primme->stats.numOuterIterations++;

      //Report timings
      if(primme->procID == 0 && (primme->stats.numOuterIterations % 100 == 0 || i % 100 == 0)){
         primme->stats.elapsedTime = primme_wTimer() - elapsed_time;  // XXX: Added this for debugging purposes - Heather
         CHKERR(print_lanczos_timings_Sprimme(i, ctx));
      }

   } /* End basis build */

   /**************************************************************************************
   * BASIS BUILD DONE.
   **************************************************************************************/

   /* If the basis build is done and not everything has converged -- move onto the eigenvalue problem */
   primme->initSize = 0;

   if(primme->projectionParams.projection == primme_proj_sketched)
   {
      /* Preparing SW matrix for the sketched Rayleigh-Ritz */
      CHKERR(ortho_Sprimme(&V[ldV*primme->maxBasisSize], ldV, &H[(primme->maxBasisSize-blockSize)*ldH + primme->maxBasisSize], ldH, 0, blockSize-1, NULL, 0, 0, primme->nLocal, primme->iseed, ctx));   /* [V_i, b_i] = qr(V_i) */
      CHKERR(ortho_Sprimme(V, ldV, NULL, 0, primme->maxBasisSize, primme->maxBasisSize+blockSize-1, NULL, 0, 0, primme->nLocal, primme->iseed, ctx));   /* Orthogonalized the last block of V against the rest of the basis */

      /* SW = SV*H */
      CHKERR(sketch_basis_Sprimme(V, ldV, SV, ldSV, T, ldT, primme->maxBasisSize, blockSize, S_rows, S_vals, ctx));
      CHKERR(Num_gemm_Sprimme("N", "N", ldSV, primme->maxBasisSize, primme->maxBasisSize+blockSize, 1.0, SV, ldSV, H, ldH, 0.0, SW, ldSW, ctx));

      /* Performing sketched Rayleigh-Ritz */
      CHKERR(sketched_RR_Sprimme(SV, ldSV, T, ldT, SW, ldSW, hVecs, ldhVecs, hVals, primme->maxBasisSize, ctx));

   } else { /* End sketched */
      solve_timer = primme_wTimer();
      CHKERR(solve_H_Sprimme(H, primme->maxBasisSize, ldH, NULL, 0, NULL, 0, NULL, 0, NULL, 0, NULL, 0, hVecs, ldhVecs, hVals, NULL, 0, ctx));
      primme->stats.timeRR += primme_wTimer() - solve_timer;
   } /* End nonsketched */

   CHKERR(Num_gemm_Sprimme("N", "N", primme->nLocal, primme->numEvals, primme->maxBasisSize, 1.0, V, ldV, hVecs, ldhVecs, 0.0, evecs, ldevecs, ctx));    /* evecs = V*hVecs */
   for(i = 0; i < primme->numEvals; i++) evals[i] = hVals[i];

   /* Ensure evecs are normal if sketching */
   if(primme->projectionParams.projection == primme_proj_sketched) {
      for(j = 0; j < primme->numEvals; j++) normalize_evecs[j] = Num_dot_Sprimme(primme->nLocal, &evecs[j*ldevecs], 1, &evecs[j*ldevecs], 1, ctx);
      CHKERR(globalSum_Rprimme(normalize_evecs, primme->numEvals, ctx));
      for(j = 0; j < primme->numEvals; j++) CHKERR(Num_scal_Sprimme(primme->nLocal, 1/sqrt(normalize_evecs[j]), &evecs[j*ldevecs], 1, ctx));
   }
   
   resid_timer = primme_wTimer();
   /* Compute residual vectors */
   CHKERR(matrixMatvec_Sprimme(evecs, primme->nLocal, ldevecs, AVhVecs, ldAVhVecs, 0, primme->numEvals, ctx)); /* AVhVecs = A*V*hVecs */
   CHKERR(Num_compute_residuals_Sprimme(primme->nLocal, primme->numEvals, evals, evecs, ldevecs, AVhVecs, ldAVhVecs, rwork, ldrwork, ctx));
   
   /* Compute residual norms */
   for(j = 0; j < primme->numEvals; j++) resNorms[j] = sqrt(Num_dot_Sprimme(ldrwork, &rwork[j*ldrwork], 1, &rwork[j*ldrwork], 1, ctx))/primme->numProcs;
   CHKERR(globalSum_Rprimme(resNorms, primme->numEvals, ctx));  
   primme->stats.timeResiduals += primme_wTimer() - resid_timer;

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
   CHKERR(Num_free_Sprimme(hVecs, ctx));

   CHKERR(Num_free_SHprimme(H, ctx));

   CHKERR(Num_free_RHprimme(hVals, ctx));

   CHKERR(Num_free_iprimme(flags, ctx));
   CHKERR(Num_free_iprimme(global_start, ctx));
   CHKERR(Num_free_iprimme(eval_perm, ctx));

   if(primme->projectionParams.projection == primme_proj_sketched)
   {
      CHKERR(Num_free_Sprimme(S_vals, ctx));
      CHKERR(Num_free_Sprimme(SV, ctx));
      CHKERR(Num_free_Sprimme(T, ctx));
      CHKERR(Num_free_Sprimme(SW, ctx));
      CHKERR(Num_free_Rprimme(normalize_evecs, ctx));
      CHKERR(Num_free_Sprimme(V_temp, ctx));
      free(S_rows);
   }

   *ret = 0;

return 0;

#else
   (void)evals;
   (void)evecs;
   (void)ldevecs;
   (void)resNorms;
   (void)ret;
   (void)numRet;
   (void)fullOrtho;
   (void)ctx;

   CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
   return 0;
#endif
}

#endif   // SUPPORT_TYPE
