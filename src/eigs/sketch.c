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
 * File: projection.c
 *
 * Purpose - This builds a random subspace embedding
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/sketch.c"
#endif

#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "sketch.h"
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
 * Subroutine apply_sketching - This routine using a random subspace embedding
 * to estimate the eigenpairs of the H matrix produced by Lanczos
 * I.e. [evecs, evals] = eig((SV)'(SAV)y = Ly)
 *
 * INPUT arrays and parameters
 * ----------------------------------
 * H        The tridiagonal matrix produced by Lanczos
 * 
 * ldH      The leading dimension of H
 *
 * V        The Lanczos basis
 * 
 * ldV      The leading dimension of V
 * 
 * ldSV     The leading dimension of SV
 *
 * ldhVecs  The leading dimension of hVecs 
 * 
 * basisSize The size of the basis
 *
 * blockSize The number of columns added to the basis per iteration in Lanczos
 * 
 * nnzPerCol Number of nonzeros per column of the sketching matrix
 * 
 * S_rows   The nonzero rows indices for all columns in the sketching matrix
 *          in CSC format
 * 
 * S_vals   The values for all nonzero spots in the sketching matrix in CSC format 
 * 
 * INPUT/OUTPUT arrays and parameters
 * ----------------------------------
 * SV       The sketched basis with dimensions ldS x ldV
 * 
 * hVecs    The eigenvectors of the sketched RR problem
 *
 * hVals    The sketched Ritz values
 * 
 * resNorms The estimated residual norms from the sketched RR
 * 
 * last_sketch The number of columns of SV that have previously been computed
 *  
 ******************************************************************************/


TEMPLATE_PLEASE
int apply_sketching_Sprimme(HSCALAR *H, PRIMME_INT ldH, SCALAR *V, PRIMME_INT ldV, SCALAR *SV, PRIMME_INT ldSV, HSCALAR *hVecs, PRIMME_INT ldhVecs, HREAL *hVals, PRIMME_INT last_sketch, PRIMME_INT basisSize, PRIMME_INT blockSize, PRIMME_INT *S_rows, SCALAR* S_vals, PRIMME_INT nnzPerCol, primme_context ctx) {
#ifdef USE_HERMITIAN
   primme_params *primme = ctx.primme;

   SCALAR *V_row;    /* Used to temporarily store a row in V to avoid numberous memory accesses */
   PRIMME_INT i, j;
   PRIMME_INT numNewCols = basisSize+blockSize-last_sketch;  

   CHKERR(Num_malloc_Sprimme(numNewCols, &V_row, ctx));
   CHKERR(Num_zero_matrix_Sprimme(&SV[last_sketch*ldSV], ldSV, numNewCols, ldSV, ctx));

   /* Sparse MM */
   for(i = 0; i < primme->nLocal; i++) /* Traverse the rows of the basis V */
   {
      CHKERR(Num_copy_Sprimme(numNewCols, &V[last_sketch*ldV+i], ldV, V_row, 1, ctx));
      for(j = 0; j < nnzPerCol; j++)
      {
         CHKERR(Num_axpy_Sprimme(numNewCols, S_vals[i*nnzPerCol+j], V_row, 1, &SV[last_sketch*ldSV + S_rows[i*nnzPerCol+j]], ldSV, ctx));
      }
   }
   CHKERR(Num_free_Sprimme(V_row, ctx)); 


   /* Find the sketched basis */
   CHKERR(globalSum_Sprimme(&SV[ldSV*last_sketch], numNewCols*ldSV, ctx));  
/*
   if(primme->procID == 0 && primme->numProcs == 1)
   {
      printf("SV for basisSize %d\n", basisSize);
      for(i = 0; i < ldSV; i++)
      {
         for(j = 0; j < basisSize+blockSize; j++)
            printf("%lf ", SV[j*ldSV+i]);
         printf("\n");
      }
      printf("\n\n");
   }
*/
   if(primme->procID == 0)
   {

      SCALAR *SW;       /* The projected sketched basis */
      SCALAR *UtSW;     /* Temporary matrix to compute left hand side of the generalized eigenvalue problem */
      SCALAR *UtSWV;    /* Left hand side of the generalized eigenvalue problem */
      SCALAR *UVecs;    /* Left singular vectors of the SVD */
      SCALAR *VVecst;   /* Right singular vectors of the SVD (transposed) */
      SCALAR *ShVals;   /* The "alpha" eigenvalues returned from LaPack's GGEV */ 
      SCALAR *hVals_b;  /* The "beta" eigenvalues returned from LaPack's GGEV */ 
      SCALAR *trunc_hVecs;  /* The stabilized eigenvectors of H */ 
      SCALAR *Sigma;    /* Right hand side of the generalized eigenvalue problem */
      SCALAR *temp_SV;  /* Temporarily store SV while performing LaPack's GESVD */
      REAL *sing_vals;  /* Returned singular values from the SVD */
      int *eval_perm;   /* To sort the eigenpairs */
      PRIMME_INT ldSW, ldVVecst, ldUVecs, ldUtSW, ldUtSWV, ldSigma, ldtrunc_hVecs;   /* Leading dimension of matrices */
      PRIMME_INT trunc_basisSize; /* The basis size after stabilization */

      ldUVecs = ldSW = ldSV;
      ldVVecst = basisSize;

      CHKERR(Num_malloc_Sprimme(ldSW*basisSize, &SW, ctx));
      CHKERR(Num_malloc_Sprimme(ldUVecs*basisSize, &UVecs, ctx));
      CHKERR(Num_malloc_Sprimme(ldVVecst*basisSize, &VVecst, ctx));
      CHKERR(Num_malloc_Sprimme(ldSV*basisSize, &temp_SV, ctx));
      CHKERR(Num_malloc_Rprimme(basisSize, &sing_vals, ctx));

      /* Project the sketched basis (SW = SV*H)*/
      CHKERR(Num_gemm_Sprimme("N", "N", ldSV, basisSize, basisSize+blockSize, 1.0, SV, ldSV, H, ldH, 0.0, SW, ldSW, ctx));

      /* Take the SVD decomposition of SV */
      CHKERR(Num_copy_matrix_Sprimme(&SV[0], ldSV, basisSize, ldSV, &temp_SV[0], ldSV, ctx));
      CHKERR(Num_gesvd_Sprimme("S", "S", ldSV, basisSize, SV, ldSV, sing_vals, UVecs, ldUVecs, VVecst, ldVVecst, ctx));
      CHKERR(Num_copy_matrix_Sprimme(&temp_SV[0], ldSV, basisSize, ldSV, &SV[0], ldSV, ctx));
      CHKERR(Num_free_Sprimme(temp_SV, ctx)); 

      /* Truncate the basis for stabilization if needed */
      trunc_basisSize = basisSize;
      for(i = basisSize-1; i > primme->numEvals; i--)
      {
         if(sing_vals[0]/sing_vals[i] >= 1/MACHINE_EPSILON)
         {
            trunc_basisSize--;
         } else {
            break;
         }
      }

      /* Build each side of the generalized eigenvalue problem after stabilization */
      ldSigma = ldUtSW = ldUtSWV = ldtrunc_hVecs = trunc_basisSize;
      CHKERR(Num_malloc_Sprimme(ldUtSW*basisSize, &UtSW, ctx));
      CHKERR(Num_malloc_Sprimme(ldUtSWV*trunc_basisSize, &UtSWV, ctx));
      CHKERR(Num_malloc_Sprimme(ldtrunc_hVecs*trunc_basisSize, &trunc_hVecs, ctx));
      CHKERR(Num_malloc_Sprimme(trunc_basisSize, &ShVals, ctx));
      CHKERR(Num_malloc_Sprimme(trunc_basisSize, &hVals_b, ctx));
      CHKERR(Num_malloc_Sprimme(ldSigma*trunc_basisSize, &Sigma, ctx));
      CHKERR(Num_malloc_iprimme(trunc_basisSize, &eval_perm, ctx));

      /* Construct Mass Matrix (Diagonal matrix with singular values as entries) */
      CHKERR(Num_zero_matrix_Sprimme(Sigma, trunc_basisSize, trunc_basisSize, ldSigma, ctx));
      for(i = 0; i < trunc_basisSize; i++) Sigma[i*ldSigma+i] = sing_vals[i];

      /* Left hand side of the Eigenvalue problem (U'(SW)Vt'), where Vt is the transpose of the right singular vectors */
      CHKERR(Num_gemm_Sprimme("C", "N", trunc_basisSize, basisSize, ldSV, 1.0, UVecs, ldUVecs, SW, ldSW, 0.0, UtSW, ldUtSW, ctx)); 
      CHKERR(Num_gemm_Sprimme("N", "C", trunc_basisSize, trunc_basisSize, basisSize, 1.0, UtSW, ldUtSW, VVecst, ldVVecst, 0.0, UtSWV, ldUtSWV, ctx)); 
      
      /* Eigenvalue problem */
      CHKERR(Num_ggev_Sprimme("N", "V", trunc_basisSize, UtSWV, ldUtSWV, Sigma, ldSigma, ShVals, NULL, hVals_b, NULL, trunc_basisSize, trunc_hVecs, ldtrunc_hVecs, ctx)); 
      CHKERR(Num_gemm_Sprimme("C", "N", basisSize, trunc_basisSize, trunc_basisSize, 1.0, VVecst, ldVVecst, trunc_hVecs, ldtrunc_hVecs, 0.0, hVecs, ldhVecs, ctx)); 

      for(i = 0; i < trunc_basisSize; i++)
      {
         hVals[i] = REAL_PART(ShVals[i]/hVals_b[i]);
         eval_perm[i] = i;
      }

      /* Sort the eigenpairs */
      for(i = 0; i < trunc_basisSize; i++) CHKERR(insertionSort_Rprimme(hVals[i], hVals, 0.0, NULL, 0, NULL, eval_perm, i, 0, ctx.primme));

      CHKERR(permute_vecs_Sprimme(hVecs, basisSize, trunc_basisSize, ldhVecs, eval_perm, ctx));

      CHKERR(Num_free_Sprimme(SW, ctx)); 
      CHKERR(Num_free_Sprimme(UVecs, ctx)); 
      CHKERR(Num_free_Sprimme(VVecst, ctx)); 
      CHKERR(Num_free_Sprimme(UtSW, ctx)); 
      CHKERR(Num_free_Sprimme(UtSWV, ctx)); 
      CHKERR(Num_free_Sprimme(trunc_hVecs, ctx)); 
      CHKERR(Num_free_Sprimme(ShVals, ctx)); 
      CHKERR(Num_free_Sprimme(hVals_b, ctx)); 
      CHKERR(Num_free_Sprimme(Sigma, ctx)); 
      CHKERR(Num_free_iprimme(eval_perm, ctx)); 
      CHKERR(Num_free_Rprimme(sing_vals, ctx)); 
   }

   //PRIMME_INT numEvals = min(primme->numEvals, basisSize);
   CHKERR(broadcast_SHprimme(hVecs, ldhVecs*primme->numEvals, ctx));
   CHKERR(broadcast_RHprimme(hVals, primme->numEvals, ctx));

   /* Form residual estimates */
/*   SCALAR *SVhVecs, *SWhVecs;
   SCALAR *sketched_residuals;
   REAL   *resNorms_normalize;
   CHKERR(Num_malloc_Sprimme(ldSV*numEvals, &sketched_residuals, ctx));
   CHKERR(Num_malloc_Sprimme(ldSV*numEvals, &SVhVecs, ctx));
   CHKERR(Num_malloc_Sprimme(ldSV*numEvals, &SWhVecs, ctx));
   CHKERR(Num_malloc_Rprimme(numEvals, &resNorms_normalize, ctx));
*/   
   /* Compute sketched Ritz vectors (SV*hVecs = SVhVecs) and projected sketched sketched Ritz vectors (SW*hVecs = SAVhVecs) */
/*   CHKERR(Num_gemm_Sprimme("N", "N", ldSV, numEvals, basisSize, 1.0, SV, ldSV, hVecs, ldhVecs, 0.0, SVhVecs, ldSV, ctx));
   CHKERR(Num_gemm_Sprimme("N", "N", ldSV, numEvals, basisSize, 1.0, SW, ldSV, hVecs, ldhVecs, 0.0, SWhVecs, ldSV, ctx));

    for(i = 0; i < numEvals; i++)
      hVals[i] = Num_dot_Sprimme(ldSV, &SVhVecs[i*ldSV], 1, &SWhVecs[i*ldSV], 1, ctx) / Num_dot_Sprimme(ldSV, &SVhVecs[i*ldSV], 1, &SVhVecs[i*ldSV], 1, ctx);

   CHKERR(Num_compute_residuals_Sprimme(ldSV, numEvals, hVals, SVhVecs, ldSV, SWhVecs, ldSV, sketched_residuals, ldSV, ctx));

   for(i = 0; i < numEvals; i++)
   {
      resNorms[i] = Num_dot_Sprimme(ldSV, &sketched_residuals[i*ldSV], 1, &sketched_residuals[i*ldSV], 1, ctx);
      resNorms_normalize[i] = Num_dot_Sprimme(ldSV, &SVhVecs[i*ldSV], 1, &SVhVecs[i*ldSV], 1, ctx);
   }
   CHKERR(globalSum_RHprimme(resNorms, numEvals, ctx));  
   CHKERR(globalSum_Rprimme(resNorms_normalize, numEvals, ctx));  

   for(i = 0; i < numEvals; i++) resNorms[i] = (sqrt(resNorms[i]) / sqrt(resNorms_normalize[i]));

   CHKERR(Num_free_Sprimme(SWhVecs, ctx)); 
   CHKERR(Num_free_Sprimme(SVhVecs, ctx)); 
   CHKERR(Num_free_Sprimme(sketched_residuals, ctx)); 
   CHKERR(Num_free_Rprimme(resNorms_normalize, ctx)); 
*/

return 0;
#else
   (void)H;
   (void)ldH;
   (void)V;
   (void)ldV;
   (void)SV;
   (void)ldSV;
   (void)hVecs;
   (void)ldhVecs;
   (void)hVals;
   (void)last_sketch;
   (void)basisSize;
   (void)blockSize;
   (void)S_rows;
   (void)S_vals;
   (void)nnzPerCol;
   (void)ctx;

   CHKERR(PRIMME_FUNCTION_UNAVAILABLE);
   return 0;
#endif
}

#endif   // SUPPORT_TYPE
