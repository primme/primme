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
 * Subroutine apply_sketching - This routine applies a random subspace embedding
 * to a basis matrix V, and W = AV using Sparse Maps
 *
 * INPUT arrays and parameters
 * ----------------------------------
 * V        Store a basis V
 * 
 * ldV      The leading dimension of V
 * 
 * basisSize: The size of the basis
 * 
 * nnzPerCol: Number of nonzeros per column of the matrix
 * 
 * INPUT/OUTPUT arrays and parameters
 * ----------------------------------
 * SV       The sketched basis with dimensions ldS x ldV
 ******************************************************************************/


TEMPLATE_PLEASE
int apply_sketching_Sprimme(SCALAR *H, PRIMME_INT ldH, SCALAR *V, PRIMME_INT ldV, SCALAR *SV, PRIMME_INT ldSV, SCALAR *hVecs, PRIMME_INT ldhVecs, REAL *hVals, PRIMME_INT *basisSize, PRIMME_INT blockSize, int *S_rows, SCALAR* S_vals, PRIMME_INT ldS, PRIMME_INT nnzPerCol, primme_context ctx) {

   primme_params *primme = ctx.primme;

   PRIMME_INT i, j;

   SCALAR *V_row;
   CHKERR(Num_malloc_Sprimme(blockSize, &V_row, ctx));

   for(i = 0; i < primme->nLocal; i++)
   {
      CHKERR(Num_copy_Sprimme(blockSize, &V[(*basisSize)*ldV+i], ldV, V_row, 1, ctx));
      for(j = 0; j < nnzPerCol; j++) CHKERR(Num_axpy_Sprimme(blockSize, S_vals[i*nnzPerCol+j], &V_row[0], 1, &SV[((*basisSize)*ldSV) + S_rows[i*nnzPerCol+j]], ldSV, ctx));
   }
   CHKERR(Num_free_Sprimme(V_row, ctx)); 

   /* Find the sketched basis */
   CHKERR(globalSum_Sprimme(&SV[ldSV*(*basisSize)], ldSV*blockSize, ctx));  

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
      REAL *sing_vals;  /* Returned singular values from the SVD */
      int *eval_perm;   /* To sort the eigenpairs */
      PRIMME_INT ldSW, ldVVecst, ldUVecs, ldUtSW, ldUtSWV, ldSigma, ldtrunc_hVecs;   /* Leading dimension of matrices */
      PRIMME_INT trunc_basisSize; /* The basis size after stabilization */

      ldSW = ldUVecs = ldS;
      ldVVecst = (*basisSize);

      CHKERR(Num_malloc_Sprimme(ldSW*(*basisSize), &SW, ctx));
      CHKERR(Num_malloc_Sprimme(ldUVecs*(*basisSize), &UVecs, ctx));
      CHKERR(Num_malloc_Sprimme(ldVVecst*(*basisSize), &VVecst, ctx));
      CHKERR(Num_malloc_Rprimme((*basisSize), &sing_vals, ctx));

      /* Project the sketched basis (SW = SV*H)*/
      CHKERR(Num_gemm_Sprimme("N", "N", ldS, (*basisSize), (*basisSize+blockSize), 1.0, SV, ldSV, H, ldH, 0.0, SW, ldSW, ctx));

      /* Take the SVD decomposition of SV */
      CHKERR(Num_gesvd_Sprimme("S", "A", ldS, (*basisSize), SV, ldSV, sing_vals, UVecs, ldUVecs, VVecst, ldVVecst, ctx));

      trunc_basisSize = (*basisSize);
      for(i = (*basisSize)-1; i >= primme->numEvals; i--)
      {
         if(sing_vals[0]/sing_vals[i] > 1/MACHINE_EPSILON)
         {
            trunc_basisSize--;
         } else {
            break;
         }
      }

      /* Build each side of the generalized eigenvalue problem after stabilization */
      ldSigma = ldUtSW = ldUtSWV = ldtrunc_hVecs = trunc_basisSize;
      CHKERR(Num_malloc_Sprimme(ldUtSW*(*basisSize), &UtSW, ctx));
      CHKERR(Num_malloc_Sprimme(ldUtSWV*trunc_basisSize, &UtSWV, ctx));
      CHKERR(Num_malloc_Sprimme(ldtrunc_hVecs*trunc_basisSize, &trunc_hVecs, ctx));
      CHKERR(Num_malloc_Sprimme(trunc_basisSize, &ShVals, ctx));
      CHKERR(Num_malloc_Sprimme(trunc_basisSize, &hVals_b, ctx));
      CHKERR(Num_malloc_Sprimme(ldSigma*trunc_basisSize, &Sigma, ctx));
      CHKERR(Num_malloc_iprimme(trunc_basisSize, &eval_perm, ctx));

      /* Right hand side */
      CHKERR(Num_zero_matrix_Sprimme(Sigma, trunc_basisSize, trunc_basisSize, ldSigma, ctx));
      for(i = 0; i < trunc_basisSize; i++) Sigma[i*ldSigma+i] = sing_vals[i];

      /* Left hand side */
      CHKERR(Num_gemm_Sprimme("C", "N", trunc_basisSize, (*basisSize), ldS, 1.0, UVecs, ldUVecs, SW, ldSW, 0.0, UtSW, ldUtSW, ctx)); /* Left side matrix */
      CHKERR(Num_gemm_Sprimme("N", "C", trunc_basisSize, trunc_basisSize, (*basisSize), 1.0, UtSW, ldUtSW, VVecst, ldVVecst, 0.0, UtSWV, ldUtSWV, ctx)); /* Left side matrix */
      
      /* Eigenvalue problem */
      CHKERR(Num_ggev_Sprimme("N", "V", trunc_basisSize, UtSWV, ldUtSWV, Sigma, ldSigma, ShVals, NULL, hVals_b, NULL, trunc_basisSize, trunc_hVecs, ldtrunc_hVecs, ctx)); /* Solve Q'SWx = RLx */
      CHKERR(Num_gemm_Sprimme("C", "N", (*basisSize), trunc_basisSize, trunc_basisSize, 1.0, VVecst, ldVVecst, trunc_hVecs, ldtrunc_hVecs, 0.0, hVecs, ldhVecs, ctx)); /* Left side matrix */
      for(i = 0; i < trunc_basisSize; i++)
      {
         hVals[i] = REAL_PART(ShVals[i]/hVals_b[i]);
         eval_perm[i] = i;
      }

      /* Sort the eigenpairs */
      for(i = 0; i < trunc_basisSize; i++) CHKERR(insertionSort_Sprimme(hVals[i], hVals, 0.0, NULL, 0, NULL, eval_perm, i, 0, ctx.primme));
      CHKERR(permute_vecs_Sprimme(hVecs, (*basisSize), trunc_basisSize, ldhVecs, eval_perm, ctx));

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

   CHKERR(broadcast_Sprimme(hVecs, ldhVecs*primme->numEvals, ctx));
   CHKERR(broadcast_Rprimme(hVals, primme->numEvals, ctx));

   (*basisSize) += blockSize; 

return 0;
}

#endif   // SUPPORT_TYPE
