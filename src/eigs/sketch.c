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

/* Swaps two integer variables */
STATIC void swap_iprimme(int *a, int *b){
    int temp = *a;
    *a = *b;
    *b = temp;
}

/* Randomly shuffle an array of integers */
STATIC void random_shuffle_iprimme(int *x, int n) {
    int i, j;
    for(i = n-1; i > 0; i--) 
    {
        j = rand() % (i+1);
        swap_iprimme(&x[i], &x[j]);
    }
}


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
 * SV       The sketched basis with dimensions sketchSize x ldV
 ******************************************************************************/


TEMPLATE_PLEASE
int apply_sketching_Sprimme(SCALAR *H, PRIMME_INT ldH, SCALAR *V, PRIMME_INT ldV, SCALAR *hVecs, PRIMME_INT ldhVecs, SCALAR *hVals, PRIMME_INT basisSize, PRIMME_INT blockSize, primme_context ctx) {

   primme_params *primme = ctx.primme;

   SCALAR *SV;       /* The sketched basis */
   SCALAR *SW;       /* The projected sketched basis */
   SCALAR *QtSW;     /* The left hand side of the generalized eigenvalue problem */
   SCALAR *Q;        /* The orthononormal matrix from the QR decomposition of SV */
   SCALAR *R;        /* The right hand side of the generalized eigenvalue problem */
   SCALAR *S_vals;   /* The nonzero values of the sketching matrix */ 
   SCALAR *hVals_b;  /* The "beta" eigenvalues returned from LaPack's GGEV */ 

   int *S_rows;     /* Store the row index of each nonzero */
   int *rand_rows;  /* Vector used to randomly select nonzero rows for each column */

   PRIMME_INT i, j;                          /* Loop variables */
   PRIMME_INT ldSV, ldSW, ldQtSW, ldQ, ldR;  /* Leading dimension of matrices */
   PRIMME_INT nLocal = primme->nLocal;       /* Number of local rows */

   PRIMME_INT sketchSize = 4*basisSize;
   PRIMME_INT nnzPerCol = ceil(2*log(basisSize+1));

   /* Allocate space */
   ldSV = ldSW = ldQ = sketchSize;
   ldQtSW = ldR = basisSize;
   CHKERR(Num_malloc_Sprimme(ldSV*(basisSize+blockSize), &SV, ctx));
   CHKERR(Num_malloc_Sprimme(ldSW*basisSize, &SW, ctx));
   CHKERR(Num_malloc_Sprimme(ldQtSW*basisSize, &QtSW, ctx));
   CHKERR(Num_malloc_Sprimme(ldQ*basisSize, &Q, ctx));
   CHKERR(Num_malloc_Sprimme(ldR*basisSize, &R, ctx));
   CHKERR(Num_malloc_Sprimme(nLocal*nnzPerCol, &S_vals, ctx));
   CHKERR(Num_malloc_Sprimme(basisSize, &hVals_b, ctx));

   CHKERR(Num_malloc_iprimme(nLocal*nnzPerCol, &S_rows, ctx));
   CHKERR(Num_malloc_iprimme(sketchSize, &rand_rows, ctx));

   /* -------------------------------------------------------------------------------
    * Build and apply sketching matrix to the basis V
    *--------------------------------------------------------------------------------*/

   /* Insert random variables into the nonzeros of the skecthing matrix (Uniform on the complex unit circle or -1 and 1) */
   CHKERR(Num_larnv_Sprimme(2, primme->iseed, primme->nLocal*nnzPerCol, &S_vals[0], ctx));
   for(i = 0; i < primme->nLocal*nnzPerCol; i++) S_vals[i] /= (sqrt(pow(REAL_PART(S_vals[i]), 2)+pow(IMAGINARY_PART(S_vals[i]), 2))*sqrt(sketchSize));

   /* Select nnzperCol random rows per column to be nonzero */
   for(i = 0; i < sketchSize; i++) rand_rows[i] = i;
   for(i = 0; i < primme->nLocal; i++)
   {
      random_shuffle_iprimme(rand_rows, sketchSize);
      for(j = 0; j < nnzPerCol; j++) S_rows[i*nnzPerCol + j] = rand_rows[j];
   }

   /* Sketching matrix built. Moving on to applying it. */
   // TODO: THIS IS TEMPORARY -------------------------------------------------------

   SCALAR *S; /* Holds sketching matrix */
   PRIMME_INT ldS; /* Leading dimension of S */
   ldS = sketchSize;  
 
   CHKERR(Num_malloc_Sprimme(sketchSize*primme->nLocal, &S, ctx));
   CHKERR(Num_zero_matrix_Sprimme(S, sketchSize, ldV, ldS, ctx));

   for(i = 0; i < nnzPerCol*primme->nLocal; i++) S[sketchSize*((int)(i/nnzPerCol)) + S_rows[i]] = S_vals[i]; 

   CHKERR(Num_gemm_Sprimme("N", "N", sketchSize, basisSize+blockSize, nLocal, 1.0, S, ldS, V, ldV, 0.0, SV, ldSV, ctx)); 

   /* Project the sketched basis (SW = SV*H)*/
   CHKERR(Num_gemm_Sprimme("N", "N", sketchSize, basisSize, basisSize+blockSize, 1.0, SV, ldSV, H, ldH, 0.0, SW, ldSW, ctx));

   /* Build the lhs and rhs matrixes for the generalized eigenvalue problem */
   CHKERR(Num_copy_Sprimme(sketchSize*basisSize, SV, 1, Q, 1, ctx));
   CHKERR(ortho_Sprimme(Q, ldQ, R, ldR, 0, basisSize-1, NULL, 0, 0, sketchSize, primme->iseed, ctx));   /* QR of the sketched basis */

   CHKERR(Num_gemm_Sprimme("C", "N", basisSize, basisSize, sketchSize, 1.0, Q, ldQ, SW, ldSW, 0.0, QtSW, ldQtSW, ctx)); /* Left side matrix */
   
   CHKERR(Num_ggev_Sprimme("N", "V", basisSize, QtSW, ldQtSW, R, ldR, hVals, NULL, hVals_b, NULL, ldhVecs, hVecs, ldhVecs, ctx)); /* Solve Q'SWx = RLx */
   for(i = 0; i < basisSize; i++) hVals[i] = hVals[i]/hVals_b[i];

   /* Cleaning up */
   CHKERR(Num_free_Sprimme(S, ctx));
   CHKERR(Num_free_Sprimme(SV, ctx));
   CHKERR(Num_free_Sprimme(SW, ctx));
   CHKERR(Num_free_Sprimme(QtSW, ctx));
   CHKERR(Num_free_Sprimme(Q, ctx));
   CHKERR(Num_free_Sprimme(R, ctx));
   CHKERR(Num_free_Sprimme(S_vals, ctx));
   CHKERR(Num_free_Sprimme(hVals_b, ctx));

   CHKERR(Num_free_iprimme(S_rows, ctx));
   CHKERR(Num_free_iprimme(rand_rows, ctx));

return 0;
}

#endif   // SUPPORT_TYPE
