/*******************************************************************************
 * Copyright (c) 2016, College of William & Mary
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
 *
 *  Example to compute the k largest singular values in Lauchli matrix.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include "primme.h"   /* header file for PRIMME SVDS too */ 

#ifndef min
#define min(A,B) ((A)<=(B)?(A):(B))
#endif
#ifndef max
#define max(A,B) ((A)>=(B)?(A):(B))
#endif


void LauchliMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                         int *transpose, primme_svds_params *primme_svds, int *ierr);
void LauchliApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                                int *mode, primme_svds_params *primme_svds, int *ierr);
void LauchliAugmentedMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);

int main (int argc, char *argv[]) {

   /* Solver arrays and parameters */
   double *svals;    /* Array with the computed singular values */
   double *rnorms;   /* Array with the computed residual norms */
   complex double *svecs;    /* Array with the computed singular vectors;
                        first right (v) vector starts in svecs[0],
                        second right (v) vector starts in svecs[primme_svd.n],
                        first left (u) vector starts in svecs[primme_svd.n*numSVals]...  */
   primme_svds_params primme_svds;
                     /* PRIMME SVDS configuration struct */

   /* Other miscellaneous items */
   int ret;
   int i;
   double mu = 1e-5;

   /* Set default values in PRIMME SVDS configuration struct */
   primme_svds_initialize(&primme_svds);

   /* Set problem matrix */
   primme_svds.matrixMatvec = LauchliMatrixMatvec;
   primme_svds.matrix = &mu;
                           /* Function that implements the matrix-vector products
                              A*x and A^t*x  */
  
   /* Set problem parameters */
   primme_svds.m = 500;
   primme_svds.n = 100; /* set problem dimension */
   primme_svds.numSvals = 4;   /* Number of wanted singular values */
   primme_svds.eps = 1e-12;     /* ||r|| <= eps * ||matrix|| */
   primme_svds.target = primme_svds_smallest;
                               /* Seeking for the largest singular values  */

   /* Set preconditioner (optional) */
   primme_svds.applyPreconditioner = LauchliApplyPreconditioner;

   /* Set method to solve the singular value problem and
      the underneath eigenvalue problem (optional) */
   primme_svds_set_method(primme_svds_hybrid, PRIMME_DYNAMIC, PRIMME_DEFAULT_MIN_TIME, &primme_svds);
   /*  Set hybrid method with PRIMME_DYNAMIC and PRIMME_DEFAULT_MIN_TIME as the underneath eigensolver configuration
       for the first and the second stage, respectively.
       PRIMME_DYNAMIC uses a runtime heuristic to choose the fastest method between
       PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can set another
       method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */

   primme_svds.printLevel = 3;

   /* Set parameters for the underneath eigensolver if you know what are you doing (optional) */
   primme_svds.primme.locking = 1; 
   primme_svds.primme.restartingParams.maxPrevRetain = 3;

    /* Display PRIMME SVDS configuration struct (optional) */
   primme_svds_display_params(primme_svds);

   /* Allocate space for converged Ritz values and residual norms */
   svals = (double*)malloc(primme_svds.numSvals*sizeof(double));
   svecs = (complex double*)malloc((primme_svds.n+primme_svds.m)
         *primme_svds.numSvals*sizeof(complex double));
   rnorms = (double*)malloc(primme_svds.numSvals*sizeof(double));

   /* Call primme_svds  */
   ret = zprimme_svds(svals, svecs, rnorms, &primme_svds);

   if (ret != 0) {
      fprintf(primme_svds.outputFile, 
         "Error: primme_svds returned with nonzero exit status: %d \n",ret);
      return -1;
   }


   /* Hybrid method second stage */
   primme_svds_set_method(primme_svds_augmented, PRIMME_DEFAULT_MIN_MATVECS,
                          PRIMME_DEFAULT_METHOD, &primme_svds);
                              /* Set DEFAULT_MIN_MATVECS as the underneath eigensolver */
   primme_svds.primme.matrixMatvec = LauchliAugmentedMatvec; 
                              /* Set a custom matrix-vector product */

   /* Reporting (optional) */
   for (i=0; i < primme_svds.initSize; i++) {
      fprintf(primme_svds.outputFile, "Sval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
         svals[i], rnorms[i]); 
   }
   fprintf(primme_svds.outputFile, " %d singular triplets converged\n", primme_svds.initSize);
   fprintf(primme_svds.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme_svds.aNorm*primme_svds.eps);
   fprintf(primme_svds.outputFile, "Iterations: %-" PRIMME_INT_P "\n", 
                                                 primme_svds.stats.numOuterIterations); 
   fprintf(primme_svds.outputFile, "Restarts  : %-" PRIMME_INT_P "\n", primme_svds.stats.numRestarts);
   fprintf(primme_svds.outputFile, "Matvecs   : %-" PRIMME_INT_P "\n", primme_svds.stats.numMatvecs);
   fprintf(primme_svds.outputFile, "Preconds  : %-" PRIMME_INT_P "\n", primme_svds.stats.numPreconds);
   if (primme_svds.primme.locking && primme_svds.primme.intWork && primme_svds.primme.intWork[0] == 1) {
      fprintf(primme_svds.outputFile, "\nA locking problem has occurred.\n"
         "Some triplets do not have a residual norm less than the tolerance.\n"
         "However, the subspace of evecs is accurate to the required tolerance.\n");
   }


   primme_svds_free(&primme_svds);
   free(svals);
   free(svecs);
   free(rnorms);

  return(0);
}

/* lauchli block matrix-vector product, y = a * x (or y = a^t * x), where

   - x, input dense matrix of size primme_svds.n (or primme_svds.m) x blocksize;
   - y, output dense matrix of size primme_svds.m (or primme_svds.n) x blocksize;
   - a, lauchli matrix of dimensions primme_svds.m x (primme_svds.m+1) with this form:

        [ 1  1  1  1  1 ...   1 ],  ei = 1 - (1 - mu)*i/(min(m,n) - 1)
        [e0  0  0  0  0 ...   0 ]
        [ 0 e1  0  0  0 ...   0 ]
         ...
        [ 0  0  0  0  0 ... en-1]
*/

void LauchliMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                         int *transpose, primme_svds_params *primme_svds, int *err) {
   
   int i;            /* vector index, from 0 to *blockSize-1 */
   int j;
   int min_m_n = min(primme_svds->m, primme_svds->n);
   complex double *xvec;     /* pointer to i-th input vector x */
   complex double *yvec;     /* pointer to i-th output vector y */
   double mu = *(double*)primme_svds->matrix;

   if (*transpose == 0) { /* Do y <- A * x */
      for (i=0; i<*blockSize; i++) { 
         xvec = (complex double *)x + (*ldx)*i;
         yvec = (complex double *)y + (*ldy)*i;
         yvec[0] = 0;
         for (j=0; j<primme_svds->n; j++) {
            yvec[0] += xvec[j];
         }
         for (j=1; j<primme_svds->m; j++) {
            yvec[j] = j-1<primme_svds->n ? xvec[j-1]*(1.0 - (1.0 - mu)*(j-1)/(min_m_n - 1)) : 0.0;
         }      
      }
   } else { /* Do y <- A^t * x */
      for (i=0; i<*blockSize; i++) {
         xvec = (complex double *)x + (*ldx)*i;
         yvec = (complex double *)y + (*ldy)*i;
         for (j=0; j<primme_svds->n; j++) {
            yvec[j] = xvec[0];
            if (j+1 < primme_svds->m) yvec[j] += xvec[j+1]*(1.0 - (1.0 - mu)*j/(min_m_n - 1));
         }
      }
   }
   *err = 0;
}

/* Example of custom product for the augmented matrix [0 A^t; A 0]. In this case
   the direct and the transpose matrix-vector products are taken from
   LauchliMatrixMatvec.
*/

void LauchliAugmentedMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {
   /* A pointer to the primme_svds_params is stored at primme.matrix */
   primme_svds_params *primme_svds = (primme_svds_params*)primme->matrix;
   complex double *x0 = (complex double*)x, *x1 = &x0[primme_svds->nLocal],
              *y0 = (complex double*)y, *y1 = &y0[primme_svds->nLocal];
   int notrans=0, trans=1;
   /* [y0; y1] <-  * [x0; x1] */

   /* y0 <- A^t * x1 */
   LauchliMatrixMatvec(x1, ldx, y0, ldy, blockSize, &trans, primme_svds, ierr);
   /* y1 <- A * x0 */
   LauchliMatrixMatvec(x0, ldx, y1, ldy, blockSize, &notrans, primme_svds, ierr);
}

/* This performs Y = M^{-1} * X, where

   - X, input dense matrix of size primme_svds.n (or primme_svds.m or m+n) x blockSize;
   - Y, output dense matrix of size primme_svds.n (or primme_svds.m or m+n) x blockSize;
   - M, preconditioner for A^t*A (or A*A^t or [0 A^t; A 0]), where A is the Lauchli matrix.
*/

void LauchliApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                                int *mode, primme_svds_params *primme_svds, int *ierr) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int j;            /* row index */
   complex double *xvec;     /* pointer to i-th input vector x */
   complex double *yvec;     /* pointer to i-th output vector y */
   int modeAtA = primme_svds_op_AtA, modeAAt = primme_svds_op_AAt;
   double mu = *(double*)primme_svds->matrix;
   complex double  *aux;
   PRIMME_INT ldaux;
   int notrans = 0, trans = 1;
   int min_m_n = min(primme_svds->m, primme_svds->n);
    
   if (*mode == primme_svds_op_AtA) {
      /* Preconditioner for A^t*A, diag(A^t*A)^{-1} */
      for (i=0; i<*blockSize; i++) { 
         xvec = (complex double *)x + (*ldx)*i;
         yvec = (complex double *)y + (*ldy)*i;
         for (j=0; j<primme_svds->n; j++) {
            double ei = j<primme_svds->m ? 1.0 - (1.0 - mu)*j/(min_m_n - 1) : 0.0;
            yvec[j] = xvec[j]/(1.0 + ei*ei);
         }      
      }
   }
   else if (*mode == primme_svds_op_AAt) {
      /* Preconditioner for A*A^t, diag(A*A^t)^{-1} */
      for (i=0; i<*blockSize; i++) {
         xvec = (complex double *)x + (*ldx)*i;
         yvec = (complex double *)y + (*ldy)*i;
         yvec[0] = xvec[0]/(complex double)primme_svds->m;
         for (j=1; j<primme_svds->m; j++) {
            double ei = j<primme_svds->n ? 1.0 - (1.0 - mu)*j/(min_m_n - 1) : 1.0;
            yvec[j] = xvec[j]/ei/ei;
         }
      }
   }
   else if (*mode == primme_svds_op_augmented) {
      /* Preconditioner for [0 A^t; A 0],
         [diag(A^t*A) 0; 0 diag(A*A^t)]^{-1}*[0 A^t; A 0] */

      /* [y0; y1] <- [0 A^t; A 0] * [x0; x1] */
      ldaux = primme_svds->n+primme_svds->m;
      aux = (complex double*)malloc(sizeof(complex double)*(*blockSize)*ldaux);
      primme_svds->matrixMatvec(x, ldx, &aux[primme_svds->n], &ldaux, blockSize, &notrans, primme_svds, ierr);
      xvec = (complex double *)x + primme_svds->n;
      primme_svds->matrixMatvec(xvec, ldx, aux, &ldaux, blockSize, &trans, primme_svds, ierr);
      /* y0 <- preconditioner for A^t*A  * y0 */
      LauchliApplyPreconditioner(aux, &ldaux, y, ldy, blockSize, &modeAtA, primme_svds, ierr);
      /* y1 <- preconditioner for A*A^t  * y1 */
      yvec = (complex double *)y + primme_svds->n;
      LauchliApplyPreconditioner(&aux[primme_svds->n], &ldaux, yvec, ldy, blockSize, &modeAAt, primme_svds, ierr);
      free(aux);
   }
   *ierr = 0;
}
