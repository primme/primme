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

#include "magma_v2.h"
#include "magmasparse.h"

#include "primme.h"   /* header file for PRIMME SVDS too */ 

#ifndef min
#define min(A,B) ((A)<=(B)?(A):(B))
#endif
#ifndef max
#define max(A,B) ((A)>=(B)?(A):(B))
#endif

void magmaSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int *trans,
                    primme_svds_params *primme_svds, int *err);
void magmaDummy(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int *mode,
                    primme_svds_params *primme_svds, int *err);

int main (int argc, char *argv[]) {

   /* Solver arrays and parameters */
   double *svals;    /* Array with the computed singular values */
   double *rnorms;   /* Array with the computed residual norms */
   double *svecs;    /* Array with the computed singular vectors;
                        first right (v) vector starts in svecs[0],
                        second right (v) vector starts in svecs[primme_svd.n],
                        first left (u) vector starts in svecs[primme_svd.n*numSVals]...  */
   primme_svds_params primme_svds;
                     /* PRIMME SVDS configuration struct */

   /* Other miscellaneous items */
   int m=100,n=500; /* problem size */
   int ret;
   int i, j, k;
   double mu = 1e-5;

   int *col, *row;
   double *val;

   /* Write the Lauchli matrix in CSR

      [ 1  1  1  1  1 ...   1 ],  ei = 1 - (1 - mu)*i/(min(m,n) - 1)
      [e0  0  0  0  0 ...   0 ]
      [ 0 e1  0  0  0 ...   0 ]
      ...
      [ 0  0  0  0  0 ... en-1]
   */

   row = (int*) calloc(m+1, sizeof(int));
   col = (int*) calloc((m>0?n:0)+min(max(0,m-1),n), sizeof(int));
   val = (double*) calloc((m>0?n:0)+min(max(0,m-1),n), sizeof(double));

   for (i = j = 0; i < m; i++) {
      row[i] = j;
      if (i == 0) {
         for (k = 0; k < n; k++) {
            col[j] = k; val[j] = 1.0; j++;
         }
      }
      else if (i-1 < n) {
         col[j] = i-1; val[j] = 1.0 - (1.0-mu)*(i-1)/(min(m,n)-1); j++;
      }
   }
   row[m] = j;

   /* Initialize MAGMA and create some LA structures */
   magma_init();
   magma_queue_t queue;
   magma_queue_create(0, &queue);

   magma_d_matrix A={Magma_CSR}, At={Magma_CSR}, dA[2]={Magma_CSR, Magma_CSR};

   /* Pass the matrix to MAGMA and copy it to the GPU */
   magma_dcsrset(m, n, row, col, val, &A, queue);
   magma_dmtransposeconjugate(A, &At, queue);
   magma_dmtransfer(A, &dA[0], Magma_CPU, Magma_DEV, queue);
   magma_dmtransfer(At, &dA[1], Magma_CPU, Magma_DEV, queue);

   /* Set default values in PRIMME SVDS configuration struct */
   primme_svds_initialize(&primme_svds);

   /* Set problem matrix */
   primme_svds.matrixMatvec = magmaSparseMatrixMatvec;
   primme_svds.matrix = &dA;
                           /* Function that implements the matrix-vector products
                              A*x and A^t*x  */
  
   /* Set problem parameters */
   primme_svds.m = m;
   primme_svds.n = n; /* set problem dimension */
   primme_svds.numSvals = 4;   /* Number of wanted eigenpairs */
   primme_svds.eps = 1e-12;     /* ||r|| <= eps * ||matrix|| */
   primme_svds.target = primme_svds_smallest;
                               /* Seeking for the largest singular values  */

   /* Set preconditioner (optional) */
   //primme_svds.applyPreconditioner = magmaDummy;
   //primme_svds.primmeStage2.projectionParams.projection = primme_proj_RR;

   /* Set method to solve the singular value problem and
      the underneath eigenvalue problem (optional) */
   primme_svds_set_method(primme_svds_default, PRIMME_DEFAULT_MIN_TIME,
                              PRIMME_DEFAULT_MIN_TIME, &primme_svds);
   /*  primme_svds_default: devs choice, now being hybrid, which first solve
       the normal equation and then the augmented problem.
       PRIMME_DEFAULT_METHOD devs choice of the solver at every stage. But other methods
       can be set such as DYNAMIC or PRIMME_LOBPCG_OrthoBasis_Window. */

   primme_svds.printLevel = 5;


   /* Set advanced parameters if you know what are you doing (optional) */
   /* Configuration for 1st stage */
   /*
   or you can do:
   primme_svds.primme.maxBasisSize = 14;
   primme_svds.primme.minRestartSize = 6;
   primme_svds.primme.maxBlockSize = 2;
   */
   /* Configuration for 2nd stage */
   /*
   primme_svds.primmeStage2.maxBasisSize = 30;
   primme_svds.primmeStage2.minRestartSize = 15;
   primme_svds.primmeStage2.maxBlockSize = 1;
   */

    /* Display PRIMME SVDS configuration struct (optional) */
   primme_svds_display_params(primme_svds);

   /* Allocate space for converged Ritz values and residual norms */
   svals = (double*)malloc(primme_svds.numSvals*sizeof(double));
   if (magma_dmalloc(&svecs,(primme_svds.n+primme_svds.m)
         *primme_svds.numSvals) != MAGMA_SUCCESS) {return -1;};
   rnorms = (double*)malloc(primme_svds.numSvals*sizeof(double));

   primme_svds.queue = &queue;

   /* Call primme_svds  */
   ret = magma_dprimme_svds(svals, svecs, rnorms, &primme_svds);

   if (ret != 0) {
      fprintf(primme_svds.outputFile, 
         "Error: primme_svds returned with nonzero exit status: %d \n",ret);
      return -1;
   }


   /* Reporting (optional) */
   for (i=0; i < primme_svds.initSize; i++) {
      fprintf(primme_svds.outputFile, "Sval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
         svals[i], rnorms[i]); 
   }
   fprintf(primme_svds.outputFile, " %d eigenpairs converged\n", primme_svds.initSize);
   fprintf(primme_svds.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme_svds.aNorm*primme_svds.eps);
   fprintf(primme_svds.outputFile, "Iterations: %-" PRIMME_INT_P "\n", 
                                                 primme_svds.stats.numOuterIterations); 
   fprintf(primme_svds.outputFile, "Restarts  : %-" PRIMME_INT_P "\n", primme_svds.stats.numRestarts);
   fprintf(primme_svds.outputFile, "Matvecs   : %-" PRIMME_INT_P "\n", primme_svds.stats.numMatvecs);
   fprintf(primme_svds.outputFile, "Preconds  : %-" PRIMME_INT_P "\n", primme_svds.stats.numPreconds);
   fprintf(primme_svds.outputFile, "Time      : (total) %g (matvec) %g (precond) %g (ortho) %g\n", primme_svds.stats.elapsedTime, primme_svds.stats.timeMatvec, primme_svds.stats.timePrecond, primme_svds.stats.timeOrtho);

   primme_svds_free(&primme_svds);
   free(svals);
   free(rnorms);

  return(0);
}

void magmaSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int *trans,
                    primme_svds_params *primme_svds, int *err) {
   int i;            /* vector index, from 0 to *blockSize-1*/
   magma_d_matrix *A = primme_svds->matrix;
 
   if (!*trans) { // y <- A*x
      for (i=0; i<*blockSize; i++) {
         magma_d_matrix vx = {Magma_CSR};  /* i-th input vector x */
         magma_d_matrix vy = {Magma_CSR};  /* i-th output vector y */

         magma_dvset_dev(primme_svds->nLocal, 1, (double *)x + *ldx*i, &vx, *(magma_queue_t*)primme_svds->primme.queue);
         magma_dvset_dev(primme_svds->mLocal, 1, (double *)y + *ldy*i, &vy, *(magma_queue_t*)primme_svds->primme.queue);

         magma_d_spmv(1.0, A[0], vx, 0.0, vy, *(magma_queue_t*)primme_svds->primme.queue);
      }
   }
   else { // y <- A'*x
      for (i=0; i<*blockSize; i++) {
         magma_d_matrix vx = {Magma_CSR};  /* i-th input vector x */
         magma_d_matrix vy = {Magma_CSR};  /* i-th output vector y */

         magma_dvset_dev(primme_svds->mLocal, 1, (double *)x + *ldx*i, &vx, *(magma_queue_t*)primme_svds->primme.queue);
         magma_dvset_dev(primme_svds->nLocal, 1, (double *)y + *ldy*i, &vy, *(magma_queue_t*)primme_svds->primme.queue);

         magma_d_spmv(1.0, A[1], vx, 0.0, vy, *(magma_queue_t*)primme_svds->primme.queue);
      }
   }
   *err = 0;
}

void magmaDummy(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int *mode,
                    primme_svds_params *primme_svds, int *err) {
 
   int nrows=0;
   if (*mode == primme_svds_op_AtA) {
     nrows = primme_svds->nLocal;
   }
   else if (*mode == primme_svds_op_AAt) {
      nrows = primme_svds->mLocal;
   }
   else if (*mode == primme_svds_op_augmented) {
      nrows = primme_svds->mLocal + primme_svds->nLocal;
   }
   else {
      *err = -1;
      return;
   }

   magma_dcopymatrix(nrows, *blockSize, (double*)x, *ldx, (double*)y, *ldy, *(magma_queue_t*)primme_svds->primme.queue);
   *err = 0;
}
