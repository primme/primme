/*******************************************************************************
 * Copyright (c) 2017, College of William & Mary
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
 *  Example to compute the k largest eigenvalues in a 1-D Laplacian matrix
 *  with MAGMA
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "magma_v2.h"
#include "magmasparse.h"

#include "primme.h"   /* header file is required to run primme */ 

#include <time.h>

void magmaSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
void magmaDummy(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);


int main (int argc, char *argv[]) {

   /* Solver arrays and parameters */
   double *evals;    /* Array with the computed eigenvalues */
   double *rnorms;   /* Array with the computed eigenpairs residual norms */
   double *evecs;    /* Array with the computed eigenvectors;
                        first vector starts in evecs[0],
                        second vector starts in evecs[primme.n],
                        third vector starts in evecs[primme.n*2]...  */
   primme_params primme;
                     /* PRIMME configuration struct */

   /* Other miscellaneous items */
   int n=1000; /* problem size */
   int ret;
   int i,j;

   int *col, *row;
   double *val;

   row = (int*) calloc(n+1, sizeof(int));
   col = (int*) calloc(n+(n>0?n-1:0)*2, sizeof(int));
   val = (double*) calloc(n+(n>0?n-1:0)*2, sizeof(double));

   for (i = j = 0; i < n; i++) {
      row[i] = j;
      if (i > 0)   {col[j] = i-1; val[j] = -1.0; j++;}
                    col[j] = i  ; val[j] =  2.0; j++;
      if (i < n-1) {col[j] = i+1; val[j] = -1.0; j++;}
   }
   row[n] = j;

   /* Initialize MAGMA and create some LA structures */
   magma_init();
   magma_queue_t queue;
   magma_queue_create(0, &queue);

   magma_d_matrix A={Magma_CSR}, dA={Magma_CSR};

   /* Pass the matrix to MAGMA and copy it to the GPU */
   magma_dcsrset(n, n, row, col, val, &A, queue);
   magma_dmtransfer(A, &dA, Magma_CPU, Magma_DEV, queue);

   /* Set default values in PRIMME configuration struct */
   primme_initialize(&primme);
 
   /* Set problem parameters */
   primme.n = n; /* set problem dimension */
   primme.numEvals = 6;   /* Number of wanted eigenpairs */
   primme.eps = 1e-12;      /* ||r|| <= eps * ||matrix|| */
   primme.target = primme_smallest;
                           /* Wanted the smallest eigenvalues */

   /* Set problem matrix */
   primme.matrixMatvec = magmaSparseMatrixMatvec;
   primme.matrix = &dA;
                           /* Function that implements the matrix-vector product
                              A*x for solving the problem A*x = l*x */
 
   /* Set preconditioner (optional) */
   primme.applyPreconditioner = magmaDummy;
   primme.correctionParams.precondition = 1;

   /* Set advanced parameters if you know what are you doing (optional) */
   /*
   primme.maxBasisSize = 14;
   primme.minRestartSize = 4;
   primme.maxBlockSize = 1;
   primme.maxMatvecs = 1000;
   */

   /* Set method to solve the problem */
   primme_set_method(PRIMME_DYNAMIC, &primme);
//   primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
//   primme_set_method(PRIMME_DEFAULT_MIN_TIME, &primme);
   /* DYNAMIC uses a runtime heuristic to choose the fastest method between
       PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
       set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */

   /* Display PRIMME configuration struct (optional) */
   primme_display_params(primme);

   /* Allocate space for converged Ritz values and residual norms */
   evals = (double*)malloc(primme.numEvals*sizeof(double));
   magma_dmalloc(&evecs, primme.n*primme.numEvals);
   rnorms = (double*)malloc(primme.numEvals*sizeof(double));

   primme.queue = &queue;

/*
   clock_t start,end;

   start = clock();
   primme.funcTime = 0;
*/
  time_t rawtime,rawtime2;
  struct tm * timeinfo,* timeinfo2;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  printf ( "Start ->Current local time and date: %s", asctime (timeinfo) );


   /* Call primme  */
   ret = magma_dprimme(evals, evecs, rnorms, &primme);

  time ( &rawtime2 );
  timeinfo2 = localtime ( &rawtime2 );
  printf ( "End   ->Current local time and date: %s", asctime (timeinfo2) );

/*
   end = clock();
   double time = (double)(end-start)/CLOCKS_PER_SEC;
*/


   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   /* Reporting (optional) */
   for (i=0; i < primme.initSize; i++) {
      fprintf(primme.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
         evals[i], rnorms[i]); 
   }
   fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
   fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme.aNorm*primme.eps);
   fprintf(primme.outputFile, "Iterations: %-" PRIMME_INT_P "\n", 
                                                 primme.stats.numOuterIterations); 
   fprintf(primme.outputFile, "Restarts  : %-" PRIMME_INT_P "\n", primme.stats.numRestarts);
   fprintf(primme.outputFile, "Matvecs   : %-" PRIMME_INT_P "\n", primme.stats.numMatvecs);
   fprintf(primme.outputFile, "Preconds  : %-" PRIMME_INT_P "\n", primme.stats.numPreconds);

   switch (primme.dynamicMethodSwitch) {
      case -1: fprintf(primme.outputFile,
            "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
      case -2: fprintf(primme.outputFile,
            "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
      case -3: fprintf(primme.outputFile,
            "Recommended method for next run: DYNAMIC (close call)\n"); break;
   }


//   printf("Time used in func:%e\n",primme.funcTime);
//  printf("Execution time: %e\n",time);

   primme_free(&primme);
   free(row);
   free(col);
   free(val);
   free(evals);
   free(rnorms);

   // Free the allocated memory...
   magma_free(evecs);
   magma_dmfree( &dA, queue );

   // and finalize MAGMA.
   magma_queue_destroy( queue );
   magma_finalize();

   return 0;
}

void magmaSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   double *xvec;     
   double *yvec;     
   magma_d_matrix *A = primme->matrix;
 
   for (i=0; i<*blockSize; i++) {
      magma_d_matrix vx = {Magma_CSR};  /* i-th input vector x */
      magma_d_matrix vy = {Magma_CSR};  /* i-th output vector y */

      magma_dvset_dev(primme->n, 1, (double *)x + *ldx*i, &vx, *(magma_queue_t*)primme->queue);
      magma_dvset_dev(primme->n, 1, (double *)y + *ldy*i, &vy, *(magma_queue_t*)primme->queue);

      magma_d_spmv(1.0, *A, vx, 0.0, vy, *(magma_queue_t*)primme->queue);
   }
   *err = 0;
}

void magmaDummy(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   magma_dcopymatrix(primme->n, *blockSize, (double*)x, *ldx, (double*)y, *ldy, *(magma_queue_t*)primme->queue);
   *err = 0;
}
