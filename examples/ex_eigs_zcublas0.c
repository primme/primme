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
 *  with CUBLAS
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <cublas_v2.h>
#include <cusparse.h>

#include "primme.h"   /* header file is required to run primme */ 

void cuSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, primme_params *primme, int *ierr);

void checkCuda(cudaError_t err) {
   if (err != cudaSuccess) {
      fprintf(stderr, "cuda call failed!\n");
      exit(-1);
   }
}

void checkCublas(cublasStatus_t err) {
   if (err != CUBLAS_STATUS_SUCCESS) {
      fprintf(stderr, "cublas call failed!\n");
      exit(-1);
   }
}

void checkCusparse(cusparseStatus_t err) {
   if (err != CUSPARSE_STATUS_SUCCESS) {
      fprintf(stderr, "cusparse call failed!\n");
      exit(-1);
   }
}

int main (int argc, char *argv[]) {

   /* Solver arrays and parameters */
   double *evals;    /* Array with the computed eigenvalues */
   double *rnorms;   /* Array with the computed eigenpairs residual norms */
   complex double *evecs;    /* Array with the computed eigenvectors;
                        first vector starts in evecs[0],
                        second vector starts in evecs[primme.n],
                        third vector starts in evecs[primme.n*2]...  */
   primme_params primme;
                     /* PRIMME configuration struct */

   /* Other miscellaneous items */
   int n=64; /* problem size */
   int ret;
   int i;

   complex double *val;

   /* Create the matrix on cpu */
   val = (complex double*) calloc(n*n, sizeof(complex double));

   for (i = 0; i < n; i++) {
      val[i+n*i] = (complex double)i;
   }

   /* Copy the matrix on the gpu */
   cublasHandle_t cublas_handle;
   checkCublas(cublasCreate(&cublas_handle));
   complex double *val_dev;
   checkCuda(cudaMalloc((void**)&val_dev, n*n*sizeof(complex double)));
   checkCublas(cublasSetVector(n*n,sizeof(complex double), val, 1, val_dev, 1));
 
   /* Set default values in PRIMME configuration struct */
   primme_initialize(&primme);
 
   /* Set problem parameters */
   primme.n = n; /* set problem dimension */
   primme.numEvals = 6;   /* Number of wanted eigenpairs */
   primme.eps = 1e-12;      /* ||r|| <= eps * ||matrix|| */
   primme.target = primme_smallest;
   primme.maxBlockSize = 2;
                           /* Wanted the smallest eigenvalues */

   /* Set problem matrix */
   primme.matrixMatvec = cuSparseMatrixMatvec;
   primme.matrix = val_dev;
                           /* Function that implements the matrix-vector product
                              A*x for solving the problem A*x = l*x */
 
   /* Set preconditioner (optional) */
   /*
   primme.applyPreconditioner = magmaDummy;
   primme.correctionParams.precondition = 1;
   */

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
   checkCuda(cudaMalloc((void**)&evecs, primme.n*primme.numEvals*sizeof(complex double)));
   rnorms = (double*)malloc(primme.numEvals*sizeof(double));
   
   primme.queue = &cublas_handle;

   /* Call primme  */
   ret = cublas_zprimme(evals, evecs, rnorms, &primme);

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


   primme_free(&primme);
   free(val);
   free(evals);
   free(rnorms);

   // Free the allocated memory...
   checkCuda(cudaFree(evecs));
   checkCuda(cudaFree(val_dev));

   // and finalize cuBLAS and cuSparse.
   checkCublas(cublasDestroy(cublas_handle));

   return 0;
}

void cuSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, primme_params *primme, int *err) {

   cublasHandle_t cublas_handle = *(cublasHandle_t *)primme->queue;
   complex double alpha = 1.0, beta = 0;
   checkCublas(cublasGemmEx(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N, primme->n,
         *blockSize, primme->n, &alpha, (const void *)primme->matrix,
         CUDA_C_64F, primme->n, x, CUDA_C_64F, *ldx, &beta, y, CUDA_C_64F, *ldy,
         CUBLAS_COMPUTE_64F, CUBLAS_GEMM_DEFAULT));
   *err = 0;
}
