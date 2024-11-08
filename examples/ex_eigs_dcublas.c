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
#include <assert.h>
#include <string.h>
#include <math.h>
#include <cublas_v2.h>
#include <cusparse.h>

#include "primme.h"   /* header file is required to run primme */ 

void cuSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, primme_params *primme, int *ierr);

typedef struct {
   cusparseHandle_t cusparse_handle;
   cusparseSpMatDescr_t desc;
   void *aux;
   size_t aux_size;
} MatrixInfo;

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
   (void)argc;
   (void)argv;

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

   /* Create the matrix on cpu */
   row = (int*) calloc(n+1, sizeof(int));
   int nnz = n+(n>0?n-1:0)*2;
   col = (int*) calloc(nnz, sizeof(int));
   val = (double*) calloc(nnz, sizeof(double));

   for (i = j = 0; i < n; i++) {
      row[i] = j;
      if (i > 0)   {col[j] = i-1; val[j] = -1.0; j++;}
                    col[j] = i  ; val[j] =  2.0; j++;
      if (i < n-1) {col[j] = i+1; val[j] = -1.0; j++;}
   }
   row[n] = j;

   /* Copy the matrix on the gpu */
   int *col_dev, *row_dev;
   double *val_dev;
   checkCuda(cudaMalloc((void**)&row_dev, (n+1)*sizeof(int)));
   checkCuda(cudaMalloc((void**)&col_dev, nnz*sizeof(int)));
   checkCuda(cudaMalloc((void**)&val_dev, nnz*sizeof(double)));
   checkCublas(cublasSetVector(n+1,sizeof(int), row, 1, row_dev, 1));
   checkCublas(cublasSetVector(nnz,sizeof(int), col, 1, col_dev, 1));
   checkCublas(cublasSetVector(nnz,sizeof(double), val, 1, val_dev, 1));
   MatrixInfo A;
   checkCusparse(cusparseCreate(&A.cusparse_handle));
   checkCusparse(cusparseCreateCsr(&A.desc, n, n, nnz, row_dev, col_dev, val_dev,
         CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO,
         CUDA_R_64F));
   A.aux = NULL;
   A.aux_size = 0;
 
   /* Set default values in PRIMME configuration struct */
   primme_initialize(&primme);
 
   /* Set problem parameters */
   primme.n = n; /* set problem dimension */
   primme.numEvals = 6;   /* Number of wanted eigenpairs */
   primme.eps = 1e-12;      /* ||r|| <= eps * ||matrix|| */
   primme.target = primme_smallest;
                           /* Wanted the smallest eigenvalues */

   /* Set problem matrix */
   primme.matrixMatvec = cuSparseMatrixMatvec;
   primme.matrix = &A;
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
   checkCuda(cudaMalloc((void**)&evecs, primme.n*primme.numEvals*sizeof(double)));
   rnorms = (double*)malloc(primme.numEvals*sizeof(double));
   
   cublasHandle_t cublas_handle;
   checkCublas(cublasCreate(&cublas_handle));
   primme.queue = &cublas_handle;

   /* Call primme  */
   ret = cublas_dprimme(evals, evecs, rnorms, &primme);

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
   checkCusparse(cusparseDestroySpMat(A.desc));
   checkCuda(cudaFree(evecs));
   checkCuda(cudaFree(row_dev));
   checkCuda(cudaFree(col_dev));
   checkCuda(cudaFree(val_dev));
   if (A.aux_size > 0) checkCuda(cudaFree(A.aux));

   // and finalize cuBLAS and cuSparse.
   checkCusparse(cusparseDestroy(A.cusparse_handle));
   checkCublas(cublasDestroy(cublas_handle));

   return 0;
}

void cuSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, primme_params *primme, int *err) {

   MatrixInfo *A = (MatrixInfo*)primme->matrix;
   cusparseDnMatDescr_t matx, maty;
   checkCusparse(cusparseCreateDnMat(&matx, primme->nLocal, *blockSize, *ldx, x, CUDA_R_64F,
         CUSPARSE_ORDER_COL));
   checkCusparse(cusparseCreateDnMat(&maty, primme->nLocal, *blockSize, *ldy, y, CUDA_R_64F,
         CUSPARSE_ORDER_COL));
   double alpha = 1.0, beta = 0;
   size_t buffer_size = 0;
   checkCusparse(cusparseSpMM_bufferSize(A->cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
         CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, A->desc, matx, &beta, maty,
         CUDA_R_64F, CUSPARSE_SPMM_ALG_DEFAULT, &buffer_size));
   if (buffer_size > A->aux_size) {
      if (A->aux) checkCuda(cudaFree(A->aux));
      checkCuda(cudaMalloc(&A->aux, buffer_size));
      A->aux_size = buffer_size;
   }
   checkCusparse(cusparseSpMM(A->cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
         CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, A->desc, matx, &beta, maty,
         CUDA_R_64F, CUSPARSE_SPMM_ALG_DEFAULT, A->aux));
   checkCusparse(cusparseDestroyDnMat(matx));
   checkCusparse(cusparseDestroyDnMat(maty));
   *err = 0;
}
