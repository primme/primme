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
 *  with HIPBLAS
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <hipblas.h>
#include <hipsparse.h>

#include "primme.h"   /* header file is required to run primme */ 

void hipSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, primme_params *primme, int *ierr);

typedef struct {
   hipsparseHandle_t hipsparse_handle;
   hipsparseSpMatDescr_t desc;
   void *aux;
   size_t aux_size;
} MatrixInfo;

void checkHip(hipError_t err) {
   if (err != hipSuccess) {
      fprintf(stderr, "hip call failed!\n");
      exit(-1);
   }
}

void checkHipblas(hipblasStatus_t err) {
   if (err != HIPBLAS_STATUS_SUCCESS) {
      fprintf(stderr, "hipblas call failed!\n");
      exit(-1);
   }
}

void checkHipsparse(hipsparseStatus_t err) {
   if (err != HIPSPARSE_STATUS_SUCCESS) {
      fprintf(stderr, "hipsparse call failed!\n");
      exit(-1);
   }
}

int main (int argc, char *argv[]) {
   (void)argc;
   (void)argv;

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
   int n=1000; /* problem size */
   int ret;
   int i,j;

   int *col, *row;
   complex double *val;

   /* Create the matrix on cpu */
   row = (int*) calloc(n+1, sizeof(int));
   int nnz = n+(n>0?n-1:0)*2;
   col = (int*) calloc(nnz, sizeof(int));
   val = (complex double*) calloc(nnz, sizeof(complex double));

   for (i = j = 0; i < n; i++) {
      row[i] = j;
      if (i > 0)   {col[j] = i-1; val[j] = -1.0; j++;}
                    col[j] = i  ; val[j] =  2.0; j++;
      if (i < n-1) {col[j] = i+1; val[j] = -1.0; j++;}
   }
   row[n] = j;

   /* Copy the matrix on the gpu */
   int *col_dev, *row_dev;
   complex double *val_dev;
   checkHip(hipMalloc((void**)&row_dev, (n+1)*sizeof(int)));
   checkHip(hipMalloc((void**)&col_dev, nnz*sizeof(int)));
   checkHip(hipMalloc((void**)&val_dev, nnz*sizeof(complex double)));
   checkHipblas(hipblasSetVector(n+1,sizeof(int), row, 1, row_dev, 1));
   checkHipblas(hipblasSetVector(nnz,sizeof(int), col, 1, col_dev, 1));
   checkHipblas(hipblasSetVector(nnz,sizeof(complex double), val, 1, val_dev, 1));
   MatrixInfo A;
   checkHipsparse(hipsparseCreate(&A.hipsparse_handle));
   checkHipsparse(hipsparseCreateCsr(&A.desc, n, n, nnz, row_dev, col_dev, val_dev,
         HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_32I, HIPSPARSE_INDEX_BASE_ZERO,
         HIP_C_64F));
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
   primme.matrixMatvec = hipSparseMatrixMatvec;
   primme.matrix = &A;
                           /* Function that implements the matrix-vector product
                              A*x for solving the problem A*x = l*x */
   primme.printLevel = 5;
 
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
   checkHip(hipMalloc((void**)&evecs, primme.n*primme.numEvals*sizeof(complex double)));
   rnorms = (double*)malloc(primme.numEvals*sizeof(double));
   
   hipblasHandle_t hipblas_handle;
   checkHipblas(hipblasCreate(&hipblas_handle));
   primme.queue = &hipblas_handle;

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


//   printf("Time used in func:%e\n",primme.funcTime);
//  printf("Execution time: %e\n",time);

   primme_free(&primme);
   free(row);
   free(col);
   free(val);
   free(evals);
   free(rnorms);

   // Free the allocated memory...
   checkHipsparse(hipsparseDestroySpMat(A.desc));
   checkHip(hipFree(evecs));
   checkHip(hipFree(row_dev));
   checkHip(hipFree(col_dev));
   checkHip(hipFree(val_dev));
   if (A.aux_size > 0) checkHip(hipFree(A.aux));

   // and finalize cuBLAS and cuSparse.
   checkHipsparse(hipsparseDestroy(A.hipsparse_handle));
   checkHipblas(hipblasDestroy(hipblas_handle));

   return 0;
}

void hipSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, primme_params *primme, int *err) {

   MatrixInfo *A = (MatrixInfo*)primme->matrix;
   hipsparseDnMatDescr_t matx, maty;
   checkHipsparse(hipsparseCreateDnMat(&matx, primme->nLocal, *blockSize, *ldx, x, HIP_C_64F,
         HIPSPARSE_ORDER_COLUMN));
   checkHipsparse(hipsparseCreateDnMat(&maty, primme->nLocal, *blockSize, *ldy, y, HIP_C_64F,
         HIPSPARSE_ORDER_COLUMN));
   double alpha = 1.0, beta = 0;
   size_t buffer_size = 0;
   checkHipsparse(hipsparseSpMM_bufferSize(A->hipsparse_handle, HIPSPARSE_OPERATION_NON_TRANSPOSE,
         HIPSPARSE_OPERATION_NON_TRANSPOSE, &alpha, A->desc, matx, &beta, maty,
         HIP_C_64F, HIPSPARSE_SPMM_ALG_DEFAULT, &buffer_size));
   if (buffer_size > A->aux_size) {
      if (A->aux) checkHip(hipFree(A->aux));
      checkHip(hipMalloc(&A->aux, buffer_size));
      A->aux_size = buffer_size;
   }
   checkHipsparse(hipsparseSpMM(A->hipsparse_handle, HIPSPARSE_OPERATION_NON_TRANSPOSE,
         HIPSPARSE_OPERATION_NON_TRANSPOSE, &alpha, A->desc, matx, &beta, maty,
         HIP_C_64F, HIPSPARSE_SPMM_ALG_DEFAULT, A->aux));
   checkHipsparse(hipsparseDestroyDnMat(matx));
   checkHipsparse(hipsparseDestroyDnMat(maty));
   *err = 0;
}
