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
 *
 *  Example to compute interior values in a normal operator.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "primme.h"   /* header file is required to run primme */
typedef PRIMME_COMPLEX_DOUBLE cdouble;

#ifdef __cplusplus
#define creal(x) std::real(x)
#define cimag(x) std::imag(x)
#define cabs(x) std::abs(x)
#define conj(x) std::conj(x)
#define I cdouble(0,1)
#endif

void LaplacianLikeMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
void LaplacianLikeApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);

int main (int argc, char *argv[]) {

   /* Solver arrays and parameters */
   cdouble *evals;    /* Array with the computed eigenvalues */
   double *rnorms;   /* Array with the computed eigenpairs residual norms */
   cdouble *evecs;    /* Array with the computed eigenvectors;
                        first vector starts in evecs[0],
                        second vector starts in evecs[primme.n],
                        third vector starts in evecs[primme.n*2]...  */
   primme_params primme;
                     /* PRIMME configuration struct */
   double targetShifts[1];

   /* Other miscellaneous items */
   int ret;
   int i;

   /* Set default values in PRIMME configuration struct */
   primme_initialize(&primme);

   /* Set problem matrix */
   primme.matrixMatvec = LaplacianLikeMatrixMatvec;
                           /* Function that implements the matrix-vector product
                              A*x for solving the problem A*x = l*x */
  
   /* Set problem parameters */
   primme.n = 100; /* set problem dimension */
   primme.numEvals = 10;   /* Number of wanted eigenpairs */
   primme.eps = 1e-9;      /* ||r|| <= eps * ||matrix|| */
   targetShifts[0] = .5;
   primme.targetShifts = targetShifts;
   primme.numTargetShifts = 1;
   primme.target = primme_closest_abs;
                           /* Wanted the smallest eigenvalues */

   /* Set preconditioner (optional) */
   primme.applyPreconditioner = LaplacianLikeApplyPreconditioner;
   primme.correctionParams.precondition = 1;

   /* Set advanced parameters if you know what are you doing (optional) */
   /*
   primme.maxBasisSize = 14;
   primme.minRestartSize = 4;
   primme.maxBlockSize = 1;
   primme.maxMatvecs = 1000;
   */

   /* Set method to solve the problem */
   /* NOTE: PRIMME_DEFAULT_MIN_TIME is not supported normal operators */
   primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
   /* You can set other methods, such as PRIMME_LOBPCG_OrthoBasis_Window */

   /* NOTE: cheap Olsen approximation is not supported for normal operators */
   primme.correctionParams.projectors.RightX = 0;

   primme.printLevel = 3;

   /* Display PRIMME configuration struct (optional) */
   primme_display_params(primme);

   /* Allocate space for converged Ritz values and residual norms */
   evals = (cdouble*)malloc(primme.numEvals*sizeof(cdouble));
   evecs = (cdouble*)malloc(primme.n*primme.numEvals*sizeof(cdouble));
   rnorms = (double*)malloc(primme.numEvals*sizeof(double));

   /* Call primme  */
   ret = zprimme_normal(evals, evecs, rnorms, &primme);

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   /* Reporting (optional) */
   for (i=0; i < primme.initSize; i++) {
      fprintf(primme.outputFile, "Eval[%d]: %-22.15E + %-22.15E I rnorm: %-22.15E\n", i+1,
         creal(evals[i]), cimag(evals[i]), rnorms[i]); 
   }
   fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
   fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme.aNorm*primme.eps);
   fprintf(primme.outputFile, "Iterations: %-" PRIMME_INT_P "\n", 
                                                 primme.stats.numOuterIterations); 
   fprintf(primme.outputFile, "Restarts  : %-" PRIMME_INT_P "\n", primme.stats.numRestarts);
   fprintf(primme.outputFile, "Matvecs   : %-" PRIMME_INT_P "\n", primme.stats.numMatvecs);
   fprintf(primme.outputFile, "Preconds  : %-" PRIMME_INT_P "\n", primme.stats.numPreconds);
   if (primme.stats.lockingIssue) {
      fprintf(primme.outputFile, "\nA locking problem has occurred.\n");
      fprintf(primme.outputFile,
         "Some eigenpairs do not have a residual norm less than the tolerance.\n");
      fprintf(primme.outputFile,
         "However, the subspace of evecs is accurate to the required tolerance.\n");
   }

   switch (primme.dynamicMethodSwitch) {
      case -1: fprintf(primme.outputFile,
            "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
      case -2: fprintf(primme.outputFile,
            "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
      case -3: fprintf(primme.outputFile,
            "Recommended method for next run: DYNAMIC (close call)\n"); break;
   }


   /* Note that d/zprimme can be called more than once before call primme_free. */
   /* Find the 5 eigenpairs closest to .5 */
   primme.numTargetShifts = 1;
   targetShifts[0] = .5;
   primme.targetShifts = targetShifts;
   primme.target = primme_closest_abs;
   primme.numEvals = 5;
   primme.initSize = 0; /* primme.initSize may be not zero after a d/zprimme;
                           so set it to zero to avoid the already converged eigenvectors
                           being used as initial vectors. */

   /* Call primme  */
   ret = zprimme_normal(evals, evecs, rnorms, &primme);

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   /* Reporting (optional) */
   for (i=0; i < primme.initSize; i++) {
      fprintf(primme.outputFile, "Eval[%d]: %-22.15E + %-22.15E I rnorm: %-22.15E\n", i+1,
         creal(evals[i]), cimag(evals[i]), rnorms[i]); 
   }
   fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
   fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme.aNorm*primme.eps);
   fprintf(primme.outputFile, "Iterations: %-" PRIMME_INT_P "\n", 
                                                 primme.stats.numOuterIterations); 
   fprintf(primme.outputFile, "Restarts  : %-" PRIMME_INT_P "\n", primme.stats.numRestarts);
   fprintf(primme.outputFile, "Matvecs   : %-" PRIMME_INT_P "\n", primme.stats.numMatvecs);
   fprintf(primme.outputFile, "Preconds  : %-" PRIMME_INT_P "\n", primme.stats.numPreconds);
   if (primme.stats.lockingIssue) {
      fprintf(primme.outputFile, "\nA locking problem has occurred.\n");
      fprintf(primme.outputFile,
         "Some eigenpairs do not have a residual norm less than the tolerance.\n");
      fprintf(primme.outputFile,
         "However, the subspace of evecs is accurate to the required tolerance.\n");
   }

   switch (primme.dynamicMethodSwitch) {
      case -1: fprintf(primme.outputFile,
            "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
      case -2: fprintf(primme.outputFile,
            "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
      case -3: fprintf(primme.outputFile,
            "Recommended method for next run: DYNAMIC (close call)\n"); break;
   }


   /* Perturb the 5 approximate eigenvectors in evecs and used them as initial solution.
      This time the solver should converge faster than the last one. */
   for (i=0; i<primme.n*5; i++)
      evecs[i] += rand()/(double)RAND_MAX*1e-4;
   primme.initSize = 5;
   primme.numEvals = 5;

   /* Call primme  */
   ret = zprimme_normal(evals, evecs, rnorms, &primme);

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   /* Reporting (optional) */
   for (i=0; i < primme.initSize; i++) {
      fprintf(primme.outputFile, "Eval[%d]: %-22.15E + %-22.15E I rnorm: %-22.15E\n", i+1,
         creal(evals[i]), cimag(evals[i]), rnorms[i]); 
   }
   fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
   fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme.aNorm*primme.eps);
   fprintf(primme.outputFile, "Iterations: %-" PRIMME_INT_P "\n", 
                                                 primme.stats.numOuterIterations); 
   fprintf(primme.outputFile, "Restarts  : %-" PRIMME_INT_P "\n", primme.stats.numRestarts);
   fprintf(primme.outputFile, "Matvecs   : %-" PRIMME_INT_P "\n", primme.stats.numMatvecs);
   fprintf(primme.outputFile, "Preconds  : %-" PRIMME_INT_P "\n", primme.stats.numPreconds);
   if (primme.stats.lockingIssue) {
      fprintf(primme.outputFile, "\nA locking problem has occurred.\n");
      fprintf(primme.outputFile,
         "Some eigenpairs do not have a residual norm less than the tolerance.\n");
      fprintf(primme.outputFile,
         "However, the subspace of evecs is accurate to the required tolerance.\n");
   }

   switch (primme.dynamicMethodSwitch) {
      case -1: fprintf(primme.outputFile,
            "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
      case -2: fprintf(primme.outputFile,
            "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
      case -3: fprintf(primme.outputFile,
            "Recommended method for next run: DYNAMIC (close call)\n"); break;
   }


   /* Find the next 5 eigenpairs closest to .5 */
   primme.initSize = 0;
   primme.numEvals = 5;
   primme.numOrthoConst = 5; /* solver will find solutions orthogonal to the already
                                5 approximate eigenvectors in evecs */

   /* Call primme  */
   ret = zprimme_normal(evals, evecs, rnorms, &primme);

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   /* Reporting (optional) */
   for (i=0; i < primme.initSize; i++) {
      fprintf(primme.outputFile, "Eval[%d]: %-22.15E + %-22.15E I rnorm: %-22.15E\n", i+1,
         creal(evals[i]), cimag(evals[i]), rnorms[i]); 
   }
   fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
   fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme.aNorm*primme.eps);
   fprintf(primme.outputFile, "Iterations: %-" PRIMME_INT_P "\n", 
                                                 primme.stats.numOuterIterations); 
   fprintf(primme.outputFile, "Restarts  : %-" PRIMME_INT_P "\n", primme.stats.numRestarts);
   fprintf(primme.outputFile, "Matvecs   : %-" PRIMME_INT_P "\n", primme.stats.numMatvecs);
   fprintf(primme.outputFile, "Preconds  : %-" PRIMME_INT_P "\n", primme.stats.numPreconds);
   if (primme.stats.lockingIssue) {
      fprintf(primme.outputFile, "\nA locking problem has occurred.\n");
      fprintf(primme.outputFile,
         "Some eigenpairs do not have a residual norm less than the tolerance.\n");
      fprintf(primme.outputFile,
         "However, the subspace of evecs is accurate to the required tolerance.\n");
   }

   switch (primme.dynamicMethodSwitch) {
      case -1: fprintf(primme.outputFile,
            "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
      case -2: fprintf(primme.outputFile,
            "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
      case -3: fprintf(primme.outputFile,
            "Recommended method for next run: DYNAMIC (close call)\n"); break;
   }

   primme_free(&primme);
   free(evals);
   free(evecs);
   free(rnorms);

  return(0);
}

/* Modified 1-D Laplacian block matrix-vector product, Y = A * X, where

   - X, input dense matrix of size primme.n x blockSize;
   - Y, output dense matrix of size primme.n x blockSize;
   - A, tridiagonal square matrix of dimension primme.n with this form:

        [ 2   -1+i  0    0    0 ... ]
        [-1+i  2   -1+i  0    0 ... ]
        [ 0   -1+i  2   -1+i  0 ... ]
         ...
*/

void LaplacianLikeMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
   cdouble *xvec;     /* pointer to i-th input vector x */
   cdouble *yvec;     /* pointer to i-th output vector y */
   
   for (i=0; i<*blockSize; i++) {
      xvec = (cdouble *)x + *ldx*i;
      yvec = (cdouble *)y + *ldy*i;
      for (row=0; row<primme->n; row++) {
         yvec[row] = 0.0;
         if (row-1 >= 0) yvec[row] += (-1.0+I)*xvec[row-1];
         yvec[row] += 2.0*xvec[row];
         if (row+1 < primme->n) yvec[row] += (-1.0+I)*xvec[row+1];
      }      
   }
   *err = 0;
}

/* This performs Y = M^{-1} * X, where

   - X, input dense matrix of size primme.n x blockSize;
   - Y, output dense matrix of size primme.n x blockSize;
   - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
*/

void LaplacianLikeApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
   cdouble *xvec;     /* pointer to i-th input vector x */
   cdouble *yvec;     /* pointer to i-th output vector y */
    
   for (i=0; i<*blockSize; i++) {
      xvec = (cdouble *)x + *ldx*i;
      yvec = (cdouble *)y + *ldy*i;
      for (row=0; row<primme->n; row++) {
         yvec[row] = xvec[row]/2.;
      }      
   }
   *ierr = 0;
}
