/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2015 College of William & Mary,
 *   James R. McCombs, Eloy Romero Alcalde, Andreas Stathopoulos, Lingfei Wu
 *
 *   This file is part of PRIMME.
 *
 *   PRIMME is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   PRIMME is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *******************************************************************************
 *
 *  Example to compute the k largest eigenvalues in a 1-D Laplacian matrix.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "primme.h"   /* header file is required to run primme */ 

void LaplacianMatrixMatvec(void *x, void *y, int *blockSize, primme_params *primme);
void LaplacianApplyPreconditioner(void *x, void *y, int *blockSize, primme_params *primme);

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
   double targetShifts[1];

   /* Other miscellaneous items */
   int ret;
   int i;

   /* Set default values in PRIMME configuration struct */
   primme_initialize(&primme);

   /* Set problem matrix */
   primme.matrixMatvec = LaplacianMatrixMatvec;
                           /* Function that implements the matrix-vector product
                              A*x for solving the problem A*x = l*x */
  
   /* Set problem parameters */
   primme.n = 100; /* set problem dimension */
   primme.numEvals = 10;   /* Number of wanted eigenpairs */
   primme.eps = 1e-9;      /* ||r|| <= eps * ||matrix|| */
   primme.target = primme_smallest;
                           /* Wanted the smallest eigenvalues */

   /* Set preconditioner (optional) */
   primme.applyPreconditioner = LaplacianApplyPreconditioner;
   primme.correctionParams.precondition = 1;

   /* Set advanced parameters if you know what are you doing (optional) */
   /*
   primme.maxBasisSize = 14;
   primme.minRestartSize = 4;
   primme.maxBlockSize = 1;
   primme.maxMatvecs = 1000;
   */

   /* Set method to solve the problem */
   primme_set_method(DYNAMIC, &primme);
   /* DYNAMIC uses a runtime heuristic to choose the fastest method between
       DEFAULT_MIN_TIME and DEFAULT_MIN_MATVECS. But you can set another
       method, such as LOBPCG_OrthoBasis_Window, directly */

   /* Display PRIMME configuration struct (optional) */
   primme_display_params(primme);

   /* Allocate space for converged Ritz values and residual norms */
   evals = (double *)primme_calloc(primme.numEvals, sizeof(double), "evals");
   evecs = (complex double *)primme_calloc(primme.n*primme.numEvals, 
                                sizeof(complex double), "evecs");
   rnorms = (double *)primme_calloc(primme.numEvals, sizeof(double), "rnorms");

   /* Call primme  */
   ret = zprimme(evals, (Complex_Z*)evecs, rnorms, &primme);

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   /* Reporting (optional) */
   primme_PrintStackTrace(primme);

   for (i=0; i < primme.initSize; i++) {
      fprintf(primme.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
         evals[i], rnorms[i]); 
   }
   fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
   fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme.aNorm*primme.eps);
   fprintf(primme.outputFile, "Iterations: %-d\n", 
                                                 primme.stats.numOuterIterations); 
   fprintf(primme.outputFile, "Restarts  : %-d\n", primme.stats.numRestarts);
   fprintf(primme.outputFile, "Matvecs   : %-d\n", primme.stats.numMatvecs);
   fprintf(primme.outputFile, "Preconds  : %-d\n", primme.stats.numPreconds);
   if (primme.locking && primme.intWork && primme.intWork[0] == 1) {
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


   /* Note that d/zprimme can be called more than once before call primme_Free. */
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
   ret = zprimme(evals, (Complex_Z*)evecs, rnorms, &primme);

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   /* Reporting (optional) */
   primme_PrintStackTrace(primme);

   for (i=0; i < primme.initSize; i++) {
      fprintf(primme.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
         evals[i], rnorms[i]); 
   }
   fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
   fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme.aNorm*primme.eps);
   fprintf(primme.outputFile, "Iterations: %-d\n", 
                                                 primme.stats.numOuterIterations); 
   fprintf(primme.outputFile, "Restarts  : %-d\n", primme.stats.numRestarts);
   fprintf(primme.outputFile, "Matvecs   : %-d\n", primme.stats.numMatvecs);
   fprintf(primme.outputFile, "Preconds  : %-d\n", primme.stats.numPreconds);
   if (primme.locking && primme.intWork && primme.intWork[0] == 1) {
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
   ret = zprimme(evals, (Complex_Z*)evecs, rnorms, &primme);

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   /* Reporting (optional) */
   primme_PrintStackTrace(primme);

   for (i=0; i < primme.initSize; i++) {
      fprintf(primme.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
         evals[i], rnorms[i]); 
   }
   fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
   fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme.aNorm*primme.eps);
   fprintf(primme.outputFile, "Iterations: %-d\n", 
                                                 primme.stats.numOuterIterations); 
   fprintf(primme.outputFile, "Restarts  : %-d\n", primme.stats.numRestarts);
   fprintf(primme.outputFile, "Matvecs   : %-d\n", primme.stats.numMatvecs);
   fprintf(primme.outputFile, "Preconds  : %-d\n", primme.stats.numPreconds);
   if (primme.locking && primme.intWork && primme.intWork[0] == 1) {
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
   ret = zprimme(evals, (Complex_Z*)evecs, rnorms, &primme);

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   /* Reporting (optional) */
   primme_PrintStackTrace(primme);

   for (i=0; i < primme.initSize; i++) {
      fprintf(primme.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
         evals[i], rnorms[i]); 
   }
   fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
   fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
                                                         primme.aNorm*primme.eps);
   fprintf(primme.outputFile, "Iterations: %-d\n", 
                                                 primme.stats.numOuterIterations); 
   fprintf(primme.outputFile, "Restarts  : %-d\n", primme.stats.numRestarts);
   fprintf(primme.outputFile, "Matvecs   : %-d\n", primme.stats.numMatvecs);
   fprintf(primme.outputFile, "Preconds  : %-d\n", primme.stats.numPreconds);
   if (primme.locking && primme.intWork && primme.intWork[0] == 1) {
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

   primme_Free(&primme);
   free(evals);
   free(evecs);
   free(rnorms);

  return(0);
}

/* 1-D Laplacian block matrix-vector product, Y = A * X, where

   - X, input dense matrix of size primme.n x blockSize;
   - Y, output dense matrix of size primme.n x blockSize;
   - A, tridiagonal square matrix of dimension primme.n with this form:

        [ 2 -1  0  0  0 ... ]
        [-1  2 -1  0  0 ... ]
        [ 0 -1  2 -1  0 ... ]
         ...
*/

void LaplacianMatrixMatvec(void *x, void *y, int *blockSize, primme_params *primme) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
   complex double *xvec;     /* pointer to i-th input vector x */
   complex double *yvec;     /* pointer to i-th output vector y */
   
   for (i=0; i<*blockSize; i++) {
      xvec = (complex double *)x + primme->n*i;
      yvec = (complex double *)y + primme->n*i;
      for (row=0; row<primme->n; row++) {
         yvec[row] = 0.0;
         if (row-1 >= 0) yvec[row] += -1.0*xvec[row-1];
         yvec[row] += 2.0*xvec[row];
         if (row+1 < primme->n) yvec[row] += -1.0*xvec[row+1];
      }      
   }
}

/* This performs Y = M^{-1} * X, where

   - X, input dense matrix of size primme.n x blockSize;
   - Y, output dense matrix of size primme.n x blockSize;
   - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
*/

void LaplacianApplyPreconditioner(void *x, void *y, int *blockSize, primme_params *primme) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
   complex double *xvec;     /* pointer to i-th input vector x */
   complex double *yvec;     /* pointer to i-th output vector y */
    
   for (i=0; i<*blockSize; i++) {
      xvec = (complex double *)x + primme->n*i;
      yvec = (complex double *)y + primme->n*i;
      for (row=0; row<primme->n; row++) {
         yvec[row] = xvec[row]/2.;
      }      
   }
}
