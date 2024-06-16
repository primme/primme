#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "primme.h"   /* header file is required to run primme */ 

void LaplacianMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
void LaplacianApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
void convtest_sm(double *eval, void *evec, double *resNorm, int *isconv, primme_params *primme, int *ierr);


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
   primme.printLevel = 3; /*  report every iteration  */
   primme.numEvals = 10;   /* Number of wanted eigenpairs */
   primme.convTestFun = convtest_sm;
   primme.target = primme_largest;
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
   primme_set_method(PRIMME_DYNAMIC, &primme);
   /* DYNAMIC uses a runtime heuristic to choose the fastest method between
       PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
       set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */

   /* Display PRIMME configuration struct (optional) */
   primme_display_params(primme);

   /* Allocate space for converged Ritz values and residual norms */
   evals = (double*)malloc(primme.numEvals*sizeof(double));
   evecs = (double*)malloc(primme.n*primme.numEvals*sizeof(double));
   rnorms = (double*)malloc(primme.numEvals*sizeof(double));

   /* Call primme  */
   ret = dprimme(evals, evecs, rnorms, &primme);

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

void LaplacianMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
   double *xvec;     /* pointer to i-th input vector x */
   double *yvec;     /* pointer to i-th output vector y */
   
   for (i=0; i<*blockSize; i++) {
      xvec = (double *)x + *ldx*i;
      yvec = (double *)y + *ldy*i;
      for (row=0; row<primme->n; row++) {
         yvec[row] = 0.0;
         if (row-1 >= 0) yvec[row] += -1.0*xvec[row-1];
         yvec[row] += 2.0*xvec[row];
         if (row+1 < primme->n) yvec[row] += -1.0*xvec[row+1];
      }      
   }
   *err = 0;
}

/* This performs Y = M^{-1} * X, where

   - X, input dense matrix of size primme.n x blockSize;
   - Y, output dense matrix of size primme.n x blockSize;
   - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
*/

void convtest_sm(double *eval, void *evec, double *resNorm, int *isconv, primme_params *primme, int *ierr){
	*isconv = (resNorm && eval && (abs(*eval) > 0.1 * (*resNorm)));
	*ierr = 0;
}

void LaplacianApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
   double *xvec;     /* pointer to i-th input vector x */
   double *yvec;     /* pointer to i-th output vector y */
    
   for (i=0; i<*blockSize; i++) {
      xvec = (double *)x + *ldx*i;
      yvec = (double *)y + *ldy*i;
      for (row=0; row<primme->n; row++) {
         yvec[row] = xvec[row]/2.;
      }      
   }
   *ierr = 0;
}
