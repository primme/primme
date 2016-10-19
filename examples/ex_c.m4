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
 *  Example to compute the k largest eigenvalues in a 1-D Laplacian matrix.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
ifdef(`USE_COMPLEX', ifdef(`USE_COMPLEX_CXX', ``#include <complex>'', ``#include <complex.h>''))
ifdef(`USE_PETSC', ``#include <petscpc.h>
#include <petscmat.h>
'')dnl
#include "primme.h"   /* header file is required to run primme */ 
define(`PRIMME_NUM', ifdef(`USE_PETSC', `PetscScalar', ifdef(`USE_COMPLEX', ifdef(`USE_COMPLEX_CXX', `std::complex<double>', `complex double'), `double')))dnl
ifdef(`USE_PETSC', `
PetscErrorCode generateLaplacian1D(int n, Mat *A);
void PETScMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
void ApplyPCPrecPETSC(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
static void par_GlobalSum(void *sendBuf, void *recvBuf, int *count,
                         primme_params *primme, int *ierr);
', `
void LaplacianMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
void LaplacianApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
')dnl

int main (int argc, char *argv[]) {

   /* Solver arrays and parameters */
   double *evals;    /* Array with the computed eigenvalues */
   double *rnorms;   /* Array with the computed eigenpairs residual norms */
   PRIMME_NUM *evecs;    /* Array with the computed eigenvectors;
                        first vector starts in evecs[0],
                        second vector starts in evecs[primme.n],
                        third vector starts in evecs[primme.n*2]...  */
   primme_params primme;
                     /* PRIMME configuration struct */
ifdef(`ADVANCED', `   double targetShifts[1];
')dnl

   /* Other miscellaneous items */
   int ret;
   int i;
ifdef(`USE_PETSC', `   Mat A; /* problem matrix */
   PC pc;            /* preconditioner */
   PetscErrorCode ierr;
   PetscInt n, nLocal;
   MPI_Comm comm;

   PetscInitialize(&argc, &argv, NULL, NULL);

')dnl

   /* Set default values in PRIMME configuration struct */
   primme_initialize(&primme);

   /* Set problem matrix */
ifdef(`USE_PETSC', `   ierr = generateLaplacian1D(100, &A); CHKERRQ(ierr);
   primme.matrix = &A;
   primme.matrixMatvec = PETScMatvec;
', `   primme.matrixMatvec = LaplacianMatrixMatvec;
')dnl
                           /* Function that implements the matrix-vector product
                              A*x for solving the problem A*x = l*x */
  
   /* Set problem parameters */
ifdef(`USE_PETSC', `   ierr = MatGetSize(A, &n, NULL); CHKERRQ(ierr);
   primme.n = (PRIMME_INT)n;',
`   primme.n = 100;') /* set problem dimension */
   primme.numEvals = 10;   /* Number of wanted eigenpairs */
   primme.eps = 1e-9;      /* ||r|| <= eps * ||matrix|| */
   primme.target = primme_smallest;
                           /* Wanted the smallest eigenvalues */

   /* Set preconditioner (optional) */
ifdef(`USE_PETSC', `   ierr = PCCreate(PETSC_COMM_WORLD, &pc); CHKERRQ(ierr);
   ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);
   ierr = PCSetOperators(pc, A, A); CHKERRQ(ierr);
   ierr = PCSetFromOptions(pc); CHKERRQ(ierr);
   ierr = PCSetUp(pc); CHKERRQ(ierr);
   primme.preconditioner = &pc;
   primme.applyPreconditioner = ApplyPCPrecPETSC;
', `   primme.applyPreconditioner = LaplacianApplyPreconditioner;
')dnl
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

ifdef(`USE_PETSC', `   /* Set parallel parameters */
   ierr = MatGetLocalSize(A, &nLocal, NULL); CHKERRQ(ierr);
   primme.nLocal = (PRIMME_INT)nLocal;
   comm = PETSC_COMM_WORLD;
   primme.commInfo = &comm;
   MPI_Comm_size(comm, &primme.numProcs);
   MPI_Comm_rank(comm, &primme.procID);
   primme.globalSumReal = par_GlobalSum;

')dnl
   /* Display PRIMME configuration struct (optional) */
ifdef(`USE_PETSC', `   if (primme.procID == 0) /* Reports process with ID 0 */
   ')   primme_display_params(primme);

   /* Allocate space for converged Ritz values and residual norms */ifdef(`USE_COMPLEX_CXX', `
   evals = new double[primme.numEvals];
   evecs = new PRIMME_NUM[primme.n*primme.numEvals];
   rnorms = new double[primme.numEvals];',`
   evals = (double*)malloc(primme.numEvals*sizeof(double));
   evecs = (PRIMME_NUM*)malloc(primme.n*primme.numEvals*sizeof(PRIMME_NUM));
   rnorms = (double*)malloc(primme.numEvals*sizeof(double));')

define(`CALL_PRIMME', `   /* Call primme  */
ifdef(`USE_PETSC', ``#if defined(PETSC_USE_COMPLEX)
   ret = zprimme(evals, evecs, rnorms, &primme);
#else
   ret = dprimme(evals, evecs, rnorms, &primme);
#endif
'',
`   ret = ifdef(`USE_COMPLEX',`z', `d')primme(evals, evecs, rnorms, &primme);
')dnl

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

ifdef(`USE_PETSC', ``   if (primme.procID == 0) { /* Reports process with ID 0 */
' define(sp, `   ')', `define(sp, `')')dnl
   sp()/* Reporting (optional) */
   sp()for (i=0; i < primme.initSize; i++) {
   sp()   fprintf(primme.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
   sp()      evals[i], rnorms[i]); 
   sp()}
   sp()fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
   sp()fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
   sp()                                                      primme.aNorm*primme.eps);
   sp()fprintf(primme.outputFile, "Iterations: %-" PRIMME_INT_P "\n", 
   sp()                                              primme.stats.numOuterIterations); 
   sp()fprintf(primme.outputFile, "Restarts  : %-" PRIMME_INT_P "\n", primme.stats.numRestarts);
   sp()fprintf(primme.outputFile, "Matvecs   : %-" PRIMME_INT_P "\n", primme.stats.numMatvecs);
   sp()fprintf(primme.outputFile, "Preconds  : %-" PRIMME_INT_P "\n", primme.stats.numPreconds);
   sp()if (primme.locking && primme.intWork && primme.intWork[0] == 1) {
   sp()   fprintf(primme.outputFile, "\nA locking problem has occurred.\n");
   sp()   fprintf(primme.outputFile,
   sp()      "Some eigenpairs do not have a residual norm less than the tolerance.\n");
   sp()   fprintf(primme.outputFile,
   sp()      "However, the subspace of evecs is accurate to the required tolerance.\n");
   sp()}

   sp()switch (primme.dynamicMethodSwitch) {
   sp()   case -1: fprintf(primme.outputFile,
   sp()         "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
   sp()   case -2: fprintf(primme.outputFile,
   sp()         "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
   sp()   case -3: fprintf(primme.outputFile,
   sp()         "Recommended method for next run: DYNAMIC (close call)\n"); break;
   sp()}
ifdef(`USE_PETSC', `   }
')dnl
')dnl end of CALL_PRIMME
CALL_PRIMME
ifdef(`ADVANCED', `
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

CALL_PRIMME

   /* Perturb the 5 approximate eigenvectors in evecs and used them as initial solution.
      This time the solver should converge faster than the last one. */
   for (i=0; i<primme.n*5; i++)
      evecs[i] += rand()/(double)RAND_MAX*1e-4;
   primme.initSize = 5;
   primme.numEvals = 5;

CALL_PRIMME

   /* Find the next 5 eigenpairs closest to .5 */
   primme.initSize = 0;
   primme.numEvals = 5;
   primme.numOrthoConst = 5; /* solver will find solutions orthogonal to the already
                                5 approximate eigenvectors in evecs */

CALL_PRIMME
')dnl
   primme_free(&primme);ifdef(`USE_COMPLEX_CXX', `
   delete [] evals;
   delete [] evecs;
   delete [] rnorms;',`
   free(evals);
   free(evecs);
   free(rnorms);')

ifdef(`USE_PETSC', `   ierr = PetscFinalize(); CHKERRQ(ierr);

')dnl
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
ifdef(`USE_PETSC', `
PetscErrorCode generateLaplacian1D(int n, Mat *A) {
   PetscScalar    value[3] = {-1.0, 2.0, -1.0};
   PetscInt       i,Istart,Iend,col[3];
   PetscBool      FirstBlock=PETSC_FALSE,LastBlock=PETSC_FALSE;
   PetscErrorCode ierr;

   PetscFunctionBegin;

   ierr = MatCreate(PETSC_COMM_WORLD, A); CHKERRQ(ierr);
   ierr = MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRQ(ierr);
   ierr = MatSetFromOptions(*A); CHKERRQ(ierr);
   ierr = MatSetUp(*A); CHKERRQ(ierr);

   ierr = MatGetOwnershipRange(*A, &Istart, &Iend); CHKERRQ(ierr);
   if (Istart == 0) FirstBlock = PETSC_TRUE;
   if (Iend == n) LastBlock = PETSC_TRUE;
   for (i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++) {
      col[0]=i-1; col[1]=i; col[2]=i+1;
      ierr = MatSetValues(*A, 1, &i, 3, col, value, INSERT_VALUES); CHKERRQ(ierr);
   }
   if (LastBlock) {
      i=n-1; col[0]=n-2; col[1]=n-1;
      ierr = MatSetValues(*A, 1, &i, 2, col, value, INSERT_VALUES); CHKERRQ(ierr);
   }
   if (FirstBlock) {
      i=0; col[0]=0; col[1]=1; value[0]=2.0; value[1]=-1.0;
      ierr = MatSetValues(*A, 1, &i, 2, col, value, INSERT_VALUES); CHKERRQ(ierr);
   }

   ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
   ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

   PetscFunctionReturn(0);
}

void PETScMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   int i;
   Mat *matrix;
   Vec xvec, yvec;
   PetscErrorCode ierr;

   matrix = (Mat *)primme->matrix;

   ierr = MatCreateVecs(*matrix, &xvec, &yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   for (i=0; i<*blockSize; i++) {
      ierr = VecPlaceArray(xvec, ((PetscScalar*)x) + *ldx*i); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecPlaceArray(yvec, ((PetscScalar*)y) + *ldy*i); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = MatMult(*matrix, xvec, yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecResetArray(xvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecResetArray(yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   }
   ierr = VecDestroy(&xvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   ierr = VecDestroy(&yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   *err = 0; 
}
', `
void LaplacianMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
   PRIMME_NUM *xvec;     /* pointer to i-th input vector x */
   PRIMME_NUM *yvec;     /* pointer to i-th output vector y */
   
   for (i=0; i<*blockSize; i++) {
      xvec = (PRIMME_NUM *)x + *ldx*i;
      yvec = (PRIMME_NUM *)y + *ldy*i;
      for (row=0; row<primme->n; row++) {
         yvec[row] = 0.0;
         if (row-1 >= 0) yvec[row] += -1.0*xvec[row-1];
         yvec[row] += 2.0*xvec[row];
         if (row+1 < primme->n) yvec[row] += -1.0*xvec[row+1];
      }      
   }
   *err = 0;
}
')dnl

/* This performs Y = M^{-1} * X, where

   - X, input dense matrix of size primme.n x blockSize;
   - Y, output dense matrix of size primme.n x blockSize;
   - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
*/
ifdef(`USE_PETSC', `
void ApplyPCPrecPETSC(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   int i;
   Mat *matrix;
   PC *pc;
   Vec xvec, yvec;
   PetscErrorCode ierr;
   
   matrix = (Mat *)primme->matrix;
   pc = (PC *)primme->preconditioner;

   ierr = MatCreateVecs(*matrix, &xvec, &yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   for (i=0; i<*blockSize; i++) {
      ierr = VecPlaceArray(xvec, ((PetscScalar*)x) + *ldx*i); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecPlaceArray(yvec, ((PetscScalar*)y) + *ldy*i); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = PCApply(*pc, xvec, yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecResetArray(xvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecResetArray(yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   }
   ierr = VecDestroy(&xvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   ierr = VecDestroy(&yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
}

static void par_GlobalSum(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme, int *ierr) {
   MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;

   if (sendBuf == recvBuf) {
     *ierr = MPI_Allreduce(MPI_IN_PLACE, recvBuf, *count, MPIU_REAL, MPI_SUM, communicator) != MPI_SUCCESS;
   } else {
     *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPIU_REAL, MPI_SUM, communicator) != MPI_SUCCESS;
   }
}
', `
void LaplacianApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
   PRIMME_NUM *xvec;     /* pointer to i-th input vector x */
   PRIMME_NUM *yvec;     /* pointer to i-th output vector y */
    
   for (i=0; i<*blockSize; i++) {
      xvec = (PRIMME_NUM *)x + *ldx*i;
      yvec = (PRIMME_NUM *)y + *ldy*i;
      for (row=0; row<primme->n; row++) {
         yvec[row] = xvec[row]/2.;
      }      
   }
   *ierr = 0;
}
')dnl
