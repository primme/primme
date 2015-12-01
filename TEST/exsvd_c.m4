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
 *  Example to compute the k largest singular values in Lauchli matrix.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
ifdef(`USE_COMPLEX', ``#include <complex.h>
'')dnl
ifdef(`USE_PETSC', ``#include <petscpc.h>
#include <petscmat.h>
'')dnl
//#include "primme.h"   /* header file for PRIMME SVDS too */ 
#include "primme_svds.h"   /* this shouldn't be in this way */ 
define(`PRIMME_NUM', ifdef(`USE_PETSC', `PetscScalar', ifdef(`USE_COMPLEX', `complex double', `double')))dnl
ifdef(`USE_PETSC', `
PetscErrorCode generateLauchli(int n, double mu, Mat *A);
void PETScMatvec(void *x, void *y, int *blockSize, primme_params *primme);
void ApplyPCPrecPETSC(void *x, void *y, int *blockSize, primme_params *primme);
static void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count,
                         primme_params *primme);
', `
void LauchliMatrixMatvec(void *x, int *ldx, void *y, int *ldy, int *blockSize,
                         int *transpose, primme_svds_params *primme_svds);
void LauchliApplyPreconditioner(void *x, int *ldx, void *y, int *ldy, int *blockSize,
                                int *mode, primme_svds_params *primme_svds);
')dnl

int main (int argc, char *argv[]) {

   /* Solver arrays and parameters */
   double *svals;    /* Array with the computed singular values */
   double *rnorms;   /* Array with the computed residual norms */
   PRIMME_NUM *svecs;    /* Array with the computed singular vectors;
                        first right (v) vector starts in svecs[0],
                        second right (v) vector starts in svecs[primme_svd.n],
                        first left (u) vector starts in svecs[primme_svd.n*numSVals]...  */
   primme_svds_params primme_svds;
   primme_params *primme;
                     /* PRIMME SVDS configuration struct */
ifdef(`ADVANCED', `   double targetShifts[1];
')dnl

   /* Other miscellaneous items */
   int ret;
   int i;
   double mu = 1.234567890123456e-7;
ifdef(`USE_PETSC', `   Mat A; /* problem matrix */
   PC pc;            /* preconditioner */
   PetscErrorCode ierr;
   MPI_Comm comm;

   PetscInitialize(&argc, &argv, NULL, NULL);

')dnl

   /* Set default values in PRIMME SVDS configuration struct */
   primme_svds_initialize(&primme_svds);

   /* Set problem matrix */
ifdef(`USE_PETSC', `   ierr = generateLauchli(100, mu, &A); CHKERRQ(ierr);
   primme_svds.matrix = &A;
   primme_svds.matrixMatvec = PETScMatvec;
', `   primme_svds.matrixMatvec = LauchliMatrixMatvec;
   primme_svds.matrix = &mu;
')dnl
                           /* Function that implements the matrix-vector products
                              A*x and A^t*x  */
  
   /* Set problem parameters */
ifdef(`USE_PETSC', `   ierr = MatGetSize(A, &primme_svds.m, &primme_svds.n); CHKERRQ(ierr);',
                   `   primme_svds.m = 100+1;
   primme_svds.n = 100;') /* set problem dimension */
   primme_svds.numSvals = 4;   /* Number of wanted eigenpairs */
   primme_svds.eps = 1e-12;     /* ||r|| <= eps * ||matrix|| */
   primme_svds.target = primme_svds_largest;
                               /* Seeking for the largest singular values  */

   /* Set preconditioner (optional) */
ifdef(`USE_PETSC', `   ierr = PCCreate(PETSC_COMM_WORLD, &pc); CHKERRQ(ierr);
   ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);
   ierr = PCSetOperators(pc, A, A); CHKERRQ(ierr);
   ierr = PCSetFromOptions(pc); CHKERRQ(ierr);
   ierr = PCSetUp(pc); CHKERRQ(ierr);
   primme.preconditioner = &pc;
   primme.applyPreconditioner = ApplyPCPrecPETSC;
', `//   primme_svds.applyPreconditioner = LauchliApplyPreconditioner;
')dnl

   /* Set method to solve the singular value problem and
      the underneath eigenvalue problem (optional) */
ifdef(`ADVANCED_HYBRID', `   primme_svds_set_method(primme_svds_hybrid, DYNAMIC, DEFAULT_MIN_TIME, &primme_svds);
   /*  Set hybrid method with DYNAMIC and DEFAULT_MIN_TIME as the underneath eigensolver configuration
       for the first and the second stage, respectively.
       DYNAMIC uses a runtime heuristic to choose the fastest method between
       DEFAULT_MIN_TIME and DEFAULT_MIN_MATVECS. But you can set another
       method, such as LOBPCG_OrthoBasis_Window, directly */
', `/*   primme_svds_set_method(primme_svds_default, &primme_svds); */
   /*  primme_svds_default: devs choice, now being hybrid, which first solve
       the normal equation and then the augmented problem.
       DYNAMIC uses a runtime heuristic to choose the fastest method between
       DEFAULT_MIN_TIME and DEFAULT_MIN_MATVECS. But you can set another
       method, such as LOBPCG_OrthoBasis_Window, directly */
')dnl
   primme_svds_set_method(primme_svds_default, &primme_svds);
   primme_set_method(GD_Olsen_plusK, &primme_svds.primme);
   primme_svds.printLevel = 3;

ifdef(`USE_PETSC', `   /* Set parallel parameters */
   ierr = MatGetLocalSize(A, &primme.nLocal, NULL); CHKERRQ(ierr);
   comm = PETSC_COMM_WORLD;
   primme.commInfo = &comm;
   MPI_Comm_size(comm, &primme.numProcs);
   MPI_Comm_rank(comm, &primme.procID);
   primme.globalSumDouble = par_GlobalSumDouble;

')dnl
ifdef(`ADVANCED_HYBRID', `   /* Set parameters for the underneath eigensolver if you know what are you doing (optional) */
   primme_svds.primme->locking = 1; 
   primme_svds.primme->restart.numPrevRetain = 3;
', `
   /* Set advanced parameters if you know what are you doing (optional) */
   /* Configuration for 1st stage */
   /*
   primme = primme_svds_get_primme_params(primme_svds);
   or you can do: primme = primme_svds.primme;
   primme->maxBasisSize = 14;
   primme->minRestartSize = 6;
   primme->maxBlockSize = 2;
   */
   /* Configuration for 2nd stage */
   /*
   primme = primme_get_stage_conf(primme_svds, 1);
   or you can do: primme = primme_svds.primme->next_conf;
   primme->maxBasisSize = 30;
   primme->minRestartSize = 15;
   primme->maxBlockSize = 1;
   */
')dnl

    /* Display PRIMME SVDS configuration struct (optional) */
ifdef(`USE_PETSC', `   if (primme.procID == 0) /* Reports process with ID 0 */
   ')   primme_svds_display_params(primme_svds);

   /* Allocate space for converged Ritz values and residual norms */
   svals = (double *)primme_calloc(primme_svds.numSvals, sizeof(double), "svals");
   svecs = (PRIMME_NUM *)primme_calloc((primme_svds.n+primme_svds.m)*primme_svds.numSvals, 
                                sizeof(PRIMME_NUM), "svecs");
   rnorms = (double *)primme_calloc(primme_svds.numSvals, sizeof(double), "rnorms");

define(`CALL_PRIMME_SVDS', `   /* Call primme_svds  */
ifdef(`USE_PETSC', ``#if defined(PETSC_USE_COMPLEX)
   ret = zprimme(evals, (Complex_Z*)evecs, rnorms, &primme);
#else
   ret = dprimme(evals, evecs, rnorms, &primme);
#endif
'',
`   ret = ifdef(`USE_COMPLEX',`z', `d')primme_svds(svals, ifdef(`USE_COMPLEX', `(Complex_Z*)')svecs, rnorms, &primme_svds);
')dnl

   if (ret != 0) {
      fprintf(primme_svds.outputFile, 
         "Error: primme_svds returned with nonzero exit status: %d \n",ret);
      return -1;
   }
')dnl end of CALL_PRIMME_SVDS

define(`REPORT_PRIMME_SVDS', `
ifdef(`USE_PETSC', ``   if (primme.procID == 0) { /* Reports process with ID 0 */
' define(sp, `   ')', `define(sp, `')')dnl
   sp()/* Reporting (optional) */
   sp()primme_PrintStackTrace(primme_svds.primme);

   sp()for (i=0; i < primme_svds.initSize; i++) {
   sp()   fprintf(primme_svds.outputFile, "Sval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
   sp()      svals[i], rnorms[i]); 
   sp()}
   sp()fprintf(primme_svds.outputFile, " %d eigenpairs converged\n", primme_svds.initSize);
   sp()fprintf(primme_svds.outputFile, "Tolerance : %-22.15E\n", 
   sp()                                                      primme_svds.aNorm*primme_svds.eps);
   sp()fprintf(primme_svds.outputFile, "Iterations: %-d\n", 
   sp()                                              primme_svds.stats.numOuterIterations); 
   sp()fprintf(primme_svds.outputFile, "Restarts  : %-d\n", primme_svds.stats.numRestarts);
   sp()fprintf(primme_svds.outputFile, "Matvecs   : %-d\n", primme_svds.stats.numMatvecs);
   sp()fprintf(primme_svds.outputFile, "Preconds  : %-d\n", primme_svds.stats.numPreconds);
   sp()if (primme_svds.primme.locking && primme_svds.primme.intWork && primme_svds.primme.intWork[0] == 1) {
   sp()   fprintf(primme_svds.outputFile, "\nA locking problem has occurred.\n"
   sp()      "Some eigenpairs do not have a residual norm less than the tolerance.\n"
   sp()      "However, the subspace of evecs is accurate to the required tolerance.\n");
   sp()}

ifdef(`USE_PETSC', `   }
')dnl
')dnl end of REPORT_PRIMME_SVDS
CALL_PRIMME_SVDS
ifdef(`ADVANCED_HYBRID', `
   /* Hybrid method second stage */
   primme_svds_set_method(primme_svds_augmented, DEFAULT_MIN_MATVECS, &primme_svds);
                              /* Set DEFAULT_MIN_MATVECS as the eigensolver for the
                                 second stage */
   primme_svds.primme.matrixMatvec = LauchliAugmentedMatvec; 
                              /* Set a custom matrix-vector product */
')dnl
REPORT_PRIMME_SVDS
ifdef(`ADVANCED', `
   /* Note that d/zprimme_svds can be called more than once before call primme_svds_Free. */

   /* Find the next 5 largest singular values */
   primme_svds.numEvals = 5;
   primme_svds.numOrthoConst = primme_svds.initSize;
                             /* the solver will find solutions orthogonal to the already
                                approximate singular vectors in svecs */
   primme_svds.initSize = 0;

CALL_PRIMME_SVDS
REPORT_PRIMME_SVDS
')dnl
   primme_svds_Free(&primme_svds);
   free(svals);
   free(svecs);
   free(rnorms);

ifdef(`USE_PETSC', `   ierr = PetscFinalize(); CHKERRQ(ierr);

')dnl
  return(0);
}

/* Lauchli block matrix-vector product, Y = A * X (or Y = A^t * X), where

   - X, input dense matrix of size primme_svds.n (or primme_svds.m) x blockSize;
   - Y, output dense matrix of size primme_svds.m (or primme_svds.n) x blockSize;
   - A, Lauchli matrix of dimensions primme_svds.m x (primme_svds.m+1) with this form:

        [ 1  1  1  1  1 ...  1 ]
        [mu  0  0  0  0 ...  0 ]
        [ 0 mu  0  0  0 ...  0 ]
         ...
        [ 0  0  0  0  0 ... mu ]
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

void PETScMatvec(void *x, void *y, int *blockSize, primme_params *primme) {
   int i;
   Mat *matrix;
   Vec xvec, yvec;
   PetscErrorCode ierr;

   matrix = (Mat *)primme->matrix;

   ierr = MatCreateVecs(*matrix, &xvec, &yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   for (i=0; i<*blockSize; i++) {
      ierr = VecPlaceArray(xvec, ((PetscScalar*)x) + primme->nLocal*i); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecPlaceArray(yvec, ((PetscScalar*)y) + primme->nLocal*i); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = MatMult(*matrix, xvec, yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecResetArray(xvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecResetArray(yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   }
   ierr = VecDestroy(&xvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   ierr = VecDestroy(&yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
}
', `
void LauchliMatrixMatvec(void *x, int *ldx, void *y, int *ldy, int *blockSize,
                         int *transpose, primme_svds_params *primme_svds) {
   
   int i;            /* vector index, from 0 to *blockSize-1 */
   int j;
   PRIMME_NUM *xvec;     /* pointer to i-th input vector x */
   PRIMME_NUM *yvec;     /* pointer to i-th output vector y */
   double mu = *(double*)primme_svds->matrix;

   if (*transpose == 0) { /* Do y <- A * x */
      for (i=0; i<*blockSize; i++) { 
         xvec = (PRIMME_NUM *)x + (*ldx)*i;
         yvec = (PRIMME_NUM *)y + (*ldy)*i;
         yvec[0] = 0;
         for (j=1; j<primme_svds->m; j++) {
            yvec[0] += xvec[j-1];
            yvec[j] = mu*xvec[j-1];
         }      
      }
   } else { /* Do y <- A^t * x */
      for (i=0; i<*blockSize; i++) {
         xvec = (PRIMME_NUM *)x + (*ldx)*i;
         yvec = (PRIMME_NUM *)y + (*ldy)*i;
         for (j=0; j<primme_svds->n; j++) {
            yvec[j] = xvec[0] + mu*xvec[j+1];
         }      
      }
   }
}
')dnl

/* This performs Y = M^{-1} * X, where

   - X, input dense matrix of size primme_svds.n (or primme_svds.m or m+n) x blockSize;
   - Y, output dense matrix of size primme_svds.n (or primme_svds.m or m+n) x blockSize;
   - M, preconditioner for A^t*A (or A*A^t or [0 A^t; A 0]), where A is the Lauchli matrix.
*/
ifdef(`USE_PETSC', `
void ApplyPCPrecPETSC(void *x, void *y, int *blockSize, primme_params *primme) {
   int i;
   Mat *matrix;
   PC *pc;
   Vec xvec, yvec;
   PetscErrorCode ierr;
   
   matrix = (Mat *)primme->matrix;
   pc = (PC *)primme->preconditioner;

   ierr = MatCreateVecs(*matrix, &xvec, &yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   for (i=0; i<*blockSize; i++) {
      ierr = VecPlaceArray(xvec, ((PetscScalar*)x) + primme->nLocal*i); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecPlaceArray(yvec, ((PetscScalar*)y) + primme->nLocal*i); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = PCApply(*pc, xvec, yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecResetArray(xvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
      ierr = VecResetArray(yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   }
   ierr = VecDestroy(&xvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   ierr = VecDestroy(&yvec); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
}

static void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme) {
   MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;

   MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator);
}
', `
void LauchliApplyPreconditioner(void *x, int *ldx, void *y, int *ldy, int *blockSize,
                                int *mode, primme_svds_params *primme_svds) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int j;            /* row index */
   PRIMME_NUM *xvec;     /* pointer to i-th input vector x */
   PRIMME_NUM *yvec;     /* pointer to i-th output vector y */
   int modeAtA = primme_svds_op_AtA, modeAAt = primme_svds_op_AAt;
   double mu = *(double*)primme_svds->matrix;
   PRIMME_NUM  *aux;
   int ldaux, notrans = 0, trans = 1;
    
   if (*mode == primme_svds_op_AtA) {
      /* Preconditioner for A^t*A, diag(A^t*A)^{-1} */
      for (i=0; i<*blockSize; i++) { 
         xvec = (PRIMME_NUM *)x + (*ldx)*i;
         yvec = (PRIMME_NUM *)y + (*ldy)*i;
         for (j=0; j<primme_svds->n; j++) {
            yvec[j] = xvec[j]/(1.0 + mu*mu);
         }      
      }
   }
   else if (*mode == primme_svds_op_AAt) {
      /* Preconditioner for A*A^t, diag(A*A^t)^{-1} */
      for (i=0; i<*blockSize; i++) {
         xvec = (PRIMME_NUM *)x + (*ldx)*i;
         yvec = (PRIMME_NUM *)y + (*ldy)*i;
         yvec[0] = xvec[0]/(PRIMME_NUM)primme_svds->m;
         for (j=1; j<primme_svds->m; j++) {
            yvec[j] = xvec[j]/mu/mu;
         }      
      }
   }
   else if (*mode == primme_svds_op_augmented) {
      /* Preconditioner for [0 A^t; A 0],
         [diag(A^t*A) 0; 0 diag(A*A^t)]^{-1}*[0 A^t; A 0] */

      /* [y0; y1] <- [0 A^t; A 0] * [x0; x1] */
      ldaux = primme_svds->n+primme_svds->m;
      aux = malloc(sizeof(PRIMME_NUM)*(*blockSize)*ldaux);
      primme_svds->matrixMatvec(x, ldx, &aux[primme_svds->n], &ldaux, blockSize, &notrans, primme_svds);
      xvec = (PRIMME_NUM *)x + primme_svds->n;
      primme_svds->matrixMatvec(xvec, ldx, aux, &ldaux, blockSize, &trans, primme_svds);
      /* y0 <- preconditioner for A^t*A  * y0 */
      LauchliApplyPreconditioner(aux, &ldaux, y, ldy, blockSize, &modeAtA, primme_svds);
      /* y1 <- preconditioner for A*A^t  * y1 */
      yvec = (PRIMME_NUM *)aux + primme_svds->n;
      LauchliApplyPreconditioner(yvec, &ldaux, yvec, ldy, blockSize, &modeAAt, primme_svds);
      free(aux);
   }
}
')dnl
