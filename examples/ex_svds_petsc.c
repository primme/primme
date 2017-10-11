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
 *  Example to compute the k largest singular values in Lauchli matrix.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <petscpc.h>
#include <petscmat.h>
#include "primme.h"   /* header file for PRIMME SVDS too */ 

#ifndef min
#define min(A,B) ((A)<=(B)?(A):(B))
#endif
#ifndef max
#define max(A,B) ((A)>=(B)?(A):(B))
#endif


PetscErrorCode generateLauchli(int m, int n, PetscReal mu, Mat *A);
void PETScMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                         int *transpose, primme_svds_params *primme_svds, int *ierr);
void ApplyPCPrecAHA(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                         int *transpose, primme_svds_params *primme_svds, int *ierr);
void par_GlobalSum(void *sendBuf, void *recvBuf, int *count,
                         primme_svds_params *primme_svds, int *ierr);

int main (int argc, char *argv[]) {

   /* Solver arrays and parameters */
   PetscReal *svals;    /* Array with the computed singular values */
   PetscReal *rnorms;   /* Array with the computed residual norms */
   PetscScalar *svecs;    /* Array with the computed singular vectors;
                        first right (v) vector starts in svecs[0],
                        second right (v) vector starts in svecs[primme_svd.n],
                        first left (u) vector starts in svecs[primme_svd.n*numSVals]...  */
   primme_svds_params primme_svds;
                     /* PRIMME SVDS configuration struct */

   /* Other miscellaneous items */
   int ret;
   int i;
   PetscReal mu = 1e-5;
   Mat A; /* problem matrix */
   Mat AHA;          /* auxiliary matrix for A^t*A */
   PC pc;            /* preconditioner */
   PetscErrorCode ierr;
   PetscInt m, n, mLocal, nLocal;
   MPI_Comm comm;

   PetscInitialize(&argc, &argv, NULL, NULL);


   /* Set default values in PRIMME SVDS configuration struct */
   primme_svds_initialize(&primme_svds);

   /* Set problem matrix */
   ierr = generateLauchli(500, 100, mu, &A); CHKERRQ(ierr);
   primme_svds.matrix = &A;
   primme_svds.matrixMatvec = PETScMatvec;
                           /* Function that implements the matrix-vector products
                              A*x and A^t*x  */
  
   /* Set problem parameters */
   ierr = MatGetSize(A, &m, &n); CHKERRQ(ierr);
   primme_svds.m = (PRIMME_INT)m;
   primme_svds.n = (PRIMME_INT)n; /* set problem dimension */
   primme_svds.numSvals = 4;   /* Number of wanted singular values */
   primme_svds.eps = 1e-6;     /* ||r|| <= eps * ||matrix|| */
   primme_svds.target = primme_svds_smallest;
                               /* Seeking for the largest singular values  */

   /* Set preconditioner (optional) */
   /* Build the Jacobi preconditioner of A^T*A, useful when m>=n */
   ierr = MatCreateNormal(A, &AHA); CHKERRQ(ierr);
   ierr = PCCreate(PETSC_COMM_WORLD, &pc); CHKERRQ(ierr);
   ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);
   ierr = PCSetOperators(pc, AHA, AHA); CHKERRQ(ierr);
   ierr = PCSetFromOptions(pc); CHKERRQ(ierr);
   ierr = PCSetUp(pc); CHKERRQ(ierr);
   primme_svds.preconditioner = &pc;
   primme_svds.applyPreconditioner = ApplyPCPrecAHA;

   /* Set method to solve the singular value problem and
      the underneath eigenvalue problem (optional) */
   primme_svds_set_method(primme_svds_default, PRIMME_DEFAULT_METHOD,
                              PRIMME_DEFAULT_METHOD, &primme_svds);
   /*  primme_svds_default: devs choice, now being hybrid, which first solve
       the normal equation and then the augmented problem.
       PRIMME_DEFAULT_METHOD devs choice of the solver at every stage. But other methods
       can be set such as DYNAMIC or PRIMME_LOBPCG_OrthoBasis_Window. */

   primme_svds.printLevel = 3;

   /* Set parallel parameters */
   ierr = MatGetLocalSize(A, &mLocal, &nLocal); CHKERRQ(ierr);
   primme_svds.mLocal = (PRIMME_INT)mLocal;
   primme_svds.nLocal = (PRIMME_INT)nLocal;
   comm = PETSC_COMM_WORLD;
   primme_svds.commInfo = &comm;
   MPI_Comm_size(comm, &primme_svds.numProcs);
   MPI_Comm_rank(comm, &primme_svds.procID);
   primme_svds.globalSumReal = par_GlobalSum;


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
   if (primme_svds.procID == 0) /* Reports process with ID 0 */
      primme_svds_display_params(primme_svds);

   /* Allocate space for converged Ritz values and residual norms */
   svals = (PetscReal*)malloc(primme_svds.numSvals*sizeof(PetscReal));
   svecs = (PetscScalar*)malloc((primme_svds.n+primme_svds.m)
         *primme_svds.numSvals*sizeof(PetscScalar));
   rnorms = (PetscReal*)malloc(primme_svds.numSvals*sizeof(PetscReal));

   /* Call primme_svds  */
#if defined(PETSC_USE_COMPLEX) && defined(PETSC_USE_REAL_SINGLE)
   ret = cprimme_svds(svals, svecs, rnorms, &primme_svds);
#elif defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_REAL_SINGLE)
   ret = zprimme_svds(svals, svecs, rnorms, &primme_svds);
#elif !defined(PETSC_USE_COMPLEX) && defined(PETSC_USE_REAL_SINGLE)
   ret = sprimme_svds(svals, svecs, rnorms, &primme_svds);
#elif !defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_REAL_SINGLE)
   ret = dprimme_svds(svals, svecs, rnorms, &primme_svds);
#endif

   if (ret != 0) {
      fprintf(primme_svds.outputFile, 
         "Error: primme_svds returned with nonzero exit status: %d \n",ret);
      return -1;
   }


   if (primme_svds.procID == 0) { /* Reports process with ID 0 */
       /* Reporting (optional) */
      for (i=0; i < primme_svds.initSize; i++) {
         fprintf(primme_svds.outputFile, "Sval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
            svals[i], rnorms[i]); 
      }
      fprintf(primme_svds.outputFile, " %d singular triplets converged\n", primme_svds.initSize);
      fprintf(primme_svds.outputFile, "Tolerance : %-22.15E\n", 
                                                            primme_svds.aNorm*primme_svds.eps);
      fprintf(primme_svds.outputFile, "Iterations: %-" PRIMME_INT_P "\n", 
                                                    primme_svds.stats.numOuterIterations); 
      fprintf(primme_svds.outputFile, "Restarts  : %-" PRIMME_INT_P "\n", primme_svds.stats.numRestarts);
      fprintf(primme_svds.outputFile, "Matvecs   : %-" PRIMME_INT_P "\n", primme_svds.stats.numMatvecs);
      fprintf(primme_svds.outputFile, "Preconds  : %-" PRIMME_INT_P "\n", primme_svds.stats.numPreconds);
      if (primme_svds.primme.locking && primme_svds.primme.intWork && primme_svds.primme.intWork[0] == 1) {
         fprintf(primme_svds.outputFile, "\nA locking problem has occurred.\n"
            "Some triplets do not have a residual norm less than the tolerance.\n"
            "However, the subspace of evecs is accurate to the required tolerance.\n");
      }

   }

   primme_svds_free(&primme_svds);
   free(svals);
   free(svecs);
   free(rnorms);

   ierr = PetscFinalize(); CHKERRQ(ierr);

  return(0);
}

/* lauchli block matrix-vector product, y = a * x (or y = a^t * x), where

   - x, input dense matrix of size primme_svds.n (or primme_svds.m) x blocksize;
   - y, output dense matrix of size primme_svds.m (or primme_svds.n) x blocksize;
   - a, lauchli matrix of dimensions primme_svds.m x (primme_svds.m+1) with this form:

        [ 1  1  1  1  1 ...   1 ],  ei = 1 - (1 - mu)*i/(min(m,n) - 1)
        [e0  0  0  0  0 ...   0 ]
        [ 0 e1  0  0  0 ...   0 ]
         ...
        [ 0  0  0  0  0 ... en-1]
*/

PetscErrorCode generateLauchli(int m, int n, PetscReal mu, Mat *A) {
   PetscInt       i,Istart,Iend;
   PetscErrorCode ierr;

   PetscFunctionBegin;

   ierr = MatCreate(PETSC_COMM_WORLD, A); CHKERRQ(ierr);
   ierr = MatSetSizes(*A, PETSC_DECIDE, PETSC_DECIDE, m, n); CHKERRQ(ierr);
   ierr = MatSetFromOptions(*A); CHKERRQ(ierr);
   ierr = MatSetUp(*A); CHKERRQ(ierr);

   ierr = MatGetOwnershipRange(*A, &Istart, &Iend); CHKERRQ(ierr);
   if (Istart == 0) {
      for (i=0; i<n; i++) {
         ierr = MatSetValue(*A, 0, i, 1.0, INSERT_VALUES); CHKERRQ(ierr);
      }
   }
   for (i=max(1,Istart); i<min(Iend,n); i++) {
      ierr = MatSetValue(*A, i, i-1, 1.0 - (1.0 - mu)*(i-1)/(min(m,n) - 1), INSERT_VALUES); CHKERRQ(ierr);
   }

   ierr = MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
   ierr = MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

   PetscFunctionReturn(0);
}

void PETScMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int *trans,
                    primme_svds_params *primme_svds, int *err) {
   int i;
   Mat *A;
   Vec xvec, yvec;
   PetscInt m, n, mLocal, nLocal;
   PetscErrorCode ierr;
   
   A = (Mat *)primme_svds->matrix;

   assert(sizeof(PetscScalar) == sizeof(PetscScalar));   
   ierr = MatGetSize(*A, &m, &n); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
   ierr = MatGetLocalSize(*A, &mLocal, &nLocal); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
   assert(m == primme_svds->m && n == primme_svds->n && mLocal == primme_svds->mLocal
         && nLocal == primme_svds->nLocal);

   #if PETSC_VERSION_LT(3,6,0)
      ierr = MatGetVecs(*A, &xvec, &yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
   #else
      ierr = MatCreateVecs(*A, &xvec, &yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
   #endif
   if (*trans) {
      Vec aux = xvec; xvec = yvec; yvec = aux;
   }
   for (i=0; i<*blockSize; i++) {
      ierr = VecPlaceArray(xvec, ((PetscScalar*)x) + (*ldx)*i); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
      ierr = VecPlaceArray(yvec, ((PetscScalar*)y) + (*ldy)*i); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
      if (*trans == 0) {
         ierr = MatMult(*A, xvec, yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
      } else {
         ierr = MatMultHermitianTranspose(*A, xvec, yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
      }
      ierr = VecResetArray(xvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
      ierr = VecResetArray(yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
   }
   ierr = VecDestroy(&xvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
   ierr = VecDestroy(&yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
   *err = 0;
}

/* This performs Y = M^{-1} * X, where

   - X, input dense matrix of size primme_svds.n (or primme_svds.m or m+n) x blockSize;
   - Y, output dense matrix of size primme_svds.n (or primme_svds.m or m+n) x blockSize;
   - M, preconditioner for A^t*A (or A*A^t or [0 A^t; A 0]), where A is the Lauchli matrix.
*/

void ApplyPCPrecAHA(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                         int *mode, primme_svds_params *primme_svds, int *err) {
   int i,j;
   Mat *matrix;
   PC *pc;
   Vec xvec, yvec;
   PetscScalar *x0 = (PetscScalar*)x, *y0 = (PetscScalar*)y;
   PetscErrorCode ierr;
   
   matrix = (Mat *)primme_svds->matrix;
   pc = (PC *)primme_svds->preconditioner;

   /* The preconditioner is only build for A^t*A; in the rest of cases y <= x */

   if (*mode == primme_svds_op_AtA) {
      ierr = MatCreateVecs(*matrix, &xvec, NULL); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
      ierr = MatCreateVecs(*matrix, &yvec, NULL); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
      for (i=0; i<*blockSize; i++) {
         ierr = VecPlaceArray(xvec, ((PetscScalar*)x) + primme_svds->nLocal*i); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
         ierr = VecPlaceArray(yvec, ((PetscScalar*)y) + primme_svds->nLocal*i); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
         ierr = PCApply(*pc, xvec, yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
         ierr = VecResetArray(xvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
         ierr = VecResetArray(yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
      }
      ierr = VecDestroy(&xvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
      ierr = VecDestroy(&yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
   }
   else if (*mode == primme_svds_op_AAt) {
      for (i=0; i<*blockSize; i++)
         for (j=0; j<primme_svds->mLocal; j++)
            y0[(*ldy)*i+j] = x0[(*ldx)*i+j];
   }
   else if (*mode == primme_svds_op_augmented) {
      for (i=0; i<*blockSize; i++)
         for (j=0; j<primme_svds->mLocal+primme_svds->nLocal; j++)
            y0[(*ldy)*i+j] = x0[(*ldx)*i+j];
   }
   *err = 0;
}

void par_GlobalSum(void *sendBuf, void *recvBuf, int *count, 
                         primme_svds_params *primme_svds, int *ierr) {
   MPI_Comm communicator = *(MPI_Comm *) primme_svds->commInfo;

   if (sendBuf == recvBuf) {
     *ierr = MPI_Allreduce(MPI_IN_PLACE, recvBuf, *count, MPIU_REAL, MPI_SUM, communicator) != MPI_SUCCESS;
   } else {
     *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPIU_REAL, MPI_SUM, communicator) != MPI_SUCCESS;
   }
}
