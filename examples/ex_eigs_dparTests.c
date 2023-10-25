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
 *  Example to compute the k largest eigenvalues in a 1-D Laplacian matrix.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <petscpc.h>
#include <petscmat.h>
#include "primme.h"   /* header file is required to run primme */ 

//PetscErrorCode generateLaplacian1D(int n, Mat *A);
PetscErrorCode CSR_to_PETSc_Matrix(char *filename, Mat *A, primme_params *primme, int *err);
void PETScMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
static void par_GlobalSum(void *sendBuf, void *recvBuf, int *count,
                         primme_params *primme, int *ierr);

void printHelp() {
   printf("Run Options:\n");
   printf("---------------------------------------\n");
   printf("--help:               Prints this help screen\n");
   printf("-basisSize <Int>:     Set maximum size of the Krylov basis   (Default: 100)\n");
   printf("-blockSize <Int>:     Set maximum size of the block          (Default:   1)\n");
   printf("-tol <Double>:        Set the convergence tolerence          (Default: 0.1)\n");
   printf("-anorm <Double>:      Set the norm of input matrix if known  (Default:   1)\n");
   printf("-numEvals <Int>:      Set number of eigenpairs searching for (Default:  10)\n");
   printf("                      NOTE: Must be <= basisSize\n");
   printf("-restartSize <Int>:   Set minimum restart size               (Default:  10)\n");
   printf("                      NOTE: Must be >= numEvals\n");
   printf("-locking <0 or 1>:    Turn on (1) or off (0) locking         (Default:   0)\n");
   printf("                      NOTE: Locking can not be used with sketching\n");
   printf("-target <Option>:     Set which end of the eigenspectra we are getting\n");
   printf("                      our eigenpair approximations from      (Default:   L)\n");
   printf("    OPTIONS: \n");
   printf("       'S' for smallest eigenpairs\n");
   printf("       'L' for largest eigenpairs\n");
   printf("-sketching <option>   Determine whether we use sketching     (Default:   0)\n");
   printf("    OPTIONS: \n");
   printf("       '0' for no sketching\n");
   printf("       '1' for sketching\n");
   printf("-expansion <option>   Determine which Krylov solver to use   (Default: default)\n");
   printf("    OPTIONS: \n");
   printf("       'fullLanczos'    for Lanczos with full orthogonalization\n");
   printf("       'partialLanczos' for 3-term recurrence only Lanczos\n");
   printf("       'davidson'       for Generalized Davidson\n");
   printf("       'default'        for default method\n");
   printf("-residual <option>    Set which method to use for finding residuals (Default: RR)\n");
   printf("                      NOTE: Can only be set for the Lanczos method\n");
   printf("    OPTIONS: \n");
   printf("       'sketched'       Finds the sketched residuals and residual norms\n");
   printf("       'RR'             Finds the residuals from the Ritz pairs\n");
   printf("       'RQ'             Finds the residuals using the Rayleigh Quotient\n");
}

int main (int argc, char *argv[]) {

   /* Solver arrays and parameters */
   PetscReal *evals;    /* Array with the computed eigenvalues */
   PetscReal *rnorms;   /* Array with the computed eigenpairs residual norms */
   PetscScalar *evecs;    /* Array with the computed eigenvectors;
                        first vector starts in evecs[0],
                        second vector starts in evecs[primme.n],
                        third vector starts in evecs[primme.n*2]...  */
   primme_params primme; /* PRIMME configuration struct */

   /* Other miscellaneous items */
   int ret;
   int i;
   Mat A; /* problem matrix */
   PetscErrorCode ierr;
   MPI_Comm comm;

   PetscInitialize(&argc, &argv, NULL, NULL);

   /* Set default values in PRIMME configuration struct */
   primme_initialize(&primme);

   comm = PETSC_COMM_WORLD;
   primme.commInfo = &comm;
   MPI_Comm_size(comm, &primme.numProcs);
   MPI_Comm_rank(comm, &primme.procID);
   primme.globalSumReal = par_GlobalSum;

   /* In this example, the matrix is distributed by rows, and the first
    * processes may have an extra row in order to distribute the remaining rows
    * n % numProcs */
   primme.globalSumReal = par_GlobalSum;

   /* Set problem matrix */
   ierr = CSR_to_PETSc_Matrix(argv[1], &A, &primme, &ierr); CHKERRQ(ierr);
   primme.matrix = &A;
   primme.matrixMatvec = PETScMatvec;
                           /* Function that implements the matrix-vector product
                              A*x for solving the problem A*x = l*x */

   /* Default Parameters */
   primme.maxBasisSize = 100;
   primme.numEvals = 10;   /* Number of wanted eigenpairs */
   primme.maxBlockSize = 1;
   primme.printLevel = 4;
   primme.minRestartSize = 10;
   primme.locking = 0;
   primme.eps = 1e-1;      /* ||r|| <= eps * ||matrix|| */
   primme.aNorm = 1.0;
   primme.target = primme_largest;
   //primme.maxMatvecs = 1000;

   primme.expansionParams.expansion = primme_expansion_davidson;
   primme.projectionParams.projection = primme_proj_default;
 
   for(i = 2; i < argc; i++)
   {

      if(strcmp(argv[i], "--help") == 0) {
         printHelp();
         return 0;
      }
      if(strcmp(argv[i], "-basisSize") == 0)    primme.maxBasisSize      = atoi(argv[i+1]); 
      if(strcmp(argv[i], "-numEvals") == 0)     primme.numEvals          = atoi(argv[i+1]); 
      if(strcmp(argv[i], "-blockSize") == 0)    primme.maxBlockSize      = atoi(argv[i+1]); 
      if(strcmp(argv[i], "-printLevel") == 0)   primme.printLevel        = atoi(argv[i+1]); 
      if(strcmp(argv[i], "-restartSize") == 0)  primme.minRestartSize    = atoi(argv[i+1]); 
      if(strcmp(argv[i], "-locking") == 0)      primme.locking           = atoi(argv[i+1]); 
      if(strcmp(argv[i], "-tol") == 0)          sscanf(argv[i+1], "%lf", &primme.eps);
      if(strcmp(argv[i], "-anorm") == 0)        sscanf(argv[i+1], "%lf", &primme.aNorm);
      if(strcmp(argv[i], "-target") == 0) 
      {
         if(strcmp(argv[i+1], "S") == 0){ 
            primme.target = primme_smallest;
         } else if(strcmp(argv[i+1], "L") == 0) {
            primme.target = primme_largest;
         }
      } /* End target */
      if(strcmp(argv[i], "-sketching") == 0)
      {
         if(atoi(argv[i+1]) == 1){
            primme.projectionParams.projection = primme_proj_sketched; 
         } else {
            primme.projectionParams.projection = primme_proj_default;
         } 
      } /* End sketching */
      if(strcmp(argv[i], "-expansion") == 0)
      {
         if(strcmp(argv[i+1], "fullLanczos") == 0) {
            primme.expansionParams.expansion = primme_expansion_fullLanczos; 
         } else if(strcmp(argv[i+1], "partialLanczos") == 0) {
            primme.expansionParams.expansion = primme_expansion_lanczos; 
         } else if(strcmp(argv[i+1], "davidson") == 0) {
            primme.expansionParams.expansion = primme_expansion_davidson; 
         } else {
            primme.expansionParams.expansion = primme_expansion_default; 
         }
      } /* End expansion */
      if(strcmp(argv[i], "-residual") == 0) {
         if(strcmp(argv[i+1], "sketched") == 0) {
            primme.residualParams.residual = primme_residual_sketched; 
         } else if (strcmp(argv[i+1], "RR") == 0) {
            primme.residualParams.residual = primme_residual_RR; 
         } else if (strcmp(argv[i+1], "RQ") == 0){
            primme.residualParams.residual = primme_residual_RQ; 
         }
      } /* End residual */
   } /* End command line argument options */

   /* Set method to solve the problem */
   primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
   //primme_set_method(PRIMME_DYNAMIC, &primme);
   /* DYNAMIC uses a runtime heuristic to choose the fastest method between
       PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
       set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */

   /* Display PRIMME configuration struct (optional) */
   if (primme.procID == 0) /* Reports process with ID 0 */
      primme_display_params(primme);

   /* Allocate space for converged Ritz values and residual norms */
   evals = (PetscReal*)malloc(primme.numEvals*sizeof(PetscReal));
   evecs = (PetscScalar*)malloc(primme.n*primme.numEvals*sizeof(PetscScalar));
   rnorms = (PetscReal*)malloc(primme.numEvals*sizeof(PetscReal));

#if defined(PETSC_USE_COMPLEX) && defined(PETSC_USE_REAL_SINGLE)
   ret = cprimme(evals, evecs, rnorms, &primme);
#elif defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_REAL_SINGLE)
   ret = zprimme(evals, evecs, rnorms, &primme);
#elif !defined(PETSC_USE_COMPLEX) && defined(PETSC_USE_REAL_SINGLE)
   ret = sprimme(evals, evecs, rnorms, &primme);
#elif !defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_REAL_SINGLE)
   ret = dprimme(evals, evecs, rnorms, &primme);
#endif

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: primme returned with nonzero exit status: %d \n",ret);
      return -1;
   }
   
   if (primme.procID == 0) { 

      for (i=0; i < primme.initSize; i++) fprintf(primme.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1, evals[i], rnorms[i]); 
      fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);
      fprintf(primme.outputFile, "Tolerance : %-22.15E\t (%4.6E tolerance x %4.6E norm)\n", primme.aNorm*primme.eps, primme.eps, primme.aNorm);
      fprintf(primme.outputFile, "Iterations: %-" PRIMME_INT_P "\n", primme.stats.numOuterIterations); 
      fprintf(primme.outputFile, "Restarts  : %-" PRIMME_INT_P "\n", primme.stats.numRestarts);
      fprintf(primme.outputFile, "Matvecs   : %-" PRIMME_INT_P "\n", primme.stats.numMatvecs);
      fprintf(primme.outputFile, "Preconds  : %-" PRIMME_INT_P "\n", primme.stats.numPreconds);
      fprintf(primme.outputFile, "Elapsed Time    : %-22.10E\n", primme.stats.elapsedTime);
      fprintf(primme.outputFile, "MatVec Time     : %-22.10E\n", primme.stats.timeMatvec);
      fprintf(primme.outputFile, "Precond Time    : %-22.10E\n", primme.stats.timePrecond);
      fprintf(primme.outputFile, "Ortho Time      : %-22.10E\n", primme.stats.timeOrtho);
      fprintf(primme.outputFile, "GlobalSum Time  : %-22.10E\n", primme.stats.timeGlobalSum);
      fprintf(primme.outputFile, "Broadcast Time  : %-22.10E\n", primme.stats.timeBroadcast);
      fprintf(primme.outputFile, "Krylov Time     : %-22.10E\n", primme.stats.timeKrylov);
      fprintf(primme.outputFile, "Sketching Time  : %-22.10E\n", primme.stats.timeSketching);
      fprintf(primme.outputFile, "Residual Time   : %-22.10E\n", primme.stats.timeResiduals);

      if (primme.stats.lockingIssue) {
         fprintf(primme.outputFile, "\nA locking problem has occurred.\n");
         fprintf(primme.outputFile, "Some eigenpairs do not have a residual norm less than the tolerance.\n");
         fprintf(primme.outputFile, "However, the subspace of evecs is accurate to the required tolerance.\n");
      }

      switch (primme.dynamicMethodSwitch) {
         case -1: fprintf(primme.outputFile,
               "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
         case -2: fprintf(primme.outputFile,
               "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
         case -3: fprintf(primme.outputFile,
               "Recommended method for next run: DYNAMIC (close call)\n"); break;
      }
   }

   primme_free(&primme);
   free(evals);
   free(evecs);
   free(rnorms);

   ierr = PetscFinalize(); CHKERRQ(ierr);

  return(0);
}


PetscErrorCode CSR_to_PETSc_Matrix(char *filename, Mat *A, primme_params *primme, int *err) {
   
   PetscInt i, N, M, nnz, start_row, end_row, nnz_Local;
   PetscInt *A_I_temp, *A_J_temp, *perm;
   PetscInt *Assign_Rows, *A_I, *A_J, *A_I_Local, *A_J_Local;
   PetscScalar *A_V_temp;
   PetscScalar *A_V, *A_V_Local;
   PetscErrorCode ierr;
   PRIMME_INT nLocal;
   
   PetscFunctionBegin;

   /* Open the file and read in CSR format */
   char* line = NULL;
   size_t len;
   ssize_t read;
   FILE *fp = fopen(filename, "r");

   while((read = getline(&line, &len, fp)) != -1) 
   {
      if(line[0] == '%')
      {
         continue;
      } else {
         sscanf(line, "%d %d %d", &N, &M, &nnz);
         break;
      }
   }

   ierr = PetscMalloc1(nnz, &A_I_temp); CHKERRQ(ierr);  
   ierr = PetscMalloc1(nnz, &A_J_temp); CHKERRQ(ierr);  
   ierr = PetscMalloc1(nnz, &A_V_temp); CHKERRQ(ierr);  
   ierr = PetscMalloc1(nnz, &perm); CHKERRQ(ierr);  
   ierr = PetscCalloc1(N+1, &A_I); CHKERRQ(ierr);  
   ierr = PetscMalloc1(nnz, &A_J); CHKERRQ(ierr);  
   ierr = PetscMalloc1(nnz, &A_V); CHKERRQ(ierr);  

   // Read in the rest of the file
   for(i = 0; i < nnz; i++)
   {
      fscanf(fp, "%d %d %lf", &A_I_temp[i], &A_J_temp[i], &A_V_temp[i]);
      A_J_temp[i]--;
      perm[i] = i;
   }
   fclose(fp);
   
   // Constructing CSR
   ierr = PetscSortIntWithArray(nnz, A_I_temp, perm);
   for(i = 0; i < nnz; i++)
   {
      A_I[A_I_temp[i]]++;
      A_J[i] = A_J_temp[perm[i]];
      A_V[i] = A_V_temp[perm[i]];
   }
   
   // Cumulative Sum
   for(i = 2; i <= N; i++) A_I[i] += A_I[i-1];

   // Set primme->n and primme->nLocal
   primme->n = N;
   nLocal = primme->n / primme->numProcs + (primme->n % primme->numProcs > primme->procID ? 1 : 0);
   primme->nLocal = nLocal; /* Number of local rows */
  
   // Which process own which rows
   ierr = PetscCalloc1(primme->numProcs+1, &Assign_Rows); CHKERRQ(ierr);  
   for(i = 0; i < primme->numProcs; i++) Assign_Rows[i+1] = primme->n / primme->numProcs + (primme->n % primme->numProcs > i ? 1 : 0);
   for(i = 2; i <= primme->numProcs; i++) Assign_Rows[i] += Assign_Rows[i-1];
 
   // Create local CSR
   start_row = Assign_Rows[primme->procID];
   end_row = Assign_Rows[primme->procID+1]-1;
   nnz_Local = A_I[end_row+1]-A_I[start_row];

   ierr = PetscCalloc1(nLocal+1, &A_I_Local); CHKERRQ(ierr); 
   ierr = PetscMalloc1(nnz_Local, &A_J_Local); CHKERRQ(ierr);  
   ierr = PetscMalloc1(nnz_Local, &A_V_Local); CHKERRQ(ierr); 
    
   for(i = 0; i <= nLocal; i++) A_I_Local[i] = A_I[start_row+i]-A_I[start_row];
   for(i = 0; i < nnz_Local; i++)
   {
      A_J_Local[i] = A_J[A_I[start_row]+i];
      A_V_Local[i] = A_V[A_I[start_row]+i];
   }
 
   MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;
   ierr = MatCreateMPIAIJWithArrays(communicator, primme->nLocal, PETSC_DECIDE, primme->n, primme->n, A_I_Local, A_J_Local, A_V_Local, A); CHKERRQ(ierr); 

   // Free arrays   
   ierr = PetscFree(A_I_temp); CHKERRQ(ierr);
   ierr = PetscFree(A_J_temp); CHKERRQ(ierr);
   ierr = PetscFree(A_V_temp); CHKERRQ(ierr);
   ierr = PetscFree(Assign_Rows); CHKERRQ(ierr);
   ierr = PetscFree(perm); CHKERRQ(ierr);
   ierr = PetscFree(A_I); CHKERRQ(ierr);
   ierr = PetscFree(A_J); CHKERRQ(ierr);
   ierr = PetscFree(A_V); CHKERRQ(ierr);
   ierr = PetscFree(A_I_Local); CHKERRQ(ierr);
   ierr = PetscFree(A_J_Local); CHKERRQ(ierr);
   ierr = PetscFree(A_V_Local); CHKERRQ(ierr);

   *err = 0;
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

static void par_GlobalSum(void *sendBuf, void *recvBuf, int *count,
        primme_params *primme, int *ierr) {
    MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;

    if (sendBuf == recvBuf) {
        *ierr = MPI_Allreduce(MPI_IN_PLACE, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator) != MPI_SUCCESS;
    } else {
        *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator) != MPI_SUCCESS;
    }
}

