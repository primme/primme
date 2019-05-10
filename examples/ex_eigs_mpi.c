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
 *  Example to compute the largest eigenvalues in a diagonal matrix in
 *  parallel using MPI.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>

#include "primme.h"   /* header file is required to run primme */ 

void DiagonalMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
static void par_GlobalSum(void *sendBuf, void *recvBuf, int *count,
                         primme_params *primme, int *ierr);

#ifndef min
#  define min(a, b) ((a) < (b) ? (a) : (b))
#endif

int main (int argc, char *argv[]) {

   /* Solver arrays and parameters */
   float *evals;    /* Array with the computed eigenvalues */
   float *rnorms;   /* Array with the computed eigenpairs residual norms */
   float *evecs;    /* Array with the computed eigenvectors;
                        first vector starts in evecs[0],
                        second vector starts in evecs[primme.n],
                        third vector starts in evecs[primme.n*2]...  */
   primme_params primme;
                     /* PRIMME configuration struct */

   /* Other miscellaneous items */
   int ret;
   int i;

   /* Initialize the infrastructure necessary for communication */
   MPI_Init(&argc, &argv);

   /* Set default values in PRIMME configuration struct */
   primme_initialize(&primme);

   /* Set problem matrix */
   primme.matrixMatvec = DiagonalMatrixMatvec;
                           /* Function that implements the matrix-vector product
                              A*x for solving the problem A*x = l*x */
  
   /* Set problem parameters */
   primme.n = 1000; /* set problem dimension */
   primme.numEvals = 1000;   /* Number of wanted eigenpairs */
   primme.eps = .1;      /* ||r|| <= eps * ||matrix|| */
   primme.target = primme_largest;
                           /* Wanted the smallest eigenvalues */

   /* Set advanced parameters if you know what are you doing (optional) */
   /*
   primme.maxBasisSize = 14;
   primme.minRestartSize = 4;
   primme.maxBlockSize = 1;
   primme.maxMatvecs = 1000;
   */

   /* Set method to solve the problem */
   primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
   /* DYNAMIC uses a runtime heuristic to choose the fastest method between
       PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
       set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */

   /* Set parallel parameters */
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_size(comm, &primme.numProcs);
   MPI_Comm_rank(comm, &primme.procID);
   primme.commInfo = &comm; /* User-defined member to pass the communicator to
                               globalSumReal and broadcastReal */
   /* In this example, the matrix is distributed by rows, and the first
    * processes may have an extra row in order to distribute the remaining rows
    * n % numProcs */
   PRIMME_INT nLocal = primme.n / primme.numProcs +
                       (primme.n % primme.numProcs > primme.procID ? 1 : 0);
   primme.nLocal = nLocal; /* Number of local rows */
   primme.globalSumReal = par_GlobalSum;

   /* Display PRIMME configuration struct (optional) */
   if (primme.procID == 0) primme_display_params(primme);

   /* Allocate space for converged Ritz values and residual norms */
   evals = (float*)malloc(primme.numEvals*sizeof(float));
   evecs = (float*)malloc(primme.n*primme.numEvals*sizeof(float));
   rnorms = (float*)malloc(primme.numEvals*sizeof(float));

   /* Call primme  */
   ret = sprimme(evals, evecs, rnorms, &primme);

   if (primme.procID == 0) {
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
     fprintf(primme.outputFile, "Orthogonalization Time : %g\n", primme.stats.timeOrtho);
     fprintf(primme.outputFile, "Matvec Time            : %g\n", primme.stats.timeMatvec);
     fprintf(primme.outputFile, "GlobalSum Time         : %g\n", primme.stats.timeGlobalSum);
     fprintf(primme.outputFile, "Broadcast Time         : %g\n", primme.stats.timeBroadcast);
     fprintf(primme.outputFile, "Total Time             : %g\n", primme.stats.elapsedTime);
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
   }

   primme_free(&primme);
   free(evals);
   free(evecs);
   free(rnorms);

   /* Tear down the communication infrastructure */
   MPI_Finalize();

   return(0);
}

/* Diagonal block matrix-vector product, Y = A * X, where

   - X, input dense matrix of size primme.n x blockSize;
   - Y, output dense matrix of size primme.n x blockSize;
   - A, tridiagonal square matrix of dimension primme.n with this form:

*/

void DiagonalMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   
   int i;            /* vector index, from 0 to *blockSize-1*/
   int row;          /* local matrix row index, from 0 to nLocal */
   /* In this example, row0 is the global index of the first local row */
   int row0 = primme->n / primme->numProcs * primme->procID +
              min(primme->n % primme->numProcs, primme->procID);
   float *xvec;     /* pointer to i-th input vector x */
   float *yvec;     /* pointer to i-th output vector y */
   
   for (i=0; i<*blockSize; i++) {
      xvec = (float *)x + *ldx*i;
      yvec = (float *)y + *ldy*i;
      for (row = 0; row < primme->nLocal; row++) {
         /* The diagonal matrix has the spectrum of a Laplacial */
         float v = sin(M_PI * (row + row0 + 1) / 2.0 / (primme->n + 1));
         yvec[row] = 4. * v * v * xvec[row];
      }
   }
   *err = 0;
}

static void par_GlobalSum(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme, int *ierr) {
   MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;

   if (sendBuf == recvBuf) {
     *ierr = MPI_Allreduce(MPI_IN_PLACE, recvBuf, *count, MPI_FLOAT, MPI_SUM, communicator) != MPI_SUCCESS;
   } else {
     *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPI_FLOAT, MPI_SUM, communicator) != MPI_SUCCESS;
   }
}
