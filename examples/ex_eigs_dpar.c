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
#include <mpi.h>
#include <assert.h>
#include "primme.h"   /* header file is required to run primme */ 

void DiagonalMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
void DiagonalApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);

static void par_GlobalSum(void *sendBuf, void *recvBuf, int *count,
        primme_params *primme, int *ierr);
static void par_Broadcast(void *sendBuf, int *count,
        primme_params *primme, int *ierr);

#ifndef min
#  define min(a, b) ((a) < (b) ? (a) : (b))
#endif

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
    primme.numEvals = 10;   /* Number of wanted eigenpairs */
    primme.eps = .1;      /* ||r|| <= eps * ||matrix|| */
    primme.target = primme_largest;
    /* Wanted the smallest eigenvalues */

    /* Set preconditioner (optional) */
    //primme.applyPreconditioner = DiagonalApplyPreconditioner;
    //primme.correctionParams.precondition = 1;

    primme.projectionParams.projection = primme_proj_sketched;

    /* Set advanced parameters if you know what are you doing (optional) */
    //primme.minRestartSize = 100;
    primme.maxBasisSize = 750;
    primme.initSize = 0;
    primme.locking = 0;
    primme.maxMatvecs = 2000;
    primme.aNorm = 1.0;
    primme.printLevel = 4;

    primme.maxBlockSize = 1;
    primme.expansionParams.expansion = primme_expansion_fullLanczos;
    /* Set method to solve the problem */
    //primme_set_method(PRIMME_DYNAMIC, &primme);
    primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
    /* DYNAMIC uses a runtime heuristic to choose the fastest method between
       PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
       set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */

    /* Display PRIMME configuration struct (optional) */
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
    primme.broadcastReal = par_Broadcast;

    /* Display PRIMME configuration struct (optional) */
    if (primme.procID == 0) primme_display_params(primme);

    /* Allocate space for converged Ritz values and residual norms */
    evals = (double*)malloc(primme.numEvals*sizeof(double));
    evecs = (double*)malloc(primme.n*primme.numEvals*sizeof(double));
    rnorms = (double*)malloc(primme.numEvals*sizeof(double));

    /* Call primme  */
    ret = dprimme(evals, evecs, rnorms, &primme);

    if(primme.procID == 0){

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

void DiagonalMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {

    int i;            /* vector index, from 0 to *blockSize-1*/
    int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
    double *xvec;     /* pointer to i-th input vector x */
    double *yvec;     /* pointer to i-th output vector y */

    for (i = 0; i < *blockSize; i++) {
        xvec = (double *)x + *ldx*i;
        yvec = (double *)y + *ldy*i;
        for (row = 0; row < primme->n; row++) {
            yvec[row] = 0.0;
            yvec[row] += (double)((row+1)*(row+1))*xvec[row];
        }      
    }
    *err = 0;
}

void DiagonalApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {

    int i;            /* vector index, from 0 to *blockSize-1*/
    int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
    double *xvec;     /* pointer to i-th input vector x */
    double *yvec;     /* pointer to i-th output vector y */

    for (i=0; i<*blockSize; i++) {
        xvec = (double *)x + *ldx*i;
        yvec = (double *)y + *ldy*i;
        for (row=0; row<primme->n; row++) {
            yvec[row] = xvec[row]/(double)((row+1)*(row+1));
        }      
    }
    *ierr = 0;
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

static void par_Broadcast(void *sendBuf, int *count, primme_params *primme, int *ierr) {
    MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;
   *ierr = MPI_Bcast(sendBuf, *count, MPI_DOUBLE, primme->procID, communicator) != MPI_SUCCESS;
}

