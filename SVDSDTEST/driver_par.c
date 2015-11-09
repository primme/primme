/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2005  James R. McCombs,  Andreas Stathopoulos
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
 * --------------------------------------------------------------------------
 *
 *  Parallel driver for dprimme. Calling format:
 *
 * 	    par_dprimme DriverConfigFileName SolverConfigFileName
 *
 *  DriverConfigFileName  includes the path and filename of the matrix
 *  			  as well as preconditioning information (eg., 
 *  			  ParaSails parameters).
 *  			  Currently, for reading the input matrix,
 *  			  full coordinate format (.mtx) and upper triangular 
 *  			  coordinate format (.U) are supported.
 *
 *         Example file:  DriverConf
 *
 *  SolverConfigFileName  includes all dprimme required information
 *  			  as stored in primme data structure.
 *
 *  	   Example files: FullConf  Full customization of primme
 *  		          LeanConf  Use a preset method and some customization
 *  		          MinConf   Provide ONLY a preset method and numEvals.
 *
 ******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include "driver_par.h"

/* Parasails library headers */
#include "ParaSails.h"

/* primme.h header file is required to run primme */
#include "primme.h"

/* wtime.h header file is included so primme's timimg functions can be used */
#include "wtime.h"

/******************************************************************************/
int main (int argc, char *argv[]) {

   /* Timing vars */
   double ut1,ut2,st1,st2,wt1,wt2;

   /* Matrix */
   int n, nLocal, nnz;
   double fnorm;
   CSRMatrix matrix; 	     // The matrix in simple SPARSKIT CSR format
   Matrix *Par_matrix;       // Pointer to the matrix in Parasails CSR format

   /* Permutation/partitioning stuff */
   int *mask, *map;
   int *fg2or, *or2fg;
   
   /* Preconditioner */
   ParaSails *Par_Factors;   // Pointer to the matrix in Parasails CSR format

   /* Files */
   char *DriverConfigFileName, *SolverConfigFileName;
   char partFileName[512];
   FILE *partFile;
   
   /* Driver and solver I/O arrays and parameters */
   double *evals, *evecs, *rnorms;
   driver_params driver;
   primme_params primme;
   primme_preset_method method;
#ifdef Cplusplus     /* C++ has a stricter type checking */
   void (*precond_function)(void *, void *, int *, primme_params *);
#else
   void *precond_function;
#endif

   /* Other miscellaneous items */
   int i,j;
   int ret, modulo;
   int rangeStart;
   int rangeEnd;
   int numProcs, procID;
   MPI_Comm comm;

   /* --------------------------------------------------------------------- */
   /* MPI INITIALIZATION						    */
   /* --------------------------------------------------------------------- */
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &procID);
   comm = MPI_COMM_WORLD;

   /* --------------------------------------------------------------------- */
   /*   Read matrix and driver setup                                        */
   /* --------------------------------------------------------------------- */

   /* ------------------------------------------------------- */
   /* Get from command line the names for the 2 config files  */
   /* ------------------------------------------------------- */

   if (argc == 3) {
      DriverConfigFileName = argv[1];
      SolverConfigFileName = argv[2];
   }
   else {
      MPI_Finalize();
      return(-1);
   }


   /* ----------------------------- */
   /* Read in the driver parameters */
   /* ----------------------------- */
   if (read_driver_params(DriverConfigFileName, &driver) < 0) {
      fprintf(stderr, "Reading driver parameters failed\n");
	fflush(stderr);
      MPI_Finalize();
      return(-1);
   }

   /* ------------------------------------------ */
   /* Read the matrix and store it in CSR format */
   /* ------------------------------------------ */
   if (procID == 0) {

      fprintf(stdout," Matrix: %s\n",driver.matrixFileName);
	fflush(stdout);

      if (!strcmp("mtx", 
		&driver.matrixFileName[strlen(driver.matrixFileName)-3])) {  
         // coordinate format storing both lower and upper triangular parts
         ret = readfullMTX(driver.matrixFileName, &matrix.AElts, &matrix.JA, 
            &matrix.IA, &n, &nnz);
         if (ret < 0) {
            fprintf(stderr, "ERROR: Could not read matrix file\n");
	    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            return(-1);
         }
      }
      else if (driver.matrixFileName[strlen(driver.matrixFileName)-1] == 'U') {
         // coordinate format storing only upper triangular part
         ret = readUpperMTX(driver.matrixFileName, &matrix.AElts, &matrix.JA,
            &matrix.IA, &n, &nnz);
         if (ret < 0) {
            fprintf(stderr, "ERROR: Could not read matrix file\n");
	    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            return(-1);
         }
      }
      else {  
         //Harwell Boeing format NOT IMPLEMENTED
         //ret = readmt()
         ret = -1;
         if (ret < 0) {
            fprintf(stderr, "ERROR: Could not read matrix file\n");
	    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            return(-1);
         }
      }
   } // if procID == 0

   /* ----------------------------------------------------------- */
   // Allocate space on other processors and broadcast the matrix
   /* ----------------------------------------------------------- */

   MPI_Bcast(&nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

   if (procID != 0) {

      matrix.AElts = (double *)primme_calloc(nnz, sizeof(double), "A");
      matrix.JA = (int *)primme_calloc(nnz, sizeof(int), "JA");
      matrix.IA = (int *)primme_calloc(n+1, sizeof(int), "IA");
   }
   else {
      // Proc 0 converts CSR to C indexing
      for (i=0; i < n+1; i++) {
         matrix.IA[i]--;
      }
      for (i=0; i < nnz; i++) {
         matrix.JA[i]--;
      }
   }

   MPI_Bcast(matrix.AElts, nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(matrix.IA, n+1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(matrix.JA, nnz, MPI_INT, 0, MPI_COMM_WORLD);

   fnorm = frobeniusNorm(n, matrix.IA, matrix.AElts);

   /* ---------------------------------------------------------------------- */
   /*  Partitioning of the matrix among the processors                       */
   /* ---------------------------------------------------------------------- */

   mask = (int *)primme_calloc(n, sizeof(int), "mask");
   map = (int *)primme_calloc(numProcs+1, sizeof(int), "map");
   fg2or = (int *)primme_calloc(n, sizeof(int), "fg2or");
   or2fg = (int *)primme_calloc(n, sizeof(int), "or2fg");
      
   /* * * * * * * * * * * * * * * * * *
    * Read the partition from a file
    * * * * * * * * * * * * * * * * * */
/*    sprintf(partFileName, "%s/%s", driver.partDir, driver.partId);
*    partFile = fopen(partFileName, "r");
* 
*    if (partFile == 0) {
*       fprintf(stderr, "ERROR: Could not open '%s'\n", partFileName);
*       MPI_Finalize();
*       return(-1);
*    }
* 
*    for (i = 0; i < n; i++) {
*       fscanf(partFile, "%d", &mask[i]);
*    }
* 
*    fclose(partFile);
*/

   /* * * * * * * * * * * * * * * * * * * * * * * */
   /* Simplistic assignment of processors to rows */
   /* * * * * * * * * * * * * * * * * * * * * * * */
   nLocal = n / numProcs;
   modulo = n % numProcs;
   rangeStart = 0;
   for (i=0; i<numProcs; i++) {
      rangeEnd = rangeStart + nLocal;
      if (i < modulo) rangeEnd = rangeEnd + 1;
      for (j = rangeStart; j< rangeEnd; j++) mask[j] = i;
      rangeStart = rangeEnd;
   }
   /* * * * * * * * * * * * * * * * * * * * * * * */

   generatePermutations(n, numProcs, mask, or2fg, fg2or, map);
   rangeStart = map[procID];
   rangeEnd = map[procID+1]-1;
   nLocal = rangeEnd - rangeStart+1;
   Par_matrix = csrToParaSails(procID, map, fg2or, or2fg, 
		               matrix.IA, matrix.JA, matrix.AElts, comm);

/* ------------------------------------------------------------------------- */
/*  Set up preconditioner if needed. For parallel programs the only choice   */
/*  in driver.PrecChoice:  (driver.PrecChoice is read in read_driver_params) */
/*     choice = 4    Parallel ParaSails preconditioners                      */
/*                   with parameters (level, threshold, filter, isymm)       */
/*                   as read in read_driver_params().                        */
/* ------------------------------------------------------------------------- */

   if (driver.PrecChoice == 4) {

      Par_Factors = generate_precond(&matrix, driver.shift, n,
	 procID, map, fg2or, or2fg, rangeStart, rangeEnd, driver.isymm, 
	 driver.level, driver.threshold, driver.filter, comm);
      precond_function = par_ApplyParasailsPrec;
   }
   else {
      Par_Factors = NULL;
      precond_function = NULL;
      // Free A as it is not further needed
      free(matrix.AElts); free(matrix.IA); free(matrix.JA);
   }

   free(mask); free(map); free(fg2or); free(or2fg);

   /* --------------------------------------------------------------------- */
   /*    Primme solver setup                                                */
   /*       primme_initialize  (not needed if ALL primme struct members set)*/
   /*       primme_set_method  (bypass it to fully customize your solver)   */
   /* --------------------------------------------------------------------- */

   /* ----------------------------- */
   /* Initialize defaults in primme */
   /* ----------------------------- */
   primme_initialize(&primme);

   /* --------------------------------------- */
   /* Read in the primme configuration file   */
   /* --------------------------------------- */
   primme.n     = n;
   primme.aNorm = fnorm; /* ||A||_frobenius. A configFile entry overwrites it */

   if (procID == 0) {
      if (read_solver_params(SolverConfigFileName, driver.outputFileName, 
			   &primme, &method) < 0) {
         fprintf(stderr, "Reading solver parameters failed\n");
	 MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
         return(-1);
      }
   }

   /* ------------------------------------------------- */
   // Send read common primme members to all processors
   // Setup the primme members local to this processor  
   /* ------------------------------------------------- */
   broadCast(&primme, &method, procID, comm);

   primme.procID = procID;
   primme.numProcs = numProcs;
   primme.nLocal = nLocal;
   primme.commInfo = &comm;
   primme.globalSumDouble = par_GlobalSumDouble;

   /* --------------------------------------- */
   /* Pick one of the default methods(if set) */
   /* --------------------------------------- */

   if (primme_set_method(method, &primme) < 0 ) {
      fprintf(primme.outputFile, "No preset method. Using custom settings\n");
   }

   /* --------------------------------------- */
   /* Optional: report memory requirements    */
   /* --------------------------------------- */

   ret = dprimme(NULL,NULL,NULL,&primme);
   fprintf(primme.outputFile,"PRIMME will allocate the following memory:\n");
   fprintf(primme.outputFile," processor %d, real workspace, %ld bytes\n",
		   		procID, primme.realWorkSize);
   fprintf(primme.outputFile," processor %d, int  workspace, %d bytes\n",
		   		procID, primme.intWorkSize);

   /* --------------------------------------- */
   /* Set up matrix vector and preconditioner */
   /* --------------------------------------- */

   primme.matrixMatvec        = par_MatrixMatvec;
   primme.applyPreconditioner = precond_function;

   /* --------------------------------------- */
   /* Optional: provide matrix/preconditioner */
   /* --------------------------------------- */
   primme.matrix         = Par_matrix;
   primme.preconditioner = Par_Factors;

   /* --------------------------------------- */
   /* Display given parameter configuration   */
   /* Place this after the dprimme() to see   */
   /* any changes dprimme() made to primme    */
   /* --------------------------------------- */

   if (procID >= 0) {
      fprintf(primme.outputFile," Matrix: %s\n",driver.matrixFileName);
      primme_display_params(primme);
   }

   /* --------------------------------------------------------------------- */
   /* 	                   Run the dprimme solver                           */
   /* --------------------------------------------------------------------- */

   /* Allocate space for converged Ritz values and residual norms */

   evals = (double *)primme_calloc(primme.numEvals, sizeof(double), "evals");
   evecs = (double *)primme_calloc(primme.nLocal*primme.numEvals, 
				sizeof(double), "evecs");
   rnorms = (double *)primme_calloc(primme.numEvals, sizeof(double), "rnorms");

   /* ------------------------ */
   /* Initial guess (optional) */
   /* ------------------------ */
    for (i=0;i<primme.nLocal;i++) evecs[i]=1/sqrt(primme.n);

   /* ------------- */
   /*  Call primme  */
   /* ------------- */

   wt1 = primme_get_wtime(); 
   primme_get_time(&ut1,&st1);

   ret = dprimme(evals, evecs, rnorms, &primme);

   wt2 = primme_get_wtime();
   primme_get_time(&ut2,&st2);

   /* --------------------------------------------------------------------- */
   /* Reporting                                                             */
   /* --------------------------------------------------------------------- */

   if (procID == 0) {
      primme_PrintStackTrace(primme);
   }

   fprintf(primme.outputFile, "Wallclock Runtime   : %-f\n", wt2-wt1);
   fprintf(primme.outputFile, "User Time           : %f seconds\n", ut2-ut1);
   fprintf(primme.outputFile, "Syst Time           : %f seconds\n", st2-st1);

   if (primme.procID == 0) {
      for (i=0; i < primme.numEvals; i++) {
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
      if (primme.locking && primme.intWork[0] == 1) {
         fprintf(primme.outputFile, "\nA locking problem has occurred.\n");
         fprintf(primme.outputFile,
            "Some eigenpairs do not have a residual norm less than the tolerance.\n");
         fprintf(primme.outputFile,
            "However, the subspace of evecs is accurate to the required tolerance.\n");
      }

      fprintf(primme.outputFile, "\n\n#,%d,%.1f\n\n", primme.stats.numMatvecs,
         wt2-wt1); 

      switch (primme.dynamicMethodSwitch) {
         case -1: fprintf(primme.outputFile,
               "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
         case -2: fprintf(primme.outputFile,
               "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
         case -3: fprintf(primme.outputFile,
               "Recommended method for next run: DYNAMIC (close call)\n"); break;
      }
   }

   if (ret != 0 && procID == 0) {
      fprintf(primme.outputFile, 
         "Error: dprimme returned with nonzero exit status: %d \n",ret);
      return -1;
   }


   fflush(primme.outputFile);
   fclose(primme.outputFile);
   primme_Free(&primme);

   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(comm);
   MPI_Finalize();

   return(0);

}
/******************************************************************************/
/* END OF MAIN DRIVER FUNCTION                                                */
/******************************************************************************/

/******************************************************************************
 * Computed the Frobenius norm of a CSR matrix 
 *
 * 	||A||_frob = sqrt( \sum_{i,j} A_ij^2 )
 *
******************************************************************************/
double frobeniusNorm(int n, int *IA, double *AElts) {

   int i, j;
   double fnorm;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   fnorm = 0.0L;

   for (i=0; i < n; i++) {
      for (j=IA[i]; j <= IA[i+1]-1; j++) {
         fnorm = fnorm + AElts[j-1]*AElts[j-1];
      }
   }

   return (sqrt(fnorm)); 
}  
   

/******************************************************************************
 * Shifts a CSR matrix by a shift
 *
 * 	A = A + shift I
 *
******************************************************************************/
void shiftCSRMatrix(double shift, int n, int *IA, int *JA, double *AElts) {

   int i, j;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   for (i=0; i < n; i++) {
      for (j=IA[i]; j <= IA[i+1]-1; j++) {

         if (JA[j-1]-1 == i) {
            AElts[j-1] = AElts[j-1] + shift;
         }

      }
   }

}

/******************************************************************************/
/* Matvec, preconditioner and other utilities                                 */

/******************************************************************************/

/******************************************************************************
 * Function to broadcast the primme data structure to all processors
 *
 * EXCEPTIONS: procID and seed[] are not copied from processor 0. 
 *             Each process creates their own.
******************************************************************************/
void broadCast(primme_params *primme, primme_preset_method *method, 
   int procID, MPI_Comm comm){

   int i;

   MPI_Bcast(&(primme->numEvals), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->target), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->numTargetShifts), 1, MPI_INT, 0, comm);

   if (primme->numTargetShifts > 0 && procID !=0) {
      primme->targetShifts = (double *)primme_calloc(
         primme->numTargetShifts, sizeof(double), "targetShifts");
   }
   for (i=0; i<primme->numTargetShifts; i++) {
      MPI_Bcast(&(primme->targetShifts[i]), 1, MPI_DOUBLE, 0, comm);
   }
   MPI_Bcast(&(primme->locking), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->dynamicMethodSwitch), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->initSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->numOrthoConst), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxBasisSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->minRestartSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxBlockSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxMatvecs), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxOuterIterations), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->aNorm), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme->eps), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme->printLevel), 1, MPI_INT, 0, comm);

   MPI_Bcast(&(primme->restartingParams.scheme), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->restartingParams.maxPrevRetain), 1, MPI_INT, 0, comm);

   MPI_Bcast(&(primme->correctionParams.precondition), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->correctionParams.robustShifts), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->correctionParams.maxInnerIterations),1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.convTest), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->correctionParams.relTolBase), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme->correctionParams.projectors.LeftQ),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.LeftX),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.RightQ), 1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.RightX), 1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.SkewQ),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.SkewX),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.SkewX),  1, MPI_INT, 0,comm);

   MPI_Bcast(method, 1, MPI_INT, 0, comm);
}

/******************************************************************************
 * Convert CSR matrix format to Parasails matrix format 
 *
******************************************************************************/
Matrix* csrToParaSails(int procID, int *map, int *fg2or, int *or2fg, int *IA,
   int *JA, double *AElts, MPI_Comm comm) {

   int i, j;
   int ncols;
   int origRow;
   int rowStart, rangeStart;
   int rowEnd, rangeEnd;
   Matrix *newMatrix;

   rangeStart = map[procID];
   rangeEnd = map[procID+1]-1;
   newMatrix = MatrixCreate(comm, rangeStart, rangeEnd);

   for (i = rangeStart; i <= rangeEnd; i++) {
      origRow = fg2or[i];
      rowStart = IA[origRow];
      rowEnd = IA[origRow+1]-1;
      ncols = rowEnd - rowStart + 1;

      for (j=rowStart; j <= rowEnd; j++) {
         JA[j] = or2fg[JA[j]];
      }

      MatrixSetRow(newMatrix, i, ncols, &JA[rowStart], &AElts[rowStart]);

      for (j=rowStart; j <= rowEnd; j++) {
         JA[j] = fg2or[JA[j]];
      }
   }

   MatrixComplete(newMatrix);

   return newMatrix;
}

/******************************************************************************
 * Generate the parallel Parasails preconditioner with parameters read from
 * the driver. 
 *
******************************************************************************/
ParaSails* generate_precond(CSRMatrix *matrix, double shift, int n, int procID,
   int *map, int *fg2or, int *or2fg, int rangeStart, int rangeEnd, int isymm, 
   int level, double threshold, double filter, MPI_Comm comm)
{
   Matrix *shiftedMatrix;    // Temporary matrix
   ParaSails *A_p;           // Pointer holding the resulting preconditioner
   double t1, t2;

   if (procID == 0) {
      fprintf(stdout, "Computing preconditioner\n");
      fprintf(stdout, "isymm: %d level: %d thresh: %f filter: %f\n", 
         isymm, level, threshold, filter);
   }
   t1 = MPI_Wtime();

   // Compute A = A-shift
   shiftCSRMatrix(-shift, n, matrix->IA, matrix->JA, matrix->AElts);

   // Change A to Parasails format
   shiftedMatrix = csrToParaSails(procID, map, fg2or, or2fg, 
		            matrix->IA, matrix->JA, matrix->AElts, comm);

   // Free A to make room for preconditioner
   free(matrix->AElts); free(matrix->IA); free(matrix->JA);
   
   // Create parasails preconditioner
   A_p = ParaSailsCreate(comm, rangeStart, rangeEnd, isymm);
   ParaSailsSetupPattern(A_p, shiftedMatrix, threshold, level);
   ParaSailsSetupValues(A_p, shiftedMatrix, filter);

   MatrixDestroy(shiftedMatrix);

   t2 = MPI_Wtime();
   if (procID == 0) {
      fprintf(stdout, "Done computing preconditioner\n");
      fprintf(stdout, "Preconditioner time: %f\n", t2-t1);
      //fprintf(report, "Preconditioner time: %f\n", t2-t1);
   }

   return(A_p);
}

/******************************************************************************
 * Parallel globalSumDouble function
 *
******************************************************************************/
void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, 
   primme_params *primme) {
   
   MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;

   MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator);
}

/******************************************************************************
 * Applies the PARALLEL matrix vector mulitplication of a block of vectors.
 * Because a block function is not available, we call blockSize times
 * the Parasails function MatrixMatvec()
 *
******************************************************************************/
void par_MatrixMatvec(void *x, void *y, int *blockSize, 
   primme_params *primme) {
   
   int i;
   double *xvec, *yvec;
   xvec = (double *)x;
   yvec = (double *)y;

   for (i=0;i<*blockSize;i++) {
      MatrixMatvec(primme->matrix, &xvec[primme->nLocal*i], 
		      		   &yvec[primme->nLocal*i]);
   }

}

/******************************************************************************
 * Apply the PARALLEL Parasails preconditioner to a block of vectors.
 * Because Parasails is not block, we apply the preconditioner for each
 * block vector.
 *
******************************************************************************/
void par_ApplyParasailsPrec(void *x, void *y, int *blockSize,
   primme_params *primme) {

   int i;
   double *xvec, *yvec;
   xvec = (double *)x;
   yvec = (double *)y;

   for (i=0;i<*blockSize;i++) {
     ParaSailsApply(primme->preconditioner, &xvec[primme->nLocal*i], 
		     			    &yvec[primme->nLocal*i]);
   }

}

/******************************************************************************
 * void generatePermutations() 
 * Given a proc array :  proc[i] = processor # where row i lies
 * it generates all other needed permutation arrays for processing in 
 * Parasails.
 *
******************************************************************************/
void generatePermutations(int n, int nParts, int *proc, int *perm,
   int *iperm, int *map) {

   int i;
   int *count;

   count = (int *)primme_calloc(nParts, sizeof(int), "counts");

   for (i=0; i < nParts; i++) {
      count[i] = 0;
   }

   for (i=0; i < n; i++) {
      count[proc[i]]++;
   }

   map[0] = 0;
   for (i=1; i <= nParts; i++) {
      map[i] = map[i-1] + count[i-1];
   }

   for (i=0; i < n; i++) {
      iperm[map[proc[i]]] = i;
      map[proc[i]]++;
   }

   for (i=0; i < n; i++) {
      perm[iperm[i]] = i;
   }

   map[0] = 0;
   for (i=1; i <= nParts; i++) {
      map[i] = map[i-1] + count[i-1];
   }

   free(count);
}
