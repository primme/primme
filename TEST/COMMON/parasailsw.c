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
 * File: parasailsw.c
 * 
 * Purpose - Parasails wrapper.
 * 
 ******************************************************************************/

#include "csr.h"
#include "ParaSails.h"
#include <mpi.h>
#include "parasailsw.h"

static void generatePermutations(int n, int nParts, int *proc, int *perm,
   int *iperm, int *map);
static Matrix* csrToParaSails(int procID, int *map, int *fg2or, int *or2fg, int *IA,
   int *JA, double *AElts, MPI_Comm comm);
static ParaSails* generate_precond(CSRMatrix *matrix, double shift, int n, int procID,
   int *map, int *fg2or, int *or2fg, int rangeStart, int rangeEnd, int isymm, 
   int level, double threshold, double filter, MPI_Comm comm);

int readMatrixAndPrecondParaSails(const char* matrixFileName, double shift,
                                  int level, double threshold, double filter,
                                  int isymm, MPI_Comm comm, double *fnorm,
                                  int *m, int *n, int *mLocal_, int *nLocal_,
                                  int *numProcs_, int *procID_,
                                  Matrix **pmatrix, ParaSails **pfactor) {
   int numProcs, procID, nLocal, modulo, rangeStart, rangeEnd;
   CSRMatrix *matrix;
   int i, j;

   /* Permutation/partitioning stuff */
   int *mask, *map;
   int *fg2or, *or2fg;
 
   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &procID);

   /* ---------------------------------------------------------------------- */
   /*  Distribution of the matrix among the processors                       */
   /* ---------------------------------------------------------------------- */

   if (procID == 0) {
      if (readMatrixNative(matrixFileName, &matrix, fnorm) !=0 )
         return -1;
   }
   else {
      matrix = (CSRMatrix *)primme_calloc(1, sizeof(CSRMatrix), "CSRMatrix");
   }
   MPI_Bcast(&matrix->nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&matrix->m, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&matrix->n, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if (fnorm) MPI_Bcast(fnorm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   if (procID != 0) {
      matrix->AElts = (double *)primme_calloc(matrix->nnz, sizeof(double), "A");
      matrix->JA = (int *)primme_calloc(matrix->nnz, sizeof(int), "JA");
      matrix->IA = (int *)primme_calloc(matrix->m+1, sizeof(int), "IA");
   }
   else {
      // Proc 0 converts CSR to C indexing
      for (i=0; i < matrix->m+1; i++) {
         matrix->IA[i]--;
      }
      for (i=0; i < matrix->nnz; i++) {
         matrix->JA[i]--;
      }
   }

   MPI_Bcast(matrix->AElts, matrix->nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Bcast(matrix->IA, matrix->m+1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(matrix->JA, matrix->nnz, MPI_INT, 0, MPI_COMM_WORLD);

   /* ---------------------------------------------------------------------- */
   /*  Partitioning of the matrix among the processors                       */
   /* ---------------------------------------------------------------------- */

   mask = (int *)primme_calloc(matrix->n, sizeof(int), "mask");
   map = (int *)primme_calloc(numProcs+1, sizeof(int), "map");
   fg2or = (int *)primme_calloc(matrix->n, sizeof(int), "fg2or");
   or2fg = (int *)primme_calloc(matrix->n, sizeof(int), "or2fg");
      
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
   mLocal = matrix->m / numProcs;
   modulo = matrix->m % numProcs;
   rangeStart = 0;
   for (i=0; i<numProcs; i++) {
      rangeEnd = rangeStart + nLocal;
      if (i < modulo) rangeEnd = rangeEnd + 1;
      for (j = rangeStart; j< rangeEnd; j++) mask[j] = i;
      rangeStart = rangeEnd;
   }
   /* * * * * * * * * * * * * * * * * * * * * * * */

   generatePermutations(matrix->n, numProcs, mask, or2fg, fg2or, map);
   rangeStart = map[procID];
   rangeEnd = map[procID+1]-1;
   *n = matrix->n;
   *numProcs_ = numProcs;
   *procID_ = procID;
   *nLocal_ = rangeEnd - rangeStart+1;
   *pmatrix = csrToParaSails(procID, map, fg2or, or2fg, 
                             matrix->IA, matrix->JA, matrix->AElts, comm);

/* ------------------------------------------------------------------------- */
/*  Set up preconditioner if needed. For parallel programs the only choice   */
/*  in driver.PrecChoice:  (driver.PrecChoice is read in read_driver_params) */
/*     choice = 4    Parallel ParaSails preconditioners                      */
/*                   with parameters (level, threshold, filter, isymm)       */
/*                   as read in read_driver_params().                        */
/* ------------------------------------------------------------------------- */

   if (pfactor) {
      *pfactor = generate_precond(matrix, shift, matrix->n, procID, map, fg2or, or2fg,
         rangeStart, rangeEnd, isymm, level, threshold, filter, comm);
   }
   else {
      // Free A as it is not further needed
      free(matrix->AElts); free(matrix->IA); free(matrix->JA);
   }

   free(mask); free(map); free(fg2or); free(or2fg); free(matrix);
   return 0;
}

/******************************************************************************
 * void generatePermutations() 
 * Given a proc array :  proc[i] = processor # where row i lies
 * it generates all other needed permutation arrays for processing in 
 * Parasails.
 *
******************************************************************************/
static void generatePermutations(int n, int nParts, int *proc, int *perm,
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

/******************************************************************************
 * Convert CSR matrix format to Parasails matrix format 
 *
******************************************************************************/
static Matrix* csrToParaSails(int procID, int *map, int *fg2or, int *or2fg, int *IA,
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
static ParaSails* generate_precond(CSRMatrix *matrix, double shift, int n, int procID,
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
   shiftCSRMatrix(-shift, matrix->n, matrix->IA, matrix->JA, matrix->AElts);

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
 * Applies the PARALLEL matrix vector mulitplication of a block of vectors.
 * Because a block function is not available, we call blockSize times
 * the Parasails function MatrixMatvec()
 *
******************************************************************************/
void ParaSailsMatrixMatvec(void *x, void *y, int *blockSize, 
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
void ApplyPrecParaSails(void *x, void *y, int *blockSize,
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
