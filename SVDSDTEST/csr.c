/*  PRIMME PReconditioned Iterative MultiMethod Eigensolver
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
 */

/* Required by qsort_r */
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "primme.h"
#include "csr.h"
#include "mmio.h"

static int readfullMTX(const char *mtfile, double **A, int **JA, int **IA, int *n, int *nnz);
static int readUpperMTX(const char *mtfile, double **A, int **JA, int **IA, int *n, int *nnz);
int ssrcsr(int *job, int *value2, int *nrow, double *a, int *ja, int *ia, 
   int *nzmax, double *ao, int *jao, int *iao, int *indu, int *iwk, int *ierr);

int readMatrixNative(const char* matrixFileName, CSRMatrix **matrix_, double *fnorm) {
   int ret;
   CSRMatrix *matrix;

   matrix = (CSRMatrix*)primme_calloc(1, sizeof(CSRMatrix), "CSRMatrix");
   if (!strcmp("mtx", &matrixFileName[strlen(matrixFileName)-3])) {  
      // coordinate format storing both lower and upper triangular parts
      ret = readfullMTX(matrixFileName, &matrix->AElts, &matrix->JA, 
         &matrix->IA, &matrix->n, &matrix->nnz);
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   else if (matrixFileName[strlen(matrixFileName)-1] == 'U') {
      // coordinate format storing only upper triangular part
      ret = readUpperMTX(matrixFileName, &matrix->AElts, &matrix->JA,
         &matrix->IA, &matrix->n, &matrix->nnz);
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   else {  
      //Harwell Boeing format NOT IMPLEMENTED
      //ret = readmt()
      ret = -1;
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   *matrix_ = matrix;
   if (fnorm)
      *fnorm = frobeniusNorm(matrix->n, matrix->IA, matrix->AElts);

   return 0;
}

static int my_comp(const void *a, const void *b, void *ctx) {
   const int ia = *(int*)a, ib = *(int*)b;
   int **p = (int **)ctx;
   return p[0][ia] != p[0][ib] ? p[0][ia] - p[0][ib] : p[1][ia] - p[1][ib];
}

static int readfullMTX(const char *mtfile, double **AA, int **JA, int **IA, int *n, int *nnz) { 
   int i,j, k, nzmax, m;
   int *I, *J, *perm;
   double *A;
   double im;
   FILE *matrixFile;
   MM_typecode type;

   matrixFile = fopen(mtfile, "r");
   if (matrixFile == NULL) {
      return(-1);  
   }

   /* first read to set matrix kind and size */
   if(mm_read_banner(matrixFile, &type) != 0) return -1;
   if (!mm_is_valid(type) || !mm_is_sparse(type) || mm_is_complex(type) ||
       mm_is_hermitian(type) || mm_is_skew(type) ||
       !(mm_is_real(type) || mm_is_pattern(type))) {
      fprintf(stderr, "Matrix format '%s' not supported!", mm_typecode_to_str(type)); 
      return -1;
   }

   if (mm_read_mtx_crd_size(matrixFile, &m, n, nnz) != 0) return -1;
   if (m != *n) return -1;

   nzmax = *nnz;
   if (mm_is_symmetric(type)) nzmax *= 2;
   A = (double *)primme_calloc(nzmax, sizeof(double), "A");
   J = (int *)primme_calloc(nzmax, sizeof(int), "J");
   I = (int *)primme_calloc(nzmax, sizeof(int), "I");

   /* Read matrix in COO */
   for (k=0, i=0; k<*nnz; k++, i++) {
      if (mm_read_mtx_crd_entry(matrixFile, &I[i], &J[i], &A[i], &im, type)!=0) return -1;
      if (mm_is_symmetric(type) && I[i] != J[i]) {
         I[i+1] = J[i]; J[i+1] = I[i]; A[i+1] = A[i]; i++;
      }
   }
   nzmax = *nnz = i;

   /* Sort COO by columns */
   perm = (int *)primme_calloc(nzmax, sizeof(int), "perm");
   for (i=0; i<nzmax; i++) perm[i] = i;
   {
      int *ctx[2] = {J, I};
      qsort_r(perm, nzmax, sizeof(int), my_comp, ctx);
   }

   /* Collapse columns */
   *IA = (int *)primme_calloc(*n+1, sizeof(int), "IA");
   (*IA)[0] = 1;
   for (i=0, j=1; i<nzmax; i++)
      while (j < J[perm[i]]) (*IA)[j++] = i+1;
   while (j <= *n) (*IA)[j++] = nzmax+1;

   /* Copy rows and values sorted */
   *JA = J;
   for (i=0; i<nzmax; i++) (*JA)[i] = I[perm[i]];
   free(I);
   *AA = (double *)primme_calloc(nzmax, sizeof(double), "AA");
   for (i=0; i<nzmax; i++) (*AA)[i] = A[perm[i]];
   free(A);

   fclose(matrixFile);

   return(0);
}

static int readUpperMTX(const char *mtfile, double **A, int **JA, int **IA, int *n, int *nnz) { 
   int i, k, nzmax;
   int job, value2;
   int row, nextRow;
   int ierror;
   int *iwk1, *iwk2;
   FILE *matrixFile;

   matrixFile = fopen(mtfile, "r");

   if (matrixFile == NULL) {
      return(-1);  
   }

   i = 0;
   nextRow = 0;
   fscanf(matrixFile, "%d %d\n", n, nnz);
   fprintf(stderr, "%d %d\n", *n, *nnz);

   nzmax = 2*(*nnz) - *n;
   *A = (double *)primme_calloc(nzmax, sizeof(double), "A");
   *JA =   (int *)primme_calloc(nzmax, sizeof(int), "JA");
   *IA = (int *)primme_calloc(*n+1, sizeof(int), "IA");

   iwk1 = (int *)primme_calloc(*n+1, sizeof(int), "iwk1");
   iwk2 = (int *)primme_calloc(*n+1, sizeof(int), "iwk2");

   for (k=1; k <= *nnz; k++) {
      int tja; double ta;
      fscanf(matrixFile, "%d %d %lf\n", &row, &tja, &ta);
      (*JA)[k-1]=tja;
      (*A)[k-1] = ta;
      if (i != row) {
         i = row;
         nextRow = nextRow + 1;
         (*IA)[nextRow-1] = k;
      }
   }

   (*IA)[*n] = (*IA)[0] + *nnz;
   fclose(matrixFile);

   job = 3;
   value2 = 1;

   ssrcsr(&job, &value2, n, *A, *JA, *IA, &nzmax, *A, *JA, *IA, iwk1, iwk2,
      &ierror);
   *nnz = 2*(*nnz) - *n;

   free(iwk1);
   free(iwk2);

   return(0);
}

/******************************************************************************
 * Computed the Frobenius norm of a CSR matrix 
 *
 *         ||A||_frob = sqrt( \sum_{i,j} A_ij^2 )
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
 *         A = A + shift I
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


