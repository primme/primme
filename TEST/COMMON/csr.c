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
 * File: csr.c
 * 
 * Purpose - Functions to read MatrixMarket format matrices and other CSR
 *           auxiliary functions.
 * 
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "mmio.h"
#include "primme.h"
#include "csr.h"

static int readfullMTX(const char *mtfile, PRIMME_NUM **A, int **JA, int **IA, int *m, int *n, int *nnz);
#ifndef USE_DOUBLECOMPLEX
static int readUpperMTX(const char *mtfile, double **A, int **JA, int **IA, int *n, int *nnz);
int ssrcsr(int *job, int *value2, int *nrow, double *a, int *ja, int *ia, 
   int *nzmax, double *ao, int *jao, int *iao, int *indu, int *iwk, int *ierr);
#endif

int readMatrixNative(const char* matrixFileName, CSRMatrix **matrix_, double *fnorm) {
   int ret;
   CSRMatrix *matrix;

   matrix = (CSRMatrix*)primme_calloc(1, sizeof(CSRMatrix), "CSRMatrix");
   if (!strcmp("mtx", &matrixFileName[strlen(matrixFileName)-3])) {  
      /* coordinate format storing both lower and upper triangular parts */
      ret = readfullMTX(matrixFileName, &matrix->AElts, &matrix->JA, 
         &matrix->IA, &matrix->m, &matrix->n, &matrix->nnz);
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   else if (matrixFileName[strlen(matrixFileName)-1] == 'U') {
      /* coordinate format storing only upper triangular part */
#ifndef USE_DOUBLECOMPLEX
      ret = readUpperMTX(matrixFileName, &matrix->AElts, &matrix->JA,
         &matrix->IA, &matrix->n, &matrix->nnz);
#else
      /* TODO: support this in complex arithmetic */
      ret = -1;
#endif
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   else {  
      /* Harwell Boeing format NOT IMPLEMENTED */
      /* ret = readmt() */
      ret = -1;
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   *matrix_ = matrix;
   if (fnorm)
      *fnorm = frobeniusNorm(matrix);

   return 0;
}

static int *my_comp_ctx[2];
static int my_comp(const void *a, const void *b)
{
   const int ia = *(int*)a, ib = *(int*)b;
   int **p = my_comp_ctx;
   return p[0][ia] != p[0][ib] ? p[0][ia] - p[0][ib] : p[1][ia] - p[1][ib];
}

static int readfullMTX(const char *mtfile, PRIMME_NUM **AA, int **JA, int **IA, int *m, int *n, int *nnz) { 
   int i,j, k, nzmax;
   int *I, *J, *perm;
   PRIMME_NUM *A;
   double re, im;
   FILE *matrixFile;
   MM_typecode type;

   matrixFile = fopen(mtfile, "r");
   if (matrixFile == NULL) {
      return(-1);  
   }

   /* first read to set matrix kind and size */
   if(mm_read_banner(matrixFile, &type) != 0) return -1;
   if (!mm_is_valid(type) || !mm_is_sparse(type) || mm_is_skew(type)
#ifndef USE_DOUBLECOMPLEX
       || mm_is_complex(type) || mm_is_hermitian(type) || !(mm_is_real(type))
#endif
      ) {
      fprintf(stderr, "Matrix format '%s' not supported!", mm_typecode_to_str(type)); 
      return -1;
   }

   if (mm_read_mtx_crd_size(matrixFile, m, n, nnz) != 0) return -1;

   nzmax = *nnz;
   if (mm_is_symmetric(type)) nzmax *= 2;
   A = (PRIMME_NUM *)primme_calloc(nzmax, sizeof(PRIMME_NUM), "A");
   J = (int *)primme_calloc(nzmax, sizeof(int), "J");
   I = (int *)primme_calloc(nzmax, sizeof(int), "I");

   /* Read matrix in COO */
   for (k=0, i=0; k<*nnz; k++, i++) {
      if (mm_read_mtx_crd_entry(matrixFile, &I[i], &J[i], &re, &im, type)!=0) return -1;
      if (mm_is_pattern(type)) A[i] = 1;
      else if (mm_is_real(type)) A[i] = re;
      else A[i] = re + IMAGINARY*im;
      if ((mm_is_symmetric(type) || mm_is_hermitian(type)) && I[i] != J[i]) {
         I[i+1] = J[i]; J[i+1] = I[i]; A[i+1] = CONJ(A[i]); i++;
      }
   }
   nzmax = *nnz = i;

   /* Sort COO by columns */
   perm = (int *)primme_calloc(nzmax, sizeof(int), "perm");
   for (i=0; i<nzmax; i++) perm[i] = i;
   my_comp_ctx[0] = I;
   my_comp_ctx[1] = J;
   qsort(perm, nzmax, sizeof(int), my_comp);

   /* Collapse columns */
   *IA = (int *)primme_calloc(*m+1, sizeof(int), "IA");
   (*IA)[0] = 1;
   for (i=0, j=1; i<nzmax; i++)
      while (j < I[perm[i]]) (*IA)[j++] = i+1;
   while (j <= *m) (*IA)[j++] = nzmax+1;

   /* Copy rows and values sorted */
   *JA = I;
   for (i=0; i<nzmax; i++) (*JA)[i] = J[perm[i]];
   free(J);
   *AA = (PRIMME_NUM *)primme_calloc(nzmax, sizeof(PRIMME_NUM), "AA");
   for (i=0; i<nzmax; i++) (*AA)[i] = A[perm[i]];
   free(A);
   free(perm);

   fclose(matrixFile);

   return 0;
}

#ifndef USE_DOUBLECOMPLEX
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
   if (fscanf(matrixFile, "%d %d\n", n, nnz) != 2) return -1;
   fprintf(stderr, "%d %d\n", *n, *nnz);

   nzmax = 2*(*nnz) - *n;
   *A = (double *)primme_calloc(nzmax, sizeof(double), "A");
   *JA =   (int *)primme_calloc(nzmax, sizeof(int), "JA");
   *IA = (int *)primme_calloc(*n+1, sizeof(int), "IA");

   iwk1 = (int *)primme_calloc(*n+1, sizeof(int), "iwk1");
   iwk2 = (int *)primme_calloc(*n+1, sizeof(int), "iwk2");

   for (k=1; k <= *nnz; k++) {
      int tja; double ta;
      if (fscanf(matrixFile, "%d %d %lf\n", &row, &tja, &ta) != 3) return -1;
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
#endif

/******************************************************************************
 * Computed the Frobenius norm of a CSR matrix 
 *
 *         ||A||_frob = sqrt( \sum_{i,j} A_ij^2 )
 *
******************************************************************************/
double frobeniusNorm(const CSRMatrix *matrix) {

   int i, j;
   double fnorm;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   fnorm = 0.0L;

   for (i=0; i < matrix->m; i++) {
      for (j=matrix->IA[i]; j <= matrix->IA[i+1]-1; j++) {
         fnorm = fnorm + REAL_PART(CONJ(matrix->AElts[j-1])*matrix->AElts[j-1]);
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
void shiftCSRMatrix(double shift, CSRMatrix *matrix) {

   int i, j, n;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   for (i=0, n=min(matrix->m, matrix->n); i < n; i++) {
      for (j=matrix->IA[i]; j <= matrix->IA[i+1]-1; j++) {

         if (matrix->JA[j-1]-1 == i) {
            matrix->AElts[j-1] = matrix->AElts[j-1] + shift;
         }
      }
   }

}


