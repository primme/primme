/* Required by qsort_r */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include "shared_utils.h"
#include "mmio.h"

int ssrcsr(int *job, int *value2, int *nrow, double *a, int *ja, int *ia,
   int *nzmax, double *ao, int *jao, int *iao, int *indu, int *iwk, int *ierr);

int readfullMTX(char *mtfile, double **A, int **JA, int **IA, int *m, int *n, int *nnz) { 

   int i,j, k, row, nzmax, tja;
   int firstFast;
   double ta;
   FILE *matrixFile;
   char ident[128];

   matrixFile = fopen(mtfile, "r");

   if (matrixFile == NULL) {
      return(-1);  
   }

   while (1) {
      if (NULL == fgets(ident, 128, matrixFile)) {
              return(-1);
      }
      else {
         fprintf(stderr,"%s",ident);
         if (ident[0] != '%') {
            sscanf(ident, "%d %d %d\n",m, n, nnz);
            fprintf(stderr,"M=%d N=%d NNZ=%d \n",*m,*n,*nnz);
            break;
         }
      }
   }

   nzmax = *nnz;
   *A = (double *)primme_calloc(nzmax, sizeof(double), "A");
   *JA =   (int *)primme_calloc(nzmax, sizeof(int), "JA");
   *IA = (int *)primme_calloc(*n+1, sizeof(int), "IA");

   /* Check in the input file which column runs fast: first or second */
   firstFast = 1;
   fscanf(matrixFile, "%d %d %lf\n", &tja, &row, &ta);
   i = row;
   j = tja;
   for (k=0;k<*nnz;k++) {
      fscanf(matrixFile, "%d %d %lf\n", &tja, &row, &ta);
      if (i == row && j != tja) {
         firstFast = 1;
         break;
      }
      if (i != row && j == tja) {
         firstFast = 0;
         break;
      }
   }
   /* Rewind the file and reread info */
   rewind(matrixFile);
   while (1) {
      if (NULL == fgets(ident, 128, matrixFile)) {
              return(-1);
      }
      else {
         if (ident[0] != '%') {
            sscanf(ident, "%d %d %d\n",n, n, nnz);
            break;
         }
      }
   }


   /* Read in the CSR format */
   i = 0;
   (*IA)[i] = 1;
   for (k=0; k < *nnz; k++) {
      if (firstFast) 
         fscanf(matrixFile, "%d %d %lf\n", &tja, &row, &ta);
      else
         fscanf(matrixFile, "%d %d %lf\n", &row, &tja, &ta);
      (*JA)[k]=tja;
      (*A)[k] = ta;
      if (i != row-1) {
	 for (j=i+1;j<row;j++) (*IA)[j]=k+1;
	 i=row-1;
      }
   }

   for (i=row;i<=*n;i++) (*IA)[i]=(*IA)[0] + *nnz;
   fclose(matrixFile);

   return(0);
}

static int my_comp(void *ctx, const void *a, const void *b) {
   const int ia = *(int*)a, ib = *(int*)b;
   int **p = (int **)ctx;
   return p[0][ia] != p[0][ib] ? p[0][ia] - p[0][ib] : p[1][ia] - p[1][ib];
}

int readfullMTXR(char *mtfile, double **AA, int **JA, int **IA, int *m, int *n, int *nnz) { 
   int i,j, k, nzmax;
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

   if (mm_read_mtx_crd_size(matrixFile, m, n, nnz) != 0) return -1;

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
      int *ctx[2] = {I, J};
      qsort_r(perm, nzmax, sizeof(int), ctx, my_comp);
   }

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
   *AA = (double *)primme_calloc(nzmax, sizeof(double), "AA");
   for (i=0; i<nzmax; i++) (*AA)[i] = A[perm[i]];
   free(A);

   fclose(matrixFile);

   return 0;
}

int readUpperMTX(char *mtfile, double **A, int **JA, int **IA, int *n, int *nnz) { 

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
