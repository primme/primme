#include <stdio.h>
#include <stdlib.h>
#include "shared_utils.h"

int ssrcsr(int *job, int *value2, int *nrow, double *a, int *ja, int *ia,
   int *nzmax, double *ao, int *jao, int *iao, int *indu, int *iwk, int *ierr);

int readfullMTX(char *mtfile, double **A, int **JA, int **IA, int *n, int *nnz) { 

   int i, k, nzmax;
   int row, nextRow;
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
            sscanf(ident, "%d %d %d\n",n, n, nnz);
            fprintf(stderr,"N=%d NNZ=%d \n",*n,*nnz);
            break;
         }
      }
   }

   i = 0;
   nextRow = 0;

   nzmax = 2*(*nnz) - *n;
   *A = (double *)primme_calloc(nzmax, sizeof(double), "A");
   *JA =   (int *)primme_calloc(nzmax, sizeof(int), "JA");
   *IA = (int *)primme_calloc(*n+1, sizeof(int), "IA");

   for (k=1; k <= *nnz; k++) {
      int tja; double ta;
      fscanf(matrixFile, "%d %d %lf\n", &tja, &row, &ta);
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


   return(0);
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
