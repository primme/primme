#include <stdio.h>
#include <stdlib.h>
#include "shared_utils.h"

int readfullMTX(char *mtfile, Complex_Z **A, int **JA, int **IA, int *n, 
		int *nnz) { 

   int i,j, k, row, nzmax, tja;
   int firstFast;
   Complex_Z ta;
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

   nzmax = *nnz;
   *A = (Complex_Z *)primme_calloc(nzmax, sizeof(Complex_Z), "A");
   *JA =   (int *)primme_calloc(nzmax, sizeof(int), "JA");
   *IA = (int *)primme_calloc(*n+1, sizeof(int), "IA");

   /* Check in the input file which column runs fast: first or second */
   firstFast = 1;
   fscanf(matrixFile, "%d %d %lf %lf\n", &tja, &row, &ta.r, &ta.i);
   i = row;
   j = tja;
   for (k=0;k<*nnz;k++) {
      fscanf(matrixFile, "%d %d %lf %lf\n", &tja, &row, &ta.r, &ta.i);
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

   i = 0;
   (*IA)[i] = 1;
   for (k=0; k < *nnz; k++) {
      if (firstFast) 
         fscanf(matrixFile, "%d %d %lf %lf\n", &tja, &row, &ta.r, &ta.i);
      else
         fscanf(matrixFile, "%d %d %lf %lf\n", &row, &tja, &ta.r, &ta.i);
      (*JA)[k]=tja;
      (*A)[k].r = ta.r;
      (*A)[k].i = ta.i;
      if (i != row-1) {
	 for (j=i+1;j<row;j++) (*IA)[j]=k+1;
	 i=row-1;
      }
   }

   for (i=row;i<=*n;i++) (*IA)[i]=(*IA)[0] + *nnz;
   fclose(matrixFile);
   return(0);
}
