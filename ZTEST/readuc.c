#include <stdio.h>
#include <stdlib.h>
#include "shared_utils.h"

int readfullMTX(char *mtfile, Complex_Z **A, int **JA, int **IA, int *n, 
		int *nnz) { 

   int i, k, nzmax;
   int row, nextRow, tja; 
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

   i = 0;
   nextRow = 0;

   nzmax = 2*(*nnz) - *n;
   *A = (Complex_Z *)primme_calloc(nzmax, sizeof(Complex_Z), "A");
   *JA =   (int *)primme_calloc(nzmax, sizeof(int), "JA");
   *IA = (int *)primme_calloc(*n+1, sizeof(int), "IA");

   for (k=1; k <= *nnz; k++) {
      fscanf(matrixFile, "%d %d %lf %lf\n", &tja, &row, &ta.r, &ta.i);
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
