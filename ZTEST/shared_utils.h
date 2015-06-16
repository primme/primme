#ifndef DRIVER_SHARED_H
#define DRIVER_SHARED_H

#include "primme.h"

#define TRUE  1
#define FALSE 0
#define ZERO 0.0L

typedef struct driver_params {

   char outputFileName[512];
   char partId[256];
   char testId[256];
   char partDir[1024];
   char matrixFileName[1024];

   int weightedPart;

   // Preconditioning paramaters for various preconditioners
   int PrecChoice;       // 0, 1, 2, or 3 for the moment (see driver)
   int isymm;
   int level;
   double threshold;
   double filter;
   double shift;
   
} driver_params;


int read_solver_params(char *configFileName, char *outputFileName,
		primme_params *primme, primme_preset_method *method);
int read_driver_params(char *configFileName, driver_params *driver);
int readfullMTX(char *mtfile, Complex_Z **A, int **JA, int **IA, int *n, 
		int *nnz);

#endif /* DRIVER_SHARED_H */
