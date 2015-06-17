#ifndef DRIVER_H
#define DRIVER_H

#include "ParaSails.h"
#include "shared_utils.h"
#include "primme.h"

typedef struct {
   int *JA;
   int *IA;
   double *AElts;
} CSRMatrix;

FILE *output;

void broadCast(primme_params *primme, primme_preset_method *method, 
   int procID, MPI_Comm comm);
ParaSails* generate_precond(CSRMatrix *matrix, double shift, int n, int procID,
   int *map, int *fg2or, int *or2fg, int rangeStart, int rangeEnd, int isymm, 
   int level, double threshold, double filter, MPI_Comm comm);
Matrix* csrToParaSails(int procID, int *map, int *fg2or, int *or2fg, int *IA, 
   int *JA, double *AElts, MPI_Comm comm);
void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, 
   primme_params *primme);
void par_MatrixMatvec(void *x, void *y, int *blockSize, 
   primme_params *primme);
void par_ApplyParasailsPrec(void *x, void *y, int *blockSize,
   primme_params *primme);
void generatePermutations(int n, int nParts, int *proc, int *perm,
   int *iperm, int *map);
double frobeniusNorm(int n, int *IA, double *AElts);
void shiftCSRMatrix(double shift, int n, int *IA, int *JA, double *AElts);

#endif /* DRIVER_H */
