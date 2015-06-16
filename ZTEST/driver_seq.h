#ifndef DRIVER_H
#define DRIVER_H

#include "shared_utils.h"
#include "primme.h"
#define max(a, b) (a > b ? a : b)

typedef struct {
   int *JA;
   int *IA;
   Complex_Z *AElts;
} CSRMatrix; 

FILE *output;

void MatrixMatvec(void *x, void *y, int *blockSize, primme_params *primme);
void Apply_Diagonal_Shifted_Prec(void *x, void *y, int *blockSize, 
		            primme_params *primme);
void Apply_Inv_Diagonal_Prec(void *x, void *y, int *blockSize, 
                            primme_params *primme);
void generate_Inv_Diagonal_Prec(int n, double shift, 
   int *IA, int *JA, Complex_Z *AElts, int *PIA, int *PJA, Complex_Z *PElts);
void generate_Diagonal_Prec(int n, int *IA, int *JA, Complex_Z *AElts, 
   int *PIA, int *PJA, Complex_Z *PElts);
double frobeniusNorm(int n, int *IA, Complex_Z *AElts);
void shiftCSRMatrix(double shift, int n, int *IA, int *JA, Complex_Z *AElts);

int create_preconditioner(CSRMatrix matrix, CSRMatrix *Factors, 
#ifdef Cplusplus     /* C++ has a stricter type checking */
   void (**precond_function)(void *, void *, int *, primme_params *),
#else
   void **precond_function,
#endif
   int n, int nnz, driver_params driver);

#ifdef Cplusplus
extern "C" {
#endif /* Cplusplus */
void zamux_(int *n, Complex_Z *x, Complex_Z *y, Complex_Z *AElts, int *JA, int *IA);
#ifdef Cplusplus
}
#endif /* Cplusplus */


#endif /* DRIVER_H */
