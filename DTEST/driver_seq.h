#ifndef DRIVER_H
#define DRIVER_H

#include "shared_utils.h"
#include "primme.h"
#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#define DOUBLE_EQ(a,b) (((a)==0. && (b)==0.) || fabs(((a)-(b))/(b)) < 1e-15)

typedef struct {
   int *JA;
   int *IA;
   double *AElts;
} CSRMatrix; 

FILE *output;

void MatrixMatvec(void *x, void *y, int *blockSize, primme_params *primme);
void Apply_Diagonal_Shifted_Prec(void *x, void *y, int *blockSize, 
                            primme_params *primme);
void Apply_Inv_Diagonal_Prec(void *x, void *y, int *blockSize, 
                            primme_params *primme);
void Apply_ILUT_Prec(void *x, void *y, int *blockSize, 
                            primme_params *primme);
void generate_Inv_Diagonal_Prec(int n, double shift, 
   int *IA, int *JA, double *AElts, int *PIA, int *PJA, double *PElts);
void generate_Diagonal_Prec(int n, int *IA, int *JA, double *AElts, 
   int *PIA, int *PJA, double *PElts);
double frobeniusNorm(int n, int *IA, double *AElts);
void shiftCSRMatrix(double shift, int n, int *IA, int *JA, double *AElts);

int create_preconditioner(CSRMatrix matrix, CSRMatrix *Factors, 
   void (**precond_function)(void *, void *, int *, primme_params *),
     int n, int nnz, driver_params driver);

int check_solution(const char *checkXFileName, primme_params *primme, double *evals,
                   double *evecs, double *rnorms, CSRMatrix *matrix);
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
void amux_(int *n, double *x, double *y, double *AElts, int *JA, int *IA);
void ilut_(int *n, double *AElts, int *JA, int *IA, int *lfil, double *tol,
   double *PAElts, int *PJA, int *PIA, int *lenFactors, double *W1, double *W2,
   int *iW1, int *iW2, int *iW3, int *ierr); 
void lusol0_(int *n, double *x, double *y, double *AElts, int *JA, int *IA);
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* DRIVER_H */
