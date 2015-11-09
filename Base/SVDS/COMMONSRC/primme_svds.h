
#ifndef PRIMME_SVDS_H
#define PRIMME_SVDS_H

#include "primme.h"

typedef enum {
    primme_svds_largest,
    primme_svds_smallest
}primme_svds_target;

typedef enum {
    primme_svds_hybrid,
    primme_svds_normalequations,
    primme_svds_augmented
} primme_svds_svdsmethod;

typedef struct primme_svds_stats {                                                                                 
    int numOuterIterations;
    int numRestarts;
    int numMatvecs;
    int numPreconds;
    double elapsedTime;
} primme_svds_stats;

typedef struct primme_svds_params {
    /* Specify the size of the rectangular matrix A */
    int m; /*row size*/ 
    int n; /*colomn size*/

    /* The user must input at least the following two arguments */
    void (*matrixMatvec) 
        (void *x,  void *y, int *blockSize, struct primme_svds_params *primme_svds, const char *transpose);
    void (*applyPreconditioner) 
        (void *x,  void *y, int *blockSize, struct primme_svds_params *primme_svds, const char *transpose);

    /* The user may provide matrix-vector and apply 
       preconditioning for normal equation matrix ATA or 
       augmented matrix B directly, not implement yet.  */
    void (*matrixATA_Matvec) 
        (void *x,  void *y, int *blockSize, struct primme_svds_params *primme_svds);
    void (*applyATA_Preconditioner) 
        (void *x,  void *y, int *blockSize, struct primme_svds_params *primme_svds);
    void (*matrixAugmented_Matvec) 
        (void *x,  void *y, int *blockSize, struct primme_svds_params *primme_svds);
    void (*applyAugmented_Preconditioner) 
        (void *x,  void *y, int *blockSize, struct primme_svds_params *primme_svds);

    /* Input for the following is only required for parallel programs */
    int numProcs;
    int procID;
    int mLocal;
    int nLocal;
    void *commInfo;
    void (*globalSumDouble)
        (void *sendBuf, void *recvBuf, int *count, struct primme_params *primme );

    /* Though primme_svds_initialize will assign defaults, most users will set these */
    int numSvals;
    primme_svds_target target;
    int numTargetShifts;    /* For primme_svds_augmented meethod, user has to */ 
    double *targetShifts;   /* make sure  at least one shift must also be set */
    primme_preset_method eigsMethod_stage1;
    primme_preset_method eigsMethod_stage2;
    primme_svds_svdsmethod svdsMethod;

    /* These pointers are not for users but for d/zprimme_svds function */
    void *primme_svds_info;
    void *realWork;
    
    /* These pointers may be used for users to provide matrix/preconditioner */
    void *matrix;
    void *preconditioner;
    
    /* The following will be given default values depending on the method */
    double aNorm;
    double eps;

    int precondition;
    int initSize;
    int maxBasisSize;
    int maxBlockSize;
    int maxMatvecs;
    int iseed[4];
    int printLevel;
    FILE *outputFile;
    struct primme_svds_stats stats;
    struct primme_params primme;

} primme_svds_params;

int dprimme_svds(double *evals, double *evecs, double *resNorms,
         primme_svds_params *primme_svds);
int zprimme_svds(double *evals, Complex_Z *evecs, double *resNorms,
         primme_svds_params *primme_svds);
void primme_svds_initialize(primme_svds_params *primme_svds);
void primme_svds_display_params(primme_svds_params primme_svds);
void primme_svds_Free(primme_svds_params *primme_svds);
void MatrixATA_Matvec(void *x, void *y, int *blockSize, primme_params *primme);
void MatrixB_Matvec(void *x, void *y, int *blockSize, primme_params *primme);
void MatrixATA_Precond(void *x, void *y, int *blockSize, primme_params *primme);
void MatrixB_Precond(void *x, void *y, int *blockSize, primme_params *primme);

#endif /* PRIMME_SVDS_H */  
