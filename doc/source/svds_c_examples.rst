.. _svdsSimple:

Simple C SVDS Example
---------------------

.. code:: c

    #include <stdlib.h>
    #include <stdio.h>
    #include <math.h>
    #include <assert.h>

    #include "primme.h"   /* header file for PRIMME SVDS too */ 

    #ifndef min
    #define min(A,B) ((A)<=(B)?(A):(B))
    #endif
    #ifndef max
    #define max(A,B) ((A)>=(B)?(A):(B))
    #endif
    int main (int argc, char *argv[]) {

    /* Solver arrays and parameters */
    double *svals;    /* Array with the computed singular values */
    double *rnorms;   /* Array with the computed residual norms */
    double *svecs;    /* Array with the computed singular vectors;
                        first right (v) vector starts in svecs[0],
                        second right (v) vector starts in svecs[primme_svd.n],
                        first left (u) vector starts in svecs[primme_svd.n*numSVals]...  */
    primme_svds_params primme_svds;
                    /* PRIMME SVDS configuration struct */

    /* Other miscellaneous items */
    int ret;
    int i;
    double mu = 1e-5;

    /* Set default values in PRIMME SVDS configuration struct */
    primme_svds_initialize(&primme_svds);

    /* Set problem matrix */
    primme_svds.matrixMatvec = LauchliMatrixMatvec;
    primme_svds.matrix = &mu;
                            /* Function that implements the matrix-vector products
                            A*x and A^t*x  */

    /* Set problem parameters */
    primme_svds.m = 500;
    primme_svds.n = 100; /* set problem dimension */
    primme_svds.numSvals = 4;   /* Number of wanted singular values */
    primme_svds.eps = 1e-6;     /* ||r|| <= eps * ||matrix|| */
    primme_svds.target = primme_svds_smallest;
                                /* Seeking for the largest singular values  */

    /* Set preconditioner (optional) */
    primme_svds.applyPreconditioner = LauchliApplyPreconditioner;

    /* Set method to solve the singular value problem and
    the underneath eigenvalue problem (optional) */
    primme_svds_set_method(primme_svds_default, PRIMME_DEFAULT_METHOD,
                            PRIMME_DEFAULT_METHOD, &primme_svds);
    /*  primme_svds_default: devs choice, now being hybrid, which first solve
        the normal equation and then the augmented problem.
        PRIMME_DEFAULT_METHOD devs choice of the solver at every stage. But other methods
        can be set such as DYNAMIC or PRIMME_LOBPCG_OrthoBasis_Window. */

    primme_svds.printLevel = 3;


    /* Set advanced parameters if you know what are you doing (optional) */
    /* Configuration for 1st stage */
    /*
    or you can do:
    primme_svds.primme.maxBasisSize = 14;
    primme_svds.primme.minRestartSize = 6;
    primme_svds.primme.maxBlockSize = 2;
    */
    /* Configuration for 2nd stage */
    /*
    primme_svds.primmeStage2.maxBasisSize = 30;
    primme_svds.primmeStage2.minRestartSize = 15;
    primme_svds.primmeStage2.maxBlockSize = 1;
    */

    /* Display PRIMME SVDS configuration struct (optional) */
    primme_svds_display_params(primme_svds);

    /* Allocate space for converged Ritz values and residual norms */
    svals = (double*)malloc(primme_svds.numSvals*sizeof(double));
    svecs = (double*)malloc((primme_svds.n+primme_svds.m)
        *primme_svds.numSvals*sizeof(double));
    rnorms = (double*)malloc(primme_svds.numSvals*sizeof(double));

    /* Call primme_svds  */
    ret = dprimme_svds(svals, svecs, rnorms, &primme_svds);

    void LauchliMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                            int *transpose, primme_svds_params *primme_svds, int *err) {
    
        int i;            /* vector index, from 0 to *blockSize-1 */
        int j;
        int min_m_n = min(primme_svds->m, primme_svds->n);
        double *xvec;     /* pointer to i-th input vector x */
        double *yvec;     /* pointer to i-th output vector y */
        double mu = *(double*)primme_svds->matrix;

        if (*transpose == 0) { /* Do y <- A * x */
            for (i=0; i<*blockSize; i++) { 
                xvec = (double *)x + (*ldx)*i;
                yvec = (double *)y + (*ldy)*i;
                yvec[0] = 0;
                for (j=0; j<primme_svds->n; j++) {
                    yvec[0] += xvec[j];
                }
                for (j=1; j<primme_svds->m; j++) {
                    yvec[j] = j-1<primme_svds->n ? xvec[j-1]*(1.0 - (1.0 - mu)*(j-1)/(min_m_n - 1)) : 0.0;
                }      
            }
        } else { /* Do y <- A^t * x */
            for (i=0; i<*blockSize; i++) {
                xvec = (double *)x + (*ldx)*i;
                yvec = (double *)y + (*ldy)*i;
                for (j=0; j<primme_svds->n; j++) {
                    yvec[j] = xvec[0];
                    if (j+1 < primme_svds->m) yvec[j] += xvec[j+1]*(1.0 - (1.0 - mu)*j/(min_m_n - 1));
                }
            }
        }
        *err = 0;
    }

    /* This performs Y = M^{-1} * X, where

    - X, input dense matrix of size primme_svds.n (or primme_svds.m or m+n) x blockSize;
    - Y, output dense matrix of size primme_svds.n (or primme_svds.m or m+n) x blockSize;
    - M, preconditioner for A^t*A (or A*A^t or [0 A^t; A 0]), where A is the Lauchli matrix.
    */

    void LauchliApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                                    int *mode, primme_svds_params *primme_svds, int *ierr) {
    
        int i;            /* vector index, from 0 to *blockSize-1*/
        int j;            /* row index */
        double *xvec;     /* pointer to i-th input vector x */
        double *yvec;     /* pointer to i-th output vector y */
        int modeAtA = primme_svds_op_AtA, modeAAt = primme_svds_op_AAt;
        double mu = *(double*)primme_svds->matrix;
        double  *aux;
        PRIMME_INT ldaux;
        int notrans = 0, trans = 1;
        int min_m_n = min(primme_svds->m, primme_svds->n);
            
        if (*mode == primme_svds_op_AtA) {
            /* Preconditioner for A^t*A, diag(A^t*A)^{-1} */
            for (i=0; i<*blockSize; i++) { 
                xvec = (double *)x + (*ldx)*i;
                yvec = (double *)y + (*ldy)*i;
                for (j=0; j<primme_svds->n; j++) {
                    double ei = j<primme_svds->m ? 1.0 - (1.0 - mu)*j/(min_m_n - 1) : 0.0;
                    yvec[j] = xvec[j]/(1.0 + ei*ei);
                }      
            }
        }
        else if (*mode == primme_svds_op_AAt) {
            /* Preconditioner for A*A^t, diag(A*A^t)^{-1} */
            for (i=0; i<*blockSize; i++) {
                xvec = (double *)x + (*ldx)*i;
                yvec = (double *)y + (*ldy)*i;
                yvec[0] = xvec[0]/(double)primme_svds->m;
                for (j=1; j<primme_svds->m; j++) {
                    double ei = j<primme_svds->n ? 1.0 - (1.0 - mu)*j/(min_m_n - 1) : 1.0;
                    yvec[j] = xvec[j]/ei/ei;
                }
            }
        }
        else if (*mode == primme_svds_op_augmented) {
            /* Preconditioner for [0 A^t; A 0],
                [diag(A^t*A) 0; 0 diag(A*A^t)]^{-1}*[0 A^t; A 0] */

            /* [y0; y1] <- [0 A^t; A 0] * [x0; x1] */
            ldaux = primme_svds->n+primme_svds->m;
            aux = (double*)malloc(sizeof(double)*(*blockSize)*ldaux);
            primme_svds->matrixMatvec(x, ldx, &aux[primme_svds->n], &ldaux, blockSize, &notrans, primme_svds, ierr);
            xvec = (double *)x + primme_svds->n;
            primme_svds->matrixMatvec(xvec, ldx, aux, &ldaux, blockSize, &trans, primme_svds, ierr);
            /* y0 <- preconditioner for A^t*A  * y0 */
            LauchliApplyPreconditioner(aux, &ldaux, y, ldy, blockSize, &modeAtA, primme_svds, ierr);
            /* y1 <- preconditioner for A*A^t  * y1 */
            yvec = (double *)y + primme_svds->n;
            LauchliApplyPreconditioner(&aux[primme_svds->n], &ldaux, yvec, ldy, blockSize, &modeAAt, primme_svds, ierr);
            free(aux);
        }
        *ierr = 0;
    }

.. _svdsParallel:

Parallel C SVDS Example
-----------------------

.. code:: c

    #include <stdlib.h>
    #include <stdio.h>
    #include <math.h>
    #include <assert.h>

    #include <petscpc.h>
    #include <petscmat.h>
    #include "primme.h"   /* header file for PRIMME SVDS too */ 

    #ifndef min
    #define min(A,B) ((A)<=(B)?(A):(B))
    #endif
    #ifndef max
    #define max(A,B) ((A)>=(B)?(A):(B))
    #endif


    PetscErrorCode generateLauchli(int m, int n, PetscReal mu, Mat *A);
    void PETScMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                            int *transpose, primme_svds_params *primme_svds, int *ierr);
    void ApplyPCPrecAHA(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                            int *transpose, primme_svds_params *primme_svds, int *ierr);
    void par_GlobalSum(void *sendBuf, void *recvBuf, int *count,
                            primme_svds_params *primme_svds, int *ierr);

    int main (int argc, char *argv[]) {

    /* Solver arrays and parameters */
    PetscReal *svals;    /* Array with the computed singular values */
    PetscReal *rnorms;   /* Array with the computed residual norms */
    PetscScalar *svecs;    /* Array with the computed singular vectors;
                            first right (v) vector starts in svecs[0],
                            second right (v) vector starts in svecs[primme_svd.n],
                            first left (u) vector starts in svecs[primme_svd.n*numSVals]...  */
    primme_svds_params primme_svds;
                        /* PRIMME SVDS configuration struct */

    /* Other miscellaneous items */
    int ret;
    int i;
    PetscReal mu = 1e-5;
    Mat A; /* problem matrix */
    Mat AHA;          /* auxiliary matrix for A^t*A */
    PC pc;            /* preconditioner */
    PetscErrorCode ierr;
    PetscInt m, n, mLocal, nLocal;
    MPI_Comm comm;

    PetscInitialize(&argc, &argv, NULL, NULL);


    /* Set default values in PRIMME SVDS configuration struct */
    primme_svds_initialize(&primme_svds);

    /* Set problem matrix */
    ierr = generateLauchli(500, 100, mu, &A); CHKERRQ(ierr);
    primme_svds.matrix = &A;
    primme_svds.matrixMatvec = PETScMatvec;
                            /* Function that implements the matrix-vector products
                                A*x and A^t*x  */
    
    /* Set problem parameters */
    ierr = MatGetSize(A, &m, &n); CHKERRQ(ierr);
    primme_svds.m = (PRIMME_INT)m;
    primme_svds.n = (PRIMME_INT)n; /* set problem dimension */
    primme_svds.numSvals = 4;   /* Number of wanted singular values */
    primme_svds.eps = 1e-6;     /* ||r|| <= eps * ||matrix|| */
    primme_svds.target = primme_svds_smallest;
                                /* Seeking for the largest singular values  */

    /* Set preconditioner (optional) */
    /* Build the Jacobi preconditioner of A^T*A, useful when m>=n */
    ierr = MatCreateNormal(A, &AHA); CHKERRQ(ierr);
    ierr = PCCreate(PETSC_COMM_WORLD, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);
    ierr = PCSetOperators(pc, AHA, AHA); CHKERRQ(ierr);
    ierr = PCSetFromOptions(pc); CHKERRQ(ierr);
    ierr = PCSetUp(pc); CHKERRQ(ierr);
    primme_svds.preconditioner = &pc;
    primme_svds.applyPreconditioner = ApplyPCPrecAHA;

    /* Set method to solve the singular value problem and
        the underneath eigenvalue problem (optional) */
    primme_svds_set_method(primme_svds_default, PRIMME_DEFAULT_METHOD,
                                PRIMME_DEFAULT_METHOD, &primme_svds);
    /*  primme_svds_default: devs choice, now being hybrid, which first solve
        the normal equation and then the augmented problem.
        PRIMME_DEFAULT_METHOD devs choice of the solver at every stage. But other methods
        can be set such as DYNAMIC or PRIMME_LOBPCG_OrthoBasis_Window. */

    primme_svds.printLevel = 3;

    /* Set parallel parameters */
    ierr = MatGetLocalSize(A, &mLocal, &nLocal); CHKERRQ(ierr);
    primme_svds.mLocal = (PRIMME_INT)mLocal;
    primme_svds.nLocal = (PRIMME_INT)nLocal;
    comm = PETSC_COMM_WORLD;
    primme_svds.commInfo = &comm;
    MPI_Comm_size(comm, &primme_svds.numProcs);
    MPI_Comm_rank(comm, &primme_svds.procID);
    primme_svds.globalSumReal = par_GlobalSum;


    /* Set advanced parameters if you know what are you doing (optional) */
    /* Configuration for 1st stage */
    /*
    or you can do:
    primme_svds.primme.maxBasisSize = 14;
    primme_svds.primme.minRestartSize = 6;
    primme_svds.primme.maxBlockSize = 2;
    */
    /* Configuration for 2nd stage */
    /*
    primme_svds.primmeStage2.maxBasisSize = 30;
    primme_svds.primmeStage2.minRestartSize = 15;
    primme_svds.primmeStage2.maxBlockSize = 1;
    */

        /* Display PRIMME SVDS configuration struct (optional) */
    if (primme_svds.procID == 0) /* Reports process with ID 0 */
        primme_svds_display_params(primme_svds);

    /* Allocate space for converged Ritz values and residual norms */
    svals = (PetscReal*)malloc(primme_svds.numSvals*sizeof(PetscReal));
    svecs = (PetscScalar*)malloc((primme_svds.n+primme_svds.m)
            *primme_svds.numSvals*sizeof(PetscScalar));
    rnorms = (PetscReal*)malloc(primme_svds.numSvals*sizeof(PetscReal));

    /* Call primme_svds  */
    #if defined(PETSC_USE_COMPLEX) && defined(PETSC_USE_REAL_SINGLE)
    ret = cprimme_svds(svals, svecs, rnorms, &primme_svds);
    #elif defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_REAL_SINGLE)
    ret = zprimme_svds(svals, svecs, rnorms, &primme_svds);
    #elif !defined(PETSC_USE_COMPLEX) && defined(PETSC_USE_REAL_SINGLE)
    ret = sprimme_svds(svals, svecs, rnorms, &primme_svds);
    #elif !defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_REAL_SINGLE)
    ret = dprimme_svds(svals, svecs, rnorms, &primme_svds);
    #endif

    void ApplyPCPrecAHA(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                            int *mode, primme_svds_params *primme_svds, int *err) {
    int i,j;
    Mat *matrix;
    PC *pc;
    Vec xvec, yvec;
    PetscScalar *x0 = (PetscScalar*)x, *y0 = (PetscScalar*)y;
    PetscErrorCode ierr;
    
    matrix = (Mat *)primme_svds->matrix;
    pc = (PC *)primme_svds->preconditioner;

    /* The preconditioner is only build for A^t*A; in the rest of cases y <= x */

    if (*mode == primme_svds_op_AtA) {
        ierr = MatCreateVecs(*matrix, &xvec, NULL); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
        ierr = MatCreateVecs(*matrix, &yvec, NULL); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
        for (i=0; i<*blockSize; i++) {
            ierr = VecPlaceArray(xvec, ((PetscScalar*)x) + primme_svds->nLocal*i); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
            ierr = VecPlaceArray(yvec, ((PetscScalar*)y) + primme_svds->nLocal*i); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
            ierr = PCApply(*pc, xvec, yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
            ierr = VecResetArray(xvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
            ierr = VecResetArray(yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
        }
        ierr = VecDestroy(&xvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
        ierr = VecDestroy(&yvec); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
    }
    else if (*mode == primme_svds_op_AAt) {
        for (i=0; i<*blockSize; i++)
            for (j=0; j<primme_svds->mLocal; j++)
                y0[(*ldy)*i+j] = x0[(*ldx)*i+j];
    }
    else if (*mode == primme_svds_op_augmented) {
        for (i=0; i<*blockSize; i++)
            for (j=0; j<primme_svds->mLocal+primme_svds->nLocal; j++)
                y0[(*ldy)*i+j] = x0[(*ldx)*i+j];
    }
    *err = 0;
    }

    void par_GlobalSum(void *sendBuf, void *recvBuf, int *count, 
                            primme_svds_params *primme_svds, int *ierr) {
    MPI_Comm communicator = *(MPI_Comm *) primme_svds->commInfo;

    if (sendBuf == recvBuf) {
        *ierr = MPI_Allreduce(MPI_IN_PLACE, recvBuf, *count, MPIU_REAL, MPI_SUM, communicator) != MPI_SUCCESS;
    } else {
        *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPIU_REAL, MPI_SUM, communicator) != MPI_SUCCESS;
    }
    }

.. _cppSvds:

C++ SVDS Example
----------------

.. code:: cpp

    #include <stdlib.h>
    #include <stdio.h>
    #include <math.h>
    #include <assert.h>
    #include <complex>
    #include "primme.h"   /* header file for PRIMME SVDS too */ 
    
    #ifndef min
    #define min(A,B) ((A)<=(B)?(A):(B))
    #endif
    #ifndef max
    #define max(A,B) ((A)>=(B)?(A):(B))
    #endif
    
    
    void LauchliMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                             int *transpose, primme_svds_params *primme_svds, int *ierr);
    void LauchliApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
                                    int *mode, primme_svds_params *primme_svds, int *ierr);
    void LauchliAugmentedMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
    
    int main (int argc, char *argv[]) {
    
       /* Solver arrays and parameters */
       double *svals;    /* Array with the computed singular values */
       double *rnorms;   /* Array with the computed residual norms */
       std::complex<double> *svecs;    /* Array with the computed singular vectors;
                            first right (v) vector starts in svecs[0],
                            second right (v) vector starts in svecs[primme_svd.n],
                            first left (u) vector starts in svecs[primme_svd.n*numSVals]...  */
       primme_svds_params primme_svds;
                         /* PRIMME SVDS configuration struct */
    
       /* Other miscellaneous items */
       int ret;
       int i;
       double mu = 1e-5;
    
       /* Set default values in PRIMME SVDS configuration struct */
       primme_svds_initialize(&primme_svds);
    
       /* Set problem matrix */
       primme_svds.matrixMatvec = LauchliMatrixMatvec;
       primme_svds.matrix = &mu;
                               /* Function that implements the matrix-vector products
                                  A*x and A^t*x  */
      
       /* Set problem parameters */
       primme_svds.m = 500;
       primme_svds.n = 100; /* set problem dimension */
       primme_svds.numSvals = 4;   /* Number of wanted singular values */
       primme_svds.eps = 1e-6;     /* ||r|| <= eps * ||matrix|| */
       primme_svds.target = primme_svds_smallest;
                                   /* Seeking for the largest singular values  */
    
       /* Set preconditioner (optional) */
       primme_svds.applyPreconditioner = LauchliApplyPreconditioner;
    
       /* Set method to solve the singular value problem and
          the underneath eigenvalue problem (optional) */
       primme_svds_set_method(primme_svds_hybrid, PRIMME_DYNAMIC, PRIMME_DEFAULT_MIN_TIME, &primme_svds);
       /*  Set hybrid method with PRIMME_DYNAMIC and PRIMME_DEFAULT_MIN_TIME as the underneath eigensolver configuration
           for the first and the second stage, respectively.
           PRIMME_DYNAMIC uses a runtime heuristic to choose the fastest method between
           PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can set another
           method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */
    
       primme_svds.printLevel = 3;
    
       /* Set parameters for the underneath eigensolver if you know what are you doing (optional) */
       primme_svds.primme.locking = 1; 
       primme_svds.primme.restartingParams.maxPrevRetain = 3;
    
        /* Display PRIMME SVDS configuration struct (optional) */
       primme_svds_display_params(primme_svds);
    
       /* Allocate space for converged Ritz values and residual norms */
       svals = new double[primme_svds.numSvals];
       svecs = new std::complex<double>[(primme_svds.n+primme_svds.m)
                              *primme_svds.numSvals];
       rnorms = new double[primme_svds.numSvals];
    
       /* Call primme_svds  */
       ret = zprimme_svds(svals, svecs, rnorms, &primme_svds);