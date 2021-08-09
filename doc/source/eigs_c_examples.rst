.. role:: ccode(code) 
   :language: c

.. highlight:: c

.. _simple:

Simple C Eigs Example
---------------------

.. code:: c

   #include <stdlib.h>
   #include <stdio.h>
   #include <math.h>
   
   #include "primme.h"   /* header file is required to run primme */ 
   
   void LaplacianMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
   void LaplacianApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
   void convtest_sm(double *eval, void *evec, double *resNorm, int *isconv, primme_params *primme, int *ierr);
   
   
   int main (int argc, char *argv[]) {
   
      /* Solver arrays and parameters */
      double *evals;    /* Array with the computed eigenvalues */
      double *rnorms;   /* Array with the computed eigenpairs residual norms */
      double *evecs;    /* Array with the computed eigenvectors;
                           first vector starts in evecs[0],
                           second vector starts in evecs[primme.n],
                           third vector starts in evecs[primme.n*2]...  */
      primme_params primme;
                        /* PRIMME configuration struct */
   
      /* Other miscellaneous items */
      int ret;
      int i;
   
      /* Set default values in PRIMME configuration struct */
      primme_initialize(&primme);
   
      /* Set problem matrix */
      primme.matrixMatvec = LaplacianMatrixMatvec;
                              /* Function that implements the matrix-vector product
                                 A*x for solving the problem A*x = l*x */
     
      /* Set problem parameters */
      primme.n = 100; /* set problem dimension */
      primme.numEvals = 10;   /* Number of wanted eigenpairs */
      primme.eps = 1e-9;      /* ||r|| <= eps * ||matrix|| */
      primme.convTestFun = convtest_sm;
      primme.target = primme_smallest;
                              /* Wanted the smallest eigenvalues */
   
      /* Set preconditioner (optional) */
      primme.applyPreconditioner = LaplacianApplyPreconditioner;
      primme.correctionParams.precondition = 1;
   
      /* Set advanced parameters if you know what are you doing (optional) */
      /*
      primme.maxBasisSize = 14;
      primme.minRestartSize = 4;
      primme.maxBlockSize = 1;
      primme.maxMatvecs = 1000;
      */
   
      /* Set method to solve the problem */
      primme_set_method(PRIMME_DYNAMIC, &primme);
      /* DYNAMIC uses a runtime heuristic to choose the fastest method between
          PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
          set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */
   
      /* Display PRIMME configuration struct (optional) */
      primme_display_params(primme);
   
      /* Allocate space for converged Ritz values and residual norms */
      evals = (double*)malloc(primme.numEvals*sizeof(double));
      evecs = (double*)malloc(primme.n*primme.numEvals*sizeof(double));
      rnorms = (double*)malloc(primme.numEvals*sizeof(double));
   
      /* Call primme  */
      ret = dprimme(evals, evecs, rnorms, &primme);
   
   }
   
   void LaplacianMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
      
      int i;            /* vector index, from 0 to *blockSize-1*/
      int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
      double *xvec;     /* pointer to i-th input vector x */
      double *yvec;     /* pointer to i-th output vector y */
      
      for (i=0; i<*blockSize; i++) {
         xvec = (double *)x + *ldx*i;
         yvec = (double *)y + *ldy*i;
         for (row=0; row<primme->n; row++) {
            yvec[row] = 0.0;
            if (row-1 >= 0) yvec[row] += -1.0*xvec[row-1];
            yvec[row] += 2.0*xvec[row];
            if (row+1 < primme->n) yvec[row] += -1.0*xvec[row+1];
         }      
      }
      *err = 0;
   }
   
   /* This performs Y = M^{-1} * X, where
   
      - X, input dense matrix of size primme.n x blockSize;
      - Y, output dense matrix of size primme.n x blockSize;
      - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
   */
   
   void convtest_sm(double *eval, void *evec, double *resNorm, int *isconv, primme_params *primme, int *ierr){
      *isconv = abs(*eval) > 0.1 * (*resNorm);
      *ierr = 0;
   }
   
   void LaplacianApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {
      
      int i;            /* vector index, from 0 to *blockSize-1*/
      int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
      double *xvec;     /* pointer to i-th input vector x */
      double *yvec;     /* pointer to i-th output vector y */
       
      for (i=0; i<*blockSize; i++) {
         xvec = (double *)x + *ldx*i;
         yvec = (double *)y + *ldy*i;
         for (row=0; row<primme->n; row++) {
            yvec[row] = xvec[row]/2.;
         }      
      }
      *ierr = 0;
   }
   

.. _parallel:

Parallel C Eigs Example
-----------------------

.. code:: c

   #include <stdlib.h>
   #include <stdio.h>
   #include <math.h>
   #include <mpi.h>
   #include <assert.h>

   #include "primme.h"   /* header file is required to run primme */ 

   void DiagonalMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
   static void par_GlobalSum(void *sendBuf, void *recvBuf, int *count,
                           primme_params *primme, int *ierr);

   #ifndef min
   #  define min(a, b) ((a) < (b) ? (a) : (b))
   #endif

   int main (int argc, char *argv[]) {

      /* Solver arrays and parameters */
      float *evals;    /* Array with the computed eigenvalues */
      float *rnorms;   /* Array with the computed eigenpairs residual norms */
      float *evecs;    /* Array with the computed eigenvectors;
                           first vector starts in evecs[0],
                           second vector starts in evecs[primme.n],
                           third vector starts in evecs[primme.n*2]...  */
      primme_params primme;
                        /* PRIMME configuration struct */

      /* Other miscellaneous items */
      int ret;
      int i;

      /* Initialize the infrastructure necessary for communication */
      MPI_Init(&argc, &argv);

      /* Set default values in PRIMME configuration struct */
      primme_initialize(&primme);

      /* Set problem matrix */
      primme.matrixMatvec = DiagonalMatrixMatvec;
                              /* Function that implements the matrix-vector product
                                 A*x for solving the problem A*x = l*x */
   
      /* Set problem parameters */
      primme.n = 1000; /* set problem dimension */
      primme.numEvals = 1000;   /* Number of wanted eigenpairs */
      primme.eps = .1;      /* ||r|| <= eps * ||matrix|| */
      primme.target = primme_largest;
                              /* Wanted the smallest eigenvalues */

      /* Set advanced parameters if you know what are you doing (optional) */
      /*
      primme.maxBasisSize = 14;
      primme.minRestartSize = 4;
      primme.maxBlockSize = 1;
      primme.maxMatvecs = 1000;
      */

      /* Set method to solve the problem */
      primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
      /* DYNAMIC uses a runtime heuristic to choose the fastest method between
         PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
         set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */

      /* Set parallel parameters */
      MPI_Comm comm = MPI_COMM_WORLD;
      MPI_Comm_size(comm, &primme.numProcs);
      MPI_Comm_rank(comm, &primme.procID);
      primme.commInfo = &comm; /* User-defined member to pass the communicator to
                                 globalSumReal and broadcastReal */
      /* In this example, the matrix is distributed by rows, and the first
      * processes may have an extra row in order to distribute the remaining rows
      * n % numProcs */
      PRIMME_INT nLocal = primme.n / primme.numProcs +
                        (primme.n % primme.numProcs > primme.procID ? 1 : 0);
      primme.nLocal = nLocal; /* Number of local rows */
      primme.globalSumReal = par_GlobalSum;

      /* Display PRIMME configuration struct (optional) */
      if (primme.procID == 0) primme_display_params(primme);

      /* Allocate space for converged Ritz values and residual norms */
      evals = (float*)malloc(primme.numEvals*sizeof(float));
      evecs = (float*)malloc(primme.n*primme.numEvals*sizeof(float));
      rnorms = (float*)malloc(primme.numEvals*sizeof(float));

      /* Call primme  */
      ret = sprimme(evals, evecs, rnorms, &primme);
   }
   void DiagonalMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   
      int i;            /* vector index, from 0 to *blockSize-1*/
      int row;          /* local matrix row index, from 0 to nLocal */
      /* In this example, row0 is the global index of the first local row */
      int row0 = primme->n / primme->numProcs * primme->procID +
               min(primme->n % primme->numProcs, primme->procID);
      float *xvec;     /* pointer to i-th input vector x */
      float *yvec;     /* pointer to i-th output vector y */
      
      for (i=0; i<*blockSize; i++) {
         xvec = (float *)x + *ldx*i;
         yvec = (float *)y + *ldy*i;
         for (row = 0; row < primme->nLocal; row++) {
            /* The diagonal matrix has the spectrum of a Laplacial */
            float v = sin(M_PI * (row + row0 + 1) / 2.0 / (primme->n + 1));
            yvec[row] = 4. * v * v * xvec[row];
         }
      }
      *err = 0;
   }

   static void par_GlobalSum(void *sendBuf, void *recvBuf, int *count, 
                           primme_params *primme, int *ierr) {
      MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;

      if (sendBuf == recvBuf) {
      *ierr = MPI_Allreduce(MPI_IN_PLACE, recvBuf, *count, MPI_FLOAT, MPI_SUM, communicator) != MPI_SUCCESS;
      } else {
      *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPI_FLOAT, MPI_SUM, communicator) != MPI_SUCCESS;
      }
   }



.. _dmagmaEigs:

dmagma Eigs Example
--------------------

.. code:: c

   #include <stdlib.h>
   #include <stdio.h>
   #include <string.h>
   #include <math.h>

   #include "magma_v2.h"
   #include "magmasparse.h"

   #include "primme.h"   /* header file is required to run primme */ 

   #include <time.h>

   void magmaSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
   void magmaDummy(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);


   int main (int argc, char *argv[]) {

      /* Solver arrays and parameters */
      double *evals;    /* Array with the computed eigenvalues */
      double *rnorms;   /* Array with the computed eigenpairs residual norms */
      double *evecs;    /* Array with the computed eigenvectors;
                           first vector starts in evecs[0],
                           second vector starts in evecs[primme.n],
                           third vector starts in evecs[primme.n*2]...  */
      primme_params primme;
                        /* PRIMME configuration struct */

      /* Other miscellaneous items */
      int n=1000; /* problem size */
      int ret;
      int i,j;

      int *col, *row;
      double *val;

      row = (int*) calloc(n+1, sizeof(int));
      col = (int*) calloc(n+(n>0?n-1:0)*2, sizeof(int));
      val = (double*) calloc(n+(n>0?n-1:0)*2, sizeof(double));

      for (i = j = 0; i < n; i++) {
         row[i] = j;
         if (i > 0)   {col[j] = i-1; val[j] = -1.0; j++;}
                     col[j] = i  ; val[j] =  2.0; j++;
         if (i < n-1) {col[j] = i+1; val[j] = -1.0; j++;}
      }
      row[n] = j;

      /* Initialize MAGMA and create some LA structures */
      magma_init();
      magma_queue_t queue;
      magma_queue_create(0, &queue);

      magma_d_matrix A={Magma_CSR}, dA={Magma_CSR};

      /* Pass the matrix to MAGMA and copy it to the GPU */
      magma_dcsrset(n, n, row, col, val, &A, queue);
      magma_dmtransfer(A, &dA, Magma_CPU, Magma_DEV, queue);

      /* Set default values in PRIMME configuration struct */
      primme_initialize(&primme);
   
      /* Set problem parameters */
      primme.n = n; /* set problem dimension */
      primme.numEvals = 6;   /* Number of wanted eigenpairs */
      primme.eps = 1e-12;      /* ||r|| <= eps * ||matrix|| */
      primme.target = primme_smallest;
                              /* Wanted the smallest eigenvalues */

      /* Set problem matrix */
      primme.matrixMatvec = magmaSparseMatrixMatvec;
      primme.matrix = &dA;
                              /* Function that implements the matrix-vector product
                                 A*x for solving the problem A*x = l*x */
   
      /* Set preconditioner (optional) */
      primme.applyPreconditioner = magmaDummy;
      primme.correctionParams.precondition = 1;

      /* Set advanced parameters if you know what are you doing (optional) */
      /*
      primme.maxBasisSize = 14;
      primme.minRestartSize = 4;
      primme.maxBlockSize = 1;
      primme.maxMatvecs = 1000;
      */

      /* Set method to solve the problem */
      primme_set_method(PRIMME_DYNAMIC, &primme);
   //   primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
   //   primme_set_method(PRIMME_DEFAULT_MIN_TIME, &primme);
      /* DYNAMIC uses a runtime heuristic to choose the fastest method between
         PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
         set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */

      /* Display PRIMME configuration struct (optional) */
      primme_display_params(primme);

      /* Allocate space for converged Ritz values and residual norms */
      evals = (double*)malloc(primme.numEvals*sizeof(double));
      magma_dmalloc(&evecs, primme.n*primme.numEvals);
      rnorms = (double*)malloc(primme.numEvals*sizeof(double));

      primme.queue = &queue;

      /*
         clock_t start,end;

         start = clock();
         primme.funcTime = 0;
      */
      time_t rawtime,rawtime2;
      struct tm * timeinfo,* timeinfo2;

      time ( &rawtime );
      timeinfo = localtime ( &rawtime );
      printf ( "Start ->Current local time and date: %s", asctime (timeinfo) );


      /* Call primme  */
      ret = magma_dprimme(evals, evecs, rnorms, &primme);
   }

   void magmaSparseMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   
      int i;            /* vector index, from 0 to *blockSize-1*/
      double *xvec;     
      double *yvec;     
      magma_d_matrix *A = primme->matrix;
   
      for (i=0; i<*blockSize; i++) {
         magma_d_matrix vx = {Magma_CSR};  /* i-th input vector x */
         magma_d_matrix vy = {Magma_CSR};  /* i-th output vector y */

         magma_dvset_dev(primme->n, 1, (double *)x + *ldx*i, &vx, *(magma_queue_t*)primme->queue);
         magma_dvset_dev(primme->n, 1, (double *)y + *ldy*i, &vy, *(magma_queue_t*)primme->queue);

         magma_d_spmv(1.0, *A, vx, 0.0, vy, *(magma_queue_t*)primme->queue);
      }
      *err = 0;
   }

   void magmaDummy(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
      magma_dcopymatrix(primme->n, *blockSize, (double*)x, *ldx, (double*)y, *ldy, *(magma_queue_t*)primme->queue);
      *err = 0;
   }

.. _normalEigs:

Normal Eigs Example
-------------------

.. code:: c

   #include <stdlib.h>
   #include <stdio.h>
   #include <math.h>
   #include <complex.h>
   #include "primme.h"   /* header file is required to run primme */ 
   
   void LaplacianLikeMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
   void LaplacianLikeApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
   
   int main (int argc, char *argv[]) {
   
      /* Solver arrays and parameters */
      complex double *evals;    /* Array with the computed eigenvalues */
      double *rnorms;   /* Array with the computed eigenpairs residual norms */
      complex double *evecs;    /* Array with the computed eigenvectors;
                           first vector starts in evecs[0],
                           second vector starts in evecs[primme.n],
                           third vector starts in evecs[primme.n*2]...  */
      primme_params primme;
                        /* PRIMME configuration struct */
      double targetShifts[1];
   
      /* Other miscellaneous items */
      int ret;
      int i;
   
      /* Set default values in PRIMME configuration struct */
      primme_initialize(&primme);
   
      /* Set problem matrix */
      primme.matrixMatvec = LaplacianLikeMatrixMatvec;
                              /* Function that implements the matrix-vector product
                                 A*x for solving the problem A*x = l*x */
     
      /* Set problem parameters */
      primme.n = 100; /* set problem dimension */
      primme.numEvals = 10;   /* Number of wanted eigenpairs */
      primme.eps = 1e-9;      /* ||r|| <= eps * ||matrix|| */
      targetShifts[0] = .5;
      primme.targetShifts = targetShifts;
      primme.numTargetShifts = 1;
      primme.target = primme_closest_abs;
                              /* Wanted the smallest eigenvalues */
   
      /* Set preconditioner (optional) */
      primme.applyPreconditioner = LaplacianLikeApplyPreconditioner;
      primme.correctionParams.precondition = 1;
   
      /* Set advanced parameters if you know what are you doing (optional) */
      /*
      primme.maxBasisSize = 14;
      primme.minRestartSize = 4;
      primme.maxBlockSize = 1;
      primme.maxMatvecs = 1000;
      */
   
      /* Set method to solve the problem */
      /* NOTE: PRIMME_DEFAULT_MIN_TIME is not supported normal operators */
      primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
      /* You can set other methods, such as PRIMME_LOBPCG_OrthoBasis_Window */
   
      /* NOTE: cheap Olsen approximation is not supported for normal operators */
      primme.correctionParams.projectors.RightX = 0;
   
      primme.printLevel = 3;
   
      /* Display PRIMME configuration struct (optional) */
      primme_display_params(primme);
   
      /* Allocate space for converged Ritz values and residual norms */
      evals = (complex double*)malloc(primme.numEvals*sizeof(complex double));
      evecs = (complex double*)malloc(primme.n*primme.numEvals*sizeof(complex double));
      rnorms = (double*)malloc(primme.numEvals*sizeof(double));
   
      /* Call primme  */
      ret = zprimme_normal(evals, evecs, rnorms, &primme);
   }

   void LaplacianLikeMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
      
      int i;            /* vector index, from 0 to *blockSize-1*/
      int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
      complex double *xvec;     /* pointer to i-th input vector x */
      complex double *yvec;     /* pointer to i-th output vector y */
      
      for (i=0; i<*blockSize; i++) {
         xvec = (complex double *)x + *ldx*i;
         yvec = (complex double *)y + *ldy*i;
         for (row=0; row<primme->n; row++) {
            yvec[row] = 0.0;
            if (row-1 >= 0) yvec[row] += (-1.0+I)*xvec[row-1];
            yvec[row] += 2.0*xvec[row];
            if (row+1 < primme->n) yvec[row] += (-1.0+I)*xvec[row+1];
         }      
      }
      *err = 0;
   }

   /* This performs Y = M^{-1} * X, where

      - X, input dense matrix of size primme.n x blockSize;
      - Y, output dense matrix of size primme.n x blockSize;
      - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
   */

   void LaplacianLikeApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {
      
      int i;            /* vector index, from 0 to *blockSize-1*/
      int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
      complex double *xvec;     /* pointer to i-th input vector x */
      complex double *yvec;     /* pointer to i-th output vector y */
      
      for (i=0; i<*blockSize; i++) {
         xvec = (complex double *)x + *ldx*i;
         yvec = (complex double *)y + *ldy*i;
         for (row=0; row<primme->n; row++) {
            yvec[row] = xvec[row]/2.;
         }      
      }
      *ierr = 0;
   }

.. _cppEigs:

C++ Eigs Example
----------------

.. code:: cpp

   #include <stdlib.h>
   #include <stdio.h>
   #include <math.h>
   #include <complex>
   #include "primme.h"   /* header file is required to run primme */ 

   void LaplacianMatrixMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   
      int i;            /* vector index, from 0 to *blockSize-1*/
      int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
      std::complex<double> *xvec;     /* pointer to i-th input vector x */
      std::complex<double> *yvec;     /* pointer to i-th output vector y */
      
      for (i=0; i<*blockSize; i++) {
         xvec = (std::complex<double> *)x + *ldx*i;
         yvec = (std::complex<double> *)y + *ldy*i;
         for (row=0; row<primme->n; row++) {
            yvec[row] = 0.0;
            if (row-1 >= 0) yvec[row] += -1.0*xvec[row-1];
            yvec[row] += 2.0*xvec[row];
            if (row+1 < primme->n) yvec[row] += -1.0*xvec[row+1];
         }      
      }
      *err = 0;
   }

   /* This performs Y = M^{-1} * X, where

      - X, input dense matrix of size primme.n x blockSize;
      - Y, output dense matrix of size primme.n x blockSize;
      - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
   */

   void LaplacianApplyPreconditioner(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {
      
      int i;            /* vector index, from 0 to *blockSize-1*/
      int row;          /* Laplacian matrix row index, from 0 to matrix dimension */
      std::complex<double> *xvec;     /* pointer to i-th input vector x */
      std::complex<double> *yvec;     /* pointer to i-th output vector y */
      
      for (i=0; i<*blockSize; i++) {
         xvec = (std::complex<double> *)x + *ldx*i;
         yvec = (std::complex<double> *)y + *ldy*i;
         for (row=0; row<primme->n; row++) {
            yvec[row] = xvec[row]/2.;
         }      
      }
      *ierr = 0;
   }
   
   int main (int argc, char *argv[]) {

      /* Solver arrays and parameters */
      double *evals;    /* Array with the computed eigenvalues */
      double *rnorms;   /* Array with the computed eigenpairs residual norms */
      std::complex<double> *evecs;    /* Array with the computed eigenvectors;
                           first vector starts in evecs[0],
                           second vector starts in evecs[primme.n],
                           third vector starts in evecs[primme.n*2]...  */
      primme_params primme;
                        /* PRIMME configuration struct */
      double targetShifts[1];

      /* Other miscellaneous items */
      int ret;
      int i;

      /* Set default values in PRIMME configuration struct */
      primme_initialize(&primme);

      /* Set problem matrix */
      primme.matrixMatvec = LaplacianMatrixMatvec;
                              /* Function that implements the matrix-vector product
                                 A*x for solving the problem A*x = l*x */
   
      /* Set problem parameters */
      primme.n = 100; /* set problem dimension */
      primme.numEvals = 10;   /* Number of wanted eigenpairs */
      primme.eps = 1e-9;      /* ||r|| <= eps * ||matrix|| */
      primme.target = primme_smallest;
                              /* Wanted the smallest eigenvalues */

      /* Set preconditioner (optional) */
      primme.applyPreconditioner = LaplacianApplyPreconditioner;
      primme.correctionParams.precondition = 1;

      /* Set advanced parameters if you know what are you doing (optional) */
      /*
      primme.maxBasisSize = 14;
      primme.minRestartSize = 4;
      primme.maxBlockSize = 1;
      primme.maxMatvecs = 1000;
      */

      /* Set method to solve the problem */
      primme_set_method(PRIMME_DYNAMIC, &primme);
      /* DYNAMIC uses a runtime heuristic to choose the fastest method between
         PRIMME_DEFAULT_MIN_TIME and PRIMME_DEFAULT_MIN_MATVECS. But you can
         set another method, such as PRIMME_LOBPCG_OrthoBasis_Window, directly */

      /* Display PRIMME configuration struct (optional) */
      primme_display_params(primme);

      /* Allocate space for converged Ritz values and residual norms */
      evals = new double[primme.numEvals];
      evecs = new std::complex<double>[primme.n*primme.numEvals];
      rnorms = new double[primme.numEvals];

      /* Call primme  */
      ret = zprimme(evals, evecs, rnorms, &primme);
   }