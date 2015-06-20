/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2005  James R. McCombs,  Andreas Stathopoulos
 *
 *   This file is part of PRIMME.
 *
 *   PRIMME is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   PRIMME is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * --------------------------------------------------------------------------
 *
 *  Sequential driver for dprimme. Calling format:
 *
 *             seq_dprimme DriverConfigFileName SolverConfigFileName
 *
 *  DriverConfigFileName  includes the path and filename of the matrix
 *                            as well as preconditioning information (eg., 
 *                            ILUT parameters).
 *                            Currently, for reading the input matrix,
 *                            full coordinate format (.mtx) and upper triangular 
 *                            coordinate format (.U) are supported.
 *
 *                Example file:  DriverConf 
 *
 *  SolverConfigFileName  includes all dprimme required information
 *                            as stored in primme data structure.
 *
 *                Example files: FullConf  Full customization of primme
 *                            LeanConf  Use a preset method and some customization
 *                            MinConf   Provide ONLY a preset method and numEvals
 *
 ******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include "../PRIMMESRC/DSRC/numerical_d.h"
#include "driver_seq.h"
#include <assert.h>

/* primme.h header file is required to run primme */
#include "primme.h"

/* wtime.h header file is included so primme's timimg functions can be used */
#include "wtime.h"

/******************************************************************************/
int main (int argc, char *argv[]) {

   /* Timing vars */
   double ut1,ut2,st1,st2,wt1,wt2;

   /* Matrix */
   int n, nnz;
   double fnorm;
   CSRMatrix matrix;

   /* Preconditioner */
   CSRMatrix Factors;

   /* Files */
   char *DriverConfigFileName=NULL, *SolverConfigFileName=NULL;
   
   /* Driver and solver I/O arrays and parameters */
   double *evals, *evecs, *rnorms;
   driver_params driver;
   primme_params primme;
   primme_preset_method method;
   void (*precond_function)(void *, void *, int *, primme_params *);

   /* Other miscellaneous items */
   int i;
   int ret, retX=0;

   /* --------------------------------------------------------------------- */
   /*   Read matrix and driver setup                                        */
   /* --------------------------------------------------------------------- */

   /* ------------------------------------------------------- */
   /* Get from command line the names for the 2 config files  */
   /* ------------------------------------------------------- */

   if (argc == 2) {
      DriverConfigFileName = argv[1];
      SolverConfigFileName = argv[1];
   } else if (argc == 3) {
      DriverConfigFileName = argv[1];
      SolverConfigFileName = argv[2];
   } else {
      fprintf(stderr, "Invalid number of arguments.\n");
      return(-1);
   }

   /* ----------------------------- */
   /* Read in the driver parameters */
   /* ----------------------------- */
   if (read_driver_params(DriverConfigFileName, &driver) < 0) {
      fprintf(stderr, "Reading driver parameters failed\n");

      return(-1);
   }

   /* ------------------------------------------ */
   /* Read the matrix and store it in CSR format */
   /* ------------------------------------------ */
   fprintf(stderr," Matrix: %s\n",driver.matrixFileName);


   if (!strcmp("mtx", &driver.matrixFileName[strlen(driver.matrixFileName)-3]))
   {  /* coordinate format storing both lower and upper triangular parts */
      ret = readfullMTX(driver.matrixFileName, &matrix.AElts, &matrix.JA, 
         &matrix.IA, &n, &nnz);
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   else if (driver.matrixFileName[strlen(driver.matrixFileName)-1] == 'U') {
      /* coordinate format storing only upper triangular part */
      ret = readUpperMTX(driver.matrixFileName, &matrix.AElts, &matrix.JA,
         &matrix.IA, &n, &nnz);
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   else {  
      /* Harwell Boeing format NOT IMPLEMENTED */
      /* ret = readmt() */
      ret = -1;
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }

   fnorm = frobeniusNorm(n, matrix.IA, matrix.AElts);

/* ------------------------------------------------------------------------- */
/*  Set up preconditioner if needed. We provide these sample options         */
/*  in driver.PrecChoice:  (driver.PrecChoice is read in read_driver_params) */
/*     choice = 0  no preconditioner                                         */
/*     choice = 1  K=Diag(A-shift),   shift provided once by user            */
/*     choice = 2  K=Diag(A-shift_i), shifts provided by primme every step   */
/*     choice = 3  K=ILUT(A-shift)  , shift provided once by user            */
/* ------------------------------------------------------------------------- */

   ret = create_preconditioner
         (matrix, &Factors, &precond_function, n, nnz, driver);
   if (ret < 0) {
      fprintf(stderr, "ERROR: Could not create requested preconditioner \n");
      return(-1);
   }

   /* --------------------------------------------------------------------- */
   /*    Primme solver setup                                                */
   /*       primme_initialize  (not needed if ALL primme struct members set)*/
   /*       primme_set_method  (bypass it to fully customize your solver)   */
   /* --------------------------------------------------------------------- */

   /* ----------------------------- */
   /* Initialize defaults in primme */
   /* ----------------------------- */
   primme_initialize(&primme);

   /* --------------------------------------- */
   /* Read in the primme configuration file   */
   /* --------------------------------------- */
   primme.n     = n;
   primme.aNorm = fnorm; /* ||A||_frobenius. A configFile entry overwrites it */

   if (read_solver_params(SolverConfigFileName, driver.outputFileName, 
                           &primme, &method) < 0) {
      fprintf(stderr, "Reading solver parameters failed\n");
      return(-1);
   }

   /* --------------------------------------- */
   /* Pick one of the default methods(if set) */
   /* --------------------------------------- */

   if (primme_set_method(method, &primme) < 0 ) {
      fprintf(primme.outputFile, "No preset method. Using custom settings\n");
   }

   /* --------------------------------------- */
   /* Optional: report memory requirements    */
   /* --------------------------------------- */

   ret = dprimme(NULL,NULL,NULL,&primme);
   fprintf(primme.outputFile,"PRIMME will allocate the following memory:\n");
   fprintf(primme.outputFile,"real workspace, %ld bytes\n",primme.realWorkSize);
   fprintf(primme.outputFile,"int  workspace, %d bytes\n",primme.intWorkSize);
   
   /* --------------------------------------- */
   /* Set up matrix vector and preconditioner */
   /* --------------------------------------- */

   primme.matrixMatvec        = MatrixMatvec;
   primme.applyPreconditioner = precond_function;

   /* --------------------------------------- */
   /* Optional: provide matrix/preconditioner */
   /* --------------------------------------- */
   primme.matrix         = &matrix;
   primme.preconditioner = &Factors;

   /* --------------------------------------- */
   /* Display given parameter configuration   */
   /* Place this after the dprimme() to see   */
   /* any changes dprimme() made to primme    */
   /* --------------------------------------- */

   driver_display_params(driver, primme.outputFile); 
   primme_display_params(primme);
   driver_display_method(method, primme.outputFile);

   /* --------------------------------------------------------------------- */
   /*                      Run the dprimme solver                           */
   /* --------------------------------------------------------------------- */

   /* Allocate space for converged Ritz values and residual norms */

   evals = (double *)primme_calloc(primme.numEvals, sizeof(double), "evals");
   evecs = (double *)primme_calloc(
                primme.n*primme.numEvals,sizeof(double), "evecs");
   rnorms = (double *)primme_calloc(primme.numEvals, sizeof(double), "rnorms");

   /* ------------------------ */
   /* Initial guess (optional) */
   /* ------------------------ */

   /* Read initial guess from a file */
   if (driver.initialGuessesFileName[0] && primme.initSize != 0) {
      FILE *f = fopen(driver.initialGuessesFileName, "rb");
      double d;
      int cols;
      assert(f);
      ret = fread(&d, sizeof(d), 1, f); assert(ret); assert(((int)d) == sizeof(d));
      fread(&d, sizeof(d), 1, f); assert(ret); assert(((int)d) == primme.n);
      fread(&d, sizeof(d), 1, f); assert(ret); cols = d;
      for (i=0; i<min(cols, primme.initSize); i++) {
         ret = fread(&evecs[primme.nLocal*i], sizeof(d), primme.nLocal, f); assert(ret == primme.nLocal);
      }

      /* Perturb the initial guesses by a vector with some norm  */
      if (driver.initialGuessesPert > 0) {
         double *r = (double *)primme_calloc(primme.nLocal,sizeof(double), "random");
         double norm;
         int j;
         for (i=0; i<min(cols, primme.initSize); i++) {
            Num_larnv_dprimme(2, primme.iseed, primme.nLocal, r);
            norm = sqrt(Num_dot_dprimme(primme.nLocal, r, 1, r, 1));
            for (j=0; j<primme.nLocal; j++) evecs[primme.nLocal*i+j] += r[j]/norm*driver.initialGuessesPert;
         }
         free(r);
      }
      Num_larnv_dprimme(2, primme.iseed, (primme.initSize-i)*primme.nLocal, &evecs[primme.nLocal*i]);
      fclose(f);
   } else if (primme.initSize > 0) {
      Num_larnv_dprimme(2, primme.iseed, primme.initSize*primme.nLocal, evecs);
   } else {
      Num_larnv_dprimme(2, primme.iseed, primme.nLocal, evecs);
   }


   /* ------------- */
   /*  Call primme  */
   /* ------------- */

   wt1 = primme_get_wtime(); 
   primme_get_time(&ut1,&st1);

   ret = dprimme(evals, evecs, rnorms, &primme);

   wt2 = primme_get_wtime();
   primme_get_time(&ut2,&st2);

   if (driver.checkXFileName[0]) {
      retX = check_solution(driver.checkXFileName, &primme, evals, evecs, rnorms,
                            &matrix);
   }

   /* --------------------------------------------------------------------- */
   /* Save evecs (optional)                                                 */
   /* --------------------------------------------------------------------- */
   if (driver.saveXFileName[0]) {
      FILE *f = fopen(driver.saveXFileName, "wb");
      double d;
      int ret;
      assert(f);
      d = sizeof(d); ret = fwrite(&d, sizeof(double), 1, f); assert(ret);
      d = primme.n; ret = fwrite(&d, sizeof(double), 1, f); assert(ret);
      d = primme.initSize; ret = fwrite(&d, sizeof(double), 1, f); assert(ret);
      for (i=0; i<primme.initSize; i++) {
         ret = fwrite(&evecs[primme.nLocal*i], sizeof(d), primme.nLocal, f); assert(ret == primme.nLocal);
      }
      ret = fwrite(&primme, sizeof(primme), 1, f); assert(ret);
      fclose(f);
   }

   /* --------------------------------------------------------------------- */
   /* Reporting                                                             */
   /* --------------------------------------------------------------------- */

   /* --------------------------------------------------------------------- */
   /* Check how PRIMME may have changed some input parameters               */
   /* primme_display_params(primme); */
   /* --------------------------------------------------------------------- */

   primme_PrintStackTrace(primme);

   fprintf(primme.outputFile, "Wallclock Runtime   : %-f\n", wt2-wt1);
   fprintf(primme.outputFile, "User Time           : %f seconds\n", ut2-ut1);
   fprintf(primme.outputFile, "Syst Time           : %f seconds\n", st2-st1);


   if (primme.procID == 0) {
      for (i=0; i < primme.numEvals; i++) {
         fprintf(primme.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
            evals[i], rnorms[i]); 
      }
      fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);

      fprintf(primme.outputFile, "Tolerance : %-22.15E\n", 
                                                            primme.aNorm*primme.eps);
      fprintf(primme.outputFile, "Iterations: %-d\n", 
                                                    primme.stats.numOuterIterations); 
      fprintf(primme.outputFile, "Restarts  : %-d\n", primme.stats.numRestarts);
      fprintf(primme.outputFile, "Matvecs   : %-d\n", primme.stats.numMatvecs);
      fprintf(primme.outputFile, "Preconds  : %-d\n", primme.stats.numPreconds);
      if (primme.locking && primme.intWork && primme.intWork[0] == 1) {
         fprintf(primme.outputFile, "\nA locking problem has occurred.\n");
         fprintf(primme.outputFile, 
            "Some eigenpairs do not have a residual norm less than the tolerance.\n");
         fprintf(primme.outputFile, 
            "However, the subspace of evecs is accurate to the required tolerance.\n");
      }


      fprintf(primme.outputFile, "\n\n#,%d,%.1f\n\n", primme.stats.numMatvecs,
         wt2-wt1); 

      switch (primme.dynamicMethodSwitch) {
         case -1: fprintf(primme.outputFile, 
               "Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
         case -2: fprintf(primme.outputFile, 
               "Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
         case -3: fprintf(primme.outputFile, 
               "Recommended method for next run: DYNAMIC (close call)\n"); break;
      }

   }

   if (ret != 0) {
      fprintf(primme.outputFile, 
         "Error: dprimme returned with nonzero exit status\n");
      return -1;
   }

   if (retX != 0) {
      fprintf(primme.outputFile, 
         "Error: found some issues in the solution return by dprimme\n");
      return -1;
   }


   fclose(primme.outputFile);
   primme_Free(&primme);

   return(0);

}
/******************************************************************************/
/* END OF MAIN DRIVER FUNCTION                                                */
/******************************************************************************/

/******************************************************************************/
/* Matvec preconditioner and other utilities                                  */

/******************************************************************************
 * Applies the matrix vector multiplication on a block of vectors.
 * Because a block function is not available, we call blockSize times
 * the SPARSKIT function amux(). Note the (void *) parameters x, y that must 
 * be cast as doubles for use in amux()
 *
******************************************************************************/
void MatrixMatvec(void *x, void *y, int *blockSize, primme_params *primme) {
   
   int i;
   double *xvec, *yvec;
   CSRMatrix *matrix;
   
   matrix = (CSRMatrix *)primme->matrix;
   xvec = (double *)x;
   yvec = (double *)y;

   for (i=0;i<*blockSize;i++) {
      amux_(&primme->n, &xvec[primme->nLocal*i], &yvec[primme->nLocal*i], 
                      matrix->AElts, matrix->JA, matrix->IA);
   }

}


/******************************************************************************
 * Creates one of 3 preconditioners depending on driver.PrecChoice:
 *
 * Our sample choices are
 *
 *   choice = 0  no preconditioner                                       
 *   choice = 1  K=Diag(A-shift),   shift provided once by user         
 *   choice = 2  K=Diag(A-shift_i), shifts provided by primme every step
 *   choice = 3  K=ILUT(A-shift)  , shift provided once by user         
 *
 * on return:
 *    Factors:                pointer to the preconditioner CSR structure
 *    precond_function: pointer to the function that applies the preconditioner
 *
 * For each choice we have:
 *
 *                   generated by function:             Applied by function
 *                 -------------------------           -------------------
 *   choice = 0             -                                   NULL
 *   choice = 1   generate_Inv_Diagonal_Prec           Apply_Inv_Diagonal_Prec
 *   choice = 2   generate_Diagonal_Prec           Apply_Diagonal_Shifted_Prec
 *   choice = 3   ilut                                   Apply_ILUT_Prec
 *
 * The user is cautioned that ILUT is a nonsymmetric preconditioner, 
 * which could potentially prevent the QMR inner solver from converging.
 * We have not noticed this in our experiments.
 *
******************************************************************************/
int create_preconditioner(CSRMatrix matrix, CSRMatrix *Factors, 
   void (**precond_function)(void *, void *, int *, primme_params *),
   int n, int nnz, driver_params driver) {

   int lenFactors;

   *precond_function = NULL;

   switch (driver.PrecChoice) {
      case 0:  /* No preconditioner created */
         break;
      case 1: case 2:
         /* Diagonal: K = diag(A-shift I), shift can be primme provided */
         Factors->AElts = (double *)primme_calloc(n, sizeof(double), 
                                                             "Factors.AElts");
         Factors->IA = (int *)primme_calloc(n+1, sizeof(int), "Factors.IA");
         Factors->JA = (int *)primme_calloc(n, sizeof(int),"Factors.JA");
         if (driver.PrecChoice == 1) {
            generate_Inv_Diagonal_Prec(n, driver.shift, matrix.IA, matrix.JA, 
               matrix.AElts, Factors->IA, Factors->JA, Factors->AElts);
            *precond_function = Apply_Inv_Diagonal_Prec;
            break;
         } 
         else {
            generate_Diagonal_Prec(n, matrix.IA, matrix.JA, matrix.AElts, 
               Factors->IA, Factors->JA, Factors->AElts);
            *precond_function = Apply_Diagonal_Shifted_Prec;
            break;
         }
      case 3: { /* ILUT(A-shift I) */
         int ierr;
         int precondlfil = driver.level;
         double precondTol = driver.threshold;
         double *W1, *W2;
         int *iW1, *iW2, *iW3;

         if (driver.shift != 0.0L) {
            shiftCSRMatrix(-driver.shift, n, matrix.IA,matrix.JA,matrix.AElts);
         }

         /* Work arrays */
         W1 = (double *)primme_calloc( n+1,  sizeof(double), "W1");
         W2 = (double *)primme_calloc( n,  sizeof(double), "W2");
         iW1 = (int *)primme_calloc( n,  sizeof(int), "iW1");
         iW2 = (int *)primme_calloc( n,  sizeof(int), "iW2");
         iW3 = (int *)primme_calloc( n,  sizeof(int), "iW2");
         /* --------------------------------------------------- */
         /* Max size of factorization                           */
         lenFactors = 9*nnz;
         /* --------------------------------------------------- */
         Factors->AElts = (double *)primme_calloc(lenFactors,
                                                          sizeof(double), "iluElts");
         Factors->JA = (int *)primme_calloc(lenFactors, sizeof(int), "Jilu");
         Factors->IA = (int *)primme_calloc(n+1, sizeof(int), "Iilu");
      
         ilut_(&n,matrix.AElts,matrix.JA,matrix.IA,&precondlfil,&precondTol,
                     Factors->AElts, Factors->JA, Factors->IA, &lenFactors, 
                     W1,W2,iW1,iW2,iW3,&ierr);
      
         if (ierr != 0)  {
            fprintf(stderr, "ILUT factorization could not be completed\n");
            return(-1);
         }

          if (driver.shift != 0.0L) {
            shiftCSRMatrix(driver.shift, n, matrix.IA,matrix.JA,matrix.AElts);
         }
         /* free workspace */
         free(W1); free(W2); free(iW1); free(iW2); free(iW3);
         }
         *precond_function = Apply_ILUT_Prec;
         break;
      case 4:  /* Parasails(A - shift I) */
         /* not implemented yet */
         break;
   }
   return 0;
}


/******************************************************************************
 * Applies a Davidson type preconditioner
 *
 *    x(i) = (Diag(A) - primme.Shifts(i) I)^(-1) * y(i),   i=1:blockSize
 *         
 * NOTE that each block vector may have its own shift provided by dprimme
 * in the array primme->ShiftsForPreconditioner
 *
 * To avoid division with too small numbers we limit how small relatively 
 * to ||A|| the denominators can be. In the absense of ||A|| we use 1e-14.
 *
******************************************************************************/

void Apply_Diagonal_Shifted_Prec(void *x, void *y, int *blockSize, 
                            primme_params *primme) {

   int i, j, index;
   double shift, denominator, minDenominator, NormEstimate;
   double *xvec, *yvec;
   CSRMatrix *diagPrecond;

   diagPrecond = (CSRMatrix *)primme->preconditioner;
   xvec = (double *)x;
   yvec = (double *)y;

   if (primme->aNorm >= 0.0L)
      NormEstimate = primme->aNorm;
   else 
      NormEstimate = 1;

   for (j=0;j<*blockSize;j++) {

      index = primme->nLocal*j;
      shift = primme->ShiftsForPreconditioner[j];

      /* Apply shifted diagonal preconditioner */
      for (i=0;i<primme->nLocal;i++) {
           denominator = diagPrecond->AElts[i] - shift;

           minDenominator = 1e-14*NormEstimate;
           if (fabs(denominator) < minDenominator) {
              if (denominator < 0) denominator = -minDenominator;
              else denominator = minDenominator;
           }
           yvec[index+i] = xvec[index+i]/denominator;
      }
   }
}


/******************************************************************************
 * Applies the (already inverted) diagonal preconditioner
 *
 *         y(i) = P*x(i), i=1:blockSize, 
 *         with P = (Diag(A)-shift)^(-1)
 *
******************************************************************************/

void Apply_Inv_Diagonal_Prec(void *x, void *y, int *blockSize, 
                                                  primme_params *primme) {
   int i;
   double *xvec, *yvec;
   CSRMatrix *prec;
   
   prec = (CSRMatrix *)primme->preconditioner;
   xvec = (double *)x;
   yvec = (double *)y;

   for (i=0;i<*blockSize;i++) {
      amux_(&primme->n, &xvec[primme->n*i], &yvec[primme->n*i],
                      prec->AElts, prec->JA, prec->IA);
   }

}

/******************************************************************************
 * Applies the ILUT preconditioner 
 *
 *         y(i) = U^(-1)*( L^(-1)*x(i)), i=1:blockSize, 
 *         with L,U = ilut(A-shift) 
 * 
 * It calls the SPARSKIT lusol0 function for each block vector.
 *
******************************************************************************/

void Apply_ILUT_Prec(void *x, void *y, int *blockSize, primme_params *primme) {

   int i;
   double *xvec, *yvec;
   CSRMatrix *prec;
   
   prec = (CSRMatrix *)primme->preconditioner;
   xvec = (double *)x;
   yvec = (double *)y;

   for (i=0;i<*blockSize;i++) {
      lusol0_(&primme->n,&xvec[primme->n*i],&yvec[primme->n*i],
                      prec->AElts, prec->JA,prec->IA);
   }

}

/******************************************************************************
 * Computed the Frobenius norm of a CSR matrix 
 *
 *         ||A||_frob = sqrt( \sum_{i,j} A_ij^2 )
 *
******************************************************************************/
double frobeniusNorm(int n, int *IA, double *AElts) {

   int i, j;
   double fnorm;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   fnorm = 0.0L;

   for (i=0; i < n; i++) {
      for (j=IA[i]; j <= IA[i+1]-1; j++) {
         fnorm = fnorm + AElts[j-1]*AElts[j-1];
      }
   }

   return (sqrt(fnorm)); 
}  
   

/******************************************************************************
 * Shifts a CSR matrix by a shift
 *
 *         A = A + shift I
 *
******************************************************************************/
void shiftCSRMatrix(double shift, int n, int *IA, int *JA, double *AElts) {

   int i, j;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   for (i=0; i < n; i++) {
      for (j=IA[i]; j <= IA[i+1]-1; j++) {

         if (JA[j-1]-1 == i) {
            AElts[j-1] = AElts[j-1] + shift;
         }

      }
   }

}

/******************************************************************************
 * Generates a shifted and inverted diagonal preconditioner
 *
 *         P = (Diag(A)-shift)^(-1)
 *
 * If A_ii - shift is too close to zero in a relative to Frob norm sense,
 * a small value is replaced. 
 * The preconditioner P is then simply multiplied to a vector.
******************************************************************************/
void generate_Inv_Diagonal_Prec(int n, double shift, 
   int *IA, int *JA, double *AElts, int *PIA, int *PJA, double *PElts) {

   int i, j;
   double temp, atemp, frobNorm;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   /* compute the frobenius norm to estimate the smallest element to consider */
   frobNorm = frobeniusNorm(n, IA, AElts);

   for (i=0; i < n; i++) {
      for (j=IA[i]; j <= IA[i+1]-1; j++) {

         if (JA[j-1]-1 == i) {
            temp = AElts[j-1]-shift;
            atemp = fabs(temp);
            atemp = max(1e-15*frobNorm, atemp);
            if (temp < 0 ) atemp = -atemp;
            PElts[i] = 1.0L/atemp;
            PIA[i] = i+1;
            PJA[i] = i+1;
         }
      }
   }

   PIA[n] = n+1;
}

/******************************************************************************
 * Generates the diagonal of A.
 *
 *         P = Diag(A)
 *
 * This will be used with solver provided shifts as (P-shift_i)^(-1) 
******************************************************************************/
void generate_Diagonal_Prec(int n, int *IA, int *JA, double *AElts, 
   int *PIA, int *PJA, double *PElts) {

   int i, j;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   for (i=0; i < n; i++) {
      for (j=IA[i]; j <= IA[i+1]-1; j++) {

         if (JA[j-1]-1 == i) {
            PElts[i] = AElts[j-1];
            PIA[i] = i+1;
            PJA[i] = i+1;
         }
      }
   }

   PIA[n] = n+1;
}

int check_solution(const char *checkXFileName, primme_params *primme, double *evals,
                   double *evecs, double *rnorms, CSRMatrix *matrix) {

   double *Ax, eval0, rnorm0, prod, *r, d, *X, *h;
   int i, j, cols, ret, retX;
   FILE *f;
   primme_params primme0;

   /* Read stored eigenvectors */
   f = fopen(checkXFileName, "rb");
   assert(f);
   ret = fread(&d, sizeof(d), 1, f); assert(ret); assert(((int)d) == sizeof(d));
   fread(&d, sizeof(d), 1, f); assert(ret); assert(((int)d) == primme->n);
   fread(&d, sizeof(d), 1, f); assert(ret); cols = d;
   X = (double *)primme_calloc(primme->n*cols, sizeof(double), "X");
   for (i=0; i<cols; i++) {
      ret = fread(&X[primme->nLocal*i], sizeof(d), primme->nLocal, f); assert(ret == primme->nLocal);
   }

   /* Read stored primme_params */
   ret = fread(&primme0, sizeof(primme0), 1, f); assert(ret);

   /* Check primme_params */
   assert(primme0.n == primme->n);
   if (primme0.numEvals == primme->numEvals &&
       primme0.target == primme->target &&
       primme0.numTargetShifts == primme->numTargetShifts &&
       primme0.dynamicMethodSwitch == primme->dynamicMethodSwitch &&
       primme0.locking == primme->locking &&
       primme0.numOrthoConst == primme->numOrthoConst &&
       primme0.maxBasisSize == primme->maxBasisSize &&
       primme0.minRestartSize == primme->minRestartSize &&
       primme0.aNorm == primme->aNorm &&
       primme0.eps == primme->eps &&
       primme0.restartingParams.scheme == primme->restartingParams.scheme &&
       primme0.restartingParams.maxPrevRetain == primme->restartingParams.maxPrevRetain &&
       primme0.correctionParams.precondition == primme->correctionParams.precondition &&
       primme0.correctionParams.robustShifts == primme->correctionParams.robustShifts &&
       primme0.correctionParams.maxInnerIterations == primme->correctionParams.maxInnerIterations &&
       primme0.correctionParams.projectors.LeftQ  == primme->correctionParams.projectors.LeftQ  &&
       primme0.correctionParams.projectors.LeftX  == primme->correctionParams.projectors.LeftX  &&
       primme0.correctionParams.projectors.RightQ == primme->correctionParams.projectors.RightQ &&
       primme0.correctionParams.projectors.RightX == primme->correctionParams.projectors.RightX &&
       primme0.correctionParams.projectors.SkewQ  == primme->correctionParams.projectors.SkewQ  &&
       primme0.correctionParams.projectors.SkewX  == primme->correctionParams.projectors.SkewX  &&
       primme0.correctionParams.convTest == primme->correctionParams.convTest &&
       primme0.correctionParams.relTolBase == primme->correctionParams.relTolBase) {
      if (abs(primme0.stats.numOuterIterations - primme->stats.numOuterIterations) > primme->stats.numOuterIterations*3/100+1) {
         fprintf(stderr, "Warning: discrepancy in numOuterIterations, %d should be close to %d\n", primme->stats.numOuterIterations, primme0.stats.numOuterIterations);
         retX = 1;
      }
      if (primme0.initSize != primme->initSize) {
         fprintf(stderr, "Warning: discrepancy in the number of converged pairs, %d should be close to %d\n", primme->initSize, primme0.initSize);
         retX = 1;
      }
   }
   else {
      fprintf(stderr, "Warning: discrepancy in some member of primme\n");
      retX = 1;
   }

   h = (double *)primme_calloc(cols, sizeof(double), "h");
   Ax = (double *)primme_calloc(primme->nLocal, sizeof(double), "Ax");
   r = (double *)primme_calloc(primme->n*primme->initSize, sizeof(double), "rwork");
   
   for (i=0; i < primme->initSize; i++) {
      /* Check |V(:,i)'A*V(:,i) - evals[i]| < |r|*|A| */
      amux_(&primme->n, &evecs[primme->nLocal*i], Ax, matrix->AElts, matrix->JA, matrix->IA);
      eval0 = Num_dot_dprimme(primme->nLocal, &evecs[primme->nLocal*i], 1, Ax, 1);
      if (fabs(evals[i] - eval0) > rnorms[i]*primme->aNorm) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E should be close to %-22.1E\n", i, evals[i], eval0);
         retX = 1;
      }
      /* Check |A*V(:,i) - (V(:,i)'A*V(:,i))*V(:,i)| < |r| */
      for (j=0; j<primme->nLocal; j++) r[j] = Ax[j] - evals[i]*evecs[primme->nLocal*i+j];
      rnorm0 = sqrt(Num_dot_dprimme(primme->nLocal, r, 1, r, 1));
      if (rnorms[i] > primme->eps*primme->aNorm*sqrt((double)(i+1))) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E, residual %5E is larger than expected %5E\n", i, evals[i], rnorms[i], primme->eps*primme->aNorm*sqrt((double)(i+1)));
         retX = 1;
      }
      if (rnorm0 > primme->eps*primme->aNorm*sqrt((double)(i+1))) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E, RR residual %5E is larger than tolerance %5E\n", i, evals[i], rnorm0, primme->eps*primme->aNorm*sqrt((double)(i+1)));
         retX = 1;
      }
      /* Check X'V(:,i) >= sqrt(1-2|r|), assuming residual of X is less than the tolerance */
      Num_gemv_dprimme("C", primme->n, cols, 1.0, X, primme->n, &evecs[primme->nLocal*i], 1, 0., h, 1);
      prod = Num_dot_dprimme(primme->nLocal, h, 1, h, 1);
      if (prod < sqrt(1.-2.*rnorms[i])) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E not found on X\n", i, evals[i]);
         retX = 1;
      }
   }
   free(h);
   free(X);
   free(r);
   free(Ax);

   return retX; 
}
