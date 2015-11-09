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
 *  Sequential driver for zprimme. Calling format:
 *
 * 	    seq_zprimme DriverConfigFileName SolverConfigFileName
 *
 *  DriverConfigFileName  includes the path and filename of the matrix
 *  			  as well as preconditioning information.
 *  			  Currently, for reading the input matrix,
 *  			  full coordinate format (.mtx) is supported only.
 *
 *     	   Example file:  DriverConf 
 *
 *  SolverConfigFileName  includes all zprimme required information
 *  			  as stored in primme data structure.
 *
 *     	   Example files: FullConf  Full customization of primme
 *  		          LeanConf  Use a preset method and some customization
 *  		          MinConf   Provide ONLY a preset method and numEvals
 *
 ******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <strings.h> //BSD
#include <string.h>  //Linux
#include <math.h>
#include "driver_seq.h"

/* primme.h header file is required to run primme */
#include "primme.h"

/* wtime.h header file is included so primme's timimg functions can be used */
#include "wtime.h"

/******************************************************************************/
int main (int argc, char *argv[]) {

   /* Timing vars */
   double wt1,wt2,ut1,ut2,st1,st2;

   /* Matrix */
   int n, nnz;
   double fnorm;
   CSRMatrix matrix;

   /* Preconditioner */
   CSRMatrix Factors;

   /* Files */
   char *DriverConfigFileName, *SolverConfigFileName;
   
   /* Driver and solver I/O arrays and parameters */
   double *evals, *rnorms;
   Complex_Z *evecs;
   driver_params driver;
   primme_params primme;
   primme_preset_method method;
#ifdef Cplusplus     /* C++ has a stricter type checking */
   void (*precond_function)(void *, void *, int *, primme_params *);
#else
   void *precond_function;
#endif

   /* Other miscellaneous items */
   int i;
   int ret;

   /* --------------------------------------------------------------------- */
   /*   Read matrix and driver setup                                        */
   /* --------------------------------------------------------------------- */

   /* ------------------------------------------------------- */
   /* Get from command line the names for the 2 config files  */
   /* ------------------------------------------------------- */

   if (argc == 3) {
      DriverConfigFileName = argv[1];
      SolverConfigFileName = argv[2];
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
   {  // coordinate format storing both lower and upper triangular parts
      ret = readfullMTX(driver.matrixFileName, &matrix.AElts, &matrix.JA, 
         &matrix.IA, &n, &nnz);
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   else {  
      //No other formats implemented for Complex
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

   ret = zprimme(NULL,NULL,NULL,&primme);
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
   /* Place this after the zprimme() to see   */
   /* any changes zprimme() made to primme    */
   /* --------------------------------------- */

   fprintf(primme.outputFile," Matrix: %s\n",driver.matrixFileName);
   primme_display_params(primme);

   /* --------------------------------------------------------------------- */
   /* 	                   Run the zprimme solver                           */
   /* --------------------------------------------------------------------- */

   /* Allocate space for converged Ritz values and residual norms */

   evals = (double *)primme_calloc(primme.numEvals, sizeof(double), "evals");
   evecs = (Complex_Z *)primme_valloc(
		primme.n*primme.numEvals*sizeof(Complex_Z), "evecs");
   rnorms = (double *)primme_calloc(primme.numEvals, sizeof(double), "rnorms");

   /* ------------------------ */
   /* Initial guess (optional) */
   /* ------------------------ */
       for (i=0;i<primme.n;i++) {
	  evecs[i].r=1/sqrt(primme.n);
	  evecs[i].i=0.0L;
       }

   /* ------------- */
   /*  Call primme  */
   /* ------------- */

   wt1 = primme_get_wtime(); 
   primme_get_time(&ut1,&st1);

   ret = zprimme(evals, evecs, rnorms, &primme);

   wt2 = primme_get_wtime();
   primme_get_time(&ut2,&st2);

   /* --------------------------------------------------------------------- */
   /* Reporting                                                             */
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
      if (primme.locking && primme.intWork[0] == 1) {
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
         "Error: zprimme returned with nonzero exit status\n");
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
 * the modified SPARSKIT function zamux(). Note the (void *) parameters x, y 
 * that must be cast as Complex_Z for use in zamux()
 *
******************************************************************************/
void MatrixMatvec(void *x, void *y, int *blockSize, primme_params *primme) {
   
   int i;
   CSRMatrix *matrix;
   Complex_Z *xvec, *yvec;
   
   matrix = (CSRMatrix *)primme->matrix;
   xvec = (Complex_Z *)x;
   yvec = (Complex_Z *)y;
   
   for (i=0;i<*blockSize;i++) {
      zamux_(&primme->n, &xvec[primme->nLocal*i], &yvec[primme->nLocal*i], 
		      matrix->AElts, matrix->JA, matrix->IA);
   }

}


/******************************************************************************
 * Creates one of 2 preconditioners depending on driver.PrecChoice:
 *
 * Our sample choices are
 *
 *   choice = 0  no preconditioner                                       
 *   choice = 1  K=Diag(A-shift),   shift provided once by user         
 *   choice = 2  K=Diag(A-shift_i), shifts provided by primme every step
 *   choice = 3,4 not implemented for Complex
 *
 * on return:
 *    Factors:		pointer to the preconditioner CSR structure
 *    precond_function: pointer to the function that applies the preconditioner
 *
 * For each choice we have:
 *
 *		   generated by function:             Applied by function
 *		 -------------------------           -------------------
 *   choice = 0  	   -				   NULL
 *   choice = 1   generate_Inv_Diagonal_Prec	   Apply_Inv_Diagonal_Prec
 *   choice = 2   generate_Diagonal_Prec           Apply_Diagonal_Shifted_Prec
 *
******************************************************************************/
int create_preconditioner(CSRMatrix matrix, CSRMatrix *Factors, 
#ifdef Cplusplus     /* C++ has a stricter type checking */
   void (**precond_function)(void *, void *, int *, primme_params *),
#else
   void **precond_function,
#endif
   int n, int nnz, driver_params driver) {

   int lenFactors;

   *precond_function = NULL;

   switch (driver.PrecChoice) {
      case 0:  // No preconditioner created
         break;
      case 1: case 2:
	 // Diagonal: K = diag(A-shift I), shift can be primme provided
         Factors->AElts = (Complex_Z *)primme_calloc(n, sizeof(Complex_Z), 
                                                             "Factors.AElts");
         Factors->IA = (int *)primme_calloc(n+1, sizeof(int), "Factors.IA");
         Factors->JA = (int *)primme_calloc(n, sizeof(int),"Factors.JA");
         printf("Generating diagonal preconditioner");
	 if (driver.PrecChoice == 1) {
	    printf(" with the user provided shift %e\n",driver.shift);
            generate_Inv_Diagonal_Prec(n, driver.shift, matrix.IA, matrix.JA, 
	       matrix.AElts, Factors->IA, Factors->JA, Factors->AElts);
            *precond_function = Apply_Inv_Diagonal_Prec;
	    break;
	 } 
	 else {
	    printf(" that will use solver provided shifts\n");
	    generate_Diagonal_Prec(n, matrix.IA, matrix.JA, matrix.AElts, 
	       Factors->IA, Factors->JA, Factors->AElts);
            *precond_function = Apply_Diagonal_Shifted_Prec;
            break;
	 }
      case 3:  // ilut
	 // not implemented yet in complex
         break;
      case 4:  // Parasails(A - shift I)
	 // not implemented yet in complex
         break;
   }
   return 0;
}


/******************************************************************************
 * Applies a Davidson type preconditioner
 *
 *    x(i) = (Diag(A) - primme.Shifts(i) I)^(-1) * y(i),   i=1:blockSize
 * 	
 * NOTE that each block vector may have its own shift provided by zprimme
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
   CSRMatrix *diagPrecond;
   Complex_Z *xvec, *yvec;
   
   diagPrecond = (CSRMatrix *)primme->preconditioner;
   xvec = (Complex_Z *)x;
   yvec = (Complex_Z *)y;

   if (primme->aNorm >= 0.0L)
      NormEstimate = primme->aNorm;
   else 
      NormEstimate = 1;

   for (j=0;j<*blockSize;j++) {

      index = primme->nLocal*j;
      shift = primme->ShiftsForPreconditioner[j];

      // Apply shifted diagonal preconditioner
      for (i=0;i<primme->nLocal;i++) {
	   denominator = diagPrecond->AElts[i].r - shift;

	   minDenominator = 1e-14*NormEstimate;
	   if (fabs(denominator) < minDenominator) {
	      if (denominator < 0) denominator = -minDenominator;
	      else denominator = minDenominator;
	   }
	   yvec[index+i].r = xvec[index+i].r/denominator;
	   yvec[index+i].i = xvec[index+i].i/denominator;
      }
   }
}


/******************************************************************************
 * Applies the (already inverted) diagonal preconditioner
 *
 * 	y(i) = P*x(i), i=1:blockSize, 
 * 	with P = (Diag(A)-shift)^(-1)
 *
******************************************************************************/

void Apply_Inv_Diagonal_Prec(void *x, void *y, int *blockSize, 
		                               primme_params *primme) {
   int i;
   CSRMatrix *prec;
   Complex_Z *xvec, *yvec;
   
   prec = (CSRMatrix *)primme->preconditioner;
   xvec = (Complex_Z *)x;
   yvec = (Complex_Z *)y;

   for (i=0;i<*blockSize;i++) {
      zamux_(&primme->n, &xvec[primme->n*i], &yvec[primme->n*i],
                      prec->AElts, prec->JA, prec->IA);
   }

}

/******************************************************************************
 * Computed the Frobenius norm of a CSR matrix 
 *
 * 	||A||_frob = sqrt( \sum_{i,j} A_ij^2 )
 *
******************************************************************************/
double frobeniusNorm(int n, int *IA, Complex_Z *AElts) {

   int i, j;
   double fnorm;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   fnorm = 0.0L;

   for (i=0; i < n; i++) {
      for (j=IA[i]; j <= IA[i+1]-1; j++) {
         fnorm = fnorm + z_abs_primme(AElts[j-1])*z_abs_primme(AElts[j-1]);
      }
   }

   return (sqrt(fnorm)); 
}  
   

/******************************************************************************
 * Shifts a CSR matrix by a shift
 *
 * 	A = A + shift I
 *
******************************************************************************/
void shiftCSRMatrix(double shift, int n, int *IA, int *JA, Complex_Z *AElts) {

   int i, j;

   /* IA and JA are indexed using C indexing, but their contents */
   /* assume Fortran indexing.  Thus, the contents of IA and JA  */
   /* must be decremented before being used in C.                */

   for (i=0; i < n; i++) {
      for (j=IA[i]; j <= IA[i+1]-1; j++) {

         if (JA[j-1]-1 == i) {
            AElts[j-1].r = AElts[j-1].r + shift;
         }

      }
   }

}

/******************************************************************************
 * Generates a shifted and inverted diagonal preconditioner
 *
 * 	P = (Diag(A)-shift)^(-1)
 *
 * If A_ii - shift is too close to zero in a relative to Frob norm sense,
 * a small value is replaced. 
 * The preconditioner P is then simply multiplied to a vector.
******************************************************************************/
void generate_Inv_Diagonal_Prec(int n, double shift, 
   int *IA, int *JA, Complex_Z *AElts, int *PIA, int *PJA, Complex_Z *PElts) {

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
	    temp = AElts[j-1].r-shift;
	    atemp = fabs(temp);
	    atemp = max(1e-15*frobNorm, atemp);
	    if (temp < 0 ) atemp = -atemp;
            PElts[i].r = 1.0L/atemp;
            PElts[i].i = 0.0L;
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
 * 	P = Diag(A)
 *
 * This will be used with solver provided shifts as (P-shift_i)^(-1) 
******************************************************************************/
void generate_Diagonal_Prec(int n, int *IA, int *JA, Complex_Z *AElts, 
   int *PIA, int *PJA, Complex_Z *PElts) {

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

