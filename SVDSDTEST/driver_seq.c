 /******************************************************************************
 *  Sequential driver for dprimme. Calling format:
 *
 * 	    seq_dprimme DriverConfigFileName SolverConfigFileName
 *
 *  DriverConfigFileName  includes the path and filename of the matrix
 *  			  as well as preconditioning information (eg., 
 *  			  ILUT parameters).
 *  			  Currently, for reading the input matrix,
 *  			  full coordinate format (.mtx) and upper triangular 
 *  			  coordinate format (.U) are supported.
 *
 *     	   Example file:  DriverConf 
 *
 *  SolverConfigFileName  includes all dprimme required information
 *  			  as stored in primme data structure.
 *
 *     	   Example files: FullConf  Full customization of primme
 *  		          LeanConf  Use a preset method and some customization
 *  		          MinConf   Provide ONLY a preset method and numEvals
 *
 ******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <math.h>
#include "driver_seq.h"

/* primme_svds.h header file is required to run primme_svds */
#include "primme_svds.h"

/* wtime.h header file is included so primme's timimg functions can be used */
#include "wtime.h"

/******************************************************************************/
int main (int argc, char *argv[]) {

   /* Timing vars */
   double ut1,ut2,st1,st2,wt1,wt2;

   /* Matrix */
   int m, n, nnz;
   double fnorm;
   CSRMatrix matrix;

   /* Preconditioner */
   CSRMatrix Factors;

   /* Files */
   char *DriverConfigFileName, *SolverConfigFileName;
   
   /* Driver and solver I/O arrays and parameters */
   double *svals, *svecs, *rnorms;
   driver_params driver;
   primme_svds_params primme_svds;
#ifdef Cplusplus     /* C++ has a stricter type checking */
   void (*precond_function)(void *, void *, int *, primme_svds_params *);
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
      ret = readfullMTXR(driver.matrixFileName, &matrix.AElts, &matrix.JA, 
         &matrix.IA, &m, &n, &nnz);
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   else if (driver.matrixFileName[strlen(driver.matrixFileName)-1] == 'U') {
      // coordinate format storing only upper triangular part
      ret = readUpperMTX(driver.matrixFileName, &matrix.AElts, &matrix.JA,
         &matrix.IA, &n, &nnz);
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }
   else {  
      //Harwell Boeing format NOT IMPLEMENTED
      //ret = readmt()
      ret = -1;
      if (ret < 0) {
         fprintf(stderr, "ERROR: Could not read matrix file\n");
         return(-1);
      }
   }

   fnorm = frobeniusNorm(m, matrix.IA, matrix.AElts);

/* ------------------------------------------------------------------------- */
/*  Set up preconditioner if needed. We provide these sample options         */
/*  in driver.PrecChoice:  (driver.PrecChoice is read in read_driver_params) */
/*     choice = 0  no preconditioner                                         */
/*     choice = 1  K=Diag(A-shift),   shift provided once by user            */
/*     choice = 2  K=Diag(A-shift_i), shifts provided by primme every step   */
/*     choice = 3  K=ILUT(A-shift)  , shift provided once by user            */
/* ------------------------------------------------------------------------- */
   if (m == n){
       ret = create_preconditioner
           (matrix, &Factors, &precond_function, n, nnz, driver);
       if (ret < 0) {
           fprintf(stderr, "ERROR: Could not create requested preconditioner \n");
           return(-1);
       }
   }

   /* --------------------------------------------------------------------- */
   /*    Primme solver setup                                                */
   /*       primme_initialize  (not needed if ALL primme struct members set)*/
   /*       primme_set_method  (bypass it to fully customize your solver)   */
   /* --------------------------------------------------------------------- */

   /* ---------------------------------- */
   /* Initialize defaults in primme_svds */
   /* ---------------------------------- */
   primme_svds_initialize(&primme_svds);

   /* -------------------------------------------- */
   /* Read in the primme_svds configuration file   */
   /* -------------------------------------------- */
   primme_svds.m     = m;
   primme_svds.n     = n;
   primme_svds.aNorm = fnorm; /* ||A||_frobenius. A configFile entry overwrites it */

   if (read_solver_params(SolverConfigFileName, driver.outputFileName, 
			   &primme_svds) < 0) {
      fprintf(stderr, "Reading solver parameters failed\n");
      return(-1);
   }
   
   /* --------------------------------------- */
   /* Set up matrix vector and preconditioner */
   /* --------------------------------------- */

   primme_svds.matrixMatvec           = MatrixMatvec;
   if (m == n){
      primme_svds.applyPreconditioner = precond_function;
   }
   /* --------------------------------------- */
   /* Optional: provide matrix/preconditioner */
   /* --------------------------------------- */
   primme_svds.matrix            = &matrix;
   if (m == n){
      primme_svds.preconditioner = &Factors;
   }

   /* --------------------------------------- */
   /* Display given parameter configuration   */
   /* Place this after the dprimme() to see   */
   /* any changes dprimme() made to primme    */
   /* --------------------------------------- */

   fprintf(primme_svds.outputFile," Matrix: %s\n",driver.matrixFileName);
   primme_svds_display_params(primme_svds);

   /* --------------------------------------------------------------------- */
   /* 	              Run the dprimme_svds solver                           */
   /* --------------------------------------------------------------------- */

   /* Allocate space for converged singular triplets and residual norms */

   svals = (double *)primme_calloc(primme_svds.numSvals, sizeof(double), "svals");
//   svecs = (double *)primme_calloc((primme_svds.m + primme_svds.n)*
//	   (primme_svds.numSvals+primme_svds.maxBlockSize), sizeof(double), "svecs");
   svecs = (double *)primme_calloc((primme_svds.m + primme_svds.n)*
	   primme_svds.numSvals, sizeof(double), "svecs");
   rnorms = (double *)primme_calloc(primme_svds.numSvals, sizeof(double), "rnorms");

   /* ------------------------ */
   /* Initial guess (optional) */
   /* ------------------------ */
   /* Initialize all vectors to avoid uninitialized iniital guesses */
//   for (i=0;i<(primme_svds.n+primme_svds.m)*primme_svds.numSvals;i++) 
//       svecs[i] = 1/sqrt(primme_svds.n);
   if (primme_svds.svdsMethod == primme_svds_augmented){
       for (i=0;i<(primme_svds.n+primme_svds.m);i++) 
           svecs[i] = 1/sqrt(primme_svds.n);
   }
   else if (primme_svds.svdsMethod == primme_svds_normalequations ||
       primme_svds.svdsMethod == primme_svds_hybrid) {
       for (i=0;i<min(primme_svds.n, primme_svds.m);i++) 
           svecs[i] = 1/sqrt(primme_svds.n);
   }

   /* ------------- */
   /*  Call primme  */
   /* ------------- */

   wt1 = primme_get_wtime(); 
   primme_get_time(&ut1,&st1);

   ret = dprimme_svds(svals, svecs, rnorms, &primme_svds);

   wt2 = primme_get_wtime();
   primme_get_time(&ut2,&st2);

   primme_svds_display_params(primme_svds);
   /* --------------------------------------- */
   /* Display given parameter configuration   */
   /* Place this after the dprimme() to see   */
   /* any changes dprimme() made to primme    */
   /* --------------------------------------- */

   fprintf(primme_svds.outputFile," Matrix: %s\n",driver.matrixFileName);

   /* --------------------------------------------------------------------- */
   /* Reporting                                                             */
   /* --------------------------------------------------------------------- */

   primme_PrintStackTrace(primme_svds.primme);

   fprintf(primme_svds.outputFile, "Wallclock Runtime   : %-f\n", wt2-wt1);
   fprintf(primme_svds.outputFile, "User Time           : %f seconds\n", ut2-ut1);
   fprintf(primme_svds.outputFile, "Syst Time           : %f seconds\n", st2-st1);
                 
   if (primme_svds.procID == 0) {
      for (i=0; i < primme_svds.numSvals; i++) {
         fprintf(primme_svds.outputFile, "Sval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
            svals[i], rnorms[i]); 
      }
      fprintf(primme_svds.outputFile, "%d singular triplets converged\n", primme_svds.primme.initSize);

      fprintf(primme_svds.outputFile, "Tolerance : %-22.15E\n", 
		      				      primme_svds.aNorm*primme_svds.eps);
      fprintf(primme_svds.outputFile, "Iterations: %-d\n", 
		      			     primme_svds.stats.numOuterIterations); 
      fprintf(primme_svds.outputFile, "Restarts  : %-d\n", primme_svds.stats.numRestarts);
      fprintf(primme_svds.outputFile, "Matvecs   : %-d\n", primme_svds.stats.numMatvecs);
      fprintf(primme_svds.outputFile, "Preconds  : %-d\n", primme_svds.stats.numPreconds);

      fprintf(primme_svds.outputFile, "\n\n#,%d,%.1f\n\n", primme_svds.stats.numMatvecs,
         wt2-wt1); 
   }

   if (ret != 0) {
      fprintf(primme_svds.outputFile, "Error: dprimme_svds returned with nonzero exit status\n");
      return -1;
   }

   fclose(primme_svds.outputFile);
   free(svals);
   free(svecs);
   free(rnorms);
//   primme_svds_Free(&primme_svds); 
/*lingfei: seq_dprimme_svds(37111,0x7fff751aa310) 
  malloc: *** error for object 0x101805e08: incorrect 
  checksum for freed object - object was probably 
  modified after being freed. *** set a breakpoint 
  in malloc_error_break to debug */

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
void MatrixMatvec(void *x, void *y, int *blockSize, 
    primme_svds_params *primme_svds, const char *transpose) {

   const char notransp[] = "notransp"; /*used for Matvec*/
   const char transp[] = "transp"; /*used for Matvec*/
   int i;
   double *xvec, *yvec;
   CSRMatrix *matrix;
   
   matrix = (CSRMatrix *)primme_svds->matrix;
   xvec = (double *)x;
   yvec = (double *)y;

   if (strcmp(transpose,notransp) == 0){
      for (i=0;i<*blockSize;i++) {
          amux_(&primme_svds->m, &xvec[primme_svds->nLocal*i], &yvec[primme_svds->nLocal*i], 
		      matrix->AElts, matrix->JA, matrix->IA);
      }
   }
   else if (strcmp(transpose,transp) == 0){
      for (i=0;i<*blockSize;i++) {
//          atmux_(&primme_svds->m, &xvec[primme_svds->nLocal*i], 
//            &yvec[primme_svds->nLocal*i], matrix->AElts, matrix->JA, matrix->IA);
          atmuxr_(&primme_svds->n, &primme_svds->m, &xvec[primme_svds->nLocal*i], 
            &yvec[primme_svds->nLocal*i], matrix->AElts, matrix->JA, matrix->IA);
      }
   }
}


/******************************************************************************
 * Creates one of 3 preconditioners depending on driver.PrecChoice:
 *
 * Our sample choices are
 *
 *   choice = 0  no preconditioner                                       
 *   choice = 1  K=Diag(A-shift),   shift provided once by user         
 *   choice = 2  K=Diag(A-shift_i), shifts provided by primme_svds every step
 *   choice = 3  K=ILUT(A-shift)  , shift provided once by user         
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
 *   choice = 3   ilut				   Apply_ILUT_Prec
 *
 * The user is cautioned that ILUT is a nonsymmetric preconditioner, 
 * which could potentially prevent the QMR inner solver from converging.
 * We have not noticed this in our experiments.
 *
******************************************************************************/
int create_preconditioner(CSRMatrix matrix, CSRMatrix *Factors, 
#ifdef Cplusplus     /* C++ has a stricter type checking */
   void (**precond_function)(void *, void *, int *, primme_svds_params *),
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
	 // Diagonal: K = diag(A-shift I), shift can be primme_svds provided
         Factors->AElts = (double *)primme_calloc(n, sizeof(double), 
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
      case 3: { // ILUT(A-shift I) 
	 int ierr;
         int precondlfil = driver.level;
         double precondTol = driver.threshold;
         double *W1, *W2;
         int *iW1, *iW2, *iW3;

	 if (driver.shift != 0.0L) {
	    shiftCSRMatrix(-driver.shift, n, matrix.IA,matrix.JA,matrix.AElts);
	 }

         // Work arrays
         W1 = (double *)primme_calloc( n+1,  sizeof(double), "W1");
         W2 = (double *)primme_calloc( n,  sizeof(double), "W2");
         iW1 = (int *)primme_calloc( n,  sizeof(int), "iW1");
         iW2 = (int *)primme_calloc( n,  sizeof(int), "iW2");
         iW3 = (int *)primme_calloc( n,  sizeof(int), "iW2");
         //---------------------------------------------------
         // Max size of factorization                       //
         lenFactors = 9*nnz;				    //
         //---------------------------------------------------
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
	 // free workspace
         free(W1); free(W2); free(iW1); free(iW2); free(iW3);
	 }
         *precond_function = Apply_ILUT_Prec;
         break;
      case 4:  // Parasails(A - shift I)
	 // not implemented yet
         break;
   }
   return 0;
}


/******************************************************************************
 * Applies a Davidson type preconditioner
 *
 *    x(i) = (Diag(A) - primme_svds.Shifts(i) I)^(-1) * y(i),   i=1:blockSize
 * 	
 * NOTE that each block vector may have its own shift provided by dprimme_svds
 * in the array primme_svds->ShiftsForPreconditioner
 *
 * To avoid division with too small numbers we limit how small relatively 
 * to ||A|| the denominators can be. In the absense of ||A|| we use 1e-14.
 *
******************************************************************************/

void Apply_Diagonal_Shifted_Prec(void *x, void *y, int *blockSize, 
		            primme_svds_params *primme_svds) {

   int i, j, index;
   double shift, denominator, minDenominator, NormEstimate;
   double *xvec, *yvec;
   CSRMatrix *diagPrecond;

   diagPrecond = (CSRMatrix *)primme_svds->preconditioner;
   xvec = (double *)x;
   yvec = (double *)y;

   if (primme_svds->aNorm >= 0.0L)
      NormEstimate = primme_svds->aNorm;
   else 
      NormEstimate = 1;

   for (j=0;j<*blockSize;j++) {

      index = primme_svds->nLocal*j;
      shift = primme_svds->primme.ShiftsForPreconditioner[j];

      // Apply shifted diagonal preconditioner
      for (i=0;i<primme_svds->nLocal;i++) {
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
 * 	y(i) = P*x(i), i=1:blockSize, 
 * 	with P = (Diag(A)-shift)^(-1)
 *
******************************************************************************/

void Apply_Inv_Diagonal_Prec(void *x, void *y, int *blockSize, 
            primme_svds_params *primme_svds, const char *transpose) {
   const char notransp[] = "notransp"; /*used for Matvec*/
   const char transp[] = "transp"; /*used for Matvec*/
   int i;
   double *xvec, *yvec;
   CSRMatrix *prec;
   
   prec = (CSRMatrix *)primme_svds->preconditioner;
   xvec = (double *)x;
   yvec = (double *)y;
   
   if (strcmp(transpose,notransp) == 0){
      for (i=0;i<*blockSize;i++) {
          amux_(&primme_svds->m, &xvec[primme_svds->nLocal*i], 
            &yvec[primme_svds->nLocal*i], prec->AElts, prec->JA, prec->IA);
      }
   }
   else if (strcmp(transpose,transp) == 0){
      for (i=0;i<*blockSize;i++) {
          atmuxr_(&primme_svds->m, &primme_svds->n, &xvec[primme_svds->nLocal*i], 
            &yvec[primme_svds->nLocal*i], prec->AElts, prec->JA, prec->IA);
      }
   }

//   for (i=0;i<*blockSize;i++) {
//      amux_(&primme_svds->n, &xvec[primme_svds->n*i], &yvec[primme_svds->n*i],
//                      prec->AElts, prec->JA, prec->IA);
//   }
//
}

/******************************************************************************
 * Applies the ILUT preconditioner 
 *
 * 	y(i) = U^(-1)*( L^(-1)*x(i)), i=1:blockSize, 
 * 	with L,U = ilut(A-shift) 
 * 
 * It calls the SPARSKIT lusol0 function for each block vector.
 *
******************************************************************************/

void Apply_ILUT_Prec(void *x, void *y, int *blockSize, primme_svds_params *primme_svds) {

   int i;
   double *xvec, *yvec;
   CSRMatrix *prec;
   
   prec = (CSRMatrix *)primme_svds->preconditioner;
   xvec = (double *)x;
   yvec = (double *)y;

   for (i=0;i<*blockSize;i++) {
      lusol0_(&primme_svds->n,&xvec[primme_svds->n*i],&yvec[primme_svds->n*i],
		      prec->AElts, prec->JA,prec->IA);
   }

}

/******************************************************************************
 * Computed the Frobenius norm of a CSR matrix 
 *
 * 	||A||_frob = sqrt( \sum_{i,j} A_ij^2 )
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
 * 	A = A + shift I
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
 * 	P = (Diag(A)-shift)^(-1)
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
 * 	P = Diag(A)
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
