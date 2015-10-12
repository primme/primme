/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2015 College of William & Mary,
 *   James R. McCombs, Eloy Romero Alcalde, Andreas Stathopoulos, Lingfei Wu
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
 *******************************************************************************
 * File: driver.c
 *
 * Purpose - driver that can read a matrix from a file and compute some
 *           eigenvalues using PRIMME.
 *
 *  Parallel driver for PRIMME. Calling format:
 *
 *             primme DriverConfigFileName SolverConfigFileName
 *
 *  DriverConfigFileName  includes the path and filename of the matrix
 *                            as well as preconditioning information (eg., 
 *                            ParaSails parameters).
 *                            Currently, for reading the input matrix,
 *                            full coordinate format (.mtx) and upper triangular 
 *                            coordinate format (.U) are supported.
 *
 *         Example file:  DriverConf
 *
 *  SolverConfigFileName  includes all d/zprimme required information
 *                            as stored in primme_params data structure.
 *
 *             Example files: FullConf  Full customization of PRIMME
 *                            LeanConf  Use a preset method and some customization
 *                            MinConf   Provide ONLY a preset method and numEvals.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>

#ifdef USE_MPI
#  include <mpi.h>
#endif
#ifdef USE_NATIVE
#  include "native.h"
#endif
#ifdef USE_PARASAILS
#  include "parasailsw.h"
#endif
#ifdef USE_PETSC
# include "petscw.h"
#endif

/* primme.h header file is required to run primme */
#include "primme.h"
#include "shared_utils.h"
/* wtime.h header file is included so primme's timimg functions can be used */
#include "wtime.h"

#define ASSERT_MSG(COND, RETURN, ...) { if (!(COND)) {fprintf(stderr, "Error in " __FUNCT__ ": " __VA_ARGS__); return (RETURN);} }

static int real_main (int argc, char *argv[]);
static int setMatrixAndPrecond(driver_params *driver, primme_params *primme, int **permutation);
#ifdef USE_MPI
static void broadCast(primme_params *primme, primme_preset_method *method, 
   driver_params *driver, int master, MPI_Comm comm);
static void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme);
#endif
static int check_solution(const char *checkXFileName, primme_params *primme, double *evals,
                          PRIMME_NUM *evecs, double *rnorms, int *perm);
static int destroyMatrixAndPrecond(driver_params *driver, primme_params *primme, int *permutation);
static int writeBinaryEvecsAndPrimmeParams(const char *fileName, PRIMME_NUM *X, int *perm,
                                           primme_params *primme);
static int readBinaryEvecsAndPrimmeParams(const char *fileName, PRIMME_NUM *X, PRIMME_NUM **Xout,
                                          int n, int Xcols, int *Xcolsout, int nLocal,
                                          int *perm, primme_params *primme);



int main (int argc, char *argv[]) {
   int ret;
#if defined(USE_PETSC)
   PetscInt ierr;
#endif

#if defined(USE_MPI) && !defined(USE_PETSC)
   MPI_Init(&argc, &argv);
#elif defined(USE_PETSC)
   PetscInitialize(&argc, &argv, NULL, NULL);
#endif

   ret = real_main(argc, argv);

#if defined(USE_MPI) && !defined(USE_PETSC)
   if (ret >= 0) {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
   }
   else {
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
   }
#elif defined(USE_PETSC)
   ierr = PetscFinalize(); CHKERRQ(ierr);
#endif

   return ret;
}

/******************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "real_main"
static int real_main (int argc, char *argv[]) {

   /* Timing vars */
   double wt1,wt2;
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
   double ut1,ut2,st1,st2;
#endif

   /* Files */
   char *DriverConfigFileName=NULL, *SolverConfigFileName=NULL;
   
   /* Driver and solver I/O arrays and parameters */
   double *evals, *rnorms;
   PRIMME_NUM *evecs;
   driver_params driver;
   primme_params primme;
   primme_preset_method method;
   int *permutation = NULL;

   /* Other miscellaneous items */
   int ret, retX=0;
   int i;
   int master = 1;
   int procID = 0;

#ifdef USE_MPI
   MPI_Comm comm;
   int numProcs;
   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &procID);

   comm = MPI_COMM_WORLD;
   master = (procID == 0);
#endif

   primme_initialize(&primme);

   if (master) {
      /* ------------------------------------------------------------ */
      /* Get from command line the names for the 1 or 2 config files  */
      /* NOTE: PETSc arguments starts with '-' and they shouldn't be  */
      /*       considered as configuration files.                     */
      /* ------------------------------------------------------------ */
   
      if (argc == 2 || (argc > 2 && argv[2][0] == '-')) {
         DriverConfigFileName = argv[1];
         SolverConfigFileName = argv[1];
      } else if (argc >= 3) {
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
         fflush(stderr);
         return(-1);
      }
   
      /* --------------------------------------- */
      /* Read in the PRIMME configuration file   */
      /* --------------------------------------- */
      if (read_solver_params(SolverConfigFileName, driver.outputFileName, 
                           &primme, &method) < 0) {
         fprintf(stderr, "Reading solver parameters failed\n");
         return(-1);
      }
   }

#ifdef USE_MPI
   /* ------------------------------------------------- */
   /* Send read common primme members to all processors */ 
   /* Setup the primme members local to this processor  */ 
   /* ------------------------------------------------- */
   broadCast(&primme, &method, &driver, master, comm);
#endif

   /* --------------------------------------- */
   /* Set up matrix vector and preconditioner */
   /* --------------------------------------- */
   if (setMatrixAndPrecond(&driver, &primme, &permutation) != 0) return -1;

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
   if (master) {
      fprintf(primme.outputFile,"PRIMME will allocate the following memory:\n");
      fprintf(primme.outputFile," processor %d, real workspace, %ld bytes\n",
                                      procID, primme.realWorkSize);
      fprintf(primme.outputFile," processor %d, int  workspace, %d bytes\n",
                                      procID, primme.intWorkSize);
   }

   /* --------------------------------------- */
   /* Display given parameter configuration   */
   /* Place this after the dprimme() to see   */
   /* any changes dprimme() made to PRIMME    */
   /* --------------------------------------- */

   if (master) {
      driver_display_params(driver, primme.outputFile); 
      primme_display_params(primme);
      driver_display_method(method, primme.outputFile);
   }

   /* --------------------------------------------------------------------- */
   /*                            Run the d/zprimme solver                   */
   /* --------------------------------------------------------------------- */

   /* Allocate space for converged Ritz values and residual norms */

   evals = (double *)primme_calloc(primme.numEvals, sizeof(double), "evals");
   evecs = (PRIMME_NUM *)primme_calloc(primme.nLocal*primme.numEvals, 
                                sizeof(PRIMME_NUM), "evecs");
   rnorms = (double *)primme_calloc(primme.numEvals, sizeof(double), "rnorms");

   /* ------------------------ */
   /* Initial guess (optional) */
   /* ------------------------ */

   /* Read initial guess from a file */
   if (driver.initialGuessesFileName[0] && primme.initSize+primme.numOrthoConst > 0) {
      int cols, i=0;
      ASSERT_MSG(readBinaryEvecsAndPrimmeParams(driver.initialGuessesFileName, evecs, NULL, primme.n,
                                                min(primme.initSize+primme.numOrthoConst, primme.numEvals),
                                                &cols, primme.nLocal, permutation, &primme) != 0, 1, "");
      primme.numOrthoConst = min(primme.numOrthoConst, cols);

      /* Perturb the initial guesses by a vector with some norm  */
      if (driver.initialGuessesPert > 0) {
         PRIMME_NUM *r = (PRIMME_NUM *)primme_calloc(primme.nLocal,sizeof(PRIMME_NUM), "random");
         double norm;
         int j;
         for (i=primme.numOrthoConst; i<min(cols, primme.initSize+primme.numOrthoConst); i++) {
            SUF(Num_larnv)(2, primme.iseed, primme.nLocal, COMPLEXZ(r));
            norm = sqrt(REAL_PARTZ(SUF(Num_dot)(primme.nLocal, COMPLEXZ(r), 1, COMPLEXZ(r), 1)));
            for (j=0; j<primme.nLocal; j++)
               evecs[primme.nLocal*i+j] += r[j]/norm*driver.initialGuessesPert;
         }
         free(r);
      }
      SUF(Num_larnv)(2, primme.iseed, (primme.initSize+primme.numOrthoConst-i)*primme.nLocal,
                     COMPLEXZ(&evecs[primme.nLocal*i]));
   } else if (primme.numOrthoConst > 0) {
      ASSERT_MSG(0, 1, "numOrthoConst > 0 but no value in initialGuessesFileName.\n");
   } else if (primme.initSize > 0) {
      SUF(Num_larnv)(2, primme.iseed, primme.initSize*primme.nLocal, COMPLEXZ(evecs));
   } else {
      SUF(Num_larnv)(2, primme.iseed, primme.nLocal, COMPLEXZ(evecs));
   }


   /* ------------- */
   /*  Call primme  */
   /* ------------- */

   wt1 = primme_get_wtime(); 
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
   primme_get_time(&ut1,&st1);
#endif

   ret = PREFIX(primme)(evals, COMPLEXZ(evecs), rnorms, &primme);

   wt2 = primme_get_wtime();
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
   primme_get_time(&ut2,&st2);
#endif

   if (driver.checkXFileName[0]) {
      retX = check_solution(driver.checkXFileName, &primme, evals, evecs, rnorms, permutation);
   }

   /* --------------------------------------------------------------------- */
   /* Save evecs and primme params  (optional)                              */
   /* --------------------------------------------------------------------- */
   if (driver.saveXFileName[0]) {
      ASSERT_MSG(writeBinaryEvecsAndPrimmeParams(driver.saveXFileName, evecs, permutation, &primme) == 0, 1, "");
   }

   /* --------------------------------------------------------------------- */
   /* Reporting                                                             */
   /* --------------------------------------------------------------------- */

   if (master) {
      primme_PrintStackTrace(primme);

      fprintf(primme.outputFile, "Wallclock Runtime   : %-f\n", wt2-wt1);
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
      fprintf(primme.outputFile, "User Time           : %f seconds\n", ut2-ut1);
      fprintf(primme.outputFile, "Syst Time           : %f seconds\n", st2-st1);
#endif

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

   fclose(primme.outputFile);
   destroyMatrixAndPrecond(&driver, &primme, permutation);
   primme_Free(&primme);
   free(evals);
   free(evecs);
   free(rnorms);

   if (ret != 0 && master) {
      fprintf(primme.outputFile, 
         "Error: dprimme returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   if (retX != 0 && master) {
      fprintf(primme.outputFile, 
         "Error: found some issues in the solution return by dprimme\n");
      return -1;
   }

  return(0);
}
/******************************************************************************/
/* END OF MAIN DRIVER FUNCTION                                                */
/******************************************************************************/

/******************************************************************************/
/* Matvec, preconditioner and other utilities                                 */

/******************************************************************************/

#ifdef USE_MPI
/******************************************************************************
 * Function to broadcast the primme data structure to all processors
 *
 * EXCEPTIONS: procID and seed[] are not copied from processor 0. 
 *             Each process creates their own.
******************************************************************************/
static void broadCast(primme_params *primme, primme_preset_method *method, 
   driver_params *driver, int master, MPI_Comm comm){

   int i;

   MPI_Bcast(driver->outputFileName, 512, MPI_CHAR, 0, comm);
   MPI_Bcast(driver->matrixFileName, 1024, MPI_CHAR, 0, comm);
   MPI_Bcast(driver->initialGuessesFileName, 1024, MPI_CHAR, 0, comm);
   MPI_Bcast(driver->saveXFileName, 1024, MPI_CHAR, 0, comm);
   MPI_Bcast(driver->checkXFileName, 1024, MPI_CHAR, 0, comm);
   MPI_Bcast(&driver->initialGuessesPert, 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&driver->matrixChoice, 1, MPI_INT, 0, comm);
   MPI_Bcast(&driver->PrecChoice, 1, MPI_INT, 0, comm);
   MPI_Bcast(&driver->isymm, 1, MPI_INT, 0, comm);
   MPI_Bcast(&driver->level, 1, MPI_INT, 0, comm);
   MPI_Bcast(&driver->threshold, 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&driver->filter, 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&driver->shift, 1, MPI_DOUBLE, 0, comm);

   MPI_Bcast(&(primme->numEvals), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->target), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->numTargetShifts), 1, MPI_INT, 0, comm);

   if (primme->numTargetShifts > 0 && !master) {
      primme->targetShifts = (double *)primme_calloc(
         primme->numTargetShifts, sizeof(double), "targetShifts");
   }
   for (i=0; i<primme->numTargetShifts; i++) {
      MPI_Bcast(&(primme->targetShifts[i]), 1, MPI_DOUBLE, 0, comm);
   }
   MPI_Bcast(&(primme->locking), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->dynamicMethodSwitch), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->initSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->numOrthoConst), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxBasisSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->minRestartSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxBlockSize), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxMatvecs), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->maxOuterIterations), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->aNorm), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme->eps), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme->printLevel), 1, MPI_INT, 0, comm);

   MPI_Bcast(&(primme->restartingParams.scheme), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->restartingParams.maxPrevRetain), 1, MPI_INT, 0, comm);

   MPI_Bcast(&(primme->correctionParams.precondition), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->correctionParams.robustShifts), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->correctionParams.maxInnerIterations),1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.convTest), 1, MPI_INT, 0, comm);
   MPI_Bcast(&(primme->correctionParams.relTolBase), 1, MPI_DOUBLE, 0, comm);
   MPI_Bcast(&(primme->correctionParams.projectors.LeftQ),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.LeftX),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.RightQ), 1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.RightX), 1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.SkewQ),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.SkewX),  1, MPI_INT, 0,comm);
   MPI_Bcast(&(primme->correctionParams.projectors.SkewX),  1, MPI_INT, 0,comm);

   MPI_Bcast(method, 1, MPI_INT, 0, comm);
}
#endif

static int setMatrixAndPrecond(driver_params *driver, primme_params *primme, int **permutation) {
   int numProcs=1;

#  if defined(USE_MPI)
   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   primme->commInfo = (MPI_Comm *)primme_calloc(1, sizeof(MPI_Comm), "MPI_Comm");
   *(MPI_Comm*)primme->commInfo = MPI_COMM_WORLD;
#  endif

   if (driver->matrixChoice == driver_default) {
      if (numProcs <= 1) {
         driver->matrixChoice = driver_native;
      } else {
#        ifdef USE_PETSC
            driver->matrixChoice = driver_petsc;
#        else
            driver->matrixChoice = driver_parasails;
#        endif
      }
   }
   switch(driver->matrixChoice) {
   case driver_default:
      assert(0);
      break;
   case driver_native:
#if !defined(USE_NATIVE)
      fprintf(stderr, "ERROR: NATIVE is needed!\n");
      return -1;
#else
#  if defined(USE_MPI)
      if (numProcs != 1) {
         fprintf(stderr, "ERROR: MPI is not supported with NATIVE, use other!\n");
         return -1;
      }
      *(MPI_Comm*)primme->commInfo = MPI_COMM_WORLD;
#  endif
      {
         CSRMatrix *matrix, *prec;
         double *diag;
         
         if (readMatrixNative(driver->matrixFileName, &matrix, &primme->aNorm) !=0 )
            return -1;
         primme->matrix = matrix;
         primme->matrixMatvec = CSRMatrixMatvec;
         primme->n = primme->nLocal = matrix->n;
         switch(driver->PrecChoice) {
         case driver_noprecond:
            primme->preconditioner = NULL;
            primme->applyPreconditioner = NULL;
            break;
         case driver_jacobi:
            createInvDiagPrecNative(matrix, driver->shift, &diag);
            primme->preconditioner = diag;
            primme->applyPreconditioner = ApplyInvDiagPrecNative;
            break;
         case driver_jacobi_i:
            createInvDavidsonDiagPrecNative(matrix, &diag);
            primme->preconditioner = diag;
            primme->applyPreconditioner = ApplyInvDavidsonDiagPrecNative;
         case driver_ilut:
            createILUTPrecNative(matrix, driver->shift, driver->level, driver->threshold,
                                 driver->filter, &prec);
            primme->preconditioner = prec;
            primme->applyPreconditioner = ApplyILUTPrecNative;
            break;
         }
      }
#endif
      break;

   case driver_parasails:
#if !defined(USE_PARASAILS)
      fprintf(stderr, "ERROR: ParaSails is needed!\n");
      return -1;
      if (driver->PrecChoice != driver_ilut) {
         fprintf(stderr, "ERROR: ParaSails only supports ILUT preconditioner!\n");
         return -1;
      }
#else
      {
         Matrix *matrix;
         ParaSails *precond=NULL;
         int m, mLocal;
         readMatrixAndPrecondParaSails(driver->matrixFileName, driver->shift, driver->level,
               driver->threshold, driver->filter, driver->isymm, MPI_COMM_WORLD, &primme->aNorm,
               &primme->n, &m, &primme->nLocal, &mLocal, &primme->numProcs, &primme->procID, &matrix,
               (driver->PrecChoice == driver_ilut) ? &precond : NULL);
         *(MPI_Comm*)primme->commInfo = MPI_COMM_WORLD;
         primme->matrix = matrix;
         primme->matrixMatvec = ParaSailsMatrixMatvec;
         primme->preconditioner = precond;
         primme->applyPreconditioner = precond ? ApplyPrecParaSails : NULL;
      }
#endif
      break;

   case driver_petsc:
#ifndef USE_PETSC
      fprintf(stderr, "ERROR: PETSc is needed!\n");
      return -1;
#else
      {
         PetscErrorCode ierr;
         Mat *matrix;
         PC *pc;
         Vec *vec;
         int m, mLocal;
         if (readMatrixPetsc(driver->matrixFileName, &primme->n, &m, &primme->nLocal, &mLocal,
                         &primme->numProcs, &primme->procID, &matrix, &primme->aNorm, permutation) != 0)
            return -1;
         *(MPI_Comm*)primme->commInfo = PETSC_COMM_WORLD;
         primme->matrix = matrix;
         primme->matrixMatvec = PETScMatvec;
         if (driver->PrecChoice == driver_noprecond) {
            primme->preconditioner = NULL;
            primme->applyPreconditioner = NULL;
         }
         else if (driver->PrecChoice != driver_jacobi_i) {
            pc = (PC *)primme_calloc(1, sizeof(PC), "pc");
            ierr = PCCreate(PETSC_COMM_WORLD, pc); CHKERRQ(ierr);
            if (driver->PrecChoice == driver_jacobi) {
               ierr = PCSetType(*pc, PCJACOBI); CHKERRQ(ierr);
            }
            else if (driver->PrecChoice == driver_ilut) {
               if (primme->numProcs <= 1) {
                  ierr = PCSetType(*pc, PCICC); CHKERRQ(ierr);
               }
               else {
                  #ifdef PETSC_HAVE_HYPRE
                     ierr = PCSetType(pc, PCHYPRE); CHKERRQ(ierr);
                     ierr = PCHYPRESetType(*pc, "parasails"); CHKERRQ(ierr);
                  #else
                     ierr = PCSetType(*pc, PCBJACOBI); CHKERRQ(ierr);
                  #endif
               }
            }

            ierr = PCSetOperators(*pc, *matrix, *matrix); CHKERRQ(ierr);
            ierr = PCSetFromOptions(*pc); CHKERRQ(ierr);
            ierr = PCSetUp(*pc); CHKERRQ(ierr);
            primme->preconditioner = pc;
            primme->applyPreconditioner = ApplyPCPrecPETSC;
         }
         else {
            vec = (Vec *)primme_calloc(1, sizeof(Vec), "Vec preconditioner");
            ierr = MatCreateVecs(*matrix, vec, NULL); CHKERRQ(ierr);
            ierr = MatGetDiagonal(*matrix, *vec); CHKERRQ(ierr);
            primme->preconditioner = vec;
            primme->applyPreconditioner = ApplyPCPrecPETSC;
         }
      }
#endif
      break;
   }

#if defined(USE_MPI)
   primme->globalSumDouble = par_GlobalSumDouble;
#endif
   return 0;
}

static int destroyMatrixAndPrecond(driver_params *driver, primme_params *primme, int *permutation) {
   switch(driver->matrixChoice) {
   case driver_default:
      assert(0);
      break;
   case driver_native:
#if !defined(USE_NATIVE)
      fprintf(stderr, "ERROR: NATIVE is needed!\n");
      return -1;
#else
      free(((CSRMatrix*)primme->matrix)->AElts);
      free(((CSRMatrix*)primme->matrix)->IA);
      free(((CSRMatrix*)primme->matrix)->JA);
      free(primme->matrix);

      switch(driver->PrecChoice) {
      case driver_noprecond:
         break;
      case driver_jacobi:
      case driver_jacobi_i:
         free(primme->preconditioner);
         break;
      case driver_ilut:
         if (primme->preconditioner) {
            free(((CSRMatrix*)primme->preconditioner)->AElts);
            free(((CSRMatrix*)primme->preconditioner)->IA);
            free(((CSRMatrix*)primme->preconditioner)->JA);
            free(primme->preconditioner);
         }
         break;
      }
#endif
      break;

   case driver_parasails:
#if !defined(USE_PARASAILS)
      fprintf(stderr, "ERROR: ParaSails is needed!\n");
      return -1;
      if (driver->PrecChoice != driver_ilut) {
         fprintf(stderr, "ERROR: ParaSails only supports ILUT preconditioner!\n");
         return -1;
      }
#else
      /* TODO: destroy ParaSail matrices */

#endif
      break;

   case driver_petsc:
#ifndef USE_PETSC
      fprintf(stderr, "ERROR: PETSc is needed!\n");
      return -1;
#else
      {
         PetscErrorCode ierr;
         ierr = MatDestroy((Mat*)primme->matrix);CHKERRQ(ierr);
         if (primme->preconditioner) {
         }
         if (driver->PrecChoice == driver_noprecond) {
         }
         else if (driver->PrecChoice != driver_jacobi_i) {
            ierr = PCDestroy((PC*)primme->preconditioner);CHKERRQ(ierr);
            free(primme->preconditioner);
         }
         else {
            ierr = VecDestroy((Vec*)primme->preconditioner);CHKERRQ(ierr);
            free(primme->preconditioner);
         }
      }
#endif
      break;
   }
#if defined(USE_MPI)
   free(primme->commInfo);
#endif
   if (permutation) free(permutation);
   return 0;
}


#ifdef USE_MPI
/******************************************************************************
 * MPI globalSumDouble function
 *
******************************************************************************/
static void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme) {
   MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;

   MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator);
}
#endif

#undef __FUNCT__
#define __FUNCT__ "check_solution"
static int check_solution(const char *checkXFileName, primme_params *primme, double *evals,
                   PRIMME_NUM *evecs, double *rnorms, int *perm) {

   double eval0, rnorm0, prod, auxd;
   PRIMME_NUM *Ax, *r, *X=NULL, *h, *h0;
   int i, j, cols, retX=0, one=1;
   primme_params primme0;

   /* Read stored eigenvectors and primme_params */
   ASSERT_MSG(readBinaryEvecsAndPrimmeParams(checkXFileName, NULL, &X, primme->n, primme->n, &cols,
                                             primme->nLocal, perm, &primme0) == 0, -1, "");
   /* Check primme_params */
#  define CHECK_PRIMME_PARAM(F) \
        if (primme0. F != primme-> F ) { \
           fprintf(stderr, "Warning: discrepancy in primme." #F ", %d should be close to %d\n", primme-> F , primme0. F ); \
           retX = 1; \
        }
#  define CHECK_PRIMME_PARAM_DOUBLE(F) \
        if (fabs(primme0. F - primme-> F) > primme-> F * 1e-14) { \
           fprintf(stderr, "Warning: discrepancy in primme." #F ", %.16e should be close to %.16e\n", primme-> F , primme0. F ); \
           retX = 1; \
        }
#  define CHECK_PRIMME_PARAM_TOL(F, T) \
        if (abs(primme0. F - primme-> F ) > primme-> F * T /100+1) { \
           fprintf(stderr, "Warning: discrepancy in primme." #F ", %d should be close to %d\n", primme-> F , primme0. F ); \
           retX = 1; \
        }

   if (primme0.n) {
      CHECK_PRIMME_PARAM(n);
      CHECK_PRIMME_PARAM(numEvals);
      CHECK_PRIMME_PARAM(target);
      CHECK_PRIMME_PARAM(numTargetShifts);
      CHECK_PRIMME_PARAM(dynamicMethodSwitch);
      CHECK_PRIMME_PARAM(locking);
      CHECK_PRIMME_PARAM(numOrthoConst);
      CHECK_PRIMME_PARAM(maxBasisSize);
      CHECK_PRIMME_PARAM(minRestartSize);
      CHECK_PRIMME_PARAM(restartingParams.scheme);
      CHECK_PRIMME_PARAM(restartingParams.maxPrevRetain);
      CHECK_PRIMME_PARAM(correctionParams.precondition);
      CHECK_PRIMME_PARAM(correctionParams.robustShifts);
      CHECK_PRIMME_PARAM(correctionParams.maxInnerIterations);
      CHECK_PRIMME_PARAM(correctionParams.projectors.LeftQ);
      CHECK_PRIMME_PARAM(correctionParams.projectors.LeftX);
      CHECK_PRIMME_PARAM(correctionParams.projectors.RightQ);
      CHECK_PRIMME_PARAM(correctionParams.projectors.RightX);
      CHECK_PRIMME_PARAM(correctionParams.projectors.SkewQ);
      CHECK_PRIMME_PARAM(correctionParams.projectors.SkewX);
      CHECK_PRIMME_PARAM(correctionParams.convTest);
      CHECK_PRIMME_PARAM_DOUBLE(aNorm);
      CHECK_PRIMME_PARAM_DOUBLE(eps);
      CHECK_PRIMME_PARAM_DOUBLE(correctionParams.relTolBase);
      CHECK_PRIMME_PARAM(initSize);
      CHECK_PRIMME_PARAM_TOL(stats.numOuterIterations, 40);
   }

   h = (PRIMME_NUM *)primme_calloc(cols*2, sizeof(PRIMME_NUM), "h"); h0 = &h[cols];
   Ax = (PRIMME_NUM *)primme_calloc(primme->nLocal, sizeof(PRIMME_NUM), "Ax");
   r = (PRIMME_NUM *)primme_calloc(primme->nLocal, sizeof(PRIMME_NUM), "r");
   
   for (i=0; i < primme->initSize; i++) {
      /* Check |V(:,i)'A*V(:,i) - evals[i]| < |r|*|A| */
      primme->matrixMatvec(&evecs[primme->nLocal*i], Ax, &one, primme);
      auxd = REAL_PARTZ(SUF(Num_dot)(primme->nLocal, COMPLEXZ(&evecs[primme->nLocal*i]), 1, COMPLEXZ(Ax), 1));
      if (primme->globalSumDouble) primme->globalSumDouble(&auxd, &eval0, &one, primme);
      else eval0 = auxd;
      if (fabs(evals[i] - eval0) > rnorms[i]*primme->aNorm && primme->procID == 0) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E should be close to %-22.1E\n", i, evals[i], eval0);
         retX = 1;
      }
      /* Check |A*V(:,i) - (V(:,i)'A*V(:,i))*V(:,i)| < |r| */
      for (j=0; j<primme->nLocal; j++) r[j] = Ax[j] - evals[i]*evecs[primme->nLocal*i+j];
      auxd = REAL_PARTZ(SUF(Num_dot)(primme->nLocal, COMPLEXZ(r), 1, COMPLEXZ(r), 1));
      if (primme->globalSumDouble) primme->globalSumDouble(&auxd, &rnorm0, &one, primme);
      else rnorm0 = auxd;
      rnorm0 = sqrt(rnorm0);
      if (fabs(rnorms[i]-rnorm0) > 4*max(primme->aNorm,fabs(evals[i]))*MACHINE_EPSILON && primme->procID == 0) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E, residual | %5E - %5E | <= %5E\n", i, evals[i], rnorms[i], rnorm0, 4*max(primme->aNorm,fabs(evals[i]))*MACHINE_EPSILON);
         retX = 1;
      }
      if (rnorm0 > primme->eps*primme->aNorm*sqrt((double)(i+1)) && primme->procID == 0) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E, RR residual %5E is larger than tolerance %5E\n", i, evals[i], rnorm0, primme->eps*primme->aNorm*sqrt((double)(i+1)));
         retX = 1;
      }
      /* Check X'V(:,i) >= sqrt(1-2|r|), assuming residual of X is less than the tolerance */
      SUF(Num_gemv)("C", primme->nLocal, cols, COMPLEXZV(1.0), COMPLEXZ(X), primme->nLocal, COMPLEXZ(&evecs[primme->nLocal*i]), 1, COMPLEXZV(0.), COMPLEXZ(h), 1);
      if (primme->globalSumDouble) {
         int cols0 = cols*sizeof(PRIMME_NUM)/sizeof(double);
         primme->globalSumDouble(h, h0, &cols0, primme);
      }
      else h0 = h;
      prod = REAL_PARTZ(SUF(Num_dot)(cols, COMPLEXZ(h0), 1, COMPLEXZ(h0), 1));
      if (prod < sqrt(1.-2.*rnorms[i]) && primme->procID == 0) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E not found on X, %5E > %5E\n", i, evals[i], prod, sqrt(1.-2.*rnorms[i]));
         retX = 1;
      }
   }
   free(h);
   free(X);
   free(r);
   free(Ax);

   return retX; 
}

#undef __FUNCT__
#define __FUNCT__ "readBinaryEvecsAndPrimmeParams"
static int readBinaryEvecsAndPrimmeParams(const char *fileName, PRIMME_NUM *X, PRIMME_NUM **Xout,
                                          int n, int Xcols, int *Xcolsout, int nLocal,
                                          int *perm, primme_params *primme_out) {

#  define FREAD(A, B, C, D) { ASSERT_MSG(fread(A, B, C, D) == (size_t)C, -1, "Unexpected end of file\n"); }

   FILE *f;
   PRIMME_NUM d;
   int i, j, cols;

   ASSERT_MSG((f = fopen(fileName, "rb")),
                  -1, "Could not open file %s\n", fileName);

   /* Check number size */
   /* NOTE: 2*IMAGINARY*IMAGINARY+1 is -1 in complex arith and 1 in real arith */
   FREAD(&d, sizeof(d), 1, f);
   ASSERT_MSG((int)(REAL_PART(d*(2.*IMAGINARY*IMAGINARY + 1.))) == (int)sizeof(d),
                  -1, "Mismatch arithmetic in file %s\n", fileName);
   /* Check matrix size */
   FREAD(&d, sizeof(d), 1, f);
   ASSERT_MSG(((int)REAL_PART(d)) == n,
                  -1, "Mismatch matrix size in file %s\n", fileName);

   /* Read X */
   FREAD(&d, sizeof(d), 1, f); cols = REAL_PART(d);
   if (Xcols > 0 && (X || Xout)) {
      if (!X) *Xout = X = (PRIMME_NUM*)malloc(sizeof(PRIMME_NUM)*min(cols, Xcols)*nLocal);
      if (Xcolsout) *Xcolsout = min(cols, Xcols);
      if (!perm) {
         for (i=0; i<min(cols, Xcols); i++) {
            fseek(f, (i*n + 3)*sizeof(d), SEEK_SET);
            FREAD(&X[nLocal*i], sizeof(d), nLocal, f);
         }
      }
      else {
         for (i=0; i<min(cols, Xcols); i++) {
            for (j=0; j<nLocal; j++) {
               fseek(f, (i*n + perm[j] + 3)*sizeof(d), SEEK_SET);
               FREAD(&X[nLocal*i+j], sizeof(d), 1, f);
            }
         }
      }
   }
   fseek(f, (cols*n + 3)*sizeof(d), SEEK_SET);

   /* Read primme_params */
   if (primme_out) {
      FREAD(&d, sizeof(d), 1, f);
      if ((int)REAL_PART(d) == (int)sizeof(*primme_out)) {
         FREAD(primme_out, sizeof(*primme_out), 1, f);
      }
      else
         primme_out->n = 0;
   }

   fclose(f);
   return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeBinaryEvecsAndPrimmeParams"
static int writeBinaryEvecsAndPrimmeParams(const char *fileName, PRIMME_NUM *X, int *perm,
                                           primme_params *primme) {

#  define FWRITE(A, B, C, D) { ASSERT_MSG(fwrite(A, B, C, D) == (size_t)C, -1, "Unexpected error writing on %s\n", fileName); }

   FILE *f;
   PRIMME_NUM d;
   int i, j;

   ASSERT_MSG((f = fopen(fileName, "wb")),
                  -1, "Could not open file %s\n", fileName);

   /* Write number size */
   if (primme->procID == 0) {
      /* NOTE: 2*IMAGINARY*IMAGINARY+1 is -1 in complex arith and 1 in real arith */
      d = (2.*IMAGINARY*IMAGINARY + 1.)*sizeof(d);
      FWRITE(&d, sizeof(d), 1, f);
      /* Write matrix size */
      d = primme->n;
      FWRITE(&d, sizeof(d), 1, f);
      /* Write number of columns */
      d = primme->initSize;
      FWRITE(&d, sizeof(d), 1, f);
   }

   /* Write X */
   if (!perm) {
      for (i=0; i<primme->initSize; i++) {
         fseek(f, (i*primme->n + 3)*sizeof(d), SEEK_SET);
         FWRITE(&X[primme->nLocal*i], sizeof(d), primme->nLocal, f);
      }
   }
   else {
      for (i=0; i<primme->initSize; i++) {
         for (j=0; j<primme->nLocal; j++) {
            fseek(f, (i*primme->n + perm[j] + 3)*sizeof(d), SEEK_SET);
            FWRITE(&X[primme->nLocal*i+j], sizeof(d), 1, f);
         }
      }
   }

   /* Write primme_params */
   if (primme->procID == 0) {
      fseek(f, sizeof(d)*(primme->n*primme->initSize + 3), SEEK_SET);
      d = sizeof(*primme);
      FWRITE(&d, sizeof(d), 1, f);
      FWRITE(primme, sizeof(*primme), 1, f);
   }

   fclose(f);
   return 0;
}
