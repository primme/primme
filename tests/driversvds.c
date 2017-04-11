/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2017 College of William & Mary,
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
 *           singular values and vectors using PRIMME.
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
 *  SolverConfigFileName  includes all d/zprimme_svds required information
 *                            as stored in primme_svds_params data structure.
 *
 *             Example files: FullConfSvds  Full customization of PRIMME
 *                            LeanConfSvds  Use a preset method and some customization
 *                            MinConfSvds   Provide ONLY a preset method and numSvals.
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
PetscLogEvent PRIMME_GLOBAL_SUM;
#endif
#ifdef USE_RSB
#  include "rsbw.h"
#endif

/* primme.h header file is required to run primme */
#include "primme_svds.h"
#include "shared_utils.h"
#include "ioandtest.h"
/* wtime.h header file is included so primme's timimg functions can be used */
#include "../src/include/wtime.h"

static int real_main (int argc, char *argv[]);
static int setMatrixAndPrecond(driver_params *driver, primme_svds_params *primme_svds, int **permutation);
static int destroyMatrixAndPrecond(driver_params *driver, primme_svds_params *primme_svds, int *permutation);



int main (int argc, char *argv[]) {
   int ret;
#if defined(USE_PETSC)
   PetscInt ierr;
#endif

#if defined(USE_MPI) && !defined(USE_PETSC)
   MPI_Init(&argc, &argv);
#elif defined(USE_PETSC)
   PetscInitialize(&argc, &argv, NULL, NULL);
   PetscLogEventRegister("PRIMME global sum", 0, &PRIMME_GLOBAL_SUM);
   #if PETSC_VERSION_LT(3,7,0)
   ierr = PetscLogBegin(); CHKERRQ(ierr);
   #else
   ierr = PetscLogDefaultBegin(); CHKERRQ(ierr);
   #endif
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
   double *svals, *rnorms;
   SCALAR *svecs;
   driver_params driver;
   primme_svds_params primme_svds;
   primme_svds_preset_method method=primme_svds_default;
   primme_preset_method primmemethod=PRIMME_DEFAULT_METHOD, primmemethodStage2=PRIMME_DEFAULT_METHOD;
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

   primme_svds_initialize(&primme_svds);

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
      /* Read in the PRIMME SVDS configuration file   */
      /* --------------------------------------- */
      if (read_solver_params_svds(SolverConfigFileName, driver.outputFileName, 
                           &primme_svds, "primme_svds.", &method, "method", &primmemethod,
                           &primmemethodStage2) < 0) {
         fprintf(stderr, "Reading solver parameters failed\n");
         return(-1);
      }
   }

#ifdef USE_MPI
   /* ------------------------------------------------- */
   /* Send read common primme members to all processors */ 
   /* Setup the primme members local to this processor  */ 
   /* ------------------------------------------------- */
   broadCast_svds(&primme_svds, &method, &primmemethod, &primmemethodStage2, &driver, master, comm);
#endif

   /* --------------------------------------- */
   /* Set up matrix vector and preconditioner */
   /* --------------------------------------- */
   if (setMatrixAndPrecond(&driver, &primme_svds, &permutation) != 0) return -1;

   /* --------------------------------------- */
   /* Pick one of the default methods(if set) */
   /* --------------------------------------- */
   primme_svds_set_method(method, primmemethod, primmemethodStage2, &primme_svds);

#ifdef NOT_USE_ALIGNMENT
   /* --------------------------------------- */
   /* Set alignment (optional)                */
   /* --------------------------------------- */
   if (primme_svds.method == primme_svds_op_AtA) {
      primme_svds.primme.ldOPs = primme_svds.nLocal;
   }
   else if (primme_svds.method == primme_svds_op_AAt) {
      primme_svds.primme.ldOPs = primme_svds.mLocal;
   }
   else {
      primme_svds.primme.ldOPs = primme_svds.mLocal + primme_svds.nLocal;
   }
   if (primme_svds.methodStage2 == primme_svds_op_augmented) {
      primme_svds.primmeStage2.ldOPs = primme_svds.mLocal + primme_svds.nLocal;
   }
#endif

   /* --------------------------------------- */
   /* Optional: report memory requirements    */
   /* --------------------------------------- */

   ret = Sprimme_svds(NULL,NULL,NULL,&primme_svds);
   if (master) {
      fprintf(primme_svds.outputFile,"PRIMME SVDS will allocate the following memory:\n");
      fprintf(primme_svds.outputFile," processor %d, real workspace, %ld bytes\n",
                                      procID, primme_svds.realWorkSize);
      fprintf(primme_svds.outputFile," processor %d, int  workspace, %d bytes\n",
                                      procID, primme_svds.intWorkSize);
   }

   /* ---------------------------------------- */
   /* Display given parameter configuration    */
   /* Place this after the dprimme() to see    */
   /* any changes primme_svds() made to PRIMME */
   /* ---------------------------------------- */

   if (master) {
      driver_display_params(driver, primme_svds.outputFile); 
      primme_svds_display_params(primme_svds);
      driver_display_methodsvd(method, "method", primme_svds.outputFile);
      driver_display_method(primmemethod, "primme.method", primme_svds.outputFile);
      driver_display_method(primmemethodStage2, "primmeStage2.method", primme_svds.outputFile);
   }

   /* --------------------------------------------------------------------- */
   /*                            Run the d/zprimme_svds solver                   */
   /* --------------------------------------------------------------------- */

   /* Allocate space for converged Ritz values and residual norms */

   svals = (double *)primme_calloc(primme_svds.numSvals, sizeof(double), "svals");
   svecs = (SCALAR *)primme_calloc((primme_svds.mLocal+primme_svds.nLocal)*
                                       primme_svds.numSvals, sizeof(SCALAR), "svecs");
   rnorms = (double *)primme_calloc(primme_svds.numSvals, sizeof(double), "rnorms");

   /* ------------------------ */
   /* Initial guess (optional) */
   /* ------------------------ */

   /* Read initial guess from a file */
   if (driver.initialGuessesFileName[0] && primme_svds.initSize+primme_svds.numOrthoConst > 0) {
      int cols, i=0, n;
      ASSERT_MSG(readBinaryEvecsAndPrimmeSvdsParams(
         driver.initialGuessesFileName, svecs, NULL, primme_svds.m, primme_svds.n,
         min(primme_svds.initSize+primme_svds.numOrthoConst,
             primme_svds.numSvals), &cols, primme_svds.mLocal, primme_svds.nLocal,
         permutation, &primme_svds) != 0, 1, "");
      primme_svds.numOrthoConst = min(primme_svds.numOrthoConst, cols);
      primme_svds.initSize = min(primme_svds.initSize, cols - primme_svds.numOrthoConst);
      n = primme_svds.initSize+primme_svds.numOrthoConst;

      /* Perturb the initial guesses by a vector with some norm  */
      if (driver.initialGuessesPert > 0) {
         SCALAR *r = (SCALAR *)primme_calloc(max(primme_svds.nLocal,primme_svds.mLocal),sizeof(SCALAR), "random");
         double norm;
         int j;
         assert(primme_svds.numProcs <= 1);
         for (i=primme_svds.numOrthoConst; i<min(cols, primme_svds.initSize+primme_svds.numOrthoConst); i++) {
            Num_larnv_Sprimme(2, primme_svds.iseed, primme_svds.mLocal, r);
            norm = sqrt(REAL_PART(Num_dot_Sprimme(primme_svds.mLocal, r, 1, r, 1)));
            for (j=0; j<primme_svds.mLocal; j++)
               svecs[primme_svds.mLocal*i+j] += r[j]/norm*driver.initialGuessesPert;
         }
         for (i=primme_svds.numOrthoConst; i<min(cols, primme_svds.initSize+primme_svds.numOrthoConst); i++) {
            Num_larnv_Sprimme(2, primme_svds.iseed, primme_svds.nLocal, r);
            norm = sqrt(REAL_PART(Num_dot_Sprimme(primme_svds.nLocal, r, 1, r, 1)));
            for (j=0; j<primme_svds.nLocal; j++)
               svecs[primme_svds.mLocal*n+primme_svds.nLocal*i+j] += r[j]/norm*driver.initialGuessesPert;
         }
         free(r);
      }
      Num_larnv_Sprimme(2, primme_svds.iseed, (primme_svds.initSize+primme_svds.numOrthoConst-i)*primme_svds.mLocal,
                     &svecs[primme_svds.mLocal*i]);
      Num_larnv_Sprimme(2, primme_svds.iseed, (primme_svds.initSize+primme_svds.numOrthoConst-i)*primme_svds.mLocal,
                     &svecs[primme_svds.mLocal*n+primme_svds.nLocal*i]);
   } else if (primme_svds.numOrthoConst > 0) {
      ASSERT_MSG(0, 1, "numOrthoConst > 0 but no value in initialGuessesFileName.\n");
   } else if (primme_svds.initSize > 0) {
      Num_larnv_Sprimme(2, primme_svds.iseed, primme_svds.initSize*(primme_svds.mLocal+primme_svds.nLocal),
                     svecs);
   }


   /* ------------------ */
   /*  Call svds_primme  */
   /* ------------------ */

   wt1 = primme_get_wtime(); 
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
   primme_get_time(&ut1,&st1);
#endif

   ret = Sprimme_svds(svals, svecs, rnorms, &primme_svds);

   wt2 = primme_get_wtime();
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
   primme_get_time(&ut2,&st2);
#endif

   if (driver.checkXFileName[0]) {
      retX = check_solution_svds(driver.checkXFileName, &primme_svds, svals, svecs, rnorms, permutation);
   }

   /* --------------------------------------------------------------------- */
   /* Save svecs and primme_svds_params  (optional)                         */
   /* --------------------------------------------------------------------- */
   if (driver.saveXFileName[0]) {
      ASSERT_MSG(writeBinaryEvecsAndPrimmeSvdsParams(driver.saveXFileName,
         svecs, permutation, &primme_svds) == 0, 1, "");
   }

   /* --------------------------------------------------------------------- */
   /* Reporting                                                             */
   /* --------------------------------------------------------------------- */

   if (master) {
      fprintf(primme_svds.outputFile, "Wallclock Runtime   : %-f\n", wt2-wt1);
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
      fprintf(primme_svds.outputFile, "User Time           : %f seconds\n", ut2-ut1);
      fprintf(primme_svds.outputFile, "Syst Time           : %f seconds\n", st2-st1);
#endif

      for (i=0; i < primme_svds.numSvals; i++) {
         fprintf(primme_svds.outputFile, "Sval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
            svals[i], rnorms[i]); 
      }
      fprintf(primme_svds.outputFile, "%d singular triplets converged\n", primme_svds.primme.initSize);

      fprintf(primme_svds.outputFile, "Tolerance : %-22.15E\n", 
                                                            primme_svds.aNorm*primme_svds.eps);
      #define PRINT_STATS(A, pre) { \
         fprintf(primme_svds.outputFile, pre "Iterations  : %-" PRIMME_INT_P "\n", (A).numOuterIterations); \
         fprintf(primme_svds.outputFile, pre "Restarts    : %-" PRIMME_INT_P "\n", (A).numRestarts);\
         fprintf(primme_svds.outputFile, pre "Matvecs     : %-" PRIMME_INT_P "\n", (A).numMatvecs);\
         fprintf(primme_svds.outputFile, pre "Preconds    : %-" PRIMME_INT_P "\n", (A).numPreconds);\
         fprintf(primme_svds.outputFile, pre "ElapsedTime : %-f\n", (A).elapsedTime);}

      if (primme_svds.methodStage2 != primme_svds_op_none) {
         PRINT_STATS(primme_svds.primme.stats, "1st ");
         PRINT_STATS(primme_svds.primmeStage2.stats, "2sd ");
      }
      PRINT_STATS(primme_svds.stats, "");
      if (primme_svds.locking && primme_svds.intWork && primme_svds.intWork[0] == 1) {
         fprintf(primme_svds.outputFile, "\nA locking problem has occurred.\n");
         fprintf(primme_svds.outputFile,
            "Some eigenpairs do not have a residual norm less than the tolerance.\n");
         fprintf(primme_svds.outputFile,
            "However, the subspace of svecs is accurate to the required tolerance.\n");
      }

      fprintf(primme_svds.outputFile, "\n\n#,%" PRIMME_INT_P ",%.1f\n\n", primme_svds.stats.numMatvecs,
         wt2-wt1); 
   }

   fclose(primme_svds.outputFile);
   destroyMatrixAndPrecond(&driver, &primme_svds, permutation);
   primme_svds_free(&primme_svds);
   free(svals);
   free(svecs);
   free(rnorms);

   if (ret != 0 && master) {
      fprintf(primme_svds.outputFile, 
         "Error: primme_svds returned with nonzero exit status: %d \n",ret);
      return -1;
   }

   if (retX != 0 && master) {
      fprintf(primme_svds.outputFile, 
         "Error: found some issues in the solution return by primme_svds\n");
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

#ifdef _OPENMP
#include <omp.h>
#endif

static int setMatrixAndPrecond(driver_params *driver,
      primme_svds_params *primme_svds, int **permutation) {
   int numProcs=1;
   double aNorm;

#  if defined(USE_MPI) || defined(USE_PETSC)
   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   primme_svds->commInfo = (MPI_Comm *)primme_calloc(1, sizeof(MPI_Comm), "MPI_Comm");
   *(MPI_Comm*)primme_svds->commInfo = MPI_COMM_WORLD;
#  endif

   if (driver->matrixChoice == driver_default) {
      if (numProcs <= 1) {
#        ifdef USE_RSB
         driver->matrixChoice = driver_rsb;
#        else
         driver->matrixChoice = driver_native;
#        endif
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
      *(MPI_Comm*)primme_svds->commInfo = MPI_COMM_WORLD;
#  endif
      {
         CSRMatrix *matrix;
         double *diag;
         /* Fix to use a single thread */
         #ifdef _OPENMP
         omp_set_num_threads(1);
         #endif
          
         if (readMatrixNative(driver->matrixFileName, &matrix, &aNorm) !=0 )
            return -1;
         primme_svds->matrix = matrix;
         primme_svds->matrixMatvec = CSRMatrixMatvecSVD;
         primme_svds->m = primme_svds->mLocal = matrix->m;
         primme_svds->n = primme_svds->nLocal = matrix->n;
         switch(driver->PrecChoice) {
         case driver_noprecond:
            primme_svds->preconditioner = NULL;
            primme_svds->applyPreconditioner = NULL;
            break;
         case driver_jacobi:
            createInvNormalPrecNative(matrix, driver->shift, &diag);
            primme_svds->preconditioner = diag;
            primme_svds->applyPreconditioner = ApplyInvNormalPrecNative;
            break;
         case driver_jacobi_i:
            createInvNormalPrecNative(matrix, 0, &diag);
            primme_svds->preconditioner = diag;
            primme_svds->applyPreconditioner = ApplyInvDavidsonNormalPrecNative;
            break;
         case driver_ilut:
            fprintf(stderr, "ERROR: ilut preconditioner is not supported with NATIVE, use other!\n");
            return -1;
         default:
            fprintf(stderr, "ERROR: preconditioner is not supported with NATIVE, use other!\n");
            return -1;
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
               driver->threshold, driver->filter, driver->isymm, MPI_COMM_WORLD, &aNorm,
               &primme_svds->n, &m, &primme_svds->nLocal, &mLocal, &primme_svds->numProcs, &primme_svds->procID, &matrix,
               (driver->PrecChoice == driver_ilut) ? &precond : NULL);
         *(MPI_Comm*)primme_svds->commInfo = MPI_COMM_WORLD;
         primme_svds->matrix = matrix;
         primme_svds->matrixMatvec = ParaSailsMatrixMatvec;
         primme_svds->preconditioner = precond;
         primme_svds->applyPreconditioner = precond ? ApplyPrecParaSails : NULL;
      }
#endif
      break;

   case driver_petsc:
#ifndef USE_PETSC
      fprintf(stderr, "ERROR: PETSc is needed!\n");
      return -1;
#else
      {
         Mat *matrix,A;
         PC *pc;
         double *diag;
         if (readMatrixPetsc(driver->matrixFileName, &primme_svds->m, &primme_svds->n, &primme_svds->mLocal,
               &primme_svds->nLocal, &primme_svds->numProcs, &primme_svds->procID, &matrix, &aNorm, permutation) != 0)
            return -1;
         *(MPI_Comm*)primme_svds->commInfo = PETSC_COMM_WORLD;
         primme_svds->matrix = matrix;
         primme_svds->matrixMatvec = PETScMatvecSVD;
         switch(driver->PrecChoice) {
         case driver_noprecond:
            primme_svds->preconditioner = NULL;
            primme_svds->applyPreconditioner = NULL;
            break;
         case driver_jacobi:
            createInvNormalPrecPETSC(*matrix, driver->shift, &diag);
            primme_svds->preconditioner = diag;
            primme_svds->applyPreconditioner = ApplyInvNormalPrecNative;
            break;
         case driver_jacobi_i:
            createInvNormalPrecPETSC(*matrix, 0, &diag);
            primme_svds->preconditioner = diag;
            primme_svds->applyPreconditioner = ApplyInvDavidsonNormalPrecNative;
            break;
         case driver_ilut:
            if (primme_svds->m == primme_svds->n) {
               PetscErrorCode ierr;
               pc = (PC *)primme_calloc(2, sizeof(PC), "pc");
               pc[1] = NULL;
               ierr = PCCreate(PETSC_COMM_WORLD, pc); CHKERRQ(ierr);
               ierr = MatDuplicate(*matrix, MAT_COPY_VALUES, &A);CHKERRQ(ierr);
               ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);
               ierr = MatShift(A, -driver->shift);CHKERRQ(ierr);
#ifdef PETSC_HAVE_HYPRE
               ierr = PCSetType(pc[0], PCHYPRE); CHKERRQ(ierr);
               ierr = PCHYPRESetType(pc[0], "boomeramg"); CHKERRQ(ierr);
               ierr = PCCreate(PETSC_COMM_WORLD, &pc[1]); CHKERRQ(ierr);
               ierr = PCSetType(pc[1], PCHYPRE); CHKERRQ(ierr);
               ierr = PCHYPRESetType(pc[1], "boomeramg"); CHKERRQ(ierr);
#else
               ierr = PCSetType(pc[0], PCBJACOBI); CHKERRQ(ierr);
#  ifdef USE_DOUBLECOMPLEX
               ierr = PCCreate(PETSC_COMM_WORLD, &pc[1]); CHKERRQ(ierr);
               ierr = PCSetType(pc[1], PCBJACOBI); CHKERRQ(ierr);
#  endif
#endif
               ierr = PCSetOperators(pc[0], A, A); CHKERRQ(ierr);
               ierr = PCSetFromOptions(pc[0]); CHKERRQ(ierr);
               ierr = PCSetUp(pc[0]); CHKERRQ(ierr);
               ierr = MatDestroy(&A);CHKERRQ(ierr);
               if (pc[1]) {
                  ierr = MatHermitianTranspose(*matrix,MAT_INITIAL_MATRIX,&A);CHKERRQ(ierr);
                  ierr = MatShift(A, -driver->shift);CHKERRQ(ierr);
                  ierr = PCSetOperators(pc[1], A, A); CHKERRQ(ierr);
                  ierr = PCSetFromOptions(pc[1]); CHKERRQ(ierr);
                  ierr = PCSetUp(pc[1]); CHKERRQ(ierr);
                  ierr = MatDestroy(&A);CHKERRQ(ierr);
               }
               primme_svds->preconditioner = pc;
               primme_svds->applyPreconditioner = ApplyPCPrecPETSCSVD;
            }
            else  {
               fprintf(stderr, "ERROR: ilut only supported for square matrices! Try bjacobi or normal\n");
               return -1;
            }
            break;
         case driver_normal:
            {
               PetscErrorCode ierr;
               Mat C;
               pc = (PC *)primme_calloc(1, sizeof(PC), "pc");
               ierr = PCCreate(PETSC_COMM_WORLD, pc); CHKERRQ(ierr);
#ifdef PETSC_HAVE_HYPRE
               ierr = PCSetType(*pc, PCHYPRE); CHKERRQ(ierr);
               ierr = PCHYPRESetType(*pc, "boomeramg"); CHKERRQ(ierr);
#else
               ierr = PCSetType(*pc, PCBJACOBI); CHKERRQ(ierr);
#endif
               ierr = MatHermitianTranspose(*matrix,MAT_INITIAL_MATRIX,&A);CHKERRQ(ierr);
               if (primme_svds->m >= primme_svds->n) {
                  ierr = MatMatMult(A,*matrix,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);
               }
               else {
                  ierr = MatMatMult(*matrix,A,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);
               }
               ierr = PCSetOperators(*pc, C, C); CHKERRQ(ierr);
               ierr = PCSetFromOptions(*pc); CHKERRQ(ierr);
               ierr = PCSetUp(*pc); CHKERRQ(ierr);
               primme_svds->preconditioner = pc;
               primme_svds->primme.applyPreconditioner = ApplyPCPrecPETSC;
               primme_svds->primme.preconditioner = pc;
               primme_svds->primme.correctionParams.precondition = 1;
               primme_svds->primme.commInfo = primme_svds->commInfo;
               ierr = MatDestroy(&A);CHKERRQ(ierr);
               ierr = MatDestroy(&C);CHKERRQ(ierr);
            }
            break;
         case driver_bjacobi:
            {
               PetscErrorCode ierr;
               Mat C, localA, localAt, At;
               IS iscols, isrows;
               PetscInt lowj, highj;
               MPI_Comm comm;
               pc = (PC *)primme_calloc(1, sizeof(PC), "pc");
               if (primme_svds->m > primme_svds->n) {
                  ierr = MatHermitianTranspose(*matrix,MAT_INITIAL_MATRIX,&A);CHKERRQ(ierr);
                  ierr = PetscObjectGetComm((PetscObject)*matrix,&comm);CHKERRQ(ierr);
                  ierr = MatGetOwnershipRangeColumn(*matrix,&lowj,&highj);CHKERRQ(ierr);  
                  ierr = ISCreateStride(comm, primme_svds->nLocal, lowj, 1, &isrows);CHKERRQ(ierr);
                  ierr = MatGetSubMatrix(A,isrows,NULL,MAT_INITIAL_MATRIX,&At);CHKERRQ(ierr);
                  ierr = MatDestroy(&A);CHKERRQ(ierr);
                  ierr = ISDestroy(&isrows);CHKERRQ(ierr);
                  ierr = MatMPIAIJGetLocalMat(At, MAT_INITIAL_MATRIX, &localA);
                  ierr = MatDestroy(&At);CHKERRQ(ierr);
               }
               else {
                  ierr = MatMPIAIJGetLocalMat(*matrix, MAT_INITIAL_MATRIX, &localA);
               }
               ierr = MatHermitianTranspose(localA, MAT_INITIAL_MATRIX, &localAt); CHKERRQ(ierr);
               ierr = MatMatMult(localA,localAt,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&C);CHKERRQ(ierr);
               ierr = MatDestroy(&localA);CHKERRQ(ierr);
               ierr = MatDestroy(&localAt);CHKERRQ(ierr);
               ierr = PetscObjectGetComm((PetscObject)C, &comm);CHKERRQ(ierr);
               ierr = PCCreate(comm, pc); CHKERRQ(ierr);
#ifdef PETSC_HAVE_HYPRE
               ierr = PCSetType(*pc, PCHYPRE); CHKERRQ(ierr);
               ierr = PCHYPRESetType(*pc, "boomeramg"); CHKERRQ(ierr);
#else
               ierr = PCSetType(*pc, PCBJACOBI); CHKERRQ(ierr);
#endif
               ierr = PCSetOperators(*pc, C, C); CHKERRQ(ierr);
               ierr = PCSetFromOptions(*pc); CHKERRQ(ierr);
               ierr = PCSetUp(*pc); CHKERRQ(ierr);
               primme_svds->preconditioner = pc;
               primme_svds->primme.applyPreconditioner = ApplyPCPrecPETSC;
               primme_svds->primme.preconditioner = pc;
               primme_svds->primme.correctionParams.precondition = 1;
               primme_svds->primme.commInfo = primme_svds->commInfo;
               ierr = MatDestroy(&C);CHKERRQ(ierr);
            }
         }
      }
#endif
      break;

   case driver_rsb:
#if !defined(USE_RSB)
      fprintf(stderr, "ERROR: RSB is needed!\n");
      return -1;
#else
#  if defined(USE_MPI)
      if (numProcs != 1) {
         fprintf(stderr, "ERROR: MPI is not supported with RSB, use other!\n");
         return -1;
      }
      *(MPI_Comm*)primme->commInfo = MPI_COMM_WORLD;
#  endif
      {
         blas_sparse_matrix *matrix;
         double *diag;
         matrix = (blas_sparse_matrix *)primme_calloc(1, sizeof(blas_sparse_matrix), "matrix");
         
         if (readMatrixRSB(driver->matrixFileName, matrix, &aNorm) !=0 )
            return -1;
         primme_svds->matrix = matrix;
         primme_svds->matrixMatvec = RSBMatvecSVD;
         primme_svds->m = primme_svds->mLocal = BLAS_usgp(*matrix, blas_num_rows);
         primme_svds->n = primme_svds->nLocal = BLAS_usgp(*matrix, blas_num_cols);
         switch(driver->PrecChoice) {
         case driver_noprecond:
            primme_svds->preconditioner = NULL;
            primme_svds->applyPreconditioner = NULL;
            break;
         case driver_jacobi:
            createInvNormalPrecRSB(*matrix, driver->shift, &diag);
            primme_svds->preconditioner = diag;
            primme_svds->applyPreconditioner = ApplyInvNormalPrecNative;
            break;
         case driver_jacobi_i:
            createInvNormalPrecRSB(*matrix, 0.0, &diag);
            primme_svds->preconditioner = diag;
            primme_svds->applyPreconditioner = ApplyInvDavidsonNormalPrecNative;
            break;
         default:
            assert(0);
            break;
         }
      }
#endif
      break;

   }

   if (primme_svds->aNorm < 0) primme_svds->aNorm = aNorm;

#if defined(USE_MPI)
   primme_svds->globalSumReal = par_GlobalSumDoubleSvds;
#endif
   return 0;
}

static int destroyMatrixAndPrecond(driver_params *driver, primme_svds_params *primme_svds, int *permutation) {
   switch(driver->matrixChoice) {
   case driver_default:
      assert(0);
      break;
   case driver_native:
#if !defined(USE_NATIVE)
      fprintf(stderr, "ERROR: NATIVE is needed!\n");
      return -1;
#else
      freeCSRMatrix((CSRMatrix*)primme_svds->matrix);

      switch(driver->PrecChoice) {
      case driver_noprecond:
         break;
      case driver_jacobi:
      case driver_jacobi_i:
         free(primme_svds->preconditioner);
         break;
      case driver_ilut:
         break;
      default:
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
         PC *pc = (PC*)primme_svds->preconditioner;
         ierr = MatDestroy((Mat*)primme_svds->matrix);CHKERRQ(ierr);
         switch(driver->PrecChoice) {
         case driver_noprecond:
         break;
         case driver_ilut:
            if (pc[1]) {
               ierr = PCDestroy(&pc[1]);CHKERRQ(ierr);
            }
         case driver_normal:
            ierr = PCDestroy(&pc[0]);CHKERRQ(ierr);
            free(primme_svds->preconditioner);
            break;
         case driver_jacobi:
         case driver_jacobi_i:
            free(primme_svds->preconditioner);
            break;
         }
      }
#endif
      break;

   case driver_rsb:
#if !defined(USE_RSB)
      fprintf(stderr, "ERROR: RSB is needed!\n");
      return -1;
#else
      {
         blas_sparse_matrix *matrix = primme_svds->matrix;
         BLAS_usds(*matrix);
         free(matrix);

         switch(driver->PrecChoice) {
         case driver_noprecond:
            break;
         case driver_jacobi:
         case driver_jacobi_i:
            free(primme_svds->preconditioner);
            break;
         default:
            assert(0);
            break;
         }
      }
#endif
      break;

   }
#if defined(USE_MPI)
   free(primme_svds->commInfo);
#endif
   if (permutation) free(permutation);
   return 0;
}
