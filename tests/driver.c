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
PetscLogEvent PRIMME_GLOBAL_SUM;
#endif
#ifdef USE_RSB
#  include "rsbw.h"
#endif

/* primme.h header file is required to run primme */
#include "primme.h"
#include "shared_utils.h"
#include "ioandtest.h"
/* wtime.h header file is included so primme's timimg functions can be used */
#include "../../src/include/wtime.h"

static int real_main (int argc, char *argv[]);
static int setMatrixAndPrecond(driver_params *driver, primme_params *primme, int **permutation);
static int destroyMatrixAndPrecond(driver_params *driver, primme_params *primme, int *permutation);



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
   double *evals, *rnorms;
   SCALAR *evecs;
   driver_params driver;
   primme_params primme;
   primme_preset_method method=PRIMME_DEFAULT_METHOD;
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
                           &primme, "primme.", &method, "method") < 0) {
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
   primme_set_method(method, &primme);

   /* --------------------------------------- */
   /* Optional: report memory requirements    */
   /* --------------------------------------- */

   ret = Sprimme(NULL,NULL,NULL,&primme);
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
      driver_display_method(method, "method", primme.outputFile);
   }

   /* --------------------------------------------------------------------- */
   /*                            Run the d/zprimme solver                   */
   /* --------------------------------------------------------------------- */

   /* Allocate space for converged Ritz values and residual norms */

   evals = (double *)primme_calloc(primme.numEvals, sizeof(double), "evals");
   evecs = (SCALAR *)primme_calloc(primme.nLocal*primme.numEvals, 
                                sizeof(SCALAR), "evecs");
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
         SCALAR *r = (SCALAR *)primme_calloc(primme.nLocal,sizeof(SCALAR), "random");
         double norm;
         int j;
         for (i=primme.numOrthoConst; i<min(cols, primme.initSize+primme.numOrthoConst); i++) {
            Num_larnv_Sprimme(2, primme.iseed, primme.nLocal, r);
            norm = sqrt(REAL_PART(Num_dot_Sprimme(primme.nLocal, r, 1, r, 1)));
            for (j=0; j<primme.nLocal; j++)
               evecs[primme.nLocal*i+j] += r[j]/norm*driver.initialGuessesPert;
         }
         free(r);
      }
      Num_larnv_Sprimme(2, primme.iseed, (primme.initSize+primme.numOrthoConst-i)*primme.nLocal,
                     &evecs[primme.nLocal*i]);
   } else if (primme.numOrthoConst > 0) {
      ASSERT_MSG(0, 1, "numOrthoConst > 0 but no value in initialGuessesFileName.\n");
   } else if (primme.initSize > 0) {
      Num_larnv_Sprimme(2, primme.iseed, primme.initSize*primme.nLocal, evecs);
   }


   /* ------------- */
   /*  Call primme  */
   /* ------------- */

   wt1 = primme_get_wtime(); 
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
   primme_get_time(&ut1,&st1);
#endif

   ret = Sprimme(evals, evecs, rnorms, &primme);

   wt2 = primme_get_wtime();
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
   primme_get_time(&ut2,&st2);
#endif

   if (driver.checkXFileName[0]) {
      retX = check_solution(driver.checkXFileName, &primme, evals, evecs, rnorms, permutation, driver.checkInterface);
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
      fprintf(primme.outputFile, "Wallclock Runtime   : %-f\n", wt2-wt1);
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
      fprintf(primme.outputFile, "User Time           : %f seconds\n", ut2-ut1);
      fprintf(primme.outputFile, "Syst Time           : %f seconds\n", st2-st1);
#endif

      for (i=0; i < primme.initSize; i++) {
         fprintf(primme.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
            evals[i], rnorms[i]); 
      }
      fprintf(primme.outputFile, " %d eigenpairs converged\n", primme.initSize);

      fprintf(primme.outputFile, "Tolerance  : %-22.15E\n", 
                                                            primme.aNorm*primme.eps);
      fprintf(primme.outputFile, "Iterations : %-" PRIMME_INT_P "\n", 
                                                    primme.stats.numOuterIterations); 
      fprintf(primme.outputFile, "Restarts   : %-" PRIMME_INT_P "\n", primme.stats.numRestarts);
      fprintf(primme.outputFile, "Matvecs    : %-" PRIMME_INT_P "\n", primme.stats.numMatvecs);
      fprintf(primme.outputFile, "Preconds   : %-" PRIMME_INT_P "\n", primme.stats.numPreconds);
      fprintf(primme.outputFile, "Time matvecs  : %f\n",  primme.stats.timeMatvec);
      fprintf(primme.outputFile, "Time precond  : %f\n",  primme.stats.timePrecond);
      fprintf(primme.outputFile, "Time ortho  : %f\n",  primme.stats.timeOrtho);
      if (primme.locking && primme.intWork && primme.intWork[0] == 1) {
         fprintf(primme.outputFile, "\nA locking problem has occurred.\n");
         fprintf(primme.outputFile,
            "Some eigenpairs do not have a residual norm less than the tolerance.\n");
         fprintf(primme.outputFile,
            "However, the subspace of evecs is accurate to the required tolerance.\n");
      }

      fprintf(primme.outputFile, "\n\n#,%" PRIMME_INT_P ",%.1f\n\n", primme.stats.numMatvecs,
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
   primme_free(&primme);
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

#ifdef _OPENMP
#include <omp.h>
#endif

static int setMatrixAndPrecond(driver_params *driver, primme_params *primme, int **permutation) {
   int numProcs=1;
   double aNorm;

#  if defined(USE_MPI) || defined(USE_PETSC)
   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   primme->commInfo = (MPI_Comm *)primme_calloc(1, sizeof(MPI_Comm), "MPI_Comm");
   *(MPI_Comm*)primme->commInfo = MPI_COMM_WORLD;
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
      *(MPI_Comm*)primme->commInfo = MPI_COMM_WORLD;
#  endif
      {
         CSRMatrix *matrix, *prec;
         double *diag;
         /* Fix to use a single thread */
         #ifdef _OPENMP
         omp_set_num_threads(1);
         #endif
          
         if (readMatrixNative(driver->matrixFileName, &matrix, &aNorm) !=0 )
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
            createInvDiagPrecNative(matrix, 0.0, &diag);
            primme->preconditioner = diag;
            primme->applyPreconditioner = ApplyInvDavidsonDiagPrecNative;
            break;
         case driver_ilut:
            createILUTPrecNative(matrix, driver->shift, driver->level, driver->threshold,
                                 driver->filter, &prec);
            primme->preconditioner = prec;
            primme->applyPreconditioner = ApplyILUTPrecNative;
            break;
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
         Mat *matrix,A;
         PC *pc;
         Vec *vec;
         PRIMME_INT m, mLocal;
         if (readMatrixPetsc(driver->matrixFileName, &primme->n, &m, &primme->nLocal, &mLocal,
                         &primme->numProcs, &primme->procID, &matrix, &aNorm, permutation) != 0)
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
                  ierr = PCSetType(*pc, PCILU); CHKERRQ(ierr);
               }
               else {
                  #ifdef PETSC_HAVE_HYPRE
                     ierr = PCSetType(*pc, PCHYPRE); CHKERRQ(ierr);
                     ierr = PCHYPRESetType(*pc, "boomeramg"); CHKERRQ(ierr);
                  #else
                     ierr = PCSetType(*pc, PCBJACOBI); CHKERRQ(ierr);
                  #endif
               }
            }
            ierr = MatDuplicate(*matrix, MAT_COPY_VALUES, &A);CHKERRQ(ierr);
            ierr = MatShift(A, -driver->shift);CHKERRQ(ierr);
            ierr = PCSetOperators(*pc, A, A); CHKERRQ(ierr);
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
         primme->matrix = matrix;
         primme->matrixMatvec = RSBMatvec;
         primme->n = primme->nLocal = BLAS_usgp(*matrix, blas_num_rows);
         switch(driver->PrecChoice) {
         case driver_noprecond:
            primme->preconditioner = NULL;
            primme->applyPreconditioner = NULL;
            break;
         case driver_jacobi:
            createInvDiagPrecRSB(*matrix, driver->shift, &diag);
            primme->preconditioner = diag;
            primme->applyPreconditioner = ApplyInvDiagPrecNative;
            break;
         case driver_jacobi_i:
            createInvDiagPrecRSB(*matrix, 0.0, &diag);
            primme->preconditioner = diag;
            primme->applyPreconditioner = ApplyInvDavidsonDiagPrecNative;
            break;
          default:
            assert(0);
            break;
         }
      }
#endif
      break;

   }

   if (primme->aNorm < 0) primme->aNorm = aNorm;

#if defined(USE_MPI)
   primme->globalSumReal = par_GlobalSumDouble;
#endif

#ifdef NOT_USE_ALIGNMENT
   primme->ldOPs = primme->nLocal ? primme->nLocal : primme->n;
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
      freeCSRMatrix((CSRMatrix*)primme->matrix);

      switch(driver->PrecChoice) {
      case driver_noprecond:
         break;
      case driver_jacobi:
      case driver_jacobi_i:
         free(primme->preconditioner);
         break;
      case driver_ilut:
         if (primme->preconditioner) {
            freeCSRMatrix((CSRMatrix*)primme->preconditioner);
         }
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

   case driver_rsb:
#if !defined(USE_RSB)
      fprintf(stderr, "ERROR: RSB is needed!\n");
      return -1;
#else
      {
         blas_sparse_matrix *matrix = primme->matrix;
         BLAS_usds(*matrix);
         free(matrix);

         switch(driver->PrecChoice) {
         case driver_noprecond:
            break;
         case driver_jacobi:
         case driver_jacobi_i:
            free(primme->preconditioner);
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
   free(primme->commInfo);
#endif
   if (permutation) free(permutation);
   return 0;
}
