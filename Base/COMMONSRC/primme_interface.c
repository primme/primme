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
 **********************************************************************
 * File: primme_interface.c
 *
 * Purpose - Contains interface functions to PRIMME named primme_* 
 *           If desired, the user can call any of these functions for 
 *           initializing parameters, setting up the method, 
 *           allocating memory, or even checking a given set of parameters. 
 *
 ******************************************************************************/

#if !(defined (__APPLE__) && defined (__MACH__))
#include <malloc.h>
#endif
#include <stdlib.h>   /* mallocs, free */
#include <stdio.h>    
#include <math.h>    
#include "primme.h"
#include "common_numerical.h"
#include "primme_interface.h"
#include "const.h"

/***************************************************************************

   Initialize the primme data structure
  
***************************************************************************/
void primme_initialize(primme_params *primme) {

   /* Essential parameters */
   primme->n                       = 0;
   primme->numEvals                = 1;
   primme->target                  = primme_smallest;
   primme->aNorm                   = 0.0L;
   primme->eps                     = 1e-12;

   /* Matvec and preconditioner */
   primme->matrixMatvec            = NULL;
   primme->applyPreconditioner     = NULL;
   primme->massMatrixMatvec        = NULL;

   /* Shifts for interior eigenvalues*/
   primme->numTargetShifts         = 0;
   primme->targetShifts            = NULL;

   /* Parallel computing parameters */
   primme->numProcs                = 1;
   primme->procID                  = 0;
   primme->nLocal                  = 0;
   primme->commInfo                = NULL;
   primme->globalSumDouble         = NULL;

   /* Initial guesses/constraints */
   primme->initSize                = 0;
   primme->numOrthoConst           = 0;

   primme->projectionParams.projection = primme_proj_default;
   primme->projectionParams.refinedScheme  = primme_ref_default;

   primme->currentEstimates.Anormest           = 0;
   primme->currentEstimates.targetRitzVal      = 0;
   primme->currentEstimates.targetRitzValNorm  = 0;
   primme->currentEstimates.targetRitzVec      = NULL;
   primme->InitBasisMode                       = primme_init_default;
   primme->ReIntroInitGuessToBasis             = 0;

   /* Eigensolver parameters (outer) */
   primme->locking                             = -1;
   primme->dynamicMethodSwitch                 = 0;
   primme->maxBasisSize                        = 0;
   primme->minRestartSize                      = 0;
   primme->maxBlockSize                        = 0;
   primme->maxMatvecs                          = INT_MAX;
   primme->maxOuterIterations                  = INT_MAX;
   primme->restartingParams.scheme             = primme_thick;
   primme->restartingParams.maxPrevRetain      = -1;

   /* correction parameters (inner) */
   primme->correctionParams.precondition       = -1;
   primme->correctionParams.robustShifts       = 0;
   primme->correctionParams.maxInnerIterations = -INT_MAX;
   primme->correctionParams.projectors.LeftQ   = 0;
   primme->correctionParams.projectors.LeftX   = 0;
   primme->correctionParams.projectors.RightQ  = 0;
   primme->correctionParams.projectors.RightX  = 0;
   primme->correctionParams.projectors.SkewQ   = 0;
   primme->correctionParams.projectors.SkewX   = 0;
   primme->correctionParams.relTolBase         = 0;
   primme->correctionParams.convTest           = primme_adaptive_ETolerance;

   /* Printing and reporting */
   primme->outputFile              = stdout;
   primme->printLevel              = 1;
   primme->stats.numOuterIterations= 0;
   primme->stats.numRestarts       = 0;
   primme->stats.numMatvecs        = 0;
   primme->stats.numPreconds       = 0;
   primme->stats.elapsedTime       = 0.0L;
   primme->stats.estimateMaxEVal   = -HUGE_VAL;
   primme->stats.estimateMinEVal   = HUGE_VAL;
   primme->stats.estimateLargestSVal = -HUGE_VAL;
   primme->stats.maxConvTol        = 0.0L;

   /* Optional user defined structures */
   primme->matrix                  = NULL;
   primme->preconditioner          = NULL;

   /* Internally used variables */
   primme->iseed[0] = -1;   /* To set iseed, we first need procID           */
   primme->iseed[1] = -1;   /* Thus we set all iseeds to -1                 */
   primme->iseed[2] = -1;   /* Unless users provide their own iseeds,       */
   primme->iseed[3] = -1;   /* PRIMME will set thse later uniquely per proc */
   primme->intWorkSize             = 0;
   primme->realWorkSize            = 0;
   primme->intWork                 = NULL;
   primme->realWork                = NULL;
   primme->stackTrace              = NULL;
   primme->ShiftsForPreconditioner = NULL;
   primme->convTestFun             = NULL;

}

/***************************************************************************
   wrappers for allocating and freeing space in a more friendly way
***************************************************************************/

/* Memory alignment at page boundary */
void *primme_valloc(size_t byteSize, const char *target) {

   void *ptr;

   if ( (ptr = malloc(byteSize)) == NULL) {
      perror("primme_alloc");
      fprintf(stderr,
         "ERROR(primme_alloc): Could not allocate %lu bytes for: %s\n",
         byteSize, target);
      fflush(stderr);
      exit(EXIT_FAILURE);
   }

   return (ptr);

}

/* Memory allocation with zeroing */
void *primme_calloc(size_t nelem, size_t elsize, const char *target) {

   void *ptr;

   if ((ptr = calloc(nelem, elsize)) == NULL) {
      perror("primme_calloc");
      fprintf(stderr, 
         "ERROR(primme_calloc): Could not allocate %lu elements of %lu bytes for: %s\n",
         nelem, elsize, target);
      fflush(stderr);
      exit(EXIT_FAILURE);
   }

   return (ptr);

}

void primme_Free(primme_params *params) {

   free(params->intWork);
   free(params->realWork);
   params->intWorkSize  = 0;
   params->realWorkSize = 0;

} /**************************************************************************/

/******************************************************************************
 * int primme_set_method(primme_preset_method method,primme_params *params)
 *
 *    Set the eigensolver parameters to implement a method requested by the user
 *    A choice of 15 preset methods is provided. These implement 12 well 
 *    known methods. The two default methods are shells for easy access to
 *    the best performing GD+k and JDQMR_ETol with expertly chosen parameters.
 *    
 *    The DYNAMIC method makes runtime measurements and switches dynamically
 *    between DEFAULT_MIN_TIME and DEFAULT_MIN_MATVEC keeping the one that
 *    performs the best. Because it usually achieves best runtime over all 
 *    methods, it is recommended for the general user.
 *
 *    For most methods the user may specify the maxBasisSize, restart size, 
 *    block size, etc. If any (or all) of these parameters are not specified, 
 *    they will be given default values that are appropriate for the method.
 *
 *    primme_set_method() will override those parameters in primme that 
 *    are needed to implement the method.
 *
 *    Note: Spectral transformations can be applied by simply providing
 *         (A-shift I)^{-1} as the matvec.
 *    
 *    For additional information see the readme file
 *
 * INPUT
 * -----
 *    method    One of the following 12 enum methods:
 *
 *        DYNAMIC,                 : Switches dynamically to the best method
 *        DEFAULT_MIN_TIME,        : Currently set at JDQMR_ETol
 *        DEFAULT_MIN_MATVECS,     : Currently set at GD+block
 *        Arnoldi,                 : obviously not an efficient choice 
 *        GD,                      : classical block Generalized Davidson 
 *        GD_plusK,                : GD+k block GD with recurrence restarting
 *        GD_Olsen_plusK,          : GD+k with approximate Olsen precond.
 *        JD_Olsen_plusK,          : GD+k, exact Olsen (two precond per step)
 *        RQI,                     : Rayleigh Quotient Iteration. Also INVIT,
 *                                 :   but for INVIT provide targetShifts
 *        JDQR,                    : Original block, Jacobi Davidson
 *        JDQMR,                   : Our block JDQMR method (similar to JDCG)
 *        JDQMR_ETol,              : Slight, but efficient JDQMR modification
 *        SUBSPACE_ITERATION,      : equiv. to GD(block,2*block)
 *        LOBPCG_OrthoBasis,       : equiv. to GD(nev,3*nev)+nev
 *        LOBPCG_OrthoBasis_Window : equiv. to GD(block,3*block)+block nev>block
 *
 *
 * INPUT/OUTPUT
 * ------
 *    params    The main structure to be used by Primme, with parameters
 *              reflecting the user's choice of method 
 *
 *
 * return value         Note: the return value is a LONG INT
 *
 *      = 0 parameters for "method" have been set.
 *
 *      < 0 no such method exists. If not defined by the user, defaults have
 *          been set for maxBasisSize, minRestartSize, and maxBlockSize.
 *
 ******************************************************************************/
int primme_set_method(primme_preset_method method, primme_params *params) {

   /* From our experience, these two methods yield the smallest matvecs/time */
   /* DYNAMIC will make some timings before it settles on one of the two     */
   if (method == DEFAULT_MIN_MATVECS) {
      method = GD_Olsen_plusK;
   }
   else if (method == DEFAULT_MIN_TIME) {
      /* JDQMR works better than JDQMR_ETol in interior problems. */
      if (params->target == primme_smallest || params->target == primme_largest) {
         method = JDQMR_ETol;
      }
      else {
         method = JDQMR;
      }
   }
   else if (method == DYNAMIC) {
      /* JDQMR works better than JDQMR_ETol in interior problems. */
      if (params->target == primme_smallest || params->target == primme_largest) {
         method = JDQMR_ETol;
      }
      else {
         method = JDQMR;
      }
      params->dynamicMethodSwitch = 1;
   }
  
   if (params->dynamicMethodSwitch < 0) {
      params->dynamicMethodSwitch = 0;
   }

   if (params->maxBlockSize == 0) {
      params->maxBlockSize = 1;
   }
   if (params->correctionParams.precondition == -1) {
      params->correctionParams.precondition = params->applyPreconditioner ? 1 : 0;
   }

   if (method == Arnoldi) {
      params->restartingParams.maxPrevRetain      = 0;
      params->correctionParams.precondition       = 0;
      params->correctionParams.maxInnerIterations = 0;
   }
   else if (method == GD) {
      params->restartingParams.maxPrevRetain      = 0;
      params->correctionParams.robustShifts       = 1;
      params->correctionParams.maxInnerIterations = 0;
      params->correctionParams.projectors.RightX  = 0;
      params->correctionParams.projectors.SkewX   = 0;
   }
   else if (method == GD_plusK) {
      if (params->restartingParams.maxPrevRetain <= 0) {
         if (params->maxBlockSize == 1 && params->numEvals > 1) {
            params->restartingParams.maxPrevRetain = 2;
         }
         else {
            params->restartingParams.maxPrevRetain = params->maxBlockSize;
         }
      }
      params->correctionParams.maxInnerIterations = 0;
      params->correctionParams.projectors.RightX  = 0;
      params->correctionParams.projectors.SkewX   = 0;
   }
   else if (method == GD_Olsen_plusK) {
      if (params->restartingParams.maxPrevRetain <= 0) {
         if (params->maxBlockSize == 1 && params->numEvals > 1) {
            params->restartingParams.maxPrevRetain = 2;
         }
         else {
            params->restartingParams.maxPrevRetain = params->maxBlockSize;
         }
      }
      params->correctionParams.maxInnerIterations = 0;
      params->correctionParams.projectors.RightX  = 1;
      params->correctionParams.projectors.SkewX   = 0;
   }
   else if (method == JD_Olsen_plusK) {
      if (params->restartingParams.maxPrevRetain <= 0) {
         if (params->maxBlockSize == 1 && params->numEvals > 1) {
            params->restartingParams.maxPrevRetain = 2;
         }
         else {
            params->restartingParams.maxPrevRetain = params->maxBlockSize;
         }
      }
      params->correctionParams.robustShifts       = 1;
      params->correctionParams.maxInnerIterations = 0;
      params->correctionParams.projectors.RightX  = 1;
      params->correctionParams.projectors.SkewX   = 1;
   }
   else if (method == RQI) {
     /* This is accelerated RQI as basisSize may be > 2                 */
     /* NOTE: If numTargetShifts > 0 and targetShifts[*] are provided,  *
      * the interior problem solved uses these shifts in the correction *
      * equation. Therefore RQI becomes INVIT in that case.             */
      params->locking                             = 1;
      params->restartingParams.maxPrevRetain      = 0;
      params->correctionParams.robustShifts       = 1;
      params->correctionParams.maxInnerIterations = -1;
      params->correctionParams.projectors.LeftQ   = 1;
      params->correctionParams.projectors.LeftX   = 1;
      params->correctionParams.projectors.RightQ  = 0;
      params->correctionParams.projectors.RightX  = 1;
      params->correctionParams.projectors.SkewQ   = 0;
      params->correctionParams.projectors.SkewX   = 0;
      params->correctionParams.convTest           = primme_full_LTolerance;
   }
   else if (method == JDQR) {
      params->locking                             = 1;
      params->restartingParams.maxPrevRetain      = 1;
      params->correctionParams.robustShifts       = 0;
      if (params->correctionParams.maxInnerIterations == -INT_MAX) {
         params->correctionParams.maxInnerIterations = 10;
      }
      params->correctionParams.projectors.LeftQ   = 0;
      params->correctionParams.projectors.LeftX   = 1;
      params->correctionParams.projectors.RightQ  = 1;
      params->correctionParams.projectors.RightX  = 1;
      params->correctionParams.projectors.SkewQ   = 1;
      params->correctionParams.projectors.SkewX   = 1;
      params->correctionParams.relTolBase         = 1.5;
      params->correctionParams.convTest           = primme_full_LTolerance;
   }
   else if (method == JDQMR) {
      if (params->restartingParams.maxPrevRetain < 0) {
         params->restartingParams.maxPrevRetain   = 1;
      }
      params->correctionParams.maxInnerIterations = -1;
      if (params->correctionParams.precondition) {
         params->correctionParams.projectors.LeftQ   = 1;
      }
      else {
         params->correctionParams.projectors.LeftQ   = 0;
      }
      params->correctionParams.projectors.LeftX   = 1;
      params->correctionParams.projectors.RightQ  = 0;
      params->correctionParams.projectors.RightX  = 0;
      params->correctionParams.projectors.SkewQ   = 0;
      params->correctionParams.projectors.SkewX   = 1;
      params->correctionParams.convTest           = primme_adaptive;
   }
   else if (method == JDQMR_ETol) {
      if (params->restartingParams.maxPrevRetain < 0) {
         params->restartingParams.maxPrevRetain   = 1;
      }
      params->correctionParams.maxInnerIterations = -1;
      if (params->correctionParams.precondition) {
         params->correctionParams.projectors.LeftQ   = 1;
      }
      else {
         params->correctionParams.projectors.LeftQ   = 0;
      }
      params->correctionParams.projectors.LeftX   = 1;
      params->correctionParams.projectors.RightQ  = 0;
      params->correctionParams.projectors.RightX  = 0;
      params->correctionParams.projectors.SkewQ   = 0;
      params->correctionParams.projectors.SkewX   = 1;
      params->correctionParams.convTest           = primme_adaptive_ETolerance;
   }
   else if (method == SUBSPACE_ITERATION) {
      params->locking                             = 1;
      params->maxBasisSize                        = params->numEvals*2;
      params->minRestartSize                      = params->numEvals;
      params->maxBlockSize                        = params->numEvals;
      params->restartingParams.scheme             = primme_thick;
      params->restartingParams.maxPrevRetain      = 0;
      params->correctionParams.robustShifts       = 0;
      params->correctionParams.maxInnerIterations = 0;
      params->correctionParams.projectors.RightX  = 1;
      params->correctionParams.projectors.SkewX   = 0;
   }
   else if (method == LOBPCG_OrthoBasis) {
      params->maxBasisSize                        = params->numEvals*3;
      params->minRestartSize                      = params->numEvals;
      params->maxBlockSize                        = params->numEvals;
      params->restartingParams.maxPrevRetain      = params->numEvals;
      params->restartingParams.scheme             = primme_thick;
      params->correctionParams.robustShifts       = 0;
      params->correctionParams.maxInnerIterations = 0;
      params->correctionParams.projectors.RightX  = 1;
      params->correctionParams.projectors.SkewX   = 0;
   }
   else if (method == LOBPCG_OrthoBasis_Window) {
      params->maxBasisSize                        = params->maxBlockSize*3;
      params->minRestartSize                      = params->maxBlockSize;
      params->restartingParams.maxPrevRetain      = params->maxBlockSize;
      params->restartingParams.scheme             = primme_thick;
      params->correctionParams.robustShifts       = 0;
      params->correctionParams.maxInnerIterations = 0;
      params->correctionParams.projectors.RightX  = 1;
      params->correctionParams.projectors.SkewX   = 0;
   }
   else {
      return -1;
   }

   primme_set_defaults(params);

   return 0;

  /**************************************************************************/
} /* end of primme_set_method */
  /**************************************************************************/

/******************************************************************************
 * Subroutine primme_set_defaults - Set valid values for options that still
 *    has the initial invalid value set in primme_initialize.
 *
 ******************************************************************************/

void primme_set_defaults(primme_params *params) {
   if (params->dynamicMethodSwitch < 0) {
      primme_set_method(DYNAMIC, params);
   }

   /* ----------------------------------------- */
   /* Set some defaults for sequential programs */
   /* ----------------------------------------- */
   if (params->numProcs <= 1) {
      params->nLocal = params->n;
      params->procID = 0;
      if (params->globalSumDouble == NULL) 
         params->globalSumDouble = primme_seq_globalSumDouble;
   }

   /* --------------------------------------------------------------------- */
   /* Decide on whether to use locking (hard locking), or not (soft locking)*/
   /* --------------------------------------------------------------------- */
   if (params->target != primme_smallest && params->target != primme_largest ) {
      /* Locking is necessary as interior Ritz values can cross shifts */
      params->locking = 1;
   }
   else if (params->numEvals > params->minRestartSize) {
      /* use locking when not enough vectors to restart with */
      params->locking = 1;
   }
   else if (params->locking == -1) {
      params->locking = 0;   
   }

   if (params->projectionParams.projection == primme_proj_default)
      params->projectionParams.projection = primme_proj_RR;
   if (params->projectionParams.refinedScheme == primme_ref_default)
      params->projectionParams.refinedScheme = primme_ref_none;
   if (params->InitBasisMode == primme_init_default)
      params->InitBasisMode = primme_init_krylov;
      
   /* Now that most of the parameters have been set, set defaults  */
   /* for basisSize, restartSize (for those methods that need it)  */
   /* For interior, larger basisSize and restartSize are advisable */
   /* If maxBlockSize is provided, assign at least 4*blocksize     */
   /* and consider also minRestartSize and maxPrevRetain           */
   if (params->maxBasisSize == 0) {
      if (params->target==primme_smallest || params->target==primme_largest)
         params->maxBasisSize   = min(params->n, max(
            max(15, 4*params->maxBlockSize+params->restartingParams.maxPrevRetain), 
            (int) 2.5*params->minRestartSize+params->restartingParams.maxPrevRetain));
      else
         params->maxBasisSize   = min(params->n, max(
            max(35, 5*params->maxBlockSize+params->restartingParams.maxPrevRetain),
            (int) 1.7*params->minRestartSize+params->restartingParams.maxPrevRetain));
   }

   if (params->minRestartSize == 0) {
      if (params->target==primme_smallest || params->target==primme_largest)
         params->minRestartSize = (int) (0.5 + 0.4*params->maxBasisSize);
      else
         params->minRestartSize = (int) (0.5 + 0.6*params->maxBasisSize);

      /* Adjust so that an integer number of blocks are added between restarts*/
      /* restart=basis-block*ceil((basis-restart-prevRetain)/block)-prevRetain*/
      if (params->maxBlockSize > 1) {
         if (params->restartingParams.maxPrevRetain > 0) 
            params->minRestartSize = params->maxBasisSize-params->maxBlockSize*
            (1 + (int) ((params->maxBasisSize - params->minRestartSize - 1 
                         -params->restartingParams.maxPrevRetain ) / (double) 
            params->maxBlockSize) ) - params->restartingParams.maxPrevRetain ;
         else 
            params->minRestartSize = params->maxBasisSize-params->maxBlockSize*
            (1 + (int) ((params->maxBasisSize - params->minRestartSize - 1)
                        /(double) params->maxBlockSize) );
      }
   }
}

/******************************************************************************
 *
 * void primme_display_params(primme_params *primme);
 *
 *    Displays the current configuration of primme data structure
 *
 *****************************************************************************/
void primme_display_params(primme_params primme) {

   fprintf(primme.outputFile,
           "// ---------------------------------------------------\n"
           "//                 primme configuration               \n"
           "// ---------------------------------------------------\n");

   primme_display_params_prefix("primme", primme);
   fflush(primme.outputFile);
}

void primme_display_params_prefix(const char* prefix, primme_params primme) {

   int i;
   FILE *outputFile = primme.outputFile;

#define PRINT(P,L) fprintf(outputFile, "%s." #P " = " #L "\n", prefix, primme. P);
#define PRINTIF(P,V) if (primme. P == V) fprintf(outputFile, "%s." #P " = " #V "\n", prefix);
#define PRINTParams(P,S,L) fprintf(outputFile, "%s." #P "." #S " = " #L "\n", \
                                    prefix, primme. P ## Params. S);
#define PRINTParamsIF(P,S,V) if (primme. P ## Params. S == V) \
                                 fprintf(outputFile, "%s." #P "." #S " = " #V "\n", prefix);

   PRINT(n, %d);
   PRINT(nLocal, %d);
   PRINT(numProcs, %d);
   PRINT(procID, %d);

   fprintf(outputFile, "\n// Output and reporting\n");
   PRINT(printLevel, %d);

   fprintf(outputFile, "\n// Solver parameters\n");
   PRINT(numEvals, %d);
   PRINT(aNorm, %e);
   PRINT(eps, %e);
   PRINT(maxBasisSize, %d);
   PRINT(minRestartSize, %d);
   PRINT(maxBlockSize, %d);
   PRINT(maxOuterIterations, %d);
   PRINT(maxMatvecs, %d);

   PRINTIF(target, primme_smallest);
   PRINTIF(target, primme_largest);
   PRINTIF(target, primme_closest_geq);
   PRINTIF(target, primme_closest_leq);
   PRINTIF(target, primme_closest_abs);

   PRINTParamsIF(projection, projection, primme_proj_default);
   PRINTParamsIF(projection, projection, primme_proj_RR);
   PRINTParamsIF(projection, projection, primme_proj_Har);

   PRINTParamsIF(projection, refinedScheme, primme_ref_default);
   PRINTParamsIF(projection, refinedScheme, primme_ref_none);
   PRINTParamsIF(projection, refinedScheme, primme_ref_MultiShifts_QR);
   PRINTParamsIF(projection, refinedScheme, primme_ref_MultiShifts_WTW);
   PRINTParamsIF(projection, refinedScheme, primme_ref_OneShift_QR);
   PRINTParamsIF(projection, refinedScheme, primme_ref_OneShift_WTW);
   PRINTParamsIF(projection, refinedScheme, primme_ref_OneAccuShift_QR);

   PRINTIF(InitBasisMode, primme_init_default);
   PRINTIF(InitBasisMode, primme_init_krylov);
   PRINTIF(InitBasisMode, primme_init_random);
   PRINTIF(InitBasisMode, primme_init_user);
   PRINT(ReIntroInitGuessToBasis, %d);

   PRINT(numTargetShifts, %d);
   if (primme.numTargetShifts > 0) {
      fprintf(outputFile, "%s.targetShifts =", prefix);
      for (i=0; i<primme.numTargetShifts;i++) {
         fprintf(outputFile, " %e",primme.targetShifts[i]);
      }
      fprintf(outputFile, "\n");
   }

   PRINT(dynamicMethodSwitch, %d);
   PRINT(locking, %d);
   PRINT(initSize, %d);
   PRINT(numOrthoConst, %d);
   fprintf(outputFile, "%s.iseed =", prefix);
   for (i=0; i<4;i++) {
      fprintf(outputFile, " %d",primme.iseed[i]);
   }
   fprintf(outputFile, "\n");

   fprintf(outputFile, "\n// Restarting\n");
   PRINTParamsIF(restarting, scheme, primme_thick);
   PRINTParamsIF(restarting, scheme, primme_dtr);

   PRINTParams(restarting, maxPrevRetain, %d);

   fprintf(outputFile, "\n// Correction parameters\n");
   PRINTParams(correction, precondition, %d);
   PRINTParams(correction, robustShifts, %d);
   PRINTParams(correction, maxInnerIterations, %d);
   PRINTParams(correction, relTolBase, %g);

   PRINTParamsIF(correction, convTest, primme_full_LTolerance);
   PRINTParamsIF(correction, convTest, primme_decreasing_LTolerance);
   PRINTParamsIF(correction, convTest, primme_adaptive_ETolerance);
   PRINTParamsIF(correction, convTest, primme_adaptive);

   fprintf(outputFile, "\n// projectors for JD cor.eq.\n");
   PRINTParams(correction, projectors.LeftQ , %d);
   PRINTParams(correction, projectors.LeftX , %d);
   PRINTParams(correction, projectors.RightQ, %d);
   PRINTParams(correction, projectors.SkewQ , %d);
   PRINTParams(correction, projectors.RightX, %d);
   PRINTParams(correction, projectors.SkewX , %d);
   fprintf(outputFile, "// ---------------------------------------------------\n");

#undef PRINT
#undef PRINTIF
#undef PRINTParams
#undef PRINTParamsIF

  /**************************************************************************/
} /* end of display params */
  /**************************************************************************/

/******************************************************************************
 * function 
 * primme_seq_globalSumDouble(void *sendBuf, double *recvBuf, int count) 
 *
 * This is the sequential default for the function globalSumDouble. 
 * If the program is parallel, the user must supply a different pointer 
 * function to the primme.globalSumDouble 
 * 
 ******************************************************************************
 *        NOTE: The count and the copying refers to double datatypes
 ******************************************************************************/

void primme_seq_globalSumDouble(void *sendBuf, void *recvBuf, int *count, 
                      primme_params *params) {

   Num_dcopy_primme(*count, (double *) sendBuf, 1, (double *) recvBuf, 1);

}


/*******************************************************************************
 * Subroutine convTestFunAbsolute - This routine implements primme_params.
 *    convTestFun and return an approximate eigenpair converged when           
 *    resNorm < eps*(aNorm != 0 ? aNorm : aNormEstimate) or
 *    resNorm is close to machineEpsilon * aNorm.          
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * evec         The approximate eigenvector
 * eval         The approximate eigenvalue 
 * rNorm        The norm of the residual vector
 * primme       Structure containing various solver parameters
 *
 * OUTPUT PARAMETERS
 * ----------------------------------
 * isConv      if it isn't zero the approximate pair is marked as converged
 ******************************************************************************/

void convTestFunAbsolute(double *eval, void *evec, double *rNorm, int *isConv,
   primme_params *primme) {

   const double machEps = Num_dlamch_primme("E");
   *isConv = *rNorm < max(
               primme->eps * (
                     primme->aNorm > 0.0 ? primme->aNorm : primme->stats.estimateLargestSVal),
               machEps * 3.16 * primme->stats.estimateLargestSVal);
}
