/*******************************************************************************
 * Copyright (c) 2016, College of William & Mary
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the College of William & Mary nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * PRIMME: https://github.com/primme/primme
 * Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
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
#  include <malloc.h>
#endif
#include <stdlib.h>   /* mallocs, free */
#include <stdio.h>    
#include <math.h>    
#include <limits.h>    
#include "template.h"
#include "primme_interface.h"
#include "const.h"

/* Only define these functions ones */
#ifdef USE_DOUBLE
#include "notemplate.h"

/*******************************************************************************
 * Subroutine primme_initialize - Set primme_params members to default values.
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 * primme  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_initialize(primme_params *primme) {

   /* Essential parameters */
   primme->n                       = 0;
   primme->numEvals                = 1;
   primme->target                  = primme_smallest;
   primme->aNorm                   = 0.0L;
   primme->eps                     = 0.0;

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
   primme->globalSumReal           = NULL;

   /* Initial guesses/constraints */
   primme->initSize                = 0;
   primme->numOrthoConst           = 0;

   primme->projectionParams.projection = primme_proj_default;

   primme->initBasisMode                       = primme_init_default;

   /* Eigensolver parameters (outer) */
   primme->locking                             = -1;
   primme->dynamicMethodSwitch                 = -1;
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
   primme->stats.volumeGlobalSum   = 0;
   primme->stats.numOrthoInnerProds= 0.0;
   primme->stats.elapsedTime       = 0.0;
   primme->stats.timeMatvec        = 0.0;
   primme->stats.timePrecond       = 0.0;
   primme->stats.timeGlobalSum     = 0.0;
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
   primme->ShiftsForPreconditioner = NULL;
   primme->convTestFun             = NULL;
   primme->ldevecs                 = 0;
   primme->ldOPs                   = 0;

}

/*******************************************************************************
 * Subroutine primme_free - Free memory resources allocated by Sprimme.
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 * primme  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_free(primme_params *primme) {

   free(primme->intWork);
   free(primme->realWork);
   primme->intWorkSize  = 0;
   primme->realWorkSize = 0;

}

/******************************************************************************
 * Subroutine primme_set_method - Set the eigensolver parameters to implement a
 *    method requested by the user.
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
 *    primme    The main structure to be used by Primme, with parameters
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
int primme_set_method(primme_preset_method method, primme_params *primme) {

   /* Set default method as DYNAMIC */
   if (method == PRIMME_DEFAULT_METHOD)
      method = PRIMME_DYNAMIC;

   /* From our experience, these two methods yield the smallest matvecs/time */
   /* DYNAMIC will make some timings before it settles on one of the two     */
   if (method == PRIMME_DEFAULT_MIN_MATVECS) {
      method = PRIMME_GD_Olsen_plusK;
   }
   else if (method == PRIMME_DEFAULT_MIN_TIME) {
      /* JDQMR works better than JDQMR_ETol in interior problems. */
      if (primme->target == primme_smallest || primme->target == primme_largest) {
         method = PRIMME_JDQMR_ETol;
      }
      else {
         method = PRIMME_JDQMR;
      }
   }
   if (method == PRIMME_DYNAMIC) {
      /* JDQMR works better than JDQMR_ETol in interior problems. */
      if (primme->target == primme_smallest || primme->target == primme_largest) {
         method = PRIMME_JDQMR_ETol;
      }
      else {
         method = PRIMME_JDQMR;
      }
      primme->dynamicMethodSwitch = 1;
   }
   else {
      primme->dynamicMethodSwitch = 0;
   }

   if (primme->maxBlockSize == 0) {
      primme->maxBlockSize = 1;
   }
   if (primme->correctionParams.precondition == -1) {
      primme->correctionParams.precondition = primme->applyPreconditioner ? 1 : 0;
   }

   if (method == PRIMME_Arnoldi) {
      primme->restartingParams.maxPrevRetain      = 0;
      primme->correctionParams.precondition       = 0;
      primme->correctionParams.maxInnerIterations = 0;
   }
   else if (method == PRIMME_GD) {
      primme->restartingParams.maxPrevRetain      = 0;
      primme->correctionParams.robustShifts       = 1;
      primme->correctionParams.maxInnerIterations = 0;
      primme->correctionParams.projectors.RightX  = 0;
      primme->correctionParams.projectors.SkewX   = 0;
   }
   else if (method == PRIMME_GD_plusK) {
      if (primme->restartingParams.maxPrevRetain <= 0) {
         if (primme->maxBlockSize == 1 && primme->numEvals > 1) {
            primme->restartingParams.maxPrevRetain = 2;
         }
         else {
            primme->restartingParams.maxPrevRetain = primme->maxBlockSize;
         }
      }
      primme->correctionParams.maxInnerIterations = 0;
      primme->correctionParams.projectors.RightX  = 0;
      primme->correctionParams.projectors.SkewX   = 0;
   }
   else if (method == PRIMME_GD_Olsen_plusK) {
      if (primme->restartingParams.maxPrevRetain <= 0) {
         if (primme->maxBlockSize == 1 && primme->numEvals > 1) {
            primme->restartingParams.maxPrevRetain = 2;
         }
         else {
            primme->restartingParams.maxPrevRetain = primme->maxBlockSize;
         }
      }
      primme->correctionParams.maxInnerIterations = 0;
      primme->correctionParams.projectors.RightX  = 1;
      primme->correctionParams.projectors.SkewX   = 0;
   }
   else if (method == PRIMME_JD_Olsen_plusK) {
      if (primme->restartingParams.maxPrevRetain <= 0) {
         if (primme->maxBlockSize == 1 && primme->numEvals > 1) {
            primme->restartingParams.maxPrevRetain = 2;
         }
         else {
            primme->restartingParams.maxPrevRetain = primme->maxBlockSize;
         }
      }
      primme->correctionParams.robustShifts       = 1;
      primme->correctionParams.maxInnerIterations = 0;
      primme->correctionParams.projectors.RightX  = 1;
      primme->correctionParams.projectors.SkewX   = 1;
   }
   else if (method == PRIMME_RQI) {
     /* This is accelerated RQI as basisSize may be > 2                 */
     /* NOTE: If numTargetShifts > 0 and targetShifts[*] are provided,  *
      * the interior problem solved uses these shifts in the correction *
      * equation. Therefore RQI becomes INVIT in that case.             */
      primme->locking                             = 1;
      primme->restartingParams.maxPrevRetain      = 0;
      primme->correctionParams.robustShifts       = 1;
      primme->correctionParams.maxInnerIterations = -1;
      primme->correctionParams.projectors.LeftQ   = 1;
      primme->correctionParams.projectors.LeftX   = 1;
      primme->correctionParams.projectors.RightQ  = 0;
      primme->correctionParams.projectors.RightX  = 1;
      primme->correctionParams.projectors.SkewQ   = 0;
      primme->correctionParams.projectors.SkewX   = 0;
      primme->correctionParams.convTest           = primme_full_LTolerance;
   }
   else if (method == PRIMME_JDQR) {
      primme->locking                             = 1;
      primme->restartingParams.maxPrevRetain      = 1;
      primme->correctionParams.robustShifts       = 0;
      if (primme->correctionParams.maxInnerIterations == -INT_MAX) {
         primme->correctionParams.maxInnerIterations = 10;
      }
      primme->correctionParams.projectors.LeftQ   = 0;
      primme->correctionParams.projectors.LeftX   = 1;
      primme->correctionParams.projectors.RightQ  = 1;
      primme->correctionParams.projectors.RightX  = 1;
      primme->correctionParams.projectors.SkewQ   = 1;
      primme->correctionParams.projectors.SkewX   = 1;
      primme->correctionParams.relTolBase         = 1.5;
      primme->correctionParams.convTest           = primme_full_LTolerance;
   }
   else if (method == PRIMME_JDQMR) {
      if (primme->restartingParams.maxPrevRetain < 0) {
         primme->restartingParams.maxPrevRetain   = 1;
      }
      primme->correctionParams.maxInnerIterations = -1;
      if (primme->correctionParams.precondition) {
         primme->correctionParams.projectors.LeftQ   = 1;
      }
      else {
         primme->correctionParams.projectors.LeftQ   = 0;
      }
      primme->correctionParams.projectors.LeftX   = 1;
      primme->correctionParams.projectors.RightQ  = 0;
      primme->correctionParams.projectors.RightX  = 0;
      primme->correctionParams.projectors.SkewQ   = 0;
      primme->correctionParams.projectors.SkewX   = 1;
      primme->correctionParams.convTest           = primme_adaptive;
   }
   else if (method == PRIMME_JDQMR_ETol) {
      if (primme->restartingParams.maxPrevRetain < 0) {
         primme->restartingParams.maxPrevRetain   = 1;
      }
      primme->correctionParams.maxInnerIterations = -1;
      if (primme->correctionParams.precondition) {
         primme->correctionParams.projectors.LeftQ   = 1;
      }
      else {
         primme->correctionParams.projectors.LeftQ   = 0;
      }
      primme->correctionParams.projectors.LeftX   = 1;
      primme->correctionParams.projectors.RightQ  = 0;
      primme->correctionParams.projectors.RightX  = 0;
      primme->correctionParams.projectors.SkewQ   = 0;
      primme->correctionParams.projectors.SkewX   = 1;
      primme->correctionParams.convTest           = primme_adaptive_ETolerance;
   }
   else if (method == PRIMME_SUBSPACE_ITERATION) {
      primme->locking                             = 1;
      primme->maxBasisSize                        = primme->numEvals*2;
      primme->minRestartSize                      = primme->numEvals;
      primme->maxBlockSize                        = primme->numEvals;
      primme->restartingParams.scheme             = primme_thick;
      primme->restartingParams.maxPrevRetain      = 0;
      primme->correctionParams.robustShifts       = 0;
      primme->correctionParams.maxInnerIterations = 0;
      primme->correctionParams.projectors.RightX  = 1;
      primme->correctionParams.projectors.SkewX   = 0;
   }
   else if (method == PRIMME_LOBPCG_OrthoBasis) {
      primme->maxBasisSize                        = primme->numEvals*3;
      primme->minRestartSize                      = primme->numEvals;
      primme->maxBlockSize                        = primme->numEvals;
      primme->restartingParams.maxPrevRetain      = primme->numEvals;
      primme->restartingParams.scheme             = primme_thick;
      primme->correctionParams.robustShifts       = 0;
      primme->correctionParams.maxInnerIterations = 0;
      primme->correctionParams.projectors.RightX  = 1;
      primme->correctionParams.projectors.SkewX   = 0;
   }
   else if (method == PRIMME_LOBPCG_OrthoBasis_Window) {
      /* Observed needing to restart with two vectors at least to converge    */
      /* in some tests like for instance testi-10-LOBPCG_OrthoBasis_Window-3- */
      /* primme_closest_geq-primme_proj_refined.F                             */
      if (primme->maxBlockSize == 1
            && (primme->target == primme_closest_leq
               || primme->target == primme_closest_geq)) {
         primme->maxBasisSize                        = 4;
         primme->minRestartSize                      = 2;
         primme->restartingParams.maxPrevRetain      = 1;
      }
      else {
         primme->maxBasisSize                        = primme->maxBlockSize*3;
         primme->minRestartSize                      = primme->maxBlockSize;
         primme->restartingParams.maxPrevRetain      = primme->maxBlockSize;
      }
      primme->restartingParams.scheme             = primme_thick;
      primme->correctionParams.robustShifts       = 0;
      primme->correctionParams.maxInnerIterations = 0;
      primme->correctionParams.projectors.RightX  = 1;
      primme->correctionParams.projectors.SkewX   = 0;
   }
   else {
      return -1;
   }

   primme_set_defaults(primme);

   return 0;

}

/******************************************************************************
 * Subroutine primme_set_defaults - Set valid values for options that still
 *    has the initial invalid value set in primme_initialize.
 *
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 * primme  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_set_defaults(primme_params *primme) {
   if (primme->dynamicMethodSwitch < 0) {
      primme_set_method(PRIMME_DYNAMIC, primme);
   }

   /* ----------------------------------------- */
   /* Set some defaults for sequential programs */
   /* ----------------------------------------- */
   if (primme->numProcs <= 1) {
      primme->nLocal = primme->n;
      primme->procID = 0;
   }

   if (primme->ldevecs == 0)
      primme->ldevecs = primme->nLocal;
   if (primme->projectionParams.projection == primme_proj_default)
      primme->projectionParams.projection = primme_proj_RR;
   if (primme->initBasisMode == primme_init_default)
      primme->initBasisMode = primme_init_krylov;

   /* If we are free to choose the leading dimension of V and W, use    */
   /* a multiple of PRIMME_BLOCK_SIZE. This may improve the performance */
   /* of Num_update_VWXR_Sprimme.                                       */

   if (primme->ldOPs == 0) {
      primme->ldOPs = min(((primme->nLocal + PRIMME_BLOCK_SIZE - 1)
               /PRIMME_BLOCK_SIZE)*PRIMME_BLOCK_SIZE, primme->nLocal);
   }
      
   /* Now that most of the parameters have been set, set defaults  */
   /* for basisSize, restartSize (for those methods that need it)  */
   /* For interior, larger basisSize and restartSize are advisable */
   /* If maxBlockSize is provided, assign at least 4*blocksize     */
   /* and consider also minRestartSize and maxPrevRetain           */
   if (primme->maxBasisSize == 0) {
      if (primme->target==primme_smallest || primme->target==primme_largest)
         primme->maxBasisSize   = min(primme->n, max(
            max(15, 4*primme->maxBlockSize+primme->restartingParams.maxPrevRetain), 
            (int) 2.5*primme->minRestartSize+primme->restartingParams.maxPrevRetain));
      else
         primme->maxBasisSize   = min(primme->n, max(
            max(35, 5*primme->maxBlockSize+primme->restartingParams.maxPrevRetain),
            (int) 1.7*primme->minRestartSize+primme->restartingParams.maxPrevRetain));
   }

   if (primme->minRestartSize == 0) {
      if (primme->target==primme_smallest || primme->target==primme_largest)
         primme->minRestartSize = (int) (0.5 + 0.4*primme->maxBasisSize);
      else
         primme->minRestartSize = (int) (0.5 + 0.6*primme->maxBasisSize);

      /* Adjust so that an integer number of blocks are added between restarts*/
      /* restart=basis-block*ceil((basis-restart-prevRetain)/block)-prevRetain*/
      if (primme->maxBlockSize > 1) {
         if (primme->restartingParams.maxPrevRetain > 0) 
            primme->minRestartSize = primme->maxBasisSize-primme->maxBlockSize*
            (1 + (int) ((primme->maxBasisSize - primme->minRestartSize - 1 
                         -primme->restartingParams.maxPrevRetain ) / (double) 
            primme->maxBlockSize) ) - primme->restartingParams.maxPrevRetain ;
         else 
            primme->minRestartSize = primme->maxBasisSize-primme->maxBlockSize*
            (1 + (int) ((primme->maxBasisSize - primme->minRestartSize - 1)
                        /(double) primme->maxBlockSize) );
      }
   }

   /* --------------------------------------------------------------------- */
   /* Decide on whether to use locking (hard locking), or not (soft locking)*/
   /* --------------------------------------------------------------------- */
   if (primme->locking >= 0) {
      /* Honor the user setup (do nothing) */
   }
   else if (primme->target != primme_smallest && primme->target != primme_largest) {
       primme->locking = 1;
   }
   else if (primme->numEvals > primme->minRestartSize) {
      /* use locking when not enough vectors to restart with */
      primme->locking = 1;
   }
   else {
      primme->locking = 0;   
   }
}

/*******************************************************************************
 * Subroutine primme_display_params - Displays the current configuration of
 *    primme data structure.
 * 
 * INPUT PARAMETERS
 * ----------------------------------
 * primme  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_display_params(primme_params primme) {

   fprintf(primme.outputFile,
           "// ---------------------------------------------------\n"
           "//                 primme configuration               \n"
           "// ---------------------------------------------------\n");

   primme_display_params_prefix("primme", primme);
   fflush(primme.outputFile);
}

/*******************************************************************************
 * Subroutine primme_params_prefix - Displays the current configuration of
 *    primme data structure. The options are printed as
 *    prefix.optionName = value.
 * 
 * INPUT PARAMETERS
 * ----------------------------------
 * prefix  Name of the structure
 * primme  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_display_params_prefix(const char* prefix, primme_params primme) {

   int i;
   FILE *outputFile = primme.outputFile;

#define PRINT(P,L) fprintf(outputFile, "%s." #P " = " #L "\n", prefix, primme. P);
#define PRINTIF(P,V) if (primme. P == V) fprintf(outputFile, "%s." #P " = " #V "\n", prefix);
#define PRINTParams(P,S,L) fprintf(outputFile, "%s." #P "." #S " = " #L "\n", \
                                    prefix, primme. P ## Params. S);
#define PRINTParamsIF(P,S,V) if (primme. P ## Params. S == V) \
                                 fprintf(outputFile, "%s." #P "." #S " = " #V "\n", prefix);
#define PRINT_PRIMME_INT(P) fprintf(outputFile, "%s." #P " = %" PRIMME_INT_P "\n", prefix, primme. P);
 
   PRINT_PRIMME_INT(n);
   PRINT_PRIMME_INT(nLocal);
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
   PRINT_PRIMME_INT(maxOuterIterations);
   PRINT_PRIMME_INT(maxMatvecs);

   PRINTIF(target, primme_smallest);
   PRINTIF(target, primme_largest);
   PRINTIF(target, primme_closest_geq);
   PRINTIF(target, primme_closest_leq);
   PRINTIF(target, primme_closest_abs);
   PRINTIF(target, primme_largest_abs);

   PRINTParamsIF(projection, projection, primme_proj_default);
   PRINTParamsIF(projection, projection, primme_proj_RR);
   PRINTParamsIF(projection, projection, primme_proj_harmonic);
   PRINTParamsIF(projection, projection, primme_proj_refined);

   PRINTIF(initBasisMode, primme_init_default);
   PRINTIF(initBasisMode, primme_init_krylov);
   PRINTIF(initBasisMode, primme_init_random);
   PRINTIF(initBasisMode, primme_init_user);

   PRINT(numTargetShifts, %d);
   if (primme.numTargetShifts > 0 && primme.targetShifts) {
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
   PRINT_PRIMME_INT(ldevecs);
   PRINT_PRIMME_INT(ldOPs);
   fprintf(outputFile, "%s.iseed =", prefix);
   for (i=0; i<4;i++) {
      fprintf(outputFile, " %" PRIMME_INT_P, primme.iseed[i]);
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
}

#endif /* USE_DOUBLE */
