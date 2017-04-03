/*******************************************************************************
 * Copyright (c) 2017, College of William & Mary
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
#include <string.h>  /* strcmp */  
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
   primme->outputFile                          = stdout;
   primme->printLevel                          = 1;
   primme->stats.numOuterIterations            = 0; 
   primme->stats.numRestarts                   = 0;
   primme->stats.numMatvecs                    = 0;
   primme->stats.numPreconds                   = 0;
   primme->stats.numGlobalSum                  = 0;
   primme->stats.volumeGlobalSum               = 0;
   primme->stats.numOrthoInnerProds            = 0.0;
   primme->stats.elapsedTime                   = 0.0;
   primme->stats.timeMatvec                    = 0.0;
   primme->stats.timePrecond                   = 0.0;
   primme->stats.timeOrtho                     = 0.0;
   primme->stats.timeGlobalSum                 = 0.0;
   primme->stats.estimateMinEVal               = -HUGE_VAL;
   primme->stats.estimateMaxEVal               = HUGE_VAL;
   primme->stats.estimateLargestSVal           = -HUGE_VAL;
   primme->stats.maxConvTol                    = 0.0;
   primme->stats.estimateResidualError         = 0.0;

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
   primme->monitorFun              = NULL;
   primme->monitor                 = NULL;
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
 *        STEEPEST_DESCENT,      : equiv. to GD(block,2*block)
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
   else if (method == PRIMME_STEEPEST_DESCENT) {
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

/*******************************************************************************
 * Subroutine primme_get_member - get the value of a parameter in primme_params
 * 
 * INPUT PARAMETERS
 * ----------------
 * primme  Structure containing various solver parameters and statistics
 * label   reference to the parameter
 *
 * OUTPUT PARAMETERS
 * -----------------
 * value   value of the parameter
 *
 * RETURN
 * ------
 * error code  zero if ok
 *
 ******************************************************************************/

int primme_get_member(primme_params *primme, primme_params_label label,
      void *value) {

   int i;
   union value_t {
      PRIMME_INT int_v;
      void (*matFunc_v) (void *,PRIMME_INT*,void *,PRIMME_INT*,int *,struct primme_params *,int*);
      void *ptr_v;
      void (*globalSumRealFunc_v) (void *,void *,int *,struct primme_params *,int*);
      void (*convTestFun_v)(double *,void*,double*,int*,struct primme_params*,int*);
      primme_target target_v;
      double double_v;
      FILE *file_v;
      primme_init init_v;
      primme_projection projection_v;
      primme_restartscheme restartscheme_v;
      primme_convergencetest convergencetest_v;
      void (*monitorFun_v)(void *basisEvals, int *basisSize, int *basisFlags,
            int *iblock, int *blockSize, void *basisNorms, int *numConverged,
            void *lockedEvals, int *numLocked, int *lockedFlags, void *lockedNorms,
            int *inner_its, void *LSRes, primme_event *event,
            struct primme_params *primme, int *err);
   } *v = (union value_t*)value;

   switch (label) {
      case PRIMME_n:
              v->int_v = primme->n;
      break;
      case PRIMME_matrixMatvec:
              v->matFunc_v = primme->matrixMatvec;
      break;
      case PRIMME_massMatrixMatvec:
              v->matFunc_v = primme->massMatrixMatvec;
      break;
      case PRIMME_applyPreconditioner:
              v->matFunc_v = primme->applyPreconditioner;
      break;
      case PRIMME_numProcs:
              v->int_v = primme->numProcs;
      break;
      case PRIMME_procID:
              v->int_v = primme->procID;
      break;
      case PRIMME_commInfo:
              v->ptr_v = primme->commInfo;
      break;
      case PRIMME_nLocal:
              v->int_v = primme->nLocal;
      break;
      case PRIMME_globalSumReal:
              v->globalSumRealFunc_v = primme->globalSumReal;
      break;
      case PRIMME_numEvals:
              v->int_v = primme->numEvals;
      break;
      case PRIMME_target:
              v->target_v = primme->target;
      break;
      case PRIMME_numTargetShifts:
              v->int_v = primme->numTargetShifts;
      break;
      case PRIMME_targetShifts:
         for (i=0; i< primme->numTargetShifts; i++) {
             (&v->double_v)[i] = primme->targetShifts[i];
         }
      break;
      case PRIMME_locking:
              v->int_v = primme->locking;
      break;
      case PRIMME_initSize:
              v->int_v = primme->initSize;
      break;
      case PRIMME_numOrthoConst:
              v->int_v = primme->numOrthoConst;
      break;
      case PRIMME_dynamicMethodSwitch:
              v->int_v = primme->dynamicMethodSwitch;
      break;
      case PRIMME_maxBasisSize:
              v->int_v = primme->maxBasisSize;
      break;
      case PRIMME_minRestartSize:
              v->int_v = primme->minRestartSize;
      break;
      case PRIMME_maxBlockSize:
              v->int_v = primme->maxBlockSize;
      break;
      case PRIMME_maxMatvecs:
              v->int_v = primme->maxMatvecs;
      break;
      case PRIMME_maxOuterIterations:
              v->int_v = primme->maxOuterIterations;
      break;
      case PRIMME_intWorkSize:
              v->int_v = primme->intWorkSize;
      break;
      case PRIMME_realWorkSize:
              v->int_v = primme->realWorkSize;
      break;
      case PRIMME_iseed:
         for (i=0; i< 4; i++) {
            (&v->int_v)[i] = primme->iseed[i];
         }
      break;
      case PRIMME_intWork:
              v->ptr_v = primme->intWork;
      break;
      case PRIMME_realWork:
              v->ptr_v = primme->realWork;
      break;
      case PRIMME_aNorm:
              v->double_v = primme->aNorm;
      break;
      case PRIMME_eps:
              v->double_v = primme->eps;
      break;
      case PRIMME_printLevel:
              v->int_v = primme->printLevel;
      break;
      case PRIMME_outputFile:
              v->file_v = primme->outputFile;
      break;
      case PRIMME_matrix:
              v->ptr_v = primme->matrix;
      break;
      case PRIMME_preconditioner:
              v->ptr_v = primme->preconditioner;
      break;
      case PRIMME_restartingParams_scheme:
              v->restartscheme_v = primme->restartingParams.scheme;
      break;
      case PRIMME_restartingParams_maxPrevRetain:
              v->int_v = primme->restartingParams.maxPrevRetain;
      break;
      case PRIMME_correctionParams_precondition:
              v->int_v = primme->correctionParams.precondition;
      break;
      case PRIMME_correctionParams_robustShifts:
              v->int_v = primme->correctionParams.robustShifts;
      break;
      case PRIMME_correctionParams_maxInnerIterations:
              v->int_v = primme->correctionParams.maxInnerIterations;
      break;
      case PRIMME_correctionParams_projectors_LeftQ:
              v->int_v = primme->correctionParams.projectors.LeftQ;
      break;
      case PRIMME_correctionParams_projectors_LeftX:
              v->int_v = primme->correctionParams.projectors.LeftX;
      break;
      case PRIMME_correctionParams_projectors_RightQ:
              v->int_v = primme->correctionParams.projectors.RightQ;
      break;
      case PRIMME_correctionParams_projectors_RightX:
              v->int_v = primme->correctionParams.projectors.RightX;
      break;
      case PRIMME_correctionParams_projectors_SkewQ:
              v->int_v = primme->correctionParams.projectors.SkewQ;
      break;
      case PRIMME_correctionParams_projectors_SkewX:
              v->int_v = primme->correctionParams.projectors.SkewX;
      break;
      case PRIMME_correctionParams_convTest:
              v->convergencetest_v = primme->correctionParams.convTest;
      break;
      case PRIMME_correctionParams_relTolBase:
              v->double_v = primme->correctionParams.relTolBase;
      break;
      case PRIMME_stats_numOuterIterations:
              v->int_v = primme->stats.numOuterIterations;
      break;
      case PRIMME_stats_numRestarts:
              v->int_v = primme->stats.numRestarts;
      break;
      case PRIMME_stats_numMatvecs:
              v->int_v = primme->stats.numMatvecs;
      break;
      case PRIMME_stats_numPreconds:
              v->int_v = primme->stats.numPreconds;
      break;
      case PRIMME_stats_numGlobalSum:
              v->int_v = primme->stats.numGlobalSum;
      break;
      case PRIMME_stats_volumeGlobalSum:
              v->int_v = primme->stats.volumeGlobalSum;
      break;
      case PRIMME_stats_numOrthoInnerProds:
              v->double_v = primme->stats.numOrthoInnerProds;
      break;
      case PRIMME_stats_elapsedTime:
              v->double_v = primme->stats.elapsedTime;
      break;
      case PRIMME_stats_timeMatvec:
              v->double_v = primme->stats.timeMatvec;
      break;
      case PRIMME_stats_timePrecond:
              v->double_v = primme->stats.timePrecond;
      break;
      case PRIMME_stats_timeOrtho:
              v->double_v = primme->stats.timeOrtho;
      break;
      case PRIMME_stats_timeGlobalSum:
              v->double_v = primme->stats.timeGlobalSum;
      break;
      case PRIMME_stats_estimateMinEVal:
              v->double_v = primme->stats.estimateMinEVal;
      break;
      case PRIMME_stats_estimateMaxEVal:
              v->double_v = primme->stats.estimateMaxEVal;
      break;
      case PRIMME_stats_estimateLargestSVal:
              v->double_v = primme->stats.estimateLargestSVal;
      break;
      case PRIMME_ldevecs:
              v->int_v = primme->ldevecs;
      break;
      case PRIMME_ldOPs:
              v->int_v = primme->ldOPs;
      break;
      case PRIMME_monitorFun:
              v->monitorFun_v = primme->monitorFun;
      break;
      case PRIMME_monitor:
              v->ptr_v = primme->monitor;
      break;
      default :
      return 1;
   }
   return 0;
}

/*******************************************************************************
 * Subroutine primme_set_member - set the value to a parameter in primme_params
 * 
 * INPUT PARAMETERS
 * ----------------
 * label   reference to the parameter
 * value   value of the parameter
 *
 * INPUT/OUTPUT PARAMETERS
 * -----------------
 * primme  Structure containing various solver parameters and statistics
 *
 * RETURN
 * ------
 * error code  zero if ok
 *
 ******************************************************************************/

int primme_set_member(primme_params *primme, primme_params_label label,
      void *value) {
   int i;

   union value_t {
      PRIMME_INT *int_v;
      void (*matFunc_v) (void *,PRIMME_INT*,void *,PRIMME_INT*,int *,struct primme_params *,int*);
      void *ptr_v;
      void (*globalSumRealFunc_v) (void *,void *,int *,struct primme_params *,int*);
      void (*convTestFun_v)(double *,void*,double*,int*,struct primme_params*,int*);
      primme_target *target_v;
      double *double_v;
      FILE *file_v;
      primme_init *init_v;
      primme_projection *projection_v;
      primme_restartscheme *restartscheme_v;
      primme_convergencetest *convergencetest_v;
      void (*monitorFun_v)(void *basisEvals, int *basisSize, int *basisFlags,
            int *iblock, int *blockSize, void *basisNorms, int *numConverged,
            void *lockedEvals, int *numLocked, int *lockedFlags, void *lockedNorms,
            int *inner_its, void *LSRes, primme_event *event,
            struct primme_params *primme, int *err);
   } v = *(union value_t*)&value;

   switch (label) {
      case PRIMME_n:
              primme->n = *v.int_v;
      break;
      case PRIMME_matrixMatvec:
              primme->matrixMatvec = v.matFunc_v;
      break;
      case PRIMME_massMatrixMatvec:
              primme->massMatrixMatvec = v.matFunc_v;
      break;
      case PRIMME_applyPreconditioner:
              primme->applyPreconditioner = v.matFunc_v;
      break;
      case PRIMME_numProcs:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->numProcs = (int)*v.int_v;
      break;
      case PRIMME_procID:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->procID = (int)*v.int_v;
      break;
      case PRIMME_commInfo:
              primme->commInfo = v.ptr_v;
      break;
      case PRIMME_nLocal:
              primme->nLocal = *v.int_v;
      break;
      case PRIMME_globalSumReal:
              primme->globalSumReal = v.globalSumRealFunc_v;
      break;
      case PRIMME_numEvals:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->numEvals = (int)*v.int_v;
      break;
      case PRIMME_target:
              primme->target = *v.target_v;
      break;
      case PRIMME_numTargetShifts:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->numTargetShifts = (int)*v.int_v;
      break;
      case PRIMME_targetShifts:
              primme->targetShifts = v.double_v;
      break;
      case PRIMME_locking:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->locking = (int)*v.int_v;
      break;
      case PRIMME_initSize:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->initSize = (int)*v.int_v;
      break;
      case PRIMME_numOrthoConst:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->numOrthoConst = (int)*v.int_v;
      break;
      case PRIMME_dynamicMethodSwitch:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->dynamicMethodSwitch = (int)*v.int_v;
      break;
      case PRIMME_maxBasisSize:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->maxBasisSize = (int)*v.int_v;
      break;
      case PRIMME_minRestartSize:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->minRestartSize = (int)*v.int_v;
      break;
      case PRIMME_maxBlockSize:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->maxBlockSize = (int)*v.int_v;
      break;
      case PRIMME_maxMatvecs:
              primme->maxMatvecs = *v.int_v;
      break;
      case PRIMME_maxOuterIterations:
              primme->maxOuterIterations = *v.int_v;
      break;
      case PRIMME_intWorkSize:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->intWorkSize = (int)*v.int_v;
      break;
      case PRIMME_realWorkSize:
              primme->realWorkSize = (size_t)*v.int_v;
      break;
      case PRIMME_iseed:
         for (i=0; i< 4; i++) {
            primme->iseed[i] = v.int_v[i];
         }
      break;
      case PRIMME_intWork:
              primme->intWork = (int*)v.int_v;
      break;
      case PRIMME_realWork:
              primme->realWork = v.ptr_v;
      break;
      case PRIMME_aNorm:
              primme->aNorm = *v.double_v;
      break;
      case PRIMME_eps:
              primme->eps = *v.double_v;
      break;
      case PRIMME_printLevel:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->printLevel = (int)*v.int_v;
      break;
      case PRIMME_outputFile:
              primme->outputFile = v.file_v;
      break;
      case PRIMME_matrix:
              primme->matrix = v.ptr_v;
      break;
      case PRIMME_preconditioner:
              primme->preconditioner = v.ptr_v;
      break;
      case PRIMME_initBasisMode:
              primme->initBasisMode = *v.init_v;
      break;
      case PRIMME_projectionParams_projection:
              primme->projectionParams.projection = *v.projection_v;
      break;
      case PRIMME_restartingParams_scheme:
              primme->restartingParams.scheme = *v.restartscheme_v;
      break;
      case PRIMME_restartingParams_maxPrevRetain:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->restartingParams.maxPrevRetain = (int)*v.int_v;
      break;
      case PRIMME_correctionParams_precondition:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->correctionParams.precondition = (int)*v.int_v;
      break;
      case PRIMME_correctionParams_robustShifts:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->correctionParams.robustShifts = (int)*v.int_v;
      break;
      case PRIMME_correctionParams_maxInnerIterations:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->correctionParams.maxInnerIterations = (int)*v.int_v;
      break;
      case PRIMME_correctionParams_projectors_LeftQ:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->correctionParams.projectors.LeftQ = (int)*v.int_v;
      break;
      case PRIMME_correctionParams_projectors_LeftX:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->correctionParams.projectors.LeftX = (int)*v.int_v;
      break;
      case PRIMME_correctionParams_projectors_RightQ:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->correctionParams.projectors.RightQ = (int)*v.int_v;
      break;
      case PRIMME_correctionParams_projectors_RightX:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->correctionParams.projectors.RightX = (int)*v.int_v;
      break;
      case PRIMME_correctionParams_projectors_SkewQ:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->correctionParams.projectors.SkewQ = (int)*v.int_v;
      break;
      case PRIMME_correctionParams_projectors_SkewX:
              if (*v.int_v > INT_MAX) return 1; else 
              primme->correctionParams.projectors.SkewX = (int)*v.int_v;
      break;
      case PRIMME_correctionParams_convTest:
              primme->correctionParams.convTest = *v.convergencetest_v;
      break;
      case PRIMME_correctionParams_relTolBase:
              primme->correctionParams.relTolBase = *v.double_v;
      break;
      case PRIMME_stats_numOuterIterations:
              primme->stats.numOuterIterations = *v.int_v;
      break;
      case PRIMME_stats_numRestarts:
              primme->stats.numRestarts = *v.int_v;
      break;
      case PRIMME_stats_numMatvecs:
              primme->stats.numMatvecs = *v.int_v;
      break;
      case PRIMME_stats_numPreconds:
              primme->stats.numPreconds = *v.int_v;
      break;
      case PRIMME_stats_volumeGlobalSum:
              primme->stats.volumeGlobalSum = *v.int_v;
      break;
      case PRIMME_stats_numOrthoInnerProds:
              primme->stats.numOrthoInnerProds = *v.double_v;
      break;
      case PRIMME_stats_elapsedTime:
              primme->stats.elapsedTime = *v.double_v;
      break;
      case PRIMME_stats_timeMatvec:
              primme->stats.timeMatvec = *v.double_v;
      break;
      case PRIMME_stats_timePrecond:
              primme->stats.timePrecond = *v.double_v;
      break;
      case PRIMME_stats_timeOrtho:
              primme->stats.timeOrtho = *v.double_v;
      break;
      case PRIMME_stats_timeGlobalSum:
              primme->stats.timeGlobalSum = *v.double_v;
      break;
      case PRIMME_stats_estimateMinEVal:
              primme->stats.estimateMinEVal = *v.double_v;
      break;
      case PRIMME_stats_estimateMaxEVal:
              primme->stats.estimateMaxEVal = *v.double_v;
      break;
      case PRIMME_stats_estimateLargestSVal:
              primme->stats.estimateLargestSVal = *v.double_v;
      break;
      case PRIMME_stats_maxConvTol:
              primme->stats.maxConvTol = *v.double_v;
      break;
      case PRIMME_convTestFun:
              primme->convTestFun = v.convTestFun_v;
      break;
      case PRIMME_ldevecs:
              primme->ldevecs = *v.int_v;
      break;
      case PRIMME_ldOPs:
              primme->ldOPs = *v.int_v;
      break;
      case PRIMME_monitorFun:
              primme->monitorFun = v.monitorFun_v;
      break;
      case PRIMME_monitor:
              primme->monitor = v.ptr_v;
      break;
      default : 
      return 1;
   }
   return 0;
}

/*******************************************************************************
 * Subroutine primme_member_info - return the label value or the label name,
 *    the type and the arity of some member of primme_params.
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------
 * label       reference to the parameter by the value in primme_params_label
 * label_name  reference by the string associated to the parameter
 *
 * OUTPUT PARAMETERS
 * -----------------
 * type        kind of value
 * arity       number of elements
 *
 * RETURN
 * ------
 * error code  zero if ok
 *
 ******************************************************************************/

int primme_member_info(primme_params_label *label_, const char** label_name_,
      primme_type *type, int *arity) {
   primme_params_label label = (primme_params_label)1000;
   const char *label_name;

   /* Quick exit when neither label nor label_name is given */

   if (label_ == NULL && (label_name_ == NULL || *label_name_ == NULL)) {
      return 1;
   }

   /* Get the label from label_name_ and label_name from label_ */

#define IF_IS(F,O) \
   if ((label_name_ && *label_name_ && strcmp(#F, *label_name_) == 0) \
         || (label_ && *label_ == PRIMME_ ## O)) { \
      label = PRIMME_ ## O; \
      label_name = #F; \
   }

   IF_IS(n                            , n);
   IF_IS(matrixMatvec                 , matrixMatvec);
   IF_IS(massMatrixMatvec             , massMatrixMatvec);
   IF_IS(applyPreconditioner          , applyPreconditioner);
   IF_IS(numProcs                     , numProcs);
   IF_IS(procID                       , procID);
   IF_IS(commInfo                     , commInfo);
   IF_IS(nLocal                       , nLocal);
   IF_IS(globalSumReal                , globalSumReal);
   IF_IS(numEvals                     , numEvals);
   IF_IS(target                       , target);
   IF_IS(numTargetShifts              , numTargetShifts);
   IF_IS(targetShifts                 , targetShifts);
   IF_IS(locking                      , locking);
   IF_IS(initSize                     , initSize);
   IF_IS(numOrthoConst                , numOrthoConst);
   IF_IS(dynamicMethodSwitch          , dynamicMethodSwitch);
   IF_IS(maxBasisSize                 , maxBasisSize);
   IF_IS(minRestartSize               , minRestartSize);
   IF_IS(maxBlockSize                 , maxBlockSize);
   IF_IS(maxMatvecs                   , maxMatvecs);
   IF_IS(maxOuterIterations           , maxOuterIterations);
   IF_IS(intWorkSize                  , intWorkSize);
   IF_IS(realWorkSize                 , realWorkSize);
   IF_IS(iseed                        , iseed);
   IF_IS(intWork                      , intWork);
   IF_IS(realWork                     , realWork);
   IF_IS(aNorm                        , aNorm);
   IF_IS(eps                          , eps);
   IF_IS(printLevel                   , printLevel);
   IF_IS(outputFile                   , outputFile);
   IF_IS(matrix                       , matrix);
   IF_IS(preconditioner               , preconditioner);
   IF_IS(initBasisMode                , initBasisMode);
   IF_IS(projection_projection        , projectionParams_projection);
   IF_IS(restarting_scheme            , restartingParams_scheme);
   IF_IS(restarting_maxPrevRetain     , restartingParams_maxPrevRetain);
   IF_IS(correction_precondition      , correctionParams_precondition);
   IF_IS(correction_robustShifts      , correctionParams_robustShifts);
   IF_IS(correction_maxInnerIterations, correctionParams_maxInnerIterations);
   IF_IS(correction_projectors_LeftQ  , correctionParams_projectors_LeftQ);
   IF_IS(correction_projectors_LeftX  , correctionParams_projectors_LeftX);
   IF_IS(correction_projectors_RightQ , correctionParams_projectors_RightQ);
   IF_IS(correction_projectors_RightX , correctionParams_projectors_RightX);
   IF_IS(correction_projectors_SkewQ  , correctionParams_projectors_SkewQ);
   IF_IS(correction_projectors_SkewX  , correctionParams_projectors_SkewX);
   IF_IS(correction_convTest          , correctionParams_convTest);
   IF_IS(correction_relTolBase        , correctionParams_relTolBase);
   IF_IS(stats_numOuterIterations     , stats_numOuterIterations);
   IF_IS(stats_numRestarts            , stats_numRestarts);
   IF_IS(stats_numMatvecs             , stats_numMatvecs);
   IF_IS(stats_numPreconds            , stats_numPreconds);
   IF_IS(stats_numGlobalSum           , stats_numGlobalSum);
   IF_IS(stats_volumeGlobalSum        , stats_volumeGlobalSum);
   IF_IS(stats_numOrthoInnerProds     , stats_numOrthoInnerProds);
   IF_IS(stats_elapsedTime            , stats_elapsedTime);
   IF_IS(stats_timeMatvec             , stats_timeMatvec);
   IF_IS(stats_timePrecond            , stats_timePrecond);
   IF_IS(stats_timeOrtho              , stats_timeOrtho);
   IF_IS(stats_timeGlobalSum          , stats_timeGlobalSum);
   IF_IS(stats_estimateMinEVal        , stats_estimateMinEVal);
   IF_IS(stats_estimateMaxEVal        , stats_estimateMaxEVal);
   IF_IS(stats_estimateLargestSVal    , stats_estimateLargestSVal);
   IF_IS(stats_maxConvTol             , stats_maxConvTol);
   IF_IS(convTestFun                  , convTestFun);
   IF_IS(ldevecs                      , ldevecs);
   IF_IS(ldOPs                        , ldOPs);
   IF_IS(monitorFun                   , monitorFun);
   IF_IS(monitor                      , monitor);
#undef IF_IS

   /* Return label/label_name */

   if (label_) *label_ = label;
   if (label_name_) *label_name_ = label_name;

   /* Return type and arity */

   switch(label) {
      /* members with type int */

      case PRIMME_n:
      case PRIMME_numEvals:
      case PRIMME_target:
      case PRIMME_locking:
      case PRIMME_initSize:
      case PRIMME_numOrthoConst:
      case PRIMME_dynamicMethodSwitch:
      case PRIMME_maxBasisSize:
      case PRIMME_minRestartSize:
      case PRIMME_maxBlockSize:
      case PRIMME_maxMatvecs:
      case PRIMME_maxOuterIterations:
      case PRIMME_initBasisMode:
      case PRIMME_projectionParams_projection:
      case PRIMME_restartingParams_scheme:
      case PRIMME_restartingParams_maxPrevRetain:
      case PRIMME_correctionParams_precondition:
      case PRIMME_correctionParams_robustShifts:
      case PRIMME_correctionParams_maxInnerIterations:
      case PRIMME_correctionParams_projectors_LeftQ:
      case PRIMME_correctionParams_projectors_LeftX:
      case PRIMME_correctionParams_projectors_RightQ:
      case PRIMME_correctionParams_projectors_RightX:
      case PRIMME_correctionParams_projectors_SkewQ:
      case PRIMME_correctionParams_projectors_SkewX:
      case PRIMME_correctionParams_convTest:
      case PRIMME_stats_numOuterIterations:
      case PRIMME_stats_numRestarts:
      case PRIMME_stats_numMatvecs:
      case PRIMME_stats_numPreconds:
      case PRIMME_stats_numGlobalSum:
      case PRIMME_stats_volumeGlobalSum:
      case PRIMME_numProcs:
      case PRIMME_procID:
      case PRIMME_nLocal:
      case PRIMME_numTargetShifts:
      case PRIMME_intWorkSize:
      case PRIMME_realWorkSize:
      case PRIMME_printLevel:
      case PRIMME_ldevecs:
      case PRIMME_ldOPs:
      if (type) *type = primme_int;
      if (arity) *arity = 1;
      break;

      case PRIMME_iseed:
      if (type) *type = primme_int;
      if (arity) *arity = 4;
      break;

      /* members with type double */

      case PRIMME_aNorm:
      case PRIMME_eps:
      case PRIMME_correctionParams_relTolBase:
      case PRIMME_stats_numOrthoInnerProds:
      case PRIMME_stats_timeMatvec:
      case PRIMME_stats_timePrecond:
      case PRIMME_stats_timeOrtho:
      case PRIMME_stats_timeGlobalSum:
      case PRIMME_stats_elapsedTime:
      case PRIMME_stats_estimateMinEVal:
      case PRIMME_stats_estimateMaxEVal:
      case PRIMME_stats_estimateLargestSVal:
      case PRIMME_stats_maxConvTol:
      if (type) *type = primme_double;
      if (arity) *arity = 1;
      break;
 
      case PRIMME_targetShifts:
      if (type) *type = primme_double;
      if (arity) *arity = 0;
      break;

      /* members with type pointer */

      case PRIMME_matrixMatvec:
      case PRIMME_applyPreconditioner:
      case PRIMME_commInfo:
      case PRIMME_globalSumReal:
      case PRIMME_intWork:
      case PRIMME_realWork:
      case PRIMME_massMatrixMatvec:
      case PRIMME_outputFile:
      case PRIMME_matrix:
      case PRIMME_preconditioner:
      case PRIMME_convTestFun:
      case PRIMME_monitorFun:
      case PRIMME_monitor:
      if (type) *type = primme_pointer;
      if (arity) *arity = 1;
      break;

      default: 
      return 1;
   }

   return 0;
}
 
/*******************************************************************************
 * Subroutine primme_constant_info - return the value of a primme enum constant.
 * 
 * INPUT PARAMETERS
 * ----------------
 * label_name  name of the constant
 *
 * OUTPUT PARAMETERS
 * -----------------
 * value       numerical value of the constant
 *
 * RETURN
 * ------
 * error code  zero if ok
 *
 ******************************************************************************/

int primme_constant_info(const char* label_name, int *value) {

#define IF_IS(F) if (strcmp(#F, label_name) == 0) {*value = (int)F; return 0;}

   /* primme_preset_method */

   IF_IS(PRIMME_DEFAULT_METHOD);
   IF_IS(PRIMME_DYNAMIC);
   IF_IS(PRIMME_DEFAULT_MIN_TIME);
   IF_IS(PRIMME_DEFAULT_MIN_MATVECS);
   IF_IS(PRIMME_Arnoldi);
   IF_IS(PRIMME_GD);
   IF_IS(PRIMME_GD_plusK);
   IF_IS(PRIMME_GD_Olsen_plusK);
   IF_IS(PRIMME_JD_Olsen_plusK);
   IF_IS(PRIMME_RQI);
   IF_IS(PRIMME_JDQR);
   IF_IS(PRIMME_JDQMR);
   IF_IS(PRIMME_JDQMR_ETol);
   IF_IS(PRIMME_STEEPEST_DESCENT);
   IF_IS(PRIMME_LOBPCG_OrthoBasis);
   IF_IS(PRIMME_LOBPCG_OrthoBasis_Window);
   
   /* enum members for targeting; restarting and innertest */
   
   IF_IS(primme_smallest);
   IF_IS(primme_largest);
   IF_IS(primme_closest_geq);
   IF_IS(primme_closest_leq);
   IF_IS(primme_closest_abs);
   IF_IS(primme_largest_abs);
   IF_IS(primme_proj_default);
   IF_IS(primme_proj_RR);
   IF_IS(primme_proj_harmonic);
   IF_IS(primme_proj_refined);
   IF_IS(primme_init_default);
   IF_IS(primme_init_krylov);
   IF_IS(primme_init_random);
   IF_IS(primme_init_user);
   IF_IS(primme_thick);
   IF_IS(primme_dtr);
   IF_IS(primme_full_LTolerance);
   IF_IS(primme_decreasing_LTolerance);
   IF_IS(primme_adaptive_ETolerance);
   IF_IS(primme_adaptive);

   /* enum member for event */

   IF_IS(primme_event_outer_iteration);
   IF_IS(primme_event_inner_iteration);
   IF_IS(primme_event_restart);
   IF_IS(primme_event_reset);
   IF_IS(primme_event_converged);
   IF_IS(primme_event_locked);
#undef IF_IS

   /* return error if label not found */

   return 1;   
}

#endif /* USE_DOUBLE */
