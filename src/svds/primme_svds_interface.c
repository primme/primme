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
 * File: primme_svds_interface.c
 *
 * Purpose - Contains interface functions to PRIMME SVDS named primme_svds_*.
 *
 ******************************************************************************/

#include <stdlib.h>   /* mallocs, free */
#include <stdio.h>    
#include <math.h>    
#include <limits.h>    
#include "numerical.h"
#include "primme_svds_interface.h"
#include "primme_interface.h"

/* Only define these functions ones */
#ifdef USE_DOUBLE
#include "notemplate.h"

static void copy_params_from_svds(primme_svds_params *primme_svds, int stage);
static void globalSumRealSvds(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme, int *ierr);

/*******************************************************************************
 * Subroutine primme__svdsinitialize - Set primme_svds_params members to default
 *    values.
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 * primme_svds  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_svds_initialize(primme_svds_params *primme_svds) {

   /* Essential parameters */
   primme_svds->m                       = 0;
   primme_svds->n                       = 0;
   primme_svds->numSvals                = 1;
   primme_svds->target                  = primme_svds_largest;
   primme_svds->method                  = primme_svds_op_none;
   primme_svds->methodStage2            = primme_svds_op_none;

   /* Shifts for primme_svds_augmented method */
   primme_svds->numTargetShifts         = 0;
   primme_svds->targetShifts            = NULL;

   /* Parallel computing parameters */
   primme_svds->numProcs                = 1;
   primme_svds->procID                  = 0;
   primme_svds->mLocal                  = 0;
   primme_svds->nLocal                  = 0;
   primme_svds->commInfo                = NULL;
   primme_svds->globalSumReal           = NULL;

   /* Use these pointers to provide matrix/preconditioner */
   primme_svds->matrix                  = NULL;
   primme_svds->preconditioner          = NULL;

   /* Matvec and preconditioner */
   primme_svds->matrixMatvec            = NULL;
   primme_svds->applyPreconditioner     = NULL;

   /* Other important parameters users may set */
   primme_svds->aNorm                   = 0.0L;
   primme_svds->eps                     = 0.0;
   primme_svds->precondition            = -1;
   primme_svds->initSize                = 0;
   primme_svds->maxBasisSize            = 0;
   primme_svds->maxBlockSize            = 0;
   primme_svds->maxMatvecs              = INT_MAX;
   primme_svds->printLevel              = 1;
   primme_svds->outputFile              = stdout;
   primme_svds->locking                 = -1;
   primme_svds->numOrthoConst           = 0;

   /* Reporting performance */
   primme_svds->stats.numOuterIterations= 0;
   primme_svds->stats.numRestarts       = 0;
   primme_svds->stats.numMatvecs        = 0;
   primme_svds->stats.numPreconds       = 0;
   primme_svds->stats.elapsedTime       = 0.0L;

   /* Internally used variables */
   primme_svds->iseed[0] = -1;   /* To set iseed, we first need procID           */ 
   primme_svds->iseed[1] = -1;   /* Thus we set all iseeds to -1                 */
   primme_svds->iseed[2] = -1;   /* Unless users provide their own iseeds,       */
   primme_svds->iseed[3] = -1;   /* PRIMME will set thse later uniquely per proc */
   primme_svds->intWorkSize             = 0;
   primme_svds->realWorkSize            = 0;
   primme_svds->intWork                 = NULL;
   primme_svds->realWork                = NULL;

   primme_initialize(&primme_svds->primme);
   primme_initialize(&primme_svds->primmeStage2);
 
}

/******************************************************************************
 * Subroutine primme_svds_set_method - Set the approach to solve the singular
 *    value problem.
 *
 * INPUT PARAMETERS
 * ----------------
 *    method   singular value method, one of:
 *
 *       primme_svds_default, currently set as primme_svds_hybrid.
 *       primme_svds_normalequations, compute the eigenvectors of A'*A or A*A'.
 *       primme_svds_augmented|, compute the eigenvectors of the augmented
 *          matrix, [zeros() A'; A zeros()].
 *       primme_svds_hybrid, start with primme_svds_normalequations; use the
 *       resulting approximate singular vectors as initial vectors for
 *       primme_svds_augmented if the required accuracy was not achieved.
 *
 *    methodStage1: preset method to compute the eigenpairs at the first stage.
 *
 *    methodStage2: preset method to compute the eigenpairs with the second
 *       stage of primme_svds_hybrid.
 *
 * INPUT/OUTPUT PARAMATERS
 * -----------------------
 *    primme_svds   parameters structure 
 *
 ******************************************************************************/

int primme_svds_set_method(primme_svds_preset_method method,
      primme_preset_method methodStage1, primme_preset_method methodStage2,
      primme_svds_params *primme_svds) {

   /* Set method and methodStage2 in primme_svds_params */
   switch(method) {
   case primme_svds_default:
   case primme_svds_hybrid:
      primme_svds->method = (primme_svds->n <= primme_svds->m ?
            primme_svds_op_AtA : primme_svds_op_AAt);
      primme_svds->methodStage2 = primme_svds_op_augmented;
      break;
   case primme_svds_normalequations:
      primme_svds->method = (primme_svds->n <= primme_svds->m ?
            primme_svds_op_AtA : primme_svds_op_AAt);
      primme_svds->methodStage2 = primme_svds_op_none;
      break;
   case primme_svds_augmented:
      primme_svds->method = primme_svds_op_augmented;
      primme_svds->methodStage2 = primme_svds_op_none;
      break;
   }

   /* Setup underneath eigensolvers based on primme_svds configuration */
   primme_svds_set_defaults(primme_svds);

   /* Set method for the first stage */
   primme_set_method(methodStage1, &primme_svds->primme);

   /* Set method for the second state */
   if (methodStage2 == PRIMME_DEFAULT_METHOD
         && primme_svds->target != primme_svds_largest)
      methodStage2 = PRIMME_JDQMR;
   if (primme_svds->methodStage2 != primme_svds_op_none)
      primme_set_method(methodStage2, &primme_svds->primmeStage2);

   return 0;
}

/******************************************************************************
 * Subroutine primme_svds_set_defaults - Set valid values for options that still
 *    has the initial invalid value set in primme_svds_initialize.
 *
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 *    primme_svds   parameters structure 
 *
 ******************************************************************************/

void primme_svds_set_defaults(primme_svds_params *primme_svds) {

   /* Set defaults for sequential programs */
   if (primme_svds->numProcs <= 1) {
      primme_svds->mLocal = primme_svds->m;
      primme_svds->nLocal = primme_svds->n;
      primme_svds->procID = 0;
      primme_svds->numProcs = 1;
   }

   /* Set svds method if none set */
   if (primme_svds->method == primme_svds_op_none) {
      primme_svds_set_method(primme_svds_default, PRIMME_DEFAULT_METHOD,
            PRIMME_DEFAULT_METHOD, primme_svds);
   }

   /* Copy values set in primme_svds to the first stage underneath eigensolver */
   copy_params_from_svds(primme_svds, 0);

   if (primme_svds->methodStage2 != primme_svds_op_none) {
      /* Copy values set in primme_svds to the second stage underneath eigensolver */
      copy_params_from_svds(primme_svds, 1);
   }
}

/*******************************************************************************
 * Subroutine copy_params_from_svds - Transfer some options in primme_svds to
 *    the primme_params for the given stage.
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 * stage        First stage (0) or second stage (1)
 * primme_svds  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

static void copy_params_from_svds(primme_svds_params *primme_svds, int stage) {
   primme_params *primme;
   primme_svds_operator method;

   primme = stage == 0 ? &primme_svds->primme : &primme_svds->primmeStage2;
   method = stage == 0 ? primme_svds->method : primme_svds->methodStage2;

   if (method == primme_svds_op_none) {
      primme->maxMatvecs = 1;
      return;
   }

   /* -----------------------------------------------*/
   /* Set important parameters for primme structure  */
   /* -----------------------------------------------*/
   primme->numEvals = primme_svds->numSvals;
   if (primme_svds->aNorm > 0.0) {
      switch(method) {
      case primme_svds_op_AtA:
      case primme_svds_op_AAt:
         primme->aNorm = primme_svds->aNorm*primme_svds->aNorm;
         break;
      case primme_svds_op_augmented:
         primme->aNorm = primme_svds->aNorm*sqrt(2.0);
         break;
      case primme_svds_op_none:
         break;
      }
   }
   primme->eps = primme_svds->eps;
   primme->initSize = primme_svds->initSize;
   if (primme_svds->maxBasisSize > 0)
      primme->maxBasisSize = primme_svds->maxBasisSize;
   if (primme_svds->maxBlockSize > 0)
      primme->maxBlockSize = primme_svds->maxBlockSize;
   primme->maxMatvecs = primme_svds->maxMatvecs;
   primme->printLevel = primme_svds->printLevel;
   primme->outputFile = primme_svds->outputFile;
   primme->numOrthoConst = primme_svds->numOrthoConst;

   /* ---------------------------------------------- */
   /* Set some parameters only for parallel programs */
   /* ---------------------------------------------- */
   if (primme_svds->numProcs > 1 && primme_svds->globalSumReal != NULL) {
      primme->procID = primme_svds->procID;
      primme->numProcs = primme_svds->numProcs;
      primme->commInfo = primme_svds->commInfo;
      primme->globalSumReal = globalSumRealSvds;
   }

   switch(method) {
   case primme_svds_op_AtA:
      primme->n = primme_svds->n;
      primme->nLocal = primme_svds->nLocal;
      break;
   case primme_svds_op_AAt:
      primme->n = primme_svds->m;
      primme->nLocal = primme_svds->mLocal;
      break;
   case primme_svds_op_augmented:
      primme->n = primme_svds->m + primme_svds->n;
      primme->nLocal = primme_svds->mLocal + primme_svds->nLocal;
      break;
   case primme_svds_op_none:
      break;
   }

   switch (primme_svds->target) {
   case primme_svds_largest:
      primme->target = primme_largest;
      break;
   case primme_svds_smallest:
      primme->target = (method == primme_svds_op_augmented) ?
         primme_closest_geq : primme_smallest;
      break;
   case primme_svds_closest_abs:
      primme->target = primme_closest_abs;
      primme->numTargetShifts = primme_svds->numTargetShifts;
      break;
   }

   if (stage == 1 && primme->initBasisMode == primme_init_default) {
      primme->initBasisMode = primme_init_user;
   }

   if (((method == primme_svds_op_augmented && 
         primme_svds->target != primme_svds_largest) ||
        primme_svds->target == primme_svds_closest_abs) &&
         primme->projectionParams.projection == primme_proj_default) {
      /* NOTE: refined extraction seems to work better than RR */
      primme->projectionParams.projection = primme_proj_refined;
   }

   if (primme_svds->locking >= 0) {
      primme->locking = primme_svds->locking;
   }

   if (primme_svds->precondition >= 0) {
      primme->correctionParams.precondition = primme_svds->precondition;
   }
   else if (primme->correctionParams.precondition < 0) {
      primme->correctionParams.precondition = primme_svds->applyPreconditioner ? 1 : 0;
   }

}

/*******************************************************************************
 * Subroutine primme_svds_display_params - Displays the current configuration of
 *    primme svds data structure.
 * 
 * INPUT PARAMETERS
 * ----------------------------------
 * primme_svds  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_svds_display_params(primme_svds_params primme_svds) {

   int i;
   FILE *outputFile = primme_svds.outputFile;

#define PRINT(P,L) fprintf(outputFile, "primme_svds." #P " = " #L "\n", primme_svds. P);
#define PRINTIF(P,V) if (primme_svds. P == V) fprintf(outputFile, "primme_svds." #P " = " #V "\n");
#define PRINTParams(P,S,L) fprintf(outputFile, "primme_svds." #P "." #S " = " #L "\n", \
                                    primme_svds. P ## Params. S);
#define PRINTParamsIF(P,S,V) if (primme_svds. P ## Params. S == V) \
                                 fprintf(outputFile, "primme_svds." #P "." #S " = " #V "\n");
#define PRINT_PRIMME_INT(P) fprintf(outputFile, "primme_svds." #P " = %" PRIMME_INT_P "\n", primme_svds. P);
 
fprintf(outputFile, "// ---------------------------------------------------\n"
                    "//            primme_svds configuration               \n"
                    "// ---------------------------------------------------\n");

   PRINT_PRIMME_INT(m);
   PRINT_PRIMME_INT(n);
   PRINT_PRIMME_INT(mLocal);
   PRINT_PRIMME_INT(nLocal);
   PRINT(numProcs, %d);
   PRINT(procID, %d);

   fprintf(outputFile, "\n// Output and reporting\n");
   PRINT(printLevel, %d);

   fprintf(outputFile, "\n// Solver parameters\n");
   PRINT(numSvals, %d);
   PRINT(aNorm, %e);
   PRINT(eps, %e);
   PRINT(maxBasisSize, %d);
   PRINT(maxBlockSize, %d);
   PRINT_PRIMME_INT(maxMatvecs);

   PRINTIF(target, primme_svds_smallest);
   PRINTIF(target, primme_svds_largest);
   PRINTIF(target, primme_svds_closest_abs);

   PRINT(numTargetShifts, %d);
   if (primme_svds.numTargetShifts > 0) {
      fprintf(outputFile, "primme_svds.targetShifts =");
      for (i=0; i<primme_svds.numTargetShifts;i++) {
         fprintf(outputFile, " %e",primme_svds.targetShifts[i]);
      }
      fprintf(outputFile, "\n");
   }

   PRINT(locking, %d);
   PRINT(initSize, %d);
   PRINT(numOrthoConst, %d);
   fprintf(outputFile, "primme_svds.iseed =");
   for (i=0; i<4;i++) {
      fprintf(outputFile, " %" PRIMME_INT_P, primme_svds.iseed[i]);
   }
   fprintf(outputFile, "\n");

   PRINT(precondition, %d);

   PRINTIF(method, primme_svds_op_none);
   PRINTIF(method, primme_svds_op_AtA);
   PRINTIF(method, primme_svds_op_AAt);
   PRINTIF(method, primme_svds_op_augmented);

   PRINTIF(methodStage2, primme_svds_op_none);
   PRINTIF(methodStage2, primme_svds_op_AtA);
   PRINTIF(methodStage2, primme_svds_op_AAt);
   PRINTIF(methodStage2, primme_svds_op_augmented);

   if (primme_svds.method != primme_svds_op_none) {
      fprintf(outputFile, "\n"
                          "// ---------------------------------------------------\n"
                          "//            1st stage primme configuration          \n"
                          "// ---------------------------------------------------\n");
      primme_svds.primme.outputFile = outputFile;
      primme_display_params_prefix("primme", primme_svds.primme);
   }

   if (primme_svds.methodStage2 != primme_svds_op_none) {
      fprintf(outputFile, "\n"
                          "// ---------------------------------------------------\n"
                          "//            2st stage primme configuration          \n"
                          "// ---------------------------------------------------\n");
      primme_svds.primmeStage2.outputFile = outputFile;
      primme_display_params_prefix("primmeStage2", primme_svds.primmeStage2);
   }
   fflush(outputFile);

}


/*******************************************************************************
 * Subroutine primme_svds_free - Free memory resources allocated by Sprimme_svds
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------------------------
 * primme_svds  Structure containing various solver parameters and statistics
 *
 ******************************************************************************/

void primme_svds_free(primme_svds_params *primme) {
    
   free(primme->intWork);
   free(primme->realWork);
   primme->intWorkSize  = 0;
   primme->realWorkSize = 0;
}

/*******************************************************************************
 * Subroutine globalSumRealSvds - implementation of primme_params' globalSumReal
 *    that uses the callback defined in primme_svds_params.
 * 
 ******************************************************************************/

static void globalSumRealSvds(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme, int *ierr) {
   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   primme_svds->globalSumReal(sendBuf, recvBuf, count, primme_svds, ierr);
}

#endif /* USE_DOUBLE */
