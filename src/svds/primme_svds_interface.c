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
 * File: primme_svds_interface.c
 *
 * Purpose - Contains interface functions to PRIMME SVDS named primme_svds_*.
 *
 ******************************************************************************/

#include <stdlib.h>   /* mallocs, free */
#include <stdio.h>    
#include <math.h>    
#include <string.h>    
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
   primme_svds->stats.numOuterIterations            = 0; 
   primme_svds->stats.numRestarts                   = 0;
   primme_svds->stats.numMatvecs                    = 0;
   primme_svds->stats.numPreconds                   = 0;
   primme_svds->stats.numGlobalSum                  = 0;
   primme_svds->stats.volumeGlobalSum               = 0;
   primme_svds->stats.numOrthoInnerProds            = 0.0;
   primme_svds->stats.elapsedTime                   = 0.0;
   primme_svds->stats.timeMatvec                    = 0.0;
   primme_svds->stats.timePrecond                   = 0.0;
   primme_svds->stats.timeOrtho                     = 0.0;
   primme_svds->stats.timeGlobalSum                 = 0.0;

   /* Internally used variables */
   primme_svds->iseed[0] = -1;   /* To set iseed, we first need procID           */ 
   primme_svds->iseed[1] = -1;   /* Thus we set all iseeds to -1                 */
   primme_svds->iseed[2] = -1;   /* Unless users provide their own iseeds,       */
   primme_svds->iseed[3] = -1;   /* PRIMME will set thse later uniquely per proc */
   primme_svds->intWorkSize             = 0;
   primme_svds->realWorkSize            = 0;
   primme_svds->intWork                 = NULL;
   primme_svds->realWork                = NULL;
   primme_svds->monitorFun              = NULL;
   primme_svds->monitor                 = NULL;

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

/*******************************************************************************
 * Subroutine primme_svds_get_member - get the value of a parameter in
 *    primme_svds_params
 * 
 * INPUT PARAMETERS
 * ----------------
 * primme_svds  Structure containing various solver parameters and statistics
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

int primme_svds_get_member(primme_svds_params *primme_svds,
      primme_svds_params_label label, void *value) {

   int i;
   union value_t {
      PRIMME_INT int_v;
      void (*matFunc_v) (void*,PRIMME_INT*,void*,PRIMME_INT*,int*,int*,struct primme_svds_params*,int*);
      void *ptr_v;
      void (*globalSumRealFunc_v) (void *,void *,int *,struct primme_svds_params*,int*);
      primme_svds_target target_v;
      primme_svds_operator operator_v;
      double double_v;
      FILE *file_v;
      void (*monitorFun_v)(void *basisSvals, int *basisSize, int *basisFlags,
            int *iblock, int *blockSize, void *basisNorms, int *numConverged,
            void *lockedSvals, int *numLocked, int *lockedFlags, void *lockedNorms,
            int *inner_its, void *LSRes, primme_event *event, int *stage,
            struct primme_svds_params *primme_svds, int *err);
   } *v = (union value_t*)value;

   switch(label) {
      case PRIMME_SVDS_primme :
         v->ptr_v = &primme_svds->primme;
         break;
      case PRIMME_SVDS_primmeStage2 :
         v->ptr_v = &primme_svds->primmeStage2;
         break;
      case PRIMME_SVDS_m :
         v->int_v = primme_svds->m;
         break;
      case PRIMME_SVDS_n :
         v->int_v = primme_svds->n;
         break;
      case PRIMME_SVDS_matrixMatvec :
         v->matFunc_v = primme_svds->matrixMatvec;
         break;
      case PRIMME_SVDS_applyPreconditioner :
         v->matFunc_v = primme_svds->applyPreconditioner;
         break;
      case PRIMME_SVDS_numProcs :
         v->int_v = primme_svds->numProcs;
         break;
      case PRIMME_SVDS_procID :
         v->int_v = primme_svds->procID;
         break;
      case PRIMME_SVDS_mLocal :
         v->int_v = primme_svds->mLocal;
         break;
      case PRIMME_SVDS_nLocal :
         v->int_v = primme_svds->nLocal;
         break;
      case PRIMME_SVDS_commInfo :
         v->ptr_v = primme_svds->commInfo;
         break;
      case PRIMME_SVDS_globalSumReal :
         v->globalSumRealFunc_v = primme_svds->globalSumReal;
         break;
      case PRIMME_SVDS_numSvals :
         v->int_v = primme_svds->numSvals;
         break;
      case PRIMME_SVDS_target :
         v->target_v = primme_svds->target;
         break;
      case PRIMME_SVDS_numTargetShifts :
         v->int_v = primme_svds->numTargetShifts;
         break;
      case PRIMME_SVDS_targetShifts :
         for (i=0; i< primme_svds->numTargetShifts; i++) {
             (&v->double_v)[i] = primme_svds->targetShifts[i];
         }
         break;
      case PRIMME_SVDS_method :
         v->operator_v = primme_svds->method;
         break;
      case PRIMME_SVDS_methodStage2 :
         v->operator_v = primme_svds->methodStage2;
         break;
      case PRIMME_SVDS_intWorkSize :
         v->int_v = primme_svds->intWorkSize;
         break;
      case PRIMME_SVDS_realWorkSize :
         v->int_v = (PRIMME_INT)primme_svds->realWorkSize;
         break;
      case PRIMME_SVDS_intWork :
         v->ptr_v = primme_svds->intWork;
         break;
      case PRIMME_SVDS_realWork :
         v->ptr_v = primme_svds->realWork;
         break;
      case PRIMME_SVDS_matrix :
         v->ptr_v = primme_svds->matrix;
         break;
      case PRIMME_SVDS_preconditioner :
         v->ptr_v = primme_svds->preconditioner;
         break;
      case PRIMME_SVDS_locking :
         v->int_v = primme_svds->locking;
         break;
      case PRIMME_SVDS_numOrthoConst :
         v->int_v = primme_svds->numOrthoConst;
         break;
      case PRIMME_SVDS_aNorm :
         v->double_v = primme_svds->aNorm;
         break;
      case PRIMME_SVDS_eps :
         v->double_v = primme_svds->eps;
         break;
      case PRIMME_SVDS_precondition :
         v->int_v = primme_svds->precondition;
         break;
      case PRIMME_SVDS_initSize :
         v->int_v = primme_svds->initSize;
         break;
      case PRIMME_SVDS_maxBasisSize :
         v->int_v = primme_svds->maxBasisSize;
         break;
      case PRIMME_SVDS_maxBlockSize :
         v->int_v = primme_svds->maxBlockSize;
         break;
      case PRIMME_SVDS_maxMatvecs :
         v->int_v = primme_svds->maxMatvecs;
         break;
      case PRIMME_SVDS_iseed :
         for (i=0; i<4; i++) {
            (&v->int_v)[i] = primme_svds->iseed[i];
         }
         break;
      case PRIMME_SVDS_printLevel :
         v->int_v = primme_svds->printLevel;
         break;
      case PRIMME_SVDS_outputFile :
         v->file_v = primme_svds->outputFile;
         break;
      case PRIMME_SVDS_stats_numOuterIterations :
         v->int_v = primme_svds->stats.numOuterIterations;
         break;
      case PRIMME_SVDS_stats_numRestarts :
         v->int_v = primme_svds->stats.numRestarts;
         break;
      case PRIMME_SVDS_stats_numMatvecs :
         v->int_v = primme_svds->stats.numMatvecs;
         break;
      case PRIMME_SVDS_stats_numPreconds :
         v->int_v = primme_svds->stats.numPreconds;
         break;
      case PRIMME_SVDS_stats_numGlobalSum:
         v->int_v = primme_svds->stats.numGlobalSum;
         break;
      case PRIMME_SVDS_stats_volumeGlobalSum:
         v->int_v = primme_svds->stats.volumeGlobalSum;
         break;
      case PRIMME_SVDS_stats_numOrthoInnerProds:
         v->double_v = primme_svds->stats.numOrthoInnerProds;
         break;
      case PRIMME_SVDS_stats_elapsedTime :
         v->double_v = primme_svds->stats.elapsedTime;
         break;
      case PRIMME_SVDS_stats_timeMatvec:
         v->double_v = primme_svds->stats.timeMatvec;
         break;
      case PRIMME_SVDS_stats_timePrecond:
         v->double_v = primme_svds->stats.timePrecond;
         break;
      case PRIMME_SVDS_stats_timeOrtho:
         v->double_v = primme_svds->stats.timeOrtho;
         break;
      case PRIMME_SVDS_stats_timeGlobalSum:
         v->double_v = primme_svds->stats.timeGlobalSum;
         break;
      case PRIMME_SVDS_monitorFun:
         v->monitorFun_v = primme_svds->monitorFun;
         break;
      case PRIMME_SVDS_monitor:
         v->ptr_v = primme_svds->monitor;
         break;
      default:
         return 1;
   }
   return 0;
}

/*******************************************************************************
 * Subroutine primme_svds_set_member - set the value to a parameter in
 *    primme_svds_params
 * 
 * INPUT PARAMETERS
 * ----------------
 * label   reference to the parameter
 * value   value of the parameter
 *
 * INPUT/OUTPUT PARAMETERS
 * -----------------
 * primme_svds  Structure containing various solver parameters and statistics
 *
 * RETURN
 * ------
 * error code  zero if ok
 *
 ******************************************************************************/

int primme_svds_set_member(primme_svds_params *primme_svds,
      primme_svds_params_label label, void *value) {
   int i;
   union value_t {
      PRIMME_INT *int_v;
      void (*matFunc_v) (void*,PRIMME_INT*,void*,PRIMME_INT*,int*,int*,struct primme_svds_params*,int*);
      void *ptr_v;
      void (*globalSumRealFunc_v) (void *,void *,int *,struct primme_svds_params*,int*);
      primme_svds_target *target_v;
      primme_svds_operator *operator_v;
      double *double_v;
      FILE *file_v;
      void (*monitorFun_v)(void *basisSvals, int *basisSize, int *basisFlags,
            int *iblock, int *blockSize, void *basisNorms, int *numConverged,
            void *lockedSvals, int *numLocked, int *lockedFlags, void *lockedNorms,
            int *inner_its, void *LSRes, primme_event *event, int *stage,
            struct primme_svds_params *primme_svds, int *err);

   } v = *(union value_t*)&value;

   switch(label) {
      case PRIMME_SVDS_primme :
         return 1;
         break;
      case PRIMME_SVDS_primmeStage2 :
         return 1;
         break;
      case PRIMME_SVDS_m :
         primme_svds->m = *v.int_v;
         break;
      case PRIMME_SVDS_n :
         primme_svds->n = *v.int_v;
         break;
      case PRIMME_SVDS_matrixMatvec :
         primme_svds->matrixMatvec = v.matFunc_v;
         break;
      case PRIMME_SVDS_applyPreconditioner :
         primme_svds->applyPreconditioner = v.matFunc_v;
         break;
      case PRIMME_SVDS_numProcs :
         if (*v.int_v > INT_MAX) return 1; else
         primme_svds->numProcs = (int)*v.int_v;
         break;
      case PRIMME_SVDS_procID :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->procID = (int)*v.int_v;
         break;
      case PRIMME_SVDS_mLocal :
         primme_svds->mLocal = *v.int_v;
         break;
      case PRIMME_SVDS_nLocal :
         primme_svds->nLocal = *v.int_v;
         break;
      case PRIMME_SVDS_commInfo :
         primme_svds->commInfo = v.ptr_v;
         break;
      case PRIMME_SVDS_globalSumReal :
         primme_svds->globalSumReal = v.globalSumRealFunc_v;
         break;
      case PRIMME_SVDS_numSvals :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->numSvals = (int)*v.int_v;
         break;
      case PRIMME_SVDS_target :
         primme_svds->target = *v.target_v;
         break;
      case PRIMME_SVDS_numTargetShifts :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->numTargetShifts = (int)*v.int_v;
         break;
      case PRIMME_SVDS_targetShifts :
         primme_svds->targetShifts = v.double_v;
         break;
      case PRIMME_SVDS_method :
         primme_svds->method = *v.operator_v;
         break;
      case PRIMME_SVDS_methodStage2 :
         primme_svds->methodStage2 = *v.operator_v;
         break;
      case PRIMME_SVDS_intWorkSize :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->intWorkSize = (int)*v.int_v;
         break;
      case PRIMME_SVDS_realWorkSize :
         primme_svds->realWorkSize = (size_t)*v.int_v;
         break;
      case PRIMME_SVDS_intWork :
         primme_svds->intWork = (int*)v.int_v;
         break;
      case PRIMME_SVDS_realWork :
         primme_svds->realWork = v.ptr_v;
         break;
      case PRIMME_SVDS_matrix :
         primme_svds->matrix = v.ptr_v;
         break;
      case PRIMME_SVDS_preconditioner :
         primme_svds->preconditioner = v.ptr_v;
         break;
      case PRIMME_SVDS_locking :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->locking = (int)*v.int_v;
         break;
      case PRIMME_SVDS_numOrthoConst :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->numOrthoConst = (int)*v.int_v;
         break;
      case PRIMME_SVDS_aNorm :
         primme_svds->aNorm = *v.double_v;
         break;
      case PRIMME_SVDS_eps :
         primme_svds->eps = *v.double_v;
         break;
      case PRIMME_SVDS_precondition :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->precondition = (int)*v.int_v;
         break;
      case PRIMME_SVDS_initSize :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->initSize = (int)*v.int_v;
         break;
      case PRIMME_SVDS_maxBasisSize :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->maxBasisSize = (int)*v.int_v;
         break;
      case PRIMME_SVDS_maxBlockSize :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->maxBlockSize = (int)*v.int_v;
         break;
      case PRIMME_SVDS_maxMatvecs :
         primme_svds->maxMatvecs = *v.int_v;
         break;
      case PRIMME_SVDS_iseed :
         for (i=0; i<4; i++) {
            primme_svds->iseed[i] = v.int_v[i];
         }
         break;
      case PRIMME_SVDS_printLevel :
         if (*v.int_v > INT_MAX) return 1; else 
         primme_svds->printLevel = (int)*v.int_v;
         break;
      case PRIMME_SVDS_outputFile :
         primme_svds->outputFile = v.file_v;
         break;
      case PRIMME_SVDS_stats_numOuterIterations :
         primme_svds->stats.numOuterIterations = *v.int_v;
         break;
      case PRIMME_SVDS_stats_numRestarts :
         primme_svds->stats.numRestarts = *v.int_v;
         break;
      case PRIMME_SVDS_stats_numMatvecs :
         primme_svds->stats.numMatvecs = *v.int_v;
         break;
      case PRIMME_SVDS_stats_numPreconds :
         primme_svds->stats.numPreconds = *v.int_v;
         break;
      case PRIMME_SVDS_stats_volumeGlobalSum:
         primme_svds->stats.volumeGlobalSum = *v.int_v;
         break;
      case PRIMME_SVDS_stats_numOrthoInnerProds:
         primme_svds->stats.numOrthoInnerProds = *v.double_v;
         break;
      case PRIMME_SVDS_stats_elapsedTime :
         primme_svds->stats.elapsedTime = *v.double_v;
         break;
      case PRIMME_SVDS_stats_timeMatvec:
         primme_svds->stats.timeMatvec = *v.double_v;
         break;
      case PRIMME_SVDS_stats_timePrecond:
         primme_svds->stats.timePrecond = *v.double_v;
         break;
      case PRIMME_SVDS_stats_timeOrtho:
         primme_svds->stats.timeOrtho = *v.double_v;
         break;
      case PRIMME_SVDS_stats_timeGlobalSum:
         primme_svds->stats.timeGlobalSum = *v.double_v;
         break;
      case PRIMME_SVDS_monitorFun:
         primme_svds->monitorFun = v.monitorFun_v;
         break;
      case PRIMME_SVDS_monitor:
         primme_svds->monitor = v.ptr_v;
      break;
      default:
         return 1;
   }
   return 0;
}

/*******************************************************************************
 * Subroutine primme_svds_member_info - return the label value or the label
 *    name, the type and the arity of some member of primme_svds_params.
 * 
 * INPUT/OUTPUT PARAMETERS
 * ----------------
 * label       reference to the parameter by the value in primme_svds_params_label
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

int primme_svds_member_info(primme_svds_params_label *label_,
      const char** label_name_, primme_type *type, int *arity) {
   primme_svds_params_label label = (primme_svds_params_label)1000;
   const char *label_name;

   /* Quick exit when neither label nor label_name is given */

   if (label_ == NULL && (label_name_ == NULL || *label_name_ == NULL)) {
      return 1;
   }

   /* Get the label from label_name_ and label_name from label_ */

#define IF_IS(F) \
   if ((label_name_ && *label_name_ && strcmp(#F, *label_name_) == 0) \
         || (label_ && *label_ == PRIMME_SVDS_ ## F)) { \
      label = PRIMME_SVDS_ ## F; \
      label_name = #F; \
   }

   IF_IS(primme);
   IF_IS(primmeStage2);
   IF_IS(m);
   IF_IS(n);
   IF_IS(matrixMatvec);
   IF_IS(applyPreconditioner);
   IF_IS(numProcs);
   IF_IS(procID);
   IF_IS(mLocal);
   IF_IS(nLocal);
   IF_IS(commInfo);
   IF_IS(globalSumReal);
   IF_IS(numSvals);
   IF_IS(target);
   IF_IS(numTargetShifts);
   IF_IS(targetShifts);
   IF_IS(method);
   IF_IS(methodStage2);
   IF_IS(intWorkSize);
   IF_IS(realWorkSize);
   IF_IS(intWork);
   IF_IS(realWork);
   IF_IS(matrix);
   IF_IS(preconditioner);
   IF_IS(locking);
   IF_IS(numOrthoConst);
   IF_IS(aNorm);
   IF_IS(eps);
   IF_IS(precondition);
   IF_IS(initSize);
   IF_IS(maxBasisSize);
   IF_IS(maxBlockSize);
   IF_IS(maxMatvecs);
   IF_IS(iseed);
   IF_IS(printLevel);
   IF_IS(outputFile);
   IF_IS(stats_numOuterIterations);
   IF_IS(stats_numRestarts);
   IF_IS(stats_numMatvecs);
   IF_IS(stats_numPreconds);
   IF_IS(stats_numGlobalSum);
   IF_IS(stats_volumeGlobalSum);
   IF_IS(stats_numOrthoInnerProds);
   IF_IS(stats_elapsedTime);
   IF_IS(stats_timeMatvec);
   IF_IS(stats_timePrecond);
   IF_IS(stats_timeOrtho);
   IF_IS(stats_timeGlobalSum);
   IF_IS(monitorFun);
   IF_IS(monitor);
#undef IF_IS

   /* Return label/label_name */

   if (label_) *label_ = label;
   if (label_name_) *label_name_ = label_name;

   /* Return type and arity */

   switch(label) {
      /* members with type int */

      case PRIMME_SVDS_m: 
      case PRIMME_SVDS_n:
      case PRIMME_SVDS_numSvals:
      case PRIMME_SVDS_target:
      case PRIMME_SVDS_method:
      case PRIMME_SVDS_methodStage2:
      case PRIMME_SVDS_locking:
      case PRIMME_SVDS_numOrthoConst:
      case PRIMME_SVDS_precondition:
      case PRIMME_SVDS_initSize:
      case PRIMME_SVDS_maxBasisSize:
      case PRIMME_SVDS_maxBlockSize:
      case PRIMME_SVDS_maxMatvecs:
      case PRIMME_SVDS_printLevel:
      case PRIMME_SVDS_stats_numOuterIterations:
      case PRIMME_SVDS_stats_numRestarts:
      case PRIMME_SVDS_stats_numMatvecs:
      case PRIMME_SVDS_stats_numPreconds:
      case PRIMME_SVDS_stats_numGlobalSum:
      case PRIMME_SVDS_stats_volumeGlobalSum:
      case PRIMME_SVDS_iseed:
      case PRIMME_SVDS_numProcs: 
      case PRIMME_SVDS_procID: 
      case PRIMME_SVDS_mLocal: 
      case PRIMME_SVDS_nLocal: 
      case PRIMME_SVDS_numTargetShifts:
      case PRIMME_SVDS_intWorkSize:
      case PRIMME_SVDS_realWorkSize:
      if (type) *type = primme_int;
      if (arity) *arity = 1;
      break;

      /* members with type double */

      case PRIMME_SVDS_aNorm:
      case PRIMME_SVDS_eps:
      case PRIMME_SVDS_stats_elapsedTime:
      case PRIMME_SVDS_stats_numOrthoInnerProds:
      case PRIMME_SVDS_stats_timeMatvec:
      case PRIMME_SVDS_stats_timePrecond:
      case PRIMME_SVDS_stats_timeOrtho:
      case PRIMME_SVDS_stats_timeGlobalSum:
      if (type) *type = primme_double;
      if (arity) *arity = 1;
      break;
  
      case PRIMME_SVDS_targetShifts:
      if (type) *type = primme_double;
      if (arity) *arity = 0;
      break;

      /* members with type pointer */

      case PRIMME_SVDS_primme:
      case PRIMME_SVDS_primmeStage2:
      case PRIMME_SVDS_matrixMatvec: 
      case PRIMME_SVDS_applyPreconditioner:
      case PRIMME_SVDS_commInfo:
      case PRIMME_SVDS_globalSumReal:
      case PRIMME_SVDS_intWork:
      case PRIMME_SVDS_realWork:
      case PRIMME_SVDS_matrix:
      case PRIMME_SVDS_preconditioner:
      case PRIMME_SVDS_outputFile:
      case PRIMME_SVDS_monitorFun:
      case PRIMME_SVDS_monitor:
      if (type) *type = primme_pointer;
      if (arity) *arity = 1;
      break;

      default: 
      return 1;
   }

   return 0;
}
 
/*******************************************************************************
 * Subroutine primme_svds_constant_info - return the value of a primme svds enum
 *    constant.
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

int primme_svds_constant_info(const char* label_name, int *value) {

#define IF_IS(F) if (strcmp(#F, label_name) == 0) {*value = (int)F; return 0;}

   /* method in primme_svds_set_method */
   
   IF_IS(primme_svds_default);
   IF_IS(primme_svds_hybrid);
   IF_IS(primme_svds_normalequations);
   IF_IS(primme_svds_augmented);
   
   /* enum members for targeting and operator */
   
   IF_IS(primme_svds_largest);
   IF_IS(primme_svds_smallest);
   IF_IS(primme_svds_closest_abs);
   IF_IS(primme_svds_op_none);
   IF_IS(primme_svds_op_AtA);
   IF_IS(primme_svds_op_AAt);
   IF_IS(primme_svds_op_augmented);
#undef IF_IS

   /* return error if label not found */

   return 1;   
}

#endif /* USE_DOUBLE */
