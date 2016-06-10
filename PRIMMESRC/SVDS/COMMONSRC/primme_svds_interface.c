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
 * File: primme_svds_interface.c
 *
 * Purpose - Contains interface functions to PRIMME SVDS named primme_svds_*.
 *
 ******************************************************************************/

#include <stdlib.h>   /* mallocs, free */
#include <stdio.h>    
#include <math.h>    
#include "primme_svds.h"
#include "primme_svds_interface.h"
#include "primme_interface.h"
#include "common_numerical.h"
#include "primme_svds_interface_private.h"

/***************************************************************************

   Initialize the primme_svds data structure
  
***************************************************************************/
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
   primme_svds->globalSumDouble         = NULL;

   /* Use these pointers to provide matrix/preconditioner */
   primme_svds->matrix                  = NULL;
   primme_svds->preconditioner          = NULL;

   /* Matvec and preconditioner */
   primme_svds->matrixMatvec            = NULL;
   primme_svds->applyPreconditioner     = NULL;

   /* Other important parameters users may set */
   primme_svds->aNorm                   = 0.0L;
   primme_svds->eps                     = 1e-12;
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

int primme_svds_set_method(primme_svds_preset_method method,
      primme_preset_method methodStage1, primme_preset_method methodStage2,
      primme_svds_params *primme_svds) {

   /* Set method and methodStage2 in primme_svds_params */
   switch(method) {
   case primme_svds_default:
   case primme_svds_hybrid:
      primme_svds->method = primme_svds->n <= primme_svds->m ? primme_svds_op_AtA : primme_svds_op_AAt;
      primme_svds->methodStage2 = primme_svds_op_augmented;
      break;
   case primme_svds_normalequations:
      primme_svds->method = primme_svds->n <= primme_svds->m ? primme_svds_op_AtA : primme_svds_op_AAt;
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
   if (methodStage2 == DEFAULT_METHOD) methodStage2 = JDQMR;
   if (primme_svds->methodStage2 != primme_svds_op_none)
      primme_set_method(methodStage2, &primme_svds->primmeStage2);

   return 0;
}

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
      primme_svds_set_method(primme_svds_default, DEFAULT_METHOD,
            DEFAULT_METHOD, primme_svds);
   }

   /* Copy values set in primme_svds to the first stage underneath eigensolver */
   copy_params_from_svds(primme_svds, 0);

   if (primme_svds->methodStage2 != primme_svds_op_none) {
      /* Copy values set in primme_svds to the second stage underneath eigensolver */
      copy_params_from_svds(primme_svds, 1);

      /* NOTE: refined extraction seems to work better than RR */
      if (primme_svds->primmeStage2.projectionParams.projection == primme_proj_default)
         primme_svds->primmeStage2.projectionParams.projection = primme_proj_refined;
   }
}

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
   if (primme_svds->numProcs > 1 && primme_svds->globalSumDouble != NULL) {
      primme->procID = primme_svds->procID;
      primme->numProcs = primme_svds->numProcs;
      primme->commInfo = primme_svds->commInfo;
      primme->globalSumDouble = globalSumDoubleSvds;
   }

   switch(method) {
   case primme_svds_op_AtA:
      primme->n = primme_svds->n;
      primme->nLocal = primme_svds->nLocal;
      primme->convTestFun = convTestFunATA;
      break;
   case primme_svds_op_AAt:
      primme->n = primme_svds->m;
      primme->nLocal = primme_svds->mLocal;
      primme->convTestFun = convTestFunATA;
      break;
   case primme_svds_op_augmented:
      primme->n = primme_svds->m + primme_svds->n;
      primme->nLocal = primme_svds->mLocal + primme_svds->nLocal;
      primme->convTestFun = convTestFunAugmented;
      break;
   case primme_svds_op_none:
      break;
   }

   if (stage == 0) {
      switch (primme_svds->target) {
      case primme_svds_largest:
         primme->target = primme_largest;
         break;
      case primme_svds_smallest:
         primme->target = primme_smallest;
         break;
      case primme_svds_closest_abs:
         primme->target = method == primme_svds_op_augmented ?
            primme_closest_geq : primme_closest_abs;
         primme->numTargetShifts = primme_svds->numTargetShifts;
         break;
      }
   }
   else {
      primme->target = primme_closest_geq;
      primme->initBasisMode = primme_init_user;
   }
      
   if (primme->target == primme_smallest || primme->target == primme_largest){
      if (primme_svds->locking >= 0)
         primme->locking = primme_svds->locking;
   }
   else {
      if (primme_svds->locking >= 0)
         primme->locking = primme_svds->locking;
      else if (primme->locking < 0)
         primme->locking = 1;
      primme->numTargetShifts = primme_svds->numSvals;
   }

   if (primme_svds->precondition >= 0)
      primme->correctionParams.precondition = primme_svds->precondition;
   else if (primme->correctionParams.precondition < 0)
      primme->correctionParams.precondition = primme_svds->applyPreconditioner ? 1 : 0;

}

/******************************************************************************
 *
 * void primme_svds_display_params(primme_params *primme);
 *
 *    Displays the current configuration of primme data structure
 *
 *****************************************************************************/
void primme_svds_display_params(primme_svds_params primme_svds) {

   int i;
   FILE *outputFile = primme_svds.outputFile;

#define PRINT(P,L) fprintf(outputFile, "primme_svds." #P " = " #L "\n", primme_svds. P);
#define PRINTIF(P,V) if (primme_svds. P == V) fprintf(outputFile, "primme_svds." #P " = " #V "\n");
#define PRINTParams(P,S,L) fprintf(outputFile, "primme_svds." #P "." #S " = " #L "\n", \
                                    primme_svds. P ## Params. S);
#define PRINTParamsIF(P,S,V) if (primme_svds. P ## Params. S == V) \
                                 fprintf(outputFile, "primme_svds." #P "." #S " = " #V "\n");

fprintf(outputFile, "// ---------------------------------------------------\n"
                    "//            primme_svds configuration               \n"
                    "// ---------------------------------------------------\n");

   PRINT(m, %d);
   PRINT(n, %d);
   PRINT(mLocal, %d);
   PRINT(nLocal, %d);
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
   PRINT(maxMatvecs, %d);

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
      fprintf(outputFile, " %d",primme_svds.iseed[i]);
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

  /**************************************************************************/
} /* end of display params */


void primme_svds_Free(primme_svds_params *params) {
    
   free(params->intWork);
   free(params->realWork);
   params->intWorkSize  = 0;
   params->realWorkSize = 0;
}

/*******************************************************************************
 * Subroutine convTestFunATA - This routine implements primme_params.
 *    convTestFun and returns an approximate eigenpair converged when           
 *    resNorm < eps * sval * primme_svds.aNorm = eps * sqrt(eval*primme.aNorm)
 *    resNorm is close to machineEpsilon * primme.aNorm.
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

static void convTestFunATA(double *eval, void *evec, double *rNorm, int *isConv,
   primme_params *primme) {

   const double machEps = Num_dlamch_primme("E");
   *isConv = *rNorm < max(
               primme->eps * sqrt(fabs(*eval * (
                  primme->aNorm > 0 ? primme->aNorm : primme->stats.estimateLargestSVal))),
               machEps * 3.16 * primme->stats.estimateLargestSVal);
}

/*******************************************************************************
 * Subroutine convTestFunAugmented - This routine implements primme_params.
 *    convTestFun and returns an approximate eigenpair converged when           
 *    resNorm < eps / sqrt(2) * primme_svds.aNorm = eps / sqrt(2) * primme.aNorm.          
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

static void convTestFunAugmented(double *eval, void *evec, double *rNorm, int *isConv,
   primme_params *primme) {

   const double machEps = Num_dlamch_primme("E");
   *isConv = *rNorm < max(
               primme->eps / sqrt(2.0) * (
                     primme->aNorm > 0.0 ? primme->aNorm : primme->stats.estimateLargestSVal),
               machEps * 3.16 * primme->stats.estimateLargestSVal);
} 

static void globalSumDoubleSvds(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme) {
   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   primme_svds->globalSumDouble(sendBuf, recvBuf, count, primme_svds);
}
