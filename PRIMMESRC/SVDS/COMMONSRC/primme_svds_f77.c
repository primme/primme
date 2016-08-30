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
 * File: primme_svds_f77.c
 *
 * Purpose - Implementation of PRIMME_SVDS F77 interface functions.
 *
 ******************************************************************************/

#include <stdlib.h>   /* free */
#include "primme_svds.h"
#include "primme_svds_f77_private.h"


void AS_FORTRAN(dprimme_svds)(double *svals, double *svecs, double *resNorms,
      primme_svds_params **primme_svds, int *ierr) {
   *ierr = dprimme_svds(svals, svecs, resNorms, *primme_svds);
}

void AS_FORTRAN(zprimme_svds)(double *svals, double complex *svecs, double *resNorms,
      primme_svds_params **primme_svds, int *ierr) {
   *ierr = zprimme_svds(svals, svecs, resNorms, *primme_svds);
}

void AS_FORTRAN(primme_svds_initialize)(primme_svds_params **primme_svds) {
   *primme_svds = (primme_svds_params*)primme_calloc(1,
         sizeof(primme_svds_params), "primme_svds");
   primme_svds_initialize(*primme_svds);
}

void AS_FORTRAN(primme_svds_set_method)(primme_svds_preset_method *method,
      primme_preset_method *methodStage1, primme_preset_method *methodStage2,
      primme_svds_params **primme_svds, int *ierr) {
   *ierr = primme_svds_set_method(*method, *methodStage1, *methodStage2,
         *primme_svds);
}

void AS_FORTRAN(primme_svds_display_params)(primme_svds_params **primme_svds) {
   primme_svds_display_params(**primme_svds);
}

void AS_FORTRAN(primme_svds_free)(primme_svds_params **primme_svds) {
   primme_svds_Free(*primme_svds);
   free(*primme_svds);
   *primme_svds = NULL;
}

void AS_FORTRAN(primme_svds_set_member)(primme_svds_params **primme_svds_, int *label,
      union f77_value v, int *ierr) {
   int i;
   primme_svds_params *primme_svds = *primme_svds_;
   *ierr = 0;

   switch(*label) {
      case PRIMMEF77_SVDS_primme :
         *ierr = 1;
         break;
      case PRIMMEF77_SVDS_primmeStage2 :
         *ierr = 1;
         break;
      case PRIMMEF77_SVDS_m :
         primme_svds->m = *v.int_v;
         break;
      case PRIMMEF77_SVDS_n :
         primme_svds->n = *v.int_v;
         break;
      case PRIMMEF77_SVDS_matrixMatvec :
         primme_svds->matrixMatvec = v.matFunc_v;
         break;
      case PRIMMEF77_SVDS_applyPreconditioner :
         primme_svds->applyPreconditioner = v.matFunc_v;
         break;
      case PRIMMEF77_SVDS_numProcs :
         primme_svds->numProcs = *v.int_v;
         break;
      case PRIMMEF77_SVDS_procID :
         primme_svds->procID = *v.int_v;
         break;
      case PRIMMEF77_SVDS_mLocal :
         primme_svds->mLocal = *v.int_v;
         break;
      case PRIMMEF77_SVDS_nLocal :
         primme_svds->nLocal = *v.int_v;
         break;
      case PRIMMEF77_SVDS_commInfo :
         primme_svds->commInfo = v.ptr_v;
         break;
      case PRIMMEF77_SVDS_globalSumDouble :
         primme_svds->globalSumDouble = v.globalSumDoubleFunc_v;
         break;
      case PRIMMEF77_SVDS_numSvals :
         primme_svds->numSvals = *v.int_v;
         break;
      case PRIMMEF77_SVDS_target :
         primme_svds->target = *v.target_v;
         break;
      case PRIMMEF77_SVDS_numTargetShifts :
         primme_svds->numTargetShifts = *v.int_v;
         break;
      case PRIMMEF77_SVDS_targetShifts :
         primme_svds->targetShifts = v.double_v;
         break;
      case PRIMMEF77_SVDS_method :
         primme_svds->method = *v.operator_v;
         break;
      case PRIMMEF77_SVDS_methodStage2 :
         primme_svds->methodStage2 = *v.operator_v;
         break;
      case PRIMMEF77_SVDS_intWorkSize :
         primme_svds->intWorkSize = *v.int_v;
         break;
      case PRIMMEF77_SVDS_realWorkSize :
         primme_svds->realWorkSize = *v.long_int_v;
         break;
      case PRIMMEF77_SVDS_intWork :
         primme_svds->intWork = v.int_v;
         break;
      case PRIMMEF77_SVDS_realWork :
         primme_svds->realWork = v.ptr_v;
         break;
      case PRIMMEF77_SVDS_matrix :
         primme_svds->matrix = v.ptr_v;
         break;
      case PRIMMEF77_SVDS_preconditioner :
         primme_svds->preconditioner = v.ptr_v;
         break;
      case PRIMMEF77_SVDS_locking :
         primme_svds->locking = *v.int_v;
         break;
      case PRIMMEF77_SVDS_numOrthoConst :
         primme_svds->numOrthoConst = *v.int_v;
         break;
      case PRIMMEF77_SVDS_aNorm :
         primme_svds->aNorm = *v.double_v;
         break;
      case PRIMMEF77_SVDS_eps :
         primme_svds->eps = *v.double_v;
         break;
      case PRIMMEF77_SVDS_precondition :
         primme_svds->precondition = *v.int_v;
         break;
      case PRIMMEF77_SVDS_initSize :
         primme_svds->initSize = *v.int_v;
         break;
      case PRIMMEF77_SVDS_maxBasisSize :
         primme_svds->maxBasisSize = *v.int_v;
         break;
      case PRIMMEF77_SVDS_maxBlockSize :
         primme_svds->maxBlockSize = *v.int_v;
         break;
      case PRIMMEF77_SVDS_maxMatvecs :
         primme_svds->maxMatvecs = *v.int_v;
         break;
      case PRIMMEF77_SVDS_iseed :
         for (i=0; i<4; i++) {
            primme_svds->iseed[i] = v.int_v[i];
         }
         break;
      case PRIMMEF77_SVDS_printLevel :
         primme_svds->printLevel = *v.int_v;
         break;
      case PRIMMEF77_SVDS_outputFile :
         primme_svds->outputFile = v.file_v;
         break;
      case PRIMMEF77_SVDS_stats_numOuterIterations :
         primme_svds->stats.numOuterIterations = *v.int_v;
         break;
      case PRIMMEF77_SVDS_stats_numRestarts :
         primme_svds->stats.numRestarts = *v.int_v;
         break;
      case PRIMMEF77_SVDS_stats_numMatvecs :
         primme_svds->stats.numMatvecs = *v.int_v;
         break;
      case PRIMMEF77_SVDS_stats_numPreconds :
         primme_svds->stats.numPreconds = *v.int_v;
         break;
      case PRIMMEF77_SVDS_stats_elapsedTime :
         primme_svds->stats.elapsedTime = *v.double_v;
         break;
      default:
         *ierr = 1;
   }
}

void AS_FORTRAN(primme_svdstop_get_member)(primme_svds_params **primme_svds_, int *label,
      union f77_value_ptr *v, int *ierr) {
   int i;
   primme_svds_params *primme_svds = *primme_svds_;
   *ierr = 0;

   switch(*label) {
      case PRIMMEF77_SVDS_primme :
         v->ptr_v = &primme_svds->primme;
         break;
      case PRIMMEF77_SVDS_primmeStage2 :
         v->ptr_v = &primme_svds->primmeStage2;
         break;
      case PRIMMEF77_SVDS_m :
         v->int_v = primme_svds->m;
         break;
      case PRIMMEF77_SVDS_n :
         v->int_v = primme_svds->n;
         break;
      case PRIMMEF77_SVDS_matrixMatvec :
         v->matFunc_v = primme_svds->matrixMatvec;
         break;
      case PRIMMEF77_SVDS_applyPreconditioner :
         v->matFunc_v = primme_svds->applyPreconditioner;
         break;
      case PRIMMEF77_SVDS_numProcs :
         v->int_v = primme_svds->numProcs;
         break;
      case PRIMMEF77_SVDS_procID :
         v->int_v = primme_svds->procID;
         break;
      case PRIMMEF77_SVDS_mLocal :
         v->int_v = primme_svds->mLocal;
         break;
      case PRIMMEF77_SVDS_nLocal :
         v->int_v = primme_svds->nLocal;
         break;
      case PRIMMEF77_SVDS_commInfo :
         v->ptr_v = primme_svds->commInfo;
         break;
      case PRIMMEF77_SVDS_globalSumDouble :
         v->globalSumDoubleFunc_v = primme_svds->globalSumDouble;
         break;
      case PRIMMEF77_SVDS_numSvals :
         v->int_v = primme_svds->numSvals;
         break;
      case PRIMMEF77_SVDS_target :
         v->target_v = primme_svds->target;
         break;
      case PRIMMEF77_SVDS_numTargetShifts :
         v->int_v = primme_svds->numTargetShifts;
         break;
      case PRIMMEF77_SVDS_targetShifts :
         for (i=0; i< primme_svds->numTargetShifts; i++) {
             (&v->double_v)[i] = primme_svds->targetShifts[i];
         }
         break;
      case PRIMMEF77_SVDS_method :
         v->operator_v = primme_svds->method;
         break;
      case PRIMMEF77_SVDS_methodStage2 :
         v->operator_v = primme_svds->methodStage2;
         break;
      case PRIMMEF77_SVDS_intWorkSize :
         v->int_v = primme_svds->intWorkSize;
         break;
      case PRIMMEF77_SVDS_realWorkSize :
         v->long_int_v = primme_svds->realWorkSize;
         break;
      case PRIMMEF77_SVDS_intWork :
         v->ptr_v = primme_svds->intWork;
         break;
      case PRIMMEF77_SVDS_realWork :
         v->ptr_v = primme_svds->realWork;
         break;
      case PRIMMEF77_SVDS_matrix :
         v->ptr_v = primme_svds->matrix;
         break;
      case PRIMMEF77_SVDS_preconditioner :
         v->ptr_v = primme_svds->preconditioner;
         break;
      case PRIMMEF77_SVDS_locking :
         v->int_v = primme_svds->locking;
         break;
      case PRIMMEF77_SVDS_numOrthoConst :
         v->int_v = primme_svds->numOrthoConst;
         break;
      case PRIMMEF77_SVDS_aNorm :
         v->double_v = primme_svds->aNorm;
         break;
      case PRIMMEF77_SVDS_eps :
         v->double_v = primme_svds->eps;
         break;
      case PRIMMEF77_SVDS_precondition :
         v->int_v = primme_svds->precondition;
         break;
      case PRIMMEF77_SVDS_initSize :
         v->int_v = primme_svds->initSize;
         break;
      case PRIMMEF77_SVDS_maxBasisSize :
         v->int_v = primme_svds->maxBasisSize;
         break;
      case PRIMMEF77_SVDS_maxBlockSize :
         v->int_v = primme_svds->maxBlockSize;
         break;
      case PRIMMEF77_SVDS_maxMatvecs :
         v->int_v = primme_svds->maxMatvecs;
         break;
      case PRIMMEF77_SVDS_iseed :
         for (i=0; i<4; i++) {
            (&v->int_v)[i] = primme_svds->iseed[i];
         }
         break;
      case PRIMMEF77_SVDS_printLevel :
         v->int_v = primme_svds->printLevel;
         break;
      case PRIMMEF77_SVDS_outputFile :
         v->file_v = primme_svds->outputFile;
         break;
      case PRIMMEF77_SVDS_stats_numOuterIterations :
         v->int_v = primme_svds->stats.numOuterIterations;
         break;
      case PRIMMEF77_SVDS_stats_numRestarts :
         v->int_v = primme_svds->stats.numRestarts;
         break;
      case PRIMMEF77_SVDS_stats_numMatvecs :
         v->int_v = primme_svds->stats.numMatvecs;
         break;
      case PRIMMEF77_SVDS_stats_numPreconds :
         v->int_v = primme_svds->stats.numPreconds;
         break;
      case PRIMMEF77_SVDS_stats_elapsedTime :
         v->double_v = primme_svds->stats.elapsedTime;
         break;
      default:
         *ierr = 1;
   }
}

void AS_FORTRAN(primme_svds_get_member)(primme_svds_params *primme_svds, int *label,
      union f77_value_ptr *v, int *ierr) {
   AS_FORTRAN(primme_svdstop_get_member)(&primme_svds, label, v, ierr);
}
