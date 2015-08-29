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
 * File: primme_f77.c
 *
 * Purpose - Implementation of PRIMME F77 interface functions.
 *
 ******************************************************************************/

#include "primme.h"
#include "primme_f77_private.h"
#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
/* The following differ only in the name mangling from C                      */
/******************************************************************************/

/*****************************************************************************
 * Wrapper for calling dprimme from Fortran 77/90. 
 * The only difference from primme: the return value passed as parameter 
 *****************************************************************************/
#ifdef F77UNDERSCORE
void dprimme_f77_(double *evals, double *evecs, double *rnorms, 
                primme_params **primme, int *ierr) {
#else
void dprimme_f77(double *evals, double *evecs, double *rnorms, 
                primme_params **primme, int *ierr) {
#endif

  *ierr = dprimme(evals, evecs, rnorms, *primme);

} /* end of dprimme_f77 wrapper for calling from Fortran */

/*****************************************************************************
 * Wrapper for calling zprimme from Fortran 77/90. 
 * The only difference from primme: the return value passed as parameter 
 *****************************************************************************/
#ifdef F77UNDERSCORE
void zprimme_f77_(double *evals, Complex_Z *evecs, double *rnorms, 
                primme_params **primme, int *ierr) {
#else
void zprimme_f77(double *evals, Complex_Z *evecs, double *rnorms, 
                primme_params **primme, int *ierr) {
#endif

  *ierr = zprimme(evals, evecs, rnorms, *primme);

} /* end of dprimme_f77 wrapper for calling from Fortran */

/*****************************************************************************
 * Initialize handles also the allocation of primme structure 
 *****************************************************************************/
#ifdef F77UNDERSCORE
void primme_initialize_f77_(primme_params **primme){
#else
void primme_initialize_f77(primme_params **primme){
#endif

   *primme = (primme_params *)primme_calloc(1,sizeof(primme_params),"primme");
   primme_initialize(*primme);
}

/*****************************************************************************
 *  * Free the internally allocated work arrays of the primme structure 
 *   *****************************************************************************/
#ifdef F77UNDERSCORE
void primme_free_f77_(primme_params **primme){
#else
void primme_free_f77(primme_params **primme){
#endif

   primme_Free(*primme);
}

/*****************************************************************************
 * Wrapper for displaying the primme parameters 
 *****************************************************************************/
#ifdef F77UNDERSCORE
void primme_display_params_f77_(primme_params **primme) {
#else
void primme_display_params_f77(primme_params **primme) {
#endif

   primme_display_params(*(*primme));
}

#ifdef F77UNDERSCORE
void primme_printstacktrace_f77_(primme_params **primme) {
#else
void primme_printstacktrace_f77(primme_params **primme) {
#endif

   primme_PrintStackTrace(*(*primme));
}

/*************************************************************************
 * set method takes "int method" as input. On return:
 *    returnValue = 0 successful completion 
 *    returnValue < 0 no such method exists. If not defined by user, defaults 
 *    have been set for maxBasisSize, minRestartSize, and maxBlockSize   
 *************************************************************************/

#ifdef F77UNDERSCORE
void primme_set_method_f77_(primme_params **primme, int *method, 
                int *returnValue)
#else
void primme_set_method_f77(primme_params **primme, int *method, 
                int *returnValue)
#endif
{
   int d;

   switch (*method) {
      case PRIMMEF77_DYNAMIC:  
              d = primme_set_method(DYNAMIC, *primme); break;
      case PRIMMEF77_DEFAULT_MIN_TIME:  
              d = primme_set_method(DEFAULT_MIN_TIME, *primme); break;
      case PRIMMEF77_DEFAULT_MIN_MATVECS:  
              d = primme_set_method(DEFAULT_MIN_MATVECS, *primme); break;
      case PRIMMEF77_Arnoldi:  
              d = primme_set_method(Arnoldi, *primme); break;
      case PRIMMEF77_GD:  
              d = primme_set_method(GD, *primme); break;
      case PRIMMEF77_GD_plusK:  
              d = primme_set_method(GD_plusK, *primme); break;
      case PRIMMEF77_GD_Olsen_plusK:  
              d = primme_set_method(GD_Olsen_plusK, *primme); break;
      case PRIMMEF77_JD_Olsen_plusK:  
              d = primme_set_method(JD_Olsen_plusK, *primme); break;
      case PRIMMEF77_RQI:  
              d = primme_set_method(RQI, *primme); break;
      case PRIMMEF77_JDQR:  
              d = primme_set_method(JDQR, *primme); break;
      case PRIMMEF77_JDQMR:  
              d = primme_set_method(JDQMR, *primme); break;
      case PRIMMEF77_JDQMR_ETol: 
              d = primme_set_method(JDQMR_ETol, *primme); break;
      case PRIMMEF77_SUBSPACE_ITERATION: 
              d = primme_set_method(SUBSPACE_ITERATION, *primme); break;
      case PRIMMEF77_LOBPCG_OrthoBasis: 
              d = primme_set_method(LOBPCG_OrthoBasis, *primme); break;
      case PRIMMEF77_LOBPCG_OrthoBasis_Window: 
              d = primme_set_method(LOBPCG_OrthoBasis_Window, *primme); break;
      default : fprintf(stderr," Using user parameter settings.\n");
              d = primme_set_method( (primme_preset_method) -1, *primme); break;
   }
   *returnValue = d;
}

/*************************************************************************
 *  Display the runtime statistics from primme (Fortran does not 
 *  have access to these.
 *************************************************************************/
#ifdef F77UNDERSCORE
void primme_display_stats_f77_(primme_params **primme)
#else
void primme_display_stats_f77(primme_params **primme)
#endif
{ 
        fprintf((*primme)->outputFile, 
                                   "--------------------------------------\n");
        fprintf((*primme)->outputFile, "Number of outer iterations: %d\n",
                                        (*primme)->stats.numOuterIterations);
        fprintf((*primme)->outputFile, "Number of Restarts: %d\n",
                                        (*primme)->stats.numRestarts);
        fprintf((*primme)->outputFile, "Number of Matrix-vector products: %d\n",
                                        (*primme)->stats.numMatvecs);
        if ((*primme)->correctionParams.precondition == 1)
        fprintf((*primme)->outputFile, "Number of Precond operations: %d\n",
                                        (*primme)->stats.numPreconds);
        fprintf((*primme)->outputFile, "Total elapsed wall clock Time: %g\n",
                                        (*primme)->stats.elapsedTime);
        fprintf((*primme)->outputFile, "--------------------------------------\n");
}

/*************************************************************************
 * primme_set_member_f77(primme_params *primme, int *label, void *ptr)
 *
 * Sets any of the members of the primme data structure from a Fortran call.
 *
 * primme is a pointer to the primme data structure.
 *
 * void *ptr generic type pointer is used for input, and is cast to the 
 * required pointer type before it is assigned to a data structure member.
 *
 * The choice of member is through an integer label. The integers describing
 * the member are prefixed with PRIMMEF77_ followed by the actual name of 
 * the member. If the member is a structure, underscores are used for dots,
 *
 * Eg the label int for:    primme->correctionParams.projectors.LeftX
 * is:                    PRIMMEF77_correctionParams_projectors_LeftX
 *
 *************************************************************************/
#ifdef F77UNDERSCORE
void primme_set_member_f77_(primme_params **primme, int *label, union f77_value v){
#else
void primme_set_member_f77(primme_params **primme, int *label, union f77_value v){
#endif
   int i;

   switch (*label) {
      case PRIMMEF77_n:
              (*primme)->n = *v.int_v;
      break;
      case PRIMMEF77_matrixMatvec:
              (*primme)->matrixMatvec = v.matFunc_v;
      break;
      case PRIMMEF77_massMatrixMatvec:
              (*primme)->massMatrixMatvec = v.matFunc_v;
      break;
      case PRIMMEF77_applyPreconditioner:
              (*primme)->applyPreconditioner = v.matFunc_v;
      break;
      case PRIMMEF77_numProcs:
              (*primme)->numProcs = *v.int_v;
      break;
      case PRIMMEF77_procID:
              (*primme)->procID = *v.int_v;
      break;
      case PRIMMEF77_commInfo:
              (*primme)->commInfo = v.ptr_v;
      break;
      case PRIMMEF77_nLocal:
              (*primme)->nLocal = *v.int_v;
      break;
      case PRIMMEF77_globalSumDouble:
              (*primme)->globalSumDouble = v.globalSumDoubleFunc_v;
      break;
      case PRIMMEF77_numEvals:
              (*primme)->numEvals = *v.int_v;
      break;
      case PRIMMEF77_target:
              (*primme)->target = *v.target_v;
      break;
      case PRIMMEF77_numTargetShifts:
              (*primme)->numTargetShifts = *v.int_v;
      break;
      case PRIMMEF77_targetShifts:
              (*primme)->targetShifts = v.double_v;
      break;
      case PRIMMEF77_locking:
              (*primme)->locking = *v.int_v;
      break;
      case PRIMMEF77_initSize:
              (*primme)->initSize = *v.int_v;
      break;
      case PRIMMEF77_numOrthoConst:
              (*primme)->numOrthoConst = *v.int_v;
      break;
      case PRIMMEF77_dynamicMethodSwitch:
              (*primme)->dynamicMethodSwitch = *v.int_v;
      break;
      case PRIMMEF77_maxBasisSize:
              (*primme)->maxBasisSize = *v.int_v;
      break;
      case PRIMMEF77_minRestartSize:
              (*primme)->minRestartSize = *v.int_v;
      break;
      case PRIMMEF77_maxBlockSize:
              (*primme)->maxBlockSize = *v.int_v;
      break;
      case PRIMMEF77_maxMatvecs:
              (*primme)->maxMatvecs = *v.int_v;
      break;
      case PRIMMEF77_maxOuterIterations:
              (*primme)->maxOuterIterations = *v.int_v;
      break;
      case PRIMMEF77_intWorkSize:
              (*primme)->intWorkSize = *v.int_v;
      break;
      case PRIMMEF77_realWorkSize:
              (*primme)->realWorkSize = *v.long_int_v;
      break;
      case PRIMMEF77_iseed:
         for (i=0; i< 4; i++) {
            (*primme)->iseed[i] = v.int_v[i];
         }
      break;
      case PRIMMEF77_intWork:
              (*primme)->intWork = v.int_v;
      break;
      case PRIMMEF77_realWork:
              (*primme)->realWork = v.ptr_v;
      break;
      case PRIMMEF77_aNorm:
              (*primme)->aNorm = *v.double_v;
      break;
      case PRIMMEF77_eps:
              (*primme)->eps = *v.double_v;
      break;
      case PRIMMEF77_printLevel:
              (*primme)->printLevel = *v.int_v;
      break;
      case PRIMMEF77_outputFile:
              (*primme)->outputFile = v.file_v;
      break;
      case PRIMMEF77_matrix:
              (*primme)->matrix = v.ptr_v;
      break;
      case PRIMMEF77_preconditioner:
              (*primme)->preconditioner = v.ptr_v;
      break;
      case PRIMMEF77_restartingParams_scheme:
              (*primme)->restartingParams.scheme = *v.restartscheme_v;
      break;
      case PRIMMEF77_restartingParams_maxPrevRetain:
              (*primme)->restartingParams.maxPrevRetain = *v.int_v;
      break;
      case PRIMMEF77_correctionParams_precondition:
              (*primme)->correctionParams.precondition = *v.int_v;
      break;
      case PRIMMEF77_correctionParams_robustShifts:
              (*primme)->correctionParams.robustShifts = *v.int_v;
      break;
      case PRIMMEF77_correctionParams_maxInnerIterations:
              (*primme)->correctionParams.maxInnerIterations = *v.int_v;
      break;
      case PRIMMEF77_correctionParams_projectors_LeftQ:
              (*primme)->correctionParams.projectors.LeftQ = *v.int_v;
      break;
      case PRIMMEF77_correctionParams_projectors_LeftX:
              (*primme)->correctionParams.projectors.LeftX = *v.int_v;
      break;
      case PRIMMEF77_correctionParams_projectors_RightQ:
              (*primme)->correctionParams.projectors.RightQ = *v.int_v;
      break;
      case PRIMMEF77_correctionParams_projectors_RightX:
              (*primme)->correctionParams.projectors.RightX = *v.int_v;
      break;
      case PRIMMEF77_correctionParams_projectors_SkewQ:
              (*primme)->correctionParams.projectors.SkewQ = *v.int_v;
      break;
      case PRIMMEF77_correctionParams_projectors_SkewX:
              (*primme)->correctionParams.projectors.SkewX = *v.int_v;
      break;
      case PRIMMEF77_correctionParams_convTest:
              (*primme)->correctionParams.convTest = *v.convergencetest_v;
      break;
      case PRIMMEF77_correctionParams_relTolBase:
              (*primme)->correctionParams.relTolBase = *v.double_v;
      break;
      case PRIMMEF77_stats_numOuterIterations:
              (*primme)->stats.numOuterIterations = *v.int_v;
      break;
      case PRIMMEF77_stats_numRestarts:
              (*primme)->stats.numRestarts = *v.int_v;
      break;
      case PRIMMEF77_stats_numMatvecs:
              (*primme)->stats.numMatvecs = *v.int_v;
      break;
      case PRIMMEF77_stats_numPreconds:
              (*primme)->stats.numPreconds = *v.int_v;
      break;
      case PRIMMEF77_stats_elapsedTime:
              (*primme)->stats.elapsedTime = *v.double_v;
      break;
      default : 
      fprintf(stderr,"Requested member (%d) does not exist: consult primme_f77.h.",*label);
      fprintf(stderr," No action taken \n");
   }

}

/*************************************************************************
 * primme_get_member_f77(primme_params *primme, int *label, void *ptr)
 *
 * Gets the value of any of the members of the primme data structure from 
 * a Fortran call.
 *
 * To be called ONLY from a f77 subroutine that is called by C PRIMME, 
 * (eg.,MV, PRECOND, GLOBALSUM) ---NOT from F77 driver.
 *
 *   F77 driver --> C PRIMME --> F77 function (MV) --> primme_get_member_f77
 *
 * For obtaining information from the driver, see below. The reason is that 
 * such f77 functions (MV, etc) receive the primme structure from PRIMME as
 * (primme_params *primme). However, the f77 driver has a pointer to **primme. 
 *
 * void *ptr generic type pointer is used for output, and is cast to the 
 * required pointer type before it is assigned the data structure member.
 *
 * The choice of member is through an integer label. The integers describing
 * the member are prefixed with PRIMMEF77_ followed by the actual name of 
 * the member. If the member is a structure, underscores are used for dots,
 *
 * Eg the label int for:    primme->correctionParams.projectors.LeftX
 * is:                    PRIMMEF77_correctionParams_projectors_LeftX
 *
 *************************************************************************/
#ifdef F77UNDERSCORE
void primme_get_member_f77_(primme_params *primme, int *label, union f77_value_ptr *v){
#else
void primme_get_member_f77(primme_params *primme, int *label, union f77_value_ptr *v){
#endif
   int i;

   switch (*label) {
      case PRIMMEF77_n:
              v->int_v = primme->n;
      break;
      case PRIMMEF77_matrixMatvec:
              v->matFunc_v = primme->matrixMatvec;
      break;
      case PRIMMEF77_massMatrixMatvec:
              v->matFunc_v = primme->massMatrixMatvec;
      break;
      case PRIMMEF77_applyPreconditioner:
              v->matFunc_v = primme->applyPreconditioner;
      break;
      case PRIMMEF77_numProcs:
              v->int_v = primme->numProcs;
      break;
      case PRIMMEF77_procID:
              v->int_v = primme->procID;
      break;
      case PRIMMEF77_commInfo:
              v->ptr_v = primme->commInfo;
      break;
      case PRIMMEF77_nLocal:
              v->int_v = primme->nLocal;
      break;
      case PRIMMEF77_globalSumDouble:
              v->globalSumDoubleFunc_v = primme->globalSumDouble;
      break;
      case PRIMMEF77_numEvals:
              v->int_v = primme->numEvals;
      break;
      case PRIMMEF77_target:
              v->target_v = primme->target;
      break;
      case PRIMMEF77_numTargetShifts:
              v->int_v = primme->numTargetShifts;
      break;
      case PRIMMEF77_targetShifts:
         for (i=0; i< primme->numTargetShifts; i++) {
             (&v->double_v)[i] = primme->targetShifts[i];
         }
      break;
      case PRIMMEF77_locking:
              v->int_v = primme->locking;
      break;
      case PRIMMEF77_initSize:
              v->int_v = primme->initSize;
      break;
      case PRIMMEF77_numOrthoConst:
              v->int_v = primme->numOrthoConst;
      break;
      case PRIMMEF77_dynamicMethodSwitch:
              v->int_v = primme->dynamicMethodSwitch;
      break;
      case PRIMMEF77_maxBasisSize:
              v->int_v = primme->maxBasisSize;
      break;
      case PRIMMEF77_minRestartSize:
              v->int_v = primme->minRestartSize;
      break;
      case PRIMMEF77_maxBlockSize:
              v->int_v = primme->maxBlockSize;
      break;
      case PRIMMEF77_maxMatvecs:
              v->int_v = primme->maxMatvecs;
      break;
      case PRIMMEF77_maxOuterIterations:
              v->int_v = primme->maxOuterIterations;
      break;
      case PRIMMEF77_intWorkSize:
              v->int_v = primme->intWorkSize;
      break;
      case PRIMMEF77_realWorkSize:
              v->long_int_v = primme->realWorkSize;
      break;
      case PRIMMEF77_iseed:
         for (i=0; i< 4; i++) {
            (&v->int_v)[i] = primme->iseed[i];
         }
      break;
      case PRIMMEF77_intWork:
              v->ptr_v = primme->intWork;
      break;
      case PRIMMEF77_realWork:
              v->ptr_v = primme->realWork;
      break;
      case PRIMMEF77_aNorm:
              v->double_v = primme->aNorm;
      break;
      case PRIMMEF77_eps:
              v->double_v = primme->eps;
      break;
      case PRIMMEF77_printLevel:
              v->int_v = primme->printLevel;
      break;
      case PRIMMEF77_outputFile:
              v->file_v = primme->outputFile;
      break;
      case PRIMMEF77_matrix:
              v->ptr_v = primme->matrix;
      break;
      case PRIMMEF77_preconditioner:
              v->ptr_v = primme->preconditioner;
      break;
      case PRIMMEF77_restartingParams_scheme:
              v->restartscheme_v = primme->restartingParams.scheme;
      break;
      case PRIMMEF77_restartingParams_maxPrevRetain:
              v->int_v = primme->restartingParams.maxPrevRetain;
      break;
      case PRIMMEF77_correctionParams_precondition:
              v->int_v = primme->correctionParams.precondition;
      break;
      case PRIMMEF77_correctionParams_robustShifts:
              v->int_v = primme->correctionParams.robustShifts;
      break;
      case PRIMMEF77_correctionParams_maxInnerIterations:
              v->int_v = primme->correctionParams.maxInnerIterations;
      break;
      case PRIMMEF77_correctionParams_projectors_LeftQ:
              v->int_v = primme->correctionParams.projectors.LeftQ;
      break;
      case PRIMMEF77_correctionParams_projectors_LeftX:
              v->int_v = primme->correctionParams.projectors.LeftX;
      break;
      case PRIMMEF77_correctionParams_projectors_RightQ:
              v->int_v = primme->correctionParams.projectors.RightQ;
      break;
      case PRIMMEF77_correctionParams_projectors_RightX:
              v->int_v = primme->correctionParams.projectors.RightX;
      break;
      case PRIMMEF77_correctionParams_projectors_SkewQ:
              v->int_v = primme->correctionParams.projectors.SkewQ;
      break;
      case PRIMMEF77_correctionParams_projectors_SkewX:
              v->int_v = primme->correctionParams.projectors.SkewX;
      break;
      case PRIMMEF77_correctionParams_convTest:
              v->convergencetest_v = primme->correctionParams.convTest;
      break;
      case PRIMMEF77_correctionParams_relTolBase:
              v->double_v = primme->correctionParams.relTolBase;
      break;
      case PRIMMEF77_stats_numOuterIterations:
              v->int_v = primme->stats.numOuterIterations;
      break;
      case PRIMMEF77_stats_numRestarts:
              v->int_v = primme->stats.numRestarts;
      break;
      case PRIMMEF77_stats_numMatvecs:
              v->int_v = primme->stats.numMatvecs;
      break;
      case PRIMMEF77_stats_numPreconds:
              v->int_v = primme->stats.numPreconds;
      break;
      case PRIMMEF77_stats_elapsedTime:
              v->double_v = primme->stats.elapsedTime;
      break;
      default :
      fprintf(stderr,"Requested member (%d) does not exist: consult primme_f77.h.",*label);
      fprintf(stderr," No action taken \n");
   }
}

/*************************************************************************
 *   If the user requires the current shift suggested by PRIMME to 
 *   be used in inverting (M-shift_i I), for the i-th eigenvalue
 *   this is at location (i-1) in C   (i-th location in Fortran)
 *   of the array primme->ShiftsForPreconditioner
 *************************************************************************/
#ifdef F77UNDERSCORE
void primme_get_prec_shift_f77_(primme_params *primme, int *i, double *shift)
#else
void primme_get_prec_shift_f77(primme_params *primme, int *i, double *shift)
#endif
{
   *shift = primme->ShiftsForPreconditioner[*i-1];
}

/*************************************************************************/
/* If the information from primme is needed from the calling driver in F77
 * the primme is passed as  (primme_params **), while any f77 function such
 * as MV or preconditioning that needs primme information it passes it as
 * (primme_params *). The reason is that PRIMME C code works with the pointer
 * to primme and that's what it calls MV/GLOBALSUM with.
 *
 *   F77 driver --> other F77 subs --> primmetop_get_member_f77
 * 
 * If the user needs to obtain structure info from the top function that 
 * initialized primme they should call the wrappers:
 *
 *      primmetop_get_member_f77
 *  or
 *      primmetop_get_prec_shift_f77   (probably not needed ever)
 *
 * From within MV, PRECONDITION, or GLOBALSUMDOUBLE call the above function
 *      primme_get_member_f77
 *
 *************************************************************************/
#ifdef F77UNDERSCORE
void primmetop_get_member_f77_(primme_params **primme, int *label, union f77_value_ptr *ptr){
   primme_get_member_f77_(*primme, label, ptr);
}
#else
void primmetop_get_member_f77(primme_params **primme, int *label, union f77_value_ptr *ptr){
   primme_get_member_f77(*primme, label, ptr);
}
#endif
/*************************************************************************/
#ifdef F77UNDERSCORE
void primmetop_get_prec_shift_f77_(primme_params **primme, int *i, double *shift)
{
   primme_get_prec_shift_f77_(*primme, i, shift);
}
#else
void primmetop_get_prec_shift_f77(primme_params **primme, int *i, double *shift)
{
   primme_get_prec_shift_f77(*primme, i, shift);
}
#endif
/*************************************************************************/

#ifdef __cplusplus
}
#endif
