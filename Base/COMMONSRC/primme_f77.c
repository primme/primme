/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2005  James R. McCombs,  Andreas Stathopoulos
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
 * --------------------------------------------------------------------------
 */
#include "primme.h"
#include "primme_f77_private.h"
#ifdef Cplusplus
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
void primme_set_member_f77_(primme_params **primme, int *label, void *ptr){
#else
void primme_set_member_f77(primme_params **primme, int *label, void *ptr){
#endif
   int i;

   switch (*label) {
      case PRIMMEF77_n:
              (*primme)->n = *(int*) ptr;
      break;
      case PRIMMEF77_matrixMatvec:
              (*primme)->matrixMatvec = (void (*)
                (void *,void *,int *,struct primme_params *) ) ptr;
      break;
      case PRIMMEF77_massMatrixMatvec:
              (*primme)->massMatrixMatvec = (void (*)
                (void *,void *,int *,struct primme_params *) ) ptr;
      break;
      case PRIMMEF77_applyPreconditioner:
              (*primme)->applyPreconditioner = (void (*)
                (void *,void *,int *,struct primme_params *) ) ptr;
      break;
      case PRIMMEF77_numProcs:
              (*primme)->numProcs = *(int*) ptr;
      break;
      case PRIMMEF77_procID:
              (*primme)->procID = *(int*) ptr;
      break;
      case PRIMMEF77_commInfo:
              (*primme)->commInfo =  ptr;
      break;
      case PRIMMEF77_nLocal:
              (*primme)->nLocal = *(int*) ptr;
      break;
      case PRIMMEF77_globalSumDouble:
              (*primme)->globalSumDouble = (void (*)
                (void *,void *,int *,struct primme_params *) ) ptr;
      break;
      case PRIMMEF77_numEvals:
              (*primme)->numEvals = *(int*) ptr;
      break;
      case PRIMMEF77_target:
              (*primme)->target = *(primme_target*) ptr;
      break;
      case PRIMMEF77_numTargetShifts:
              (*primme)->numTargetShifts = *(int*) ptr;
      break;
      case PRIMMEF77_targetShifts:
              (*primme)->targetShifts = (double*) ptr;
      break;
      case PRIMMEF77_locking:
              (*primme)->locking = *(int*) ptr;
      break;
      case PRIMMEF77_initSize:
              (*primme)->initSize = *(int*) ptr;
      break;
      case PRIMMEF77_numOrthoConst:
              (*primme)->numOrthoConst = *(int*) ptr;
      break;
      case PRIMMEF77_dynamicMethodSwitch:
              (*primme)->dynamicMethodSwitch = *(int*) ptr;
      break;
      case PRIMMEF77_maxBasisSize:
              (*primme)->maxBasisSize = *(int*) ptr;
      break;
      case PRIMMEF77_minRestartSize:
              (*primme)->minRestartSize = *(int*) ptr;
      break;
      case PRIMMEF77_maxBlockSize:
              (*primme)->maxBlockSize = *(int*) ptr;
      break;
      case PRIMMEF77_maxMatvecs:
              (*primme)->maxMatvecs = *(int*) ptr;
      break;
      case PRIMMEF77_maxOuterIterations:
              (*primme)->maxOuterIterations = *(int*) ptr;
      break;
      case PRIMMEF77_intWorkSize:
              (*primme)->intWorkSize = *(int*) ptr;
      break;
      case PRIMMEF77_realWorkSize:
              (*primme)->realWorkSize = *(long int*) ptr;
      break;
      case PRIMMEF77_iseed:
         for (i=0; i< 4; i++) {
            (*primme)->iseed[i] = ((int *) ptr)[i];
         }
      break;
      case PRIMMEF77_intWork:
              (*primme)->intWork = (int*) ptr;
      break;
      case PRIMMEF77_realWork:
              (*primme)->realWork = (void*) ptr;
      break;
      case PRIMMEF77_aNorm:
              (*primme)->aNorm = *(double*) ptr;
      break;
      case PRIMMEF77_eps:
              (*primme)->eps = *(double*) ptr;
      break;
      case PRIMMEF77_printLevel:
              (*primme)->printLevel = *(int*) ptr;
      break;
      case PRIMMEF77_outputFile:
              (*primme)->outputFile = (FILE*) ptr;
      break;
      case PRIMMEF77_matrix:
              (*primme)->matrix = ptr;
      break;
      case PRIMMEF77_preconditioner:
              (*primme)->preconditioner = ptr;
      break;
      case PRIMMEF77_restartingParams_scheme:
              (*primme)->restartingParams.scheme = *(primme_restartscheme*) ptr;
      break;
      case PRIMMEF77_restartingParams_maxPrevRetain:
              (*primme)->restartingParams.maxPrevRetain = *(int*) ptr;
      break;
      case PRIMMEF77_correctionParams_precondition:
              (*primme)->correctionParams.precondition = *(int*) ptr;
      break;
      case PRIMMEF77_correctionParams_robustShifts:
              (*primme)->correctionParams.robustShifts = *(int*) ptr;
      break;
      case PRIMMEF77_correctionParams_maxInnerIterations:
              (*primme)->correctionParams.maxInnerIterations = *(int*) ptr;
      break;
      case PRIMMEF77_correctionParams_projectors_LeftQ:
              (*primme)->correctionParams.projectors.LeftQ = *(int*) ptr;
      break;
      case PRIMMEF77_correctionParams_projectors_LeftX:
              (*primme)->correctionParams.projectors.LeftX = *(int*) ptr;
      break;
      case PRIMMEF77_correctionParams_projectors_RightQ:
              (*primme)->correctionParams.projectors.RightQ = *(int*) ptr;
      break;
      case PRIMMEF77_correctionParams_projectors_RightX:
              (*primme)->correctionParams.projectors.RightX = *(int*) ptr;
      break;
      case PRIMMEF77_correctionParams_projectors_SkewQ:
              (*primme)->correctionParams.projectors.SkewQ = *(int*) ptr;
      break;
      case PRIMMEF77_correctionParams_projectors_SkewX:
              (*primme)->correctionParams.projectors.SkewX = *(int*) ptr;
      break;
      case PRIMMEF77_correctionParams_convTest:
              (*primme)->correctionParams.convTest = 
                                        *(primme_convergencetest*) ptr;
      break;
      case PRIMMEF77_correctionParams_relTolBase:
              (*primme)->correctionParams.relTolBase = *(double*) ptr;
      break;
      case PRIMMEF77_stats_numOuterIterations:
              (*primme)->stats.numOuterIterations = *(int*) ptr;
      break;
      case PRIMMEF77_stats_numRestarts:
              (*primme)->stats.numRestarts = *(int*) ptr;
      break;
      case PRIMMEF77_stats_numMatvecs:
              (*primme)->stats.numMatvecs = *(int*) ptr;
      break;
      case PRIMMEF77_stats_numPreconds:
              (*primme)->stats.numPreconds = *(int*) ptr;
      break;
      case PRIMMEF77_stats_elapsedTime:
              (*primme)->stats.elapsedTime = *(double*) ptr;
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
void primme_get_member_f77_(primme_params *primme, int *label, void *ptr){
#else
void primme_get_member_f77(primme_params *primme, int *label, void *ptr){
#endif
   int i;

   switch (*label) {
      case PRIMMEF77_n:
              *(int*) ptr = primme->n;
      break;
      case PRIMMEF77_matrixMatvec:
              *(void **) ptr = (void *) primme->matrixMatvec;
      break;
      case PRIMMEF77_massMatrixMatvec:
              *(void **) ptr = (void *) primme->massMatrixMatvec;
      break;
      case PRIMMEF77_applyPreconditioner:
              *(void **) ptr = (void *) primme->applyPreconditioner;
      break;
      case PRIMMEF77_numProcs:
              *(int*) ptr = primme->numProcs;
      break;
      case PRIMMEF77_procID:
              *(int*) ptr = primme->procID;
      break;
      case PRIMMEF77_commInfo:
              *(void **) ptr = primme->commInfo;
      break;
      case PRIMMEF77_nLocal:
              *(int*) ptr = primme->nLocal;
      break;
      case PRIMMEF77_globalSumDouble:
              *(void **) ptr = (void *) primme->globalSumDouble;
      break;
      case PRIMMEF77_numEvals:
              *(int*) ptr = primme->numEvals;
      break;
      case PRIMMEF77_target:
              *(primme_target*) ptr = primme->target;
      break;
      case PRIMMEF77_numTargetShifts:
              *(int*) ptr = primme->numTargetShifts;
      break;
      case PRIMMEF77_targetShifts:
         for (i=0; i< primme->numTargetShifts; i++) {
             ((double*) ptr)[i] = primme->targetShifts[i];
         }
      break;
      case PRIMMEF77_locking:
              *(int*) ptr = primme->locking;
      break;
      case PRIMMEF77_initSize:
              *(int*) ptr = primme->initSize;
      break;
      case PRIMMEF77_numOrthoConst:
              *(int*) ptr = primme->numOrthoConst;
      break;
      case PRIMMEF77_dynamicMethodSwitch:
              *(int*) ptr = primme->dynamicMethodSwitch;
      break;
      case PRIMMEF77_maxBasisSize:
              *(int*) ptr = primme->maxBasisSize;
      break;
      case PRIMMEF77_minRestartSize:
              *(int*) ptr = primme->minRestartSize;
      break;
      case PRIMMEF77_maxBlockSize:
              *(int*) ptr = primme->maxBlockSize;
      break;
      case PRIMMEF77_maxMatvecs:
              *(int*) ptr = primme->maxMatvecs;
      break;
      case PRIMMEF77_maxOuterIterations:
              *(int*) ptr = primme->maxOuterIterations;
      break;
      case PRIMMEF77_intWorkSize:
              *(int*) ptr = primme->intWorkSize;
      break;
      case PRIMMEF77_realWorkSize:
              *(long int*) ptr = primme->realWorkSize;
      break;
      case PRIMMEF77_iseed:
         for (i=0; i< 4; i++) {
            ((double*) ptr)[i] = primme->iseed[i];
         }
      break;
      case PRIMMEF77_intWork:
              *(void **) ptr = (void *) primme->intWork;
      break;
      case PRIMMEF77_realWork:
              *(void **) ptr = (void *) primme->realWork;
      break;
      case PRIMMEF77_aNorm:
              *(double*) ptr = primme->aNorm;
      break;
      case PRIMMEF77_eps:
              *(double*) ptr = primme->eps;
      break;
      case PRIMMEF77_printLevel:
              *(int*) ptr = primme->printLevel;
      break;
      case PRIMMEF77_outputFile:
              *(void **) ptr = (void *) primme->outputFile;
      break;
      case PRIMMEF77_matrix:
              *(void **) ptr = primme->matrix;
      break;
      case PRIMMEF77_preconditioner:
              *(void **) ptr = primme->preconditioner;
      break;
      case PRIMMEF77_restartingParams_scheme:
              *(primme_restartscheme*) ptr = primme->restartingParams.scheme;
      break;
      case PRIMMEF77_restartingParams_maxPrevRetain:
              *(int*) ptr = primme->restartingParams.maxPrevRetain;
      break;
      case PRIMMEF77_correctionParams_precondition:
              *(int*) ptr = primme->correctionParams.precondition;
      break;
      case PRIMMEF77_correctionParams_robustShifts:
              *(int*) ptr = primme->correctionParams.robustShifts;
      break;
      case PRIMMEF77_correctionParams_maxInnerIterations:
              *(int*) ptr = primme->correctionParams.maxInnerIterations;
      break;
      case PRIMMEF77_correctionParams_projectors_LeftQ:
              *(int*) ptr = primme->correctionParams.projectors.LeftQ;
      break;
      case PRIMMEF77_correctionParams_projectors_LeftX:
              *(int*) ptr = primme->correctionParams.projectors.LeftX;
      break;
      case PRIMMEF77_correctionParams_projectors_RightQ:
              *(int*) ptr = primme->correctionParams.projectors.RightQ;
      break;
      case PRIMMEF77_correctionParams_projectors_RightX:
              *(int*) ptr = primme->correctionParams.projectors.RightX;
      break;
      case PRIMMEF77_correctionParams_projectors_SkewQ:
              *(int*) ptr = primme->correctionParams.projectors.SkewQ;
      break;
      case PRIMMEF77_correctionParams_projectors_SkewX:
              *(int*) ptr = primme->correctionParams.projectors.SkewX;
      break;
      case PRIMMEF77_correctionParams_convTest:
             *(primme_convergencetest*) ptr = primme->correctionParams.convTest;
      break;
      case PRIMMEF77_correctionParams_relTolBase:
              *(double*) ptr = primme->correctionParams.relTolBase;
      break;
      case PRIMMEF77_stats_numOuterIterations:
              *(int*) ptr = primme->stats.numOuterIterations;
      break;
      case PRIMMEF77_stats_numRestarts:
              *(int*) ptr = primme->stats.numRestarts;
      break;
      case PRIMMEF77_stats_numMatvecs:
              *(int*) ptr = primme->stats.numMatvecs;
      break;
      case PRIMMEF77_stats_numPreconds:
              *(int*) ptr = primme->stats.numPreconds;
      break;
      case PRIMMEF77_stats_elapsedTime:
              *(double*) ptr = primme->stats.elapsedTime;
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
void primmetop_get_member_f77_(primme_params **primme, int *label, void *ptr){
   primme_get_member_f77_(*primme, label, ptr);
}
#else
void primmetop_get_member_f77(primme_params **primme, int *label, void *ptr){
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

#ifdef Cplusplus
}
#endif
