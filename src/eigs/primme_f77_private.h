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
 * File: primme_f77_private.h
 *
 * Purpose - Definitions used exclusively by primme_f77.c
 *
 ******************************************************************************/

#ifndef PRIMME_F77_PRIVATE_H
#define PRIMME_F77_PRIVATE_H

/*-------------------------------------------------------*/
/*     Defining easy to remember labels for setting the  */
/*     members of the primme structure from Fortran      */
/*-------------------------------------------------------*/
#define PRIMMEF77_n  0
#define PRIMMEF77_matrixMatvec  1
#define PRIMMEF77_applyPreconditioner  2
#define PRIMMEF77_numProcs  3
#define PRIMMEF77_procID  4
#define PRIMMEF77_commInfo  5
#define PRIMMEF77_nLocal  6
#define PRIMMEF77_globalSumReal  7
#define PRIMMEF77_numEvals  8
#define PRIMMEF77_target  9
#define PRIMMEF77_numTargetShifts  10
#define PRIMMEF77_targetShifts  11
#define PRIMMEF77_locking  12
#define PRIMMEF77_initSize  13
#define PRIMMEF77_numOrthoConst  14
#define PRIMMEF77_maxBasisSize  15
#define PRIMMEF77_minRestartSize  16
#define PRIMMEF77_maxBlockSize  17
#define PRIMMEF77_maxMatvecs  18
#define PRIMMEF77_maxOuterIterations  19
#define PRIMMEF77_intWorkSize  20
#define PRIMMEF77_realWorkSize  21
#define PRIMMEF77_iseed  22
#define PRIMMEF77_intWork  23
#define PRIMMEF77_realWork  24
#define PRIMMEF77_aNorm  25
#define PRIMMEF77_eps  26
#define PRIMMEF77_printLevel  27
#define PRIMMEF77_outputFile  28
#define PRIMMEF77_matrix  29
#define PRIMMEF77_preconditioner  30
#define PRIMMEF77_initBasisMode   301
#define PRIMMEF77_projectionParams_projection  302
#define PRIMMEF77_restartingParams_scheme  31
#define PRIMMEF77_restartingParams_maxPrevRetain  32
#define PRIMMEF77_correctionParams_precondition  33
#define PRIMMEF77_correctionParams_robustShifts  34
#define PRIMMEF77_correctionParams_maxInnerIterations  35
#define PRIMMEF77_correctionParams_projectors_LeftQ  36
#define PRIMMEF77_correctionParams_projectors_LeftX  37
#define PRIMMEF77_correctionParams_projectors_RightQ  38
#define PRIMMEF77_correctionParams_projectors_RightX  39
#define PRIMMEF77_correctionParams_projectors_SkewQ  40
#define PRIMMEF77_correctionParams_projectors_SkewX  41
#define PRIMMEF77_correctionParams_convTest  42
#define PRIMMEF77_correctionParams_relTolBase  43
#define PRIMMEF77_stats_numOuterIterations  44
#define PRIMMEF77_stats_numRestarts  45
#define PRIMMEF77_stats_numMatvecs  46
#define PRIMMEF77_stats_numPreconds  47
#define PRIMMEF77_stats_elapsedTime  48
#define PRIMMEF77_stats_estimateMinEVal  481
#define PRIMMEF77_stats_estimateMaxEVal  482
#define PRIMMEF77_stats_estimateLargestSVal  483
#define PRIMMEF77_stats_maxConvTol  484
#define PRIMMEF77_dynamicMethodSwitch 49
#define PRIMMEF77_massMatrixMatvec  50
#define PRIMMEF77_convTestFun  51
#define PRIMMEF77_ldevecs  52
#define PRIMMEF77_ldOPs  53

#include "template.h"
#include "primme_interface.h" /* for Sprimme */

/* Prototypes for Fortran-C interface */

#ifdef __cplusplus
extern "C" {
#endif

union f77_value {
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
};
union f77_value_ptr {
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
};

#define AS_FORTRAN(X) AS_FORTRANX(X)
#define AS_FORTRANX(X) FORTRAN_FUNCTION(X ## _f77)

void AS_FORTRAN(Sprimme)(REAL *evals, SCALAR *evecs,
      REAL *rnorms, primme_params **primme, int *ierr);

/* Only define these functions ones */
#ifdef USE_DOUBLE
void AS_FORTRAN(primme_initialize)(primme_params **primme);
void AS_FORTRAN(primme_free)(primme_params **primme);
void AS_FORTRAN(primme_display_params)(primme_params **primme);
void AS_FORTRAN(primme_printstacktrace)(primme_params **primme);
void AS_FORTRAN(primme_set_method)(primme_params **primme, primme_preset_method *method, int *returnValue);
void AS_FORTRAN(primme_set_member)(primme_params **primme, int *label, union f77_value ptr, int *ierr);
void AS_FORTRAN(primme_get_prec_shift)(primme_params *primme, int *i, double *shift);
void AS_FORTRAN(primme_get_member)(primme_params *primme, int *label, union f77_value_ptr *ptr, int *ierr);
void AS_FORTRAN(primmetop_get_member)(primme_params **primme, int *label, union f77_value_ptr *ptr, int *ierr);
void AS_FORTRAN(primmetop_get_prec_shift)(primme_params **primme, int *i, double *shift);
#endif

#ifdef __cplusplus
}
#endif

#endif /* PRIMME_F77_PRIVATE_H */
