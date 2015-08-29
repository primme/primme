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
/*     method in primme_set_method from Fortran          */
/*-------------------------------------------------------*/
#define PRIMMEF77_DYNAMIC  0
#define PRIMMEF77_DEFAULT_MIN_TIME  1
#define PRIMMEF77_DEFAULT_MIN_MATVECS  2
#define PRIMMEF77_Arnoldi  3
#define PRIMMEF77_GD  4
#define PRIMMEF77_GD_plusK  5
#define PRIMMEF77_GD_Olsen_plusK  6
#define PRIMMEF77_JD_Olsen_plusK  7
#define PRIMMEF77_RQI  8
#define PRIMMEF77_JDQR  9
#define PRIMMEF77_JDQMR  10
#define PRIMMEF77_JDQMR_ETol  11
#define PRIMMEF77_SUBSPACE_ITERATION  12
#define PRIMMEF77_LOBPCG_OrthoBasis  13
#define PRIMMEF77_LOBPCG_OrthoBasis_Window  14
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
#define PRIMMEF77_globalSumDouble  7
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
#define PRIMMEF77_dynamicMethodSwitch 49
#define PRIMMEF77_massMatrixMatvec  50

/*----------------------------------------------------------*/
/*     Defining easy to remember labels for setting the     */
/*     enum members for targeting, restarting and innertest */
/*----------------------------------------------------------*/
#define PRIMMEF77_smallest  0
#define PRIMMEF77_largest  1
#define PRIMMEF77_closest_geq  2
#define PRIMMEF77_closest_leq  3
#define PRIMMEF77_closest_abs  4
/*-------------------------------------------------------*/
#define PRIMMEF77_thick  0
#define PRIMMEF77_dt  1
/*-------------------------------------------------------*/
#define PRIMMEF77_full_LTolerance  0
#define PRIMMEF77_decreasing_LTolerance  1
#define PRIMMEF77_adaptive_ETolerance  2
#define PRIMMEF77_adaptive  3


/* Prototypes for Fortran-C interface */

#ifdef __cplusplus
extern "C" {
#endif

union f77_value {
   int *int_v;
   void (*matFunc_v) (void *,void *,int *,struct primme_params *);
   void *ptr_v;
   void (*globalSumDoubleFunc_v) (void *,void *,int *,struct primme_params *);
   primme_target *target_v;
   double *double_v;
   long int *long_int_v;
   FILE *file_v;
   primme_restartscheme *restartscheme_v;
   primme_convergencetest *convergencetest_v;
};
union f77_value_ptr {
   int int_v;
   void (*matFunc_v) (void *,void *,int *,struct primme_params *);
   void *ptr_v;
   void (*globalSumDoubleFunc_v) (void *,void *,int *,struct primme_params *);
   primme_target target_v;
   double double_v;
   long int long_int_v;
   FILE *file_v;
   primme_restartscheme restartscheme_v;
   primme_convergencetest convergencetest_v;
};


#ifdef F77UNDERSCORE
void dprimme_f77_(double *evals, double *evecs, double *rnorms, 
                  primme_params **primme, int *ierr);
void zprimme_f77_(double *evals, Complex_Z *evecs, double *rnorms, 
                  primme_params **primme, int *ierr);
void primme_initialize_f77_(primme_params **primme);
void primme_free_f77_(primme_params **primme);
void primme_display_params_f77_(primme_params **primme);
void primme_printstacktrace_f77_(primme_params **primme);
void primme_set_method_f77_(primme_params **primme, int *method,int *returnVal);
void primme_set_member_f77_(primme_params **primme, int *label, union f77_value ptr);

void primme_get_prec_shift_f77_(primme_params *primme, int *i, double *shift);
void primme_get_member_f77_(primme_params *primme, int *label, union f77_value_ptr *ptr);
void primmetop_get_member_f77_(primme_params **primme, int *label, union f77_value_ptr *ptr);
void primmetop_get_prec_shift_f77_(primme_params **primme, int *i, 
                                   double *shift);
#else
void dprimme_f77(double *evals, double *evecs, double *rnorms, 
                 primme_params **primme, int *ierr);
void zprimme_f77(double *evals, Complex_Z *evecs, double *rnorms, 
                  primme_params **primme, int *ierr);
void primme_initialize_f77(primme_params **primme);
void primme_free_f77(primme_params **primme);
void primme_display_params_f77(primme_params **primme);
void primme_printstacktrace_f77(primme_params **primme);
void primme_set_method_f77(primme_params **primme, int *method, int *returnVal);
void primme_set_member_f77(primme_params **primme, int *label, union f77_value ptr);

void primme_get_prec_shift_f77(primme_params *primme, int *i, double *shift);
void primme_get_member_f77(primme_params *primme, int *label, union f77_value_ptr *ptr);
void primmetop_get_member_f77(primme_params **primme, int *label, union f77_value_ptr *ptr);
void primmetop_get_prec_shift_f77(primme_params **primme, int *i, 
                                  double *shift);
#endif

#ifdef __cplusplus
}
#endif

#endif /* PRIMME_F77_PRIVATE_H */
