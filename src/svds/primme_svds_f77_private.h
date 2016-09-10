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
 * File: primme_svds_f77_private.h
 *
 * Purpose - Definitions used exclusively by primme_svds_f77.c
 *
 ******************************************************************************/

#ifndef PRIMME_SVDS_F77_PRIVATE_H
#define PRIMME_SVDS_F77_PRIVATE_H

/*-------------------------------------------------------*/
/*     Defining easy to remember labels for setting the  */
/*     method in primme_svds_set_method from Fortran     */
/*-------------------------------------------------------*/
#define PRIMMEF77_SVDS_DEFAULT 0
#define PRIMMEF77_SVDS_HYBRID 1
#define PRIMMEF77_SVDS_NORMALEQUATIONS 2
#define PRIMMEF77_SVDS_AUGMENTED 3

/*-------------------------------------------------------*/
/*     Defining easy to remember labels for setting the  */
/*     members of the primme_svds structure from Fortran */
/*-------------------------------------------------------*/
#define PRIMMEF77_SVDS_primme 0
#define PRIMMEF77_SVDS_primmeStage2 1
#define PRIMMEF77_SVDS_m 2
#define PRIMMEF77_SVDS_n 3
#define PRIMMEF77_SVDS_matrixMatvec 4
#define PRIMMEF77_SVDS_applyPreconditioner 5
#define PRIMMEF77_SVDS_numProcs 6
#define PRIMMEF77_SVDS_procID 7
#define PRIMMEF77_SVDS_mLocal 8
#define PRIMMEF77_SVDS_nLocal 9
#define PRIMMEF77_SVDS_commInfo 10
#define PRIMMEF77_SVDS_globalSumReal 11
#define PRIMMEF77_SVDS_numSvals 12
#define PRIMMEF77_SVDS_target 13
#define PRIMMEF77_SVDS_numTargetShifts 14
#define PRIMMEF77_SVDS_targetShifts 15
#define PRIMMEF77_SVDS_method 16
#define PRIMMEF77_SVDS_methodStage2 17
#define PRIMMEF77_SVDS_intWorkSize 18
#define PRIMMEF77_SVDS_realWorkSize 19
#define PRIMMEF77_SVDS_intWork 20
#define PRIMMEF77_SVDS_realWork 21
#define PRIMMEF77_SVDS_matrix 22
#define PRIMMEF77_SVDS_preconditioner 23
#define PRIMMEF77_SVDS_locking 24
#define PRIMMEF77_SVDS_numOrthoConst 25
#define PRIMMEF77_SVDS_aNorm 26
#define PRIMMEF77_SVDS_eps 27
#define PRIMMEF77_SVDS_precondition 28
#define PRIMMEF77_SVDS_initSize 29
#define PRIMMEF77_SVDS_maxBasisSize 30
#define PRIMMEF77_SVDS_maxBlockSize 31
#define PRIMMEF77_SVDS_maxMatvecs 32
#define PRIMMEF77_SVDS_iseed 33
#define PRIMMEF77_SVDS_printLevel 34
#define PRIMMEF77_SVDS_outputFile 35
#define PRIMMEF77_SVDS_stats_numOuterIterations 36
#define PRIMMEF77_SVDS_stats_numRestarts 37
#define PRIMMEF77_SVDS_stats_numMatvecs 38
#define PRIMMEF77_SVDS_stats_numPreconds 39
#define PRIMMEF77_SVDS_stats_elapsedTime 40

/*-------------------------------------------------------*/
/*    Defining easy to remember labels for setting the   */
/*    enum members for targeting and operator            */
/*-------------------------------------------------------*/
#define PRIMMEF77_SVDS_largest 0
#define PRIMMEF77_SVDS_smallest 1
#define PRIMMEF77_SVDS_closest_abs 2
/*-------------------------------------------------------*/
#define PRIMMEF77_SVDS_op_none 0
#define PRIMMEF77_SVDS_op_AtA 1
#define PRIMMEF77_SVDS_op_AAt 2
#define PRIMMEF77_SVDS_op_augmented 3
 
#include "template.h"
#include "primme_svds_interface.h"

/* Prototypes for Fortran-C interface */

#ifdef __cplusplus
extern "C" {
#endif

union f77_value {
   PRIMME_INT *int_v;
   void (*matFunc_v) (void*,PRIMME_INT*,void*,PRIMME_INT*,int*,int*,struct primme_svds_params*,int*);
   void *ptr_v;
   void (*globalSumRealFunc_v) (void *,void *,int *,struct primme_svds_params*,int*);
   primme_svds_target *target_v;
   primme_svds_operator *operator_v;
   double *double_v;
   FILE *file_v;
};
union f77_value_ptr {
   PRIMME_INT int_v;
   void (*matFunc_v) (void*,PRIMME_INT*,void*,PRIMME_INT*,int*,int*,struct primme_svds_params*,int*);
   void *ptr_v;
   void (*globalSumRealFunc_v) (void *,void *,int *,struct primme_svds_params*,int*);
   primme_svds_target target_v;
   primme_svds_operator operator_v;
   double double_v;
   FILE *file_v;
};

#define AS_FORTRAN(X) AS_FORTRANX(X)
#define AS_FORTRANX(X) FORTRAN_FUNCTION(X ## _f77)

void AS_FORTRAN(Sprimme_svds)(REAL *svals, SCALAR *svecs,
      REAL *resNorms, primme_svds_params **primme_svds, int *ierr);

/* Only define these functions ones */
#ifdef USE_DOUBLE
void AS_FORTRAN(primme_svds_initialize)(primme_svds_params **primme_svds);
void AS_FORTRAN(primme_svds_set_method)(primme_svds_preset_method *method,
      primme_preset_method *methodStage1, primme_preset_method *methodStage2,
      primme_svds_params **primme_svds, int *ierr);
void AS_FORTRAN(primme_svds_display_params)(primme_svds_params **primme_svds);
void AS_FORTRAN(primme_svds_free)(primme_svds_params **primme_svds);
void AS_FORTRAN(primme_svds_set_member)(primme_svds_params **primme_svds, int *label, union f77_value ptr, int *ierr);
void AS_FORTRAN(primme_svdstop_get_member)(primme_svds_params **primme_svds, int *label, union f77_value_ptr *ptr, int *ierr);
void AS_FORTRAN(primme_svds_get_member)(primme_svds_params *primme_svds, int *label, union f77_value_ptr *ptr, int *ierr);
#endif

#ifdef __cplusplus
}
#endif

#endif /* PRIMME_SVDS_F77_PRIVATE_H */
