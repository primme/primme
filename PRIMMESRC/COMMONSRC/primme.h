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
 * File: primme.h
 * 
 * Purpose - Main header with the PRIMME C interface functions.
 * 
 ******************************************************************************/

#ifndef PRIMME_H
#define PRIMME_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <sys/types.h>
#include <limits.h>
#include "Complexz.h"

#define PRIMME_MAX_NAME_LENGTH 128

typedef enum {
   Primme_dprimme,
   Primme_zprimme,
   Primme_check_input,
   Primme_allocate_workspace,
   Primme_main_iter,
   Primme_init_basis,
   Primme_init_block_krylov,
   Primme_init_krylov,
   Primme_ortho,
   Primme_solve_h,
   Primme_restart,
   Primme_restart_h,
   Primme_insert_submatrix,
   Primme_lock_vectors,
   Primme_num_dsyev,
   Primme_num_zheev,
   Primme_num_dspev,
   Primme_num_zhpev,
   Primme_ududecompose,
   Primme_udusolve,
   Primme_apply_projected_preconditioner,
   Primme_apply_skew_projector,
   Primme_inner_solve,
   Primme_solve_correction,
   Primme_fopen,
   Primme_malloc
} primme_function;


typedef enum {
   primme_smallest,
   primme_largest,
   primme_closest_geq,
   primme_closest_leq,
   primme_closest_abs
} primme_target;


typedef enum {
   primme_thick,
   primme_dtr
} primme_restartscheme;


typedef enum {
   primme_full_LTolerance,
   primme_decreasing_LTolerance,
   primme_adaptive_ETolerance,
   primme_adaptive
} primme_convergencetest;


typedef struct stackTraceNode {
   primme_function callingFunction;
   primme_function failedFunction;
   int errorCode;
   int lineNumber;
   char fileName[PRIMME_MAX_NAME_LENGTH];
   struct stackTraceNode *nextNode;
} stackTraceNode;


typedef struct primme_stats {
   int numOuterIterations;
   int numRestarts;
   int numMatvecs;
   int numPreconds;
   double elapsedTime; 
} primme_stats;
   
typedef struct JD_projectors {
   int LeftQ;
   int LeftX;
   int RightQ;
   int RightX;
   int SkewQ;
   int SkewX;
} JD_projectors;

typedef struct correction_params {
   int precondition;
   int robustShifts;
   int maxInnerIterations;
   struct JD_projectors projectors;
   primme_convergencetest convTest;
   double relTolBase;
} correction_params;


typedef struct restarting_params {
   primme_restartscheme scheme;
   int maxPrevRetain;
} restarting_params;
   

/*--------------------------------------------------------------------------*/
typedef struct primme_params {

   /* The user must input at least the following two arguments */
   int n;
   void (*matrixMatvec)
      ( void *x,  void *y, int *blockSize, struct primme_params *primme);

   /* Preconditioner applied on block of vectors (if available) */
   void (*applyPreconditioner)
      ( void *x,  void *y, int *blockSize, struct primme_params *primme);

   /* Matrix times a multivector for mass matrix B for generalized Ax = xBl */
   void (*massMatrixMatvec)
      ( void *x,  void *y, int *blockSize, struct primme_params *primme);

   /* input for the following is only required for parallel programs */
   int numProcs;
   int procID;
   int nLocal;
   void *commInfo;
   void (*globalSumDouble)
      (void *sendBuf, void *recvBuf, int *count, struct primme_params *primme );

   /*Though primme_initialize will assign defaults, most users will set these */
   int numEvals;          
   primme_target target; 
   int numTargetShifts;              /* For targeting interior epairs,      */
   double *targetShifts;             /* at least one shift must also be set */

   /* the following will be given default values depending on the method */
   int dynamicMethodSwitch;
   int locking;
   int initSize;
   int numOrthoConst;
   int maxBasisSize;
   int minRestartSize;
   int maxBlockSize;
   int maxMatvecs;
   int maxOuterIterations;
   int intWorkSize;
   long int realWorkSize;
   int iseed[4];
   int *intWork;
   void *realWork;
   double aNorm;
   double eps;

   int printLevel;
   FILE *outputFile;
   
   void *matrix;
   void *preconditioner;
   double *ShiftsForPreconditioner;

   struct restarting_params restartingParams;
   struct correction_params correctionParams;
   struct primme_stats stats;
   struct stackTraceNode *stackTrace;
   
} primme_params;
/*---------------------------------------------------------------------------*/

typedef enum {
   DYNAMIC,
   DEFAULT_MIN_TIME,
   DEFAULT_MIN_MATVECS,
   Arnoldi,
   GD,
   GD_plusK,
   GD_Olsen_plusK,
   JD_Olsen_plusK,
   RQI,
   JDQR,
   JDQMR,
   JDQMR_ETol,
   SUBSPACE_ITERATION,
   LOBPCG_OrthoBasis,
   LOBPCG_OrthoBasis_Window
} primme_preset_method;



int dprimme(double *evals, double *evecs, double *resNorms, 
            primme_params *primme);
int zprimme(double *evals, Complex_Z *evecs, double *resNorms, 
            primme_params *primme);
void primme_initialize(primme_params *primme);
int  primme_set_method(primme_preset_method method, primme_params *params);
void primme_display_params(primme_params primme);
void *primme_valloc(size_t byteSize, const char *target);
void *primme_calloc(size_t nelem, size_t elsize, const char *target);
void primme_Free(primme_params *primme);
void primme_seq_globalSumDouble(void *sendBuf, void *recvBuf, int *count,
                                                   primme_params *params);
void primme_PushErrorMessage(const primme_function callingFunction, 
     const primme_function failedFunction, const int errorCode, 
     const char *fileName, const int lineNumber, primme_params *primme);
void primme_PrintStackTrace(const primme_params primme);
void primme_DeleteStackTrace(primme_params *primme);

#ifdef __cplusplus
}
#endif


#endif /* PRIMME_H */
