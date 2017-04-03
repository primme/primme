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
 *******************************************************************************
 * File: primme_eigs.h
 * 
 * Purpose - Main header with the PRIMME EIGS C interface functions.
 * 
 ******************************************************************************/

#ifndef PRIMME_EIGS_H
#define PRIMME_EIGS_H

#include <stdio.h>

#include "primme.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
   primme_smallest,        /* leftmost eigenvalues */
   primme_largest,         /* rightmost eigenvalues */
   primme_closest_geq,     /* leftmost but greater than the target shift */
   primme_closest_leq,     /* rightmost but less than the target shift */
   primme_closest_abs,     /* the closest to the target shift */
   primme_largest_abs      /* the farthest to the target shift */
} primme_target;

/* projection methods for extraction */
typedef enum {
   primme_proj_default,
   primme_proj_RR,          /* Rayleigh-Ritz */
   primme_proj_harmonic,    /* Harmonic Rayleigh-Ritz */
   primme_proj_refined      /* refined with fixed target */
} primme_projection;

typedef enum {         /* Initially fill up the search subspace with: */
   primme_init_default,
   primme_init_krylov, /* a) Krylov with the last vector provided by the user or random */
   primme_init_random, /* b) just random vectors */
   primme_init_user    /* c) provided vectors or a single random vector */
} primme_init;


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


/* Identifies the type of event for which monitor is being called */
typedef enum {
   primme_event_outer_iteration,    /* report at every outer iteration        */
   primme_event_inner_iteration,    /* report at every QMR iteration          */
   primme_event_restart,            /* report at every basis restart          */
   primme_event_reset,              /* event launch if basis reset            */
   primme_event_converged,          /* report new pair marked as converged    */
   primme_event_locked              /* report new pair marked as locked       */
} primme_event;

typedef struct primme_stats {
   PRIMME_INT numOuterIterations;
   PRIMME_INT numRestarts;
   PRIMME_INT numMatvecs;
   PRIMME_INT numPreconds;
   PRIMME_INT numGlobalSum;         /* times called globalSumReal */
   PRIMME_INT volumeGlobalSum;      /* number of SCALARs reduced by globalSumReal */
   double numOrthoInnerProds;       /* number of inner prods done by Ortho */
   double elapsedTime; 
   double timeMatvec;               /* time expend by matrixMatvec */
   double timePrecond;              /* time expend by applyPreconditioner */
   double timeOrtho;                /* time expend by ortho  */
   double timeGlobalSum;            /* time expend by globalSumReal  */
   double estimateMinEVal;          /* the leftmost Ritz value seen */
   double estimateMaxEVal;          /* the rightmost Ritz value seen */
   double estimateLargestSVal;      /* absolute value of the farthest to zero Ritz value seen */
   double maxConvTol;               /* largest norm residual of a locked eigenpair */
   double estimateResidualError;    /* accumulated error in V and W */
} primme_stats;

typedef struct JD_projectors {
   int LeftQ;
   int LeftX;
   int RightQ;
   int RightX;
   int SkewQ;
   int SkewX;
} JD_projectors;

typedef struct projection_params {
   primme_projection projection;
} projection_params;

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
   PRIMME_INT n;
   void (*matrixMatvec)
      ( void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
        struct primme_params *primme, int *ierr);

   /* Preconditioner applied on block of vectors (if available) */
   void (*applyPreconditioner)
      ( void *x, PRIMME_INT *ldx,  void *y, PRIMME_INT *ldy, int *blockSize,
        struct primme_params *primme, int *ierr);

   /* Matrix times a multivector for mass matrix B for generalized Ax = xBl */
   void (*massMatrixMatvec)
      ( void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,
        struct primme_params *primme, int *ierr);

   /* input for the following is only required for parallel programs */
   int numProcs;
   int procID;
   PRIMME_INT nLocal;
   void *commInfo;
   void (*globalSumReal)
      (void *sendBuf, void *recvBuf, int *count, struct primme_params *primme,
       int *ierr );

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
   PRIMME_INT maxMatvecs;
   PRIMME_INT maxOuterIterations;
   int intWorkSize;
   size_t realWorkSize;
   PRIMME_INT iseed[4];
   int *intWork;
   void *realWork;
   double aNorm;
   double eps;

   int printLevel;
   FILE *outputFile;

   void *matrix;
   void *preconditioner;
   double *ShiftsForPreconditioner;
   primme_init initBasisMode;
   PRIMME_INT ldevecs;
   PRIMME_INT ldOPs;

   struct projection_params projectionParams; 
   struct restarting_params restartingParams;
   struct correction_params correctionParams;
   struct primme_stats stats;

   void (*convTestFun)(double *eval, void *evec, double *rNorm, int *isconv, 
         struct primme_params *primme, int *ierr);
   void *convtest;
   void (*monitorFun)(void *basisEvals, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms, int *numConverged,
      void *lockedEvals, int *numLocked, int *lockedFlags, void *lockedNorms,
      int *inner_its, void *LSRes, primme_event *event,
      struct primme_params *primme, int *err);
   void *monitor;
} primme_params;
/*---------------------------------------------------------------------------*/

typedef enum {
   PRIMME_DEFAULT_METHOD,
   PRIMME_DYNAMIC,
   PRIMME_DEFAULT_MIN_TIME,
   PRIMME_DEFAULT_MIN_MATVECS,
   PRIMME_Arnoldi,
   PRIMME_GD,
   PRIMME_GD_plusK,
   PRIMME_GD_Olsen_plusK,
   PRIMME_JD_Olsen_plusK,
   PRIMME_RQI,
   PRIMME_JDQR,
   PRIMME_JDQMR,
   PRIMME_JDQMR_ETol,
   PRIMME_STEEPEST_DESCENT,
   PRIMME_LOBPCG_OrthoBasis,
   PRIMME_LOBPCG_OrthoBasis_Window
} primme_preset_method;

typedef enum {
   primme_int,
   primme_double,
   primme_pointer
} primme_type;

typedef enum {
   PRIMME_n =  0,
   PRIMME_matrixMatvec =  1,
   PRIMME_applyPreconditioner =  2,
   PRIMME_numProcs =  3,
   PRIMME_procID =  4,
   PRIMME_commInfo =  5,
   PRIMME_nLocal =  6,
   PRIMME_globalSumReal =  7,
   PRIMME_numEvals =  8,
   PRIMME_target =  9,
   PRIMME_numTargetShifts =  10,
   PRIMME_targetShifts =  11,
   PRIMME_locking =  12,
   PRIMME_initSize =  13,
   PRIMME_numOrthoConst =  14,
   PRIMME_maxBasisSize =  15,
   PRIMME_minRestartSize =  16,
   PRIMME_maxBlockSize =  17,
   PRIMME_maxMatvecs =  18,
   PRIMME_maxOuterIterations =  19,
   PRIMME_intWorkSize =  20,
   PRIMME_realWorkSize =  21,
   PRIMME_iseed =  22,
   PRIMME_intWork =  23,
   PRIMME_realWork =  24,
   PRIMME_aNorm =  25,
   PRIMME_eps =  26,
   PRIMME_printLevel =  27,
   PRIMME_outputFile =  28,
   PRIMME_matrix =  29,
   PRIMME_preconditioner =  30,
   PRIMME_initBasisMode =   301,
   PRIMME_projectionParams_projection =  302,
   PRIMME_restartingParams_scheme =  31,
   PRIMME_restartingParams_maxPrevRetain =  32,
   PRIMME_correctionParams_precondition =  33,
   PRIMME_correctionParams_robustShifts =  34,
   PRIMME_correctionParams_maxInnerIterations =  35,
   PRIMME_correctionParams_projectors_LeftQ =  36,
   PRIMME_correctionParams_projectors_LeftX =  37,
   PRIMME_correctionParams_projectors_RightQ =  38,
   PRIMME_correctionParams_projectors_RightX =  39,
   PRIMME_correctionParams_projectors_SkewQ =  40,
   PRIMME_correctionParams_projectors_SkewX =  41,
   PRIMME_correctionParams_convTest =  42,
   PRIMME_correctionParams_relTolBase =  43,
   PRIMME_stats_numOuterIterations =  44,
   PRIMME_stats_numRestarts =  45,
   PRIMME_stats_numMatvecs =  46,
   PRIMME_stats_numPreconds =  47,
   PRIMME_stats_numGlobalSum =  471,
   PRIMME_stats_volumeGlobalSum =  472,
   PRIMME_stats_numOrthoInnerProds =  473,
   PRIMME_stats_elapsedTime =  48,
   PRIMME_stats_timeMatvec =  4801,
   PRIMME_stats_timePrecond =  4802,
   PRIMME_stats_timeOrtho =  4803,
   PRIMME_stats_timeGlobalSum =  4804,
   PRIMME_stats_estimateMinEVal =  481,
   PRIMME_stats_estimateMaxEVal =  482,
   PRIMME_stats_estimateLargestSVal =  483,
   PRIMME_stats_maxConvTol =  484,
   PRIMME_dynamicMethodSwitch = 49,
   PRIMME_massMatrixMatvec =  50,
   PRIMME_convTestFun =  51,
   PRIMME_ldevecs =  52,
   PRIMME_ldOPs =  53,
   PRIMME_monitorFun = 54,
   PRIMME_monitor = 55
} primme_params_label;

int sprimme(float *evals, float *evecs, float *resNorms, 
      primme_params *primme);
int cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, 
      primme_params *primme);
int dprimme(double *evals, double *evecs, double *resNorms, 
      primme_params *primme);
int zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, 
      primme_params *primme);
void primme_initialize(primme_params *primme);
int  primme_set_method(primme_preset_method method, primme_params *params);
void primme_display_params(primme_params primme);
void primme_free(primme_params *primme);
int primme_get_member(primme_params *primme, primme_params_label label,
      void *value);
int primme_set_member(primme_params *primme, primme_params_label label,
      void *value);
int primme_member_info(primme_params_label *label, const char** label_name,
      primme_type *type, int *arity);
int primme_constant_info(const char* label_name, int *value);

#ifdef __cplusplus
}
#endif

#endif /* PRIMME_EIGS_H */
