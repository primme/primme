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

#include <stdio.h>

#ifdef __cplusplus
#  include <complex>
#  define __PRIMME_COMPLEX_DOUBLE__ std::complex<double>
extern "C" {
#else
#  include <complex.h>
#  define __PRIMME_COMPLEX_DOUBLE__ double complex
#endif

#if !defined(PRIMME_INT_SIZE) || PRIMME_INT_SIZE == 64
#  include <stdint.h>
#  include <inttypes.h>
#  define PRIMME_INT int64_t
#  define PRIMME_INT_P PRId64
#  define PRIMME_INT_MAX INT64_MAX
#elif PRIMME_INT_SIZE == 0
#  include <limits.h>
#  define PRIMME_INT int
#  define PRIMME_INT_P "d"
#  define PRIMME_INT_MAX INT_MAX
#elif PRIMME_INT == 32
#  include <stdint.h>
#  include <inttypes.h>
#  define PRIMME_INT int32_t
#  define PRIMME_INT_P PRId32
#  define PRIMME_INT_MAX INT32_MAX
#else
#  define PRIMME_INT PRIMME_INT_SIZE
#  define PRIMME_INT_P "d"
#  define PRIMME_INT_MAX INT_MAX
#endif

#define PRIMME_MAX_NAME_LENGTH 128

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


typedef struct primme_stats {
   PRIMME_INT numOuterIterations;
   PRIMME_INT numRestarts;
   PRIMME_INT numMatvecs;
   PRIMME_INT numPreconds;
   double elapsedTime; 
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
   PRIMME_INT nLocal;
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

   struct projection_params projectionParams; 
   struct restarting_params restartingParams;
   struct correction_params correctionParams;
   struct primme_stats stats;

   void (*convTestFun)(double *eval, void *evec, double *rNorm, int *isconv, 
         struct primme_params *primme);
} primme_params;
/*---------------------------------------------------------------------------*/

typedef enum {
   DEFAULT_METHOD,
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
int zprimme(double *evals, __PRIMME_COMPLEX_DOUBLE__ *evecs, double *resNorms, 
      primme_params *primme);
void primme_initialize(primme_params *primme);
int  primme_set_method(primme_preset_method method, primme_params *params);
void primme_display_params(primme_params primme);
void *primme_valloc(size_t byteSize, const char *target);
void *primme_calloc(size_t nelem, size_t elsize, const char *target);
void primme_Free(primme_params *primme);

#ifdef __cplusplus
}
#endif

#endif /* PRIMME_H */
