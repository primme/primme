C*******************************************************************************
C    PRIMME PReconditioned Iterative MultiMethod Eigensolver
C    Copyright (C) 2015 College of William & Mary,
C    James R. McCombs, Eloy Romero Alcalde, Andreas Stathopoulos, Lingfei Wu
C 
C    This file is part of PRIMME.
C 
C    PRIMME is free software; you can redistribute it and/or
C    modify it under the terms of the GNU Lesser General Public
C    License as published by the Free Software Foundation; either
C    version 2.1 of the License, or (at your option) any later version.
C 
C    PRIMME is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C    Lesser General Public License for more details.
C 
C    You should have received a copy of the GNU Lesser General Public
C    License along with this library; if not, write to the Free Software
C    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
C 
C*******************************************************************************
C  File: primme_eigs_f77.h
C  
C  Purpose - Main header with the PRIMME EIGS F77 interface functions.
C  
C******************************************************************************

C-------------------------------------------------------
C     Defining easy to remember labels for setting the 
C     method in primme_set_method from Fortran
C-------------------------------------------------------
      integer
     : PRIMME_DEFAULT_METHOD,
     : PRIMME_DYNAMIC,
     : PRIMME_DEFAULT_MIN_TIME,
     : PRIMME_DEFAULT_MIN_MATVECS,
     : PRIMME_Arnoldi,
     : PRIMME_GD,
     : PRIMME_GD_plusK,
     : PRIMME_GD_Olsen_plusK,
     : PRIMME_JD_Olsen_plusK,
     : PRIMME_RQI,
     : PRIMME_JDQR,
     : PRIMME_JDQMR,
     : PRIMME_JDQMR_ETol,
     : PRIMME_SUBSPACE_ITERATION,
     : PRIMME_LOBPCG_OrthoBasis,
     : PRIMME_LOBPCG_OrthoBasis_Window

      parameter(
     : PRIMME_DEFAULT_METHOD = 0,
     : PRIMME_DYNAMIC = 1,
     : PRIMME_DEFAULT_MIN_TIME = 2,
     : PRIMME_DEFAULT_MIN_MATVECS = 3,
     : PRIMME_Arnoldi = 4,
     : PRIMME_GD = 5,
     : PRIMME_GD_plusK = 6,
     : PRIMME_GD_Olsen_plusK = 7,
     : PRIMME_JD_Olsen_plusK = 8,
     : PRIMME_RQI = 9,
     : PRIMME_JDQR = 10,
     : PRIMME_JDQMR = 11,
     : PRIMME_JDQMR_ETol = 12,
     : PRIMME_SUBSPACE_ITERATION = 13,
     : PRIMME_LOBPCG_OrthoBasis = 14,
     : PRIMME_LOBPCG_OrthoBasis_Window = 15
     :)

C-------------------------------------------------------
C     Defining easy to remember labels for setting the 
C     members of the primme structure from Fortran
C-------------------------------------------------------
      INTEGER
     : PRIMME_n,
     : PRIMME_matrixMatvec,
     : PRIMME_applyPreconditioner,
     : PRIMME_numProcs,
     : PRIMME_procID,
     : PRIMME_commInfo,
     : PRIMME_nLocal,
     : PRIMME_globalSumReal,
     : PRIMME_numEvals,
     : PRIMME_target,
     : PRIMME_numTargetShifts,
     : PRIMME_targetShifts,
     : PRIMME_locking,
     : PRIMME_initSize,
     : PRIMME_numOrthoConst,
     : PRIMME_maxBasisSize,
     : PRIMME_minRestartSize,
     : PRIMME_maxBlockSize,
     : PRIMME_maxMatvecs,
     : PRIMME_maxOuterIterations,
     : PRIMME_intWorkSize,
     : PRIMME_realWorkSize,
     : PRIMME_iseed,
     : PRIMME_intWork,
     : PRIMME_realWork,
     : PRIMME_aNorm,
     : PRIMME_eps,
     : PRIMME_printLevel,
     : PRIMME_outputFile,
     : PRIMME_matrix,
     : PRIMME_preconditioner,
     : PRIMME_initBasisMode,
     : PRIMME_projectionParams_projection,
     : PRIMME_restartingParams_scheme,
     : PRIMME_restartingParams_maxPrevRetain,
     : PRIMME_correctionParams_precondition,
     : PRIMME_correctionParams_robustShifts,
     : PRIMME_correctionParams_maxInnerIterations,
     : PRIMME_correctionParams_projectors_LeftQ,
     : PRIMME_correctionParams_projectors_LeftX,
     : PRIMME_correctionParams_projectors_RightQ,
     : PRIMME_correctionParams_projectors_RightX,
     : PRIMME_correctionParams_projectors_SkewQ,
     : PRIMME_correctionParams_projectors_SkewX,
     : PRIMME_correctionParams_convTest,
     : PRIMME_correctionParams_relTolBase,
     : PRIMME_stats_numOuterIterations,
     : PRIMME_stats_numRestarts,
     : PRIMME_stats_numMatvecs,
     : PRIMME_stats_numPreconds,
     : PRIMME_stats_elapsedTime,
     : PRIMME_stats_estimateMinEVal,
     : PRIMME_stats_estimateMaxEVal,
     : PRIMME_stats_estimateLargestSVal,
     : PRIMME_stats_maxConvTol,
     : PRIMME_dynamicMethodSwitch,
     : PRIMME_massMatrixMatvec,
     : PRIMME_convTestFun

      parameter(
     : PRIMME_n = 0,
     : PRIMME_matrixMatvec = 1, 
     : PRIMME_applyPreconditioner = 2,
     : PRIMME_numProcs = 3,
     : PRIMME_procID = 4,
     : PRIMME_commInfo = 5,
     : PRIMME_nLocal = 6,
     : PRIMME_globalSumReal = 7,
     : PRIMME_numEvals = 8,
     : PRIMME_target = 9,
     : PRIMME_numTargetShifts = 10,
     : PRIMME_targetShifts = 11,
     : PRIMME_locking = 12,
     : PRIMME_initSize = 13,
     : PRIMME_numOrthoConst = 14,
     : PRIMME_maxBasisSize = 15,
     : PRIMME_minRestartSize = 16,
     : PRIMME_maxBlockSize = 17,
     : PRIMME_maxMatvecs = 18,
     : PRIMME_maxOuterIterations = 19,
     : PRIMME_intWorkSize = 20,
     : PRIMME_realWorkSize = 21,
     : PRIMME_iseed = 22,
     : PRIMME_intWork = 23,
     : PRIMME_realWork = 24,
     : PRIMME_aNorm = 25,
     : PRIMME_eps = 26,
     : PRIMME_printLevel = 27,
     : PRIMME_outputFile = 28,
     : PRIMME_matrix = 29,
     : PRIMME_preconditioner = 30,
     : PRIMME_initBasisMode = 301,
     : PRIMME_projectionParams_projection = 302,
     : PRIMME_restartingParams_scheme = 31,
     : PRIMME_restartingParams_maxPrevRetain = 32,
     : PRIMME_correctionParams_precondition = 33,
     : PRIMME_correctionParams_robustShifts = 34,
     : PRIMME_correctionParams_maxInnerIterations = 35,
     : PRIMME_correctionParams_projectors_LeftQ = 36,
     : PRIMME_correctionParams_projectors_LeftX = 37,
     : PRIMME_correctionParams_projectors_RightQ = 38,
     : PRIMME_correctionParams_projectors_RightX = 39,
     : PRIMME_correctionParams_projectors_SkewQ = 40,
     : PRIMME_correctionParams_projectors_SkewX = 41,
     : PRIMME_correctionParams_convTest = 42,
     : PRIMME_correctionParams_relTolBase = 43,
     : PRIMME_stats_numOuterIterations = 44,
     : PRIMME_stats_numRestarts = 45,
     : PRIMME_stats_numMatvecs = 46,
     : PRIMME_stats_numPreconds = 47,
     : PRIMME_stats_elapsedTime = 48,
     : PRIMME_stats_estimateMinEVal = 481,
     : PRIMME_stats_estimateMaxEVal = 482,
     : PRIMME_stats_estimateLargestSVal = 483,
     : PRIMME_stats_maxConvTol = 484,
     : PRIMME_dynamicMethodSwitch = 49,
     : PRIMME_massMatrixMatvec = 50,
     : PRIMME_convTestFun = 51
     : )

C-------------------------------------------------------
C    Defining easy to remember labels for setting the 
C    enum members for targeting, restarting and innertest
C-------------------------------------------------------

      integer 
     : primme_smallest,
     : primme_largest,
     : primme_closest_geq,
     : primme_closest_leq,
     : primme_closest_abs,
     : primme_largest_abs,
     : primme_proj_default,
     : primme_proj_RR,
     : primme_proj_harmonic,
     : primme_proj_refined,
     : primme_init_default,
     : primme_init_krylov,
     : primme_init_random,
     : primme_init_user,
     : primme_thick,
     : primme_dtr,
     : primme_full_LTolerance,
     : primme_decreasing_LTolerance,
     : primme_adaptive_ETolerance,
     : primme_adaptive

      parameter(
     : primme_smallest = 0,
     : primme_largest = 1,
     : primme_closest_geq = 2,
     : primme_closest_leq = 3,
     : primme_closest_abs = 4,
     : primme_largest_abs = 5,
     : primme_proj_default = 0,
     : primme_proj_RR = 1,
     : primme_proj_harmonic = 2,
     : primme_proj_refined = 3,
     : primme_init_default = 0,
     : primme_init_krylov = 1,
     : primme_init_random = 2,
     : primme_init_user = 3,
     : primme_thick = 0,
     : primme_dtr = 1,
     : primme_full_LTolerance = 0,
     : primme_decreasing_LTolerance = 1,
     : primme_adaptive_ETolerance = 2,
     : primme_adaptive = 3
     : )
