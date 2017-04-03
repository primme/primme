C*******************************************************************************
C  Copyright (c) 2017, College of William & Mary                                   
C  All rights reserved.                                                            
C                                                                                  
C  Redistribution and use in source and binary forms, with or without
C  modification, are permitted provided that the following conditions are met:     
C      * Redistributions of source code must retain the above copyright
C        notice, this list of conditions and the following disclaimer.             
C      * Redistributions in binary form must reproduce the above copyright         
C        notice, this list of conditions and the following disclaimer in the       
C        documentation and/or other materials provided with the distribution.      
C      * Neither the name of the College of William & Mary nor the
C        names of its contributors may be used to endorse or promote products      
C        derived from this software without specific prior written permission.     
C                                                                                  
C  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
C  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
C  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          
C  DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY       
C  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      
C  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
C  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
C  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
C  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C  
C  PRIMME: https://github.com/primme/primme
C  Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
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
     : PRIMME_STEEPEST_DESCENT,
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
     : PRIMME_STEEPEST_DESCENT = 13,
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
     : PRIMME_stats_numGlobalSum,
     : PRIMME_stats_volumeGlobalSum,
     : PRIMME_stats_numOrthoInnerProds,
     : PRIMME_stats_elapsedTime,
     : PRIMME_stats_timeMatvec,
     : PRIMME_stats_timePrecond,
     : PRIMME_stats_timeOrtho,
     : PRIMME_stats_timeGlobalSum,
     : PRIMME_stats_estimateMinEVal,
     : PRIMME_stats_estimateMaxEVal,
     : PRIMME_stats_estimateLargestSVal,
     : PRIMME_stats_maxConvTol,
     : PRIMME_dynamicMethodSwitch,
     : PRIMME_massMatrixMatvec,
     : PRIMME_convTestFun,
     : PRIMME_ldevecs,
     : PRIMME_ldOPs,
     : PRIMME_monitorFun,
     : PRIMME_monitor

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
     : PRIMME_stats_numGlobalSum =  471,
     : PRIMME_stats_volumeGlobalSum =  472,
     : PRIMME_stats_numOrthoInnerProds =  473,
     : PRIMME_stats_elapsedTime = 48,
     : PRIMME_stats_timeMatvec =  4801,
     : PRIMME_stats_timePrecond =  4802,
     : PRIMME_stats_timeOrtho =  4803,
     : PRIMME_stats_timeGlobalSum =  4804,
     : PRIMME_stats_estimateMinEVal = 481,
     : PRIMME_stats_estimateMaxEVal = 482,
     : PRIMME_stats_estimateLargestSVal = 483,
     : PRIMME_stats_maxConvTol = 484,
     : PRIMME_dynamicMethodSwitch = 49,
     : PRIMME_massMatrixMatvec = 50,
     : PRIMME_convTestFun = 51,
     : PRIMME_ldevecs = 52,
     : PRIMME_ldOPs = 53,
     : PRIMME_monitorFun = 54,
     : PRIMME_monitor = 55
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
     : primme_adaptive,
     : primme_event_outer_iteration,
     : primme_event_inner_iteration,
     : primme_event_restart,
     : primme_event_reset,
     : primme_event_converged,
     : primme_event_locked

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
     : primme_adaptive = 3,
     : primme_event_outer_iteration = 0,
     : primme_event_inner_iteration = 1,
     : primme_event_restart = 2,
     : primme_event_reset = 3,
     : primme_event_converged = 4,
     : primme_event_locked = 5
     : )
