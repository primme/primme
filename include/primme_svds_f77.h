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
C  File: primme_svds_f77.h
C  
C  Purpose - Main header with the PRIMME SVDS F77 interface functions.
C  
C******************************************************************************

C-------------------------------------------------------
C     Defining easy to remember labels for setting the 
C     method in primme_svds_set_method from Fortran
C-------------------------------------------------------
      integer
     : primme_svds_default,
     : primme_svds_hybrid,
     : primme_svds_normalequations,
     : primme_svds_augmented

      parameter(
     : primme_svds_default = 0,
     : primme_svds_hybrid = 1,
     : primme_svds_normalequations = 2,
     : primme_svds_augmented = 3
     :)

C-------------------------------------------------------
C     Defining easy to remember labels for setting the 
C     members of the primme_svds structure from Fortran
C-------------------------------------------------------
 
      INTEGER
     : PRIMME_SVDS_primme,
     : PRIMME_SVDS_primmeStage2,
     : PRIMME_SVDS_m,
     : PRIMME_SVDS_n,
     : PRIMME_SVDS_matrixMatvec ,
     : PRIMME_SVDS_applyPreconditioner,
     : PRIMME_SVDS_numProcs,
     : PRIMME_SVDS_procID,
     : PRIMME_SVDS_mLocal,
     : PRIMME_SVDS_nLocal,
     : PRIMME_SVDS_commInfo,
     : PRIMME_SVDS_globalSumReal,
     : PRIMME_SVDS_numSvals,
     : PRIMME_SVDS_target,
     : PRIMME_SVDS_numTargetShifts,
     : PRIMME_SVDS_targetShifts,
     : PRIMME_SVDS_method,
     : PRIMME_SVDS_methodStage2,
     : PRIMME_SVDS_intWorkSize,
     : PRIMME_SVDS_realWorkSize,
     : PRIMME_SVDS_intWork,
     : PRIMME_SVDS_realWork,
     : PRIMME_SVDS_matrix,
     : PRIMME_SVDS_preconditioner,
     : PRIMME_SVDS_locking,
     : PRIMME_SVDS_numOrthoConst,
     : PRIMME_SVDS_aNorm,
     : PRIMME_SVDS_eps,
     : PRIMME_SVDS_precondition,
     : PRIMME_SVDS_initSize,
     : PRIMME_SVDS_maxBasisSize,
     : PRIMME_SVDS_maxBlockSize,
     : PRIMME_SVDS_maxMatvecs,
     : PRIMME_SVDS_iseed,
     : PRIMME_SVDS_printLevel,
     : PRIMME_SVDS_outputFile,
     : PRIMME_SVDS_stats_numOuterIterations, 
     : PRIMME_SVDS_stats_numRestarts,
     : PRIMME_SVDS_stats_numMatvecs,
     : PRIMME_SVDS_stats_numPreconds,
     : PRIMME_SVDS_stats_elapsedTime

      parameter(
     : PRIMME_SVDS_primme = 0,
     : PRIMME_SVDS_primmeStage2 = 1,
     : PRIMME_SVDS_m = 2,
     : PRIMME_SVDS_n = 3,
     : PRIMME_SVDS_matrixMatvec = 4,
     : PRIMME_SVDS_applyPreconditioner = 5,
     : PRIMME_SVDS_numProcs = 6,
     : PRIMME_SVDS_procID = 7,
     : PRIMME_SVDS_mLocal = 8,
     : PRIMME_SVDS_nLocal = 9,
     : PRIMME_SVDS_commInfo = 10,
     : PRIMME_SVDS_globalSumReal = 11,
     : PRIMME_SVDS_numSvals = 12,
     : PRIMME_SVDS_target = 13,
     : PRIMME_SVDS_numTargetShifts = 14,
     : PRIMME_SVDS_targetShifts = 15,
     : PRIMME_SVDS_method = 16,
     : PRIMME_SVDS_methodStage2 = 17,
     : PRIMME_SVDS_intWorkSize = 18,
     : PRIMME_SVDS_realWorkSize = 19,
     : PRIMME_SVDS_intWork = 20,
     : PRIMME_SVDS_realWork = 21,
     : PRIMME_SVDS_matrix = 22,
     : PRIMME_SVDS_preconditioner = 23,
     : PRIMME_SVDS_locking = 24,
     : PRIMME_SVDS_numOrthoConst = 25,
     : PRIMME_SVDS_aNorm = 26,
     : PRIMME_SVDS_eps = 27,
     : PRIMME_SVDS_precondition = 28,
     : PRIMME_SVDS_initSize = 29,
     : PRIMME_SVDS_maxBasisSize = 30,
     : PRIMME_SVDS_maxBlockSize = 31,
     : PRIMME_SVDS_maxMatvecs = 32,
     : PRIMME_SVDS_iseed = 33,
     : PRIMME_SVDS_printLevel = 34,
     : PRIMME_SVDS_outputFile = 35,
     : PRIMME_SVDS_stats_numOuterIterations = 36,
     : PRIMME_SVDS_stats_numRestarts = 37,
     : PRIMME_SVDS_stats_numMatvecs = 38,
     : PRIMME_SVDS_stats_numPreconds = 39,
     : PRIMME_SVDS_stats_elapsedTime = 40
     :)

C-------------------------------------------------------
C    Defining easy to remember labels for setting the 
C    enum members for targeting and operator
C-------------------------------------------------------

      integer 
     : primme_svds_largest,
     : primme_svds_smallest,
     : primme_svds_closest_abs,
     : primme_svds_op_none,
     : primme_svds_op_AtA,
     : primme_svds_op_AAt,
     : primme_svds_op_augmented

      parameter(
     : primme_svds_largest = 0,
     : primme_svds_smallest = 1,
     : primme_svds_closest_abs = 2,
     : primme_svds_op_none = 0,
     : primme_svds_op_AtA = 1,
     : primme_svds_op_AAt = 2,
     : primme_svds_op_augmented = 3
     :)
