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
     : PRIMMEF77_SVDS_default,
     : PRIMMEF77_SVDS_hybrid,
     : PRIMMEF77_SVDS_normalequations,
     : PRIMMEF77_SVDS_augmented

      parameter(
     : PRIMMEF77_SVDS_default = 0,
     : PRIMMEF77_SVDS_hybrid = 1,
     : PRIMMEF77_SVDS_normalequations = 2,
     : PRIMMEF77_SVDS_augmented = 3
     :)

C-------------------------------------------------------
C     Defining easy to remember labels for setting the 
C     members of the primme_svds structure from Fortran
C-------------------------------------------------------
 
      INTEGER
     : PRIMMEF77_SVDS_primme,
     : PRIMMEF77_SVDS_primmeStage2,
     : PRIMMEF77_SVDS_m,
     : PRIMMEF77_SVDS_n,
     : PRIMMEF77_SVDS_matrixMatvec ,
     : PRIMMEF77_SVDS_applyPreconditioner,
     : PRIMMEF77_SVDS_numProcs,
     : PRIMMEF77_SVDS_procID,
     : PRIMMEF77_SVDS_mLocal,
     : PRIMMEF77_SVDS_nLocal,
     : PRIMMEF77_SVDS_commInfo,
     : PRIMMEF77_SVDS_globalSumDouble,
     : PRIMMEF77_SVDS_numSvals,
     : PRIMMEF77_SVDS_target,
     : PRIMMEF77_SVDS_numTargetShifts,
     : PRIMMEF77_SVDS_targetShifts,
     : PRIMMEF77_SVDS_method,
     : PRIMMEF77_SVDS_methodStage2,
     : PRIMMEF77_SVDS_intWorkSize,
     : PRIMMEF77_SVDS_realWorkSize,
     : PRIMMEF77_SVDS_intWork,
     : PRIMMEF77_SVDS_realWork,
     : PRIMMEF77_SVDS_matrix,
     : PRIMMEF77_SVDS_preconditioner,
     : PRIMMEF77_SVDS_locking,
     : PRIMMEF77_SVDS_numOrthoConst,
     : PRIMMEF77_SVDS_aNorm,
     : PRIMMEF77_SVDS_eps,
     : PRIMMEF77_SVDS_precondition,
     : PRIMMEF77_SVDS_initSize,
     : PRIMMEF77_SVDS_maxBasisSize,
     : PRIMMEF77_SVDS_maxBlockSize,
     : PRIMMEF77_SVDS_maxMatvecs,
     : PRIMMEF77_SVDS_iseed,
     : PRIMMEF77_SVDS_printLevel,
     : PRIMMEF77_SVDS_outputFile,
     : PRIMMEF77_SVDS_stats_numOuterIterations, 
     : PRIMMEF77_SVDS_stats_numRestarts,
     : PRIMMEF77_SVDS_stats_numMatvecs,
     : PRIMMEF77_SVDS_stats_numPreconds,
     : PRIMMEF77_SVDS_stats_elapsedTime

      parameter(
     : PRIMMEF77_SVDS_primme = 0,
     : PRIMMEF77_SVDS_primmeStage2 = 1,
     : PRIMMEF77_SVDS_m = 2,
     : PRIMMEF77_SVDS_n = 3,
     : PRIMMEF77_SVDS_matrixMatvec = 4,
     : PRIMMEF77_SVDS_applyPreconditioner = 5,
     : PRIMMEF77_SVDS_numProcs = 6,
     : PRIMMEF77_SVDS_procID = 7,
     : PRIMMEF77_SVDS_mLocal = 8,
     : PRIMMEF77_SVDS_nLocal = 9,
     : PRIMMEF77_SVDS_commInfo = 10,
     : PRIMMEF77_SVDS_globalSumDouble = 11,
     : PRIMMEF77_SVDS_numSvals = 12,
     : PRIMMEF77_SVDS_target = 13,
     : PRIMMEF77_SVDS_numTargetShifts = 14,
     : PRIMMEF77_SVDS_targetShifts = 15,
     : PRIMMEF77_SVDS_method = 16,
     : PRIMMEF77_SVDS_methodStage2 = 17,
     : PRIMMEF77_SVDS_intWorkSize = 18,
     : PRIMMEF77_SVDS_realWorkSize = 19,
     : PRIMMEF77_SVDS_intWork = 20,
     : PRIMMEF77_SVDS_realWork = 21,
     : PRIMMEF77_SVDS_matrix = 22,
     : PRIMMEF77_SVDS_preconditioner = 23,
     : PRIMMEF77_SVDS_locking = 24,
     : PRIMMEF77_SVDS_numOrthoConst = 25,
     : PRIMMEF77_SVDS_aNorm = 26,
     : PRIMMEF77_SVDS_eps = 27,
     : PRIMMEF77_SVDS_precondition = 28,
     : PRIMMEF77_SVDS_initSize = 29,
     : PRIMMEF77_SVDS_maxBasisSize = 30,
     : PRIMMEF77_SVDS_maxBlockSize = 31,
     : PRIMMEF77_SVDS_maxMatvecs = 32,
     : PRIMMEF77_SVDS_iseed = 33,
     : PRIMMEF77_SVDS_printLevel = 34,
     : PRIMMEF77_SVDS_outputFile = 35,
     : PRIMMEF77_SVDS_stats_numOuterIterations = 36,
     : PRIMMEF77_SVDS_stats_numRestarts = 37,
     : PRIMMEF77_SVDS_stats_numMatvecs = 38,
     : PRIMMEF77_SVDS_stats_numPreconds = 39,
     : PRIMMEF77_SVDS_stats_elapsedTime = 40
     :)

C-------------------------------------------------------
C    Defining easy to remember labels for setting the 
C    enum members for targeting and operator
C-------------------------------------------------------

      integer 
     : PRIMMEF77_SVDS_largest,
     : PRIMMEF77_SVDS_smallest,
     : PRIMMEF77_SVDS_closest_abs,
     : PRIMMEF77_SVDS_op_none,
     : PRIMMEF77_SVDS_op_AtA,
     : PRIMMEF77_SVDS_op_AAt,
     : PRIMMEF77_SVDS_op_augmented

      parameter(
     : PRIMMEF77_SVDS_largest = 0,
     : PRIMMEF77_SVDS_smallest = 1,
     : PRIMMEF77_SVDS_closest_abs = 2,
     : PRIMMEF77_SVDS_op_none = 0,
     : PRIMMEF77_SVDS_op_AtA = 1,
     : PRIMMEF77_SVDS_op_AAt = 2,
     : PRIMMEF77_SVDS_op_augmented = 3
     :)
