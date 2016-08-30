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
 * File: restart_private.h
 *
 * Purpose - Definitions used exclusively by restart.c
 *
 ******************************************************************************/

#ifndef RESTART_PRIVATE_H
#define RESTART_PRIVATE_H

#define INSERT_SUBMATRIX_FAILURE -1
#define RESTART_H_FAILURE        -2
#define NUM_DSYEV_FAILURE        -3
#define UDUDECOMPOSE_FAILURE     -4
#define PSEUDOLOCK_FAILURE       -5

static int restart_soft_locking_zprimme(int *restartSize, __PRIMME_COMPLEX_DOUBLE__ *V,
       __PRIMME_COMPLEX_DOUBLE__ *W, int nLocal, int basisSize, int ldV, __PRIMME_COMPLEX_DOUBLE__ **X,
       __PRIMME_COMPLEX_DOUBLE__ **R, __PRIMME_COMPLEX_DOUBLE__ *hVecs, int ldhVecs, int *restartPerm,
       double *hVals, int *flags, int *iev, int *ievSize, double *blockNorms,
       __PRIMME_COMPLEX_DOUBLE__ *evecs, double *evals, double *resNorms, __PRIMME_COMPLEX_DOUBLE__ *evecsHat,
       int ldevecsHat, __PRIMME_COMPLEX_DOUBLE__ *M, int ldM, int *numConverged,
       int *numConvergedStored, int numPrevRetained, int *indexOfPreviousVecs,
       int *hVecsPerm, int reset, double machEps, __PRIMME_COMPLEX_DOUBLE__ *rwork, int rworkSize,
       int *iwork, primme_params *primme);

static int restart_projection_zprimme(__PRIMME_COMPLEX_DOUBLE__ *V, int ldV, __PRIMME_COMPLEX_DOUBLE__ *W,
      int ldW, __PRIMME_COMPLEX_DOUBLE__ *H, int ldH, __PRIMME_COMPLEX_DOUBLE__ *Q, int nLocal, int ldQ,
      __PRIMME_COMPLEX_DOUBLE__ *R, int ldR, __PRIMME_COMPLEX_DOUBLE__ *QtV, int ldQtV, __PRIMME_COMPLEX_DOUBLE__ *hU, int ldhU,
      int newldhU, int indexOfPreviousVecsBeforeRestart, __PRIMME_COMPLEX_DOUBLE__ *hVecs,
      int ldhVecs, int newldhVecs, double *hVals, double *hSVals,
      int *restartPerm, int *hVecsPerm, int restartSize, int basisSize,
      int numPrevRetained, int indexOfPreviousVecs, __PRIMME_COMPLEX_DOUBLE__ *evecs,
      int *evecsSize, int ldevecs, __PRIMME_COMPLEX_DOUBLE__ *evecsHat, int ldevecsHat,
      __PRIMME_COMPLEX_DOUBLE__ *M, int ldM, __PRIMME_COMPLEX_DOUBLE__ *UDU, int ldUDU, int *ipivot,
      int *targetShiftIndex, int numConverged, int *numArbitraryVecs,
      __PRIMME_COMPLEX_DOUBLE__ *hVecsRot, int ldhVecsRot,
      int rworkSize, __PRIMME_COMPLEX_DOUBLE__ *rwork, int *iwork,
      double machEps, primme_params *primme);

static int dtr_zprimme(int numLocked, __PRIMME_COMPLEX_DOUBLE__ *hVecs, double *hVals, int *flags, 
  int basisSize, int numFree, int *iev, __PRIMME_COMPLEX_DOUBLE__ *rwork, primme_params *primme);

static int restart_RR(__PRIMME_COMPLEX_DOUBLE__ *H, int ldH, __PRIMME_COMPLEX_DOUBLE__ *hVecs, int ldhVecs,
      int newldhVecs, double *hVals, int restartSize, int basisSize, int numLocked,
      int numPrevRetained, int indexOfPreviousVecs, int *hVecsPerm,
      double machEps, int rworkSize, __PRIMME_COMPLEX_DOUBLE__ *rwork, int *iwork,
      primme_params *primme);

static int restart_refined(__PRIMME_COMPLEX_DOUBLE__ *V, int ldV, __PRIMME_COMPLEX_DOUBLE__ *W, int ldW, __PRIMME_COMPLEX_DOUBLE__ *H,
      int ldH, __PRIMME_COMPLEX_DOUBLE__ *Q, int nLocal, int ldQ, __PRIMME_COMPLEX_DOUBLE__ *R, int ldR,
      __PRIMME_COMPLEX_DOUBLE__ *hU, int ldhU, int newldhU, int indexOfPreviousVecsBeforeRestart,
      __PRIMME_COMPLEX_DOUBLE__ *hVecs, int ldhVecs, int newldhVecs, double *hVals,
      double *hSVals, int *restartPerm, int *hVecsPerm, int restartSize,
      int basisSize, int numPrevRetained, int indexOfPreviousVecs,
      int *targetShiftIndex, int numConverged, int *numArbitraryVecs,
      __PRIMME_COMPLEX_DOUBLE__ *hVecsRot, int ldhVecsRot, int rworkSize,
      __PRIMME_COMPLEX_DOUBLE__ *rwork, int *iwork, double machEps, primme_params *primme);

static int restart_harmonic(__PRIMME_COMPLEX_DOUBLE__ *V, int ldV, __PRIMME_COMPLEX_DOUBLE__ *W, int ldW, __PRIMME_COMPLEX_DOUBLE__ *H,
   int ldH, __PRIMME_COMPLEX_DOUBLE__ *Q, int nLocal, int ldQ, __PRIMME_COMPLEX_DOUBLE__ *R, int ldR, __PRIMME_COMPLEX_DOUBLE__ *QtV,
   int ldQtV, __PRIMME_COMPLEX_DOUBLE__ *hU, int ldhU, int newldhU, __PRIMME_COMPLEX_DOUBLE__ *hVecs, int ldhVecs,
   int newldhVecs, double *hVals, double *hSVals, int *restartPerm, int *hVecsPerm,
   int restartSize, int basisSize, int numPrevRetained, int indexOfPreviousVecs,
   int *targetShiftIndex, int numConverged, int *numArbitraryVecs, __PRIMME_COMPLEX_DOUBLE__ *hVecsRot,
   int ldhVecsRot, int rworkSize, __PRIMME_COMPLEX_DOUBLE__ *rwork, int *iwork,
   double machEps, primme_params *primme);

static int ortho_coefficient_vectors_zprimme(__PRIMME_COMPLEX_DOUBLE__ *hVecs, int basisSize,
      int ldhVecs, int indexOfPreviousVecs, __PRIMME_COMPLEX_DOUBLE__ *hU, int ldhU, __PRIMME_COMPLEX_DOUBLE__ *R,
      int ldR, int *numPrevRetained, __PRIMME_COMPLEX_DOUBLE__ *outR, int ldoutR, double machEps,
      __PRIMME_COMPLEX_DOUBLE__ *rwork, int rworkSize, primme_params *primme);

#endif /* RESTART_PRIVATE_H */
