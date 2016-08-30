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

static int restart_soft_locking_zprimme(int *restartSize, complex double *V,
       complex double *W, int nLocal, int basisSize, int ldV, complex double **X,
       complex double **R, complex double *hVecs, int ldhVecs, int *restartPerm,
       double *hVals, int *flags, int *iev, int *ievSize, double *blockNorms,
       complex double *evecs, double *evals, double *resNorms, complex double *evecsHat,
       int ldevecsHat, complex double *M, int ldM, int *numConverged,
       int *numConvergedStored, int numPrevRetained, int *indexOfPreviousVecs,
       int *hVecsPerm, int reset, double machEps, complex double *rwork, int rworkSize,
       int *iwork, primme_params *primme);

static int restart_projection_zprimme(complex double *V, int ldV, complex double *W,
      int ldW, complex double *H, int ldH, complex double *Q, int nLocal, int ldQ,
      complex double *R, int ldR, complex double *QtV, int ldQtV, complex double *hU, int ldhU,
      int newldhU, int indexOfPreviousVecsBeforeRestart, complex double *hVecs,
      int ldhVecs, int newldhVecs, double *hVals, double *hSVals,
      int *restartPerm, int *hVecsPerm, int restartSize, int basisSize,
      int numPrevRetained, int indexOfPreviousVecs, complex double *evecs,
      int *evecsSize, int ldevecs, complex double *evecsHat, int ldevecsHat,
      complex double *M, int ldM, complex double *UDU, int ldUDU, int *ipivot,
      int *targetShiftIndex, int numConverged, int *numArbitraryVecs,
      complex double *hVecsRot, int ldhVecsRot,
      int rworkSize, complex double *rwork, int *iwork,
      double machEps, primme_params *primme);

static int dtr_zprimme(int numLocked, complex double *hVecs, double *hVals, int *flags, 
  int basisSize, int numFree, int *iev, complex double *rwork, primme_params *primme);

static int restart_RR(complex double *H, int ldH, complex double *hVecs, int ldhVecs,
      int newldhVecs, double *hVals, int restartSize, int basisSize, int numLocked,
      int numPrevRetained, int indexOfPreviousVecs, int *hVecsPerm,
      double machEps, int rworkSize, complex double *rwork, int *iwork,
      primme_params *primme);

static int restart_refined(complex double *V, int ldV, complex double *W, int ldW, complex double *H,
      int ldH, complex double *Q, int nLocal, int ldQ, complex double *R, int ldR,
      complex double *hU, int ldhU, int newldhU, int indexOfPreviousVecsBeforeRestart,
      complex double *hVecs, int ldhVecs, int newldhVecs, double *hVals,
      double *hSVals, int *restartPerm, int *hVecsPerm, int restartSize,
      int basisSize, int numPrevRetained, int indexOfPreviousVecs,
      int *targetShiftIndex, int numConverged, int *numArbitraryVecs,
      complex double *hVecsRot, int ldhVecsRot, int rworkSize,
      complex double *rwork, int *iwork, double machEps, primme_params *primme);

static int restart_harmonic(complex double *V, int ldV, complex double *W, int ldW, complex double *H,
   int ldH, complex double *Q, int nLocal, int ldQ, complex double *R, int ldR, complex double *QtV,
   int ldQtV, complex double *hU, int ldhU, int newldhU, complex double *hVecs, int ldhVecs,
   int newldhVecs, double *hVals, double *hSVals, int *restartPerm, int *hVecsPerm,
   int restartSize, int basisSize, int numPrevRetained, int indexOfPreviousVecs,
   int *targetShiftIndex, int numConverged, int *numArbitraryVecs, complex double *hVecsRot,
   int ldhVecsRot, int rworkSize, complex double *rwork, int *iwork,
   double machEps, primme_params *primme);

static int ortho_coefficient_vectors_zprimme(complex double *hVecs, int basisSize,
      int ldhVecs, int indexOfPreviousVecs, complex double *hU, int ldhU, complex double *R,
      int ldR, int *numPrevRetained, complex double *outR, int ldoutR, double machEps,
      complex double *rwork, int rworkSize, primme_params *primme);

#endif /* RESTART_PRIVATE_H */
