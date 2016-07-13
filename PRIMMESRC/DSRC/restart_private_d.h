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

static int restart_soft_locking_dprimme(int *restartSize, double *V,
       double *W, int nLocal, int basisSize, int ldV, double **X,
       double **R, double *hVecs, int ldhVecs, int *restartPerm,
       double *hVals, int *flags, int *iev, int *ievSize, double *blockNorms,
       double *evecs, double *evals, double *resNorms, double *evecsHat,
       int ldevecsHat, double *M, int ldM, int *numConverged,
       int *numConvergedStored, int numPrevRetained, int *indexOfPreviousVecs,
       int *hVecsPerm, int reset, double machEps, double *rwork, int rworkSize,
       int *iwork, primme_params *primme);

static int restart_projection_dprimme(double *V, int ldV, double *W,
      int ldW, double *H, int ldH, double *Q, int nLocal, int ldQ,
      double *R, int ldR, double *QtV, int ldQtV, double *hU, int ldhU,
      int newldhU, double *hVecs, int ldhVecs, int newldhVecs, double *hVals,
      double *hSVals, int *restartPerm, int *hVecsPerm, int restartSize,
      int basisSize, int numPrevRetained, int indexOfPreviousVecs,
      double *evecs, int *evecsSize, int ldevecs, double *evecsHat,
      int ldevecsHat, double *M, int ldM, double *UDU, int ldUDU,
      int *ipivot, int *targetShiftIndex, int numConverged,
      int numArbitraryVecs,  double *hVecsRot, int ldhVecsRot,
      double *previousHU, int ldpreviousHU, double *prevhSvals,
      int numRecentlyLocked, int rworkSize, double *rwork, int *iwork,
      double machEps, primme_params *primme);

static int dtr_dprimme(int numLocked, double *hVecs, double *hVals, int *flags, 
  int basisSize, int numFree, int *iev, double *rwork, primme_params *primme);

static int restart_RR(double *H, int ldH, double *hVecs, int ldhVecs,
      int newldhVecs, double *hVals, int restartSize, int basisSize, int numLocked,
      int numPrevRetained, int indexOfPreviousVecs, int *hVecsPerm,
      double machEps, int rworkSize, double *rwork, int *iwork,
      primme_params *primme);

static int restart_refined(double *V, int ldV, double *W, int ldW, double *H,
   int ldH, double *Q, int nLocal, int ldQ, double *R, int ldR,
   double *hU, int ldhU, int newldhU, double *hVecs, int ldhVecs,
   int newldhVecs, double *hVals, double *hSVals, int *restartPerm, int *hVecsPerm,
   int restartSize, int basisSize, int numPrevRetained, int indexOfPreviousVecs,
   int *targetShiftIndex, int numConverged, int numArbitraryVecs, double *hVecsRot,
   int ldhVecsRot, double *previousHU, int ldpreviousHU, double *prevhSvals,
   int numRecentlyLocked, int rworkSize, double *rwork, int *iwork, double machEps, 
   primme_params *primme);

static int restart_harmonic(double *V, int ldV, double *W, int ldW, double *H,
   int ldH, double *Q, int nLocal, int ldQ, double *R, int ldR, double *QtV,
   int ldQtV, double *hU, int ldhU, int newldhU, double *hVecs, int ldhVecs,
   int newldhVecs, double *hVals, double *hSVals, int *restartPerm, int *hVecsPerm,
   int restartSize, int basisSize, int numPrevRetained, int indexOfPreviousVecs,
   int *targetShiftIndex, int numConverged, int numArbitraryVecs, double *hVecsRot,
   int ldhVecsRot, double *previousHU, int ldpreviousHU, double *prevhSvals,
   int numRecentlyLocked, int rworkSize, double *rwork, int *iwork, double machEps, 
   primme_params *primme);

static int ortho_coefficient_vectors_dprimme(double *hVecs, int basisSize,
      int ldhVecs, int indexOfPreviousVecs, double *hU, int ldhU, double *R,
      int ldR, int *numPrevRetained, double *prevhSvals, double *previousHU,
      int ldpreviousHU, double machEps, double *rwork, int rworkSize,
      primme_params *primme);

#endif /* RESTART_PRIVATE_H */
