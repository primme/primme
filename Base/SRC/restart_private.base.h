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

static int restart_soft_locking_@(pre)primme(int *restartSize, @(type) *V, @(type) *W, int nLocal,
   int basisSize, int ldV, @(type) **X, @(type) **R, @(type) *hVecs, int ldhVecs,
   int *restartPerm, double *hVals, int *flags, int *iev, int *ievSize,
   double *blockNorms, @(type) *evecs, double *evals, double *resNorms,
   @(type) *evecsHat, int ldevecsHat, @(type) *M, int ldM, int *numConverged,
   int *numConvergedStored, @(type) *previousHVecs, int *numPrevRetained,
   int ldpreviousHVecs, int *indexOfPreviousVecs, int *hVecsPerm, double machEps,
   @(type) *rwork, int rworkSize, int *iwork, primme_params *primme);

static int restart_projection_@(pre)primme(@(type) *V, int ldV, @(type) *W, int ldW,
   @(type) *H, int ldH, @(type) *Q, int nLocal, int ldQ, @(type) *R, int ldR,
   @(type) *hU, int ldhU, int newldhU, @(type) *hVecs, int ldhVecs, int newldhVecs,
   double *hVals, double *hSVals, int *restartPerm, int *hVecsPerm,
   int restartSize, int basisSize, int numPrevRetained,
   int indexOfPreviousVecs, @(type) *evecs, int *evecsSize,
   int ldevecs, @(type) *evecsHat, int ldevecsHat, @(type) *M, int ldM, @(type) *UDU,
   int ldUDU, int *ipivot, int *targetShiftIndex, int numConverged,
   int rworkSize, @(type) *rwork, int *iwork, double machEps, primme_params *primme);

static int dtr_@(pre)primme(int numLocked, @(type) *hVecs, double *hVals, int *flags, 
  int basisSize, int numFree, int *iev, @(type) *rwork, primme_params *primme);

static int restart_RR(@(type) *H, int ldH, @(type) *hVecs, int ldhVecs,
   int newldhVecs, double *hVals, int restartSize, int basisSize, int numLocked,
   int numPrevRetained, int indexOfPreviousVecs, int *hVecsPerm,
   int rworkSize, @(type) *rwork, int *iwork, primme_params *primme);

static int restart_ref(@(type) *V, int ldV, @(type) *W, int ldW, @(type) *H,
   int ldH, @(type) *Q, int nLocal, int ldQ, @(type) *R, int ldR, @(type) *hU,
   int ldhU, int newldhU, @(type) *hVecs, int ldhVecs, int newldhVecs,
   double *hVals, double *hSVals, int *restartPerm, int *hVecsPerm, int restartSize,
   int basisSize, int numPrevRetained, int indexOfPreviousVecs,
   int *targetShiftIndex, int numConverged, int rworkSize,
   @(type) *rwork, int *iwork, double machEps, primme_params *primme);

static int compute_submatrix(@(type) *X, int nX, int ldX, 
   @(type) *H, int nH, int ldH, @(type) *R, int ldR,
   @(type) *rwork, int lrwork);

#endif /* RESTART_PRIVATE_H */
