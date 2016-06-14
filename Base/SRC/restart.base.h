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
 * File: restart.h
 *
 * Purpose - Contains prototypes for functions related to restarting.
 *
 ******************************************************************************/

#ifndef RESTART_H
#define RESTART_H

void reset_flags_@(pre)primme(int *flags, int first, int last);

int restart_@(pre)primme(@(type) *V, @(type) *W, int nLocal, int basisSize, int ldV,
   double *hVals, double *hSVals, int *flags, int *iev, int *ievSize,
   double *blockNorms, @(type) *evecs, int *evecsPerm, double *evals, double *resNorms,
   @(type) *evecsHat, int ldevecsHat, @(type) *M, int ldM, @(type) *UDU,
   int ldUDU, int *ipivot, int *numConverged, int *numLocked,
   int *numConvergedStored, @(type) *previousHVecs, int *numPrevRetained,
   int ldpreviousHVecs, int numGuesses, double *prevRitzVals, int *numPrevRitzVals,
   @(type) *H, int ldH, @(type) *Q, int ldQ, @(type) *R, int ldR, @(type)* QtV, int ldQtV,
   @(type) *hU, int ldhU, int newldhU, @(type) *hVecs, int ldhVecs, int newldhVecs,
   int *restartSizeOutput, int *targetShiftIndex, int numArbitraryVecs,
   int *restartsSinceReset, int *reset, double machEps,
   @(type) *rwork, int rworkSize, int *iwork, primme_params *primme);

int ortho_coefficient_vectors_@(pre)primme(@(type) *hVecs, int basisSize, int ldhVecs,
   int indexOfPreviousVecs, int newBasisSize, int *perm, @(type) *hU, int ldhU,
   @(type) *R, int ldR, int numPrevRetained, double machEps, int *iwork,
   @(type) *rwork, int rworkSize, primme_params *primme);

int Num_reset_update_VWXR_@(pre)primme(@(type) *V, @(type) *W, int mV, int nV, int ldV,
   @(type) *h, int nh, int ldh, double *hVals,
   @(type) *X0, int nX0b, int nX0e, int ldX0,
   @(type) *X1, int nX1b, int nX1e, int ldX1,
   @(type) *evecs, int evecsSize, int nX2b, int nX2e, int ldevecs,
   @(type) *Wo, int nWob, int nWoe, int ldWo,
   @(type) *R, int nRb, int nRe, int ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   int reset, double machEps, @(type) *rwork, int lrwork, primme_params *primme);

#endif
