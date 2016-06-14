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

void reset_flags_zprimme(int *flags, int first, int last);

int restart_zprimme(Complex_Z *V, Complex_Z *W, int nLocal, int basisSize, int ldV,
   double *hVals, double *hSVals, int *flags, int *iev, int *ievSize,
   double *blockNorms, Complex_Z *evecs, int *evecsPerm, double *evals, double *resNorms,
   Complex_Z *evecsHat, int ldevecsHat, Complex_Z *M, int ldM, Complex_Z *UDU,
   int ldUDU, int *ipivot, int *numConverged, int *numLocked,
   int *numConvergedStored, Complex_Z *previousHVecs, int *numPrevRetained,
   int ldpreviousHVecs, int numGuesses, double *prevRitzVals, int *numPrevRitzVals,
   Complex_Z *H, int ldH, Complex_Z *Q, int ldQ, Complex_Z *R, int ldR, Complex_Z* QtV, int ldQtV,
   Complex_Z *hU, int ldhU, int newldhU, Complex_Z *hVecs, int ldhVecs, int newldhVecs,
   int *restartSizeOutput, int *targetShiftIndex, int numArbitraryVecs,
   int *restartsSinceReset, int *reset, double machEps,
   Complex_Z *rwork, int rworkSize, int *iwork, primme_params *primme);

int ortho_coefficient_vectors_zprimme(Complex_Z *hVecs, int basisSize, int ldhVecs,
   int indexOfPreviousVecs, int newBasisSize, int *perm, Complex_Z *hU, int ldhU,
   Complex_Z *R, int ldR, int numPrevRetained, double machEps, int *iwork,
   Complex_Z *rwork, int rworkSize, primme_params *primme);

int Num_reset_update_VWXR_zprimme(Complex_Z *V, Complex_Z *W, int mV, int nV, int ldV,
   Complex_Z *h, int nh, int ldh, double *hVals,
   Complex_Z *X0, int nX0b, int nX0e, int ldX0,
   Complex_Z *X1, int nX1b, int nX1e, int ldX1,
   Complex_Z *evecs, int evecsSize, int nX2b, int nX2e, int ldevecs,
   Complex_Z *Wo, int nWob, int nWoe, int ldWo,
   Complex_Z *R, int nRb, int nRe, int ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   int reset, double machEps, Complex_Z *rwork, int lrwork, primme_params *primme);

#endif
