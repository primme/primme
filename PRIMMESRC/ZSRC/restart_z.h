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

int restart_zprimme(complex double *V, complex double *W, int nLocal, int basisSize,
       int ldV, double *hVals, double *hSVals, int *flags, int *iev,
       int *ievSize, double *blockNorms, complex double *evecs, int *evecsPerm,
       double *evals, double *resNorms, complex double *evecsHat, int ldevecsHat,
       complex double *M, int ldM, complex double *UDU, int ldUDU, int *ipivot,
       int *numConverged, int *numLocked, int *numConvergedStored,
       complex double *previousHVecs, int *numPrevRetained, int ldpreviousHVecs,
       int numGuesses, double *prevRitzVals, int *numPrevRitzVals, complex double *H,
       int ldH, complex double *Q, int ldQ, complex double *R, int ldR, complex double* QtV,
       int ldQtV, complex double *hU, int ldhU, int newldhU, complex double *hVecs,
       int ldhVecs, int newldhVecs, int *restartSizeOutput,
       int *targetShiftIndex, int *numArbitraryVecs, complex double *hVecsRot,
       int ldhVecsRot, int *restartsSinceReset, int *reset,
       double machEps, complex double *rwork, int rworkSize, int *iwork,
       primme_params *primme);

int Num_reset_update_VWXR_zprimme(complex double *V, complex double *W, int mV, int nV, int ldV,
   complex double *h, int nh, int ldh, double *hVals,
   complex double *X0, int nX0b, int nX0e, int ldX0,
   complex double *X1, int nX1b, int nX1e, int ldX1,
   complex double *evecs, int evecsSize, int nX2b, int nX2e, int ldevecs,
   complex double *Wo, int nWob, int nWoe, int ldWo,
   complex double *R, int nRb, int nRe, int ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   int reset, double machEps, complex double *rwork, int lrwork, primme_params *primme);

int retain_previous_coefficients_zprimme(complex double *hVecs, int ldhVecs,
   complex double *hU, int ldhU, double *hSVals, complex double *previousHVecs,
   int ldpreviousHVecs, int mprevious, int basisSize, int *iev, int blockSize,
   int *flags, int *numPrevRetained, int *iwork, complex double *rwork, int rworkSize,
   primme_params *primme);

#endif
