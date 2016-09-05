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

int restart_@(pre)primme(SCALAR *V, SCALAR *W, PRIMME_INT nLocal, int basisSize,
       PRIMME_INT ldV, double *hVals, double *hSVals, int *flags, int *iev,
       int *ievSize, double *blockNorms, SCALAR *evecs, int *evecsPerm,
       double *evals, double *resNorms, SCALAR *evecsHat, PRIMME_INT ldevecsHat,
       SCALAR *M, int ldM, SCALAR *UDU, int ldUDU, int *ipivot,
       int *numConverged, int *numLocked, int *numConvergedStored,
       SCALAR *previousHVecs, int *numPrevRetained, int ldpreviousHVecs,
       int numGuesses, double *prevRitzVals, int *numPrevRitzVals, SCALAR *H,
       int ldH, SCALAR *Q, PRIMME_INT ldQ, SCALAR *R, int ldR, SCALAR* QtV,
       int ldQtV, SCALAR *hU, int ldhU, int newldhU, SCALAR *hVecs,
       int ldhVecs, int newldhVecs, int *restartSizeOutput,
       int *targetShiftIndex, int *numArbitraryVecs, SCALAR *hVecsRot,
       int ldhVecsRot, int *restartsSinceReset, int *reset,
       double machEps, SCALAR *rwork, size_t *rworkSize, int *iwork,
       int iworkSize, primme_params *primme);

int Num_reset_update_VWXR_@(pre)primme(SCALAR *V, SCALAR *W, PRIMME_INT mV,
   int nV, PRIMME_INT ldV,
   SCALAR *h, int nh, int ldh, double *hVals,
   SCALAR *X0, int nX0b, int nX0e, PRIMME_INT ldX0,
   SCALAR *X1, int nX1b, int nX1e, PRIMME_INT ldX1,
   SCALAR *evecs, int evecsSize, int nX2b, int nX2e, PRIMME_INT ldevecs,
   SCALAR *Wo, int nWob, int nWoe, PRIMME_INT ldWo,
   SCALAR *R, int nRb, int nRe, PRIMME_INT ldR, double *Rnorms,
   double *rnorms, int nrb, int nre,
   int reset, double machEps, SCALAR *rwork, size_t *lrwork,
   primme_params *primme);

int retain_previous_coefficients_@(pre)primme(SCALAR *hVecs, int ldhVecs,
   SCALAR *hU, int ldhU, SCALAR *previousHVecs, int ldpreviousHVecs,
   int mprevious, int basisSize, int *iev, int blockSize, int *flags,
   int *numPrevRetained, int *iwork, int iworkSize, primme_params *primme);

#endif
