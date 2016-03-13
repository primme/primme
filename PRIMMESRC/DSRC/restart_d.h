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

void reset_flags_dprimme(int *flags, int first, int last);

int restart_dprimme(double *V, double *W, int nLocal, int basisSize, int ldV,
   double *hVals, double *hSVals, int *flags, int *iev, int *ievSize,
   double *blockNorms, double *evecs, int *evecsPerm, double *evals, double *resNorms,
   double *evecsHat, int ldevecsHat, double *M, int ldM, double *UDU,
   int ldUDU, int *ipivot, int *numConverged, int *numLocked,
   int *numConvergedStored, double *previousHVecs, int *numPrevRetained,
   int ldpreviousHVecs, int numGuesses, double *prevRitzVals, int *numPrevRitzVals,
   double *H, int ldH, double *Q, int ldQ, double *R, int ldR, double* QV, int ldQV,
   double *hU, int ldhU, int newldhU, double *hVecs, int ldhVecs, int newldhVecs,
   int *restartSizeOutput, int *targetShiftIndex, double machEps,
   double *rwork, int rworkSize, int *iwork, primme_params *primme);

int ortho_coefficient_vectors_dprimme(double *hVecs, int basisSize, int ldhVecs,
   int indexOfPreviousVecs, int newBasisSize, int *perm, double *hU, int ldhU,
   double *R, int ldR, int numPrevRetained, int machEps, int *iwork,
   double *rwork, int rworkSize, primme_params *primme);

#endif
