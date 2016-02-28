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
 * File: locking.h
 *
 * Purpose - Contains prototypes for locking routines.
 *
 ******************************************************************************/


#ifndef LOCKING_H
#define LOCKING_H

int restart_locking_@(pre)primme(int *restartSize, @(type) *V, @(type) *W,
   int nLocal, @(type) *hR, int ldhR, @(type) *hU, int ldhU,
   int basisSize, int ldV, @(type) **X, @(type) **R, @(type) *hVecs,
   int ldhVecs, int *restartPerm, double *hVals, int *flags, int *iev, int *ievSize,
   double *blockNorms, @(type) *evecs, double *evals, int *numConverged, int *numLocked,
   double *resNorms, int *evecsperm, int numGuesses, @(type) *previousHVecs,
   int *numPrevRetained, int ldpreviousHVecs, int *indexOfPreviousVecs, int *hVecsPerm,
   double machEps, @(type) *rwork, int rworkSize, int *iwork, primme_params *primme);

#endif
