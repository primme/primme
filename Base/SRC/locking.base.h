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

int lock_vectors_@(pre)primme(double tol, double *aNormEstimate, double *maxConvTol, 
   int *basisSize, int *numLocked, int *numGuesses, int *nextGuess,
   @(type) *V, @(type) *W, @(type) *H, @(type) *evecsHat, @(type) *M, 
   @(type) *UDU, int *ipivot, double *hVals, @(type) *hVecs, 
   @(type) *evecs, double *evals, int *perm, double machEps, 
   double *resNorms, int *numPrevRitzVals, double *prevRitzVals, 
   int *flag, @(type) *rwork, int rworkSize, int *iwork, 
   int *LockingProblem, primme_params *primme);

#endif
