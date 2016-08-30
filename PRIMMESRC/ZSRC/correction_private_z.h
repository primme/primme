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
 * File: correction_private.h
 *
 * Purpose - Definitions used exclusively by correction.c
 *
 ******************************************************************************/

#ifndef CORRECTION_PRIVATE_H
#define CORRECTION_PRIVATE_H

#define INNER_SOLVE_FAILURE -1

static double computeRobustShift(int blockIndex, double resNorm, 
   double *prevRitzVals, int numPrevRitzVals, double *sortedRitzVals, 
   double *approxOlsenShift, int numSorted, int *ilev, primme_params *primme);

static void mergeSort(double *lockedEvals, int numLocked, double *ritzVals, 
   int *flag, int basisSize, double *sortedEvals, int *ilev, int blockSize,
   primme_params *primme);

static void Olsen_preconditioner_block(__PRIMME_COMPLEX_DOUBLE__ *r, __PRIMME_COMPLEX_DOUBLE__ *x,
   int blockSize, __PRIMME_COMPLEX_DOUBLE__ *rwork, primme_params *primme) ;

static void apply_preconditioner_block(__PRIMME_COMPLEX_DOUBLE__ *v, __PRIMME_COMPLEX_DOUBLE__ *result,
   int blockSize, primme_params *primme);

static void setup_JD_projectors(__PRIMME_COMPLEX_DOUBLE__ *x, __PRIMME_COMPLEX_DOUBLE__ *r, __PRIMME_COMPLEX_DOUBLE__ *evecs,
   __PRIMME_COMPLEX_DOUBLE__ *evecsHat, __PRIMME_COMPLEX_DOUBLE__ *Kinvx, __PRIMME_COMPLEX_DOUBLE__ *xKinvx,
   __PRIMME_COMPLEX_DOUBLE__ **Lprojector, __PRIMME_COMPLEX_DOUBLE__ **RprojectorQ, __PRIMME_COMPLEX_DOUBLE__ **RprojectorX,
   int *sizeLprojector, int *sizeRprojectorQ, int *sizeRprojectorX,
   int numLocked, int numConverged, primme_params *primme);


#endif
