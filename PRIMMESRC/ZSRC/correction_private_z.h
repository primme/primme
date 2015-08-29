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

static void Olsen_preconditioner_block(Complex_Z *r, Complex_Z *x,
   int blockSize, Complex_Z *rwork, primme_params *primme) ;

static void apply_preconditioner_block(Complex_Z *v, Complex_Z *result,
   int blockSize, primme_params *primme);

static void setup_JD_projectors(Complex_Z *x, Complex_Z *r, Complex_Z *evecs,
   Complex_Z *evecsHat, Complex_Z *Kinvx, Complex_Z *xKinvx,
   Complex_Z **Lprojector, Complex_Z **RprojectorQ, Complex_Z **RprojectorX,
   int *sizeLprojector, int *sizeRprojectorQ, int *sizeRprojectorX,
   int numLocked, int numConverged, primme_params *primme);


#endif
