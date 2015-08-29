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

void reset_flags_@(pre)primme(int *flag, int first, int last);

int restart_@(pre)primme(@(type) *V, @(type) *W, @(type) *H, @(type) *hVecs, 
   double *hVals, int *flags, int *iev, @(type) *evecs, @(type) *evecsHat, 
   @(type) *M, @(type) *UDU, int *ipivot, int basisSize, int numConverged, 
   int *numConvergedStored, int numLocked, int numGuesses, 
   @(type) *previousHVecs, int numPrevRetained, double machEps, 
   @(type) *rwork, int rworkSize, primme_params *primme);

#endif
