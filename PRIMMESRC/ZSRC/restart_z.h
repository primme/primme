/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2005  James R. McCombs,  Andreas Stathopoulos
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
 * File: restart.h
 *
 * Purpose - Contains prototypes for functions related to restarting.
 *
 * Module name      : %M%
 * SID              : %I%
 * Date             : %G%
 ******************************************************************************/

#ifndef RESTART_H
#define RESTART_H

void reset_flags_zprimme(int *flag, int first, int last);

int restart_zprimme(Complex_Z *V, Complex_Z *W, Complex_Z *H, Complex_Z *hVecs, 
   double *hVals, int *flags, int *iev, Complex_Z *evecs, Complex_Z *evecsHat, 
   Complex_Z *M, Complex_Z *UDU, int *ipivot, int basisSize, int numConverged, 
   int *numConvergedStored, int numLocked, int numGuesses, 
   Complex_Z *previousHVecs, int numPrevRetained, double machEps, 
   Complex_Z *rwork, int rworkSize, primme_params *primme);

#endif
