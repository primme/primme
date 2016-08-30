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
 * File: init.h
 *
 * Purpose - Function headers for basis initialization
 *
 ******************************************************************************/

#ifndef INIT_H
#define INIT_H

int init_basis_zprimme(__PRIMME_COMPLEX_DOUBLE__ *V, int nLocal, int ldV, __PRIMME_COMPLEX_DOUBLE__ *W, int ldW,
   __PRIMME_COMPLEX_DOUBLE__ *evecs, int ldevecs, __PRIMME_COMPLEX_DOUBLE__ *evecsHat, int ldevecsHat, __PRIMME_COMPLEX_DOUBLE__ *M,
   int ldM, __PRIMME_COMPLEX_DOUBLE__ *UDU, int ldUDU, int *ipivot, double machEps, __PRIMME_COMPLEX_DOUBLE__ *rwork,
   int rworkSize, int *basisSize, int *nextGuess, int *numGuesses,
   double *timeForMV, primme_params *primme);

#endif
