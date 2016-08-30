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
 * File: inner_solve.h
 *
 * Purpose - Contains prototypes for calling the inner solve. 
 *
 ******************************************************************************/

#ifndef INNER_SOLVE_H
#define INNER_SOLVE_H

int inner_solve_zprimme(__PRIMME_COMPLEX_DOUBLE__ *x, __PRIMME_COMPLEX_DOUBLE__ *r, double *rnorm, __PRIMME_COMPLEX_DOUBLE__ *evecs,
   __PRIMME_COMPLEX_DOUBLE__ *evecsHat, __PRIMME_COMPLEX_DOUBLE__ *UDU, int *ipivot, __PRIMME_COMPLEX_DOUBLE__ *xKinvx, 
   __PRIMME_COMPLEX_DOUBLE__ *Lprojector, __PRIMME_COMPLEX_DOUBLE__ *RprojectorQ, __PRIMME_COMPLEX_DOUBLE__ *RprojectorX,
   int sizeLprojector, int sizeRprojectorQ, int sizeRprojectorX,
   __PRIMME_COMPLEX_DOUBLE__ *sol, double eval, double shift, double eresTol, 
   double aNormEstimate, double machEps, __PRIMME_COMPLEX_DOUBLE__ *rwork, int rworkSize, 
   primme_params *primme);

#endif
