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
 * File: ortho.h
 *
 * Purpose - Header file containing prototypes for in-core and 
 *           out-of-core orthogonalization routines.  This header
 *           is to be included in user/developer source files that
 *           need to call these routines.
 *
 ******************************************************************************/

#ifndef ORTHO_H
#define ORTHO_H

int ortho_zprimme(__PRIMME_COMPLEX_DOUBLE__ *basis, int ldBasis, __PRIMME_COMPLEX_DOUBLE__ *R, int ldR,
   int b1, int b2, __PRIMME_COMPLEX_DOUBLE__ *locked, int ldLocked, int numLocked,
   int nLocal, int *iseed, double machEps, __PRIMME_COMPLEX_DOUBLE__ *rwork, int rworkSize,
   primme_params *primme);

int ortho_single_iteration_zprimme(__PRIMME_COMPLEX_DOUBLE__ *Q, int mQ, int nQ, int ldQ, __PRIMME_COMPLEX_DOUBLE__ *X,
   int *inX, int nX, int ldX, double *overlaps, double *norms, __PRIMME_COMPLEX_DOUBLE__ *rwork, int lrwork,
   primme_params *primme);

#endif
