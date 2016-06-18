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

int ortho_@(pre)primme(@(type) *basis, int ldBasis, @(type) *R, int ldR,
   int b1, int b2, @(type) *locked, int ldLocked, int numLocked,
   int nLocal, int *iseed, double machEps, @(type) *rwork, int rworkSize,
   primme_params *primme);

int ortho_single_iteration_@(pre)primme(@(type) *Q, int mQ, int nQ, int ldQ, @(type) *X,
   int *inX, int nX, int ldX, double *overlaps, double *norms, @(type) *rwork, int lrwork,
   primme_params *primme);

#endif
