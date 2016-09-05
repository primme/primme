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

int ortho_zprimme(SCALAR *basis, PRIMME_INT ldBasis, SCALAR *R,
      PRIMME_INT ldR, int b1, int b2, SCALAR *locked, PRIMME_INT ldLocked,
      int numLocked, PRIMME_INT nLocal, PRIMME_INT *iseed, double machEps,
      SCALAR *rwork, size_t *rworkSize, primme_params *primme);

 int ortho_single_iteration_zprimme(SCALAR *Q, PRIMME_INT mQ, PRIMME_INT nQ,
      PRIMME_INT ldQ, SCALAR *X, int *inX, int nX, PRIMME_INT ldX,
      REAL *overlaps, REAL *norms, SCALAR *rwork, size_t *lrwork,
      primme_params *primme);

#endif
