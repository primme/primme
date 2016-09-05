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
 * File: solve_H.h
 *
 * Purpose - Contains prototype for solving the projected eigenproblem.
 *
 ******************************************************************************/

#ifndef SOLVE_H_H
#define SOLVE_H_H

int solve_H_@(pre)primme(SCALAR *H, int basisSize, int ldH, SCALAR *R, int ldR,
   SCALAR *QtV, int ldQtV, SCALAR *hU, int ldhU, SCALAR *hVecs, int ldhVecs,
   double *hVals, double *hSVals, int numConverged, double machEps, size_t *lrwork,
   SCALAR *rwork, int liwork, int *iwork, primme_params *primme);

int solve_H_RR_@(pre)primme(SCALAR *H, int ldH, SCALAR *hVecs,
   int ldhVecs, double *hVals, int basisSize, int numConverged, size_t *lrwork,
   SCALAR *rwork, int liwork, int *iwork, primme_params *primme);

int prepare_vecs_@(pre)primme(int basisSize, int i0, int blockSize,
      SCALAR *H, int ldH, double *hVals, double *hSVals, SCALAR *hVecs,
      int ldhVecs, int targetShiftIndex, int *arbitraryVecs,
      double smallestResNorm, int *flags, int RRForAll, SCALAR *hVecsRot,
      int ldhVecsRot, double machEps, size_t *rworkSize, SCALAR *rwork,
      int iworkSize, int *iwork, primme_params *primme);

#endif
