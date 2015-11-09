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
 * File: solve_H.h
 *
 * Purpose - Contains prototype for solving the projected eigenproblem.
 *
 * Module name      : %M%
 * SID              : %I%
 * Date             : %G%
 ******************************************************************************/

#ifndef SOLVE_H_H
#define SOLVE_H_H

int solve_H_@(pre)primme(@(type) *H, @(type) *hVecs, @(type) *Q, @(type) *R, double *hVals, 
   int basisSize, int maxBasisSize, double *largestEval, int numLocked, double machEps,
   int lrwork, @(type) *rwork, int *perm, primme_params *primme, 
   @(type) *V, @(type) *W, int recentlyConverged);

#ifdefarithm L_DEFCPLX
void permute_evecs_zprimme(double *evecs, int elemSize, int *perm, 
   double *rwork, int nev, int nLocal);
#endifarithm
#ifdefarithm L_DEFREAL
void permute_evecs_dprimme(double *evecs, int *perm, double *rwork, 
   int nev, int nLocal);
#endifarithm


#endif
