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

int solve_H_zprimme(Complex_Z *H, int basisSize, int ldH, Complex_Z *R, int ldR,
   Complex_Z *hU, int ldhU, Complex_Z *hVecs, int ldhVecs, double *hVals, double *hSVals,
   int numConverged, int lrwork, Complex_Z *rwork, int *iwork, primme_params *primme);

int solve_H_RR_zprimme(Complex_Z *H, int maxBasisSize, Complex_Z *hVecs,
   int ldhVecs, double *hVals, int basisSize, int numLocked, int lrwork,
   Complex_Z *rwork, int *iwork, primme_params *primme);

int solve_H_Ref_zprimme(Complex_Z *H, int ldH, Complex_Z *hVecs,
   int ldhVecs, Complex_Z *hU, int ldhU, double *hSVals, Complex_Z *R, int ldR,
   double *hVals, int basisSize, int lrwork, Complex_Z *rwork, primme_params *primme);

#endif
