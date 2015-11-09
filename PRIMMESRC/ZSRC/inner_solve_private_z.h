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
 * File: inner_solve_private.h
 *
 * Purpose - Prototypes used exclusively by inner_solve.c
 *
 ******************************************************************************/

#ifndef INNER_SOLVE_PRIVATE_H
#define INNER_SOLVE_PRIVATE_H

#define APPLYPROJECTEDPRECONDITIONER_FAILURE -1
#define APPLYSKEWPROJECTOR_FAILURE           -2
#define UDUSOLVE_FAILURE                     -3

static int apply_projected_preconditioner(Complex_Z *v, Complex_Z *Q, 
   Complex_Z *RprojectorQ, Complex_Z *x, Complex_Z *RprojectorX, 
   int sizeRprojectorQ, int sizeRprojectorX, Complex_Z *xKinvx, 
   Complex_Z *UDU, int *ipivot, Complex_Z *result, Complex_Z *rwork, 
   primme_params *primme);

static int apply_skew_projector(Complex_Z *Q, Complex_Z *Qhat, Complex_Z *UDU, 
   int *ipivot, int numCols, Complex_Z *v, Complex_Z *rwork, 
   primme_params *primme);

static void apply_projected_matrix(Complex_Z *v, double shift, Complex_Z *Q, 
   int dimQ, Complex_Z *result, Complex_Z *rwork, primme_params *primme);

static void apply_projector(Complex_Z *Q, int numCols, Complex_Z *v, 
   Complex_Z *rwork, primme_params *primme);

static Complex_Z dist_dot(Complex_Z *x, int incx,
   Complex_Z *y, int incy, primme_params *primme);

#endif
