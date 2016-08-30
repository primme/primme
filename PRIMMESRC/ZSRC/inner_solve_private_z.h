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

static int apply_projected_preconditioner(complex double *v, complex double *Q, 
   complex double *RprojectorQ, complex double *x, complex double *RprojectorX, 
   int sizeRprojectorQ, int sizeRprojectorX, complex double *xKinvx, 
   complex double *UDU, int *ipivot, complex double *result, complex double *rwork, 
   primme_params *primme);

static int apply_skew_projector(complex double *Q, complex double *Qhat, complex double *UDU, 
   int *ipivot, int numCols, complex double *v, complex double *rwork, 
   primme_params *primme);

static void apply_projected_matrix(complex double *v, double shift, complex double *Q, 
   int dimQ, complex double *result, complex double *rwork, primme_params *primme);

static void apply_projector(complex double *Q, int numCols, complex double *v, 
   complex double *rwork, primme_params *primme);

static complex double dist_dot(complex double *x, int incx,
   complex double *y, int incy, primme_params *primme);

#endif
