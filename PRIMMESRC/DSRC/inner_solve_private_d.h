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

static int apply_projected_preconditioner(double *v, double *Q, 
   double *RprojectorQ, double *x, double *RprojectorX, 
   int sizeRprojectorQ, int sizeRprojectorX, double *xKinvx, 
   double *UDU, int *ipivot, double *result, double *rwork, 
   primme_params *primme);

static int apply_skew_projector(double *Q, double *Qhat, double *UDU, 
   int *ipivot, int numCols, double *v, double *rwork, 
   primme_params *primme);

static void apply_projected_matrix(double *v, double shift, double *Q, 
   int dimQ, double *result, double *rwork, primme_params *primme);

static void apply_projector(double *Q, int numCols, double *v, 
   double *rwork, primme_params *primme);

static double dist_dot(double *x, int incx,
   double *y, int incy, primme_params *primme);

#endif
