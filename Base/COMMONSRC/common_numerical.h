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
 * File: common_numerical.h
 *
 * Purpose - Contains prototypes for common numerical functions.
 *
 ******************************************************************************/

#ifndef COMMON_NUMERICAL_H
#define COMMON_NUMERICAL_H

#define max(a, b) (a > b ? a : b)
#define min(a, b) (a < b ? a : b)

void Num_dcopy_primme(int n, double *x, int incx, double *y, int incy);
double Num_dlamch_primme(char *cmach);
int Num_imax_primme(int numArgs, int val1, int val2, ...);
double Num_fmin_primme(int numArgs, double val1, double val2, ...);
double Num_fmax_primme(int numArgs, double val1, double val2, ...);

#endif


