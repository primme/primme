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
 * File: common_numerical_private.h
 *
 * Purpose - Contains definitions and prototypes for exclusive use with 
 *           common_numerical.c. There are various definitions for use 
 *           with Sun, IBM, and Cray.
 *
 ******************************************************************************/

#if !defined(NUM_SUN) && !defined(NUM_IBM) && !defined(NUM_CRAY)
#define NUM_SUN
#endif

#ifdef NUM_SUN

#define DCOPY  dcopy_
#define DLAMCH dlamch_

#elif defined(NUM_IBM)

#define DCOPY  dcopy
#define DLAMCH dlamch

#ifdef NUM_ESSL
#include <essl.h>
#endif

#elif defined(NUM_CRAY)
#include <fortran.h>
#include <string.h>

#define DCOPY  SCOPY
#define DLAMCH SLAMCH

#endif

#ifdef Cplusplus
extern "C"
{
#endif /* Cplusplus */

#ifndef NUM_CRAY

void DCOPY(int *n, double *x, int *incx, double *y, int *incy);
double DLAMCH(char *cmach);

#ifdef NUM_ESSL
#endif

#else

void DCOPY(int *n, double *x, int *incx, double *y, int *incy);
double DLAMCH(_fcd cmach_fcd);

#endif

#ifdef Cplusplus
}
#endif /* Cplusplus */
