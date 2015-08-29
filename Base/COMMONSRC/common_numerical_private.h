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
 * File: common_numerical_private.h
 *
 * Purpose - Contains definitions and prototypes for exclusive use with 
 *           common_numerical.c. There are various definitions for use 
 *           with Sun, IBM, and Cray.
 *
 ******************************************************************************/

#ifndef COMMON_NUMERICAL_PRIVATE_H
#define COMMON_NUMERICAL_PRIVATE_H

#include "common_numerical.h"

#ifdef F77UNDERSCORE
#define FORTRAN_FUNCTION(X) X ## _
#else
#define FORTRAN_FUNCTION(X) X
#endif

#ifndef NUM_CRAY

#define DCOPY  FORTRAN_FUNCTION(dcopy)
#define DLAMCH FORTRAN_FUNCTION(dlamch)

#ifdef NUM_ESSL
#include <essl.h>
#endif

#else /* NUM_CRAY */

#include <fortran.h>
#include <string.h>

#define DCOPY  SCOPY
#define DLAMCH SLAMCH

#endif /* NUM_CRAY */

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#ifndef NUM_CRAY

void DCOPY(PRIMME_BLASINT *n, double *x, PRIMME_BLASINT *incx, double *y, PRIMME_BLASINT *incy);
double DLAMCH(const char *cmach);

#else

void DCOPY(PRIMME_BLASINT *n, double *x, PRIMME_BLASINT *incx, double *y, PRIMME_BLASINT *incy);
double DLAMCH(_fcd cmach_fcd);

#endif

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* COMMON_NUMERICAL_PRIVATE_H */
