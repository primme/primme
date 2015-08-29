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
 * File: num.h
 * 
 * Purpose - Definitions used by the driver.
 * 
 ******************************************************************************/

#ifndef NUM_H
#define NUM_H

#ifdef USE_DOUBLECOMPLEX
#  include "../../PRIMMESRC/ZSRC/numerical_z.h"
#  include <complex.h>
#  undef I
#  define IMAGINARY _Complex_I
#  define PRIMME_NUM complex double
#  define REAL_PART(x) (creal(x))
#  define REAL_PARTZ(x) ((x).r)
#  define CONJ(x) (conj(x))
#  define SUF(NAME) NAME ## _zprimme
#  define SUFX(NAME) Num_z ## NAME ## _zprimme
#  define PREFIX(NAME) z ## NAME
#  define COMPLEXZ(X) ((Complex_Z*)(X))
   static inline Complex_Z COMPLEXZV(PRIMME_NUM x) {Complex_Z z={creal(x), cimag(x)}; return z;}
#else
#  include "../../PRIMMESRC/DSRC/numerical_d.h"
#  define IMAGINARY 0.0
#  define PRIMME_NUM double
#  define REAL_PART(x) (x)
#  define REAL_PARTZ(x) (x)
#  define CONJ(x) (x)
#  define SUF(NAME) NAME ## _dprimme
#  define SUFX(NAME) Num_d ## NAME ## _dprimme
#  define PREFIX(NAME) d ## NAME
#  define COMPLEXZ(X) (X)
#  define COMPLEXZV(X) (X)
#endif
#define MACHINE_EPSILON 1.11e-16

#ifdef F77UNDERSCORE
#define FORTRAN_FUNCTION(X) X ## _
#else
#define FORTRAN_FUNCTION(X) X
#endif

#ifndef max
#  define max(a, b) ((a) > (b) ? (a) : (b))
#  define min(a, b) ((a) < (b) ? (a) : (b))
#endif

#endif /* NUM_H */
