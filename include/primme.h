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
 * File: primme.h
 * 
 * Purpose - Main header with the PRIMME C interface functions.
 * 
 ******************************************************************************/

#ifndef PRIMME_H
#define PRIMME_H

/* A C99 code with complex type is not a valid C++ code. However C++          */
/* compilers usually can take it. Nevertheless in order to avoid the warnings */
/* while compiling in pedantic mode, we use the proper complex type for C99   */
/* (complex double and complex float) and C++ (std::complex<double> and       */
/* std::complex<float>). Of course both complex types are binary compatible.  */

#ifdef __cplusplus
#  include <complex>
#  define PRIMME_COMPLEX_FLOAT std::complex<float>
#  define PRIMME_COMPLEX_DOUBLE std::complex<double>
extern "C" {
#else
#  include <complex.h>
#  define PRIMME_COMPLEX_FLOAT float complex
#  define PRIMME_COMPLEX_DOUBLE double complex
#endif

#if !defined(PRIMME_INT_SIZE) || PRIMME_INT_SIZE == 64
#  include <stdint.h>
#  include <inttypes.h>
#  define PRIMME_INT int64_t
#  define PRIMME_INT_P PRId64
#  define PRIMME_INT_MAX INT64_MAX
#elif PRIMME_INT_SIZE == 0
#  include <limits.h>
#  define PRIMME_INT int
#  define PRIMME_INT_P "d"
#  define PRIMME_INT_MAX INT_MAX
#elif PRIMME_INT == 32
#  include <stdint.h>
#  include <inttypes.h>
#  define PRIMME_INT int32_t
#  define PRIMME_INT_P PRId32
#  define PRIMME_INT_MAX INT32_MAX
#else
#  define PRIMME_INT PRIMME_INT_SIZE
#  define PRIMME_INT_P "d"
#  define PRIMME_INT_MAX INT_MAX
#endif

#include "primme_eigs.h"
#include "primme_svds.h"

#endif /* PRIMME_H */
