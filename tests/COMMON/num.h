/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2017 College of William & Mary,
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

#include "../../src/include/template.h"
#include "../../src/include/blaslapack.h"
#ifdef USE_COMPLEX
#  ifndef __cplusplus
#     define IMAGINARY _Complex_I
#  else
#     define IMAGINARY std::complex<SCALAR>(0.0, 1.0)
#  endif
#else
#  define IMAGINARY 0.0
#endif
#define Sprimme CONCAT(SCALAR_PRE,primme)
#define Sprimme_svds CONCAT(SCALAR_PRE,primme_svds)
#if !(defined (__APPLE__) && defined (__MACH__))
#  include <malloc.h> /* malloc */
#endif
#include <stdlib.h>   /* malloc, free */
#define primme_calloc(N,S,D) (malloc((N)*(S)))
#define ASSERT_MSG(COND, RETURN, ...) { if (!(COND)) {fprintf(stderr, "Error in " __FUNCT__ ": " __VA_ARGS__); return (RETURN);} }

#endif /* NUM_H */
