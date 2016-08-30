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
#  define PRIMME_NUM complex double
#  define SUF(NAME) NAME ## _zprimme
#  define PREFIX(NAME) z ## NAME
#else
#  include "../../PRIMMESRC/DSRC/numerical_d.h"
#  define PRIMME_NUM double
#  define SUF(NAME) NAME ## _dprimme
#  define PREFIX(NAME) d ## NAME
#endif

#define ASSERT_MSG(COND, RETURN, ...) { if (!(COND)) {fprintf(stderr, "Error in " __FUNCT__ ": " __VA_ARGS__); return (RETURN);} }

#endif /* NUM_H */
