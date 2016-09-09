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
 * File: notemplate.h
 *
 * Purpose - Protect functions that are built only once and shouldn't use
 *             SCALAR nor REAL nor derivates.
 *
 ******************************************************************************/

#ifndef NOTEMPLATE_H
#define NOTEMPLATE_H

#ifdef SCALAR
#  undef SCALAR
#endif

#ifdef REAL
#  undef REAL
#endif

#ifdef SCALAR_SUF
#  undef SCALAR_SUF
#endif

#ifdef REAL_SUF
#  undef REAL_SUF
#endif

#ifdef REAL_PART
#  undef REAL_PART
#endif

#ifdef ABS
#  undef ABS
#endif

#ifdef CONJ
#  undef CONJ
#endif

#ifdef TEMPLATE_H
#  undef TEMPLATE_H
#endif

#endif
