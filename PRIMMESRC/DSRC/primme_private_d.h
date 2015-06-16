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
 * File: primme_private.h
 *
 * Purpose - Contains definitions for exclusive use by primme.c
 *
 * Module name      : %M%
 * SID              : %I%
 * Date             : %G%
 ******************************************************************************/


#ifndef DPRIMME_PRIVATE_H
#define DPRIMME_PRIVATE_H

#define ALLOCATE_WORKSPACE_FAILURE -1
#define MALLOC_FAILURE             -2
#define MAIN_ITER_FAILURE          -3

static int allocate_workspace(primme_params *primme, int allocate);
static int check_input(double *evals, double *evecs, double *resNorms,
		       primme_params *primme);
#endif
