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
 * File: primme_svds_private.h
 *
 * Purpose - Contains definitions for exclusive use by primme_svds.c
 *
 ******************************************************************************/


#ifndef DPRIMME_SVDS_PRIVATE_H
#define DPRIMME_SVDS_PRIVATE_H

#include "primme_svds.h"

#define ALLOCATE_WORKSPACE_FAILURE -1
#define MALLOC_FAILURE             -3

static int primme_svds_check_input(double *svals, double *svecs, 
        double *resNorms, primme_svds_params *primme_svds);
static double* copy_last_params_from_svds(primme_svds_params *primme_svds, int stage,
		double *svals, double *svecs, double *rnorms);
static void copy_last_params_to_svds(primme_svds_params *primme_svds, int stage,
                                     double *svals, double *svecs, double *rnorms);
static void matrixMatvecSVDS(void *x_, void *y_, int *blockSize, primme_params *primme);
static void applyPreconditionerSVDS(void *x, void *y, int *blockSize, primme_params *primme);
static void Num_scalInv_dmatrix(double *x, int m, int n, int ldx, double *factors,
                                       primme_svds_params *primme_svds);
static int allocate_workspace_svds(primme_svds_params *primme_svds, int allocate);
#endif
