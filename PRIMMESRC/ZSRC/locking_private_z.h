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
 * File: locking_private.h
 *
 * Purpose - Definitions used exclusively by locking.c
 *
 ******************************************************************************/

#ifndef LOCKING_PRIVATE_H
#define LOCKING_PRIVATE_H

#define ORTHO_FAILURE             -1 
#define INIT_BLOCK_KRYLOV_FAILURE -2
#define SOLVE_H_FAILURE           -3
#define UDUDECOMPOSE_FAILURE      -4

static int swap_flagVecs_toEnd(int basisSize, int flagValue, Complex_Z *V, 
   Complex_Z *W, Complex_Z *H, double *hVals, int *flag, primme_params *primme);

static void insertionSort(double newVal, double *evals, double newNorm,
   double *resNorms, int *perm, int numConverged, primme_params *primme);

#endif
