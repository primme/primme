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
 * File: restart_private.h
 *
 * Purpose - Definitions used exclusively by restart.c
 *
 ******************************************************************************/

#ifndef RESTART_PRIVATE_H
#define RESTART_PRIVATE_H

#define INSERT_SUBMATRIX_FAILURE -1
#define RESTART_H_FAILURE        -2
#define NUM_DSYEV_FAILURE        -3
#define UDUDECOMPOSE_FAILURE     -4
#define PSEUDOLOCK_FAILURE       -5

static void restart_X(@(type) *X, @(type) *hVecs, int nLocal, 
   int basisSize, int restartSize, @(type) *rwork, int rworkSize);

static int restart_H(@(type) *H, @(type) *hVecs, double *hVals, 
   int restartSize, int basisSize, @(type) *previousHVecs, 
   int numPrevRetained, int indexOfPreviousVecs, int rworkSize, 
   @(type) *rwork, primme_params *primme);

static int dtr(int numLocked, @(type) *hVecs, double *hVals, int *flags, 
 int basisSize, int numFree, int *iev, @(type) *rwork, primme_params *primme);

static int pack_converged_coefficients(int *restartSize, int basisSize, 
   int *numPrevRetained, int numLocked, int numGuesses, @(type) *hVecs, 
   double *hVals, int *flag, primme_params *primme);

static int combine_retained_vectors(double *hVals, int *flags, @(type) *hVecs, 
   int basisSize, int *restartSize, int numPacked, @(type) *previousHVecs, 
   int *numPrevRetained, double machEps, @(type) *rwork, primme_params *primme);

static void compute_submatrix(@(type) *previousHVecs, int numPrevRetained, 
   @(type) *H, int basisSize, int maxBasisSize, @(type) *subMatrix, 
   @(type) *rwork);

static int insert_submatrix(@(type) *H, double *hVals, @(type) *hVecs, 
   int restartSize, @(type) *subMatrix, int numPrevRetained, 
   int indexOfPreviousVecs, int rworkSize, @(type) *rwork, 
   primme_params *primme);

static void apply_preconditioner_block(@(type) *v, @(type) *result,
   int blockSize, primme_params *primme);

#endif /* RESTART_PRIVATE_H */
