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

static void restart_X(Complex_Z *X, Complex_Z *hVecs, int nLocal, 
   int basisSize, int restartSize, Complex_Z *rwork, int rworkSize);

static int restart_H(Complex_Z *H, Complex_Z *hVecs, double *hVals, 
   int restartSize, int basisSize, Complex_Z *previousHVecs, 
   int numPrevRetained, int indexOfPreviousVecs, int rworkSize, 
   Complex_Z *rwork, primme_params *primme);

static int dtr(int numLocked, Complex_Z *hVecs, double *hVals, int *flags, 
 int basisSize, int numFree, int *iev, Complex_Z *rwork, primme_params *primme);

static int pack_converged_coefficients(int *restartSize, int basisSize, 
   int *numPrevRetained, int numLocked, int numGuesses, Complex_Z *hVecs, 
   double *hVals, int *flag, primme_params *primme);

static int combine_retained_vectors(double *hVals, int *flags, Complex_Z *hVecs, 
   int basisSize, int *restartSize, int numPacked, Complex_Z *previousHVecs, 
   int *numPrevRetained, double machEps, Complex_Z *rwork, primme_params *primme);

static void compute_submatrix(Complex_Z *previousHVecs, int numPrevRetained, 
   Complex_Z *H, int basisSize, int maxBasisSize, Complex_Z *subMatrix, 
   Complex_Z *rwork);

static int insert_submatrix(Complex_Z *H, double *hVals, Complex_Z *hVecs, 
   int restartSize, Complex_Z *subMatrix, int numPrevRetained, 
   int indexOfPreviousVecs, int rworkSize, Complex_Z *rwork, 
   primme_params *primme);

static void apply_preconditioner_block(Complex_Z *v, Complex_Z *result,
   int blockSize, primme_params *primme);

#endif /* RESTART_PRIVATE_H */
