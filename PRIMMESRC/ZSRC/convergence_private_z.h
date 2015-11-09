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
 * File: convergence_private.h
 *
 * Purpose - Definitions used exclusively by convergence.c
 *
 ******************************************************************************/

#ifndef CONVERGENCE_PRIVATE_H
#define CONVERGENCE_PRIVATE_H

static void compute_resnorms(Complex_Z *V, Complex_Z *W, Complex_Z *hVecs, 
   double *hVals, int basisSize, double *blockNorms, int *iev, int left,
   int right, void *rwork, primme_params *primme);

static void print_residuals(double *ritzValues, double *blockNorms, 
   int numConverged, int numLocked, int *iev, int left, int right, 
   primme_params *primme);

static void swap_UnconvVecs(Complex_Z *V, Complex_Z *W, int nLocal, 
   int basisSize, int *iev, int *flag, double *blockNorms, int dimEvecs, 
   int blockSize, int left);

static void replace_vectors(int *iev, int *flag, int blockSize, int basisSize,
   int numVacancies, int *left, int *right, int *ievMax);

static void check_practical_convergence(Complex_Z *V, Complex_Z *W, 
  Complex_Z *evecs, int numLocked, int basisSize, int blockSize, int start, 
  int numToProject, int *iev, int *flags, double *blockNorms, double tol, 
  int *recentlyConverged, int *numVacancies, Complex_Z *rwork, 
  primme_params *primme);

#endif /* CONVERGENCE_PRIVATE_H */
