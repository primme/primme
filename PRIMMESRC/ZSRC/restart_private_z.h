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

static int restart_RR(Complex_Z *H, int ldH, Complex_Z *hVecs, int ldhVecs,
   int newldhVecs, double *hVals, int restartSize, int basisSize, int numLocked,
   int numPrevRetained, int indexOfPreviousVecs, int *hVecsPerm,
   int rworkSize, Complex_Z *rwork, int *iwork, primme_params *primme);

static int restart_ref(Complex_Z *V, int ldV, Complex_Z *W, int ldW, Complex_Z *H,
   int ldH, Complex_Z *Q, int nLocal, int ldQ, Complex_Z *R, int ldR, Complex_Z *hU,
   int ldhU, int newldhU, Complex_Z *hVecs, int ldhVecs, int newldhVecs,
   double *hVals, double *hSVals, int *restartPerm, int *hVecsPerm, int restartSize,
   int basisSize, int numPrevRetained, int indexOfPreviousVecs,
   int numConvergedBeforeRestart, int numConverged, int rworkSize,
   Complex_Z *rwork, int *iwork, double machEps, primme_params *primme);

static int compute_submatrix(Complex_Z *X, int nX, int ldX, 
   Complex_Z *H, int nH, int ldH, Complex_Z *R, int ldR,
   Complex_Z *rwork, int lrwork);

#endif /* RESTART_PRIVATE_H */
