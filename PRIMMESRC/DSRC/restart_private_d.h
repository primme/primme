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

static int restart_RR(double *H, int ldH, double *hVecs, int ldhVecs,
   int newldhVecs, double *hVals, int restartSize, int basisSize, int numLocked,
   int numPrevRetained, int indexOfPreviousVecs, int *hVecsPerm,
   int rworkSize, double *rwork, int *iwork, primme_params *primme);

static int restart_qr(double *V, int ldV, double *W, int ldW, double *H,
   int ldH, double *Q, int nLocal, int ldQ, double *R, int ldR, double *QV,
   int ldQV, double *hU, int ldhU, int newldhU, double *hVecs, int ldhVecs,
   int newldhVecs, double *hVals, double *hSVals, int *restartPerm, int *hVecsPerm,
   int restartSize, int basisSize, int numPrevRetained, int indexOfPreviousVecs,
   int numConvergedBeforeRestart, int numConverged, int rworkSize,
   double *rwork, int *iwork, double machEps, primme_params *primme);

static int compute_submatrix(double *X, int nX, int ldX, 
   double *H, int nH, int ldH, double *R, int ldR,
   double *rwork, int lrwork);

#endif /* RESTART_PRIVATE_H */
