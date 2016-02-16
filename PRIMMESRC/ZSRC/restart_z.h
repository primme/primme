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
 * File: restart.h
 *
 * Purpose - Contains prototypes for functions related to restarting.
 *
 ******************************************************************************/

#ifndef RESTART_H
#define RESTART_H

int restart_zprimme(int *restartSize, Complex_Z *V, Complex_Z *W, int nLocal,
   int basisSize, int ldV, Complex_Z **X, Complex_Z **R, Complex_Z *hVecs, int ldhVecs,
   int *hVecsperm, double *hVals, int *flags, int *iev, int *ievSize,
   double *blockNorms, Complex_Z *evecs, double *evals, double *resNorms,
   Complex_Z *evecsHat, int ldevecsHat, Complex_Z *M, int ldM, int *numConverged,
   int *numConvergedStored, Complex_Z *previousHVecs, int *numPrevRetained,
   int ldpreviousHVecs, int *indexOfPreviousVecs, int *Vperm, double machEps,
   Complex_Z *rwork, int rworkSize, int *iwork, primme_params *primme);

int after_restart_zprimme(Complex_Z *V, int ldV, Complex_Z *W, int ldW,
   Complex_Z *H, int ldH, Complex_Z *Q, int nLocal, int ldQ, Complex_Z *R, int ldR,
   Complex_Z *hU, int ldhU, int newldhU, Complex_Z *hVecs, int ldhVecs, int newldhVecs,
   double *hVals, double *hSVals, int *hVecsperm, int *Vperm,
   int restartSize, int basisSize, int numPrevRetained,
   int indexOfPreviousVecs, Complex_Z *evecs, int *evecsSize,
   int ldevecs, Complex_Z *evecsHat, int ldevecsHat, Complex_Z *M, int ldM, Complex_Z *UDU,
   int ldUDU, int *ipivot, int numConvergedBeforeRestart, int numConverged,
   int rworkSize, Complex_Z *rwork, int *iwork, double machEps, primme_params *primme);

void reset_flags_zprimme(int *flags, int first, int last);

int dtr_z(int numLocked, Complex_Z *hVecs, double *hVals, int *flags, 
  int basisSize, int numFree, int *iev, Complex_Z *rwork, primme_params *primme);

#endif
