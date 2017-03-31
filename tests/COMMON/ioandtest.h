/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2017 College of William & Mary,
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
 * File: ioandtest.h
 *
 * Purpose - Definitions in ioandtest.h.
 * 
 ******************************************************************************/

#ifndef IOANDTEST_H
#define IOANDTEST_H

#include "primme_svds.h"
#include "shared_utils.h"
 
int check_solution(const char *checkXFileName, primme_params *primme, double *evals,
                   SCALAR *evecs, double *rnorms, int *perm, int checkInterface);
int writeBinaryEvecsAndPrimmeParams(const char *fileName, SCALAR *X, int *perm,
                                    primme_params *primme);
int readBinaryEvecsAndPrimmeParams(const char *fileName, SCALAR *X, SCALAR **Xout,
                                   int n, int Xcols, int *Xcolsout, int nLocal,
                                   int *perm, primme_params *primme);
int check_solution_svds(const char *checkXFileName, primme_svds_params *primme_svds, double *svals,
                        SCALAR *svecs, double *rnorms, int *perm);
int readBinaryEvecsAndPrimmeSvdsParams(const char *fileName, SCALAR *X, SCALAR **Xout,
                                       int m, int n, int Xcols, int *Xcolsout, int mLocal, int nLocal,
                                       int *perm, primme_svds_params *primme_svds_out);
int writeBinaryEvecsAndPrimmeSvdsParams(const char *fileName, SCALAR *X, int *perm,
                                    primme_svds_params *primme_svds);
#endif
