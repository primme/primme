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
 * File: rsbw.h
 * 
 * Purpose - librsb wrapper.
 * 
 ******************************************************************************/

#ifndef RSBW_H
#define RSBW_H

#include <rsb.h>
#include <blas_sparse.h>
#include "primme_svds.h"

int readMatrixRSB(const char* matrixFileName, blas_sparse_matrix *matrix, double *fnorm);
void RSBMatvec(void *x, void *y, int *blockSize, primme_params *primme);
void RSBMatvecSVD(void *x, int *ldx, void *y, int *ldy, int *blockSize,
                  int *trans, primme_svds_params *primme_svds);
int createInvDiagPrecRSB(blas_sparse_matrix matrix, double shift, double **prec);
int createInvNormalPrecRSB(blas_sparse_matrix, double shift, double **prec);

#endif

