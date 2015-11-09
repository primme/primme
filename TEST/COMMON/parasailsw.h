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
 * File: parasailsw.h
 * 
 * Purpose - Definitions of the Parasails wrapper used by the driver.
 * 
 ******************************************************************************/

#ifndef PARASAILSW_H

#include "ParaSails.h"
#include <mpi.h>
#include "primme.h"
#include "csr.h"

int readMatrixAndPrecondParaSails(const char* matrixFileName, double shift,
                                  int level, double threshold, double filter,
                                  int isymm, MPI_Comm comm, double *fnorm,
                                  int *n, int *nLocal_, int *numProc_, int *procID_,
                                  Matrix **pmatrix, ParaSails **pfactor);
void ParaSailsMatrixMatvec(void *x, void *y, int *blockSize, 
                           primme_params *primme);
void ApplyPrecParaSails(void *x, void *y, int *blockSize,
   primme_params *primme);

#define PARASAILS_H
#endif
