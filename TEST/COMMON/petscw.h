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
 * File: petscw.h
 * 
 * Purpose - Definitions of the PETSc wrapper used by the driver.
 * 
 ******************************************************************************/

#ifndef PETSCW_H

#include <petscpc.h>
#include <petscmat.h>
#include "primme.h"
#include "num.h"

int readMatrixPetsc(const char* matrixFileName, int *m, int *n, int *mLocal, int *nLocal,
                    int *numProcs, int *procID, Mat **matrix, double *fnorm_, int **perm);
void PETScMatvec(void *x, void *y, int *blockSize, primme_params *primme);
void ApplyPCPrecPETSC(void *x, void *y, int *blockSize, primme_params *primme);

#define PETSCW_H
#endif

