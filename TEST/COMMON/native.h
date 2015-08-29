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
 * File: native.h
 * 
 * Purpose - Definitions of CSR matrix functions used by the driver.
 * 
 ******************************************************************************/

#ifndef NATIVE_H
#define NATIVE_H

#include "csr.h"
#include "primme.h"

void CSRMatrixMatvec(void *x, void *y, int *blockSize, primme_params *primme);
int createInvDiagPrecNative(const CSRMatrix *matrix, double shift, double **prec);
void ApplyInvDiagPrecNative(void *x, void *y, int *blockSize, 
                                        primme_params *primme);
int createInvDavidsonDiagPrecNative(const CSRMatrix *matrix, double **prec);
void ApplyInvDavidsonDiagPrecNative(void *x, void *y, int *blockSize, 
                                        primme_params *primme);
int createILUTPrecNative(const CSRMatrix *matrix, double shift, int level,
                         double threshold, double filter, CSRMatrix **prec);
void ApplyILUTPrecNative(void *x, void *y, int *blockSize, primme_params *primme);

#endif

