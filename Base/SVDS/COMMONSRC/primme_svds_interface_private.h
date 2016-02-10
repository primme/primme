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
 **********************************************************************
 * File: primme_interface_private.h
 *
 * Purpose - Contains definitions and prototypes for exclusive use with 
 *           primme_svds_interface.c.
 *
 ******************************************************************************/

#ifndef PRIMME_SVDS_PRIVATE_H
#define PRIMME_SVDS_PRIVATE_H

static void convTestFunAugmented(double *eval, void *evec, double *rNorm, int *isConv,
   primme_params *primme);
static void convTestFunATA(double *eval, void *evec, double *rNorm, int *isConv,
   primme_params *primme);
static void copy_params_from_svds(primme_svds_params *primme_svds, int stage);
static void globalSumDoubleSvds(void *sendBuf, void *recvBuf, int *count, 
                         primme_params *primme);

#endif
