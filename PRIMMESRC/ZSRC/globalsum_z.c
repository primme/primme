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
 * File: globalsum.c
 *
 * Purpose - Wrappers around primme->globalSumDouble
 *
 ******************************************************************************/

#include "primme.h"
#include "numerical_z.h"

void globalSum_zprimme(complex double *sendBuf, complex double *recvBuf, int count, 
      primme_params *primme) {

   count *= 2;

   if (primme && primme->globalSumDouble) {
      primme->globalSumDouble(sendBuf, recvBuf, &count, primme);
   }
   else {
      Num_copy_dprimme(count, (double*)sendBuf, 1, (double*)recvBuf, 1);
   }
}
