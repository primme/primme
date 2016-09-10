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

#include "numerical.h"
#include "globalsum.h"
#include "wtime.h"

TEMPLATE_PLEASE
int globalSum_Sprimme(SCALAR *sendBuf, SCALAR *recvBuf, int count, 
      primme_params *primme) {

   int ierr;
   double t0=0.0;

   if (primme && primme->globalSumReal) {
      t0 = primme_wTimer(0);

      /* If it is a complex type, count real and imaginary part */
#ifdef USE_COMPLEX
      count *= 2;
#endif
      CHKERRM((primme->globalSumReal(sendBuf, recvBuf, &count, primme, &ierr),
               ierr), -1,
            "Error returned by 'globalSumReal' %d", ierr);

      primme->stats.timeGlobalSum += primme_wTimer(0) - t0;
      primme->stats.volumeGlobalSum += count;
   }
   else {
      Num_copy_Sprimme(count, sendBuf, 1, recvBuf, 1);
   }

   return 0;
}
