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

TEMPLATE_PLEASE
int globalSum_Sprimme(SCALAR *sendBuf, SCALAR *recvBuf, int count, 
      primme_params *primme) {

   if (primme && primme->globalSumDouble) {
      /* If it is a complex type, count real and imaginary part */
      #ifdef USE_COMPLEX
         count *= 2;
      #endif
      #if defined(USE_DOUBLE) || defined(USE_DOUBLECOMPLEX)
         primme->globalSumDouble((double*)sendBuf, (double*)recvBuf, &count,
               primme);
      #else
         double *sendBufd, *recvBufd;
         int i;
         sendBufd = (double*)malloc(sizeof(double)*count*2);
         recvBufd = sendBufd + count;
         for(i=0; i<count; i++) {
            sendBufd[i] = ((REAL*)sendBuf)[i];
         }
         primme->globalSumDouble(sendBufd, recvBufd, &count, primme);
         for(i=0; i<count; i++) {
            ((REAL*)recvBuf)[i] = recvBufd[i];
         }
         free(sendBufd);
      #endif
   }
   else {
      Num_copy_Sprimme(count, sendBuf, 1, recvBuf, 1);
   }

   return 0;
}
