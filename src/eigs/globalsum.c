/*******************************************************************************
 * Copyright (c) 2017, College of William & Mary
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the College of William & Mary nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * PRIMME: https://github.com/primme/primme
 * Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
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
