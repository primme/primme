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
 * File: numerical.c
 *
 * Purpose - This file contains C wrapper routines for certain numerical 
 *           routines to be used by programs of either precision.
 *
 ******************************************************************************/

#include <stdarg.h>
#include "common_numerical_private.h"
#include "common_numerical.h"

/******************************************************************************/
void Num_dcopy_primme(int n, double *x, int incx, double *y, int incy) {
   PRIMME_BLASINT ln = n;
   PRIMME_BLASINT lincx = incx;
   PRIMME_BLASINT lincy = incy;

   DCOPY(&ln, x, &lincx, y, &lincy);
}

/******************************************************************************/
double Num_dlamch_primme(const char *cmach) {
#ifdef NUM_CRAY
   _fcd cmach_fcd;

   cmach_fcd = _cptofcd(cmach, strlen(cmach));
   return (DLAMCH(cmach_fcd));
#else
   return (DLAMCH(cmach));
#endif

}

/******************************************************************************/

int Num_imax_primme(int numArgs, int val1, int val2, ...) {

   int maxVal, nextVal;
   va_list argPtr;

   va_start(argPtr, val2);

   if (numArgs < 2) {
      numArgs = 2;
   }

   maxVal = (val1 > val2 ? val1 : val2);

   for (numArgs -= 2; numArgs; numArgs--) {
      nextVal = va_arg(argPtr, int);

      maxVal = (nextVal > maxVal ? nextVal : maxVal);
   }

   va_end(argPtr);

   return maxVal;

}


double Num_fmin_primme(int numArgs, double val1, double val2, ...) {

   double minVal, nextVal;
   va_list argPtr;

   va_start(argPtr, val2);

   if (numArgs < 2) {
      numArgs = 2;
   }

   minVal = (val1 < val2 ? val1 : val2);

   for (numArgs -= 2; numArgs; numArgs--) {
      nextVal = va_arg(argPtr, double);

      minVal = (nextVal < minVal ? nextVal : minVal);
   }

   va_end(argPtr);

   return minVal;

}


double Num_fmax_primme(int numArgs, double val1, double val2, ...) {

   double maxVal, nextVal;
   va_list argPtr;

   va_start(argPtr, val2);

   if (numArgs < 2) {
      numArgs = 2;
   }

   maxVal = (val1 > val2 ? val1 : val2);

   for (numArgs -= 2; numArgs; numArgs--) {
      nextVal = va_arg(argPtr, double);

      maxVal = (nextVal > maxVal ? nextVal : maxVal);
   }

   va_end(argPtr);

   return maxVal;

}

