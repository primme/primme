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
 *******************************************************************************
 * File: primme_f77_private.h
 *
 * Purpose - Definitions used exclusively by primme_f77.c
 *
 ******************************************************************************/

#ifndef PRIMME_F77_PRIVATE_H
#define PRIMME_F77_PRIVATE_H

#include "template.h"
#include "primme_interface.h" /* for Sprimme */

/* Prototypes for Fortran-C interface */

#ifdef __cplusplus
extern "C" {
#endif

#define AS_FORTRAN(X) AS_FORTRANX(X)
#define AS_FORTRANX(X) FORTRAN_FUNCTION(X ## _f77)

void AS_FORTRAN(Sprimme)(REAL *evals, SCALAR *evecs,
      REAL *rnorms, primme_params **primme, int *ierr);

/* Only define these functions ones */
#ifdef USE_DOUBLE
void AS_FORTRAN(primme_initialize)(primme_params **primme);
void AS_FORTRAN(primme_free)(primme_params **primme);
void AS_FORTRAN(primme_display_params)(primme_params **primme);
void AS_FORTRAN(primme_printstacktrace)(primme_params **primme);
void AS_FORTRAN(primme_set_method)(primme_params **primme, primme_preset_method *method, int *returnValue);
void AS_FORTRAN(primme_set_member)(primme_params **primme, int *label, void *ptr, int *ierr);
void AS_FORTRAN(primme_get_prec_shift)(primme_params *primme, int *i, double *shift);
void AS_FORTRAN(primme_get_member)(primme_params *primme, int *label, void *ptr, int *ierr);
void AS_FORTRAN(primmetop_get_member)(primme_params **primme, int *label, void *ptr, int *ierr);
void AS_FORTRAN(primmetop_get_prec_shift)(primme_params **primme, int *i, double *shift);
#endif

#ifdef __cplusplus
}
#endif

#endif /* PRIMME_F77_PRIVATE_H */
