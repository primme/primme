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
 * File: template.h
 *
 * Purpose - Contains definitions of macros used along PRIMME.
 *
 ******************************************************************************/

#ifndef TEMPLATE_H
#define TEMPLATE_H

#include <limits.h>    
#include "primme.h"

/*****************************************************************************/
/* Arithmetic                                                                */
/*****************************************************************************/

/**********************************************************************
 * Macros USE_FLOAT, USE_FLOATCOMPLEX, USE_DOUBLE and USE_DOUBLECOMPLEX -
 *    only one of them is defined at the same time, and identifies the
 *    type of SCALAR, one of float, complex float, double or complex double.
 *
 * Macro SCALAR - type of the matrices problem and the eigenvectors.
 *
 * Macro REAL - type of the real part of SCALAR, and the type of the
 *    eigenvalues.
 *
 * Macro SCALAR_PRE - prefix used in Xprimme C interface function (dprimme,
 *    zprimme...) and in BLAS/LAPACK functions.
 *
 * Macro REAL_PRE - prefix used in BLAS/LAPACK functions.
 *
 * Macro SCALAR_SUF - suffix appended in functions whose prototypes have
 *    SCALAR and REAL arguments and is called from other files.
 *
 * Macro REAL_SUF - suffix appended to call the REAL version of the
 *    function.
 *
 * Macro USE_COMPLEX - only defined when SCALAR is a complex type.
 **********************************************************************/

#if defined(USE_DOUBLE)
#  define SCALAR_PRE d
#  define REAL_PRE d
#  define SCALAR_SUF dprimme
#  define REAL_SUF dprimme
#  define SCALAR double
#  define REAL double
#elif defined(USE_DOUBLECOMPLEX)
#  define SCALAR_PRE z
#  define REAL_PRE d
#  define SCALAR_SUF zprimme
#  define REAL_SUF dprimme
#  define USE_COMPLEX
#  define REAL double
#elif defined(USE_FLOAT)
#  define SCALAR_PRE s
#  define REAL_PRE s
#  define SCALAR_SUF sprimme
#  define REAL_SUF sprimme
#  define SCALAR float
#  define REAL float
#elif defined(USE_FLOATCOMPLEX)
#  define SCALAR_PRE c
#  define REAL_PRE s
#  define SCALAR_SUF cprimme
#  define REAL_SUF sprimme
#  define USE_COMPLEX
#  define REAL float
#endif

/* A C99 code with complex type is not a valid C++ code. However C++          */
/* compilers usually can take it. Nevertheless in order to avoid the warnings */
/* while compiling in pedantic mode, we use the proper complex type for C99   */
/* (complex double and complex float) and C++ (std::complex<double> and       */
/* std::complex<float>). Of course both complex types are binary compatible.  */

#ifdef USE_COMPLEX
#  ifndef __cplusplus
#     include <complex.h> /* definition of creal, cabs, conj */
#     define SCALAR complex REAL
#     define REAL_PART(x) (creal(x))
#     define ABS(x) (cabs(x))
#     define CONJ(x) (conj(x))
#  else
#     include <complex> /* definition of real, abs, conj */
#     define SCALAR std::complex<REAL>
#     define REAL_PART(x) (std::real(x))
#     define ABS(x) (std::abs(x))
#     define CONJ(x) (std::conj(x))
#  endif
#else
#  define REAL_PART(x) (x)
#  define ABS(x) (fabs(x))
#  define CONJ(x) (x)
#endif

/* complex.h may be defined in primme.h or here; so undefine I */
#ifdef I
#   undef I
#endif

#ifndef __cplusplus
#include <tgmath.h>   /* select proper function abs from fabs, cabs... */
#endif

#define MACHINE_EPSILON 1.11e-16

/* TEMPLATE_PLEASE tags the functions whose prototypes depends on macros and  */
/* are used in other files. The macro has value only when the tool ctemplate  */
/* will inspect the source files, which happens when the macro CHECK_TEMPLATE */
/* is defined. See Makefile and tools/ctemplate.                              */
/*                                                                            */
/* When SCALAR is not a complex type (e.g., float and double) the function it */
/* will be referred as _Sprimme and _Rprimme. Otherwise it will be referred   */
/* only as _Sprimme. The term TEMPLATE_PLEASE should prefix every function    */
/* that will be instantiated with different values for SCALAR and REAL.       */

#ifdef CHECK_TEMPLATE
#  ifdef USE_DOUBLE
#     define TEMPLATE_PLEASE \
        APPEND_FUNC(Sprimme,SCALAR_SUF) USE(Sprimme,"SCALAR_SUF") USE(Rprimme,"REAL_SUF")
#  else
#     define TEMPLATE_PLEASE \
        APPEND_FUNC(Sprimme,SCALAR_SUF)
#  endif
#else
#  define TEMPLATE_PLEASE
#endif

/*****************************************************************************/
/* Error management                                                          */
/*****************************************************************************/

/* The error code is checked with an assert. Then by default the code will   */
/* abort in case of error, unless it is compiled with the macro NDEBUG       */
/* defined. In that case it is printed out the file and the line in the      */
/* source that caused the error in the file descriptor primme->outputFile    */
/* or primme_svds->outputFile. Use CHKERR and CHKERRM in functions that has  */
/* primme_params* primme, and CHKERRS and CHKERRMS in functions that has     */
/* primme_svds_params* primme_svds.                                          */

#include <assert.h>

/**********************************************************************
 * Macro CHKERR - assert that ERRN == 0. If not it is printed out in
 *    primme->outputFile the file and the line of the caller, and
 *    force the caller function to return RETURN.
 *
 *    ERRN is only evaluated once.
 *
 * INPUT PARAMETERS
 * ----------------
 * ERRN    Expression that returns an error code
 * RETURN  Value that the caller function will return in case of error
 *
 **********************************************************************/

#define CHKERR(ERRN, RETURN) { \
   int err = (ERRN); assert(err==0);\
   if (err) {\
      if (primme->outputFile) \
         fprintf(primme->outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", err, __LINE__, #ERRN );\
      return (RETURN);\
   }\
}

/**********************************************************************
 * Macro CHKERRM - assert that ERRN == 0. If not it is printed out in
 *    primme->outputFile the file and the line of the caller, in
 *    addition to an error message passed as arguments. As CHKERR
 *    the caller function is forced to return RETURN in case of error.
 *
 *    ERRN is only evaluated once.
 *
 * INPUT PARAMETERS
 * ----------------
 * ERRN    Expression that returns an error code
 * RETURN  Value that the caller function will return in case of error
 *
 * EXAMPLE
 * -------
 *   CHKERRM((ptr = malloc(n*sizeof(double))) == NULL, -1,
 *        "malloc could not allocate %d doubles\n", n);
 *
 **********************************************************************/

#define CHKERRM(ERRN, RETURN, ...) { \
   int err = (ERRN); assert(err==0);\
   if (err) {\
      if (primme->outputFile) {\
         fprintf(primme->outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", err, __LINE__, #ERRN );\
         fprintf(primme->outputFile, "PRIMME: " __VA_ARGS__);\
      }\
      return (RETURN);\
   }\
}

/**********************************************************************
 * Macro CHKERRS - assert that ERRN == 0. If not it is printed out in
 *    primme_svds->outputFile the file and the line of the caller, and
 *    force the caller function to return RETURN.
 *
 *    ERRN is only evaluated once.
 *
 * INPUT PARAMETERS
 * ----------------
 * ERRN    Expression that returns an error code
 * RETURN  Value that the caller function will return in case of error
 *
 **********************************************************************/

#define CHKERRS(ERRN, RETURN) { \
   int err = (ERRN); assert(err==0);\
   if (err) {\
      fprintf(primme_svds->outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", err, __LINE__, #ERRN );\
      return (RETURN);\
   }\
}

/**********************************************************************
 * Macro CHKERRMS - assert that ERRN == 0. If not it is printed out in
 *    primme->outputFile the file and the line of the caller, in
 *    addition to an error message passed as arguments. As CHKERRS
 *    the caller function is forced to return RETURN in case of error.
 *
 *    ERRN is only evaluated once.
 *
 * INPUT PARAMETERS
 * ----------------
 * ERRN    Expression that returns an error code
 * RETURN  Value that the caller function will return in case of error
 *
 * EXAMPLE
 * -------
 *   CHKERRMS((ptr = malloc(n*sizeof(double))) == NULL, -1,
 *        "malloc could not allocate %d doubles\n", n);
 *
 **********************************************************************/

#define CHKERRMS(ERRN, RETURN, ...) { \
   int err = (ERRN); assert(err==0);\
   if (err) {\
      fprintf(primme_svds->outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", err, __LINE__, #ERRN );\
      fprintf(primme_svds->outputFile, "PRIMME: " __VA_ARGS__);\
      return (RETURN);\
   }\
}

/*****************************************************************************/
/* Memory management                                                         */
/*****************************************************************************/

/**********************************************************************
 * Macro MALLOC_PRIMME - malloc NELEM of type **X and write down the
 *    pointer in *X.
 *
 * INPUT PARAMETERS
 * ----------------
 * NELEM   Number of sizeof(**X) to allocate
 * X       Where to store the pointer
 *
 * RETURN VALUE
 * ------------
 * error code
 *
 * EXAMPLE
 * -------
 *   double *values;
 *   CHKERRM(MALLOC_PRIMME(n, &values), -1,
 *        "malloc could not allocate %d doubles\n", n);
 *
 **********************************************************************/

#if !(defined (__APPLE__) && defined (__MACH__))
#  include <malloc.h> /* malloc */
#endif
#include <stdlib.h>   /* malloc, free */

#define MALLOC_PRIMME(NELEM, X) (*((void**)X) = malloc((NELEM)*sizeof(**(X))), *(X) == NULL)

/*****************************************************************************/
/* Miscellanea                                                               */
/*****************************************************************************/

/* Used by macros emitted by ctemplate */
#define CONCAT(a, b) CONCATX(a, b)
#define CONCATX(a, b) a ## b

#define TO_INT(X) ((X) < INT_MAX ? (X) : INT_MAX)

#ifdef F77UNDERSCORE
#define FORTRAN_FUNCTION(X) X ## _
#else
#define FORTRAN_FUNCTION(X) X
#endif

#ifndef max
#  define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#  define min(a, b) ((a) < (b) ? (a) : (b))
#endif

#endif
