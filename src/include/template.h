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
 * File: template.h
 *
 * Purpose - Contains definitions of macros used along PRIMME.
 *
 ******************************************************************************/

#ifndef TEMPLATE_H
#define TEMPLATE_H

#include <limits.h>    
#include <float.h>
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
#  define MACHINE_EPSILON DBL_EPSILON
#elif defined(USE_DOUBLECOMPLEX)
#  define SCALAR_PRE z
#  define REAL_PRE d
#  define SCALAR_SUF zprimme
#  define REAL_SUF dprimme
#  define USE_COMPLEX
#  define SCALAR PRIMME_COMPLEX_DOUBLE
#  define REAL double
#  define MACHINE_EPSILON DBL_EPSILON
#elif defined(USE_FLOAT)
#  define SCALAR_PRE s
#  define REAL_PRE s
#  define SCALAR_SUF sprimme
#  define REAL_SUF sprimme
#  define SCALAR float
#  define REAL float
#  define MACHINE_EPSILON FLT_EPSILON
#elif defined(USE_FLOATCOMPLEX)
#  define SCALAR_PRE c
#  define REAL_PRE s
#  define SCALAR_SUF cprimme
#  define REAL_SUF sprimme
#  define USE_COMPLEX
#  define SCALAR PRIMME_COMPLEX_FLOAT
#  define REAL float
#  define MACHINE_EPSILON FLT_EPSILON
#else
#  error "An arithmetic should be selected, please define one of USE_DOUBLE, USE_DOUBLECOMPLEX, USE_FLOAT or USE_FLOATCOMPLEX."
#endif

/* A C99 code with complex type is not a valid C++ code. However C++          */
/* compilers usually can take it. Nevertheless in order to avoid the warnings */
/* while compiling in pedantic mode, we use the proper complex type for C99   */
/* (complex double and complex float) and C++ (std::complex<double> and       */
/* std::complex<float>). Of course both complex types are binary compatible.  */

#ifdef USE_COMPLEX
#  ifndef __cplusplus
#     define REAL_PART(x) (creal(x))
#     define ABS(x) (cabs(x))
#     define CONJ(x) (conj(x))
#  else
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
   /* Avoid to use the final type for complex in generated headers file.      */
   /* Instead use PRIMME_COMPLEX_FLOAT and _DOUBLE.                           */
#  undef PRIMME_COMPLEX_FLOAT
#  undef PRIMME_COMPLEX_DOUBLE
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
   int __err = (ERRN); assert(__err==0);\
   if (__err) {\
      if (primme->printLevel > 0 && primme->outputFile) \
         fprintf(primme->outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", __err, __LINE__, #ERRN );\
      return (RETURN);\
   }\
}

/**********************************************************************
 * Macro CHKERRNOABORT - If ERRN != 0, it is printed out in
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

#define CHKERRNOABORT(ERRN, RETURN) { \
   int __err = (ERRN);\
   if (__err) {\
      if (primme->printLevel > 0 && primme->outputFile) \
         fprintf(primme->outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", __err, __LINE__, #ERRN );\
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
   int __err = (ERRN); assert(__err==0);\
   if (__err) {\
      if (primme->printLevel > 0 && primme->outputFile) {\
         fprintf(primme->outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", __err, __LINE__, #ERRN );\
         fprintf(primme->outputFile, "PRIMME: " __VA_ARGS__);\
         fprintf(primme->outputFile, "\n");\
      }\
      return (RETURN);\
   }\
}

/**********************************************************************
 * Macro CHKERRNOABORTM - If ERRN == 0, it is printed out in
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
 *   CHKERRNOABORTM((ptr = malloc(n*sizeof(double))) == NULL, -1,
 *        "malloc could not allocate %d doubles\n", n);
 *
 **********************************************************************/

#define CHKERRNOABORTM(ERRN, RETURN, ...) { \
   int __err = (ERRN);\
   if (__err) {\
      if (primme->printLevel > 0 && primme->outputFile) {\
         fprintf(primme->outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", __err, __LINE__, #ERRN );\
         fprintf(primme->outputFile, "PRIMME: " __VA_ARGS__);\
         fprintf(primme->outputFile, "\n");\
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
   int __err = (ERRN); assert(__err==0);\
   if (__err) {\
      if (primme_svds->printLevel > 0 && primme_svds->outputFile) {\
         fprintf(primme_svds->outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", __err, __LINE__, #ERRN );\
      }\
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
   int __err = (ERRN); assert(__err==0);\
   if (__err) {\
      if (primme_svds->printLevel > 0 && primme_svds->outputFile) {\
         fprintf(primme_svds->outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", __err, __LINE__, #ERRN );\
         fprintf(primme_svds->outputFile, "PRIMME: " __VA_ARGS__);\
         fprintf(primme_svds->outputFile, "\n");\
      }\
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


/**********************************************************************
 * Macro WRKSP_MALLOC_PRIMME - borrow NELEM of type **X from workspace *rwork;
 *    on return *X is an aligned pointer to the borrowed space and workspace
 *    *rwork points to the first free address.
 *
 * INPUT/OUTPUT PARAMETERS
 * -----------------------
 * NELEM   Number of sizeof(**X) to allocate
 * X       Where to store the pointer of the borrowed space
 * RWORK   Reference to the workspace from take space; *RWORK is updated
 *         with the first free address
 * LRWORK  Reference to the number of elements in *RWORK; *LRWORK is updated
 *         with the left number of elements in *RWORK
 *
 * RETURN VALUE
 * ------------
 * error code
 *
 * EXAMPLE
 * -------
 *   double *values, *rwork; size_t *rworkSize;
 *   CHKERR(WRKSP_MALLOC_PRIMME(n, &values, &rwork, &rworkSize), -1);
 *
 **********************************************************************/

#define WRKSP_MALLOC_PRIMME(NELEM, X, RWORK, LRWORK) (\
      *(uintptr_t*)(X) = ALIGN_BY_SIZE(*(RWORK), sizeof(**(X))), \
      *(uintptr_t*)(RWORK) = ALIGN_BY_SIZE(*(X)+(NELEM), sizeof(**(RWORK))), \
      /* Check that there are enough elements in *RWORK */ \
      /* NOTE: the check is pessimistic */ \
      (sizeof(**(X))*((NELEM)+1) + sizeof(**(RWORK)) - 2 \
        <= *(LRWORK)*sizeof(**(RWORK))) \
        /* If there is, subtract the used number of elements */ \
        ? (*(LRWORK) -= (sizeof(**(X))*((NELEM)+1) + sizeof(**(RWORK)) - 2) \
                           / sizeof(**(RWORK)), 0 /* return success */)\
        /* Else, return error */ \
        : -1)

#define ALIGN_BY_SIZE(ptr,size_of_T) (((uintptr_t)(ptr)+(size_of_T)-1) & -(size_of_T))

#define ALIGN(ptr, T) (T*) ALIGN_BY_SIZE(ptr, sizeof(T))


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
