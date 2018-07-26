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

/* Including MAGMA headers before C99 complex.h avoid compiler issues */
#if !defined(CHECK_TEMPLATE) && (defined(USE_FLOAT_MAGMA) || defined(USE_FLOATCOMPLEX_MAGMA) || defined(USE_DOUBLE_MAGMA) || defined(USE_DOUBLECOMPLEX_MAGMA))
#  include <magma_v2.h>
#endif

#include <limits.h>    
#include <float.h>
#include <stdint.h>
#include "primme.h"

/*****************************************************************************/
/* Arithmetic                                                                */
/*****************************************************************************/

/* Fake types for MAGMA */

typedef double magma_double;
typedef float magma_float;
typedef PRIMME_COMPLEX_DOUBLE magma_complex_double;
typedef PRIMME_COMPLEX_FLOAT magma_complex_float;

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
 * Macro HOST_SCALAR_SUF and HOST_REAL_SUF - suffix appended to call the
 *    CPU version.
 *
 * Macro HSCALAR and HREAL - cpu versions of the types.
 *
 * Macro USE_COMPLEX - only defined when SCALAR is a complex type.
 *
 * Macro USE_HOST - CPU version
 *
 * Macro USE_MAGMA - MAGMA version
 **********************************************************************/

#if defined(USE_DOUBLE)
#  define SCALAR_PRE d
#  define REAL_PRE d
#  define SCALAR_SUF dprimme
#  define REAL_SUF dprimme
#  define SCALAR double
#  define REAL double
#  define MACHINE_EPSILON DBL_EPSILON
#  define USE_HOST
#elif defined(USE_DOUBLECOMPLEX)
#  define SCALAR_PRE z
#  define REAL_PRE d
#  define SCALAR_SUF zprimme
#  define REAL_SUF dprimme
#  define USE_COMPLEX
#  define SCALAR PRIMME_COMPLEX_DOUBLE
#  define REAL double
#  define MACHINE_EPSILON DBL_EPSILON
#  define USE_HOST
#elif defined(USE_FLOAT)
#  define SCALAR_PRE s
#  define REAL_PRE s
#  define SCALAR_SUF sprimme
#  define REAL_SUF sprimme
#  define SCALAR float
#  define REAL float
#  define MACHINE_EPSILON FLT_EPSILON
#  define USE_HOST
#elif defined(USE_FLOATCOMPLEX)
#  define SCALAR_PRE c
#  define REAL_PRE s
#  define SCALAR_SUF cprimme
#  define REAL_SUF sprimme
#  define USE_COMPLEX
#  define SCALAR PRIMME_COMPLEX_FLOAT
#  define REAL float
#  define MACHINE_EPSILON FLT_EPSILON
#  define USE_HOST

/* MAGMA types */

#elif defined(USE_DOUBLE_MAGMA)
#  define SCALAR_PRE magma_d
#  define REAL_PRE magma_d
#  define SCALAR_SUF dmagmaprimme
#  define REAL_SUF dmagmaprimme
#  define HOST_SCALAR_SUF dprimme
#  define HOST_REAL_SUF dprimme
#  define SCALAR magma_double
#  define REAL magma_double
#  define HSCALAR double
#  define HREAL double
#  define MAGMA_SCALAR double
#  define MAGMA_REAL double
#  define MACHINE_EPSILON DBL_EPSILON
#  define USE_MAGMA
#elif defined(USE_DOUBLECOMPLEX_MAGMA)
#  define SCALAR_PRE magma_z
#  define REAL_PRE magma_d
#  define SCALAR_SUF zmagmaprimme
#  define REAL_SUF dmagmaprimme
#  define HOST_SCALAR_SUF zprimme
#  define HOST_REAL_SUF dprimme
#  define USE_COMPLEX
#  define SCALAR magma_complex_double
#  define REAL magma_double
#  define HSCALAR PRIMME_COMPLEX_DOUBLE
#  define HREAL double
#  define MAGMA_SCALAR magmaDoubleComplex
#  define MAGMA_REAL double
#  define MACHINE_EPSILON DBL_EPSILON
#  define USE_MAGMA
#elif defined(USE_FLOAT_MAGMA)
#  define SCALAR_PRE magma_s
#  define REAL_PRE magma_s
#  define SCALAR_SUF smagmaprimme
#  define REAL_SUF smagmaprimme
#  define HOST_SCALAR_SUF sprimme
#  define HOST_REAL_SUF sprimme
#  define SCALAR magma_float
#  define REAL magma_float
#  define HSCALAR float
#  define HREAL float
#  define MAGMA_SCALAR float
#  define MAGMA_REAL float
#  define MACHINE_EPSILON FLT_EPSILON
#  define USE_MAGMA
#elif defined(USE_FLOATCOMPLEX_MAGMA)
#  define SCALAR_PRE magma_c
#  define REAL_PRE magma_s
#  define SCALAR_SUF cmagmaprimme
#  define REAL_SUF smagmaprimme
#  define HOST_SCALAR_SUF cprimme
#  define HOST_REAL_SUF sprimme
#  define USE_COMPLEX
#  define SCALAR magma_complex_float
#  define REAL magma_float
#  define HSCALAR PRIMME_COMPLEX_FLOAT
#  define HREAL float
#  define MAGMA_SCALAR magmaFloatComplex
#  define MAGMA_REAL float
#  define MACHINE_EPSILON FLT_EPSILON
#  define USE_MAGMA
#else
#  error "An arithmetic should be selected, please define one of USE_DOUBLE, USE_DOUBLECOMPLEX, USE_FLOAT or USE_FLOATCOMPLEX."
#endif

#ifndef HOST_SCALAR_SUF
#  define HOST_SCALAR_SUF SCALAR_SUF
#endif
#ifndef HOST_REAL_SUF
#  define HOST_REAL_SUF REAL_SUF
#endif
#ifndef HSCALAR
#  define HSCALAR SCALAR
#endif
#ifndef HREAL
#  define HREAL REAL
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

#ifndef __cplusplus
#  define ISFINITE isfinite
#else
#  define ISFINITE std::isfinite
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
        APPEND_FUNC(Sprimme,SCALAR_SUF) USE(Sprimme,"SCALAR_SUF") USE(Rprimme,"REAL_SUF") USE(SHprimme,"HOST_SCALAR_SUF") USE(RHprimme,"HOST_REAL_SUF")
#  else
#     define TEMPLATE_PLEASE \
        APPEND_FUNC(Sprimme,SCALAR_SUF)
#  endif
   /* Avoid to use the final type for integers and complex in generated       */
   /* headers file. Instead use PRIMME_COMPLEX_FLOAT and _DOUBLE.             */
#  undef PRIMME_COMPLEX_FLOAT
#  undef PRIMME_COMPLEX_DOUBLE
#  undef PRIMME_INT
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

#define CHKERR(ERRN) { \
   int __err = (ERRN); assert(__err==0);\
   if (__err) {\
      if (ctx.printLevel > 0 && ctx.outputFile) \
         fprintf(ctx.outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", __err, __LINE__, #ERRN );\
      return __err;\
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
      if (ctx.printLevel > 0 && ctx.outputFile) {\
         fprintf(ctx.outputFile, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", __err, __LINE__, #ERRN );\
         fprintf(ctx.outputFile, "PRIMME: " __VA_ARGS__);\
         fprintf(ctx.outputFile, "\n");\
      }\
      return (RETURN);\
   }\
}

/*****************************************************************************/
/* Memory, error and params  management                                      */
/*****************************************************************************/

typedef struct {
   /* For PRIMME */
   primme_params *primme;
   primme_svds_params *primme_svds;

   /* For output */
   int printLevel;
   FILE *outputFile;

   /* for MPI */
   int numProcs;     /* number of processes */
   int procID;       /* process id */
   void *mpicomm;    /* MPI communicator */

   /* For MAGMA */
   void *queue;   	/* magma device queue (magma_queue_t*) */
} primme_context;

/*****************************************************************************/
/* Miscellanea                                                               */
/*****************************************************************************/

/* Used by macros emitted by ctemplate */
#define CONCAT(a, b) CONCATX(a, b)
#define CONCATX(a, b) a ## b

#define TO_INT(X) ((X) < INT_MAX ? (X) : INT_MAX)

#ifdef F77UNDERSCORE
#define FORTRAN_FUNCTION(X) CONCAT(X,_)
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
