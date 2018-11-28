/*******************************************************************************
 * Copyright (c) 2018, College of William & Mary
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
 *    In short source files, *.c, are compiled several times, every time for a
 *    different type, referred as SCALAR. Examples of types are float, double,
 *    complex float, complex double, and the corresponding GPU versions. All
 *    types are described in section "Arithmetic". Other macros are defined to
 *    refer derived types and functions. An example is REAL, which is
 *    defined as the real (non-complex) version of SCALAR. For instance,
 *    when SCALAR is complex double, REAL is double. Similarly, it is possible
 *    to call the real version of a function. For instance,
 *    permute_vecs_Sprimme permutes vectors with type SCALAR and
 *    permute_vecs_Rprimme permutes vectors with type REAL.
 *
 *    When SCALAR is a GPU type, the pointers SCALAR* are supposed to point out
 *    memory allocated on GPUs, also called devices. For instance,
 *    Num_malloc_Sprimme allocates GPU memory when SCALAR is a GPU type. Use
 *    HSCALAR and HREAL as the non-GPU, also called host, versions of SCALAR and
 *    REAL. Also to use the non-GPU version use the suffices _SHprimme and
 *    _RHprimme. For instance,
 *    Num_malloc_SHprimme allocates memory on the host, and
 *    permute_vecs_RHprimme permute REAL vectors on the host.
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
#include <stdlib.h>   /* malloc, free */
#include <string.h>   /* strlen */
#include <stdio.h>    /* snprintf */
#include "primme.h"

/*****************************************************************************/
/* Arithmetic                                                                */
/*****************************************************************************/

/* Fake types for MAGMA. Compiler will complain when assigning a value and  */
/* add, subtract, multiply and divide with a GPU type. Also when a GPU type */
/* argument is passed when a non-GPU is expected.                           */

typedef struct {double a;}  magma_double;
typedef struct {float a;}  magma_float;
typedef struct {PRIMME_COMPLEX_DOUBLE a;} magma_complex_double;
typedef struct {PRIMME_COMPLEX_FLOAT a;} magma_complex_float;

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

/* exp and log macros cause warnings on some systems. PRIMME only uses  */
/* the double version of these functions.                               */

#ifdef pow
#  undef pow
#endif
#ifdef exp
#  undef exp
#endif
#ifdef log
#  undef log
#endif

#ifndef __cplusplus
#  define ISFINITE isfinite
#else
#  define ISFINITE std::isfinite
#endif

/* TEMPLATE_PLEASE tags the functions whose prototypes depends on macros and  */
/* are used in other files. The macro has value only when the tool ctemplate  */
/* is inspecting the source files, which is indicated by the macro            */
/* CHECK_TEMPLATE begin defined. See Makefile and tools/ctemplate.            */
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
/* Memory management                                                         */
/*****************************************************************************/

/* Recent versions of PRIMME are using dynamic memory to manage the memory.  */
/* In general the use of dynamic memory simplifies the code by not needing   */
/* to take care of providing enough working space for all functions in       */
/* PRIMME. The small drawback of dynamic memory is to mingle with error      */
/* management. C++ shines in this situation. But we are restricted to C,     */
/* probably, not for good reasons. The goal is to avoid writing specific code*/
/* to free allocated memory in case of an error happening in the body of a   */
/* function. The next macros together with the functions in memman.c         */
/* provide an interface to track memory allocations and free them in case    */
/* of error. These macros and functions are going to be used mostly by       */
/* error management macros, eg CHKERR, and memory allocation functions, eg   */
/* Num_malloc_Sprimme. The error manager is going to assume that every       */
/* allocation is going to be freed at the end of the function. If            */
/* that is not the case, call Mem_keep_frame(ctx). Examples of this are the  */
/* functions Num_malloc_Sprimme.                                             */

#include "memman.h"

/**********************************************************************
 * Macro MEM_PUSH_FRAME - Set a new frame in which new allocations are going
 *    to be registered.
 *
 * NOTE: only one call is allowed in the same block of code (anything between
 *       curly brackets).
 **********************************************************************/

#define MEM_PUSH_FRAME \
   primme_frame __frame = {NULL, 0, ctx.mm}; \
   ctx.mm = &__frame;

/**********************************************************************
 * Macro MEM_POP_FRAME(ERRN) - If ERRN is 0, it just removes all registered
 *    allocations. Otherwise it also frees all allocations.
 *
 * NOTE: only one call is allowed in the same block of code (anything between
 *       curly brackets).
 **********************************************************************/

#ifndef NDEBUG
#define MEM_POP_FRAME(ERRN) { \
   if (ERRN) {\
      Mem_pop_clean_frame(ctx);\
   } else {\
      Mem_debug_frame(__FILE__ ": " STR(__LINE__), ctx);\
      Mem_pop_frame(&ctx); \
   }\
}
#else
#define MEM_POP_FRAME(ERRN) { \
   if (ERRN) {\
      Mem_pop_clean_frame(ctx);\
   } else {\
      Mem_pop_frame(&ctx); \
   }\
}
#endif

/*****************************************************************************/
/* Reporting                                                                 */
/*****************************************************************************/

/**********************************************************************
 * Macro PRINTFALLCTX - report some information. It invokes ctx.report if the
 *    level >= CTX.printLevel regardless of the processor's id.
 *
 * INPUT PARAMETERS
 * ----------------
 * L     print level
 * ...   printf arguments
 *
 * EXAMPLE
 * -------
 *   PRINTFALLCTX(ctx, 2, "Alpha=%f is to close to zero", alpha);
 *
 **********************************************************************/

#define PRINTFALLCTX(CTX, L, ...) { \
   if ((CTX).report && (L) <= (CTX).printLevel) { \
      int len = snprintf(NULL, 0, __VA_ARGS__)+1; \
      char *str = malloc(len); \
      snprintf(str, len, __VA_ARGS__); \
      (CTX).report(str, -1.0, (CTX)); \
      free(str); \
   }\
}

/**********************************************************************
 * Macro PRINTFALL - report some information. It invokes ctx.report if the
 *    level >= ctx.printLevel regardless of the processor's id.
 *
 * INPUT PARAMETERS
 * ----------------
 * L     print level
 * ...   printf arguments
 *
 * EXAMPLE
 * -------
 *   PRINTFALL(2, "Alpha=%f is to close to zero", alpha);
 *
 **********************************************************************/

#define PRINTFALL(L, ...) PRINTFALLCTX(ctx, (L), __VA_ARGS__)

/**********************************************************************
 * Macro PRINTF - report some information. It invokes ctx.report if the
 *    level >= ctx.printLevel and it is processor 0.
 *
 * INPUT PARAMETERS
 * ----------------
 * L     print level
 * ...   printf arguments
 *
 * EXAMPLE
 * -------
 *   PRINTF(2, "Alpha=%f is to close to zero", alpha);
 *
 **********************************************************************/

#define PRINTF(L, ...) {if (ctx.procID == 0) PRINTFALL((L), __VA_ARGS__);}

/*****************************************************************************/
/* Profiling                                                                 */
/*****************************************************************************/

#include "wtime.h"

#ifdef PRIMME_PROFILE

#include <regex.h>

static inline const char *__compose_function_name(const char *path,
      const char *call, const char *filename, const char *line) {

   int len = strlen(path) + 1 + 45 + strlen(filename) + 1 + strlen(line) + 10;
   char *s = (char*)malloc(len);
   snprintf(s, len-1, "%s~%.40s@%s:%s", path, call, filename, line);
   s[len-1] = 0;
   return s;
}

#define PROFILE_BEGIN(CALL) \
   double ___t0 = 0; \
   const char *old_path = ctx.path; \
   if (ctx.path && ctx.report) { \
      ctx.path = __compose_function_name(ctx.path, CALL, __FILE__, STR(__LINE__)); \
      if (regexec(&ctx.profile, ctx.path, 0, NULL, 0) == 0) { \
            ctx.report(ctx.path, -.5, ctx); \
         ___t0 = primme_wTimer() - *ctx.timeoff; \
      } \
   }

#define PROFILE_END \
   if (ctx.path) { \
      if (___t0 > 0) { \
         double ___t1 = primme_wTimer(); \
         ___t0 = ___t1 - ___t0 - *ctx.timeoff; \
         if (___t0 > 0 && ctx.report) { \
            ctx.report(ctx.path, ___t0, ctx); \
            *ctx.timeoff += primme_wTimer() - ___t1; \
         } \
      } \
      free((void*)ctx.path); \
      ctx.path = old_path; \
   }
 
#else
#define PROFILE_BEGIN(CALL)
#define PROFILE_END
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
 * Macro CHKERR - assert that ERRN == 0. Otherwise it is printed out in
 *    ctx.outputFile the file and the line of the caller, and
 *    forced the caller function to return ERRN.
 *
 *    ERRN is only evaluated once.
 *
 * INPUT PARAMETERS
 * ----------------
 * ERRN    Expression that returns an error code
 *
 **********************************************************************/

#define CHKERR(ERRN) { \
   MEM_PUSH_FRAME; \
   PROFILE_BEGIN(STR(ERRN)); \
   int __err = (ERRN); assert(__err==0);\
   PROFILE_END; \
   MEM_POP_FRAME(__err); \
   if (__err) {\
      PRINTFALL(1, "PRIMME: Error %d in (" __FILE__ ":%d): %s", __err, __LINE__, #ERRN );\
      return __err;\
   }\
}

/**********************************************************************
 * Macro CHKERRM - assert that ERRN == 0. Otherwise it is printed out in
 *    ctx.outputFile the file and the line of the caller, in
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
   MEM_PUSH_FRAME; \
   PROFILE_BEGIN(STR(ERRN)); \
   int __err = (ERRN); assert(__err==0);\
   PROFILE_END; \
   MEM_POP_FRAME(__err); \
   if (__err) {\
      PRINTFALL(1, "PRIMME: Error %d in (" __FILE__ ":%d): %s", __err, __LINE__, #ERRN );\
      PRINTFALL(1, "PRIMME: " __VA_ARGS__);\
      return (RETURN);\
   }\
}


/**********************************************************************
 * Macro CHKERRA - assert that ERRN == 0. Otherwise it is printed out in
 *    ctx.outputFile the file and the line of the caller, and
 *    forced the caller function to execute ACTION.
 *
 *    ERRN is only evaluated once.
 *
 * INPUT PARAMETERS
 * ----------------
 * ERRN    Expression that returns an error code
 * ACTION  Expression that the caller function will execute in case of error
 *
 **********************************************************************/

#define CHKERRA(ERRN, ACTION) { \
   MEM_PUSH_FRAME; \
   PROFILE_BEGIN(STR(ERRN)); \
   int __err = (ERRN); assert(__err==0);\
   PROFILE_END; \
   MEM_POP_FRAME(__err); \
   if (__err) {\
      PRINTFALL(1, "PRIMME: Error %d in (" __FILE__ ":%d): %s\n", __err, __LINE__, #ERRN );\
      ACTION;\
      return;\
   } \
}

/*****************************************************************************/
/* Memory, error and params  management                                      */
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

#define MALLOC_PRIMME(NELEM, X)                                                \
   (*((void **)X) = malloc((NELEM) * sizeof(**(X))),                           \
         *(X) == NULL ? PRIMME_MALLOC_FAILURE : 0)

typedef struct primme_context_str {
   /* For PRIMME */
   primme_params *primme;
   primme_svds_params *primme_svds;

   /* For output */
   int printLevel;
   FILE *outputFile;
   int (*report)(const char *fun, double time, struct primme_context_str ctx);

   /* For memory management */
   primme_frame *mm;

   /* for MPI */
   int numProcs;     /* number of processes */
   int procID;       /* process id */
   void *mpicomm;    /* MPI communicator */

   /* For MAGMA */
   void *queue;      /* magma device queue (magma_queue_t*) */

   #ifdef PRIMME_PROFILE
   /* For profiling */
   regex_t profile;  /* Pattern of the functions to profile */
   const char *path; /* Path of the current function */
   double *timeoff;  /* Accumulated overhead of profiling and reporting */
   #endif
} primme_context;

/*****************************************************************************/
/* Miscellanea                                                               */
/*****************************************************************************/

/* Used by macros emitted by ctemplate */
#define CONCAT(a, b) CONCATX(a, b)
#define CONCATX(a, b) a ## b
#define STR(X) STR0(X)
#define STR0(X) #X

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
