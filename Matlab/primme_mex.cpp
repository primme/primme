/*******************************************************************************
 * Copyright (c) 2016, College of William & Mary
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
 * File: primme_mex.c
 * 
 * Purpose - PRIMME MEX interface.
 * 
 * Currently we don't support distributed computing. The PRIMME MEX
 * reads the inputs from MATLAB, constructs the necessary structures
 * and then calls PRIMME. The desired results are then returned to MATLAB.
 * 
 * Matrix-vector and preconditioning functions are performed by callbacks
 * to MATLAB functions.
 *
 * For details about PRIMME parameters, methods, and settings see ../readme.txt
 *
 ******************************************************************************/

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <complex>
#include <cassert>
#include "mex.h"
#include "primme.h"
#ifdef USE_GPUARRAY
#  include "gpu/mxGPUArray.h"
#  include <cuda_runtime.h>
#  include "magma_v2.h"
#endif

// Attempt to capture ctrl+c
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__)) || defined (__FreeBSD__)
#include <signal.h>
#if (defined (__APPLE__) && defined (__MACH__)) || defined (__FreeBSD__)
# define SIGHANDLER_T sig_t
#else
# define SIGHANDLER_T sighandler_t
#endif

static volatile int keepRunning = 1;
void interrumptHandler(int sig) {
   keepRunning = 0;
}
#endif

#ifndef macro_max
#define macro_max(a, b) ((a) > (b) ? (a) : (b))
#endif

////////////////////////////////////////////////////////////////////////////////
//
// Auxiliary functions
//

// Allow 1 extra argument as printf, e.g.,
// mexErrMsgTxtPrintf1("The argument %d is empty", argumentNumber)

#define mexErrMsgTxtPrintf1(A, B) \
   {\
      char buffer[500]; \
      snprintf(buffer, 500, A, B); \
      mexErrMsgTxt(buffer); \
   }

// Allow 2 extra arguments as printf, e.g.,
// mexErrMsgTxtPrintf2("The matrix should be of size %dx%d ", m, n)

#define mexErrMsgTxtPrintf2(A, B, C) \
   {\
      char buffer[500]; \
      snprintf(buffer, 500, A, B, C); \
      mexErrMsgTxt(buffer); \
   }

// Return the corresponding REAL type of a given type, for instance
// Real<std::complex<double> >::type is double and
// Real<float>::type is float.

template <typename T>
struct Real { typedef T type; };
template <>
struct Real<std::complex<double> > { typedef double type; };
template <>
struct Real<std::complex<float> > { typedef float type; };

// Return whether the type is complex, for instance
// isComplex<double>() is false and
// isComplex<std::complex<float> >() is true.

template <typename T>
static inline bool isComplex (void) { return false; }
template <>
bool isComplex<std::complex<float> > (void) { return true; }
template <>
bool isComplex<std::complex<double> > (void) { return true; }

// Return the mxClassID corresponding to a C type

template <typename T>
static inline mxClassID toClassID (void);
template <>
mxClassID toClassID<float>(void) { return mxSINGLE_CLASS; }
template <>
mxClassID toClassID<std::complex<float> >(void) { return mxSINGLE_CLASS; }
template <>
mxClassID toClassID<double>(void) { return mxDOUBLE_CLASS; }
template <>
mxClassID toClassID<std::complex<double> >(void) { return mxDOUBLE_CLASS; }
template <>
mxClassID toClassID<PRIMME_INT>(void) { return mxINT64_CLASS; }
template <>
mxClassID toClassID<int>(void) { return mxINT32_CLASS; }

// Return the size of an element of mxClassID

size_t getSizeMxClassId(mxClassID id) {
   switch(id) {
   case mxDOUBLE_CLASS: return sizeof(double);
   case mxSINGLE_CLASS: return sizeof(float);
   case mxINT8_CLASS: return (size_t)1;
   case mxUINT8_CLASS: return (size_t)1;
   case mxINT16_CLASS: return (size_t)2;
   case mxUINT16_CLASS: return (size_t)2;
   case mxINT32_CLASS: return (size_t)4;
   case mxUINT32_CLASS: return (size_t)4;
   case mxINT64_CLASS: return (size_t)8;
   case mxUINT64_CLASS: return (size_t)8;
   case mxLOGICAL_CLASS: return sizeof(mxLogical);
   default:
      mexErrMsgTxt("Unsupported matrix type");
   }
   return 0;
}

// Template cpu/gpu

struct CPU {
   CPU() {}
};
struct GPU {
   GPU() {}
};
static bool isCPU(CPU) { return true; }
static bool isGPU(CPU) { return false; }
static bool isCPU(GPU) { return false; }
static bool isGPU(GPU) { return true; }


// Auxiliary function for copy_mxArray, copy the content of a mxArray with
// compatible C type TX to a C array of type TY.
// Arguments:
// - x: MATLAB array from copy the values
// - y: C type array to copy the values
// - m: number of rows of matrix y and x
// - n: number of columns of matrix y and x
// - ldy: leading dimension of y

template <typename TX, typename TY, typename I>
static void copy_mxArray_kernel(const mxArray *x, TY *y, I m, I n, I ldy) {
   if (!mxIsComplex(x)) {
      TX *px = (TX*)mxGetData(x);
      for (I i=0; i<n; i++) for (I j=0; j<m; j++) y[ldy*i+j] = px[m*i+j];
   }
   else {
      mexErrMsgTxt("Unsupported matrix type");
   }
}

template <typename TX, typename TY, typename I>
static void copy_mxArray_kernel(const mxArray *x, std::complex<TY> *y, I m, I n, I ldy) {
   if (!mxIsComplex(x)) {
      TX *px = (TX*)mxGetData(x);
      for (I i=0; i<n; i++) for (I j=0; j<m; j++) y[ldy*i+j] = px[m*i+j];
   }
   else {
      TX *pxr = (TX*)mxGetData(x), *pxi = (TX*)mxGetImagData(x);
      for (I i=0; i<n; i++) for (I j=0; j<m; j++)
         y[ldy*i+j] = std::complex<TY>(pxr[m*i+j], pxi[m*i+j]);
   }
}

// Copy the content of a mxArray to a C array
// Arguments:
// - x: MATLAB array from copy the values
// - y: C type array to copy the values
// - m: number of rows of matrix y and x
// - n: number of columns of matrix y and x
// - ldy: leading dimension of y

template <typename TY, typename I>
static void copy_mxArray(const mxArray *x, TY *y, I m, I n, I ldy, CPU) {

   /* Check dimensions */

   const mwSize *dims = mxGetDimensions(x);
   mwSize ndims = mxGetNumberOfDimensions(x);
   assert(ldy >= m);
   if (!((ndims == 1 && dims[0] == m && n == 1)
            || (ndims == 2 && dims[0] == 1 && dims[1] == m && n == 1)
            || (ndims == 2 && dims[0] == m && dims[1] == n))) {
      mexErrMsgTxtPrintf2("Unsupported matrix dimension; it should be %dx%d",
            (int)m, (int)n);
      return;
   }

   /* Do the copy */

   switch(mxGetClassID(x)) {
   case mxDOUBLE_CLASS:
      copy_mxArray_kernel<double,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   case mxSINGLE_CLASS:
      copy_mxArray_kernel<float,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   case mxINT8_CLASS:
      copy_mxArray_kernel<char,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   case mxUINT8_CLASS:
      copy_mxArray_kernel<unsigned char,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   case mxINT16_CLASS:
      copy_mxArray_kernel<short,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   case mxUINT16_CLASS:
      copy_mxArray_kernel<unsigned short,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   case mxINT32_CLASS:
      copy_mxArray_kernel<int,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   case mxUINT32_CLASS:
      copy_mxArray_kernel<unsigned int,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   case mxINT64_CLASS:
      copy_mxArray_kernel<int64_T,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   case mxUINT64_CLASS:
      copy_mxArray_kernel<uint64_T,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   case mxLOGICAL_CLASS:
      copy_mxArray_kernel<mxLogical,typename Real<TY>::type,I>(x, y, m, n, ldy);
      break;
   default:
      mexErrMsgTxt("Unsupported matrix type");
   }      
}

#ifdef USE_GPUARRAY

template <typename TY, typename I>
static void copy_mxArray(const mxArray *x, TY *y, I m, I n, I ldy, GPU) {

   /* Get gpuArray */

   mxGPUArray const *xgpu = mxGPUCreateFromMxArray(x);

   /* Check dimensions */

   const mwSize *dims = mxGPUGetDimensions(xgpu);
   mwSize ndims = mxGPUGetNumberOfDimensions(xgpu);
   assert(ldy >= m);
   if (!((ndims == 1 && dims[0] == m && n == 1)
            || (ndims == 2 && dims[0] == 1 && dims[1] == m && n == 1)
            || (ndims == 2 && dims[0] == m && dims[1] == n))) {
      mexErrMsgTxtPrintf2("Unsupported matrix dimension; it should be %dx%d",
            (int)m, (int)n);
      return;
   }

   /* Check that the size of x and y are the same */

   if (getSizeMxClassId(mxGPUGetClassID(xgpu)) *
               (mxGPUGetComplexity(xgpu) == mxREAL ? 1 : 2) !=
         sizeof(TY)) {
      mexErrMsgTxt("Not supported to copy matrices of different types");
   }

   /* Do the copy */

   void const *xdata = mxGPUGetDataReadOnly(xgpu);
   if (cudaMemcpy2D(y, ldy * sizeof(TY), xdata, m * sizeof(TY), m * sizeof(TY),
             n, cudaMemcpyDeviceToDevice) != cudaSuccess) {
      mexErrMsgTxt("Error copying data");
   }

   /* Destroy gpuArray */

   mxGPUDestroyGPUArray(xgpu);
}
#endif /* USE_GPUARRAY */

// Creates a mxArray with the content of a C array
// Arguments:
// - y: C type array from to get the values
// - m: number of rows of matrix y and output mxArray
// - n: number of columns of matrix y and output mxArray
// - ldy: leading dimension of y
// - avoidCopy: if true, try to use y as the data of the output mxArray

template <typename T, typename I>
static mxArray* create_mxArray(T *y, I m, I n, I ldy, CPU, bool avoidCopy=false) {

   if (isComplex<T>())
         mexErrMsgTxt("This should not happen");

   // We avoid to copy when the T isn't complex and the leading dimension of
   // y is the number of row. This trick only works in Octave, MATLAB requires
   // the data in mxArray being created with mxMalloc, etc.

#ifdef HAVE_OCTAVE
   if (avoidCopy && m == ldy) {
      // Create an empty mxArray to avoid MATLAB to allocate space

      mxArray *x = mxCreateNumericMatrix(0, 0, toClassID<T>(), mxREAL);

      // Set data and dimensions of the new mxArray

      mxSetData(x, y);
      mxSetM(x, (mwSize)macro_max(m, 0));
      mxSetN(x, (mwSize)macro_max(n, 0));

      return x;
   }
   else
#endif
   {
      // Create mxArray with dimensions m x n and proper complexity

      mxArray *x = mxCreateNumericMatrix((mwSize)macro_max(m, 0),
            (mwSize)macro_max(n, 0), toClassID<T>(), mxREAL);

      // Copy the content of y into the mxArray
      if (y) {
         T *px = (T *)mxGetData(x);
         for (I i = 0; i < n; i++)
            for (I j = 0; j < m; j++) px[m * i + j] = y[ldy * i + j];
      }

      return x;
   }
}

template <typename T, typename I>
static mxArray* create_mxArray(std::complex<T> *y, I m, I n, I ldy,
      CPU, bool avoidCopy=false) {

   // When y is complex isn't possible to save the copy
   (void)avoidCopy;

   // Create mxArray with dimensions m x n and proper complexity

   mxArray *x = mxCreateNumericMatrix((mwSize)macro_max(m, 0),
         (mwSize)macro_max(n, 0), toClassID<T>(), mxCOMPLEX);

   // Copy the content of y into the mxArray

   if (y) {
      T *pxr = (T *)mxGetData(x), *pxi = (T *)mxGetImagData(x);
      for (I i = 0; i < n; i++)
         for (I j = 0; j < m; j++) {
            pxr[m * i + j] = std::real(y[ldy * i + j]);
            pxi[m * i + j] = std::imag(y[ldy * i + j]);
         }
   }

   return x;
}

#ifdef USE_GPUARRAY

template <typename T, typename I>
static mxArray* create_mxArray(T *y, I m, I n, I ldy, GPU, bool avoidCopy=false) {

   m = macro_max(m, 0);
   n = macro_max(n, 0);

   // Create mxArray with dimensions m x n and proper complexity

   mwSize dims[2] = {(mwSize)m, (mwSize)n};
   mxGPUArray *xgpu = mxGPUCreateGPUArray(2, dims, toClassID<T>(),
         isComplex<T>() ? mxCOMPLEX : mxREAL, MX_GPU_DO_NOT_INITIALIZE);

   // Copy the content of y into the gpuArray

   if (y) {
      void *xdata = mxGPUGetData(xgpu);
      if (cudaMemcpy2D(xdata, m * sizeof(T), y, ldy * sizeof(T), m * sizeof(T),
                n, cudaMemcpyDeviceToDevice) != cudaSuccess) {
         mexErrMsgTxt("Error copying data");
      }
   }

   mxArray* x = mxGPUCreateMxArrayOnGPU(xgpu);
   mxGPUDestroyGPUArray(xgpu);

   return x;
}

#endif /* GPU_ARRAY */

static mxArray *create_mxArray(
      const char *y, bool avoidCopy) {

   (void)avoidCopy;
   return mxCreateString(y ? y : "");
}

// Template version of sprimme, cprimme, dprimme and zprimme

static int tprimme(float *evals, float *evecs, float *resNorms, primme_params *primme, CPU) {
      return sprimme(evals, evecs, resNorms, primme);
}
static int tprimme(float *evals, std::complex<float> *evecs, float *resNorms, primme_params *primme, CPU) {
      return cprimme(evals, evecs, resNorms, primme);
}
static int tprimme(double *evals, double *evecs, double *resNorms, primme_params *primme, CPU) {
      return dprimme(evals, evecs, resNorms, primme);
}
static int tprimme(double *evals, std::complex<double> *evecs, double *resNorms, primme_params *primme, CPU) {
      return zprimme(evals, evecs, resNorms, primme);
}
static int tprimme(std::complex<float> *evals, std::complex<float> *evecs, float *resNorms, primme_params *primme, CPU) {
      return cprimme_normal(evals, evecs, resNorms, primme);
}
static int tprimme(std::complex<double> *evals, std::complex<double> *evecs, double *resNorms, primme_params *primme, CPU) {
      return zprimme_normal(evals, evecs, resNorms, primme);
}
#ifdef USE_GPUARRAY
static int tprimme(float *evals, float *evecs, float *resNorms, primme_params *primme, GPU) {
      return magma_sprimme(evals, evecs, resNorms, primme);
}
static int tprimme(float *evals, std::complex<float> *evecs, float *resNorms, primme_params *primme, GPU) {
      return magma_cprimme(evals, evecs, resNorms, primme);
}
static int tprimme(double *evals, double *evecs, double *resNorms, primme_params *primme, GPU) {
      return magma_dprimme(evals, evecs, resNorms, primme);
}
static int tprimme(double *evals, std::complex<double> *evecs, double *resNorms, primme_params *primme, GPU) {
      return magma_zprimme(evals, evecs, resNorms, primme);
}
static int tprimme(std::complex<float> *evals, std::complex<float> *evecs, float *resNorms, primme_params *primme, GPU) {
      return magma_cprimme_normal(evals, evecs, resNorms, primme);
}
static int tprimme(std::complex<double> *evals, std::complex<double> *evecs, double *resNorms, primme_params *primme, GPU) {
      return magma_zprimme_normal(evals, evecs, resNorms, primme);
}
#endif /* GPU_ARRAY */


// Template version of sprimme_svds, cprimme_svds, dprimme_svds and zprimme_svds

static int tprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds, CPU) { 
   return sprimme_svds(svals, svecs, resNorms, primme_svds);
}
static int tprimme_svds(float *svals, std::complex<float> *svecs, float *resNorms, primme_svds_params *primme_svds, CPU) {
   return cprimme_svds(svals, svecs, resNorms, primme_svds);
}
static int tprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds, CPU) { 
   return dprimme_svds(svals, svecs, resNorms, primme_svds);
}
static int tprimme_svds(double *svals, std::complex<double> *svecs, double *resNorms, primme_svds_params *primme_svds, CPU) {
   return zprimme_svds(svals, svecs, resNorms, primme_svds);
}
#ifdef USE_GPUARRAY
static int tprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds, GPU) { 
   return magma_sprimme_svds(svals, svecs, resNorms, primme_svds);
}
static int tprimme_svds(float *svals, std::complex<float> *svecs, float *resNorms, primme_svds_params *primme_svds, GPU) {
   return magma_cprimme_svds(svals, svecs, resNorms, primme_svds);
}
static int tprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds, GPU) { 
   return magma_dprimme_svds(svals, svecs, resNorms, primme_svds);
}
static int tprimme_svds(double *svals, std::complex<double> *svecs, double *resNorms, primme_svds_params *primme_svds, GPU) {
   return magma_zprimme_svds(svals, svecs, resNorms, primme_svds);
}
#endif /* GPU_ARRAY */


// Select the function based on primme_op_datatype

template <typename T, typename F, typename G,
      template <typename, typename, typename> class S,
      typename R = typename S<T, F, G>::t>
typename S<T, F, G>::t select_fun(primme_op_datatype t) {
   if (!isComplex<T>()) {
      switch (t) {
      case primme_op_default: return S<T,F,G>::f;
      case primme_op_double:  return S<double,F,G>::f;
      case primme_op_float:   return S<float,F,G>::f;
      default: return nullptr;
      }
   } else {
      switch (t) {
      case primme_op_default: return S<T,F,G>::f;
      case primme_op_double:  return S<std::complex<double>,F,G>::f;
      case primme_op_float:   return S<std::complex<float>,F,G>::f;
      default: return nullptr;
      }
   }
}


// Check that there are N input arguments in a MATLAB function

#define ASSERT_NUMARGSIN(N) \
   if (nrhs != (N)) { \
      mexErrMsgTxtPrintf1("Required %d arguments\n", (N)+1); \
   }

// Check that there are N output arguments in a MATLAB function

#define ASSERT_NUMARGSOUT(N) \
   if (nlhs != (N)) { \
      mexErrMsgTxtPrintf1("Required %d output arguments", (N)+1); \
   }

// Check that there are N output arguments at least in a MATLAB function

#define ASSERT_NUMARGSOUTGE(N) \
   if (nlhs <= (N)) { \
      mexErrMsgTxtPrintf1("Required %d output arguments at least", (N)+1); \
   }

// Check that argument NARG is compatible with a pointer in a MATLAB function

#define ASSERT_POINTER(NARG) \
   if (sizeof(void*)==4) { \
      if (!mxIsUint32(prhs[(NARG)])) { \
         mexErrMsgTxtPrintf1("Argument %d should be uint32", (NARG)+2); \
      } \
   } \
   if (sizeof(void*)==8) { \
      if (!mxIsUint64(prhs[(NARG)])) { \
         mexErrMsgTxtPrintf1("Argument %d should be uint64", (NARG)+2); \
      } \
   }

// Check that argument NARG is compatible with a number in a MATLAB function

#define ASSERT_NUMERIC(NARG) \
   if (!mxIsNumeric(prhs[(NARG)])) { \
      mexErrMsgTxtPrintf1("Argument %d should be numeric", (NARG)+2); \
   }

// Check that argument NARG is a function handler in a MATLAB function

#define ASSERT_FUNCTION(NARG) \
   if (mxGetClassID(prhs[(NARG)]) != mxFUNCTION_CLASS) { \
      mexErrMsgTxtPrintf1("Argument %d should be function handler", (NARG)+2); \
   }

// Check that argument NARG is compatible with a number/string in a MATLAB function

#define ASSERT_NUMERIC_OR_CHAR(NARG) \
   if (!mxIsNumeric(prhs[(NARG)]) && !mxIsChar(prhs[(NARG)])) { \
      mexErrMsgTxtPrintf1("Argument %d should be numeric or string", (NARG)+2); \
   }

// Check that C returns zero, if not print error

#define CHKERR(C) {\
   int __ierr = (C); \
   if (__ierr != 0) { \
      mexErrMsgTxtPrintf1("Call '" #C "' returned %d", __ierr); \
   } \
}

// Return a mxArray wrapping a pointer up

static mxArray* mxArrayFromPointer(void *p) {
   mxArray *a = mxCreateNumericMatrix(1, 1,
         sizeof(void*)==4?mxUINT32_CLASS:mxUINT64_CLASS, mxREAL);
   *(void**)mxGetData(a) = p;
   return a;
}

// Return the pointer wrapped up in a mxArray

static void* mxArrayToPointer(const mxArray *a) {
   return *(void**)mxGetData(a);
}

// Return a label

static primme_params_label mxArrayToLabel(const mxArray *a) {
   if (mxIsChar(a)) {
      const char *label_name = mxArrayToString(a);
      primme_params_label label = (primme_params_label)-1;
      CHKERR(primme_member_info(&label, &label_name, NULL, NULL));
      return label;
   }
   else if (mxIsNumeric(a)) {
      return (primme_params_label)mxGetScalar(a);
   }
   else {
      mexErrMsgTxt("Not valid value for a label");
   }
   return (primme_params_label)0;
}

static primme_svds_params_label mxArrayToLabelSvds(const mxArray *a) {
   if (mxIsChar(a)) {
      const char *label_name = mxArrayToString(a);
      primme_svds_params_label label = (primme_svds_params_label)-1;
      CHKERR(primme_svds_member_info(&label, &label_name, NULL, NULL));
      return label;
   }
   else if (mxIsNumeric(a)) {
      return (primme_svds_params_label)mxGetScalar(a);
   }
   else {
      mexErrMsgTxt("Not valid value for a label");
   }
   return (primme_svds_params_label)0;
}

// Return a constant

static int mxArrayToConstant(const mxArray *a) {
   if (mxIsChar(a)) {
      const char *constant_name = mxArrayToString(a);
      int constant;
      if(primme_constant_info(constant_name, &constant)) {
         mexErrMsgTxt("Not valid constant");
         return -1;
      }
      return constant;
   }
   else if (mxIsNumeric(a)) {
      return (int)mxGetScalar(a);
   }
   else {
      mexErrMsgTxt("Not valid constant");
   }
   return -1;
}

static int mxArrayToConstantSvds(const mxArray *a) {
   if (mxIsChar(a)) {
      const char *constant_name = mxArrayToString(a);
      int constant;
      if(primme_svds_constant_info(constant_name, &constant)) {
         mexErrMsgTxt("Not valid constant");
         return -1;
      }
      return constant;
   }
   else if (mxIsNumeric(a)) {
      return (int)mxGetScalar(a);
   }
   else {
      mexErrMsgTxt("Not valid constant");
   }
   return -1;
}

////////////////////////////////////////////////////////////////////////////////
//
// MATLAB wrappers around functions in PRIMME interface
//

// Wrapper around primme_initialize; prototype:
// [primme] = mexFunction_primme_initialize()

static void mexFunction_primme_initialize(int nlhs, mxArray *plhs[], int nrhs,
      const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(0);
   ASSERT_NUMARGSOUT(1);

   primme_params *primme = new primme_params;
   primme_initialize(primme);
   plhs[0] = mxArrayFromPointer(primme);
}

// Wrapper around primme_set_method; prototype:
// mexFunction_primme_set_method(primme_preset_method, primme)

static void mexFunction_primme_set_method(int nlhs, mxArray *plhs[], int nrhs,
      const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(2);
   ASSERT_NUMARGSOUT(0);
   ASSERT_NUMERIC_OR_CHAR(0);
   ASSERT_POINTER(1);

   primme_preset_method method = (primme_preset_method)mxArrayToConstant(prhs[0]);
   primme_params *primme = (primme_params*)mxArrayToPointer(prhs[1]);
   CHKERR(primme_set_method(method, primme));
}

// Wrapper around primme_free; prototype:
// mexFunction_primme_free(primme)

static void mexFunction_primme_free(int nlhs, mxArray *plhs[], int nrhs,
      const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(1);
   ASSERT_NUMARGSOUT(0);
   ASSERT_POINTER(0);

   primme_params *primme = (primme_params*)mxArrayToPointer(prhs[0]);
   if (primme->targetShifts) delete [] primme->targetShifts;
   if (primme->matrix) mxDestroyArray((mxArray*)primme->matrix);
   if (primme->massMatrix) mxDestroyArray((mxArray*)primme->massMatrix);
   if (primme->preconditioner) mxDestroyArray((mxArray*)primme->preconditioner);
   if (primme->convtest) mxDestroyArray((mxArray*)primme->convtest);
   if (primme->monitor) mxDestroyArray((mxArray*)primme->monitor);
   if (primme->commInfo) mxDestroyArray((mxArray*)primme->commInfo);
#ifdef USE_GPUARRAY
   if (primme->queue) {
      // and finalize MAGMA.
      magma_queue_destroy(*(magma_queue_t*)primme->queue);
      free(primme->queue);
      magma_finalize();
   }
#endif

   primme_free(primme);
   delete primme;
}

// Wrapper around primme_set_member; prototype:
// mexFunction_primme_set_member(primme, label, value)

static void mexFunction_primme_set_member(int nlhs, mxArray *plhs[], int nrhs,
      const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(3);
   ASSERT_NUMARGSOUT(0);
   ASSERT_POINTER(0);
   ASSERT_NUMERIC_OR_CHAR(1);

   primme_params *primme = (primme_params*)mxArrayToPointer(prhs[0]);
   primme_params_label label = mxArrayToLabel(prhs[1]);

   switch(label) {
      // Set members with arity > 1

      case PRIMME_iseed:
      {
         ASSERT_NUMERIC(2);
         copy_mxArray(prhs[2], primme->iseed, 4, 1, 4, CPU());
         break;
      }

      case PRIMME_targetShifts:
      {
         ASSERT_NUMERIC(2);
         if (primme->targetShifts) delete [] primme->targetShifts;
         int n = (int)mxGetNumberOfElements(prhs[2]);
         primme->targetShifts = new double[n];
         primme->numTargetShifts = n;
         copy_mxArray(prhs[2], primme->targetShifts, n, 1, n, CPU());
         break;
      }

      // The function handlers are stored in the user data fields in
      // primme_params, e.g., in matrix for matrixMatvec and preconditioner
      // for applyPreconditioner

      case PRIMME_matrixMatvec:
      {
         ASSERT_FUNCTION(2);
         if (primme->matrix) mxDestroyArray((mxArray*)primme->matrix);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme->matrix = (void*)a;
         break;
      }
      case PRIMME_applyPreconditioner:
      {
         ASSERT_FUNCTION(2);
         if (primme->preconditioner) mxDestroyArray((mxArray*)primme->preconditioner);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme->preconditioner = (void*)a;
         break;
      }
      case PRIMME_massMatrixMatvec:
      {
         ASSERT_FUNCTION(2);
         if (primme->massMatrix) mxDestroyArray((mxArray*)primme->massMatrix);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme->massMatrix = (void*)a;
         break;
      }
      case PRIMME_convTestFun:
      {
         ASSERT_FUNCTION(2);
         if (primme->convtest) mxDestroyArray((mxArray*)primme->convtest);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme->convtest = (void*)a;
         break;
      }
      case PRIMME_monitorFun:
      {
         ASSERT_FUNCTION(2);
         if (primme->monitor) mxDestroyArray((mxArray*)primme->monitor);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme->monitor = (void*)a;
         break;
      }

      case PRIMME_commInfo:
      {
         ASSERT_NUMERIC_OR_CHAR(2);
         if (primme->commInfo) mxDestroyArray((mxArray*)primme->commInfo);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme->commInfo = (void*)a;
         break;
      }


      // Forbidden members
 
      case PRIMME_numProcs:
      case PRIMME_procID:
      case PRIMME_nLocal:
      case PRIMME_globalSumReal:
      case PRIMME_numTargetShifts:
      case PRIMME_outputFile:
      case PRIMME_matrix:
      case PRIMME_preconditioner:
      case PRIMME_ldevecs:
      case PRIMME_ldOPs:
      case PRIMME_monitor:
      case PRIMME_convtest:
         mexErrMsgTxt("Unsupported to set this option");
         break;

      default : 
      {
         primme_type ptype;
         int arity;

         // Get type of the option
         int ierr = primme_member_info(&label, NULL, &ptype, &arity);
         if (ierr != 0) {
            mexErrMsgTxt("Unknown option");
         }
         CHKERR(arity != 1);

         // Set members with type int

         if (ptype == primme_int) {
            PRIMME_INT v = mxArrayToConstant(prhs[2]);
            CHKERR(primme_set_member(primme, label, &v));
         }

         // Set members with type double

         else if (ptype == primme_double) {
            ASSERT_NUMERIC(2);
            double v;
            copy_mxArray(prhs[2], &v, 1, 1, 1, CPU());
            CHKERR(primme_set_member(primme, label, &v));
         }

         // Set members with type string

         else if (ptype == primme_string) {
            const char *str = mxArrayToString(prhs[2]);
            CHKERR(primme_set_member(primme, label, (void*)str));
         }

         else {
            /* This shouldn't happen */
            CHKERR(1);
         }
      }
   }
}

// Wrapper around primme_get_member; prototype:
// mexFunction_primme_get_member(primme, label, value)

static void mexFunction_primme_get_member(int nlhs, mxArray *plhs[], int nrhs,
      const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(2);
   ASSERT_NUMARGSOUT(1);
   ASSERT_POINTER(0);
   ASSERT_NUMERIC_OR_CHAR(1);

   primme_params *primme = (primme_params*)mxArrayToPointer(prhs[0]);
   primme_params_label label = mxArrayToLabel(prhs[1]);

   switch(label) {
      // Set members with arity > 1

      case PRIMME_iseed:
      {
         plhs[0] = create_mxArray(primme->iseed, 4, 1, 4, CPU());
         break;
      }

      case PRIMME_targetShifts:
      {
         plhs[0] = create_mxArray(primme->targetShifts, primme->numTargetShifts,
               1, primme->numTargetShifts, CPU());
         break;
      }

      // The function handlers are stored in the user data fields in
      // primme_params, e.g., in matrix for matrixMatvec and preconditioner
      // for applyPreconditioner

      case PRIMME_matrixMatvec:
      {
         plhs[0] = (mxArray*)primme->matrix;
         break;
      }
      case PRIMME_applyPreconditioner:
      {
         plhs[0] = (mxArray*)primme->preconditioner;
         break;
      }
      case PRIMME_massMatrixMatvec:
      {
         plhs[0] = (mxArray*)primme->massMatrix;
         break;
      }
      case PRIMME_convTestFun:
      {
         plhs[0] = (mxArray*)primme->convtest;
         break;
      }
      case PRIMME_monitorFun:
      {
         plhs[0] = (mxArray*)primme->monitor;
         break;
      }
  
      // Forbidden members

      case PRIMME_numProcs:
      case PRIMME_procID:
      case PRIMME_commInfo:
      case PRIMME_nLocal:
      case PRIMME_globalSumReal:
      case PRIMME_numTargetShifts:
      case PRIMME_outputFile:
      case PRIMME_matrix:
      case PRIMME_preconditioner:
      case PRIMME_convtest:
      case PRIMME_ldevecs:
      case PRIMME_ldOPs:
      case PRIMME_monitor:
         mexErrMsgTxt("Unsupported to set this option");
         break;

      default : 
      {
         primme_type ptype;
         int arity;

         // Get type of the option
         int ierr = primme_member_info(&label, NULL, &ptype, &arity);
         if (ierr != 0) {
            mexErrMsgTxt("Unknown option");
         }
         CHKERR(arity != 1);

         // Get members with type int

         if (ptype == primme_int) {
            PRIMME_INT v;
            CHKERR(primme_get_member(primme, label, &v));
            plhs[0] = create_mxArray(&v, 1, 1, 1, CPU());
         }

         // Get members with type double

         else if (ptype == primme_double) {
            double v;
            CHKERR(primme_get_member(primme, label, &v));
            plhs[0] = create_mxArray(&v, 1, 1, 1, CPU());
         }

         // Get members with type string

         else if (ptype == primme_string) {
            const char *v;
            CHKERR(primme_get_member(primme, label, &v));
            plhs[0] = mxCreateString(v);
         }

         else {
            /* This shouldn't happen */
            CHKERR(1);
         }
      }
   }
}

// Auxiliary functions for mexFunction_xprimme; they return the value of
// some option in primme_params

struct getMatrixField {
   static void *get(primme_params *primme) {
      return primme->matrix;
   }
};

struct getPreconditionerField {
   static void* get(primme_params *primme) {
      return primme->preconditioner;
   }
};

struct getMassMatrixField {
   static void *get(primme_params *primme) {
      return primme->massMatrix;
   }
};

// Auxiliary function for mexFunction_xprimme; PRIMME wrapper around
// matrixMatvec, massMatrixMatvec and applyPreconditioner. Create a mxArray
// from input vector x, call the function handler returned by F(primme) and
// copy the content of its returned mxArray into the output vector y.

template <typename T, typename F, typename CPUGPU>
struct matrixMatvecEigs {
   typedef void (*t)(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, struct primme_params *primme, int *ierr);
   static void f(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, struct primme_params *primme, int *ierr)
   {  
      mxArray *prhs[2], *plhs[1];

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
      // Check interrupt handler
      if (!keepRunning) {
         *ierr = 1;
      return;
      }
#endif

      if (*blockSize <= 0) {*ierr = 0; return;}

      // Create input vector x (avoid copy if possible)

      prhs[1] = create_mxArray(
            (T *)x, primme->n, (PRIMME_INT)*blockSize, *ldx, CPUGPU(), true);

      // Call the callback

      prhs[0] = (mxArray*)F::get(primme);
      *ierr = mexCallMATLAB(1, plhs, 2, prhs, "feval");

      // Copy lhs[0] to y and destroy it

      if (plhs[0]) {
         copy_mxArray(plhs[0], (T *)y, primme->n, (PRIMME_INT)*blockSize, *ldy,
               CPUGPU());
         mxDestroyArray(plhs[0]);
      }

      // Destroy prhs[1]

      if (isCPU(CPUGPU()) && mxGetData(prhs[1]) == x) mxSetData(prhs[1], NULL);
      mxDestroyArray(prhs[1]); 
   }
};

template <typename T, typename CPUGPU, typename EVAL>
static void convTestFunEigs(double *eval, void *evec, double *rNorm, int *isConv, 
         struct primme_params *primme, int *ierr)
{  
   mxArray *prhs[4], *plhs[1];

   // Create input vectors (avoid copy if possible)

   prhs[1] = create_mxArray((EVAL *)eval, eval ? 1 : 0, 1, eval ? 1 : 0, CPU());
   prhs[2] = create_mxArray((T *)evec, evec ? primme->nLocal : (PRIMME_INT)0,
         (PRIMME_INT)1, evec ? primme->nLocal : (PRIMME_INT)0, CPUGPU(), true);
   prhs[3] = create_mxArray<double,int>(rNorm, rNorm?1:0, 1, rNorm?1:0, CPU());

   // Call the callback

   prhs[0] = (mxArray*)primme->convtest;
   *ierr = mexCallMATLAB(1, plhs, 4, prhs, "feval");

   // Copy plhs[0] to isConv and destroy it

   if (plhs[0]) {
      copy_mxArray(plhs[0], isConv, 1, 1, 1, CPU());
      mxDestroyArray(plhs[0]);
   }

   // Destroy prhs[1..3]

   if (isCPU(CPUGPU()) && mxGetData(prhs[2]) == evec) mxSetData(prhs[2], NULL);
   for (int i=1; i<4; i++) mxDestroyArray(prhs[i]); 
}

template <typename T, typename EVAL>
static void monitorFunEigs(void *basisEvals, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms, int *numConverged,
      void *lockedEvals, int *numLocked, int *lockedFlags, void *lockedNorms,
      int *inner_its, void *LSRes, const char *msg, double *time,
      primme_event *event, struct primme_params *primme, int *ierr) {

   mxArray *prhs[14];

   // Create input vectors (avoid copy if possible)

   typedef typename Real<T>::type R;
   prhs[1] = create_mxArray<R, int>((EVAL *)basisEvals, basisSize ? *basisSize : 0,
         1, basisSize ? *basisSize : 0, CPU(), true);
   prhs[2] = create_mxArray<int, int>(basisFlags, basisFlags ? *basisSize : 0,
         1, basisFlags ? *basisSize : 0, CPU(), true);
   prhs[3] = create_mxArray<int, int>(iblock, blockSize ? *blockSize : 0, 1,
         blockSize ? *blockSize : 0, CPU(), true);
   prhs[4] = create_mxArray<R, int>((R *)basisNorms, basisSize ? *basisSize : 0,
         1, basisSize ? *basisSize : 0, CPU(), true);
   prhs[5] = create_mxArray<int, int>(
         numConverged, numConverged ? 1 : 0, 1, 1, CPU(), true);
   int numLocked0 = numLocked && *numLocked > 0 ? *numLocked : 0;
   prhs[6] = create_mxArray<R, int>(
         (EVAL *)lockedEvals, numLocked0, 1, numLocked0, CPU(), true);
   prhs[7] = create_mxArray<int, int>(
         lockedFlags, numLocked0, 1, numLocked0, CPU(), true);
   prhs[8] = create_mxArray<R, int>(
         (R *)lockedNorms, numLocked0, 1, numLocked0, CPU(), true);
   prhs[9] = create_mxArray<int, int>(inner_its, inner_its ? 1 : 0, 1, 1, CPU(), true);
   prhs[10] = create_mxArray<R, int>((R *)LSRes, LSRes ? 1 : 0, 1, 1, CPU(), true);
   prhs[11] = create_mxArray(msg, true);
   prhs[12] = create_mxArray<double, int>(time, time ? 1 : 0, 1, 1, CPU(), true);
   prhs[13] = create_mxArray<int, int>((int *)event, event ? 1 : 0, 1, 1, CPU(), true);

   // Call the callback

   prhs[0] = (mxArray*)primme->monitor;
   *ierr = mexCallMATLAB(0, NULL, 14, prhs, "feval");

   // Destroy prhs[1..11]

   if (mxGetData(prhs[1]) == basisEvals)   mxSetData(prhs[1], NULL);
   if (mxGetData(prhs[2]) == basisFlags)   mxSetData(prhs[2], NULL);
   if (mxGetData(prhs[3]) == iblock)       mxSetData(prhs[3], NULL);
   if (mxGetData(prhs[4]) == basisNorms)   mxSetData(prhs[4], NULL);
   if (mxGetData(prhs[5]) == numConverged) mxSetData(prhs[5], NULL);
   if (mxGetData(prhs[6]) == lockedEvals)  mxSetData(prhs[6], NULL);
   if (mxGetData(prhs[7]) == lockedFlags)  mxSetData(prhs[7], NULL);
   if (mxGetData(prhs[8]) == lockedNorms)  mxSetData(prhs[8], NULL);
   if (mxGetData(prhs[9]) == inner_its)    mxSetData(prhs[9], NULL);
   if (mxGetData(prhs[10]) == LSRes)       mxSetData(prhs[10], NULL);
   if (mxGetData(prhs[11]) == msg)         mxSetData(prhs[11], NULL);
   if (mxGetData(prhs[12]) == time)        mxSetData(prhs[12], NULL);
   if (mxGetData(prhs[13]) == event)       mxSetData(prhs[13], NULL);
   for (int i=1; i<14; i++) mxDestroyArray(prhs[i]); 
}


// Wrapper around xprimme; prototype:
// [ret, evals, rnorms, evecs] = mexFunction_xprimme(init_guesses, primme)

template<typename T, typename CPUGPU, typename EVAL>
static void mexFunction_xprimme(int nlhs, mxArray *plhs[], int nrhs,
      const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(2);
   ASSERT_NUMARGSOUTGE(1);
   ASSERT_POINTER(1);

   primme_params *primme = (primme_params*)mxArrayToPointer(prhs[1]);

#ifdef USE_GPUARRAY
   if (isGPU(CPUGPU())) {
      if (!primme->commInfo)
         mexErrMsgTxt("Set primme.commInfo with the GPU index");

      // Initialize the MathWorks GPU API
      mxInitGPU();

      // Initialize MAGMA
      magma_init();

      // Create a context
      int gpuDevice;
      copy_mxArray((mxArray*)primme->commInfo, &gpuDevice, 1, 1, 1, CPU());
      if (primme->queue)
         magma_queue_destroy(*(magma_queue_t *)primme->queue);
      else
         primme->queue = malloc(sizeof(magma_queue_t));
      magma_queue_create(gpuDevice, (magma_queue_t*)primme->queue);
   }
#endif

   // Allocate evals, rnorms and evecs; if possible create the mxArray and use
   // its data

   EVAL *evals;
   typename Real<T>::type *rnorms;
   T *evecs;
   mxArray *mxEvals, *mxRnorms, *mxEvecs = nullptr;

   if (nlhs <= 0 || isComplex<EVAL>()) {
      evals = new EVAL[primme->numEvals];
      mxEvals = NULL;
   }
   else {
      mxEvals = mxCreateNumericMatrix(primme->numEvals, 1, toClassID<EVAL>(),
            mxREAL);
      evals = (EVAL *)mxGetData(mxEvals);
   }

   if (nlhs <= 1) {
      rnorms = new typename Real<T>::type[primme->numEvals];
      mxRnorms = NULL;
   }
   else {
      mxRnorms = mxCreateNumericMatrix(primme->numEvals, 1,
            toClassID<typename Real<T>::type>(), mxREAL);
      rnorms = (typename Real<T>::type*)mxGetData(mxRnorms);
   }

   // Number of columns in matrix evecs
   int nevecs =
         primme->numOrthoConst + macro_max(primme->numEvals, primme->initSize);

#ifdef USE_GPUARRAY
   mxGPUArray *mxgpuEvecs  = nullptr; 
   if (isGPU(CPUGPU())) {
      mwSize dims[2] = {
            (mwSize)macro_max(primme->n, 0), (mwSize)macro_max(nevecs, 0)};
      mxgpuEvecs = mxGPUCreateGPUArray(2, dims, toClassID<T>(),
            isComplex<T>() ? mxCOMPLEX : mxREAL, MX_GPU_DO_NOT_INITIALIZE);
      evecs = (T*)mxGPUGetData(mxgpuEvecs);
   } else 
#endif
   if (nlhs <= 3 || isComplex<T>() || primme->numOrthoConst > 0) {
      evecs = new T[primme->n * nevecs];
      mxEvecs = NULL;
   }
   else {
      mxEvecs =
            mxCreateNumericMatrix(primme->n, nevecs, toClassID<T>(), mxREAL);
      evecs = (T*)mxGetData(mxEvecs);
   }

   // Copy initial vectors

   if (primme->numOrthoConst + primme->initSize > 0) {
      ASSERT_NUMERIC(0);
      copy_mxArray(prhs[0], evecs, primme->n,
            (PRIMME_INT)primme->numOrthoConst + primme->initSize, primme->n,
            CPUGPU());
   }

   // Set matvec and preconditioner and monitorFun and convTestFun

   primme->matrixMatvec = select_fun<T, getMatrixField, CPUGPU, matrixMatvecEigs>(
         primme->matrixMatvec_type);
   if (primme->massMatrix) {
      primme->massMatrixMatvec =
            select_fun<T, getMassMatrixField, CPUGPU, matrixMatvecEigs>(
                  primme->massMatrixMatvec_type);
   }
   if (primme->correctionParams.precondition) {
      primme->applyPreconditioner =
            select_fun<T, getPreconditionerField, CPUGPU, matrixMatvecEigs>(
                  primme->applyPreconditioner_type);
   }
   if (primme->monitor) {
      primme->monitorFun = monitorFunEigs<T, EVAL>;
   }
   if (primme->convtest) {
      primme->convTestFun = convTestFunEigs<T, CPUGPU, EVAL>;
   }


#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__)) || defined (__FreeBSD__)
   // Set ctrl+c handler

   keepRunning = 1;
   SIGHANDLER_T prev_handler = signal(SIGINT, interrumptHandler);
   if (prev_handler == interrumptHandler) prev_handler = NULL;
#endif

   // Call xprimme

   int ret = tprimme(evals, evecs, rnorms, primme, CPUGPU());

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__)) || defined (__FreeBSD__)
   // Unset ctrl+c handler

   signal(SIGINT, prev_handler);
#endif

   // Return error code

   plhs[0] = create_mxArray(&ret, 1, 1, 1, CPU());

   // Return evals

   if (nlhs >= 2) {
      if (!mxEvals) {
         mxEvals = create_mxArray(evals, (int)primme->initSize, 1,
               (int)primme->initSize, CPUGPU());
         delete [] evals;
      }
      mxSetM(mxEvals, primme->initSize);
      plhs[1] = mxEvals;
   }
   else {
      delete [] evals;
   }

   // Return rnorms

   if (nlhs >= 3) {
      mxSetM(mxRnorms, primme->initSize);
      plhs[2] = mxRnorms;
   }
   else {
      delete [] rnorms;
   }

   // Return evecs

   if (nlhs >= 4) {
      if (primme->numOrthoConst > 0 || (isCPU(CPUGPU()) && !mxEvecs)) {
         mxEvecs = create_mxArray(&evecs[primme->n * primme->numOrthoConst],
               primme->n, (PRIMME_INT)primme->initSize, primme->n, CPUGPU());
         if (isCPU(CPUGPU()) && evecs) {
            delete[] evecs;
         }
      }
#ifdef USE_GPUARRAY
      else if (!mxEvecs && isGPU(CPUGPU())) {
         mxEvecs = mxGPUCreateMxArrayOnGPU(mxgpuEvecs);
      }
#endif
      plhs[3] = mxEvecs;
   }
   else if (isCPU(CPUGPU())) {
      if (evecs) delete [] evecs;
   }

#ifdef USE_GPUARRAY
   if (isGPU(CPUGPU())) {
      if (mxgpuEvecs) mxGPUDestroyGPUArray(mxgpuEvecs);
   }
#endif
}

// Wrapper around primme_svds_initialize; prototype:
// [primme_svds] = mexFunction_primme_svds_initialize()

static void mexFunction_primme_svds_initialize(int nlhs, mxArray *plhs[],
      int nrhs, const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(0);
   ASSERT_NUMARGSOUT(1);

   primme_svds_params *primme_svds = new primme_svds_params;
   primme_svds_initialize(primme_svds);
   plhs[0] = mxArrayFromPointer(primme_svds);
}

// Wrapper around primme_svds_set_method; prototype:
// mexFunction_primme_svds_set_method(primme_svds_preset_method method,
//    primme_preset_method methodStage1, primme_preset_method methodStage2,
//    primme_svds)

static void mexFunction_primme_svds_set_method(int nlhs, mxArray *plhs[],
      int nrhs, const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(4);
   ASSERT_NUMARGSOUT(0);
   ASSERT_NUMERIC_OR_CHAR(0);
   ASSERT_NUMERIC_OR_CHAR(1);
   ASSERT_NUMERIC_OR_CHAR(2);
   ASSERT_POINTER(3);

   primme_svds_preset_method method = (primme_svds_preset_method)mxArrayToConstantSvds(prhs[0]);
   primme_preset_method methodStage1 = (primme_preset_method)mxArrayToConstant(prhs[1]);
   primme_preset_method methodStage2 = (primme_preset_method)mxArrayToConstant(prhs[2]);
   primme_svds_params *primme_svds = (primme_svds_params*)mxArrayToPointer(prhs[3]);
   CHKERR(primme_svds_set_method(method, methodStage1, methodStage2,
            primme_svds));
}

// Wrapper around primme_svds_free; prototype:
// mexFunction_primme_svds_free(primme_svds)

static void mexFunction_primme_svds_free(int nlhs, mxArray *plhs[], int nrhs,
      const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(1);
   ASSERT_NUMARGSOUT(0);
   ASSERT_POINTER(0);

   primme_svds_params *primme_svds = (primme_svds_params*)mxArrayToPointer(prhs[0]);
   if (primme_svds->targetShifts) delete [] primme_svds->targetShifts;
   if (primme_svds->matrix) mxDestroyArray((mxArray*)primme_svds->matrix);
   if (primme_svds->preconditioner) mxDestroyArray((mxArray*)primme_svds->preconditioner);
   if (primme_svds->monitor) mxDestroyArray((mxArray*)primme_svds->monitor);
   if (primme_svds->convtest) mxDestroyArray((mxArray*)primme_svds->convtest);
   if (primme_svds->commInfo) mxDestroyArray((mxArray*)primme_svds->commInfo);
#ifdef USE_GPUARRAY
   if (primme_svds->queue) {
      // and finalize MAGMA.
      magma_queue_destroy(*(magma_queue_t*)primme_svds->queue);
      free(primme_svds->queue);
      magma_finalize();
   }
#endif

   primme_svds_free(primme_svds);
   delete primme_svds;
}

// Wrapper around primme_svds_set_member; prototype:
// mexFunction_primme_svds_set_member(primme_svds, label, value)

static void mexFunction_primme_svds_set_member(int nlhs, mxArray *plhs[],
      int nrhs, const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(3);
   ASSERT_NUMARGSOUT(0);
   ASSERT_POINTER(0);
   ASSERT_NUMERIC_OR_CHAR(1);

   primme_svds_params *primme_svds = (primme_svds_params*)mxArrayToPointer(prhs[0]);
   primme_svds_params_label label = mxArrayToLabelSvds(prhs[1]);
   switch(label) {
      // Set members with arity > 1
 
      case PRIMME_SVDS_iseed:
      {
         ASSERT_NUMERIC(2);
         copy_mxArray(prhs[2], primme_svds->iseed, 4, 1, 4, CPU());
         break;
      }


      case PRIMME_SVDS_targetShifts:
      {
         ASSERT_NUMERIC(2);
         if (primme_svds->targetShifts) delete [] primme_svds->targetShifts;
         int n = (int)mxGetNumberOfElements(prhs[2]);
         primme_svds->targetShifts = new double[n];
         primme_svds->numTargetShifts = n;
         copy_mxArray(prhs[2], primme_svds->targetShifts, n, 1, n, CPU());
         break;
      }

      // The function handlers are stored in the user data fields in
      // primme_params, e.g., in matrix for matrixMatvec and preconditioner
      // for applyPreconditioner

      case PRIMME_SVDS_matrixMatvec: 
      {
         ASSERT_FUNCTION(2);
         if (primme_svds->matrix) mxDestroyArray((mxArray*)primme_svds->matrix);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme_svds->matrix = (void*)a;
         break;
      }
      case PRIMME_SVDS_applyPreconditioner:
      {
         ASSERT_FUNCTION(2);
         if (primme_svds->preconditioner) mxDestroyArray((mxArray*)primme_svds->preconditioner);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme_svds->preconditioner = (void*)a;
         break;
      }
      case PRIMME_SVDS_convTestFun:
      {
         ASSERT_FUNCTION(2);
         if (primme_svds->convtest) mxDestroyArray((mxArray*)primme_svds->convtest);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme_svds->convtest = (void*)a;
         break;
      }
      case PRIMME_SVDS_monitorFun:
      {
         ASSERT_FUNCTION(2);
         if (primme_svds->monitor) mxDestroyArray((mxArray*)primme_svds->monitor);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme_svds->monitor = (void*)a;
         break;
      }
      case PRIMME_SVDS_commInfo:
      {
         ASSERT_NUMERIC_OR_CHAR(2);
         if (primme_svds->commInfo) mxDestroyArray((mxArray*)primme_svds->commInfo);
         mxArray *a = mxDuplicateArray(prhs[2]);
         mexMakeArrayPersistent(a);
         primme_svds->commInfo = (void*)a;
         break;
      }


      // Forbidden members
 
      case PRIMME_SVDS_primme: 
      case PRIMME_SVDS_primmeStage2:
      case PRIMME_SVDS_numProcs: 
      case PRIMME_SVDS_procID: 
      case PRIMME_SVDS_mLocal: 
      case PRIMME_SVDS_nLocal: 
      case PRIMME_SVDS_globalSumReal:
      case PRIMME_SVDS_numTargetShifts:
      case PRIMME_SVDS_matrix:
      case PRIMME_SVDS_preconditioner:
      case PRIMME_SVDS_outputFile:
      case PRIMME_SVDS_convtest:
      case PRIMME_SVDS_monitor:
         mexErrMsgTxt("Unsupported to set this option");
         break;

      default : 
      {
         primme_type ptype;
         int arity;

         // Get type of the option
         int ierr = primme_svds_member_info(&label, NULL, &ptype, &arity);
         if (ierr != 0) {
            mexErrMsgTxt("Unknown option");
         }
         CHKERR(arity != 1);

         // Set members with type int

         if (ptype == primme_int) {
            PRIMME_INT v = mxArrayToConstantSvds(prhs[2]);
            CHKERR(primme_svds_set_member(primme_svds, label, &v));
         }

         // Set members with type double

         else if (ptype == primme_double) {
            ASSERT_NUMERIC(2);
            double v;
            copy_mxArray(prhs[2], &v, 1, 1, 1, CPU());
            CHKERR(primme_svds_set_member(primme_svds, label, &v));
         }

         // Set members with type string

         else if (ptype == primme_string) {
            const char *str = mxArrayToString(prhs[2]);
            CHKERR(primme_svds_set_member(primme_svds, label, (void*)str));
         }

         else {
            /* This shouldn't happen */
            CHKERR(1);
         }
      }
   }
}

// Wrapper around primme_svds_get_member; prototype:
// mexFunction_primme_svds_get_member(primme_svds, label, value)

static void mexFunction_primme_svds_get_member(int nlhs, mxArray *plhs[],
      int nrhs, const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(2);
   ASSERT_NUMARGSOUT(1);
   ASSERT_POINTER(0);
   ASSERT_NUMERIC_OR_CHAR(1);

   primme_svds_params *primme_svds = (primme_svds_params*)mxArrayToPointer(prhs[0]);
   primme_svds_params_label label = mxArrayToLabelSvds(prhs[1]);
   switch(label) {
      // Set members with arity > 1

      case PRIMME_SVDS_iseed:
      {
         plhs[0] = create_mxArray(primme_svds->iseed, 4, 1, 4, CPU());
         break;
      }

      case PRIMME_SVDS_targetShifts:
      {
         plhs[0] = create_mxArray(primme_svds->targetShifts,
               primme_svds->numTargetShifts, 1, primme_svds->numTargetShifts, CPU());
         break;
      }

      // The function handlers are stored in the user data fields in
      // primme_params, e.g., in matrix for matrixMatvec and preconditioner
      // for applyPreconditioner

      case PRIMME_SVDS_matrixMatvec: 
      {
         plhs[0] = (mxArray*)primme_svds->matrix;
         break;
      }
      case PRIMME_SVDS_applyPreconditioner:
      {
         plhs[0] = (mxArray*)primme_svds->preconditioner;
         break;
      }
      case PRIMME_SVDS_convTestFun:
      {
         plhs[0] = (mxArray*)primme_svds->convtest;
         break;
      }
      case PRIMME_SVDS_monitorFun:
      {
         plhs[0] = (mxArray*)primme_svds->monitor;
         break;
      }

      // Get primme_params references

      case PRIMME_SVDS_primme:
      { 
         plhs[0] = mxArrayFromPointer(&primme_svds->primme);
         break;
      }
      case PRIMME_SVDS_primmeStage2:
      { 
         plhs[0] = mxArrayFromPointer(&primme_svds->primmeStage2);
         break;
      }

      // Forbidden members
 
      case PRIMME_SVDS_numProcs: 
      case PRIMME_SVDS_procID: 
      case PRIMME_SVDS_mLocal: 
      case PRIMME_SVDS_nLocal: 
      case PRIMME_SVDS_commInfo:
      case PRIMME_SVDS_globalSumReal:
      case PRIMME_SVDS_numTargetShifts:
      case PRIMME_SVDS_matrix:
      case PRIMME_SVDS_preconditioner:
      case PRIMME_SVDS_outputFile:
      case PRIMME_SVDS_convtest:
      case PRIMME_SVDS_monitor:
         mexErrMsgTxt("Unsupported to set this option");
         break;

      default : 
      {
         primme_type ptype;
         int arity;

         // Get type of the option
         int ierr = primme_svds_member_info(&label, NULL, &ptype, &arity);
         if (ierr != 0) {
            mexErrMsgTxt("Unknown option");
         }
         CHKERR(arity != 1);

         // Get members with type int

         if (ptype == primme_int) {
            PRIMME_INT v;
            CHKERR(primme_svds_get_member(primme_svds, label, &v));
            plhs[0] = create_mxArray(&v, 1, 1, 1, CPU());
         }

         // Set members with type double

         else if (ptype == primme_double) {
            double v;
            CHKERR(primme_svds_get_member(primme_svds, label, &v));
            plhs[0] = create_mxArray(&v, 1, 1, 1, CPU());
         }

         // Get members with type string

         else if (ptype == primme_string) {
            const char *v;
            CHKERR(primme_svds_get_member(primme_svds, label, &v));
            plhs[0] = mxCreateString(v);
         }

         else {
            /* This shouldn't happen */
            CHKERR(1);
         }
      }
   }
}

// Auxiliary functions for mexFunction_xprimme_svds; they return the next
// information based on the mode and primme_svds_params:
// Arguments:
// - mode: the corresponding input value in primme_svds.matrixMatvec and
//         primme_svds.applyPreconditioner.
// - primme_svds: primme_svds_params
// - mx: return the number of rows of input vectors x
// - my: return the number of rows of output vectors y
// - str: return the corresponding string for mode (notransp/transp or
//        AHA/AAH/aug).

struct getSvdsForMatrix {
   static void get(int transpose, primme_svds_params *primme_svds,
         PRIMME_INT *mx, PRIMME_INT *my, mxArray **AFUN, const char **str) {
      *AFUN = (mxArray*)primme_svds->matrix;
      if (transpose == 0) { /* Doing y <- A * x */
         *mx = primme_svds->n;
         *my = primme_svds->m;
         *str = "notransp";
      }
      else { /* Doing y <- A' * x */
         *mx = primme_svds->m;
         *my = primme_svds->n;
         *str = "transp";
      }
   }
};

struct getSvdsForPreconditioner {
   static void get(int mode, primme_svds_params *primme_svds,
         PRIMME_INT *mx, PRIMME_INT *my, mxArray **AFUN, const char **str) {
      *AFUN = (mxArray*)primme_svds->preconditioner;
      if (mode == primme_svds_op_AtA) {
         /* Preconditioner for A^t*A */
         *mx = *my = primme_svds->n;
         *str = "AHA";
      }
      else if (mode == primme_svds_op_AAt) {
         /* Preconditioner for A*A^t */
         *mx = *my = primme_svds->m;
         *str = "AAH";
      }
      else if (mode == primme_svds_op_augmented) {
         /* Preconditioner for [0 A^t; A 0] */
         *mx = *my = primme_svds->m + primme_svds->n;
         *str = "aug";
      }
      else {
         mexErrMsgTxt("Unsupported preconditioner type");
         *mx = *my = -1;
         *str = "";
      }
   }
};


// Auxiliary function for mexFunction_xprimme_svds; PRIMME wrapper around
// matrixMatvec and applyPreconditioner. Create a mxArray
// from input vector x, call the function handler returned by F and
// copy the content of its returned mxArray into the output vector y. The
// functor F returns also the number of rows in x and y and the string
// passed in callback depending on mode.

template <typename T, typename F, typename CPUGPU>
struct matrixMatvecSvds {
   typedef void (*t)(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, int *mode, struct primme_svds_params *primme_svds,
      int *ierr);
   static void f(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, int *mode, struct primme_svds_params *primme_svds,
      int *ierr)
   {  
      mxArray *prhs[3], *plhs[1];

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
      // Check interrupt handler
      if (!keepRunning) {
         *ierr = 1;
         return;
      }
#endif

      if (*blockSize <= 0) {*ierr = 0; return;}

      // Get numbers of rows of x and y
      PRIMME_INT mx, my;
      const char *str;
      F::get(*mode, primme_svds, &mx, &my, &prhs[0], &str);
      assert(mx > 0);

      // Create input vector x (avoid copy if possible)

      prhs[1] = create_mxArray(
            (T *)x, mx, (PRIMME_INT)*blockSize, *ldx, CPUGPU(), true);
      prhs[2] = mxCreateString(str);

      // Call the callback

      *ierr = mexCallMATLAB(1, plhs, 3, prhs, "feval");

      // Copy lhs[0] to y and destroy it

      if (plhs[0]) {
         copy_mxArray(
               plhs[0], (T *)y, my, (PRIMME_INT)*blockSize, *ldy, CPUGPU());
         mxDestroyArray(plhs[0]);
      }

      // Destroy prhs[*]

      if (isCPU(CPUGPU()) && mxGetData(prhs[1]) == x) mxSetData(prhs[1], NULL);
      mxDestroyArray(prhs[1]); 
      mxDestroyArray(prhs[2]); 
   }
};

template <typename T, typename CPUGPU>
static void convTestFunSvds(double *sval, void *leftsvec, void *rightsvec,
      double *rNorm, int *method, int *isConv,
      struct primme_svds_params *primme_svds, int *ierr) {
   (void)method;
 
   mxArray *prhs[5], *plhs[1];

   // Create input vectors (avoid copy if possible)

   prhs[1] = create_mxArray<double,int>(sval, sval?1:0, 1, sval?1:0, CPU());
   prhs[2] = create_mxArray((T *)leftsvec,
         leftsvec ? primme_svds->mLocal : (PRIMME_INT)0, (PRIMME_INT)1,
         leftsvec ? primme_svds->mLocal : (PRIMME_INT)0, CPUGPU(), true);
   prhs[3] = create_mxArray((T *)rightsvec,
         rightsvec ? primme_svds->nLocal : (PRIMME_INT)0, (PRIMME_INT)1,
         rightsvec ? primme_svds->nLocal : (PRIMME_INT)0, CPUGPU(), true);
   prhs[4] = create_mxArray<double,int>(rNorm, rNorm?1:0, 1, rNorm?1:0, CPU());

   // Call the callback

   prhs[0] = (mxArray*)primme_svds->convtest;
   *ierr = mexCallMATLAB(1, plhs, 5, prhs, "feval");

   // Copy plhs[0] to isConv and destroy it

   if (plhs[0]) {
      copy_mxArray(plhs[0], isConv, 1, 1, 1, CPU());
      mxDestroyArray(plhs[0]);
   }

   // Destroy prhs[1..4]

   if (isCPU(CPUGPU()) && mxGetData(prhs[2]) == leftsvec)
      mxSetData(prhs[2], NULL);
   if (isCPU(CPUGPU()) && mxGetData(prhs[3]) == rightsvec)
      mxSetData(prhs[3], NULL);
   for (int i=1; i<5; i++) mxDestroyArray(prhs[i]); 
}

template <typename T>
static void monitorFunSvds(void *basisSvals, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms, int *numConverged,
      void *lockedSvals, int *numLocked, int *lockedFlags, void *lockedNorms,
      int *inner_its, void *LSRes, const char *msg, double *time,
      primme_event *event, int *stage, struct primme_svds_params *primme_svds,
      int *ierr) {

   mxArray *prhs[15];

   // Create input vectors (avoid copy if possible)

   typedef typename Real<T>::type R;
   prhs[1] = create_mxArray<R, int>((R *)basisSvals, basisSize ? *basisSize : 0,
         1, basisSize ? *basisSize : 0, CPU(), true);
   prhs[2] = create_mxArray<int, int>(basisFlags, basisFlags ? *basisSize : 0,
         1, basisFlags ? *basisSize : 0, CPU(), true);
   prhs[3] = create_mxArray<int, int>(iblock, blockSize ? *blockSize : 0, 1,
         blockSize ? *blockSize : 0, CPU(), true);
   prhs[4] = create_mxArray<R, int>((R *)basisNorms, basisSize ? *basisSize : 0,
         1, basisSize ? *basisSize : 0, CPU(), true);
   prhs[5] = create_mxArray<int, int>(
         numConverged, numConverged ? 1 : 0, 1, 1, CPU(), true);
   int numLocked0 = numLocked && *numLocked > 0 ? *numLocked : 0;
   prhs[6] = create_mxArray<R, int>(
         (R *)lockedSvals, numLocked0, 1, numLocked0, CPU(), true);
   prhs[7] = create_mxArray<int, int>(
         lockedFlags, numLocked0, 1, numLocked0, CPU(), true);
   prhs[8] = create_mxArray<R, int>(
         (R *)lockedNorms, numLocked0, 1, numLocked0, CPU(), true);
   prhs[9] = create_mxArray<int, int>(inner_its, inner_its ? 1 : 0, 1, 1, CPU(), true);
   prhs[10] = create_mxArray<R, int>((R *)LSRes, LSRes ? 1 : 0, 1, 1, CPU(), true);
   prhs[11] = create_mxArray(msg, true);
   prhs[12] = create_mxArray<double, int>(time, time ? 1 : 0, 1, 1, CPU(), true);
   prhs[13] = create_mxArray<int, int>((int *)event, event ? 1 : 0, 1, 1, CPU(), true);
   prhs[14] = create_mxArray<int, int>(stage, stage ? 1 : 0, 1, 1, CPU(), true);

   // Call the callback

   prhs[0] = (mxArray*)primme_svds->monitor;
   *ierr = mexCallMATLAB(0, NULL, 15, prhs, "feval");

   // Destroy prhs[1..12]

   if (mxGetData(prhs[1]) == basisSvals)   mxSetData(prhs[1], NULL);
   if (mxGetData(prhs[2]) == basisFlags)   mxSetData(prhs[2], NULL);
   if (mxGetData(prhs[3]) == iblock)       mxSetData(prhs[3], NULL);
   if (mxGetData(prhs[4]) == basisNorms)   mxSetData(prhs[4], NULL);
   if (mxGetData(prhs[5]) == numConverged) mxSetData(prhs[5], NULL);
   if (mxGetData(prhs[6]) == lockedSvals)  mxSetData(prhs[6], NULL);
   if (mxGetData(prhs[7]) == lockedFlags)  mxSetData(prhs[7], NULL);
   if (mxGetData(prhs[8]) == lockedNorms)  mxSetData(prhs[8], NULL);
   if (mxGetData(prhs[9]) == inner_its)    mxSetData(prhs[9], NULL);
   if (mxGetData(prhs[10]) == LSRes)       mxSetData(prhs[10], NULL);
   if (mxGetData(prhs[11]) == msg)         mxSetData(prhs[11], NULL);
   if (mxGetData(prhs[12]) == time)        mxSetData(prhs[12], NULL);
   if (mxGetData(prhs[13]) == event)       mxSetData(prhs[13], NULL);
   if (mxGetData(prhs[14]) == stage)       mxSetData(prhs[14], NULL);
   for (int i=1; i<15; i++) mxDestroyArray(prhs[i]); 
}


// Wrapper around xprimme_svds; prototype:
// [ret, evals, rnorms, evecs] = mexFunction_xprimme_svds(...
//                            init_guesses_left, init_guesses_right, primme_svds)

template<typename T, typename CPUGPU>
static void mexFunction_xprimme_svds(int nlhs, mxArray *plhs[], int nrhs,
      const mxArray *prhs[])
{
   ASSERT_NUMARGSIN(3);
   ASSERT_NUMARGSOUTGE(1);
   ASSERT_POINTER(2);

   primme_svds_params *primme_svds = (primme_svds_params*)mxArrayToPointer(prhs[2]);

#ifdef USE_GPUARRAY
   if (isGPU(CPUGPU())) {
      if (!primme_svds->commInfo)
         mexErrMsgTxt("Set primme_svds.commInfo with the GPU index");

      // Initialize the MathWorks GPU API
      mxInitGPU();

      // Initialize MAGMA
      magma_init();

      // Create a context
      int gpuDevice;
      copy_mxArray((mxArray*)primme_svds->commInfo, &gpuDevice, 1, 1, 1, CPU());
      if (primme_svds->queue)
         magma_queue_destroy(*(magma_queue_t *)primme_svds->queue);
      else
         primme_svds->queue = malloc(sizeof(magma_queue_t));
      magma_queue_create(gpuDevice, (magma_queue_t*)primme_svds->queue);
   }
#endif


   // Allocate svals, rnorms and svecs; if possible create the mxArray and use
   // its data

   typename Real<T>::type *svals, *rnorms;
   T *svecs;
   mxArray *mxSvals, *mxRnorms;

   if (nlhs <= 0) {
      svals = new typename Real<T>::type[primme_svds->numSvals];
      mxSvals = NULL;
   }
   else {
      mxSvals = mxCreateNumericMatrix(primme_svds->numSvals, 1,
            toClassID<typename Real<T>::type>(), mxREAL);
      svals = (typename Real<T>::type*)mxGetData(mxSvals);
   }

   if (nlhs <= 1) {
      rnorms = new typename Real<T>::type[primme_svds->numSvals];
      mxRnorms = NULL;
   }
   else {
      mxRnorms = mxCreateNumericMatrix(primme_svds->numSvals, 1,
            toClassID<typename Real<T>::type>(), mxREAL);
      rnorms = (typename Real<T>::type*)mxGetData(mxRnorms);
   }

   PRIMME_INT n = primme_svds->numOrthoConst
      + macro_max(primme_svds->initSize, primme_svds->numSvals);
#ifdef USE_GPUARRAY
   mxGPUArray *mxgpuSvecs  = nullptr; 
   if (isGPU(CPUGPU())) {
      mwSize dims[2] = {(mwSize)primme_svds->m+primme_svds->n, (mwSize)n};
      mxgpuSvecs = mxGPUCreateGPUArray(2, dims, toClassID<T>(),
            isComplex<T>() ? mxCOMPLEX : mxREAL, MX_GPU_DO_NOT_INITIALIZE);
      svecs = (T*)mxGPUGetData(mxgpuSvecs);
   } else 
#endif
   {
      svecs = new T[n*(primme_svds->m+primme_svds->n)];
   }

   // Copy initial vectors

   if (primme_svds->numOrthoConst + primme_svds->initSize > 0) {
      ASSERT_NUMERIC(0);
      int ninit = primme_svds->numOrthoConst+primme_svds->initSize;
      copy_mxArray(prhs[0], svecs, primme_svds->m, (PRIMME_INT)ninit,
            primme_svds->m, CPUGPU());
      copy_mxArray(prhs[1], &svecs[primme_svds->m*ninit], primme_svds->n,
            (PRIMME_INT)ninit, primme_svds->n, CPUGPU());
   }

   // Set matvec and preconditioner and monitorFun and convTestFun

   primme_svds->matrixMatvec =
         select_fun<T, getSvdsForMatrix, CPUGPU, matrixMatvecSvds>(
               primme_svds->matrixMatvec_type);
   if (primme_svds->preconditioner) {
      primme_svds->applyPreconditioner =
            select_fun<T, getSvdsForPreconditioner, CPUGPU, matrixMatvecSvds>(
                  primme_svds->applyPreconditioner_type);
   }
   if (primme_svds->monitor) {
      primme_svds->monitorFun = monitorFunSvds<T>;
   }
   if (primme_svds->convtest) {
      primme_svds->convTestFun = convTestFunSvds<T, CPUGPU>;
   }

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__)) || defined (__FreeBSD__)
   // Set ctrl+c handler

   keepRunning = 1;
   SIGHANDLER_T prev_handler = signal(SIGINT, interrumptHandler);
   if (prev_handler == interrumptHandler) prev_handler = NULL;
#endif

   // Call xprimme_svds

   int ret = tprimme_svds(svals, svecs, rnorms, primme_svds, CPUGPU());

#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__)) || defined (__FreeBSD__)
   // Unset ctrl+c handler

   signal(SIGINT, prev_handler);
#endif

   // Return error code

   plhs[0] = create_mxArray(&ret, 1, 1, 1, CPU());

   // Return svals

   if (nlhs >= 2) {
      mxSetM(mxSvals, primme_svds->initSize);
      plhs[1] = mxSvals;
   }
   else {
      delete [] svals;
   }

   // Return rnorms

   if (nlhs >= 3) {
      mxSetM(mxRnorms, primme_svds->initSize);
      plhs[2] = mxRnorms;
   }
   else {
      delete [] rnorms;
   }

   // Return svecs

   if (nlhs >= 4) {
      plhs[3] = create_mxArray(
            &svecs[primme_svds->m * primme_svds->numOrthoConst], primme_svds->m,
            (PRIMME_INT)primme_svds->initSize, primme_svds->m, CPUGPU());
   }

   if (nlhs >= 5) {
      plhs[4] = create_mxArray(
            &svecs[primme_svds->m *
                         (primme_svds->numOrthoConst + primme_svds->initSize) +
                   primme_svds->n * primme_svds->numOrthoConst],
            primme_svds->n, (PRIMME_INT)primme_svds->initSize, primme_svds->n,
            CPUGPU());
   }

#ifdef USE_GPUARRAY
   if (isGPU(CPUGPU())) {
      if (mxgpuSvecs) mxGPUDestroyGPUArray(mxgpuSvecs);
   } else 
#endif
      if (svecs) delete [] svecs;
}

// Main function: dispatch the call to the proper mexFunction_* function

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 0 || !mxIsChar(prhs[0])) {
      mexErrMsgTxt("The first argument should be char");
   }

   char *function_name = mxArrayToString(prhs[0]);
   bool called = false;

   #define PRIMME_TRY_CALL(F) \
      if (strcmp(#F, function_name) == 0) called=true,mexFunction_ ## F (nlhs, plhs, nrhs-1, &prhs[1]);
   #define PRIMME_TRY_CALL_T(F, FT, ...) \
      if (strcmp(#F, function_name) == 0) called=true,mexFunction_ ## FT < __VA_ARGS__ > (nlhs, plhs, nrhs-1, &prhs[1]);

   PRIMME_TRY_CALL_T(sprimme, xprimme, float, CPU, float);
   PRIMME_TRY_CALL_T(cprimme, xprimme, std::complex<float>, CPU, float);
   PRIMME_TRY_CALL_T(dprimme, xprimme, double, CPU, double);
   PRIMME_TRY_CALL_T(zprimme, xprimme, std::complex<double>, CPU, double);
   PRIMME_TRY_CALL_T(cprimme_normal, xprimme, std::complex<float>, CPU,  std::complex<float>);
   PRIMME_TRY_CALL_T(zprimme_normal, xprimme, std::complex<double>, CPU, std::complex<double>);
#ifdef USE_GPUARRAY
   PRIMME_TRY_CALL_T(magma_sprimme, xprimme, float, GPU, float);
   PRIMME_TRY_CALL_T(magma_cprimme, xprimme, std::complex<float>, GPU, float);
   PRIMME_TRY_CALL_T(magma_dprimme, xprimme, double, GPU, double);
   PRIMME_TRY_CALL_T(magma_zprimme, xprimme, std::complex<double>, GPU, double);
   PRIMME_TRY_CALL_T(magma_cprimme_normal, xprimme, std::complex<float>, GPU,  std::complex<float>);
   PRIMME_TRY_CALL_T(magma_zprimme_normal, xprimme, std::complex<double>, GPU, std::complex<double>);
#endif
   PRIMME_TRY_CALL(primme_initialize);
   PRIMME_TRY_CALL(primme_set_method);
   PRIMME_TRY_CALL(primme_get_member);
   PRIMME_TRY_CALL(primme_set_member);
   PRIMME_TRY_CALL(primme_free);
   PRIMME_TRY_CALL_T(sprimme_svds, xprimme_svds, float, CPU);
   PRIMME_TRY_CALL_T(cprimme_svds, xprimme_svds, std::complex<float>, CPU);
   PRIMME_TRY_CALL_T(dprimme_svds, xprimme_svds, double, CPU);
   PRIMME_TRY_CALL_T(zprimme_svds, xprimme_svds, std::complex<double>, CPU);
#ifdef USE_GPUARRAY
   PRIMME_TRY_CALL_T(magma_sprimme_svds, xprimme_svds, float, GPU);
   PRIMME_TRY_CALL_T(magma_cprimme_svds, xprimme_svds, std::complex<float>, GPU);
   PRIMME_TRY_CALL_T(magma_dprimme_svds, xprimme_svds, double, GPU);
   PRIMME_TRY_CALL_T(magma_zprimme_svds, xprimme_svds, std::complex<double>, GPU);
#endif
   PRIMME_TRY_CALL(primme_svds_initialize);
   PRIMME_TRY_CALL(primme_svds_set_method);
   PRIMME_TRY_CALL(primme_svds_set_member);
   PRIMME_TRY_CALL(primme_svds_get_member);
   PRIMME_TRY_CALL(primme_svds_free);

   #undef PRIMME_TRY_CALL
   #undef PRIMME_TRY_CALL_T

   if (!called)
      mexErrMsgTxtPrintf1("Unavailable function: %s", function_name);
}
