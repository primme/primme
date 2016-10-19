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
 * File: primme.i
 * 
 * Purpose - SWIG description of PYTHON interface to PRIMME.
 * 
 ******************************************************************************/

%define DOCSTRING
"Find a few eigenvalues and eigenvectors, and also singular triplets of a matrix.
Uses PRIMME: https://github.com/primme/primme
"
%enddef

%module(docstring=DOCSTRING,directors="1") Primme

%pythoncode %{
__all__ = ['PrimmeParams', 'sprimme', 'cprimme', 'dprimme', 'zprimme', 'eigsh', 'PrimmeError', 'PRIMME_Arnoldi', 'PRIMME_DEFAULT_METHOD', 'PRIMME_DEFAULT_MIN_MATVECS', 'PRIMME_DEFAULT_MIN_TIME', 'PRIMME_DYNAMIC', 'PRIMME_GD', 'PRIMME_GD_Olsen_plusK', 'PRIMME_GD_plusK', 'PRIMME_JDQMR', 'PRIMME_JDQMR_ETol', 'PRIMME_JDQR', 'PRIMME_JD_Olsen_plusK', 'PRIMME_LOBPCG_OrthoBasis', 'PRIMME_LOBPCG_OrthoBasis_Window', 'PRIMME_RQI', 'PRIMME_SUBSPACE_ITERATION', 'primme_adaptive', 'primme_adaptive_ETolerance', 'primme_closest_abs', 'primme_closest_geq', 'primme_closest_leq', 'primme_decreasing_LTolerance', 'primme_dtr', 'primme_full_LTolerance', 'primme_init_default', 'primme_init_krylov', 'primme_init_random', 'primme_init_user', 'primme_largest', 'primme_largest_abs', 'primme_proj_RR', 'primme_proj_default', 'primme_proj_harmonic', 'primme_proj_refined', 'primme_smallest', 'primme_thick', 'PrimmeSvdsParams', 'svds', 'primme_svds_augmented', 'primme_svds_closest_abs', 'primme_svds_default', 'primme_svds_hybrid', 'primme_svds_largest', 'primme_svds_normalequations', 'primme_svds_op_AAt', 'primme_svds_op_AtA', 'primme_svds_op_augmented', 'primme_svds_op_none', 'primme_svds_smallest', 'sprimme_svds', 'cprimme_svds', 'dprimme_svds', 'zprimme_svds', 'PrimmeSvdsError']
%}
// Support PRIMME_INT for int64_t
%include "stdint.i"
%{
#define SWIG_FILE_WITH_INIT
#include "primmew.h"
%}

// Get the NumPy typemaps
%include "numpy.i"

// Handle standard exceptions
%include "exception.i"
%exception
{
  try
  {
    $action
  }
  catch (const std::invalid_argument& e)
  {
    SWIG_exception(SWIG_ValueError, e.what());
  }
  catch (const std::out_of_range& e)
  {
    SWIG_exception(SWIG_IndexError, e.what());
  }
  catch (Swig::DirectorException &e)
  {
     SWIG_fail;
  }
  if (PyErr_Occurred()) SWIG_fail;
}
%init %{
  import_array();
%}

// Global ignores
%ignore primme_display_params;
%ignore primme_set_method;
%ignore primme_initialize;
%ignore primme_free;

%ignore primme_svds_display_params;
%ignore primme_svds_set_method;
%ignore primme_svds_initialize;
%ignore primme_svds_free;

%ignore PRId64;

%ignore tprimme;
%ignore tprimme_svds;

%ignore PrimmeParams::matrixMatvec;
%ignore PrimmeParams::massMatrixMatvec;
%ignore PrimmeParams::applyPreconditioner;
%ignore PrimmeParams::convTestFun;
%ignore PrimmeParams::targetShifts;
%ignore PrimmeParams::numTargetShifts;
%ignore PrimmeParams::globalSumDouble;
%ignore PrimmeParams::commInfo;
%ignore PrimmeParams::intWork;
%ignore PrimmeParams::realWork;
%ignore PrimmeParams::outputFile;
%ignore PrimmeParams::matrix;
%ignore PrimmeParams::preconditioner;
%ignore PrimmeParams::ShiftsForPreconditioner;
%ignore primme_params::matrixMatvec;
%ignore primme_params::massMatrixMatvec;
%ignore primme_params::applyPreconditioner;
%ignore primme_params::convTestFun;
%ignore primme_params::targetShifts;
%ignore primme_params::numTargetShifts;
%ignore primme_params::globalSumDouble;
%ignore primme_params::commInfo;
%ignore primme_params::intWork;
%ignore primme_params::realWork;
%ignore primme_params::outputFile;
%ignore primme_params::matrix;
%ignore primme_params::preconditioner;
%ignore primme_params::ShiftsForPreconditioner;
%ignore PrimmeSvdsParams::matrixMatvec;
%ignore PrimmeSvdsParams::applyPreconditioner;
%ignore PrimmeSvdsParams::convTestFun;
%ignore PrimmeSvdsParams::targetShifts;
%ignore PrimmeSvdsParams::numTargetShifts;
%ignore PrimmeSvdsParams::commInfo;
%ignore PrimmeSvdsParams::globalSumDouble;
%ignore PrimmeSvdsParams::intWork;
%ignore PrimmeSvdsParams::realWork;
%ignore PrimmeSvdsParams::outputFile;
%ignore PrimmeSvdsParams::matrix;
%ignore PrimmeSvdsParams::preconditioner;
%ignore PrimmeSvdsParams::ShiftsForPreconditioner;
%ignore PrimmeSvdsParams::primme;
%ignore PrimmeSvdsParams::primmeStage2;
%ignore primme_svds_params::matrixMatvec;
%ignore primme_svds_params::applyPreconditioner;
%ignore primme_svds_params::convTestFun;
%ignore primme_svds_params::targetShifts;
%ignore primme_svds_params::numTargetShifts;
%ignore primme_svds_params::globalSumDouble;
%ignore primme_svds_params::commInfo;
%ignore primme_svds_params::intWork;
%ignore primme_svds_params::realWork;
%ignore primme_svds_params::outputFile;
%ignore primme_svds_params::matrix;
%ignore primme_svds_params::preconditioner;
%ignore primme_svds_params::ShiftsForPreconditioner;
%ignore primme_svds_params::primme;
%ignore primme_svds_params::primmeStage2;


%fragment("NumPy_Array_Requirements_extra",
          "header",
          fragment="NumPy_Array_Requirements")
{
  /* Require the given PyArrayObject to to be Fortran ordered.  If the
   * the PyArrayObject is already Fortran ordered, do nothing.  Else,
   * set the Fortran ordering flag and recompute the strides.
   * NOTE: based on require_fortran in numpy.i
   */
  int require_fortran2(PyArrayObject* ary, int ld)
  {
    int success = 1;
    if (array_numdims(ary) != 2) return 0;
    int single_dim = (array_size(ary, 0) == 1 || array_size(ary, 1) == 1);
    npy_intp * strides = array_strides(ary);
    if (!array_is_fortran(ary)) {
      strides[0] = strides[1];
      strides[1] = single_dim ? strides[0] : strides[0]*ld;
    } else {
      strides[1] = single_dim ? strides[0] : strides[0]*ld;
    }
    /* Set the Fortran ordered flag */
    /* Note that this should be done after strides change */
    PyArray_UpdateFlags(ary, NPY_ARRAY_FARRAY);
    return success;
  }
}
 
%define %numpy_typemaps_ext(DATA_TYPE, DATA_TYPECODE, DIM_TYPE)

/* Typemap suite for (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* IN_FARRAY2D)
   See description of ARGOUTVIEW_FARRAY2 in numpy.i
 */
%typemap(directorin,
         fragment="NumPy_Backward_Compatibility,NumPy_Array_Requirements_extra,NumPy_Fragments")
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DIM_TYPE LD, DATA_TYPE* IN_FARRAY2D)
{
  npy_intp dims[2] = { $1, $2 };
  PyObject* obj = PyArray_SimpleNewFromData(2, dims, DATA_TYPECODE, (void*)($4));
  PyArrayObject* array = (PyArrayObject*) obj;

  if (!array || !require_fortran2(array, $3))
        throw Swig::DirectorMethodException();
  $input = obj;
}

/* Typemap suite for (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* INPLACE_ARRAY2)
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DIM_TYPE LD, DATA_TYPE* OUT_FARRAY2D)
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input),
                                                 DATA_TYPECODE);
}
%typemap(in,numinputs=0)
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DIM_TYPE LD, DATA_TYPE* OUT_FARRAY2D)
{ $1 = $2 = $3 = 0; }

%typemap(directorargout,
         fragment="NumPy_Fragments")
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DIM_TYPE LD, DATA_TYPE* OUT_FARRAY2D)
  (PyArrayObject* array=NULL, PyObject* o=NULL)
{
  o = $result;
  array = obj_to_array_no_conversion(o, DATA_TYPECODE);
  if (!array || !require_dimensions(array,2) || !require_native(array) ||
        !require_c_or_f_contiguous(array))
     Swig::DirectorMethodException::raise("No valid type for object returned by $symname");
  if (($1) != (DIM_TYPE) array_size(array,0) ||
      ($2) != (DIM_TYPE) array_size(array,1))
          {Swig::DirectorMethodException::raise("No valid dimensions for object returned by $symname");}
  npy_intp * strides = array_strides(array);
  if (array_is_fortran(array)) {
    copy_matrix((DATA_TYPE*)array_data(array), ($1), ($2), (DIM_TYPE)(strides[1]/strides[0]), ($4), ($3));
  } else {
      DATA_TYPE *x = (DATA_TYPE*)array_data(array);
      npy_intp ldx = strides[0]/strides[1];
      for (int i=0; i<($1); i++)
         for (int j=0; j<($2); j++)
            ($4)[i+j*($3)] = x[i*ldx+j];
  }
}

/* Typemap suite for (DIM_TYPE DIM, DATA_TYPE* IN_ARRAY1D)
   See description of ARGOUTVIEW_FARRAY2 in numpy.i
 */
%typemap(directorin,
         fragment="NumPy_Backward_Compatibility,NumPy_Array_Requirements_extra,NumPy_Fragments")
  (DIM_TYPE DIM1, DATA_TYPE* IN_ARRAY1D)
{
  npy_intp dims[1] = { $1 };
  PyObject* obj = PyArray_SimpleNewFromData(1, dims, DATA_TYPECODE, (void*)($2));
  PyArrayObject* array = (PyArrayObject*) obj;

  if (!array || !require_c_or_f_contiguous(array))
        throw Swig::DirectorMethodException();
  $input = obj;
}

/* Typemap suite for (DIM_TYPE DIM, DATA_TYPE* OUT_ARRAY1D)
   See description of IN_ARRAY1 in numpy.i
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (DIM_TYPE DIM1, DATA_TYPE* OUT_ARRAY1D)
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input),
                                                 DATA_TYPECODE);
}
%typemap(in,numinputs=0)
  (DIM_TYPE DIM1, DATA_TYPE* OUT_ARRAY1D)
{ $1 = 0; }

%typemap(directorargout,
         fragment="NumPy_Fragments")
  (DIM_TYPE DIM1, DATA_TYPE* OUT_ARRAY1D)
  (PyArrayObject* array=NULL, PyObject* o=NULL)
{
  o = $result;
  array = obj_to_array_no_conversion(o, DATA_TYPECODE);
  if (!array || !require_dimensions(array,1) || !require_native(array) ||
        !require_c_or_f_contiguous(array))
     Swig::DirectorMethodException::raise("No valid type for object returned by $symname");
  if (($1) != (DIM_TYPE) array_size(array,0))
          {Swig::DirectorMethodException::raise("No valid dimensions for object returned by $symname");}
  copy_matrix((DATA_TYPE*)array_data(array), 1, ($1), 1, ($2), 1);
}
%enddef    /* %numpy_typemaps_ext() macro */

%numpy_typemaps_ext(double              , NPY_DOUBLE , int)
%numpy_typemaps_ext(std::complex<double>, NPY_CDOUBLE, int)
%numpy_typemaps_ext(float               , NPY_FLOAT  , int)
%numpy_typemaps_ext(std::complex<float> , NPY_CFLOAT , int)



%apply (int DIM1, double* INPLACE_ARRAY1) {
   (int lenEvals, double* evals),
   (int lenSvals, double* svals),
   (int lenResNorms, double* resNorms)};
%apply (int DIM1, float* INPLACE_ARRAY1) {
   (int lenEvals, float* evals),
   (int lenSvals, float* svals),
   (int lenResNorms, float* resNorms)};
%apply (int DIM1, int DIM2, double* INPLACE_FARRAY2) {
   (int len1Evecs, int len2Evecs, double* evecs),
   (int len1SvecsLeft, int len2SvecsLeft, double* svecsLeft),
   (int len1SvecsRight, int len2SvecsRight, double* svecsRight)};
%apply (int DIM1, int DIM2, float* INPLACE_FARRAY2) {
   (int len1Evecs, int len2Evecs, float* evecs),
   (int len1SvecsLeft, int len2SvecsLeft, float* svecsLeft),
   (int len1SvecsRight, int len2SvecsRight, float* svecsRight)};
%apply (int DIM1, int DIM2, std::complex<double>* INPLACE_FARRAY2) {
   (int len1Evecs, int len2Evecs, std::complex<double>* evecs),
   (int len1SvecsLeft, int len2SvecsLeft, std::complex<double>* svecsLeft),
   (int len1SvecsRight, int len2SvecsRight, std::complex<double>* svecsRight)};
%apply (int DIM1, int DIM2, std::complex<float>* INPLACE_FARRAY2) {
   (int len1Evecs, int len2Evecs, std::complex<float>* evecs),
   (int len1SvecsLeft, int len2SvecsLeft, std::complex<float>* svecsLeft),
   (int len1SvecsRight, int len2SvecsRight, std::complex<float>* svecsRight)};

%apply (int DIM1, int DIM2, int LD, double* IN_FARRAY2D) {
   (int len1YD, int len2YD, int ldYD, double* yd)};
%apply (int DIM1, int DIM2, int LD, float* IN_FARRAY2D) {
   (int len1YD, int len2YD, int ldYD, float* yd)};
%apply (int DIM1, int DIM2, int LD, double* OUT_FARRAY2D) {
   (int len1XD, int len2XD, int ldXD, double* xd)};
%apply (int DIM1, int DIM2, int LD, float* OUT_FARRAY2D) {
   (int len1XD, int len2XD, int ldXD, float* xd)};
%apply (int DIM1, int DIM2, int LD, std::complex<double>* IN_FARRAY2D) {
   (int len1YD, int len2YD, int ldYD, std::complex<double>* yd)};
%apply (int DIM1, int DIM2, int LD, std::complex<float>* IN_FARRAY2D) {
   (int len1YD, int len2YD, int ldYD, std::complex<float>* yd)};
%apply (int DIM1, int DIM2, int LD, std::complex<double>* OUT_FARRAY2D) {
   (int len1XD, int len2XD, int ldXD, std::complex<double>* xd)};
%apply (int DIM1, int DIM2, int LD, std::complex<float>* OUT_FARRAY2D) {
   (int len1XD, int len2XD, int ldXD, std::complex<float>* xd)};

/* typemaps for targetShift and numTargetShifts */
 
%apply (double* IN_ARRAY1, int DIM1) {
   (double *targetShifts, int n)};
%apply (double **ARGOUTVIEW_ARRAY1, int *DIM1) {
   (double **targetShifts, int *n)};

/* typemaps for globalSumReal */
%apply (int DIM1, double* IN_ARRAY1D) {
   (int lenYD, double *yd)};
%apply (int DIM1, float* IN_ARRAY1D) {
   (int lenYD, float *yd)};
%apply (int DIM1, double *OUT_ARRAY1D) {
   (int lenXD, double *xd)};
%apply (int DIM1, float *OUT_ARRAY1D) {
   (int lenXD, float *xd)};

%inline %{
template <typename T, typename I, typename J>
static void copy_matrix(T *x, I m, J n, I ldx, T *y, I ldy) {
   I i; J j;

   assert(ldx >= m && ldy >= m);

   /* Do nothing if x and y are the same matrix */
   if (x == y && ldx == ldy) return;

   /* Copy a contiguous memory region */
   if (ldx == ldy && ldx == m) {
      memmove(y, x, sizeof(T)*(size_t)m*(size_t)n);
   }

   /* Copy matrix some rows down or up */
   else if (ldx == ldy && (y > x ? y-x : x-y) < ldx) {
      for (i=0; i<n; i++)
         memmove(&y[i*ldy], &x[i*ldx], sizeof(T)*(size_t)m);
   }

   /* Copy matrix some columns forward */
   else if (ldx == ldy && y > x && y-x > ldx) {
      for (i=n-1; i>=0; i--)
         for (j=0; j<m; j++)
            y[i*ldy+j] = x[i*ldx+j];
   }

   /* Copy matrix some columns backward and the general case */
   else {
      /* TODO: assert x and y don't overlap */
      for (i=0; i<n; i++)
         for (j=0; j<m; j++)
            y[i*ldy+j] = x[i*ldx+j];
   }

}


static int tprimme(float *evals, float *evecs, float *resNorms, primme_params *primme) {
      return sprimme(evals, evecs, resNorms, primme);
}
static int tprimme(float *evals, std::complex<float> *evecs, float *resNorms, primme_params *primme) {
      return cprimme(evals, evecs, resNorms, primme);
}
static int tprimme(double *evals, double *evecs, double *resNorms, primme_params *primme) {
      return dprimme(evals, evecs, resNorms, primme);
}
static int tprimme(double *evals, std::complex<double> *evecs, double *resNorms, primme_params *primme) {
      return zprimme(evals, evecs, resNorms, primme);
}

template <typename T>
static void mymatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, struct primme_params *primme, int *ierr) {
    PrimmeParams *pp = static_cast<PrimmeParams*>(primme);
    pp->matvec((int)primme->nLocal, *blockSize, (int)*ldx, (T*)x, (int)primme->nLocal, *blockSize, (int)*ldy, (T*)y);
    *ierr = 0; 
}

template <typename T>
static void myprevec(void *x, PRIMME_INT *ldx,  void *y, PRIMME_INT *ldy, int *blockSize, struct primme_params *primme, int *ierr) {
    PrimmeParams *pp = static_cast<PrimmeParams*>(primme);
    pp->prevec((int)primme->nLocal, *blockSize, (int)*ldx, (T*)x, (int)primme->nLocal, *blockSize, (int)*ldy, (T*)y);
    *ierr = 0; 
}

template <typename T>
struct Real { typedef T type;};
template <>
struct Real<std::complex<double> > { typedef double type;};
template <>
struct Real<std::complex<float> > { typedef float type;};

template <typename T>
static void myglobalSum(void *sendBuf, void *recvBuf, int *count, struct primme_params *primme, int *ierr) {
    PrimmeParams *pp = static_cast<PrimmeParams*>(primme);
    pp->globalSum(*count, static_cast<typename Real<T>::type*>(sendBuf), *count, static_cast<typename Real<T>::type*>(recvBuf));
    *ierr = 0;
}

template <typename T, typename R>
int my_primme(int lenEvals, R *evals,
            int len1Evecs, int len2Evecs, T *evecs,
            int lenResNorms, R *resNorms, 
            PrimmeParams *primme) {
   if (lenEvals < primme->numEvals) {
        PyErr_Format(PyExc_ValueError,
                     "Length of `evals' should be at least %d",
                     primme->numEvals);
        return -30;
   }
   if (primme->nLocal == 0)
        primme->nLocal = primme->n;
   if (len1Evecs < primme->nLocal || len2Evecs < primme->numEvals) {
        PyErr_Format(PyExc_ValueError,
                     "Size of `evecs' should be at least (%" PRIMME_INT_P ", %d)",
                     primme->nLocal, primme->numEvals);
        return -31;
   }
   if (lenResNorms < primme->numEvals) {
        PyErr_Format(PyExc_ValueError,
                     "Length of `resNorms' should be at least %d",
                     primme->numEvals);
        return -32;
   }
   primme->matrixMatvec = mymatvec<T>;
   if (primme->correctionParams.precondition) 
      primme->applyPreconditioner = myprevec<T>;
   if (primme->globalSum_set)
      primme->globalSumReal = myglobalSum<T>;
   int ret = tprimme(evals, evecs, resNorms, static_cast<primme_params*>(primme));
   return ret;
}

static int tprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds) { 
   return sprimme_svds(svals, svecs, resNorms, primme_svds);
}
static int tprimme_svds(float *svals, std::complex<float> *svecs, float *resNorms, primme_svds_params *primme_svds) {
   return cprimme_svds(svals, svecs, resNorms, primme_svds);
}
static int tprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds) { 
   return dprimme_svds(svals, svecs, resNorms, primme_svds);
}
static int tprimme_svds(double *svals, std::complex<double> *svecs, double *resNorms, primme_svds_params *primme_svds) {
   return zprimme_svds(svals, svecs, resNorms, primme_svds);
}

template <typename T>
static void myglobalSum_svds(void *sendBuf, void *recvBuf, int *count, struct primme_svds_params *primme_svds, int *ierr) {
    PrimmeSvdsParams *pp = static_cast<PrimmeSvdsParams*>(primme_svds);
    pp->globalSum(*count, static_cast<typename Real<T>::type*>(sendBuf), *count, static_cast<typename Real<T>::type*>(recvBuf));
    *ierr = 0;
}

template <typename T>
static void mymatvec_svds(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int *transpose, struct primme_svds_params *primme_svds, int *ierr) {
   PrimmeSvdsParams *pp = static_cast<PrimmeSvdsParams*>(primme_svds);
   PRIMME_INT m, n;
   if (*transpose == 0) {
      m = primme_svds->mLocal;
      n = primme_svds->nLocal;
   }
   else {
      n = primme_svds->mLocal;
      m = primme_svds->nLocal;
   }
   pp->matvec((int)n, *blockSize, (int)*ldx, (T*)x, (int)m, *blockSize, (int)*ldy, (T*)y, *transpose);
   *ierr = 0;
}

template <typename T>
static void myprevec_svds(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int *mode, struct primme_svds_params *primme_svds, int *ierr) {
   PrimmeSvdsParams *pp = static_cast<PrimmeSvdsParams*>(primme_svds);
   PRIMME_INT m=0;
   if (*mode == primme_svds_op_AtA) {
      m = primme_svds->nLocal;
   } else if (*mode ==  primme_svds_op_AAt) {
      m = primme_svds->mLocal;
   } else if (*mode == primme_svds_op_augmented) {
      m = primme_svds->mLocal + primme_svds->nLocal;
   }
   pp->prevec((int)m, *blockSize, (int)*ldx, (T*)x, (int)m, *blockSize, (int)*ldy, (T*)y, *mode);
   *ierr = 0;
}

template <typename T, typename R>
int my_primme_svds(int lenSvals, R *svals,
            int len1SvecsLeft, int len2SvecsLeft, T *svecsLeft,
            int len1SvecsRight, int len2SvecsRight, T *svecsRight,
            int lenResNorms, R *resNorms, 
            PrimmeSvdsParams *primme_svds) {
   if (lenSvals < primme_svds->numSvals) {
        PyErr_Format(PyExc_ValueError,
                     "Length of `svals' should be at least %d",
                     primme_svds->numSvals);
        return -30;
   }
   if (primme_svds->mLocal == 0)
        primme_svds->mLocal = primme_svds->m;
   if (primme_svds->nLocal == 0)
        primme_svds->nLocal = primme_svds->n;
   if (len1SvecsLeft < primme_svds->mLocal || len2SvecsLeft < primme_svds->numSvals) {
        PyErr_Format(PyExc_ValueError,
                     "Size of `svecsleft' should be at least (%" PRIMME_INT_P ", %d)",
                     primme_svds->mLocal, primme_svds->numSvals);
        return -31;
   }
   if (len1SvecsRight < primme_svds->nLocal || len2SvecsRight < primme_svds->numSvals) {
        PyErr_Format(PyExc_ValueError,
                     "Size of `svecsright' should be at least (%" PRIMME_INT_P ", %d)",
                     primme_svds->nLocal, primme_svds->numSvals);
        return -31;
   }
   if (lenResNorms < primme_svds->numSvals) {
        PyErr_Format(PyExc_ValueError,
                     "Length of `resNorms' should be at least %d",
                     primme_svds->numSvals);
        return -32;
   }
   primme_svds->matrixMatvec = mymatvec_svds<T>;
   if (primme_svds->precondition) 
      primme_svds->applyPreconditioner = myprevec_svds<T>;
   if (primme_svds->globalSum_set)
      primme_svds->globalSumReal = myglobalSum_svds<T>;
   T *svecs = new T[(primme_svds->nLocal+primme_svds->mLocal)*(primme_svds->numOrthoConst+primme_svds->numSvals)];
   copy_matrix(svecsLeft, primme_svds->mLocal, primme_svds->numOrthoConst,
         (PRIMME_INT)len1SvecsLeft, svecs, primme_svds->mLocal);
   copy_matrix(svecsRight, primme_svds->nLocal, primme_svds->numOrthoConst,
         (PRIMME_INT)len1SvecsRight, &svecs[primme_svds->numOrthoConst*primme_svds->mLocal],
         primme_svds->nLocal);
   int ret = tprimme_svds(svals, svecs, resNorms, static_cast<primme_svds_params*>(primme_svds));
   copy_matrix(&svecs[primme_svds->mLocal*primme_svds->numOrthoConst],
         primme_svds->mLocal, primme_svds->numSvals,
         primme_svds->mLocal, &svecsLeft[len1SvecsLeft*primme_svds->numOrthoConst], (PRIMME_INT)len1SvecsLeft);
   copy_matrix(&svecs[primme_svds->mLocal*(primme_svds->numOrthoConst
            +primme_svds->initSize)
         + primme_svds->nLocal*primme_svds->numOrthoConst], primme_svds->nLocal,
         primme_svds->initSize, primme_svds->nLocal, &svecsRight[len1SvecsRight*primme_svds->numOrthoConst], (PRIMME_INT)len1SvecsRight);
   delete [] svecs;

   return ret;
}
%}

%template (sprimme) my_primme<float,float>;
%template (cprimme) my_primme<std::complex<float>,float>;
%template (dprimme) my_primme<double,double>;
%template (zprimme) my_primme<std::complex<double>,double>;
%template (sprimme_svds) my_primme_svds<float,float>;
%template (cprimme_svds) my_primme_svds<std::complex<float>,float>;
%template (dprimme_svds) my_primme_svds<double,double>;
%template (zprimme_svds) my_primme_svds<std::complex<double>,double>;


%feature("director") PrimmeParams;
%feature("director") PrimmeSvdsParams;
%feature("director:except") {
    if ($error != NULL) {
        throw Swig::DirectorMethodException();
    }
}

%define DOCSTRING_PrimmeParams
"Abstract class to specify the eigenvalue problem and the options for calling
dprimme and zprimme.

Example
-------
>>> import Primme, scipy.sparse, numpy as np
>>> A = scipy.sparse.spdiags(range(100), [0], 100, 100) # sparse diag. matrix
>>> class PP(Primme.PrimmeParams):
... 	def __init__(self):
... 		Primme.PrimmeParams.__init__(self)
... 	def matvec(self, X):
... 		return A*X
>>> pp = PP()
>>> pp.n = A.shape[0] # set problem dimension
>>> pp.numEvals = 3   # set number of eigenvalues
>>> pp.target = Primme.primme_largest # find the largest eigenvalues
>>> pp.eps = 1e-6     # residual norm tolerance
>>> evals = np.zeros(pp.numEvals)                    # store eigenvalues
>>> evecs = np.zeros((pp.n, pp.numEvals), order='F') # store eigenvectors
>>> norms = np.zeros(pp.numEvals)                    # store residual norms
>>> ret = Primme.dprimme(evals, evecs, norms, pp) # call the solver
>>> ret  # error code, 0 is success!
0
>>> pp.initSize # number of converged eigenpairs
3
>>> evals[0:pp.initSize] # converged values 
array([ 99.,  98.,  97.])
>>> # Time in seconds and A*v times that took
>>> pp.stats.elapsedTime, pp.stats.numMatvecs # doctest: +SKIP
0.5, 110
"
%enddef
%feature("docstring") PrimmeParams DOCSTRING_PrimmeParams;

%define DOCSTRING_PrimmeSvdsParams
"Abstract class to specify the eigenvalue problem and the options for calling
dprimme_svds and zprimme_svds.

Example
-------
>>> import Primme, scipy.sparse, numpy as np
>>> A = scipy.sparse.spdiags(range(1,11), [0], 100, 10) # sparse diag. rect. matrix
>>> class PSP(Primme.PrimmeSvdsParams):
... 	def __init__(self):
... 		Primme.PrimmeSvdsParams.__init__(self)
... 	def matvec(self, X, transpose):
... 		return A*X if transpose == 0 else A.H*X
>>> pp = PSP()
>>> pp.m, pp.n = A.shape # set problem dimensions
>>> pp.numSvals = 3   # set number of singular values to seek
>>> pp.target = Primme.primme_svds_smallest # find the smallest singular values
>>> pp.eps = 1e-6     # residual norm tolerance
>>> svals = np.zeros(pp.numSvals)                     # store singular values
>>> svecsl = np.zeros((pp.m, pp.numSvals), order='F') # store left singular vectors
>>> svecsr = np.zeros((pp.n, pp.numSvals), order='F') # store right singular vectors
>>> norms = np.zeros(pp.numSvals)                     # store residual norms
>>> ret = Primme.dprimme_svds(svals, svecsl, svecsr, norms, pp) # call the solver
>>> ret  # error code, 0 is success!
0
>>> pp.initSize # number of converged singular pairs
3
>>> svals[0:pp.initSize] # converged singular values 
array([ 1.,  2.,  3.])
>>> # Time in seconds and A*v times that took
>>> pp.stats.elapsedTime, pp.stats.numMatvecs # doctest: +SKIP
0.7, 94
"
%enddef
%feature("docstring") PrimmeSvdsParams DOCSTRING_PrimmeSvdsParams;


%include "../include/primme.h"
%include "../include/primme_eigs.h"
%include "../include/primme_svds.h"
%include "primmew.h"

%pythoncode "wrappers.py"
%pythoncode %{
PrimmeParams.targetShifts = property(PrimmeParams._get_targetShifts, PrimmeParams._set_targetShifts)
PrimmeSvdsParams.targetShifts = property(PrimmeSvdsParams._get_targetShifts, PrimmeSvdsParams._set_targetShifts)
%}
