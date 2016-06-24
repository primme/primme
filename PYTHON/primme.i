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
 * File: primme.i
 * 
 * Purpose - SWIG description of PYTHON interface to PRIMME.
 * 
 ******************************************************************************/

%define DOCSTRING
"Find a few eigenvectors and eigenvalues of a matrix.
Uses PRIMME: https://github.com/primme/primme

Example
-------
>>> import Primme, numpy as np
>>> from scipy.sparse import *
>>> 
>>> # A = [ 2  1  0 ...
>>> #      -1  2 -1 0 ...
>>> #       0 -1  2 -1 0 ... ]
>>> a = np.ones(10)
>>> A = spdiags(np.array([a*(-1.), a*2., a*(-1.)]), np.array([-1, 0, 1]), 10, 10)
>>> 
>>> class PPd(Primme.primme_params_w):
... 	def __init__(self):
... 		Primme.primme_params_w.__init__(self)
... 	def matvec(self):
... 		self.setY(A*self.getX())
>>> pp = PPd()
>>> pp.n = A.shape[0]
>>> pp.maxBasisSize = 3
>>> pp.minRestartSize = 1
>>> pp.numEvals = 3
>>> pp.restartingParams.maxPrevRetain = 1
>>> pp.set_method(Primme.DYNAMIC)
>>> pp.display()
>>> evals = np.zeros(pp.numEvals)
>>> evecs = np.zeros((pp.n, pp.numEvals))
>>> norms = np.zeros(pp.numEvals)
>>> print Primme.dprimme(evals, evecs, norms, pp)
>>> print pp.initSize, evals, norms
>>> 
>>> class PPz(Primme.primme_params_w):
... 	def __init__(self, matrix=None):
... 		Primme.primme_params_w.__init__(self)
... 		self.mymatrix = matrix
... 	def matvec(self):
... 		self.setY(self.mymatrix*self.getX())
>>> 
>>> a = np.ones(10, complex)
>>> A = spdiags(np.array([a*(-1.), a*2., a*(-1.)]), np.array([-1, 0, 1]), 10, 10)
>>> pp = PPz(A)
>>> pp.n = A.shape[0]
>>> pp.maxBasisSize = 3
>>> pp.minRestartSize = 1
>>> pp.numEvals = 3
>>> pp.set_method(Primme.DYNAMIC)
>>> pp.display()
>>> evals = np.zeros(pp.numEvals)
>>> evecs = np.zeros((pp.n, pp.numEvals), complex)
>>> norms = np.zeros(pp.numEvals)
>>> print Primme.zprimme(evals, evecs, norms, pp)
>>> print pp.initSize, evals, norms, pp.stats.numMatvecs"
%enddef

%module(docstring=DOCSTRING,directors="1") Primme

%pythoncode %{
__all__ = ['primme_params_w', 'dprimme', 'zprimme', 'eigsh', 'PRIMMEError', 'PRIMMENoConvergence', 'Arnoldi', 'DEFAULT_METHOD', 'DEFAULT_MIN_MATVECS', 'DEFAULT_MIN_TIME', 'DYNAMIC', 'GD', 'GD_Olsen_plusK', 'GD_plusK', 'JDQMR', 'JDQMR_ETol', 'JDQR', 'JD_Olsen_plusK', 'LOBPCG_OrthoBasis', 'LOBPCG_OrthoBasis_Window', 'RQI', 'SUBSPACE_ITERATION', 'primme_adaptive', 'primme_adaptive_ETolerance', 'primme_closest_abs', 'primme_closest_geq', 'primme_closest_leq', 'primme_decreasing_LTolerance', 'primme_dtr', 'primme_full_LTolerance', 'primme_init_default', 'primme_init_krylov', 'primme_init_random', 'primme_init_user', 'primme_largest', 'primme_largest_abs', 'primme_params', 'primme_proj_RR', 'primme_proj_default', 'primme_proj_harmonic', 'primme_proj_refined', 'primme_smallest', 'primme_thick']
%}

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
}
%init %{
  import_array();
%}

// Global ignores
%ignore PRIMME_MAX_NAME_LENGTH;
%ignore stackTraceNode;
%ignore primme_valloc;
%ignore primme_calloc;
%ignore primme_malloc;
%ignore primme_display_params;
%ignore primme_set_method;
%ignore primme_initialize;
%ignore primme_seq_globalSumDouble;
%ignore primme_PushErrorMessage;
%ignore primme_PrintStackTrace;
%ignore primme_DeleteStackTrace;
%ignore primme_Free;

%rename (dprimme) my_dprimme;
%exception my_dprimme {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (zprimme) my_zprimme;
%exception my_zprimme {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%define %numpy_typemaps_ext(DATA_TYPE, DATA_TYPECODE, DIM_TYPE)

/* Typemap suite for (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* IN_FARRAY2D)
 */
%typemap(directorin,
         fragment="NumPy_Fragments")
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* IN_FARRAY2D)
  (PyArrayObject* array=NULL, int is_new_object=1)
{
  npy_intp dims[2] = { $1, $2 };
  PyObject* obj = PyArray_SimpleNewFromData(2, dims, DATA_TYPECODE, (void*)($3));
  //$input = Py_BuildValue("(O)", obj);
  $input = obj;
}

%typemap(directorfreearg)
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* IN_FARRAY2D)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}


/* Typemap suite for (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* INPLACE_ARRAY2)
 */
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY,
           fragment="NumPy_Macros")
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* OUT_FARRAY2D)
{
  $1 = is_array($input) && PyArray_EquivTypenums(array_type($input),
                                                 DATA_TYPECODE);
}
%typemap(in,numinputs=0)
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* OUT_FARRAY2D)
{}

%typemap(directorargout,
         fragment="NumPy_Fragments")
  (DIM_TYPE DIM1, DIM_TYPE DIM2, DATA_TYPE* OUT_FARRAY2D)
  (PyArrayObject* array=NULL, PyObject* o=NULL)
{
  o = $result;
  if (!is_array(o) || !PyArray_EquivTypenums(array_type(o), DATA_TYPECODE))
     Swig::DirectorMethodException::raise("No valid type for $3_name");
  array = obj_to_array_no_conversion(o, DATA_TYPECODE);
  if (!array || !require_dimensions(array,2) || !require_contiguous(array) ||
      !require_native(array))
          Swig::DirectorMethodException::raise("No valid type for $3_name");
  if (($1) != (DIM_TYPE) array_size(array,0) ||
      ($2) != (DIM_TYPE) array_size(array,1))
          {Swig::DirectorMethodException::raise("No valid dimensions for $3_name");}
  memcpy(($3), array_data(array), sizeof(DATA_TYPE)*($1)*($2));
}
%enddef    /* %numpy_typemaps_ext() macro */

%numpy_typemaps_ext(double            , NPY_DOUBLE   , int)
%numpy_typemaps_ext(std::complex<double>, NPY_CDOUBLE, int)


%apply (int DIM1, double* INPLACE_ARRAY1) {
   (int lenEvals, double* evals),
   (int lenResNorms, double* resNorms)};
%apply (int DIM1, int DIM2, double* INPLACE_FARRAY2) {
   (int len1Evecs, int len2Evecs, double* evecs)};
%apply (int DIM1, int DIM2, std::complex<double>* INPLACE_FARRAY2) {
   (int len1Evecs, int len2Evecs, std::complex<double>* evecs)};
%apply (int* DIM1, int* DIM2, double** ARGOUTVIEW_FARRAY2) {
   (int *len1X, int *len2X, double **x)};
%apply (int* DIM1, int* DIM2, std::complex<double>** ARGOUTVIEW_FARRAY2) {
   (int *len1X, int *len2X, std::complex<double> **x)};
%apply (int DIM1, int DIM2, double* IN_FARRAY2) {
   (int len1Y, int len2Y, double* y)};
%apply (int DIM1, int DIM2, std::complex<double>* IN_FARRAY2) {
   (int len1Y, int len2Y, std::complex<double>* y)};

%apply (int DIM1, int DIM2, double* IN_FARRAY2D) {
   (int len1YD, int len2YD, double* yd)};
%apply (int DIM1, int DIM2, double* OUT_FARRAY2D) {
   (int len1XD, int len2XD, double* xd)};
%apply (int DIM1, int DIM2, std::complex<double>* IN_FARRAY2D) {
   (int len1YD, int len2YD, std::complex<double>* yd)};
%apply (int DIM1, int DIM2, std::complex<double>* OUT_FARRAY2D) {
   (int len1XD, int len2XD, std::complex<double>* xd)};


%inline %{
int my_dprimme(int lenEvals, double *evals,
            int len1Evecs, int len2Evecs, double *evecs,
            int lenResNorms, double *resNorms, 
            primme_params_w *primme) {
   if (lenEvals < primme->numEvals) {
        PyErr_Format(PyExc_ValueError,
                     "Length of `evals' should be at least %d",
                     primme->numEvals);
        return -30;
   }
   if (len1Evecs < primme->nLocal || len2Evecs < primme->numEvals) {
        PyErr_Format(PyExc_ValueError,
                     "Size of `evecs' should be at least (%d, %d)",
                     primme->nLocal, primme->numEvals);
        return -31;
   }
   if (lenResNorms < primme->numEvals) {
        PyErr_Format(PyExc_ValueError,
                     "Length of `resNorms' should be at least %d",
                     primme->numEvals);
        return -32;
   }
   primme->__kind = 1;
   int ret = dprimme(evals, evecs, resNorms, static_cast<primme_params*>(primme));
   primme_Free(static_cast<primme_params*>(primme));
   return ret;
}
%}

%inline %{
int my_zprimme(int lenEvals, double *evals,
            int len1Evecs, int len2Evecs, std::complex<double> *evecs,
            int lenResNorms, double *resNorms, 
            primme_params_w *primme) {
   if (lenEvals < primme->numEvals) {
        PyErr_Format(PyExc_ValueError,
                     "Length of `evals' should be at least %d",
                     primme->numEvals);
        return -30;
   }
   if (len1Evecs < primme->nLocal || len2Evecs < primme->numEvals) {
        PyErr_Format(PyExc_ValueError,
                     "Size of `evecs' should be at least (%d, %d)",
                     primme->nLocal, primme->numEvals);
        return -31;
   }
   if (lenResNorms < primme->numEvals) {
        PyErr_Format(PyExc_ValueError,
                     "Length of `resNorms' should be at least %d",
                     primme->numEvals);
        return -32;
   }
   primme->__kind = 3;
   int ret = zprimme(evals, (Complex_Z*)evecs, resNorms, static_cast<primme_params*>(primme));
   primme_Free(static_cast<primme_params*>(primme));
   return ret;
}
%}

%feature("director") primme_params_w;
%include "../PRIMMESRC/COMMONSRC/primme.h"
%include "primmew.h"

%pythoncode "wrappers.py"
