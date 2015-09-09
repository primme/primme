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


%module(docstring="",directors="1") Primme

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


%feature("shadow") primme_params_w::getX() %{
def getX(self):
   if self.__kind == 1:
      return self.getXd()
   elif self.__kind == 3:
      return self.getXz()
%}

%feature("shadow") primme_params_w::setY() %{
def setY(self, Y):
   if self.__kind == 1:
      self.setYd(Y)
   elif self.__kind == 3:
      self.setYz(Y)
%}

%extend primme_params_w {
   void getX() {}
   void setY() {}
}

%feature("director") primme_params_w;
%include "../PRIMMESRC/COMMONSRC/primme.h"
%include "primmew.h"


