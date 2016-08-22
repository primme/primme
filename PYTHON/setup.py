#! /usr/bin/env python
from __future__ import division, print_function

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.
numpy_include = numpy.get_include()

# Array extension module
_Primme = Extension("_Primme",
                   ["primme_wrap.cxx", "primmew.cxx"],
                   include_dirs = [numpy_include],
                   library_dirs = [".."],
                   libraries = ["primme", "blas", "lapack"],
                   #extra_compile_args = ["-g", "-O0"]
                   )

# NumyTypemapTests setup
setup(name        = "NumpyPrimme",
      description = "PRIMME wrapper for Numpy",
      author      = "Eloy Romero",
      py_modules  = ["Primme"],
      ext_modules = [_Primme]
      )
