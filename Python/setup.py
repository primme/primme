#! /usr/bin/env python
from __future__ import division, print_function

# System imports
from distutils.core import *

# Third-party modules - we depend on numpy for everything
import numpy
try:
    from numpy.distutils.system_info import get_info
except:
    from numpy.__config__ import get_info

from codecs import open
from os import path

# Obtain the numpy include directory
numpy_include = numpy.get_include()

# PRIMME source files
primme_source_files = []
for f in (
"eigs/auxiliary_eigs",
"eigs/main_iter",
"eigs/ortho",
"eigs/solve_projection",
"eigs/init",
"eigs/restart",
"eigs/locking",
"eigs/update_projection",
"eigs/correction",
"eigs/primme_interface",
"eigs/globalsum",
"eigs/primme",
"eigs/update_W",
"eigs/convergence",
"eigs/inner_solve",
"eigs/primme_f77",
"eigs/factorize",
"linalg/blaslapack",
"linalg/wtime",
"linalg/auxiliary",
"svds/primme_svds",
"svds/primme_svds_f77",
"svds/primme_svds_interface"):
    for v in ("double", "doublecomplex", "float", "floatcomplex"):
        primme_source_files.append("primme/src/%s%s.c" % (f, v))


# Obtain BLAS/LAPACK linking options
lapack_info = get_info('lapack_opt')
blas_info = get_info('blas_opt')
blaslapack_libraries = lapack_info.get('libraries', []) + blas_info.get('libraries', [])
blaslapack_library_dirs = lapack_info.get('library_dirs', []) + blas_info.get('library_dirs', [])
blaslapack_extra_link_args = lapack_info.get('extra_link_args', []) + blas_info.get('extra_link_args', [])
if not blaslapack_libraries and not blaslapack_extra_link_args:
    blaslapack_libraries = ['lapack', 'blas']


# Get the long description from the README file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


# Array extension module
_Primme = Extension("_Primme",
                   ["primme_wrap.cxx", "primmew.cxx"] + primme_source_files,
                   include_dirs = [numpy_include, "primme/include", "primme/src/include"],
                   library_dirs = blaslapack_library_dirs,
                   libraries = blaslapack_libraries,
                   extra_link_args = blaslapack_extra_link_args
                   #extra_compile_args = ["-g", "-O0", "-Wall", "-Wextra"]
                   )

# NumyTypemapTests setup
setup(name        = "primme",
      version     = "2.1.1",
      description = "PRIMME wrapper for Python",
      long_description = long_description,
      url         = "https://github.com/primme/primme",
      author      = "Eloy Romero Alcalde, Andreas Stathopoulos and Lingfei Wu",
      author_email = "eloy@cs.wm.edu",
      license     = "BSD",
      classifiers=[
# How mature is this project? Common values are
#   3 - Alpha
#   4 - Beta
#   5 - Production/Stable
      'Development Status :: 4 - Beta',

# Indicate who your project is intended for
      'Intended Audience :: Science/Research',
      'Intended Audience :: Developers',

# Pick your license as you wish (should match "license" above)
      'License :: OSI Approved :: BSD License',

      'Topic :: Software Development',
      'Topic :: Scientific/Engineering',
      'Intended Audience :: Information Technology',
      'Operating System :: Microsoft :: Windows',
      'Operating System :: POSIX',
      'Operating System :: Unix',
      'Operating System :: MacOS',

# Specify the Python versions you support here. In particular, ensure
# that you indicate whether you support Python 2, Python 3 or both.
      'Programming Language :: C',
      'Programming Language :: Python :: 2',
      'Programming Language :: Python :: 2.6',
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3',
      'Programming Language :: Python :: 3.2',
      'Programming Language :: Python :: 3.3',
      'Programming Language :: Python :: 3.4',
      ],
      keywords = "eigenvalues singular values Davidson-type high-performance large-scale matrix",
      setup_requires = ['numpy', 'scipy'],
      install_requires = ['numpy', 'scipy'],
      py_modules  = ["Primme"],
      ext_modules = [_Primme]
      )
