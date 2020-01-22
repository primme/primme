#! /usr/bin/env python

from codecs import open
from os import path
import sys

primme_source_files = [
   "primme/src/eigs/auxiliary_eigs.cpp",
   "primme/src/eigs/auxiliary_eigs_normal.cpp",
   "primme/src/eigs/convergence.cpp",
   "primme/src/eigs/correction.cpp",
   "primme/src/eigs/factorize.cpp",
   "primme/src/eigs/init.cpp",
   "primme/src/eigs/inner_solve.cpp",
   "primme/src/eigs/main_iter.cpp",
   "primme/src/eigs/ortho.cpp",
   "primme/src/eigs/primme_c.cpp",
   "primme/src/eigs/primme_f77.cpp",
   "primme/src/eigs/primme_interface.cpp",
   "primme/src/eigs/restart.cpp",
   "primme/src/eigs/solve_projection.cpp",
   "primme/src/eigs/update_projection.cpp",
   "primme/src/eigs/update_W.cpp",
   "primme/src/linalg/auxiliary.cpp",
   "primme/src/linalg/blaslapack.cpp",
   "primme/src/linalg/magma_wrapper.cpp",
   "primme/src/linalg/memman.cpp",
   "primme/src/linalg/wtime.cpp",
   "primme/src/svds/primme_svds_c.cpp",
   "primme/src/svds/primme_svds_f77.cpp",
   "primme/src/svds/primme_svds_interface.cpp"
]

def get_numpy_options():
   # Third-party modules - we depend on numpy for everything
   import numpy
   try:
       from numpy.distutils.system_info import get_info
   except:
       from numpy.__config__ import get_info
   
   # Obtain the numpy include directory
   numpy_include = numpy.get_include()

   # Obtain BLAS/LAPACK linking options
   lapack_info = get_info('lapack_opt')
   blas_info = get_info('blas_opt')
   using_atlas = False
   using_f77blas = False
   using_lapack = False
   for l in lapack_info.get('libraries', []) + blas_info.get('libraries', []):
      if "atlas" in l: using_atlas = True
      if "f77blas" in l: using_f77blas = True
      if "lapack" in l: using_lapack = True
   if using_atlas and (not using_f77blas or not using_lapack):
      lapack_info = get_info('atlas')
      # ATLAS notices an incomplete LAPACK by not setting language to f77
      complete_lapack = lapack_info.get('language', "") == "f77"
      if complete_lapack:
         blas_info = {}
      else:
         # If ATLAS has an incomplete LAPACK, use a regular one
         blas_info = get_info('atlas_blas')
         lapack_info = get_info('lapack')
   
   blaslapack_libraries = lapack_info.get('libraries', []) + blas_info.get('libraries', [])
   blaslapack_library_dirs = lapack_info.get('library_dirs', []) + blas_info.get('library_dirs', [])
   blaslapack_extra_link_args = lapack_info.get('extra_link_args', []) + blas_info.get('extra_link_args', [])
   if not blaslapack_libraries and not blaslapack_extra_link_args:
       blaslapack_libraries = ['lapack', 'blas']

   return dict(
                   include_dirs = [numpy_include, "primme/include", "primme/src/include"],
                   library_dirs = blaslapack_library_dirs,
                   libraries = blaslapack_libraries,
                   extra_link_args = blaslapack_extra_link_args,
                   define_macros = [("NDEBUG", None), ("F77UNDERSCORE", None)]
   )

def get_basic_options():
   return dict(
                   include_dirs = ["primme/include", "primme/src/include"],
                   define_macros = [("NDEBUG", None), ("F77UNDERSCORE", None)]
   )

def setup_package():
   import sys
   from distutils.core import setup
   from distutils.extension import Extension
   
   try:
      import numpy
   except:
      if ('bdist_wheel' in sys.argv[1:] or 'install' in sys.argv[1:] or
          'build_ext' in sys.argv[1:]):
         raise Exception("numpy not installed; please, install numpy before primme")
      extra_options = get_basic_options()
   else:
      extra_options = get_numpy_options()
   
   # Get the long description from the README file
   here = path.abspath(path.dirname(__file__))
   with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
       long_description = f.read()
   
   
   # Array extension module
   _Primme = Extension("primme",
                      ["primme.cpp"] + primme_source_files,
                      **extra_options
                      #extra_compile_args = ["-g", "-O0", "-Wall", "-Wextra"]
                      )
   
   # NumyTypemapTests setup
   setup(name        = "primme",
         version     = "3.0.3",
         description = "PRIMME wrapper for Python",
         long_description = long_description,
         long_description_content_type = "text/x-rst",
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
         'Programming Language :: Python :: 3'
         ],
         keywords = "eigenvalues singular values generalized Hermitian symmetric Davidson-type high-performance large-scale matrix",
         setup_requires = ['numpy', 'scipy'],
         install_requires = ['future', 'numpy', 'scipy'],
         ext_modules = [_Primme]
         )

if __name__ == '__main__':
   setup_package()
