import sys

def get_numpy_options():
   # Obtain the numpy include directory
   import numpy
   numpy_include = numpy.get_include()

   # Check for blas/lapack libraries
   r = {}
   find_lapack = False
   try:
      import pkgconfig
      for pkg_name in [ prefix + lib for prefix in ["", "lib"] for lib in ['openblas', 'mkl', 'blis', 'lapack'] ]:
         if not pkgconfig.exists(pkg_name): continue
         r = pkgconfig.parse(pkg_name)
         find_lapack = True
         break
   except:
      pass
   if not find_lapack:
      try:
         import scipy_openblas32
         r = dict(
            library_dirs = [scipy_openblas32.get_lib_dir()],
            libraries = [scipy_openblas32.get_library()]
         )
         find_lapack = True
      except:
         pass
   if not find_lapack:
      try:
         import scipy_openblas64
         r = dict(
            define_macros = [("PRIMME_BLASINT_SIZE","64")],
            library_dirs = [scipy_openblas64.get_lib_dir()],
            libraries = [scipy_openblas64.get_library()]
         )
         find_lapack = True
      except:
         pass
   if not find_lapack:
      r = dict( libraries = ['lapack', 'blas'] )
          
   r['define_macros'] = r.get('define_macros', []) + [("NDEBUG", None), ("F77UNDERSCORE", None)]
   r['include_dirs'] = [numpy_include, "src/primme/include", "src/primme/svds", "src/primme/eigs", "src/primme/linalg"]

   return r

from setuptools import setup, Extension

try:
   import numpy
   extra_options = get_numpy_options()
except:
   extra_options = {}

_Primme = Extension(name="primme",
                   sources=["primme.pyx",
                            "src/primme/svds/primme_svds_f77.c",
                            "src/primme/svds/primme_svds_interface.c",
                            "src/primme/svds/primme_svds_c.c",
                            "src/primme/eigs/primme_f77.c",
                            "src/primme/eigs/auxiliary_eigs.c",
                            "src/primme/eigs/convergence.c",
                            "src/primme/eigs/restart.c",
                            "src/primme/eigs/ortho.c",
                            "src/primme/eigs/update_projection.c",
                            "src/primme/eigs/factorize.c",
                            "src/primme/eigs/primme_interface.c",
                            "src/primme/eigs/inner_solve.c",
                            "src/primme/eigs/update_W.c",
                            "src/primme/eigs/init.c",
                            "src/primme/eigs/correction.c",
                            "src/primme/eigs/main_iter.c",
                            "src/primme/eigs/auxiliary_eigs_normal.c",
                            "src/primme/eigs/solve_projection.c",
                            "src/primme/eigs/primme_c.c",
                            "src/primme/linalg/wtime.c",
                            "src/primme/linalg/magma_wrapper.c",
                            "src/primme/linalg/memman.c",
                            "src/primme/linalg/cublas_wrapper.c",
                            "src/primme/linalg/blaslapack.c",
                            "src/primme/linalg/auxiliary.c",
                   ],
                   **extra_options
                   #extra_compile_args = ["-g", "-O0", "-Wall", "-Wextra"]
                   )

setup(ext_modules = [_Primme])
