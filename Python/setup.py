#! /usr/bin/env python

from codecs import open
from os import path
import sys

def get_numpy_options():
   # Obtain the numpy include directory
   import numpy
   numpy_include = numpy.get_include()

   r = dict(
                   include_dirs = [numpy_include, "primme/include", "primme/src/include"],
                   libraries = ['lapack', 'blas'],
   )

   # Link dynamically on Windows and statically otherwise
   if sys.platform == 'win32':
      r['libraries'] = ['primme'] + r['libraries']
   else:
      r['extra_objects'] = ['../lib/libprimme.a']

   return r

def setup_package():
   import sys
   from setuptools import setup, Extension
   from Cython.Build import cythonize
   
   try:
      import numpy
   except:
      raise Exception("numpy not installed; please, install numpy and scipy before primme")
   else:
      extra_options = get_numpy_options()
   
   # Get the long description from the README file
   here = path.abspath(path.dirname(__file__))
   with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
       long_description = f.read()
   
   
   # Array extension module
   _Primme = Extension("primme",
                      ["primme.pyx"],
                      **extra_options
                      #extra_compile_args = ["-g", "-O0", "-Wall", "-Wextra"]
                      )
   
   # NumyTypemapTests setup
   setup(name        = "primme",
         version     = "3.1.0",
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
         install_requires = ['future', 'numpy', 'scipy'],
         ext_modules = cythonize([_Primme])
         )

if __name__ == '__main__':
   setup_package()
