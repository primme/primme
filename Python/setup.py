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
   
   try:
      import numpy
   except:
      raise Exception("numpy not installed; please, install numpy and scipy before primme")
   else:
      extra_options = get_numpy_options()
   
   _Primme = Extension(name="primme",
                      sources=["primme.pyx"],
                      **extra_options
                      #extra_compile_args = ["-g", "-O0", "-Wall", "-Wextra"]
                      )
   
   setup(ext_modules = [_Primme])

if __name__ == '__main__':
   setup_package()
