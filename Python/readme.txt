-----------------------------------------------------------------------------
                 Primme.py: A Python Interface for PRIMME
   
Copyright (c) 2016, College of William & Mary
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the College of William & Mary nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

PRIMME: https://github.com/primme/primme
Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
-----------------------------------------------------------------------------

Primme.py is a Python interface for the free software PRIMME (PReconditioned
Iterative MultiMethod Eigensolver), which finds a number of eigenvalues and 
their corresponding eigenvectors of a real symmetric, or complex hermitian
matrix A. And also recently it support singular values. It is a useful tool
for both non-experts and experts to easily call PRIMME. Largest, smallest and
interior eigenvalues and singular values are supported.  Preconditioning can
be used to accelerate convergence. 

-----------------------------------------------------------------------------
1. Directory Structure 
-----------------------------------------------------------------------------
PRIMME/PYTHON/

> ls 
readme.txt            <- this file
Makefile              <- make with actions to build, test and clean
numpy.i and
pyfragments.swg       <- SWIG files distributed by Numpy (in dir tools/swig)
primme.i              <- SWIG file with PRIMME interface description
primmew.h and
primmew.cxx           <- C++ class encapsulating PRIMME and PRIMME_SVDS
wrappers.py           <- implementation of eigs and svds
primme_wrap.h and
primme_wrap.cxx and
Python.py             <- files generated from primme.i
_Primme.so (generated)<- shared library with PRIMME and PYTHON interface
setup.py              <- disutils script to build _Primme.so
tests.py              <- tests for the python interface
examples.py           <- few examples with eigs and svds

-----------------------------------------------------------------------------
2. _Primme.so compilation 
-----------------------------------------------------------------------------

The python interface is composed by Primme.py and _Primme.so. Follow this
instruction to build _Primme.so.

Previously, make sure that libprimme.a is suitable for dynamic linking (i.e.,
it has been compiled with -fPIC). For instance do the following on the root
directory of PRIMME:

  make clean lib CFLAGS="-fPIC -O2"

After that, execute one of the following commands.

a) To generate _Primme.so in the current working path do:

   make all

or

   python setup.py build_ext -i

b) To install Primme with other python packages do:

   python setup.py install

To verify the installation try to run the tests:

   make test

or execute the examples:

   python examples.py

-----------------------------------------------------------------------------
3. Interface description
-----------------------------------------------------------------------------

Primme module offers eigs and svds as high level functions to compute
eigenpairs and singular triples. Too see their description, please refer
to the documentation or use help from an interactive python:

$ python
>>> import Primme
>>> help(Primme.eigs)
>>> help(Primme.svds)

The abstract classes PrimmeParams and PrimmeSvdsParams in Primme module
encapsulates primme_params and primme_svds_params respectively from the C
interface. To see how they works, look for the definitions of eigs and svds
in wrappers.py.
