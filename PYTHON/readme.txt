-----------------------------------------------------------------------------
                 Primme.py: A Python Interface for PRIMME
   
                Copyright (C) 2015 College of William & Mary,
   James R. McCombs, Eloy Romero Alcalde, Andreas Stathopoulos, Lingfei Wu
-----------------------------------------------------------------------------
 
   This file is part of PRIMME.

   PRIMME is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   PRIMME is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


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

   make tests

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
