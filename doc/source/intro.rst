
PRIMME: PReconditioned Iterative MultiMethod Eigensolver
--------------------------------------------------------

PRIMME, pronounced as *prime*, finds a number of eigenvalues and their corresponding eigenvectors of a 
real symmetric, or complex hermitian matrix A. Largest, smallest and interior 
eigenvalues are supported. Preconditioning can be used to accelerate 
convergence. 
PRIMME is written in C99, but complete interfaces are provided for Fortran 77 and MATLAB.
  
Changelog
---------

Changes in PRIMME 1.2.1 (released on ???):

* Added MATLAB interface to full PRIMME functionality.

* Support for BLAS_/LAPACK_ with 64bits integers (``-DPRIMME_BLASINT_SIZE=64``).

* Simplified configuration of Make_flags and Make_links (removed ``TOP`` variable 
  and replaced defines ``NUM_SUM`` and ``NUM_IBM`` by ``F77UNDERSCORE``).

* Replaced directories ``DTEST`` and ``ZTEST`` by ``TEST``, that has:

  * ``driver.c``: read matrices in MatrixMarket format and PETSc binary and
    call PRIMME with the parameters specified in a file; support
    complex arithmetic and MPI and can use PETSc preconditioners.
  * ``ex*.c`` and ``ex*.f``: small, didactic examples of usage in C and Fortran
    and in parallel (with PETSc).

* Fixed a few minor bugs and improved documentation (specially the F77 interface).

* Using Sphinx_ to manage documentation.

Changes in PRIMME 1.2 (released on December 21, 2014):

* A Fortran compiler is no longer required for building the PRIMME library.
  Fortran programs can still be linked to PRIMME's F77 interface.

* Fixed some uncommon issues with the F77 interface.

* PRIMME can be called now multiple times from the same program.

* Performance improvements in the QMR inner solver, especially for 
  complex arithmetic.

* Fixed a couple of bugs with the locking functionality. 

  * In certain extreme cases where all eigenvalues of a matrix were needed.
  * The order of selecting interior eigenvalues.

  The above fixes have improved robustness and performance. 

* PRIMME now assigns unique random seeds per parallel process 
  for up to 4096^3  (140 trillion processes).

* For the DYNAMIC method, fixed issues with initialization and 
  synchronization decisions across multiple processes.

* Fixed uncommon library interface bugs, coordinated better
  setting the method and the user setting of parameters, and improved 
  the interface in the sample programs and makefiles.

* Other performance and documentation improvements.


Citing this code 
---------------- 

Please cite:

.. [r1] A. Stathopoulos and J. R. McCombs PRIMME: *PReconditioned Iterative
   MultiMethod Eigensolver: Methods and software description*, ACM
   Transaction on Mathematical Software Vol. 37, No. 2, (2010),
   21:1-21:30.

More information on the algorithms and research that led to this
software can be found in the rest of the papers. The work has been
supported by a number of grants from the National Science Foundation.

.. [r2] A. Stathopoulos, *Nearly optimal preconditioned methods for hermitian
   eigenproblems under limited memory. Part I: Seeking one eigenvalue*, SIAM
   J. Sci. Comput., Vol. 29, No. 2, (2007), 481--514.

.. [r3] A. Stathopoulos and J. R. McCombs, *Nearly optimal preconditioned
   methods for hermitian eigenproblems under limited memory. Part II:
   Seeking many eigenvalues*, SIAM J. Sci. Comput., Vol. 29, No. 5, (2007),
   2162-2188.

.. [r4] J. R. McCombs and A. Stathopoulos, *Iterative Validation of
   Eigensolvers: A Scheme for Improving the Reliability of Hermitian
   Eigenvalue Solvers*, SIAM J. Sci. Comput., Vol. 28, No. 6, (2006),
   2337-2358.

.. [r5] A. Stathopoulos, *Locking issues for finding a large number of eigenvectors
   of hermitian matrices*, Tech Report: WM-CS-2005-03, July, 2005.

License Information
-------------------

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


Contact Information 
-------------------

For reporting bugs or questions about functionality contact `Andreas Stathopoulos`_ by
email, `andreas` at `cs.wm.edu`. See further information in
the webpage http://www.cs.wm.edu/~andreas .


Directory Structure
-------------------

The next directories and files should be available:

* ``COPYING.txt``, LGPL License;
* ``Make_flags``,  flags to be used by makefiles to compile library and tests;
* ``Link_flags``,  flags needed in making and linking the test programs;
* ``PRIMMESRC/``,  directory with source code in the following subdirectories:

   * ``COMMONSRC/``, interface and common functions used by all precision versions;
   * ``DSRC/``,      the source code for the double precision :c:func:`dprimme`;
   * ``ZSRC/``,      the source code for the double complex precision :c:func:`zprimme`;

* ``MEX/``,          MATLAB interface for PRIMME;
* ``TEST/``,         driver and samples in C and F77, both sequential and parallel;
* ``libprimme.a``,   the PRIMME library (to be made);
* ``makefile``       main make file;
* ``readme.txt``     text version of the documentation;
* ``doc/``           directory with the HTML and PDF versions of the documentation.


Making and Linking
------------------

``Make_flags`` has the flags and compilers used to make ``libprimme.a``:

* `CC`, compiler program such as ``gcc``, ``clang`` or ``icc``.
* `CFLAGS`, compiler options such as ``-g`` or ``-O3``. Also include some of these
  options if the BLAS_ and LAPACK_ that will be linked:

  * ``-DF77UNDERSCORE``, the Fortran function names has appended an underscore
    (usually they does).
  * ``-DPRIMME_BLASINT_SIZE=64``, integers are 64-bit integer (``kind=8``) type
    (usually they doesn't).

.. note::

   When ``-DPRIMME_BLASINT_SIZE=64`` is set the code uses the type ``int64_t``
   supported by the starndard C99. In case the compiler doesn't honor the
   standard, replace the next lines in :file:`PRIMMESRC/COMMONSRC/common_numerical.h`::

      #if !defined(PRIMME_BLASINT_SIZE)
      #  define PRIMME_BLASINT int
      #else
      #  include <stdint.h>
      #  define GENERIC_INT(N) int ## N ## _t
      #  define XGENERIC_INT(N) GENERIC_INT(N)
      #  define PRIMME_BLASINT XGENERIC_INT(PRIMME_BLASINT_SIZE)
      #endif

   by the next macro definition with the proper type for an ``int`` of 64 bits::

      #define PRIMME_BLASINT __int64


After customizing ``Make_flags``, type this to generate ``libprimme.a``::

    make lib

Making can be also done at the command line::

    make lib CC=clang CFLAGS='-O3'

``Link_flags`` has the flags for linking with external libraries and making the executables
located in ``TEST``:

* `LDFLAGS`, linker flags such as ``-framework Accelerate``.
* `LIBS`, flags to link with libraries (BLAS_ and LAPACK_ are required), such as ``-lprimme -llapack -lblas -lgfortran -lm``.

After that type this to compile and execute a simple test::

    make test 

If it worked, try with other examples in ``TEST`` (see ``README`` in ``TEST`` for more
information about how to compile the driver and the examples).

Full description of actions that `make` can take:

* `make lib`, builds ``libprimme.a``; alternatively:
* `make libd`, if only :c:func:`dprimme` is of interest, build ``libdprimme.a``:
* `make libz`, if only :c:func:`zprimme` is of interest, build ``libzprimme.a``;
* `make test`, build and execute a simple example; 
* `make clean`, removes all ``*.o``, ``a.out``, and core files from all directories.

Tested Systems
--------------

PRIMME is primary developed with GNU gcc, g++ and gfortran (versions 4.8 and later).
Many users have reported builds on several other platforms/compilers:

* SUSE 13.1 & 13.2
* CentOS 6.6
* Ubuntu 14.04
* MacOS X 10.9 & 10.10 
* Cygwin & MinGW
* Cray XC30
* SunOS 5.9, quad processor Sun-Fire-280R, and several other UltraSparcs
* AIX 5.2 IBM SP POWER 3+, 16-way SMP, 375 MHz nodes (seaborg at nersc.gov)

.. include:: epilog.inc
