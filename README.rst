
PRIMME: PReconditioned Iterative MultiMethod Eigensolver
========================================================

PRIMME, pronounced as *prime*, finds a number of eigenvalues and their corresponding eigenvectors of a 
real symmetric, or complex hermitian matrix A. Largest, smallest and interior 
eigenvalues are supported. Preconditioning can be used to accelerate 
convergence. 
PRIMME is written in C, but complete interfaces are provided for Fortran 77 and MATLAB.
  
Making and Linking
------------------

`Make_flags` has the flags and compilers used to make `libprimme.a`:

* `CC`, compiler program such as ``gcc``, ``clang`` or ``icc``.
* `CFLAGS`, compiler options such as ``-g`` or ``-O3``.

After customizing `Make_flags`, type this to generate `libprimme.a`::

    make lib

Making can be also done at the command line::

    make lib CC=clang CFLAGS='-O3'


C Library Interface
-------------------

To solve real symmetric standard eigenproblems call::

    int dprimme(double *evals, double *evecs, double *resNorms, 
                primme_params *primme);

To solve Hermitian standard eigenproblems call::

    int zprimme(double *evals, Complex_Z *evecs, double *resNorms, 
                primme_params *primme);

The call arguments are:

* `evals`, array to return the found eigenvalues;
* `evecs`, array to return the found eigenvectors;
* `resNorms`, array to return the residual norms of the found eigenpairs; and
* `primme`, structure that specify the matrix problem, which eigenvalues are wanted and several method options.

See documentation in `readme.txt` file and in ``doc`` directory.

Citing this code 
----------------

Please cite:

* A. Stathopoulos and J. R. McCombs PRIMME: *PReconditioned Iterative
  MultiMethod Eigensolver: Methods and software description*, ACM
  Transaction on Mathematical Software Vol. 37, No. 2, (2010),
  21:1-21:30.

More information on the algorithms and research that led to this
software can be found in the rest of the papers. The work has been
supported by a number of grants from the National Science Foundation.

* A. Stathopoulos, *Nearly optimal preconditioned methods for hermitian
  eigenproblems under limited memory. Part I: Seeking one eigenvalue*, SIAM
  J. Sci. Comput., Vol. 29, No. 2, (2007), 481--514.

* A. Stathopoulos and J. R. McCombs, *Nearly optimal preconditioned
  methods for hermitian eigenproblems under limited memory. Part II:
  Seeking many eigenvalues*, SIAM J. Sci. Comput., Vol. 29, No. 5, (2007),
  2162-2188.

* J. R. McCombs and A. Stathopoulos, *Iterative Validation of
  Eigensolvers: A Scheme for Improving the Reliability of Hermitian
  Eigenvalue Solvers*, SIAM J. Sci. Comput., Vol. 28, No. 6, (2006),
  2337-2358.

* A. Stathopoulos, *Locking issues for finding a large number of eigenvectors
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
the webpage http://www.cs.wm.edu/~andreas/software .

.. _`Andreas Stathopoulos`: http://www.cs.wm.edu/~andreas/software
