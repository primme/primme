
PRIMME: PReconditioned Iterative MultiMethod Eigensolver
========================================================

PRIMME, pronounced as *prime*, computes
a few eigenvalues and their corresponding eigenvectors of a real symmetric or complex Hermitian matrix. 
It can also compute singular values and vectors of a square or rectangular matrix. 
It can find largest, smallest, or interior singular/eigenvalues and can use preconditioning to accelerate convergence. 
It is especially optimized for large, difficult problems, and can be a useful tool for both non-experts and experts.
PRIMME is written in C99, but complete interfaces are provided for Fortran 77, MATLAB, Python, and R.

Making and Linking
------------------

`Make_flags` has the flags and compilers used to make `libprimme.a`:

* `CC`, compiler program such as ``gcc``, ``clang`` or ``icc``.
* `CFLAGS`, compiler options such as ``-g`` or ``-O3``.

After customizing `Make_flags`, type this to generate `libprimme.a`::

    make lib

Making can be also done at the command line::

    make lib CC=clang CFLAGS='-O3'

Optionally for building some of the external interfaces just do::

    make matlab
    make octave
    make python
    make R_install

Alternatively to install the development version of PRIMME on R::

    library(devtools)
    install_github("primme/primme", subdir="R")

C Library Interface
-------------------

To compute few eigenvalues and eigenvectors from a real symmetric matrix call::

    int dprimme(double *evals, double *evecs, double *resNorms, 
                primme_params *primme);

The call arguments are:

* `evals`, array to return the found eigenvalues;
* `evecs`, array to return the found eigenvectors;
* `resNorms`, array to return the residual norms of the found eigenpairs; and
* `primme`, structure that specify the matrix problem, which eigenvalues are wanted and several method options.

To compute few singular values and vectors from a matrix call::

    int dprimme_svds(double *svals, double *svecs, double *resNorms, 
                primme_svds_params *primme_svds);

The call arguments are:

* `svals`, array to return the found singular values;
* `svecs`, array to return the found vectors;
* `resNorms`, array to return the residual norms of the triplets; and
* `primme_svds`, structure that specify the matrix problem, which values are wanted and several method options.

There are available versions for complex double, float and float double.
See documentation in `readme.txt` file and in ``doc`` directory; also it is online at doc_.
The `examples` directory is plenty of self-contained examples in C, C++ and F77, and some of them using PETSc_.

Citing this code 
----------------

Please cite (bibtex_):

* A. Stathopoulos and J. R. McCombs *PRIMME: PReconditioned Iterative
  MultiMethod Eigensolver: Methods and software description*, ACM
  Transaction on Mathematical Software Vol. 37, No. 2, (2010),
  21:1-21:30.

* L. Wu, E. Romero and A. Stathopoulos, *PRIMME_SVDS: A High-Performance
  Preconditioned SVD Solver for Accurate Large-Scale Computations*,
  arXiv:1607.01404

More information on the algorithms and research that led to this
software can be found in the rest of the papers. The work has been
supported by a number of grants from the National Science Foundation.

* A. Stathopoulos, *Nearly optimal preconditioned methods for Hermitian
  eigenproblems under limited memory. Part I: Seeking one eigenvalue*, SIAM
  J. Sci. Comput., Vol. 29, No. 2, (2007), 481--514.

* A. Stathopoulos and J. R. McCombs, *Nearly optimal preconditioned
  methods for Hermitian eigenproblems under limited memory. Part II:
  Seeking many eigenvalues*, SIAM J. Sci. Comput., Vol. 29, No. 5, (2007),
  2162-2188.

* J. R. McCombs and A. Stathopoulos, *Iterative Validation of
  Eigensolvers: A Scheme for Improving the Reliability of Hermitian
  Eigenvalue Solvers*, SIAM J. Sci. Comput., Vol. 28, No. 6, (2006),
  2337-2358.

* A. Stathopoulos, *Locking issues for finding a large number of eigenvectors
  of Hermitian matrices*, Tech Report: WM-CS-2005-03, July, 2005.

* L. Wu and A. Stathopoulos, *A Preconditioned Hybrid SVD Method for Computing
  Accurately Singular Triplets of Large Matrices*, SIAM J. Sci. Comput. 37-5(2015),
  pp. S365-S388.

License Information
-------------------

PRIMME is licensed under the 3-clause license BSD.
Python and Matlab interfaces have BSD-compatible licenses.
Source code under `tests` is compatible with LGPLv3.
Details can be taken from COPYING.txt.

Contact Information 
-------------------

For reporting bugs or questions about functionality contact `Andreas Stathopoulos`_ by
email, `andreas` at `cs.wm.edu`. See further information in
the webpage http://www.cs.wm.edu/~andreas/software.

.. _`Andreas Stathopoulos`: http://www.cs.wm.edu/~andreas/software
.. _`github`: https://github.com/primme/primme
.. _`doc`: http://www.cs.wm.edu/~andreas/software/doc/readme.html
.. _PETSc : http://www.mcs.anl.gov/petsc/
.. _`bibtex`: https://raw.githubusercontent.com/primme/primme/master/doc/primme.bib
