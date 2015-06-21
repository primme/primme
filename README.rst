
PRIMME: PReconditioned Iterative MultiMethod Eigensolver
========================================================

PRIMME, pronounced as *prime*, finds a number of eigenvalues and their corresponding eigenvectors of a 
real symmetric, or complex hermitian matrix A. Largest, smallest and interior 
eigenvalues are supported. Preconditioning can be used to accelerate 
convergence. 
PRIMME is written in C, but a complete Fortran 77 interface is also provided.
  
Making & Linking
----------------

`Make_flags` has the flags and compilers used to make `libprimme.a`. Set at minimum:

* `TOP`, path where the PRIMME directory is located.
* `CC`, compiler program like `gcc`, `clang` or `icc`. 

Then do this to generate `libprimme.a`::

    make lib

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
* `rNorms`, array to return the residual of the found eigenpairs; and
* `primme`, structure that specify the matrix problem, which eigenvalues are wanted and several method options.

See documentation in `readme.txt` and `doc.pdf`.

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

Contact Information 
-------------------

For reporting bugs or questions about functionality contact `Andreas Stathopoulos`_

.. _`Andreas Stathopoulos`: http://www.cs.wm.edu/~andreas/
