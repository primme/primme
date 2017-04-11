
PRIMME: PReconditioned Iterative MultiMethod Eigensolver
========================================================

`Primme` is a Python interface to PRIMME_, a C library for computing a few
eigenvalues and their corresponding eigenvectors of a real symmetric or complex
Hermitian matrix. It can also compute singular values and vectors of a square
or rectangular matrix. It can find largest, smallest, or interior
singular/eigenvalues and can use preconditioning to accelerate convergence. It
is especially optimized for large, difficult problems, and can be a useful tool
for both non-experts and experts.

Install
-------

You can install the latest version with `pip`::

    pip install numpy   # if numpy is not installed yet
    pip install scipy   # if scipy is not installed yet
    pip install primme

Optionally for building the development version do::

    git clone https://github.com/primme/primme
    cd primme
    make python_install

Usage
-----

In the following examples it is computed few eigenvalues and eigenvectors from a real symmetric matrix::

    >>> import Primme, scipy.sparse
    >>> A = scipy.sparse.spdiags(range(100), [0], 100, 100) # sparse diag. matrix
    >>> evals, evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA')
    >>> evals # the three largest eigenvalues of A
    array([ 99.,  98.,  97.])

    >>> new_evals, new_evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA', ortho=evecs)
    >>> new_evals # the next three largest eigenvalues
    array([ 96.,  95.,  94.])

In the following examples it is computed few singular values and vectors::

    >>> import Primme, scipy.sparse
    >>> A = scipy.sparse.spdiags(range(1, 11), [0], 100, 10) # sparse diag. rect. matrix
    >>> svecs_left, svals, svecs_right = Primme.svds(A, 3, tol=1e-6, which='SM')
    >>> svals # the three smallest singular values of A
    array([ 1.,  2.,  3.])

    >>> A = scipy.sparse.rand(10000, 100, random_state=10)
    >>> prec = scipy.sparse.spdiags(np.reciprocal(A.multiply(A).sum(axis=0)),
    ...           [0], 100, 100) # square diag. preconditioner
    >>> svecs_left, svals, svecs_right = Primme.svds(A, 3, which=6.0, tol=1e-6,
    ...           precAHA=prec)
    >>> ["%.5f" % x for x in svals.flat] # the three closest singular values of A to 0.5
    ['5.99871', '5.99057', '6.01065']

Check further examples_ and the documentation of eigsh_ and svds_.

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

License Information
-------------------

PRIMME and this interface is licensed under the 3-clause license BSD.

Contact Information 
-------------------

For reporting bugs or questions about functionality contact `Andreas Stathopoulos`_ by
email, `andreas` at `cs.wm.edu`. See further information in
the webpage http://www.cs.wm.edu/~andreas/software.

.. _PRIMME: https://github.com/primme/primme
.. _`Andreas Stathopoulos`: http://www.cs.wm.edu/~andreas/software
.. _`github`: https://github.com/primme/primme
.. _`doc`: http://www.cs.wm.edu/~andreas/software/doc/readme.html
.. _PETSc : http://www.mcs.anl.gov/petsc/
.. _`bibtex`: https://raw.githubusercontent.com/primme/primme/master/doc/primme.bib
.. _eigsh: http://www.cs.wm.edu/~andreas/software/doc/pyeigsh.html
.. _svds: http://www.cs.wm.edu/~andreas/software/doc/pysvds.html
.. _examples: https://github.com/primme/primme/blob/master/Python/examples.py
