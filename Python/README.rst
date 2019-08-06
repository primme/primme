
PRIMME: PReconditioned Iterative MultiMethod Eigensolver
========================================================

`primme` is a Python interface to PRIMME_, a high-performance library for computing a few eigenvalues/eigenvectors, and singular values/vectors.
PRIMME is especially optimized for large, difficult problems.
Real symmetric and complex Hermitian problems, standard `A x = \lambda x` and generalized `A x = \lambda B x`, are supported.
It can find largest, smallest, or interior singular/eigenvalues, and can use preconditioning to accelerate convergence.

The main contributors to PRIMME are James R. McCombs, Eloy Romero Alcalde, Andreas Stathopoulos and Lingfei Wu.

Install
-------

You can install the latest version with `pip`::

    pip install numpy   # if numpy is not installed yet
    pip install scipy   # if scipy is not installed yet
    pip install future  # if using python 2
    conda install mkl-devel # if using Anaconda Python distribution
    pip install primme

Optionally for building the development version do::

    git clone https://github.com/primme/primme
    cd primme
    make python_install

Usage
-----

The following examples compute a few eigenvalues and eigenvectors from a real symmetric matrix::

    >>> import Primme, scipy.sparse
    >>> A = scipy.sparse.spdiags(range(100), [0], 100, 100) # sparse diag. matrix
    >>> evals, evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA')
    >>> evals # the three largest eigenvalues of A
    array([ 99.,  98.,  97.])

    >>> new_evals, new_evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA', ortho=evecs)
    >>> new_evals # the next three largest eigenvalues
    array([ 96.,  95.,  94.])

    >>> evals, evecs = primme.eigsh(A, 3, tol=1e-6, which=50.1)
    >>> evals # the three closest eigenvalues to 50.1
    array([ 50.,  51.,  49.])


The following examples compute a few eigenvalues and eigenvectors from a generalized Hermitian problem, without factorizing or inverting `B`::

    >>> import Primme, scipy.sparse
    >>> A = scipy.sparse.spdiags(range(100), [0], 100, 100) # sparse diag. matrix
    >>> M = scipy.sparse.spdiags(np.asarray(range(99,-1,-1)), [0], 100, 100)
    >>> evals, evecs = primme.eigsh(A, 3, M=M, tol=1e-6, which='SA')
    >>> evals
    array([1.0035e-07, 1.0204e-02, 2.0618e-02])

The following examples compute a few singular values and vectors::

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

Further examples_.

Documentation of eigsh_ and svds_.

Citing this code 
----------------

Please cite (bibtex_):

* A. Stathopoulos and J. R. McCombs *PRIMME: PReconditioned Iterative
  MultiMethod Eigensolver: Methods and software description*, ACM
  Transaction on Mathematical Software Vol. 37, No. 2, (2010),
  21:1-21:30.

* L. Wu, E. Romero and A. Stathopoulos, *PRIMME_SVDS: A High-Performance
  Preconditioned SVD Solver for Accurate Large-Scale Computations*,
  J. Sci. Comput., Vol. 39, No. 5, (2017), S248--S271.

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
.. _PETSc: http://www.mcs.anl.gov/petsc/
.. _`bibtex`: https://raw.githubusercontent.com/primme/primme/master/doc/primme.bib
.. _eigsh: http://www.cs.wm.edu/~andreas/software/doc/pyeigsh.html
.. _svds: http://www.cs.wm.edu/~andreas/software/doc/pysvds.html
.. _examples: https://github.com/primme/primme/blob/master/Python/examples.py
