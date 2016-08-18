import numpy as np
from scipy.sparse.linalg.interface import aslinearoperator

__docformat__ = "restructuredtext en"

__PRIMMEErrors = {
0: "success",
1: "reported only amount of required memory",
-1: "failed in allocating int or real workspace",
-2: "malloc failed in allocating a permutation integer array",
-3: "main_iter() encountered problem; the calling stack of the functions where the error occurred was printed in 'stderr'",
-4: "if argument 'primme' is NULL",
-5: "if 'n' <= 0 or 'nLocal' <= 0",
-6: "if 'numProcs' < 1",
-7: "if 'matrixMatvec' is NULL",
-8: "if 'applyPreconditioner' is NULL and 'precondition' is not NULL",
-9: "if 'globalSumDouble' is NULL",
-10: "if 'numEvals' > 'n'",
-11: "if 'numEvals' < 0",
-12: "if 'eps' > 0 and 'eps' < machine precision",
-13: "if 'target' is not properly defined",
-14: "if 'target' is one of 'primme_closest_geq', 'primme_closest_leq' or 'primme_closest_abs' but 'numTargetShifts' <= 0 (no shifts)",
-15: "if 'target' is one of 'primme_closest_geq', 'primme_closest_leq' or 'primme_closest_abs' but 'targetShifts' is NULL  (no shifts array)",
-16: "if 'numOrthoConst' < 0 or 'numOrthoConst' >= 'n'. (no free dimensions left)",
-17: "if 'maxBasisSize' < 2",
-18: "if 'minRestartSize' <= 0",
-19: "if 'maxBlockSize' <= 0",
-20: "if 'maxPrevRetain' < 0",
-21: "if 'scheme' is not one of *primme_thick* or *primme_dtr*",
-22: "if 'initSize' < 0",
-23: "if not 'locking' and 'initSize' > 'maxBasisSize'",
-24: "if 'locking' and 'initSize' > 'numEvals'",
-25: "if 'maxPrevRetain' + 'minRestartSize' >= 'maxBasisSize'",
-26: "if 'minRestartSize' >= 'n'",
-27: "if 'printLevel' < 0 or 'printLevel' > 5",
-28: "if 'convTest' is not one of 'primme_full_LTolerance', 'primme_decreasing_LTolerance', 'primme_adaptive_ETolerance' or 'primme_adaptive'",
-29: "if 'convTest' == 'primme_decreasing_LTolerance' and 'relTolBase' <= 1",
-30: "if 'evals' is NULL, but not 'evecs' and 'resNorms'",
-31: "if 'evecs' is NULL, but not 'evals' and 'resNorms'",
-32: "if 'resNorms' is NULL, but not 'evecs' and 'evals'",
-33: "if not 'locking' and 'minRestartSize' < 'numEvals'"
}


class PrimmeError(RuntimeError):
    """
    PRIMME error
    """
    def __init__(self, err):
        self.err = err
        RuntimeError.__init__(self, "PRIMME error %d: %s" % (err, __PRIMMEErrors[err]))

class PrimmeSvdsError(RuntimeError):
    """
    PRIMME SVDS error
    """
    def __init__(self, err):
        self.err = err
        RuntimeError.__init__(self, "PRIMME SVDS error %d: %s" % (err, __PRIMMEErrors[err]))


def eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None,
          ncv=None, maxiter=None, tol=0, return_eigenvectors=True,
          Minv=None, OPinv=None, mode='normal'):
    """
    Find k eigenvalues and eigenvectors of the real symmetric square matrix
    or complex hermitian matrix A.

    Solves ``A * x[i] = w[i] * x[i]``, the standard eigenvalue problem for
    w[i] eigenvalues with corresponding eigenvectors x[i].

    If M is specified, solves ``A * x[i] = w[i] * M * x[i]``, the
    generalized eigenvalue problem for w[i] eigenvalues
    with corresponding eigenvectors x[i]

    Parameters
    ----------
    A : An N x N matrix, array, sparse matrix, or LinearOperator representing
        the operation A * x, where A is a real symmetric matrix
        For buckling mode (see below) A must additionally be positive-definite
    k : int, optional
        The number of eigenvalues and eigenvectors desired.
        `k` must be smaller than N. It is not possible to compute all
        eigenvectors of a matrix.

    Returns
    -------
    w : array
        Array of k eigenvalues
    v : array
        An array representing the `k` eigenvectors.  The column ``v[:, i]`` is
        the eigenvector corresponding to the eigenvalue ``w[i]``.

    Other Parameters
    ----------------
    M : An N x N matrix, array, sparse matrix, or linear operator representing
        the operation M * x for the generalized eigenvalue problem

            A * x = w * M * x.

        M must represent a real, symmetric matrix if A is real, and must
        represent a complex, hermitian matrix if A is complex. For best
        results, the data type of M should be the same as that of A.
        Additionally:

            If sigma is None, M is symmetric positive definite

            If sigma is specified, M is symmetric positive semi-definite

            In buckling mode, M is symmetric indefinite.

        If sigma is None, eigsh requires an operator to compute the solution
        of the linear equation ``M * x = b``. This is done internally via a
        (sparse) LU decomposition for an explicit matrix M, or via an
        iterative solver for a general linear operator.  Alternatively,
        the user can supply the matrix or operator Minv, which gives
        ``x = Minv * b = M^-1 * b``.
    sigma : real
        Find eigenvalues near sigma using shift-invert mode.  This requires
        an operator to compute the solution of the linear system
        `[A - sigma * M] x = b`, where M is the identity matrix if
        unspecified.  This is computed internally via a (sparse) LU
        decomposition for explicit matrices A & M, or via an iterative
        solver if either A or M is a general linear operator.
        Alternatively, the user can supply the matrix or operator OPinv,
        which gives ``x = OPinv * b = [A - sigma * M]^-1 * b``.
        Note that when sigma is specified, the keyword 'which' refers to
        the shifted eigenvalues ``w'[i]`` where:

            if mode == 'normal', ``w'[i] = 1 / (w[i] - sigma)``.

            if mode == 'cayley', ``w'[i] = (w[i] + sigma) / (w[i] - sigma)``.

            if mode == 'buckling', ``w'[i] = w[i] / (w[i] - sigma)``.

        (see further discussion in 'mode' below)
    v0 : ndarray, optional
        Starting vector for iteration.
        Default: random
    ncv : int, optional
        The number of Lanczos vectors generated ncv must be greater than k and
        smaller than n; it is recommended that ``ncv > 2*k``.
        Default: ``min(n, 2*k + 1)``
    which : str ['LM' | 'SM' | 'LA' | 'SA' | 'BE']
        If A is a complex hermitian matrix, 'BE' is invalid.
        Which `k` eigenvectors and eigenvalues to find:

            'LM' : Largest (in magnitude) eigenvalues

            'SM' : Smallest (in magnitude) eigenvalues

            'LA' : Largest (algebraic) eigenvalues

            'SA' : Smallest (algebraic) eigenvalues

            'BE' : Half (k/2) from each end of the spectrum (not supported)

        When sigma != None, 'which' refers to the shifted eigenvalues ``w'[i]``
        (see discussion in 'sigma', above).
    maxiter : int, optional
        Maximum number of restarts update iterations allowed
        Default: ``n*10``
    tol : float
        Accuracy for eigenvalues (stopping criterion).
        The default value is sqrt of machine precision.
    Minv : N x N matrix, array, sparse matrix, or LinearOperator
        See notes in M, above
    OPinv : N x N matrix, array, sparse matrix, or LinearOperator
        See notes in sigma, above.
    return_eigenvectors : bool
        Return eigenvectors (True) in addition to eigenvalues
    mode : string ['normal' | 'buckling' | 'cayley']
        Only 'normal' mode is supported.

    Raises
    ------
    PrimmeError
        When the requested convergence is not obtained.

        The PRIMME error code can be found as ``err`` attribute of the exception
        object.

    See Also
    --------
    eigs : eigenvalues and eigenvectors for a general (nonsymmetric) matrix A
    svds : singular value decomposition for a matrix A

    Notes
    -----
    This function is a wrapper to PRIMME functions to find the eigenvalues and
    eigenvectors [2]_.

    References
    ----------
    .. [1] ARPACK Software, http://www.caam.rice.edu/software/ARPACK/
    .. [2] R. B. Lehoucq, D. C. Sorensen, and C. Yang,  ARPACK USERS GUIDE:
       Solution of Large Scale Eigenvalue Problems by Implicitly Restarted
       Arnoldi Methods. SIAM, Philadelphia, PA, 1998.

    Examples
    --------
    >>> import Primme
    >>> import numpy as np
    >>> from scipy.sparse import spdiags
    >>> a = np.ones(3)
    >>> A  = spdiags(np.array([a*(-1.), a*2., a*(-1.)]), np.array([-1, 0, 1]), 10, 10)
    >>> r = Primme.eigsh(A, which='LA')
    >>> r['w'] % values 
    """

    A = aslinearoperator(A)
    if A.shape[0] != A.shape[1]:
        raise ValueError('expected square matrix (shape=%s)' % (A.shape,))

    if M is not None:
        raise ValueError('generalized problems (M != None) are not supported')

    if OPinv != None:
        OPinv = aslinearoperator(OPinv)
        if OPinv.shape[0] != OPinv.shape[1] or OPinv.shape[0] != A.shape[0]:
            raise ValueError('expected square matrix with same shape as A (shape=%s)' % (OPinv.shape,))

    class PP(PrimmeParams):
        def __init__(self):
            PrimmeParams.__init__(self)
        def matvec(self, X):
            return A.matmat(X)
        def prevec(self, X):
            return OPinv.matmat(X)

    pp = PP()
 
    pp.n = A.shape[0]

    if k <= 0 or k >= pp.n:
        raise ValueError("k must be between 1 and the order of the "
                         "square input matrix.")
    pp.numEvals = k
    pp.correctionParams.precondition = 1 if OPinv != None else 0

    if which == 'LA' and sigma == None:
        pp.target = primme_largest
    elif which == 'SA' and sigma == None:
        pp.target = primme_smallest
    elif which == 'SM':
        pp.target = primme_closest_abs
        pp.numTargetShifts = 1
        if sigma != None:
            sigma = 0.0
        pp.targetShifts = np.array([sigma], dtype=np.dtype('d'))
    else:
        raise ValueError("which value '%s' and sigma value '%s' not supported" % (which, sigma))

    pp.eps = tol

    if ncv != None:
        pp.maxBasisSize = ncv

    if maxiter != None:
        pp.maxMatvecs = maxiter

    pp.set_method(DYNAMIC)

    evals = np.zeros(pp.numEvals)
    norms = np.zeros(pp.numEvals)
    evecs = np.zeros((pp.n, pp.numEvals), A.dtype)

    if v0 != None:
        pp.initSize = v0.shape[1]
        np.copyto(evecs[:, 0:pp.initSize], v0)

    if A.dtype is np.dtype(np.complex128):
        err = zprimme(evals, evecs, norms, pp)
    elif A.dtype is np.dtype('d'):
        err = dprimme(evals, evecs, norms, pp)
    else:
        raise ValueError("dtype of A not supported")

    if err != 0:
        raise PrimmeError(err)

    return { 'w': evals, 'v': evecs }


def svds(A, k=6, ncv=None, tol=0, which='LM', v0=None,
         maxiter=None, return_singular_vectors=True):
    """Compute the largest k singular values/vectors for a sparse matrix.
    Parameters
    ----------
    A : {sparse matrix, LinearOperator}
        Array to compute the SVD on, of shape (M, N)
    k : int, optional
        Number of singular values and vectors to compute.
        Must be 1 <= k < min(A.shape).
    ncv : int, optional
        The maximum size of the basis
    tol : float, optional
        Tolerance for singular values. Zero (default) means machine precision.
    which : str, ['LM' | 'SM'], optional
        Which `k` singular values to find:
            - 'LM' : largest singular values
            - 'SM' : smallest singular values
    v0 : ndarray, optional
        Starting vectors for iteration, of length min(A.shape). Should be
        (approximate) left singular vectors if N > M and a right singular
        vectors otherwise.
    maxiter : int, optional
        Maximum number of iterations.
    Returns
    -------
    u : ndarray, shape=(M, k)
        Unitary matrix having left singular vectors as columns.
        If `return_singular_vectors` is "vh", this variable is not computed,
        and None is returned instead.
    s : ndarray, shape=(k,)
        The singular values.
    vt : ndarray, shape=(k, N)
        Unitary matrix having right singular vectors as rows.
        If `return_singular_vectors` is "u", this variable is not computed,
        and None is returned instead.
    """

    A = aslinearoperator(A)

    n, m = A.shape

    if k <= 0 or k >= min(n, m):
        raise ValueError("k must be between 1 and min(A.shape), k=%d" % k)

    class PSP(PrimmeSvdsParams):
        def __init__(self):
            PrimmeSvdsParams.__init__(self)

        def matvec(self, X, transpose):
            if transpose == 0:
                return A.matmat(X)
            else:
                return A.H.matmat(X) 

    pp = PSP()

    pp.m = A.shape[0]
    pp.n = A.shape[1]

    pp.numSvals = k

    if which == 'LM':
        pp.target = primme_svds_largest
    elif which == 'SM':
        pp.target = primme_svds_smallest
    else:
        raise ValueError("which must be either 'LM' or 'SM'.")

    pp.eps = tol

    if v0 != None:
        pp.initSize = v0.shape[1]

    if ncv:
        pp.maxBasisSize = ncv

    if maxiter != None:
        pp.maxMatvecs = maxiter

    svals = np.zeros(pp.numSvals)
    svecsl = np.zeros((pp.m, pp.numSvals), A.dtype)
    svecsr = np.zeros((pp.n, pp.numSvals), A.dtype)
    norms = np.zeros(pp.numSvals)

    if v0 != None:
        pp.initSize = v0.shape[1]
        np.copyto(evecs[:, 0:pp.initSize], v0)

    if A.dtype is np.dtype('d'):
        err = dprimme_svds(svals, svecsl, svecsr, norms, pp)
    elif A.dtype is np.dtype(np.complex128):
        err = zprimme_svds(svals, svecsl, svecsr, norms, pp)
    else:
        raise ValueError("dtype of A not supported")

    if err != 0:
        raise PrimmeSvdsError(err)

    if not return_singular_vectors:
        return svals

    # Transpose conjugate svecsr
    svecsr = svecsr.T.conj()

    return svecsl, svals, svecsr

# if __name__ == '__main__':
#     from scipy.sparse import spdiags
#     a = np.ones(10)
#     A  = spdiags(np.array([a*(-1.), a*2., a*(-1.)]), np.array([-1, 0, 1]), 10, 10)
#     print(eigsh_primme(A, which='LA'))
