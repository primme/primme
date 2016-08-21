import numpy as np
from scipy.sparse.linalg.interface import aslinearoperator

__docformat__ = "restructuredtext en"

_PRIMMEErrors = {
0: "success",
1: "reported only amount of required memory",
-1: "failed in allocating int or real workspace",
-2: "malloc failed in allocating a permutation integer array",
-3: "main_iter() encountered problem; the calling stack of the functions where the error occurred was printed in 'stderr'",
-4: "argument 'primme' is NULL",
-5: "'n' <= 0 or 'nLocal' <= 0",
-6: "'numProcs' < 1",
-7: "'matrixMatvec' is NULL",
-8: "'applyPreconditioner' is NULL and 'precondition' is not NULL",
-9: "'globalSumDouble' is NULL",
-10: "'numEvals' > 'n'",
-11: "'numEvals' < 0",
-12: "'eps' > 0 and 'eps' < machine precision",
-13: "'target' is not properly defined",
-14: "'target' is one of 'primme_largest_abs', 'primme_closest_geq', 'primme_closest_leq' or 'primme_closest_abs' but 'numTargetShifts' <= 0 (no shifts)",
-15: "'target' is one of 'primme_largest_abs', 'primme_closest_geq', 'primme_closest_leq' or 'primme_closest_abs' but 'targetShifts' is NULL  (no shifts array)",
-16: "'numOrthoConst' < 0 or 'numOrthoConst' >= 'n'. (no free dimensions left)",
-17: "'maxBasisSize' < 2",
-18: "'minRestartSize' <= 0",
-19: "'maxBlockSize' <= 0",
-20: "'maxPrevRetain' < 0",
-21: "'scheme' is not one of *primme_thick* or *primme_dtr*",
-22: "'initSize' < 0",
-23: "not 'locking' and 'initSize' > 'maxBasisSize'",
-24: "'locking' and 'initSize' > 'numEvals'",
-25: "'maxPrevRetain' + 'minRestartSize' >= 'maxBasisSize'",
-26: "'minRestartSize' >= 'n'",
-27: "'printLevel' < 0 or 'printLevel' > 5",
-28: "'convTest' is not one of 'primme_full_LTolerance', 'primme_decreasing_LTolerance', 'primme_adaptive_ETolerance' or 'primme_adaptive'",
-29: "'convTest' == 'primme_decreasing_LTolerance' and 'relTolBase' <= 1",
-30: "'evals' is NULL, but not 'evecs' and 'resNorms'",
-31: "'evecs' is NULL, but not 'evals' and 'resNorms'",
-32: "'resNorms' is NULL, but not 'evecs' and 'evals'",
-33: "not 'locking' and 'minRestartSize' < 'numEvals'"
}

_PRIMMESvdsErrors = {
0   : "success",
1   : "reported only amount of required memory",
-1  : "failed in allocating int or real workspace",
-2  : "malloc failed in allocating a permutation integer array",
-3  : "main_iter() encountered problem; the calling stack of the functions where the error occurred was printed in 'stderr'",
-4  : "primme_svds is NULL",
-5  : "Wrong value for m or n",
-6  : "Wrong value for numProcs",
-7  : "matrixMatvec is not set",
-8  : "applyPreconditioner is not set but precondition == 1 ",
-9  : "numProcs >1 but globalSumDouble is not set",
-10 : "Wrong value for numSvals, it's larger than min(m, n)",
-11 : "Wrong value for numSvals, it's smaller than 1",
-13 : "Wrong value for target",
-14 : "Wrong value for method",
-15 : "Not supported combination of method and methodStage2",
-16 : "Wrong value for printLevel",
-17 : "svals is not set",
-18 : "svecs is not set",
-19 : "resNorms is not set"
}


class PrimmeError(RuntimeError):
    """
    PRIMME error
    """
    def __init__(self, err):
        self.err = err
        RuntimeError.__init__(self, "PRIMME error %d: %s" % (err, _PRIMMEErrors[err]))

class PrimmeSvdsError(RuntimeError):
    """
    PRIMME SVDS error
    """
    def __init__(self, err):
        self.err = err
        if err < 100:
            msg = _PRIMMESvdsErrors[err]
        elif err < 200:
            msg = "Error from PRIMME first stage: " + _PRIMMEErrors[err+100]
        elif err < 300:
            msg = "Error from PRIMME second stage: " + _PRIMMEErrors[err+200]
        RuntimeError.__init__(self, "PRIMME SVDS error %d: %s" % (err, msg))


def eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None,
          ncv=None, maxiter=None, tol=0, return_eigenvectors=True,
          Minv=None, OPinv=None, mode='normal', lock=None):
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
    sigma : real, optional
        Find eigenvalues near sigma.
    v0 : N x i, ndarray, optional
        Starting vectors for iteration.
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
    maxiter : int, optional
        Maximum number of restarts update iterations allowed
        Default: ``n*10``
    tol : float
        Accuracy for eigenvalues (stopping criterion).
        The default value is sqrt of machine precision.
    Minv : (not supported)
    OPinv : N x N matrix, array, sparse matrix, or LinearOperator
        Preconditioner to accelerate the convergence. Usually it is an
        approximation of the inverse of (A - sigma*M).
    return_eigenvectors : bool
        Return eigenvectors (True) in addition to eigenvalues
    mode : string ['normal' | 'buckling' | 'cayley']
        Only 'normal' mode is supported.
    lock : N x i, ndarray, optional
        Seek the eigenvectors orthogonal to these ones. The provided
        vectors *should* be orthonormal. Useful to not converge some already
        computed solutions.

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
    eigenvectors [1]_.

    References
    ----------
    .. [1] PRIMME Software, https://github.com/primme/primme

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
        raise ValueError('A: expected square matrix (shape=%s)' % (A.shape,))

    if M is not None:
        raise ValueError('generalized problems (M != None) are not supported')

    if OPinv is not None:
        OPinv = aslinearoperator(OPinv)
        if OPinv.shape[0] != OPinv.shape[1] or OPinv.shape[0] != A.shape[0]:
            raise ValueError('OPinv: expected square matrix with same shape as A (shape=%s)' % (OPinv.shape,))

    class PP(PrimmeParams):
        def __init__(self):
            PrimmeParams.__init__(self)
        def matvec(self, X):
            return A.matmat(X)
        def prevec(self, X):
            return OPinv.matmat(X)

    pp = PP()
 
    pp.n = A.shape[0]

    if k <= 0 or k > pp.n:
        raise ValueError("k=%d must be between 1 and %d, the order of the "
                         "square input matrix." % (k, pp.n))
    pp.numEvals = k
    pp.correctionParams.precondition = 0 if OPinv is None else 1

    if which == 'LM':
        pp.target = primme_largest_abs
        if sigma is None:
            sigma = 0.0
    elif which == 'LA':
        pp.target = primme_largest
        sigma = None
    elif which == 'SA':
        pp.target = primme_smallest
        sigma = None
    elif which == 'SM':
        pp.target = primme_closest_abs
        if sigma is None:
            sigma = 0.0
    else:
        raise ValueError("which='%s' not supported" % which)

    if sigma is not None:
        pp.targetShifts = np.array([sigma], dtype=np.dtype('d'))

    pp.eps = tol

    if ncv is not None:
        pp.maxBasisSize = ncv

    if maxiter is not None:
        pp.maxMatvecs = maxiter

    if OPinv is not None:
        pp.precondition = 1

    if lock is not None:
        if lock.shape[0] != n:
            raise ValueError('lock: expected matrix with the same columns as A (shape=%s)' % (lock.shape,))
        pp.numOrthoConst = min(v0.shape[1], n)

    evals = np.zeros(pp.numEvals)
    norms = np.zeros(pp.numEvals)
    evecs = np.zeros((pp.n, pp.numOrthoConst+pp.numEvals), A.dtype, order='F')

    if lock is not None:
        np.copyto(evecs[:, 0:pp.numOrthoConst], lock[:, 0:pp.numOrthoConst])

    if v0 is not None:
        pp.initSize = min(v0.shape[1], pp.numEvals)
        np.copyto(evecs[:, pp.numOrthoConst:pp.numOrthoConst+pp.initSize],
            v0[:, 0:pp.initSize])

    if A.dtype is np.dtype(np.complex128):
        err = zprimme(evals, evecs, norms, pp)
    elif A.dtype is np.dtype('d'):
        err = dprimme(evals, evecs, norms, pp)
    else:
        raise ValueError("dtype of A not supported")

    if err != 0:
        raise PrimmeError(err)

    evecs = evecs[:, pp.numOrthoConst:]
    return evals, evecs


def svds(A, k=6, ncv=None, tol=0, which='LM', v0=None,
         maxiter=None, return_singular_vectors=True,
         precAHA=None, precAAH=None, precAug=None,
         u0=None, locku0=None, lockv0=None):
    """Compute k singular values/vectors for a sparse matrix.
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
    which : str ['LM' | 'SM'] or number, optional
        Which `k` singular values to find:
            - 'LM' : largest singular values
            - 'SM' : smallest singular values
            - number : closest singular values to (referred as sigma later)
    u0, v0 : ndarray, optional
        Starting vectors for the iterations. Should be approximate left singular
        vectors and right singular vectors respectively. If only u0 or v0 is
        provided, the other is computed.
    maxiter : int, optional
        Maximum number of iterations.
    precAHA : {N x N matrix, array, sparse matrix, LinearOperator}, optional
        Approximate inverse of (A.H*A - sigma*I). If provided and M>N, it
        usually accelerates the convergence.
    precAAH : {M x M matrix, array, sparse matrix, LinearOperator}, optional
        Approximate inverse of (A*A.H - sigma*I). If provided and M<N, it
        usually accelerates the convergence.
    precAug : {(M+N) x (M+N) matrix, array, sparse matrix, LinearOperator}, optional
        Approximate inverse of ([zeros() A.H; zeros() A] - sigma*I). It usually
        accelerates the convergence if tol<dtype.eps**.5.
    locku0, lockv0 : ndarray, optional
        Seek singular triplets orthogonal to these ones. The provided vectors
        *should* be orthonormal. If only locku0 or lockv0 is provided, the other
        is computed. Useful to not converge some already computed solutions.

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

    m, n = A.shape

    if k <= 0 or k > min(n, m):
        raise ValueError("k=%d must be between 1 and min(A.shape)=%d" % (k, min(n, m)))

    if precAHA is not None:
        precAHA = aslinearoperator(precAHA)
        if precAHA.shape[0] != precAHA.shape[1] or precAHA.shape[0] != n:
            raise ValueError('precAHA: expected square matrix with size %d' % n)

    if precAAH is not None:
        precAAH = aslinearoperator(precAAH)
        if precAAH.shape[0] != precAAH.shape[1] or precAAH.shape[0] != m:
            raise ValueError('precAAH: expected square matrix with size %d' % m)

    if precAug is not None:
        precAug = aslinearoperator(precAug)
        if precAug.shape[0] != precAug.shape[1] or precAug.shape[0] != m+n:
            raise ValueError('precAug: expected square matrix with size %d' % (m+n))

    class PSP(PrimmeSvdsParams):
        def __init__(self):
            PrimmeSvdsParams.__init__(self)

        def matvec(self, X, transpose):
            if transpose == 0:
                return A.matmat(X)
            else:
                return A.H.matmat(X) 

        def prevec(self, X, mode):
            if mode == primme_svds_op_AtA and precAHA is not None:
                return precAHA.matmat(X)
            elif mode == primme_svds_op_AAt and precAAH is not None:
                return precAAH.matmat(X) 
            elif mode == primme_svds_op_augmented and precAug is not None:
                return precAug.matmat(X) 
            else:
                raise ValueError('Not expected mode')
            return X

    pp = PSP()

    pp.m = A.shape[0]
    pp.n = A.shape[1]

    pp.numSvals = k

    if which == 'LM':
        pp.target = primme_svds_largest
    elif which == 'SM':
        pp.target = primme_svds_smallest
    else:
        try:
            which = float(which)
        except:
            raise ValueError("which must be either 'LM', 'SM' or a number.")
        pp.target = primme_svds_closest_abs
        pp.targetShifts = np.array([which], dtype='d')

    pp.eps = tol

    if ncv:
        pp.maxBasisSize = ncv

    if maxiter:
        pp.maxMatvecs = maxiter

    def check_pair(u, v, var_names):
        if ((u is not None and u.shape[0] != m) or
                (v is not None and v.shape[0] != n)):
            aux = v; v = u; u = aux

        if ((u is not None and u.shape[0] != m) or
                (v is not None and v.shape[0] != n)):
            aux = v; v = u; u = aux
            raise ValueError("%s don't have the expected number of rows." % var_names)

        if u is not None and v is not None and u.shape[1] != v.shape[1]:
            raise ValueError("%s don't have the same number of columns." % var_names)

        if u is not None and v is None:
            v, _ = np.linalg.qr(A.H.matmult(u))

        if v is not None and u is None:
            u, _ = np.linalg.qr(A.matmult(v))

        return u, v

    locku0, lockv0 = check_pair(locku0, lockv0, "lockv0 or locku0")

    if locku0 is not None:
        pp.numOrthoConst = min(locku0.shape[1], min(m,n))

    svals = np.zeros(pp.numSvals)
    svecsl = np.zeros((pp.m, pp.numOrthoConst+pp.numSvals), A.dtype, order='F')
    svecsr = np.zeros((pp.n, pp.numOrthoConst+pp.numSvals), A.dtype, order='F')
    norms = np.zeros(pp.numSvals)

    if locku0 is not None:
        np.copyto(svecsl[:, 0:pp.numOrthoConst], locku0[:, 0:pp.numOrthoConst])
        np.copyto(svecsr[:, 0:pp.numOrthoConst], lockv0[:, 0:pp.numOrthoConst])

    u0, v0 = check_pair(u0, v0, "v0 or u0")
    
    if v0 is not None:
        pp.initSize = min(v0.shape[1], pp.numSvals)
        np.copyto(svecsl[:, pp.numOrthoConst:pp.numOrthoConst+pp.initSize], u0[:, 0:pp.initSize])
        np.copyto(svecsr[:, pp.numOrthoConst:pp.numOrthoConst+pp.initSize], v0[:, 0:pp.initSize])

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

    svecsl = svecsl[:, pp.numOrthoConst:]
    svecsr = svecsr[:, pp.numOrthoConst:]

    # Transpose conjugate svecsr
    svecsr = svecsr.T.conj()

    return svecsl, svals, svecsr

# if __name__ == '__main__':
#     from scipy.sparse import spdiags
#     a = np.ones(10)
#     A  = spdiags(np.array([a*(-1.), a*2., a*(-1.)]), np.array([-1, 0, 1]), 10, 10)
#     print(eigsh_primme(A, which='LA'))
