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
-5: "'n' < 0 or 'nLocal' < 0 or 'nLocal' > 'n'",
-6: "'numProcs' < 1",
-7: "'matrixMatvec' is NULL",
-8: "'applyPreconditioner' is NULL and 'precondition' is not NULL",
-9: "'not used",
-10: "'numEvals' > 'n'",
-11: "'numEvals' < 0",
-12: "'eps' > 0 and 'eps' < machine precision",
-13: "'target' is not properly defined",
-14: "'target' is one of 'primme_largest_abs', 'primme_closest_geq', 'primme_closest_leq' or 'primme_closest_abs' but 'numTargetShifts' <= 0 (no shifts)",
-15: "'target' is one of 'primme_largest_abs', 'primme_closest_geq', 'primme_closest_leq' or 'primme_closest_abs' but 'targetShifts' is NULL  (no shifts array)",
-16: "'numOrthoConst' < 0 or 'numOrthoConst' > 'n'. (no free dimensions left)",
-17: "'maxBasisSize' < 2",
-18: "'minRestartSize' < 0 or 'minRestartSize' shouldn't be zero",
-19: "'maxBlockSize' < 0 or 'maxBlockSize' shouldn't be zero",
-20: "'maxPrevRetain' < 0",
-21: "'scheme' is not one of *primme_thick* or *primme_dtr*",
-22: "'initSize' < 0",
-23: "'locking' == 0 and 'initSize' > 'maxBasisSize'",
-24: "'locking' and 'initSize' > 'numEvals'",
-25: "'maxPrevRetain' + 'minRestartSize' >= 'maxBasisSize'",
-26: "'minRestartSize' >= 'n'",
-27: "'printLevel' < 0 or 'printLevel' > 5",
-28: "'convTest' is not one of 'primme_full_LTolerance', 'primme_decreasing_LTolerance', 'primme_adaptive_ETolerance' or 'primme_adaptive'",
-29: "'convTest' == 'primme_decreasing_LTolerance' and 'relTolBase' <= 1",
-30: "'evals' is NULL, but not 'evecs' and 'resNorms'",
-31: "'evecs' is NULL, but not 'evals' and 'resNorms'",
-32: "'resNorms' is NULL, but not 'evecs' and 'evals'",
-33: "'locking' == 0 and 'minRestartSize' < 'numEvals'",
-34: "'ldevecs' is less than 'nLocal'",
-35: "'ldOPs' is non-zero and less than 'nLocal'",
-36 : "not enough memory for realWork",
-37 : "not enough memory for intWork",
-38 : "'locking' == 0 and 'target' is 'primme_closest_leq' or 'primme_closet_geq'"
}

_PRIMMESvdsErrors = {
0   : "success",
1   : "reported only amount of required memory",
-1  : "failed in allocating int or real workspace",
-2  : "malloc failed in allocating a permutation integer array",
-3  : "main_iter() encountered problem; the calling stack of the functions where the error occurred was printed in 'stderr'",
-4  : "primme_svds is NULL",
-5  : "Wrong value for m or n or mLocal or nLocal",
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
-19 : "resNorms is not set",
-20 : "not enough memory for realWork",
-21 : "not enough memory for intWork"
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
        msg = ""
        if err >= -100:
            msg = _PRIMMESvdsErrors[err]
        elif err >= -200:
            msg = "Error from PRIMME first stage: " + _PRIMMEErrors[err+100]
        elif err >= -300:
            msg = "Error from PRIMME second stage: " + _PRIMMEErrors[err+200]
        RuntimeError.__init__(self, "PRIMME SVDS error %d: %s" % (err, msg))


def eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None,
          ncv=None, maxiter=None, tol=0, return_eigenvectors=True,
          Minv=None, OPinv=None, mode='normal', lock=None,
          return_stats=False, maxBlockSize=0, minRestartSize=0,
          maxPrevRetain=0, method=None, return_history=False, **kargs):
    """
    Find k eigenvalues and eigenvectors of the real symmetric square matrix
    or complex Hermitian matrix A.

    Solves ``A * x[i] = w[i] * x[i]``, the standard eigenvalue problem for
    w[i] eigenvalues with corresponding eigenvectors x[i].

    If M is specified, solves ``A * x[i] = w[i] * M * x[i]``, the
    generalized eigenvalue problem for w[i] eigenvalues
    with corresponding eigenvectors x[i]

    Parameters
    ----------
    A : An N x N matrix, array, sparse matrix, or LinearOperator
        the operation A * x, where A is a real symmetric matrix or complex
        Hermitian.
    k : int, optional
        The number of eigenvalues and eigenvectors desired.
    M : An N x N matrix, array, sparse matrix, or LinearOperator
        (not supported yet)
        the operation M * x for the generalized eigenvalue problem

            A * x = w * M * x.

        M must represent a real, symmetric matrix if A is real, and must
        represent a complex, hermitian matrix if A is complex. For best
        results, the data type of M should be the same as that of A.
    sigma : real, optional
        Find eigenvalues near sigma.
    v0 : N x i, ndarray, optional
        Starting vectors for iteration.
    ncv : int, optional
        The maximum size of the basis
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
        Maximum number of iterations.
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
    maxBlockSize : int, optional
        Maximum number of vectors added at every iteration.
    minRestartSize : int, optional
        Number of approximate eigenvectors kept from last iteration in restart.
    maxPrevRetain: int, optional
        Number of approximate eigenvectors kept from previous iteration in
        restart. Also referred as +k vectors in GD+k.
    method : int, optional
        Preset method, one of:

        - DEFAULT_MIN_TIME : a variant of JDQMR,
        - DEFAULT_MIN_MATVECS : GD+k
        - DYNAMIC : choose dynamically between both previous methods.

        See a detailed description of the methods and other possible values
        in [2]_.
        
    return_stats : bool, optional
        If True, it is also returned extra information from PRIMME.
    return_history: bool, optional
        If True, it is also returned performance information at every iteration.

    Returns
    -------
    w : array
        Array of k eigenvalues
    v : array
        An array representing the `k` eigenvectors.  The column ``v[:, i]`` is
        the eigenvector corresponding to the eigenvalue ``w[i]``.
    stats : dict, optional (if return_stats)
        Extra information reported by PRIMME:

        - "numOuterIterations": number of outer iterations
        - "numRestarts": number of restarts
        - "numMatvecs": number of A*v
        - "numPreconds": number of OPinv*v
        - "elapsedTime": time that took 
        - "estimateMinEVal": the leftmost Ritz value seen
        - "estimateMaxEVal": the rightmost Ritz value seen
        - "estimateLargestSVal": the largest singular value seen
        - "rnorms" : ||A*x[i] - x[i]*w[i]||
        - "hist" : (if return_history) report at every outer iteration of:

          - "elapsedTime": time spent up to now
          - "numMatvecs": number of A*v spent up to now
          - "nconv": number of converged pair
          - "eval": eigenvalue of the first unconverged pair
          - "resNorm": residual norm of the first unconverged pair

    Raises
    ------
    PrimmeError
        When the requested convergence is not obtained.

        The PRIMME error code can be found as ``err`` attribute of the exception
        object.

    See Also
    --------
    scipy.sparse.linalg.eigs : eigenvalues and eigenvectors for a general (nonsymmetric) matrix A
    Primme.svds : singular value decomposition for a matrix A

    Notes
    -----
    This function is a wrapper to PRIMME functions to find the eigenvalues and
    eigenvectors [1]_.

    References
    ----------
    .. [1] PRIMME Software, https://github.com/primme/primme
    .. [2] Preset Methods, http://www.cs.wm.edu/~andreas/software/doc/readme.html#preset-methods

    Examples
    --------
    >>> import Primme, scipy.sparse
    >>> A = scipy.sparse.spdiags(range(100), [0], 100, 100) # sparse diag. matrix
    >>> evals, evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA')
    >>> evals # the three largest eigenvalues of A
    array([ 99.,  98.,  97.])
    >>> evals, evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA', lock=evecs)
    >>> evals # the next three largest eigenvalues
    array([ 96.,  95.,  94.])
    """

    A = aslinearoperator(A)
    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('A: expected square matrix (shape=%s)' % (A.shape,))

    if M is not None:
        raise ValueError('generalized problems (M != None) are not supported')

    hist = {"numMatvecs": [], "elapsedTime": [], "nconv": [],
            "eval": [], "resNorm": []}

    class PP(PrimmeParams):
        def __init__(self):
            PrimmeParams.__init__(self)
        def matvec(self, X):
            return A.matmat(X)
        def prevec(self, X):
            return OPinv.matmat(X)
        def mon(self, basisEvals, basisFlags, iblock, basisNorms, numConverged,
                    lockedEvals, lockedFlags, lockedNorms, inner_its, LSRes, event):
            if event == 0 and len(iblock)>0: # event iteration
                hist["numMatvecs"].append(self.stats.numMatvecs)
                hist["elapsedTime"].append(self.stats.elapsedTime)
                hist["nconv"].append(numConverged)
                hist["eval"].append(basisEvals[iblock[0]])
                hist["resNorm"].append(basisNorms[iblock[0]])

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
        OPinv = aslinearoperator(OPinv)
        if OPinv.shape[0] != OPinv.shape[1] or OPinv.shape[0] != A.shape[0]:
            raise ValueError('OPinv: expected square matrix with same shape as A (shape=%s)' % (OPinv.shape,))
        pp.correctionParams.precondition = 1

    if lock is not None:
        if lock.shape[0] != pp.n:
            raise ValueError('lock: expected matrix with the same columns as A (shape=%s)' % (lock.shape,))
        pp.numOrthoConst = min(lock.shape[1], pp.n)

    if return_history and return_stats:
        pp.monitor_set = 1

    # Set other parameters
    for dk, dv in kargs.items():
      setattr(pp, dk, dv)

    if A.dtype.kind in frozenset(["b", "i", "u"]) or A.dtype.type is np.double:
        dtype = np.dtype("d")
    else:
        dtype = A.dtype

    if dtype.type is np.complex64:
        Xprimme = cprimme
        rtype = np.dtype(np.float32)
    elif dtype.type is np.float32:
        Xprimme = sprimme
        rtype = np.dtype(np.float32)
    elif dtype.type is np.float64:
        Xprimme = dprimme
        rtype = np.dtype(np.float64)
    else:
        Xprimme = zprimme
        rtype = np.dtype(np.float64)

    evals = np.zeros(pp.numEvals, rtype)
    norms = np.zeros(pp.numEvals, rtype)
    evecs = np.zeros((pp.n, pp.numOrthoConst+pp.numEvals), dtype, order='F')

    if lock is not None:
        np.copyto(evecs[:, 0:pp.numOrthoConst], lock[:, 0:pp.numOrthoConst])

    if v0 is not None:
        pp.initSize = min(v0.shape[1], pp.numEvals)
        np.copyto(evecs[:, pp.numOrthoConst:pp.numOrthoConst+pp.initSize],
            v0[:, 0:pp.initSize])

    if maxBlockSize:
        pp.maxBlockSize = maxBlockSize

    if minRestartSize:
        pp.minRestartSize = minRestartSize

    if maxPrevRetain:
        pp.restartingParams.maxPrevRetain = maxPrevRetain

    if method is not None:
        pp.set_method(method)
 
    err = Xprimme(evals, evecs, norms, pp)

    if err != 0:
        raise PrimmeError(err)

    evals = evals[0:pp.initSize]
    norms = norms[0:pp.initSize]
    evecs = evecs[:, pp.numOrthoConst:pp.numOrthoConst+pp.initSize]

    if return_stats:
        stats = dict((f, getattr(pp.stats, f)) for f in [
            "numOuterIterations", "numRestarts", "numMatvecs",
            "numPreconds", "elapsedTime", "estimateMinEVal",
            "estimateMaxEVal", "estimateLargestSVal"])
        stats['rnorms'] = norms
        if return_history:
            stats["hist"] = hist
        return evals, evecs, stats
    else:
        return evals, evecs


def svds(A, k=6, ncv=None, tol=0, which='LM', v0=None,
         maxiter=None, return_singular_vectors=True,
         precAHA=None, precAAH=None, precAug=None,
         u0=None, locku0=None, lockv0=None,
         return_stats=False, maxBlockSize=0,
         method=None, methodStage1=None, methodStage2=None,
         return_history=False, **kargs):
    """
    Compute k singular values and vectors for a sparse matrix.

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

    u0 : ndarray, optional
        Left starting vectors for the iterations.

        Should be approximate left singular vectors. If only u0 or v0 is
        provided, the other is computed.
    v0 : ndarray, optional
        Right starting vectors for the iterations.
    maxiter : int, optional
        Maximum number of iterations.
    precAHA : {N x N matrix, array, sparse matrix, LinearOperator}, optional
        Approximate inverse of (A.H*A - sigma**2*I). If provided and M>N, it
        usually accelerates the convergence.
    precAAH : {M x M matrix, array, sparse matrix, LinearOperator}, optional
        Approximate inverse of (A*A.H - sigma**2*I). If provided and M<N, it
        usually accelerates the convergence.
    precAug : {(M+N) x (M+N) matrix, array, sparse matrix, LinearOperator}, optional
        Approximate inverse of ([zeros() A.H; zeros() A] - sigma*I). It usually
        accelerates the convergence if tol<dtype.eps**.5.
    locku0 : ndarray, optional
        Left orthogonal vector constrain.

        Seek singular triplets orthogonal to locku0 and lockv0. The provided vectors
        *should* be orthonormal. If only locku0 or lockv0 is provided, the other
        is computed. Useful to not converge some already computed solutions.
    lockv0 : ndarray, optional
        Right orthogonal vector constrain. See locku0.
    maxBlockSize : int, optional
        Maximum number of vectors added at every iteration.
    return_stats : bool, optional
        If True, it is also returned extra information from PRIMME.
    return_history: bool, optional
        If True, it is also returned performance information at every iteration.

    Returns
    -------
    u : ndarray, shape=(M, k), optional
        Unitary matrix having left singular vectors as columns.
        Returned if `return_singular_vectors` is True.
    s : ndarray, shape=(k,)
        The singular values.
    vt : ndarray, shape=(k, N), optional
        Unitary matrix having right singular vectors as rows.
        Returned if `return_singular_vectors` is True.
    stats : dict, optional (if return_stats)
        Extra information reported by PRIMME:

        - "numOuterIterations": number of outer iterations
        - "numRestarts": number of restarts
        - "numMatvecs": number of A*v
        - "numPreconds": number of OPinv*v
        - "elapsedTime": time that took 
        - "rnorms" : ||A*v[i] - u[i]*s[i]||
        - "hist" : (if return_history) report at every outer iteration of:

          - "elapsedTime": time spent up to now
          - "numMatvecs": number of A*v spent up to now
          - "nconv": number of converged pair
          - "eval": eigenvalue of the first unconverged pair
          - "resNorm": residual norm of the first unconverged pair

    See Also
    --------
    Primme.eigsh : eigenvalue decomposition for a sparse symmetrix/complex Hermitian matrix A
    scipy.sparse.linalg.eigs : eigenvalues and eigenvectors for a general (nonsymmetric) matrix A

    Examples
    --------
    >>> import Primme, scipy.sparse
    >>> A = scipy.sparse.spdiags(range(1, 11), [0], 100, 10) # sparse diag. rect. matrix
    >>> svecs_left, svals, svecs_right = Primme.svds(A, 3, tol=1e-6, which='SM')
    >>> svals # the three smallest singular values of A
    array([ 1.,  2.,  3.])

    >>> import Primme, scipy.sparse
    >>> A = scipy.sparse.rand(10000, 100, random_state=10)
    >>> prec = scipy.sparse.spdiags(np.reciprocal(A.multiply(A).sum(axis=0)),
    ...           [0], 100, 100) # square diag. preconditioner
    >>> svecs_left, svals, svecs_right = Primme.svds(A, 3, which=6.0, tol=1e-6, precAHA=prec)
    >>> ["%.5f" % x for x in svals.flat] # the three closest singular values of A to 0.5
    ['5.99871', '5.99057', '6.01065']
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

    hist = {"numMatvecs": [], "elapsedTime": [], "nconv": [],
            "sval": [], "resNorm": []}

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
            return X

        def mon(self, basisSvals, basisFlags, iblock, basisNorms, numConverged,
                    lockedSvals, lockedFlags, lockedNorms, inner_its, LSRes,
                    event, stage):
            if event == 0 and len(iblock)>0: # event iteration
                hist["numMatvecs"].append(self.stats.numMatvecs)
                hist["elapsedTime"].append(self.stats.elapsedTime)
                hist["nconv"].append(numConverged)
                hist["sval"].append(basisSvals[iblock[0]])
                hist["resNorm"].append(basisNorms[iblock[0]])

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
        # NOTE: every eigensolver iteration spend two matvecs*blockSize
        pp.maxMatvecs = maxiter*(maxBlockSize if maxBlockSize else 1)/2

    if maxBlockSize:
        pp.maxBlockSize = maxBlockSize

    if precAHA is not None or precAAH is not None or precAug is not None:
        pp.precondition = 1

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

    if return_history and return_stats:
        pp.monitor_set = 1

    # Set other parameters
    for dk, dv in kargs.items():
      setattr(pp, dk, dv)

    if A.dtype.kind in frozenset(["b", "i", "u"]) or A.dtype.type is np.double:
        dtype = np.dtype("d")
    else:
        dtype = A.dtype

    if dtype.type is np.complex64:
        Xprimme_svds = cprimme_svds
        rtype = np.dtype(np.float32)
    elif dtype.type is np.float32:
        Xprimme_svds = sprimme_svds
        rtype = np.dtype(np.float32)
    elif dtype.type is np.float64:
        Xprimme_svds = dprimme_svds
        rtype = np.dtype(np.float64)
    else:
        Xprimme_svds = zprimme_svds
        rtype = np.dtype(np.float64)

    svals = np.zeros(pp.numSvals, rtype)
    svecsl = np.zeros((pp.m, pp.numOrthoConst+pp.numSvals), dtype, order='F')
    svecsr = np.zeros((pp.n, pp.numOrthoConst+pp.numSvals), dtype, order='F')
    norms = np.zeros(pp.numSvals, rtype)

    if locku0 is not None:
        np.copyto(svecsl[:, 0:pp.numOrthoConst], locku0[:, 0:pp.numOrthoConst])
        np.copyto(svecsr[:, 0:pp.numOrthoConst], lockv0[:, 0:pp.numOrthoConst])

    u0, v0 = check_pair(u0, v0, "v0 or u0")
    
    if v0 is not None:
        pp.initSize = min(v0.shape[1], pp.numSvals)
        np.copyto(svecsl[:, pp.numOrthoConst:pp.numOrthoConst+pp.initSize], u0[:, 0:pp.initSize])
        np.copyto(svecsr[:, pp.numOrthoConst:pp.numOrthoConst+pp.initSize], v0[:, 0:pp.initSize])

    # Set method
    if method is not None or methodStage1 is not None or methodStage2 is not None:
        if method is None: method = primme_svds_default
        if methodStage1 is None: methodStage1 = PRIMME_DEFAULT_METHOD
        if methodStage2 is None: methodStage2 = PRIMME_DEFAULT_METHOD
        pp.set_method(method, methodStage1, methodStage2)

    err = Xprimme_svds(svals, svecsl, svecsr, norms, pp)

    if err != 0:
        raise PrimmeSvdsError(err)

    if return_stats:
        stats = dict((f, getattr(pp.stats, f)) for f in [
            "numOuterIterations", "numRestarts", "numMatvecs",
            "numPreconds", "elapsedTime"])
        stats["rnorms"] = norms
        if return_history:
            stats["hist"] = hist
 
    if not return_singular_vectors:
        return svals if not return_stats else (svals, stats)

    svals = svals[0:pp.initSize]
    norms = norms[0:pp.initSize]
    svecsl = svecsl[:, pp.numOrthoConst:pp.numOrthoConst+pp.initSize]
    svecsr = svecsr[:, pp.numOrthoConst:pp.numOrthoConst+pp.initSize]

    # Transpose conjugate svecsr
    svecsr = svecsr.T.conj()

    if not return_stats:
        return svecsl, svals, svecsr
    else:
        return svecsl, svals, svecsr, stats
