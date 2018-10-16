
import numpy as np
cimport numpy as np
from scipy.sparse.linalg.interface import aslinearoperator
cimport cython
from cython cimport view

cdef extern from "../include/primme.h":
    struct primme_params:
        pass
    ctypedef int primme_preset_method
    ctypedef enum primme_type:
        primme_int, primme_double, primme_pointer
    ctypedef int primme_params_label
    ctypedef int primme_event
    int sprimme(float *evals, float *evecs, float *resNorms, primme_params *primme)
    int cprimme(float *evals, np.complex64_t *evecs, float *resNorms, primme_params *primme)
    int dprimme(double *evals, double *evecs, double *resNorms, primme_params *primme)
    int zprimme(double *evals, np.complex128_t *evecs, double *resNorms, primme_params *primme)
    int magma_sprimme(float *evals, float *evecs, float *resNorms, primme_params *primme)
    int magma_cprimme(float *evals, np.complex64_t *evecs, float *resNorms, primme_params *primme)
    int magma_dprimme(double *evals, double *evecs, double *resNorms, primme_params *primme)
    int magma_zprimme(double *evals, np.complex128_t *evecs, double *resNorms, primme_params *primme)
    primme_params* primme_params_create()
    int primme_params_destroy(primme_params *primme)
    void primme_initialize(primme_params *primme)
    int  primme_set_method(primme_preset_method method, primme_params *params)
    void primme_free(primme_params *primme)
    int primme_get_member(primme_params *primme, primme_params_label label, void *value)
    int primme_set_member(primme_params *primme, primme_params_label label, void *value)
    int primme_member_info(primme_params_label *label, const char** label_name, primme_type *t, int *arity)
    int primme_constant_info(const char* label_name, int *value)

ctypedef fused numerics:
    float
    double
    np.complex64_t
    np.complex128_t

ctypedef fused numerics_real:
    float
    double

cdef class PrimmeParams:
    cpdef primme_params *pp
    def __cinit__(self):
        self.pp = primme_params_create()
        if self.pp is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self.pp is not NULL:
            primme_params_destroy(self.pp)
    
def primme_params_get(PrimmeParams pp_, field_):
    cdef primme_params *primme = <primme_params*>(pp_.pp)
    cdef const char* field = <const char *>field_
    cdef primme_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_member_info(&l, <const char **>&field, &t, &arity)
    if not l >= 0 or arity != 1:
        raise "Invalid field '%s'" % field_
    cdef np.int64_t v_int
    cdef double v_double
    cdef void *v_pvoid
    if t == primme_int:
        primme_get_member(primme, l, &v_int)
        return v_int
    elif t == primme_double:
        primme_get_member(primme, l, &v_double)
        return v_double
    elif t == primme_pointer:
        primme_get_member(primme, l, &v_pvoid)
        return <object>v_pvoid
    else:
        raise ValueError("Not supported type for member '%s'" % field)

cdef object primme_params_get_object(primme_params *primme, cython.p_char field):
    cdef primme_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_member_info(&l, <const char **>&field, &t, &arity)
    assert l >= 0 and arity == 1 and t == primme_pointer, "Invalid field '%s'" % <bytes>field
    cdef void *v_pvoid
    primme_get_member(primme, l, &v_pvoid)
    return <object>v_pvoid

cdef np.int64_t primme_params_get_int(primme_params *primme, cython.p_char field):
    cdef primme_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_member_info(&l, <const char **>&field, &t, &arity)
    assert l >= 0 and arity == 1 and t == primme_int, "Invalid field '%s'" % <bytes>field
    cdef np.int64_t v_int
    primme_get_member(primme, l, &v_int)
    return v_int

def primme_params_set(PrimmeParams pp_, field_, value):
    cdef primme_params *primme = <primme_params*>(pp_.pp)
    cdef const char* field = <const char *>field_
    cdef primme_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_member_info(&l, <const char **>&field, &t, &arity)
    if l< 0 or arity != 1:
        raise "Invalid field '%s'" % field_
    cdef np.int64_t v_int
    cdef double v_double
    cdef int i
    if t == primme_pointer:
        primme_set_member(primme, l, <void*>value)
    elif t == primme_int:
        if isinstance(value, str):
            primme_constant_info(<const char*>value, &i)
            value = i
        v_int = value
        primme_set_member(primme, l, &v_int)
    elif t == primme_double:
        v_double = value
        primme_set_member(primme, l, &v_double)
    else:
        raise ValueError("Not supported type for member '%s'" % field)
   
cdef void primme_params_set_pointer(primme_params *primme, cython.p_char field, void* value):
    cdef primme_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_member_info(&l, <const char **>&field, &t, &arity)
    assert(l >= 0 and arity == 1 and t == primme_pointer)
    primme_set_member(primme, l, value)

cdef void primme_params_set_doubles(primme_params *primme, cython.p_char field, double *value):
    cdef primme_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_member_info(&l, <const char **>&field, &t, &arity)
    assert(l >= 0 and arity == 0 and t == primme_double)
    primme_set_member(primme, l, value)


cdef np.dtype get_np_type(numerics *p):
    cdef double *pd
    cdef float *pf
    cdef np.complex128_t *pcd
    cdef np.complex64_t *pcf
    if cython.typeof(p) == cython.typeof(pd):
        return np.dtype(np.double)
    elif cython.typeof(p) == cython.typeof(pf):
        return np.dtype(np.float32)
    elif cython.typeof(p) == cython.typeof(pcd):
        return np.dtype(np.complex128)
    elif cython.typeof(p) == cython.typeof(pcf):
        return np.dtype(np.complex64)

cdef void c_matvec_numpy(numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
    if blockSize[0] <= 0:
        ierr[0] = 0
        return
    ierr[0] = 1
    cdef object matvec = primme_params_get_object(primme, 'matrix')
    n = primme_params_get_int(primme, "nLocal")
    x_py = np.asarray(<numerics[:ldx[0]:1, :blockSize[0]]>x)[0:n,:]
    cdef np.ndarray[numerics, ndim=2, mode="strided"] y_py
    y_py = matvec(x_py)
    (<numerics[:ldy[0]:1, :blockSize[0]]>y)[0:n,:] = y_py[:,:]
    ierr[0] = 0

cdef void c_massmatvec_numpy(numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
    if blockSize[0] <= 0:
        ierr[0] = 0
        return
    ierr[0] = 1
    cdef object massmatvec = primme_params_get_object(primme, 'massMatrix')
    n = primme_params_get_int(primme, "nLocal")
    x_py = np.asarray(<numerics[:ldx[0]:1, :blockSize[0]]>x)[0:n,:]
    cdef np.ndarray[numerics, ndim=2, mode="strided"] y_py
    y_py = massmatvec(x_py)
    (<numerics[:ldy[0]:1, :blockSize[0]]>y)[0:n,:] = y_py[:,:]
    ierr[0] = 0

cdef void c_precond_numpy(numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
    if blockSize[0] <= 0:
        ierr[0] = 0
        return
    ierr[0] = 1
    cdef object precond = primme_params_get_object(primme, 'preconditioner')
    n = primme_params_get_int(primme, "nLocal")
    x_py = np.asarray(<numerics[:ldx[0]:1, :blockSize[0]]>x)[0:n,:]
    cdef np.ndarray[numerics, ndim=2, mode="strided"] y_py
    y_py = precond(x_py)
    (<numerics[:ldy[0]:1, :blockSize[0]]>y)[0:n,:] = y_py[:,:]
    ierr[0] = 0

cdef void c_monitor(numerics_real *basisEvals, int *basisSize, int *basisFlags, int *iblock, int *blockSize, numerics_real *basisNorms, int *numConverged, numerics_real *lockedEvals, int *numLocked, int *lockedFlags, numerics_real *lockedNorms, int *inner_its, numerics_real *LSRes, primme_event *event, primme_params *primme, int *ierr):
    ierr[0] = 1
    cdef object monitor = primme_params_get_object(primme, 'monitor')
    cdef int bs = basisSize[0] if basisSize is not NULL else 0
    cdef int blks = blockSize[0] if blockSize is not NULL else 0
    cdef int nLocked = numLocked[0] if numLocked is not NULL else 0
    monitor(
        <numerics_real[:bs]>basisEvals if basisEvals is not NULL and bs > 0 else None,
        <int[:bs]>basisFlags if basisFlags is not NULL and bs > 0 else None,
        <int[:blks]>iblock if iblock is not NULL and blks > 0 else None,
        <numerics_real[:blks]>basisNorms if basisNorms is not NULL and blks > 0 else None,
        numConverged[0] if numConverged is not NULL else None,
        <numerics_real[:nLocked]>lockedEvals if lockedEvals is not NULL and nLocked > 0 else None,
        <int[:nLocked]>lockedFlags if lockedFlags is not NULL and nLocked > 0 else None,
        <numerics_real[:nLocked]>lockedNorms if lockedNorms is not NULL and nLocked > 0 else None,
        inner_its[0] if inner_its is not NULL else None,
        LSRes[0] if LSRes is not NULL else None,
        event[0])
    ierr[0] = 0


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
        The number of eigenvalues and eigenvectors to be computed. Must be
        1 <= k < min(A.shape).
    M : An N x N matrix, array, sparse matrix, or LinearOperator
        (not supported yet)
        the operation M * x for the generalized eigenvalue problem

            A * x = w * M * x.

        M must represent a real, symmetric matrix if A is real, and must
        represent a complex, Hermitian matrix if A is complex. For best
        results, the data type of M should be the same as that of A.
    sigma : real, optional
        Find eigenvalues near sigma.
    v0 : N x i, ndarray, optional
        Initial guesses to the eigenvectors.
    ncv : int, optional
        The maximum size of the basis
    which : str ['LM' | 'SM' | 'LA' | 'SA']
        Which `k` eigenvectors and eigenvalues to find:

            'LM' : Largest in magnitude eigenvalues; the farthest from sigma

            'SM' : Smallest in magnitude eigenvalues; the closest to sigma

            'LA' : Largest algebraic eigenvalues

            'SA' : Smallest algebraic eigenvalues

            'CLT' : closest but left to sigma

            'CGT' : closest but greater than sigma

        When sigma == None, 'LM', 'SM', 'CLT', and 'CGT' treat sigma as zero. 
    maxiter : int, optional
        Maximum number of iterations.
    tol : float
        Required accuracy for eigenpairs (stopping criterion).
        The default value is sqrt of machine precision.
    Minv : (not supported yet)
        The inverse of M in the generalized eigenproblem.
    OPinv : N x N matrix, array, sparse matrix, or LinearOperator, optional
        Preconditioner to accelerate the convergence. Usually it is an
        approximation of the inverse of (A - sigma*M).
    return_eigenvectors : bool, optional
        Return eigenvectors (True) in addition to eigenvalues
    mode : string ['normal' | 'buckling' | 'cayley']
        Only 'normal' mode is supported.
    lock : N x i, ndarray, optional
        Seek the eigenvectors orthogonal to these ones. The provided
        vectors *should* be orthonormal. Useful to avoid converging to previously
        computed solutions.
    maxBlockSize : int, optional
        Maximum number of vectors added at every iteration.
    minRestartSize : int, optional
        Number of approximate eigenvectors kept during restart.
    maxPrevRetain: int, optional
        Number of approximate eigenvectors kept from previous iteration in
        restart. Also referred as +k vectors in GD+k.
    method : int, optional
        Preset method, one of:

        - DEFAULT_MIN_TIME : a variant of JDQMR,
        - DEFAULT_MIN_MATVECS : GD+k
        - DYNAMIC : choose dynamically between these previous methods.

        See a detailed description of the methods and other possible values
        in [2]_.
        
    return_stats : bool, optional
        If True, the function returns extra information (see stats in Returns).
    return_history: bool, optional
        If True, the function returns performance information at every iteration
        (see hist in Returns).

    Returns
    -------
    w : array
        Array of k eigenvalues ordered to best satisfy "which".
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
    primme.svds : singular value decomposition for a matrix A

    Notes
    -----
    This function is a wrapper to PRIMME functions to find the eigenvalues and
    eigenvectors [1]_.

    References
    ----------
    .. [1] PRIMME Software, https://github.com/primme/primme
    .. [2] Preset Methods, http://www.cs.wm.edu/~andreas/software/doc/readme.html#preset-methods
    .. [3] A. Stathopoulos and J. R. McCombs PRIMME: PReconditioned
           Iterative MultiMethod Eigensolver: Methods and software
           description, ACM Transaction on Mathematical Software Vol. 37,
           No. 2, (2010), 21:1-21:30.

    Examples
    --------
    >>> import primme, scipy.sparse
    >>> A = scipy.sparse.spdiags(range(100), [0], 100, 100) # sparse diag. matrix
    >>> evals, evecs = primme.eigsh(A, 3, tol=1e-6, which='LA')
    >>> evals # the three largest eigenvalues of A
    array([99., 98., 97.])
    >>> new_evals, new_evecs = primme.eigsh(A, 3, tol=1e-6, which='LA', lock=evecs)
    >>> new_evals # the next three largest eigenvalues
    array([96., 95., 94.])
    """

    PP = PrimmeParams()
    cdef primme_params *pp = PP.pp
 
    A = aslinearoperator(A)
    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('A: expected square matrix (shape=%s)' % A.shape)

    primme_params_set(PP, "matrix", A)
    n = A.shape[0]
    primme_params_set(PP, "n", n)

    if M is not None:
        M = aslinearoperator(M)
        if len(M.shape) != 2 or A.shape[0] != M.shape[0]:
            raise ValueError('M: expected square matrix (shape=%s)' % A.shape)
        primme_params_set(PP, "massMatrix", M)

    if k <= 0 or k > n:
        raise ValueError("k=%d must be between 1 and %d, the order of the "
                         "square input matrix." % (k, n))
    primme_params_set(PP, "numEvals", k)

    if which == 'LM':
        primme_params_set(PP, "target", "primme_largest_abs")
        if sigma is None:
            sigma = 0.0
    elif which == 'LA':
        primme_params_set(PP, "target", "primme_largest")
        sigma = None
    elif which == 'SA':
        primme_params_set(PP, "target", "primme_smallest")
        sigma = None
    elif which == 'SM':
        primme_params_set(PP, "target", "primme_closest_abs")
        if sigma is None:
            sigma = 0.0
    elif which == 'CLT':
        primme_params_set(PP, "target", "primme_closest_leq")
        if sigma is None:
            sigma = 0.0
    elif which == 'CGT':
        primme_params_set(PP, "target", "primme_closest_geq")
        if sigma is None:
            sigma = 0.0
    else:
        raise ValueError("which='%s' not supported" % which)

    cdef double sigma_c
    if sigma is not None:
        sigma_c = float(sigma)
        primme_params_set(PP, "numTargetShifts", 1)
        primme_params_set_doubles(pp, "targetShifts", &sigma_c)

    primme_params_set(PP, "eps", tol)

    if ncv is not None:
        primme_params_set(PP, "maxBasisSize", ncv)

    if maxiter is not None:
        primme_params_set(PP, "maxMatvecs", maxiter)

    if OPinv is not None:
        OPinv = aslinearoperator(OPinv)
        if OPinv.shape[0] != OPinv.shape[1] or OPinv.shape[0] != A.shape[0]:
            raise ValueError('OPinv: expected square matrix with same shape as A (shape=%s)' % (OPinv.shape,n))
        primme_params_set(PP, "correction_precondition", 1)
        primme_params_set(PP, "preconditioner", <object>OPinv)
    else:
        primme_params_set(PP, "correction_precondition", 0)

    numOrthoConst = 0
    if lock is not None:
        if lock.shape[0] != n:
            raise ValueError('lock: expected matrix with the same columns as A (shape=%s)' % (lock.shape,n))
        numOrthoConst = min(lock.shape[1], n)
        primme_params_set(PP, "numOrthoConst", numOrthoConst)

    # Set other parameters
    for dk, dv in kargs.items():
      try:
        primme_params_set(PP, dk, dv)
      except:
        raise ValueError("Invalid option '%s' with value '%s'" % (dk, dv))

    hist = {"numMatvecs": [], "elapsedTime": [], "nconv": [],
            "eval": [], "resNorm": []}

    def mon(basisEvals, basisFlags, iblock, basisNorms, numConverged,
            lockedEvals, lockedFlags, lockedNorms, inner_its, LSRes, event):
        
        if event == 0 and iblock and len(iblock)>0: # event iteration
            hist["numMatvecs"].append(primme_params_get(PP, 'stats_numMatvecs'))
            hist["elapsedTime"].append(primme_params_get(PP, 'stats_elapsedTime'))
            hist["nconv"].append(numConverged)
            hist["eval"].append(basisEvals[iblock[0]])
            hist["resNorm"].append(basisNorms[0])

    if return_history:
        primme_params_set(PP, 'monitor', mon)

    if A.dtype.kind in frozenset(["b", "i", "u"]) or A.dtype.type is np.double:
        dtype = np.dtype("d")
    else:
        dtype = A.dtype

    if dtype.type is np.complex64:
        rtype = np.dtype(np.float32)
        primme_params_set_pointer(pp, "matrixMatvec", <void*>c_matvec_numpy[np.complex64_t])
        if M: 
            primme_params_set_pointer(pp, "massMatrixMatvec", <void*>c_massmatvec_numpy[np.complex64_t])
        primme_params_set_pointer(pp, "applyPreconditioner", <void*>c_precond_numpy[np.complex64_t])
        if return_history:
            primme_params_set_pointer(pp, "monitorFun", <void*>c_monitor[float])
    elif dtype.type is np.float32:
        rtype = np.dtype(np.float32)
        primme_params_set_pointer(pp, "matrixMatvec", <void*>c_matvec_numpy[float])
        if M: 
            primme_params_set_pointer(pp, "massMatrixMatvec", <void*>c_massmatvec_numpy[float])
        primme_params_set_pointer(pp, "applyPreconditioner", <void*>c_precond_numpy[float])
        if return_history:
            primme_params_set_pointer(pp, "monitorFun", <void*>c_monitor[float])
    elif dtype.type is np.float64:
        rtype = np.dtype(np.float64)
        primme_params_set_pointer(pp, "matrixMatvec", <void*>c_matvec_numpy[double])
        if M: 
            primme_params_set_pointer(pp, "massMatrixMatvec", <void*>c_massmatvec_numpy[double])
        primme_params_set_pointer(pp, "applyPreconditioner", <void*>c_precond_numpy[double])
        if return_history:
            primme_params_set_pointer(pp, "monitorFun", <void*>c_monitor[double])
    else:
        rtype = np.dtype(np.float64)
        primme_params_set_pointer(pp, "matrixMatvec", <void*>c_matvec_numpy[np.complex128_t])
        if M: 
            primme_params_set_pointer(pp, "massMatrixMatvec", <void*>c_massmatvec_numpy[np.complex128_t])
        primme_params_set_pointer(pp, "applyPreconditioner", <void*>c_precond_numpy[np.complex128_t])
        if return_history:
            primme_params_set_pointer(pp, "monitorFun", <void*>c_monitor[double])

    evals = np.zeros(k, rtype)
    norms = np.zeros(k, rtype)
    evecs = np.zeros((n, numOrthoConst+k), dtype, order='F')

    if lock is not None:
        np.copyto(evecs[:, 0:numOrthoConst], lock[:, 0:numOrthoConst])

    if v0 is not None:
        initSize = min(v0.shape[1], k)
        primme_params_set(PP, "initSize", initSize)
        np.copyto(evecs[:, numOrthoConst:numOrthoConst+initSize],
            v0[:, 0:initSize])

    if maxBlockSize:
        primme_params_set(PP, "maxBlockSize", maxBlockSize)

    if minRestartSize:
        primme_params_set(PP, "minRestartSize", minRestartSize)

    if maxPrevRetain:
        primme_params_set(PP, "restarting_maxPrevRetain", maxPrevRetain)

    cdef int method_int = -1;
    if method is not None:
        primme_constant_info(<const char *>method, &method_int)
        if method_int < 0:
            raise ValueError('Not valid "method": %s' % method)
        primme_set_method(method_int, pp)
 
    cdef np.ndarray[double, ndim=1, mode="c"] evals_d, norms_d
    cdef np.ndarray[float, ndim=1, mode="c"] evals_f, norms_f
    cdef np.ndarray[float, ndim=2, mode="fortran"] evecs_f
    cdef np.ndarray[double, ndim=2, mode="fortran"] evecs_d
    cdef np.ndarray[np.complex64_t, ndim=2, mode="fortran"] evecs_c
    cdef np.ndarray[np.complex128_t, ndim=2, mode="fortran"] evecs_z

    if dtype.type is np.complex64:
        evals_f, norms_f, evecs_c = evals, norms, evecs
        err = cprimme(&evals_f[0], &evecs_c[0,0], &norms_f[0], pp)
    elif dtype.type is np.float32:
        evals_f, norms_f, evecs_f = evals, norms, evecs
        err = sprimme(&evals_f[0], &evecs_f[0,0], &norms_f[0], pp)
    elif dtype.type is np.float64:
        evals_d, norms_d, evecs_d = evals, norms, evecs
        err = dprimme(&evals_d[0], &evecs_d[0,0], &norms_d[0], pp)
    else:
        evals_d, norms_d, evecs_z = evals, norms, evecs
        err = zprimme(&evals_d[0], &evecs_z[0,0], &norms_d[0], pp)

    if err != 0:
        raise PrimmeError(err)

    initSize = primme_params_get(PP, "initSize")
    evals = evals[0:initSize]
    norms = norms[0:initSize]
    evecs = evecs[:, numOrthoConst:numOrthoConst+initSize]

    if return_stats:
        stats = dict((f, primme_params_get(PP, "stats_" + f)) for f in [
            "numOuterIterations", "numRestarts", "numMatvecs",
            "numPreconds", "elapsedTime", "estimateMinEVal",
            "estimateMaxEVal", "estimateLargestSVal"])
        stats['rnorms'] = norms
        if return_history:
            stats["hist"] = hist
        return evals, evecs, stats
    else:
        return evals, evecs



cdef extern from "../include/primme.h":
    struct primme_svds_params:
        pass
    ctypedef int primme_svds_preset_method
    ctypedef int primme_svds_params_label
    ctypedef enum primme_svds_operator:
        primme_svds_op_none,
        primme_svds_op_AtA,
        primme_svds_op_AAt,
        primme_svds_op_augmented
    int sprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds)
    int cprimme_svds(float *svals, np.complex64_t *svecs, float *resNorms, primme_svds_params *primme_svds)
    int dprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds)
    int zprimme_svds(double *svals, np.complex128_t *svecs, double *resNorms, primme_svds_params *primme_svds)
    int magma_sprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds)
    int magma_cprimme_svds(float *svals,  np.complex64_t *svecs, float *resNorms, primme_svds_params *primme_svds)
    int magma_dprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds)
    int magma_zprimme_svds(double *svals,  np.complex128_t *svecs, double *resNorms, primme_svds_params *primme_svds)
    primme_svds_params* primme_svds_params_create()
    int primme_svds_params_destroy(primme_svds_params *primme_svds)
    void primme_svds_initialize(primme_svds_params *primme_svds)
    int primme_svds_set_method(primme_svds_preset_method method, primme_preset_method methodStage1, primme_preset_method methodStage2, primme_svds_params *primme_svds)
    void primme_svds_free(primme_svds_params *primme_svds)
    int primme_svds_get_member(primme_svds_params *primme_svds, primme_svds_params_label label, void *value)
    int primme_svds_set_member(primme_svds_params *primme_svds, primme_svds_params_label label, void *value)
    int primme_svds_member_info(primme_svds_params_label *label, const char** label_name, primme_type *t, int *arity)
    int primme_svds_constant_info(const char* label_name, int *value)


cdef class PrimmeSvdsParams:
    cpdef primme_svds_params *pp
    def __cinit__(self):
        self.pp = primme_svds_params_create()
        if self.pp is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self.pp is not NULL:
            primme_svds_params_destroy(self.pp)
    
def primme_svds_params_get(PrimmeSvdsParams pp_, field_):
    cdef primme_svds_params *primme_svds = <primme_svds_params*>(pp_.pp)
    cdef const char* field = <const char *>field_
    cdef primme_svds_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_svds_member_info(&l, <const char **>&field, &t, &arity)
    if not l >= 0 or arity != 1:
        raise "Invalid field '%s'" % field_
    cdef np.int64_t v_int
    cdef double v_double
    cdef void *v_pvoid
    if t == primme_int:
        primme_svds_get_member(primme_svds, l, &v_int)
        return v_int
    elif t == primme_double:
        primme_svds_get_member(primme_svds, l, &v_double)
        return v_double
    elif t == primme_pointer:
        primme_svds_get_member(primme_svds, l, &v_pvoid)
        return <object>v_pvoid
    else:
        raise ValueError("Not supported type for member '%s'" % field)

cdef object primme_svds_params_get_object(primme_svds_params *primme_svds, char *field):
    cdef primme_svds_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_svds_member_info(&l, <const char **>&field, &t, &arity)
    assert l >= 0 and arity == 1 and t == primme_pointer, "Invalid field '%s'" % <bytes>field
    cdef void *v_pvoid
    primme_svds_get_member(primme_svds, l, &v_pvoid)
    return <object>v_pvoid

cdef np.int64_t primme_svds_params_get_int(primme_svds_params *primme_svds, cython.p_char field):
    cdef primme_svds_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_svds_member_info(&l, <const char **>&field, &t, &arity)
    assert l >= 0 and arity == 1 and t == primme_int, "Invalid field '%s'" % <bytes>field
    cdef np.int64_t v_int
    primme_svds_get_member(primme_svds, l, &v_int)
    return v_int

def primme_svds_params_set(PrimmeSvdsParams pp_, field_, value):
    cdef primme_svds_params *primme_svds = <primme_svds_params*>(pp_.pp)
    cdef const char* field = <const char *>field_
    cdef primme_svds_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_svds_member_info(&l, <const char **>&field, &t, &arity)
    if l< 0 or arity != 1:
        raise "Invalid field '%s'" % field_
    cdef np.int64_t v_int
    cdef double v_double
    cdef int i
    if t == primme_pointer:
        primme_svds_set_member(primme_svds, l, <void*>value)
    elif t == primme_int:
        if isinstance(value, str):
            primme_svds_constant_info(<const char*>value, &i)
            value = i
        v_int = value
        primme_svds_set_member(primme_svds, l, &v_int)
    elif t == primme_double:
        v_double = value
        primme_svds_set_member(primme_svds, l, &v_double)
    else:
        raise ValueError("Not supported type for member '%s'" % field)
   
cdef void primme_svds_params_set_pointer(primme_svds_params *primme_svds, cython.p_char field, void* value):
    cdef primme_svds_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_svds_member_info(&l, <const char **>&field, &t, &arity)
    assert(l >= 0 and arity == 1 and t == primme_pointer)
    primme_svds_set_member(primme_svds, l, value)

cdef void primme_svds_params_set_doubles(primme_svds_params *primme_svds, cython.p_char field, double *value):
    cdef primme_svds_params_label l = -1
    cdef primme_type t
    cdef int arity
    primme_svds_member_info(&l, <const char **>&field, &t, &arity)
    assert(l >= 0 and arity == 0 and t == primme_double)
    primme_svds_set_member(primme_svds, l, value)


cdef void c_svds_matvec_numpy(numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, int *transpose, primme_svds_params *primme_svds, int *ierr):
    if blockSize[0] <= 0:
        ierr[0] = 0
        return
    ierr[0] = 1
    cdef object A = primme_svds_params_get_object(primme_svds, 'matrix')
    m = primme_svds_params_get_int(primme_svds, "mLocal")
    n = primme_svds_params_get_int(primme_svds, "nLocal")
    cdef np.ndarray[numerics, ndim=2, mode="strided"] y_py
    if transpose[0] == 0:
        x_py = np.asarray(<numerics[:ldx[0]:1, :blockSize[0]]>x)[0:n,:]
        y_py = A.matmat(x_py)
        (<numerics[:ldy[0]:1, :blockSize[0]]>y)[0:m,:] = y_py[:,:]
    else:
        x_py = np.asarray(<numerics[:ldx[0]:1, :blockSize[0]]>x)[0:m,:]
        y_py = A.H.matmat(x_py)
        (<numerics[:ldy[0]:1, :blockSize[0]]>y)[0:n,:] = y_py[:,:]
    ierr[0] = 0

cdef void c_svds_precond_numpy(numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, primme_svds_operator *mode, primme_svds_params *primme_svds, int *ierr):
    if blockSize[0] <= 0:
        ierr[0] = 0
        return
    ierr[0] = 1
    cdef object precond = primme_svds_params_get_object(primme_svds, 'preconditioner')
    m = primme_svds_params_get_int(primme_svds, "mLocal")
    n = primme_svds_params_get_int(primme_svds, "nLocal")
    cdef np.ndarray[numerics, ndim=2, mode="strided"] y_py
    if mode[0] == primme_svds_op_AtA:
        x_py = np.asarray(<numerics[:ldx[0]:1, :blockSize[0]]>x)[0:n,:]
        y_py = precond(x_py, mode[0])
        (<numerics[:ldy[0]:1, :blockSize[0]]>y)[0:n,:] = y_py[:,:]
    elif mode[0] == primme_svds_op_AAt:
        x_py = np.asarray(<numerics[:ldx[0]:1, :blockSize[0]]>x)[0:m,:]
        y_py = precond(x_py, mode[0])
        (<numerics[:ldy[0]:1, :blockSize[0]]>y)[0:m,:] = y_py[:,:]
    elif mode[0] == primme_svds_op_augmented:
        x_py = np.asarray(<numerics[:ldx[0]:1, :blockSize[0]]>x)[0:m+n,:]
        y_py = precond(x_py, mode[0])
        (<numerics[:ldy[0]:1, :blockSize[0]]>y)[0:m+n,:] = y_py[:,:]
    else:
        return
    ierr[0] = 0

cdef void c_svds_monitor(numerics_real *basisSvals, int *basisSize, int *basisFlags, int *iblock, int *blockSize,
      numerics_real *basisNorms, int *numConverged, numerics_real *lockedSvals, int *numLocked, int *lockedFlags, numerics_real *lockedNorms,
      int *inner_its, numerics_real *LSRes, primme_event *event, int *stage, primme_svds_params *primme_svds, int *ierr):
    ierr[0] = 1
    cdef object monitor = primme_svds_params_get_object(primme_svds, 'monitor')
    cdef int blks = blockSize[0] if blockSize is not NULL else 0
    cdef int bs = basisSize[0] if basisSize is not NULL else 0
    cdef int nLocked = numLocked[0] if numLocked is not NULL else 0
    monitor(
        <numerics_real[:bs]>basisSvals if basisSvals is not NULL and bs > 0 else None,
        <int[:bs]>basisFlags if basisFlags is not NULL and bs > 0 else None,
        <int[:blks]>iblock if iblock is not NULL and blks > 0 else None,
        <numerics_real[:blks]>basisNorms if basisNorms is not NULL and blks > 0 else None,
        numConverged[0] if numConverged is not NULL else None,
        <numerics_real[:nLocked]>lockedSvals if lockedSvals is not NULL and nLocked > 0 else None,
        <int[:nLocked]>lockedFlags if lockedFlags is not NULL and nLocked > 0 else None,
        <numerics_real[:nLocked]>lockedNorms if lockedNorms is not NULL and nLocked > 0 else None,
        inner_its[0] if inner_its is not NULL else None,
        LSRes[0] if LSRes is not NULL else None,
        event[0], stage[0])
    ierr[0] = 0


def svds(A, k=6, ncv=None, tol=0, which='LM', v0=None,
         maxiter=None, return_singular_vectors=True,
         precAHA=None, precAAH=None, precAug=None,
         u0=None, orthou0=None, orthov0=None,
         return_stats=False, maxBlockSize=0,
         method=None, methodStage1=None, methodStage2=None,
         return_history=False, **kargs):
    """
    Compute k singular values and vectors of the matrix A.

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
        Tolerance for singular values. Zero (default) means 10**4 times the machine precision.

        A triplet ``(u,sigma,v)`` is marked as converged when
        (||A*v - sigma*u||**2 + ||A.H*u - sigma*v||**2)**.5
        is less than "tol" * ||A||, or close to the minimum tolerance that
        the method can achieve. See the note.
    which : str ['LM' | 'SM'] or number, optional
        Which `k` singular values to find:

            - 'LM' : largest singular values
            - 'SM' : smallest singular values
            - number : closest singular values to (referred as sigma later)

    u0 : ndarray, optional
        Initial guesses for the left singular vectors.

        If only u0 or v0 is provided, the other is computed. If both are
        provided, u0 and v0 should have the same number of columns.
    v0 : ndarray, optional
        Initial guesses for the right singular vectors.
    maxiter : int, optional
        Maximum number of matvecs with A and A.H. 
    precAHA : {N x N matrix, array, sparse matrix, LinearOperator}, optional
        Approximate inverse of (A.H*A - sigma**2*I). If provided and M>=N, it
        usually accelerates the convergence.
    precAAH : {M x M matrix, array, sparse matrix, LinearOperator}, optional
        Approximate inverse of (A*A.H - sigma**2*I). If provided and M<N, it
        usually accelerates the convergence.
    precAug : {(M+N) x (M+N) matrix, array, sparse matrix, LinearOperator}, optional
        Approximate inverse of ([zeros() A.H; zeros() A] - sigma*I).
    orthou0 : ndarray, optional
        Left orthogonal vector constrain.

        Seek singular triplets orthogonal to orthou0 and orthov0. The provided vectors
        *should* be orthonormal. If only orthou0 or orthov0 is provided, the other
        is computed. Useful to avoid converging to previously computed solutions.
    orthov0 : ndarray, optional
        Right orthogonal vector constrain. See orthou0.
    maxBlockSize : int, optional
        Maximum number of vectors added at every iteration.
    return_stats : bool, optional
        If True, the function returns extra information (see stats in Returns).
    return_history: bool, optional
        If True, the function returns performance information at every iteration

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
        - "numMatvecs": number of matvecs with A and A.H
        - "numPreconds": cumulative number of applications of precAHA, precAAH
          and precAug
        - "elapsedTime": time that took 
        - "rnorms" : (||A*v[:,i] - sigma[i]*u[:,i]||**2 + ||A.H*u[:,i] - sigma[i]*v[:,i]||**2)**.5
        - "hist" : (if return_history) report at every outer iteration of:

          - "elapsedTime": time spent up to now
          - "numMatvecs": number of A*v and A.H*v spent up to now
          - "nconv": number of converged triplets
          - "sval": singular value of the first unconverged triplet
          - "resNorm": residual norm of the first unconverged triplet

    Notes
    -----
    The default method used is the hybrid method, which first solves the
    equivalent eigenvalue problem A.H*A or A*A.H (normal equations) and then
    refines the solution solving the augmented problem. The minimum tolerance
    that this method can achieve is ||A||*epsilon, where epsilon is the
    machine precision. However it may not return triplets with singular values
    smaller than ||A||*epsilon if "tol" is smaller than ||A||*epsilon/sigma.

    This function is a wrapper to PRIMME functions to find singular values and
    vectors [1]_.

    References
    ----------
    .. [1] PRIMME Software, https://github.com/primme/primme

    .. [2] L. Wu, E. Romero and A. Stathopoulos, PRIMME_SVDS: A High-
           Performance Preconditioned SVD Solver for Accurate Large-Scale
           Computations. https://arxiv.org/abs/1607.01404

    See Also
    --------
    primme.eigsh : eigenvalue decomposition for a sparse symmetrix/complex Hermitian matrix A
    scipy.sparse.linalg.eigs : eigenvalues and eigenvectors for a general (nonsymmetric) matrix A

    Examples
    --------
    >>> import primme, scipy.sparse
    >>> A = scipy.sparse.spdiags(range(1, 11), [0], 100, 10) # sparse diag. rect. matrix
    >>> svecs_left, svals, svecs_right = primme.svds(A, 3, tol=1e-6, which='SM')
    >>> svals # the three smallest singular values of A
    array([1., 2., 3.])

    >>> import primme, scipy.sparse
    >>> A = scipy.sparse.rand(10000, 100, random_state=10)
    >>> prec = scipy.sparse.spdiags(np.reciprocal(A.multiply(A).sum(axis=0)),
    ...           [0], 100, 100) # square diag. preconditioner
    >>> svecs_left, svals, svecs_right = primme.svds(A, 3, which=6.0, tol=1e-6, precAHA=prec)
    >>> ["%.5f" % x for x in svals.flat] # the three closest singular values of A to 0.5
    ['5.99871', '5.99057', '6.01065']
    """
    PP = PrimmeSvdsParams()
    cdef primme_svds_params *pp = PP.pp
 
    A = aslinearoperator(A)

    cdef int m, n
    m, n = A.shape
    primme_svds_params_set(PP, "matrix", A)
    primme_svds_params_set(PP, "m", m)
    primme_svds_params_set(PP, "n", n)

    if k <= 0 or k > min(n, m):
        raise ValueError("k=%d must be between 1 and min(A.shape)=%d" % (k, min(n, m)))
    primme_svds_params_set(PP, "numSvals", k)

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

    def prevec(X, mode):
        if mode == primme_svds_op_AtA and precAHA is not None:
            return precAHA.matmat(X)
        elif mode == primme_svds_op_AAt and precAAH is not None:
            return precAAH.matmat(X) 
        elif mode == primme_svds_op_augmented and precAug is not None:
            return precAug.matmat(X) 
        return X

    if precAHA is not None or precAAH is not None or precAug is not None:
        primme_svds_params_set(PP, "precondition", 1)
        primme_svds_params_set(PP, "preconditioner", prevec)
    else:
        primme_svds_params_set(PP, "precondition", 0)

    hist = {"numMatvecs": [], "elapsedTime": [], "nconv": [],
            "sval": [], "resNorm": []}

    def mon(basisSvals, basisFlags, iblock, basisNorms, numConverged,
            lockedSvals, lockedFlags, lockedNorms, inner_its, LSRes,
            event, stage):
        if event == 0 and iblock and len(iblock)>0: # event iteration
            hist["numMatvecs"].append(primme_svds_params_get(PP, 'stats_numMatvecs'))
            hist["elapsedTime"].append(primme_svds_params_get(PP, 'stats_elapsedTime'))
            hist["nconv"].append(numConverged)
            hist["sval"].append(basisSvals[iblock[0]])
            hist["resNorm"].append(basisNorms[0])

    if return_history:
        primme_svds_params_set(PP, 'monitor', mon)

    cdef double sigma_c
    if which == 'LM':
        primme_svds_params_set(PP, "target", "primme_svds_largest")
    elif which == 'SM':
        primme_svds_params_set(PP, "target", "primme_svds_smallest")
    else:
        try:
            sigma_c = float(which)
        except:
            raise ValueError("which must be either 'LM', 'SM' or a number.")
        primme_svds_params_set(PP, "numTargetShifts", 1)
        primme_svds_params_set_doubles(pp, "targetShifts", &sigma_c)
        primme_svds_params_set(PP, "target", "primme_svds_closest_abs")

    primme_svds_params_set(PP, "eps", tol)

    if ncv:
        primme_svds_params_set(PP, "maxBasisSize", ncv)

    if maxiter:
        # NOTE: every eigensolver iteration spend two matvecs*blockSize
        primme_svds_params_set(PP, "maxMatvecs", maxiter*(maxBlockSize if maxBlockSize else 1)/2)

    if maxBlockSize:
        primme_svds_params_set(PP, "maxBlockSize", maxBlockSize)

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

    orthou0, orthov0 = check_pair(orthou0, orthov0, "orthov0 or orthou0")

    cdef int numOrthoConst = 0
    if orthou0 is not None:
        numOrthoConst = min(orthou0.shape[1], min(m,n))
        primme_svds_params_set(PP, "numOrthoConst", numOrthoConst)

    # Set other parameters
    for dk, dv in kargs.items():
      try:
        primme_svds_params_set(PP, dk, dv)
      except:
        raise ValueError("Invalid option '%s' with value '%s'" % (dk, dv))

    if A.dtype.kind in frozenset(["b", "i", "u"]) or A.dtype.type is np.double:
        dtype = np.dtype("d")
    else:
        dtype = A.dtype

    if dtype.type is np.complex64:
        rtype = np.dtype(np.float32)
        primme_svds_params_set_pointer(pp, "matrixMatvec", <void*>c_svds_matvec_numpy[np.complex64_t])
        primme_svds_params_set_pointer(pp, "applyPreconditioner", <void*>c_svds_precond_numpy[np.complex64_t])
        if return_history:
            primme_svds_params_set_pointer(pp, "monitorFun", <void*>c_svds_monitor[float])
    elif dtype.type is np.float32:
        rtype = np.dtype(np.float32)
        primme_svds_params_set_pointer(pp, "matrixMatvec", <void*>c_svds_matvec_numpy[float])
        primme_svds_params_set_pointer(pp, "applyPreconditioner", <void*>c_svds_precond_numpy[float])
        if return_history:
            primme_svds_params_set_pointer(pp, "monitorFun", <void*>c_svds_monitor[float])
    elif dtype.type is np.float64:
        rtype = np.dtype(np.float64)
        primme_svds_params_set_pointer(pp, "matrixMatvec", <void*>c_svds_matvec_numpy[double])
        primme_svds_params_set_pointer(pp, "applyPreconditioner", <void*>c_svds_precond_numpy[double])
        if return_history:
            primme_svds_params_set_pointer(pp, "monitorFun", <void*>c_svds_monitor[double])
    else:
        rtype = np.dtype(np.float64)
        primme_svds_params_set_pointer(pp, "matrixMatvec", <void*>c_svds_matvec_numpy[np.complex128_t])
        primme_svds_params_set_pointer(pp, "applyPreconditioner", <void*>c_svds_precond_numpy[np.complex128_t])
        if return_history:
            primme_svds_params_set_pointer(pp, "monitorFun", <void*>c_svds_monitor[double])

    svals = np.zeros(k, rtype)
    svecs = np.empty(((m+n)*(numOrthoConst+k),), dtype)
    norms = np.zeros(k, rtype)

    u0, v0 = check_pair(u0, v0, "v0 or u0")
    
    # Set method
    cdef int method_int = -1, methodStage1_int = -1, methodStage2_int = -1;
    if method is not None or methodStage1 is not None or methodStage2 is not None:
        primme_svds_constant_info(method if method is not None else "primme_svds_default", &method_int)
        if method_int < 0:
            raise ValueError('Not valid "method": %s' % method)
        primme_constant_info(methodStage1 if methodStage1 is not None else "PRIMME_DEFAULT_METHOD", &methodStage1_int)
        if methodStage1_int < 0:
            raise ValueError('Not valid "methodStage1": %s' % methodStage1)
        primme_constant_info(methodStage2 if methodStage2 is not None else "PRIMME_DEFAULT_METHOD", &methodStage2_int)
        if methodStage2_int < 0:
            raise ValueError('Not valid "methodStage2": %s' % methodStage2)
        primme_svds_set_method(method_int, methodStage1_int, methodStage2_int, pp)

    cdef np.ndarray[double, ndim=1, mode="c"] svals_d, norms_d
    cdef np.ndarray[float, ndim=1, mode="c"] svals_f, norms_f
    cdef np.ndarray[float, ndim=1, mode="c"] svecs_f
    cdef np.ndarray[double, ndim=1, mode="c"] svecs_d
    cdef np.ndarray[np.complex64_t, ndim=1, mode="c"] svecs_c
    cdef np.ndarray[np.complex128_t, ndim=1, mode="c"] svecs_z

    if dtype.type is np.complex64:
        svals_f, norms_f, svecs_c = svals, norms, svecs
    elif dtype.type is np.float32:
        svals_f, norms_f, svecs_f = svals, norms, svecs
    elif dtype.type is np.float64:
        svals_d, norms_d, svecs_d = svals, norms, svecs
    else:
        svals_d, norms_d, svecs_z = svals, norms, svecs

    cdef int initSize = 0
    if v0 is not None:
        initSize = min(v0.shape[1], k)
        primme_svds_params_set(PP, "initSize", initSize)
        if dtype.type is np.complex64:
            (<np.complex64_t[:m:1, :initSize]>&svecs_c[m*numOrthoConst])[:,:] = u0[:, 0:initSize]
            (<np.complex64_t[:n:1, :initSize]>&svecs_c[m*(numOrthoConst+initSize)+n*numOrthoConst])[:,:] = v0[:, 0:initSize]
        elif dtype.type is np.float32:
            (<float[:m:1, :initSize]>&svecs_f[m*numOrthoConst])[:,:] = u0[:, 0:initSize]
            (<float[:n:1, :initSize]>&svecs_f[m*(numOrthoConst+initSize)+n*numOrthoConst])[:,:] = v0[:, 0:initSize]
        elif dtype.type is np.float64:
            (<double[:m:1, :initSize]>&svecs_d[m*numOrthoConst])[:,:] = u0[:, 0:initSize]
            (<double[:n:1, :initSize]>&svecs_d[m*(numOrthoConst+initSize)+n*numOrthoConst])[:,:] = v0[:, 0:initSize]
        else:
            (<np.complex128_t[:m:1, :initSize]>&svecs_z[m*numOrthoConst])[:,:] = u0[:, 0:initSize]
            (<np.complex128_t[:n:1, :initSize]>&svecs_z[m*(numOrthoConst+initSize)+n*numOrthoConst])[:,:] = v0[:, 0:initSize]

    if orthou0 is not None:
        if dtype.type is np.complex64:
            (<np.complex64_t[:m:1, :numOrthoConst]>&svecs_c[0])[:,:] = orthou0[:, 0:numOrthoConst]
            (<np.complex64_t[:n:1, :numOrthoConst]>&svecs_c[m*(numOrthoConst+initSize)])[:,:] = orthov0[:, 0:numOrthoConst]
        elif dtype.type is np.float32:
            (<float[:m:1, :numOrthoConst]>&svecs_f[0])[:,:] = orthou0[:, 0:numOrthoConst]
            (<float[:n:1, :numOrthoConst]>&svecs_f[m*(numOrthoConst+initSize)])[:,:] = orthov0[:, 0:numOrthoConst]
        elif dtype.type is np.float64:
            (<double[:m:1, :numOrthoConst]>&svecs_d[0])[:,:] = orthou0[:, 0:numOrthoConst]
            (<double[:n:1, :numOrthoConst]>&svecs_d[m*(numOrthoConst+initSize)])[:,:] = orthov0[:, 0:numOrthoConst]
        else:
            (<np.complex128_t[:m:1, :numOrthoConst]>&svecs_z[0])[:,:] = orthou0[:, 0:numOrthoConst]
            (<np.complex128_t[:n:1, :numOrthoConst]>&svecs_z[m*(numOrthoConst+initSize)])[:,:] = orthov0[:, 0:numOrthoConst]

    if dtype.type is np.complex64:
        err = cprimme_svds(&svals_f[0], &svecs_c[0], &norms_f[0], pp)
    elif dtype.type is np.float32:
        err = sprimme_svds(&svals_f[0], &svecs_f[0], &norms_f[0], pp)
    elif dtype.type is np.float64:
        err = dprimme_svds(&svals_d[0], &svecs_d[0], &norms_d[0], pp)
    else:
        err = zprimme_svds(&svals_d[0], &svecs_z[0], &norms_d[0], pp)

    if err != 0:
        raise PrimmeSvdsError(err)

    if return_stats:
        stats = dict((f, primme_svds_params_get(PP, 'stats_' + f)) for f in [
            "numOuterIterations", "numRestarts", "numMatvecs",
            "numPreconds", "elapsedTime"])
        stats["rnorms"] = norms
        if return_history:
            stats["hist"] = hist
 
    initSize = primme_svds_params_get(PP, "initSize")
    svals = svals[0:initSize]
    if not return_singular_vectors:
        return svals if not return_stats else (svals, stats)

    numOrthoConst = primme_svds_params_get(PP, "numOrthoConst")
    norms = norms[0:initSize]
    if dtype.type is np.complex64:
        svecsl = np.asarray(<np.complex64_t[:m:1, :initSize]>&svecs_c[m*numOrthoConst])
        svecsr = np.asarray(<np.complex64_t[:n:1, :initSize]>&svecs_c[m*(numOrthoConst+initSize)+n*numOrthoConst])
    elif dtype.type is np.float32:
        svecsl = np.asarray(<float[:m:1, :initSize]>&svecs_f[m*numOrthoConst])
        svecsr = np.asarray(<float[:n:1, :initSize]>&svecs_f[m*(numOrthoConst+initSize)+n*numOrthoConst])
    elif dtype.type is np.float64:
        svecsl = np.asarray(<double[:m:1, :initSize]>&svecs_d[m*numOrthoConst])
        svecsr = np.asarray(<double[:n:1, :initSize]>&svecs_d[m*(numOrthoConst+initSize)+n*numOrthoConst])
    else:
        svecsl = np.asarray(<np.complex128_t[:m:1, :initSize]>&svecs_z[m*numOrthoConst])
        svecsr = np.asarray(<np.complex128_t[:n:1, :initSize]>&svecs_z[m*(numOrthoConst+initSize)+n*numOrthoConst])

    # Make copies and transpose conjugate svecsr
    svecsl = svecsl.copy()
    svecsr = svecsr.T.conj().copy()

    if not return_stats:
        return svecsl, svals, svecsr
    else:
        return svecsl, svals, svecsr, stats

#
# PRIMME Errors
#


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



