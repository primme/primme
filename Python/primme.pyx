# cython: language_level=2, c_string_type=bytes, c_string_encoding=ascii, embedsignature=True

import traceback
import numpy as np
cimport numpy as np
import scipy.sparse.linalg
cimport cython
from cython cimport view

try:
    import cupy, cupyx
    __with_cupy = True
except Exception as e:
    cupy_exception = e
    __with_cupy = False

try:
    from inspect import signature
    def __getnumargs(f):
            return len(signature(f).parameters)
except:
    from inspect import getargspec
    def __getnumargs(f):
            return len(getargspec(f).args)

def __get_method(method):
    if method is None: return PRIMME_DEFAULT_METHOD
    cdef int method_int = -1;
    method0 = method
    if not method0.startswith('PRIMME_'):
        method0 = 'PRIMME_' + method0;
    method0 = bytesp23(method0, 'ASCII')
    primme_constant_info(<const char *>method0, &method_int)
    if method_int < 0:
        raise ValueError('Not valid "method": %s' % method)
    return method_int

IF BUILD_WITH_MAGMA: 
    cdef extern from "magma_v2.h":
        struct magma_queue:
            pass
        ctypedef magma_queue* magma_queue_t
        int magma_init()
        int magma_finalize()
        void magma_queue_create(int device, magma_queue_t *queue)
        void magma_queue_destroy(magma_queue_t queue)

try:
    from builtins import bytes as bytesp23 # bytes compatibility Py2/3
except Exception as e:
    raise Exception('package future is not installed') from e

ctypedef fused numerics:
    float
    double
    np.complex64_t
    np.complex128_t

ctypedef fused numerics_real:
    float
    double

# Initialize MAGMA
# NOTE: see documentation for set_device

IF BUILD_WITH_MAGMA: 
    __device_id = -1
    cdef magma_queue_t __queue = NULL
    
    def set_device(device):
        """
        Set the GPU device that the functions eigsh and svds will use in case the
        input matrix is on GPU.
    
        WARNING: setting a new device may invalidate previous CUDA handlers.
    
        Parameters
        ----------
        device : int, positive or zero
            GPU device index
        """
    
        device = int(device)
        if device < 0:
            raise ValueError("device should be >= 0")
    
        # Avoid creating new queues
        global __device_id
        if device == __device_id:
            return 
    
        if __queue is not NULL:
            magma_queue_destroy(__queue)
        magma_queue_create(device, &__queue)
        __device_id = device
    
    if __with_cupy:
        if magma_init() != 0:
            raise RuntimeError('MAGMA initialization failed')
        set_device(cupy.cuda.runtime.getDevice())

# Exception captured in user-defined functions
cdef Exception __user_function_exception = None

cdef np.dtype get_np_type(numerics *p):
    if numerics == double: return np.dtype(np.double)
    elif numerics == float: return np.dtype(np.float32)
    elif numerics == np.complex64_t: return np.dtype(np.complex64)
    elif numerics == np.complex128_t: return np.dtype(np.complex128)

cdef np.dtype get_np_type0(numerics *p):
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
    else:
        return None

def __get_real_dtype(dtype):
    if dtype.type is np.complex64 or dtype.type is np.float32:
        return np.dtype(np.float32)
    else:
        return np.dtype(np.float64)


cdef extern from "../include/primme.h":
    struct primme_params:
        pass
    struct PRIMME_COMPLEX_FLOAT:
        pass
    struct PRIMME_COMPLEX_DOUBLE:
        pass
    ctypedef enum primme_preset_method:
        PRIMME_DEFAULT_METHOD # We only use this value
    ctypedef enum primme_type:
        primme_int, primme_double, primme_pointer # We only use these values
    ctypedef enum primme_params_label: PRIMME_invalid_label
    ctypedef int primme_event
    int sprimme(float *evals, void *evecs, float *resNorms, primme_params *primme)
    int cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, primme_params *primme)
    int dprimme(double *evals, void *evecs, double *resNorms, primme_params *primme)
    int zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, primme_params *primme)
    int magma_sprimme(float *evals, void *evecs, float *resNorms, primme_params *primme)
    int magma_cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, primme_params *primme)
    int magma_dprimme(double *evals, void *evecs, double *resNorms, primme_params *primme)
    int magma_zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, primme_params *primme)
    void primme_display_params(primme_params primme)
    primme_params* primme_params_create()
    int primme_params_destroy(primme_params *primme)
    void primme_initialize(primme_params *primme)
    int  primme_set_method(primme_preset_method method, primme_params *params)
    void primme_free(primme_params *primme)
    int primme_get_member(primme_params *primme, primme_params_label label, void *value)
    int primme_set_member(primme_params *primme, primme_params_label label, void *value)
    int primme_member_info(primme_params_label *label, const char** label_name, primme_type *t, int *arity)
    int primme_constant_info(const char* label_name, int *value)

# Block size passed in a callback such as matrixMatvec, massMatrixMatvec, applyPreconditioner
cdef int __current_block_size  = -1
# Current primme_params running
__current_primme_params = None

cdef class PrimmeParams:
    cdef primme_params *pp
    cdef int _ownpp
    def __cinit__(self, create=True):
        if create:	
            self.pp = primme_params_create()
            if self.pp is NULL:
                raise MemoryError()
            IF BUILD_WITH_MAGMA: 
                primme_params_set_pointer(self.pp, "queue", <void*>&__queue)
            self._ownpp = 1
        else:
            self.pp = NULL
            self._ownpp = 0

    @staticmethod
    cdef from_ptr(void* pp):
        PP = PrimmeParams(False)
        PP.pp = <primme_params*>pp
        PP._ownpp = 0
        return PP

    def __dealloc__(self):
        if self.pp is not NULL and self._ownpp == 1:
            primme_params_destroy(self.pp)
            self.pp = NULL
            self._ownpp = 0
    
def __primme_params_get(PrimmeParams pp_, field_):
    field_ = bytesp23(field_, 'ASCII')
    cdef primme_params *primme = <primme_params*>(pp_.pp)
    cdef const char* field = <const char *>field_
    cdef primme_params_label l = PRIMME_invalid_label
    cdef primme_type t
    cdef int arity, r
    cdef void *v_pvoid

    r = primme_member_info(&l, <const char **>&field, &t, &arity)
    if field_ == bytesp23(b'ShiftsForPreconditioner'):
        if r != 0:
            raise ValueError("Invalid field '%s'" % field_)
        global __current_block_size
        if __current_block_size < 0: raise ValueError('Getting ShiftsForPreconditioner from an invalid callback')
        if __current_block_size == 0: return []
        primme_get_member(primme, l, &v_pvoid)
        return <double[:__current_block_size]> <double *>v_pvoid

    if r != 0 or arity != 1:
        raise ValueError("Invalid field '%s'" % field_)
    cdef np.int64_t v_int
    cdef double v_double
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
    cdef primme_params_label l = PRIMME_invalid_label
    cdef primme_type t
    cdef int arity, r
    cdef void *v_pvoid
    try:
        r = primme_member_info(&l, <const char **>&field, &t, &arity)
        assert r == 0 and arity == 1 and t == primme_pointer, "Invalid field '%s'" % <bytes>field
        r = primme_get_member(primme, l, &v_pvoid)
        assert r == 0, "Invalid field '%s'" % <bytes>field
        if v_pvoid is NULL: return None
        return <object>v_pvoid
    except:
        return None

cdef np.int64_t primme_params_get_int(primme_params *primme, cython.p_char field):
    cdef primme_params_label l = PRIMME_invalid_label
    cdef primme_type t
    cdef int arity, r
    cdef np.int64_t v_int
    try:
        r = primme_member_info(&l, <const char **>&field, &t, &arity)
        assert r == 0 and arity == 1 and t == primme_int, "Invalid field '%s'" % <bytes>field
        r = primme_get_member(primme, l, &v_int)
        assert r == 0, "Invalid field '%s'" % <bytes>field
        return v_int
    except:
        return -1

def __primme_params_set(PrimmeParams pp_, field_, value):
    field_ = bytesp23(field_, 'ASCII')
    cdef primme_params *primme = <primme_params*>(pp_.pp)
    cdef const char* field = <const char*>field_
    cdef primme_params_label l = PRIMME_invalid_label
    cdef primme_type t
    cdef int arity, r
    r = primme_member_info(&l, <const char **>&field, &t, &arity)
    if r != 0  or arity != 1:
        raise ValueError("Invalid field '%s'" % field_)
    cdef np.int64_t v_int
    cdef double v_double
    cdef int i
    if t == primme_pointer:
        r = primme_set_member(primme, l, <void*>value)
        if r != 0: raise Exception("Something went wrong setting the field '%s'" % field_)
    elif t == primme_int:
        if isinstance(value, (bytesp23,str)):
            value = bytesp23(value, 'ASCII')
            r = primme_constant_info(<const char*>value, &i)
            if r != 0: raise ValueError("Invalid value '%s' for field '%s'" % (value, field_))
            value = i
        v_int = value
        r = primme_set_member(primme, l, &v_int)
        if r != 0: raise Exception("Something went wrong setting the field '%s'" % field_)
    elif t == primme_double:
        v_double = value
        r = primme_set_member(primme, l, &v_double)
        if r != 0: raise Exception("Something went wrong setting the field '%s'" % field_)
    else:
        raise ValueError("Not supported type for member '%s'" % field_)
   
cdef void primme_params_set_pointer(primme_params *primme, cython.p_char field, void* value) except *:
    cdef primme_params_label l = PRIMME_invalid_label
    cdef primme_type t
    cdef int arity, r
    r = primme_member_info(&l, <const char **>&field, &t, &arity)
    assert(r == 0 and arity == 1 and t == primme_pointer, "Invalid field '%s'" % <bytes>field)
    r = primme_set_member(primme, l, value)
    assert(r == 0, "Invalid field '%s'" % <bytes>field)

cdef void primme_params_set_doubles(primme_params *primme, cython.p_char field, double *value) except *:
    cdef primme_params_label l = PRIMME_invalid_label
    cdef primme_type t
    cdef int arity, r
    r = primme_member_info(&l, <const char **>&field, &t, &arity)
    assert(r == 0 and arity == 0 and t == primme_double, "Invalid field '%s'" % <bytes>field)
    r = primme_set_member(primme, l, value)
    assert(r == 0, "Invalid field '%s'" % <bytes>field)


cdef void c_matvec_gen_numpy(cython.p_char operator, numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
    global __current_block_size
    if blockSize[0] <= 0:
        ierr[0] = 0
        __current_block_size = -1
        return
    ierr[0] = 1
    cdef object matvec = primme_params_get_object(primme, operator)
    n = primme_params_get_int(primme, "nLocal")
    global __user_function_exception
    try:
        __current_block_size = blockSize[0]
        matvec((<numerics[:ldx[0]:1,:blockSize[0]]>x)[:n,:], (<numerics[:ldy[0]:1,:blockSize[0]]>y)[:n,:])
        ierr[0] = 0
        __current_block_size = -1
    except Exception as e:
        __user_function_exception = e
        __current_block_size = -1

cdef void c_matvec_numpy(numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
    c_matvec_gen_numpy("matrix", x, ldx, y, ldy, blockSize, primme, ierr)

cdef void c_massmatvec_numpy(numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
    c_matvec_gen_numpy("massMatrix", x, ldx, y, ldy, blockSize, primme, ierr)

cdef void c_precond_numpy(numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
    c_matvec_gen_numpy("preconditioner", x, ldx, y, ldy, blockSize, primme, ierr)

cdef void c_monitor(numerics_real *basisEvals, int *basisSize, int *basisFlags, int *iblock, int *blockSize, numerics_real *basisNorms, int *numConverged, numerics_real *lockedEvals, int *numLocked, int *lockedFlags, numerics_real *lockedNorms, int *inner_its, numerics_real *LSRes, const char *msg, double *time, primme_event *event, primme_params *primme, int *ierr):
    ierr[0] = 1
    cdef object monitor = primme_params_get_object(primme, 'monitor')
    if monitor is None: return
    cdef int bs = basisSize[0] if basisSize is not NULL else 0
    cdef int blks = blockSize[0] if blockSize is not NULL else 0
    cdef int nLocked = numLocked[0] if numLocked is not NULL else 0
    global __user_function_exception
    try:
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
            event[0] if event is not NULL else None)
        ierr[0] = 0
    except Exception as e:
        __user_function_exception = e

cdef void c_convtest_numpy(double *eval, numerics *evec, double *resNorm, int *isconv, primme_params *primme, int *ierr):
    ierr[0] = 1
    cdef object convtest = primme_params_get_object(primme, 'convtest')
    if convtest is None: return
    global __user_function_exception
    try:
        n = primme_params_get_int(primme, "nLocal")
        isconv[0] = 1 if convtest(eval[0] if eval is not NULL else None,
            <numerics[:n]>evec if evec is not NULL else None,
            resNorm[0] if resNorm is not NULL else None) else 0
        ierr[0] = 0
    except Exception as e:
        __user_function_exception = e
 
IF BUILD_WITH_MAGMA: 
    cdef void c_matvec_gen_cupy(cython.p_char operator, void *x, np.int64_t *ldx, void *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
        if blockSize[0] <= 0:
            ierr[0] = 0
            return
        ierr[0] = 1
        cdef object matvec
        global __user_function_exception
        try:
            matvec = primme_params_get_object(primme, operator)
            matvec(matvec.__toObj(<size_t>x, blockSize[0], ldx[0]), matvec.__toObj(<size_t>y, blockSize[0], ldy[0]))
            ierr[0] = 0
        except Exception as e:
            __user_function_exception = e
    
    cdef void c_matvec_cupy(void *x, np.int64_t *ldx, void *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
        c_matvec_gen_cupy("matrix", x, ldx, y, ldy, blockSize, primme, ierr)
    
    cdef void c_massmatvec_cupy(void *x, np.int64_t *ldx, void *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
        c_matvec_gen_cupy("massMatrix", x, ldx, y, ldy, blockSize, primme, ierr)
    
    cdef void c_precond_cupy(void *x, np.int64_t *ldx, void *y, np.int64_t *ldy, int *blockSize, primme_params *primme, int *ierr):
        c_matvec_gen_cupy("preconditioner", x, ldx, y, ldy, blockSize, primme, ierr)
    
    cdef void c_convtest_cupy(double *eval, void *evec, double *resNorm, int *isconv, primme_params *primme, int *ierr):
        ierr[0] = 1
        cdef object convtest = primme_params_get_object(primme, 'convtest')
        cdef object matvec = primme_params_get_object(primme, 'matrix')
        global __user_function_exception
        try:
            n = primme_params_get_int(primme, "nLocal")
            isconv[0] = 1 if convtest(eval[0] if eval is not NULL else None,
                matvec.__toObj(<size_t>evec, 1, n) if evec is not NULL else None,
                resNorm[0] if resNorm is not NULL else None) else 0
            ierr[0] = 0
        except Exception as e:
            __user_function_exception = e

def __to_eigsh_matvec(A, driver, argname):
    if isinstance(A, (np.ndarray,scipy.sparse.spmatrix,scipy.sparse.linalg.LinearOperator)) or (__with_cupy and isinstance(A, (cupy.ndarray,cupyx.scipy.sparse.spmatrix))):
        def matvec(x, y): y[:,:] = A.dot(x)
    elif __getnumargs(A) == 1:
        if driver is None: raise Exception('Please set `driver` when the valude of `%s` is a function' % argname)
        def matvec(x, y): y[:,:] = A(x)
    elif __getnumargs(A) == 2:
        if driver is None: raise Exception('Please set `driver` when the valude of `%s` is a function' % argname)
        matvec = A
    else:
        raise Exception('Unsupported value of `%s` as matrix or function' % argname)

    n = A.shape[0]
    dtype = A.dtype
    itemsize = dtype.itemsize
    if isinstance(A, (np.ndarray,scipy.sparse.spmatrix,scipy.sparse.linalg.LinearOperator)) or driver == 'numpy':
        if driver != 'numpy' and driver is not None:
            raise Exception('Incompatible value of `driver` with `%s`' % argname)
        def toObj(ptr, blockSize, ld):
            return np.ndarray((n,blockSize), dtype, ptr, strides=(itemsize,itemsize*ld), order='F')
        matvec.__toObj = toObj
    elif __with_cupy and (isinstance(A, (cupy.ndarray,cupyx.scipy.sparse.spmatrix)) or driver == 'cupy'):
        if driver != 'cupy' and driver is not None:
            raise Exception('Incompatible value of `driver` with `%s`' % argname)
        def toObj(ptr, blockSize, ld):
            return cupy.ndarray((n,blockSize), dtype, cupy.cuda.MemoryPointer(cupy.cuda.UnownedMemory(ptr,0,None),0), strides=(itemsize,itemsize*ld), order='F')
        matvec.__toObj = toObj

    return matvec

def eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None,
          ncv=None, maxiter=None, tol=None, reltol=None, return_eigenvectors=True,
          Minv=None, OPinv=None, mode='normal', lock=None, driver=None,
          return_stats=False, maxBlockSize=0, minRestartSize=0,
          maxPrevRetain=0, method=None, return_unconverged=False,
          return_history=False, convtest=None, raise_for_unconverged=True, **kargs):
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
    A : matrix (numpy.ndarray, scipy.sparse.spmatrix, scipy.sparse.linalg.LinearOperator, cupy.ndarray,cupyx.scipy.sparse.spmatrix) or function with one argument
        the operation A * x, where A is a real symmetric matrix or complex
        Hermitian. If it is a function, A(x) returns A*x and A.shape and
        A.dtype give the shape and the dtype of the matrix.
    k : int, optional
        The number of eigenvalues and eigenvectors to be computed. Must be
        1 <= k < min(A.shape).
    M : matrix (numpy.ndarray, scipy.sparse.spmatrix, scipy.sparse.linalg.LinearOperator, cupy.ndarray,cupyx.scipy.sparse.spmatrix) or function with one argument, optional
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
    which : str ['LM' | 'SM' | 'LA' | 'SA' | number]
        Which `k` eigenvectors and eigenvalues to find:

            'LM' : Largest in magnitude eigenvalues; the farthest from sigma

            'SM' : Smallest in magnitude eigenvalues; the closest to sigma

            'LA' : Largest algebraic eigenvalues

            'SA' : Smallest algebraic eigenvalues

            'CLT' : closest but left to sigma

            'CGT' : closest but greater than sigma

            number : the closest to which

        When sigma == None, 'LM', 'SM', 'CLT', and 'CGT' treat sigma as zero. 
    maxiter : int, optional
        Maximum number of iterations.
    tol : float
        Tolerance for eigenpairs (stopping criterion). The default value is sqrt of machine precision.

        An eigenpair ``(lamba,v)`` is marked as converged when ||A*v - lambda*B*v|| < max(|eig(A,B)|)*tol.

        An error is raised if convtest is also provided.
    reltol : float
        Relative tolerance for eigenpairs (stopping criterion). The default value is None.

        An eigenpair ``(lamba,v)`` is marked as converged when ||A*v - lambda*B*v|| < max(|eig(A,B)|)*tol (if given)
        and ||A*v - lambda*B*v|| < |lambda|*reltol.

        An error is raised if convtest is also provided.
    Minv : (not supported yet)
        The inverse of M in the generalized eigenproblem.
    OPinv : matrix, scipy.sparse.linalg.interface.LinearOperator, pycuda.sparse.operator.OperatorBase, or function, optional
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
    driver : string ['numpy' | 'cupy' ]
        Kind of object passed to A, M and OPinv, and returned as evecs. If 'numpy' or 'cupy' is given, then
        the passed input and the expected returned object is from class `ndarray' in
        numpy or cupy respectively.
        
        If no value is given, and `A' is scipy.sparse.spmatrix or cupyx.scipy.sparse.spmatrix,
        then it is set to 'numpy' and 'cupy' respectively.
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
    return_unconverged: bool, optional
        If True, return all requested eigenvalues and vectors regardless of
        being marked as converged. The default is False.
    convtest : callable
        User-defined function to mark an approximate eigenpair as converged.

        The function is called as convtest(eval, evec, resNorm) and returns
        True if the eigenpair with value `eval`, vector `evec` and residual
        norm `resNorm` is considered converged.
    raise_for_unconverged: bool, optional
        If True and return_unconverged is False, raise an exception when not
        all requested eigenvalues and vectors converged. The default is True.

    return_stats : bool, optional
        If True, the function returns extra information (see stats in Returns).
    return_history: bool, optional
        If True, the function returns performance information at every iteration
        (see hist in Returns).

    Returns
    -------
    w : array
        Array of k eigenvalues ordered to best satisfy "which".
    v : array, optional (if return_eigenvectors)
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
    >>> evals, evecs = primme.eigsh(A, 3, tol=1e-6, which=50.1)
    >>> evals # the three closest eigenvalues to 50.1
    array([50.,  51.,  49.])
    >>> M = scipy.sparse.spdiags(np.asarray(range(99,-1,-1)), [0], 100, 100)
    >>> # the smallest eigenvalues of the eigenproblem (A,M)
    >>> evals, evecs = primme.eigsh(A, 3, M=M, tol=1e-6, which='SA')
    >>> evals # doctest: +SKIP
    array([1.0035e-07, 1.0204e-02, 2.0618e-02])

    >>> # Giving the matvec as a function
    >>> import primme, scipy.sparse, numpy as np
    >>> Adiag = np.arange(0, 100).reshape((100,1))
    >>> def Amatmat(x):
    ...    if len(x.shape) == 1: x = x.reshape((100,1))
    ...    return Adiag * x   # equivalent to diag(Adiag).dot(x)
    ...
    >>> A = scipy.sparse.linalg.LinearOperator((100,100), matvec=Amatmat, matmat=Amatmat)
    >>> evals, evecs = primme.eigsh(A, 3, tol=1e-6, which='LA')
    >>> evals
    array([99., 98., 97.])

    >>> import primme, scipy.sparse
    >>> import pycuda.sparse.packeted
    >>> A = scipy.sparse.spdiags(range(100), [0], 100, 100) # sparse diag. matrix
    >>> Agpu = pycuda.sparse.packeted.PacketedSpMV(A, True, A.dtype)
    >>> evals, evecs = primme.eigsh(Agpu, 3, tol=1e-6, which='LA')
    >>> evals # the three largest eigenvalues of A
"""

    PP = PrimmeParams()
    cdef primme_params *pp = PP.pp
 
    shape = A.shape
        
    if len(shape) != 2 or shape[0] != shape[1]:
        raise ValueError('A: expected square matrix (shape=%s)' % shape)

    if driver is None:
        if isinstance(A, (np.ndarray,scipy.sparse.spmatrix,scipy.sparse.linalg.LinearOperator)):
            driver = 'numpy'
        elif __with_cupy and isinstance(A, (cupy.ndarray,cupyx.scipy.sparse.spmatrix)):
            driver = 'cupy'
        else:
            raise Exception('Please set `driver` when giving a function as `A`')

    Af = __to_eigsh_matvec(A, driver, 'A')
    __primme_params_set(PP, "matrix", Af)

    n = shape[0]
    __primme_params_set(PP, "n", n)

    if M is not None:
        if len(M.shape) != 2 or shape[0] != M.shape[0]:
            raise ValueError('M: expected square matrix (shape=%s)' % A.shape)
        Mf = __to_eigsh_matvec(M, driver, 'M')
        __primme_params_set(PP, "massMatrix", Mf)

    if k <= 0 or k > n:
        raise ValueError("k=%d must be between 1 and %d, the order of the "
                         "square input matrix." % (k, n))
    __primme_params_set(PP, "numEvals", k)

    if which == 'LM':
        __primme_params_set(PP, "target", "primme_largest_abs")
        if sigma is None:
            sigma = 0.0
    elif which == 'LA':
        __primme_params_set(PP, "target", "primme_largest")
        sigma = None
    elif which == 'SA':
        __primme_params_set(PP, "target", "primme_smallest")
        sigma = None
    elif which == 'SM':
        __primme_params_set(PP, "target", "primme_closest_abs")
        if sigma is None:
            sigma = 0.0
    elif which == 'CLT':
        __primme_params_set(PP, "target", "primme_closest_leq")
        if sigma is None:
            sigma = 0.0
    elif which == 'CGT':
        __primme_params_set(PP, "target", "primme_closest_geq")
        if sigma is None:
            sigma = 0.0
    else:
        try:
            sigma0 = float(which)
        except:
            raise ValueError("which='%s' not supported. It should be 'LM', 'LA', 'SA', 'SM', 'CLT', 'CGT' or a number" % which)
        if sigma is not None:
            raise ValueError("Giving a numeric value in `which`, and also giving `sigma`. Set only one of those.")
        sigma = sigma0
        __primme_params_set(PP, "target", "primme_closest_abs")

    cdef double sigma_c
    if sigma is not None:
        sigma_c = float(sigma)
        __primme_params_set(PP, "numTargetShifts", 1)
        primme_params_set_doubles(pp, "targetShifts", &sigma_c)

    if tol is not None:
        if convtest is not None:
            raise ValueError('Giving tol and convtest is not allowed.')
        __primme_params_set(PP, "eps", tol)

    if reltol is not None:
        if convtest is not None:
            raise ValueError('Giving reltol and convtest is not allowed.')
        def check_convergence(eval, evec, resNorm):
            if resNorm is None: return True
            if tol is not None:
                aNorm = __primme_params_get(PP, "stats_estimateLargestSVal")
                if resNorm > aNorm*tol: return False
            if eval is not None and resNorm > abs(eval)*reltol:
                return False
            return True

        convtest = check_convergence

    if ncv is not None:
        __primme_params_set(PP, "maxBasisSize", ncv)

    if maxiter is not None:
        __primme_params_set(PP, "maxMatvecs", maxiter)

    if OPinv is not None:
        if OPinv.shape[0] != OPinv.shape[1] or OPinv.shape[0] != A.shape[0]:
            raise ValueError('OPinv: expected square matrix with same shape as A (shape=%s)' % (OPinv.shape,n))
        OPinvf = __to_eigsh_matvec(OPinv, driver, 'OPinv')
        __primme_params_set(PP, "correction_precondition", 1)
        __primme_params_set(PP, "preconditioner", OPinvf)
    else:
        __primme_params_set(PP, "correction_precondition", 0)

    numOrthoConst = 0
    if lock is not None:
        if lock.shape[0] != n:
            raise ValueError('lock: expected matrix with the same columns as A (shape=%s)' % (lock.shape,n))
        numOrthoConst = min(lock.shape[1], n)
        __primme_params_set(PP, "numOrthoConst", numOrthoConst)

    if convtest is not None:
        __primme_params_set(PP, "convtest", convtest)

    # Set other parameters
    for dk, dv in kargs.items():
      try:
        __primme_params_set(PP, dk, dv)
      except:
        raise ValueError("Invalid option '%s' with value '%s'" % (dk, dv))

    hist = {"numMatvecs": [], "elapsedTime": [], "nconv": [],
            "eval": [], "resNorm": []}

    def mon(basisEvals, basisFlags, iblock, basisNorms, numConverged,
            lockedEvals, lockedFlags, lockedNorms, inner_its, LSRes, event):
        
        if event == 0 and iblock and len(iblock)>0: # event iteration
            hist["numMatvecs"].append(__primme_params_get(PP, 'stats_numMatvecs'))
            hist["elapsedTime"].append(__primme_params_get(PP, 'stats_elapsedTime'))
            hist["nconv"].append(numConverged)
            hist["eval"].append(basisEvals[iblock[0]])
            hist["resNorm"].append(basisNorms[0])

    if return_history:
        __primme_params_set(PP, 'monitor', mon)

    if A.dtype.kind in frozenset(["b", "i", "u"]) or A.dtype.type is np.double:
        dtype = np.dtype("d")
    else:
        dtype = A.dtype


    if driver == 'numpy':
        if dtype.type is np.complex64:
            primme_params_set_pointer(pp, "matrixMatvec", <void*>c_matvec_numpy[np.complex64_t])
            if M is not None: 
                primme_params_set_pointer(pp, "massMatrixMatvec", <void*>c_massmatvec_numpy[np.complex64_t])
            primme_params_set_pointer(pp, "applyPreconditioner", <void*>c_precond_numpy[np.complex64_t])
            if convtest:
                primme_params_set_pointer(pp, "convTestFun", <void*>c_convtest_numpy[np.complex64_t])
        elif dtype.type is np.float32:
            primme_params_set_pointer(pp, "matrixMatvec", <void*>c_matvec_numpy[float])
            if M is not None: 
                primme_params_set_pointer(pp, "massMatrixMatvec", <void*>c_massmatvec_numpy[float])
            primme_params_set_pointer(pp, "applyPreconditioner", <void*>c_precond_numpy[float])
            if convtest:
                primme_params_set_pointer(pp, "convTestFun", <void*>c_convtest_numpy[float])
        elif dtype.type is np.float64:
            primme_params_set_pointer(pp, "matrixMatvec", <void*>c_matvec_numpy[double])
            if M is not None: 
                primme_params_set_pointer(pp, "massMatrixMatvec", <void*>c_massmatvec_numpy[double])
            primme_params_set_pointer(pp, "applyPreconditioner", <void*>c_precond_numpy[double])
            if convtest:
                primme_params_set_pointer(pp, "convTestFun", <void*>c_convtest_numpy[double])
        else:
            primme_params_set_pointer(pp, "matrixMatvec", <void*>c_matvec_numpy[np.complex128_t])
            if M is not None: 
                primme_params_set_pointer(pp, "massMatrixMatvec", <void*>c_massmatvec_numpy[np.complex128_t])
            primme_params_set_pointer(pp, "applyPreconditioner", <void*>c_precond_numpy[np.complex128_t])
            if convtest:
                primme_params_set_pointer(pp, "convTestFun", <void*>c_convtest_numpy[np.complex128_t])
    elif driver == 'cupy':
        IF BUILD_WITH_MAGMA: 
            primme_params_set_pointer(pp, "matrixMatvec", <void*>c_matvec_cupy)
            if M is not None: 
                primme_params_set_pointer(pp, "massMatrixMatvec", <void*>c_massmatvec_cupy)
            primme_params_set_pointer(pp, "applyPreconditioner", <void*>c_precond_cupy)
            if convtest:
                primme_params_set_pointer(pp, "convTestFun", <void*>c_convtest_cupy)
        ELSE:
            raise Exception('primme is not compiled with cupy support')

    if return_history:
        primme_params_set_pointer(pp, "monitorFun", <void*>c_monitor[float] if dtype.type is np.float32 or dtype.type is np.complex64 else <void*>c_monitor[double])

    cdef double[::1] evals_d, norms_d
    cdef float[::1] evals_s, norms_s
    cdef float[::1, :] evecs_s
    cdef double[::1, :] evecs_d
    cdef np.complex64_t[::1, :] evecs_c
    cdef np.complex128_t[::1, :] evecs_z
    cdef void *evecs_p = <void*>0

    rtype = __get_real_dtype(dtype);
    evals = np.zeros(k, rtype)
    norms = np.zeros(k, rtype)
    if rtype.type is np.float32:
        evals_s, norms_s = evals, norms
    else:
        evals_d, norms_d = evals, norms
    
    if driver == 'numpy':
        evecs = np.zeros((n, numOrthoConst+k), dtype, order='F')
        if dtype.type is np.float64:
            evecs_d = evecs
            evecs_p = <void*>&evecs_d[0,0]
        elif dtype.type is np.float32:
            evecs_s = evecs
            evecs_p = <void*>&evecs_s[0,0]
        elif dtype.type is np.complex64:
            evecs_c = evecs
            evecs_p = <void*>&evecs_c[0,0]
        else:
            evecs_z = evecs
            evecs_p = <void*>&evecs_z[0,0]
    else:
        evecs = cupy.zeros((n, numOrthoConst+k), dtype, order='F')
        evecs_p = <void*>(<size_t>evecs.data.ptr)

    if lock is not None:
        evecs[:, 0:numOrthoConst] = lock[:, 0:numOrthoConst]

    if v0 is not None:
        initSize = min(v0.shape[1], k)
        __primme_params_set(PP, "initSize", initSize)
        evecs[:, numOrthoConst:numOrthoConst+initSize] = v0[:, 0:initSize]

    if maxBlockSize:
        __primme_params_set(PP, "maxBlockSize", maxBlockSize)

    if minRestartSize:
        __primme_params_set(PP, "minRestartSize", minRestartSize)

    if maxPrevRetain:
        __primme_params_set(PP, "restarting_maxPrevRetain", maxPrevRetain)

    cdef int method_int = -1;
    if method is not None:
        method_int = __get_method(method)
        primme_set_method(<primme_preset_method>method_int, pp)

        # Set other parameters (again)
        for dk, dv in kargs.items():
          try:
            __primme_params_set(PP, dk, dv)
          except:
            raise ValueError("Invalid option '%s' with value '%s'" % (dk, dv))


    if __primme_params_get(PP, "printLevel") >= 5:
        primme_display_params(pp[0]);

    global __user_function_exception
    __user_function_exception = None
    global __current_primme_params
    __current_primme_params = PP
    err = 0 # avoid warning
    if driver == 'numpy':
        if dtype.type is np.complex64:
            err = cprimme(&evals_s[0], <PRIMME_COMPLEX_FLOAT*>evecs_p, &norms_s[0], pp)
        elif dtype.type is np.float32:
            err = sprimme(&evals_s[0], <float*>evecs_p, &norms_s[0], pp)
        elif dtype.type is np.float64:
            err = dprimme(&evals_d[0], <double*>evecs_p, &norms_d[0], pp)
        else:
            err = zprimme(&evals_d[0], <PRIMME_COMPLEX_DOUBLE*>evecs_p, &norms_d[0], pp)
    elif driver == 'cupy':
        if dtype.type is np.complex64:
            err = magma_cprimme(&evals_s[0], <PRIMME_COMPLEX_FLOAT*>evecs_p, &norms_s[0], pp)
        elif dtype.type is np.float32:
            err = magma_sprimme(&evals_s[0], <float*>evecs_p, &norms_s[0], pp)
        elif dtype.type is np.float64:
            err = magma_dprimme(&evals_d[0], <double*>evecs_p, &norms_d[0], pp)
        else:
            err = magma_zprimme(&evals_d[0], <PRIMME_COMPLEX_DOUBLE*>evecs_p, &norms_d[0], pp)
    __current_primme_params = None

    if err != 0:
        if __user_function_exception is not None:
            raise PrimmeError(err) from __user_function_exception
        elif err != -3 or (not return_unconverged and raise_for_unconverged):
            raise PrimmeError(err)

    initSize = k if return_unconverged else __primme_params_get(PP, "initSize")

    # Always return eigenvalues
    evals = evals[0:initSize]
    result = [evals]

    if return_eigenvectors:
        evecs = evecs[:, numOrthoConst:numOrthoConst+initSize]
        result.append(evecs)

    if return_stats:
        norms = norms[0:initSize]
        stats = dict((f, __primme_params_get(PP, "stats_" + f)) for f in [
            "numOuterIterations", "numRestarts", "numMatvecs",
            "numPreconds", "elapsedTime", "estimateMinEVal",
            "estimateMaxEVal", "estimateLargestSVal"])
        stats['rnorms'] = norms
        if return_history:
            stats["hist"] = hist
        result.append(stats)

    if len(result) == 1:
        # Don't return tuple if just eigenvalues
        return result[0]

    return tuple(result)

def get_eigsh_param(field):
    """
    Get the value of a PRIMME parameter within a callback function.

    Parameters
    ----------
    field : parameter name
      'n'
      'numEvals'
      'target'
      'locking'
      'initSize'
      'numOrthoConst'
      'maxBasisSize'
      'minRestartSize'
      'maxBlockSize'
      'maxMatvecs'
      'maxOuterIterations'
      'aNorm'
      'BNorm'
      'invBNorm'
      'eps'
      'orth'
      'internalPrecision'
      'printLevel'
      'matrix'
      'massMatrix'
      'preconditioner'
      'initBasisMode'
      'projectionParams_projection'
      'restartingParams_maxPrevRetain'
      'correctionParams_precondition'
      'correctionParams_robustShifts'
      'correctionParams_maxInnerIterations'
      'correctionParams_projectors_LeftQ'
      'correctionParams_projectors_LeftX'
      'correctionParams_projectors_RightQ'
      'correctionParams_projectors_RightX'
      'correctionParams_projectors_SkewQ'
      'correctionParams_projectors_SkewX'
      'correctionParams_convTest'
      'correctionParams_relTolBase'
      'stats_numOuterIterations'
      'stats_numRestarts'
      'stats_numMatvecs'
      'stats_numPreconds'
      'stats_numGlobalSum'
      'stats_volumeGlobalSum'
      'stats_numBroadcast'
      'stats_volumeBroadcast'
      'stats_flopsDense'
      'stats_numOrthoInnerProds'
      'stats_elapsedTime'
      'stats_timeMatvec'
      'stats_timePrecond'
      'stats_timeOrtho'
      'stats_timeGlobalSum'
      'stats_timeBroadcast'
      'stats_timeDense'
      'stats_estimateMinEVal'
      'stats_estimateMaxEVal'
      'stats_estimateLargestSVal'
      'stats_estimateBNorm'
      'stats_estimateInvBNorm'
      'stats_maxConvTol'
      'stats_lockingIssue'
      'dynamicMethodSwitch'
      'convtest'
      'ldevecs'
      'ldOPs'
      'monitor'

    See http://www.cs.wm.edu/~andreas/software/doc/appendix.html#primme-params for a
    description of each parameter.
 
    Examples
    --------
    >>> import primme, scipy.sparse
    >>> A = scipy.sparse.spdiags(range(100), [0], 100, 100) # sparse diag. matrix
    >>> def convtest(eval, evec, rnorm):
    ...   estimateAnorm = primme.get_eigsh_param('stats_estimateLargestSVal')
    ...   return rnorm <= estimateAnorm * 1e-3
    ...
    >>> evals, evecs = primme.eigsh(A, 3, which='LA', convtest=convtest)
 
    >>> import primme, scipy.sparse
    >>> A = scipy.sparse.spdiags(range(100), [0], 100, 100) # sparse diag. matrix
    >>> def P(x):
    ...    # The scipy.sparse.linalg.LinearOperator constructor may call this function giving a vector
    ...    # as input; detect that case and return whatever
    ...    if x.ndim == 1:
    ...       return x / A.diagonal()
    ...    shifts = primme.get_eigsh_param('ShiftsForPreconditioner')
    ...    y = np.copy(x)
    ...    for i in range(x.shape[1]): y[:,i] = x[:,i] / (A.diagonal() - shifts[i])
    ...    return y
    ...
    >>> Pop = scipy.sparse.linalg.LinearOperator(A.shape, matvec=P, matmat=P)
    >>> evals, evecs = primme.eigsh(A, 3, OPinv=Pop, tol=1e-3, which='LA')
    """ 
    if __current_primme_params is None:
        raise RuntimeError("no eigsh running; call this function within a callback function while `eigsh` is running")

    return __primme_params_get(__current_primme_params, field)
 
cdef extern from "../include/primme.h":
    struct primme_svds_params:
        pass
    ctypedef enum primme_svds_preset_method:
        primme_svds_default # We only use this value
    ctypedef enum primme_svds_params_label: PRIMME_SVDS_invalid_label
    ctypedef enum primme_svds_operator:
        primme_svds_op_none,
        primme_svds_op_AtA,
        primme_svds_op_AAt,
        primme_svds_op_augmented
    int sprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds)
    int cprimme_svds(float *svals, PRIMME_COMPLEX_FLOAT *svecs, float *resNorms, primme_svds_params *primme_svds)
    int dprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds)
    int zprimme_svds(double *svals, PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms, primme_svds_params *primme_svds)
    int magma_sprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds)
    int magma_cprimme_svds(float *svals,  PRIMME_COMPLEX_FLOAT *svecs, float *resNorms, primme_svds_params *primme_svds)
    int magma_dprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds)
    int magma_zprimme_svds(double *svals,  PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms, primme_svds_params *primme_svds)
    primme_svds_params* primme_svds_params_create()
    int primme_svds_params_destroy(primme_svds_params *primme_svds)
    void primme_svds_initialize(primme_svds_params *primme_svds)
    void primme_svds_display_params(primme_svds_params primme_svds)
    int primme_svds_set_method(primme_svds_preset_method method, primme_preset_method methodStage1, primme_preset_method methodStage2, primme_svds_params *primme_svds)
    void primme_svds_free(primme_svds_params *primme_svds)
    int primme_svds_get_member(primme_svds_params *primme_svds, primme_svds_params_label label, void *value)
    int primme_svds_set_member(primme_svds_params *primme_svds, primme_svds_params_label label, void *value)
    int primme_svds_member_info(primme_svds_params_label *label, const char** label_name, primme_type *t, int *arity)
    int primme_svds_constant_info(const char* label_name, int *value)


cdef class PrimmeSvdsParams:
    cdef primme_svds_params *pp
    def __cinit__(self):
        self.pp = primme_svds_params_create()
        if self.pp is NULL:
            raise MemoryError()
        IF BUILD_WITH_MAGMA: 
            primme_svds_params_set_pointer(self.pp, "queue", <void*>&__queue)

    def __dealloc__(self):
        if self.pp is not NULL:
            primme_svds_params_destroy(self.pp)
    
def __primme_svds_params_get(PrimmeSvdsParams pp_, field_):
    field_ = bytesp23(field_, 'ASCII')
    cdef primme_svds_params *primme_svds = <primme_svds_params*>(pp_.pp)
    cdef const char* field = <const char *>field_
    cdef primme_svds_params_label l = PRIMME_SVDS_invalid_label
    cdef primme_type t
    cdef int arity, r
    r = primme_svds_member_info(&l, <const char **>&field, &t, &arity)
    if r != 0 or arity != 1:
        raise ValueError("Invalid field '%s'" % field_)
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

cdef object primme_svds_params_get_object(primme_svds_params *primme_svds, cython.p_char field):
    cdef primme_svds_params_label l = PRIMME_SVDS_invalid_label
    cdef primme_type t
    cdef int arity, r
    cdef void *v_pvoid
    try:
        r = primme_svds_member_info(&l, <const char **>&field, &t, &arity)
        assert r == 0 and arity == 1 and t == primme_pointer, "Invalid field '%s'" % <bytes>field
        r = primme_svds_get_member(primme_svds, l, &v_pvoid)
        assert r == 0, "Invalid field '%s'" % <bytes>field
        if v_pvoid is NULL: return None
        return <object>v_pvoid
    except:
        return None

cdef void* primme_svds_params_get_pointer(primme_svds_params *primme_svds, const char* field):
    cdef primme_svds_params_label l = <primme_svds_params_label>0
    cdef primme_type t
    cdef int arity, r
    cdef void *v_pvoid
    try:
        r = primme_svds_member_info(&l, <const char **>&field, &t, &arity)
        assert r == 0 and l >= 0 and l < 1000 and arity == 1 and t == primme_pointer, "Invalid field '%s'" % <bytes>field
        r = primme_svds_get_member(primme_svds, l, &v_pvoid)
        assert r == 0, "Invalid field '%s'" % <bytes>field
        return v_pvoid
    except:
        return NULL

cdef np.int64_t primme_svds_params_get_int(primme_svds_params *primme_svds, cython.p_char field):
    cdef primme_svds_params_label l = PRIMME_SVDS_invalid_label
    cdef primme_type t
    cdef int arity, r
    cdef np.int64_t v_int
    try:
        r = primme_svds_member_info(&l, <const char **>&field, &t, &arity)
        assert r == 0 and arity == 1 and t == primme_int, "Invalid field '%s'" % <bytes>field
        r = primme_svds_get_member(primme_svds, l, &v_int)
        assert r == 0, "Invalid field '%s'" % <bytes>field
        return v_int
    except:
        return -1

def __primme_svds_params_set(PrimmeSvdsParams pp_, field_, value):
    field0 = field_
    fieldStartsWithPrimme = field_.startswith('primme')
    field_ = bytesp23(field_, 'ASCII')
    cdef primme_svds_params *primme_svds = <primme_svds_params*>(pp_.pp)
    cdef const char* field = <const char *>field_
    cdef primme_svds_params_label l = PRIMME_SVDS_invalid_label
    cdef primme_type t
    cdef int arity, r
    r = primme_svds_member_info(&l, <const char **>&field, &t, &arity)
    if r != 0 or arity != 1:
        raise ValueError("Invalid field '%s'" % field_)
    cdef np.int64_t v_int
    cdef double v_double
    cdef int i
    cdef void* primme
    if t == primme_pointer and fieldStartsWithPrimme:
        if not isinstance(value, dict):
            raise Exception("Invalid value for the field '%s': it should be a dictionary" % field_)
        primme = primme_svds_params_get_pointer(primme_svds, field)
        PP = PrimmeParams.from_ptr(primme)
        # Set the parameters in value
        for dk, dv in value.items():
          try:
            __primme_params_set(PP, dk, dv)
          except Exception as e:
            raise ValueError("Invalid option '%s.%s' with value '%s':\n%s" % (field_, dk, dv, e))
    elif t == primme_pointer:
        r = primme_svds_set_member(primme_svds, l, <void*>value)
        if r != 0: raise Exception("Something went wrong setting the field '%s'" % field_)
    elif t == primme_int:
        if isinstance(value, (bytesp23,str)):
            value = bytesp23(value, 'ASCII')
            r = primme_svds_constant_info(<const char*>value, &i)
            if r != 0: raise ValueError("Invalid value '%s' for field '%s'" % (value, field_))
            value = i
        v_int = value
        r = primme_svds_set_member(primme_svds, l, &v_int)
        if r != 0: raise Exception("Something went wrong setting the field '%s'" % field_)
    elif t == primme_double:
        v_double = value
        r = primme_svds_set_member(primme_svds, l, &v_double)
        if r != 0: raise Exception("Something went wrong setting the field '%s'" % field_)
    else:
        raise ValueError("Not supported type for member '%s'" % field_)
   
cdef void primme_svds_params_set_pointer(primme_svds_params *primme_svds, cython.p_char field, void* value) except *:
    cdef primme_svds_params_label l = PRIMME_SVDS_invalid_label
    cdef primme_type t
    cdef int arity, r
    r = primme_svds_member_info(&l, <const char **>&field, &t, &arity)
    assert(r == 0 and arity == 1 and t == primme_pointer)
    r = primme_svds_set_member(primme_svds, l, value)
    assert(r == 0, "Invalid field '%s'" % <bytes>field)

cdef void primme_svds_params_set_doubles(primme_svds_params *primme_svds, cython.p_char field, double *value) except *:
    cdef primme_svds_params_label l = PRIMME_SVDS_invalid_label
    cdef primme_type t
    cdef int arity, r
    r = primme_svds_member_info(&l, <const char **>&field, &t, &arity)
    assert(r == 0 and arity == 0 and t == primme_double)
    r = primme_svds_set_member(primme_svds, l, value)
    assert(r == 0, "Invalid field '%s'" % <bytes>field)


cdef void c_svds_matvec_numpy(numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, int *transpose, primme_svds_params *primme_svds, int *ierr):
    if blockSize[0] <= 0:
        ierr[0] = 0
        return
    ierr[0] = 1
    cdef object matvec = primme_svds_params_get_object(primme_svds, 'matrix')
    m = primme_svds_params_get_int(primme_svds, "mLocal")
    n = primme_svds_params_get_int(primme_svds, "nLocal")
    if transpose[0] != 0: m, n = n, m
    global __user_function_exception
    try:
        matvec((<numerics[:ldx[0]:1,:blockSize[0]]>x)[:n,:], (<numerics[:ldy[0]:1,:blockSize[0]]>y)[:m,:], transpose[0] != 0)
        ierr[0] = 0
    except Exception as e:
        __user_function_exception = e


cdef void c_svds_precond_numpy(numerics *x, np.int64_t *ldx, numerics *y, np.int64_t *ldy, int *blockSize, primme_svds_operator *mode, primme_svds_params *primme_svds, int *ierr):
    if blockSize[0] <= 0:
        ierr[0] = 0
        return
    ierr[0] = 1
    cdef object precond = primme_svds_params_get_object(primme_svds, 'preconditioner')
    m = primme_svds_params_get_int(primme_svds, "mLocal")
    n = primme_svds_params_get_int(primme_svds, "nLocal")
    global __user_function_exception
    try:
        if mode[0] == primme_svds_op_AtA:
            precond((<numerics[:ldx[0]:1,:blockSize[0]]>x)[:n,:], (<numerics[:ldy[0]:1,:blockSize[0]]>y)[:n,:], mode[0])
        elif mode[0] == primme_svds_op_AAt:
            precond((<numerics[:ldx[0]:1,:blockSize[0]]>x)[:m,:], (<numerics[:ldy[0]:1,:blockSize[0]]>y)[:m,:], mode[0])
        elif mode[0] == primme_svds_op_augmented:
            precond((<numerics[:ldx[0]:1,:blockSize[0]]>x)[:m+n,:], (<numerics[:ldy[0]:1,:blockSize[0]]>y)[:m+n,:], mode[0])
        else:
            return
        ierr[0] = 0
    except Exception as e:
        __user_function_exception = e

cdef void c_svds_monitor(numerics_real *basisSvals, int *basisSize, int *basisFlags, int *iblock, int *blockSize,
      numerics_real *basisNorms, int *numConverged, numerics_real *lockedSvals, int *numLocked, int *lockedFlags, numerics_real *lockedNorms,
      int *inner_its, numerics_real *LSRes, const char *msg, double *time, primme_event *event, int *stage, primme_svds_params *primme_svds, int *ierr):
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
        event[0] if event is not NULL else None,
        stage[0] if stage is not NULL else None)
    ierr[0] = 0

cdef void c_svds_convtest_numpy(double *sval, numerics *svecleft, numerics *svecright, double *resNorm, int *method, int *isconv, primme_svds_params *primme_svds, int *ierr):
    ierr[0] = 1
    cdef object convtest = primme_svds_params_get_object(primme_svds, 'convtest')
    if convtest is None: return
    global __user_function_exception
    try:
        m = primme_svds_params_get_int(primme_svds, "mLocal")
        n = primme_svds_params_get_int(primme_svds, "nLocal")
        isconv[0] = 1 if convtest(sval[0] if sval is not NULL else None,
            <numerics[:m]>svecleft if svecleft is not NULL else None,
            <numerics[:n]>svecright if svecright is not NULL else None,
            resNorm[0] if resNorm is not NULL else None) else 0
        ierr[0] = 0
    except Exception as e:
        __user_function_exception = e
 
IF BUILD_WITH_MAGMA: 
    cdef void c_svds_matvec_cupy(void *x, np.int64_t *ldx, void *y, np.int64_t *ldy, int *blockSize, int *transpose, primme_svds_params *primme_svds, int *ierr):
        if blockSize[0] <= 0:
            ierr[0] = 0
            return
        ierr[0] = 1
        global __user_function_exception
        cdef object A = primme_svds_params_get_object(primme_svds, 'matrix')
        try:
            m = primme_svds_params_get_int(primme_svds, "mLocal")
            n = primme_svds_params_get_int(primme_svds, "nLocal")
            x_shape = (n if transpose[0] == 0 else m, blockSize[0])
            y_shape = (m if transpose[0] == 0 else n, blockSize[0])
            A(A.__toObj(<size_t>x, x_shape, ldx[0]), A.__toObj(<size_t>y, y_shape, ldy[0]), transpose[0] == 0)
            ierr[0] = 0
        except Exception as e:
            __user_function_exception = e
    
    cdef void c_svds_precond_cupy(void *x, np.int64_t *ldx, void *y, np.int64_t *ldy, int *blockSize, primme_svds_operator *mode, primme_svds_params *primme_svds, int *ierr):
        if blockSize[0] <= 0:
            ierr[0] = 0
            return
        ierr[0] = 1
        global __user_function_exception
        cdef object precond = primme_svds_params_get_object(primme_svds, 'preconditioner')
        cdef object A = primme_svds_params_get_object(primme_svds, 'matrix')
        try:
            m = primme_svds_params_get_int(primme_svds, "mLocal")
            n = primme_svds_params_get_int(primme_svds, "nLocal")
            if mode[0] == primme_svds_op_AtA:
                shape = (n, blockSize[0])
            elif mode[0] == primme_svds_op_AAt:
                shape = (m, blockSize[0])
            elif mode[0] == primme_svds_op_augmented:
                shape = (m+n, blockSize[0])
            else:
                return
            precond(A.__toObj(<size_t>x, shape, ldx[0]), A.__toObj(<size_t>y, shape, ldy[0]), mode[0])
            ierr[0] = 0
        except Exception as e:
            __user_function_exception = e
    
    cdef void c_svds_convtest_cupy(double *sval, void *svecleft, void *svecright, double *resNorm, int *method, int *isconv, primme_svds_params *primme_svds, int *ierr):
        ierr[0] = 1
        cdef object convtest = primme_svds_params_get_object(primme_svds, 'convtest')
        cdef object A = primme_svds_params_get_object(primme_svds, 'matrix')
        if convtest is None: return
        global __user_function_exception
        try:
            m = primme_svds_params_get_int(primme_svds, "mLocal")
            n = primme_svds_params_get_int(primme_svds, "nLocal")
            isconv[0] = 1 if convtest(sval[0] if sval is not NULL else None,
                A.__toObj(<size_t>svecleft, (m,1), m) if svecleft is not NULL else None,
                A.__toObj(<size_t>svecright, (n,1), n) if svecright is not NULL else None,
                resNorm[0] if resNorm is not NULL else None) else 0
            ierr[0] = 0
        except Exception as e:
            __user_function_exception = e

def __get_svds_method(method):
    if method is None: return primme_svds_default
    cdef int method_int = -1;
    method0 = method
    if not method0.startswith('primme_svds_'):
        method0 = 'primme_svds' + method0;
    method0 = bytesp23(method0, 'ASCII')
    primme_svds_constant_info(<const char *>method0, &method_int)
    if method_int < 0:
        raise ValueError('Not valid "method": %s' % method)
    return method_int

def __to_svds_A_matvec(A, driver):
    m = A.shape[0]
    n = A.shape[1]
    dtype = A.dtype
    itemsize = dtype.itemsize
    if isinstance(A, (np.ndarray,scipy.sparse.spmatrix,scipy.sparse.linalg.LinearOperator)) or (__with_cupy and isinstance(A, (cupy.ndarray,cupyx.scipy.sparse.spmatrix))):
        def matvec(x, y, transpose): y[:,:] = A.dot(x) if not transpose else A.H.dot(x)
    elif __getnumargs(A) == 2:
        if driver is None: raise Exception('Please set `driver` when `A` is a function')
        def matvec(x, y, transpose): y[:,:] = A(x, transpose)
    elif __getnumargs(A) == 3:
        if driver is None: raise Exception('Please set `driver` when `A` is a function')
        def matvec(x, y, transpose): A(x, transpose, y)
    else:
        raise Exception('Unsupported value of `A` as matrix or function')

    if __with_cupy and (isinstance(A, (cupy.ndarray,cupyx.scipy.sparse.spmatrix)) or driver == 'cupy'):
        if driver == 'numpy':
            raise Exception('Incompatible value of `driver` with `A`')
        def toObj(ptr, shape, ld):
            return cupy.ndarray(shape, A.dtype, cupy.cuda.MemoryPointer(cupy.cuda.UnownedMemory(ptr,0,None),0), strides=(A.dtype.itemsize,A.dtype.itemsize*ld), order='F')
        matvec.__toObj = toObj

    return matvec

def svds(A, k=6, ncv=None, tol=0, which='LM', v0=None,
         maxiter=None, return_singular_vectors=True,
         precAHA=None, precAAH=None, precAug=None,
         u0=None, orthou0=None, orthov0=None, reltol=None, driver=None,
         return_stats=False, maxBlockSize=0,
         method=None, methodStage1=None, methodStage2=None,
         return_history=False, convtest=None, raise_for_unconverged=True, **kargs):
    """
    Compute k singular values and vectors of the matrix A.

    Parameters
    ----------
    A : matrix (numpy.ndarray, scipy.sparse.spmatrix, scipy.sparse.linalg.LinearOperator, cupy.ndarray,cupyx.scipy.sparse.spmatrix) or function with one argument
        Matrix to compute the SVD on, of shape (M, N).

        If it is a function, A(x) returns A*x and A.shape and
        A.dtype give the shape and the dtype of the matrix.

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
        the method can achieve.

        An error is raised if convtest is also provided.
    reltol : float, optional
        Relative tolerance for singular values. The default value is None.

        A triplet ``(u,sigma,v)`` is marked as converged when
        (||A*v - sigma*u||**2 + ||A.H*u - sigma*v||**2)**.5
        is less than "tol" * ||A|| if "tol" is given, and
        (||A*v - sigma*u||**2 + ||A.H*u - sigma*v||**2)**.5
        is less than "reltol" * sigma.

        An error is raised if convtest is also provided.
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
    precAHA : matrix (numpy.ndarray, scipy.sparse.spmatrix, scipy.sparse.linalg.LinearOperator, cupy.ndarray,cupyx.scipy.sparse.spmatrix) or function with one argument, optional
        Approximate inverse of (A.H*A - sigma**2*I). If provided and M>=N, it
        usually accelerates the convergence.
    precAAH : matrix (numpy.ndarray, scipy.sparse.spmatrix, scipy.sparse.linalg.LinearOperator, cupy.ndarray,cupyx.scipy.sparse.spmatrix) or function with one argument, optional
        Approximate inverse of (A*A.H - sigma**2*I). If provided and M<N, it
        usually accelerates the convergence.
    precAug : matrix (numpy.ndarray, scipy.sparse.spmatrix, scipy.sparse.linalg.LinearOperator, cupy.ndarray,cupyx.scipy.sparse.spmatrix) or function with one argument, optional
        Approximate inverse of ([zeros() A.H; zeros() A] - sigma*I).
    orthou0 : ndarray, optional
        Left orthogonal vector constrain.

        Seek singular triplets orthogonal to orthou0 and orthov0. The provided vectors
        *should* be orthonormal. If only orthou0 or orthov0 is provided, the other
        is computed. Useful to avoid converging to previously computed solutions.
    orthov0 : ndarray, optional
        Right orthogonal vector constrain. See orthou0.
    driver : string ['numpy' | 'cupy' ]
        Kind of object passed to A, precAHA, precAAH, precAug, and returned as svecs. If 'numpy' or 'cupy' is given, then
        the passed input and the expected returned object is from class `ndarray' in
        numpy or cupy respectively.
        
        If no value is given, and `A' is scipy.sparse.spmatrix or cupyx.scipy.sparse.spmatrix,
        then it is set to 'numpy' and 'cupy' respectively.
    maxBlockSize : int, optional
        Maximum number of vectors added at every iteration.
    convtest : callable
        User-defined function to mark an approximate singular triplet as converged.

        The function is called as convtest(sval, svecleft, svecright, resNorm)
        and returns True if the triplet with value `sval`, left vector `svecleft`,
        right vector `svecright`, and residual norm `resNorm` is considered converged.
    raise_for_unconverged: bool, optional
        If True, raise an exception when not all requested singular values
        converged. The default is True.

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
    >>> svecs_left, svals, svecs_right = primme.svds(A, 3, tol=1e-6, which='LM')
    >>> svals # the three largest singular values of A
    array([10., 9., 8.])

    >>> import primme, scipy.sparse, numpy as np
    >>> A = scipy.sparse.rand(10000, 100, random_state=10)
    >>> prec = scipy.sparse.spdiags(np.reciprocal(A.multiply(A).sum(axis=0)),
    ...           [0], 100, 100) # square diag. preconditioner
    >>> # the three smallest singular values of A, using preconditioning
    >>> svecs_left, svals, svecs_right = primme.svds(A, 3, which='SM', tol=1e-6, precAHA=prec)
    >>> ["%.5f" % x for x in svals.flat] # doctest: +SKIP
    ['4.57263', '4.78752', '4.82229']

    >>> # Giving the matvecs as functions
    >>> import primme, scipy.sparse, numpy as np
    >>> Bdiag = np.arange(0, 100).reshape((100,1))
    >>> Bdiagr = np.concatenate((np.arange(0, 100).reshape((100,1)).astype(np.float32), np.zeros((100,1), dtype=np.float32)), axis=None).reshape((200,1))
    >>> def Bmatmat(x):
    ...    if len(x.shape) == 1: x = x.reshape((100,1))
    ...    return np.vstack((Bdiag * x, np.zeros((100, x.shape[1]), dtype=np.float32)))
    ...
    >>> def Brmatmat(x):
    ...    if len(x.shape) == 1: x = x.reshape((200,1))
    ...    return (Bdiagr * x)[0:100,:]
    ...
    >>> B = scipy.sparse.linalg.LinearOperator((200,100), matvec=Bmatmat, matmat=Bmatmat, rmatvec=Brmatmat, dtype=np.float32)
    >>> svecs_left, svals, svecs_right = primme.svds(B, 5, which='LM', tol=1e-6)
    >>> svals # doctest: +SKIP
    array([99., 98., 97., 96., 95.])

    >>> import primme, scipy.sparse
    >>> import pycuda.sparse.packeted
    >>> A = scipy.sparse.rand(10000, 100, random_state=10)
    >>> Agpu = pycuda.sparse.packeted.PacketedSpMV(A, False, A.dtype)
    >>> evals, evecs = primme.eigsh(Agpu, 3, tol=1e-6, which='LA')
    >>> evals # the three largest eigenvalues of A
    """
    PP = PrimmeSvdsParams()
    cdef primme_svds_params *pp = PP.pp
 
    if driver is None:
        if isinstance(A, (np.ndarray,scipy.sparse.spmatrix,scipy.sparse.linalg.LinearOperator)):
            driver = 'numpy'
        elif __with_cupy and isinstance(A, (cupy.ndarray,cupyx.scipy.sparse.spmatrix)):
            driver = 'cupy'
        else:
            raise Exception('Please set `driver` when giving a function as `A`')

    try:
        Af = __to_svds_A_matvec(A, driver)
    except Exception as e:
        raise Exception('Invalid value for argument `A`)') from e
    __primme_svds_params_set(PP, "matrix", Af)

    cdef int m, n
    m, n = A.shape
    __primme_svds_params_set(PP, "m", m)
    __primme_svds_params_set(PP, "n", n)

    if k <= 0 or k > min(n, m):
        raise ValueError("k=%d must be between 1 and min(A.shape)=%d" % (k, min(n, m)))
    __primme_svds_params_set(PP, "numSvals", k)

    precAHAf = precAAHf = precAugf = None
    if precAHA is not None:
        if precAHA.shape[0] != precAHA.shape[1] or precAHA.shape[0] != n:
            raise ValueError('precAHA: expected square matrix with size %d' % n)
        precAHAf = __to_eigsh_matvec(precAHA, driver, 'precAHA')

    if precAAH is not None:
        if precAAH.shape[0] != precAAH.shape[1] or precAAH.shape[0] != m:
            raise ValueError('precAAH: expected square matrix with size %d' % m)
        precAAHf = __to_eigsh_matvec(precAAH, driver, 'precAAH')

    if precAug is not None:
        if precAug.shape[0] != precAug.shape[1] or precAug.shape[0] != m+n:
            raise ValueError('precAug: expected square matrix with size %d' % (m+n))
        precAugf = __to_eigsh_matvec(precAug, driver, 'precAug')

    def prevec(X, Y, mode):
        if mode == primme_svds_op_AtA and precAHAf is not None:
            precAHAf(X, Y)
        elif mode == primme_svds_op_AAt and precAAHf is not None:
            precAAHf(X, Y) 
        elif mode == primme_svds_op_augmented and precAugf is not None:
            precAugf(X, Y) 
        else:
            Y[:,:] = X[:,:]

    if precAHA is not None or precAAH is not None or precAug is not None:
        __primme_svds_params_set(PP, "precondition", 1)
        __primme_svds_params_set(PP, "preconditioner", prevec)
    else:
        __primme_svds_params_set(PP, "precondition", 0)

    hist = {"numMatvecs": [], "elapsedTime": [], "nconv": [],
            "sval": [], "resNorm": []}

    def mon(basisSvals, basisFlags, iblock, basisNorms, numConverged,
            lockedSvals, lockedFlags, lockedNorms, inner_its, LSRes,
            event, stage):
        if event == 0 and iblock and len(iblock)>0: # event iteration
            hist["numMatvecs"].append(__primme_svds_params_get(PP, 'stats_numMatvecs'))
            hist["elapsedTime"].append(__primme_svds_params_get(PP, 'stats_elapsedTime'))
            hist["nconv"].append(numConverged)
            hist["sval"].append(basisSvals[iblock[0]])
            hist["resNorm"].append(basisNorms[0])

    if return_history:
        __primme_svds_params_set(PP, 'monitor', mon)

    cdef double sigma_c
    if which == 'LM':
        __primme_svds_params_set(PP, "target", "primme_svds_largest")
    elif which == 'SM':
        __primme_svds_params_set(PP, "target", "primme_svds_smallest")
    else:
        try:
            sigma_c = float(which)
        except:
            raise ValueError("which must be either 'LM', 'SM' or a number.")
        __primme_svds_params_set(PP, "numTargetShifts", 1)
        primme_svds_params_set_doubles(pp, "targetShifts", &sigma_c)
        __primme_svds_params_set(PP, "target", "primme_svds_closest_abs")

    if tol:
        if convtest is not None:
            raise ValueError('Giving tol and convtest is not allowed.')
        __primme_svds_params_set(PP, "eps", tol)

    if ncv:
        __primme_svds_params_set(PP, "maxBasisSize", ncv)

    if maxiter:
        # NOTE: every eigensolver iteration spend two matvecs*blockSize
        __primme_svds_params_set(PP, "maxMatvecs", maxiter*(maxBlockSize if maxBlockSize else 1)/2)

    if reltol is not None:
        if convtest is not None:
            raise ValueError('Giving reltol and convtest is not allowed.')
        def check_convergence(sval, lsvec, rsvec, resNorm):
            if resNorm is None: return True
            if tol is not None:
                aNorm = __primme_svds_params_get(PP, "aNorm")
                if resNorm > aNorm*tol: return False
            if sval is not None and resNorm > sval*reltol:
                return False
            return True

        convtest = check_convergence


    if maxBlockSize:
        __primme_svds_params_set(PP, "maxBlockSize", maxBlockSize)

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
            v, _ = np.linalg.qr(A.rmatmult(u))

        if v is not None and u is None:
            u, _ = np.linalg.qr(A.matmult(v))

        return u, v

    orthou0, orthov0 = check_pair(orthou0, orthov0, "orthov0 or orthou0")

    cdef int numOrthoConst = 0
    if orthou0 is not None:
        numOrthoConst = min(orthou0.shape[1], min(m,n))
        __primme_svds_params_set(PP, "numOrthoConst", numOrthoConst)

    if convtest is not None:
        __primme_svds_params_set(PP, "convtest", convtest)

    # Set other parameters
    for dk, dv in kargs.items():
        __primme_svds_params_set(PP, dk, dv)

    if A.dtype.kind in frozenset(["b", "i", "u"]) or A.dtype.type is np.double:
        dtype = np.dtype("d")
    else:
        dtype = A.dtype


    cdef void *svecs_p = <void*>0
    if driver == 'numpy':
        if dtype.type is np.complex64:
            primme_svds_params_set_pointer(pp, "matrixMatvec", <void*>c_svds_matvec_numpy[np.complex64_t])
            primme_svds_params_set_pointer(pp, "applyPreconditioner", <void*>c_svds_precond_numpy[np.complex64_t])
            if convtest:
                primme_svds_params_set_pointer(pp, "convTestFun", <void*>c_svds_convtest_numpy[np.complex64_t])
        elif dtype.type is np.float32:
            primme_svds_params_set_pointer(pp, "matrixMatvec", <void*>c_svds_matvec_numpy[float])
            primme_svds_params_set_pointer(pp, "applyPreconditioner", <void*>c_svds_precond_numpy[float])
            if convtest:
                primme_svds_params_set_pointer(pp, "convTestFun", <void*>c_svds_convtest_numpy[float])
        elif dtype.type is np.float64:
            primme_svds_params_set_pointer(pp, "matrixMatvec", <void*>c_svds_matvec_numpy[double])
            primme_svds_params_set_pointer(pp, "applyPreconditioner", <void*>c_svds_precond_numpy[double])
            if convtest:
                primme_svds_params_set_pointer(pp, "convTestFun", <void*>c_svds_convtest_numpy[double])
        else:
            primme_svds_params_set_pointer(pp, "matrixMatvec", <void*>c_svds_matvec_numpy[np.complex128_t])
            primme_svds_params_set_pointer(pp, "applyPreconditioner", <void*>c_svds_precond_numpy[np.complex128_t])
            if convtest:
                primme_svds_params_set_pointer(pp, "convTestFun", <void*>c_svds_convtest_numpy[np.complex128_t])
    else:
        IF BUILD_WITH_MAGMA: 
            primme_svds_params_set_pointer(pp, "matrixMatvec", <void*>c_svds_matvec_cupy)
            primme_svds_params_set_pointer(pp, "applyPreconditioner", <void*>c_svds_precond_cupy)
            if convtest:
                primme_svds_params_set_pointer(pp, "convTestFun", <void*>c_svds_convtest_cupy)
        ELSE:
            raise Exception('primme is not compiled with cupy support')

    if return_history:
        primme_svds_params_set_pointer(pp, "monitorFun", <void*>c_svds_monitor[float] if dtype.type is np.float32 or dtype.type is np.complex64 else <void*>c_svds_monitor[double])

    cdef double[::1] svals_d, norms_d
    cdef float[::1] svals_s, norms_s
    cdef float[::1] svecs_s
    cdef double[::1] svecs_d
    cdef np.complex64_t[::1] svecs_c
    cdef np.complex128_t[::1] svecs_z

    rtype = __get_real_dtype(dtype);
    svals = np.zeros(k, rtype)
    norms = np.zeros(k, rtype)
    if rtype.type is np.float64:
        svals_d, norms_d = svals, norms
    else:
        svals_s, norms_s = svals, norms
 
    if driver == 'numpy':
        svecs = np.empty(((m+n)*(numOrthoConst+k),), dtype)
        if dtype.type is np.float64:
            svecs_d = svecs
            svecs_p = <void*>&svecs_d[0]
        elif dtype.type is np.float32:
            svecs_s = svecs
            svecs_p = <void*>&svecs_s[0]
        elif dtype.type is np.complex64:
            svecs_c = svecs
            svecs_p = <void*>&svecs_c[0]
        elif dtype.type is np.complex128:
            svecs_z = svecs
            svecs_p = <void*>&svecs_z[0]
    elif driver == 'cupy':
        svecs = cupy.zeros(((m+n)*(numOrthoConst+k),), dtype, order='F')
        svecs_p = <void*>(<size_t>svecs.data.ptr)

    u0, v0 = check_pair(u0, v0, "v0 or u0")
    
    # Set method
    cdef int method_int = -1, methodStage1_int = -1, methodStage2_int = -1;
    if method is not None or methodStage1 is not None or methodStage2 is not None:
        method_int = __get_svds_method(method)
        methodStage1_int = __get_method(methodStage1)
        methodStage2_int = __get_method(methodStage2)
        primme_svds_set_method(<primme_svds_preset_method>method_int, <primme_preset_method>methodStage1_int, <primme_preset_method>methodStage2_int, pp)
        # Set other parameters (again)
        for dk, dv in kargs.items():
            try:
                __primme_svds_params_set(PP, dk, dv)
            except Exception as e:
                raise ValueError("Invalid value in field '%s': %s\n%s" % (dk, dv, e))


    cdef int initSize = 0
    if v0 is not None:
        initSize = min(v0.shape[1], k)
        __primme_svds_params_set(PP, "initSize", initSize)
        svecs[m*numOrthoConst:m*(numOrthoConst+initSize)].reshape((m,initSize), order='F')[:,:] = u0[:,:initSize]
        svecs[m*(numOrthoConst+initSize)+n*numOrthoConst:(m+n)*(numOrthoConst+initSize)].reshape((n,initSize), order='F')[:,:] = v0[:,:initSize]

    if orthou0 is not None:
        svecs[0:m*numOrthoConst].reshape((m,numOrthoConst), order='F')[:,:] = orthou0[:,0:numOrthoConst]
        svecs[m*(numOrthoConst+initSize):m*(numOrthoConst+initSize)+n*numOrthoConst].reshape((n,numOrthoConst), order='F')[:,:] = orthov0[:,:initSize]

    if __primme_svds_params_get(PP, "printLevel") >= 5:
        primme_svds_display_params(pp[0]);

    global __user_function_exception
    __user_function_exception = None
    err = 0 # avoid warning
    if driver == 'numpy':
        if dtype.type is np.complex64:
            err = cprimme_svds(&svals_s[0], <PRIMME_COMPLEX_FLOAT*>svecs_p, &norms_s[0], pp)
        elif dtype.type is np.float32:
            err = sprimme_svds(&svals_s[0], <float*>svecs_p, &norms_s[0], pp)
        elif dtype.type is np.float64:
            err = dprimme_svds(&svals_d[0], <double*>svecs_p, &norms_d[0], pp)
        else:
            err = zprimme_svds(&svals_d[0], <PRIMME_COMPLEX_DOUBLE*>svecs_p, &norms_d[0], pp)
    elif driver == 'cupy':
        if dtype.type is np.complex64:
            err = magma_cprimme_svds(&svals_s[0], <PRIMME_COMPLEX_FLOAT*>svecs_p, &norms_s[0], pp)
        elif dtype.type is np.float32:
            err = magma_sprimme_svds(&svals_s[0], <float*>svecs_p, &norms_s[0], pp)
        elif dtype.type is np.float64:
            err = magma_dprimme_svds(&svals_d[0], <double*>svecs_p, &norms_d[0], pp)
        else:
            err = magma_zprimme_svds(&svals_d[0], <PRIMME_COMPLEX_DOUBLE*>svecs_p, &norms_d[0], pp)

    if err != 0:
        if __user_function_exception is not None:
            raise PrimmeSvdsError(err) from __user_function_exception
        elif err == -3 and raise_for_unconverged:
            raise PrimmeSvdsError(err)

    if return_stats:
        stats = dict((f, __primme_svds_params_get(PP, 'stats_' + f)) for f in [
            "numOuterIterations", "numRestarts", "numMatvecs",
            "numPreconds", "elapsedTime"])
        stats["rnorms"] = norms
        if return_history:
            stats["hist"] = hist
 
    initSize = __primme_svds_params_get(PP, "initSize")
    svals = svals[0:initSize]
    if not return_singular_vectors:
        return svals if not return_stats else (svals, stats)

    numOrthoConst = __primme_svds_params_get(PP, "numOrthoConst")
    norms = norms[0:initSize]

    # Make copies and transpose conjugate svecsr
    svecsl = svecs[m*numOrthoConst:m*(numOrthoConst+initSize)].reshape((m,initSize), order='F').T.copy().T
    svecsr = svecs[m*(numOrthoConst+initSize)+n*numOrthoConst:(m+n)*(numOrthoConst+initSize)].reshape((n,initSize), order='F').T.conj().copy()

    if not return_stats:
        return svecsl, svals, svecsr
    else:
        return svecsl, svals, svecsr, stats

#
# PRIMME Errors
#


_PRIMMEErrors = {
0: "success",
-1: "unexpected internal error; please consider to set 'printLevel' to a value larger than 0 to see the call stack and to report these errors because they may be bugs",
-2: "memory allocation failure",
-3: "maximum iterations or matvecs reached",
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
-30: "'evals' is NULL",
-31: "'evecs' is NULL",
-32: "'resNorms' is NULL",
-33: "'locking' == 0 and 'minRestartSize' < 'numEvals'",
-34: "'ldevecs' is less than 'nLocal'",
-35: "'ldOPs' is non-zero and less than 'nLocal'",
-38: "'locking' == 0 and 'target' is 'primme_closest_leq' or 'primme_closet_geq'",
-40: "some LAPACK function performing a factorization returned an error code; set 'printLevel' > 0 to see the error code and the call stack",
-41: "error happened at the matvec or applying the preconditioner",
-42: "the matrix provided in 'lock' is not full rank",
-43: "parallel failure",
-44: "unavailable functionality; PRIMME was not compiled with support for the requesting precision or for GPUs"
}

_PRIMMESvdsErrors = {
0   : "success",
-1  : "unexpected internal error; please consider to set 'printLevel' to a value larger than 0 to see the call stack and to report these errors because they may be bugs",
-2  : "memory allocation failure",
-3  : "maximum iterations or matvecs reached",
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
-40 : "some LAPACK function performing a factorization returned an error code; set 'printLevel' > 0 to see the error code and the call stack",
-41 : "error happened at the matvec or applying the preconditioner",
-42 : "the matrix provided in 'lock' is not full rank",
-43 : "parallel failure",
-44 : "unavailable functionality; PRIMME was not compiled with support for the requesting precision or for GPUs"
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



