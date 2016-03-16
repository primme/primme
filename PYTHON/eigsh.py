import Primme
import numpy as np
from scipy.sparse import spdiags

def eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None, ncv=None, maxiter=None, tol=0, return_eigenvectors=True, Minv=None, OPinv=None, mode='normal'):

    class PP(Primme.primme_params_w):
        def __init__(self):
            Primme.primme_params_w.__init__(self)
        def matvec(self):
            self.setY(A*self.getX())

    pp = PP()
    pp.n = A.shape[0]
    pp.numEvals = k

    if sigma != None:
        pp.numTargetShifts = 1
        pp.targetShifts = sigma

    if which == 'LA':
        pp.target = Primme.primme_largest
    elif which == 'SA':
        pp.target = Primme.primme_smallest
    else:
        print('not supported') # THROW ERROR HERE
        return

    if v0 != None:
        pp.initSize = len(v0)

    if ncv != None:
        pp.maxBasisSize = ncv

    if maxiter != None:
        pp.maxMatvecs = maxiter

    pp.set_method(Primme.DYNAMIC)

    evals = np.zeros(pp.numEvals)
    norms = np.zeros(pp.numEvals)

    if A.dtype is np.dtype(np.complex64):
        evecs = np.zeros((pp.n, pp.numEvals), complex)
        Primme.zprimme(evals, evecs, norms, pp)
    elif A.dtype is np.dtype('d'):
        evecs = np.zeros((pp.n, pp.numEvals))
        print(Primme.dprimme(evals, evecs, norms, pp))

    return { 'w': evals, 'v': evecs }


if __name__ == '__main__':
    a = np.ones(10)
    A  = spdiags(np.array([a*(-1.), a*2., a*(-1.)]), np.array([-1, 0, 1]), 10, 10)
    print(eigsh(A, which='LA'))

