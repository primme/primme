#!/usr/bin/env python
"""
Test Functions for primme eigsh and svds.
Based on test_lobpcg.py and test_arpack.py in scipy,
see https://github.com/scipy/scipy
"""

import warnings
import numpy as np
from numpy.testing import run_module_suite, assert_allclose
from scipy import ones, r_, diag
from scipy.sparse.linalg import aslinearoperator
from scipy.sparse import csr_matrix
import math
import primme
from primme import eigsh, svds
from compare import stats as st
from builtins import str

#
# Collection of problems
#

def ElasticRod(n, dtype=np.dtype("d")):
    # Fixed-free elastic rod
    L = 1.0
    le = L/n
    rho = 7.85e3
    S = 1.e-4
    E = 2.1e11
    mass = rho*S*le/6.
    k = E*S/le
    c = math.sqrt(k)
    A = c*(diag(r_[2.*ones(n-1, dtype=dtype),1])-diag(ones(n-1, dtype=dtype),1)-diag(ones(n-1, dtype=dtype),-1))
    B = mass/c*(diag(r_[4.*ones(n-1, dtype=dtype),2])+diag(ones(n-1, dtype=dtype),1)+diag(ones(n-1, dtype=dtype),-1))
    return A, B

def MikotaPair(n, dtype=np.dtype("d")):
    # Mikota pair acts as a nice test since the eigenvalues
    # are the squares of the integers n, n=1,2,...
    x = np.arange(1,n+1, dtype=dtype)
    B = diag(1./x)
    y = np.arange(n-1,0,-1, dtype=dtype)
    z = np.arange(2*n-1,0,-2, dtype=dtype)
    A = diag(z)-diag(y,-1)-diag(y,1)
    return A, B

def diagonal(n, dtype=np.dtype("d")):
    # Diagonal problem
    vals = np.arange(1, n+1, dtype=dtype)
    A = diag(vals)
    return A, None

def Lauchli_like(m, n=None, c=1e-5, dtype=float):
   k = min(m-1, n)
   diagonal = diag(np.array(np.linspace(1.0, 5*c, k), dtype=dtype))
   d = np.c_[diagonal, np.zeros((k, n-k), dtype=dtype)]
   return np.r_[ones((1,n), dtype), d, np.zeros((m-k-1,n), dtype=dtype)]

#
# Tools
#

def toStandardProblem(p):
   """
   Return a symmetric/Hermitian problem if the original was generalized.
   """

   A, B = p
   if B is None:
      return A
   L = np.linalg.cholesky(B)
   # A = L\(A/L')
   A = np.linalg.solve(L, np.linalg.solve(L, A.T.conj()).T.conj())
   return A

def jacobi_prec(A, sigma):
   """
   Return the Jacobi preconditioner.
   """

   if sigma is None: sigma = 0
   with warnings.catch_warnings():
      warnings.simplefilter("error")
      while True:
         try:
            return np.diag(np.reciprocal(np.diag(A)) - sigma)
         except RuntimeWarning:
            sigma = sigma*1.001 if sigma != 0.0 else 0.1

def sqr_diagonal_prec(A, sigma):
   """
   Return the square diagonal preconditioner.
   """

   if sigma is None: sigma = 0
   with warnings.catch_warnings():
      warnings.simplefilter("error")
      while True:
         try:
            return {'precAHA': np.diag(np.reciprocal(np.diag(A.T.conj().dot(A)) - sigma)),
                    'precAAH': np.diag(np.reciprocal(np.diag(A.dot(A.T.conj())) - sigma))}
         except RuntimeWarning:
            sigma = sigma*1.001 if sigma != 0.0 else 0.1

#
# Tests
#
         
def select_pairs_eigsh(k, sigma, which, evals):
   """
   Return the k values that eigsh should return for that sigma and which.
   """

   if sigma == None:
      sigma = 0
   if which == 'LM':
      f = lambda x : -abs(x - sigma)
   elif which == 'SM':
      f = lambda x : abs(x - sigma)
   elif which == 'LA':
      f = lambda x : -x
   elif which == 'SA':
      f = lambda x : x
   else:
      raise ValueError("which='%s' not supported" % which)

   n = max(evals.shape)
   return np.array(sorted(evals, key=f)[0:k])

def to_primme_datatype(dtype):
   if isinstance(dtype, (str,bytes)):
      return dtype
   if dtype == np.float16:
      return 'primme_op_half'
   elif dtype == np.float32 or dtype == np.complex64:
      return 'primme_op_float'
   else:
      return 'primme_op_double'

def dtype_to_str(precision, complexity):
    if complexity == np.float64:
        return str(precision)
    if precision == np.float16:
        return "<class 'numpy.complex32'>"
    if precision == np.float32:
        return "<class 'numpy.complex64'>"
    if precision == np.float64:
        return "<class 'numpy.complex128'>"
    
def eigsh_check(eigsh_solver, A, B, normInvB, k, M, which, sigma, tol,
                exact_evals, dtype, case_desc, add_stats=True):
   """
   Test eigsh
   """

   try:
      evals, evecs, stats = eigsh_solver(A, k, B, sigma, which, tol=tol, OPinv=M,
            maxMatvecs=70000, return_stats=True, internalPrecision=to_primme_datatype(dtype))
   except Exception as e:
      raise Exception("Ups! Case %s\n%s" % (case_desc, e))
   sol_evals = select_pairs_eigsh(k, sigma, which, exact_evals)

   # Check eigenvalues are close enough to the exact ones
   ANorm = np.amax(np.fabs(exact_evals))
   assert_allclose(evals, sol_evals, atol=ANorm*tol*normInvB, rtol=1, err_msg=case_desc)

   # Check the residual norm associated to the returned pairs
   if B is None: 
      R = A.dot(evecs) - evecs.dot(np.diag(evals))
   else:
      R = A.dot(evecs) - B.dot(evecs.dot(np.diag(evals)))
   Rnorms = np.linalg.norm(R, axis=0)
   eps = np.finfo(dtype).eps
   assert_allclose(Rnorms, np.zeros(k), atol=ANorm*(tol+eps)*(k**.5), rtol=1, err_msg=case_desc + (" |A|=%s tol=%f" % (ANorm, tol)))

   # Add stats
   if add_stats:
      st.add("eigsh: " + case_desc, ("eigsh",), mv=stats['numMatvecs'], time=stats['elapsedTime'])

def test_primme_eigsh():
   """
   Test cases for primme.eighs for standard problems.
   """

   for n in (2, 3, 5, 10, 100):
      for gen in (ElasticRod, MikotaPair, diagonal):
         evals = np.linalg.eigvalsh(toStandardProblem(gen(n)))
         for complexity in (np.float64, np.complex128):
            A = toStandardProblem(gen(n, dtype=complexity))
            sigma0 = evals[0]*.51 + evals[-1]*.49
            nA = np.max(np.fabs(evals))
            for precision in (np.float16, np.float32, np.float64):
               tol = np.finfo(precision).eps**.5 * 0.1
               for which, sigma in [(w, None) for w in ('LM', 'SM', 'LA', 'SA')] + [('SM', sigma0)] :
                  if gen.__name__ != "ElasticRod" and sigma is not None:
                     precs = (None, jacobi_prec(A, sigma))
                  else:
                     precs = (None,)
                  if precision is np.float16 and which != 'LM':
                     continue
                  if precision is np.float16 and gen.__name__ == "ElasticRod":
                     A0, evals0 = A/nA, evals/nA
                  else:
                     A0, evals0 = A, evals
                  for prec in precs:
                     for k in (1, 2, 3, 5, 10, 70):
                        if k > n: continue
                        case_desc = ("A=%s(%d, %s), k=%d, M=%s, which=%s, sigma=%s" %
                              (gen.__name__, n, dtype_to_str(precision, complexity), k, prec is not None, which, sigma))
                        yield (eigsh_check, eigsh, A0, None, 1, k, prec, which, sigma, tol, evals0, precision, case_desc)

def test_primme_eigsh_gen():
   """
   Test cases for primme.eighs for generalized problems.
   """

   for n in (2, 3, 5, 10, 100):
      for gen in (ElasticRod, MikotaPair):
         A, B = gen(n)
         evals = np.linalg.eigvalsh(toStandardProblem((A, B)))
         normInvB = 1./min(np.linalg.eigvalsh(B))
         nA = np.max(np.fabs(evals))
         sigma0 = evals[0]*.51 + evals[-1]*.49
         for complexity in (np.float64, np.complex128):
            A, B = gen(n, dtype=complexity)
            stdP = toStandardProblem((A,B))
            for precision in (np.float16, np.float32, np.float64):
               tol = np.finfo(precision).eps**.5 * 0.1
               for which, sigma in [(w, None) for w in ('LM', 'SM', 'LA', 'SA')] + [('SM', sigma0)] :
                  if gen.__name__ != "ElasticRod" and sigma is not None:
                     precs = (None, jacobi_prec(stdP, sigma))
                  else:
                     precs = (None,)
                  if precision is np.float16 and which != 'LM':
                     continue
                  if precision is np.float16 and gen.__name__ == "ElasticRod":
                     A0, evals0 = A/nA, evals/nA
                  else:
                     A0, evals0 = A, evals
                  for prec in precs:
                     for k in (1, 2, 3, 5, 10, 50):
                        if k > n: continue
                        case_desc = ("A,B=%s(%d, %s), k=%d, M=%s, which=%s, sigma=%s" %
                              (gen.__name__, n, dtype_to_str(precision, complexity), k, prec is not None, which, sigma))
                        yield (eigsh_check, eigsh, A0, B, normInvB, k, prec, which, sigma, tol, evals0, precision, case_desc)

def test_primme_eigsh_matrix_types():
   """
   Test cases for primme.eighs with csr and LinearOperator matrix types.
   """
   n = 10
   for dtype in (np.float64, np.complex64):
      A = toStandardProblem(MikotaPair(n, dtype=dtype))
      evals, evecs = np.linalg.eigh(A)
      sigma0 = evals[0]*.51 + evals[-1]*.49
      for op in ((lambda x : x), csr_matrix, aslinearoperator): 
         which, sigma = 'SM', sigma0
         prec = jacobi_prec(A, sigma)
         k = 5
         M = op(prec) if prec is not None else None
         case_desc = ("A=%s(%d, %s), k=%d, M=%s, which=%s, sigma=%s" %
                      (MikotaPair.__name__, n, dtype, k, prec is None, which, sigma))
         yield (eigsh_check, eigsh, op(A), None, 1, k, M, which, sigma, 1e-6, evals, dtype, case_desc, False)


def select_pairs_svds(k, which, svals):
   """
   Return the k values that svds should return for that which.
   """

   if which == 'LM':
      f = lambda x : -x
   elif which == 'SM':
      f = lambda x : x
   else:
      sigma = float(which)
      f = lambda x : abs(x - sigma)

   n = max(svals.shape)
   return np.array(sorted(svals, key=f)[0:k])
 
def svds_check(svds_solver, A, k, M, which, tol, exact_svals, dtype, case_desc, add_stats=True):
   """
   Test svds
   """

   try:
      svl, sva, svr, stats = svds_solver(A, k, None, which=which, tol=tol,
            maxMatvecs=30000, return_stats=True, internalPrecision=to_primme_datatype(dtype), **M)
   except Exception as e:
      raise Exception("Ups! Case %s\n%s" % (case_desc, e))
   sol_svals = select_pairs_svds(k, which, exact_svals)
   svr = svr.T.conj()

   # Check singular values are close enough to the exact ones
   ANorm = np.amax(exact_svals)
   assert_allclose(sorted(sva), sorted(sol_svals), atol=ANorm*tol, rtol=1, err_msg=case_desc)

   # Check the residual norm associated to the returned pairs
   R = A.dot(svr) - svl.dot(np.diag(sva))
   Rnorms = np.linalg.norm(R, axis=0)
   eps = np.finfo(dtype).eps
   assert_allclose(Rnorms, np.zeros(k), atol=ANorm*(tol+eps)*(k**.5), rtol=1, err_msg=case_desc)

   # Add stats
   if add_stats:
      st.add("svds: " + case_desc, ("svds",), mv=stats['numMatvecs'], time=stats['elapsedTime'])

def test_primme_svds():
   """
   Generate all test cases for primme.svds.
   """

   for n in (2, 3, 5, 10, 50, 100):
      for precision in (np.float16, np.float32, np.float64):
         c = np.finfo(precision).eps**.333
         for gen_name, gen in (("MikotaPair", (lambda n, d: toStandardProblem(MikotaPair(n, dtype=d)))),
                            ("Lauchli_like_vert", (lambda n, d: Lauchli_like(n*2, n, c, dtype=d))),
                            ("Lauchli_like_hori", (lambda n, d: Lauchli_like(n, n*2, c, dtype=d)))):
            sva = np.linalg.svd(gen(n, np.float64), full_matrices=False, compute_uv=False)
            sigma0 = sva[0]*.51 + sva[-1]*.49
            for complexity in (np.float64, np.complex128):
               A = gen(n, complexity)
               tol = np.finfo(precision).eps**.5 * 0.1
               for which, sigma in [('LM', 0), ('SM', 0), (sigma0, sigma0)]:
                  for prec in (({},) if which == 'LM' else ({}, sqr_diagonal_prec(A, sigma))):
                     # If the condition number is too large, the first stage may end
                     # with approximations of larger values than the actual smallest
                     if (gen_name == "MikotaPair" and n > 50
                              and (precision is np.float32 or precision is np.complex64)
                              and which != 'LM'):
                        continue
                     if precision is np.float16 and which != 'LM':
                        continue
                     if (n >= 50 and
                         precision is np.float16 and which == 'LM'):
                        continue
                     for k in (1, 2, 3, 5, 10, 15):
                        if k > n: continue
                        case_desc = ("A=%s(%d, %s), k=%d, M=%s, which=%s, tol=%g" %
                              (gen_name, n, dtype_to_str(precision, complexity), k, bool(prec), which, tol))
                        yield (svds_check, svds, A, k, prec, which, tol, sva, precision, case_desc)

def test_primme_svds_matrix_types():
   """
   Generate all test cases for primme.svds with csr and LinearOperator matrix types.
   """

   n = 10
   for dtype in (np.float64, np.complex64):
      A = Lauchli_like(n*2, n, dtype=dtype)
      svl, sva, svr = np.linalg.svd(A, full_matrices=False)
      sigma0 = sva[0]*.51 + sva[-1]*.49
      for op in ((lambda x : x), csr_matrix, aslinearoperator): 
         which, sigma = 'SM', 0
         prec = sqr_diagonal_prec(A, sigma)
         k = 2
         case_desc = ("A=%s(%d, %s), k=%d, M=%s, which=%s" %
                      ("Lauchli_like_vert", n, dtype, k, bool(prec), which))
         yield (svds_check, svds, op(A), k, prec, which, 1e-5, sva, dtype, case_desc, False)

def test_examples_from_doc():
   import doctest
   doctest.testmod(primme, raise_on_error=True, optionflags=doctest.NORMALIZE_WHITESPACE)

def test_return_stats():
    A, _ = diagonal(100)
    evals, evecs, stats = primme.eigsh(A, 3, tol=1e-6, which='LA',
            return_stats=True, return_history=True)
    assert(stats["hist"]["numMatvecs"])

    svecs_left, svals, svecs_right, stats = primme.svds(A, 3, tol=1e-6,
            which='SM', return_stats=True, return_history=True)
    assert(stats["hist"]["numMatvecs"])


if __name__ == "__main__":
    run_module_suite()
    st.dump('tests.json')
