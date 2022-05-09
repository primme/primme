#!/usr/bin/env python
"""
Test Functions for primme eigsh and svds.
Based on test_lobpcg.py and test_arpack.py in scipy,
see https://github.com/scipy/scipy
"""

import warnings
import numpy as np
from numpy.testing import run_module_suite, assert_allclose, assert_
#from nose.tools import nottest
from scipy import ones, r_, diag
from scipy.sparse.linalg import aslinearoperator
from scipy.sparse import csr_matrix
import math
import primme
from primme import eigsh, svds
from compare import stats as st
from builtins import str

TEST_THIS_CASE = None
#TEST_THIS_CASE = "A=ElasticRod(100, <class 'numpy.float16'>), k=2, M=False, which=SA, sigma=None, tol=0.00390625, bs=3, method=DEFAULT_MIN_MATVECS, with_gpu=False"
#TEST_THIS_CASE = "A=ElasticRod(100, <class 'numpy.complex128'>), k=70, M=False, which=SM, sigma=1572736246213.813, tol=0.00371088, bs=2, method=DEFAULT_MIN_MATVECS, with_gpu=False"

#
# GPU stuff
#

try:
    import cupy
    test_gpu = True
    raise Exception("caca")
except Exception:
    test_gpu = False
    print("Not testing GPU interface")

def togpu(A, dtype=None):
   if A is None: return A
   return cupy.asarray(A, dtype)

def tocpu(A):
   if isinstance(A, np.ndarray): return A
   return A.get()

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
   sorted_evals = sorted(evals, key=f)
   sel_evals = np.array(sorted_evals[0:k])
   if k < n:
      gap = f(sorted_evals[k]) - f(sorted_evals[k-1])
   else:
      gap = 0
   return sel_evals, gap

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

def get_tol_eigsh(k, sigma, which, evals, precision):
    """Return a tolerance considering the spectrum"""

    # Get the smallest eigenvalue (in absolute value) that is going to converge
    # and the smallest distance between the converged values and the rest of
    # the spectrum, gap
    sel_evals, gap = select_pairs_eigsh(k, sigma, which, evals)
    smallest_eval = np.fabs(sel_evals).min()

    # A good tolerance is smaller than the smallest eigenvalue and the gap
    nA = np.fabs(evals).max()
    tol = np.min([smallest_eval*.1, gap])/nA

    # Limit the tolerance to be larger than the used precision's machine epsilon 
    return np.max([tol, np.finfo(precision).eps**.8])
    
def eigsh_check(eigsh_solver, A, B, normInvB, k, M, which, sigma, tol, bs, method,
                exact_evals, dtype, case_desc, with_gpu=False, add_stats=True):
   """
   Test eigsh
   """

   if TEST_THIS_CASE and TEST_THIS_CASE != case_desc:
      return

   extra = {'printLevel': 5 } if TEST_THIS_CASE == case_desc else {}

   import sys
   sys.stdout.write("Case %s\n" % case_desc)
   sys.stdout.flush()

   f = togpu if with_gpu else lambda x:x 
   try:
      evals, evecs, stats = eigsh_solver(f(A), k, f(B), sigma, which, tol=tol, reltol=.5,
            OPinv=f(M), maxBlockSize=bs, method=method, maxMatvecs=70000 if B is None else 700000,
            return_stats=True, internalPrecision=to_primme_datatype(dtype),
            driver="numpy" if not with_gpu else "cupy", **extra)
   except Exception as e:
      raise Exception("Ups! Case %s\n%s" % (case_desc, e))
   evecs = tocpu(evecs)
   sol_evals, _ = select_pairs_eigsh(k, sigma, which, exact_evals)

   # Check eigenvalues are close enough to the exact ones
   ANorm = np.amax(np.fabs(exact_evals))
   assert_allclose(np.sort(evals), np.sort(sol_evals), atol=ANorm*tol*normInvB, rtol=1, err_msg=case_desc)

   # Check the residual norm associated to the returned pairs
   if B is None: 
      R = A.dot(evecs) - evecs.dot(np.diag(evals))
   else:
      R = A.dot(evecs) - B.dot(evecs.dot(np.diag(evals)))
   Rnorms = np.linalg.norm(R, axis=0)
   eps = np.finfo(dtype).eps
   for i in range(k):
      if Rnorms[i] >= ANorm*(tol+eps)*(k**.5):
         if B is None:
            res_defl_i = np.linalg.norm(R[:,i] - evecs.dot(evecs.T.conj().dot(R[:,i])), axis=0)
         else:
            res_defl_i = np.linalg.norm(R[:,i] - B.dot(evecs.dot(evecs.T.conj().dot(R[:,i]))), axis=0)
         assert_(res_defl_i < ANorm*(tol+eps),
            "The pair %d failed the error checking for case_desc %s:\n|r|=%g  |(I-XX')r|=%g that should be <= %g"
               % (i, case_desc, Rnorms[i], res_defl_i, ANorm*tol))
         
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
            cA = nA/np.min(np.fabs(evals)[np.fabs(evals) > 0])
            for with_gpu in ((False, True) if test_gpu else (False,)):
               for precision in (np.float16, np.float32, np.float64):
                  for which, sigma in [(w, None) for w in ('LM', 'SM', 'LA', 'SA')] + [('SM', sigma0)] :
                     # CUBLAS does not support complex half with gpus
                     if with_gpu and complexity is np.complex128 and precision is np.float16:
                        continue
                     # Normalize the problem in half precision to avoid overflows
                     if precision is np.float16:
                        A0, evals0, nA0 = A/nA, evals/nA, 1
                        sigma = sigma/nA if sigma is not None else None
                     else:
                        A0, evals0, nA0 = A, evals, nA
                     # Jacobi preconditioner is terrible for ElasticRod
                     if gen.__name__ != "ElasticRod" and sigma is not None:
                        precs = (None, jacobi_prec(A, sigma))
                     else:
                        precs = (None,)
                     for method in ("DEFAULT_MIN_MATVECS", "DEFAULT_MIN_TIME"):
                        for prec in precs:
                           for k in (1, 2, 3, 5, 10, 70):
                              if k > n: continue
                              tol = get_tol_eigsh(k, sigma, which, evals0, precision)
                              for bs in (1, 2, 3, 5):
                                 if bs*3 > n: continue
                                 case_desc = ("A=%s(%d, %s), k=%d, M=%s, which=%s, sigma=%s, tol=%g, bs=%d, method=%s, with_gpu=%s" %
                                       (gen.__name__, n, dtype_to_str(precision, complexity), k, prec is not None, which, sigma, tol, bs, method, with_gpu))
                                 yield (eigsh_check, eigsh, A0, None, 1, k, prec, which, sigma, tol, bs, method, evals0, precision, case_desc, with_gpu)

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
            for with_gpu in ((False, True) if test_gpu else (False,)):
               for precision in (np.float16, np.float32, np.float64):
                  for which, sigma in [(w, None) for w in ('LM', 'SM', 'LA', 'SA')] + [('SM', sigma0)] :
                     # CUBLAS does not support complex half with gpus
                     if with_gpu and complexity is np.complex128 and precision is np.float16:
                        continue
                     # Normalize the problem in half precision to avoid overflows
                     if precision is np.float16:
                        A0, evals0, nA0 = A/nA, evals/nA, 1
                        sigma = sigma/nA if sigma is not None else None
                     else:
                        A0, evals0, nA0 = A, evals, nA
                     # Jacobi preconditioner is terrible for ElasticRod
                     if gen.__name__ != "ElasticRod" and sigma is not None:
                        precs = (None, jacobi_prec(stdP/nA0, sigma))
                     else:
                        precs = (None,)
                     for method in ("DEFAULT_MIN_MATVECS", "DEFAULT_MIN_TIME"):
                        for prec in precs:
                           for k in (1, 2, 3, 5, 10, 50):
                              if k > n: continue
                              tol = get_tol_eigsh(k, sigma, which, evals0, precision)
                              for bs in (1, 2, 3, 5, 10, 70):
                                 if bs*3 > n: continue
                                 case_desc = ("A,B=%s(%d, %s), k=%d, M=%s, which=%s, sigma=%s, bs=%d, method=%s, with_gpu=%s" %
                                       (gen.__name__, n, dtype_to_str(precision, complexity), k, prec is not None, which, sigma, bs, method, with_gpu))
                                 yield (eigsh_check, eigsh, A0, B, normInvB, k, prec, which, sigma, tol, bs, method, evals0, precision, case_desc, with_gpu)

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
         case_desc = ("A=%s(%d, %s), k=%d, M=%s, which=%s, sigma=%s, bs=%d, method=%s" %
                      (MikotaPair.__name__, n, dtype, k, prec is None, which, sigma, 1, 'DYNAMIC'))
         yield (eigsh_check, eigsh, op(A), None, 1, k, M, which, sigma, 1e-6, 1, 'DYNAMIC', evals, dtype, case_desc, False, False)


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
   sorted_svals = sorted(svals, key=f)
   sel_svals = np.array(sorted_svals[0:k])
   if k < n:
      gap = f(sorted_svals[k]) - f(sorted_svals[k-1])
   else:
      gap = 0
   return sel_svals, gap


def get_tol_svds(k, which, svals, precision):
    """Return a tolerance considering the spectrum"""

    # Get the smallest singular value (in absolute value) that is going to converge
    # and the smallest distance between the converged values and the rest of
    # the spectrum, gap
    sel_svals, gap = select_pairs_svds(k, which, svals)
    smallest_sval = np.fabs(sel_svals).min()

    # A good tolerance is smaller than the smallest singular value and the gap
    nA = np.fabs(svals).max()
    tol = np.min([smallest_sval*.1, gap])/nA

    # Limit the tolerance to be larger than the used precision's machine epsilon 
    return np.max([tol, np.finfo(precision).eps**.8])
 
def svds_check(svds_solver, A, k, M, which, tol, bs, methodStage1, exact_svals, dtype,
               case_desc, with_gpu=False, add_stats=True):
   """
   Test svds
   """

   if TEST_THIS_CASE and TEST_THIS_CASE != case_desc:
      return

   if TEST_THIS_CASE:
      if M is None: M = {}
      M['printLevel'] = 5

   import sys
   sys.stdout.write("Case %s\n" % case_desc)
   sys.stdout.flush()

   f = togpu if with_gpu else lambda x:x 

   if with_gpu and M:
      M = M.copy()
      if 'precAHA' in M:
         M['precAHA'] = togpu(M['precAHA'])
      if 'precAAH' in M:
         M['precAAH'] = togpu(M['precAAH'])

   try:
      svl, sva, svr, stats = svds_solver(f(A), k, None, which=which, tol=tol, reltol=.1,
            maxBlockSize=bs, methodStage1=methodStage1,
            maxMatvecs=3000000, return_stats=True,
            internalPrecision=to_primme_datatype(dtype),
            driver="numpy" if not with_gpu else "cupy", **M)
   except Exception as e:
      raise Exception("Ups! Case %s\n%s" % (case_desc, e))
   sol_svals, _ = select_pairs_svds(k, which, exact_svals)
   svl = tocpu(svl)
   svr = tocpu(svr).T.conj()

   # Check the residual norm associated to the returned pairs
   ANorm = np.amax(exact_svals)
   Rnorms = np.sqrt(
      np.linalg.norm(A.dot(svr) - svl.dot(np.diag(sva)), axis=0)**2 +
      np.linalg.norm(A.T.conj().dot(svl) - svr.dot(np.diag(sva)), axis=0)**2)

   eps = np.finfo(dtype).eps
   assert_allclose(Rnorms, np.zeros(k), atol=ANorm*(tol+eps)*max(k**.5,4), rtol=1, err_msg=case_desc)

   # Check singular values are close enough to the exact ones
   assert_allclose(sorted(sva), sorted(sol_svals), atol=ANorm*tol, rtol=1, err_msg=case_desc)

   # Add stats
   if add_stats:
      st.add("svds: " + case_desc, ("svds",), mv=stats['numMatvecs'], time=stats['elapsedTime'])

def test_primme_svds():
   """
   Generate all test cases for primme.svds.
   """

   for n in (2, 3, 5, 10, 50, 100, 200):
      for precision in (np.float64, np.float32, np.float16):
         c = np.finfo(precision).eps**.333
         for gen_name, gen in (("MikotaPair", (lambda n, d: toStandardProblem(MikotaPair(n, dtype=d)))),
                            ("Lauchli_like_vert", (lambda n, d: Lauchli_like(n*2, n, c, dtype=d))),
                            ("Lauchli_like_hori", (lambda n, d: Lauchli_like(n, n*2, c, dtype=d)))):
            sva = np.linalg.svd(gen(n, np.float64), full_matrices=False, compute_uv=False)
            nA = np.max(np.fabs(sva))
            sigma0 = sva[0]*.51 + sva[-1]*.49
            for complexity in (np.float64, np.complex128):
               A = gen(n, complexity)
               #tol = np.finfo(precision).eps**.5 * 0.1
               for with_gpu in ((False, True) if test_gpu else (False,)):
                  if with_gpu and complexity is np.complex128 and precision is np.float16:
                      continue
                  for which, sigma in [('LM', None), ('SM', 0), (sigma0, sigma0)]:
                      if precision is np.float16:
                          A0, sva0 = A/nA, sva/nA
                          if sigma is not None:
                             sigma = sigma/nA
                          if which not in frozenset(('LM', 'SM')):
                             which = which/nA
                      else:
                          A0, sva0 = A, sva
                      for prec in (({},) if which == 'LM' else ({}, sqr_diagonal_prec(A0, sigma))):
                        for k in (1, 2, 3, 5, 10, 15):
                           if k > n: continue
                           tol = get_tol_svds(k, which, sva0, precision)
                           for method in ("DEFAULT_MIN_MATVECS", "DEFAULT_MIN_TIME"):
                              for bs in (1, 2, 3, 5, 10, 70):
                                 prec0 = prec.copy()
                                 if gen_name == "MikotaPair" and n >= 50 and precision is np.float16 and which != 'LM':
                                    if bs >= 2: continue
                                    #if prec0:
                                    #   prec0['primmeStage2'] = {'projection_projection': 'primme_proj_RR'}
                                 if which != 'LM' and which != 'SM' and bs > 1:
                                    continue
                                 # if gen_name == "Lauchli_like_hori" and n >= 200 and not bool(prec) and which == sigma and k >= 15 and bs >= 2:
                                 #    continue
                                 if bs*3 > n: break
                                 case_desc = ("A=%s(%d, %s), k=%d, M=%s, which=%s, tol=%s, bs=%d, method=%s, with_gpu=%s" %
                                       (gen_name, n, dtype_to_str(precision, complexity), k, bool(prec), which, tol, bs, method, with_gpu))
                                 yield (svds_check, svds, A0, k, prec0, which, tol, bs, method, sva0, precision, case_desc, with_gpu)
                                 if bs > k: break

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
         yield (svds_check, svds, op(A), k, prec, which, 1e-5, 1, 'DYNAMIC', sva, dtype, case_desc, False, False)

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
