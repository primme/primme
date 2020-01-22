#!/usr/bin/env python

#  Copyright (c) 2016, College of William & Mary
#  All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met:
#      * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.
#      * Redistributions in binary form must reproduce the above copyright
#        notice, this list of conditions and the following disclaimer in the
#        documentation and/or other materials provided with the distribution.
#      * Neither the name of College of William & Mary nor the
#        names of its contributors may be used to endorse or promote products
#        derived from this software without specific prior written permission.
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
#  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
#  PRIMME: https://github.com/primme/primme
#  Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u

import numpy as np
from numpy.testing import assert_allclose
import scipy.sparse
import primme


# Sparse diagonal matrix of size 100
A = scipy.sparse.spdiags(np.asarray(range(100), dtype=np.float32), [0], 100, 100)

# Compute the three largest eigenvalues of A with a residual norm tolerance of 1e-6
evals, evecs = primme.eigsh(A, 3, tol=1e-6, which='LA')
assert_allclose(evals, [ 99.,  98.,  97.], atol=1e-6*100)
print(evals) # [ 99.,  98.,  97.]

# Compute the three largest eigenvalues of A orthogonal to the previous computed
# eigenvectors, i.e., the next three eigenvalues
evals, evecs = primme.eigsh(A, 3, tol=1e-6, which='LA', lock=evecs)
assert_allclose(evals, [ 96.,  95.,  94.], atol=1e-6*100)
print(evals) # [ 96.,  95.,  94.]

# Compute the three closest eigenvalues to 50.1
evals, evecs = primme.eigsh(A, 3, tol=1e-6, which=50.1)
assert_allclose(evals, [ 50.,  51.,  49.], atol=1e-6*100)
print(evals) # [ 50.,  51.,  49.]

# Estimation of the largest eigenvalue in magnitude
def convtest_lm(eval, evecl, rnorm):
   return np.abs(eval) > 0.1 * rnorm
eval, evec = primme.eigsh(A, 1, which='LM', convtest=convtest_lm)
assert_allclose(eval, [ 99.], atol=.1)

# User-defined matvec: implicit diagonal matrix
Adiag = np.arange(0, 100).reshape((100,1))
def Amatmat(x):
   if len(x.shape) == 1: x = x.reshape((100,1))
   return Adiag * x   # equivalent to diag(Adiag).dot(x)
A = scipy.sparse.linalg.LinearOperator((100,100), matvec=Amatmat, matmat=Amatmat)
evals, evecs = primme.eigsh(A, 3, tol=1e-6, which='LA')
assert_allclose(evals, [ 99.,  98.,  97.], atol=1e-6*100)

# Sparse singular mass matrix
A = scipy.sparse.spdiags(np.asarray(range(100), dtype=np.float32), [0], 100, 100)
M = scipy.sparse.spdiags(np.asarray(range(99,-1,-1), dtype=np.float32), [0], 100, 100)
evals, evecs = primme.eigsh(A, 3, M=M, tol=1e-6, which='SA')
assert_allclose(evals, [ 0./99.,  1./98.,  2./97.], atol=1e-6*100)
print(evals)


# Sparse rectangular matrix 100x10 with non-zeros on the main diagonal
A = scipy.sparse.spdiags(range(10), [0], 100, 10)

# Compute the three closest to 4.1 singular values and the left and right corresponding
# singular vectors
svecs_left, svals, svecs_right = primme.svds(A, 3, tol=1e-6, which=4.1)
assert_allclose(sorted(svals), [ 3.,  4.,  5.], atol=1e-6*10)
print(svals) # [ 4.,  5.,  3.]

# Sparse random rectangular matrix 10^5x100
A = scipy.sparse.rand(10000, 100, density=0.001, random_state=10)

# Compute the three closest singular values to 6.0 with a tolerance of 1e-6
svecs_left, svals, svecs_right, stats = primme.svds(A, 3, which='SM', tol=1e-6,
                                                    return_stats=True)
A_svals = svals
print(svals)
print(stats["elapsedTime"], stats["numMatvecs"])

# Compute the square diagonal preconditioner
prec = scipy.sparse.spdiags(np.reciprocal(A.multiply(A).sum(axis=0)),
          [0], 100, 100)

# Recompute the singular values but using the preconditioner
svecs_left, svals, svecs_right, stats = primme.svds(A, 3, which='SM', tol=1e-6,
                        precAHA=prec, return_stats=True)
assert_allclose(svals, A_svals, atol=1e-6*100)
print(stats["elapsedTime"], stats["numMatvecs"])

# Estimation of the smallest singular value
def convtest_sm(sval, svecl, svecr, rnorm):
   return sval > 0.1 * rnorm
svec_left, sval, svec_right, stats = primme.svds(A, 1, which='SM',
                        convtest=convtest_sm, return_stats=True)
assert_allclose(sval, [ 1.], atol=.1)

# User-defined matvec: implicit rectangular matrix with nonzero elements on the diagonal only
Bdiag = np.arange(0, 100).reshape((100,1))
Bdiagr = np.concatenate((np.arange(0, 100).reshape((100,1)).astype(np.float32), np.zeros((100,1), dtype=np.float32)), axis=None).reshape((200,1))
def Bmatmat(x):
   if len(x.shape) == 1: x = x.reshape((100,1))
   return np.vstack((Bdiag * x, np.zeros((100, x.shape[1]), dtype=np.float32)))
def Brmatmat(x):
   if len(x.shape) == 1: x = x.reshape((200,1))
   return (Bdiagr * x)[0:100,:]

B = scipy.sparse.linalg.LinearOperator((200,100), matvec=Bmatmat, matmat=Bmatmat, rmatvec=Brmatmat, dtype=np.float32)
svecs_left, svals, svecs_right = primme.svds(B, 3, which='LM', tol=1e-6)
assert_allclose(svals, [ 99.,  98.,  97.], atol=1e-6*100)
