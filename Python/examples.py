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
import scipy.sparse
import Primme


# Sparse diagonal matrix of size 100
A = scipy.sparse.spdiags(range(100), [0], 100, 100)

# Compute the three largest eigenvalues of A with a residual norm tolerance of 1e-6
evals, evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA')
print(evals) # [ 99.,  98.,  97.]

# Compute the three largest eigenvalues of A orthogonal to the previous computed
# eigenvectors, i.e., the next three eigenvalues
evals, evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA', lock=evecs)
print(evals) # [ 96.,  95.,  94.]

# Sparse rectangular matrix 100x10 with non-zeros on the main diagonal
A = scipy.sparse.spdiags(range(10), [0], 100, 10)

# Compute the three closest to 4.1 singular values and the left and right corresponding
# singular vectors
svecs_left, svals, svecs_right = Primme.svds(A, 3, tol=1e-6, which=4.1)
print(svals) # [ 4.,  5.,  3.]

# Sparse random rectangular matrix 10^5x100
A = scipy.sparse.rand(10000, 100, density=0.001, random_state=10)

# Compute the three closest singular values to 6.0 with a tolerance of 1e-6
svecs_left, svals, svecs_right, stats = Primme.svds(A, 3, which='SM', tol=1e-6,
                                                    return_stats=True)
print(svals) # [ 0.79488437  0.85890809  0.87174328]
print(stats["elapsedTime"], stats["numMatvecs"]) # it took that seconds and 101 matvecs

# Compute the square diagonal preconditioner
prec = scipy.sparse.spdiags(np.reciprocal(A.multiply(A).sum(axis=0)),
          [0], 100, 100)

# Recompute the singular values but using the preconditioner
svecs_left, svals, svecs_right, stats = Primme.svds(A, 3, which='SM', tol=1e-6,
                        precAHA=prec, return_stats=True)
print(stats["elapsedTime"], stats["numMatvecs"]) # it took that seconds and 45 matvecs
