# ******************************************************************************
#    PRIMME PReconditioned Iterative MultiMethod Eigensolver
#    Copyright (C) 2015 College of William & Mary,
#    James R. McCombs, Eloy Romero Alcalde, Andreas Stathopoulos, Lingfei Wu
# 
#    This file is part of PRIMME.
# 
#    PRIMME is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
# 
#    PRIMME is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
# 
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# 
# ******************************************************************************
#  File: examples.py
#  
#  Purpose - Simple examples using the PYTHON interface to PRIMME.
#  
# ******************************************************************************

import numpy as np
import scipy.sparse
import Primme


# Sparse diagonal matrix of size 100
A = scipy.sparse.spdiags(range(100), [0], 100, 100)

# Compute the three largest eigenvalues of A with a residual norm tolerance of 1e-6
evals, evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA')
print evals # [ 99.,  98.,  97.]

# Compute the three largest eigenvalues of A orthogonal to the previous computed
# eigenvectors, i.e., the next three eigenvalues
evals, evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA', lock=evecs)
print evals # [ 96.,  95.,  94.]

# Sparse rectangular matrix 100x10 with non-zeros on the main diagonal
A = scipy.sparse.spdiags(range(10), [0], 100, 10)

# Compute the three closest to 4.1 singular values and the left and right corresponding
# singular vectors
svecs_left, svals, svecs_right = Primme.svds(A, 3, tol=1e-6, which=4.1)
print svals # [ 4.,  5.,  3.]

# Sparse random rectangular matrix 10^5x100
A = scipy.sparse.rand(10000, 100, density=0.001, random_state=10)

# Compute the three closest singular values to 6.0 with a tolerance of 1e-6
svecs_left, svals, svecs_right, stats = Primme.svds(A, 3, which='SM', tol=1e-6,
                                                    return_stats=True)
print svals # [ 0.79488437  0.85890809  0.87174328]
print stats["elapsedTime"], stats["numMatvecs"] # it took that seconds and 101 matvecs

# Compute the square diagonal preconditioner
prec = scipy.sparse.spdiags(np.reciprocal(A.multiply(A).sum(axis=0)),
          [0], 100, 100)

# Recompute the singular values but using the preconditioner
svecs_left, svals, svecs_right, stats = Primme.svds(A, 3, which='SM', tol=1e-6,
                        precAHA=prec, return_stats=True)
print stats["elapsedTime"], stats["numMatvecs"] # it took that seconds and 45 matvecs
