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
#  File: test1.py
#  
#  Purpose - Simple test for PYTHON interface to PRIMME.
#  
# ******************************************************************************

import Primme, numpy as np
from scipy.sparse import *

# A = [ 2  1  0 ...
#      -1  2 -1 0 ...
#       0 -1  2 -1 0 ... ]
a = np.ones(10)
A = spdiags(np.array([a*(-1.), a*2., a*(-1.)]), np.array([-1, 0, 1]), 10, 10)

class PP(Primme.primme_params_w):
	def __init__(self):
		Primme.primme_params_w.__init__(self)
	def matvec(self, X):
		return A*X
pp = PP()
pp.n = A.shape[0]
#pp.maxBasisSize = 3
#pp.minRestartSize = 1
pp.numEvals = 3
#pp.restartingParams.maxPrevRetain = 1

pp.set_method(Primme.DYNAMIC)
pp.display()
evals = np.zeros(pp.numEvals)
evecs = np.zeros((pp.n, pp.numEvals))
norms = np.zeros(pp.numEvals)
print Primme.dprimme(evals, evecs, norms, pp)
print pp.initSize, evals, norms

class PPs(Primme.primme_params_w):
	def __init__(self, matrix=None):
		Primme.primme_params_w.__init__(self)
		self.mymatrix = matrix
	def matvec(self, X):
		return self.mymatrix*X

a = np.ones(10, complex)
A = spdiags(np.array([a*(-1.), a*2., a*(-1.)]), np.array([-1, 0, 1]), 10, 10)
pp = PPs(A)
pp.n = A.shape[0]
#pp.maxBasisSize = 3
#pp.minRestartSize = 1
pp.numEvals = 3
pp.set_method(Primme.DYNAMIC)
#pp.display()
evals = np.zeros(pp.numEvals)
evecs = np.zeros((pp.n, pp.numEvals), complex)
norms = np.zeros(pp.numEvals)
print Primme.zprimme(evals, evecs, norms, pp)
print pp.initSize, evals, norms, pp.stats.numMatvecs
