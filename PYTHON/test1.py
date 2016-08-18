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

a = np.ones(10, complex)
Az = spdiags(np.array([a*(-1.), a*2., a*(-1.)]), np.array([-1, 0, 1]), 10, 10)

#
# High level tests
#

print Primme.eigsh(A, k=3, tol=1e-6, which='SA')
print Primme.eigsh(Az, k=3, tol=1e-6, which='SA')
print Primme.svds(A, k=3, tol=1e-6, which='SM')
print Primme.svds(Az, k=3, tol=1e-6, which='SM')

#
# Low level tests
#

class PP(Primme.PrimmeParams):
	def __init__(self):
		Primme.PrimmeParams.__init__(self)
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

class PPc(Primme.PrimmeParams):
	def __init__(self, matrix=None):
		Primme.PrimmeParams.__init__(self)
		self.mymatrix = matrix
	def matvec(self, X):
		return self.mymatrix*X

pp = PPc(Az)
pp.n = Az.shape[0]
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

class PSP(Primme.PrimmeSvdsParams):
	def __init__(self):
		Primme.PrimmeSvdsParams.__init__(self)
		self._At = A.T
	def matvec(self, X, transpose):
		if not transpose:
			return A*X
		else:
			return self._At*X

pp = PSP()
pp.m = A.shape[0]
pp.n = A.shape[1]
#pp.maxBasisSize = 3
#pp.minRestartSize = 1
pp.numSvals = 3
pp.eps = 1e-6
pp.target = Primme.primme_svds_largest

pp.set_method(Primme.primme_svds_default, Primme.DEFAULT_METHOD, Primme.DEFAULT_METHOD)
pp.display()
svals = np.zeros(pp.numSvals)
svecsl = np.zeros((pp.m, pp.numSvals))
svecsr = np.zeros((pp.n, pp.numSvals))
norms = np.zeros(pp.numSvals)
print Primme.dprimme_svds(svals, svecsl, svecsr, norms, pp)
print pp.initSize, svals, norms

class PSPc(Primme.PrimmeSvdsParams):
	def __init__(self):
		Primme.PrimmeSvdsParams.__init__(self)
		self._At = Az.T.conj()
	def matvec(self, X, transpose):
		if not transpose:
			return Az*X
		else:
			return self._At*X
pp = PSPc()
pp.m = Az.shape[0]
pp.n = Az.shape[1]
#pp.maxBasisSize = 3
#pp.minRestartSize = 1
pp.numSvals = 3
pp.eps = 1e-6

pp.set_method(Primme.primme_svds_default, Primme.DEFAULT_METHOD, Primme.DEFAULT_METHOD)

svals = np.zeros(pp.numSvals)
svecsl = np.zeros((pp.m, pp.numSvals), Az.dtype)
svecsr = np.zeros((pp.n, pp.numSvals), Az.dtype)
norms = np.zeros(pp.numSvals)
print Primme.zprimme_svds(svals, svecsl, svecsr, norms, pp)
print pp.initSize, svals, norms
