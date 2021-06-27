# !pip install scipy numpy primme sklearn
# !pip install -e git+https://github.com/bwlewis/irlbpy.git#egg=irlbpy
# !pip install -e git+https://github.com/jakevdp/pypropack.git#egg=pypropack

import primme, irlb, timeit, sklearn.decomposition, scipy.sparse, numpy as np
#from pypropack import svdp

#------------------------
#EigenValue Tests
#------------------------

A = np.diag(np.asarray(range(1,12001), dtype=np.float32), 0)
for i in range(1,11999):
  A[i,i+1] = 1
  A[i+1,i] = 1

normA = 12000

#PRIMME performance test
start = timeit.default_timer()
evals, evecs, stats = primme.eigsh(A, 2, tol=1e-5, return_stats =True, return_history=True)
end = timeit.default_timer()
print("Primme:", end - start, "rnorm:", stats['rnorms'], "Eigenvalues:", evals)

#scipy.sparse.linalg.lobpcg performance test
y = np.eye(12000, 2)
rng = np.random.default_rng()
X = rng.random((12000, 2))
start = timeit.default_timer()
evals, evecs = scipy.sparse.linalg.lobpcg(A, X=X, tol=1e-5*normA, largest=True, maxiter=1000)
end = timeit.default_timer()
rnorm = np.linalg.norm(np.matmul(A, evecs) - (evecs.dot(np.diag(evals))), axis=0)
print("scipy.sparse.linalg.lobpcg:", end - start, "rnorm:", rnorm, "Eigenvalues:", evals, "\n")

#------------------------
#Singular Value Tests
#------------------------

A = np.random.normal(size=(6000,6000))

#PRIMME performance test
start = timeit.default_timer()
u, sdvals, vt, stats = primme.svds(A, 2, tol=1e-5, return_stats=True)
end = timeit.default_timer()
print("Primme:", end - start, "rnorm:", stats['rnorms'], "Eigenvalues:", sdvals)
#print(np.linalg.norm(A.dot(vt.T.conj()) - u.dot(np.diag(sdvals)), axis=0))

#scipy.sparse.linalg.svds performance test
start = timeit.default_timer()
u, sdvals, vt = scipy.sparse.linalg.svds(A, 2, tol=1e-5, return_singular_vectors=True)
end = timeit.default_timer()
rnorm = np.linalg.norm(A.dot(vt.T.conj()) - u.dot(np.diag(sdvals)), axis=0)
print("scipy.sparse.linalg.svds:", end - start, "rnorm:", rnorm, "Eigenvalues:", sdvals, "\n")

#irlb performance test
start = timeit.default_timer()
sdvals = irlb.irlb(A, 2, tol=1e-5)
end = timeit.default_timer()
u = sdvals[0]
vt = sdvals[2]
sdvals = sdvals[1]
rnorm = np.linalg.norm(A.dot(vt.conj()) - u.dot(np.diag(sdvals)), axis=0)
print("irlb:", end - start, "rnorm:", rnorm, "Eigenvalues:", sdvals, "\n")

#sklearn.decomposition permormance test
start = timeit.default_timer()
svdals = sklearn.decomposition.TruncatedSVD(n_components= 2, tol=1e-5, n_iter = 100)
svdals.fit(A)
end = timeit.default_timer()
#rnorm = np.linalg.norm(A.dot(vt.T.conj()) - u.dot(np.diag(sdvals)), axis=0)
print("sklearn:", end - start, "rnorm:", '---', "Eigenvalues:", sdvals, "\n")


