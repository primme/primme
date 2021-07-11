# !pip install scipy numpy primme sklearn
# !pip install -e git+https://github.com/bwlewis/irlbpy.git#egg=irlbpy

import primme, irlb, timeit, sklearn.decomposition, scipy.sparse.linalg, numpy as np
from tabulate import tabulate

#------------------------
#EigenValue Tests
#------------------------

table = list()
tolerance = 1e-5
mvCounter = 0
A = np.diag(np.asarray(range(1,12001), dtype=np.float32), 0)
for i in range(1,11999):
  A[i,i+1] = 1
  A[i+1,i] = 1

def eigMatVec(v):
    global mvCounter
    mvCounter +=1
    return np.matmul(A,v)

lin = scipy.sparse.linalg.LinearOperator((12000, 12000), matvec=eigMatVec, matmat=eigMatVec)
normA = 12000
headers = ["test", "time", "matvec", "rnorm"]

#PRIMME performance test
start = timeit.default_timer()
evals, evecs, stats = primme.eigsh(A, 2, tol=tolerance, return_stats =True, return_history=True)
end = timeit.default_timer()
table.append(["Primme", end - start, stats["numMatvecs"], max(stats['rnorms'])])

mvCounter = 0
#scipy.sparse.linalg.lobpcg performance test
y = np.eye(12000, 2)
rng = np.random.default_rng()
X = rng.random((12000, 2))
start = timeit.default_timer()
evals, evecs = scipy.sparse.linalg.lobpcg(A, X=X, tol=tolerance*normA, largest=True, maxiter=1000)
end = timeit.default_timer()
rnorm = np.linalg.norm(np.matmul(A, evecs) - (evecs.dot(np.diag(evals))), axis=0)
table.append(["sp.sparse.linalg.lobpcg", end - start, "---", max(rnorm)])
print(tabulate(table, headers=headers, tablefmt="rst"), "\n")

#------------------------------------------
#EigenValue Preconditioner comparison Tests
#------------------------------------------

table = list()
#PRIMME performance test
start = timeit.default_timer()
evals, evecs, stats = primme.eigsh(A, 2, which= 'SM', tol=tolerance, return_stats =True, return_history=True)
end = timeit.default_timer()
table.append(["Primme", end - start, stats["numMatvecs"], max(stats['rnorms'])])


diagInv = np.diag(np.reciprocal(np.asarray(range(1,12001), dtype=np.float32)), 0)
#PRIMME Preconditioned performance test
start = timeit.default_timer()
evals, evecs, stats = primme.eigsh(A, 2, OPinv= diagInv, tol=tolerance, return_stats =True, which= 'SM', return_history=True)
end = timeit.default_timer()
table.append(["Primme Prec", end - start, stats["numMatvecs"], max(stats['rnorms'])])
print(tabulate(table, headers=headers, tablefmt= "rst"), "\n")

#------------------------
#Singular Value Tests
#------------------------

A = np.random.normal(size=(6000,6000))
normA = 154.8
lin = scipy.sparse.linalg.LinearOperator(shape = (6000, 6000), matvec=eigMatVec, matmat=eigMatVec)
k = [2,100]
for n in k:
  table = list()
  #PRIMME performance test
  start = timeit.default_timer()
  u, sdvals, vt, stats = primme.svds(A, n, tol=tolerance, return_stats=True)
  end = timeit.default_timer()
  table.append(["Primme", end - start, stats["numMatvecs"], max(stats['rnorms'])])
  
  #scipy.sparse.linalg.svds performance test
  start = timeit.default_timer()
  u, sdvals, vt = scipy.sparse.linalg.svds(A, n, tol=tolerance*normA, return_singular_vectors=True)
  end = timeit.default_timer()
  rnorm = np.linalg.norm(A.dot(vt.T.conj()) - u.dot(np.diag(sdvals)),axis=0)**2+np.linalg.norm( u.T.conj().dot(A).T.conj() -(vt.T.conj()).dot(np.diag(sdvals)), axis=0)**2
  table.append(["sp.sparse.linalg.svds", end - start, "---", max(rnorm)])
  
  #irlb performance test
  start = timeit.default_timer()
  u, sdvals, v, it, mprod  = irlb.irlb(A, n, tol=tolerance, maxit = 1000)
  end = timeit.default_timer()
  rnorm = np.linalg.norm(A.dot(v) - u.dot(np.diag(sdvals)),axis=0)**2+np.linalg.norm( u.T.conj().dot(A).T.conj() -v.dot(np.diag(sdvals)), axis=0)**2
  table.append(["irlb", end - start, "---", max(np.sqrt(rnorm))])
  
  if n == 2:
    #sklearn.decomposition permormance test
    start = timeit.default_timer()
    svdals = sklearn.decomposition.TruncatedSVD(n_components= n, tol=tolerance, n_iter = 1000)
    svdals.fit(A)
    end = timeit.default_timer()
    table.append(["sklearn", end - start, '---', '---'])
  print(tabulate(table, headers=headers, tablefmt="rst"), "\n")


