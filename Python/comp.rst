PRIMME Python comparison
========================
All tests were run on a Ryzen 2700X with 32GB of RAM available

Eigensolver Comparison
----------------------

PRIMME eigsh is based on Davidson-type methods and they may be faster than Lanczos/Arnoldi based method in difficult problems that eigenpairs take many iterations to converge or an efficient preconditioner is available::
    
    A = np.diag(np.asarray(range(1,12001), dtype=np.float32), 0)
    for i in range(1,11999):
        A[i,i+1] = 1
        A[i+1,i] = 1    
    evals, evecs, stats = primme.eigsh(A, 2, tol=tolerance, return_stats =True, return_history=True)
    evals, evecs = scipy.sparse.linalg.lobpcg(A, X=X, tol=tolerance*normA, largest=True, maxiter=1000)

=======================  ========  ========  ========
test                         time  matvec       rnorm
=======================  ========  ========  ========
Primme                    50.8575  503       0.118076
sp.sparse.linalg.lobpcg  124.044   ---       0.119619
=======================  ========  ========  ========

PRIMME Preconditioned Comparison
--------------------------------

By default PRIMME tries to guess the best configuration, but a little hint can help sometimes. Convergence can be accelerated by providing PRIMME with a preconditioner::
    
    diagInv = np.diag(np.reciprocal(np.asarray(range(1,12001), dtype=np.float32)), 0)
    evals, evecs, stats = primme.eigsh(A, 2, OPinv= diagInv, tol=tolerance, return_stats =True, which= 'SM', return_history=True)

===========  ========  ========  ========
test             time    matvec     rnorm
===========  ========  ========  ========
Primme       37.7227        450  0.11965
Primme Prec   1.89826        29  0.108924
===========  ========  ========  ======== 

Small Singular Value Comparison
-------------------------------
For SVD problems, the package provides a similar interface::

    u, sdvals, vt, stats = primme.svds(A, 2, tol=tolerance, return_stats=True)
    
PRIMME svds may perform as good as similar methods in the packages scipy.linalg.sparse and irlb in solving few singular values. PRIMME can take advantage of a light matrix-vector product::
    
    A = np.random.normal(size=(6000,6000))    
    u, sdvals, vt, stats = primme.svds(A, 2, tol=tolerance, maxBlockSize = 2, return_stats=True)
    u, sdvals, vt = scipy.sparse.linalg.svds(A, 2, tol=tolerance*normA, return_singular_vectors=True)
    u, sdvals, v, it, mprod  = irlb.irlb(A, 2, tol=tolerance, maxit = 1000)
    svdals = sklearn.decomposition.TruncatedSVD(n_components=2, tol=tolerance, n_iter = 1000)
    svdals.fit(A)

=====================  =========  ========  =====================
test                        time  matvec    rnorm
=====================  =========  ========  =====================
Primme                   5.21035  320       0.0014542283007373514
sp.sparse.linalg.svds    8.17256  ---       1.473394835276466e-09
irlb                     6.68747  ---       0.0011971646391146895
sklearn                141.192    ---       ---
=====================  =========  ========  ===================== 

Large Singular Value Comparison
-------------------------------

PRIMME performs similarly for larger problems::
    
    A = np.random.normal(size=(6000,6000))    
    u, sdvals, vt, stats = primme.svds(A, 100, tol=tolerance, maxBloackSize = 100, return_stats=True)
    u, sdvals, vt = scipy.sparse.linalg.svds(A, 100, tol=tolerance*normA, return_singular_vectors=True)
    u, sdvals, v, it, mprod  = irlb.irlb(A, 100, tol=tolerance, maxit = 1000)

=====================  =======  ========  ==========
test                      time  matvec         rnorm
=====================  =======  ========  ==========
Primme                 26.944   5672      0.00212836
sp.sparse.linalg.svds  33.9405  ---       7.8088e-11
irlb                   41.0081  ---       0.00132166
=====================  =======  ========  ========== 
