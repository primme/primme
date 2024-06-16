PRIMME Matlab Comparison
========================
All tests run on a core i7-9700k with 32GB of RAM available

Eigensolver Comparison
----------------------

PRIMME eigsh is based on Davidson-type methods and they may be faster than Lanczos/Arnoldi based method in difficult problems that eigenpairs take many iterations to converge or an efficient preconditioner is available::

    A = sparse(diag(1:12000) + diag(ones(11999,1),1) + diag(ones(11999,1), -1));
    ops = struct();
    ops.tol =1e-5;
    [evecs, evals, rnorms, stats] = primme_eigs(A, 2, 'LA', ops);
    [evecs, evals] = eigs(A, 2, 'largestreal', ops);

======  =======  ========  ========
test       time    matvec     rnorm
======  =======  ========  ========
Primme    0.142       526     0.117 
eigs      0.215       ---     0.116
======  =======  ========  ========

PRIMME Preconditioned Comparison
--------------------------------

By default PRIMME tries to guess the best configuration, but a little hint can help sometimes. Convergence can be accelerated by providing PRIMME with a preconditioner::
    
    [evecs, evals, rnorms, stats] = primme_eigs(A, 2, 'SM', ops);
    
    diagInv = diag((1./(1:12000)), 0);
    [evecs, evals, rnorms, stats] = primme_eigs(A, 2, 'SM', ops, [], diagInv);

===========  ========  ========  ========
test             time    matvec     rnorm
===========  ========  ========  ========
Primme          0.184       376     0.118
Primme Prec     0.452        29     0.045
===========  ========  ========  ========

Small Singular Value Comparison
-------------------------------

For SVD problems, the package provides a similar interface
PRIMME svds may perform as good as similar methods in other packages. PRIMME can take advantage of a light matrix-vector product::
 
    A = randn(6000, 6000);   
    [u, svals, v, rnorms] = primme_svds(A, 2, 'LA', ops);
    [u, svals, v] = svds(A, 2, 'largest', ops);

======  =======  ========  ========
test       time    matvec     rnorm
======  =======  ========  ========
Primme    2.096       246     0.002 
svds      4.259       ---     0.001
======  =======  ========  ========

Large Singular Value Comparison
-------------------------------

PRIMME performs similarly for larger problems::

    [u, svals, v, rnorms] = primme_svds(A, 100, 'LA', ops);
    [u, svals, v] = svds(A, 100, 'largest', ops);


======  =======  ========  ========
test       time    matvec     rnorm
======  =======  ========  ========
Primme    6.291      6612     0.002 
svds     20.919       ---     0.000
======  =======  ========  ========
