.. highlight:: matlab

MATLAB Interface
----------------

.. mat:function:: function [varargout] = primme_eigs(varargin)

   :mat:func:`primme_eigs` finds a few eigenvalues and their corresponding eigenvectors 
   of a real symmetric or Hermitian matrix, ``A``, by calling PRIMME_.

   ``D = primme_eigs(A)`` returns a vector of ``A``'s 6 largest magnitude eigenvalues.

   ``D = primme_eigs(Afun,dim)`` accepts a function ``Afun`` instead of a matrix. ``Afun``
   is a function handle and ``y = Afun(x)`` returns the matrix-vector product ``A*x``.
   In all the following syntaxes, ``A`` can be replaced by ``Afun, dim``.

   ``D = primme_eigs(A,k)`` finds the ``k`` largest magnitude eigenvalues. ``k`` must be
   less than the dimension of the matrix ``A``.

   ``D = primme_eigs(A,k,target)`` returns ``k`` eigenvalues such that: 
   If ``target`` is a real number, it finds the closest eigenvalues to ``target``.
   If ``target`` is

      * ``'LA'`` or ``'SA'``, eigenvalues with the largest or smallest algebraic value.

      * ``'LM'`` or ``'SM'``, eigenvalues with the largest or smallest magnitude if
        ``OPTS.targetShifts`` is empty. If ``target`` is a real or complex 
        scalar including 0, :mat:func:`primme_eigs` finds the eigenvalues closest 
        to ``target``.

        In addition, if some values are provided in ``OPTS.targetShifts``,
        it finds eigenvalues that are farthest (``'LM'``) or closest (``'SM'``) in 
        absolute value from the given values.

        Examples:

        ``k=1``, ``'LM'``, ``OPTS.targetShifts=[]`` returns the largest magnitude ``eig(A)``.
        ``k=1``, ``'SM'``, ``OPTS.targetShifts=[]`` returns the smallest magnitude ``eig(A)``.
        ``k=3``, ``'SM'``, ``OPTS.targetShifts=[2, 5]`` returns the closest eigenvalue in 
        absolute sense to 2, and the two closest eigenvalues to 5.

      * ``'CLT'`` or ``'CGT'``, find eigenvalues closest to but less or greater than
        the given values in ``OPTS.targetShifts``.

   ``D = primme_eigs(A,k,target,OPTS)`` specifies extra solver parameters. Some
   default values are indicated in brackets {}:

      * |aNorm|: the estimated 2-norm of A {0.0 (estimate the norm internally)}
      * ``tol``: convergence tolerance: ``NORM(A*X(:,i)-X(:,i)*D(i,i)) < tol*NORM(A)``
        (see |eps|) {:math:`10^4` times the machine precision}
      * |maxBlockSize|: maximum block size (useful for high multiplicities) {1}
      * ``disp``: different level reporting (0-3) (see HIST) {no output 0}
      * ``isreal``: whether A represented by ``Afun`` is real or complex {false}
      * |targetShifts|: shifts for interior eigenvalues (see ``target``) {[]}
      * ``v0``: any number of initial guesses to the eigenvectors (see |initSize| {[]}
      * ``orthoConst``: external orthogonalization constraints (see |numOrthoConst| {[]}
      * |locking|: 1, hard locking; 0, soft locking
      * ``p``: maximum size of the search subspace (see |maxBasisSize|)
      * |minRestartSize|: minimum Ritz vectors to keep in restarting
      * |maxMatvecs|: maximum number of matrix vector multiplications {Inf}
      * ``maxit``: maximum number of outer iterations (see |maxOuterIterations|) {Inf}
      * |scheme|: the restart scheme {'primme_thick'}
      * |maxPrevRetain|: number of Ritz vectors from previous iteration that are kept after restart {typically >0}
      * |robustShifts|: setting to true may avoid stagnation or misconvergence 
      * |maxInnerIterations|: maximum number of inner solver iterations
      * |LeftQ|: use the locked vectors in the left projector
      * |LeftX|: use the approx. eigenvector in the left projector
      * |RightQ|: use the locked vectors in the right projector
      * |RightX|: use the approx. eigenvector in the right projector
      * |SkewQ|: use the preconditioned locked vectors in the right projector
      * |SkewX|: use the preconditioned approx. eigenvector in the right projector
      * |relTolBase|: a legacy from classical JDQR (not recommended)
      * |convTest|: how to stop the inner QMR Method
      * |iseed|: random seed

   ``D = primme_eigs(A,k,target,OPTS,METHOD)`` specifies the eigensolver method.
   METHOD can be one of the next strings:

      * '|DYNAMIC|', (default)        switches dynamically to the best method
      * '|DEFAULT_MIN_TIME|',         best method for low-cost matrix-vector product
      * '|DEFAULT_MIN_MATVECS|',      best method for heavy matvec/preconditioner
      * '|Arnoldi|',                  Arnoldi not implemented efficiently
      * '|GD|',                       classical block Generalized Davidson 
      * '|GD_plusK|',                 GD+k block GD with recurrence restarting
      * '|GD_Olsen_plusK|',           GD+k with approximate Olsen precond.
      * '|JD_Olsen_plusK|',           GD+k, exact Olsen (two precond per step)
      * '|RQI|',                      Rayleigh Quotient Iteration. Also INVIT, but for INVIT provide OPTS.targetShifts
      * '|JDQR|',                     Original block, Jacobi Davidson
      * '|JDQMR|',                    Our block JDQMR method (similar to JDCG)
      * '|JDQMR_ETol|',               Slight, but efficient JDQMR modification
      * '|STEEPEST_DESCENT|',         equivalent to GD(block,2*block)
      * '|LOBPCG_OrthoBasis|',        equivalent to GD(nev,3*nev)+nev
      * '|LOBPCG_OrthoBasis_Window|'  equivalent to GD(block,3*block)+block nev>block

   ``D = primme_eigs(A,k,target,OPTS,METHOD,P)``

   ``D = primme_eigs(A,k,target,OPTS,METHOD,P1,P2)`` uses preconditioner ``P`` or
   ``P = P1*P2`` to accelerate convergence of the method. Applying ``P\x`` should
   approximate ``(A-sigma*eye(N))\x``, for ``sigma`` near the wanted eigenvalue(s).
   If ``P`` is ``[]`` then a preconditioner is not applied. ``P`` may be a function 
   handle ``PFUN`` such that ``PFUN(x)`` returns ``P\x``.

   ``[X,D] = primme_eigs(...)`` returns a diagonal matrix ``D`` with the eigenvalues
   and a matrix ``X`` whose columns are the corresponding eigenvectors. 
 
   ``[X,D,R] = primme_eigs(...)`` also returns an array of the residual norms of
   the computed eigenpairs.

   ``[X,D,R,STATS] = primme_eigs(...)`` returns a ``struct`` to report statistical
   information about number of matvecs, elapsed time, and estimates for the
   largest and smallest algebraic eigenvalues of ``A``.

   ``[X,D,R,STATS,HIST] = primme_eigs(...)`` it returns the convergence history,
   instead of printing it. Every row is a record, and the columns report:
  
      * ``HIST(:,1)``: number of matvecs
      * ``HIST(:,2)``: time
      * ``HIST(:,3)``: number of converged/locked pairs
      * ``HIST(:,4)``: block index
      * ``HIST(:,5)``: approximate eigenvalue
      * ``HIST(:,6)``: residual norm
      * ``HIST(:,7)``: QMR residual norm

   ``OPTS.disp`` controls the granularity of the record. If ``OPTS.disp == 1``, ``HIST``
   has one row per converged eigenpair and only the first three columns are
   reported; if ``OPTS.disp == 2``, ``HIST`` has one row per outer iteration and only
   the first six columns are reported; and otherwise ``HIST`` has one row per QMR
   iteration and all columns are reported.

   Examples::

      A = diag(1:100);

      d = primme_eigs(A,10) % the 10 largest magnitude eigenvalues

      d = primme_eigs(A,10,'SM') % the 10 smallest magnitude eigenvalues

      d = primme_eigs(A,10,25.0) % the 10 closest eigenvalues to 25.0

      opts.targetShifts = [2 20];
      d = primme_eigs(A,10,'SM',opts) % 1 eigenvalue closest to 2 and 
                                      % 9 eigenvalues closest to 20

      opts = struct();
      opts.tol = 1e-4; % set tolerance
      opts.maxBlockSize = 2; % set block size
      [x,d] = primme_eigs(A,10,'SA',opts,'DEFAULT_MIN_TIME')

      opts.orthoConst = x;  
      [d,rnorms] = primme_eigs(A,10,'SA',opts) % find another 10 with the default method

      % Compute the 6 eigenvalues closest to 30.5 using ILU(0) as a preconditioner
      % by passing the matrices L and U.
      A = sparse(diag(1:50) + diag(ones(49,1), 1) + diag(ones(49,1), -1));
      [L,U] = ilu(A, struct('type', 'nofill'));
      d = primme_eigs(A, k, 30.5, [], [], L, U);

      % Compute the 6 eigenvalues closest to 30.5 using Jacobi preconditioner
      % by passing a function.
      Pfun = @(x)(diag(A) - 30.5)\x;
      d = primme_eigs(A,6,30.5,[],[],Pfun) % find the closest 5 to 30.5

   See also: `MATLAB eigs`_, :mat:func:`primme_svds`

.. include:: epilog.inc
