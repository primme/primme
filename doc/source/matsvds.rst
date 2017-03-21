.. highlight:: matlab

MATLAB Interface
----------------

.. mat:function:: function [varargout] = primme_svds(varargin)

   :mat:func:`primme_svds` find a few singular values and vectors of large,
   sparse matrices. by calling PRIMME_.

   ``S = primme_svds(A)`` computes the 6 largest singular values of ``A``.

   ``S = primme_svds(AFUN,M,N)`` accepts the function handle ``AFUN`` to perform
   the matrix vector products with an M-by-N matrix ``A``. 
   ``AFUN(X,'notransp')`` returns ``A*X`` while ``AFUN(X,'transp')`` returns ``A’*X``.
   In all the following, ``A`` can be replaced by ``AFUN,M,N``.
 
   ``S = primme_svds(A,K)`` computes the ``K`` largest singular values of ``A``.

   ``S = primme_svds(A,K,SIGMA)`` computes the K singular values closest to the
   scalar shift ``SIGMA``.

      * If ``SIGMA`` is a vector, find a singular value closest to each ``SIGMA(i)``
      * If ``SIGMA`` is ``'L'``, it computes the largest singular values.
      * if ``SIGMA`` is ``'S'``, it computes the smallest singular values.

   ``S = primme_svds(A,K,SIGMA,OPTIONS)`` specifies extra solver parameters.
   Some default values are indicated in brackets {}:

      * |SaNorm|:    estimation of the 2-norm A                    
      * ``tol``:     convergence tolerance ``NORM([A*V-U*S;A'*U-V*S]) <= tol * NORM(A)`` (see |Seps|) { ``1e-10``}
      * ``maxit``:   maximum number of iterations (see |SmaxMatvecs|)  {inf}
      * ``p``:       maximum basis size (see |SmaxBasisSize|)
      * ``disp``:    level of reporting 0-3 (see HIST) {0: no output}
      * ``isreal``:  if 0, the matrix is complex; else it's real {1: complex}
      * ``isdouble``: if 0, the matrix is single; else it's double {1: double}
      * ``method``:  which equivalent eigenproblem to solve

         * '|primme_svds_normalequations|': ``A'*A`` or ``A*A'``
         * '|primme_svds_augmented|': ``[0 A';A 0]``
         * '|primme_svds_hybrid|': first normal equations and then augmented

      * ``u0``:       approximate left singular vectors (see |SinitSize|) {[]}
      * ``v0``:       approximate right singular vectors {[]}
      * ``orthoConst``: external orthogonalization constraints (see |SnumOrthoConst|) {[]}
      * |Slocking|:  1, hard locking; 0, soft locking 
      * |SmaxBlockSize|: maximum block size
      * |Siseed|:    random seed
      * |Sprimme|:   options for first stage solver
      * |SprimmeStage2|: options for second stage solver

   The available options for ``OPTIONS.primme`` and ``primmeStage2`` are
   the same as :mat:func:`primme_eigs`, plus the option ``'method'``.

   ``S = primme_svds(A,K,SIGMA,OPTIONS,P)``
   ``S = primme_svds(A,K,SIGMA,OPTIONS,P1,P2)`` makes use of a preconditioner,
   applying ``P\X`` or ``(P1*P2)\X``. If ``P`` is ``[]`` then a preconditioner is not
   applied. ``P`` may be a function handle ``PFUN`` such that ``PFUN(X,'AHA')``
   returns an approximation of ``(A'*A)\X``, ``PFUN(X,'AAH')``, of ``(A*A')\X`` and
   ``PFUN(X,'aug')``, of ``[zeros(N,N) A';A zeros(M,M)]\X``.

   ``[U,S,V] = primme_svds(...)`` returns also the corresponding singular vectors.
   If ``A`` is M-by-N and ``K`` singular triplets are computed, then ``U`` is M-by-K
   with orthonormal columns, ``S`` is K-by-K diagonal, and ``V`` is N-by-K with
   orthonormal columns.

   ``[S,R] = primme_svds(...)``
   ``[U,S,V,R] = primme_svds(...)`` returns upper bounds of the residual norm
   of each ``K`` triplet, ``NORM([A*V(:,i)-S(i,i)*U(:,i); A'*U(:,i)-S(i,i)*V(:,i)])``.

   ``[U,S,V,R,STATS] = primme_svds(...)`` returns how many times ``A`` and ``P`` were
   used and elapsed time. The application of ``A`` is counted independently from
   the application of ``A'``.

   ``[U,S,V,R,STATS,HIST] = primme_svds(...)`` instead of printing the convergence
   history, it is returned. Every row is a record, and the columns report:
  
      * ``HIST(:,1)``: number of matvecs
      * ``HIST(:,2)``: time
      * ``HIST(:,3)``: number of converged/locked triplets
      * ``HIST(:,4)``: stage
      * ``HIST(:,5)``: block index
      * ``HIST(:,6)``: approximate singular value
      * ``HIST(:,7)``: residual norm
      * ``HIST(:,8)``: QMR residual norm

   ``OPTS.disp`` controls the granularity of the record. If ``OPTS.disp == 1``, ``HIST``
   has one row per converged triplet and only the first four columns are
   reported; if ``OPTS.disp == 2``, ``HIST`` has one row per outer iteration and only
   the first seven columns are reported; and otherwise ``HIST`` has one row per QMR
   iteration and all columns are reported.

   Examples::

      A = diag(1:50); A(200,1) = 0; % rectangular matrix of size 200x50

      s = primme_svds(A,10) % the 10 largest singular values

      s = primme_svds(A,10,'S') % the 10 smallest singular values

      s = primme_svds(A,10,25) % the 10 closest singular values to 25

      opts = struct();
      opts.tol = 1e-4; % set tolerance
      opts.method = 'primme_svds_normalequations' % set svd solver method
      opts.primme.method = 'DEFAULT_MIN_TIME' % set first stage eigensolver method
      opts.primme.maxBlockSize = 2; % set block size for first stage
      [u,s,v] = primme_svds(A,10,'S',opts); % find 10 smallest svd triplets

      opts.orthoConst = {u,v};  
      [s,rnorms] = primme_svds(A,10,'S',opts) % find another 10

      % Compute the 5 smallest singular values of a rectangular matrix using
      % Jacobi preconditioner on (A'*A)
      A = sparse(diag(1:50) + diag(ones(49,1), 1));
      A(200,50) = 1;  % size(A)=[200 50]
      Pstruct = struct('AHA', diag(A'*A),...
                       'AAH', ones(200,1), 'aug', ones(250,1));
      Pfun = @(x,mode)Pstruct.(mode).\x;
      s = primme_svds(A,5,'S',[],Pfun) % find the 5 smallest values

   See also: `MATLAB svds`_, :mat:func:`primme_eigs`

.. include:: epilog.inc
