function [varargout] = primme_eigs(varargin)
%PRIMME_EIGS  Find a few eigenvalues/vectors of large, sparse Hermitian matrices
%
%   D = PRIMME_EIGS(A) returns a vector of A's 6 largest magnitude eigenvalues.
%
%   D = PRIMME_EIGS(A,B) returns a vector of the 6 largest magnitude eigenvalues
%   of the generalized eigenproblem (A,B).
%
%   D = PRIMME_EIGS(AFUN,DIM)
%   D = PRIMME_EIGS(AFUN,BFUN,DIM) accepts the functions AFUN and BFUN instead
%   of matrices. AFUN and BFUN are function handles. AFUN(x) and BFUN(x) return
%   the matrix-vector product A*x and B*x.
%
%   D = PRIMME_EIGS(...,K) finds the K largest magnitude eigenvalues. K must be
%   less than the dimension of the matrix A.
%
%   D = PRIMME_EIGS(...,K,TARGET) returns K eigenvalues such that: 
%     If TARGET is a real number, it finds the closest eigenvalues to TARGET.
%     If TARGET is
%       'LA' or 'SA', eigenvalues with the largest or smallest algebraic value
%       'LM' or 'SM', eigenvalues with the largest or smallest magnitude if
%                 OPTS.targetShifts is empty. If TARGET is a real or complex 
%                 scalar including 0, PRIMME_EIGS finds the eigenvalues closest 
%                 to TARGET.
%                 In addition, if m values are provided in OPTS.targetShifts, 
%                 find eigenvalues that are farthest (LM) or closest (SM) in 
%                 absolute value from the given values. 
%                 Examples: 
%                 k=1, 'LM', OPTS.targetShifts=[] returns the largest magnitude lambda(A).
%                 k=1, 'SM', OPTS.targetShifts=[] returns the smallest magnitude lambda(A).
%                 k=3, 'SM', OPTS.targetShifts=[2, 5] returns the closest eigenvalue in 
%                 absolute sense to 2, and the two closest eigenvalues to 5.
%       'CLT' or 'CGT', find eigenvalues closest to but less or greater than
%                 the given values in OPTS.targetShifts.
%
%   D = PRIMME_EIGS(...,K,TARGET,OPTS) specifies extra solver parameters. Some
%     default values are indicated in brackets {}:
%
%     OPTS.aNorm: the estimated 2-norm of A {0.0 (estimate the norm internally)}
%     OPTS.tol: convergence tolerance:                      {eps*1e4}
%                NORM(A*X(:,i)-X(:,i)*D(i,i)) < tol*NORM(A)
%     OPTS.maxBlockSize: maximum block size (useful for high multiplicities) {1}
%     OPTS.disp: different level reporting (0-3) (see HIST) {no output 0}
%     OPTS.display: toggle information display (see HIST)
%     OPTS.isreal: whether A represented by AFUN is real or complex {false}
%     OPTS.isdouble: whether the class of in/out vectors in AFUN are
%          double or single {false}
%     OPTS.isgpu: whether the class of in/out vectors in AFUN are gpuArray {false}
%     OPTS.ishermitian: whether A is Hermitian; otherwise it is considered normal {true}
%     OPTS.targetShifts: shifts for interior eigenvalues (see TARGET) {[]}
%     OPTS.v0: any number of initial guesses to the eigenvectors {[]}
%     OPTS.orthoConst: external orthogonalization constraints {[]}
%     OPTS.locking: 1, hard locking; 0, soft locking
%     OPTS.p: maximum size of the search subspace
%     OPTS.minRestartSize: minimum Ritz vectors to keep in restarting
%     OPTS.maxMatvecs: maximum number of matrix vector multiplications {Inf}
%     OPTS.maxit: maximum number of outer iterations {Inf}
%     OPTS.scheme: the restart scheme {'primme_thick'}
%     OPTS.maxPrevRetain: number of Ritz vectors from previous iteration
%          that are kept after restart {typically >0, see PRIMME doc}
%     OPTS.robustShifts: setting to true may avoid stagnation or misconvergence 
%     OPTS.maxInnerIterations: maximum number of inner solver iterations
%     OPTS.LeftQ: use the locked vectors in the left projector
%     OPTS.LeftX: use the approx. eigenvector in the left projector
%     OPTS.RightQ: use the locked vectors in the right projector
%     OPTS.RightX: use the approx. eigenvector in the right projector
%     OPTS.SkewQ: use the preconditioned locked vectors in the right projector
%     OPTS.SkewX: use the preconditioned approx. eigenvector in the right projector
%     OPTS.relTolBase: a legacy from classical JDQR (not recommended)
%     OPTS.convTest: how to stop the inner QMR Method
%     OPTS.convTestFun: function handler with an alternative convergence criterion.
%          If FUN(EVAL,EVEC,RNORM) returns a nonzero value, the pair (EVAL,EVEC)
%          with residual norm RNORM is considered converged.
%     OPTS.iseed: random seed
%     OPTS.profiler: return times from selected PRIMME's internal functions.
%          If 1, STATS returns times for the main functions. If it is a cell,
%          STATS returns times for those functions only. For instance,
%          {'Bortho'} returns times for all function calls containing 'Bortho'
%          on main_iteration. {'??/Bortho'} returns all calls regardless the
%          caller. {'**/Borth'} returns all calls containing 'Bortho' grouped by
%          callers. And {'++/Borth'} grouped by invocations. {'Bortho/*'}
%          returns times taken by functions called at Bortho.
%
%   For detailed descriptions of the above options, visit:
%   http://www.cs.wm.edu/~andreas/software/doc/primmec.html#parameters-guide
%
%   D = PRIMME_EIGS(...,K,TARGET,OPTS,METHOD) specifies the eigensolver method:
%     'DYNAMIC', (default)        switches dynamically to the best method
%     'DEFAULT_MIN_TIME',         best method for low-cost matrix-vector product
%     'DEFAULT_MIN_MATVECS',      best method for heavy matvec/preconditioner
%     'Arnoldi',                  Arnoldi not implemented efficiently
%     'GD',                       classical block Generalized Davidson 
%     'GD_plusK',                 GD+k block GD with recurrence restarting
%     'GD_Olsen_plusK',           GD+k with approximate Olsen precond.
%     'JD_Olsen_plusK',           GD+k, exact Olsen (two precond per step)
%     'RQI',                      Rayleigh Quotient Iteration. Also INVIT,
%                                 but for INVIT provide OPTS.targetShifts
%     'JDQR',                     Original block, Jacobi Davidson
%     'JDQMR',                    Our block JDQMR method (similar to JDCG)
%     'JDQMR_ETol',               Slight, but efficient JDQMR modification
%     'STEEPEST_DESCENT',         equiv. to GD(block,2*block)
%     'LOBPCG_OrthoBasis',        equiv. to GD(nev,3*nev)+nev
%     'LOBPCG_OrthoBasis_Window'  equiv. to GD(block,3*block)+block nev>block
%
%   For further description of the method visit:
%   http://www.cs.wm.edu/~andreas/software/doc/appendix.html#preset-methods
%
%   D = PRIMME_EIGS(...,K,TARGET,OPTS,METHOD,P) 
%   D = PRIMME_EIGS(...,K,TARGET,OPTS,METHOD,P1,P2) uses preconditioner P or
%   P = P1*P2 to accelerate convergence of the method. Applying P\x should
%   approximate (A-sigma*eye(N))\x, for sigma near the wanted eigenvalue(s).
%   If P is [] then a preconditioner is not applied. P may be a function 
%   handle PFUN such that PFUN(x) returns P\x.
%
%   [X,D] = PRIMME_EIGS(...) returns a diagonal matrix D with the eigenvalues
%   and a matrix X whose columns are the corresponding eigenvectors. 
% 
%   [X,D,R] = PRIMME_EIGS(...) also returns an array of the residual norms of
%   the computed eigenpairs.
%
%   [X,D,R,STATS] = PRIMME_EIGS(...) returns a struct to report statistical
%   information about number of matvecs, elapsed time, and estimates for the
%   largest and smallest algebraic eigenvalues of A, and functions selected
%   OPTS.profile.
%
%   [X,D,R,STATS,HIST] = PRIMME_EIGS(...) it returns the convergence history,
%   instead of printing it. Every row is a record, and the columns report:
%  
%   HIST(:,1): number of matvecs
%   HIST(:,2): time
%   HIST(:,3): number of converged/locked pairs
%   HIST(:,4): block index
%   HIST(:,5): approximate eigenvalue
%   HIST(:,6): residual norm
%   HIST(:,7): QMR residual norm
%
%   OPTS.disp controls the granularity of the record. If OPTS.disp == 1, HIST
%   has one row per converged eigenpair and only the first three columns
%   together with the fifth and the sixth are reported. If OPTS.disp == 2, HIST
%   has one row per outer iteration and converged value, and only the first six
%   columns are reported. Otherwise HIST has one row per QMR iteration, outer
%   iteration and converged value, and all columns are reported.
%
%   The convergence history is displayed if OPTS.disp > 0 and either HIST is
%   not returned or OPTS.display == 1.
%  
%   Examples:
%      A = diag(1:100);
%
%      d = primme_eigs(A,10) % the 10 largest magnitude eigenvalues
%
%      d = primme_eigs(A,10,'SM') % the 10 smallest magnitude eigenvalues
%
%      d = primme_eigs(A,10,25.0) % the 10 closest eigenvalues to 25.0
%
%      opts.targetShifts = [2 20];
%      d = primme_eigs(A,10,'SM',opts) % 1 eigenvalue closest to 2 and 
%                                      % 9 eigenvalues closest to 20
%      B = diag(100:-1:1);
%      d = primme_eigs(A,B,10,'SM') % the 10 smallest magnitude eigenvalues
%
%      opts = struct();
%      opts.tol = 1e-4; % set tolerance
%      opts.maxBlockSize = 2; % set block size
%      [x,d] = primme_eigs(A,10,'SA',opts,'DEFAULT_MIN_TIME')
%
%      opts.orthoConst = x;  
%      [d,rnorms] = primme_eigs(A,10,'SA',opts) % find another 10
%
%      % Compute the 6 eigenvalues closest to 30.5 using ILU(0) as a precond.
%      % by passing the matrices L and U.
%      A = sparse(diag(1:50) + diag(ones(49,1), 1) + diag(ones(49,1), -1));
%      [L,U] = ilu(A, struct('type', 'nofill'));
%      d = primme_eigs(A, k, 30.5, [], [], L, U);
%
%      % Compute the 6 eigenvalues closest to 30.5 using Jacobi preconditioner
%      % by passing a function.
%      Pfun = @(x)(diag(A) - 30.5)\x;
%      d = primme_eigs(A,6,30.5,[],[],Pfun);
%
%   For more details see PRIMME documentation at
%   http://www.cs.wm.edu/~andreas/software/doc/readme.html 
%
%   See also PRIMME_SVDS, EIGS.

   % Check primme_mex exists
   if ~ exist('primme_mex')
      warning 'primme_mex is not available. Building PRIMME...'
      make
   end

   % Check arity of input and output arguments
   minInputs = 1;
   maxInputs = 9;
   narginchk(minInputs,maxInputs);

   minOutputs = 0;
   maxOutputs = 5;
   nargoutchk(minOutputs,maxOutputs);

   % Check input arguments
   opts = struct();
   A = varargin{1};
   nextArg = 2;
   isgeneralized = 0;
   Acomplex = true;
   Adouble = true;
   Agpu = false;
   Aherm = true;
   if isnumeric(A)
      % Check matrix is Hermitian and get matrix dimension
      [m, n] = size(A);
      opts.n = n;
      opts.matrixMatvec = @(x)A*x;

      % Get type and complexity
      Acomplex = ~isreal(A);
      Agpu = strcmp(class(A), 'gpuArray');
      if Agpu
         Adouble = strcmp(classUnderlying(A), 'double');
      else
         Adouble = strcmp(class(A), 'double');
      end
      ABfun = 0;
   else
      opts.matrixMatvec = fcnchk_gen(A); % get the function handle of user's function
      ABfun = 1;
   end

   if nargin >= nextArg && (~isnumeric(varargin{nextArg}) || ~isscalar(varargin{nextArg}))
      B = varargin{nextArg};
      if isnumeric(B)
         % Check matrix is Hermitian and get matrix dimension
         [m, n] = size(B);
         if m ~= n || m < 1e4 && ~ishermitian(B)
            error('Input matrix must be real symmetric or complex Hermitian');
         elseif ~ABfun && m ~= opts.n
            error('Input matrices A and B should have the same dimensions');
         end
         opts.massMatrixMatvec = @(x)B*x;

         % Get type and complexity
         Acomplex = Acomplex || ~isreal(B);
         Agpu = Agpu || strcmp(class(B), 'gpuArray');
         if strcmp(class(B), 'gpuArray')
            Adouble = Adouble || strcmp(classUnderlying(B), 'double');
         else
            Adouble = Adouble || strcmp(class(B), 'double');
         end
         isgeneralized = 1;
      elseif ~isempty(B)
         opts.massMatrixMatvec = fcnchk_gen(B); % get the function handle of user's function
         ABfun = 1;
         isgeneralized = 1;
      end
      nextArg = nextArg + 1;
   end

   if ABfun
      n = varargin{nextArg};
      if ~isscalar(n) || ~isnumeric(n) || (n<0) || ~isfinite(n)
         error(message('The size of input matrices must be a positive integer'));
      end
      n = round(n);
      opts.n = n;
      nextArg = nextArg + 1;
   end

   if nargin >= nextArg
      opts.numEvals = varargin{nextArg};
      if ~isscalar(opts.numEvals) || ~isnumeric(opts.numEvals) || (opts.numEvals<0) || ~isfinite(opts.numEvals)
         error(message('The argument numEvals must be a positive integer'));
      end
      opts.numEvals = round(opts.numEvals);
      nextArg = nextArg + 1;
   else
      opts.numEvals = min(6, opts.n);
   end

   if nargin >= nextArg
      target = varargin{nextArg};
      if isnumeric(target)
         opts.target = 'primme_closest_abs';
         opts.targetShifts = target;
      elseif ischar(target)
         targets = struct('LA', 'primme_largest', ...
                          'LM', 'primme_largest_abs', ...
                          'SA', 'primme_smallest', ...
                          'CGT', 'primme_closest_geq', ...
                          'CLT', 'primme_closest_leq', ...
                          'SM', 'primme_closest_abs');
         if ~isfield(targets, target)
            error('target must be LA, SA, LM, SM, CGT or CLT');
         end
         opts.target = getfield(targets, target);
         if (strcmp(target, 'SM') || strcmp(target, 'LM')) && ~isfield(opts, 'targetShifts')
            opts.targetShifts = 0;
         end
      else
         error('target must be a number or a string');
      end
      nextArg = nextArg + 1;
   else
      opts.target = 'primme_largest_abs';
      opts.targetShifts = 0;
   end

   if nargin >= nextArg
      if ~isempty(varargin{nextArg})
         opts0 = varargin{nextArg};
         if ~isstruct(opts0)
            error('opts must be a struct');
         end
         opts0_names = fieldnames(opts0);
         for i=1:numel(opts0_names)
            opts.(opts0_names{i}) = opts0.(opts0_names{i});
         end
      end
      nextArg = nextArg + 1;
   end

   method = 'PRIMME_DEFAULT_METHOD';
   if nargin >= nextArg
      if ~isempty(varargin{nextArg})
         method = varargin{nextArg};
         if ischar(method)
            method = ['PRIMME_' method];
         end
      end
      nextArg = nextArg + 1;
   end

   if nargin >= nextArg
      P = varargin{nextArg};
      if isnumeric(P)
         P = @(x)P\x;
      else
         P = fcnchk_gen(P); % get the function handle of user's function
      end
      nextArg = nextArg + 1;
   else
      P = [];
   end

   if nargin >= nextArg
      P2 = varargin{nextArg};
      if isnumeric(P2)
         P2 = @(x)P2\x;
      else
         P2 = fcnchk_gen(P2); % get the function handle of user's function
      end
      P = @(x)P2(P(x));
   end
   if ~isempty(P)
      opts.applyPreconditioner = P;
      opts.correction.precondition = 1;
   end
 
   % Process 'isreal' in opts
   if isfield(opts, 'isreal')
      Acomplex = ~opts.isreal;
      opts = rmfield(opts, 'isreal');
   end

   % Process 'isdouble' in opts
   if isfield(opts, 'isdouble')
      Adouble = opts.isdouble;
      opts = rmfield(opts, 'isdouble');
   end

   % Process 'ishermitian' in opts
   if isfield(opts, 'ishermitian')
      Aherm = opts.ishermitian;
      opts = rmfield(opts, 'ishermitian');
   end
   if isnumeric(A) && m < 1e4 && Aherm && ~ishermitian(A)
      error('Input matrix must be real symmetric or complex Hermitian, or set OPTS.ishermitian to false');
   end

   % Process 'isgpu' in opts
   if isfield(opts, 'isgpu')
      Agpu = opts.isgpu;
      opts = rmfield(opts, 'isgpu');
   end
   if Adouble
      Aclass = 'double';
   else
      Aclass = 'single';
   end
   if Agpu
      d = gpuDevice;
      opts.commInfo = d.Index - 1;
   end
   if isnumeric(A) && issparse(A) && strcmp(Aclass, 'single')
      opts.matrixMatvec_type = 'primme_op_double';
      Aclass = 'double';
   end

   % Test whether the given matrix and preconditioner are valid
   try
      if ~Agpu
         test_x = ones(opts.n, 1, Aclass);
      else
         test_x = ones(opts.n, 1, Aclass, 'gpuArray');
      end
      x = opts.matrixMatvec(test_x);
      if isfield(opts, 'applyPreconditioner')
         x = opts.applyPreconditioner(test_x);
      end
      clear test_x;
      clear x;
   catch ME
      rethrow(ME);
   end

   % Process 'display' in opts
   showHist = [];
   dispLevel = 0;
   if isfield(opts, 'display')
      showHist = opts.display;
      if numel(showHist) ~= 1 || (showHist ~= 0 && showHist ~= 1)
         error('Invalid value in opts.display; it should be 0 or 1');
      end
      opts = rmfield(opts, 'display');
      if showHist
         dispLevel = 1;
      end
   elseif nargout >= 5
      showHist = false;
      dispLevel = 1;
   end

   % Process 'disp' in opts
   if isfield(opts, 'disp')
      dispLevel = opts.disp;
      if dispLevel > 3 || dispLevel < 0
         error('Invalid value in opts.disp; it should be 0, 1, 2 or 3');
      end
      opts = rmfield(opts, 'disp');
   elseif nargout >= 5 || (~isempty(showHist) && showHist)
      dispLevel = 1;
   end
   if isempty(showHist)
      showHist = dispLevel > 0;
   end

   % Process profile
   profile0 = {};
   if isfield(opts, 'profile')
      if isnumeric(opts.profile) && numel(opts.profile) == 1 && opts.profile == 1
         opts.profile = {'init','update_Q','update_projection','solve_H','check_convergence','prepare_candidates','Bortho','restart'};
      end
      profile0 = opts.profile;
      opts.profile = get_regex_from_cell(opts.profile);
   end

   % Rename tol, maxit and p as eps, maxOuterIterations and maxBasisSize.
   %  Also move options that are outside of primme_params' hierarchy.
   changes = {{'tol', 'eps'}, {'maxit', 'maxOuterIterations'}, {'p', 'maxBasisSize'}, ...
              {'projection',         'projection_projection'}, ...
              {'scheme',             'restarting_scheme'}, ...
              {'maxPrevRetain',      'restarting_maxPrevRetain'}, ...
              {'precondition',       'correction_precondition'}, ...
              {'robustShifts',       'correction_robustShifts'}, ...
              {'maxInnerIterations', 'correction_maxInnerIterations'}, ...
              {'LeftQ',              'correction_projectors_LeftQ'}, ...
              {'LeftX',              'correction_projectors_LeftX'}, ...
              {'RightQ',             'correction_projectors_RightQ'}, ...
              {'RightX',             'correction_projectors_RightX'}, ...
              {'SkewQ',              'correction_projectors_SkewQ'}, ...
              {'SkewX',              'correction_projectors_SkewX'}, ...
              {'convTest',           'correction_convTest'}, ...
              {'relTolBase',         'correction_relTolBase'}};

   for i=1:numel(changes)
      if isfield(opts, changes{i}{1})
         opts.(changes{i}{2}) = opts.(changes{i}{1});
         opts = rmfield(opts, changes{i}{1});
      end
   end

   % Prepare numOrthoConst and initSize
   if isfield(opts, 'orthoConst')
      init = opts.orthoConst;
      if size(init, 1) ~= opts.n
         error('Invalid matrix dimensions in opts.orthoConst');
      end
      opts = rmfield(opts, 'orthoConst');
      opts.numOrthoConst = size(init, 2);
   else
      init = [];
   end

   if isfield(opts, 'v0')
      init0 = opts.v0;
      if size(init0, 1) ~= opts.n
         error('Invalid matrix dimensions in opts.v0');
      end
      opts = rmfield(opts, 'v0');
      opts.initSize = size(init0, 2);
      init = [init init0];
   end

   % Set default tol
   if ~isfield(opts, 'eps')
      opts.eps = eps(Aclass)*1e4;
   end 

   % Create primme_params
   primme = primme_mex('primme_initialize');

   % This long try-catch make sure that primme_free is called
   try
      % Set other options in primme_params
      primme_set_members(opts, primme);

      % Set method
      primme_mex('primme_set_method', method, primme);

      % Set monitor and shared variables with the monitor
      hist = [];
      histSize = 0;
      prof = struct();
      locking = primme_mex('primme_get_member', primme, 'locking');
      nconv = [];
      return_hist = nargout >= 5;
      return_prof = nargout >= 4;

      if dispLevel > 0 || return_prof
         % NOTE: Octave doesn't support function handler for nested functions
         primme_mex('primme_set_member', primme, 'monitorFun', ...
               @(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)record_history(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12));
      end
      if showHist
         if dispLevel == 1
            fprintf('#  MV\t Time\t NConv\t Value\t  Res\n');
         elseif dispLevel == 2
            fprintf('#  MV\t Time\t NConv\t  Idx\t Value\t  Res\n');
         elseif dispLevel == 3
            fprintf('#  MV\t Time\t NConv\t  Idx\t Value\t  Res\t  QMR_Res\n');
         end
      end

      % Select solver
      if Adouble
         if Acomplex
            type = 'z';
         else
            type = 'd';
         end
      else
         if Acomplex
            type = 'c';
         else
            type = 's';
         end
      end
      if Agpu
         type = ['magma_' type];
      end
      xprimme = [type 'primme'];
      if ~Aherm
         xprimme = [xprimme '_normal'];
      end

      % Call xprimme
      [ierr, evals, norms, evecs] = primme_mex(xprimme, init, primme); 

      % Process error code and return the required arguments
      if ierr == -3
         warning([xprimme ' returned ' num2str(ierr) ': ' primme_error_msg(ierr)]);
      elseif ierr ~= 0
         error([xprimme ' returned ' num2str(ierr) ': ' primme_error_msg(ierr)]);
      end
      
      % Return interior eigenvalues in descending order
      if ~strcmp(opts.target,'primme_largest') ...
            && ~strcmp(opts.target,'primme_smallest') ...
            && ~strcmp(opts.target,'primme_largest_abs')
         [evals,ind] = sort(evals,'descend');
         evecs = evecs(:,ind);
      end

      if (nargout <= 1)
         varargout{1} = evals;
      end
      if (nargout >= 2)
         varargout{1} = evecs;
         varargout{2} = diag(evals);
      end
      if (nargout >= 3)
         varargout{3} = norms;
      end
      if (nargout >= 4)
         if return_prof
            stats = make_nice_profile(prof, profile0);
         else
            stats = struct();
         end
         stats.numMatvecs = primme_mex('primme_get_member', primme, 'stats_numMatvecs');
         stats.timeMatvecs = primme_mex('primme_get_member', primme, 'stats_timeMatvec');
         stats.numPreconds = primme_mex('primme_get_member', primme, 'stats_numPreconds');
         stats.timeOrtho = primme_mex('primme_get_member', primme, 'stats_timeOrtho');
         stats.numOrthoInnerProds = primme_mex('primme_get_member', primme, 'stats_numOrthoInnerProds');
         stats.elapsedTime = primme_mex('primme_get_member', primme, 'stats_elapsedTime');
         stats.estimateMinEVal = primme_mex('primme_get_member', primme, 'stats_estimateMinEVal');
         stats.estimateMaxEVal = primme_mex('primme_get_member', primme, 'stats_estimateMaxEVal');
         stats.estimateLargestSVal = primme_mex('primme_get_member', primme, 'stats_estimateLargestSVal');
         varargout{4} = stats;
      end
      if (nargout >= 5)
         varargout{5} = hist(1:histSize,:);
      end
   catch ME
      primme_mex('primme_free', primme);
      rethrow(ME);
   end
   primme_mex('primme_free', primme);

   function record_history(basisEvals, basisFlags, iblock, basisNorms, ...
         numConverged, lockedEvals, lockedFlags, lockedNorms, inner_its, ...
         LSRes, msg, time, event)

      if event == 6 % primme_event_message
         warning(['PRIMME: ' msg]);
         return;
      elseif event == 7 % primme_event_profiler
         if return_prof
            if ~isfield(prof, msg)
               prof.(msg) = [];
            end
            prof.(msg) = [prof.(msg) time];
         end
         return;
      end

      numMatvecs = double(primme_mex('primme_get_member', primme, 'stats_numMatvecs'));
      maxInnerIterations = primme_mex('primme_get_member', primme, 'correction_maxInnerIterations');
      elapsedTime = primme_mex('primme_get_member', primme, 'stats_elapsedTime');
      histline = [];
      if event == 0 || (event == 4 && ~locking) || event == 5
         if ~locking
            nconv = double(numConverged);
         else
            nconv = numel(lockedEvals);
         end
      end
      if dispLevel == 0
         % Do nothing
      elseif dispLevel == 1
         if event == 4 && ~locking
            for i=1:numel(iblock)
               histline = [numMatvecs elapsedTime nconv basisEvals(iblock(i)+1) basisNorms(iblock(i)+1)];
            end
         elseif event == 5
            histline = [histline; numMatvecs elapsedTime nconv lockedEvals(end) lockedNorms(end)];
         end
      elseif dispLevel == 2
         if (event == 4 && ~locking) || event == 0
            for i=1:numel(iblock)
               histline = [histline; numMatvecs elapsedTime nconv i basisEvals(iblock(i)+1) basisNorms(iblock(i)+1)];
            end
         elseif event == 5
               histline = [histline; numMatvecs elapsedTime nconv 1 lockedEvals(end) lockedNorms(end)];
         end
      elseif dispLevel == 3
         if event == 1
            if ~isempty(basisEvals)
               value = basisEvals(iblock(1)+1);
               resNorm = basisNorms(iblock(1)+1);
            else
               value = nan;
               resNorm = nan;
            end
            histline = [histline; numMatvecs elapsedTime nconv nan value resNorm  LSRes];
         elseif (maxInnerIterations == 0 || nconv == opts.numEvals) && (event == 0 || (event == 4 && ~locking))
            for i=1:numel(iblock)
               histline = [histline; numMatvecs elapsedTime nconv i basisEvals(iblock(i)+1) basisNorms(iblock(i)+1) nan];
            end
         elseif (maxInnerIterations == 0 || nconv == opts.numEvals) && event == 5
               histline = [histline; numMatvecs elapsedTime nconv 1 lockedEvals(end) lockedNorms(end) nan];
         end
      end
      if showHist && size(histline,1) > 0
         template{1} = '%7d\t%-5.g\t%7d\t%s\t%-5.1e\n';
         template{2} = '%7d\t%-5.g\t%7d\t%7d\t%s\t%5.1e\n';
         template{3} = '%7d\t%-5.g\t%7d\t%7d\t%s\t%5.1e\t%5.1e\n';
         for i=1:size(histline,1)
            a = num2cell(histline(i,:));
            if dispLevel == 1, ieval = 4; else ieval = 5; end
            a{ieval} = num2str(a{ieval}, '%-5.1e');
            fprintf(template{dispLevel}, a{:});
         end
      end
      if return_hist
         if size(hist,1) < histSize + size(histline,1)
            l = max(histSize*2, histSize + size(histline,1));
            hist(l,size(histline,2)) = 0;
         end
         hist(histSize+1:histSize+size(histline,1),:) = histline;
         histSize = histSize + size(histline,1);
      end
   end
end

function [f] = fcnchk_gen(x)
   if exist('fcnchk', 'var')
      f = fcnchk(x);
   else 
      f = x;
   end
end

function primme_set_members(opts, primme, prefix)
%PRIMME_SET_MEMBERS  Set options in primme_params
%   PRIMME_SET_MEMBERS(S, P) sets the options in struct S into the primme_params
%   reference P.
%
%   Example:
%     primme = primme_mex('primme_initialize');
%     ops.n = 10;
%     ops.target = 'primme_largest';
%     primme_set_members(ops, primme);

   % NOTE: Expensive Mathworks' MATLAB doesn't support default values in function
   %       declaration, Octave does.
   if nargin < 3, prefix = ''; end
   
   fields = fieldnames(opts);
   for i=1:numel(fields)
      value = getfield(opts, fields{i});
      label = fields{i};
      if isstruct(value)
         primme_set_members(value, primme, [prefix label '_']);
      else
        try
      	  primme_mex('primme_set_member', primme, [prefix label], value);
        catch ME
          if isnumeric(value)
            error(['Error setting the option ' prefix label ' to value ' num2str(value)]);
          else
            error(['Error setting the option ' prefix label ' to value ' value]);
          end
        end
      end
   end
end

function s = primme_error_msg(errorCode)

   msg = {};
   msg{45+  0} = 'success';
   msg{45+  1} = 'reported only amount of required memory';
   msg{45+ -1} = 'unexpected failure';
   msg{45+ -2} = 'memory allocation failure';
   msg{45+ -3} = 'iteration error; usually maximum iterations or matvecs reached';
   msg{45+ -4} = 'argument primme is NULL';
   msg{45+ -5} = 'n < 0 or nLocal < 0 or nLocal > n';
   msg{45+ -6} = 'numProcs' < 1';
   msg{45+ -7} = 'matrixMatvec is NULL';
   msg{45+ -8} = 'applyPreconditioner is NULL and precondition is not NULL';
   msg{45+ -9} = 'not used';
   msg{45+-10} = 'numEvals > n';
   msg{45+-11} = 'numEvals < 0';
   msg{45+-12} = 'eps > 0 and eps < machine precision';
   msg{45+-13} = 'target is not properly defined';
   msg{45+-14} = 'target is one of primme_largest_abs, primme_closest_geq, primme_closest_leq or primme_closest_abs but numTargetShifts <= 0 (no shifts)';
   msg{45+-15} = 'target is one of primme_largest_abs primme_closest_geq primme_closest_leq or primme_closest_abs but targetShifts is NULL  (no shifts array)';
   msg{45+-16} = 'numOrthoConst < 0 or numOrthoConst > n (no free dimensions left)';
   msg{45+-17} = 'maxBasisSize < 2';
   msg{45+-18} = 'minRestartSize < 0 or minRestartSize shouldn''t be zero';
   msg{45+-19} = 'maxBlockSize < 0 or maxBlockSize shouldn''t be zero';
   msg{45+-20} = 'maxPrevRetain < 0';
   msg{45+-21} = 'scheme is not one of *primme_thick* or *primme_dtr*';
   msg{45+-22} = 'initSize < 0';
   msg{45+-23} = 'locking == 0 and initSize > maxBasisSize';
   msg{45+-24} = 'locking and initSize > numEvals';
   msg{45+-25} = 'maxPrevRetain + minRestartSize >= maxBasisSize';
   msg{45+-26} = 'minRestartSize >= n';
   msg{45+-27} = 'printLevel < 0 or printLevel > 5';
   msg{45+-28} = 'convTest is not one of primme_full_LTolerance primme_decreasing_LTolerance primme_adaptive_ETolerance or primme_adaptive';
   msg{45+-29} = 'convTest == primme_decreasing_LTolerance and relTolBase <= 1';
   msg{45+-30} = 'evals is NULL, but not evecs and resNorms';
   msg{45+-31} = 'evecs is NULL, but not evals and resNorms';
   msg{45+-32} = 'resNorms is NULL, but not evecs and evals';
   msg{45+-33} = 'locking == 0 and minRestartSize < numEvals';
   msg{45+-34} = 'ldevecs is less than nLocal';
   msg{45+-35} = 'ldOPs is non-zero and less than nLocal';
   msg{45+-36} = 'not enough memory for realWork';
   msg{45+-37} = 'not enough memory for intWork';
   msg{45+-38} = 'locking == 0 and target is primme_closest_leq or primme_closet_geq';
   msg{45+-40} = 'factorization failure';
   msg{45+-41} = 'user cancelled execution';
   msg{45+-42} = 'orthogonalization failure';
   msg{45+-43} = 'parallel failure';
   msg{45+-44} = 'unavailable functionality';

   errorCode = errorCode + 45;
   if errorCode > 0 && errorCode <= numel(msg)
      s = msg{errorCode};
   else
      s = 'Unknown error code';
   end
end

function x = replace_globs(x)
   x = strrep(strrep(strrep(x, '??', '%'), '++', '%'), '**', '%');
   x = strrep(strrep(strrep(x, '*', '[^~]*'), '+', '[^~]*'), '?', '[^~]*');
   x = strrep(x, '%', '.*');
end

function r = get_regex(c)
   if ~ischar(c)
      error('Not valid regex');
      return;
   end
  
   r = ['^~[^~]*' strjoin(cellfun(@(x)sprintf('~%s[^~]*', replace_globs(x)), strsplit(c, '/'), 'UniformOutput',false), '') '$'];
end

function r = get_regex_from_cell(c)
   if ischar(c)
      r = c;
   elseif iscell(c)
      r = strjoin(cellfun(@(x)['\(' get_regex(x) '\)'],c,'UniformOutput',false),'\\|');
   else
      error('Not valid profile');
   end
end

function r = get_group_name(f, s)
   s = strsplit(['?/' s], '/');
   r = {};
   for si = 1:numel(s)
      p = sprintf('^~%s[^~]*', replace_globs(s{si}));
      [~,l] = regexp(f, p);
      d = f(2:l);
      if strfind(s{si}, '?')
         r{end+1} = s{si};
      elseif strfind(s{si}, '**')
         r{end+1} = strjoin(cellfun(@(x)regexp(x,'[^(]+','match','once'), strsplit(d,'~'), 'UniformOutput',false), '/');
      elseif strfind(s{si}, '+')
         r{end+1} = d;
      elseif strfind(s{si}, '*')
         r{end+1} = regexp(d,'[^(]+','match','once');
      else
         r{end+1} = s{si};
      end
      f = f(l+1:end);
   end
   r = strjoin(r(2:end), '/');
end

function r = make_nice_profile(prof, groups)
   if ~iscell(groups)
      r = prof;
      return;
   end

   r = struct();
   func_names = fieldnames(prof);
   for fi = 1:numel(func_names);
      f = func_names{fi};
      for s = groups
         s = s{1};
         if ~isempty(regexp(f, get_regex(s)))
            g = get_group_name(f, s);
            if ~isfield(r, g)
               r.(g) = [];
            end
            r.(g) = [r.(g) prof.(f)];
            continue
         end
      end
   end
end
