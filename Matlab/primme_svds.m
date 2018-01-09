function [varargout] = primme_svds(varargin)
%PRIMME_SVDS  Find a few singular values and vectors of large, sparse matrices
%
%   S = PRIMME_SVDS(A) returns a vector with the 6 largest singular values of A.
%
%   S = PRIMME_SVDS(AFUN,M,N) accepts the function handle AFUN to perform
%   the matrix vector products with an M-by-N matrix A. 
%   AFUN(X,'notransp') returns A*X while AFUN(X,'transp') returns A’*X.
%   In all the following, A can be replaced by AFUN,M,N.
% 
%   S = PRIMME_SVDS(A,K) computes the K largest singular values of A.
%
%   S = PRIMME_SVDS(A,K,SIGMA) computes the K singular values closest to the
%   scalar shift SIGMA. 
%   If SIGMA is a vector, find a singular value closest to each SIGMA(i)
%   If SIGMA is 'L', it computes the largest singular values.
%   if SIGMA is 'S', it computes the smallest singular values.
%
%   S = PRIMME_SVDS(A,K,SIGMA,OPTIONS) specifies extra solver parameters:
%   (for some parameters we refer to PRIMME_EIGS)
%
%   Field name       Parameter                               Default
%
%   OPTIONS.aNorm    estimation of the 2-norm A                  0.0
%   OPTIONS.tol      convergence tolerance (see eps):    1e-10 (1e-3 for single)
%                    NORM([A*V-U*S;A'*U-V*S]) <= tol * NORM(A).
%   OPTIONS.maxit    maximum number of matvecs  (see maxMatvecs) inf
%   OPTIONS.p        maximum basis size (see maxBasisSize)         -
%   OPTIONS.disp     level of reporting 0-3 (see HIST)             0
%   OPTIONS.isreal   if 0, the matrix is complex; else it's real   0
%   OPTIONS.isdouble if 0, the matrix is single; else it's double  1
%   OPTIONS.method   which equivalent eigenproblem to solve
%                    - 'primme_svds_normalequations': A'*A or A*A'
%                    - 'primme_svds_augmented': [0 A';A 0]
%                    - 'primme_svds_hybrid':               (default)
%                       first normal equations and then augmented
%   OPTIONS.u0       approx. left singular vectors                []
%   OPTIONS.v0       approx. right singular vectors               []
%   OPTIONS.orthoConst external orthogonalization constraints     [] 
%   OPTIONS.locking  1, hard locking; 0, soft locking              -
%   OPTIONS.maxBlockSize maximum block size                        1
%   OPTIONS.iseed    random seed
%   OPTIONS.primme   options for first stage solver                -
%   OPTIONS.primmeStage2 options for second stage solver           -
%   OPTIONS.convTestFun  alternative convergence criterion         -
%
%   If OPTIONS.convTestFun(SVAL,LSVEC,RSVEC,RNORM) returns a nonzero
%   value, the triplet (SVAL,LSVEC,RSVEC) with residual norm RNORM
%   is considered converged.
%
%   The available options for OPTIONS.primme and primmeStage2 are
%   the same as PRIMME_EIGS, plus the option 'method'. For detailed
%   descriptions of the above options, visit:
%   http://www.cs.wm.edu/~andreas/software/doc/svdsc.html#parameters-guide
%   and for further descriptions of the methods visit:
%   http://www.cs.wm.edu/~andreas/software/doc/appendixsvds.html#preset-methods
%
%   S = PRIMME_SVDS(A,K,SIGMA,OPTIONS,P) applies a preconditioner P as follows.
%   If P is a matrix it applies P\X and P'\X to approximate A\X and A'\X.
%   If P is a function handle, PFUN, PFUN(X,'notransp') returns P\X and
%   PFUN(X,'transp') returns P’\X, approximating A\X and A'\X respectively.
%   If P is a struct, it can have one or more of the following fields:
%     P.AHA\X or P.AHA(X) returns an approximation of (A'*A)\X, 
%     P.AAH\X or P.AAH(X) returns an approximation of (A*A')\X,
%     P.aug\X or P.aug(X) returns an approximation of [zeros(N,N) A';A zeros(M,M)]\X.
%   If P is [] then no preconditioner is applied.
%
%   S = PRIMME_SVDS(A,K,SIGMA,OPTIONS,P1,P2) applies a factorized preconditioner.
%   If both P1,P2 are nonempty, apply (P1*P2)\X to approximate A\X. 
%   If P1 is [] and P2 is nonempty, then (P2'*P2)\X approximates A'*A. 
%   P2 can be the R factor of an (incomplete) QR factorization of A or the L
%   factor of an (incomplete) LL' factorization of A'*A (RIF).
%   If both P1 and P2 are [] then no preconditioner is applied.
%
%   [U,S,V] = PRIMME_SVDS(...) returns also the corresponding singular vectors.
%   If A is M-by-N and K singular triplets are computed, then U is M-by-K
%   with orthonormal columns, S is K-by-K diagonal, and V is N-by-K with
%   orthonormal columns.
%
%   [S,R] = PRIMME_SVDS(...)
%   [U,S,V,R] = PRIMME_SVDS(...) returns the residual norm of each K triplet,
%   NORM([A*V(:,i)-S(i,i)*U(:,i); A'*U(:,i)-S(i,i)*V(:,i)]).
%
%   [U,S,V,R,STATS] = PRIMME_SVDS(...) returns how many times A and P were
%   used and elapsed time. The application of A is counted independently from
%   the application of A'.
%
%   [U,S,V,R,STATS,HIST] = PRIMME_SVDS(...) returns the convergence history,
%   instead of printing it. Every row is a record, and the columns report:
%  
%   HIST(:,1): number of matvecs
%   HIST(:,2): time
%   HIST(:,3): number of converged/locked triplets
%   HIST(:,4): stage
%   HIST(:,5): block index
%   HIST(:,6): approximate singular value
%   HIST(:,7): residual norm
%   HIST(:,8): QMR residual norm
%
%   OPTS.disp controls the granularity of the record. If OPTS.disp == 1, HIST
%   has one row per converged triplet and only the first four columns are
%   reported; if OPTS.disp == 2, HIST has one row per outer iteration and only
%   the first seven columns are reported; and otherwise HIST has one row per QMR
%   iteration and all columns are reported.
%
%   Examples:
%      A = diag(1:50); A(200,1) = 0; % rectangular matrix of size 200x50
%
%      s = primme_svds(A,10) % the 10 largest singular values
%
%      s = primme_svds(A,10,'S') % the 10 smallest singular values
%
%      s = primme_svds(A,10,25) % the 10 closest singular values to 25
%
%      opts = struct();
%      opts.tol = 1e-4; % set tolerance
%      opts.method = 'primme_svds_normalequations' % set svd solver method
%      opts.primme.method = 'DEFAULT_MIN_TIME' % set first stage eigensolver method
%      opts.primme.maxBlockSize = 2; % set block size for first stage
%      [u,s,v] = primme_svds(A,10,'S',opts); % find 10 smallest svd triplets
%
%      opts.orthoConst = {u,v};  
%      [s,rnorms] = primme_svds(A,10,'S',opts) % find another 10
%
%      % Compute the 5 smallest singular values of a square matrix using ILU(0)
%      % as a preconditioner
%      A = sparse(diag(1:50) + diag(ones(49,1), 1));
%      [L,U] = ilu(A, struct('type', 'nofill'));
%      svals = primme_svds(A, 5, 'S', [], L, U);
%      
%      % Compute the 5 smallest singular values of a rectangular matrix using
%      % Jacobi preconditioner on (A'*A)
%      A = sparse(diag(1:50) + diag(ones(49,1), 1));
%      A(200,50) = 1;  % size(A)=[200 50]
%      P = diag(sum(abs(A).^2));
%      precond.AHA = @(x)P\x;
%      s = primme_svds(A,5,'S',[],precond) % find the 5 smallest values
%
%      % Estimation of the smallest singular value
%      A = diag([1 repmat(2,1,1000) 3:100]);
%      [~,sval,~,rnorm] = primme_svds(A,1,'S',struct('convTestFun',@(s,u,v,r)r<s*.1));
%      sval - rnorm % approximate smallest singular value
%
%   For more details see PRIMME documentation at
%   http://www.cs.wm.edu/~andreas/software/doc/readme.html 
%
%   See also PRIMME_EIGS, SVDS.

   % Check primme_mex exists
   if ~ exist('primme_mex')
      error 'primme_mex is not available. Try to recompile the MATLAB/Octave''s PRIMME module'
   end

   % Check arity of input and output arguments
   minInputs = 1;
   maxInputs = 8;
   narginchk(minInputs,maxInputs);

   minOutputs = 0;
   maxOutputs = 6;
   nargoutchk(minOutputs,maxOutputs);

   % Check input arguments
   opts = struct();
   A = varargin{1};
   nextArg = 2;
   if isnumeric(A)
      % Check matrix is Hermitian and get matrix dimension
      [m, n] = size(A);
      opts.m = m;
      opts.n = n;
      opts.matrixMatvec = @(x,mode)matvecsvds(A,x,mode);

      % Get type and complexity
      Acomplex = ~isreal(A);
      Adouble = strcmp(class(A), 'double');
   else
      opts.matrixMatvec = fcnchk_gen(A); % get the function handle of user's function
      m = round(varargin{nextArg});
      n = round(varargin{nextArg+1});
      if ~isscalar(m) || ~isreal(m) || (m<0) || ~isfinite(m) || ...
         ~isscalar(n) || ~isreal(n) || (n<0) || ~isfinite(n)
         error(message('The size of input matrix A must be an positive integer'));
      end
      opts.m = m;
      opts.n = n;
      nextArg = nextArg + 2;

      % Assume complex double matrix
      Acomplex = 1;
      Adouble = 1;
   end

   if nargin >= nextArg
      opts.numSvals = round(varargin{nextArg});
      nextArg = nextArg + 1;
   else
      opts.numSvals = min([6 opts.m opts.n]);
   end

   if nargin >= nextArg
      target = varargin{nextArg};
      if ischar(target)
         targets = struct('L', 'primme_svds_largest', ...
                          'S', 'primme_svds_smallest');
         if ~isfield(targets, target(1))
            error('target must be L, S or real non-negative numbers');
         end
         opts.target = getfield(targets, target(1));
      elseif isnumeric(target) && all(target == 0)
         opts.target = 'primme_svds_smallest';
      elseif isnumeric(target) && all(target >= 0)
         opts.targetShifts = target;
         opts.target = 'primme_svds_closest_abs';
      else
         error('target must be L, S or real non-negative numbers');
      end
      nextArg = nextArg + 1;
   else
      opts.target = 'primme_svds_largest';
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

   if nargin == nextArg || (nargin > nextArg && isempty(varargin{nextArg+1}))
      P = varargin{nextArg};
      if isnumeric(P)
         if ~isempty(P)
            P = @(x,mode)precondsvds_Pmat(P,x,mode);
         end
      elseif isstruct(P)
         P = @(x,mode)precondsvds_Pstruct(P,x,mode);
      else
         P = fcnchk_gen(P); % get the function handle of user's function
         P = @(x,mode)precondsvds_Pfun(P,x,mode,opts.m);
      end
      if ~isempty(P)
         opts.applyPreconditioner = P;
         opts.precondition = 1;
      end
   elseif nargin >= nextArg
      P1 = varargin{nextArg};
      P2 = varargin{nextArg+1};
      if (~isempty(P1) && ~isnumeric(P1)) || ~isnumeric(P2)
         error('P1 and P2 must be matrices');
      end
      P = @(x,mode)precondsvds_P1P2(P1, P2, x, mode);
      opts.applyPreconditioner = P;
      opts.precondition = 1;
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
   if Adouble
      Aclass = 'double';
   else
      Aclass = 'single';
   end

   % Test whether the given matrix and preconditioner are valid
   try
      x = opts.matrixMatvec(ones(opts.n, 1, Aclass), 'notransp');
      x = opts.matrixMatvec(ones(opts.m, 1, Aclass), 'transp');
      if isfield(opts, 'applyPreconditioner')
         x = opts.applyPreconditioner(ones(opts.n, 1, Aclass), 'AHA');
         x = opts.applyPreconditioner(ones(opts.m, 1, Aclass), 'AAH');
         x = opts.applyPreconditioner(ones(opts.m+opts.n, 1, Aclass), 'aug');
      end
      clear x;
   catch ME
      rethrow(ME);
   end

   % Process 'disp' in opts
   if isfield(opts, 'disp')
      dispLevel = opts.disp;
      if dispLevel > 3 || dispLevel < 0
         error('Invalid value in opts.disp; it should be 0, 1, 2 or 3');
      end
      opts = rmfield(opts, 'disp');
   elseif nargout >= 6
      dispLevel = 1;
   else
      dispLevel = 0;
   end

   % Rename tol, maxit and p as eps, maxMatvecs and maxBasisSize
   changes = {{'tol', 'eps'}, {'maxit', 'maxMatvecs'}, {'p', 'maxBasisSize'}};
   for i=1:numel(changes)
      if isfield(opts, changes{i}{1})
         opts.(changes{i}{2}) = opts.(changes{i}{1});
         opts = rmfield(opts, changes{i}{1});
      end
   end

   % Set default tol
   if ~isfield(opts, 'eps')
      if Adouble
         opts.eps = 1e-10;
      else
         opts.eps = eps(Aclass)*1e4;
      end
   end 

   % Move options that are outside of primme_parms' hierarchy
   changes = {{'projection',         'projection_projection'}, ...
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
   primme_fields = {'primme', 'primmeStage2'};
   for j=1:numel(primme_fields)
      if isfield(opts, primme_fields{j})
         opts0 = opts.(primme_fields{j});
         for i=1:numel(changes)
            if isfield(opts0, changes{i}{1})
               opts0.(changes{i}{2}) = opts0.(changes{i}{1});
               opts0 = rmfield(opts0, changes{i}{1});
            end
         end
         opts.(primme_fields{j}) = opts0;
      end
   end

   % Process method, primme.method and primmeStage2.method
   if isfield(opts, 'method')
      method = opts.method;
      opts = rmfield(opts, 'method');
   else
      method = 'primme_svds_default';
   end
   if isfield(opts, 'primme') && isfield(opts.primme, 'method')
      primmeStage0method = opts.primme.method;
      opts.primme = rmfield(opts.primme, 'method');
      if ischar(primmeStage0method)
         primmeStage0method = ['PRIMME_' primmeStage0method];
      end
   else
      primmeStage0method = 'PRIMME_DEFAULT_METHOD';
   end
   if isfield(opts, 'primmeStage2') && isfield(opts.primmeStage2, 'method')
      primmeStage1method = opts.primmeStage2.method;
      opts.primmeStage2 = rmfield(opts.primmeStage2, 'method');
      if ischar(primmeStage1method)
         primmeStage1method = ['PRIMME_' primmeStage1method];
      end
   else
      primmeStage1method = 'PRIMME_DEFAULT_METHOD';
   end
     
   % Prepare numOrthoConst and initSize
   if isfield(opts, 'orthoConst')
      init = opts.orthoConst;
      if ~iscell(init) || numel(init) ~= 2 || (isempty(init{1}) && isempty(init{2}))
         error('opts.orthoConst should be {left_vectors, right_vectors}');
      end
      if isempty(init{1})
         init{1} = opts.matrixMatvec(init{2}, 'notransp');
      elseif isempty(init{2})
         init{2} = opts.matrixMatvec(init{1}, 'transp');
      end
      if size(init{1}, 1) ~= opts.m || size(init{2}, 1) ~= opts.n || ...
         size(init{1}, 2) ~= size(init{2}, 2)
         error('Invalid matrix dimensions in opts.orthoConst');
      end
      opts = rmfield(opts, 'orthoConst');
      opts.numOrthoConst = size(init{1}, 2);
   else
      init = {[],[]};
   end

   if isfield(opts, 'v0') || isfield(opts, 'u0')
      if ~isfield(opts, 'v0'), opts.v0 = []; end
      if ~isfield(opts, 'u0'), opts.u0 = []; end
      init0 = {opts.v0, opts.u0};
      if isempty(init0{1})
         init0{1} = opts.matrixMatvec(init0{2}, 'notransp');
      elseif isempty(init{2})
         init0{2} = opts.matrixMatvec(init0{1}, 'transp');
      end
      if size(init0{1}, 1) ~= opts.m || size(init0{2}, 1) ~= opts.n || ...
         size(init0{1}, 2) ~= size(init0{2}, 2)
         error('Invalid matrix dimensions in opts.init');
      end
      opts = rmfield(opts, 'u0');
      opts = rmfield(opts, 'v0');
      opts.initSize = size(init0{1}, 2);
      init = {[init{1} init0{1}], [init{2} init0{2}]};
   end

   % Create primme_params
   primme_svds = primme_mex('primme_svds_initialize');

   % Set other options in primme_svds_params
   primme_svds_set_members(opts, primme_svds);

   % Set method in primme_svds_params
   primme_mex('primme_svds_set_method', method, primmeStage0method, ...
                                        primmeStage1method, primme_svds);

   % Set monitor and shared variables with the monitor
   hist = [];
   %locking = primme_mex('primme_get_member', primme, 'locking');
   nconv = 0;
   return_hist = 0;
   if dispLevel > 0
      % NOTE: Octave doesn't support function handler for nested functions
      primme_mex('primme_svds_set_member', primme_svds, 'monitorFun', ...
            @(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11)record_history(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11));
   end
   if nargout >= 5
      return_hist = 1;
   elseif dispLevel == 1
      fprintf('#MV\tTime\t\tNConv\tStage\n');
   elseif dispLevel == 2
      fprintf('#MV\tTime\t\tNConv\tStage\tIdx\tValue\tRes\n');
   elseif dispLevel == 3
      fprintf('#MV\tTime\t\tNConv\tStage\tIdx\tValue\tRes\tQMR_Res\n');
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
   xprimme_svds = [type 'primme_svds'];

   % Call xprimme_svds
   [ierr, svals, norms, svecsl, svecsr] = primme_mex(xprimme_svds, init{1}, ...
               init{2}, primme_svds); 

   % Process error code and return the required arguments
   if mod(ierr, -100) == -3 % if it is -3, -103 or -203
      warning([xprimme_svds ' returned ' num2str(ierr) ': ' primme_svds_error_msg(ierr)]);
   elseif ierr ~= 0
      error([xprimme_svds ' returned ' num2str(ierr) ': ' primme_svds_error_msg(ierr)]);
   end
   
   % Return smallest or interior singular triplets in descending order
   if strcmp(opts.target,'primme_svds_smallest') || strcmp(opts.target,'primme_svds_closest_abs')
      [svals,ind] = sort(svals,'descend');
      svecsl = svecsl(:,ind);
      svecsr = svecsr(:,ind);
   end

   if nargout <= 1
      varargout{1} = svals;
   elseif nargout == 2
      varargout{1} = svals;
      varargout{2} = norms;
   elseif nargout >= 3
      varargout{1} = svecsl;
      varargout{2} = diag(svals);
      varargout{3} = svecsr;
   end
   if (nargout >= 4)
      varargout{4} = norms;
   end
   if (nargout >= 5)
      stats = struct();
      stats.numMatvecs = primme_mex('primme_svds_get_member', primme_svds, 'stats_numMatvecs');
      stats.elapsedTime = primme_mex('primme_svds_get_member', primme_svds, 'stats_elapsedTime');
      stats.aNorm = primme_mex('primme_svds_get_member', primme_svds, 'aNorm');
      varargout{5} = stats;
   end
   if (nargout >= 6)
      varargout{6} = hist;
   end

   function record_history(basisSvals, basisFlags, iblock, basisNorms, ...
         numConverged, lockedSvals, lockedFlags, lockedNorms, inner_its, ...
         LSRes, event, stage)

      numMatvecs = double(primme_mex('primme_svds_get_member', primme_svds, 'stats_numMatvecs'));
      methodStage2 = double(primme_mex('primme_svds_get_member', primme_svds, 'methodStage2'));
      if stage == 0
         primme = primme_mex('primme_svds_get_member', primme_svds, 'primme');
      else
         primme = primme_mex('primme_svds_get_member', primme_svds, 'primmeStage2');
      end
      if stage == 0 && methodStage2 ~= 0
         locking = 1;
      else
         locking = primme_mex('primme_get_member', primme, 'locking');
      end
      maxInnerIterations = primme_mex('primme_get_member', primme, 'correction_maxInnerIterations');
      elapsedTime = primme_mex('primme_svds_get_member', primme_svds, 'stats_elapsedTime');
      hist_rows = size(hist, 1);
      if event == 0 || (event == 4 && ~locking) || event == 5
         if ~locking && ~isempty(numConverged)
            nconv = double(numConverged);
         elseif locking && ~isempty(lockedSvals)
            nconv = numel(lockedSvals);
         end
      end
      stage = double(stage) + 1;
      if dispLevel == 0
      elseif dispLevel == 1
         if (event == 4 && ~locking) || event == 5
            hist = [hist; numMatvecs elapsedTime nconv stage];
         end
      elseif dispLevel == 2
         if event == 0 || (nconv == opts.numSvals && ((event == 4 && ~locking) || event == 5))
            for i=1:numel(iblock)
               hist = [hist; numMatvecs elapsedTime nconv stage i basisSvals(iblock(i)+1) basisNorms(iblock(i)+1)];
            end
         end
      elseif dispLevel == 3
         if event == 1
            if ~isempty(basisSvals)
               value = basisSvals(iblock(1)+1);
               resNorm = basisNorms(iblock(1)+1);
            else
               value = nan;
               resNorm = nan;
            end
            hist = [hist; numMatvecs elapsedTime nconv stage nan value resNorm  LSRes];
         elseif (maxInnerIterations == 0 || nconv == opts.numSvals) && (event == 0 || ((event == 4 && ~locking) || event == 5))
            for i=1:numel(iblock)
               hist = [hist; numMatvecs elapsedTime nconv stage i basisSvals(iblock(i)+1) basisNorms(iblock(i)+1) nan];
            end
         end
      end
      if ~return_hist && size(hist,1) > hist_rows
         template{1} = '%d\t%f\t%d\t%d\n';
         template{2} = '%d\t%f\t%d\t%d\t%d\t%g\t%e\n';
         template{3} = '%d\t%f\t%d\t%d\t%d\t%g\t%e\t%e\n';
         for i=hist_rows+1:size(hist,1)
            a = num2cell(hist(i,:));
            fprintf(template{dispLevel}, a{:});
         end
         hist = [];
      end
   end
end

function [y] = matvecsvds(A, x, mode)
   if mode(1) == 'n'
      y = A*x;
   else
      y = A'*x;
   end
end

function [y] = precondsvds_Pmat(P, x, mode)
   if strcmp(mode, 'AHA')
      y = P\(P'\x);
   elseif strcmp(mode, 'AAH')
      y = P'\(P\x);
   else
      y = [P\x(size(P,1)+1:end,:); P'\x(1:size(P,1),:)];
   end
end

function [y] = precondsvds_Pfun(P, x, mode, m)
   if strcmp(mode, 'AHA')
      y = P(P(x, 'transp'), 'notransp');
   elseif strcmp(mode, 'AAH')
      y = P(P(x, 'notransp'), 'transp');
   else
      y = [P(x(m+1:end,:), 'notransp'); P(x(1:m,:), 'transp')];
   end
end

function [y] = precondsvds_P1P2(P1, P2, x, mode)
   if ~isempty(P1)
      if strcmp(mode, 'AHA')
         y = P2\(P1\(P1'\(P2'\x)));
      elseif strcmp(mode, 'AAH')
         y = P1'\(P2'\(P2\(P1\x)));
      else
         y = [P2\(P1\x(size(P1,1)+1:end,:)); P1'\(P2'\x(1:size(P1,1),:))];
      end
   else
      if strcmp(mode, 'AHA')
         y = P2\(P2'\x);
      elseif strcmp(mode, 'AAH')
         y = P2'\(P2\x);
      else
         y = x;
      end
   end
end

function [y] = precondsvds_Pstruct(P, x, mode)
   if isfield(P, mode)
      M = P.(mode);
      if isnumeric(M)
         y = M\x;
      else
         y = M(x);
      end
   else
      y = x;
   end
end


function [f] = fcnchk_gen(x)
   if exist('fcnchk', 'var')
      f = fcnchk(x);
   else 
      f = x;
   end
end

function primme_svds_set_members(opts, primme_svds, f, prefix)
%PRIMME_SVDS_SET_MEMBERS  Set options in primme_svds_params
%   PRIMME_SVDS_SET_MEMBERS(S, P) sets the options in struct S into the
%   primme_svds_params reference P.
%
%   Example:
%     primme_svds = primme_mex('primme_svds_initialize');
%     ops.n = 10;
%     ops.target = 'primme_svds_largest';
%     primme_svds_set_members(ops, primme_svds);

   % NOTE: MATLAB doesn't support default values in function
   %       declaration, Octave does.
   if nargin < 3, f = 'primme_svds_set_member'; end
   if nargin < 4, prefix = ''; end
   
   fields = fieldnames(opts);
   for i=1:numel(fields)
      value = getfield(opts, fields{i});
      label = fields{i};
      if isstruct(value) && ~strcmp('primme', label) && ~strcmp('primmeStage2', label)
         primme_svds_set_members(value, primme_svds, f, [prefix label '_']);
      elseif isstruct(value)
         primme0 = primme_mex('primme_svds_get_member', primme_svds, [prefix label]);
         primme_svds_set_members(value, primme0, 'primme_set_member');
      else
         try
            primme_mex(f, primme_svds, [prefix label], value);
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
   msg{39+  0} = 'success';
   msg{39+  1} = 'reported only amount of required memory';
   msg{39+ -1} = 'failed in allocating int or real workspace';
   msg{39+ -2} = 'malloc failed in allocating a permutation integer array';
   msg{39+ -3} = 'main_iter() encountered problem; the calling stack of the functions where the error occurred was printed in stderr';
   msg{39+ -4} = 'argument primme is NULL';
   msg{39+ -5} = 'n < 0 or nLocal < 0 or nLocal > n';
   msg{39+ -6} = 'numProcs' < 1';
   msg{39+ -7} = 'matrixMatvec is NULL';
   msg{39+ -8} = 'applyPreconditioner is NULL and precondition is not NULL';
   msg{39+ -9} = 'not used';
   msg{39+-10} = 'numEvals > n';
   msg{39+-11} = 'numEvals < 0';
   msg{39+-12} = 'eps > 0 and eps < machine precision';
   msg{39+-13} = 'target is not properly defined';
   msg{39+-14} = 'target is one of primme_largest_abs, primme_closest_geq, primme_closest_leq or primme_closest_abs but numTargetShifts <= 0 (no shifts)';
   msg{39+-15} = 'target is one of primme_largest_abs primme_closest_geq primme_closest_leq or primme_closest_abs but targetShifts is NULL  (no shifts array)';
   msg{39+-16} = 'numOrthoConst < 0 or numOrthoConst > n (no free dimensions left)';
   msg{39+-17} = 'maxBasisSize < 2';
   msg{39+-18} = 'minRestartSize < 0 or minRestartSize shouldn''t be zero';
   msg{39+-19} = 'maxBlockSize < 0 or maxBlockSize shouldn''t be zero';
   msg{39+-20} = 'maxPrevRetain < 0';
   msg{39+-21} = 'scheme is not one of *primme_thick* or *primme_dtr*';
   msg{39+-22} = 'initSize < 0';
   msg{39+-23} = 'locking == 0 and initSize > maxBasisSize';
   msg{39+-24} = 'locking and initSize > numEvals';
   msg{39+-25} = 'maxPrevRetain + minRestartSize >= maxBasisSize';
   msg{39+-26} = 'minRestartSize >= n';
   msg{39+-27} = 'printLevel < 0 or printLevel > 5';
   msg{39+-28} = 'convTest is not one of primme_full_LTolerance primme_decreasing_LTolerance primme_adaptive_ETolerance or primme_adaptive';
   msg{39+-29} = 'convTest == primme_decreasing_LTolerance and relTolBase <= 1';
   msg{39+-30} = 'evals is NULL, but not evecs and resNorms';
   msg{39+-31} = 'evecs is NULL, but not evals and resNorms';
   msg{39+-32} = 'resNorms is NULL, but not evecs and evals';
   msg{39+-33} = 'locking == 0 and minRestartSize < numEvals';
   msg{39+-34} = 'ldevecs is less than nLocal';
   msg{39+-35} = 'ldOPs is non-zero and less than nLocal';
   msg{39+-36} = 'not enough memory for realWork';
   msg{39+-37} = 'not enough memory for intWork';
   msg{39+-38} = 'locking == 0 and target is primme_closest_leq or primme_closet_geq';

   errorCode = errorCode + 39;
   if errorCode > 0 && errorCode <= numel(msg)
      s = msg{errorCode};
   else
      s = 'Unknown error code';
   end
end

function s = primme_svds_error_msg(errorCode)
   msg = {};
   msg{22+  0} = 'success';
   msg{22+  1} = 'reported only amount of required memory';
   msg{22+ -1} = 'failed in allocating int or real workspace';
   msg{22+ -2} = 'malloc failed in allocating a permutation integer array';
   msg{22+ -3} = 'main_iter() encountered problem; the calling stack of the functions where the error occurred was printed in stderr';
   msg{22+ -4} = 'primme_svds is NULL';
   msg{22+ -5} = 'Wrong value for m or n or mLocal or nLocal';
   msg{22+ -6} = 'Wrong value for numProcs';
   msg{22+ -7} = 'matrixMatvec is not set';
   msg{22+ -8} = 'applyPreconditioner is not set but precondition == 1 ';
   msg{22+ -9} = 'numProcs >1 but globalSumDouble is not set';
   msg{22+-10} = 'Wrong value for numSvals, it''s larger than min(m, n)';
   msg{22+-11} = 'Wrong value for numSvals, it''s smaller than 1';
   msg{22+-13} = 'Wrong value for target';
   msg{22+-14} = 'Wrong value for method';
   msg{22+-15} = 'Not supported combination of method and methodStage2';
   msg{22+-16} = 'Wrong value for printLevel';
   msg{22+-17} = 'svals is not set';
   msg{22+-18} = 'svecs is not set';
   msg{22+-19} = 'resNorms is not set';
   msg{22+-20} = 'not enough memory for realWork';
   msg{22+-21} = 'not enough memory for intWork';

   if errorCode >= -100
      errorCode = errorCode + 22;
      if errorCode > 0 && errorCode < numel(msg)
         s = msg{errorCode};
      else
         s = 'Unknown error code';
      end
   elseif errorCode >= -200
      s = ['Error from first stage: ' primme_error_msg(errorCode+100)];
   else
      s = ['Error from second stage: ' primme_error_msg(errorCode+200)];
   end
end
