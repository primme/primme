function [varargout] = primme_svds(varargin)
%PRIMME_SVDS   Find few singular values and vectors on large, sparse matrices
%   S = PRIMME_SVDS(A) computes the 6 largest singular values of A.
%
%   S = PRIMME_SVDS(AFUN,M,N) accepts the matvec function handle AFUN instead of
%   the matrix A. AFUN(X,'notransp') returns A*X while AFUN(X,'transp') returns
%   A'*X. The matrix A is M-by-N.
% 
%   S = PRIMME_SVDS(...,K) computes the K largest singular values of A.
%
%   S = PRIMME_SVDS(...,K,SIGMA) computes the K singular values closest to the
%   scalar shift SIGMA. If SIGMA is a vector, find the closest singular value to
%   each element in SIGMA. If SIGMA is 'L', it computes the largest singular
%   values; if SIGMA is 'S', it computes the smallest.
%
%   S = PRIMME_SVDS(A,K,SIGMA,OPTIONS) specifies extra solver parameters:
%
%   Field name       Parameter                               Default
%
%   OPTIONS.aNorm    estimation of the 2-norm A                    -
%   OPTIONS.tol      convergence tolerance (see eps):          1e-10
%                    NORM([A*V-U*S;A'*U-V*S]) <= tol * NORM(A).
%   OPTIONS.maxit    maximum number of iterat. (see maxMatvecs)  inf
%   OPTIONS.p        maximum basis size (see maxBasisSize)         -
%   OPTIONS.disp     level of message reporting (see printLevel)   0
%   OPTIONS.isreal   if 0, the matrix is complex; else it's real   1
%   OPTIONS.isdouble if 0, the matrix is single; else it's double  1
%   OPTIONS.method   which equivalent eigenproblemto to solve
%                    - 'primme_svds_normalequation': A'*A or A*A'
%                    - 'primme_svds_augmented': [0 A';A 0]
%                    - 'primme_svds_hybrid': first normal eq. and
%                      then augmented (default).                   
%   OPTIONS.v0       approx. left and right singular vectors {[],[]}
%   OPTIONS.orthoConst external orthogonalization constraints     [] 
%   OPTIONS.locking  1, hard locking; 0, soft locking              -
%   OPTIONS.maxBlockSize maximum block size                        1
%   OPTIONS.iseed    random seed
%   OPTIONS.primme   options for first stage solver                -
%   OPTIONS.primmeStage2 options for second stage solver           -
%
%   The available options for OPTIONS.primme and primmeStage1 are
%   the same as PRIMME_EIGS, plus the option 'method'. For detailed
%   descriptions of the above options, visit:
%   http://www.cs.wm.edu/~andreas/software/doc/svdsc.html#parameters-guide
%   and for further descriptions of the methods visit:
%   http://www.cs.wm.edu/~andreas/software/doc/appendixsvds.html#preset-methods
%
%   S = PRIMME_SVDS(...,K,SIGMA,OPTIONS,P)
%   S = PRIMME_SVDS(...,K,SIGMA,OPTIONS,P1,P2) makes use of a preconditioner,
%   applying P\X or (P1*P2)\X. If P is [] then a preconditioner is not
%   applied. P may be a function handle PFUN such that PFUN(X,'AHA')
%   returns an approximation of (A'*A)\X, PFUN(X,'AAH'), of (A*A')\X and
%   PFUN(X,'aug'), of [zeros(M,N) A;A' zeros(N,M)]\X.
%
%   [U,S,V] = PRIMME_SVDS(...) returns the singular vectors as well.
%   If A is M-by-N and K singular triplets are computed, then U is M-by-K
%   with orthonormal columns, S is K-by-K diagonal, and V is N-by-K with
%   orthonormal columns.
%
%   [S,R] = PRIMME_SVDS(...)
%   [U,S,V,R] = PRIMME_SVDS(...) returns upper bounds of the residual norm
%   of each K triplet, NORM([A*V(:,i)-S(i,i)*U(:,i); A'*U(:,i)-S(i,i)*V(:,i)]).
%
%   [U,S,V,R,STATS] = PRIMME_SVDS(...) returns how many times A and P were
%   used and elapsed time.
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
%      opts.method = 'primme_svds_normalequations' % set solver method
%      opts.primme.method = 'DEFAULT_MIN_TIME' % set first stage solver method
%      opts.primme.maxBlockSize = 2; % set block size for first stage
%      [u,s,v] = primme_svds(A,10,'S',opts)
%
%      opts.orthoConst = {u,v};  
%      [s,rnorms] = primme_svds(A,10,'S',opts) % find another 10
%
%      % Define a preconditioner only for first stage (A'*A)
%      Pstruct = struct('AHA', diag(A'*A),...
%                       'AAH', ones(200,1), 'aug', ones(250,1));
%      Pfun = @(x,mode)Pstruct.(mode).\x;
%      s = primme_svds(A,5,'S',[],Pfun) % find the 5 smallest values
%
%   For more details see PRIMME documentation at
%   http://www.cs.wm.edu/~andreas/software/doc/readme.html 
%
%   See also PRIMME_EIGS, SVDS.

   % Check primme_mex exists
   if ~ exist('primme_mex')
      error 'primme_mex is not available. Try to recompile the MATLAB''s PRIMME module'
   end

   % Check arity of input and output arguments
   minInputs = 1;
   maxInputs = 8;
   narginchk(minInputs,maxInputs);

   minOutputs = 0;
   maxOutputs = 5;
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
      Adouble = class(A) == 'double';
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
            error('target must be L, S or C');
         end
         opts.target = getfield(targets, target(1));
      elseif isnumeric(target)
         opts.targetShifts = target;
         opts.target = 'primme_svds_closest_abs';
      else
         error('target must be L, S or a real number');
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

   if nargin >= nextArg
      P = varargin{nextArg};
      if isnumeric(P)
         P = @(x,mode)precondsvds(P,x,mode);
      else
         P = fcnchk_gen(P); % get the function handle of user's function
      end
      opts.applyPreconditioner = P;
      opts.precondition = 1;
      nextArg = nextArg + 1;
   end
 
   % Test whether the given matrix and preconditioner are valid
   try
      x = opts.matrixMatvec(ones(opts.n, 1), 'notransp');
      x = opts.matrixMatvec(ones(opts.m, 1), 'transp');
      if isfield(opts, 'applyPreconditioner')
         x = opts.applyPreconditioner(ones(opts.n, 1), 'AHA');
         x = opts.applyPreconditioner(ones(opts.m, 1), 'AAH');
         x = opts.applyPreconditioner(ones(opts.m+opts.n, 1), 'aug');
      end
      clear x;
   catch ME
      rethrow(ME);
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

   % Rename tol, maxit, p and disp as eps, maxMatvecs, maxBasisSize and
   % printLevel
   changes = {{'tol', 'eps'}, {'maxit', 'maxMatvecs'}, {'p', 'maxBasisSize'}, ...
              {'disp', 'printLevel'}};
   for i=1:numel(changes)
      if isfield(opts, changes{i}{1})
         opts.(changes{i}{2}) = opts.(changes{i}{1});
         opts = rmfield(opts, changes{i}{1});
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
   primme_fields = {'primme', 'primmeStage1'};
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

   % Process method, primme.method and primmeStage1.method
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
   if isfield(opts, 'primmeStage1') && isfield(opts.primmeStage1, 'method')
      primmeStage1method = opts.primmeStage1.method;
      opts.primmeStage1 = rmfield(opts.primmeStage1, 'method');
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

   if isfield(opts, 'v0')
      init0 = opts.v0;
      if ~iscell(init0) || numel(init0) ~= 2 || (isempty(init0{1}) && isempty(init0{2}))
         error('opts.v0 should be {left_vectors, right_vectors}');
      end
      if isempty(init0{1})
         init0{1} = opts.matrixMatvec(init0{2}, 'notransp');
      elseif isempty(init{2})
         init0{2} = opts.matrixMatvec(init0{1}, 'transp');
      end
      if size(init0{1}, 1) ~= opts.m || size(init0{2}, 1) ~= opts.n || ...
         size(init0{1}, 2) ~= size(init0{2}, 2)
         error('Invalid matrix dimensions in opts.init');
      end
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
   if ierr ~= 0
      error([xprimme_svds ' returned ' num2str(ierr)]);
   end

   if nargout == 1
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
end

function [y] = matvecsvds(A, x, mode)
   if mode(1) == 'n'
      y = A*x;
   else
      y = A'*x;
   end
end

function [y] = precondsvds(P, x, mode)
   if mode == 'AHA'
      y = P\(P'\x);
   elseif mode == 'AAH'
      y = P'\(P\x);
   else
      y = [P'\(P\x(1:size(P,1),:)); P\(P'\x(size(P,1):end,:))];
   end
end

function [f] = fcnchk_gen(x, n)
   if exist('fcnchk')
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

   % NOTE: Expensive Mathworks' MATLAB doesn't support default values in function
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
