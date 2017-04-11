function [varargout] = primme_eigs(varargin)
%PRIMME_EIGS  Find a few eigenvalues/vectors of large, sparse Hermitian matrices
%
%   D = PRIMME_EIGS(A) returns a vector of A's 6 largest magnitude eigenvalues.
%
%   D = PRIMME_EIGS(AFUN,DIM) accepts a function AFUN instead of a matrix. AFUN
%   is a function handle and y = AFUN(x) returns the matrix-vector product A*x.
%   In all the following syntaxes, A can be replaced by AFUN,DIM.
%
%   D = PRIMME_EIGS(A,K) finds the K largest magnitude eigenvalues. K must be
%   less than the dimension of the matrix A.
%
%   D = PRIMME_EIGS(A,K,TARGET) returns K eigenvalues such that: 
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
%   D = PRIMME_EIGS(A,K,TARGET,OPTS) specifies extra solver parameters. Some
%     default values are indicated in brackets {}:
%
%     OPTS.aNorm: the estimated 2-norm of A {0.0 (estimate the norm internally)}
%     OPTS.tol: convergence tolerance:                      {eps*1e4}
%                NORM(A*X(:,i)-X(:,i)*D(i,i)) < tol*NORM(A)
%     OPTS.maxBlockSize: maximum block size (useful for high multiplicities) {1}
%     OPTS.disp: different level reporting (0-3) (see HIST) {no output 0}
%     OPTS.isreal: whether A represented by AFUN is real or complex {false}
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
%     OPTS.iseed: random seed
%
%   For detailed descriptions of the above options, visit:
%   http://www.cs.wm.edu/~andreas/software/doc/primmec.html#parameters-guide
%
%   D = PRIMME_EIGS(A,K,TARGET,OPTS,METHOD) specifies the eigensolver method:
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
%   D = PRIMME_EIGS(A,K,TARGET,OPTS,METHOD,P) 
%   D = PRIMME_EIGS(A,K,TARGET,OPTS,METHOD,P1,P2) uses preconditioner P or
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
%   largest and smallest algebraic eigenvalues of A.
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
%   has one row per converged eigenpair and only the first three columns are
%   reported; if OPTS.disp == 2, HIST has one row per outer iteration and only
%   the first six columns are reported; and otherwise HIST has one row per QMR
%   iteration and all columns are reported.
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
      error 'primme_mex is not available. Try to recompile the MATLAB/Octave''s PRIMME module'
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
      if ~ishermitian(A)
         error('Input matrix must be real symmetric or complex Hermitian');
      end
      opts.n = n;
      opts.matrixMatvec = @(x)A*x;

      % Get type and complexity
      Acomplex = ~isreal(A);
      Adouble = strcmp(class(A), 'double');
   else
      opts.matrixMatvec = fcnchk_gen(A); % get the function handle of user's function
      n = round(varargin{nextArg});
      if ~isscalar(n) || ~isreal(n) || (n<0) || ~isfinite(n)
         error(message('The size of input matrix A must be an positive integer'));
      end
      opts.n = n;
      nextArg = nextArg + 1;

      % Assume complex double matrix
      Acomplex = 1;
      Adouble = 1;
   end

   if nargin >= nextArg
      opts.numEvals = round(varargin{nextArg});
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
 
   % Test whether the given matrix and preconditioner are valid
   try
      x = opts.matrixMatvec(ones(opts.n, 1));
      if isfield(opts, 'applyPreconditioner')
         x = opts.applyPreconditioner(ones(opts.n, 1));
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

   % Process 'disp' in opts
   if isfield(opts, 'disp')
      dispLevel = opts.disp;
      if dispLevel > 3 || dispLevel < 0
         error('Invalid value in opts.disp; it should be 0, 1, 2 or 3');
      end
      opts = rmfield(opts, 'disp');
   elseif nargout >= 5
      dispLevel = 1;
   else
      dispLevel = 0;
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
      if Adouble
         opts.eps = eps*1e4;
      else
         opts.eps = sqrt(eps)*1e4;
      end
   end 

   % Create primme_params
   primme = primme_mex('primme_initialize');

   % Set other options in primme_params
   primme_set_members(opts, primme);

   % Set method
   primme_mex('primme_set_method', method, primme);

   % Set monitor and shared variables with the monitor
   hist = [];
   locking = primme_mex('primme_get_member', primme, 'locking');
   nconv = [];
   return_hist = 0;
   if dispLevel > 0
      % NOTE: Octave doesn't support function handler for nested functions
      primme_mex('primme_set_member', primme, 'monitorFun', ...
            @(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)record_history(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10));
   end
   if nargout >= 5
      return_hist = 1;
   elseif dispLevel == 1
      fprintf('#MV\tTime\t\tNConv\n');
   elseif dispLevel == 2
      fprintf('#MV\tTime\t\tNConv\tIdx\tValue\tRes\n');
   elseif dispLevel == 3
      fprintf('#MV\tTime\t\tNConv\tIdx\tValue\tRes\tQMR_Res\n');
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
   xprimme = [type 'primme'];

   % Call xprimme
   [ierr, evals, norms, evecs] = primme_mex(xprimme, init, primme); 

   % Process error code and return the required arguments
   if ierr ~= 0
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
      stats = struct();
      stats.numMatvecs = primme_mex('primme_get_member', primme, 'stats_numMatvecs');
      stats.elapsedTime = primme_mex('primme_get_member', primme, 'stats_elapsedTime');
      stats.estimateMinEVal = primme_mex('primme_get_member', primme, 'stats_estimateMinEVal');
      stats.estimateMaxEVal = primme_mex('primme_get_member', primme, 'stats_estimateMaxEVal');
      stats.estimateAnorm = primme_mex('primme_get_member', primme, 'stats_estimateLargestSVal');
      varargout{4} = stats;
   end
   if (nargout >= 5)
      varargout{5} = hist;
   end

   function record_history(basisEvals, basisFlags, iblock, basisNorms, ...
         numConverged, lockedEvals, lockedFlags, lockedNorms, inner_its, ...
         LSRes, event)

      numMatvecs = double(primme_mex('primme_get_member', primme, 'stats_numMatvecs'));
      maxInnerIterations = primme_mex('primme_get_member', primme, 'correction_maxInnerIterations');
      elapsedTime = primme_mex('primme_get_member', primme, 'stats_elapsedTime');
      hist_rows = size(hist, 1);
      if event == 0 || (event == 4 && ~locking) || event == 5
         if ~locking
            nconv = double(numConverged);
         else
            nconv = numel(lockedEvals);
         end
      end
      if dispLevel == 0
      elseif dispLevel == 1
         if (event == 4 && ~locking) || event == 5
            hist = [hist; numMatvecs elapsedTime nconv];
         end
      elseif dispLevel == 2
         if event == 0 || (nconv == opts.numEvals && ((event == 4 && ~locking) || event == 5))
            for i=1:numel(iblock)
               hist = [hist; numMatvecs elapsedTime nconv i basisEvals(iblock(i)+1) basisNorms(iblock(i)+1)];
            end
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
            hist = [hist; numMatvecs elapsedTime nconv nan value resNorm  LSRes];
         elseif (maxInnerIterations == 0 || nconv == opts.numEvals) && (event == 0 || ((event == 4 && ~locking) || event == 5))
            for i=1:numel(iblock)
               hist = [hist; numMatvecs elapsedTime nconv i basisEvals(iblock(i)+1) basisNorms(iblock(i)+1) nan];
            end
         end
      end
      if ~return_hist && size(hist,1) > hist_rows
         template{1} = '%d\t%f\t%d\n';
         template{2} = '%d\t%f\t%d\t%d\t%g\t%e\n';
         template{3} = '%d\t%f\t%d\t%d\t%g\t%e\t%e\n';
         for i=hist_rows+1:size(hist,1)
            a = num2cell(hist(i,:));
            fprintf(template{dispLevel}, a{:});
         end
         hist = [];
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
   msg{39+-38} = '"locking == 0 and target is primme_closest_leq or primme_closet_geq';

   errorCode = errorCode + 39;
   if errorCode > 0 && errorCode <= numel(msg)
      s = msg{errorCode};
   else
      s = 'Unknown error code';
   end
end
