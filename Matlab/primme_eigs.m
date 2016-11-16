function [varargout] = primme_eigs(varargin)
%PRIMME_EIGS   Find few eigenvalues and vectors on large, sparse Hermitian matrices
%   PRIMME_EIGS finds few eigenvalues and eigenvectors of a real symmetric
%   or Hermitian matrix, A, by calling the PRIMME driver PRIMME_MEX.
%   Almost full PRIMME functionality is supported.
%
%   D = PRIMME_EIGS(A) returns a vector of A's 6 largest algebraic eigenvalues.
%
%   D = PRIMME_EIGS(AFUN,DIM) accepts a function AFUN instead of a matrix. AFUN
%   is a function handle and y = AFUN(x) returns the matrix-vector product A*x.
%
%   D = PRIMME_EIGS(...,K) finds the K largest magnitude eigenvalues. K must be
%   less than the dimension of the matrix A.
%
%   D = PRIMME_EIGS(...,K,TARGET) returns K eigenvalues from different parts of
%   the spectrum. If TARGET is
%       'LM' or 'SM', find eigenvalues D with the largest or smallest ABS(D-S)
%       'LA' or 'SA', find eigenvalues with the largest or smallest value
%       'CLT' or 'CGT', find eigenvalues closest to S but less or greater than S.
%   S is each value in OPTS.targetShifts.
%   If TARGET is a real number, it finds the closest eigenvalues to TARGET.
%
%   D = PRIMME_EIGS(...,K,TARGET,OPTS) specifies extra solver parameters. Some
%     default values are indicated in brackets {}:
%
%     OPTS.aNorm: the estimate 2-norm of matrix A {estimate the norm}
%     OPTS.tol: convergence tolerance: NORM(A*X(:,i)-X(:,i)*D(i,i)) < tol*NORM(A)
%     OPTS.maxBlockSize: maximum block size {1}
%     OPTS.disp: different level reporting(0-5) {no output 0}
%     OPTS.isreal: the complexity of A represented by AFUN {false}
%     OPTS.targetShifts: shifts for interior eigenvalues (see target)
%     OPTS.v0: approximate eigenvectors {[]}
%     OPTS.orthoConst: external orthogonalization constraints {[]}
%     OPTS.locking: 1, hard locking; 0, soft locking
%     OPTS.p: maximum size of the search subspace
%     OPTS.minRestartSize: minimum Ritz vectors to keep in restarting
%     OPTS.maxMatvecs: maximum number of matrix vector multiplications {Inf}
%     OPTS.maxit: maximum number of outer iterations {Inf}
%     OPTS.scheme: the restart scheme {'primme_thick'}
%     OPTS.maxPrevRetain: number of approximate eigenvectors retained from
%       previous iteration, that are kept after restart.
%     OPTS.robustShifts: set to true to avoid stagnation.
%     OPTS.maxInnerIterations: maximum number of inner QMR iterations
%     OPTS.LeftQ: use the locked vectors in the left projector
%     OPTS.LeftX: use the approx. eigenvector in the left projector
%     OPTS.RightQ: use the locked vectors in the right projector
%     OPTS.RightX: use the approx. eigenvector in the right projector
%     OPTS.SkewQ: use the preconditioned locked vectors in the right projector
%     OPTS.SkewX: use the preconditioned approx. eigenvector in the right projector
%     OPTS.relTolBase: a legacy from classical JDQR (recommend not use)
%     OPTS.convTest: how to stop the inner QMR Method
%     OPTS.iseed: random seed
%
%   For detailed descriptions of the above options, visit:
%   http://www.cs.wm.edu/~andreas/software/doc/primmec.html#parameters-guide
%
%   D = PRIMME_EIGS(...,K,TARGET,OPTS,METHOD) specifies the eigensolver method:
%     'DYNAMIC', (default)        switches dynamically to the best method
%     'DEFAULT_MIN_TIME',         best method for light matrix-vector product
%     'DEFAULT_MIN_MATVECS',      best method for heavy matvec/preconditioner
%     'Arnoldi',                  Arnoldi not implemented efficiently
%     'GD',                       classical block Generalized Davidson 
%     'GD_plusK',                 GD+k block GD with recurrence restarting
%     'GD_Olsen_plusK',           GD+k with approximate Olsen precond.
%     'JD_Olsen_plusK',           GD+k, exact Olsen (two precond per step)
%     'RQI',                      Rayleigh Quotient Iteration. Also INVIT,
%                                 but for INVIT provide targetShifts
%     'JDQR',                     Original block, Jacobi Davidson
%     'JDQMR',                    Our block JDQMR method (similar to JDCG)
%     'JDQMR_ETol',               Slight, but efficient JDQMR modification
%     'SUBSPACE_ITERATION',       equiv. to GD(block,2*block)
%     'LOBPCG_OrthoBasis',        equiv. to GD(nev,3*nev)+nev
%     'LOBPCG_OrthoBasis_Window'  equiv. to GD(block,3*block)+block nev>block
%
%   For further description of the method visit:
%   http://www.cs.wm.edu/~andreas/software/doc/appendix.html#preset-methods
%
%   D = PRIMME_EIGS(...,K,TARGET,OPTS,METHOD,P) 
%   D = PRIMME_EIGS(...,K,TARGET,OPTS,METHOD,P1,P2) uses preconditioner P or
%   P = P1*P2 to accelerate convergence of the method.
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
%   largest and the smallest algebraic eigenvalues on A.
%
%   Examples:
%      A = diag(1:100);
%
%      d = primme_eigs(A,10) % the 10 largest magnitude eigenvalues
%
%      d = primme_eigs(A,10,'SM') % the 10 smallest magnitude eigenvalues
%
%      d = primme_eigs(A,10,25) % the 10 closest eigenvalues to 25
%
%      opts = struct();
%      opts.tol = 1e-4; % set tolerance
%      opts.maxBlockSize = 2; % set block size
%      [x,d] = primme_eigs(A,10,'S',opts,'DEFAULT_MIN_TIME')
%
%      opts.orthoConst = x;  
%      [d,rnorms] = primme_eigs(A,10,'S',opts) % find another 10
%
%      % Build a Jacobi preconditioner (too convenient for a diagonal matrix!)
%      Pfun = @(x)(diag(A) - 30.5)\x;
%      d = primme_eigs(A,5,30.5,[],[],Pfun) % find the closest 5 to 30.5
%
%   For more details see PRIMME documentation at
%   http://www.cs.wm.edu/~andreas/software/doc/readme.html 
%
%   See also PRIMME_SVDS, EIGS.

   % Check primme_mex exists
   if ~ exist('primme_mex')
      error 'primme_mex is not available. Try to recompile the MATLAB''s PRIMME module'
   end

   % Check arity of input and output arguments
   minInputs = 1;
   maxInputs = 8;
   narginchk(minInputs,maxInputs);

   minOutputs = 0;
   maxOutputs = 4;
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
      Adouble = class(A) == 'double';
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
      nextArg = nextArg + 1;
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

   % Rename tol, maxit, p and disp as eps, maxOuterIterations, maxBasisSize and
   % printLevel. Also move options that are outside of primme_params' hierarchy.
   changes = {{'tol', 'eps'}, {'maxit', 'maxOuterIterations'}, {'p', 'maxBasisSize'}, ...
              {'disp', 'printLevel'}, ...
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

   if isfield(opts, 'init')
      init0 = opts.init;
      if size(init0, 1) ~= opts.n
         error('Invalid matrix dimensions in opts.init');
      end
      opts = rmfield(opts, 'init');
      opts.initSize = size(init0, 2);
      init = [init init0];
   end

   % Default printLevel is 0
   if ~isfield(opts, 'printLevel')
      opts.printLevel = 0;
   end

   % Create primme_params
   primme = primme_mex('primme_initialize');

   % Set other options in primme_params
   primme_set_members(opts, primme);

   % Set method
   primme_mex('primme_set_method', method, primme);

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
      error([xprimme ' returned ' num2str(ierr)]);
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
end

function [f] = fcnchk_gen(x, n)
   if exist('fcnchk')
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
