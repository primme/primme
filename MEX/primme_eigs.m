function [varargout] = primme_eigs(varargin)

%  primme_eigs finds a few eigenvalues and eigenvectors of a real symmetric
%  or Hermitian matrix, A, by calling the function PRIMME_mex(flag,dim,...).
%  This in turn calls PRIMME. Full PRIMME functionality is supported.
%  
%  Input: [A, numEvals, target, opts, eigsMethod, P]
% 
%  Output: [evals, evecs, norms, primmeout]
%  
%  We provide different levels of function calls, similarly to MATLAB eigs():
%
%   primme_eigs(A)
%   primme_eigs(A, numEvals)
%   primme_eigs(A, numEvals, target)
%   primme_eigs(A, numEvals, target, opts)
%   primme_eigs(A, numEvals, target, opts, eigsMethod)
%   primme_eigs(A, numEvals, target, opts, eigsMethod, P)
%   primme_eigs(A, numEvals, target, opts, eigsMethod, P1,P2)
%   primme_eigs(A, numEvals, target, opts, eigsMethod, Pfun)
%   primme_eigs(Afun, dim,...)
%
%  primme_eigs(A) returns a vector of A's 6 largest algebraic eigenvalues.
%  A must be real symmetric or complex Hermitian and should be large and sparse.
%
%  primme_eigs(Afun, dim) accepts a function AFUN instead of a matrix. AFUN is
%  a function handle and y = Afun(x) returns the matrix-vector product A*x.
%  primme_eigs(A,...) could be replaced by primme_eigs(Afun, dim,...) in any of
%  above levels of function calls. Examples are given in PRIMME_MEX_Readme.txt
%  in the root directory of PRIMME_MEX folder.
%
%  [V, D] = primme_eigs(A) returns a diagonal matrix D, of A's 6 largest 
%  algebraic eigenvalues and a matrix V whose columns are the corresponding
%  eigenvectors. 
% 
%  [V, D, norms, primmeout] = primme_eigs(A) also returns an array of the 
%  residual norms of the computed eigenpairs, and a struct to report statistical
%  information about numOuterIterations, numRestarts, numMatvecs and numPreconds.
%
%  primme_eigs(A, numEvals) finds the numEvals largest algebraic eigenvalues.
%  numEvals must be less than the dimension of the matrix A.
%
%  primme_eigs(A, numEvals, target) returns numEvals target eigenvalues.
%  target could be a string like below:
%     'LA' ------ primme_largest (default)
%     'SA' ------ primme_smallest    
%     'CGT'------ primme_closest_geq  
%     'CLT'------ primme_closest_leq  
%     'CT' ------ primme_closest_abs   
%
%  primme_eigs(A, numEvals, target, opts, eigsMethod) specifies any of a 
%  set of possible options as explained below in the opts structure. 
%
%  eigsMethod is an integer specifying one of the preset methods in PRIMME:
%  eigsMethod   corresponding PRIMME method
%    0:    DYNAMIC, (default)        Switches dynamically to the best method
%    1:    DEFAULT_MIN_TIME,         Currently set at JDQMR_ETol
%    2:    DEFAULT_MIN_MATVECS,      Currently set at GD+block
%    3:    Arnoldi,                  obviously not an efficient choice 
%    4:    GD,                       classical block Generalized Davidson 
%    5:    GD_plusK,                 GD+k block GD with recurrence restarting
%    6:    GD_Olsen_plusK,           GD+k with approximate Olsen precond.
%    7:    JD_Olsen_plusK,           GD+k, exact Olsen (two precond per step)
%    8:    RQI,                      Rayleigh Quotient Iteration. Also INVIT,
%                                    but for INVIT provide targetShifts
%    9:    JDQR,                     Original block, Jacobi Davidson
%    10:   JDQMR,                    Our block JDQMR method (similar to JDCG)
%    11:   JDQMR_ETol,               Slight, but efficient JDQMR modification
%    12:   SUBSPACE_ITERATION,       equiv. to GD(block,2*block)
%    13:   LOBPCG_OrthoBasis,        equiv. to GD(nev,3*nev)+nev
%    14:   LOBPCG_OrthoBasis_Window  equiv. to GD(block,3*block)+block nev>block
%
%  primme_eigs(A, numEvals, target, opts, eigsMethod, P) 
%  primme_eigs(A, numEvals, target, opts, eigsMethod, P1, P2) uses 
%  preconditioner P or P = P1*P2 to accelerate convergence of the methods.
%  If P is [] then a preconditioner is not applied. P may be a function 
%  handle Pfun such that Pfun(x) returns P\x.
%
%  opts is an option structure which contain following parameters:         
%  opts.aNorm: the estimate norm value of matrix A [{0.0}|scaler]
%  opts.eps: desired computing accuracy [{1e-12}|scaler]
%  opts.maxBlockSize: maximum block size the PRIMME uses [{1}|scaler]
%  opts.printLevel: different level reporting(0-5) [{1}|scaler]
%  opts.outputFile: output file name where user wants to save results
%  opts.precondition: set to 1 if use preconditioner [{0}|1]
%  opts.isreal: the complexity of A represented by AFUN [{ture}|false]
%  opts.numTargetShifts: number of shifts for interior eigenvalues [{0}|scaler]
%  opts.targetShifts: shifts for interior eigenvalues [{}|vector]
%  opts.initSize: On INPUT, the number of initial guesses provided in evecs 
%  array. ON OUTPUT, the number of converged eigenpairs [{0}|scaler]
%  opts.numOrthoConst: Number of external orthogonalization constraints 
%  provided in the first numOrthoConst vectors of evecs [{0}|scaler]
%  opts.locking: If set to 1, hard locking will be used, otherwise the code
%  will try to use soft locking [{0}|1]
%  opts.maxBasisSize: maximum basis size allowed in the main iteration
%  opts.minRestartSize: minimum Ritz vectors to restart
%  opts.maxMatvecs: maximum number of matrix vector multiplications
%  [{INT_MAX}|scaler]
%  opts.maxOuterIterations: maximum number of outer iterations
%  [{INT_MAX}|scaler]
%  opts.restartingParams.scheme: the restart scheme [{primme_thick}| primme_dtr]
%  opts.restartingParams.maxPrevRetain: number of approximations from 
%  previous iteration to be retained after restart [{1}|scaler]
%  opts.robustShifts: set to 1 if use robustShifting to help avoid
%  stagnation and misconverge [{0}|1]
%  opts.maxInnerIterations: number of inner QMR iterations [{0}|scaler]
%  opts.LeftQ: a projector with Q must be applied on the left [{0}|1]
%  opts.LeftX: a projector with X must be applied on the left [{0}|1]
%  opts.RightQ: a projector with Q must be applied on the right [{0}|1]
%  opts.RightX: a projector with X must be applied on the right [{0}|1]
%  opts.SkewQ: the Q right projector must be skew [{0}|1]
%  opts.SkewX: the X right projector must be skew [{0}|1]
%  opts.relTolBase: a legacy from calssical JDQR (recommend not use)
%  opts.convTest: how to stop the inner QMR Method
%  opts.iseed: set iseed value for initialization
%  opts.intWorkSize: memory size for integer workspace
%  opts.realWorkSize: memory size for real or complex workspace
%
%  See also PRIMME_MEX_Readme.txt
%  For more details on PRIMME's functionality see PRIMMEtopdir/readme.txt 

   clear global eigsFunCallFlag;
   clear functions;

   minInputs = 1;
   maxInputs = 8;
   narginchk(minInputs,maxInputs);
   minOutputs = 0;
   maxOutputs = 4;
   nargoutchk(minOutputs,maxOutputs);

   global eigsFunCallFlag; % flag if user called primme_eigs function
   if isempty(eigsFunCallFlag) 
      eigsFunCallFlag = true;
   end

    global primmeA;  
    global Amatrix;   % flag if A is a matrix (otherwise a matrix function)
    global P1;  
    global P1matrix;  % flag if P1 preconditioner is a matrix (or a function)
    global P2;        % if P2 second preconditioner is given it must be a matrix

    if isfloat(varargin{1}) % check if the first input is matrix or matrix funtion
        primmeA = varargin{1};
        Amatrix = true;
    else
        % By checking the function A with fcnchk, we can now use direct
        % function evaluation on the result, without resorting to feval
        primmeA = fcnchk(varargin{1}); % get the function handle of user's funciton
        Amatrix = false;
    end

   if (Amatrix) % A is an matrix
       [M, N] = size(primmeA); % get the dimension of input matrix A
       if isequal(primmeA,primmeA')
          dim = M;
       else
          error('Input matrix must be real symmetric or complex hermitian');
       end

       if isreal(primmeA)   % mark input matrix A is real or complex
          flag = 1;
       else
          flag = 1+1i;
       end

       if (nargin >= 2)
            numEvals = varargin{2};
       else
            numEvals = min(6,dim);
       end
       if (nargin >= 3)
           if (ischar(varargin{3}))
                target = varargin{3};
           else
                error('target must be a string');
           end
       else
           target = 'LA';
       end             
       if (nargin >= 4)
           if (isstruct(varargin{4}))
               opts = varargin{4};
           else
               error('opts must be a struct');
           end
       else
           opts = struct();
       end                   
       if (nargin >= 5)
           if (isreal(varargin{5}) && 0<= varargin{5} <=14)
               eigsMethod = varargin{5};                       
           else
               error('eigsMethod must be a real number between 0 and 14');
           end
       else
           eigsMethod = 0; 
       end

       if (nargin >= 6) % check if the sixth input is matrix or matrix function      
            if isfloat(varargin{6}) % the sixth input is matrix
                P1 = varargin{6};
                P1matrix = true;
                if (nargin == 7)
                    if isfloat(varargin{7}) % check if the seventh input is matrix
                        P2 = varargin{7};                                 
                    else
                        error('The second preconditioner must be a matrix');
                    end
                end
            else% the sixth input is matrix function
                % By checking the function P1 with fcnchk, we can now use direct
                % function evaluation on the result, without resorting to feval
                P1 = fcnchk(varargin{6}); % get the function handle of user's function
                P1matrix = false;
            end
       end
   else  % A is a matrix function
       dim = varargin{2};
       if ~isscalar(dim) || ~isreal(dim) || (dim<0) || ~isfinite(dim)
            error(message('The size of input matrix A must be an positive integer'));
       end            
       % Test whether the given matvec function is valid
       try
          xvec = randn(dim,1);
          yvev = primmeA(xvec);
          clear xvec, yvec;
       catch ME
          if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
             disp('Matvec function AFUN does not exist.');
             rethrow(ME);
          end
       end
       if (nargin >= 3)
            numEvals = varargin{3};
       else
            numEvals = min(6,dim);
       end
       if (nargin >= 4)
           if (ischar(varargin{4}))
                target = varargin{4};
           else
                error('target must be a string');
           end
       else
           target = 'LA';
       end             
       if (nargin >= 5)
           if (isstruct(varargin{5}))
               opts = varargin{5};
           else
               error('opts must be a struct');
           end
       else
           opts = struct();
       end                   
       if (nargin >= 6)
           if (isreal(varargin{6}) && 0<= varargin{6} <=14)
               eigsMethod = varargin{6};                       
           else
               error('eigsMethod must be a real number between 0 and 14');
           end
       else
           eigsMethod = 0; 
       end
       if (nargin >= 7) % check if the seventh input is matrix or matrix funtion
            if isfloat(varargin{7}) % the seventh input is matrix 
                P1 = varargin{7};
                P1matrix = true;
                if (nargin == 8)
                    if isfloat(varargin{8}) % check if the eighth input is matrix
                        P2 = varargin{8};
                    else
                        error('The second preconditioner must be a matrix');
                    end
                end
            else% the seventh input is matrix function
                % By checking the function P1 with fcnchk, we can now use direct
                % function evaluation on the result, without resorting to feval
                % get the function handle of user's funciton
                P1 = fcnchk(varargin{7}); 
                P1matrix = false;
            end
       end

       % if the user does not specifies the field "opts.isreal"
       if (isfield(opts,'isreal') == 0) 
           flag = 1; % The complexity of matrix function is default real.
       else
            if (opts.isreal)
                flag = 1;
            else
                flag = 1+1i;
            end
       end
   end

   if (~isempty(P1) && ~P1matrix)
       % Test whether the given preconditioner function is valid
       try
          xvec = randn(dim,1);
          yvev = P1(xvec);
          clear xvec, yvec;
       catch ME
          if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
             disp('Preconditioner P1 function does not exist.');
             rethrow(ME);
          end
       end
   end

   [evals,evecs,norms,primmeout]= PRIMME_mex(flag, dim, numEvals, target, eigsMethod, opts); 

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
       varargout{4} = primmeout;
   end

end
