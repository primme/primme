function [varargout] = primme_eigs(varargin)

%  primme_eigs finds a few eigenvalues and eigenvectors. If A is M-by-N, 
%  a few eigenvalues and eigenvectors of A are returned by PRIMME_mex.
%  (flag,dim,...), where input matrix is A.
% 
%  It is an interface between user function call and PRIMME_mex function. 
%  The primme_eigs function requires at least one input (input matrix A) 
%  and output at least eigenvalues of A.
%  
%  Like eigs() function in the matlab, we provide different level function
%  calls to satisfy users' demands:
%  
%  Input: [A, numEvals, target, opts, eigsMethod, P]
% 
%  Output: [evals, evecs, norms, primmeout]
%  
%  Function call:  
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
%  primme_D = primme_eigs(A) returns a vector of A's 6 largest algebraic
%  eigenvalues. A must be real symmetric or complex hermitian and should
%  be large and sparse. The primme_eigs(Afun, dim) accepts a function 
%  AFUN instead of the matrix A. AFUN is a function handle and y = Afun(x) 
%  returns the matrix-vector product A*x. In all these primme_eigs
%  function syntaxes, primme_eigs(A,...) could be replaced by
%  primme_eigs(Afun, dim,...). More examples about how to use function Afun
%  are presented in the Primme_eigsTest6, 7, 8, and 9 in the root directory
%  of PRIMME_MEX folder.
%
%  [primme_V, primme_D] = primme_eigs(A) returns a diagonal matrix of
%  primme_D of A's 6 largest algebraic eigenvalues and a matrix primme_V
%  whose columns are the corresponding eigenvectors. 
% 
%  [primme_V, primme_D, norms, primmeout] = primme_eigs(A, numEvals)
%  returns a diagonal matrix of primme_D of A's numEvals largest algebraic 
%  eigenvalues, a matrix primme_V whose columns are the corresponding 
%  eigenvectors, a double array of the residual norm of eigenvalues  and
%  a struct to report statistical information about numOuterIterations,
%  numRestarts, numMatvecs and numPreconds. numEvals is the number of 
%  eigenvalues that users want to find. It must be less than dimention of 
%  the matrix A.
%
%  primme_eigs(A, numEvals, target) returns numEvals target eigenvlaues.
%  target could be a string like below:
%     'LA' ------ primme_largest (default)
%     'SA' ------ primme_smallest    
%     'CGT'------ primme_closest_geq  
%     'CLT'------ primme_closest_leq  
%     'CT' ------ primme_closest_abs   
%
%  primme_eigs(A, numEvals, target, opts, eigsMethod) specify options that 
%  are listed and explained in the last. eigsMethod is the solver method
%  the PRIMME uses to find eigenvalues and eigenvectors. eigsMethod could
%  be:
%  typedef enum{
%  DYNAMIC, (default)        ---0: Switches dynamically to the best method
%  DEFAULT_MIN_TIME,         ---1: Currently set at JDQMR_ETol
%  DEFAULT_MIN_MATVECS,      ---2: Currently set at GD_Olsen_plusK
%  Arnoldi,                  ---3: obviously not an efficient choice 
%  GD,                       ---4: classical block Generalized Davidson 
%  GD_plusK,                 ---5: GD+k block GD with recurrence restarting
%  GD_Olsen_plusK,           ---6: GD+k with approximate Olsen precond.
%  JD_Olsen_plusK,           ---7: GD+k, exact Olsen (two precond per step)
%  RQI,                      ---8: Rayleigh Quotient Iteration. Also INVIT,
%                                :   but for INVIT provide targetShifts
%  JDQR,                     ---9: Original block, Jacobi Davidson
%  JDQMR,                   ---10: Our block JDQMR method (similar to JDCG)
%  JDQMR_ETol,              ---11: Slight, but efficient JDQMR modification
%  SUBSPACE_ITERATION,      ---12: equiv. to GD(block,2*block)
%  LOBPCG_OrthoBasis,       ---13: equiv. to GD(nev,3*nev)+nev
%  LOBPCG_OrthoBasis_Window ---14: equiv. to GD(block,3*block)+block nev>block
%
%  primme_eigs(A, numEvals, target, opts, eigsMethod, P) uses preconditioner 
%  P or P = P1*P2 to solve eigenvalue problem for large sparse matrix.
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
%  opts.isreal: the complexity of A represented by afun [{ture}|false]
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

%  clear global;
   clear global eigsFunCallFlag;
   clear functions;

   minInputs = 1;
   maxInputs = 8;
   narginchk(minInputs,maxInputs);
   minOutputs = 0;
   maxOutputs = 4;
   nargoutchk(minOutputs,maxOutputs);

   global eigsFunCallFlag; % mark if user called primme_eigs function
   if isempty(eigsFunCallFlag) 
      eigsFunCallFlag = true;
   end

    global primmeA;  
    global Amatrix;   % mark if A is a marix or matrix function 
    global P1;  
    global P1matrix;  % mark if the first preconditioner is a marix or a function
    global P2;  % P2 is the second preconditioner that must be a matrix

    if isfloat(varargin{1}) % check if the first input is matrix or matrix funtion
        primmeA = varargin{1};
        Amatrix = true;
    else
        % By checking the function A with fcnchk, we can now use direct
        % function evaluation on the result, without resorting to feval
        primmeA = fcnchk(varargin{1}); % get the function handle of user's funciton
        Amatrix = false;
    end

   if (Amatrix == 1) % A is an matrix
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
           if(ischar(varargin{3}))
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

       if (nargin >= 6) % check if the sixth input is matrix or matrix funtion      
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
            else % the sixth input is matrix function
                % By checking the function P1 with fcnchk, we can now use direct
                % function evaluation on the result, without resorting to feval
                P1 = fcnchk(varargin{6}); % get the function handle of user's funciton
                P1matrix = false;
            end
       end
   else  % A is a matrix function
       dim = varargin{2};
       if ~isscalar(dim) || ~isreal(dim) || (dim<0) || ~isfinite(dim)
            error(message('The size of input matrix A must be an positive integer'));
       end            
       if (nargin >= 3)
            numEvals = varargin{3};
       else
            numEvals = min(6,dim);
       end
       if (nargin >= 4)
           if(ischar(varargin{4}))
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
            else % the seventh input is matrix function
                % By checking the function P1 with fcnchk, we can now use direct
                % function evaluation on the result, without resorting to feval
                % get the function handle of user's funciton
                P1 = fcnchk(varargin{7}); 
                P1matrix = false;
            end
       end

       % if use does not specifies the field "opts.isreal"
       if(isfield(opts,'isreal') == 0) 
           flag = 1; % The complexity of matrix function is default real.
       else
            if(opts.isreal == 1)
                flag = 1;
            else
                flag = 1+1i;
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
      
function y = getMatvecHandle(x)
%  PRIMME_MEX calls this function to get corresponding MATVEC function 
%  handle, then passes the blocked vector x to right MATVEC function and 
%  return the result y of MATRIX-VECTOR operations to PRIMME_MEX.  
%  Detailed explanation goes here

    global primmeA; % primmeA is a matrix or a matrix function
    global Amatrix; % mark if A is  a matrix or a matrix function
    global eigsFunCallFlag; % mark if primme_egis is called by users
    
    % check user calls primme_eigs 
    if  ~isempty(eigsFunCallFlag) 
        if Amatrix % primmeA is a matrix for primme_eigs
            y = primmeA * x;
        else % primmeA is a matrix function for primme_eigs
            y = primmeA(x);
        end
    end
end

function y = getPrecondHandle(x)
%  PRIMME_MEX calls this function to get corresponding precondition function 
%  handle, then passes the blocked vector x to right preconditioning  
%  function and return the result y of MATRIX-VECTOR operations to PRIMME_MEX.  
%  Detailed explanation goes here

    global P1;        % P1 is the first preconditioner matrix or function
    global P1matrix;  
    %     P1matrix = 0, First preconditioner is a matrix funciton
    %     P1matrix = 1, First preconditioner is matrix for A
    %     P1matrix = 2, First preconditioner is matrix directly for ATA or OAAO                
    global P2;        % P2 is the second preconditioner matrix
    global eigsFunCallFlag; % mark if primme_egis is called by users
        
    % check user calls primme_eigs
    if  ~isempty(eigsFunCallFlag)
        if P1matrix % P1 is a matrix for primme_eigs
            y = P1\x;
            if  ~isempty(P2) % P2 exists and is a matrix for primme_eigs
                y = P2\y;
            end
        else % P1 is a matrix function for primme_eigs
            y = P1(x);
        end
    end
end
