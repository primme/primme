.. highlight:: matlab

MATLAB Interface
----------------

.. mat:function:: function [varargout] = primme_eigs(varargin)

   :mat:func:`primme_eigs` finds a few eigenvalues and eigenvectors of a real symmetric
   or Hermitian matrix, A, by calling the function ``PRIMME_mex`` (flag,dim,...).
   This in turn calls PRIMME. Full PRIMME functionality is supported.
   
   Input: [A, numEvals, target, opts, eigsMethod, P]
  
   Output: [evals, evecs, norms, primmeout]
   
   We provide different levels of function calls, similarly to MATLAB eigs()::
 
      primme_eigs(A)
      primme_eigs(A, numEvals)
      primme_eigs(A, numEvals, target)
      primme_eigs(A, numEvals, target, opts)
      primme_eigs(A, numEvals, target, opts, eigsMethod)
      primme_eigs(A, numEvals, target, opts, eigsMethod, P)
      primme_eigs(A, numEvals, target, opts, eigsMethod, P1,P2)
      primme_eigs(A, numEvals, target, opts, eigsMethod, Pfun)
      primme_eigs(Afun, dim,...)
 
   ``primme_eigs(A)`` returns a vector of A's 6 largest algebraic eigenvalues.
   A must be real symmetric or complex Hermitian and should be large and sparse.
 
   ``primme_eigs(Afun, dim)`` accepts a function AFUN instead of a matrix. AFUN is
   a function handle and ``y = Afun(x)`` returns the matrix-vector product A*x.
   ``primme_eigs(A,...)`` could be replaced by primme_eigs(Afun, dim,...) in any of
   above levels of function calls. Examples are given in PRIMME_MEX_Readme.txt
   in the root directory of PRIMME_MEX folder.
 
   ``[V, D] = primme_eigs(A)`` returns a diagonal matrix D, of A's 6 largest 
   algebraic eigenvalues and a matrix V whose columns are the corresponding
   eigenvectors. 
  
   ``[V, D, norms, primmeout] = primme_eigs(A)`` also returns an array of the 
   residual norms of the computed eigenpairs, and a struct to report statistical
   information about |numOuterIterations|, |numRestarts|, |numMatvecs| and |numPreconds|.
 
   ``primme_eigs(A, numEvals)`` finds the |numEvals| largest algebraic eigenvalues.
   numEvals must be less than the dimension of the matrix A.
 
   ``primme_eigs(A, numEvals, target)`` returns numEvals target eigenvalues.
   ``target`` could be a string like below:

   * 'LA' : |primme_largest| (default)
   * 'SA' : |primme_smallest|    
   * 'CGT': |primme_closest_geq|
   * 'CLT': |primme_closest_leq|  
   * 'CT' : |primme_closest_abs|   
 
   ``primme_eigs(A, numEvals, target, opts, eigsMethod)`` specifies any of a 
   set of possible options as explained below in the opts structure. 
 
   ``eigsMethod`` is an integer specifying one of the preset methods in PRIMME:

   * 0:    |DYNAMIC|, (default)        Switches dynamically to the best method
   * 1:    |DEFAULT_MIN_TIME|,         Currently set at JDQMR_ETol
   * 2:    |DEFAULT_MIN_MATVECS|,      Currently set at GD+block
   * 3:    |Arnoldi|,                  obviously not an efficient choice 
   * 4:    |GD|,                       classical block Generalized Davidson 
   * 5:    |GD_plusK|,                 GD+k block GD with recurrence restarting
   * 6:    |GD_Olsen_plusK|,           GD+k with approximate Olsen precond.
   * 7:    |JD_Olsen_plusK|,           GD+k, exact Olsen (two precond per step)
   * 8:    |RQI|,                      Rayleigh Quotient Iteration. Also INVIT, but for INVIT provide targetShifts
   * 9:    |JDQR|,                     Original block, Jacobi Davidson
   * 10:   |JDQMR|,                    Our block JDQMR method (similar to JDCG)
   * 11:   |JDQMR_ETol|,               Slight, but efficient JDQMR modification
   * 12:   |SUBSPACE_ITERATION|,       equiv. to GD(block,2*block)
   * 13:   |LOBPCG_OrthoBasis|,        equiv. to GD(nev,3*nev)+nev
   * 14:   |LOBPCG_OrthoBasis_Window|  equiv. to GD(block,3*block)+block nev>block
 
   ``primme_eigs(A, numEvals, target, opts, eigsMethod, P)``

   ``primme_eigs(A, numEvals, target, opts, eigsMethod, P1, P2)`` uses 
   preconditioner P or P = P1*P2 to accelerate convergence of the methods.
   If P is [] then a preconditioner is not applied. P may be a function 
   handle Pfun such that Pfun(x) returns P\x.
 
   ``opts`` is an option structure which contain following parameters:         

   * |aNorm|: the estimate norm value of matrix A [{0.0}|scaler]
   * |eps|: desired computing accuracy [{1e-12}|scaler]
   * |maxBlockSize|: maximum block size the PRIMME uses [{1}|scaler]
   * |printLevel|: different level reporting(0-5) [{1}|scaler]
   * |outputFile|: output file name where user wants to save results
   * |precondition|: set to 1 if use preconditioner [{0}|1]
   * isreal: the complexity of A represented by AFUN [{ture}|false]
   * |numTargetShifts|: number of shifts for interior eigenvalues [{0}|scaler]
   * |targetShifts|: shifts for interior eigenvalues [{}|vector]
   * |initSize|: On INPUT, the number of initial guesses provided in evecs 
     array. ON OUTPUT, the number of converged eigenpairs [{0}|scaler]
   * |numOrthoConst|: Number of external orthogonalization constraints 
     provided in the first numOrthoConst vectors of evecs [{0}|scaler]
   * locking: If set to 1, hard locking will be used, otherwise the code
     will try to use soft locking [{0}|1]
   * |maxBasisSize|: maximum basis size allowed in the main iteration
   * |minRestartSize|: minimum Ritz vectors to restart
   * |maxMatvecs|: maximum number of matrix vector multiplications
     [{INT_MAX}|scaler]
   * |maxOuterIterations|: maximum number of outer iterations
     [{INT_MAX}|scaler]
   * restartingParams. |scheme|: the restart scheme [{primme_thick}| primme_dtr]
   * restartingParams. |maxPrevRetain|: number of approximations from 
     previous iteration to be retained after restart [{1}|scaler]
   * |robustShifts|: set to 1 if use robustShifting to help avoid
     stagnation and misconverge [{0}|1]
   * |maxInnerIterations|: number of inner QMR iterations [{0}|scaler]
   * |LeftQ|: a projector with Q must be applied on the left [{0}|1]
   * |LeftX|: a projector with X must be applied on the left [{0}|1]
   * |RightQ|: a projector with Q must be applied on the right [{0}|1]
   * |RightX|: a projector with X must be applied on the right [{0}|1]
   * |SkewQ|: the Q right projector must be skew [{0}|1]
   * |SkewX|: the X right projector must be skew [{0}|1]
   * |relTolBase|: a legacy from calssical JDQR (recommend not use)
   * |convTest|: how to stop the inner QMR Method
   * |iseed|: set iseed value for initialization
   * |intWorkSize|: memory size for integer workspace
   * |realWorkSize|: memory size for real or complex workspace
 
   See also :file:`Matlab/readme.txt`.

.. include:: epilog.inc
