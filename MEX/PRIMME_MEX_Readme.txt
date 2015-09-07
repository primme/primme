-----------------------------------------------------------------------------
                 PRIMME MEX: A MATLAB Interface for PRIMME
   
                Copyright (C) 2015 College of William & Mary,
   James R. McCombs, Eloy Romero Alcalde, Andreas Stathopoulos, Lingfei Wu
-----------------------------------------------------------------------------
 
   This file is part of PRIMME.

   PRIMME is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   PRIMME is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


-----------------------------------------------------------------------------
PRIMME MEX is a MATLAB interface for free software PRIMME (PReconditioned
Iterative MultiMethod Eigensolver) which finds a number of eigenvalues and 
their corresponding eigenvectors of a real symmetric, or complex hermitian
matrix A. It is a useful tool for both non-experts and experts to easily 
call PRIMME. Largest, smallest and interior eigenvalues are supported. 
Preconditioning can be used to accelerate convergence. 

-----------------------------------------------------------------------------
	Contents 
-----------------------------------------------------------------------------

1.  Directory Structure
2.  PRIMME Making & Linking
3.  PRIMME MEX Compilation
4.  Workflow in PRIMME MEX
5.  MATLAB Function primme_eigs
6.  Input and Output Descriptions of PRIMME MEX-file
7.  Examples
	
-----------------------------------------------------------------------------
1.	Directory Structure 
-----------------------------------------------------------------------------
PRIMME_MEX/

> ls 
PIRMME_mex.c          <- C language source MEX-file of PRIMME MEX
PIRMME_mex.mex???     <- PRIMME MEX interface (generated)
primme_eigs.m         <- MATLAB function for solving eigenpair problems
getMatvecHandle.m     <- perform matvec operations or get user's matvec function handle
getPrecondHandle.m    <- perform preconditioning or get user's preconditioner 
primme_eigs_example.m <- sample code using primme_eigs
PRIMME_MEX_Readme.txt <- this file

-----------------------------------------------------------------------------
2.	PRIMME Making & Linking 
-----------------------------------------------------------------------------
Users first must generate the libprimme.a library in the PRIMME root directory.
For more information about Making and Linking in PRIMME, please refer to the 
readme, Make_flags, Link_flags, makefile files in the PRIMME root directory.

PRIMME and PRIMME MEX interface both require the BLAS and LAPACK libraries.
MATLAB 2012a provides the mwlapack and mwblas libraries in 
   MATLABroot/extern/lib
These only support 64-bit integers and pointers (ILP64) so to work with PRIMME
on LP64 systems, the compilation flags in Make_flags in the top directory of 
PRIMME distribution must include the following definition:

   CFLAGS += -DPRIMME_BLASINT_SIZE=64

-----------------------------------------------------------------------------
3.	PRIMME MEX Compilation 
-----------------------------------------------------------------------------
For general information about building a MATLAB MEX file, please refer to 
www.mathworks.com/help/MATLAB/MATLAB_external/building-mex-files.html.
There are two steps to build a MATLAB MEX file.
 
1. confirm that your compiler is supported by your current MATLAB version.
   http://www.mathworks.com/support/compilers/R2015a/
   http://www.mathworks.com/support/compilers/R2014a/
   http://www.mathworks.com/support/compilers/R2012a/
   http://www.mathworks.com/support/compilers/R2010a/
   We have tested PRIMME MEX on SUSE Linux with R2012a, R2010a and R2015a;
   and on MAC OSX 10.10 with R2015a. 

2. build the PRIMME_mex.c using the makefile in the root directory of 
   PRIMME MEX. Alternatively, if $TOP is the path to the PRIMME installation, 
   in the MATLAB command prompt type:

   mex -v -O -largeArrayDims PRIMME_mex.c -L$TOP/PRIMME \
	-I$TOP/PRIMME/PRIMMESRC/COMMONSRC -lprimme -lm  -lmwlapack -lmwblas 


-----------------------------------------------------------------------------
4.	Workflow in PRIMME MEX 
-----------------------------------------------------------------------------
Although the generated PRIMME_mex binary could be called directly as an M-file
in MATLAB, users typically call the primme_eigs M-file. This presents a user 
interface that extends the interface of the built-in eigs function, performs 
initializations, sets up the MATLAB matvec/preconditioning functions from 
the user input, and calls PRIMME_mex. 

The PRIMME_mex MEX function converts the MATLAB inputs into C inputs as
expected by the PRIMME library and calls PRIMME. Moreover, it sets up
the matvec/preconditioning C functions to be passed to PRIMME. When these 
functions are called, they propagate the information back to their MATLAB
counterparts so that the operations are performed in MATLAB.
After PRIMME returns, PRIMME_mex returns the C output to MATLAB primme_eigs,
which returns it to the user.

-----------------------------------------------------------------------------
5.   MATLAB Function primme_eigs
-----------------------------------------------------------------------------  
Typically users will call the primme_eigs MATLAB function. 
The matrix A must be real symmetric or complex Hermitian.

Input: [A, numEvals, target, opts, eigsMethod, P]

Output: [evals, evecs, norms, primmeout]

Accepted function calls:  
primme_eigs(A)
primme_eigs(A, numEvals)
primme_eigs(A, numEvals, target)
primme_eigs(A, numEvals, target, opts)
primme_eigs(A, numEvals, target, opts, eigsMethod)
primme_eigs(A, numEvals, target, opts, eigsMethod, P)
primme_eigs(A, numEvals, target, opts, eigsMethod, P1,P2)
primme_eigs(A, numEvals, target, opts, eigsMethod, Pfun)

primme_eigs(Afun, dim,...)
Instead of a matrix, the user can provide a function handle Afun()
where Afun(x) returns the matrix-vector product y=A*x, where x,y are vectors
of length dim. See examples Primme_eigsTest6, 7, 8, and 9 later in this file.

primme_D = primme_eigs(A) returns a vector of the 6 largest algebraic
eigenvalues of A for the parameter defaults shown below.

[primme_V, primme_D] = primme_eigs(A) returns a diagonal matrix primme_D of 
the 6 largest algebraic eigenvalues and a matrix primme_V whose columns 
are the corresponding eigenvectors. 

[primme_V, primme_D, norms, primmeout] = primme_eigs(A, numEvals)
returns in addition a double array of the residual norms of the eigenpairs,
and a struct to report statistical information about numOuterIterations,
numRestarts, numMatvecs and numPreconds. 

numEvals is the number of eigenvalues that the user wants to compute.
It must be less than the dimension of the matrix A.

primme_eigs(A, numEvals, target) returns numEvals target eigenvalues.
target should be one of the strings below. See PRIMME's readme for details.
 'LA' ------ primme_largest (default)
 'SA' ------ primme_smallest    
 'CGT'------ primme_closest_geq  
 'CLT'------ primme_closest_leq  
 'CT' ------ primme_closest_abs   

primme_eigs(A, numEvals, target, opts) specifies various options
for PRIMME parameters (see below).

primme_eigs(A, numEvals, target, opts, eigsMethod) specifies in addition 
the eigensolver method in PRIMME. The integer eigsMethod should take one 
of the values:

eigsMethod	corresponding PRIMME method
  0:	   DYNAMIC, (default)        Switches dynamically to the best method
  1: 	   DEFAULT_MIN_TIME,         Currently set at JDQMR_ETol
  2:	   DEFAULT_MIN_MATVECS,      Currently set at GD+block
  3:	   Arnoldi,                  obviously not an efficient choice 
  4:	   GD,                       classical block Generalized Davidson 
  5:	   GD_plusK,                 GD+k block GD with recurrence restarting
  6:	   GD_Olsen_plusK,           GD+k with approximate Olsen precond.
  7:	   JD_Olsen_plusK,           GD+k, exact Olsen (two precond per step)
  8:	   RQI,                      Rayleigh Quotient Iteration. Also INVIT,
                                     but for INVIT provide targetShifts
  9:	   JDQR,                     Original block, Jacobi Davidson
  10:	   JDQMR,                    Our block JDQMR method (similar to JDCG)
  11:	   JDQMR_ETol,               Slight, but efficient JDQMR modification
  12:	   SUBSPACE_ITERATION,       equiv. to GD(block,2*block)
  13:	   LOBPCG_OrthoBasis,        equiv. to GD(nev,3*nev)+nev
  14:	   LOBPCG_OrthoBasis_Window  equiv. to GD(block,3*block)+block nev>block

primme_eigs(A, numEvals, target, opts, eigsMethod, P) uses preconditioner P
or factorized P = P1*P2. If P is [] or omitted, no preconditioner is applied. 
P may be a function handle Pfun such that Pfun(x) returns P\x.

The opts contains 31 different parameters which can be used to fine tune
the execution of PRIMME as described in PRIMME's readme file.
The most common options are:

opts.aNorm: the estimate norm value of matrix A [{0.0}|scalar]
opts.eps: desired computing accuracy [{1e-12}|scalar]
opts.maxBlockSize: maximum block size in PRIMME [{1}|integer]
opts.printLevel: different level reporting(0-5) [{1}|integer]
opts.outputFile: file name where the user wants to save the results
opts.precondition: set to 1 if preconditioner is to be used [{0}|1]
opts.isreal: Whether afun represents a real/complex matrix [{true}|false]
opts.numTargetShifts: number of shifts for interior eigenvalues [{0}|integer]
opts.targetShifts: shifts for interior eigenvalues [{}|scalar vector]
 
-----------------------------------------------------------------------------
6.	Input and Output Descriptions of PRIMME MEX-file
-----------------------------------------------------------------------------
If PRIMME_mex(flag,dim,...) is called directly, it returns a few eigenvalues
and eigenvectors of the current Hermitian matrix A. Its interface is closer to
PRIMME functions d/zprimme.

 Syntax:
 
 [evals, evecs, norms, primmeout] = PRIMME_mex(flag, dim, Numeigs, target, method, opts)
  
  evals = PRIMME_mex(flag, dim) 
  evals = PRIMME_mex(flag, dim, numEvals)
  evals = RRIMME_mex(flag, dim, numEvals, target)
  evals = PRIMME_mex(flag, dim, numEvals, target, method)
  evals = PRIMME_mex(flag, dim, numEvals, target, method, opts)
 [evals,evecs] = PRIMME_mex(flag, dim, numEvals, target, method, opts)
 [evals, evecs, norms] = PRIMME_mex(flag, dim, numEvals, target, method, opts)
 [evals, evecs, norms, primmeout] = PRIMME_mex(flag, dim, numEvals, target, method, opts)

 Output desciption: [evals, evecs, norms, primmeout]

  evals: a double array of the target eigenvalues returned to MATLAB

  evecs: a double array of the corresponding eigenvectors returned to MATLAB

  norms: a double array of the residual norms of the eigenvalues

  primmeout: a struct to report statistics which contains four members:
     numOuterIterations
     numRestarts
     numMatvecs
     numPreconds

 Input description: (flag, dim, numEvals, target, method, opts) 

  flag: mark the input matrix is real or complex
 
  dim: the dimension of large and sparse symmetric and hermitian matrix.
 
  numEvals: number of eigenvalues required

  target: Which eigenpairs to find.  target can be any of the following enum:
 
    primme_smallest --- 'SA':Smallest algebraic eigenvalues. Target shifts ignored
    primme_largest  --- 'LA':Largest  algebraic eigenvalues. Target shifts ignored
    primme_closest_geq -'CGT':Closest to, but greater or equal than a set of shifts
    primme_closest_leq -'CLT':Closest to, but less or equal than a set of shifts
    primme_closest_abs - 'CT':Closest in absolute value to a set of shifts

  method: which solver method to choose. method is an enum type below:
 
  method	corresponding PRIMME method
    0:	   DYNAMIC, (default)        Switches dynamically to the best method
    1: 	   DEFAULT_MIN_TIME,         Currently set at JDQMR_ETol
    2:	   DEFAULT_MIN_MATVECS,      Currently set at GD+block
    3:	   Arnoldi,                  obviously not an efficient choice 
    4:	   GD,                       classical block Generalized Davidson 
    5:	   GD_plusK,                 GD+k block GD with recurrence restarting
    6:	   GD_Olsen_plusK,           GD+k with approximate Olsen precond.
    7:	   JD_Olsen_plusK,           GD+k, exact Olsen (two precond per step)
    8:	   RQI,                      Rayleigh Quotient Iteration. Also INVIT,
                                       but for INVIT provide targetShifts
    9:	   JDQR,                     Original block, Jacobi Davidson
    10:	   JDQMR,                    Our block JDQMR method (similar to JDCG)
    11:	   JDQMR_ETol,               Slight, but efficient JDQMR modification
    12:	   SUBSPACE_ITERATION,       equiv. to GD(block,2*block)
    13:	   LOBPCG_OrthoBasis,        equiv. to GD(nev,3*nev)+nev
    14:	   LOBPCG_OrthoBasis_Window  equiv. to GD(block,3*block)+block nev>block

  'opts' is an option structure which contain following parameters of the 
     primme_params structure: 

    opts.aNorm: the estimate norm value of matrix A [default = 0]
    opts.eps: desired computing accuracy [default = 1e-12]
    opts.numTargetShifts: shifts for interior eigenvalues [default = 0]
    opts.targetShifts: pointer to get each shift for interior eigenvalues
    opts.initSize: initial guesses/constraints [default = 0]
    opts.numOrthoConst: [default = 0]
    opts.locking: 0 or 1 
    opts.dynamicMethodSwitch: from -3 to 1
    opts.maxBasisSize: maximum basis size allowed in the main iteration
    opts.minRestartSize: minimum Ritz vectors to restart
    opts.maxBlockSize:  [default = 1]
    opts.maxMatvecs: [default = INT_MAX]
    opts.maxOuterIterations: [default = INT_MAX]
    opts.restartingParams.scheme: [default = primme_thick]
    opts.restartingParams.maxPrevRetain: [default = 1]
    opts.precondition: set to 1 if preconditioning is to be performed,
        make sure the applyPreconditioner is not NULL [default = 0]
    opts.robustShifts: set to 1 to use robustShifting
    opts.maxInnerIterations: 0 (no inner iterations), = k (perform 
        at most k inner iterations per outer step)
    opts.LeftQ: 0 or 1
    opts.LeftX: 0 or 1
    opts.RightQ: 0 or 1
    opts.RightX: 0 or 1
    opts.SkewQ: 0 or 1
    opts.SkewX: 0 or 1
    opts.relTolBase: a legacy from calssical JDQR
    opts.convTest: how to stop the inner QMR method
    opts.printLevel: 0-5 (different level reporting) [default = 1]
    opts.outputFile: output file name where user wants to save results
    opts.iseed: set iseed value for initialization
    opts.intWorkSize: memory size for int workspace
    opts.realWorkSize: memory size for real workspace        

-----------------------------------------------------------------------------
7.	Examples on how to call primme_eigs()
-----------------------------------------------------------------------------

Example 1: Primme_eigsTest.m
%
% Simple test case for dprimme
% Function call: primme_eigs(A, numEvals, target)
%
       clear all
       A = delsq(numgrid('C',300));
       
       numEvals = 5;
       target = 'LA';
       
       primme_start = tic;
       [primme_V,primme_D,norms,primmeout] = primme_eigs(A, numEvals, target);
       primme_telapsed = toc(primme_start)
       primme_evalues = primme_D
       primme_numMatvec = primmeout(3)
       
       global Matvec_counter;
       Matvec_counter = 0;
       n = size(A,1);
       opts.tol = 1e-12;
       opts.issym = 1;
       opts.isreal = 1;
       opts.maxit = 65536;
       opts.v0 = (1/sqrt(n))*ones(n,1);
       eigs_start = tic;
       [V ,D] = eigs(A , numEvals, 'SA', opts);
       %[V ,D] = eigs(@(x)test_eigs(x), n , numEvals, 'LA', opts);
       eigs_evalues = D
       eigs_telapsed = toc(eigs_start)
       eigs_numMatvec = Matvec_counter
       
Results:

primme_telapsed =

   14.1366


primme_evalues =

    7.9997         0         0         0         0
         0    7.9994         0         0         0
         0         0    7.9992         0         0
         0         0         0    7.9989         0
         0         0         0         0    7.9987


primme_numMatvec =

        4541


eigs_evalues =

    7.9997         0         0         0         0
         0    7.9994         0         0         0
         0         0    7.9992         0         0
         0         0         0    7.9989         0
         0         0         0         0    7.9987


eigs_telapsed =

   33.1543


eigs_numMatvec =

        6000


Example 2: Primme_eigsTest2.m
%
% Senior test case for dprimme
% Function call: primme_eigs(A, numEvals, target, opts, eigsMethod)
%
       clear all
       load cfd1.mat;
       
       if (exist('A','var') == 0)
          A = Problem.A;
       end

       numEvals = 5;
       target = 'SA';
       eigsMethod = 1;
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'printLevel', 2, 'outputFileName', {'sampleout'}); 

       primme_start = tic;
       [primme_V,primme_D,norms,primmeout]= primme_eigs(A, numEvals, target, opts, eigsMethod);
       primme_telapsed = toc(primme_start)
       primme_evalues = primme_D
       primme_numMatvec = primmeout(3)
 
        
       global Matvec_counter;
       Matvec_counter = 0;
       n = size(A,1);
       opts.tol = 1e-10;
       opts.issym = 1;
       opts.isreal = 1;
       opts.maxit = 65536;
       opts.v0 = (1/sqrt(n))*ones(n,1);
       eigs_start = tic;
       [V,D] = eigs(@(x)test_eigs(x), n , numEvals, 'SA', opts);
       eigs_evalues = D
       eigs_telapsed = toc(eigs_start)
       eigs_numMatvec = Matvec_counter
       
Results:

primme_telapsed =

   36.1040


primme_evalues =

   1.0e-03 *

    0.0200         0         0         0         0
         0    0.0886         0         0         0
         0         0    0.1357         0         0
         0         0         0    0.1815         0
         0         0         0         0    0.2158


primme_numMatvec =

        6201


eigs_evalues =

   1.0e-03 *

    0.0200         0         0         0         0
         0    0.0886         0         0         0
         0         0    0.1357         0         0
         0         0         0    0.1815         0
         0         0         0         0    0.2158


eigs_telapsed =

  235.2509


eigs_numMatvec =

       28569


Example 3: Primme_eigsTest3.m

%
% Simple test case for zprimme. The input is complex hermitian matrix.
% Function call: primme_eigs(A, numEvals, target)
%
       clear all    
       load HA.mat;

       numEvals = 5;
       target = 'LA';
       
       primme_start = tic;
       [primme_V,primme_D,norms,primmeout]= primme_eigs(A, numEvals, target);
       primme_telapsed = toc(primme_start)
       primme_evalues = primme_D
       primme_numMatvec = primmeout(3)
       
       global Matvec_counter;
       Matvec_counter = 0;
       n = size(A,1);
       opts.tol = 1e-12;
       opts.issym = 1;
       opts.isreal = 0;
       opts.maxit = 65536;
       vecr = (1/sqrt(n))*ones(n,1);
       veci = zeros(n,1);
       opts.v0 = complex(vecr,veci);
       eigs_start = tic;
       [V,D] = eigs(@(x)test_eigsZ(x), n , numEvals, 'Lr', opts);
       eigs_evalues = D
       eigs_telapsed = toc(eigs_start)
       eigs_numMatvec = Matvec_counter
       
results:

primme_telapsed =

    0.5496


primme_evalues =

    3.3077         0         0         0         0
         0    3.3076         0         0         0
         0         0    3.3073         0         0
         0         0         0    3.3073         0
         0         0         0         0    3.3072


primme_numMatvec =

        1715


eigs_evalues =

    3.3077         0         0         0         0
         0    3.3076         0         0         0
         0         0    3.3072         0         0
         0         0         0    3.3073         0
         0         0         0         0    3.3072


eigs_telapsed =

    0.3097


eigs_numMatvec =

   343


Example 4: Primme_eigsTest4.m
%
% Senior test case for zprimme. The input is complex hermitian matrix.
% Function call: primme_eigs(A, numEvals, target, opts, eigsMethod)
%
       clear all
       load 3Dspectralwave2.mat;
       
       if (exist('A','var') == 0)
          A = Problem.A;
       end
       
       numEvals = 5;
       target = 'SA';
       eigsMethod = 2;
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'printLevel', 2, 'outputFileName', {'sampleout'}); 
  
       primme_start = tic;
       [primme_V,primme_D,norms,primmeout]= primme_eigs(A, numEvals, target, opts, eigsMethod);
       primme_telapsed = toc(primme_start)
       primme_evalues = primme_D
       Primme_numMatvec = primmeout(3)
       
       global Matvec_counter;
       Matvec_counter = 0;
       n = size(A,1);
       opts.tol = 1e-10;
       opts.issym = 1;
       opts.isreal = 0;
       opts.maxit = 65536;
       vecr = (1/sqrt(n))*ones(n,1);
       veci = zeros(n,1);
       opts.v0 = complex(vecr,veci);
       eigs_start = tic;
       [V,D] = eigs(@(x)test_eigsZ(x), n , numEvals, 'Sr', opts);
       eigs_evalues = D
       eigs_telapsed = toc(eigs_start)
       eigs_numMatvec = Matvec_counter

Results:

primme_telapsed =

   37.8866


primme_evalues =

  -42.7384         0         0         0         0
         0  -42.0612         0         0         0
         0         0  -41.4534         0         0
         0         0         0  -40.9778         0
         0         0         0         0  -40.9267


Primme_numMatvec =

   367


eigs_evalues =

  -42.7384         0         0         0         0
         0  -42.0612         0         0         0
         0         0  -41.4534         0         0
         0         0         0  -40.9267         0
         0         0         0         0  -40.9778


eigs_telapsed =

   32.0614


eigs_numMatvec =

   259


Example 5: Primme_eigsTest5.m
%
% Senior test case for user to provide matrix function instead of the 
% matrix A to find a few of eigenvalues and eigenvectors. The eigsMatvec  
% is a matrix function handle and Y = Afun(x) should return the result 
% y = A*x. 
% Function call: primme_eigs(Afun, dim, numEvals, target, opts, eigsMethod)
%
       clear all
       k = 1000;
       on = ones(k,1); 

       numEvals = 5;
       target = 'SA';
       eigsMethod = 1;
       dim = k;
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'precondition', 0,'printLevel', 2, 'outputFileName', {'sampleout'}); 

       primme_start = tic;
       [primme_V,primme_D,norms,primmeout]= primme_eigs(@(x)eigsMatvec(x,on,k), dim, numEvals, target, opts, eigsMethod);
       primme_telapsed = toc(primme_start)
       primme_evalues = primme_D
       primme_numMatvec = primmeout(3)
 
        
       global Matvec_counter;
       Matvec_counter = 0;
       n = k;
       opts.tol = 1e-10;
       opts.issym = 1;
       opts.isreal = 1;
       opts.maxit = 65536;
       opts.v0 = (1/sqrt(n))*ones(n,1);
       eigs_start = tic;
       [V,D] = eigs(@(x)test_eigs(x), n , numEvals, 'SA', opts);
       eigs_evalues = D
       eigs_telapsed = toc(eigs_start)
       eigs_numMatvec = Matvec_counter

Results:

primme_telapsed =

    1.9573


primme_evalues =

   1.0e-03 *

    0.0197         0         0         0         0
         0    0.0788         0         0         0
         0         0    0.1773         0         0
         0         0         0    0.3152         0
         0         0         0         0    0.4925


primme_numMatvec =

        3509


eigs_evalues =

   1.0e-03 *

    0.0197         0         0         0         0
         0    0.0788         0         0         0
         0         0    0.1773         0         0
         0         0         0    0.3152         0
         0         0         0         0    0.4925


eigs_telapsed =

    4.9927


eigs_numMatvec =

        6675


Example 6: Primme_eigsTest6.m
%
% Senior test case for user to provide their preconditioner to speed up the
% whole process. The input A is a matrix and preconditioner P is also a
% matrix.
% Function call: primme_eigs(A, numEvals, target, opts, eigsMethod, P)
%
       clear all
       k = 1000;
       on = ones(k,1); 
       A = spdiags([-2*on 4*on -2*on],-1:1,k,k);
       

       numEvals = 5;
       target = 'SA';
       eigsMethod = 1;
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'precondition', 1,'printLevel', 2, 'outputFileName', {'sampleout'}); 
       
       P = spdiags(4*on,0,k,k);

       primme_start = tic;
       [primme_V,primme_D,norms,primmeout]= primme_eigs(A, numEvals, target, opts, eigsMethod, P);
       primme_telapsed = toc(primme_start)
       primme_evalues = primme_D
       primme_numMatvec = primmeout(3)
 
        
       global Matvec_counter;
       Matvec_counter = 0;
       n = size(A,1);
       opts.tol = 1e-10;
       opts.issym = 1;
       opts.isreal = 1;
       opts.maxit = 65536;
       opts.v0 = (1/sqrt(n))*ones(n,1);
       eigs_start = tic;
       [V,D] = eigs(@(x)test_eigs(x), n , numEvals, 'SA', opts);
       eigs_evalues = D
       eigs_telapsed = toc(eigs_start)
       eigs_numMatvec = Matvec_counter

Results:

primme_telapsed =

    0.5341


primme_evalues =

   1.0e-03 *

    0.0197         0         0         0         0
         0    0.0788         0         0         0
         0         0    0.1773         0         0
         0         0         0    0.3152         0
         0         0         0         0    0.4925


primme_numMatvec =

        3509


eigs_evalues =

   1.0e-03 *

    0.0197         0         0         0         0
         0    0.0788         0         0         0
         0         0    0.1773         0         0
         0         0         0    0.3152         0
         0         0         0         0    0.4925


eigs_telapsed =

    1.9208


eigs_numMatvec =

        6675

Example 7: Primme_eigsTest7.m
%
% Senior test case for user to provide their preconditioner to speed up the
% whole process. The input A is a matrix and preconditioner P1 and P2 are
% both preconditioners.
% Function call: primme_eigs(A, numEvals, target, opts, eigsMethod, P1,P2)
%
       clear all
       k = 1000;
       on = ones(k,1); 
       A = spdiags([-2*on 4*on -2*on],-1:1,k,k);
       

       numEvals = 5;
       target = 'SA';
       eigsMethod = 2;
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'precondition', 1,'printLevel', 2, 'outputFileName', {'sampleout'}); 
       
       P1 = spdiags([on/(-2) on],-1:0,k,k); 
       P2 = spdiags([4*on -on],0:1,k,k);

       primme_start = tic;
       [primme_V,primme_D,norms,primmeout]= primme_eigs(A, numEvals, target, opts, eigsMethod, P1, P2);
       primme_telapsed = toc(primme_start)
       primme_evalues = primme_D
       primme_numMatvec = primmeout(3)
 
        
       global Matvec_counter;
       Matvec_counter = 0;
       n = size(A,1);
       opts.tol = 1e-10;
       opts.issym = 1;
       opts.isreal = 1;
       opts.maxit = 65536;
       opts.v0 = (1/sqrt(n))*ones(n,1);
       eigs_start = tic;
       [V,D] = eigs(@(x)test_eigs(x), n , numEvals, 'SA', opts);
       eigs_evalues = D
       eigs_telapsed = toc(eigs_start)
       eigs_numMatvec = Matvec_counter

Results:

primme_telapsed =

    0.4476


primme_evalues =

   1.0e-03 *

    0.0197         0         0         0         0
         0    0.0788         0         0         0
         0         0    0.1773         0         0
         0         0         0    0.3152         0
         0         0         0         0    0.4925


primme_numMatvec =

        1801


eigs_evalues =

   1.0e-03 *

    0.0197         0         0         0         0
         0    0.0788         0         0         0
         0         0    0.1773         0         0
         0         0         0    0.3152         0
         0         0         0         0    0.4925


eigs_telapsed =

    1.9122


eigs_numMatvec =

        6675

Example 8: Primme_eigsTest8.m
%
% Senior test case for user to provide function handle instread of the 
% matrix A to find a few of eigenvalues and eigenvectors. The eigsMatvec  
% is a matrix function handle and Y = Afun(x) should return the result 
% y = A*x. Also, the matrix P is provided as preconditioner.
% Function call: primme_eigs(Afun, dim, numEvals, target, opts, eigsMethod, P)
%
       clear all
       k = 1000;
       on = ones(k,1); 
    
       numEvals = 5;
       target = 'SA';
       eigsMethod = 1;
       dim = k;
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'precondition', 1,'printLevel', 2, 'outputFileName', {'sampleout'}); 
       
       P = spdiags(4*on,0,k,k);

       primme_start = tic;
       [primme_V,primme_D,norms,primmeout]= primme_eigs(@(x)eigsMatvec(x,on,k), dim, numEvals, target, opts, eigsMethod, P);
       primme_telapsed = toc(primme_start)
       primme_evalues = primme_D
       primme_numMatvec = primmeout(3)
 
        
       global Matvec_counter;
       Matvec_counter = 0;
       n = k;
       opts.tol = 1e-10;
       opts.issym = 1;
       opts.isreal = 1;
       opts.maxit = 65536;
       opts.v0 = (1/sqrt(n))*ones(n,1);
       eigs_start = tic;
       [V, D] = eigs(@(x)test_eigs(x), n , numEvals, 'SA', opts);
       eigs_evalues = D
       eigs_telapsed = toc(eigs_start)
       eigs_numMatvec = Matvec_counter
       
Results:

primme_telapsed =

    2.2001


primme_evalues =

   1.0e-03 *

    0.0197         0         0         0         0
         0    0.0788         0         0         0
         0         0    0.1773         0         0
         0         0         0    0.3152         0
         0         0         0         0    0.4925


primme_numMatvec =

        3509


eigs_evalues =

   1.0e-03 *

    0.0197         0         0         0         0
         0    0.0788         0         0         0
         0         0    0.1773         0         0
         0         0         0    0.3152         0
         0         0         0         0    0.4925


eigs_telapsed =

    4.9958


eigs_numMatvec =

        6675


Example 9: Primme_eigsTest9.m
%
% Senior test case for user to provide function handle instead of the 
% matrix A to find a few of eigenvalues and eigenvectors. The eigsMatvec  
% is a matrix function handle and Y = Afun(x) should return the result 
% y = A*x. Also, the eigsPrecond is a function handle as preconditioner. Y
% = Pfun(X) should return the result y = P\x.
% Function call: primme_eigs(Afun, dim, numEvals, target, opts, eigsMethod, Pfun)
%
       clear all
       k = 1000;
       on = ones(k,1); 
    
       numEvals = 5;
       target = 'SA';
       eigsMethod = 2;
       dim = k;
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'precondition', 1,'printLevel', 2, 'outputFileName', {'sampleout'});

       primme_start = tic;
       [primme_V,primme_D,norms,primmeout]= primme_eigs(@(x)eigsMatvec(x,on,k), dim, numEvals, target, opts, eigsMethod, @(x)eigsPrecond(x,on,k));
       primme_telapsed = toc(primme_start)
       primme_evalues = primme_D
       primme_numMatvec = primmeout(3)
 
        
       global Matvec_counter;
       Matvec_counter = 0;
       n = k;
       opts.tol = 1e-10;
       opts.issym = 1;
       opts.isreal = 1;
       opts.maxit = 65536;
       opts.v0 = (1/sqrt(n))*ones(n,1);
       eigs_start = tic;
       [V,D] = eigs(@(x)test_eigs(x), n , numEvals, 'SA', opts);
       eigs_evalues = D
       eigs_telapsed = toc(eigs_start)
       eigs_numMatvec = Matvec_counter

Results:

primme_telapsed =

    2.1459


primme_evalues =

   1.0e-03 *

    0.0197         0         0         0         0
         0    0.0788         0         0         0
         0         0    0.1773         0         0
         0         0         0    0.3152         0
         0         0         0         0    0.4925


primme_numMatvec =

        2487


eigs_evalues =

   1.0e-03 *

    0.0197         0         0         0         0
         0    0.0788         0         0         0
         0         0    0.1773         0         0
         0         0         0    0.3152         0
         0         0         0         0    0.4925


eigs_telapsed =

    4.9730


eigs_numMatvec =

        6675
