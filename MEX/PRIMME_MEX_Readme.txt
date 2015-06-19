
-----------------------------------------------------------------------------
        PRIMME MEX: A MATLAB Interface for PRIMME
   
      Copyright (C) 2012  Lingfei Wu,  Andreas Stathopoulos
-----------------------------------------------------------------------------
PRIMME MEX is a MATLAB interface for free software PRIMME (PReconditioned
Iterative MultiMethod Eigensolver) which finds a number of eigenvalues and 
their corresponding eigenvectors of a real symmetric, or complex hermitian
matrix A. It is a useful tool for both non-experts and experts to easily 
call PRIMME. Largest, smallest and interior eigenvalues are supported. 
Preconditioning can be used to accelerate convergence. PRIMME is written 
in C, but a complete Fortran77 interface is provided.

-----------------------------------------------------------------------------
	Contents 
-----------------------------------------------------------------------------

1.	Directory Structure
2. 	PRIMME Changes for Supporting MATLAB MEX Interface
3.	PRIMME Making & Linking
4.	PRIMME MEX Compilation
5.	Workflow in PRIMME MEX
6.  MATLAB Function Call for solving eigenpair problems
7.  MATLAB Function Call for solving singular value problems
8.  Input and Output Descriptions of PRIMME MEX-file
9.	Examples
	
-----------------------------------------------------------------------------
1.	Directory Structure 
-----------------------------------------------------------------------------
PRIMME_MEX/

> ls 
PIRMME_mex.c          <- C language source MEX-file of PRIMME MEX
PIRMME_mex.mexa64     <- executable mex file of PRIMME MEX in linux SUSE
primme_eigs.m         <- MATLAB function for solving eigenpair problems
primme_svds.m         <- MATLAB function for solving singluar value problems
getMatvecHandle.m     <- perform matvec operations or get user's matvec function handle
getPrecondHandle.m    <- perform preconditioning or get user's preconditioner 
eigsMatvec.m	      <- test case of user's matvec function for primme_eigs 
eigsPrecond.m         <- test case of user's preconditioner for primme_eigs
svdsMatvec.m	      <- test case of user's matvec function for primme_svds 
svdsPrecond.m         <- test case of user's preconditioner for primme_svds
Primme_eigsTest1-9.m  <- total 9 test cases for primme_eigs
Primme_svdsTest1-9.m  <- total 9 test cases for primme_svds
PRIMME_MEX_Readme.txt <- this file

-----------------------------------------------------------------------------
2.	PRIMME Changes for Supporting MATLAB MEX Interface 
-----------------------------------------------------------------------------
PRIMME and PRIMME MEX interface both utilize BLAS and LAPACK libraries. MATLAB
2012a provides the mwlapack and mwblas libraries in MATLABroot/extern/lib, which 
only support 64-bit integers for matrix dimensions. In order to call mwlapck
and mwblas libraries in the PRIMME, several functions in the common_numerical.c,
numerical_d.c, and numerical_z.c need to be changed. In these files, int type 
variables or pointer variables have changed to long long int type ones. The 
related functions are as follows:

   void Num_dcopy_primme(int n, double *x, int incx, double *y, int incy);   
   int Num_dspev_dprimme(int iopt, double *ap, double *w, double *z, int ldz, 
   int n, double *aux, int naux);
   void Num_dsyev_dprimme(char *jobz, char *uplo, int n, double *a, int lda, 
   double *w, double *work, int ldwork, int *info);
   void Num_dsytrf_dprimme(char *uplo, int n, double *a, int lda, int *ipivot, 
   double *work, int ldwork, int *info);
   void Num_dsytrs_dprimme(char *uplo, int n, int nrhs, double *a, int lda, 
   int *ipivot, double *b, int ldb, int *info);
   void Num_dcopy_dprimme(int n, double *x, int incx, double *y, int incy);
   double Num_dot_dprimme(int n, double *x, int incx, double *y, int incy);
   void Num_gemm_dprimme(char *transa, char *transb, int m, int n, int k, 
   double alpha, double *a, int lda, double *b, int ldb, 
   double beta, double *c, int ldc);
   void Num_symm_dprimme(char *side, char *uplo, int m, int n, double alpha, 
   double *a, int lda, double *b, int ldb, double beta, 
   double *c, int ldc);
   void Num_axpy_dprimme(int n, double alpha, double *x, int incx, 
   double *y, int incy);
   void Num_gemv_dprimme(char *transa, int m, int n, double alpha, double *a,
   int lda, double *x, int incx, double beta, double *y, int incy);
   void Num_larnv_dprimme(int idist, int *iseed, int length, double *x);
   void Num_scal_dprimme(int n, double alpha, double *x, int incx);
   void Num_swap_dprimme(int n, double *x, int incx, double *y, int incy);
   
   int Num_zhpev_zprimme(int iopt, Complex_Z *ap, double *w, Complex_Z *z, int ldz, 
   int n, Complex_Z *aux, double *rwork, int naux);
   void Num_zheev_zprimme(char *jobz, char *uplo, int n, Complex_Z *a, int lda, 
   double *w, Complex_Z *work, int ldwork, double *rwork, int *info);
   void Num_zhetrf_zprimme(char *uplo, int n, Complex_Z *a, int lda, int *ipivot,
   Complex_Z *work, int ldwork, int *info);
   void Num_zhetrs_zprimme(char *uplo, int n, int nrhs, Complex_Z *a, int lda, 
   int *ipivot, Complex_Z *b, int ldb, int *info);
   void Num_zcopy_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy);
   Complex_Z Num_dot_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy);
   void Num_gemm_zprimme(char *transa, char *transb, int m, int n, int k, 
   Complex_Z alpha, Complex_Z *a, int lda, Complex_Z *b, int ldb, 
   Complex_Z beta, Complex_Z *c, int ldc);
   void Num_symm_zprimme(char *side, char *uplo, int m, int n, Complex_Z alpha, 
   Complex_Z *a, int lda, Complex_Z *b, int ldb, Complex_Z beta, 
   Complex_Z *c, int ldc);
   void Num_axpy_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx, 
   Complex_Z *y, int incy);
   void Num_gemv_zprimme(char *transa, int m, int n, Complex_Z alpha, Complex_Z *a,
   int lda, Complex_Z *x, int incx, Complex_Z beta, Complex_Z *y, int incy);
   void Num_larnv_zprimme(int idist, int *iseed, int length, Complex_Z *x);
   void Num_scal_zprimme(int n, Complex_Z alpha, Complex_Z *x, int incx);
   void Num_swap_zprimme(int n, Complex_Z *x, int incx, Complex_Z *y, int incy);

In all functions listed above, "int" is replaced by "long long int" so that 
PRIMME MEX can work correctly. Otherwise, it is suffered by different kinds
of memory violation reports when running PRIMME MEX. 
   
In the future version of PRIMME, the more elegant solution may be defining a 
generic int type such as primme_integer type which represents 4 Byte int type 
or 8 Byte int type depending on different systems and MATLAB versions.

-----------------------------------------------------------------------------
3.	PRIMME Making & Linking 
-----------------------------------------------------------------------------
Users must customize Make_flags to create the library. 
Users may customize Link_flags to create the test programs.

To support PRIMME MEX interface, the first thing is to change all single 
line comments with "//..." by using "/*...*/" in primme.h file since MATLAB
MEX-file compilation only supports this comment method. For customizing the 
Make_flags and Link_flags files, please refer to readme file in the PRIMME.
Overall, you may change several important things in the Make_flags file and
the Link_flags file depending on your systems. 

Firstly, modify the path of the installed PRIMME directory. Secondly, for 
Make_flags file, you need change "F77 = f77" to "F77 = gfortran" if you 
OS is Suse Linux. You also need change "CFLAGS =" to "CFLAGS = -fPIC". 
In addtion, "FFLAGS = -O -fno-second-underscore" is replaced to "FFLAGS 
= -O -fno-second-underscore -fPIC". For Link_flags file, change "LDR = 
f77" to "LDR = gfortran" if your OS is Suse Linux.

makefile can perform the following functions:

make all  		builds: lib depends seqs pars
 make lib 		builds libprimme.a in PRIMME/. Alternatively:
 make libd 		  if only dprimme is of interest build libdprimme.a
 make libz 		  if only zprimme is of interest build libzprimme.a
 make depends		builds the dependencies files in (D)ZTEST/
   make ddepends_seq       builds only the ddependencies_seq in DTEST
   make zdepends_seq       builds only the zdependencies_seq in ZTEST
   make ddepends_par       builds only the ddependencies_par in DTEST
 make seqs		builds all sequential executables in (D)ZTEST
   make seq_dprimme 	   builds only the seq real C executable in DTEST
   make seqf77_dprimme 	   builds only the seq real F77 executable in DTEST
   make seq_zprimme 	   builds only the seq herm C executable in ZTEST
   make seqf77_zprimme	   builds only the seq herm F77 executable in ZTEST
 make pars		builds all the parallel executables in DTEST
   make par_dprimme	   currently the only parallel C executable in DTEST
 make clean             removes all *.o, a.out, and core files from all dirs
 make backup		makes a tar.gz dated archive of entire PRIMME directory

To build PRIMME MEX interface, you need generate libprimme.a library in the 
PRIMME root directory. For more information about Making and Linking in PRIMME, 
please refer to readme file in the PRIMME root directory.

-----------------------------------------------------------------------------
4.	PRIMME MEX Compilation 
-----------------------------------------------------------------------------
For general information about building a MATLAB MEX file, please refer to 
www.mathworks.com/help/MATLAB/MATLAB_external/building-mex-files.html.

It is not necessary to build PRIMME MEX file if you don't have any special
modifications. For experts, if the PRIMME_mex.c file is modified or other
files are changed, Generally, there are two steps to build a MATLAB MEX file.
 
Firstly, make sure that proper gcc versions supported by current MATLAB version
are installed in your operating system. At present, PRIMME MEX has been tested
in SUSE linux. MATLAB 2012a supports gcc 4.4.6 version while MATLAB 2010a 
supports gcc 4.2.5 version. However, if you have older gcc version in your 
linux system, MATLAB should support them too. It is recommended to install 
the required gcc version in your linux system.

Secondly, after right gcc version is installed, you can build PRIMME_mex.c 
mex-file in the terminal at the root directory of the PRIMME MEX or in the 
MATLAB command prompt. Recommended command is as follows:

mex -v -O -largeArrayDims PRIMME_mex.c -L/$TOP/PRIMME -I/$TOP/PRIMME/PRIMMESRC/
COMMONSRC -lprimme -lm  -lmwlapack -lmwblas 

In the mex command, "$TOP" is your path that PRIMME is installed. For 64 bit 
Suse linux, it generates PRIMME_mex binary MEX-file which is like M-file in 
the MATLAB. 

-----------------------------------------------------------------------------
5.	Workflow in PRIMME MEX 
-----------------------------------------------------------------------------
Our MATLAB MEX interface is just as simple as MATLAB build-in eigs function
and svds function. Thus, it is very easy to use PRIMME MEX which is a
powerful tool to solve eigenpair problems and singular value problems. 
There are currently three-layers function call in MATLAB as shown below:

1) MATLAB users' function call

2) primme_eigs or primme_svds MATLAB function calls

3) PRIMME_mex MEX function call 

Frist layer is user's MATLAB script or a MATLAB function in which there is 
at least a large sparse symmetric or hermitian matirx as input argument. Also, 
users can specify more input arguments in this layer so as to obtain desired 
results from PRIMME MEX. You can refer to our examples to see how to set input 
arguments and call primme_eigs or primme_svds functions. 

Second layer is primme_eigs or primme_svds MATLAB functions which serve 
as MATLAB interface between PRIMME_mex function and users' function, which 
calls PRIMME_mex function and returns the outputs to users's function. These
functions perform the initializations, specify the matvec function, call 
PRIMME_mex MEX function and return the outputs to users. It is recommended 
to call primme_eigs and primme_svd functions to solve the eigenpair 
problems and singular value problems.
 
Third layer is PRIMME_mex MEX function. It receives the inputs from the 
primme_eigs or primme_svds MATLAB function calls and then transforms the 
MATLAB inputs into C inputs in order to call PRIMME library. After PRIMME 
returns the results to PRIMME_mex MEX function, it transforms C outputs 
into MATLAB outputs so as to return them to MATLAB. 

Its classical workflow can be illustrated simply as follows:

1) In a MATLAB function call or MATLAB script, user speficy appropriate input
arguments such as Matrix, number of eigenvalues, target, method and other opts.
Then call primme_eigs or primme_svds functions.

2) The primme_eigs or primme_svds functions receives the inputs from
users' function call, performs the initialization, specifies the matvec
function, and then calls PRIMME_mex MEX function. 

3) The PRIMME_mex MEX function receives the inputs from the primme_eigs 
or primme_svds function call and then transforms the MATLAB inputs into 
C inputs in order to call dprimme function or zprimme function in PIRMME.

4) The dprimme function or zprimme function calls other functions in PRIMME.
When PRIMME needs blocked matrix-vector calculation or applies preconditioner,
it calls back PRIMME_mex MEX function to call subroutines like MatrixMatvec_d, 
MatrixMatvec_z, Preconditioner_d, and Preconditioner_z.

5) In these subroutines, C inputs are translated into MATLAB inputs and then
call corresponding default or user-defined MATLAB matvec function or user-
defined preconditioner function. After finishing execution, it returned MATLAB
outputs back to these subrotines, where MATLAB outputs are translated into 
C inputs to return intermediate results to PRIMME. Depending to the dimensions
of Matrix and number of eigenvalues, it usually performs hundreds of thousands
of times matvec and precondition operations.

6) When dprimme function or zprimme function returns outputs back to 
PRIMME_mex MEX function, it transforms C outputs into MATLAB outputs so 
as to return them to primme_eigs or primme_svds functions. At the same
time, if user specifies outputFileName, it saves the primme configuration and
corresponding results into this file.

7) When primme_eigs or primme_svds functions received the MATLAB outputs
from PRIMME_egis MEX function, they return these results to user's function. 

8) When MATLAB function call or MATLAB script receives the results returned
from primme_eigs or primme_svds function, processes them and saves them 
into the same text or binary file. (option)

-----------------------------------------------------------------------------
6.   MATLAB Function Call for solving eigenpair problems 
-----------------------------------------------------------------------------  
PRIMME MEX contains two parts: 1) primme_eigs: seek a few of eigenvalues
and eigenvectors; 2) primme_svds: find a few of singular values and vectors.

The primme_eigs is a MATLAB function call which is served as MATLAB inter-
face between PRIMME_mex amd users' function calls, which calls PRIMME_mex
and returns the outputs to users. If A is M-by-N, a few eigenvalues and 
eigenvectors of A are returned by PRIMME_mex(flag,dim,...), where implicit 
input matrix is A. 

Like eigs() function in the matlab, we provide different level function
calls to satisfy users' demands:

Input: [A, numEvals, target, opts, eigsMethod, P]

Output: [evals, evecs, norms, primmeout]

Function call:  
primme_eigs(A)
primme_eigs(A, numEvals)
primme_eigs(A, numEvals, target)
primme_eigs(A, numEvals, target, opts)
primme_eigs(A, numEvals, target, opts, eigsMethod)
primme_eigs(A, numEvals, target, opts, eigsMethod, P)
primme_eigs(A, numEvals, target, opts, eigsMethod, P1,P2)
primme_eigs(A, numEvals, target, opts, eigsMethod, Pfun)
primme_eigs(Afun, dim,...)

primme_D = primme_eigs(A) returns a vector of A's 6 largest algebraic
eigenvalues. A must be real symmetric or complex hermitian and should
be large and sparse. The primme_eigs(Afun, dim) accepts a function 
AFUN instead of the matrix A. AFUN is a function handle and y = Afun(x) 
returns the matrix-vector product A*x. In all these primme_eigs
function syntaxes, primme_eigs(A,...) could be replaced by
primme_eigs(Afun, dim,...). More examples about how to use function Afun
are presented in the Primme_eigsTest6, 7, 8, and 9 in the root directory
of PRIMME_MEX folder.

[primme_V, primme_D] = primme_eigs(A) returns a diagonal matrix of
primme_D of A's 6 largest algebraic eigenvalues and a matrix primme_V
whose columns are the corresponding eigenvectors. 

[primme_V, primme_D, norms, primmeout] = primme_eigs(A, numEvals)
returns a diagonal matrix of primme_D of A's numEvals largest algebraic 
eigenvalues, a matrix primme_V whose columns are the corresponding 
eigenvectors, a double array of the residual norm of eigenvalues  and
a struct to report statistical information about numOuterIterations,
numRestarts, numMatvecs and numPreconds. numEvals is the number of 
eigenvalues that users want to find. It must be less than dimention of 
the matrix A.

primme_eigs(A, numEvals, target) returns numEvals target eigenvlaues.
target could be a string like below:
 'LA' ------ primme_largest (default)
 'SA' ------ primme_smallest    
 'CGT'------ primme_closest_geq  
 'CLT'------ primme_closest_leq  
 'CT' ------ primme_closest_abs   

primme_eigs(A, numEvals, target, opts, eigsMethod) specify options that 
are listed and explained in the last. eigsMethod is the solver method
the PRIMME uses to find eigenvalues and eigenvectors. eigsMethod could
be:
typedef enum{
DYNAMIC, (default)        ---0: Switches dynamically to the best method
DEFAULT_MIN_TIME,         ---1: Currently set at JDQMR_ETol
DEFAULT_MIN_MATVECS,      ---2: Currently set at GD+block
Arnoldi,                  ---3: obviously not an efficient choice 
GD,                       ---4: classical block Generalized Davidson 
GD_plusK,                 ---5: GD+k block GD with recurrence restarting
GD_Olsen_plusK,           ---6: GD+k with approximate Olsen precond.
JD_Olsen_plusK,           ---7: GD+k, exact Olsen (two precond per step)
RQI,                      ---8: Rayleigh Quotient Iteration. Also INVIT,
                            :   but for INVIT provide targetShifts
JDQR,                     ---9: Original block, Jacobi Davidson
JDQMR,                   ---10: Our block JDQMR method (similar to JDCG)
JDQMR_ETol,              ---11: Slight, but efficient JDQMR modification
SUBSPACE_ITERATION,      ---12: equiv. to GD(block,2*block)
LOBPCG_OrthoBasis,       ---13: equiv. to GD(nev,3*nev)+nev
LOBPCG_OrthoBasis_Window ---14: equiv. to GD(block,3*block)+block nev>block

primme_eigs(A, numEvals, target, opts, eigsMethod, P) uses preconditioner 
P or P = P1*P2 to solve eigenvalue problem for large sparse matrix.
If P is [] then a preconditioner is not applied. P may be a function 
handle Pfun such that Pfun(x) returns P\x.

The opts contains 31 different options which can be used by experts to
achieve the best performance. The most used options are:
opts.aNorm: the estimate norm value of matrix A [{0.0}|scaler]
opts.eps: desired computing accuracy [{1e-12}|scaler]
opts.maxBlockSize: maximum block size the PRIMME uses [{1}|scaler]
opts.printLevel: different level reporting(0-5) [{1}|scaler]
opts.outputFile: output file name where user wants to save results
opts.precondition: set to 1 if use preconditioner [{0}|1]
opts.isreal: the complexity of A represented by afun [{ture}|false]
opts.numTargetShifts: number of shifts for interior eigenvalues [{0}|scaler]
opts.targetShifts: shifts for interior eigenvalues [{}|vector]
 
-----------------------------------------------------------------------------
7.   MATLAB Function Call for solving singular value problems 
-----------------------------------------------------------------------------      
PRIMME_mex MEX file can also be used to solve singular value problem. There
are two different shemes in primme_svds: 1) primme_svds_ATA; 2) primme_svds_OAAO. 
The primme_svds_ATA is a faster scheme while the primme_svds_OAAO is more 
accurate one. Users can specify primme_svds scheme by setting svdsMethod.  

The primme_svds_OAAO is a more robust scheme that finds a few singular 
values and vectors. If A is M-by-N, a few singular values and vectors of A 
are obtained by seeking a few eigenvalues and eigenvectors returned by 
PRIMME_mex(flag,dim,...), where implicit input matrix = [sparse(N,N) A'; 
A sparse(M,M)].

The primme_svds_ATA is a more faster scheme that finds a few singular 
values and vectors. If A is M-by-N, a few singular values and vectors of A 
are obtained by seeking a few eigenvalues and eigenvectors returned by 
PRIMME_mex(flag,dim,...), where implicit input matrix = A'*A. This 
implementation assumes M > N. If M < N, simply seek the singular values 
of transpose of A and then obtain corresponding singular values of A.

primme_svds is an interface between user function call and PRIMME_mex 
function. The primme_svds function requires at least one input (input 
matrix A) and output at least singular values of A.

Like svds() function in the Matlab, we provide different level function 
calls to satisfy users' demands:

Input: [A, numSvds, target, opts, eigsMethod, svdsMethod, P]

Output: [PRIMME_U, PRIMME_S, PRIMME_V , norms, primmeout]

Function call:  
primme_svds(A)
primme_svds(A, numSvds)
primme_svds(A, numSvds, target)
primme_svds(A, numSvds, target, opts)
primme_svds(A, numSvds, target, opts, eigsMethod)
primme_svds(A, numSvds, target, opts, eigsMethod, svdsMethod)
primme_svds(A, numSvds, target, opts, eigsMethod, svdsMethod, P)
primme_svds(A, numSvds, target, opts, eigsMethod, svdsMethod, P1,P2)
primme_svds(A, numSvds, target, opts, eigsMethod, svdsMethod, Pfun)       
primme_svds(Afun, M, N,...)

primme_S = primme_svds(A) returns the 6 largest singular values of A.
A could be any matrix and should be large and sparse. The 
primme_svds(Afun, M, N) accepts a function AFUN instead of the matrix A. 
AFUN(X,'notransp') accepts a vector input X and returns the matrix-vector 
product A*X while AFUN(X,'transp') returns A'*X. In all these primme_svds
function syntaxes, primme_svds(A,...) could be replaced by 
primme_svds(Afun, M, N,...). More examples about how to use function Afun
are presented in the Primme_svdsTest6, 7, 8, and 9 in the root directory
of PRIMME_MEX folder.

[primme_U, primme_S, primme_V] = primme_svds(A,...) seeks the singular
vectors as well. If A is M-by-N and numSvds singular values are sought,
then the left singular vecotr primme_U is M-by-numSvds with orthonormal 
columns, primme_S is numSvds-by-numSvds diagonal matrix with singular 
values in decreasing order, and the right singular vector primme_V is 
N-by-numSvds with orthonormal columns. 

[primme_U, primme_S, primme_V, norms, primmeout] = primme_svds(A,numSvds)
also returns a double array of the residual norm of singular values and
a struct to report statistical information about numOuterIterations,
numRestarts, numMatvecs and numPreconds. numSvds is the number of 
singular values that usrs want. It must be less than minimum dimention 
of the matrix A.

primme_svds(A, numSvds, target) returns numSvds target singular values.
target could be a string like below:
 'LA' ------ primme_largest (default)
 'SA' ------ primme_smallest    
 'CGT'------ primme_closest_geq  
 'CLT'------ primme_closest_leq  
 'CT' ------ primme_closest_abs   

primme_svds(A, numSvds, target, opts, eigsMethod) specify options that 
are listed and explained in the last. eigsMethod is the solver method
the PRIMME uses to find eigenvalues and eigenvectors for implicit input
matrix A'*A. eigsMethod could
be:
typedef enum{
DYNAMIC, (default)        ---0: Switches dynamically to the best method
DEFAULT_MIN_TIME,         ---1: Currently set at JDQMR_ETol
DEFAULT_MIN_MATVECS,      ---2: Currently set at GD+block
Arnoldi,                  ---3: obviously not an efficient choice 
GD,                       ---4: classical block Generalized Davidson 
GD_plusK,                 ---5: GD+k block GD with recurrence restarting
GD_Olsen_plusK,           ---6: GD+k with approximate Olsen precond.
JD_Olsen_plusK,           ---7: GD+k, exact Olsen (two precond per step)
RQI,                      ---8: Rayleigh Quotient Iteration. Also INVIT,
                            :   but for INVIT provide targetShifts
JDQR,                     ---9: Original block, Jacobi Davidson
JDQMR,                   ---10: Our block JDQMR method (similar to JDCG)
JDQMR_ETol,              ---11: Slight, but efficient JDQMR modification
SUBSPACE_ITERATION,      ---12: equiv. to GD(block,2*block)
LOBPCG_OrthoBasis,       ---13: equiv. to GD(nev,3*nev)+nev
LOBPCG_OrthoBasis_Window ---14: equiv. to GD(block,3*block)+block nev>block

primme_svds(A, numSvds, target, opts, eigsMethod, svdsMethod) also
choose which scheme to solve singular triplet problem. Currently,
svdsMethod could be a sting like below:
1) svdsMethod = 'ATA' : choose primme_svds_ATA scheme (default)
2) svdsMethod = 'OAAO': choose primme_svds_OAAO scheme

primme_svds(A, numSvds, target, opts, eigsMethod, svdsMethod, P) uses 
preconditioner P or P = P1*P2 to solve singular triplet problem for 
large sparse matrix. If P is [] then a preconditioner is not applied. 
P may be a function handle Pfun such that Pfun(X,'notransp')returns P\X 
and Pfun(X,'transp') returns P'\X.             

The opts contains 32 different options which can be used by experts to
get the best performance. The most used options are:
opts.aNorm: the estimate norm value of matrix A [{0.0}|scaler]
opts.eps: desired computing accuracy [{1e-12}|scaler]
opts.maxBlockSize: maximum block size the PRIMME uses [{1}|scaler]
opts.printLevel: different level reporting(0-5) [{1}|scaler]
opts.outputFile: output file name where user wants to save results
opts.precondition: set to 1 if use preconditioner [{0}|1]
opts.isreal: the complexity of A represented by afun [{ture}|false]
opts.svdsPrecondType: specify if preconditioner is provided directly for 
ATA or OAAO. [true|{false}]
opts.numTargetShifts: number of shifts for interior eigenvalues [{0}|scaler]
opts.targetShifts: shifts for interior eigenvalues [{}|vector]
 
-----------------------------------------------------------------------------
 8.	Input and Output Descriptions of PRIMME MEX-file
-----------------------------------------------------------------------------
 PRIMME MEX-file provides different level function calls to satisfy users' demands:

 Syntax:
 
 [evals, evecs, norms, primmeout] = PRIMME_mex(flag, dim, Numeigs, target, method, opts)
  
  evals = PRIMME_mex(flag, dim) 
  evals = PRIMME_mex(flag, dim, numEvals)
  evals = RRIMME_eigs(flag, dim, numEvals, target)
  evals = PRIMME_mex(flag, dim, numEvals, target, method)
  evals = PRIMME_mex(flag, dim, numEvals, target, method, opts)
 [evals,evecs] = PRIMME_mex(flag, dim, numEvals, target, method, opts)
 [evals, evecs, norms] = PRIMME_mex(flag, dim, numEvals, target, method, opts)
 [evals, evecs, norms, primmeout] = PRIMME_mex(flag, dim, numEvals, target, method, opts)

 Output desciption: [evals, evecs, norms, primmeout]

 [1] evals: a double array of the target eigenvalues returned to MATLAB

 [2] evecs: a double array of the corresponding eigenvectors returned to MATLAB

 [3] norms: a double array of the residual norm of the eigenvalues

 [4] primmeout: a struct to report statistics which contains four members:
    (1) numOuterIterations
    (2) numRestarts
    (3) numMatvecs
    (4) numPreconds

 Input description: (flag, dim, numEvals, target, method, opts) 

 [1] flag: mark the input matrix is real or complex
 
 [2] dim: the dimension of large and sparse symmetric and hermitian matrix.
 
 [3] numEvals: number of eigenvalues required

 [4] target: Which eigenpairs to find.  target can be any of the following enum:
 
    primme_smallest --- 'SA':Smallest algebraic eigenvalues. Target shifts ignored
    primme_largest  --- 'LA':Largest  algebraic eigenvalues. Target shifts ignored
    primme_closest_geq -'CGT':Closest to, but greater or equal than a set of shifts
    primme_closest_leq -'CLT':Closest to, but less or equal than a set of shifts
    primme_closest_abs - 'CT':Closest in absolute value to a set of shifts

 [5] method: which solver method to choose. method is an enum type below:
 
    typedef enum{
    DYNAMIC,                  ---0: Switches dynamically to the best method
    DEFAULT_MIN_TIME,         ---1: Currently set at JDQMR_ETol
    DEFAULT_MIN_MATVECS,      ---2: Currently set at GD+block
    Arnoldi,                  ---3: obviously not an efficient choice 
    GD,		              ---4: classical block Generalized Davidson 
    GD_plusK,		      ---5: GD+k block GD with recurrence restarting
    GD_Olsen_plusK,           ---6: GD+k with approximate Olsen precond.
    JD_Olsen_plusK,           ---7: GD+k, exact Olsen (two precond per step)
    RQI,          	      ---8: Rayleigh Quotient Iteration. Also INVIT,
     			          :   but for INVIT provide targetShifts
    JDQR,          	      ---9: Original block, Jacobi Davidson
    JDQMR,          	      ---10: Our block JDQMR method (similar to JDCG)
    JDQMR_ETol,               ---11: Slight, but efficient JDQMR modification
    SUBSPACE_ITERATION,       ---12: equiv. to GD(block,2*block)
    LOBPCG_OrthoBasis,        ---13: equiv. to GD(nev,3*nev)+nev
    LOBPCG_OrthoBasis_Window  ---14: equiv. to GD(block,3*block)+block nev>block
    } primme_preset_method; 

 [6] 'opts' is an option structure which contain following parameters in the 
     primme_params structure: 
    (0) opts.aNorm: the estimate norm value of matrix A [default = 0]
    (1) opts.eps: desired computing accuracy [default = 1e-12]
    (2) opts.numTargetShifts: shifts for interior eigenvalues [default = 0]
    (3) opts.targetShifts: pointer to get each shift for interior eigenvalues
    (4) opts.initSize: initial guesses/constraints [default = 0]
    (5) opts.numOrthoConst: [default = 0]
    (6) opts.locking: 0 or 1 
    (7) opts.dynamicMethodSwitch: from -3 to 1
    (8) opts.maxBasisSize: maximum basis size allowed in the main iteration
    (9) opts.minRestartSize: minimum Ritz vectors to restart
   (10) opts.maxBlockSize:  [default = 1]
   (11) opts.maxMatvecs: [default = INT_MAX]
   (12) opts.maxOuterIterations: [default = INT_MAX]
   (13) opts.restartingParams.scheme: [default = primme_thick]
   (14) opts.restartingParams.maxPrevRetain: [default = 1]
   (15) opts.precondition: set to 1 if preconditioning is to be performed,
        make sure the applyPreconditioner is not NULL [default = 0]
   (16) opts.robustShifts: set to 1 to use robustShifting
   (17) opts.maxInnerIterations: 0 (no inner iterations), = k (perform 
        at most k inner iterations per outer step)
   (18) opts.LeftQ: 0 or 1
   (19) opts.LeftX: 0 or 1
   (20) opts.RightQ: 0 or 1
   (21) opts.RightX: 0 or 1
   (22) opts.SkewQ: 0 or 1
   (23) opts.SkewX: 0 or 1
   (24) opts.relTolBase: a legacy from calssical JDQR
   (25) opts.convTest: how to stop the inner QMR method
   (26) opts.printLevel: 0-5 (different level reporting) [default = 1]
   (27) opts.outputFile: output file name where user wants to save results
   (28) opts.iseed: set iseed value for initialization
   (29) opts.intWorkSize: memory size for int workspace
   (30) opts.realWorkSize: memory size for real workspace        

-----------------------------------------------------------------------------
9.	Examples 
-----------------------------------------------------------------------------
In order to illustrate how to use PRIMME MEX interface, we provide eighteen 
examples in two parts. The first part is to use primme_eigs() function to seek
a few of eigenvalues and eigenvectors for eigenpair problems. Nine test case 
examples are provided.

-----------------------------------------------------------------------------
 Part I: sovling eigenpair problems
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
       [V ,D] = eigs(@(x)test_eigs(x), n , numEvals, 'LA', opts);
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
% Senior test case for user to provide matrix function instread of the 
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
% Senior test case for user to provide function handle instread of the 
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


-----------------------------------------------------------------------------
 Part II: sovling singular value problems
-----------------------------------------------------------------------------
The second part is to use primme_svds_ATA() or primme_svds_OAAO() function to 
seek a few of singular values and vectors for singular value problems. Nine test 
case examples are provided.

Example 10: Primme_svdsTest.m

%
% Simple test case for primme_svds (using default ATA eigsMethod) to seek
% largest singular triplets
% Function call: primme_svds(A, numSvds, target)
%
       clear all
       load Andrews.mat;
       
       if (exist('A','var') == 0)
          A = Problem.A;
       end
      
       primme_start = tic; 
       numSvds = 5;
       target = 'LA';
      
       [primme_U, primme_S, primme_V]= primme_svds(A, numSvds, target);  
       primme_S
       primme_svds_telapsed = toc(primme_start)

       svds_start = tic;
       OPTIONS.tol = 1e-12;
       OPTIONS.maxit = 65536;
       [U, S, V] = svds(A, numSvds, 'L', OPTIONS);
       S
       svds_telapsed = toc(svds_start)      
       
Results:

primme_S =

   36.4853         0         0         0         0
         0   36.4492         0         0         0
         0         0   36.4179         0         0
         0         0         0   36.1741         0
         0         0         0         0   36.0190


primme_svds_telapsed =

    2.4846


S =

   36.4853         0         0         0         0
         0   36.4492         0         0         0
         0         0   36.4179         0         0
         0         0         0   36.1741         0
         0         0         0         0   36.0190


svds_telapsed =

    3.9968


Example 11: Primme_svdsTest2.m

%
% Simple test case for primme_svds using OAAO eigsMethod to seek smallest singular
% triplets
% Function call: primme_svds(A, numSvds, target, opts, eigsMethod, svdsMethod)
%
       
       clear all
       load Andrews.mat; 
       if (exist('A') == 0)
           A = Problem.A;
       end

       primme_start = tic;
       numSvds = 5;
       target = 'SA';
       eigsMethod = 0;
       svdsMethod = 'OAAO';
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'printLevel', 2, 'outputFileName', {'sampleout'}); 
     
       [primme_U, primme_S, primme_V]= primme_svds(A, numSvds, target, opts, eigsMethod, svdsMethod);
       primme_S
       primme_svds_telapsed = toc(primme_start)
       

       svds_start = tic;
       OPTIONS.tol = 1e-10;
       OPTIONS.maxit = 65536;
       [U, S, V] = svds(A, numSvds, 0, OPTIONS);
       S
       svds_telapsed = toc(svds_start)
       
Results:

primme_S =

   1.1839e-01            0            0            0            0
            0   1.1192e-01            0            0            0
            0            0   8.8282e-02            0            0
            0            0            0   6.6916e-02            0
            0            0            0            0   1.0570e-15


primme_svds_telapsed =

  595.9585

Warning: NORMEST did not converge for 100 iterations with tolerance 1e-06 
> In normest at 41
  In svds at 146
  In Primme_svdsTest2 at 35 

S =

   1.1839e-01            0            0            0            0
            0   1.1192e-01            0            0            0
            0            0   8.8282e-02            0            0
            0            0            0   6.6916e-02            0
            0            0            0            0   1.7463e-14


svds_telapsed =

  235.2963


Example 12: Primme_svdsTest3.m
%
% Senior test case for user to provide matrix function instread of the 
% matrix A using ATA eigsMethod to seek singular triplets. The primme_svds  
% accepts a function handle AFUN instead of the matrix A. AFUN(X,'notransp')
% accepts a vector input X and returns the matrix-vector product A*X 
% while AFUN(X,'transp') returns A'*X. A Could be any matrix.
% Function call: primme_svds(Afun, M, N, numSvds, target, opts, eigsMethod)
%

       clear all
       k = 1000;
       on = ones(k,1); 
       
       primme_start = tic; 
       m = k;
       n = k;
       numSvds = 5;
       target = 'LA';
       eigsMethod = 2;
       
       opts = struct('aNorm',0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'printLevel', 2, 'outputFileName', {'sampleout'});
      
       [primme_U, primme_S, primme_V]= primme_svds(@(x,tflag)svdsMatvec(x,on,k,tflag), m, n, numSvds, target, opts, eigsMethod);     
       primme_S
       primme_svds_telapsed = toc(primme_start)
       
       A = spdiags([-2*on 4*on -2*on],-1:1,k,k);
       svds_start = tic;
       OPTIONS.tol = 1e-10;
       OPTIONS.maxit = 65536;
       [U, S, V] = svds(A, numSvds, 'L', OPTIONS);
        S
       svds_telapsed = toc(svds_start)

Results:

primme_S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


primme_svds_telapsed =

   2.1293e+00


S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


svds_telapsed =

   4.2368e+00

Example 13: Primme_svdsTest4.m
%
% Senior test case for user to provide matrix function instread of the 
% matrix A using OAAO eigsMethod to seek singular triplets. The primme_svds 
% accepts a function handle AFUN instead of the matrix A. AFUN(X,'notransp')
% accepts a vector input X and returns the matrix-vector product A*X 
% while AFUN(X,'transp') returns A'*X. A Could be any matrix.
% Function call: primme_svds(Afun, M, N, numSvds, target, opts, eigsMethod, svdsMethod)
%

       clear all
       k = 1000;
       on = ones(k,1); 
       
       primme_start = tic; 
       m = k;
       n = k;
       numSvds = 5;
       target = 'LA';
       eigsMethod = 2;
       svdsMethod = 'OAAO';
       
       opts = struct('aNorm',0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'printLevel', 2, 'outputFileName', {'sampleout'});
      
       [primme_U, primme_S, primme_V]= primme_svds(@(x,tflag)svdsMatvec(x,on,k,tflag), m, n, numSvds, target, opts, eigsMethod, svdsMethod); 
       primme_S
       primme_svds_telapsed = toc(primme_start)
       
       A = spdiags([-2*on 4*on -2*on],-1:1,k,k);
       svds_start = tic;
       OPTIONS.tol = 1e-10;
       OPTIONS.maxit = 65536;
       [U, S, V] = svds(A, numSvds, 'L', OPTIONS);
       S
       svds_telapsed = toc(svds_start)

Results:

primme_S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


primme_svds_telapsed =

   4.5547e+00


S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


svds_telapsed =

   4.4275e+00


Example 14: Primme_svdsTest5.m
%
% Simple test case for primme_svds using ATA eigsMethod to seek singlar
% triplets of complex hermitian matrix 3Dspectralwave2 from university of
% florida. If the input is function handle instead of a matrix, the usrs
% must specify opts.isreal feild. Set to 1 if the input is real; set to 0
% if the input is complex. If the user does not specify opts.isreal, the
% default case is the real input matrix or function handle.
% Function call: primme_svds(A, numSvds, target, opts)
%
       
       clear all
       load 3Dspectralwave2.mat;
       
       if (exist('A','var') == 0)
          A = Problem.A;
       end

       primme_start = tic;
       numSvds = 5;
       target = 'LA';
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'printLevel', 2, 'outputFileName', {'sampleout'}); 
     
       [primme_U, primme_S, primme_V]= primme_svds(A, numSvds, target, opts);
       primme_S
       primme_svds_telapsed = toc(primme_start)
       

       svds_start = tic;
       OPTIONS.tol = 1e-10;
       OPTIONS.maxit = 65536;
       [U, S, V] = svds(A, numSvds, 'L', OPTIONS);
       S
       svds_telapsed = toc(svds_start)

Results:

primme_S =

   6.9065e+01            0            0            0            0
            0   6.8789e+01            0            0            0
            0            0   6.7348e+01            0            0
            0            0            0   6.7217e+01            0
            0            0            0            0   6.7111e+01


primme_svds_telapsed =

   5.2051e+01


S =

   6.9065e+01            0            0            0            0
            0   6.8789e+01            0            0            0
            0            0   6.7348e+01            0            0
            0            0            0   6.7217e+01            0
            0            0            0            0   6.7111e+01


svds_telapsed =

   1.0423e+02


Example 15: Primme_svdsTest6.m
%
% Senior test case for user to provide their preconditioner to speed up the
% whole process. The input A is a matrix and preconditioner P is also a
% matrix.
% Function call: primme_svds(A, numSvds, target, opts, eigsMethod, P)
%
       clear all
       k = 1000;
       on = ones(k,1); 
       A = spdiags([-2*on 4*on -2*on],-1:1,k,k);
       

       numSvds = 5;
       target = 'LA';
       eigsMethod = 1;
       svdsMethod = 'ATA';
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'precondition', 1,'printLevel', 2, 'outputFileName', {'sampleout'}); 
       
       P = spdiags(4*on,0,k,k);

       primme_start = tic;
       [primme_U, primme_S, primme_V]= primme_svds(A, numSvds, target, opts, eigsMethod, svdsMethod, P);
       primme_S
       primme_svds_telapsed = toc(primme_start)
    
 
        
       svds_start = tic;
       OPTIONS.tol = 1e-10;
       OPTIONS.maxit = 65536;
       [U, S, V] = svds(A, numSvds, 'L', OPTIONS);
       S
       svds_telapsed = toc(svds_start)

Results:

primme_S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


primme_svds_telapsed =

   6.5965e-01


S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


svds_telapsed =

   5.2035e+00


Example 16: Primme_svdsTest7.m
%
% Senior test case for user to provide their preconditioner to speed up the
% whole process. The input A is a matrix and preconditioner P1 and P2 are
% both preconditioners.
% Function call: primme_svds(A, numSvds, target, opts, eigsMethod, P1,P2)
%
       clear all
       k = 1000;
       on = ones(k,1); 
       A = spdiags([-2*on 4*on -2*on],-1:1,k,k);
       

       numSvds = 5;
       target = 'LA';
       eigsMethod = 2;
       svdsMethod = 'ATA';
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'precondition', 1,'printLevel', 2, 'outputFileName', {'sampleout'}); 
       
       P1 = spdiags([on/(-2) on],-1:0,k,k); 
       P2 = spdiags([4*on -on],0:1,k,k);

       primme_start = tic;
       [primme_U, primme_S, primme_V]= primme_svds(A, numSvds, target, opts, eigsMethod, svdsMethod, P1, P2);
       primme_S
       primme_svds_telapsed = toc(primme_start)
    
 
        
       svds_start = tic;
       OPTIONS.tol = 1e-10;
       OPTIONS.maxit = 65536;
       [U, S, V] = svds(A, numSvds, 'L', OPTIONS);
       S
       svds_telapsed = toc(svds_start)

Results:

primme_S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


primme_svds_telapsed =

   3.1477e+00


S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


svds_telapsed =

   4.6079e+00


Example 17: Primme_svdsTest8.m
%
% Senior test case for user to provide their preconditioner to speed up the
% whole process. The matrix P is provided as preconditioner. The primme_svds 
% accepts a function handle AFUN instead of the matrix A. AFUN(X,'notransp')
% accepts a vector input X and returns the matrix-vector product A*X 
% while AFUN(X,'transp') returns A'*X.
% Function call: primme_svds(Afun, M, N, numSvds, target, opts, eigsMethod, P)
%
       clear all
       k = 1000;
       on = ones(k,1); 

       m = k;
       n = k;
       numSvds = 5;
       target = 'LA';
       eigsMethod = 2;
       svdsMethod = 'ATA';
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'precondition', 1,'printLevel', 2, 'outputFileName', {'sampleout'}); 
       
       P = spdiags(4*on,0,k,k);

       primme_start = tic;
       [primme_U, primme_S, primme_V]= primme_svds(@(x,tflag)svdsMatvec(x,on,k,tflag), m, n, numSvds, target, opts, eigsMethod, svdsMethod, P);
       primme_S
       primme_svds_telapsed = toc(primme_start)
    
 
       A = spdiags([-2*on 4*on -2*on],-1:1,k,k);
       svds_start = tic;
       OPTIONS.tol = 1e-10;
       OPTIONS.maxit = 65536;
       [U, S, V] = svds(A, numSvds, 'L', OPTIONS);
       S
       svds_telapsed = toc(svds_start)

Results:

primme_S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


primme_svds_telapsed =

   2.3180e+00


S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


svds_telapsed =

   4.2569e+00


Example 18: Primme_svdsTest9.m
%
% Senior test case for user to provide their preconditioner to speed up the
% whole process. The primme_svds accepts a function handle AFUN instead of 
% the matrix A. AFUN(X,'notransp') accepts a vector input X and returns the 
% matrix-vector product A*X while AFUN(X,'transp') returns A'*X. The
% primme_svds also accepts a function handle Pfun as a preconditioner
% instead of the matrix P. PFUN(X,'notransp') accepts a vector input X and 
% returns P\X while AFUN(X,'transp') returns P'\X.
% Function call: primme_svds(Afun, M, N, numSvds, target, opts, eigsMethod, Pfun)
%
       clear all
       k = 1000;
       on = ones(k,1); 

       m = k;
       n = k;
       numSvds = 5;
       target = 'LA';
       eigsMethod = 2;
       svdsMethod = 'ATA';
       
       opts = struct('aNorm', 0.0, 'eps', 1e-10, 'maxBlockSize', 1, 'initSize', 1, 'iseed', {[-1, 0, 1, 327221]}, 'precondition', 1,'printLevel', 2, 'outputFileName', {'sampleout'}); 

       primme_start = tic;
       [primme_U, primme_S, primme_V]= primme_svds(@(x,tflag)svdsMatvec(x,on,k,tflag), m, n, numSvds, target, opts, eigsMethod, svdsMethod, @(x,tflag)svdsPrecond(x,on,k,tflag));
       primme_S
       primme_svds_telapsed = toc(primme_start)
    
 
       A = spdiags([-2*on 4*on -2*on],-1:1,k,k);
       svds_start = tic;
       OPTIONS.tol = 1e-10;
       OPTIONS.maxit = 65536;
       [U, S, V] = svds(A, numSvds, 'L', OPTIONS);
       S
       svds_telapsed = toc(svds_start)

Results:

primme_S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


primme_svds_telapsed =

   2.9495e+00


S =

   8.0000e+00            0            0            0            0
            0   7.9999e+00            0            0            0
            0            0   7.9998e+00            0            0
            0            0            0   7.9997e+00            0
            0            0            0            0   7.9995e+00


svds_telapsed =

   4.6180e+00





