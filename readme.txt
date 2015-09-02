
Welcome to PRIMME's documentation!
**********************************

Table of Contents:

* PRIMME: PReconditioned Iterative MultiMethod Eigensolver

* Changelog

* Citing this code

* License Information

* Contact Information

* Directory Structure

* Making and Linking

* Tested Systems

* C Library Interface

* FORTRAN Library Interface


PRIMME: PReconditioned Iterative MultiMethod Eigensolver
********************************************************

PRIMME, pronounced as *prime*, finds a number of eigenvalues and their
corresponding eigenvectors of a real symmetric, or complex hermitian
matrix A. Largest, smallest and interior eigenvalues are supported.
Preconditioning can be used to accelerate convergence. PRIMME is
written in C99, but complete interfaces are provided for Fortran 77
and MATLAB.


Changelog
*********

Changes in PRIMME 1.2.1:

* Added MATLAB interface to full PRIMME functionality.

* Support for BLAS/LAPACK with 64bits integers
  ("-DPRIMME_BLASINT_SIZE=64").

* Simplified configuration of Make_flags and Make_links (removed
  "TOP" variable and replaced defines "NUM_SUM" and "NUM_IBM" by
  "F77UNDERSCORE").

* Replaced directories "DTEST" and "ZTEST" by "TEST", that has:

  * "driver.c": read matrices in MatrixMarket format and PETSc
    binary and call PRIMME with the parameters specified in a file;
    support complex arithmetic and MPI and can use PETSc
    preconditioners.

  * "ex*.c" and "ex*.f": small, didactic examples of usage in C and
    Fortran and in parallel (with PETSc).

* Using Sphinx to manage documentation. Detailed Fortran 77
  interface.

Changes in PRIMME 1.2:

* A Fortran compiler is no longer required for building the PRIMME
  library. Fortran programs can still be linked to PRIMME's F77
  interface.

* Fixed some uncommon issues with the F77 interface.

* PRIMME can be called now multiple times from the same program.

* Performance improvements in the QMR inner solver, especially for
  complex arithmetic.

* Fixed a couple of bugs with the locking functionality.

  * In certain extreme cases where all eigenvalues of a matrix were
    needed.

  * The order of selecting interior eigenvalues.

  The above fixes have improved robustness and performance.

* PRIMME now assigns unique random seeds per parallel process for up
  to 4096^3  (140 trillion processes).

* For the DYNAMIC method, fixed issues with initialization and
  synchronization decisions across multiple processes.

* Fixed uncommon library interface bugs, coordinated better setting
  the method and the user setting of parameters, and improved the
  interface in the sample programs and makefiles.

* Other performance and documentation improvements.


Citing this code
****************

Please cite:

[r1] A. Stathopoulos and J. R. McCombs PRIMME: *PReconditioned
     Iterative MultiMethod Eigensolver: Methods and software
     description*, ACM Transaction on Mathematical Software Vol. 37,
     No. 2, (2010), 21:1-21:30.

More information on the algorithms and research that led to this
software can be found in the rest of the papers. The work has been
supported by a number of grants from the National Science Foundation.

[r2] A. Stathopoulos, *Nearly optimal preconditioned methods for
     hermitian eigenproblems under limited memory. Part I: Seeking one
     eigenvalue*, SIAM J. Sci. Comput., Vol. 29, No. 2, (2007), 481--
     514.

[r3] A. Stathopoulos and J. R. McCombs, *Nearly optimal
     preconditioned methods for hermitian eigenproblems under limited
     memory. Part II: Seeking many eigenvalues*, SIAM J. Sci. Comput.,
     Vol. 29, No. 5, (2007), 2162-2188.

[r4] J. R. McCombs and A. Stathopoulos, *Iterative Validation of
     Eigensolvers: A Scheme for Improving the Reliability of Hermitian
     Eigenvalue Solvers*, SIAM J. Sci. Comput., Vol. 28, No. 6,
     (2006), 2337-2358.

[r5] A. Stathopoulos, *Locking issues for finding a large number
     of eigenvectors of hermitian matrices*, Tech Report: WM-
     CS-2005-03, July, 2005.


License Information
*******************

PRIMME is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 2.1 of the License, or
(at your option) any later version.

PRIMME is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
USA


Contact Information
*******************

For reporting bugs or questions about functionality contact Andreas
Stathopoulos


Directory Structure
*******************

The next directories and files should be available:

* "COPYING.txt", LGPL License;

* "Make_flags",  flags to be used by makefiles to compile library
  and tests;

* "Link_flags",  flags needed in making and linking the test
  programs;

* "PRIMMESRC/",  directory with source code in the following
  subdirectories:

     * "COMMONSRC/", interface and common functions used by all
       precision versions;

     * "DSRC/",      the source code for the double precision
       "dprimme()";

     * "ZSRC/",      the source code for the double complex
       precision "zprimme()";

* "MEX/",          MATLAB interface for PRIMME;

* "TEST/",         driver and samples in C and F77, both sequential
  and parallel;

* "libprimme.a",   the PRIMME library (to be made);

* "makefile"       main make file;

* "readme.txt"     text version of the documentation;

* "doc/"           directory with the HTML and PDF versions of the
  documentation.


Making and Linking
******************

"Make_flags" has the flags and compilers used to make "libprimme.a":

* *CC*, compiler program such as "gcc", "clang" or "icc".

* *CFLAGS*, compiler options such as "-g" or "-O3". Also include
  some of these options if the BLAS and LAPACK that will be linked:

  * "-DF77UNDERSCORE", the Fortran function names has appended an
    underscore (usually they does).

  * "-DPRIMME_BLASINT_SIZE=64", integers are 64-bit integer
    ("kind=8") type (usually they doesn't).

After customizing "Make_flags", type this to generate "libprimme.a":

   make lib

Making can be also done at the command line:

   make lib CC=clang CFLAGS='-O3'

"Link_flags" has the flags for linking with external libraries and
making the executables located in "TEST":

* *LDFLAGS*, linker flags such as "-framework Accelerate".

* *LIBS*, flags to link with libraries (BLAS and LAPACK are
  required), such as "-lprimme -llapack -lblas -lgfortran -lm".

After that type this to compile and execute a simple test:

   make test

If it worked, try with other examples in "TEST" (see "README" in
"TEST" for more information about how to compile the driver and the
examples).

Full description of actions that *make* can take:

* *make lib*, builds "libprimme.a"; alternatively:

* *make libd*, if only "dprimme()" is of interest, build
  "libdprimme.a":

* *make libz*, if only "zprimme()" is of interest, build
  "libzprimme.a";

* *make test*, build and execute a simple example;

* *make clean*, removes all "*.o", "a.out", and core files from all
  directories.


Tested Systems
**************

PRIMME is primary developed with GNU gcc, g++ and gfortran (versions
4.8 and later). Many users have reported builds on several other
platforms/compilers:

* SUSE 13.1 & 13.2

* CentOS 6.6

* Ubuntu 14.04

* MacOS X 10.9 & 10.10

* Cray XC30

* SunOS 5.9, quad processor Sun-Fire-280R, and several other
  UltraSparcs

* AIX 5.2 IBM SP POWER 3+, 16-way SMP, 375 MHz nodes (seaborg at
  nersc.gov)

C Library Interface
*******************

The next enumerations and functions are declared in "primme.h".

primme_preset_method

   Enumeration of preset configurations.

   DYNAMIC

      Switches to the best method dynamically; currently, between
      "JDQMR_ETol" and "GD_Olsen_plusK". Set "dynamicMethodSwitch" to
      1.

   DEFAULT_MIN_TIME

      Currently set as "JDQMR_ETol"; this method is usually the
      fastest if the cost of the matrix vector product is the order of
      the matrix dimension.

   DEFAULT_MIN_MATVECS

      Currently set as "GD_Olsen_plusK"; this method usually spent
      less matrix vector products than the others, so it's a good
      choice when this operation is expensive.

   Arnoldi

      Arnoldi implemented à la Generalized Davidson.

   GD

      Generalized Davidson.

   GD_plusK

      GD with locally optimal restarting. See "maxPrevRetain".

   GD_Olsen_plusK

      GD+k and the cheap Olsen's Method (set "RightX" to 1 and "SkewX"
      to 0). See "SkewX".

   JD_Olsen_plusK

      GD+k and Olsen's Method (set "RightX" to 1 and "SkewX" to 1).
      See "SkewX".

   RQI

      (Accelerated) Rayleigh Quotient Iteration.

   JDQR

      Jacobi-Davidson with fixed number of inner steps. See
      "maxInnerIterations".

   JDQMR

      Jacobi-Davidson with adaptive stopping criterion for inner Quasi
      Minimum Residual (QMR). See "convTest".

   JDQMR_ETol

      JDQMR but QMR stops after residual norm reduces by a 0.1 factor.
      See "convTest".

   SUBSPACE_ITERATION

      Subspace iteration.

   LOBPCG_OrthoBasis

      LOBPCG, the basis size is set to the number of wanted
      eigenvalues "numEvals".

   LOBPCG_OrthoBasis_Window

      LOBPCG with sliding window of "maxBlockSize" < "numEvals".

void primme_initialize(primme_params *primme)

   Set PRIMME parameters structure to the default values.

   Parameters:
      * **primme** (*primme_params**) -- parameters structure.

int primme_set_method(primme_preset_method method, primme_params *primme)

   Set PRIMME parameters to one of the preset configurations.

   Parameters:
      * **method** (*primme_preset_method*) -- preset configuration.

      * **primme** (*primme_params**) -- parameters structure.

   Depending on the method some fields in "primme" are changed:

   * "DEFAULT_MIN_TIME" is like "JDQMR_ETol".

   * "DEFAULT_MIN_MATVECS" is like "GD_Olsen_plusK".

   * "DYNAMIC" is like "JDQMR_ETol" and changes
     "dynamicMethodSwitch" to 1.

   * "Arnoldi" changes:

     * "locking" = 0;

     * "maxPrevRetain" = 0;

     * "precondition" = 0;

     * "maxInnerIterations" = 0.

   * "GD" changes:

     * "locking" = 0;

     * "maxPrevRetain" = 0;

     * "robustShifts" = 1;

     * "maxInnerIterations" = 0;

     * "RightX" = 0;

     * "SkewX" = 0.

   * "GD_plusK"  changes:

     * "maxPrevRetain" to 2 if "maxBlockSize" is 1 and "numEvals" >
       1; otherwise set "maxPrevRetain" to "maxBlockSize".

     * "locking" = 0;

     * "maxInnerIterations" = 0;

     * "RightX" = 0;

     * "SkewX" = 0.

   * "GD_Olsen_plusK" is like "GD_plusK" and changes "RightX" to 1.

   * "JD_Olsen_plusK" is like "GD_plusK" and changes:

     * "robustShifts" = 1;

     * "RightX" to 1;

     * "SkewX" to 1;

   * "RQI" changes:

     * "locking" = 1;

     * "maxPrevRetain" = 0;

     * "robustShifts"  = 1;

     * "maxInnerIterations" = -1;

     * "LeftQ"   = 1;

     * "LeftX"   = 1;

     * "RightQ"  = 0;

     * "RightX"  = 1;

     * "SkewQ"   = 0;

     * "SkewX"   = 0;

     * "convTest" = "primme_full_LTolerance".

   * "JDQR" changes:

     * "locking"     = 1;

     * "maxPrevRetain"      = 1;

     * "robustShifts"       = 0;

     * "maxInnerIterations" = 10 if it is 0;

     * "LeftQ"   = 0;

     * "LeftX"   = 1;

     * "RightQ"  = 1;

     * "RightX"  = 1;

     * "SkewQ"   = 1;

     * "SkewX"   = 1;

     * "relTolBase" = 1.5;

     * "convTest" = "primme_full_LTolerance".

   * "JDQMR" changes:

     * "locking" = 0;

     * "maxPrevRetain" = 1 if it is 0

     * "maxInnerIterations" = -1;

     * "LeftQ"   = "precondition";

     * "LeftX"   = 1;

     * "RightQ"  = 0;

     * "RightX"  = 0;

     * "SkewQ"   = 0;

     * "SkewX"   = 1;

     * "convTest"  = "primme_adaptive".

   * "JDQMR_ETol" is like "JDQMR" and changes "convTest" =
     "primme_adaptive_ETolerance".

   * "SUBSPACE_ITERATION" changes:

     * "locking"    = 1;

     * "maxBasisSize" = "numEvals" *** 2;

     * "minRestartSize" = "numEvals";

     * "maxBlockSize" = "numEvals";

     * "scheme"  = "primme_thick";

     * "maxPrevRetain"      = 0;

     * "robustShifts"       = 0;

     * "maxInnerIterations" = 0;

     * "RightX"  = 1;

     * "SkewX"   = 0.

   * "LOBPCG_OrthoBasis" changes:

     * "locking"    = 0;

     * "maxBasisSize" = "numEvals" *** 3;

     * "minRestartSize" = "numEvals";

     * "maxBlockSize" = "numEvals";

     * "scheme"  = "primme_thick";

     * "maxPrevRetain"      = "numEvals";

     * "robustShifts"       = 0;

     * "maxInnerIterations" = 0;

     * "RightX"  = 1;

     * "SkewX"   = 0.

   * "LOBPCG_OrthoBasis_Window" changes:

     * "locking"    = 0;

     * "maxBasisSize" = "maxBlockSize" *** 3;

     * "minRestartSize" = "maxBlockSize";

     * "maxBlockSize" = "numEvals";

     * "scheme"  = "primme_thick";

     * "maxPrevRetain"      = "maxBlockSize";

     * "robustShifts"       = 0;

     * "maxInnerIterations" = 0;

     * "RightX"  = 1;

     * "SkewX"   = 0.

   Returns:
      if 0, successful; if negative, something went wrong.

void primme_Free(primme_params *primme)

   Free memory allocated by PRIMME.

   Parameters:
      * **primme** (*primme_params**) -- parameters structure.

int dprimme(double *evals, double *evecs, double *resNorms, primme_params *primme)

   Solve a real symmetric standard eigenproblems.

   Parameters:
      * **evals** (*double**) -- array at least of size "numEvals"
        to store the computed eigenvalues; all parallel calls return
        the same values in this array.

      * **resNorms** (*double**) -- array at least of size
        "numEvals" to store the residual norms of the computed
        eigenpairs; all parallel calls return the same value in this
        array.

      * **evecs** (*double**) -- array at least of size "nLocal"
        times "numEvals" to store columnwise the (local part of the)
        computed eigenvectors.

      * **primme** (*primme_params**) -- parameters structure.

   Returns:
      error indicator:

      * 0: success.

      * 1: reported only amount of required memory.

      * -1: failed in allocating int or real workspace.

      * -2: malloc failed in allocating a permutation integer array.

      * -3: main_iter() encountered problem; the calling stack of
        the functions where the error occurred was printed in
        "stderr".

      * -4: if argument "primme" is NULL.

      * -5: if "n" <= 0 or "nLocal" <= 0.

      * -6: if "numProcs" < 1.

      * -7: if "matrixMatvec" is NULL.

      * -8: if "applyPreconditioner" is NULL and "precondition" is
        not NULL.

      * -9: if "globalSumDouble" is NULL.

      * -10: if "numEvals" > "n".

      * -11: if "numEvals" < 0.

      * -12: if "eps" > 0 and "eps" < machine precision.

      * -13: if "target" is not properly defined.

      * -14: if "target" is one of "primme_closest_geq",
        "primme_closest_leq" or "primme_closest_abs" but
        "numTargetShifts" <= 0 (no shifts).

      * -15: if "target" is one of "primme_closest_geq",
        "primme_closest_leq" or "primme_closest_abs" but
        "targetShifts" is NULL  (no shifts array).

      * -16: if "numOrthoConst" < 0 or "numOrthoConst" >= "n". (no
        free dimensions left).

      * -17: if "maxBasisSize" < 2.

      * -18: if "minRestartSize" <= 0.

      * -19: if "maxBlockSize" <= 0.

      * -20: if "maxPrevRetain" < 0.

      * -21: if "scheme" is not one of *primme_thick* or
        *primme_dtr*.

      * -22: if "initSize" < 0.

      * -23: if not "locking" and "initSize" > "maxBasisSize".

      * -24: if "locking" and "initSize" > "numEvals".

      * -25: if "maxPrevRetain" + "minRestartSize" >=
        "maxBasisSize".

      * -26: if "minRestartSize" >= "n".

      * -27: if "printLevel" < 0 or "printLevel" > 5.

      * -28: if "convTest" is not one of "primme_full_LTolerance",
        "primme_decreasing_LTolerance", "primme_adaptive_ETolerance"
        or "primme_adaptive".

      * -29: if "convTest" == "primme_decreasing_LTolerance" and
        "relTolBase" <= 1.

      * -30: if "evals" is NULL, but not "evecs" and "resNorms".

      * -31: if "evecs" is NULL, but not "evals" and "resNorms".

      * -32: if "resNorms" is NULL, but not "evecs" and "evals".

zprimme(double *evals, Complex_Z *evecs, double *resNorms, primme_params *primme)

   Solve a Hermitian standard eigenproblems; see function "dprimme()".

primme_params

   Structure to set the problem matrices and eigensolver options.

   int n

      Dimension of the matrix.

   void (*matrixMatvec)(void *x, void *y, int *blockSize, primme_params *primme)

      Block matrix-multivector multiplication, y = A x in solving A x
      = \lambda x or A x = \lambda B x.

      Parameters:
         * **x** (*void**) --

         * **y** (*void**) -- one dimensional array containing the
           "blockSize" vectors packed one after the other (i.e., the
           leading dimension is the vector size), each of size
           "nLocal". The real type is "double*" and "Complex_Z*" when
           called from "dprimme()" and "zprimme()" respectively.

         * **blockSize** (*int**) -- number of vectors in x and y.

         * **primme** (*primme_params**) -- parameters structure.

      Note: Argument "blockSize" is passed by reference to make
        easier the interface to other languages (like Fortran).

   void (*applyPreconditioner)(void *x, void *y, int *blockSize, struct primme_params *primme)

      Block preconditioner-multivector application, y = M^{-1}x where
      M is usually an approximation of A - \sigma I or A - \sigma B
      for finding eigenvalues close to \sigma. The function follows
      the convention of "matrixMatvec".

   void (*massMatrixMatvec)(void *x, void *y, int *blockSize, struct primme_params *primme)

      Block matrix-multivector multiplication, y = B x in solving A x
      = \lambda B x. The function follows the convention of
      "matrixMatvec".

      Warning: Generalized eigenproblems not implemented in current
        version. This member is included for future compatibility.

   int numProcs

      Number of processes calling in parallel to "dprimme()" or
      "zprimme()". The default value is 1.

   int procID

      The identity of the process that is calling in parallel to
      "dprimme()" or "zprimme()". Only the process with id 0 prints
      information. The default value is 0.

   int nLocal

      Number of local rows on this process. The default value is "n"
      if "numProcs" is 1.

   void *commInfo

      A pointer to whatever parallel environment structures needed.
      For example, with MPI, it could be a pointer to the MPI
      communicator. PRIMME does not use this. It is available for
      possible use in user functions defined in "matrixMatvec",
      "applyPreconditioner", "massMatrixMatvec" and "globalSumDouble".
      The default values is NULL.

   void (*globalSumDouble)(double *sendBuf, double *recvBuf, int *count, primme_params *primme)

      Global sum reduction function.

      Parameters:
         * **sendBuf** (*double**) -- array of size count with the
           input local values.

         * **recvBuf** (*double**) -- array of size count with the
           output global values so that i-th element of recvBuf is the
           sum over all processes of the i-th element of sendBuf.

         * **count** (*int**) -- array size of sendBuf and recvBuf.

         * **primme** (*primme_params**) -- parameters structure.

      The default value is NULL if "numProcs" is 1. When MPI this can
      be a simply wrapper to MPI_Allreduce().

         void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count,
                                  primme_params *primme) {
            MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;
            MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM,
                          communicator);
         }

      Note: Argument "count" is passed by reference to make easier
        the interface to other languages (like Fortran).

      Note: The arguments "sendBuf" and "recvBuf" are always double
        arrays and "count" is always the number of double elements in
        both arrays, even for "zprimme()".

   int numEvals

      Number of eigenvalues wanted. The default value is 1.

   primme_target target

      Which eigenpairs to find:

      "primme_smallest"
         Smallest algebraic eigenvalues; "targetShifts" is ignored.

      "primme_largest"
         Largest algebraic eigenvalues; "targetShifts" is ignored.

      "primme_closest_geq"
         Closest to, but greater or equal than the shifts in
         "targetShifts".

      "primme_closest_leq"
         Closest to, but less or equal than the shifts in
         "targetShifts".

      "primme_closest_abs"
         Closest in absolute value to than the shifts in
         "targetShifts".

      The default value is "primme_smallest".

      Note: * If some shift is close to the lower (higher) end of
        the

          spectrum, use either "primme_closest_geq"
          ("primme_closest_leq") or "primme_closest_abs".

        * "primme_closest_leq" and "primme_closest_geq" are more
          efficient than "primme_closest_abs".

   int numTargetShifts

      Size of the array "targetShifts". Used only when "target" is
      "primme_closest_geq", "primme_closest_leq" or
      "primme_closest_abs". The default values is 0.

   double *targetShifts

      Array of shifts, at least of size "numTargetShifts". Used only
      when "target" is "primme_closest_geq", "primme_closest_leq" or
      "primme_closest_abs". The default values is NULL.

      The i-th shift (or the last one, if it is not given) is taken
      into account in finding the i-th eigenvalue.

      Note: For code efficiency and robustness, the shifts should be
        ordered. Order them in ascending (descending) order for shifts
        closer to the lower (higher) end of the spectrum.

      int printLevel

         The level of message reporting from the code. One of:

         * 0: silent.

         * 1: print some error messages when these occur.

         * 2: as 1, and info about targeted eigenpairs when they are
           marked as converged:

              #Converged $1 eval[ $2 ]= $3 norm $4 Mvecs $5 Time $7

           or locked:

              #Lock epair[ $1 ]= $3 norm $4 Mvecs $5 Time $7

         * 3: as 2, and info about targeted eigenpairs every outer
           iteration:

              OUT $6 conv $1 blk $8 MV $5 Sec $7 EV $3 |r| $4

           Also, if it is used the dynamic method, show JDQMR/GDk
           performance ratio and the current method in use.

         * 4: as 3, and info about targeted eigenpairs every inner
           iteration:

              INN MV $5 Sec $7 Eval $3 Lin|r| $9 EV|r| $4

         * 5: as 4, and verbose info about certain choices of the
           algorithm.

         Output key:

         * $1: Number of converged pairs up to now.

         * $2: The index of the pair currently converged.

         * $3: The eigenvalue.

         * $4: Its residual norm.

         * $5: The current number of matrix-vector products.

         * $6: The current number of outer iterations.

         * $7: The current elapsed time.

         * $8: Index within the block of the targeted pair .

         * $9: QMR norm of the linear system residual.

         In parallel programs, output is produced in call with
         "procID" 0 when "printLevel" is from 0 to 4. If "printLevel"
         is 5 output can be produced in any of the parallel calls.

      Note: Convergence history for plotting may be produced simply
        by:

           grep OUT outpufile | awk '{print $8" "$14}' > out
           grep INN outpufile | awk '{print $3" "$11}' > inn

        Then in Matlab:

           plot(out(:,1),out(:,2),'bo');hold; plot(inn(:,1),inn(:,2),'r');

        Or in gnuplot:

           plot 'out' w lp, 'inn' w lp

   double aNorm

      An estimate of norm of the matrix A that is used in the
      convergence criterion (see "eps"). If it is less or equal to 0,
      it is used the largest absolute Ritz value seen. And on return,
      it is replaced with that value.

      The default value is 0.

   double eps

      An eigenpairs is marked as converged when the 2-norm of the
      residual is less than "eps" times "aNorm". The residual vector
      is A x - \lambda x or A x - \lambda B x.

      The default value is 10^{-12}.

   FILE *outputFile

      Opened file to write down the output.

      The default value is the standard output.

   int dynamicMethodSwitch

      If this value is 1, it alternates dynamically between
      "DEFAULT_MIN_TIME" and "DEFAULT_MIN_MATVECS", trying to identify
      the fastest method.

      On exit, it holds a recommended method for future runs on this
      problem:

      * -1: use "DEFAULT_MIN_MATVECS" next time.

      * -2: use "DEFAULT_MIN_TIME" next time.

      * -3: close call, use "DYNAMIC" next time again.

      Note: Even for expert users we do not recommend setting
        "dynamicMethodSwitch" directly, but through
        "primme_set_method()".

      Note: The code obtains timings by the "gettimeofday" Unix
        utility. If a cheaper, more accurate timer is available,
        modify the "PRIMMESRC/COMMONSRC/wtime.c"

   int locking

      If set to 1, hard locking will be used (locking converged
      eigenvectors out of the search basis). Otherwise the code will
      try to use soft locking (à la ARPACK), when large enough
      "minRestartSize" is available.

      The default depends on the method and the value of some options.

   int initSize

      On input, the number of initial vector guesses provided in
      "evecs" argument in "dprimme()" or "zprimme()". On output, the
      number of converged eigenpairs. During execution, it holds the
      current number of converged eigenpairs. If in addition locking
      is used, these are accessible in "evals" and "evecs".

      The default value is 0.

   int numOrthoConst

      Number of external orthogonalization constraint vectors provided
      "evecs" argument in "dprimme()" or "zprimme()".

      Then eigenvectors are found orthogonal to those constraints
      (equivalent to solving the problem with (I-YY^*)A(I-YY^*) and
      (I-YY^*)B(I-YY^*) instead where Y are the given constraint
      vectors). This is a handy feature if some eigenvectors are
      already known, or for finding more eigenvalues after a call to
      "dprimme()" or "zprimme()".

      The default value is 0.

   int maxBasisSize

      The maximum basis size allowed in the main iteration. This has
      memory implications.

      The default depends on method.

      Note: For interior eigenvalues use a larger value than usual.

   int minRestartSize

      Maximum Ritz vectors kept after restarting the basis.

      The default depends on "maxBasisSize", "maxBlockSize" and
      method.

   int maxBlockSize

      The maximum block size the code will try to use.

      The user should set this based on the architecture specifics of
      the target computer, as well as any a priori knowledge of
      multiplicities. The code does *not* require to be greater than 1
      to find multiple eigenvalues. For some methods, keeping to 1
      yields the best overall performance.

      The default value is 1.

      Note: Inner iterations of QMR are not performed in a block
        fashion. Every correction equation from a block is solved
        independently.

   int maxMatvecs

      Maximum number of matrix vector multiplications (approximately
      equal to the number of preconditioning operations) that the code
      is allowed to perform before it exits.

      The default value is "INT_MAX".

   int maxOuterIterations

      Maximum number of outer iterations that the code is allowed to
      perform before it exits.

      The default value is "INT_MAX".

   int intWorkSize

      If "dprimme()" or "zprimme()" are called with all arguments as
      NULL but "primme_params" then it has the size *in bytes* of the
      integer workspace that is required.

      Otherwise if not 0, it is the size of the integer work array *in
      bytes* that the user provides in "intWork". If it is 0, the code
      will allocate the required space and should be freed by calling
      "primme_Free()".

      The default value is 0.

   long int realWorkSize

      If "dprimme()" or "zprimme()" are called with all arguments as
      NULL but "primme_params" then it has the size *in bytes* of the
      real workspace that is required.

      Otherwise if not 0, it is the size of the real work array *in
      bytes* that the user provides in "realWork". If it is 0, the
      code will allocate the required space and should be freed by
      calling "primme_Free()".

      The default value is 0.

   int *intWork

      Integer work array.

      If NULL, the code will allocate its own workspace. If the
      provided space is not enough, the code will free it and allocate
      a new space.

      On exit, the first element shows if a locking problem has
      occurred. Using locking for large "numEvals" may, in some rare
      cases, cause some pairs to be practically converged, in the
      sense that their components are in the basis of "evecs". If this
      is the case, a Rayleigh Ritz on returned "evecs" would provide
      the accurate eigenvectors (see [r4]).

      The default value is NULL.

   void *realWork

      Real work array.

      If NULL, the code will allocate its own workspace. If the
      provided space is not enough, the code will free it and allocate
      a new space.

      The default value is NULL.

   int iseed

      The "int iseed[4]" is an array with the seeds needed by the
      LAPACK dlarnv and zlarnv.

      The default value is an array with values 1, 2, 3 and 5.

   void *matrix

      This field may be used to pass any required information in the
      matrix-vector product "matrixMatvec".

      The default value is NULL.

   void *preconditioner

      This field may be used to pass any required information in the
      matrix-vector product "applyPreconditioner".

      The default value is NULL.

   double *ShiftsForPreconditioner

      Array of size "blockSize" provided during execution of
      :c:member:dprimme and :c:member:zprimme holding the shifts to be
      used (if needed) in the preconditioning operation.

      For example if the block size is 3, there will be an array of
      three shifts in "ShiftsForPreconditioner". Then the user can
      invert a shifted preconditioner for each of the block vectors
      (M-ShiftsForPreconditioner_i)^{-1} x_i. Classical Davidson
      (diagonal) preconditioning is an example of this.

   primme_restartscheme restartingParams.scheme

      Select a restarting strategy:

      * "primme_thick", Thick restarting. This is the most efficient
        and robust in the general case.

      * "primme_dtr", Dynamic thick restarting. Helpful without
        preconditioning but it is expensive to implement.

      The default value is "primme_thick".

   int restartingParams.maxPrevRetain

      Number of approximations from previous iteration to be retained
      after restart (see [r2]). The restart size is "minRestartSize"
      plus "maxPrevRetain".

      The default value is 1.

   int correctionParams.precondition

      Set to 1 to use preconditioning. Make sure "applyPreconditioner"
      is not NULL then!

      The default value is 0.

   int correctionParams.robustShifts

      Set to 1 to use robust shifting. It tries to avoid stagnation
      and misconvergence by providing as shifts in
      "ShiftsForPreconditioner" the Ritz values displaced by an
      approximation of the eigenvalue error.

      The default value depends on method.

   int correctionParams.maxInnerIterations

      Control the maximum number of inner QMR iterations:

      * 0:  no inner iterations;

      * >0: perform at most that number of inner iterations per
        outer step;

      * <0: perform at most the rest of the remaining matrix-vector
        products up to reach "maxMatvecs".

      The default value depends on method.

      See also "convTest".

   double correctionParams.relTolBase

      Parameter used when "convTest" is
      "primme_decreasing_LTolerance".

   primme_convergencetest correctionParams.convTest

      Set how to stop the inner QMR method:

      * "primme_full_LTolerance": stop by iterations only;

      * "primme_decreasing_LTolerance", stop when
        \text{relTolBase}^{-\text{outIts}} where outIts is the number
        of outer iterations and retTolBase is set in "relTolBase";
        This is a legacy option from classical JDQR and we recommend
        **strongly** against its use.

      * "primme_adaptive", stop when the estimated eigenvalue
        residual has reached the required tolerance (based on Notay's
        JDCG).

      * "primme_adaptive_ETolerance", as "primme_adaptive" but also
        stopping when the estimated eigenvalue residual has reduced 10
        times.

      The default value depends on method.

      Note: Avoid to set "maxInnerIterations" to -1 and "convTest"
        to "primme_full_LTolerance".

      See also "maxInnerIterations".

   int correctionParams.projectors.LeftQ

   int correctionParams.projectors.LeftX

   int correctionParams.projectors.RightQ

   int correctionParams.projectors.RightX

   int correctionParams.projectors.SkewQ

   int correctionParams.projectors.SkewX

      Control the projectors involved in the computation of the
      correction appended to the basis every (outer) iteration.

      Given the current selected Ritz value \Lambda and vectors X, the
      residual associated vectors R=AX-X\Lambda, the previous locked
      vectors Q and the preconditioner M^{-1}. The correction D
      appended to the basis in GD (when "maxInnerIterations" is 0) is:

      +----------+---------+--------------------------------------------------------------------+
      | RightX   | SkewX   | D                                                                  |
      +==========+=========+====================================================================+
      | 0        | 0       | M^{-1}R (Classic GD)                                               |
      +----------+---------+--------------------------------------------------------------------+
      | 1        | 0       | M^{-1}(R-\Delta X) (cheap Olsen's Method)                          |
      +----------+---------+--------------------------------------------------------------------+
      | 1        | 1       | (I- M^{-1}X(X^*M^{-1}X)^{-1}X^*)M^{-1}R (Olsen's Method)           |
      +----------+---------+--------------------------------------------------------------------+
      | 0        | 1       | error                                                              |
      +----------+---------+--------------------------------------------------------------------+

      Where \Delta is a diagonal matrix that \Delta_{i,i} holds an
      estimation of the error of the approximate eigenvalue
      \Lambda_{i,i}.

      The values of "RightQ", "SkewQ", "LeftX" and "LeftQ" are
      ignored.

      The correction D in JD (when "maxInnerIterations" isn't 0)
      results from solving:

         P_Q^l P_X^l (A-\sigma I) P_X^r P_Q^r M^{-1} D' = -R, \ \ \  D
         = P_X^r P_Q^l M^{-1}D'.

      For "LeftQ" (and similarly for "LeftX"):

      * 0: P_Q^l = I;

      * 1: P_Q^l = I - QQ^*.

      For "RightQ" and "SkewQ" (and similarly for "RightX" and
      "SkewX"):

      +----------+---------+---------------------------------+
      | RightQ   | SkewQ   | P_Q^r                           |
      +==========+=========+=================================+
      | 0        | 0       | I                               |
      +----------+---------+---------------------------------+
      | 1        | 0       | I - QQ^*                        |
      +----------+---------+---------------------------------+
      | 1        | 1       | I - KQ(Q^*KQ)^{-1}Q^*           |
      +----------+---------+---------------------------------+
      | 0        | 1       | error                           |
      +----------+---------+---------------------------------+

      The default value depends on method.

      See [r3] for a study about different projector configuration in
      JD.

   int stats.numOuterIterations

      Hold the number of outer iterations. The value is available
      during execution and at the end.

   int stats.numRestarts

      Hold the number of restarts during execution and at the end.

   int stats.numMatvecs

      Hold how many vectors the operator in "matrixMatvec" has been
      applied on. The value is available during execution and at the
      end.

   int stats.numPreconds

      Hold how many vectors the operator in "applyPreconditioner" has
      been applied on. The value is available during execution and at
      the end.

   int stats.elapsedTime

      Hold the wall clock time spent by the call to "dprimme()" or
      "zprimme()". The value is available at the end of the execution.

FORTRAN Library Interface
*************************

The next enumerations and functions are declared in "primme_f77.h".

ptr

   Fortran datatype with the same size as a pointer. Use "integer*4"
   when compiling in 32 bits and "integer*8" in 64 bits.

primme_initialize_f77(primme)

   Set PRIMME parameters structure to the default values.

   Parameters:
      * **primme** (*ptr*) -- (output) parameters structure.

primme_set_method_f77(method, primme, ierr)

   Set PRIMME parameters to one of the preset configurations.

   Parameters:
      * **method** (*integer*) --

        (input) preset configuration. One of:

        * "PRIMMEF77_DYNAMIC"

        * "PRIMMEF77_DEFAULT_MIN_TIME"

        * "PRIMMEF77_DEFAULT_MIN_MATVECS"

        * "PRIMMEF77_Arnoldi"

        * "PRIMMEF77_GD"

        * "PRIMMEF77_GD_plusK"

        * "PRIMMEF77_GD_Olsen_plusK"

        * "PRIMMEF77_JD_Olsen_plusK"

        * "PRIMMEF77_RQI"

        * "PRIMMEF77_JDQR"

        * "PRIMMEF77_JDQMR"

        * "PRIMMEF77_JDQMR_ETol"

        * "PRIMMEF77_SUBSPACE_ITERATION"

        * "PRIMMEF77_LOBPCG_OrthoBasis"

        * "PRIMMEF77_LOBPCG_OrthoBasis_Window"

        See "primme_preset_method".

      * **primme** (*ptr*) -- (input) parameters structure.

      * **ierr** (*integer*) -- (output) if 0, successful; if
        negative, something went wrong.

primme_Free_f77(primme)

   Free memory allocated by PRIMME.

   Parameters:
      * **primme** (*ptr*) -- parameters structure.

dprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a real symmetric standard eigenproblems.

   Parameters:
      * **evals(*)** (*double precision*) -- (output) array at least
        of size "numEvals" to store the computed eigenvalues; all
        parallel calls return the same value in this array.

      * **resNorms(*)** (*double precision*) -- (output) array at
        least of size "numEvals" to store the residual norms of the
        computed eigenpairs; all parallel calls return the same value
        in this array.

      * **evecs(*)** (*double precision*) -- (input/output) array at
        least of size "nLocal" times "numEvals" to store columnwise
        the (local part of the) computed eigenvectors.

      * **primme** (*ptr*) -- parameters structure.

      * **ierr** (*integer*) -- (output) error indicator; see the
        returned value of function "dprimme()".

zprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblems. The arguments have the
   same meaning like in function "dprimme_f77()".

   Parameters:
      * **evals(*)** (*double precision*) -- (output)

      * **resNorms(*)** (*double precision*) -- (output)

      * **evecs(*)** (*complex double precision*) -- (input/output)

      * **primme** (*ptr*) -- (input) parameters structure.

      * **ierr** (*integer*) -- (output) error indicator.

primmetop_set_member_f77(primme, label, value)

   Set a value in some field of the parameter structure.

   Parameters:
      * **primme** (*ptr*) -- (input) parameters structure.

      * **label** (*integer*) --

        field where to set value. One of:

        * "PRIMMEF77_n",                                     in
          field "primme_params.n".

        * "PRIMMEF77_matrixMatvec",                          in
          field "primme_params.matrixMatvec".

        * "PRIMMEF77_applyPreconditioner",                   in
          field "primme_params.applyPreconditioner".

        * "PRIMMEF77_numProcs",                              in
          field "primme_params.numProcs".

        * "PRIMMEF77_procID",                                in
          field "primme_params.procID".

        * "PRIMMEF77_commInfo",                              in
          field "primme_params.commInfo".

        * "PRIMMEF77_nLocal",                                in
          field "primme_params.nLocal".

        * "PRIMMEF77_globalSumDouble",                       in
          field "primme_params.globalSumDouble".

        * "PRIMMEF77_numEvals",                              in
          field "primme_params.numEvals".

        * "PRIMMEF77_target",                                in
          field "primme_params.target".

        * "PRIMMEF77_numTargetShifts",                       in
          field "primme_params.numTargetShifts".

        * "PRIMMEF77_targetShifts",                          in
          field "primme_params.targetShifts".

        * "PRIMMEF77_locking",                               in
          field "primme_params.locking".

        * "PRIMMEF77_initSize",                              in
          field "primme_params.initSize".

        * "PRIMMEF77_numOrthoConst",                         in
          field "primme_params.numOrthoConst".

        * "PRIMMEF77_maxBasisSize",                          in
          field "primme_params.maxBasisSize".

        * "PRIMMEF77_minRestartSize",                        in
          field "primme_params.minRestartSize".

        * "PRIMMEF77_maxBlockSize",                          in
          field "primme_params.maxBlockSize".

        * "PRIMMEF77_maxMatvecs",                            in
          field "primme_params.maxMatvecs".

        * "PRIMMEF77_maxOuterIterations",                    in
          field "primme_params.maxOuterIterations".

        * "PRIMMEF77_intWorkSize",                           in
          field "primme_params.intWorkSize".

        * "PRIMMEF77_realWorkSize",                          in
          field "primme_params.realWorkSize".

        * "PRIMMEF77_iseed",                                 in
          field "primme_params.iseed".

        * "PRIMMEF77_intWork",                               in
          field "primme_params.intWork".

        * "PRIMMEF77_realWork",                              in
          field "primme_params.realWork".

        * "PRIMMEF77_aNorm",                                 in
          field "primme_params.aNorm".

        * "PRIMMEF77_eps",                                   in
          field "primme_params.eps".

        * "PRIMMEF77_printLevel",                            in
          field "primme_params.printLevel".

        * "PRIMMEF77_outputFile",                            in
          field "primme_params.outputFile".

        * "PRIMMEF77_matrix",                                in
          field "primme_params.matrix".

        * "PRIMMEF77_preconditioner",                        in
          field "primme_params.preconditioner".

        * "PRIMMEF77_restartingParams_scheme",               in
          field "primme_params.restartingParams.scheme".

        * "PRIMMEF77_restartingParams_maxPrevRetain",        in
          field "primme_params.restartingParams.maxPrevRetain".

        * "PRIMMEF77_correctionParams_precondition",         in
          field "primme_params.correctionParams.precondition".

        * "PRIMMEF77_correctionParams_robustShifts",         in
          field "primme_params.correctionParams.robustShifts".

        * "PRIMMEF77_correctionParams_maxInnerIterations",   in
          field "primme_params.correctionParams.maxInnerIterations".

        * "PRIMMEF77_correctionParams_projectors_LeftQ",     in
          field "primme_params.correctionParams.projectors.LeftQ".

        * "PRIMMEF77_correctionParams_projectors_LeftX",     in
          field "primme_params.correctionParams.projectors.LeftX".

        * "PRIMMEF77_correctionParams_projectors_RightQ",    in
          field "primme_params.correctionParams.projectors.RightQ".

        * "PRIMMEF77_correctionParams_projectors_RightX",    in
          field "primme_params.correctionParams.projectors.RightX".

        * "PRIMMEF77_correctionParams_projectors_SkewQ",     in
          field "primme_params.correctionParams.projectors.SkewQ".

        * "PRIMMEF77_correctionParams_projectors_SkewX",     in
          field "primme_params.correctionParams.projectors.SkewX".

        * "PRIMMEF77_correctionParams_convTest",             in
          field "primme_params.correctionParams.convTest".

        * "PRIMMEF77_correctionParams_relTolBase",           in
          field "primme_params.correctionParams.relTolBase".

        * "PRIMMEF77_stats_numOuterIterations",              in
          field "primme_params.stats.numOuterIterations".

        * "PRIMMEF77_stats_numRestarts",                     in
          field "primme_params.stats.numRestarts".

        * "PRIMMEF77_stats_numMatvecs",                      in
          field "primme_params.stats.numMatvecs".

        * "PRIMMEF77_stats_numPreconds",                     in
          field "primme_params.stats.numPreconds".

        * "PRIMMEF77_stats_elapsedTime",                     in
          field "primme_params.stats.elapsedTime".

        * "PRIMMEF77_dynamicMethodSwitch",                   in
          field "primme_params.dynamicMethodSwitch".

        * "PRIMMEF77_massMatrixMatvec",                      in
          field "primme_params.massMatrixMatvec".

      * **value** -- (input) value to set.

primmetop_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   Parameters:
      * **primme** (*ptr*) -- (input) parameters structure.

      * **label** (*integer*) -- (input) field where to get value.
        One of the detailed in function "primmetop_set_member_f77()".

      * **value** -- (output) value of the field.

primmetop_get_prec_shift_f77(primme, index, value)

   Get the value in some position of the array
   "ShiftsForPreconditioner".

   Parameters:
      * **primme** (*ptr*) -- (input) parameters structure.

      * **index** (*integer*) -- (input) position of the array; the
        first position is 1.

      * **value** -- (output) value of the array at that position.

primme_set_member_f77(primme, label, value)

   Set a value in some field of the parameter structure.

   Parameters:
      * **primme** (*ptr*) -- (input) parameters structure.

      * **label** (*integer*) --

        field where to set value. One of:

        * "PRIMMEF77_n",                                     in
          field "primme_params.n".

        * "PRIMMEF77_matrixMatvec",                          in
          field "primme_params.matrixMatvec".

        * "PRIMMEF77_applyPreconditioner",                   in
          field "primme_params.applyPreconditioner".

        * "PRIMMEF77_numProcs",                              in
          field "primme_params.numProcs".

        * "PRIMMEF77_procID",                                in
          field "primme_params.procID".

        * "PRIMMEF77_commInfo",                              in
          field "primme_params.commInfo".

        * "PRIMMEF77_nLocal",                                in
          field "primme_params.nLocal".

        * "PRIMMEF77_globalSumDouble",                       in
          field "primme_params.globalSumDouble".

        * "PRIMMEF77_numEvals",                              in
          field "primme_params.numEvals".

        * "PRIMMEF77_target",                                in
          field "primme_params.target".

        * "PRIMMEF77_numTargetShifts",                       in
          field "primme_params.numTargetShifts".

        * "PRIMMEF77_targetShifts",                          in
          field "primme_params.targetShifts".

        * "PRIMMEF77_locking",                               in
          field "primme_params.locking".

        * "PRIMMEF77_initSize",                              in
          field "primme_params.initSize".

        * "PRIMMEF77_numOrthoConst",                         in
          field "primme_params.numOrthoConst".

        * "PRIMMEF77_maxBasisSize",                          in
          field "primme_params.maxBasisSize".

        * "PRIMMEF77_minRestartSize",                        in
          field "primme_params.minRestartSize".

        * "PRIMMEF77_maxBlockSize",                          in
          field "primme_params.maxBlockSize".

        * "PRIMMEF77_maxMatvecs",                            in
          field "primme_params.maxMatvecs".

        * "PRIMMEF77_maxOuterIterations",                    in
          field "primme_params.maxOuterIterations".

        * "PRIMMEF77_intWorkSize",                           in
          field "primme_params.intWorkSize".

        * "PRIMMEF77_realWorkSize",                          in
          field "primme_params.realWorkSize".

        * "PRIMMEF77_iseed",                                 in
          field "primme_params.iseed".

        * "PRIMMEF77_intWork",                               in
          field "primme_params.intWork".

        * "PRIMMEF77_realWork",                              in
          field "primme_params.realWork".

        * "PRIMMEF77_aNorm",                                 in
          field "primme_params.aNorm".

        * "PRIMMEF77_eps",                                   in
          field "primme_params.eps".

        * "PRIMMEF77_printLevel",                            in
          field "primme_params.printLevel".

        * "PRIMMEF77_outputFile",                            in
          field "primme_params.outputFile".

        * "PRIMMEF77_matrix",                                in
          field "primme_params.matrix".

        * "PRIMMEF77_preconditioner",                        in
          field "primme_params.preconditioner".

        * "PRIMMEF77_restartingParams_scheme",               in
          field "primme_params.restartingParams.scheme".

        * "PRIMMEF77_restartingParams_maxPrevRetain",        in
          field "primme_params.restartingParams.maxPrevRetain".

        * "PRIMMEF77_correctionParams_precondition",         in
          field "primme_params.correctionParams.precondition".

        * "PRIMMEF77_correctionParams_robustShifts",         in
          field "primme_params.correctionParams.robustShifts".

        * "PRIMMEF77_correctionParams_maxInnerIterations",   in
          field "primme_params.correctionParams.maxInnerIterations".

        * "PRIMMEF77_correctionParams_projectors_LeftQ",     in
          field "primme_params.correctionParams.projectors.LeftQ".

        * "PRIMMEF77_correctionParams_projectors_LeftX",     in
          field "primme_params.correctionParams.projectors.LeftX".

        * "PRIMMEF77_correctionParams_projectors_RightQ",    in
          field "primme_params.correctionParams.projectors.RightQ".

        * "PRIMMEF77_correctionParams_projectors_RightX",    in
          field "primme_params.correctionParams.projectors.RightX".

        * "PRIMMEF77_correctionParams_projectors_SkewQ",     in
          field "primme_params.correctionParams.projectors.SkewQ".

        * "PRIMMEF77_correctionParams_projectors_SkewX",     in
          field "primme_params.correctionParams.projectors.SkewX".

        * "PRIMMEF77_correctionParams_convTest",             in
          field "primme_params.correctionParams.convTest".

        * "PRIMMEF77_correctionParams_relTolBase",           in
          field "primme_params.correctionParams.relTolBase".

        * "PRIMMEF77_stats_numOuterIterations",              in
          field "primme_params.stats.numOuterIterations".

        * "PRIMMEF77_stats_numRestarts",                     in
          field "primme_params.stats.numRestarts".

        * "PRIMMEF77_stats_numMatvecs",                      in
          field "primme_params.stats.numMatvecs".

        * "PRIMMEF77_stats_numPreconds",                     in
          field "primme_params.stats.numPreconds".

        * "PRIMMEF77_stats_elapsedTime",                     in
          field "primme_params.stats.elapsedTime".

        * "PRIMMEF77_dynamicMethodSwitch",                   in
          field "primme_params.dynamicMethodSwitch".

        * "PRIMMEF77_massMatrixMatvec",                      in
          field "primme_params.massMatrixMatvec".

      * **value** -- (input) value to set.

   Note: Use this function exclusively inside the function
     "matrixMatvec", "massMatrixMatvec", or "applyPreconditioner".
     Otherwise use the function "primmetop_set_member_f77()".

primme_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   Parameters:
      * **primme** (*ptr*) -- (input) parameters structure.

      * **label** (*integer*) -- (input) field where to get value.
        One of the detailed in function "primmetop_set_member_f77()".

      * **value** -- (output) value of the field.

   Note: Use this function exclusively inside the function
     "matrixMatvec", "massMatrixMatvec", or "applyPreconditioner".
     Otherwise use the function "primmetop_get_member_f77()".

primme_get_prec_shift_f77(primme, index, value)

   Get the value in some position of the array
   "ShiftsForPreconditioner".

   Parameters:
      * **primme** (*ptr*) -- (input) parameters structure.

      * **index** (*integer*) -- (input) position of the array; the
        first position is 1.

      * **value** -- (output) value of the array at that position.

   Note: Use this function exclusively inside the function
     "matrixMatvec", "massMatrixMatvec", or "applyPreconditioner".
     Otherwise use the function "primmetop_get_prec_shift_f77()".
