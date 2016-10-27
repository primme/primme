Note: For hyperlinked html and pdf versions of this document see
  directory "doc".


PRIMME Documentation
********************

Table Of Contents:

* PRIMME: PReconditioned Iterative MultiMethod Eigensolver

  * Incompatibilities

  * Changelog

  * License Information

  * Citing the code

  * Contact Information

  * Directory Structure

  * Making and Linking

  * Tested Systems

  * Main Contributors

* Eigenvalue Problems

  * C Library Interface

  * FORTRAN Library Interface

  * Python Interface

  * MATLAB Interface

  * Appendix

* Singular Value Problems

  * C Library Interface

  * FORTRAN Library Interface

  * Python Interface

  * Appendix

PRIMME: PReconditioned Iterative MultiMethod Eigensolver
********************************************************

PRIMME, pronounced as *prime*, finds a number of eigenvalues and their
corresponding eigenvectors of a real symmetric, or Hermitian matrix.
Also singular values and vectors can be computed. Largest, smallest
and interior eigenvalues and singular values are supported.
Preconditioning can be used to accelerate convergence. PRIMME is
written in C99, but complete interfaces are provided for Fortran 77,
MATLAB and Python.


Incompatibilities
=================

From PRIMME 1.x to 2.0:

* Prototype of callbacks has changed: "matrixMatvec",
  "applyPreconditioner", "massMatrixMatvec" and "globalSumReal".

* The next parameters are "PRIMME_INT": "n", "nLocal", "maxMatvecs",
  "iseed", "numOuterIterations", "numRestarts", "numMatvecs" and
  "numMatvecs"; use the macro "PRIMME_INT_P" to print the values.

* Rename the values of the enum "primme_preset_method".

* Rename "primme_Free" to "primme_free()".

* Integer parameters in Fortran functions are of the same size as
  "PRIMME_INT", which is "integer*8" by default.

* Extra parameter in many Fortran functions to return the error
  code.

* Removed "primme_display_stats_f77".


Changelog
=========

Changes in PRIMME 2.1 (released on XXX):

* Support Octave.

Changes in PRIMME 2.0 (released on September 19, 2016):

* Changed license to BSD 3-clause.

* New support for singular value problems; see "dprimme_svds()".

* New support for "float" and "complex float" arithmetic.

* Support for problem dimensions larger than 2^31, without requiring
  BLAS and LAPACK compiled with 64-bits integers.

* Improve robustness and performance for interior problems;
  implemented advanced refined and harmonic-Ritz extractions.

* Python interface compatible with NumPy and SciPy Library.

* Added parameter to indicate the leading dimension of the
  input/output matrices and to return an error code in callbacks
  "matrixMatvec", "applyPreconditioner", "massMatrixMatvec" and
  "globalSumReal".

* Changed to type "PRIMME_INT" the options "n", "nLocal",
  "maxMatvecs" and "iseed", and the stats counters
  "numOuterIterations", "numRestarts", "numMatvecs", "numPreconds".
  Also changed "realWorkSize" to "size_t". Fortran interface functions
  will expect an "interger" of size compatible with "PRIMME_INT" for
  all parameters with integer type: "int", "PRIMME_INT" and "size_t";
  see also parameter "value" in functions
  "primmetop_set_member_f77()", "primmetop_get_member_f77()",
  "primme_set_member_f77()" and "primme_get_member_f77()".

* Added parameter to return an error code in Fortran interface
  functions: "primmetop_set_member_f77()",
  "primmetop_get_member_f77()", "primme_set_member_f77()" and
  "primme_get_member_f77()".

* Added leading dimension for "evecs" "ldevecs" and preferred
  leading dimension for the operators "ldOPs", such as "matrixMatvec".

* Optional user-defined convergence function, "convTestFun".

* Prefixed methods with "PRIMME_". Rename Fortran constants from
  "PRIMMEF77_" to "PRIMME_".

* Removed "primme_display_stats_f77".

Changes in PRIMME 1.2.2 (released on October 13, 2015):

* Fixed wrong symbols in "libdprimme.a" and "libzprimme.a".

* "primme_set_method()" sets "PRIMME_JDQMR" instead of
  "PRIMME_JDQMR_ETol" for preset methods "PRIMME_DEFAULT_MIN_TIME" and
  "PRIMME_DYNAMIC" when seeking interior values.

* Fixed compilation of driver with a PETSc installation without
  HYPRE.

* Included the content of the environment variable "INCLUDE" for
  compiling the driver.

Changes in PRIMME 1.2.1 (released on September 7, 2015):

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

* Fixed a few minor bugs and improved documentation (especially the
  F77 interface).

* Using Sphinx to manage documentation.

Changes in PRIMME 1.2 (released on December 21, 2014):

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
  to 4096^3  (140 trillion) processes.

* For the "PRIMME_DYNAMIC" method, fixed issues with initialization
  and synchronization decisions across multiple processes.

* Fixed uncommon library interface bugs, coordinated better setting
  the method and the user setting of parameters, and improved the
  interface in the sample programs and makefiles.

* Other performance and documentation improvements.


License Information
===================

PRIMME is licensed under the 3-clause license BSD. Python and Matlab
interfaces have BSD-compatible licenses. Source code under
file:*tests* is compatible with LGPLv3. Details can be taken from
COPYING.txt.


Citing the code
===============

Please cite:

[r1] A. Stathopoulos and J. R. McCombs PRIMME: *PReconditioned
     Iterative MultiMethod Eigensolver: Methods and software
     description*, ACM Transaction on Mathematical Software Vol. 37,
     No. 2, (2010), 21:1-21:30.

[r6] L. Wu, E. Romero and A. Stathopoulos, *PRIMME_SVDS: A High-
     Performance Preconditioned SVD Solver for Accurate Large-Scale
     Computations*, arXiv:1607.01404

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

[r7] L. Wu and A. Stathopoulos, *A Preconditioned Hybrid SVD
     Method for Computing Accurately Singular Triplets of Large
     Matrices*, SIAM J. Sci. Comput. 37-5(2015), pp. S365-S388.


Contact Information
===================

For reporting bugs or questions about functionality contact Andreas
Stathopoulos by email, *andreas* at *cs.wm.edu*. See further
information in the webpage http://www.cs.wm.edu/~andreas/software and
on github.


Directory Structure
===================

The next directories and files should be available:

* "COPYING.txt", license;

* "Make_flags",  flags to be used by makefiles to compile library
  and tests;

* "Link_flags",  flags needed in making and linking the test
  programs;

* "include/",    directory with headers files;

* "src/",        directory with the source code for "libprimme":

     * "include/",   common headers;

     * "eigs/",      eigenvalue interface and implementation;

     * "svds/",      singular value interface and implementation;

     * "tools/",     tools used to generated some headers;

* "Matlab/",       Matlab interface;

* "PYTHON/",       Python interface;

* "examples/",     sample programs in C, C++ and F77, both
  sequential and parallel;

* "tests/",        drivers for testing purpose and test cases;

* "lib/libprimme.a",   the PRIMME library (to be made);

* "makefile"       main make file;

* "readme.txt"     text version of the documentation;

* "doc/"           directory with the HTML and PDF versions of the
  documentation.


Making and Linking
==================

"Make_flags" has the flags and compilers used to make "libprimme.a":

* *CC*, compiler program such as "gcc", "clang" or "icc".

* *CFLAGS*, compiler options such as "-g" or "-O3" and macro
  definitions

     like the ones described next.

Compiler flags for the BLAS and LAPACK libraries:

* "-DF77UNDERSCORE", if Fortran appends an underscore to function
  names (usually it does).

* "-DPRIMME_BLASINT_SIZE=64", if the library integers are 64-bit
  integer ("kind=8") type, aka ILP64 interface; usually integers are
  32-bits even in 64-bit architectures (aka LP64 interface).

By default PRIMME sets the integer type for matrix dimensions and
counters ("PRIMME_INT") to 64 bits integer "int64_t". This can be
changed by setting the macro "PRIMME_INT_SIZE" to one of the following
values:

* "0": use the regular "int" of your compiler.

* "32": use C99 "int32_t".

* "64": use C99 "int64_t".

Note: When "-DPRIMME_BLASINT_SIZE=64" is set the code uses the type
  "int64_t" supported by the C99 standard. In case the compiler
  doesn't honor the standard, you can set the corresponding type name
  supported, for instance "-DPRIMME_BLASINT_SIZE=__int64".

After customizing "Make_flags", type this to generate "libprimme.a":

   make lib

Making can be also done at the command line:

   make lib CC=clang CFLAGS='-O3'

"Link_flags" has the flags for linking with external libraries and
making the executables located in "examples" and "tests":

* *LDFLAGS*, linker flags such as "-framework Accelerate".

* *LIBS*, flags to link with libraries (BLAS and LAPACK are
  required), such as "-lprimme -llapack -lblas -lgfortran -lm".

After that, type this to compile and execute a simple test:

   $ make test
   ...
   Test passed!
   ...
   Test passed!

In case of linking problems check flags in *LDFLAGS* and *LIBS* and
consider to add/remove "-DF77UNDERSCORE" from *CFLAGS*. If the
execution fails consider to add/remove "-DPRIMME_BLASINT_SIZE=64" from
*CFLAGS*.

Full description of actions that *make* can take:

* *make lib*, builds the static library "libprimme.a".

* *make solib*, builds the shared library "libprimme.so".

* *make test*, build and execute simple examples.

* *make clean*, removes all "*.o", "a.out", and core files from
  "src".


Considerations using an IDE
---------------------------

PRIMME can be built in other environments such as Anjuta, Eclipse,
KDevelop, Qt Creator, Visual Studio and XCode. To build the PRIMME
library do the following:

1. Create a new project and include the source files under the
   directory "src".

2. Add the directories "include" and "src/include" as include
   directories.

To build an example code using PRIMME make sure:

* to add a reference for PRIMME, BLAS and LAPACK libraries;

* to add the directory "include" as an include directory.


Tested Systems
==============

PRIMME is primary developed with GNU gcc, g++ and gfortran (versions
4.8 and later). Many users have reported builds on several other
platforms/compilers:

* SUSE 13.1 & 13.2

* CentOS 6.6

* Ubuntu 14.04

* MacOS X 10.9 & 10.10

* Cygwin & MinGW

* Cray XC30

* SunOS 5.9, quad processor Sun-Fire-280R, and several other
  UltraSparcs

* AIX 5.2 IBM SP POWER 3+, 16-way SMP, 375 MHz nodes (seaborg at
  nersc.gov)


Main Contributors
=================

* James R. McCombs

* Eloy Romero Alcalde

* Andreas Stathopoulos

* Lingfei Wu

Eigenvalue Problems
*******************

* C Library Interface

* FORTRAN Library Interface

* Python Interface

* MATLAB Interface

* Appendix

C Library Interface
*******************

The PRIMME interface is composed of the following functions. To solve
real symmetric and Hermitian standard eigenproblems call respectively:

   int sprimme(float *evals, float *evecs, float *resNorms,
               primme_params *primme);

   int cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms,
               primme_params *primme);

   int dprimme(double *evals, double *evecs, double *resNorms,
               primme_params *primme);

   int zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms,
               primme_params *primme);

Other useful functions:

   void primme_initialize(primme_params *primme);
   int primme_set_method(primme_preset_method method,
                                        primme_params *params);
   void primme_display_params(primme_params primme);
   void primme_free(primme_params primme);

PRIMME stores its data on the structure "primme_params". See
*Parameters Guide* for an introduction about its fields.


Running
=======

To use PRIMME, follow these basic steps.

1. Include:

      #include "primme.h"   /* header file is required to run primme */

2. Initialize a PRIMME parameters structure for default settings:

      primme_params primme;

      primme_initialize(&primme);

3. Set problem parameters (see also *Parameters Guide*), and,
   optionally, set one of the "preset methods":

      primme.matrixMatvec = LaplacianMatrixMatvec; /* MV product */
      primme.n = 100;                   /* set problem dimension */
      primme.numEvals = 10;       /* Number of wanted eigenpairs */
      ret = primme_set_method(method, &primme);
      ...

4. Then to solve real symmetric standard eigenproblems call:

      ret = dprimme(evals, evecs, resNorms, &primme);

   The previous is the double precision call. There is available calls
   for complex double, single and complex single; check it out
   "zprimme()", "sprimme()" and "cprimme()".

   The call arguments are:

   * *evals*, array to return the found eigenvalues;

   * *evecs*, array to return the found eigenvectors;

   * *resNorms*, array to return the residual norms of the found
     eigenpairs; and

   * *ret*, returned error code.

5. To free the work arrays in PRIMME:

      primme_free(&primme);


Parameters Guide
================

PRIMME stores the data on the structure "primme_params", which has the
next fields:

   /* Basic */
   PRIMME_INT n;                                      // matrix dimension
   void (*matrixMatvec)(...);             // matrix-vector product
   int numEvals;                    // how many eigenpairs to find
   primme_target target;              // which eigenvalues to find
   int numTargetShifts;       // for targeting interior eigenpairs
   double *targetShifts;
   double eps;            // tolerance of the converged eigenpairs

   /* For parallel programs */
   int numProcs;           // number of processes
   int procID;             // rank of this process
   PRIMME_INT nLocal;      // number of rows stored in this process
   void (*globalSumReal)(...); // sum reduction among processes

   /* Accelerate the convergence */
   void (*applyPreconditioner)(...);     // precond-vector product
   int initSize;       // initial vectors as approximate solutions
   int maxBasisSize;
   int minRestartSize;
   int maxBlockSize;

   /* User data */
   void *commInfo;
   void *matrix;
   void *preconditioner;

   /* Advanced options */
   PRIMME_INT ldevecs; // leading dimension of the evecs
   int numOrthoConst; // orthogonal constrains to the eigenvectors
   int dynamicMethodSwitch;
   int locking;
   PRIMME_INT maxMatvecs;
   PRIMME_INT maxOuterIterations;
   int intWorkSize;
   size_t realWorkSize;
   PRIMME_INT iseed[4];
   int *intWork;
   void *realWork;
   double aNorm;
   int printLevel;
   FILE *outputFile;
   double *ShiftsForPreconditioner;
   primme_init initBasisMode;
   struct projection_params projectionParams;
   struct restarting_params restartingParams;
   struct correction_params correctionParams;
   struct primme_stats stats;
   void (*convTestFun)(...);
   PRIMME_INT ldOPS;   // leading dimension to use in matrixMatvec...

PRIMME requires the user to set at least the dimension of the matrix
("n") and the matrix-vector product ("matrixMatvec"), as they define
the problem to be solved. For parallel programs, "nLocal", "procID"
and "globalSumReal" are also required.

In addition, most users would want to specify how many eigenpairs to
find, and provide a preconditioner (if available).

It is useful to have set all these before calling
"primme_set_method()". Also, if users have a preference on
"maxBasisSize", "maxBlockSize", etc, they should also provide them
into "primme_params" prior to the "primme_set_method()" call. This
helps "primme_set_method()" make the right choice on other parameters.
It is sometimes useful to check the actual parameters that PRIMME is
going to use (before calling it) or used (on return) by printing them
with "primme_display_params()".


Interface Description
=====================

The next enumerations and functions are declared in "primme.h".


sprimme
-------

int sprimme(float *evals, float *evecs, float *resNorms, primme_params *primme)

   Solve a real symmetric standard eigenproblem.

   Parameters:
      * **evals** -- array at least of size "numEvals" to store the
        computed eigenvalues; all processes in a parallel run return
        this local array with the same values.

      * **resNorms** -- array at least of size "numEvals" to store
        the residual norms of the computed eigenpairs; all processes
        in parallel run return this local array with the same values.

      * **evecs** -- array at least of size "nLocal" times
        "numEvals" to store columnwise the (local part of the)
        computed eigenvectors.

      * **primme** -- parameters structure.

   Returns:
      error indicator; see *Error Codes*.


dprimme
-------

int dprimme(double *evals, double *evecs, double *resNorms, primme_params *primme)

   Solve a real symmetric standard eigenproblem.

   Parameters:
      * **evals** -- array at least of size "numEvals" to store the
        computed eigenvalues; all processes in a parallel run return
        this local array with the same values.

      * **resNorms** -- array at least of size "numEvals" to store
        the residual norms of the computed eigenpairs; all processes
        in parallel run return this local array with the same values.

      * **evecs** -- array at least of size "nLocal" times
        "numEvals" to store columnwise the (local part of the)
        computed eigenvectors.

      * **primme** -- parameters structure.

   Returns:
      error indicator; see *Error Codes*.


cprimme
-------

int cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, primme_params *primme)

   Solve a Hermitian standard eigenproblem; see function "sprimme()".


zprimme
-------

int zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, primme_params *primme)

   Solve a Hermitian standard eigenproblem; see function "dprimme()".


primme_initialize
-----------------

void primme_initialize(primme_params *primme)

   Set PRIMME parameters structure to the default values.

   Parameters:
      * **primme** -- parameters structure.


primme_set_method
-----------------

int primme_set_method(primme_preset_method method, primme_params *primme)

   Set PRIMME parameters to one of the preset configurations.

   Parameters:
      * **method** --

        preset configuration; one of

           "PRIMME_DYNAMIC"
           "PRIMME_DEFAULT_MIN_TIME"
           "PRIMME_DEFAULT_MIN_MATVECS"
           "PRIMME_Arnoldi"
           "PRIMME_GD"
           "PRIMME_GD_plusK"
           "PRIMME_GD_Olsen_plusK"
           "PRIMME_JD_Olsen_plusK"
           "PRIMME_RQI"
           "PRIMME_JDQR"
           "PRIMME_JDQMR"
           "PRIMME_JDQMR_ETol"
           "PRIMME_SUBSPACE_ITERATION"
           "PRIMME_LOBPCG_OrthoBasis"
           "PRIMME_LOBPCG_OrthoBasis_Window"

      * **primme** -- parameters structure.

   See also *Preset Methods*.


primme_display_params
---------------------

void primme_display_params(primme_params primme)

   Display all printable settings of "primme" into the file descriptor
   "outputFile".

   Parameters:
      * **primme** -- parameters structure.


primme_free
-----------

void primme_free(primme_params *primme)

   Free memory allocated by PRIMME.

   Parameters:
      * **primme** -- parameters structure.

FORTRAN Library Interface
*************************

The next enumerations and functions are declared in "primme_f77.h".

ptr

   Fortran datatype with the same size as a pointer. Use "integer*4"
   when compiling in 32 bits and "integer*8" in 64 bits.


primme_initialize_f77
=====================

primme_initialize_f77(primme)

   Set PRIMME parameters structure to the default values.

   Parameters:
      * **primme** (*ptr*) -- (output) parameters structure.


primme_set_method_f77
=====================

primme_set_method_f77(method, primme, ierr)

   Set PRIMME parameters to one of the preset configurations.

   Parameters:
      * **method** (*integer*) --

        (input) preset configuration. One of:

           "PRIMME_DYNAMIC"
           "PRIMME_DEFAULT_MIN_TIME"
           "PRIMME_DEFAULT_MIN_MATVECS"
           "PRIMME_Arnoldi"
           "PRIMME_GD"
           "PRIMME_GD_plusK"
           "PRIMME_GD_Olsen_plusK"
           "PRIMME_JD_Olsen_plusK"
           "PRIMME_RQI"
           "PRIMME_JDQR"
           "PRIMME_JDQMR"
           "PRIMME_JDQMR_ETol"
           "PRIMME_SUBSPACE_ITERATION"
           "PRIMME_LOBPCG_OrthoBasis"
           "PRIMME_LOBPCG_OrthoBasis_Window"

        See "primme_preset_method".

      * **primme** (*ptr*) -- (input) parameters structure.

      * **ierr** (*integer*) -- (output) if 0, successful; if
        negative, something went wrong.


primme_free_f77
===============

primme_free_f77(primme)

   Free memory allocated by PRIMME and delete all values set.

   Parameters:
      * **primme** (*ptr*) -- (input/output) parameters structure.


sprimme_f77
===========

sprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a real symmetric standard eigenproblem using single
   precision.

   Parameters:
      * **evals(*)** (*real*) -- (output) array at least of size
        "numEvals" to store the computed eigenvalues; all parallel
        calls return the same value in this array.

      * **resNorms(*)** (*real*) -- (output) array at least of size
        "numEvals" to store the residual norms of the computed
        eigenpairs; all parallel calls return the same value in this
        array.

      * **evecs(*)** (*real*) -- (input/output) array at least of
        size "nLocal" times "numEvals" to store columnwise the (local
        part of the) computed eigenvectors.

      * **primme** (*ptr*) -- parameters structure.

      * **ierr** (*integer*) -- (output) error indicator; see *Error
        Codes*.


cprimme_f77
===========

cprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblem. The arguments have the
   same meaning as in function "sprimme_f77()".

   Parameters:
      * **evals(*)** (*real*) -- (output)

      * **resNorms(*)** (*real*) -- (output)

      * **evecs(*)** (*complex real*) -- (input/output)

      * **primme** (*ptr*) -- (input) parameters structure.

      * **ierr** (*integer*) -- (output) error indicator; see *Error
        Codes*.


dprimme_f77
===========

dprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a real symmetric standard eigenproblem using double
   precision.

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

      * **ierr** (*integer*) -- (output) error indicator; see *Error
        Codes*.


zprimme_f77
===========

zprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblem. The arguments have the
   same meaning as in function "dprimme_f77()".

   Parameters:
      * **evals(*)** (*double precision*) -- (output)

      * **resNorms(*)** (*double precision*) -- (output)

      * **evecs(*)** (*complex double precision*) -- (input/output)

      * **primme** (*ptr*) -- (input) parameters structure.

      * **ierr** (*integer*) -- (output) error indicator; see *Error
        Codes*.


primme_set_member_f77
=====================

primme_set_member_f77(primme, label, value)

   Set a value in some field of the parameter structure.

   Parameters:
      * **primme** (*ptr*) -- (input) parameters structure.

      * **label** (*integer*) --

        field where to set value. One of:

           "PRIMME_n"
           "PRIMME_matrixMatvec"
           "PRIMME_applyPreconditioner"
           "PRIMME_numProcs"
           "PRIMME_procID"
           "PRIMME_commInfo"
           "PRIMME_nLocal"
           "PRIMME_globalSumReal"
           "PRIMME_numEvals"
           "PRIMME_target"
           "PRIMME_numTargetShifts"
           "PRIMME_targetShifts"
           "PRIMME_locking"
           "PRIMME_initSize"
           "PRIMME_numOrthoConst"
           "PRIMME_maxBasisSize"
           "PRIMME_minRestartSize"
           "PRIMME_maxBlockSize"
           "PRIMME_maxMatvecs"
           "PRIMME_maxOuterIterations"
           "PRIMME_intWorkSize"
           "PRIMME_realWorkSize"
           "PRIMME_iseed"
           "PRIMME_intWork"
           "PRIMME_realWork"
           "PRIMME_aNorm"
           "PRIMME_eps"
           "PRIMME_printLevel"
           "PRIMME_outputFile"
           "PRIMME_matrix"
           "PRIMME_preconditioner"
           "PRIMME_restartingParams_scheme".
           "PRIMME_restartingParams_maxPrevRetain"
           "PRIMME_correctionParams_precondition"
           "PRIMME_correctionParams_robustShifts"
           "PRIMME_correctionParams_maxInnerIterations"
           "PRIMME_correctionParams_projectors_LeftQ"
           "PRIMME_correctionParams_projectors_LeftX"
           "PRIMME_correctionParams_projectors_RightQ"
           "PRIMME_correctionParams_projectors_RightX"
           "PRIMME_correctionParams_projectors_SkewQ"
           "PRIMME_correctionParams_projectors_SkewX"
           "PRIMME_correctionParams_convTest"
           "PRIMME_correctionParams_relTolBase"
           "PRIMME_stats_numOuterIterations"
           "PRIMME_stats_numRestarts"
           "PRIMME_stats_numMatvecs"
           "PRIMME_stats_numPreconds"
           "PRIMME_stats_elapsedTime"
           "PRIMME_dynamicMethodSwitch"
           "PRIMME_massMatrixMatvec"

      * **value** --

        (input) value to set.

        If the type of the option is integer ("int", "PRIMME_INT",
        "size_t"), the type of "value" should be as long as
        "PRIMME_INT", which is "integer*8" by default.

   Note: **Don't use** this function inside PRIMME's callback
     functions, e.g., "matrixMatvec" or "applyPreconditioner", or in
     functions called by these functions.


primmetop_get_member_f77
========================

primmetop_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   Parameters:
      * **primme** (*ptr*) -- (input) parameters structure.

      * **label** (*integer*) -- (input) field where to get value.
        One of the detailed in function "primmetop_set_member_f77()".

      * **value** --

        (output) value of the field.

        If the type of the option is integer ("int", "PRIMME_INT",
        "size_t"), the type of "value" should be as long as
        "PRIMME_INT", which is "integer*8" by default.

   Note: **Don't use** this function inside PRIMME's callback
     functions, e.g., "matrixMatvec" or "applyPreconditioner", or in
     functions called by these functions. In those cases use
     "primme_get_member_f77()".

   Note: When "label" is one of "PRIMME_matrixMatvec",
     "PRIMME_applyPreconditioner", "PRIMME_commInfo",
     "PRIMME_intWork", "PRIMME_realWork", "PRIMME_matrix" and
     "PRIMME_preconditioner", the returned "value" is a C pointer
     ("void*"). Use Fortran pointer or other extensions to deal with
     it. For instance:

        use iso_c_binding
        MPI_Comm comm

        comm = MPI_COMM_WORLD
        call primme_set_member_f77(primme, PRIMME_commInfo, comm)
        ...
        subroutine par_GlobalSumDouble(x,y,k,primme)
        use iso_c_binding
        implicit none
        ...
        MPI_Comm, pointer :: comm
        type(c_ptr) :: pcomm

        call primme_get_member_f77(primme, PRIMME_commInfo, pcomm)
        call c_f_pointer(pcomm, comm)
        call MPI_Allreduce(x,y,k,MPI_DOUBLE,MPI_SUM,comm,ierr)

     Most users would not need to retrieve these pointers in their
     programs.


primmetop_get_prec_shift_f77
============================

primmetop_get_prec_shift_f77(primme, index, value)

   Get the value in some position of the array
   "ShiftsForPreconditioner".

   Parameters:
      * **primme** (*ptr*) -- (input) parameters structure.

      * **index** (*integer*) -- (input) position of the array; the
        first position is 1.

      * **value** -- (output) value of the array at that position.


primme_get_member_f77
=====================

primme_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   Parameters:
      * **primme** (*ptr*) -- (input) parameters structure.

      * **label** (*integer*) -- (input) field where to get value.
        One of the detailed in function "primmetop_set_member_f77()".

      * **value** --

        (output) value of the field.

        If the type of the option is integer ("int", "PRIMME_INT",
        "size_t"), the type of "value" should be as long as
        "PRIMME_INT", which is "integer*8" by default.

   Note: Use this function exclusively inside PRIMME's callback
     functions, e.g., "matrixMatvec" or "applyPreconditioner", or in
     functions called by these functions. Otherwise, e.g., from the
     main program, use the function "primmetop_get_member_f77()".

   Note: When "label" is one of "PRIMME_matrixMatvec",
     "PRIMME_applyPreconditioner", "PRIMME_commInfo",
     "PRIMME_intWork", "PRIMME_realWork", "PRIMME_matrix" and
     "PRIMME_preconditioner", the returned "value" is a C pointer
     ("void*"). Use Fortran pointer or other extensions to deal with
     it. For instance:

        use iso_c_binding
        MPI_Comm comm

        comm = MPI_COMM_WORLD
        call primme_set_member_f77(primme, PRIMME_commInfo, comm)
        ...
        subroutine par_GlobalSumDouble(x,y,k,primme)
        use iso_c_binding
        implicit none
        ...
        MPI_Comm, pointer :: comm
        type(c_ptr) :: pcomm

        call primme_get_member_f77(primme, PRIMME_commInfo, pcomm)
        call c_f_pointer(pcomm, comm)
        call MPI_Allreduce(x,y,k,MPI_DOUBLE,MPI_SUM,comm,ierr)

     Most users would not need to retrieve these pointers in their
     programs.


primme_get_prec_shift_f77
=========================

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

Appendix
********


Types
=====

The following data types are macros used in PRIMME as followed.

PRIMME_INT

   Integer type used in matrix dimensions (such as "n" and "nLocal")
   and counters (such as "numMatvecs").

   The integer size is controlled by the compilation flag
   "PRIMME_INT_SIZE", see *Making and Linking*.

PRIMME_COMPLEX_FLOAT

   Macro that is "complex float" in C and "std::complex<float>" in
   C++.

PRIMME_COMPLEX_DOUBLE

   Macro that is "complex double" in C and "std::complex<double>" in
   C++.


primme_params
=============

primme_params

   Structure to set the problem matrices and eigensolver options.

   PRIMME_INT n

      Dimension of the matrix.

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read by "dprimme()".

   void (*matrixMatvec)(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr)

      Block matrix-multivector multiplication, y = A x in solving A x
      = \lambda x or A x = \lambda B x.

      Parameters:
         * **x** -- matrix of size "nLocal" x "blockSize" in column-
           major order with leading dimension "ldx".

         * **ldx** -- the leading dimension of the array "x".

         * **y** -- matrix of size "nLocal" x "blockSize" in column-
           major order with leading dimension "ldy".

         * **ldy** -- the leading dimension of the array "y".

         * **blockSize** -- number of columns in "x" and "y".

         * **primme** -- parameters structure.

         * **ierr** -- output error code; if it is set to non-zero,
           the current call to PRIMME will stop.

      The actual type of "x" and "y" depends on which function is
      being calling. For "dprimme()", it is "double", for "zprimme()"
      it is "PRIMME_COMPLEX_DOUBLE", for "sprimme()" it is "float" and
      for "cprimme()" it is "PRIMME_COMPLEX_FLOAT".

      Input/output:

            "primme_initialize()" sets this field to NULL;
            this field is read by "dprimme()".

   Note: If you have performance issues with leading dimension
     different from "nLocal", set "ldOPs" to "nLocal".

   void (*applyPreconditioner)(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr)

      Block preconditioner-multivector application, y = M^{-1}x where
      M is usually an approximation of A - \sigma I or A - \sigma B
      for finding eigenvalues close to \sigma. The function follows
      the convention of "matrixMatvec".

      Input/output:

            "primme_initialize()" sets this field to NULL;
            this field is read by "dprimme()".

   void (*massMatrixMatvec)(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr)

      Block matrix-multivector multiplication, y = B x in solving A x
      = \lambda B x. The function follows the convention of
      "matrixMatvec".

      Input/output:

            "primme_initialize()" sets this field to NULL;
            this field is read by "dprimme()".

      Warning: Generalized eigenproblems not implemented in current
        version. This member is included for future compatibility.

   int numProcs

      Number of processes calling "dprimme()" or "zprimme()" in
      parallel.

      Input/output:

            "primme_initialize()" sets this field to 1;
            this field is read by "dprimme()".

   int procID

      The identity of the local process within a parallel execution
      calling "dprimme()" or "zprimme()". Only the process with id 0
      prints information.

      Input/output:

            "primme_initialize()" sets this field to 0;
            "dprimme()" sets this field to 0 if "numProcs" is 1;
            this field is read by "dprimme()".

   int nLocal

      Number of local rows on this process.

      Input/output:

            "primme_initialize()" sets this field to 0;
            "dprimme()" sets this field to "n" if "numProcs" is 1;
            this field is read by "dprimme()".

   void *commInfo

      A pointer to whatever parallel environment structures needed.
      For example, with MPI, it could be a pointer to the MPI
      communicator. PRIMME does not use this. It is available for
      possible use in user functions defined in "matrixMatvec",
      "applyPreconditioner", "massMatrixMatvec" and "globalSumReal".

      Input/output:

            "primme_initialize()" sets this field to NULL;

   void (*globalSumReal)(void *sendBuf, void *recvBuf, int *count, primme_params *primme, int *ierr)

      Global sum reduction function. No need to set for sequential
      programs.

      Parameters:
         * **sendBuf** -- array of size "count" with the local input
           values.

         * **recvBuf** -- array of size "count" with the global
           output values so that the i-th element of recvBuf is the
           sum over all processes of the i-th element of "sendBuf".

         * **count** -- array size of "sendBuf" and "recvBuf".

         * **primme** -- parameters structure.

         * **ierr** -- output error code; if it is set to non-zero,
           the current call to PRIMME will stop.

      The actual type of "sendBuf" and "recvBuf" depends on which
      function is being calling. For "dprimme()" and "zprimme()" it is
      "double", and for "sprimme()" and  "cprimme()" it is "float".
      Note that "count" is the number of values of the actual type.

      Input/output:

            "primme_initialize()" sets this field to an internal function;
            "dprimme()" sets this field to an internal function if "numProcs" is 1 and "globalSumReal" is NULL;
            this field is read by "dprimme()".

      When MPI is used, this can be a simply wrapper to
      MPI_Allreduce() as shown below:

         void par_GlobalSumForDouble(void *sendBuf, void *recvBuf, int *count,
                                  primme_params *primme, int *ierr) {
            MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;
            if(MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM,
                          communicator) == MPI_SUCCESS) {
               *ierr = 0;
            } else {
               *ierr = 1;
            }
         }

      }

      When calling "sprimme()" and "cprimme()" replace "MPI_DOUBLE" by
      "`MPI_FLOAT".

   int numEvals

      Number of eigenvalues wanted.

      Input/output:

            "primme_initialize()" sets this field to 1;
            this field is read by "primme_set_method()" (see *Preset Methods*) and "dprimme()".

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
         Closest in absolute value to the shifts in "targetShifts".

      "primme_largest_abs"
         Furthest in absolute value to the shifts in "targetShifts".

      Input/output:

            "primme_initialize()" sets this field to "primme_smallest";
            this field is read by "dprimme()".

   int numTargetShifts

      Size of the array "targetShifts". Used only when "target" is
      "primme_closest_geq", "primme_closest_leq", "primme_closest_abs"
      or "primme_largest_abs". The default values is 0.

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read by "dprimme()".

   double *targetShifts

      Array of shifts, at least of size "numTargetShifts". Used only
      when "target" is "primme_closest_geq", "primme_closest_leq",
      "primme_closest_abs" or "primme_largest_abs".

      Eigenvalues are computed in order so that the i-th eigenvalue is
      the closest (or closest but left or closest but right, see
      "target") to the i-th shift. If "numTargetShifts" < "numEvals",
      the last shift given is used for all the remaining i's.

      Input/output:

            "primme_initialize()" sets this field to NULL;
            this field is read by "dprimme()".

      Note: Considerations for interior problems:

        * PRIMME will try to compute the eigenvalues in the order
          given in the "targetShifts". However, for code efficiency
          and robustness, the shifts should be ordered. Order them in
          ascending (descending) order for shifts closer to the lower
          (higher) end of the spectrum.

        * If some shift is close to the lower (higher) end of the
          spectrum, use either "primme_closest_geq"
          ("primme_closest_leq") or "primme_closest_abs".

        * "primme_closest_leq" and "primme_closest_geq" are more
          efficient than "primme_closest_abs".

        * For interior eigenvalues larger "maxBasisSize" is usually
          more robust.

        * To find the largest magnitude eigenvalues set "target" to
          "primme_largest_abs", "numTargetShifts" to 1 and
          "targetShifts" to an array with a zero value.

   int printLevel

      The level of message reporting from the code. All output is
      writen in "outputFile".

      One of:

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

         $1: Number of converged pairs up to now.
         $2: The index of the pair currently converged.
         $3: The eigenvalue.
         $4: Its residual norm.
         $5: The current number of matrix-vector products.
         $6: The current number of outer iterations.
         $7: The current elapsed time.
         $8: Index within the block of the targeted pair .
         $9: QMR norm of the linear system residual.

      In parallel programs, output is produced in call with "procID" 0
      when "printLevel" is from 0 to 4. If "printLevel" is 5 output
      can be produced in any of the parallel calls.

      Input/output:

            "primme_initialize()" sets this field to 1;
            this field is read by "dprimme()".

   Note: Convergence history for plotting may be produced simply by:

        grep OUT outpufile | awk '{print $8" "$14}' > out
        grep INN outpufile | awk '{print $3" "$11}' > inn

     Then in Matlab:

        plot(out(:,1),out(:,2),'bo');hold; plot(inn(:,1),inn(:,2),'r');

     Or in gnuplot:

        plot 'out' w lp, 'inn' w lp

   double aNorm

      An estimate of the norm of A, which is used in the default
      convergence criterion (see "eps").

      If "aNorm" is less than or equal to 0, the code uses the largest
      absolute Ritz value seen. On return, "aNorm" is then replaced
      with that value.

      Input/output:

            "primme_initialize()" sets this field to 0.0;
            this field is read and written by "dprimme()".

   double eps

      If "convTestFun" is NULL, an eigenpairs is marked as converged
      when the 2-norm of the residual vector is less than "eps" *
      "aNorm". The residual vector is A x - \lambda x or A x - \lambda
      B x.

      The default value is machine precision times 10^4.

      Input/output:

            "primme_initialize()" sets this field to 0.0;
            this field is read and written by "dprimme()".

   FILE *outputFile

      Opened file to write down the output.

      Input/output:

            "primme_initialize()" sets this field to the standard output;
            this field is read by "dprimme()" and "primme_display_params()".

   int dynamicMethodSwitch

      If this value is 1, it alternates dynamically between
      "PRIMME_DEFAULT_MIN_TIME" and "PRIMME_DEFAULT_MIN_MATVECS",
      trying to identify the fastest method.

      On exit, it holds a recommended method for future runs on this
      problem:

            -1: use "PRIMME_DEFAULT_MIN_MATVECS" next time.
            -2: use "PRIMME_DEFAULT_MIN_TIME" next time.
            -3: close call, use "PRIMME_DYNAMIC" next time again.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

      Note: Even for expert users we do not recommend setting
        "dynamicMethodSwitch" directly, but through
        "primme_set_method()".

      Note: The code obtains timings by the "gettimeofday" Unix
        utility. If a cheaper, more accurate timer is available,
        modify the "PRIMMESRC/COMMONSRC/wtime.c"

   int locking

      If set to 1, hard locking will be used (locking converged
      eigenvectors out of the search basis). If set to 0, the code
      will try to use soft locking (à la ARPACK), when large enough
      "minRestartSize" is available.

      Input/output:

            "primme_initialize()" sets this field to -1;
            written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

   int initSize

      On input, the number of initial vector guesses provided in
      "evecs" argument in "dprimme()" or "zprimme()".

      On output, "initSize" holds the number of converged eigenpairs.
      Without "locking" all "numEvals" approximations are in "evecs"
      but only the "initSize" ones are converged.

      During execution, it holds the current number of converged
      eigenpairs. In addition, if locking is used, these are
      accessible in "evals" and "evecs".

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read and written by "dprimme()".

   PRIMME_INT ldevecs

      The leading dimension of "evecs". The default is "nLocal".

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read by "dprimme()".

   int numOrthoConst

      Number of vectors to be used as external orthogonalization
      constraints. These vectors are provided in the first
      "numOrthoConst" positions of the "evecs" argument in "dprimme()"
      or "zprimme()" and must be orthonormal.

      PRIMME finds new eigenvectors orthogonal to these constraints
      (equivalent to solving the problem with (I-YY^*)A(I-YY^*) and
      (I-YY^*)B(I-YY^*) matrices where Y are the given constraint
      vectors). This is a handy feature if some eigenvectors are
      already known, or for finding more eigenvalues after a call to
      "dprimme()" or "zprimme()", possibly with different parameters
      (see an example in "TEST/ex_zseq.c").

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read by "dprimme()".

   int maxBasisSize

      The maximum basis size allowed in the main iteration. This has
      memory implications.

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read and written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

   int minRestartSize

      Maximum Ritz vectors kept after restarting the basis.

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read and written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

   int maxBlockSize

      The maximum block size the code will try to use.

      The user should set this based on the architecture specifics of
      the target computer, as well as any a priori knowledge of
      multiplicities. The code does *not* require that "maxBlockSize"
      > 1 to find multiple eigenvalues. For some methods, keeping to 1
      yields the best overall performance.

      Input/output:

            "primme_initialize()" sets this field to 1;
            this field is read and written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

      Note: Inner iterations of QMR are not performed in a block
        fashion. Every correction equation from a block is solved
        independently.

   PRIMME_INT maxMatvecs

      Maximum number of matrix vector multiplications (approximately
      equal to the number of preconditioning operations) that the code
      is allowed to perform before it exits.

      Input/output:

            "primme_initialize()" sets this field to "INT_MAX";
            this field is read by "dprimme()".

   PRIMME_INT maxOuterIterations

      Maximum number of outer iterations that the code is allowed to
      perform before it exits.

      Input/output:

            "primme_initialize()" sets this field to "INT_MAX";
            this field is read by "dprimme()".

   int intWorkSize

      If "dprimme()" or "zprimme()" is called with all arguments as
      NULL except for "primme_params" then PRIMME returns immediately
      with "intWorkSize" containing the size *in bytes* of the integer
      workspace that will be required by the parameters set in PRIMME.

      Otherwise if "intWorkSize" is not 0, it should be the size of
      the integer work array *in bytes* that the user provides in
      "intWork". If "intWorkSize" is 0, the code will allocate the
      required space, which can be freed later by calling
      "primme_free()".

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read and written by "dprimme()".

   size_t realWorkSize

      If "dprimme()" or "zprimme()" is called with all arguments as
      NULL except for "primme_params" then PRIMME returns immediately
      with "realWorkSize" containing the size *in bytes* of the real
      workspace that will be required by the parameters set in PRIMME.

      Otherwise if "realWorkSize" is not 0, it should be the size of
      the real work array *in bytes* that the user provides in
      "realWork". If "realWorkSize" is 0, the code will allocate the
      required space, which can be freed later by calling
      "primme_free()".

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read and written by "dprimme()".

   int *intWork

      Integer work array.

      If NULL, the code will allocate its own workspace. If the
      provided space is not enough, the code will return the error
      code "-37".

      On exit, the first element shows if a locking problem has
      occurred. Using locking for large "numEvals" may, in some rare
      cases, cause some pairs to be practically converged, in the
      sense that their components are in the basis of "evecs". If this
      is the case, a Rayleigh Ritz on returned "evecs" would provide
      the accurate eigenvectors (see [r4]).

      Input/output:

            "primme_initialize()" sets this field to NULL;
            this field is read and written by "dprimme()".

   void *realWork

      Real work array.

      If NULL, the code will allocate its own workspace. If the
      provided space is not enough, the code will return the error
      code "-36".

      Input/output:

            "primme_initialize()" sets this field to NULL;
            this field is read and written by "dprimme()".

   PRIMME_INT iseed

      The "PRIMME_INT iseed[4]" is an array with the seeds needed by
      the LAPACK dlarnv and zlarnv.

      The default value is an array with values -1, -1, -1 and -1. In
      that case, "iseed" is set based on the value of "procID" to
      avoid every parallel process generating the same sequence of
      pseudorandom numbers.

      Input/output:

            "primme_initialize()" sets this field to "[-1, -1, -1, -1]";
            this field is read and written by "dprimme()".

   void *matrix

      This field may be used to pass any required information in the
      matrix-vector product "matrixMatvec".

      Input/output:

            "primme_initialize()" sets this field to NULL;

   void *preconditioner

      This field may be used to pass any required information in the
      preconditioner function "applyPreconditioner".

      Input/output:

            "primme_initialize()" sets this field to NULL;

   double *ShiftsForPreconditioner

      Array of size "blockSize" provided during execution of
      "dprimme()" and "zprimme()" holding the shifts to be used (if
      needed) in the preconditioning operation.

      For example if the block size is 3, there will be an array of
      three shifts in "ShiftsForPreconditioner". Then the user can
      invert a shifted preconditioner for each of the block vectors
      (M-ShiftsForPreconditioner_i)^{-1} x_i. Classical Davidson
      (diagonal) preconditioning is an example of this.

         this field is read and written by "dprimme()".

   primme_init initBasisMode

      Select how the search subspace basis is initialized up to
      "minRestartSize" vectors if not enough initial vectors are
      provided (see "initSize"):

      * "primme_init_krylov", with a block Krylov subspace generated
        by the matrix problem and the last initial vectors if given or
        a random vector otherwise; the size of the block is
        "maxBlockSize".

      * "primme_init_random", with random vectors.

      * "primme_init_user", the initial basis will have only initial
        vectors if given, or a single random vector.

      Input/output:

            "primme_initialize()" sets this field to "primme_init_krylov";
            this field is read by "dprimme()".

   primme_projection projectionParams.projection

      Select the extraction technique, i.e., how the approximate
      eigenvectors x_i and eigenvalues \lambda_i are computed from the
      search subspace \mathcal V:

      * "primme_proj_RR", Rayleigh-Ritz, Ax_i - Bx_i\lambda_i \perp
        \mathcal V.

      * "primme_proj_harmonic", Harmonic Rayleigh-Ritz, Ax_i -
        Bx_i\lambda_i \perp (A-\tau B)\mathcal V, where \tau is the
        current target shift (see "targetShifts").

      * "primme_proj_refined", refined extraction, compute ||x_i||=1
        so that minimizes ||(A-\tau B)x_i||; the eigenvalues are
        computed as the Rayleigh quotients,
        \lambda_i=\frac{x_i^*Ax_i}{x_i^*Bx_i}.

      Input/output:

            "primme_initialize()" sets this field to "primme_proj_default";
            "primme_set_method()" and "dprimme()" sets it to "primme_proj_RR" if it is "primme_proj_default".

   primme_restartscheme restartingParams.scheme

      Select a restarting strategy:

      * "primme_thick", Thick restarting. This is the most efficient
        and robust in the general case.

      * "primme_dtr", Dynamic thick restarting. Helpful without
        preconditioning but it is expensive to implement.

      Input/output:

            "primme_initialize()" sets this field to "primme_thick";
            written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

   int restartingParams.maxPrevRetain

      Number of approximations from previous iteration to be retained
      after restart (this is the locally optimal restarting, see
      [r2]). The restart size is "minRestartSize" plus
      "maxPrevRetain".

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read and written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

   int correctionParams.precondition

      Set to 1 to use preconditioning. Make sure "applyPreconditioner"
      is not NULL then!

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read and written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

   int correctionParams.robustShifts

      Set to 1 to use robust shifting. It tries to avoid stagnation
      and misconvergence by providing as shifts in
      "ShiftsForPreconditioner" the Ritz values displaced by an
      approximation of the eigenvalue error.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

   int correctionParams.maxInnerIterations

      Control the maximum number of inner QMR iterations:

      * 0:  no inner iterations;

      * >0: perform at most that number of inner iterations per
        outer step;

      * <0: perform at most the rest of the remaining matrix-vector
        products up to reach "maxMatvecs".

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read and written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

      See also "convTest".

   double correctionParams.relTolBase

      Parameter used when "convTest" is
      "primme_decreasing_LTolerance".

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

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

      Input/output:

            "primme_initialize()" sets this field to "primme_adaptive_ETolerance";
            written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

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

      Consider the current selected Ritz value \Lambda and vectors X,
      the residual associated vectors R=AX-X\Lambda, the previous
      locked vectors Q, and the preconditioner M^{-1}.

      When "maxInnerIterations" is 0, the correction D appended to the
      basis in GD is:

      +----------+---------+------------------------------------------------------------+
      | RightX   | SkewX   | D                                                          |
      +==========+=========+============================================================+
      | 0        | 0       | M^{-1}R (Classic GD)                                       |
      +----------+---------+------------------------------------------------------------+
      | 1        | 0       | M^{-1}(R-Delta X) (cheap Olsen's Method)                   |
      +----------+---------+------------------------------------------------------------+
      | 1        | 1       | (I- M^{-1}X(X^*M^{-1}X)^{-1}X^*)M^{-1}R (Olsen's Method)   |
      +----------+---------+------------------------------------------------------------+
      | 0        | 1       | error                                                      |
      +----------+---------+------------------------------------------------------------+

      Where \Delta is a diagonal matrix that \Delta_{i,i} holds an
      estimation of the error of the approximate eigenvalue
      \Lambda_{i,i}.

      The values of "RightQ", "SkewQ", "LeftX" and "LeftQ" are
      ignored.

      When "maxInnerIterations" is not 0, the correction D in Jacobi-
      Davidson results from solving:

         P_Q^l P_X^l (A-\sigma I) P_X^r P_Q^r M^{-1} D' = -R, \ \ \  D
         = P_X^r P_Q^l M^{-1}D'.

      For "LeftQ":

            0: P_Q^l = I;
            1: P_Q^l = I - QQ^*.

      For "LeftX":

            0: P_X^l = I;
            1: P_X^l = I - XX^*.

      For "RightQ" and "SkewQ":

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

      For "RightX" and "SkewX":

      +----------+---------+---------------------------------+
      | RightX   | SkewX   | P_X^r                           |
      +==========+=========+=================================+
      | 0        | 0       | I                               |
      +----------+---------+---------------------------------+
      | 1        | 0       | I - XX^*                        |
      +----------+---------+---------------------------------+
      | 1        | 1       | I - KX(X^*KX)^{-1}X^*           |
      +----------+---------+---------------------------------+
      | 0        | 1       | error                           |
      +----------+---------+---------------------------------+

      Input/output:

            "primme_initialize()" sets all of them to 0;
            this field is written by "primme_set_method()" (see *Preset Methods*);
            this field is read by "dprimme()".

      See [r3] for a study about different projector configurations in
      JD.

   PRIMME_INT ldOPs

      Recommended leading dimension to be used in "matrixMatvec",
      "applyPreconditioner" and "massMatrixMatvec". The default value
      is zero, which means no user recommendation. In that case,
      PRIMME computes ldOPs internally to get better memory
      performance.

      Input/output:

            "primme_initialize()" sets this field to 0;
            this field is read by "dprimme()".

   PRIMME_INT stats.numOuterIterations

      Hold the number of outer iterations. The value is available
      during execution and at the end.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   PRIMME_INT stats.numRestarts

      Hold the number of restarts during execution and at the end.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   PRIMME_INT stats.numMatvecs

      Hold how many vectors the operator in "matrixMatvec" has been
      applied on. The value is available during execution and at the
      end.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   PRIMME_INT stats.numPreconds

      Hold how many vectors the operator in "applyPreconditioner" has
      been applied on. The value is available during execution and at
      the end.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   PRIMME_INT stats.numGlobalSum

      Hold how many times "globalSumReal" has been called. The value
      is available during execution and at the end.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   double stats.volumeGlobalSum

      Hold how many "REAL" have been reduced by "globalSumReal". The
      value is available during execution and at the end.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   double stats.elapsedTime

      Hold the wall clock time spent by the call to "dprimme()" or
      "zprimme()". The value is available at the end of the execution.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   double stats.timeMatvec

      Hold the wall clock time spent by "matrixMatvec". The value is
      available at the end of the execution.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   double stats.timePrecond

      Hold the wall clock time spent by "applyPreconditioner". The
      value is available at the end of the execution.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   double stats.timeOrtho

      Hold the wall clock time spent by orthogonalization. The value
      is available at the end of the execution.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   double stats.timeGlobalSum

      Hold the wall clock time spent by "globalSumReal". The value is
      available at the end of the execution.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   double stats.estimateMinEVal

      Hold the estimation of the smallest eigenvalue for the current
      eigenproblem. The value is available during execution and at the
      end.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   double stats.estimateMaxEVal

      Hold the estimation of the largest eigenvalue for the current
      eigenproblem. The value is available during execution and at the
      end.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   double stats.estimateLargestSVal

      Hold the estimation of the largest singular value (i.e., the
      absolute value of the eigenvalue with largest absolute value)
      for the current eigenproblem. The value is available during
      execution and at the end.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   double stats.maxConvTol

      Hold the maximum residual norm of the converged eigenvectors.
      The value is available during execution and at the end.

      Input/output:

            "primme_initialize()" sets this field to 0;
            written by "dprimme()".

   void (*convTestFun)(double *eval, void *evecs, double *resNorm, int *isconv, primme_params *primme, int *ierr)

      Function that evaluates if the approximate eigenpair has
      converged. If NULL, it is used the default convergence criteria
      (see "eps").

      Parameters:
         * **eval** -- the approximate value to evaluate.

         * **x** -- one dimensional array of size "nLocal"
           containing the approximate vector; it can be NULL. The
           actual type depends on which function is being calling. For
           "dprimme()", it is "double", for "zprimme()" it is
           "PRIMME_COMPLEX_DOUBLE", for "sprimme()" it is "float" and
           for for "cprimme()" it is "PRIMME_COMPLEX_FLOAT".

         * **resNorm** -- the norm of residual vector.

         * **isconv** -- (output) the function sets zero if the pair
           is not converged and non zero otherwise.

         * **primme** -- parameters structure.

         * **ierr** -- output error code; if it is set to non-zero,
           the current call to PRIMME will stop.

      Input/output:

            "primme_initialize()" sets this field to NULL;
            this field is read by "dprimme()".


Error Codes
===========

The functions "dprimme()" and "zprimme()" return one of the next
values:

* 0: success.

* 1: reported only amount of required memory.

* -1: failed in allocating int or real workspace.

* -2: malloc failed in allocating a permutation integer array.

* -3: main_iter() encountered problem; the calling stack of the
  functions where the error occurred was printed in "stderr".

* -4: if argument "primme" is NULL.

* -5: if "n" < 0 or "nLocal" < 0 or "nLocal" > "n".

* -6: if "numProcs" < 1.

* -7: if "matrixMatvec" is NULL.

* -8: if "applyPreconditioner" is NULL and "precondition" > 0.

* -10: if "numEvals" > "n".

* -11: if "numEvals" < 0.

* -12: if "eps" > 0 and "eps" < machine precision.

* -13: if "target" is not properly defined.

* -14: if "target" is one of "primme_closest_geq",
  "primme_closest_leq", "primme_closest_abs" or "primme_largest_abs"
  but "numTargetShifts" <= 0 (no shifts).

* -15: if "target" is one of "primme_closest_geq",
  "primme_closest_leq", "primme_closest_abs" or "primme_largest_abs"
  but "targetShifts" is NULL  (no shifts array).

* -16: if "numOrthoConst" < 0 or "numOrthoConst" > "n". (no free
  dimensions left).

* -17: if "maxBasisSize" < 2.

* -18: if "minRestartSize" < 0 or "minRestartSize" shouldn't be
  zero.

* -19: if "maxBlockSize" < 0 or "maxBlockSize" shouldn't be zero.

* -20: if "maxPrevRetain" < 0.

* -21: if "scheme" is not one of *primme_thick* or *primme_dtr*.

* -22: if "initSize" < 0.

* -23: if "locking" == 0 and "initSize" > "maxBasisSize".

* -24: if "locking" and "initSize" > "numEvals".

* -25: if "maxPrevRetain" + "minRestartSize" >= "maxBasisSize".

* -26: if "minRestartSize" >= "n".

* -27: if "printLevel" < 0 or "printLevel" > 5.

* -28: if "convTest" is not one of "primme_full_LTolerance",
  "primme_decreasing_LTolerance", "primme_adaptive_ETolerance" or
  "primme_adaptive".

* -29: if "convTest" == "primme_decreasing_LTolerance" and
  "relTolBase" <= 1.

* -30: if "evals" is NULL, but not "evecs" and "resNorms".

* -31: if "evecs" is NULL, but not "evals" and "resNorms".

* -32: if "resNorms" is NULL, but not "evecs" and "evals".

* -33: if "locking" == 0 and "minRestartSize" < "numEvals".

* -34: if "ldevecs" < "nLocal".

* -35: if "ldOPs" is not zero and less than "nLocal".

* -36: not enough memory for "realWork".

* -37: not enough memory for "intWork".

* -38: if "locking" == 0 and "target" is "primme_closest_leq" or
  "primme_closest_geq".


Preset Methods
==============

primme_preset_method

   PRIMME_DEFAULT_MIN_TIME

      Set as "PRIMME_JDQMR_ETol" when "target" is either
      "primme_smallest" or "primme_largest", and as "PRIMME_JDQMR"
      otherwise. This method is usually the fastest if the cost of the
      matrix vector product is inexpensive.

   PRIMME_DEFAULT_MIN_MATVECS

      Currently set as "PRIMME_GD_Olsen_plusK"; this method usually
      performs fewer matrix vector products than other methods, so
      it's a good choice when this operation is expensive.

   PRIMME_DYNAMIC

      Switches to the best method dynamically; currently, between
      methods "PRIMME_DEFAULT_MIN_TIME" and
      "PRIMME_DEFAULT_MIN_MATVECS".

      With "PRIMME_DYNAMIC" "primme_set_method()" sets
      "dynamicMethodSwitch" = 1 and makes the same changes as for
      method "PRIMME_DEFAULT_MIN_TIME".

   PRIMME_Arnoldi

      Arnoldi implemented à la Generalized Davidson.

      With "PRIMME_Arnoldi" "primme_set_method()" sets:

      * "locking" = 0;

      * "maxPrevRetain" = 0;

      * "precondition" = 0;

      * "maxInnerIterations" = 0.

   PRIMME_GD

      Generalized Davidson.

      With "PRIMME_GD" "primme_set_method()" sets:

      * "locking" = 0;

      * "maxPrevRetain" = 0;

      * "robustShifts" = 1;

      * "maxInnerIterations" = 0;

      * "RightX" = 0;

      * "SkewX" = 0.

   PRIMME_GD_plusK

      GD with locally optimal restarting.

      With "PRIMME_GD_plusK" "primme_set_method()" sets
      "maxPrevRetain" = 2 if "maxBlockSize" is 1 and "numEvals" > 1;
      otherwise it sets "maxPrevRetain" to "maxBlockSize". Also:

      * "locking" = 0;

      * "maxInnerIterations" = 0;

      * "RightX" = 0;

      * "SkewX" = 0.

   PRIMME_GD_Olsen_plusK

      GD+k and the cheap Olsen's Method.

      With "PRIMME_GD_Olsen_plusK" "primme_set_method()" makes the
      same changes as for method "PRIMME_GD_plusK" and sets "RightX" =
      1.

   PRIMME_JD_Olsen_plusK

      GD+k and Olsen's Method.

      With "PRIMME_JD_Olsen_plusK" "primme_set_method()" makes the
      same changes as for method "PRIMME_GD_plusK" and also sets
      "robustShifts" = 1, "RightX" to 1, and "SkewX" to 1.

   PRIMME_RQI

      (Accelerated) Rayleigh Quotient Iteration.

      With "PRIMME_RQI" "primme_set_method()" sets:

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

      Note: If "numTargetShifts" > 0 and "targetShifts" are
        provided, the interior problem solved uses these shifts in the
        correction equation. Therefore RQI becomes INVIT (inverse
        iteration) in that case.

   PRIMME_JDQR

      Jacobi-Davidson with fixed number of inner steps.

      With "PRIMME_JDQR" "primme_set_method()" sets:

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

   PRIMME_JDQMR

      Jacobi-Davidson with adaptive stopping criterion for inner Quasi
      Minimum Residual (QMR).

      With "PRIMME_JDQMR" "primme_set_method()" sets:

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

   PRIMME_JDQMR_ETol

      JDQMR but QMR stops after residual norm reduces by a 0.1 factor.

      With "PRIMME_JDQMR_ETol" "primme_set_method()" makes the same
      changes as for the method "PRIMME_JDQMR" and sets "convTest" =
      "primme_adaptive_ETolerance".

   PRIMME_SUBSPACE_ITERATION

      Subspace iteration.

      With "PRIMME_SUBSPACE_ITERATION" "primme_set_method()" sets:

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

   PRIMME_LOBPCG_OrthoBasis

      LOBPCG with orthogonal basis.

      With "PRIMME_LOBPCG_OrthoBasis" "primme_set_method()" sets:

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

   PRIMME_LOBPCG_OrthoBasis_Window

      LOBPCG with sliding window of "maxBlockSize" < 3 *** "numEvals".

      With "PRIMME_LOBPCG_OrthoBasis_Window" "primme_set_method()"
      sets:

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

Python Interface
****************

Primme.eigsh(A, k=6, M=None, sigma=None, which='LM', v0=None, ncv=None, maxiter=None, tol=0, return_eigenvectors=True, Minv=None, OPinv=None, mode='normal', lock=None, return_stats=False, maxBlockSize=0, minRestartSize=0, maxPrevRetain=0, method=None, **kargs)

   Find k eigenvalues and eigenvectors of the real symmetric square
   matrix or complex Hermitian matrix A.

   Solves "A * x[i] = w[i] * x[i]", the standard eigenvalue problem
   for w[i] eigenvalues with corresponding eigenvectors x[i].

   If M is specified, solves "A * x[i] = w[i] * M * x[i]", the
   generalized eigenvalue problem for w[i] eigenvalues with
   corresponding eigenvectors x[i]

   Parameters:
      * **A** (*An N x N matrix, array, sparse matrix, or
        LinearOperator*) -- the operation A * x, where A is a real
        symmetric matrix or complex Hermitian.

      * **k** (*int, optional*) -- The number of eigenvalues and
        eigenvectors desired.

      * **M** (*An N x N matrix, array, sparse matrix, or
        LinearOperator*) --

        (not supported yet) the operation M * x for the generalized
        eigenvalue problem

           A * x = w * M * x.

        M must represent a real, symmetric matrix if A is real, and
        must represent a complex, hermitian matrix if A is complex.
        For best results, the data type of M should be the same as
        that of A.

      * **sigma** (*real, optional*) -- Find eigenvalues near sigma.

      * **v0** (*N x i, ndarray, optional*) -- Starting vectors for
        iteration.

      * **ncv** (*int, optional*) -- The maximum size of the basis

      * **which** (*str ['LM' | 'SM' | 'LA' | 'SA' | 'BE']*) --

        If A is a complex hermitian matrix, 'BE' is invalid. Which *k*
        eigenvectors and eigenvalues to find:

           'LM' : Largest (in magnitude) eigenvalues

           'SM' : Smallest (in magnitude) eigenvalues

           'LA' : Largest (algebraic) eigenvalues

           'SA' : Smallest (algebraic) eigenvalues

           'BE' : Half (k/2) from each end of the spectrum (not
           supported)

        When sigma != None, 'which' refers to the shifted eigenvalues
        "w'[i]"

      * **maxiter** (*int, optional*) -- Maximum number of
        iterations.

      * **tol** (*float*) -- Accuracy for eigenvalues (stopping
        criterion). The default value is sqrt of machine precision.

      * **Minv** (*(not supported)*) --

      * **OPinv** (*N x N matrix, array, sparse matrix, or
        LinearOperator*) -- Preconditioner to accelerate the
        convergence. Usually it is an approximation of the inverse of
        (A - sigma*M).

      * **return_eigenvectors** (*bool*) -- Return eigenvectors
        (True) in addition to eigenvalues

      * **mode** (*string ['normal' | 'buckling' | 'cayley']*) --
        Only 'normal' mode is supported.

      * **lock** (*N x i, ndarray, optional*) -- Seek the
        eigenvectors orthogonal to these ones. The provided vectors
        *should* be orthonormal. Useful to not converge some already
        computed solutions.

      * **maxBlockSize** (*int, optional*) -- Maximum number of
        vectors added at every iteration.

      * **minRestartSize** (*int, optional*) -- Number of
        approximate eigenvectors kept from last iteration in restart.

      * **maxPrevRetain** (*int, optional*) -- Number of approximate
        eigenvectors kept from previous iteration in restart. Also
        referred as +k vectors in GD+k.

      * **method** (*int, optional*) --

        Preset method, one of:

        * DEFAULT_MIN_TIME : a variant of JDQMR,

        * DEFAULT_MIN_MATVECS : GD+k

        * DYNAMIC : choose dynamically between both previous
          methods.

        See a detailed description of the methods and other possible
        values in [2].

      * **report_stats** (*bool, optional*) -- If True, it is also
        returned extra information from PRIMME.

   Returns:
      * **w** (*array*) -- Array of k eigenvalues

      * **v** (*array*) -- An array representing the *k*
        eigenvectors. The column "v[:, i]" is the eigenvector
        corresponding to the eigenvalue "w[i]".

      * **stats** (*dict, optional (if return_stats)*) -- Extra
        information reported by PRIMME:

        * "numOuterIterations": number of outer iterations

        * "numRestarts": number of restarts

        * "numMatvecs": number of A*v

        * "numPreconds": number of OPinv*v

        * "elapsedTime": time that took

        * "estimateMinEVal": the leftmost Ritz value seen

        * "estimateMaxEVal": the rightmost Ritz value seen

        * "estimateLargestSVal": the largest singular value seen

        * "rnorms" : ||A*x[i] - x[i]*w[i]||

   Raises:
      "PrimmeError" -- When the requested convergence is not obtained.

      The PRIMME error code can be found as "err" attribute of the
      exception object.

   See also:

     "scipy.sparse.linalg.eigs()"
        eigenvalues and eigenvectors for a general (nonsymmetric)
        matrix A

     "Primme.svds()"
        singular value decomposition for a matrix A

   -[ Notes ]-

   This function is a wrapper to PRIMME functions to find the
   eigenvalues and eigenvectors [1].

   -[ References ]-

   [1] PRIMME Software, https://github.com/primme/primme

   [2] Preset Methods,
       http://www.cs.wm.edu/~andreas/software/doc/readme.html#preset-
       methods

   -[ Examples ]-

   >>> import Primme, scipy.sparse
   >>> A = scipy.sparse.spdiags(range(100), [0], 100, 100) # sparse diag. matrix
   >>> evals, evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA')
   >>> evals # the three largest eigenvalues of A
   array([ 99.,  98.,  97.])
   >>> evals, evecs = Primme.eigsh(A, 3, tol=1e-6, which='LA', lock=evecs)
   >>> evals # the next three largest eigenvalues
   array([ 96.,  95.,  94.])

MATLAB Interface
****************

function [varargout] = primme_eigs(varargin)

   "primme_eigs()" finds a few eigenvalues and eigenvectors of a real
   symmetric or Hermitian matrix, A, by calling the function
   "PRIMME_mex" (flag,dim,...). This in turn calls PRIMME. Full PRIMME
   functionality is supported.

   Input: [A, numEvals, target, opts, eigsMethod, P]

   Output: [evals, evecs, norms, primmeout]

   We provide different levels of function calls, similarly to MATLAB
   eigs():

      primme_eigs(A)
      primme_eigs(A, numEvals)
      primme_eigs(A, numEvals, target)
      primme_eigs(A, numEvals, target, opts)
      primme_eigs(A, numEvals, target, opts, eigsMethod)
      primme_eigs(A, numEvals, target, opts, eigsMethod, P)
      primme_eigs(A, numEvals, target, opts, eigsMethod, P1,P2)
      primme_eigs(A, numEvals, target, opts, eigsMethod, Pfun)
      primme_eigs(Afun, dim,...)

   "primme_eigs(A)" returns a vector of A's 6 largest algebraic
   eigenvalues. A must be real symmetric or complex Hermitian and
   should be large and sparse.

   "primme_eigs(Afun, dim)" accepts a function AFUN instead of a
   matrix. AFUN is a function handle and "y = Afun(x)" returns the
   matrix-vector product A*x. "primme_eigs(A,...)" could be replaced
   by primme_eigs(Afun, dim,...) in any of above levels of function
   calls. Examples are given in PRIMME_MEX_Readme.txt in the root
   directory of PRIMME_MEX folder.

   "[V, D] = primme_eigs(A)" returns a diagonal matrix D, of A's 6
   largest algebraic eigenvalues and a matrix V whose columns are the
   corresponding eigenvectors.

   "[V, D, norms, primmeout] = primme_eigs(A)" also returns an array
   of the residual norms of the computed eigenpairs, and a struct to
   report statistical information about "numOuterIterations",
   "numRestarts", "numMatvecs" and "numPreconds".

   "primme_eigs(A, numEvals)" finds the "numEvals" largest algebraic
   eigenvalues. numEvals must be less than the dimension of the matrix
   A.

   "primme_eigs(A, numEvals, target)" returns numEvals target
   eigenvalues. "target" could be a string like below:

   * 'LA' : "primme_largest" (default)

   * 'SA' : "primme_smallest"

   * 'CGT': "primme_closest_geq"

   * 'CLT': "primme_closest_leq"

   * 'CT' : "primme_closest_abs"

   "primme_eigs(A, numEvals, target, opts, eigsMethod)" specifies any
   of a set of possible options as explained below in the opts
   structure.

   "eigsMethod" is an integer specifying one of the preset methods in
   PRIMME:

   * 0:    "PRIMME_DYNAMIC", (default)        Switches dynamically
     to the best method

   * 1:    "PRIMME_DEFAULT_MIN_TIME",         Currently set at
     JDQMR_ETol

   * 2:    "PRIMME_DEFAULT_MIN_MATVECS",      Currently set at
     GD+block

   * 3:    "PRIMME_Arnoldi",                  obviously not an
     efficient choice

   * 4:    "PRIMME_GD",                       classical block
     Generalized Davidson

   * 5:    "PRIMME_GD_plusK",                 GD+k block GD with
     recurrence restarting

   * 6:    "PRIMME_GD_Olsen_plusK",           GD+k with approximate
     Olsen precond.

   * 7:    "PRIMME_JD_Olsen_plusK",           GD+k, exact Olsen (two
     precond per step)

   * 8:    "PRIMME_RQI",                      Rayleigh Quotient
     Iteration. Also INVIT, but for INVIT provide targetShifts

   * 9:    "PRIMME_JDQR",                     Original block, Jacobi
     Davidson

   * 10:   "PRIMME_JDQMR",                    Our block JDQMR method
     (similar to JDCG)

   * 11:   "PRIMME_JDQMR_ETol",               Slight, but efficient
     JDQMR modification

   * 12:   "PRIMME_SUBSPACE_ITERATION",       equiv. to
     GD(block,2*block)

   * 13:   "PRIMME_LOBPCG_OrthoBasis",        equiv. to
     GD(nev,3*nev)+nev

   * 14:   "PRIMME_LOBPCG_OrthoBasis_Window"  equiv. to
     GD(block,3*block)+block nev>block

   "primme_eigs(A, numEvals, target, opts, eigsMethod, P)"

   "primme_eigs(A, numEvals, target, opts, eigsMethod, P1, P2)" uses
   preconditioner P or P = P1*P2 to accelerate convergence of the
   methods. If P is [] then a preconditioner is not applied. P may be
   a function handle Pfun such that Pfun(x) returns Px.

   "opts" is an option structure which contain following parameters:

   * "aNorm": the estimate norm value of matrix A [{0.0}|scaler]

   * "eps": desired computing accuracy [{1e-12}|scaler]

   * "maxBlockSize": maximum block size the PRIMME uses [{1}|scaler]

   * "printLevel": different level reporting(0-5) [{1}|scaler]

   * "outputFile": output file name where user wants to save results

   * "precondition": set to 1 if use preconditioner [{0}|1]

   * isreal: the complexity of A represented by AFUN [{ture}|false]

   * "numTargetShifts": number of shifts for interior eigenvalues
     [{0}|scaler]

   * "targetShifts": shifts for interior eigenvalues [{}|vector]

   * "initSize": On INPUT, the number of initial guesses provided in
     evecs array. ON OUTPUT, the number of converged eigenpairs
     [{0}|scaler]

   * "numOrthoConst": Number of external orthogonalization
     constraints provided in the first numOrthoConst vectors of evecs
     [{0}|scaler]

   * locking: If set to 1, hard locking will be used, otherwise the
     code will try to use soft locking [{0}|1]

   * "maxBasisSize": maximum basis size allowed in the main
     iteration

   * "minRestartSize": minimum Ritz vectors to restart

   * "maxMatvecs": maximum number of matrix vector multiplications
     [{INT_MAX}|scaler]

   * "maxOuterIterations": maximum number of outer iterations
     [{INT_MAX}|scaler]

   * restartingParams. "scheme": the restart scheme [{primme_thick}|
     primme_dtr]

   * restartingParams. "maxPrevRetain": number of approximations
     from previous iteration to be retained after restart [{1}|scaler]

   * "robustShifts": set to 1 if use robustShifting to help avoid
     stagnation and misconverge [{0}|1]

   * "maxInnerIterations": number of inner QMR iterations
     [{0}|scaler]

   * "LeftQ": a projector with Q must be applied on the left [{0}|1]

   * "LeftX": a projector with X must be applied on the left [{0}|1]

   * "RightQ": a projector with Q must be applied on the right
     [{0}|1]

   * "RightX": a projector with X must be applied on the right
     [{0}|1]

   * "SkewQ": the Q right projector must be skew [{0}|1]

   * "SkewX": the X right projector must be skew [{0}|1]

   * "relTolBase": a legacy from calssical JDQR (recommend not use)

   * "convTest": how to stop the inner QMR Method

   * "iseed": set iseed value for initialization

   * "intWorkSize": memory size for integer workspace

   * "realWorkSize": memory size for real or complex workspace

   See also "Matlab/readme.txt".

Singular Value Problems
***********************

* C Library Interface

* FORTRAN Library Interface

* Python Interface

* Appendix

C Library Interface
*******************

The PRIMME SVDS interface is composed of the following functions. To
solve real and complex singular value problems call respectively:

   int sprimme_svds(float *svals, float *svecs, float *resNorms,
               primme_svds_params *primme_svds);

   int cprimme_svds(float *svals, PRIMME_COMPLEX_FLOAT *svecs, float *resNorms,
               primme_svds_params\*primme_svds);

   int dprimme_svds(double *svals, double *svecs, double *resNorms,
               primme_svds_params *primme);

   int zprimme_svds(double *svals, PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms,
               primme_svds_params *primme);

Other useful functions:

   void primme_svds_initialize(primme_svds_params *primme_svds);
   int primme_svds_set_method(primme_svds_preset_method method,
         primme_preset_method methodStage1, primme_preset_method methodStage2,
         primme_svds_params *primme_svds);
   void primme_svds_display_params(primme_svds_params primme_svds);
   void primme_svds_Free(primme_svds_params *primme_svds);

PRIMME SVDS stores its data on the structure "primme_svds_params". See
*Parameters Guide* for an introduction about its fields.


Running
=======

To use PRIMME SVDS, follow these basic steps.

1. Include:

      #include "primme.h"   /* header file is required to run primme */

2. Initialize a PRIMME SVDS parameters structure for default
   settings:

      primme_svds_params primme_svds;

      primme_svds_initialize(&primme_svds);

3. Set problem parameters (see also *Parameters Guide*), and,
   optionally, set one of the "preset methods":

      primme_svds.matrixMatvec = matrixMatvec; /* MV product */
      primme_svds.m = 1000;                    /* set problem dimension */
      primme_svds.n = 100;
      primme_svds.numSvals = 10;    /* Number of wanted singular values */
      primme_svds_set_method(primme_svds_hybrid, PRIMME_DEFAULT_METHOD,
                                PRIMME_DEFAULT_METHOD, &primme_svds);
      ...

4. Then to solve a real singular value problem call:

      ret = dprimme_svds(svals, svecs, resNorms, &primme_svds);

   The previous is the double precision call. There is available calls
   for complex double, single and complex single; check it out
   "zprimme_svds()", "sprimme_svds()" and "cprimme_svds()".

   To solve complex singular value problems call:

      ret = zprimme_svds(svals, svecs, resNorms, &primme_svds);

   The call arguments are:

   * *svals*, array to return the found singular values;

   * *svecs*, array to return the found left and right singular
     vectors;

   * *resNorms*, array to return the residual norms of the found
     triplets; and

   * *ret*, returned error code.

5. To free the work arrays in PRIMME SVDS:

      primme_svds_free(&primme_svds);


Parameters Guide
================

PRIMME SVDS stores the data on the structure "primme_svds_params",
which has the next fields:

   /* Basic */
   PRIMME_INT m;                    // number of rows of the matrix
   PRIMME_INT n;                 // number of columns of the matrix
   void (*matrixMatvec)(...);              // matrix-vector product
   int numSvals;              // how many singular triplets to find
   primme_svds_target target;      // which singular values to find
   double eps;               // tolerance of the converged triplets

   /* For parallel programs */
   int numProcs;          // number of processes
   int procID;            // rank of this process
   PRIMME_INT mLocal;     // number of rows stored in this process
   PRIMME_INT nLocal;     // number of columns stored in this process
   void (*globalSumReal)(...); // sum reduction among processes

   /* Accelerate the convergence */
   void (*applyPreconditioner)(...); // preconditioner-vector product
   int initSize;        // initial vectors as approximate solutions
   int maxBasisSize;
   int minRestartSize;
   int maxBlockSize;

   /* User data */
   void *commInfo;
   void *matrix;
   void *preconditioner;

   /* Advanced options */
   int numTargetShifts;        // for targeting interior values
   double *targetShifts;
   int numOrthoConst;   // orthogonal constrains to the vectors
   int locking;
   PRIMME_INT maxMatvecs;
   int intWorkSize;
   size_t realWorkSize;
   PRIMME_INT iseed[4];
   int *intWork;
   void *realWork;
   double aNorm;
   int printLevel;
   FILE * outputFile;
   primme_svds_operator method;
   primme_svds_operator methodStage2;
   primme_params primme;
   primme_params primmeStage2;

PRIMME SVDS requires the user to set at least the matrix dimensions
("m" x "n") and the matrix-vector product ("matrixMatvec"), as they
define the problem to be solved. For parallel programs, "mLocal",
"nLocal", "procID" and "globalSumReal" are also required.

In addition, most users would want to specify how many singular
triplets to find, and provide a preconditioner (if available).

It is useful to have set all these before calling
"primme_svds_set_method()". Also, if users have a preference on
"maxBasisSize", "maxBlockSize", etc, they should also provide them
into "primme_svds_params" prior to the "primme_svds_set_method()"
call. This helps "primme_svds_set_method()" make the right choice on
other parameters. It is sometimes useful to check the actual
parameters that PRIMME SVDS is going to use (before calling it) or
used (on return) by printing them with "primme_svds_display_params()".


Interface Description
=====================

The next enumerations and functions are declared in "primme.h".


sprimme_svds
------------

int sprimme_svds(float *svals, float *svecs, float *resNorms, primme_svds_params *primme_svds)

   Solve a real singular value problem.

   Parameters:
      * **svals** -- array at least of size "numSvals" to store the
        computed singular values; all processes in a parallel run
        return this local array with the same values.

      * **resNorms** -- array at least of size "numSvals" to store
        the residual norms of the computed triplets; all processes in
        parallel run return this local array with the same values.

      * **svecs** -- array at least of size ("mLocal" + "nLocal")
        times "numSvals" to store columnwise the (local part of the)
        computed left singular vectors and the right singular vectors.

      * **primme_svds** -- parameters structure.

   Returns:
      error indicator; see *Error Codes*.

   On input, "svecs" should start with the content of the
   "numOrthoConst" left vectors, followed by the "initSize" left
   vectors, followed by the "numOrthoConst" right vectors and followed
   by the "initSize" right vectors. The i-th left vector starts at
   svecs[i* "mLocal" ]. The i-th right vector starts at svecs[(
   "numOrthoConst" + "initSize" )* "mLocal" + i* "nLocal" ].

   On return, the i-th left singular vector starts at svecs[(
   "numOrthoConst" +i)* "mLocal" ]. The i-th right singular vector
   starts at svecs[( "numOrthoConst" + "initSize" )* "mLocal" + (
   "numOrthoConst" +i)* "nLocal" ]. The first vector has i=0.


dprimme_svds
------------

int dprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds)

   Solve a real singular value problem.

   Parameters:
      * **svals** -- array at least of size "numSvals" to store the
        computed singular values; all processes in a parallel run
        return this local array with the same values.

      * **resNorms** -- array at least of size "numSvals" to store
        the residual norms of the computed triplets; all processes in
        parallel run return this local array with the same values.

      * **svecs** -- array at least of size ("mLocal" + "nLocal")
        times "numSvals" to store columnwise the (local part of the)
        computed left singular vectors and the right singular vectors.

      * **primme_svds** -- parameters structure.

   Returns:
      error indicator; see *Error Codes*.

   On input, "svecs" should start with the content of the
   "numOrthoConst" left vectors, followed by the "initSize" left
   vectors, followed by the "numOrthoConst" right vectors and followed
   by the "initSize" right vectors. The i-th left vector starts at
   svecs[i* "mLocal" ]. The i-th right vector starts at svecs[(
   "numOrthoConst" + "initSize" )* "mLocal" + i* "nLocal" ].

   On return, the i-th left singular vector starts at svecs[(
   "numOrthoConst" +i)* "mLocal" ]. The i-th right singular vector
   starts at svecs[( "numOrthoConst" + "initSize" )* "mLocal" + (
   "numOrthoConst" +i)* "nLocal" ]. The first vector has i=0.


cprimme_svds
------------

int cprimme_svds(float *svals, PRIMME_COMPLEX_FLOAT *svecs, float *resNorms, primme_svds_params *primme_svds)

   Solve a complex singular value problem; see function
   "dprimme_svds()".


zprimme_svds
------------

int zprimme_svds(double *svals, PRIMME_COMPLEX_DOUBLE *svecs, double *resNorms, primme_svds_params *primme_svds)

   Solve a complex singular value problem; see function
   "dprimme_svds()".


primme_svds_initialize
----------------------

void primme_svds_initialize(primme_svds_params *primme_svds)

   Set PRIMME SVDS parameters structure to the default values.

   Parameters:
      * **primme_svds** -- parameters structure.


primme_svds_set_method
----------------------

int primme_svds_set_method(primme_svds_preset_method method, primme_preset_method methodStage1, primme_preset_method methodStage2, primme_svds_params *primme_svds)

   Set PRIMME SVDS parameters to one of the preset configurations.

   Parameters:
      * **method** --

        preset method to compute the singular triplets; one of

        * "primme_svds_default", currently set as
          "primme_svds_hybrid".

        * "primme_svds_normalequations", compute the eigenvectors of
          A^*A or A A^*.

        * "primme_svds_augmented", compute the eigenvectors of the
          augmented matrix, \left(\begin{array}{cc} 0 & A^* \\ A & 0
          \end{array}\right).

        * "primme_svds_hybrid", start with
          "primme_svds_normalequations"; use the resulting approximate
          singular vectors as initial vectors for
          "primme_svds_augmented" if the required accuracy was not
          achieved.

      * **methodStage1** -- preset method to compute the eigenpairs
        at the first stage; see available values at
        "primme_set_method()".

      * **methodStage2** -- preset method to compute the eigenpairs
        with the second stage of "primme_svds_hybrid"; see available
        values at "primme_set_method()".

      * **primme_svds** -- parameters structure.

   See also *Preset Methods*.


primme_svds_display_params
--------------------------

void primme_svds_display_params(primme_svds_params primme_svds)

   Display all printable settings of "primme_svds" into the file
   descriptor "outputFile".

   Parameters:
      * **primme_svds** -- parameters structure.


primme_svds_free
----------------

void primme_svds_free(primme_svds_params *primme_svds)

   Free memory allocated by PRIMME SVDS.

   Parameters:
      * **primme_svds** -- parameters structure.

FORTRAN Library Interface
*************************

The next enumerations and functions are declared in
"primme_svds_f77.h".


sprimme_svds_f77
================

sprimme_svds_f77(svals, svecs, resNorms, primme_svds)

   Solve a real singular value problem using single precision.

   Parameters:
      * **svals(*)** (*real*) -- (output) array at least of size
        "numSvals" to store the computed singular values; all
        processes in a parallel run return this local array with the
        same values.

      * **resNorms(*)** (*real*) -- array at least of size
        "numSvals" to store the residual norms of the computed
        triplets; all processes in parallel run return this local
        array with the same values.

      * **svecs(*)** (*real*) -- array at least of size ("mLocal" +
        "nLocal") times "numSvals" to store columnwise the (local part
        of the) computed left singular vectors and the right singular
        vectors.

      * **primme_svds** (*ptr*) -- parameters structure.

   Returns:
      error indicator; see *Error Codes*.


cprimme_svds_f77
================

cprimme_svds_f77(svals, svecs, resNorms, primme_svds)

   Solve a complex singular value problem using single precision.

   Parameters:
      * **svals(*)** (*real*) -- (output) array at least of size
        "numSvals" to store the computed singular values; all
        processes in a parallel run return this local array with the
        same values.

      * **resNorms(*)** (*real*) -- array at least of size
        "numSvals" to store the residual norms of the computed
        triplets; all processes in parallel run return this local
        array with the same values.

      * **svecs(*)** (*complex*) -- array at least of size ("mLocal"
        + "nLocal") times "numSvals" to store columnwise the (local
        part of the) computed left singular vectors and the right
        singular vectors.

      * **primme_svds** (*ptr*) -- parameters structure.

   Returns:
      error indicator; see *Error Codes*.


dprimme_svds_f77
================

dprimme_svds_f77(svals, svecs, resNorms, primme_svds)

   Solve a real singular value problem using double precision.

   Parameters:
      * **svals(*)** (*double precision*) -- (output) array at least
        of size "numSvals" to store the computed singular values; all
        processes in a parallel run return this local array with the
        same values.

      * **resNorms(*)** (*double precision*) -- array at least of
        size "numSvals" to store the residual norms of the computed
        triplets; all processes in parallel run return this local
        array with the same values.

      * **svecs(*)** (*double precision*) -- array at least of size
        ("mLocal" + "nLocal") times "numSvals" to store columnwise the
        (local part of the) computed left singular vectors and the
        right singular vectors.

      * **primme_svds** (*ptr*) -- parameters structure.

   Returns:
      error indicator; see *Error Codes*.


zprimme_svds_f77
================

zprimme_svds_f77(svals, svecs, resNorms, primme_svds)

   Solve a complex singular value problem using double precision.

   Parameters:
      * **svals(*)** (*double precision*) -- (output) array at least
        of size "numSvals" to store the computed singular values; all
        processes in a parallel run return this local array with the
        same values.

      * **resNorms(*)** (*double precision*) -- array at least of
        size "numSvals" to store the residual norms of the computed
        triplets; all processes in parallel run return this local
        array with the same values.

      * **svecs(*)** (*complex*16*) -- array at least of size
        ("mLocal" + "nLocal") times "numSvals" to store columnwise the
        (local part of the) computed left singular vectors and the
        right singular vectors.

      * **primme_svds** (*ptr*) -- parameters structure.

   Returns:
      error indicator; see *Error Codes*.


primme_svds_initialize_f77
==========================

primme_svds_initialize_f77(primme_svds)

   Set PRIMME SVDS parameters structure to the default values.

   Parameters:
      * **primme_svds** (*ptr*) -- (output) parameters structure.


primme_svds_set_method_f77
==========================

primme_svds_set_method_f77(method, methodStage1, methodStage2, primme_svds, ierr)

   Set PRIMME SVDS parameters to one of the preset configurations.

   Parameters:
      * **method** (*integer*) --

        (input) preset configuration to compute the singular triplets;
        one of

        * "PRIMME_SVDS_default", currently set as
          "PRIMME_SVDS_hybrid".

        * "PRIMME_SVDS_normalequations", compute the eigenvectors of
          A^*A or A A^*.

        * "PRIMME_SVDS_augmented", compute the eigenvectors of the
          augmented matrix, \left(\begin{array}{cc} 0 & A^* \\ A & 0
          \end{array}\right).

        * "PRIMME_SVDS_hybrid", start with
          "PRIMME_SVDS_normalequations"; use the resulting approximate
          singular vectors as initial vectors for
          "PRIMME_SVDS_augmented" if the required accuracy was not
          achieved.

      * **methodStage1** (*primme_preset_method*) -- (input) preset
        method to compute the eigenpairs at the first stage; see
        available values at "primme_set_method_f77()".

      * **methodStage2** (*primme_preset_method*) -- (input) preset
        method to compute the eigenpairs with the second stage of
        "PRIMME_SVDS_hybrid"; see available values at
        "primme_set_method_f77()".

      * **primme_svds** (*ptr*) -- (input/output) parameters
        structure.

      * **ierr** (*integer*) -- (output) if 0, successful; if
        negative, something went wrong.


primme_svds_display_params_f77
==============================

primme_svds_display_params_f77(primme_svds)

   Display all printable settings of "primme_svds" into the file
   descriptor "outputFile".

   Parameters:
      * **primme_svds** (*ptr*) -- (input) parameters structure.


primme_svds_free_f77
====================

primme_svds_free_f77(primme_svds)

   Free memory allocated by PRIMME SVDS and delete all values set.

   Parameters:
      * **primme_svds** (*ptr*) -- (input/output) parameters
        structure.


primme_svds_set_member_f77
==========================

primme_svds_set_member_f77(primme_svds, label, value)

   Set a value in some field of the parameter structure.

   Parameters:
      * **primme_svds** (*ptr*) -- (input) parameters structure.

      * **label** (*integer*) --

        field where to set value. One of:

           "PRIMME_SVDS_primme"
           "PRIMME_SVDS_primmeStage2"
           "PRIMME_SVDS_m"
           "PRIMME_SVDS_n"
           "PRIMME_SVDS_matrixMatvec"
           "PRIMME_SVDS_applyPreconditioner"
           "PRIMME_SVDS_numProcs"
           "PRIMME_SVDS_procID"
           "PRIMME_SVDS_mLocal"
           "PRIMME_SVDS_nLocal"
           "PRIMME_SVDS_commInfo"
           "PRIMME_SVDS_globalSumReal"
           "PRIMME_SVDS_numSvals"
           "PRIMME_SVDS_target"
           "PRIMME_SVDS_numTargetShifts"
           "PRIMME_SVDS_targetShifts"
           "PRIMME_SVDS_method"
           "PRIMME_SVDS_methodStage2"
           "PRIMME_SVDS_intWorkSize"
           "PRIMME_SVDS_realWorkSize"
           "PRIMME_SVDS_intWork"
           "PRIMME_SVDS_realWork"
           "PRIMME_SVDS_matrix"
           "PRIMME_SVDS_preconditioner"
           "PRIMME_SVDS_locking"
           "PRIMME_SVDS_numOrthoConst"
           "PRIMME_SVDS_aNorm"
           "PRIMME_SVDS_eps"
           "PRIMME_SVDS_precondition"
           "PRIMME_SVDS_initSize"
           "PRIMME_SVDS_maxBasisSize"
           "PRIMME_SVDS_maxBlockSize"
           "PRIMME_SVDS_maxMatvecs"
           "PRIMME_SVDS_iseed"
           "PRIMME_SVDS_printLevel"
           "PRIMME_SVDS_outputFile"
           "PRIMME_SVDS_stats_numOuterIterations"
           "PRIMME_SVDS_stats_numRestarts"
           "PRIMME_SVDS_stats_numMatvecs"
           "PRIMME_SVDS_stats_numPreconds"
           "PRIMME_SVDS_stats_elapsedTime"

      * **value** -- (input) value to set.

   Note: **Don't use** this function inside PRIMME SVDS's callback
     functions, e.g., "matrixMatvec" or "applyPreconditioner", or in
     functions called by these functions.


primme_svdstop_get_member_f77
=============================

primme_svdstop_get_member_f77(primme_svds, label, value)

   Get the value in some field of the parameter structure.

   Parameters:
      * **primme_svds** (*ptr*) -- (input) parameters structure.

      * **label** (*integer*) -- (input) field where to get value.
        One of the detailed in function
        "primmesvds_top_set_member_f77()".

      * **value** -- (output) value of the field.

   Note: **Don't use** this function inside PRIMME SVDS's callback
     functions, e.g., "matrixMatvec" or "applyPreconditioner", or in
     functions called by these functions. In those cases use
     "primme_svds_get_member_f77()".

   Note: When "label" is one of "PRIMME_SVDS_matrixMatvec",
     "PRIMME_SVDS_applyPreconditioner", "PRIMME_SVDS_commInfo",
     "PRIMME_SVDS_intWork", "PRIMME_SVDS_realWork",
     "PRIMME_SVDS_matrix" and "PRIMME_SVDS_preconditioner", the
     returned "value" is a C pointer ("void*"). Use Fortran pointer or
     other extensions to deal with it. For instance:

        use iso_c_binding
        MPI_Comm comm

        comm = MPI_COMM_WORLD
        call primme_svds_set_member_f77(primme_svds, PRIMME_SVDS_commInfo, comm)
        ...
        subroutine par_GlobalSumDouble(x,y,k,primme_svds)
        use iso_c_binding
        implicit none
        ...
        MPI_Comm, pointer :: comm
        type(c_ptr) :: pcomm

        call primme_svds_get_member_f77(primme_svds, PRIMME_SVDS_commInfo, pcomm)
        call c_f_pointer(pcomm, comm)
        call MPI_Allreduce(x,y,k,MPI_DOUBLE,MPI_SUM,comm,ierr)

     Most users would not need to retrieve these pointers in their
     programs.


primme_svds_get_member_f77
==========================

primme_svds_get_member_f77(primme_svds, label, value)

   Get the value in some field of the parameter structure.

   Parameters:
      * **primme_svds** (*ptr*) -- (input) parameters structure.

      * **label** (*integer*) -- (input) field where to get value.
        One of the detailed in function
        "primme_svdstop_set_member_f77()".

      * **value** -- (output) value of the field.

   Note: Use this function exclusively inside PRIMME SVDS's callback
     functions, e.g., "matrixMatvec" or "applyPreconditioner", or in
     functions called by these functions. Otherwise, e.g., from the
     main program, use the function "primme_svdstop_get_member_f77()".

   Note: When "label" is one of "PRIMME_SVDS_matrixMatvec",
     "PRIMME_SVDS_applyPreconditioner", "PRIMME_SVDS_commInfo",
     "PRIMME_SVDS_intWork", "PRIMME_SVDS_realWork",
     "PRIMME_SVDS_matrix" and "PRIMME_SVDS_preconditioner", the
     returned "value" is a C pointer ("void*"). Use Fortran pointer or
     other extensions to deal with it. For instance:

        use iso_c_binding
        MPI_Comm comm

        comm = MPI_COMM_WORLD
        call primme_svds_set_member_f77(primme_svds, PRIMME_SVDS_commInfo, comm)
        ...
        subroutine par_GlobalSumDouble(x,y,k,primme_svds)
        use iso_c_binding
        implicit none
        ...
        MPI_Comm, pointer :: comm
        type(c_ptr) :: pcomm

        call primme_svds_get_member_f77(primme_svds, PRIMME_SVDS_commInfo, pcomm)
        call c_f_pointer(pcomm, comm)
        call MPI_Allreduce(x,y,k,MPI_DOUBLE,MPI_SUM,comm,ierr)

     Most users would not need to retrieve these pointers in their
     programs.

Appendix
********


primme_svds_params
==================

primme_svds_params

      Structure to set the problem matrix and the solver options.

      PRIMME_INT m

         Number of rows of the matrix.

         Input/output:

               "primme_initialize()" sets this field to 0;
               this field is read by "dprimme()".

      PRIMME_INT n

         Number of columns of the matrix.

         Input/output:

               "primme_initialize()" sets this field to 0;
               this field is read by "dprimme()".

      void (*matrixMatvec)(void *x, PRIMME_INT ldx, void *y, PRIMME_INT ldy, int *blockSize, int *transpose, primme_svds_params *primme_svds, int *ierr)

         Block matrix-multivector multiplication, y = A x if
         "transpose" is zero, and y = A^*x otherwise.

         Parameters:
            * **x** -- input array.

            * **ldx** -- leading dimension of "x".

            * **y** -- output array.

            * **ldy** -- leading dimension of "y".

            * **blockSize** -- number of columns in "x" and "y".

            * **transpose** -- if non-zero, the transpose A should
              be applied.

            * **primme_svds** -- parameters structure.

            * **ierr** -- output error code; if it is set to non-
              zero, the current call to PRIMME will stop.

         If "transpose" is zero, then "x" and "y" are arrays of
         dimensions "nLocal" x "blockSize" and "mLocal" x "blockSize"
         respectively. Elsewhere they have dimensions "mLocal" x
         "blockSize" and "nLocal" x "blockSize". Both arrays are
         column-major (consecutive rows are consecutive in memory).

         The actual type of "x" and "y" depends on which function is
         being calling. For "dprimme_svds()", it is "double", for
         "zprimme_svds()" it is "PRIMME_COMPLEX_DOUBLE", for
         "sprimme_svds()" it is "float" and for "cprimme_svds()" it is
         "PRIMME_COMPLEX_FLOAT".

         Input/output:

               "primme_initialize()" sets this field to NULL;
               this field is read by "dprimme_svds()" and "zprimme_svds()".

         Note: Integer arguments are passed by reference to make
           easier the interface to other languages (like Fortran).

      void (*applyPreconditioner)(void *x, PRIMME_INT ldx, void *y, PRIMME_INT ldy, int *blockSize, int *mode, primme_svds_params *primme_svds, int *ierr)

         Block preconditioner-multivector application. Depending on
         "mode" it is expected an approximation of the inverse of

         * "primme_svds_op_AtA": y = A^*Ax - \sigma^2 I,

         * "primme_svds_op_AAt": y = AA^*x - \sigma^2 I,

         * "primme_svds_op_augmented": \left(\begin{array}{cc} 0 &
           A^* \\ A & 0 \end{array}\right) - \sigma I.

         Where \sigma is the current target (see "targetShifts") (for
         finding the smallest \sigma is zero).

         Parameters:
            * **x** -- input array.

            * **ldx** -- leading dimension of "x".

            * **y** -- output array.

            * **ldy** -- leading dimension of "y".

            * **blockSize** -- number of columns in "x" and "y".

            * **mode** -- one of "primme_svds_op_AtA",
              "primme_svds_op_AAt" or "primme_svds_op_augmented".

            * **primme_svds** -- parameters structure.

            * **ierr** -- output error code; if it is set to non-
              zero, the current call to PRIMME will stop.

         If "mode" is "primme_svds_op_AtA", then "x" and "y" are
         arrays of dimensions "nLocal" x "blockSize"; if mode is
         "primme_svds_op_AAt", they are "mLocal" x "blockSize"; and
         otherwise they are ("mLocal" + "nLocal") x "blockSize". Both
         arrays are column-major (consecutive rows are consecutive in
         memory).

         The actual type of "x" and "y" depends on which function is
         being calling. For "dprimme_svds()", it is "double", for
         "zprimme_svds()" it is "PRIMME_COMPLEX_DOUBLE", for
         "sprimme_svds()" it is "float" and for "cprimme_svds()" it is
         "PRIMME_COMPLEX_FLOAT".

         Input/output:

               "primme_initialize()" sets this field to NULL;
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      int numProcs

         Number of processes calling "dprimme_svds()" or
         "zprimme_svds()" in parallel.

         Input/output:

               "primme_initialize()" sets this field to 1;
               this field is read by "dprimme()" and "zprimme_svds()".

      int procID

         The identity of the local process within a parallel execution
         calling "dprimme_svds()" or "zprimme_svds()". Only the
         process with id 0 prints information.

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               "dprimme_svds()" sets this field to 0 if "numProcs" is 1;
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      PRIMME_INT mLocal

         Number of local rows on this process. The value depends on
         how the matrix and preconditioner is distributed along the
         processes.

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               "dprimme_svds()" sets this field to "m" if "numProcs" is 1;
               this field is read by "dprimme_svds()" and "zprimme_svds()".

         See also: "matrixMatvec" and "applyPreconditioner".

      PRIMME_INT nLocal

         Number of local columns on this process. The value depends on
         how the matrix and preconditioner is distributed along the
         processes.

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               "dprimme_svds()" sets this field to to "n" if "numProcs" is 1;
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      void *commInfo

         A pointer to whatever parallel environment structures needed.
         For example, with MPI, it could be a pointer to the MPI
         communicator. PRIMME does not use this. It is available for
         possible use in user functions defined in "matrixMatvec",
         "applyPreconditioner" and "globalSumReal".

         Input/output:

               "primme_svds_initialize()" sets this field to NULL;

      void (*globalSumReal)(double *sendBuf, double *recvBuf, int *count, primme_svds_params *primme_svds, int *ierr)

         Global sum reduction function. No need to set for sequential
         programs.

         Parameters:
            * **sendBuf** -- array of size "count" with the local
              input values.

            * **recvBuf** -- array of size "count" with the global
              output values so that the i-th element of recvBuf is the
              sum over all processes of the i-th element of "sendBuf".

            * **count** -- array size of "sendBuf" and "recvBuf".

            * **primme_svds** -- parameters structure.

            * **ierr** -- output error code; if it is set to non-
              zero, the current call to PRIMME will stop.

         The actual type of "sendBuf" and "recvBuf" depends on which
         function is being calling. For "dprimme_svds()" and
         "zprimme_svds()" it is "double", and for "sprimme_svds()" and
         "cprimme_svds()" it is "float". Note that "count" is the
         number of values of the actual type.

         Input/output:

               "primme_svds_initialize()" sets this field to an internal function;
               "dprimme_svds()" sets this field to an internal function if "numProcs" is 1 and "globalSumReal" is NULL;
               this field is read by "dprimme_svds()" and "zprimme_svds()".

         When MPI is used, this can be a simply wrapper to
         MPI_Allreduce() as shown below:

            void par_GlobalSumForDouble(void *sendBuf, void *recvBuf, int *count,
                                     primme_svds_params *primme_svds, int *ierr) {
               MPI_Comm communicator = *(MPI_Comm *) primme_svds->commInfo;
               if (MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM,
                             communicator) == MPI_SUCCESS) {
                  *ierr = 0;
               } else {
                  *ierr = 1;
               }
            }

            When calling :c:func:`sprimme_svds` and :c:func:`cprimme_svds` replace ``MPI_DOUBLE`` by ```MPI_FLOAT``.

      int numSvals

         Number of singular triplets wanted.

         Input/output:

               "primme_svds_initialize()" sets this field to 1;
               this field is read by "primme_svds_set_method()" (see *Preset Methods*) and "dprimme_svds()".

      primme_svds_target target

         Which singular values to find:

         "primme_svds_smallest"
            Smallest singular values; "targetShifts" is ignored.

         "primme_svds_largest"
            Largest singular values; "targetShifts" is ignored.

         "primme_svds_closest_abs"
            Closest in absolute value to the shifts in "targetShifts".

         Input/output:

               "primme_svds_initialize()" sets this field to "primme_svds_smallest";
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      int numTargetShifts

         Size of the array "targetShifts". Used only when "target" is
         "primme_svds_closest_abs". The default values is 0.

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      double *targetShifts

         Array of shifts, at least of size "numTargetShifts". Used
         only when "target" is "primme_svds_closest_abs".

         Singular values are computed in order so that the i-th
         singular value is the closest to the i-th shift. If
         "numTargetShifts" < "numSvals", the last shift given is used
         for all the remaining i's.

         Input/output:

               "primme_svds_initialize()" sets this field to NULL;
               this field is read by "dprimme_svds()" and "zprimme_svds()".

         Note: Eventually this is used by  "dprimme()" and
           "zprimme()". Please see considerations of "targetShifts".

      int printLevel

         The level of message reporting from the code. For now it
         controls the reporting level of the underneath eigensolvers.
         See "printLevel" in primme_params.

         All output is writen in "outputFile".

         Input/output:

               "primme_svds_initialize()" sets this field to 1;
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      double aNorm

         An estimate of the 2-norm of A, which is used in the default
         convergence criterion (see "eps").

         If "aNorm" is less than or equal to 0, the code uses the
         largest absolute Ritz value seen. On return, "aNorm" is then
         replaced with that value.

         Input/output:

               "primme_svds_initialize()" sets this field to 0.0;
               this field is read and written by "dprimme_svds()" and "zprimme_svds()".

      double eps

         A triplet is marked as converged when the 2-norm of the
         residual vectors is less than "eps" * "aNorm". The residual
         vectors are A v - \sigma u and A^* u - \sigma v for the
         triplet (u,\sigma,v).

         The default value is machine precision times 10^4.

         Input/output:

               "primme_svds_initialize()" sets this field to 0.0;
               this field is read and written by "dprimme_svds()" and "zprimme_svds()".

      FILE *outputFile

         Opened file to write down the output.

         Input/output:

               "primme_svds_initialize()" sets this field to the standard output;
               this field is read by "dprimme_svds()", "zprimme_svds()" and "primme_svds_display_params()"

      int locking

         If set to 1, the underneath eigensolvers will use hard
         locking. See "locking".

         Input/output:

               "primme_svds_initialize()" sets this field to -1;
               written by "primme_svds_set_method()" (see *Preset Methods*);
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      int initSize

            On input, the number of initial vector guesses provided in
            "svecs" argument in "dprimme_svds()" and "zprimme_svds()".

            On output, "initSize" holds the number of converged
            triplets. Without "locking" all "numSvals" approximations
            are in "svecs" but only the first "initSize" are
            converged.

            During execution, it holds the current number of converged
            triplets.

            Input/output:

                  "primme_svds_initialize()" sets this field to 0;
                  this field is read and written by "dprimme_svds()" and "zprimme_svds()".

         int numOrthoConst

            Number of vectors to be used as external orthogonalization
            constraints. The left and the right vector constraints are
            provided as input of the "svecs" argument in
            "sprimme_svds()" or other variant, and must be
            orthonormal.

            PRIMME SVDS finds new triplets orthogonal to these
            constraints (equivalent to solving the problem
            (I-UU^*)A(I-VV^*) where U and V are the given left and
            right constraint vectors). This is a handy feature if some
            singular triplets are already known, or for finding more
            triplets after a call to "dprimme_svds()" or
            "zprimme_svds()", possibly with different parameters (see
            an example in "TEST/exsvd_zseq.c").

            Input/output:

                  "primme_svds_initialize()" sets this field to 0;
                  this field is read by "dprimme_svds()" and "zprimme_svds()".

      int maxBasisSize

         The maximum basis size allowed in the main iteration. This
         has memory implications.

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               this field is read and written by "primme_svds_set_method()" (see *Preset Methods*);
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      int maxBlockSize

         The maximum block size the code will try to use.

         The user should set this based on the architecture specifics
         of the target computer, as well as any a priori knowledge of
         multiplicities. The code does *not* require that
         "maxBlockSize" > 1 to find multiple triplets. For some
         methods, keeping to 1 yields the best overall performance.

         Input/output:

               "primme_svds_initialize()" sets this field to 1;
               this field is read and written by "primme_svds_set_method()" (see *Preset Methods*);
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      PRIMME_INT maxMatvecs

         Maximum number of matrix vector multiplications
         (approximately half the number of preconditioning operations)
         that the code is allowed to perform before it exits.

         Input/output:

               "primme_svds_initialize()" sets this field to "INT_MAX";
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      int intWorkSize

         If "dprimme_svds()" or "zprimme_svds()" is called with all
         arguments as NULL except for "primme_svds_params" then it
         returns immediately with "intWorkSize" containing the size
         *in bytes* of the integer workspace that will be required by
         the parameters set.

         Otherwise if "intWorkSize" is not 0, it should be the size of
         the integer work array *in bytes* that the user provides in
         "intWork". If "intWorkSize" is 0, the code will allocate the
         required space, which can be freed later by calling
         "primme_svds_free()".

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               this field is read and written by "dprimme_svds()" and "zprimme_svds()".

      size_t realWorkSize

         If "dprimme_svds()" or "zprimme_svds()" is called with all
         arguments as NULL except for "primme_svds_params" then it
         returns immediately with "realWorkSize" containing the size
         *in bytes* of the real workspace that will be required by the
         parameters set.

         Otherwise if "realWorkSize" is not 0, it should be the size
         of the real work array *in bytes* that the user provides in
         "realWork". If "realWorkSize" is 0, the code will allocate
         the required space, which can be freed later by calling
         "primme_svds_free()".

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               this field is read and written by "dprimme_svds()" and "zprimme_svds()".

      int *intWork

         Integer work array.

         If NULL, the code will allocate its own workspace. If the
         provided space is not enough, the code will return the error
         code "-21".

         Input/output:

               "primme_svds_initialize()" sets this field to NULL;
               this field is read and written by "dprimme_svds()" and "zprimme_svds()".

      void *realWork

         Real work array.

         If NULL, the code will allocate its own workspace. If the
         provided space is not enough, the code will return the error
         code "-20".

         Input/output:

               "primme_svds_initialize()" sets this field to NULL;
               this field is read and written by "dprimme_svds()" and "zprimme_svds()".

      PRIMME_INT iseed

         The "PRIMME_INT iseed[4]" is an array with the seeds needed
         by the LAPACK dlarnv and zlarnv.

         The default value is an array with values -1, -1, -1 and -1.
         In that case, "iseed" is set based on the value of "procID"
         to avoid every parallel process generating the same sequence
         of pseudorandom numbers.

         Input/output:

               "primme_svds_initialize()" sets this field to "[-1, -1, -1, -1]";
               this field is read and written by "dprimme_svds()" and "zprimme_svds()".

      void *matrix

         This field may be used to pass any required information in
         the matrix-vector product "matrixMatvec".

         Input/output:

               "primme_svds_initialize()" sets this field to NULL;

      void *preconditioner

         This field may be used to pass any required information in
         the preconditioner function "applyPreconditioner".

         Input/output:

               "primme_svds_initialize()" sets this field to NULL;

      int precondition

         Set to 1 to use preconditioning. Make sure
         "applyPreconditioner" is not NULL then!

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               this field is read and written by "primme_svds_set_method()" (see *Preset Methods*);
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      primme_svds_op_operator method

         Select the equivalent eigenvalue problem that will be solved:

         * "primme_svds_op_AtA": A^*Ax = \sigma^2 x,

         * "primme_svds_op_AAt": AA^*x = \sigma^2 x,

         * "primme_svds_op_augmented": \left(\begin{array}{cc} 0 &
           A^* \\ A & 0 \end{array}\right) x = \sigma x.

         The options for this solver are stored in "primme".

         Input/output:

               "primme_svds_initialize()" sets this field to "primme_svds_op_none";
               this field is read and written by "primme_svds_set_method()" (see *Preset Methods*);
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      primme_svds_op_operator methodStage2

         Select the equivalent eigenvalue problem that will be solved
         to refine the solution. The allowed options are
         "primme_svds_op_none" to not refine the solution and
         "primme_svds_op_augmented" to refine the solution by solving
         the augmented problem with the current solution as the
         initial vectors. See "method".

         The options for this solver are stored in "primmeStage2".

         Input/output:

               "primme_svds_initialize()" sets this field to "primme_svds_op_none";
               this field is read and written by "primme_svds_set_method()" (see *Preset Methods*);
               this field is read by "dprimme_svds()" and "zprimme_svds()".

      primme_params primme

         Parameter structure storing the options for underneath
         eigensolver that will be called at the first stage. See
         "method".

         Input/output:

               "primme_svds_initialize()" initialize this structure;
               this field is read and written by "primme_svds_set_method()" (see *Preset Methods*);
               this field is read and written by "dprimme_svds()" and "zprimme_svds()".

      primme_params primmeStage2

         Parameter structure storing the options for underneath
         eigensolver that will be called at the second stage. See
         "methodStage2".

         Input/output:

               "primme_svds_initialize()" initialize this structure;
               this field is read and written by "primme_svds_set_method()" (see *Preset Methods*);
               this field is read and written by "dprimme_svds()" and "zprimme_svds()".

      PRIMME_INT stats.numOuterIterations

         Hold the number of outer iterations.

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               written by "dprimme_svds()" and "zprimme_svds()".

      PRIMME_INT stats.numRestarts

         Hold the number of restarts.

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               written by "dprimme_svds()" and "zprimme_svds()".

      PRIMME_INT stats.numMatvecs

         Hold how many vectors the operator in "matrixMatvec" has been
         applied on.

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               written by "dprimme_svds()" and "zprimme_svds()".

      PRIMME_INT stats.numPreconds

         Hold how many vectors the operator in "applyPreconditioner"
         has been applied on.

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               written by "dprimme_svds()" and "zprimme_svds()".

      double stats.elapsedTime

         Hold the wall clock time spent by the call to
         "dprimme_svds()" or "zprimme_svds()".

         Input/output:

               "primme_svds_initialize()" sets this field to 0;
               written by "dprimme_svds()" and "zprimme_svds()".


Error Codes
===========

The functions "dprimme_svds()" and "zprimme_svds()" return one of the
next values:

* 0: success,

* 1: reported only amount of required memory,

* -1: failed in allocating int or real workspace,

* -2: malloc failed in allocating a permutation integer array,

* -3: main_iter() encountered problem; the calling stack of the
  functions where the error occurred was printed in 'stderr',

* -4: "primme_svds" is NULL,

* -5: Wrong value for "m" or "n" or "mLocal" or "nLocal",

* -6: Wrong value for "numProcs",

* -7: "matrixMatvec" is not set,

* -8: "applyPreconditioner" is not set but "precondition" == 1 ,

* -9: "numProcs" >1 but "globalSumReal" is not set,

* -10: Wrong value for "numSvals", it's larger than min("m", "n"),

* -11: Wrong value for "numSvals", it's smaller than 1,

* -13: Wrong value for "target",

* -14: Wrong value for "method",

* -15: Not supported combination of method and "methodStage2",

* -16: Wrong value for "printLevel",

* -17: "svals" is not set,

* -18: "svecs" is not set,

* -19: "resNorms" is not set

* -20: not enough memory for "realWork"

* -21: not enough memory for "intWork"

* -100 up to -199: eigensolver error from first stage; see the value
  plus 100 in *Error Codes*.

* -200 up to -299: eigensolver error from second stage; see the
  value plus 200 in *Error Codes*.


Preset Methods
==============

primme_svds_preset_method

   primme_svds_default

      Set as "primme_svds_hybrid".

   primme_svds_normalequations

      Solve the equivalent eigenvalue problem A^*A V = \Sigma^2 V and
      computes U by normalizing the vectors AV. If "m" is smaller than
      "n", AA^* is solved instead.

      With "primme_svds_normalequations" "primme_svds_set_method()"
      sets "method" to "primme_svds_op_AtA" if "m" is larger or equal
      than "n", and to "primme_svds_op_AAt" otherwise; and
      "methodStage2" is set to "primme_svds_op_none".

   primme_svds_augmented

      Solve the equivalent eigenvalue problem \left(\begin{array}{cc}
      0 & A^* \\ A & 0 \end{array}\right) X = \sigma X with X =
      \left(\begin{array}{cc}V\\U\end{array}\right).

      With "primme_svds_augmented" "primme_svds_set_method()" sets
      "method" to "primme_svds_op_augmented" and "methodStage2" to
      "primme_svds_op_none".

   primme_svds_hybrid

      First solve the equivalent normal equations (see
      "primme_svds_normalequations") and then refine the solution
      solving the augmented problem (see "primme_svds_augmented").

      With "primme_svds_normalequations" "primme_svds_set_method()"
      sets "method" to "primme_svds_op_AtA" if "m" is larger or equal
      than "n", and to "primme_svds_op_AAt" otherwise; and
      "methodStage2" is set to "primme_svds_op_augmented".

Python Interface
****************

Primme.svds(A, k=6, ncv=None, tol=0, which='LM', v0=None, maxiter=None, return_singular_vectors=True, precAHA=None, precAAH=None, precAug=None, u0=None, locku0=None, lockv0=None, return_stats=False, maxBlockSize=0, method=None, methodStage1=None, methodStage2=None, **kargs)

   Compute k singular values and vectors for a sparse matrix.

   Parameters:
      * **A** (*{sparse matrix, LinearOperator}*) -- Array to
        compute the SVD on, of shape (M, N)

      * **k** (*int, optional*) -- Number of singular values and
        vectors to compute. Must be 1 <= k < min(A.shape).

      * **ncv** (*int, optional*) -- The maximum size of the basis

      * **tol** (*float, optional*) -- Tolerance for singular
        values. Zero (default) means machine precision.

      * **which** (*str ['LM' | 'SM'] or number, optional*) --

        Which *k* singular values to find:

           * 'LM' : largest singular values

           * 'SM' : smallest singular values

           * number : closest singular values to (referred as sigma
             later)

      * **u0** (*ndarray, optional*) --

        Left starting vectors for the iterations.

        Should be approximate left singular vectors. If only u0 or v0
        is provided, the other is computed.

      * **v0** (*ndarray, optional*) -- Right starting vectors for
        the iterations.

      * **maxiter** (*int, optional*) -- Maximum number of
        iterations.

      * **precAHA** (*{N x N matrix, array, sparse matrix,
        LinearOperator}, optional*) -- Approximate inverse of (A.H*A -
        sigma**2*I). If provided and M>N, it usually accelerates the
        convergence.

      * **precAAH** (*{M x M matrix, array, sparse matrix,
        LinearOperator}, optional*) -- Approximate inverse of (A*A.H -
        sigma**2*I). If provided and M<N, it usually accelerates the
        convergence.

      * **precAug** (*{(M+N) x (M+N) matrix, array, sparse matrix,
        LinearOperator}, optional*) -- Approximate inverse of
        ([zeros() A.H; zeros() A] - sigma*I). It usually accelerates
        the convergence if tol<dtype.eps**.5.

      * **locku0** (*ndarray, optional*) --

        Left orthogonal vector constrain.

        Seek singular triplets orthogonal to locku0 and lockv0. The
        provided vectors *should* be orthonormal. If only locku0 or
        lockv0 is provided, the other is computed. Useful to not
        converge some already computed solutions.

      * **lockv0** (*ndarray, optional*) -- Right orthogonal vector
        constrain. See locku0.

      * **maxBlockSize** (*int, optional*) -- Maximum number of
        vectors added at every iteration.

      * **report_stats** (*bool, optional*) -- If True, it is also
        returned extra information from PRIMME.

   Returns:
      * **u** (*ndarray, shape=(M, k), optional*) -- Unitary matrix
        having left singular vectors as columns. Returned if
        *return_singular_vectors* is True.

      * **s** (*ndarray, shape=(k,)*) -- The singular values.

      * **vt** (*ndarray, shape=(k, N), optional*) -- Unitary matrix
        having right singular vectors as rows. Returned if
        *return_singular_vectors* is True.

      * **stats** (*dict, optional (if return_stats)*) -- Extra
        information reported by PRIMME:

        * "numOuterIterations": number of outer iterations

        * "numRestarts": number of restarts

        * "numMatvecs": number of A*v

        * "numPreconds": number of OPinv*v

        * "elapsedTime": time that took

        * "rnorms" : ||A*v[i] - u[i]*s[i]||

        Returned if *return_stats* is True.

   See also:

     "Primme.eigsh()"
        eigenvalue decomposition for a sparse symmetrix/complex
        Hermitian matrix A

     "scipy.sparse.linalg.eigs()"
        eigenvalues and eigenvectors for a general (nonsymmetric)
        matrix A

   -[ Examples ]-

   >>> import Primme, scipy.sparse
   >>> A = scipy.sparse.spdiags(range(1, 11), [0], 100, 10) # sparse diag. rect. matrix
   >>> svecs_left, svals, svecs_right = Primme.svds(A, 3, tol=1e-6, which='SM')
   >>> svals # the three smallest singular values of A
   array([ 1.,  2.,  3.])

   >>> import Primme, scipy.sparse
   >>> A = scipy.sparse.rand(10000, 100, random_state=10)
   >>> prec = scipy.sparse.spdiags(np.reciprocal(A.multiply(A).sum(axis=0)),
   ...           [0], 100, 100) # square diag. preconditioner
   >>> svecs_left, svals, svecs_right = Primme.svds(A, 3, which=6.0, tol=1e-6, precAHA=prec)
   >>> ["%.5f" % x for x in svals.flat] # the three closest singular values of A to 0.5
   ['5.99871', '5.99057', '6.01065']
