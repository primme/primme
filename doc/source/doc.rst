
PRIMME: PReconditioned Iterative MultiMethod Eigensolver
========================================================

PRIMME, pronounced as *prime*, finds a number of eigenvalues and their corresponding eigenvectors of a 
real symmetric, or complex hermitian matrix A. Largest, smallest and interior 
eigenvalues are supported. Preconditioning can be used to accelerate 
convergence. 
PRIMME is written in C, but a complete Fortran 77 interface is also provided.
  
Making & Linking
----------------

`Make_flags` has the flags and compilers used to make `libprimme.a`. Set at minimum:

* `TOP`, path where the PRIMME directory is located.
* `CC`, compiler program like `gcc`, `clang` or `icc`. 

Then do this to generate `libprimme.a`::

    make lib

.. role:: ccode(code) 
   :language: c

C Library Interface
-------------------

The next functions are declared in ``primme.h``.

.. c:function:: void primme_initialize(primme_params *primme)

   Set PRIMME parameters structure to the default values.

   :param primme_params* primme: parameters structure.


.. c:function:: int primme_set_method(primme_preset_method method, primme_params *primme)

   Set PRIMME parameters to one of the preset configurations.

   :param primme_preset_method method: preset configuration.

   :param primme_params* primme: parameters structure.

   :return: if 0, successful; if negative, something went wrong.

.. c:type:: primme_preset_method

   Enumeration of preset configurations.

   .. c:member:: DYNAMIC

      Switches to the best method dynamically; currently, between
      :c:member:`primme_preset_method.JDQMR_ETol` and 
      :c:member:`primme_preset_method.GD_Olsen_plusK`.
      Set :c:member:`primme_params.dynamicMethodSwitch` to 1.

   .. c:member:: DEFAULT_MIN_TIME

      Currently set as :c:member:`primme_preset_method.JDQMR_ETol`; this method is usually the fastest if the cost of the matrix vector product is the order of the matrix dimension.

   .. c:member:: DEFAULT_MIN_MATVECS

      Currently set as :c:member:`primme_preset_method.GD_Olsen_plusK`; this method usually spent less matrix vector products than the others, so it's a good choice when this operation is expensive.

   .. c:member:: Arnoldi

      Arnoldi implemented à la Generalized Davidson.

   .. c:member:: GD

      Generalized Davidson.

   .. c:member:: GD_plusK

      GD with locally optimal restarting. See :c:member:`primme_params.restartingParams.maxPrevRetain`.

   .. c:member:: GD_Olsen_plusK

      GD+k and a cheap Olsen's Method (set
      :c:member:`primme_params.correctionParams.projectors.RightX` to 1 and
      :c:member:`primme_params.correctionParams.projectors.SkewX` to 0).
      See :c:member:`primme_params.correctionParams.projectors.SkewX`.

   .. c:member:: JD_Olsen_plusK

      GD+k and Olsen's Method (set
      :c:member:`primme_params.correctionParams.projectors.RightX` to 1 and
      :c:member:`primme_params.correctionParams.projectors.SkewX` to 1).
      See :c:member:`primme_params.correctionParams.projectors.SkewX`.

   .. c:member:: RQI

      (Accelerated) Rayleigh Quotient Iteration.

   .. c:member:: JDQR

      Jacobi-Davidson with fixed number of inner steps.
      See :c:member:`primme_params.correctionParams.maxInnerIterations`.

   .. c:member:: JDQMR

      Jacobi-Davidson with adaptive stopping criterion for inner Quasi Minimum Residual (QMR).
      See :c:member:`primme_params.correctionParams.convTest`.

   .. c:member:: JDQMR_ETol

      JDQMR but QMR stops after residual norm reduces by a 0.1 factor.
      See :c:member:`primme_params.correctionParams.convTest`.

   .. c:member:: SUBSPACE_ITERATION

      Subspace iteration.

   .. c:member:: LOBPCG_OrthoBasis

      LOBPCG, the basis size is set to the number of wanted eigenvalues
      :c:member:`primme_params.numEvals`.

   .. c:member:: LOBPCG_OrthoBasis_Window

      LOBPCG with sliding window of 
      :c:member:`primme_params.maxBlockSize` <
      :c:member:`primme_params.numEvals`.

.. c:function:: void primme_Free(primme_params *primme)

   Free memory allocated by PRIMME.

   :param primme_params* primme: parameters structure.

.. c:function:: int dprimme(double *evals, double *evecs, double *resNorms, primme_params *primme)

   Solve a real symmetric standard eigenproblems.

   :param double* evals: array at least of size :c:member:`primme_params.numEvals` to store the
      computed eigenvalues; all parallel calls return the same value in this array.

   :param double* resNorms: array at least of size :c:member:`primme_params.numEvals` to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.

   :param double* evecs: array at least of size :c:member:`primme_params.nLocal` times :c:member:`primme_params.numEvals`
      to store columnwise the (local part of the) computed eigenvectors.

   :param primme_params* primme: parameters structure.

   :return: error indicator:

      *  0: success.
      *  1: reported only amount of required memory.
      * -1: failed in allocating int or real workspace.
      * -2: malloc failed in allocating a permutation integer array.
      * -3: main_iter() encountered problem; the calling stack of thei
         functions where the error occurred was printed in stderr.
      * -4: if argument primme is NULL.
      * -5: if :c:member:`primme_params.n` <= 0 or :c:member:`primme_params.nLocal` <= 0.
      * -6: if :c:member:`primme_params.numProcs` < 1.
      * -7: if :c:member:`primme_params.matrixMatvec` is NULL.
      * -8: if :c:member:`primme_params.applyPreconditioner` is NULL and 
         :c:member:`primme_params.correctionParams.precondition` is not NULL.
      * -9: if :c:member:`primme_params.globalSumDouble` is NULL.
      * -10: if :c:member:`primme_params.numEvals` > :c:member:`primme_params.n`.
      * -11: if :c:member:`primme_params.numEvals` < 0.
      * -12: if :c:member:`primme_params.eps` > 0 and
         :c:member:`primme_params.eps` < machine precision.
      * -13: if :c:member:`primme_params.target` is not properly defined.
      * -14: if :c:member:`primme_params.target` is one of ``primme_closest_geq``,
         ``primme_closest_leq`` or ``primme_closest_abs`` but
         :c:member:`primme_params.numTargetShifts` <= 0 (no shifts).
      * -15: if :c:member:`primme_params.target` is one of ``primme_closest_geq``,
         ``primme_closest_leq`` or ``primme_closest_abs`` but
         :c:member:`primme_params.targetShifts` is NULL  (no shifts array).
      * -16: if :c:member:`primme_params.numOrthoConst` < 0 or
         :c:member:`primme_params.numOrthoConst` >= :c:member:`primme_params.n`.
         (no free dimensions left).
      * -17: if :c:member:`primme_params.maxBasisSize` < 2.
      * -18: if :c:member:`primme_params.minRestartSize` <= 0.
      * -19: if :c:member:`primme_params.maxBlockSize` <= 0.
      * -20: if :c:member:`primme_params.restartingParams.maxPrevRetain` < 0.
      * -21: if :c:member:`primme_params.restartingParams.scheme` is not one of
         `primme_thick` or `primme_dtr`.
      * -22: if :c:member:`primme_params.initSize` < 0
      * -23: if not :c:member:`primme_params.locking` and
         :c:member:`primme_params.initSize` > :c:member:`primme_params.maxBasisSize`.
      * -24: if :c:member:`primme_params.locking` and
         :c:member:`primme_params.initSize` > :c:member:`primme_params.numEvals`.
      * -25: if :c:member:`primme_params.restartingParams.maxPrevRetain` +
         :c:member:`primme_params.minRestartSize` >= :c:member:`primme_params.maxBasisSize`.
      * -26: if :c:member:`primme_params.minRestartSize` >= :c:member:`primme_params.n`.
      * -27: if :c:member:`primme_params.printLevel` < 0 or
         :c:member:`primme_params.printLevel` > 5.
      * -28: if :c:member:`primme_params.correctionParams.convTest` is not one of
         ``primme_full_LTolerance``, ``primme_decreasing_LTolerance``,
         ``primme_adaptive_ETolerance`` or ``primme_adaptive``.
      * -29: if :c:member:`primme_params.correctionParams.convTest` ==
         ``primme_decreasing_LTolerance`` and
      *   :c:member:`primme_params.correctionParams.relTolBase` <= 1.
      * -30: if evals is NULL, but not evecs and resNorms.
      * -31: if evecs is NULL, but not evals and resNorms.
      * -32: if resNorms is NULL, but not evecs and evals.

.. c:function:: zprimme(double *evals, Complex_Z *evecs, double *resNorms, primme_params *primme)

   Solve a Hermitian standard eigenproblems; see function :c:func:`dprimme`.

.. c:type:: primme_params

   Structure to set the problem matrices and eigensolver options.

   .. c:member:: int n

      Dimension of the matrix.

   .. c:member:: void (*matrixMatvec) (void *x, void *y, int *blockSize, primme_params *primme)

      Block matrix-multivector multiplication, :math:`y = A x` in solving :math:`A x = \lambda x` or :math:`A x = \lambda B x`.
   
      :param void* x:
      :param void* y: one dimensional array containing the ``blockSize`` vectors 
         packed one after the other (i.e., the leading dimension is the vector size), each of size :c:member:`primme_params.nLocal`.
         The real type is ``double*`` and ``Complex_Z*`` when called from :c:func:`dprimme` and :c:func:`zprimme` respectively.
      :param int* blockSize: number of vectors in x and y.
      :param primme_params* primme: parameters structure.

      .. note::

         Argument ``blockSize`` is passed by reference to make easier the interface to other
         languages (like Fortran).

   .. c:member:: void (*applyPreconditioner)(void *x, void *y, int *blockSize, struct primme_params *primme)

      Block preconditioner-multivector application, :math:`y = M^{-1}x` where :math:`M` is usually an approximation of :math:`A - \sigma I` or :math:`A - \sigma B` for finding eigenvalues close to :math:`\sigma`.
      The function follows the convention of :c:member:`primme_params.matrixMatvec`.

 
   .. c:member:: void (*massMatrixMatvec)(void *x, void *y, int *blockSize, struct primme_params *primme)

      Block matrix-multivector multiplication, :math:`y = B x` in solving :math:`A x = \lambda B x`.
      The function follows the convention of :c:member:`primme_params.matrixMatvec`.

      .. warning::

         Generalized eigenproblems not implemented in current version.
         This member is included for future compatibility.

   .. c:member:: int numProcs

      Number of processes calling in parallel to :c:func:`dprimme` or :c:func:`zprimme`.
      The default value is 1.

   .. c:member:: int procID

      The identity of the process that is calling in parallel to :c:func:`dprimme` or
      :c:func:`zprimme`.
      Only the process with id 0 prints information.
      The default value is 0.

   .. c:member:: int nLocal

      Number of local rows on this process.
      The default value is :c:member:`primme_params.n` if :c:member:`primme_params.numProcs` is 1.

   .. c:member:: void *commInfo

      A pointer to whatever parallel environment structures needed.
      For example, with MPI, it could be a pointer to the MPI communicator.
      PRIMME does not use this. It is available for possible use in 
      user functions defined in :c:member:`primme_params.matrixMatvec`,
      :c:member:`primme_params.applyPreconditioner`, :c:member:`primme_params.massMatrixMatvec` and
      :c:member:`primme_params.globalSumDouble`.
      The default values is NULL.

   .. c:member:: void (*globalSumDouble)(double *sendBuf, double *recvBuf, int *count, primme_params *primme)

      Global sum reduction function. 

      :param double* sendBuf: array of size count with the input local values.
      :param double* recvBuf: array of size count with the output global values
         so that i-th element of recvBuf is the sum over all processes of the i-th element
         of sendBuf.
      :param int* count: array size of sendBuf and recvBuf.
      :param primme_params* primme: parameters structure.

      The default value is NULL if :c:member:`primme_params.numProcs` is 1.
      When MPI this can be a simply wrapper to MPI_Allreduce().

      .. code:: c

         void par_GlobalSumDouble(void *sendBuf, void *recvBuf, int *count, 
                                  primme_params *primme) {
            MPI_Comm communicator = *(MPI_Comm *) primme->commInfo;
            MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM,
                          communicator);
         }

      .. note::

         Argument ``count`` is passed by reference to make easier the interface to other
         languages (like Fortran).

      .. note::

         The arguments ``sendBuf`` and ``recvBuf`` are always double arrays and ``count``
         is always the number of double elements in both arrays, even for :c:func:`zprimme`.


   .. c:member:: int numEvals

      Number of eigenvalues wanted.
      The default value is 1.

   .. c:member:: primme_target target

      Which eigenpairs to find.
      The default value is ``primme_smallest``.

      ``primme_smallest``
         smallest algebraic eigenvalues; :c:member:`primme_params.targetShifts` is ignored.

      ``primme_largest``
         largest algebraic eigenvalues; :c:member:`primme_params.targetShifts` is ignored.

      ``primme_closest_geq``
         closest to, but greater or equal than the shifts in :c:member:`primme_params.targetShifts`.

      ``primme_closest_leq``
         closest to, but less or equal than the shifts in :c:member:`primme_params.targetShifts`.

      ``primme_closest_abs``
         closest in absolute value to than the shifts in :c:member:`primme_params.targetShifts`.

      .. note::
         * If some shift is close to the lower (higher) end of the spectrum,
           use either ``primme_closest_geq`` (``primme_closest_leq``) or
           ``primme_closest_abs``.
         * ``primme_closest_leq`` and ``primme_closest_geq`` are more efficient
           than ``primme_closest_abs``.
 
   .. c:member:: int numTargetShifts
 
      Size of the array :c:member:`primme_params.targetShifts`.
      Used only when :c:member:`primme_params.target` is ``primme_closest_geq``,
      ``primme_closest_leq`` or ``primme_closest_abs``.
      The default values is 0.

   .. c:member:: double *targetShifts

      Array of shifts, at least of size :c:member:`primme_params.numTrargetShifts`.
      Used only when :c:member:`primme_params.target` is ``primme_closest_geq``,
      ``primme_closest_leq`` or ``primme_closest_abs``.
      The default values is NULL.

      The i-th shift (or the last one, if it is not given) is taken into account in
      finding the i-th eigenvalue.

      .. note::

         For code efficiency and robustness, the shifts should be ordered.
         Order them in ascending (descending) order for shifts closer
         to the lower (higher) end of the spectrum.

      .. c:member:: int printLevel

         The level of message reporting from the code.
         One of:
 
         * 0: silent.
         * 1: print some error messages when these occur.
         * 2: as 1, and info about targeted eigenpairs when they are marked as converged::
         
               #Converged $1 eval[ $2 ]= $3 norm $4 Mvecs $5 Time $7

           or locked::

               #Lock epair[ $1 ]= $3 norm $4 Mvecs $5 Time $7

         * 3: as 2, and info about targeted eigenpairs every outer iteration::
         
               OUT $6 conv $1 blk $8 MV $5 Sec $7 EV $3 |r| $4

           Also, if it is used the dynamic method, show JDQMR/GDk performance ratio and
           the current method in use.
         * 4: as 3, and info about targeted eigenpairs every inner iteration::
         
               INN MV $5 Sec $7 Eval $3 Lin|r| $9 EV|r| $4
         
         * 5: as 4, and verbose info about certain choices of the algorithm.
         
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
         :c:member:`primme_params.procID` 0 when :c:member:`primme_params.printLevel`
         is from 0 to 4.
         If :c:member:`primme_params.printLevel` is 5 output can be produced in any of
         the parallel calls.

      .. note::

         Convergence history for plotting may be produced simply by::

            grep OUT outpufile | awk '{print $8" "$14}' > out
            grep INN outpufile | awk '{print $3" "$11}' > inn

         Then in Matlab::

            plot(out(:,1),out(:,2),'bo');hold; plot(inn(:,1),inn(:,2),'r');

         Or in gnuplot::

            plot 'out' w lp, 'inn' w lp

   .. c:member:: double aNorm

      An estimate of norm of the matrix A that is used in the convergence criterion
      (see :c:member:`primme_params.eps`).
      If it is less or equal to 0, it is used the largest absolute Ritz value seen.
      And on return, it is replaced with that value.

      The default value is 0.

   .. c:member:: double eps

      An eigenpairs is marked as converged when the 2-norm of the residual is less
      than :c:member:`primme_params.eps` times :c:member:`primme_params.aNorm`.
      The residual vector is :math:`A x - \lambda x` or :math:`A x - \lambda B x`.

      The default value is :math:`10^{-12}`.
 
   .. c:member:: FILE *outputFile

      Opened file to write down the output.

      The default value is the standard output.

   .. c:member:: int dynamicMethodSwitch

      If this value is 1, it alternates dynamically between :c:member:`primme_preset_method.DEFAULT_MIN_TIME`
      and :c:member:`primme_preset_method.DEFAULT_MIN_MATVECS`, trying to identify the fastest method.

      On exit, it holds a recommended method for future runs on this problem:

      * -1: use :c:member:`primme_preset_method.DEFAULT_MIN_MATVECS` next time.
      * -2: use :c:member:`primme_preset_method.DEFAULT_MIN_TIME` next time.
      * -3: close call, use :c:member:`primme_preset_method.DYNAMIC` next time again.
      
      .. note::

         Even for expert users we do not recommend setting :c:member:`primme_params.dynamicMethodSwitch`
         directly, but through :c:func:`primme_set_method`.

      .. note::

         The code obtains timings by the ``gettimeofday`` Unix utility. If a cheaper, more
         accurate timer is available, modify the ``PRIMMESRC/COMMONSRC/wtime.c``

   .. c:member:: int locking

      If set to 1, hard locking will be used (locking converged eigenvectors
      out of the search basis). Otherwise the code will try to use soft
      locking (à la ARPACK), when large enough :c:member:`primme_params.minRestartSize` is available.

      The default depends on the method and the value of some options.

   .. c:member:: int initSize
 
      On input, the number of initial vector guesses provided in ``evecs`` argument in :c:func:`dprimme`
      or :c:func:`zprimme`.
      On output, the number of converged eigenpairs.
      During execution, it holds the current number of converged eigenpairs.
      If in addition locking is used, these are accessible in ``evals`` and ``evecs``.

      The default value is 0.
      
   .. c:member:: int numOrthoConst

      Number of external orthogonalization constraint vectors provided ``evecs`` argument in
      :c:func:`dprimme` or :c:func:`zprimme`.

      Then eigenvectors are found orthogonal to those constraints (equivalent to solving
      the problem with :math:`(I-YY^*)A(I-YY^*)` and :math:`(I-YY^*)B(I-YY^*)` instead
      where :math:`Y` are the given constraint vectors).
      This is a handy feature if some eigenvectors are already known, or 
      for finding more eigenvalues after a call to :c:func:`dprimme` or :c:func:`zprimme`.

      The default value is 0.

   .. c:member:: int maxBasisSize

      The maximum basis size allowed in the main iteration. This has memory
      implications.

      The default depends on method.

      .. note::

         For interior eigenvalues use a larger value than usual.

   .. c:member:: int minRestartSize

      Maximum Ritz vectors kept after restarting the basis.

      The default depends on :c:member:`primme_params.maxBasisSize`,
      :c:member:`primme_params.blockSize` and method.

   .. c:member:: int maxBlockSize
 
      The maximum block size the code will try to use.

      The user should set
      this based on the architecture specifics of the target computer, 
      as well as any a priori knowledge of multiplicities. The code does 
      *not* require to be greater than 1 to find multiple eigenvalues. For some 
      methods, keeping to 1 yields the best overall performance.

      The default value is 1.

      .. note::

         Inner iterations of QMR are not performed in a block fashion.
         Every correction equation from a block is solved independently.

   .. c:member:: int maxMatvecs

      Maximum number of matrix vector multiplications (approximately equal to 
      the number of preconditioning operations) that the code is allowed to 
      perform before it exits.

      The default value is ``INT_MAX``. 

   .. c:member:: int maxOuterIterations

      Maximum number of outer iterations that the code is allowed to perform 
      before it exits.

      The default value is ``INT_MAX``. 

   .. c:member:: int intWorkSize

      If :c:func:`dprimme` or :c:func:`zprimme` are called with all arguments as NULL
      but :c:type:`primme_params` then it has the size *in bytes* of the integer
      workspace that is required.

      Otherwise if not 0, it is the size of the integer work array *in bytes* that
      the user provides in :c:member:`primme_params.intWork`. If it is 0, the code
      will allocate the required space and should be freed by calling :c:func:`primme_Free`.

      The default value is 0.

   .. c:member:: long int realWorkSize

      If :c:func:`dprimme` or :c:func:`zprimme` are called with all arguments as NULL
      but :c:type:`primme_params` then it has the size *in bytes* of the real
      workspace that is required.

      Otherwise if not 0, it is the size of the real work array *in bytes* that
      the user provides in :c:member:`primme_params.realWork`. If it is 0, the code
      will allocate the required space and should be freed by calling :c:func:`primme_Free`.

      The default value is 0.

   .. c:member:: int *intWork

      Integer work array.

      If NULL, the code will allocate its own workspace. If the provided space is not
      enough, the code will free it and allocate a new space.

      On exit, the first element shows if a locking problem has occurred.
      Using locking for large :c:member:`primme_params.numEvals` may, in some rare cases,
      cause some pairs to be practically converged, in the sense that their components 
      are in the basis of ``evecs``. If this is the case, a Rayleigh Ritz on returned
      ``evecs`` would provide the accurate eigenvectors (see [4]_).

      The default value is NULL. 

   .. c:member:: void *realWork

      Real work array.

      If NULL, the code will allocate its own workspace. If the provided space is not
      enough, the code will free it and allocate a new space.

      The default value is NULL. 

   .. c:member:: int iseed

      The ``int iseed[4]`` is an array with the seeds needed by the LAPACK_ dlarnv and zlarnv.

      The default value is an array with values 1, 2, 3 and 5.

   .. c:member:: void *matrix

      This field may be used to pass any required information 
      in the matrix-vector product :c:member:`primme_params.matrixMatvec`.

      The default value is NULL.
      
   .. c:member:: void *preconditioner

      This field may be used to pass any required information 
      in the matrix-vector product :c:member:`primme_params.applyPreconditioner`.

      The default value is NULL.

   .. c:member:: double *ShiftsForPreconditioner

      Array of size ``blockSize`` provided during execution of :c:member:dprimme and :c:member:zprimme holding
      the shifts to be used (if needed) in the preconditioning operation.

      For example if the block size is 3,
      there will be an array of three shifts in :c:member:`primme_params.ShiftsForPreconditioner`.
      Then the user can invert a shifted preconditioner for each of the 
      block vectors :math:`(M-ShiftsForPreconditioner_i)^{-1} x_i`.
      Classical Davidson (diagonal) preconditioning is an example of this.
   
   .. c:member:: primme_restartscheme restartingParams.scheme

      Select a restarting strategy:

      * ``primme_thick``, Thick restarting. This is the most efficient and robust
        in the general case.
      * ``primme_dtr``, Dynamic thick restarting. Helpful without 
        preconditioning but it is expensive to implement.

      The default value is ``primme_thick``.

   .. c:member:: int restartingParams.maxPrevRetain

      Number of approximations from previous iteration to be retained
      after restart (see [2]_). The restart size is :c:member:`primme_params.minRestartSize`
      plus :c:member:`primme_params.restartingParams.maxPrevRetain`.

      The default value is 1.

   .. c:member:: int correctionParams.precondition

      Set to 1 to use preconditioning.
      Make sure :c:member:`primme_params.applyPreconditioner` is not NULL then!

      The default value is 0.

   .. c:member:: int correctionParams.robustShifts

      Set to 1 to use robust shifting. It tries to avoid stagnation and 
      missconvergence by providing as shifts in :c:member:`primme_params.ShiftsForPreconditioner`
      the Ritz values displaced by an approximation of the eigenvalue error.

      The default value depends on method.

   .. c:member:: int correctionParams.maxInnerIterations

      Control the maximum number of inner QMR iterations:

      * 0:  no inner iterations;
      * >0: perform at most that number of inner iterations per outer step;
      * <0: perform at most the rest of the remaining matrix-vector products
        up to reach :c:member:`primme_params.maxMatvecs`.

      The default value depends on method.

      See also :c:member:`primme_params.correctionParams.convTest`.

   .. c:member:: double correctionParams.relTolBase

      Parameter used when :c:member:`primme_params.correctionParams.convTest`
      is ``primme_decreasing_LTolerance``.

   .. c:member:: primme_convergencetest correctionParams.convTest

      Set how to stop the inner QMR method:

      * ``primme_full_LTolerance``: stop by iterations only;

      * ``primme_decreasing_LTolerance``, stop when
        :math:`\text{relTolBase}^{-\text{outIts}}` where outIts
        is the number of outer iterations and retTolBase is set in
        :c:member:`primme_params.correctionParams.relTolBase`;
        This is a legacy option from classical JDQR and we recommend
        **strongly** against its use.

      * ``primme_adaptive``, stop when the estimated eigenvalue residual
        has reached the required tolerance (based on Notay's JDCG).

      * ``primme_adaptive_ETolerance``, as ``primme_adaptive`` but also
        stopping when the estimated eigenvalue residual has reduced 10
        times.

      The default value depends on method.

      .. note::

         Avoid to set :c:member:`primme_params.correctionParams.maxInnerIterations` to -1 and :c:member:`primme_params.correctionParams.convTest` to ``primme_full_LTolerance``.

      See also :c:member:`primme_params.correctionParams.maxInnerIterations`.

   .. c:member:: int correctionParams.projectors.LeftQ
   .. c:member:: int correctionParams.projectors.LeftX
   .. c:member:: int correctionParams.projectors.RightQ
   .. c:member:: int correctionParams.projectors.RightX
   .. c:member:: int correctionParams.projectors.SkewQ
   .. c:member:: int correctionParams.projectors.SkewX

      Control the projectors involved in the computation of the correction
      appended to the basis every (outer) iteration.

      Given the current selected Ritz value :math:`\Lambda` and vectors
      :math:`X`, the residual associated vectors :math:`R=AX-X\Lambda`, the previous locked
      vectors :math:`Q` and the preconditioner :math:`M^{-1}`.
      The correction :math:`D` appended to the basis in GD
      (when :c:member:`primme_params.correctionParams.maxInnerIterations` is 0) is:

      +--------+-------+-----------------------------------------------------------------------+
      | RightX | SkewX |        :math:`D`                                                      |
      +========+=======+=======================================================================+
      |    0   |   0   | :math:`M^{-1}R` (Classic GD)                                          |
      +--------+-------+-----------------------------------------------------------------------+
      |    1   |   0   | :math:`M^{-1}(R-\Delta X)` (cheap Olsen's Method)                     |
      +--------+-------+-----------------------------------------------------------------------+
      |    1   |   1   | :math:`(I- M^{-1}X(X^*M^{-1}X)^{-1}X^*)M^{-1}R` (Olsen's Method)      |
      +--------+-------+-----------------------------------------------------------------------+
      |    0   |   1   | error                                                                 |
      +--------+-------+-----------------------------------------------------------------------+

      Where :math:`\Delta` is a diagonal matrix that :math:`\Delta_{i,i}` holds an estimation
      of the error of the approximate eigenvalue :math:`\Lambda_{i,i}`.
 
      The values of ``RightQ``, ``SkewQ``, ``LeftX`` and ``LeftQ`` are ignored.

      The correction :math:`D` in JD
      (when :c:member:`primme_params.correctionParams.maxInnerIterations` isn't 0)
      results from solving:

      .. math::

         P_Q^l P_X^l (A-\sigma I) P_X^r P_Q^r M^{-1} D' = -R, \ \ \  D = P_X^r P_Q^l M^{-1}D'.

      For ``LeftQ`` (and similarly for ``LeftX``):

      * 0: :math:`P_Q^l = I`;
      * 1: :math:`P_Q^l = I - QQ^*`.

      For ``RightQ`` and ``SkewQ`` (and similarly for ``RightX`` and ``SkewX``):

      +--------+-------+-------------------------------+
      | RightQ | SkewQ |        :math:`P_Q^r`          |
      +========+=======+===============================+
      |    0   |   0   | :math:`I`                     |
      +--------+-------+-------------------------------+
      |    1   |   0   | :math:`I - QQ^*`              |
      +--------+-------+-------------------------------+
      |    1   |   1   | :math:`I - KQ(Q^*KQ)^{-1}Q^*` |
      +--------+-------+-------------------------------+
      |    0   |   1   | error                         |
      +--------+-------+-------------------------------+

      The default value depends on method.

      See [3]_ for a study about different projector configuration in JD.

   .. c:member: int stats.numOuterIterations

      Hold the number of outer iterations. The value is available during execution and at the end.

   .. c:member: int stats.numRestarts

      Hold the number of restarts during execution and at the end.

   .. c:member: int stats.numMatvecs

      Hold how many times :c:member:`primme_params.matrixMatvec has been called.
      The value is available during execution and at the end.

   .. c:member: int stats.elapsedTime

      Hold the wall clock time spent by the call to :c:func:`dprimme` or :c:func:`zprimme`.
      The value is available at the end of the execution.


..   struct stackTraceNode *stackTrace;            (OUTPUT)
..
.. Struct with the following members. If an error occurs the function
.. primme_PrintStackTrace(primme) prints the calling stack from top to the 
..       function that caused the error. Nothing to set.
..
.. int callingFunction;
..    int failedFunction;
..    int errorCode;
..    int lineNumber;
..    char fileName[PRIMME_MAX_NAME_LENGTH];
..    struct stackTraceNode *nextNode;
..
..

FORTRAN Library Interface
-------------------------

.. c:type:: ptr

   Fortran datatype with the same size as a pointer.
   Use ``integer*4`` when compiling in 32 bits and ``integer*8`` in 64 bits.

.. c:function:: primme_initialize_f77(primme)

   Set PRIMME parameters structure to the default values.

   :param ptr primme: (output) parameters structure.

.. c:function:: primme_set_method_f77(method, primme, ierr)

   Set PRIMME parameters to one of the preset configurations.

   :param integer method: (input) preset configuration. One of:

      * ``PRIMMEF77_DYNAMIC``
      * ``PRIMMEF77_DEFAULT_MIN_TIME``
      * ``PRIMMEF77_DEFAULT_MIN_MATVECS``
      * ``PRIMMEF77_Arnoldi``
      * ``PRIMMEF77_GD``
      * ``PRIMMEF77_GD_plusK``
      * ``PRIMMEF77_GD_Olsen_plusK``
      * ``PRIMMEF77_JD_Olsen_plusK``
      * ``PRIMMEF77_RQI``
      * ``PRIMMEF77_JDQR``
      * ``PRIMMEF77_JDQMR``
      * ``PRIMMEF77_JDQMR_ETol``
      * ``PRIMMEF77_SUBSPACE_ITERATION``
      * ``PRIMMEF77_LOBPCG_OrthoBasis``
      * ``PRIMMEF77_LOBPCG_OrthoBasis_Window``

      See :c:type:`primme_preset_method`.

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) if 0, successful; if negative, something went wrong.

.. c:function:: primme_Free_f77(primme)

   Free memory allocated by PRIMME.

   :param ptr primme: parameters structure.

.. c:function:: dprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a real symmetric standard eigenproblems.

   :param double precision evals(*): (output) array at least of size :c:member:`primme_params.numEvals` to store the
      computed eigenvalues; all parallel calls return the same value in this array.

   :param double precision resNorms(*): (output) array at least of size :c:member:`primme_params.numEvals` to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.

   :param double precision evecs(*): (input/output) array at least of size :c:member:`primme_params.nLocal` times :c:member:`primme_params.numEvals`
      to store columnwise the (local part of the) computed eigenvectors.

   :param ptr primme: parameters structure.

   :param integer ierr: (output) error indicator:

      *  0: success.
      *  1: reported only amount of required memory.
      * -1: failed in allocating int or real workspace.
      * -2: malloc failed in allocating a permutation integer array.
      * -3: main_iter() encountered problem; the calling stack of thei
         functions where the error occurred was printed in stderr.
      * -4: if argument primme is NULL.
      * -5: if :c:member:`primme_params.n` <= 0 or :c:member:`primme_params.nLocal` <= 0.
      * -6: if :c:member:`primme_params.numProcs` < 1.
      * -7: if :c:member:`primme_params.matrixMatvec` is NULL.
      * -8: if :c:member:`primme_params.applyPreconditioner` is NULL and 
         :c:member:`primme_params.correctionParams.precondition` is not NULL.
      * -9: if :c:member:`primme_params.globalSumDouble` is NULL.
      * -10: if :c:member:`primme_params.numEvals` > :c:member:`primme_params.n`.
      * -11: if :c:member:`primme_params.numEvals` < 0.
      * -12: if :c:member:`primme_params.eps` > 0 and
         :c:member:`primme_params.eps` < machine precision.
      * -13: if :c:member:`primme_params.target` is not properly defined.
      * -14: if :c:member:`primme_params.target` is one of ``primme_closest_geq``,
         ``primme_closest_leq`` or ``primme_closest_abs`` but
         :c:member:`primme_params.numTargetShifts` <= 0 (no shifts).
      * -15: if :c:member:`primme_params.target` is one of ``primme_closest_geq``,
         ``primme_closest_leq`` or ``primme_closest_abs`` but
         :c:member:`primme_params.targetShifts` is NULL  (no shifts array).
      * -16: if :c:member:`primme_params.numOrthoConst` < 0 or
         :c:member:`primme_params.numOrthoConst` >= :c:member:`primme_params.n`.
         (no free dimensions left).
      * -17: if :c:member:`primme_params.maxBasisSize` < 2.
      * -18: if :c:member:`primme_params.minRestartSize` <= 0.
      * -19: if :c:member:`primme_params.maxBlockSize` <= 0.
      * -20: if :c:member:`primme_params.restartingParams.maxPrevRetain` < 0.
      * -21: if :c:member:`primme_params.restartingParams.scheme` is not one of
         `primme_thick` or `primme_dtr`.
      * -22: if :c:member:`primme_params.initSize` < 0
      * -23: if not :c:member:`primme_params.locking` and
         :c:member:`primme_params.initSize` > :c:member:`primme_params.maxBasisSize`.
      * -24: if :c:member:`primme_params.locking` and
         :c:member:`primme_params.initSize` > :c:member:`primme_params.numEvals`.
      * -25: if :c:member:`primme_params.restartingParams.maxPrevRetain` +
         :c:member:`primme_params.minRestartSize` >= :c:member:`primme_params.maxBasisSize`.
      * -26: if :c:member:`primme_params.minRestartSize` >= :c:member:`primme_params.n`.
      * -27: if :c:member:`primme_params.printLevel` < 0 or
         :c:member:`primme_params.printLevel` > 5.
      * -28: if :c:member:`primme_params.correctionParams.convTest` is not one of
         ``primme_full_LTolerance``, ``primme_decreasing_LTolerance``,
         ``primme_adaptive_ETolerance`` or ``primme_adaptive``.
      * -29: if :c:member:`primme_params.correctionParams.convTest` ==
         ``primme_decreasing_LTolerance`` and
      *   :c:member:`primme_params.correctionParams.relTolBase` <= 1.
      * -30: if evals is NULL, but not evecs and resNorms.
      * -31: if evecs is NULL, but not evals and resNorms.
      * -32: if resNorms is NULL, but not evecs and evals.

.. c:function:: zprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblems. The arguments have the
   same meaning like in function :c:func:`dprimme_f77`.

   :param double precision evals(*): (output) 

   :param double precision resNorms(*): (output)

   :param complex double precision evecs(*): (input/output) 

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) error indicator.

.. c:function:: primmetop_set_member_f77(primme, label, value)

   Set a value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: field where to set value. One of:

      * ``PRIMMEF77_n``,                                     in field :c:member:`primme_params.n`.
      * ``PRIMMEF77_matrixMatvec``,                          in field :c:member:`primme_params.matrixMatvec`.
      * ``PRIMMEF77_applyPreconditioner``,                   in field :c:member:`primme_params.applyPreconditioner`.
      * ``PRIMMEF77_numProcs``,                              in field :c:member:`primme_params.numProcs`.
      * ``PRIMMEF77_procID``,                                in field :c:member:`primme_params.procID`.
      * ``PRIMMEF77_commInfo``,                              in field :c:member:`primme_params.commInfo`.
      * ``PRIMMEF77_nLocal``,                                in field :c:member:`primme_params.nLocal`.
      * ``PRIMMEF77_globalSumDouble``,                       in field :c:member:`primme_params.globalSumDouble`.
      * ``PRIMMEF77_numEvals``,                              in field :c:member:`primme_params.numEvals`.
      * ``PRIMMEF77_target``,                                in field :c:member:`primme_params.target`.
      * ``PRIMMEF77_numTargetShifts``,                       in field :c:member:`primme_params.numTargetShifts`.
      * ``PRIMMEF77_targetShifts``,                          in field :c:member:`primme_params.targetShifts`.
      * ``PRIMMEF77_locking``,                               in field :c:member:`primme_params.locking`.
      * ``PRIMMEF77_initSize``,                              in field :c:member:`primme_params.initSize`.
      * ``PRIMMEF77_numOrthoConst``,                         in field :c:member:`primme_params.numOrthoConst`.
      * ``PRIMMEF77_maxBasisSize``,                          in field :c:member:`primme_params.maxBasisSize`.
      * ``PRIMMEF77_minRestartSize``,                        in field :c:member:`primme_params.minRestartSize`.
      * ``PRIMMEF77_maxBlockSize``,                          in field :c:member:`primme_params.maxBlockSize`.
      * ``PRIMMEF77_maxMatvecs``,                            in field :c:member:`primme_params.maxMatvecs`.
      * ``PRIMMEF77_maxOuterIterations``,                    in field :c:member:`primme_params.maxOuterIterations`.
      * ``PRIMMEF77_intWorkSize``,                           in field :c:member:`primme_params.intWorkSize`.
      * ``PRIMMEF77_realWorkSize``,                          in field :c:member:`primme_params.realWorkSize`.
      * ``PRIMMEF77_iseed``,                                 in field :c:member:`primme_params.iseed`.
      * ``PRIMMEF77_intWork``,                               in field :c:member:`primme_params.intWork`.
      * ``PRIMMEF77_realWork``,                              in field :c:member:`primme_params.realWork`.
      * ``PRIMMEF77_aNorm``,                                 in field :c:member:`primme_params.aNorm`.
      * ``PRIMMEF77_eps``,                                   in field :c:member:`primme_params.eps`.
      * ``PRIMMEF77_printLevel``,                            in field :c:member:`primme_params.printLevel`.
      * ``PRIMMEF77_outputFile``,                            in field :c:member:`primme_params.outputFile`.
      * ``PRIMMEF77_matrix``,                                in field :c:member:`primme_params.matrix`.
      * ``PRIMMEF77_preconditioner``,                        in field :c:member:`primme_params.preconditioner`.
      * ``PRIMMEF77_restartingParams_scheme``,               in field :c:member:`primme_params.restartingParams.scheme`.
      * ``PRIMMEF77_restartingParams_maxPrevRetain``,        in field :c:member:`primme_params.restartingParams.maxPrevRetain`.
      * ``PRIMMEF77_correctionParams_precondition``,         in field :c:member:`primme_params.correctionParams.precondition`.
      * ``PRIMMEF77_correctionParams_robustShifts``,         in field :c:member:`primme_params.correctionParams.robustShifts`.
      * ``PRIMMEF77_correctionParams_maxInnerIterations``,   in field :c:member:`primme_params.correctionParams.maxInnerIterations`.
      * ``PRIMMEF77_correctionParams_projectors_LeftQ``,     in field :c:member:`primme_params.correctionParams.projectors_LeftQ`.
      * ``PRIMMEF77_correctionParams_projectors_LeftX``,     in field :c:member:`primme_params.correctionParams.projectors_LeftX`.
      * ``PRIMMEF77_correctionParams_projectors_RightQ``,    in field :c:member:`primme_params.correctionParams.projectors_RightQ`.
      * ``PRIMMEF77_correctionParams_projectors_RightX``,    in field :c:member:`primme_params.correctionParams.projectors_RightX`.
      * ``PRIMMEF77_correctionParams_projectors_SkewQ``,     in field :c:member:`primme_params.correctionParams.projectors_SkewQ`.
      * ``PRIMMEF77_correctionParams_projectors_SkewX``,     in field :c:member:`primme_params.correctionParams.projectors_SkewX`.
      * ``PRIMMEF77_correctionParams_convTest``,             in field :c:member:`primme_params.correctionParams.convTest`.
      * ``PRIMMEF77_correctionParams_relTolBase``,           in field :c:member:`primme_params.correctionParams.relTolBase`.
      * ``PRIMMEF77_stats_numOuterIterations``,              in field :c:member:`primme_params.stats.numOuterIterations`.
      * ``PRIMMEF77_stats_numRestarts``,                     in field :c:member:`primme_params.stats.numRestarts`.
      * ``PRIMMEF77_stats_numMatvecs``,                      in field :c:member:`primme_params.stats.numMatvecs`.
      * ``PRIMMEF77_stats_numPreconds``,                     in field :c:member:`primme_params.stats.numPreconds`.
      * ``PRIMMEF77_stats_elapsedTime``,                     in field :c:member:`primme_params.stats.elapsedTime`.
      * ``PRIMMEF77_dynamicMethodSwitch``,                   in field :c:member:`primme_params.dynamicMethodSwitch`.
      * ``PRIMMEF77_massMatrixMatvec``,                      in field :c:member:`primme_params.massMatrixMatvec`.

   :param value: (input) value to set.

.. c:function:: primmetop_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: (input) field where to get value. One of
      the detailed in function :c:func:`primmetop_set_member_f77`.

   :param value: (output) value of the field.

.. c:function:: primmetop_get_prec_shift_f77(primme, index, value)

   Get the value in some position of the array :c:member:`primme_params.ShiftsForPreconditioner`.

   :param ptr primme: (input) parameters structure.

   :param integer index: (input) position of the array; the first position is 1.

   :param value: (output) value of the array at that position.

.. c:function:: primme_set_member_f77(primme, label, value)

   Set a value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: field where to set value. One of:

      * ``PRIMMEF77_n``,                                     in field :c:member:`primme_params.n`.
      * ``PRIMMEF77_matrixMatvec``,                          in field :c:member:`primme_params.matrixMatvec`.
      * ``PRIMMEF77_applyPreconditioner``,                   in field :c:member:`primme_params.applyPreconditioner`.
      * ``PRIMMEF77_numProcs``,                              in field :c:member:`primme_params.numProcs`.
      * ``PRIMMEF77_procID``,                                in field :c:member:`primme_params.procID`.
      * ``PRIMMEF77_commInfo``,                              in field :c:member:`primme_params.commInfo`.
      * ``PRIMMEF77_nLocal``,                                in field :c:member:`primme_params.nLocal`.
      * ``PRIMMEF77_globalSumDouble``,                       in field :c:member:`primme_params.globalSumDouble`.
      * ``PRIMMEF77_numEvals``,                              in field :c:member:`primme_params.numEvals`.
      * ``PRIMMEF77_target``,                                in field :c:member:`primme_params.target`.
      * ``PRIMMEF77_numTargetShifts``,                       in field :c:member:`primme_params.numTargetShifts`.
      * ``PRIMMEF77_targetShifts``,                          in field :c:member:`primme_params.targetShifts`.
      * ``PRIMMEF77_locking``,                               in field :c:member:`primme_params.locking`.
      * ``PRIMMEF77_initSize``,                              in field :c:member:`primme_params.initSize`.
      * ``PRIMMEF77_numOrthoConst``,                         in field :c:member:`primme_params.numOrthoConst`.
      * ``PRIMMEF77_maxBasisSize``,                          in field :c:member:`primme_params.maxBasisSize`.
      * ``PRIMMEF77_minRestartSize``,                        in field :c:member:`primme_params.minRestartSize`.
      * ``PRIMMEF77_maxBlockSize``,                          in field :c:member:`primme_params.maxBlockSize`.
      * ``PRIMMEF77_maxMatvecs``,                            in field :c:member:`primme_params.maxMatvecs`.
      * ``PRIMMEF77_maxOuterIterations``,                    in field :c:member:`primme_params.maxOuterIterations`.
      * ``PRIMMEF77_intWorkSize``,                           in field :c:member:`primme_params.intWorkSize`.
      * ``PRIMMEF77_realWorkSize``,                          in field :c:member:`primme_params.realWorkSize`.
      * ``PRIMMEF77_iseed``,                                 in field :c:member:`primme_params.iseed`.
      * ``PRIMMEF77_intWork``,                               in field :c:member:`primme_params.intWork`.
      * ``PRIMMEF77_realWork``,                              in field :c:member:`primme_params.realWork`.
      * ``PRIMMEF77_aNorm``,                                 in field :c:member:`primme_params.aNorm`.
      * ``PRIMMEF77_eps``,                                   in field :c:member:`primme_params.eps`.
      * ``PRIMMEF77_printLevel``,                            in field :c:member:`primme_params.printLevel`.
      * ``PRIMMEF77_outputFile``,                            in field :c:member:`primme_params.outputFile`.
      * ``PRIMMEF77_matrix``,                                in field :c:member:`primme_params.matrix`.
      * ``PRIMMEF77_preconditioner``,                        in field :c:member:`primme_params.preconditioner`.
      * ``PRIMMEF77_restartingParams_scheme``,               in field :c:member:`primme_params.restartingParams.scheme`.
      * ``PRIMMEF77_restartingParams_maxPrevRetain``,        in field :c:member:`primme_params.restartingParams.maxPrevRetain`.
      * ``PRIMMEF77_correctionParams_precondition``,         in field :c:member:`primme_params.correctionParams.precondition`.
      * ``PRIMMEF77_correctionParams_robustShifts``,         in field :c:member:`primme_params.correctionParams.robustShifts`.
      * ``PRIMMEF77_correctionParams_maxInnerIterations``,   in field :c:member:`primme_params.correctionParams.maxInnerIterations`.
      * ``PRIMMEF77_correctionParams_projectors_LeftQ``,     in field :c:member:`primme_params.correctionParams.projectors_LeftQ`.
      * ``PRIMMEF77_correctionParams_projectors_LeftX``,     in field :c:member:`primme_params.correctionParams.projectors_LeftX`.
      * ``PRIMMEF77_correctionParams_projectors_RightQ``,    in field :c:member:`primme_params.correctionParams.projectors_RightQ`.
      * ``PRIMMEF77_correctionParams_projectors_RightX``,    in field :c:member:`primme_params.correctionParams.projectors_RightX`.
      * ``PRIMMEF77_correctionParams_projectors_SkewQ``,     in field :c:member:`primme_params.correctionParams.projectors_SkewQ`.
      * ``PRIMMEF77_correctionParams_projectors_SkewX``,     in field :c:member:`primme_params.correctionParams.projectors_SkewX`.
      * ``PRIMMEF77_correctionParams_convTest``,             in field :c:member:`primme_params.correctionParams.convTest`.
      * ``PRIMMEF77_correctionParams_relTolBase``,           in field :c:member:`primme_params.correctionParams.relTolBase`.
      * ``PRIMMEF77_stats_numOuterIterations``,              in field :c:member:`primme_params.stats.numOuterIterations`.
      * ``PRIMMEF77_stats_numRestarts``,                     in field :c:member:`primme_params.stats.numRestarts`.
      * ``PRIMMEF77_stats_numMatvecs``,                      in field :c:member:`primme_params.stats.numMatvecs`.
      * ``PRIMMEF77_stats_numPreconds``,                     in field :c:member:`primme_params.stats.numPreconds`.
      * ``PRIMMEF77_stats_elapsedTime``,                     in field :c:member:`primme_params.stats.elapsedTime`.
      * ``PRIMMEF77_dynamicMethodSwitch``,                   in field :c:member:`primme_params.dynamicMethodSwitch`.
      * ``PRIMMEF77_massMatrixMatvec``,                      in field :c:member:`primme_params.massMatrixMatvec`.

   :param value: (input) value to set.

   .. note::

      Use this function exclusively inside the function
      :c:member:`primme_params.matrixMatvec`,
      :c:member:`primme_params.massMatrixMatvec`, or
      :c:member:`primme_params.applyPreconditioner`.
      Otherwise use the function :c:func:`primmetop_set_member_f77`.

.. c:function:: primme_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: (input) field where to get value. One of
      the detailed in function :c:func:`primmetop_set_member_f77`.

   :param value: (output) value of the field.

   .. note::

      Use this function exclusively inside the function
      :c:member:`primme_params.matrixMatvec`,
      :c:member:`primme_params.massMatrixMatvec`, or
      :c:member:`primme_params.applyPreconditioner`.
      Otherwise use the function :c:func:`primmetop_get_member_f77`.

.. c:function:: primme_get_prec_shift_f77(primme, index, value)

   Get the value in some position of the array :c:member:`primme_params.ShiftsForPreconditioner`.

   :param ptr primme: (input) parameters structure.

   :param integer index: (input) position of the array; the first position is 1.

   :param value: (output) value of the array at that position.

   .. note::

      Use this function exclusively inside the function
      :c:member:`primme_params.matrixMatvec`,
      :c:member:`primme_params.massMatrixMatvec`, or
      :c:member:`primme_params.applyPreconditioner`.
      Otherwise use the function :c:func:`primmetop_get_prec_shift_f77`.

Citing this code 
---------------- 

Please cite:

.. [1] A. Stathopoulos and J. R. McCombs PRIMME: *PReconditioned Iterative
   MultiMethod Eigensolver: Methods and software description*, ACM
   Transaction on Mathematical Software Vol. 37, No. 2, (2010),
   21:1-21:30.

More information on the algorithms and research that led to this
software can be found in the rest of the papers. The work has been
supported by a number of grants from the National Science Foundation.

.. [2] A. Stathopoulos, *Nearly optimal preconditioned methods for hermitian
   eigenproblems under limited memory. Part I: Seeking one eigenvalue*, SIAM
   J. Sci. Comput., Vol. 29, No. 2, (2007), 481--514.

.. [3] A. Stathopoulos and J. R. McCombs, *Nearly optimal preconditioned
   methods for hermitian eigenproblems under limited memory. Part II:
   Seeking many eigenvalues*, SIAM J. Sci. Comput., Vol. 29, No. 5, (2007),
   2162-2188.

.. [4] J. R. McCombs and A. Stathopoulos, *Iterative Validation of
   Eigensolvers: A Scheme for Improving the Reliability of Hermitian
   Eigenvalue Solvers*, SIAM J. Sci. Comput., Vol. 28, No. 6, (2006),
   2337-2358.

.. [5] A. Stathopoulos, *Locking issues for finding a large number of eigenvectors
   of hermitian matrices*, Tech Report: WM-CS-2005-03, July, 2005.

License Information
-------------------

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


Contact Information 
-------------------

For reporting bugs or questions about functionality contact `Andreas Stathopoulos`_

.. _`Andreas Stathopoulos`: http://www.cs.wm.edu/~andreas/
.. _LAPACK: http://www.netlib.org/lapack/
