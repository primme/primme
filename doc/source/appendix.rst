
Appendix
--------

primme_params
"""""""""""""

.. c:type:: primme_params

   Structure to set the problem matrices and eigensolver options.

   .. c:member:: int n

      Dimension of the matrix.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: void (*matrixMatvec) (void *x, void *y, int *blockSize, primme_params *primme)

      Block matrix-multivector multiplication, :math:`y = A x` in solving :math:`A x = \lambda x` or :math:`A x = \lambda B x`.
   
      :param x: one dimensional array containing the ``blockSize`` vectors 
         packed one after the other (i.e., the leading dimension is the vector size), each of size |nLocal|.
         The real type is ``double*`` and ``Complex_Z*`` when called from :c:func:`dprimme` and :c:func:`zprimme` respectively.
      :param y: one dimensional array containing the ``blockSize`` vectors 
         packed one after the other (i.e., the leading dimension is the vector size), each of size |nLocal|.
         The real type is ``double*`` and ``Complex_Z*`` when called from :c:func:`dprimme` and :c:func:`zprimme` respectively.
      :param blockSize: number of vectors in x and y.
      :param primme: parameters structure.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme`.

      .. note::

         Argument ``blockSize`` is passed by reference to make easier the interface to other
         languages (like Fortran).

   .. c:member:: void (*applyPreconditioner)(void *x, void *y, int *blockSize, struct primme_params *primme)

      Block preconditioner-multivector application, :math:`y = M^{-1}x` where :math:`M` is usually an approximation of :math:`A - \sigma I` or :math:`A - \sigma B` for finding eigenvalues close to :math:`\sigma`.
      The function follows the convention of |matrixMatvec|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme`.
 
   .. c:member:: void (*massMatrixMatvec)(void *x, void *y, int *blockSize, struct primme_params *primme)

      Block matrix-multivector multiplication, :math:`y = B x` in solving :math:`A x = \lambda B x`.
      The function follows the convention of |matrixMatvec|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme`.

      .. warning::

         Generalized eigenproblems not implemented in current version.
         This member is included for future compatibility.

   .. c:member:: int numProcs

      Number of processes calling :c:func:`dprimme` or :c:func:`zprimme` in parallel.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 1;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int procID

      The identity of the local process within a parallel execution calling :c:func:`dprimme` or
      :c:func:`zprimme`.
      Only the process with id 0 prints information.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | :c:func:`dprimme` sets this field to 0 if |numProcs| is 1;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int nLocal

      Number of local rows on this process.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | :c:func:`dprimme` sets this field to to |n| if |numProcs| is 1;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: void *commInfo

      A pointer to whatever parallel environment structures needed.
      For example, with MPI, it could be a pointer to the MPI communicator.
      PRIMME does not use this. It is available for possible use in 
      user functions defined in |matrixMatvec|,
      |applyPreconditioner|, |massMatrixMatvec| and
      |globalSumDouble|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;

   .. c:member:: void (*globalSumDouble)(double *sendBuf, double *recvBuf, int *count, primme_params *primme)

      Global sum reduction function. No need to set for sequential programs.

      :param sendBuf: array of size ``count`` with the local input values.
      :param recvBuf: array of size ``count`` with the global output values
         so that the i-th element of recvBuf is the sum over all processes of the i-th element
         of ``sendBuf``.
      :param count: array size of ``sendBuf`` and ``recvBuf``.
      :param primme: parameters structure.

      Input/output:

         | :c:func:`primme_initialize` sets this field to an internal function;
         | :c:func:`dprimme` sets this field to an internal function if |numProcs| is 1 and |globalSumDouble| is NULL;
         | this field is read by :c:func:`dprimme`.

      When MPI is used this can be a simply wrapper to MPI_Allreduce().

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

      Input/output:

         | :c:func:`primme_initialize` sets this field to 1;
         | this field is read by :c:func:`primme_set_method` (see :ref:`methods`) and :c:func:`dprimme`.

   .. index:: interior problem

   .. c:member:: primme_target target

      Which eigenpairs to find:

      ``primme_smallest``
         Smallest algebraic eigenvalues; |targetShifts| is ignored.

      ``primme_largest``
         Largest algebraic eigenvalues; |targetShifts| is ignored.

      ``primme_closest_geq``
         Closest to, but greater or equal than the shifts in |targetShifts|.

      ``primme_closest_leq``
         Closest to, but less or equal than the shifts in |targetShifts|.

      ``primme_closest_abs``
         Closest in absolute value to than the shifts in |targetShifts|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to |primme_smallest|;
         | this field is read by :c:func:`dprimme`.

   .. index:: interior problem

   .. c:member:: int numTargetShifts
 
      Size of the array |targetShifts|.
      Used only when |target| is |primme_closest_geq|,
      |primme_closest_leq| or |primme_closest_abs|.
      The default values is 0.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read by :c:func:`dprimme`.

   .. index:: interior problem

   .. c:member:: double *targetShifts

      Array of shifts, at least of size |numTargetShifts|.
      Used only when |target| is |primme_closest_geq|,
      |primme_closest_leq| or |primme_closest_abs|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme`.

      The i-th shift (or the last one, if it is not given) is taken into account in
      finding the i-th eigenvalue.

      .. note::

         Considerations for interior problems:

         * PRIMME will try to compute the eigenvalues in the order given in the
           |targetShifts|. However, for code efficiency and robustness, the shifts
           should be ordered. Order them in ascending (descending) order for shifts
           closer to the lower (higher) end of the spectrum.
         * If some shift is close to the lower (higher) end of the spectrum,
           use either |primme_closest_geq| (|primme_closest_leq|) or
           |primme_closest_abs|.
         * |primme_closest_leq| and |primme_closest_geq| are more efficient
           than |primme_closest_abs|.
         * For interior eigenvalues larger |maxBasisSize| is usually more robust.
 
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

      | $1: Number of converged pairs up to now.
      | $2: The index of the pair currently converged.
      | $3: The eigenvalue.
      | $4: Its residual norm.
      | $5: The current number of matrix-vector products.
      | $6: The current number of outer iterations.
      | $7: The current elapsed time.
      | $8: Index within the block of the targeted pair .
      | $9: QMR norm of the linear system residual.

      In parallel programs, output is produced in call with
      |procID| 0 when |printLevel|
      is from 0 to 4.
      If |printLevel| is 5 output can be produced in any of
      the parallel calls.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 1;
         | this field is read by :c:func:`dprimme`.

   .. note::

      Convergence history for plotting may be produced simply by::

         grep OUT outpufile | awk '{print $8" "$14}' > out
         grep INN outpufile | awk '{print $3" "$11}' > inn

      Then in Matlab::

         plot(out(:,1),out(:,2),'bo');hold; plot(inn(:,1),inn(:,2),'r');

      Or in gnuplot::

         plot 'out' w lp, 'inn' w lp

   .. c:member:: double aNorm

      An estimate of the norm of :math:`A`, which is used in the convergence
      criterion (see |eps|).

      If |aNorm| is less than or equal to 0, the code uses the largest absolute
      Ritz value seen. On return, |aNorm| is then replaced with that value.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0.0;
         | this field is read and written by :c:func:`dprimme`.

   .. c:member:: double eps

      An eigenpairs is marked as converged when the 2-norm of the residual is less
      than |eps| \* |aNorm|.
      The residual vector is :math:`A x - \lambda x` or :math:`A x - \lambda B x`.

      Input/output:

         | :c:func:`primme_initialize` sets this field to :math:`10^{-12}`;
         | this field is read by :c:func:`dprimme`.
 
   .. c:member:: FILE *outputFile

      Opened file to write down the output.

      Input/output:

         | :c:func:`primme_initialize` sets this field to the standard output;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int dynamicMethodSwitch

      If this value is 1, it alternates dynamically between |DEFAULT_MIN_TIME|
      and |DEFAULT_MIN_MATVECS|, trying to identify the fastest method.

      On exit, it holds a recommended method for future runs on this problem:

         | -1: use |DEFAULT_MIN_MATVECS| next time.
         | -2: use |DEFAULT_MIN_TIME| next time.
         | -3: close call, use |DYNAMIC| next time again.
      
      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

      .. note::

         Even for expert users we do not recommend setting |dynamicMethodSwitch|
         directly, but through :c:func:`primme_set_method`.

      .. note::

         The code obtains timings by the ``gettimeofday`` Unix utility. If a cheaper, more
         accurate timer is available, modify the ``PRIMMESRC/COMMONSRC/wtime.c``

   .. c:member:: int locking

      If set to 1, hard locking will be used (locking converged eigenvectors
      out of the search basis). Otherwise the code will try to use soft
      locking (à la ARPACK), when large enough |minRestartSize| is available.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int initSize
 
      On input, the number of initial vector guesses provided in ``evecs`` argument in
      :c:func:`dprimme` or :c:func:`zprimme`.

      On output, |initSize| holds the number of converged eigenpairs. Without |locking| all
      |numEvals| approximations are in ``evecs`` but only the |initSize| ones are
      converged.

      During execution, it holds the current number of converged eigenpairs.
      In addition, if locking is used, these are accessible in ``evals`` and ``evecs``.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`dprimme`.
      
   .. c:member:: int numOrthoConst

      Number of vectors to be used as external orthogonalization constraints.
      These vectors are provided in the first |numOrthoConst| positions of
      the ``evecs`` argument in :c:func:`dprimme` or :c:func:`zprimme` and must be orthonormal.

      PRIMME finds new eigenvectors orthogonal to these constraints (equivalent to solving
      the problem with :math:`(I-YY^*)A(I-YY^*)` and :math:`(I-YY^*)B(I-YY^*)` matrices
      where :math:`Y` are the given constraint vectors).
      This is a handy feature if some eigenvectors are already known, or 
      for finding more eigenvalues after a call to :c:func:`dprimme` or :c:func:`zprimme`,
      possibly with different parameters (see an example in :file:`TEST/ex_zseq.c`).

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int maxBasisSize

      The maximum basis size allowed in the main iteration. This has memory
      implications.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int minRestartSize

      Maximum Ritz vectors kept after restarting the basis.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int maxBlockSize
 
      The maximum block size the code will try to use.

      The user should set
      this based on the architecture specifics of the target computer, 
      as well as any a priori knowledge of multiplicities. The code does 
      *not* require that |maxBlockSize| > 1 to find multiple eigenvalues. For some 
      methods, keeping to 1 yields the best overall performance.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 1;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

      .. note::

         Inner iterations of QMR are not performed in a block fashion.
         Every correction equation from a block is solved independently.

   .. index:: stopping criterion

   .. c:member:: int maxMatvecs

      Maximum number of matrix vector multiplications (approximately equal to 
      the number of preconditioning operations) that the code is allowed to 
      perform before it exits.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``INT_MAX``;
         | this field is read by :c:func:`dprimme`.

   .. index:: stopping criterion

   .. c:member:: int maxOuterIterations

      Maximum number of outer iterations that the code is allowed to perform 
      before it exits.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``INT_MAX``;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int intWorkSize

      If :c:func:`dprimme` or :c:func:`zprimme` is called with all arguments as NULL
      except for :c:type:`primme_params` then PRIMME returns immediately with |intWorkSize|
      containing the size *in bytes* of the integer workspace that will be required by the
      parameters set in PRIMME.

      Otherwise if |intWorkSize| is not 0, it should be the size of the integer work array
      *in bytes* that the user provides in |intWork|. If |intWorkSize| is 0, the code
      will allocate the required space, which can be freed later by calling :c:func:`primme_Free`.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`dprimme`.

   .. c:member:: long int realWorkSize

      If :c:func:`dprimme` or :c:func:`zprimme` is called with all arguments as NULL
      except for :c:type:`primme_params` then PRIMME returns immediately with |realWorkSize|
      containing the size *in bytes* of the real workspace that will be required by the
      parameters set in PRIMME.

      Otherwise if |realWorkSize| is not 0, it should be the size of the real work array
      *in bytes* that the user provides in |realWork|. If |realWorkSize| is 0, the code
      will allocate the required space, which can be freed later by calling :c:func:`primme_Free`.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`dprimme`.

   .. c:member:: int *intWork

      Integer work array.

      If NULL, the code will allocate its own workspace. If the provided space is not
      enough, the code will free it and allocate a new space.

      On exit, the first element shows if a locking problem has occurred.
      Using locking for large |numEvals| may, in some rare cases,
      cause some pairs to be practically converged, in the sense that their components 
      are in the basis of ``evecs``. If this is the case, a Rayleigh Ritz on returned
      ``evecs`` would provide the accurate eigenvectors (see [r4]_).

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read and written by :c:func:`dprimme`.

   .. c:member:: void *realWork

      Real work array.

      If NULL, the code will allocate its own workspace. If the provided space is not
      enough, the code will free it and allocate a new space.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read and written by :c:func:`dprimme`.

   .. c:member:: int iseed

      The ``int iseed[4]`` is an array with the seeds needed by the LAPACK_ dlarnv and zlarnv.

      The default value is an array with values -1, -1, -1 and -1. In that case, ``iseed``
      is set based on the value of |procID| to avoid every parallel process generating
      the same sequence of pseudorandom numbers.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``[-1, -1, -1, -1]``;
         | this field is read and written by :c:func:`dprimme`.

   .. c:member:: void *matrix

      This field may be used to pass any required information 
      in the matrix-vector product |matrixMatvec|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
      
   .. c:member:: void *preconditioner

      This field may be used to pass any required information 
      in the preconditioner function |applyPreconditioner|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;

   .. c:member:: double *ShiftsForPreconditioner

      Array of size ``blockSize`` provided during execution of :c:func:`dprimme` and
      :c:func:`zprimme` holding
      the shifts to be used (if needed) in the preconditioning operation.

      For example if the block size is 3,
      there will be an array of three shifts in |ShiftsForPreconditioner|.
      Then the user can invert a shifted preconditioner for each of the 
      block vectors :math:`(M-ShiftsForPreconditioner_i)^{-1} x_i`.
      Classical Davidson (diagonal) preconditioning is an example of this.
   
      | this field is read and written by :c:func:`dprimme`.

   .. c:member:: primme_restartscheme restartingParams.scheme

      Select a restarting strategy:

      * ``primme_thick``, Thick restarting. This is the most efficient and robust
        in the general case.
      * ``primme_dtr``, Dynamic thick restarting. Helpful without 
        preconditioning but it is expensive to implement.

      Input/output:

         | :c:func:`primme_initialize` sets this field to |primme_thick|;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int restartingParams.maxPrevRetain

      Number of approximations from previous iteration to be retained
      after restart (this is the locally optimal restarting, see [r2]_).
      The restart size is |minRestartSize| plus |maxPrevRetain|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int correctionParams.precondition

      Set to 1 to use preconditioning.
      Make sure |applyPreconditioner| is not NULL then!

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int correctionParams.robustShifts

      Set to 1 to use robust shifting. It tries to avoid stagnation and 
      misconvergence by providing as shifts in |ShiftsForPreconditioner|
      the Ritz values displaced by an approximation of the eigenvalue error.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: int correctionParams.maxInnerIterations

      Control the maximum number of inner QMR iterations:

      * 0:  no inner iterations;
      * >0: perform at most that number of inner iterations per outer step;
      * <0: perform at most the rest of the remaining matrix-vector products
        up to reach |maxMatvecs|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

      See also |convTest|.

   .. c:member:: double correctionParams.relTolBase

      Parameter used when |convTest|
      is |primme_decreasing_LTolerance|.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

   .. c:member:: primme_convergencetest correctionParams.convTest

      Set how to stop the inner QMR method:

      * ``primme_full_LTolerance``: stop by iterations only;

      * ``primme_decreasing_LTolerance``, stop when
        :math:`\text{relTolBase}^{-\text{outIts}}` where outIts
        is the number of outer iterations and retTolBase is set in
        |relTolBase|;
        This is a legacy option from classical JDQR and we recommend
        **strongly** against its use.

      * ``primme_adaptive``, stop when the estimated eigenvalue residual
        has reached the required tolerance (based on Notay's JDCG).

      * ``primme_adaptive_ETolerance``, as |primme_adaptive| but also
        stopping when the estimated eigenvalue residual has reduced 10
        times.

      Input/output:

         | :c:func:`primme_initialize` sets this field to ``primme_adaptive_ETolerance``;
         | written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

      .. note::

         Avoid to set |maxInnerIterations| to -1 and |convTest| to |primme_full_LTolerance|.

      See also |maxInnerIterations|.

   .. c:member:: int correctionParams.projectors.LeftQ
   .. c:member:: int correctionParams.projectors.LeftX
   .. c:member:: int correctionParams.projectors.RightQ
   .. c:member:: int correctionParams.projectors.RightX
   .. c:member:: int correctionParams.projectors.SkewQ
   .. c:member:: int correctionParams.projectors.SkewX

      Control the projectors involved in the computation of the correction
      appended to the basis every (outer) iteration.

      Consider the current selected Ritz value :math:`\Lambda` and vectors :math:`X`,
      the residual associated vectors :math:`R=AX-X\Lambda`, the previous locked vectors
      :math:`Q`, and the preconditioner :math:`M^{-1}`.

      When |maxInnerIterations| is 0, the correction :math:`D` appended to the basis
      in GD is:

      .. only:: not text

         +--------+-------+------------------------------------------------------------------+
         | RightX | SkewX |        :math:`D`                                                 |
         +========+=======+==================================================================+
         |    0   |   0   | :math:`M^{-1}R` (Classic GD)                                     |
         +--------+-------+------------------------------------------------------------------+
         |    1   |   0   | :math:`M^{-1}(R-\Delta X)` (cheap Olsen's Method)                |
         +--------+-------+------------------------------------------------------------------+
         |    1   |   1   | :math:`(I- M^{-1}X(X^*M^{-1}X)^{-1}X^*)M^{-1}R` (Olsen's Method) |
         +--------+-------+------------------------------------------------------------------+
         |    0   |   1   | error                                                            |
         +--------+-------+------------------------------------------------------------------+

      .. only:: text

         +--------+-------+----------------------------------------------------------+
         | RightX | SkewX | D                                                        |
         +========+=======+==========================================================+
         |    0   |   0   | M^{-1}R (Classic GD)                                     |
         +--------+-------+----------------------------------------------------------+
         |    1   |   0   | M^{-1}(R-\Delta X) (cheap Olsen's Method)                |
         +--------+-------+----------------------------------------------------------+
         |    1   |   1   | (I- M^{-1}X(X^*M^{-1}X)^{-1}X^*)M^{-1}R (Olsen's Method) |
         +--------+-------+----------------------------------------------------------+
         |    0   |   1   | error                                                    |
         +--------+-------+----------------------------------------------------------+


      Where :math:`\Delta` is a diagonal matrix that :math:`\Delta_{i,i}` holds an estimation
      of the error of the approximate eigenvalue :math:`\Lambda_{i,i}`.
 
      The values of ``RightQ``, ``SkewQ``, ``LeftX`` and ``LeftQ`` are ignored.

      When |maxInnerIterations| is not 0, the correction :math:`D` in Jacobi-Davidson results
      from solving:

      .. math::

         P_Q^l P_X^l (A-\sigma I) P_X^r P_Q^r M^{-1} D' = -R, \ \ \  D = P_X^r P_Q^l M^{-1}D'.

      For ``LeftQ``:

         | 0: :math:`P_Q^l = I`;
         | 1: :math:`P_Q^l = I - QQ^*`.

      For ``LeftX``:

         | 0: :math:`P_X^l = I`;
         | 1: :math:`P_X^l = I - XX^*`.

      For ``RightQ`` and ``SkewQ``:

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

      For ``RightX`` and ``SkewX``:

      +--------+-------+-------------------------------+
      | RightX | SkewX |        :math:`P_X^r`          |
      +========+=======+===============================+
      |    0   |   0   | :math:`I`                     |
      +--------+-------+-------------------------------+
      |    1   |   0   | :math:`I - XX^*`              |
      +--------+-------+-------------------------------+
      |    1   |   1   | :math:`I - KX(X^*KX)^{-1}X^*` |
      +--------+-------+-------------------------------+
      |    0   |   1   | error                         |
      +--------+-------+-------------------------------+

      Input/output:

         | :c:func:`primme_initialize` sets all of them to 0;
         | this field is written by :c:func:`primme_set_method` (see :ref:`methods`);
         | this field is read by :c:func:`dprimme`.

      See [r3]_ for a study about different projector configurations in JD.

   .. c:member:: int stats.numOuterIterations

      Hold the number of outer iterations. The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: int stats.numRestarts

      Hold the number of restarts during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: int stats.numMatvecs

      Hold how many vectors the operator in |matrixMatvec| has been applied on.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: int stats.numPreconds

      Hold how many vectors the operator in |applyPreconditioner| has been applied on.
      The value is available during execution and at the end.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

   .. c:member:: int stats.elapsedTime

      Hold the wall clock time spent by the call to :c:func:`dprimme` or :c:func:`zprimme`.
      The value is available at the end of the execution.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | written by :c:func:`dprimme`.

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


.. _error-codes:

Error Codes
"""""""""""

The functions :c:func:`dprimme` and :c:func:`zprimme` return one of the next values:

*  0: success.
*  1: reported only amount of required memory.
* -1: failed in allocating int or real workspace.
* -2: malloc failed in allocating a permutation integer array.
* -3: main_iter() encountered problem; the calling stack of the
  functions where the error occurred was printed in ``stderr``.
* -4: if argument ``primme`` is NULL.
* -5: if |n| <= 0 or |nLocal| <= 0.
* -6: if |numProcs| < 1.
* -7: if |matrixMatvec| is NULL.
* -8: if |applyPreconditioner| is NULL and 
  |precondition| is not NULL.
* -9: if |globalSumDouble| is NULL.
* -10: if |numEvals| > |n|.
* -11: if |numEvals| < 0.
* -12: if |eps| > 0 and |eps| < machine precision.
* -13: if |target| is not properly defined.
* -14: if |target| is one of |primme_closest_geq|,
  |primme_closest_leq| or |primme_closest_abs| but
  |numTargetShifts| <= 0 (no shifts).
* -15: if |target| is one of |primme_closest_geq|,
  |primme_closest_leq| or |primme_closest_abs| but
  |targetShifts| is NULL  (no shifts array).
* -16: if |numOrthoConst| < 0 or
  |numOrthoConst| >= |n|.
  (no free dimensions left).
* -17: if |maxBasisSize| < 2.
* -18: if |minRestartSize| <= 0.
* -19: if |maxBlockSize| <= 0.
* -20: if |maxPrevRetain| < 0.
* -21: if |scheme| is not one of `primme_thick` or `primme_dtr`.
* -22: if |initSize| < 0.
* -23: if not |locking| and |initSize| > |maxBasisSize|.
* -24: if |locking| and |initSize| > |numEvals|.
* -25: if |maxPrevRetain| + |minRestartSize| >= |maxBasisSize|.
* -26: if |minRestartSize| >= |n|.
* -27: if |printLevel| < 0 or |printLevel| > 5.
* -28: if |convTest| is not one of
  |primme_full_LTolerance|, |primme_decreasing_LTolerance|,
  |primme_adaptive_ETolerance| or |primme_adaptive|.
* -29: if |convTest| == |primme_decreasing_LTolerance| and |relTolBase| <= 1.
* -30: if ``evals`` is NULL, but not ``evecs`` and ``resNorms``.
* -31: if ``evecs`` is NULL, but not ``evals`` and ``resNorms``.
* -32: if ``resNorms`` is NULL, but not ``evecs`` and ``evals``.

.. _methods:

Preset Methods
""""""""""""""

.. c:type:: primme_preset_method

   .. c:member:: DEFAULT_MIN_TIME

      Set as |JDQMR_ETol|; this method is usually the fastest if
      the cost of the matrix vector product is inexpensive.

   .. c:member:: DEFAULT_MIN_MATVECS

      Currently set as |GD_Olsen_plusK|; this method usually performs
      fewer matrix vector products than other methods, so it's a good
      choice when this operation is expensive.

   .. c:member:: DYNAMIC

      Switches to the best method dynamically; currently, between
      methods |DEFAULT_MIN_TIME| and |DEFAULT_MIN_MATVECS|.

      With |DYNAMIC| :c:func:`primme_set_method` sets
      |dynamicMethodSwitch| = 1 and makes the same changes as
      for method |JDQMR_ETol|.

   .. c:member:: Arnoldi

      Arnoldi implemented à la Generalized Davidson.

      With |Arnoldi| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking| = 0;
         * |maxPrevRetain| = 0;
         * |precondition| = 0;
         * |maxInnerIterations| = 0.

   .. c:member:: GD

      Generalized Davidson.

      With |GD| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking| = 0;
         * |maxPrevRetain| = 0;
         * |robustShifts| = 1;
         * |maxInnerIterations| = 0;
         * |RightX| = 0;
         * |SkewX| = 0.

   .. c:member:: GD_plusK

      GD with locally optimal restarting. 

      With |GD_plusK| :c:func:`primme_set_method` sets
      |maxPrevRetain| = 2 if |maxBlockSize| is 1 and |numEvals| > 1;
      otherwise it sets |maxPrevRetain| to |maxBlockSize|. Also:

      .. hlist::

         * |locking| = 0;
         * |maxInnerIterations| = 0;
         * |RightX| = 0;
         * |SkewX| = 0.

   .. c:member:: GD_Olsen_plusK

      GD+k and the cheap Olsen's Method.

      With |GD_Olsen_plusK| :c:func:`primme_set_method` makes the
      same changes as for method |GD_plusK| and sets |RightX| = 1.

   .. c:member:: JD_Olsen_plusK

      GD+k and Olsen's Method.

      With |JD_Olsen_plusK| :c:func:`primme_set_method` makes the
      same changes as for method |GD_plusK| and also sets
      |robustShifts| = 1, |RightX| to 1, and |SkewX| to 1.

   .. c:member:: RQI

      (Accelerated) Rayleigh Quotient Iteration.

      With |RQI| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking| = 1;
         * |maxPrevRetain| = 0;
         * |robustShifts|  = 1;
         * |maxInnerIterations| = -1;
         * |LeftQ|   = 1;
         * |LeftX|   = 1;
         * |RightQ|  = 0;
         * |RightX|  = 1;
         * |SkewQ|   = 0;
         * |SkewX|   = 0;
         * |convTest| = |primme_full_LTolerance|.

      .. note::

         If |numTargetShifts| > 0 and |targetShifts| are provided, the interior problem
         solved uses these shifts in the correction equation. Therefore RQI becomes
         INVIT (inverse iteration) in that case.

   .. c:member:: JDQR

      Jacobi-Davidson with fixed number of inner steps.

      With |JDQR| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking|     = 1;
         * |maxPrevRetain|      = 1;
         * |robustShifts|       = 0;
         * |maxInnerIterations| = 10 if it is 0;
         * |LeftQ|   = 0;
         * |LeftX|   = 1;
         * |RightQ|  = 1;
         * |RightX|  = 1;
         * |SkewQ|   = 1;
         * |SkewX|   = 1;
         * |relTolBase| = 1.5;
         * |convTest| = |primme_full_LTolerance|.

   .. c:member:: JDQMR

      Jacobi-Davidson with adaptive stopping criterion for inner Quasi Minimum Residual (QMR).

      With |JDQMR| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking| = 0;
         * |maxPrevRetain| = 1 if it is 0
         * |maxInnerIterations| = -1;
         * |LeftQ|   = |precondition|;
         * |LeftX|   = 1;
         * |RightQ|  = 0;
         * |RightX|  = 0;
         * |SkewQ|   = 0;
         * |SkewX|   = 1;
         * |convTest|  = |primme_adaptive|.

   .. c:member:: JDQMR_ETol

      JDQMR but QMR stops after residual norm reduces by a 0.1 factor.

      With |JDQMR_ETol| :c:func:`primme_set_method` makes the same
      changes as for the method |JDQMR| and sets
      |convTest| = |primme_adaptive_ETolerance|.

   .. c:member:: SUBSPACE_ITERATION

      Subspace iteration.

      With |SUBSPACE_ITERATION| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking|    = 1;
         * |maxBasisSize| = |numEvals| `*` 2;
         * |minRestartSize| = |numEvals|;
         * |maxBlockSize| = |numEvals|;
         * |scheme|  = |primme_thick|;
         * |maxPrevRetain|      = 0;
         * |robustShifts|       = 0;
         * |maxInnerIterations| = 0;
         * |RightX|  = 1;
         * |SkewX|   = 0.

   .. c:member:: LOBPCG_OrthoBasis

      LOBPCG with orthogonal basis.

      With |LOBPCG_OrthoBasis| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking|    = 0;
         * |maxBasisSize| = |numEvals| `*` 3;
         * |minRestartSize| = |numEvals|;
         * |maxBlockSize| = |numEvals|;
         * |scheme|  = |primme_thick|;
         * |maxPrevRetain|      = |numEvals|;
         * |robustShifts|       = 0;
         * |maxInnerIterations| = 0;
         * |RightX|  = 1;
         * |SkewX|   = 0.

   .. c:member:: LOBPCG_OrthoBasis_Window

      LOBPCG with sliding window of |maxBlockSize| < 3 `*` |numEvals|.

      With |LOBPCG_OrthoBasis_Window| :c:func:`primme_set_method` sets:

      .. hlist::

         * |locking|    = 0;
         * |maxBasisSize| = |maxBlockSize| `*` 3;
         * |minRestartSize| = |maxBlockSize|;
         * |maxBlockSize| = |numEvals|;
         * |scheme|  = |primme_thick|;
         * |maxPrevRetain|      = |maxBlockSize|;
         * |robustShifts|       = 0;
         * |maxInnerIterations| = 0;
         * |RightX|  = 1;
         * |SkewX|   = 0.

.. include:: epilog.inc
