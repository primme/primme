.. highlight:: c

Parameter Description
---------------------

primme_svds_params
""""""""""""""""""

.. c:type:: primme_svds_params

   Structure to set the problem matrix and the solver options.

   .. c:member:: PRIMME_INT m

      Number of rows of the matrix.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: PRIMME_INT n

      Number of columns of the matrix.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 0;
         | this field is read by :c:func:`dprimme`.

   .. c:member:: void (*matrixMatvec) (void *x, PRIMME_INT ldx, void *y, PRIMME_INT ldy, int *blockSize, int *transpose, primme_svds_params *primme_svds, int *ierr)

      Block matrix-multivector multiplication, :math:`y = A x` if ``transpose`` is zero, and :math:`y = A^*x` otherwise.

      :param x: input array.
      :param ldx: leading dimension of ``x``.
      :param y: output array.
      :param ldy: leading dimension of ``y``.
      :param blockSize: number of columns in ``x`` and ``y``.
      :param transpose: if non-zero, the transpose A should be applied.
      :param primme_svds: parameters structure.
      :param ierr: output error code; if it is set to non-zero, the current call to PRIMME will stop.

      If ``transpose`` is zero, then ``x`` and ``y`` are arrays of dimensions |SnLocal| x ``blockSize`` and |SmLocal| x ``blockSize``
      respectively. Elsewhere they have dimensions |SmLocal| x ``blockSize`` and |SnLocal| x ``blockSize``.
      Both arrays are in column-major_ order (elements in the same column with consecutive row indices are consecutive in memory).

      The actual type of ``x`` and ``y`` depends on which function is being calling. For :c:func:`dprimme_svds`, it is ``double``,
      for :c:func:`zprimme_svds` it is :c:type:`PRIMME_COMPLEX_DOUBLE`, for :c:func:`sprimme_svds` it is ``float`` and
      for :c:func:`cprimme_svds` it is :c:type:`PRIMME_COMPLEX_FLOAT`.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

      .. note::

         Integer arguments are passed by reference to make easier the interface to other
         languages (like Fortran).

   .. c:member:: void (*applyPreconditioner)(void *x, PRIMME_INT ldx, void *y, PRIMME_INT ldy, int *blockSize, int *mode, primme_svds_params *primme_svds, int *ierr)

      Block preconditioner-multivector application, :math:`y = M^{-1}x` for finding singular values close to :math:`\sigma`.
      Depending on ``mode``, :math:`M` is expected to be an approximation of the following operators:

      * ``primme_svds_op_AtA``: :math:`M \approx A^*Ax - \sigma^2 I`,
      * ``primme_svds_op_AAt``: :math:`M \approx AA^*x - \sigma^2 I`,
      * ``primme_svds_op_augmented``: :math:`M \approx \left(\begin{array}{cc} 0 & A^* \\ A & 0 \end{array}\right) - \sigma I`.

      :param x: input array.
      :param ldx: leading dimension of ``x``.
      :param y: output array.
      :param ldy: leading dimension of ``y``.
      :param blockSize: number of columns in ``x`` and ``y``.
      :param mode: one of ``primme_svds_op_AtA``, ``primme_svds_op_AAt`` or ``primme_svds_op_augmented``.
      :param primme_svds: parameters structure.
      :param ierr: output error code; if it is set to non-zero, the current call to PRIMME will stop.

      If ``mode`` is ``primme_svds_op_AtA``, then ``x`` and ``y`` are arrays of dimensions |SnLocal| x ``blockSize``; if mode is
      ``primme_svds_op_AAt``, they are |SmLocal| x ``blockSize``; and otherwise they are (|SmLocal| + |SnLocal|) x ``blockSize``.
      Both arrays are in column-major_ order (elements in the same column with consecutive row indices are consecutive in memory).

      The actual type of ``x`` and ``y`` depends on which function is being calling. For :c:func:`dprimme_svds`, it is ``double``,
      for :c:func:`zprimme_svds` it is :c:type:`PRIMME_COMPLEX_DOUBLE`, for :c:func:`sprimme_svds` it is ``float`` and
      for :c:func:`cprimme_svds` it is :c:type:`PRIMME_COMPLEX_FLOAT`.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: int numProcs

      Number of processes calling :c:func:`dprimme_svds` or :c:func:`zprimme_svds` in parallel.

      Input/output:

         | :c:func:`primme_initialize` sets this field to 1;
         | this field is read by :c:func:`dprimme` and :c:func:`zprimme_svds`.

   .. c:member:: int procID

      The identity of the local process within a parallel execution calling :c:func:`dprimme_svds` or
      :c:func:`zprimme_svds`.
      Only the process with id 0 prints information.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | :c:func:`dprimme_svds` sets this field to 0 if |SnumProcs| is 1;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: PRIMME_INT mLocal

      Number of local rows on this process. The value depends on how the matrix and
      preconditioner is distributed along the processes.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | :c:func:`dprimme_svds` sets this field to |Sm| if |SnumProcs| is 1;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

      See also: |SmatrixMatvec| and |SapplyPreconditioner|.

   .. c:member:: PRIMME_INT nLocal

      Number of local columns on this process. The value depends on how the matrix and
      preconditioner is distributed along the processes.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | :c:func:`dprimme_svds` sets this field to to |n| if |SnumProcs| is 1;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: void *commInfo

      A pointer to whatever parallel environment structures needed.
      For example, with MPI, it could be a pointer to the MPI communicator.
      PRIMME does not use this. It is available for possible use in 
      user functions defined in |SmatrixMatvec|, |SapplyPreconditioner| and
      |SglobalSumReal|.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to NULL;

   .. c:member:: void (*globalSumReal)(double *sendBuf, double *recvBuf, int *count, primme_svds_params *primme_svds, int *ierr)

      Global sum reduction function. No need to set for sequential programs.

      :param sendBuf: array of size ``count`` with the local input values.
      :param recvBuf: array of size ``count`` with the global output values
         so that the i-th element of recvBuf is the sum over all processes of the i-th element
         of ``sendBuf``.
      :param count: array size of ``sendBuf`` and ``recvBuf``.
      :param primme_svds: parameters structure.
      :param ierr: output error code; if it is set to non-zero, the current call to PRIMME will stop.

      The actual type of ``sendBuf`` and ``recvBuf`` depends on which function is being calling. For :c:func:`dprimme_svds`
      and :c:func:`zprimme_svds` it is ``double``, and for :c:func:`sprimme_svds` and  :c:func:`cprimme_svds` it is ``float``.
      Note that ``count`` is the number of values of the actual type.
 
      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to an internal function;
         | :c:func:`dprimme_svds` sets this field to an internal function if |SnumProcs| is 1 and |SglobalSumReal| is NULL;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

      When MPI is used, this can be a simply wrapper to MPI_Allreduce() as shown below::

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

   .. c:member:: int numSvals

      Number of singular triplets wanted.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 1;
         | this field is read by :c:func:`primme_svds_set_method` (see :ref:`methods_svds`) and :c:func:`dprimme_svds`.


   .. index:: interior problem

   .. c:member:: primme_svds_target target

      Which singular values to find:

      ``primme_svds_smallest``
         Smallest singular values; |StargetShifts| is ignored.

      ``primme_svds_largest``
         Largest singular values; |StargetShifts| is ignored.

      ``primme_svds_closest_abs``
         Closest in absolute value to the shifts in |StargetShifts|.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to |primme_svds_smallest|;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. index:: interior problem

   .. c:member:: int numTargetShifts
 
      Size of the array |StargetShifts|.
      Used only when |Starget| is |primme_svds_closest_abs|.
      The default values is 0.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. index:: interior problem

   .. c:member:: double *targetShifts

      Array of shifts, at least of size |SnumTargetShifts|.
      Used only when |Starget| is |primme_svds_closest_abs|.

      Singular values are computed in order so that the
      i-th singular value is the closest to the i-th shift. If |SnumTargetShifts| < |SnumSvals|, the last shift given
      is used for all the remaining i's.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to NULL;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

      .. note::

         Eventually this is used by  :c:func:`dprimme` and :c:func:`zprimme`. Please
         see considerations of |targetShifts|.

   .. c:member:: int printLevel

      The level of message reporting from the code. All output is written in |SoutputFile|.

      One of:
 
      * 0: silent.
      * 1: print some error messages when these occur.
      * 2: as in 1, and info about targeted singular triplets when they are marked as converged::
      
            #Converged $1 sval[ $2 ]= $3 norm $4 Mvecs $5 Time $7 stage $10

        or locked::

            #Lock striplet[ $1 ]= $3 norm $4 Mvecs $5 Time $7 stage $10

      * 3: as in 2, and info about targeted singular triplets every outer iteration::
      
            OUT $6 conv $1 blk $8 MV $5 Sec $7 SV $3 |r| $4 stage $10

        Also, if using |DYNAMIC|, show JDQMR/GD+k performance ratio and
        the current method in use.
      * 4: as in 3, and info about targeted singular triplets every inner iteration::
      
            INN MV $5 Sec $7 Sval $3 Lin|r| $9 SV|r| $4 stage $10
      
      * 5: as in 4, and verbose info about certain choices of the algorithm.
      
      Output key:

      | $1: Number of converged triplets up to now.
      | $2: The index of the triplet currently converged.
      | $3: The singular value.
      | $4: Its residual norm.
      | $5: The current number of matrix-vector products.
      | $6: The current number of outer iterations.
      | $7: The current elapsed time.
      | $8: Index within the block of the targeted triplet.
      | $9: QMR norm of the linear system residual.
      | $10: stage (1 or 2)

      In parallel programs, when |SprintLevel| is 0 to 4 only |SprocID| 0 produces output.
      For |SprintLevel| 5 output can be produced in any of the parallel calls.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 1;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

    .. note::

      Convergence history for plotting may be produced simply by:

      .. code-block:: bash

         grep OUT outpufile | awk '{print $8" "$14}' > out
         grep INN outpufile | awk '{print $3" "$11}' > inn

      Or in gnuplot:

      .. code-block:: gnuplot

         plot 'out' w lp, 'inn' w lp


   .. c:member:: double aNorm

      An estimate of the 2-norm of :math:`A`, which is used in the default convergence
      criterion (see |Seps|).

      If |aNorm| is less than or equal to 0, the code uses the largest absolute
      Ritz value seen. On return, |SaNorm| is then replaced with that value.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0.0;
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: double eps

      A triplet :math:`(u,\sigma,v)` is marked as converged when
      :math:`\sqrt{\|A v - \sigma u\|^2 + \|A^* u - \sigma v\|^2}`
      is less than |Seps| \* |SaNorm|, or close to the minimum tolerance that
      the selected method can achieve in the given machine precision. See :ref:`methods_svds`.

      The default value is machine precision times :math:`10^4`.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0.0;
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.
 
   .. c:member:: FILE *outputFile

      Opened file to write down the output.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to the standard output;
         | this field is read by :c:func:`dprimme_svds`, :c:func:`zprimme_svds` and :c:func:`primme_svds_display_params`

   .. c:member:: int locking

      If set to 1, the underneath eigensolvers will use hard locking. See |locking|.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to -1;
         | written by :c:func:`primme_svds_set_method` (see :ref:`methods_svds`);
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: int initSize
 
      On input, the number of initial vector guesses provided in ``svecs`` argument in
      :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

      On output, |SinitSize| holds the number of converged triplets. Without |Slocking| all
      |SnumSvals| approximations are in ``svecs`` but only the first |SinitSize| are
      converged.

      During execution, it holds the current number of converged triplets.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

    .. c:member:: int numOrthoConst

      Number of vectors to be used as external orthogonalization constraints.
      The left and the right vector constraints are provided as input of
      the ``svecs`` argument in :c:func:`sprimme_svds` or other variant, and must be orthonormal.

      PRIMME SVDS finds new triplets orthogonal to these constraints (equivalent to solving
      the problem :math:`(I-UU^*)A(I-VV^*)` where :math:`U` and :math:`V` are the given left and right
      constraint vectors).
      This is a handy feature if some singular triplets are already known, or 
      for finding more triplets after a call to :c:func:`dprimme_svds` or :c:func:`zprimme_svds`,
      possibly with different parameters (see an example in :file:`TEST/exsvd_zseq.c`).

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: int maxBasisSize

      The maximum basis size allowed in the main iteration. This has memory
      implications.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_svds_set_method` (see :ref:`methods_svds`);
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: int maxBlockSize
 
      The maximum block size the code will try to use.

      The user should set
      this based on the architecture specifics of the target computer, 
      as well as any a priori knowledge of multiplicities. The code does 
      *not* require that |maxBlockSize| > 1 to find multiple triplets. For some 
      methods, keeping to 1 yields the best overall performance.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 1;
         | this field is read and written by :c:func:`primme_svds_set_method` (see :ref:`methods_svds`);
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. index:: stopping criterion

   .. c:member:: PRIMME_INT maxMatvecs

      Maximum number of matrix vector multiplications (approximately half 
      the number of preconditioning operations) that the code is allowed to 
      perform before it exits.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to ``INT_MAX``;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: int intWorkSize

      If :c:func:`dprimme_svds` or :c:func:`zprimme_svds` is called with all arguments as NULL
      except for :c:type:`primme_svds_params` then it returns immediately with |SintWorkSize|
      containing the size *in bytes* of the integer workspace that will be required by the
      parameters set.

      Otherwise if |SintWorkSize| is not 0, it should be the size of the integer work array
      *in bytes* that the user provides in |SintWork|. If |SintWorkSize| is 0, the code
      will allocate the required space, which can be freed later by calling :c:func:`primme_svds_free`.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: size_t realWorkSize

      If :c:func:`dprimme_svds` or :c:func:`zprimme_svds` is called with all arguments as NULL
      except for :c:type:`primme_svds_params` then it returns immediately with |SrealWorkSize|
      containing the size *in bytes* of the real workspace that will be required by the
      parameters set.

      Otherwise if |SrealWorkSize| is not 0, it should be the size of the real work array
      *in bytes* that the user provides in |SrealWork|. If |SrealWorkSize| is 0, the code
      will allocate the required space, which can be freed later by calling :c:func:`primme_svds_free`.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: int *intWork

      Integer work array.

      If NULL, the code will allocate its own workspace. If the provided space is not
      enough, the code will return the error code ``-21``.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to NULL;
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: void *realWork

      Real work array.

      If NULL, the code will allocate its own workspace. If the provided space is not
      enough, the code will return the error code ``-20``.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to NULL;
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: PRIMME_INT iseed

      The ``PRIMME_INT iseed[4]`` is an array with the seeds needed by the LAPACK_ dlarnv and zlarnv.

      The default value is an array with values -1, -1, -1 and -1. In that case, ``iseed``
      is set based on the value of |SprocID| to avoid every parallel process generating
      the same sequence of pseudorandom numbers.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to ``[-1, -1, -1, -1]``;
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: void *matrix

      This field may be used to pass any required information 
      in the matrix-vector product |SmatrixMatvec|.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to NULL;
  
   .. c:member:: void *preconditioner

      This field may be used to pass any required information 
      in the preconditioner function |SapplyPreconditioner|.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to NULL;

   .. c:member:: int precondition

      Set to 1 to use preconditioning.
      Make sure |SapplyPreconditioner| is not NULL then!

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | this field is read and written by :c:func:`primme_svds_set_method` (see :ref:`methods_svds`);
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: primme_svds_op_operator method

      Select the equivalent eigenvalue problem that will be solved:

      * ``primme_svds_op_AtA``: :math:`A^*Ax = \sigma^2 x`,
      * ``primme_svds_op_AAt``: :math:`AA^*x = \sigma^2 x`,
      * ``primme_svds_op_augmented``: :math:`\left(\begin{array}{cc} 0 & A^* \\ A & 0 \end{array}\right) x = \sigma x`.

      The options for this solver are stored in |Sprimme|.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to ``primme_svds_op_none``;
         | this field is read and written by :c:func:`primme_svds_set_method` (see :ref:`methods_svds`);
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: primme_svds_op_operator methodStage2

      Select the equivalent eigenvalue problem that will be solved to refine the solution. The allowed options
      are ``primme_svds_op_none`` to not refine the solution and ``primme_svds_op_augmented`` to refine the
      solution by solving the augmented problem with the current solution as the initial vectors. See |Smethod|.

      The options for this solver are stored in |SprimmeStage2|.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to ``primme_svds_op_none``;
         | this field is read and written by :c:func:`primme_svds_set_method` (see :ref:`methods_svds`);
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: primme_params primme

      Parameter structure storing the options for underneath eigensolver that will be called at the first stage.
      See |Smethod|.

      Input/output:

         | :c:func:`primme_svds_initialize` initialize this structure;
         | this field is read and written by :c:func:`primme_svds_set_method` (see :ref:`methods_svds`);
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: primme_params primmeStage2

      Parameter structure storing the options for underneath eigensolver that will be called at the second stage.
      See |SmethodStage2|.

      Input/output:

         | :c:func:`primme_svds_initialize` initialize this structure;
         | this field is read and written by :c:func:`primme_svds_set_method` (see :ref:`methods_svds`);
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: void (*monitorFun)(void *basisSvals, int *basisSize, int *basisFlags, int *iblock, int *blockSize, void *basisNorms, int *numConverged, void *lockedSvals, int *numLocked, int *lockedFlags, void *lockedNorms, int *inner_its, void *LSRes, primme_event *event, int *stage, struct primme_svds_params *primme_svds, int *ierr)


      Convergence monitor. Used to customize how to report solver 
      information during execution (stage, iteration number, matvecs, time, 
      residual norms, targets, etc).

      :param basisSvals:   array with approximate singular values of the basis.
      :param basisSize:    size of the arrays ``basisSvals``, ``basisFlags`` and ``basisNorms``.
      :param basisFlags:   state of every approximate triplet in the basis.
      :param iblock:       indices of the approximate triplet in the block.
      :param blockSize:    size of array ``iblock``.
      :param basisNorms:   array with residual norms of the triplets in the basis.
      :param numConverged: number of triplets converged in the basis plus the number of the locked triplets (note that this value isn't monotonic).
      :param lockedSvals:  array with the locked triplets.
      :param numLocked:    size of the arrays ``lockedSvals``, ``lockedFlags`` and ``lockedNorms``.
      :param lockedFlags:  state of each locked triplets.
      :param lockedNorms:  array with residual norms of the locked triplets.
      :param inner_its:    number of performed QMR iterations in the current correction equation.
      :param LSRes:        residual norm of the linear system at the current QMR iteration.
      :param event:        event reported.
      :param stage:        ``0`` for first stage, ``1`` for second stage.
      :param primme_svds:  parameters structure; the counter in ``stats`` are updated with the current number of matrix-vector products, iterations, elapsed time, etc., since start.
      :param ierr:         output error code; if it is set to non-zero, the current call to PRIMME will stop.

      This function is called at the next events:

      * ``*event == primme_event_outer_iteration``: every outer iterations.

        It is provided ``basisSvals``, ``basisSize``, ``basisFlags``, ``iblock`` and ``blockSize``.

        ``basisNorms[iblock[i]]`` has the residual norms for the selected triplets in the block.
        PRIMME avoids computing the residual of soft-locked triplets, ``basisNorms[i]`` for ``i<iblock[0]``.
        So those values may correspond to previous iterations. The values ``basisNorms[i]`` for ``i>iblock[blockSize-1]``
        are not valid.

        If |Slocking| is enabled, ``lockedSvals``, ``numLocked``, ``lockedFlags`` and ``lockedNorms`` are also provided.

        ``inner_its`` and  ``LSRes`` are not provided.

      * ``*event == primme_event_inner_iteration``: every QMR iteration.

        ``basisSvals[0]`` and ``basisNorms[0]`` provides the approximate singular value and the residual norm
        of the triplet which is improved in the current correction equation. If |convTest| is |primme_adaptive|
        or |primme_adaptive_ETolerance|, ``basisSvals[0]``, and ``basisNorms[0]`` are updated every QMR iteration.

        ``inner_its`` and  ``LSRes`` are also provided.

        ``lockedSvals``, ``numLocked``, ``lockedFlags``, and ``lockedNorms`` may not be provided.

      * ``*event == primme_event_convergence``: a new triplet in the basis passed the convergence criterion

        ``iblock[0]`` is the index of the newly converged triplet in the basis which will be locked or soft locked.
        The following are provided: ``basisSvals``, ``basisSize``, ``basisFlags`` and ``blockSize[0]==1``.

        ``lockedSvals``, ``numLocked``, ``lockedFlags`` and ``lockedNorms`` may not be provided.

        ``inner_its`` and  ``LSRes`` are not provided.

      * ``*event == primme_event_locked``: a new triplet added to the locked singular vectors.

        ``lockedSvals``, ``numLocked``, ``lockedFlags`` and ``lockedNorms`` are provided.
        The last element of ``lockedSvals``, ``lockedFlags`` and ``lockedNorms`` corresponds
        to the recent locked triplet.
 
        ``basisSvals``, ``numConverged``, ``basisFlags`` and ``basisNorms`` may not be provided.

        ``inner_its`` and  ``LSRes`` are not be provided.

      The values of ``basisFlags`` and ``lockedFlags`` are:

      * ``0``: unconverged.
      * ``1``: internal use; only in ``basisFlags``.
      * ``2``: passed convergence test (see |Seps|).
      * ``3``: converged because the solver may not be able to reduce the residual norm further.

      Input/output:

         | :c:func:`primme_initialize` sets this field to NULL;
         | :c:func:`dprimme_svds` sets this field to an internal function if it is NULL;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.


   .. c:member:: PRIMME_INT stats.numOuterIterations

      Hold the number of outer iterations.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: PRIMME_INT stats.numRestarts

      Hold the number of restarts.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: PRIMME_INT stats.numMatvecs

      Hold how many vectors the operator in |SmatrixMatvec| has been applied on.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: PRIMME_INT stats.numPreconds

      Hold how many vectors the operator in |SapplyPreconditioner| has been applied on.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: double stats.elapsedTime

      Hold the wall clock time spent by the call to :c:func:`dprimme_svds` or :c:func:`zprimme_svds`.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0;
         | written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

.. _methods_svds:

Preset Methods
--------------

.. c:type:: primme_svds_preset_method

   .. c:member:: primme_svds_default

      Set as :c:member:`primme_svds_hybrid`.

   .. c:member:: primme_svds_normalequations

      Solve the equivalent eigenvalue problem :math:`A^*A V = \Sigma^2 V` and computes :math:`U` by normalizing
      the vectors :math:`AV`. If |Sm| is smaller than |Sn|, :math:`AA^*` is solved instead.
  
      With :c:member:`primme_svds_normalequations` :c:func:`primme_svds_set_method` sets
      |Smethod| to ``primme_svds_op_AtA`` if |Sm| is larger or equal than |Sn|, and to ``primme_svds_op_AAt``
      otherwise; and |SmethodStage2| is set to ``primme_svds_op_none``.

      The minimum residual norm that this method can achieve is :math:`\|A\|\epsilon\sigma^{-1}`,
      where :math:`\epsilon` is the machine precision and :math:`\sigma` the required singular value.
 
   .. c:member:: primme_svds_augmented

      Solve the equivalent eigenvalue problem :math:`\left(\begin{array}{cc} 0 & A^* \\ A & 0 \end{array}\right) X = \sigma X`
      with :math:`X = \left(\begin{array}{cc}V\\U\end{array}\right)`.
  
      With :c:member:`primme_svds_augmented` :c:func:`primme_svds_set_method` sets
      |Smethod| to ``primme_svds_op_augmented`` and |SmethodStage2| to ``primme_svds_op_none``.
 
      The minimum residual norm that this method can achieve is :math:`\|A\|\epsilon`,
      where :math:`\epsilon` is the machine precision.
      However it may not return triplets with singular values smaller than :math:`\|A\|\epsilon`.
 
   .. c:member:: primme_svds_hybrid

      First solve the equivalent normal equations (see :c:member:`primme_svds_normalequations`) and then
      refine the solution solving the augmented problem (see :c:member:`primme_svds_augmented`).
  
      With :c:member:`primme_svds_normalequations` :c:func:`primme_svds_set_method` sets
      |Smethod| to ``primme_svds_op_AtA`` if |Sm| is larger or equal than |Sn|, and to ``primme_svds_op_AAt``
      otherwise; and |SmethodStage2| is set to ``primme_svds_op_augmented``.
 
      The minimum residual norm that this method can achieve is :math:`\|A\|\epsilon`,
      where :math:`\epsilon` is the machine precision.
      However it may not return triplets with singular values smaller than :math:`\|A\|\epsilon`
      if |Seps| is smaller than :math:`\|A\|\epsilon\sigma^{-1}`.

 .. _error-codes-svds:

Error Codes
-----------

The functions :c:func:`dprimme_svds` and :c:func:`zprimme_svds` return one of the next values:

*  0: success,
*  1: reported only amount of required memory,
* -1: failed in allocating int or real workspace,
* -2: malloc failed in allocating a permutation integer array,
* -3: main_iter() encountered problem; the calling stack of the functions where the error occurred was printed in 'stderr',
* -4: ``primme_svds`` is NULL,
* -5: Wrong value for |Sm| or |Sn| or |SmLocal| or |SnLocal|,
* -6: Wrong value for |SnumProcs|,
* -7: |SmatrixMatvec| is not set,
* -8: |SapplyPreconditioner| is not set but |Sprecondition| == 1 ,
* -9: |SnumProcs| >1 but |SglobalSumReal| is not set,
* -10: Wrong value for |SnumSvals|, it's larger than min(|Sm|, |Sn|),
* -11: Wrong value for |SnumSvals|, it's smaller than 1,
* -13: Wrong value for |Starget|,
* -14: Wrong value for |Smethod|,
* -15: Not supported combination of method and |SmethodStage2|,
* -16: Wrong value for |SprintLevel|,
* -17: ``svals`` is not set,
* -18: ``svecs`` is not set,
* -19: ``resNorms`` is not set
* -20: not enough memory for |SrealWork|
* -21: not enough memory for |SintWork|
* -100 up to -199: eigensolver error from first stage; see the value plus 100 in :ref:`error-codes`.
* -200 up to -299: eigensolver error from second stage; see the value plus 200 in :ref:`error-codes`.

.. include:: epilog.inc
