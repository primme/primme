Appendix
--------

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
      respectively. Elsewhere they have dimensions |SmLocal| x ``blockSize`` and |SnLocal| x ``blockSize``. Both arrays are column-major
      (consecutive rows are consecutive in memory).

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

      Block preconditioner-multivector application. Depending on ``mode`` it is expected an approximation of the inverse of

      * ``primme_svds_op_AtA``: :math:`y = A^*Ax - \sigma^2 I`,
      * ``primme_svds_op_AAt``: :math:`y = AA^*x - \sigma^2 I`,
      * ``primme_svds_op_augmented``: :math:`\left(\begin{array}{cc} 0 & A^* \\ A & 0 \end{array}\right) - \sigma I`.

      Where :math:`\sigma` is the current target (see |targetShifts|) (for finding the smallest :math:`\sigma` is zero).

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
      Both arrays are column-major (consecutive rows are consecutive in memory).

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

      When MPI is used, this can be a simply wrapper to MPI_Allreduce() as shown below:

      .. code:: c

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

      The level of message reporting from the code. For now it controls the reporting
      level of the underneath eigensolvers. See |printLevel| in primme_params.

      All output is writen in |SoutputFile|.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 1;
         | this field is read by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: double aNorm

      An estimate of the 2-norm of :math:`A`, which is used in the default convergence
      criterion (see |Seps|).

      If |aNorm| is less than or equal to 0, the code uses the largest absolute
      Ritz value seen. On return, |SaNorm| is then replaced with that value.

      Input/output:

         | :c:func:`primme_svds_initialize` sets this field to 0.0;
         | this field is read and written by :c:func:`dprimme_svds` and :c:func:`zprimme_svds`.

   .. c:member:: double eps

      A triplet is marked as converged when the 2-norm
      of the residual vectors is less than |Seps| \* |SaNorm|.
      The residual vectors are :math:`A v - \sigma u` and :math:`A^* u - \sigma v` for the
      triplet :math:`(u,\sigma,v)`.

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

 .. _error-codes-svds:

Error Codes
"""""""""""

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

.. _methods_svds:

Preset Methods
""""""""""""""

.. c:type:: primme_svds_preset_method

   .. c:member:: primme_svds_default

      Set as :c:member:`primme_svds_hybrid`.

   .. c:member:: primme_svds_normalequations

      Solve the equivalent eigenvalue problem :math:`A^*A V = \Sigma^2 V` and computes :math:`U` by normalizing
      the vectors :math:`AV`. If |Sm| is smaller than |Sn|, :math:`AA^*` is solved instead.
  
      With :c:member:`primme_svds_normalequations` :c:func:`primme_svds_set_method` sets
      |Smethod| to ``primme_svds_op_AtA`` if |Sm| is larger or equal than |Sn|, and to ``primme_svds_op_AAt``
      otherwise; and |SmethodStage2| is set to ``primme_svds_op_none``.
 
   .. c:member:: primme_svds_augmented

      Solve the equivalent eigenvalue problem :math:`\left(\begin{array}{cc} 0 & A^* \\ A & 0 \end{array}\right) X = \sigma X`
      with :math:`X = \left(\begin{array}{cc}V\\U\end{array}\right)`.
  
      With :c:member:`primme_svds_augmented` :c:func:`primme_svds_set_method` sets
      |Smethod| to ``primme_svds_op_augmented`` and |SmethodStage2| to ``primme_svds_op_none``.
 
   .. c:member:: primme_svds_hybrid

      First solve the equivalent normal equations (see :c:member:`primme_svds_normalequations`) and then
      refine the solution solving the augmented problem (see :c:member:`primme_svds_augmented`).
  
      With :c:member:`primme_svds_normalequations` :c:func:`primme_svds_set_method` sets
      |Smethod| to ``primme_svds_op_AtA`` if |Sm| is larger or equal than |Sn|, and to ``primme_svds_op_AAt``
      otherwise; and |SmethodStage2| is set to ``primme_svds_op_augmented``.
 
.. include:: epilog.inc
