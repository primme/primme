
Appendix
--------

primme_params
"""""""""""""

.. c:type:: primme_params

   Structure to set the problem matrices and eigensolver options.

   .. c:member:: int n

      Dimension of the matrix.

   .. c:member:: void (*matrixMatvec) (void *x, void *y, int *blockSize, primme_params *primme)

      Block matrix-multivector multiplication, :math:`y = A x` in solving :math:`A x = \lambda x` or :math:`A x = \lambda B x`.
   
      :param x:
      :param y: one dimensional array containing the ``blockSize`` vectors 
         packed one after the other (i.e., the leading dimension is the vector size), each of size |nLocal|.
         The real type is ``double*`` and ``Complex_Z*`` when called from :c:func:`dprimme` and :c:func:`zprimme` respectively.
      :param blockSize: number of vectors in x and y.
      :param primme: parameters structure.

      .. note::

         Argument ``blockSize`` is passed by reference to make easier the interface to other
         languages (like Fortran).

   .. c:member:: void (*applyPreconditioner)(void *x, void *y, int *blockSize, struct primme_params *primme)

      Block preconditioner-multivector application, :math:`y = M^{-1}x` where :math:`M` is usually an approximation of :math:`A - \sigma I` or :math:`A - \sigma B` for finding eigenvalues close to :math:`\sigma`.
      The function follows the convention of |matrixMatvec|.

 
   .. c:member:: void (*massMatrixMatvec)(void *x, void *y, int *blockSize, struct primme_params *primme)

      Block matrix-multivector multiplication, :math:`y = B x` in solving :math:`A x = \lambda B x`.
      The function follows the convention of |matrixMatvec|.

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
      The default value is |n| if |numProcs| is 1.

   .. c:member:: void *commInfo

      A pointer to whatever parallel environment structures needed.
      For example, with MPI, it could be a pointer to the MPI communicator.
      PRIMME does not use this. It is available for possible use in 
      user functions defined in |matrixMatvec|,
      |applyPreconditioner|, |massMatrixMatvec| and
      |globalSumDouble|.
      The default values is NULL.

   .. c:member:: void (*globalSumDouble)(double *sendBuf, double *recvBuf, int *count, primme_params *primme)

      Global sum reduction function. 

      :param sendBuf: array of size count with the input local values.
      :param recvBuf: array of size count with the output global values
         so that i-th element of recvBuf is the sum over all processes of the i-th element
         of sendBuf.
      :param count: array size of sendBuf and recvBuf.
      :param primme: parameters structure.

      The default value is NULL if |numProcs| is 1.
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

      The default value is |primme_smallest|.

      .. note::
         * If some shift is close to the lower (higher) end of the spectrum,
           use either |primme_closest_geq| (|primme_closest_leq|) or
           |primme_closest_abs|.
         * |primme_closest_leq| and |primme_closest_geq| are more efficient
           than |primme_closest_abs|.
 
   .. c:member:: int numTargetShifts
 
      Size of the array |targetShifts|.
      Used only when |target| is |primme_closest_geq|,
      |primme_closest_leq| or |primme_closest_abs|.
      The default values is 0.

   .. c:member:: double *targetShifts

      Array of shifts, at least of size |numTargetShifts|.
      Used only when |target| is |primme_closest_geq|,
      |primme_closest_leq| or |primme_closest_abs|.
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
         |procID| 0 when |printLevel|
         is from 0 to 4.
         If |printLevel| is 5 output can be produced in any of
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
      (see |eps|).
      If it is less or equal to 0, it is used the largest absolute Ritz value seen.
      And on return, it is replaced with that value.

      The default value is 0.

   .. c:member:: double eps

      An eigenpairs is marked as converged when the 2-norm of the residual is less
      than |eps| times |aNorm|.
      The residual vector is :math:`A x - \lambda x` or :math:`A x - \lambda B x`.

      The default value is :math:`10^{-12}`.
 
   .. c:member:: FILE *outputFile

      Opened file to write down the output.

      The default value is the standard output.

   .. c:member:: int dynamicMethodSwitch

      If this value is 1, it alternates dynamically between |DEFAULT_MIN_TIME|
      and |DEFAULT_MIN_MATVECS|, trying to identify the fastest method.

      On exit, it holds a recommended method for future runs on this problem:

      * -1: use |DEFAULT_MIN_MATVECS| next time.
      * -2: use |DEFAULT_MIN_TIME| next time.
      * -3: close call, use |DYNAMIC| next time again.
      
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

      The default depends on the method and the value of some options.

   .. c:member:: int initSize
 
      On input, the number of initial vector guesses provided in ``evecs`` argument in :c:func:`dprimme`
      or :c:func:`zprimme`.
      On output, the number of converged eigenpairs.
      During execution, it holds the current number of converged eigenpairs.
      If in addition locking is used, these are accessible in ``evals`` and ``evecs``.

      The default value is 0.
      
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

      The default value is 0.

   .. c:member:: int maxBasisSize

      The maximum basis size allowed in the main iteration. This has memory
      implications.

      The default depends on method.

      .. note::

         For interior eigenvalues use a larger value than usual.

   .. c:member:: int minRestartSize

      Maximum Ritz vectors kept after restarting the basis.

      The default depends on |maxBasisSize|,
      |maxBlockSize| and method.

   .. c:member:: int maxBlockSize
 
      The maximum block size the code will try to use.

      The user should set
      this based on the architecture specifics of the target computer, 
      as well as any a priori knowledge of multiplicities. The code does 
      *not* require that |maxBlockSize| > 1 to find multiple eigenvalues. For some 
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

      If :c:func:`dprimme` or :c:func:`zprimme` is called with all arguments as NULL
      except for :c:type:`primme_params` then PRIMME returns immediately with |intWorkSize|
      containing the size *in bytes* of the integer workspace that will be required by the
      parameters set in PRIMME.

      Otherwise if |intWorkSize| is not 0, it should be the size of the integer work array
      *in bytes* that the user provides in |intWork|. If |intWorkSize| is 0, the code
      will allocate the required space, which can be freed later by calling :c:func:`primme_Free`.

      The default value is 0.

   .. c:member:: long int realWorkSize

      If :c:func:`dprimme` or :c:func:`zprimme` is called with all arguments as NULL
      except for :c:type:`primme_params` then PRIMME returns immediately with |realWorkSize|
      containing the size *in bytes* of the real workspace that will be required by the
      parameters set in PRIMME.

      Otherwise if |realWorkSize| is not 0, it should be the size of the real work array
      *in bytes* that the user provides in |realWork|. If |realWorkSize| is 0, the code
      will allocate the required space, which can be freed later by calling :c:func:`primme_Free`.

      The default value is 0.

   .. c:member:: int *intWork

      Integer work array.

      If NULL, the code will allocate its own workspace. If the provided space is not
      enough, the code will free it and allocate a new space.

      On exit, the first element shows if a locking problem has occurred.
      Using locking for large |numEvals| may, in some rare cases,
      cause some pairs to be practically converged, in the sense that their components 
      are in the basis of ``evecs``. If this is the case, a Rayleigh Ritz on returned
      ``evecs`` would provide the accurate eigenvectors (see [r4]_).

      The default value is NULL. 

   .. c:member:: void *realWork

      Real work array.

      If NULL, the code will allocate its own workspace. If the provided space is not
      enough, the code will free it and allocate a new space.

      The default value is NULL. 

   .. c:member:: int iseed

      The ``int iseed[4]`` is an array with the seeds needed by the LAPACK_ dlarnv and zlarnv.

      The default value is an array with values -1, -1, -1 and -1. In that case, ``iseed``
      is set based on the value of |procID| to avoid every process generating the same
      sequence of pseudorandom numbers.

   .. c:member:: void *matrix

      This field may be used to pass any required information 
      in the matrix-vector product |matrixMatvec|.

      The default value is NULL.
      
   .. c:member:: void *preconditioner

      This field may be used to pass any required information 
      in the matrix-vector product |applyPreconditioner|.

      The default value is NULL.

   .. c:member:: double *ShiftsForPreconditioner

      Array of size ``blockSize`` provided during execution of :c:member:dprimme and :c:member:zprimme holding
      the shifts to be used (if needed) in the preconditioning operation.

      For example if the block size is 3,
      there will be an array of three shifts in |ShiftsForPreconditioner|.
      Then the user can invert a shifted preconditioner for each of the 
      block vectors :math:`(M-ShiftsForPreconditioner_i)^{-1} x_i`.
      Classical Davidson (diagonal) preconditioning is an example of this.
   
   .. c:member:: primme_restartscheme restartingParams.scheme

      Select a restarting strategy:

      * ``primme_thick``, Thick restarting. This is the most efficient and robust
        in the general case.
      * ``primme_dtr``, Dynamic thick restarting. Helpful without 
        preconditioning but it is expensive to implement.

      The default value is |primme_thick|.

   .. c:member:: int restartingParams.maxPrevRetain

      Number of approximations from previous iteration to be retained
      after restart (see [r2]_). The restart size is |minRestartSize|
      plus |maxPrevRetain|.

      The default value is 1.

   .. c:member:: int correctionParams.precondition

      Set to 1 to use preconditioning.
      Make sure |applyPreconditioner| is not NULL then!

      The default value is 0.

   .. c:member:: int correctionParams.robustShifts

      Set to 1 to use robust shifting. It tries to avoid stagnation and 
      misconvergence by providing as shifts in |ShiftsForPreconditioner|
      the Ritz values displaced by an approximation of the eigenvalue error.

      The default value depends on method.

   .. c:member:: int correctionParams.maxInnerIterations

      Control the maximum number of inner QMR iterations:

      * 0:  no inner iterations;
      * >0: perform at most that number of inner iterations per outer step;
      * <0: perform at most the rest of the remaining matrix-vector products
        up to reach |maxMatvecs|.

      The default value depends on method.

      See also |convTest|.

   .. c:member:: double correctionParams.relTolBase

      Parameter used when |convTest|
      is |primme_decreasing_LTolerance|.

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

      The default value depends on method.

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

      Given the current selected Ritz value :math:`\Lambda` and vectors
      :math:`X`, the residual associated vectors :math:`R=AX-X\Lambda`, the previous locked
      vectors :math:`Q` and the preconditioner :math:`M^{-1}`.
      The correction :math:`D` appended to the basis in GD
      (when |maxInnerIterations| is 0) is:

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

      Where :math:`\Delta` is a diagonal matrix that :math:`\Delta_{i,i}` holds an estimation
      of the error of the approximate eigenvalue :math:`\Lambda_{i,i}`.
 
      The values of ``RightQ``, ``SkewQ``, ``LeftX`` and ``LeftQ`` are ignored.

      The correction :math:`D` in JD
      (when |maxInnerIterations| isn't 0)
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

      See [r3]_ for a study about different projector configuration in JD.

   .. c:member:: int stats.numOuterIterations

      Hold the number of outer iterations. The value is available during execution and at the end.

   .. c:member:: int stats.numRestarts

      Hold the number of restarts during execution and at the end.

   .. c:member:: int stats.numMatvecs

      Hold how many vectors the operator in |matrixMatvec| has been applied on.
      The value is available during execution and at the end.

   .. c:member:: int stats.numPreconds

      Hold how many vectors the operator in |applyPreconditioner| has been applied on.
      The value is available during execution and at the end.

   .. c:member:: int stats.elapsedTime

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
      |dynamicMethodSwitch| = 1 and makes the sames changes as
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
