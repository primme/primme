
FORTRAN Library Interface
-------------------------

The next enumerations and functions are declared in ``primme_f77.h``.

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

   :param evals(*): (output) array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): double precision

   :param resNorms(*): (output) array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): double precision

   :param evecs(*): (input/output) array at least of size |nLocal| times |numEvals|
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): double precision

   :param ptr primme: parameters structure.

   :param integer ierr: (output) error indicator; see the returned value of function :c:func:`dprimme`.

.. c:function:: zprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblems. The arguments have the
   same meaning like in function :c:func:`dprimme_f77`.

   :param evals(*): (output) 
   :type evals(*): double precision

   :param resNorms(*): (output)
   :type resNorms(*): double precision

   :param evecs(*): (input/output) 
   :type evecs(*): complex double precision

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
      * ``PRIMMEF77_correctionParams_projectors_LeftQ``,     in field :c:member:`primme_params.correctionParams.projectors.LeftQ`.
      * ``PRIMMEF77_correctionParams_projectors_LeftX``,     in field :c:member:`primme_params.correctionParams.projectors.LeftX`.
      * ``PRIMMEF77_correctionParams_projectors_RightQ``,    in field :c:member:`primme_params.correctionParams.projectors.RightQ`.
      * ``PRIMMEF77_correctionParams_projectors_RightX``,    in field :c:member:`primme_params.correctionParams.projectors.RightX`.
      * ``PRIMMEF77_correctionParams_projectors_SkewQ``,     in field :c:member:`primme_params.correctionParams.projectors.SkewQ`.
      * ``PRIMMEF77_correctionParams_projectors_SkewX``,     in field :c:member:`primme_params.correctionParams.projectors.SkewX`.
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


   .. note::

      When ``label`` is one of ``PRIMMEF77_matrixMatvec``, ``PRIMMEF77_applyPreconditioner``, ``PRIMMEF77_commInfo``,
      ``PRIMMEF77_intWork``, ``PRIMMEF77_realWork``, ``PRIMMEF77_matrix`` and ``PRIMMEF77_preconditioner``,
      the returned ``value`` is a C pointer (``void*``). Use Fortran pointer or other extensions to deal with it.
      For instance::

         use iso_c_binding
         MPI_Comm comm

         comm = MPI_COMM_WORLD
         call primme_set_member_f77(primme, PRIMMEF77_commInfo, comm)
         ...
         subroutine par_GlobalSumDouble(x,y,k,primme)
         use iso_c_binding
         implicit none
         ...
         MPI_Comm, pointer :: comm
         type(c_ptr) :: pcomm

         call primme_get_member_f77(primme, PRIMMEF77_commInfo, pcomm)
         call c_f_pointer(pcomm, comm)
         call MPI_Allreduce(x,y,k,MPI_DOUBLE,MPI_SUM,comm,ierr)
          
         
.. c:function:: primmetop_get_prec_shift_f77(primme, index, value)

   Get the value in some position of the array |ShiftsForPreconditioner|.

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
      * ``PRIMMEF77_correctionParams_projectors_LeftQ``,     in field :c:member:`primme_params.correctionParams.projectors.LeftQ`.
      * ``PRIMMEF77_correctionParams_projectors_LeftX``,     in field :c:member:`primme_params.correctionParams.projectors.LeftX`.
      * ``PRIMMEF77_correctionParams_projectors_RightQ``,    in field :c:member:`primme_params.correctionParams.projectors.RightQ`.
      * ``PRIMMEF77_correctionParams_projectors_RightX``,    in field :c:member:`primme_params.correctionParams.projectors.RightX`.
      * ``PRIMMEF77_correctionParams_projectors_SkewQ``,     in field :c:member:`primme_params.correctionParams.projectors.SkewQ`.
      * ``PRIMMEF77_correctionParams_projectors_SkewX``,     in field :c:member:`primme_params.correctionParams.projectors.SkewX`.
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

      Use this function exclusively inside the function |matrixMatvec|, |massMatrixMatvec|, or |applyPreconditioner|.
      Otherwise use the function :c:func:`primmetop_set_member_f77`.

.. c:function:: primme_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: (input) field where to get value. One of
      the detailed in function :c:func:`primmetop_set_member_f77`.

   :param value: (output) value of the field.

   .. note::

      Use this function exclusively inside the function |matrixMatvec|, |massMatrixMatvec|, or |applyPreconditioner|.
      Otherwise use the function :c:func:`primmetop_get_member_f77`.

   .. note::

      When ``label`` is one of ``PRIMMEF77_matrixMatvec``, ``PRIMMEF77_applyPreconditioner``, ``PRIMMEF77_commInfo``,
      ``PRIMMEF77_intWork``, ``PRIMMEF77_realWork``, ``PRIMMEF77_matrix`` and ``PRIMMEF77_preconditioner``,
      the returned ``value`` is a C pointer (``void*``). Use Fortran pointer or other extensions to deal with it.
      For instance::

         use iso_c_binding
         MPI_Comm comm

         comm = MPI_COMM_WORLD
         call primme_set_member_f77(primme, PRIMMEF77_commInfo, comm)
         ...
         subroutine par_GlobalSumDouble(x,y,k,primme)
         use iso_c_binding
         implicit none
         ...
         MPI_Comm, pointer :: comm
         type(c_ptr) :: pcomm

         call primme_get_member_f77(primme, PRIMMEF77_commInfo, pcomm)
         call c_f_pointer(pcomm, comm)
         call MPI_Allreduce(x,y,k,MPI_DOUBLE,MPI_SUM,comm,ierr)
 
.. c:function:: primme_get_prec_shift_f77(primme, index, value)

   Get the value in some position of the array |ShiftsForPreconditioner|.

   :param ptr primme: (input) parameters structure.

   :param integer index: (input) position of the array; the first position is 1.

   :param value: (output) value of the array at that position.

   .. note::

      Use this function exclusively inside the function |matrixMatvec|, |massMatrixMatvec|, or |applyPreconditioner|.
      Otherwise use the function :c:func:`primmetop_get_prec_shift_f77`.

.. include:: epilog.inc
