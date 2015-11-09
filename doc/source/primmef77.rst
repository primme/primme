
FORTRAN Library Interface
-------------------------

The next enumerations and functions are declared in ``primme_f77.h``.

.. c:type:: ptr

   Fortran datatype with the same size as a pointer.
   Use ``integer*4`` when compiling in 32 bits and ``integer*8`` in 64 bits.

primme_initialize_f77
"""""""""""""""""""""

.. c:function:: primme_initialize_f77(primme)

   Set PRIMME parameters structure to the default values.

   :param ptr primme: (output) parameters structure.

primme_set_method_f77
"""""""""""""""""""""

.. c:function:: primme_set_method_f77(method, primme, ierr)

   Set PRIMME parameters to one of the preset configurations.

   :param integer method: (input) preset configuration. One of:

      | ``PRIMMEF77_DYNAMIC``
      | ``PRIMMEF77_DEFAULT_MIN_TIME``
      | ``PRIMMEF77_DEFAULT_MIN_MATVECS``
      | ``PRIMMEF77_Arnoldi``
      | ``PRIMMEF77_GD``
      | ``PRIMMEF77_GD_plusK``
      | ``PRIMMEF77_GD_Olsen_plusK``
      | ``PRIMMEF77_JD_Olsen_plusK``
      | ``PRIMMEF77_RQI``
      | ``PRIMMEF77_JDQR``
      | ``PRIMMEF77_JDQMR``
      | ``PRIMMEF77_JDQMR_ETol``
      | ``PRIMMEF77_SUBSPACE_ITERATION``
      | ``PRIMMEF77_LOBPCG_OrthoBasis``
      | ``PRIMMEF77_LOBPCG_OrthoBasis_Window``

      See :c:type:`primme_preset_method`.

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) if 0, successful; if negative, something went wrong.

primme_Free_f77
"""""""""""""""

.. c:function:: primme_Free_f77(primme)

   Free memory allocated by PRIMME.

   :param ptr primme: parameters structure.

dprimme_f77
"""""""""""

.. c:function:: dprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a real symmetric standard eigenproblem.

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

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

zprimme_f77
"""""""""""

.. c:function:: zprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblem. The arguments have the
   same meaning as in function :c:func:`dprimme_f77`.

   :param evals(*): (output) 
   :type evals(*): double precision

   :param resNorms(*): (output)
   :type resNorms(*): double precision

   :param evecs(*): (input/output) 
   :type evecs(*): complex double precision

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

primmetop_set_member_f77
""""""""""""""""""""""""

.. c:function:: primmetop_set_member_f77(primme, label, value)

   Set a value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: field where to set value. One of:

      | :c:member:`PRIMMEF77_n                                   <primme_params.n>`
      | :c:member:`PRIMMEF77_matrixMatvec                        <primme_params.matrixMatvec>`
      | :c:member:`PRIMMEF77_applyPreconditioner                 <primme_params.applyPreconditioner>`
      | :c:member:`PRIMMEF77_numProcs                            <primme_params.numProcs>`
      | :c:member:`PRIMMEF77_procID                              <primme_params.procID>`
      | :c:member:`PRIMMEF77_commInfo                            <primme_params.commInfo>`
      | :c:member:`PRIMMEF77_nLocal                              <primme_params.nLocal>`
      | :c:member:`PRIMMEF77_globalSumDouble                     <primme_params.globalSumDouble>`
      | :c:member:`PRIMMEF77_numEvals                            <primme_params.numEvals>`
      | :c:member:`PRIMMEF77_target                              <primme_params.target>`
      | :c:member:`PRIMMEF77_numTargetShifts                     <primme_params.numTargetShifts>`
      | :c:member:`PRIMMEF77_targetShifts                        <primme_params.targetShifts>`
      | :c:member:`PRIMMEF77_locking                             <primme_params.locking>`
      | :c:member:`PRIMMEF77_initSize                            <primme_params.initSize>`
      | :c:member:`PRIMMEF77_numOrthoConst                       <primme_params.numOrthoConst>`
      | :c:member:`PRIMMEF77_maxBasisSize                        <primme_params.maxBasisSize>`
      | :c:member:`PRIMMEF77_minRestartSize                      <primme_params.minRestartSize>`
      | :c:member:`PRIMMEF77_maxBlockSize                        <primme_params.maxBlockSize>`
      | :c:member:`PRIMMEF77_maxMatvecs                          <primme_params.maxMatvecs>`
      | :c:member:`PRIMMEF77_maxOuterIterations                  <primme_params.maxOuterIterations>`
      | :c:member:`PRIMMEF77_intWorkSize                         <primme_params.intWorkSize>`
      | :c:member:`PRIMMEF77_realWorkSize                        <primme_params.realWorkSize>`
      | :c:member:`PRIMMEF77_iseed                               <primme_params.iseed>`
      | :c:member:`PRIMMEF77_intWork                             <primme_params.intWork>`
      | :c:member:`PRIMMEF77_realWork                            <primme_params.realWork>`
      | :c:member:`PRIMMEF77_aNorm                               <primme_params.aNorm>`
      | :c:member:`PRIMMEF77_eps                                 <primme_params.eps>`
      | :c:member:`PRIMMEF77_printLevel                          <primme_params.printLevel>`
      | :c:member:`PRIMMEF77_outputFile                          <primme_params.outputFile>`
      | :c:member:`PRIMMEF77_matrix                              <primme_params.matrix>`
      | :c:member:`PRIMMEF77_preconditioner                      <primme_params.preconditioner>`
      | :c:member:`PRIMMEF77_restartingParams_scheme             <primme_params.restartingParams.scheme>`.
      | :c:member:`PRIMMEF77_restartingParams_maxPrevRetain      <primme_params.restartingParams.maxPrevRetain>`
      | :c:member:`PRIMMEF77_correctionParams_precondition       <primme_params.correctionParams.precondition>`
      | :c:member:`PRIMMEF77_correctionParams_robustShifts       <primme_params.correctionParams.robustShifts>`
      | :c:member:`PRIMMEF77_correctionParams_maxInnerIterations <primme_params.correctionParams.maxInnerIterations>`
      | :c:member:`PRIMMEF77_correctionParams_projectors_LeftQ   <primme_params.correctionParams.projectors.LeftQ>`
      | :c:member:`PRIMMEF77_correctionParams_projectors_LeftX   <primme_params.correctionParams.projectors.LeftX>`
      | :c:member:`PRIMMEF77_correctionParams_projectors_RightQ  <primme_params.correctionParams.projectors.RightQ>`
      | :c:member:`PRIMMEF77_correctionParams_projectors_RightX  <primme_params.correctionParams.projectors.RightX>`
      | :c:member:`PRIMMEF77_correctionParams_projectors_SkewQ   <primme_params.correctionParams.projectors.SkewQ>`
      | :c:member:`PRIMMEF77_correctionParams_projectors_SkewX   <primme_params.correctionParams.projectors.SkewX>`
      | :c:member:`PRIMMEF77_correctionParams_convTest           <primme_params.correctionParams.convTest>`
      | :c:member:`PRIMMEF77_correctionParams_relTolBase         <primme_params.correctionParams.relTolBase>`
      | :c:member:`PRIMMEF77_stats_numOuterIterations            <primme_params.stats.numOuterIterations>`
      | :c:member:`PRIMMEF77_stats_numRestarts                   <primme_params.stats.numRestarts>`
      | :c:member:`PRIMMEF77_stats_numMatvecs                    <primme_params.stats.numMatvecs>`
      | :c:member:`PRIMMEF77_stats_numPreconds                   <primme_params.stats.numPreconds>`
      | :c:member:`PRIMMEF77_stats_elapsedTime                   <primme_params.stats.elapsedTime>`
      | :c:member:`PRIMMEF77_dynamicMethodSwitch                 <primme_params.dynamicMethodSwitch>`
      | :c:member:`PRIMMEF77_massMatrixMatvec                    <primme_params.massMatrixMatvec>`

   :param value: (input) value to set.

   .. note::

      **Don't use** this function inside PRIMME's callback functions, e.g., |matrixMatvec| or
      |applyPreconditioner|, or in functions called by these functions. In those cases use
      :c:func:`primme_set_member_f77`.


primmetop_get_member_f77
""""""""""""""""""""""""

.. c:function:: primmetop_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: (input) field where to get value. One of
      the detailed in function :c:func:`primmetop_set_member_f77`.

   :param value: (output) value of the field.

   .. note::

      **Don't use** this function inside PRIMME's callback functions, e.g., |matrixMatvec| or
      |applyPreconditioner|, or in functions called by these functions. In those cases use
      :c:func:`primme_get_member_f77`.

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

      Most users would not need to retrieve these pointers in their programs.
    
primmetop_get_prec_shift_f77
""""""""""""""""""""""""""""
         
.. c:function:: primmetop_get_prec_shift_f77(primme, index, value)

   Get the value in some position of the array |ShiftsForPreconditioner|.

   :param ptr primme: (input) parameters structure.

   :param integer index: (input) position of the array; the first position is 1.

   :param value: (output) value of the array at that position.

primme_set_member_f77
"""""""""""""""""""""

.. c:function:: primme_set_member_f77(primme, label, value)

   Set a value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: field where to set value. One of the vales defined
       in :c:func:`primmetop_set_member_f77`.

   :param value: (input) value to set.

   .. note::

      Use this function exclusively inside PRIMME's callback functions, e.g., |matrixMatvec|
      or |applyPreconditioner|, or in functions called by these functions. Otherwise, e.g.,
      from the main program, use the function :c:func:`primmetop_set_member_f77`.

primme_get_member_f77
"""""""""""""""""""""

.. c:function:: primme_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: (input) field where to get value. One of
      the detailed in function :c:func:`primmetop_set_member_f77`.

   :param value: (output) value of the field.

   .. note::

      Use this function exclusively inside PRIMME's callback functions, e.g., |matrixMatvec|
      or |applyPreconditioner|, or in functions called by these functions. Otherwise, e.g.,
      from the main program, use the function :c:func:`primmetop_get_member_f77`.

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

      Most users would not need to retrieve these pointers in their programs.

primme_get_prec_shift_f77
"""""""""""""""""""""""""
 
.. c:function:: primme_get_prec_shift_f77(primme, index, value)

   Get the value in some position of the array |ShiftsForPreconditioner|.

   :param ptr primme: (input) parameters structure.

   :param integer index: (input) position of the array; the first position is 1.

   :param value: (output) value of the array at that position.

   .. note::

      Use this function exclusively inside the function |matrixMatvec|, |massMatrixMatvec|, or |applyPreconditioner|.
      Otherwise use the function :c:func:`primmetop_get_prec_shift_f77`.

.. include:: epilog.inc
