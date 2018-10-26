
.. highlight:: fortran

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

      | ``PRIMME_DYNAMIC``
      | ``PRIMME_DEFAULT_MIN_TIME``
      | ``PRIMME_DEFAULT_MIN_MATVECS``
      | ``PRIMME_Arnoldi``
      | ``PRIMME_GD``
      | ``PRIMME_GD_plusK``
      | ``PRIMME_GD_Olsen_plusK``
      | ``PRIMME_JD_Olsen_plusK``
      | ``PRIMME_RQI``
      | ``PRIMME_JDQR``
      | ``PRIMME_JDQMR``
      | ``PRIMME_JDQMR_ETol``
      | ``PRIMME_STEEPEST_DESCENT``
      | ``PRIMME_LOBPCG_OrthoBasis``
      | ``PRIMME_LOBPCG_OrthoBasis_Window``

      See :c:type:`primme_preset_method`.

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) if 0, successful; if negative, something went wrong.

primme_free_f77
"""""""""""""""

.. c:function:: primme_free_f77(primme)

   Free memory allocated by PRIMME and delete all values set.

   :param ptr primme: (input/output) parameters structure.

sprimme_f77
"""""""""""

.. c:function:: sprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a real symmetric standard eigenproblem using single precision.

   :param evals(*): (output) array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): real

   :param resNorms(*): (output) array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): real

   :param evecs(*): (input/output) array at least of size |nLocal| times |numEvals|
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): real

   :param ptr primme: parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

cprimme_f77
"""""""""""

.. c:function:: cprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblem. The arguments have the
   same meaning as in function :c:func:`sprimme_f77`.

   :param evals(*): (output) 
   :type evals(*): real

   :param resNorms(*): (output)
   :type resNorms(*): real

   :param evecs(*): (input/output) 
   :type evecs(*): complex real

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

dprimme_f77
"""""""""""

.. c:function:: dprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a real symmetric standard eigenproblem using double precision.

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

primme_set_member_f77
""""""""""""""""""""""""

.. c:function:: primme_set_member_f77(primme, label, value)

   Set a value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: field where to set value. One of:

      | :c:member:`PRIMME_n                                   <primme_params.n>`
      | :c:member:`PRIMME_matrixMatvec                        <primme_params.matrixMatvec>`
      | :c:member:`PRIMME_applyPreconditioner                 <primme_params.applyPreconditioner>`
      | :c:member:`PRIMME_numProcs                            <primme_params.numProcs>`
      | :c:member:`PRIMME_procID                              <primme_params.procID>`
      | :c:member:`PRIMME_commInfo                            <primme_params.commInfo>`
      | :c:member:`PRIMME_nLocal                              <primme_params.nLocal>`
      | :c:member:`PRIMME_globalSumReal                       <primme_params.globalSumReal>`
      | :c:member:`PRIMME_numEvals                            <primme_params.numEvals>`
      | :c:member:`PRIMME_target                              <primme_params.target>`
      | :c:member:`PRIMME_numTargetShifts                     <primme_params.numTargetShifts>`
      | :c:member:`PRIMME_targetShifts                        <primme_params.targetShifts>`
      | :c:member:`PRIMME_locking                             <primme_params.locking>`
      | :c:member:`PRIMME_initSize                            <primme_params.initSize>`
      | :c:member:`PRIMME_numOrthoConst                       <primme_params.numOrthoConst>`
      | :c:member:`PRIMME_maxBasisSize                        <primme_params.maxBasisSize>`
      | :c:member:`PRIMME_minRestartSize                      <primme_params.minRestartSize>`
      | :c:member:`PRIMME_maxBlockSize                        <primme_params.maxBlockSize>`
      | :c:member:`PRIMME_maxMatvecs                          <primme_params.maxMatvecs>`
      | :c:member:`PRIMME_maxOuterIterations                  <primme_params.maxOuterIterations>`
      | :c:member:`PRIMME_intWorkSize                         <primme_params.intWorkSize>`
      | :c:member:`PRIMME_realWorkSize                        <primme_params.realWorkSize>`
      | :c:member:`PRIMME_iseed                               <primme_params.iseed>`
      | :c:member:`PRIMME_intWork                             <primme_params.intWork>`
      | :c:member:`PRIMME_realWork                            <primme_params.realWork>`
      | :c:member:`PRIMME_aNorm                               <primme_params.aNorm>`
      | :c:member:`PRIMME_eps                                 <primme_params.eps>`
      | :c:member:`PRIMME_printLevel                          <primme_params.printLevel>`
      | :c:member:`PRIMME_outputFile                          <primme_params.outputFile>`
      | :c:member:`PRIMME_matrix                              <primme_params.matrix>`
      | :c:member:`PRIMME_preconditioner                      <primme_params.preconditioner>`
      | :c:member:`PRIMME_restartingParams_scheme             <primme_params.restartingParams.scheme>`.
      | :c:member:`PRIMME_restartingParams_maxPrevRetain      <primme_params.restartingParams.maxPrevRetain>`
      | :c:member:`PRIMME_correctionParams_precondition       <primme_params.correctionParams.precondition>`
      | :c:member:`PRIMME_correctionParams_robustShifts       <primme_params.correctionParams.robustShifts>`
      | :c:member:`PRIMME_correctionParams_maxInnerIterations <primme_params.correctionParams.maxInnerIterations>`
      | :c:member:`PRIMME_correctionParams_projectors_LeftQ   <primme_params.correctionParams.projectors.LeftQ>`
      | :c:member:`PRIMME_correctionParams_projectors_LeftX   <primme_params.correctionParams.projectors.LeftX>`
      | :c:member:`PRIMME_correctionParams_projectors_RightQ  <primme_params.correctionParams.projectors.RightQ>`
      | :c:member:`PRIMME_correctionParams_projectors_RightX  <primme_params.correctionParams.projectors.RightX>`
      | :c:member:`PRIMME_correctionParams_projectors_SkewQ   <primme_params.correctionParams.projectors.SkewQ>`
      | :c:member:`PRIMME_correctionParams_projectors_SkewX   <primme_params.correctionParams.projectors.SkewX>`
      | :c:member:`PRIMME_correctionParams_convTest           <primme_params.correctionParams.convTest>`
      | :c:member:`PRIMME_correctionParams_relTolBase         <primme_params.correctionParams.relTolBase>`
      | :c:member:`PRIMME_stats_numOuterIterations            <primme_params.stats.numOuterIterations>`
      | :c:member:`PRIMME_stats_numRestarts                   <primme_params.stats.numRestarts>`
      | :c:member:`PRIMME_stats_numMatvecs                    <primme_params.stats.numMatvecs>`
      | :c:member:`PRIMME_stats_numPreconds                   <primme_params.stats.numPreconds>`
      | :c:member:`PRIMME_stats_elapsedTime                   <primme_params.stats.elapsedTime>`
      | :c:member:`PRIMME_dynamicMethodSwitch                 <primme_params.dynamicMethodSwitch>`
      | :c:member:`PRIMME_massMatrixMatvec                    <primme_params.massMatrixMatvec>`

   :param value: (input) value to set.

      If the type of the option is integer (``int``, :c:type:`PRIMME_INT`, ``size_t``), the
      type of ``value`` should be as long as :c:type:`PRIMME_INT`, which is ``integer*8`` by default.

   .. note::

      **Don't use** this function inside PRIMME's callback functions, e.g., |matrixMatvec| or
      |applyPreconditioner|, or in functions called by these functions.


primmetop_get_member_f77
""""""""""""""""""""""""

.. c:function:: primmetop_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: (input) field where to get value. One of
      the detailed in function :c:func:`primmetop_set_member_f77`.

   :param value: (output) value of the field.

      If the type of the option is integer (``int``, :c:type:`PRIMME_INT`, ``size_t``), the
      type of ``value`` should be as long as :c:type:`PRIMME_INT`, which is ``integer*8`` by default.

   .. note::

      **Don't use** this function inside PRIMME's callback functions, e.g., |matrixMatvec| or
      |applyPreconditioner|, or in functions called by these functions. In those cases use
      :c:func:`primme_get_member_f77`.

   .. note::

      When ``label`` is one of ``PRIMME_matrixMatvec``, ``PRIMME_applyPreconditioner``, ``PRIMME_commInfo``,
      ``PRIMME_intWork``, ``PRIMME_realWork``, ``PRIMME_matrix`` and ``PRIMME_preconditioner``,
      the returned ``value`` is a C pointer (``void*``). Use Fortran pointer or other extensions to deal with it.
      For instance::

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

      Most users would not need to retrieve these pointers in their programs.
    
primmetop_get_prec_shift_f77
""""""""""""""""""""""""""""
         
.. c:function:: primmetop_get_prec_shift_f77(primme, index, value)

   Get the value in some position of the array |ShiftsForPreconditioner|.

   :param ptr primme: (input) parameters structure.

   :param integer index: (input) position of the array; the first position is 1.

   :param value: (output) value of the array at that position.

primme_get_member_f77
"""""""""""""""""""""

.. c:function:: primme_get_member_f77(primme, label, value)

   Get the value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: (input) field where to get value. One of
      the detailed in function :c:func:`primmetop_set_member_f77`.

   :param value: (output) value of the field.

      If the type of the option is integer (``int``, :c:type:`PRIMME_INT`, ``size_t``), the
      type of ``value`` should be as long as :c:type:`PRIMME_INT`, which is ``integer*8`` by default.

   .. note::

      Use this function exclusively inside PRIMME's callback functions, e.g., |matrixMatvec|
      or |applyPreconditioner|, or in functions called by these functions. Otherwise, e.g.,
      from the main program, use the function :c:func:`primmetop_get_member_f77`.

   .. note::

      When ``label`` is one of ``PRIMME_matrixMatvec``, ``PRIMME_applyPreconditioner``, ``PRIMME_commInfo``,
      ``PRIMME_intWork``, ``PRIMME_realWork``, ``PRIMME_matrix`` and ``PRIMME_preconditioner``,
      the returned ``value`` is a C pointer (``void*``). Use Fortran pointer or other extensions to deal with it.
      For instance::

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
