
.. highlight:: fortran

FORTRAN Library Interface
-------------------------

The next enumerations and functions are declared in ``primme_svds_f77.h``.

sprimme_svds_f77
""""""""""""""""

.. c:function:: sprimme_svds_f77(svals, svecs, resNorms, primme_svds)

   Solve a real singular value problem using single precision.

   :param svals(*): (output) array at least of size |SnumSvals| to store the
      computed singular values; all processes in a parallel run return this local array with the same values.
   :type svals(*): real

   :param resNorms(*): array at least of size |SnumSvals| to store the
      residual norms of the computed triplets; all processes in parallel run return this local array with
      the same values.
   :type resNorms(*): real

   :param svecs(*): array at least of size (|SmLocal| + |SnLocal|) times |SnumSvals|
      to store columnwise the (local part of the) computed left singular vectors
      and the right singular vectors.
   :type svecs(*): real

   :param ptr primme_svds: parameters structure.

   :return: error indicator; see :ref:`error-codes-svds`.

cprimme_svds_f77
""""""""""""""""

.. c:function:: cprimme_svds_f77(svals, svecs, resNorms, primme_svds)

   Solve a complex singular value problem using single precision.

   :param svals(*): (output) array at least of size |SnumSvals| to store the
      computed singular values; all processes in a parallel run return this local array with the same values.
   :type svals(*): real

   :param resNorms(*): array at least of size |SnumSvals| to store the
      residual norms of the computed triplets; all processes in parallel run return this local array with
      the same values.
   :type resNorms(*): real

   :param svecs(*): array at least of size (|SmLocal| + |SnLocal|) times |SnumSvals|
      to store columnwise the (local part of the) computed left singular vectors
      and the right singular vectors.
   :type svecs(*): complex

   :param ptr primme_svds: parameters structure.

   :return: error indicator; see :ref:`error-codes-svds`.

dprimme_svds_f77
""""""""""""""""

.. c:function:: dprimme_svds_f77(svals, svecs, resNorms, primme_svds)

   Solve a real singular value problem using double precision.

   :param svals(*): (output) array at least of size |SnumSvals| to store the
      computed singular values; all processes in a parallel run return this local array with the same values.
   :type svals(*): double precision

   :param resNorms(*): array at least of size |SnumSvals| to store the
      residual norms of the computed triplets; all processes in parallel run return this local array with
      the same values.
   :type resNorms(*): double precision

   :param svecs(*): array at least of size (|SmLocal| + |SnLocal|) times |SnumSvals|
      to store columnwise the (local part of the) computed left singular vectors
      and the right singular vectors.
   :type svecs(*): double precision

   :param ptr primme_svds: parameters structure.

   :return: error indicator; see :ref:`error-codes-svds`.

zprimme_svds_f77
""""""""""""""""

.. c:function:: zprimme_svds_f77(svals, svecs, resNorms, primme_svds)

   Solve a complex singular value problem using double precision.

   :param svals(*): (output) array at least of size |SnumSvals| to store the
      computed singular values; all processes in a parallel run return this local array with the same values.
   :type svals(*): double precision

   :param resNorms(*): array at least of size |SnumSvals| to store the
      residual norms of the computed triplets; all processes in parallel run return this local array with
      the same values.
   :type resNorms(*): double precision

   :param svecs(*): array at least of size (|SmLocal| + |SnLocal|) times |SnumSvals|
      to store columnwise the (local part of the) computed left singular vectors
      and the right singular vectors.
   :type svecs(*): complex*16

   :param ptr primme_svds: parameters structure.

   :return: error indicator; see :ref:`error-codes-svds`.

primme_svds_initialize_f77
""""""""""""""""""""""""""

.. c:function:: primme_svds_initialize_f77(primme_svds)

   Set PRIMME SVDS parameters structure to the default values.

   :param ptr primme_svds: (output) parameters structure.

primme_svds_set_method_f77
""""""""""""""""""""""""""

.. c:function:: primme_svds_set_method_f77(method, methodStage1, methodStage2, primme_svds, ierr)

   Set PRIMME SVDS parameters to one of the preset configurations.

   :param integer method: (input) preset configuration to compute the singular triplets; one of

      * |PRIMME_SVDS_default|, currently set as |PRIMME_SVDS_hybrid|.
      * |PRIMME_SVDS_normalequations|, compute the eigenvectors of :math:`A^*A` or :math:`A A^*`.
      * |PRIMME_SVDS_augmented|, compute the eigenvectors of the augmented matrix, :math:`\left(\begin{array}{cc} 0 & A^* \\ A & 0 \end{array}\right)`.
      * |PRIMME_SVDS_hybrid|, start with |PRIMME_SVDS_normalequations|; use the
        resulting approximate singular vectors as initial vectors for
        |PRIMME_SVDS_augmented| if the required accuracy was not achieved.

   :param primme_preset_method methodStage1: (input) preset method to compute the eigenpairs at the first stage; see available values at :c:func:`primme_set_method_f77`.

   :param primme_preset_method methodStage2: (input) preset method to compute the eigenpairs with
      the second stage of |PRIMME_SVDS_hybrid|; see available values at :c:func:`primme_set_method_f77`.

   :param ptr primme_svds: (input/output) parameters structure.

   :param integer ierr: (output) if 0, successful; if negative, something went wrong.

primme_svds_display_params_f77
""""""""""""""""""""""""""""""

.. c:function:: primme_svds_display_params_f77(primme_svds)

   Display all printable settings of ``primme_svds`` into the file descriptor |SoutputFile|.

   :param ptr primme_svds: (input) parameters structure.

primme_svds_free_f77
""""""""""""""""""""

.. c:function:: primme_svds_free_f77(primme_svds)

   Free memory allocated by PRIMME SVDS and delete all values set.

   :param ptr primme_svds: (input/output) parameters structure.

primme_svds_set_member_f77
"""""""""""""""""""""""""""""

.. c:function:: primme_svds_set_member_f77(primme_svds, label, value)

   Set a value in some field of the parameter structure.

   :param ptr primme_svds: (input) parameters structure.

   :param integer label: field where to set value. One of:

     | :c:member:`PRIMME_SVDS_primme                         <primme_svds_params.primme>`
     | :c:member:`PRIMME_SVDS_primmeStage2                   <primme_svds_params.primmeStage2>`
     | :c:member:`PRIMME_SVDS_m                              <primme_svds_params.m>`
     | :c:member:`PRIMME_SVDS_n                              <primme_svds_params.n>`
     | :c:member:`PRIMME_SVDS_matrixMatvec                   <primme_svds_params.matrixMatvec>`
     | :c:member:`PRIMME_SVDS_applyPreconditioner            <primme_svds_params.applyPreconditioner>`
     | :c:member:`PRIMME_SVDS_numProcs                       <primme_svds_params.numProcs>`
     | :c:member:`PRIMME_SVDS_procID                         <primme_svds_params.procID>`
     | :c:member:`PRIMME_SVDS_mLocal                         <primme_svds_params.mLocal>`
     | :c:member:`PRIMME_SVDS_nLocal                         <primme_svds_params.nLocal>`
     | :c:member:`PRIMME_SVDS_commInfo                       <primme_svds_params.commInfo>`
     | :c:member:`PRIMME_SVDS_globalSumReal                  <primme_svds_params.globalSumReal>`
     | :c:member:`PRIMME_SVDS_numSvals                       <primme_svds_params.numSvals>`
     | :c:member:`PRIMME_SVDS_target                         <primme_svds_params.target>`
     | :c:member:`PRIMME_SVDS_numTargetShifts                <primme_svds_params.numTargetShifts>`
     | :c:member:`PRIMME_SVDS_targetShifts                   <primme_svds_params.targetShifts>`
     | :c:member:`PRIMME_SVDS_method                         <primme_svds_params.method>`
     | :c:member:`PRIMME_SVDS_methodStage2                   <primme_svds_params.methodStage2>`
     | :c:member:`PRIMME_SVDS_intWorkSize                    <primme_svds_params.intWorkSize>`
     | :c:member:`PRIMME_SVDS_realWorkSize                   <primme_svds_params.realWorkSize>`
     | :c:member:`PRIMME_SVDS_intWork                        <primme_svds_params.intWork>`
     | :c:member:`PRIMME_SVDS_realWork                       <primme_svds_params.realWork>`
     | :c:member:`PRIMME_SVDS_matrix                         <primme_svds_params.matrix>`
     | :c:member:`PRIMME_SVDS_preconditioner                 <primme_svds_params.preconditioner>`
     | :c:member:`PRIMME_SVDS_locking                        <primme_svds_params.locking>`
     | :c:member:`PRIMME_SVDS_numOrthoConst                  <primme_svds_params.numOrthoConst>`
     | :c:member:`PRIMME_SVDS_aNorm                          <primme_svds_params.aNorm>`
     | :c:member:`PRIMME_SVDS_eps                            <primme_svds_params.eps>`
     | :c:member:`PRIMME_SVDS_precondition                   <primme_svds_params.precondition>`
     | :c:member:`PRIMME_SVDS_initSize                       <primme_svds_params.initSize>`
     | :c:member:`PRIMME_SVDS_maxBasisSize                   <primme_svds_params.maxBasisSize>`
     | :c:member:`PRIMME_SVDS_maxBlockSize                   <primme_svds_params.maxBlockSize>`
     | :c:member:`PRIMME_SVDS_maxMatvecs                     <primme_svds_params.maxMatvecs>`
     | :c:member:`PRIMME_SVDS_iseed                          <primme_svds_params.iseed>`
     | :c:member:`PRIMME_SVDS_printLevel                     <primme_svds_params.printLevel>`
     | :c:member:`PRIMME_SVDS_outputFile                     <primme_svds_params.outputFile>`
     | :c:member:`PRIMME_SVDS_stats_numOuterIterations       <primme_svds_params.stats_numOuterIterations>`
     | :c:member:`PRIMME_SVDS_stats_numRestarts              <primme_svds_params.stats_numRestarts>`
     | :c:member:`PRIMME_SVDS_stats_numMatvecs               <primme_svds_params.stats_numMatvecs>`
     | :c:member:`PRIMME_SVDS_stats_numPreconds              <primme_svds_params.stats_numPreconds>`
     | :c:member:`PRIMME_SVDS_stats_elapsedTime              <primme_svds_params.stats_elapsedTime>`

   :param value: (input) value to set.

   .. note::

      **Don't use** this function inside PRIMME SVDS's callback functions, e.g., |SmatrixMatvec| or
      |SapplyPreconditioner|, or in functions called by these functions.

primme_svdstop_get_member_f77
"""""""""""""""""""""""""""""

.. c:function:: primme_svdstop_get_member_f77(primme_svds, label, value)

   Get the value in some field of the parameter structure.

   :param ptr primme_svds: (input) parameters structure.

   :param integer label: (input) field where to get value. One of
      the detailed in function :c:func:`primmesvds_top_set_member_f77`.

   :param value: (output) value of the field.

   .. note::

      **Don't use** this function inside PRIMME SVDS's callback functions, e.g., |SmatrixMatvec| or
      |SapplyPreconditioner|, or in functions called by these functions. In those cases use
      :c:func:`primme_svds_get_member_f77`.

   .. note::

      When ``label`` is one of ``PRIMME_SVDS_matrixMatvec``, ``PRIMME_SVDS_applyPreconditioner``, ``PRIMME_SVDS_commInfo``,
      ``PRIMME_SVDS_intWork``, ``PRIMME_SVDS_realWork``, ``PRIMME_SVDS_matrix`` and ``PRIMME_SVDS_preconditioner``,
      the returned ``value`` is a C pointer (``void*``). Use Fortran pointer or other extensions to deal with it.
      For instance::

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

      Most users would not need to retrieve these pointers in their programs.

primme_svds_get_member_f77
""""""""""""""""""""""""""

.. c:function:: primme_svds_get_member_f77(primme_svds, label, value)

   Get the value in some field of the parameter structure.

   :param ptr primme_svds: (input) parameters structure.

   :param integer label: (input) field where to get value. One of
      the detailed in function :c:func:`primme_svdstop_set_member_f77`.

   :param value: (output) value of the field.

   .. note::

      Use this function exclusively inside PRIMME SVDS's callback functions, e.g., |SmatrixMatvec|
      or |SapplyPreconditioner|, or in functions called by these functions. Otherwise, e.g.,
      from the main program, use the function :c:func:`primme_svdstop_get_member_f77`.

   .. note::

      When ``label`` is one of ``PRIMME_SVDS_matrixMatvec``, ``PRIMME_SVDS_applyPreconditioner``, ``PRIMME_SVDS_commInfo``,
      ``PRIMME_SVDS_intWork``, ``PRIMME_SVDS_realWork``, ``PRIMME_SVDS_matrix`` and ``PRIMME_SVDS_preconditioner``,
      the returned ``value`` is a C pointer (``void*``). Use Fortran pointer or other extensions to deal with it.
      For instance::

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

      Most users would not need to retrieve these pointers in their programs.

.. include:: epilog.inc
