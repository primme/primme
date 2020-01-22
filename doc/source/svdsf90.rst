
.. highlight:: Fortran

.. f:currentmodule::

FORTRAN 90 Library Interface
----------------------------

.. versionadded:: 3.0

The next enumerations and functions are declared in ``primme_f90.inc``.

.. f:subroutine:: primme_svds_matvec(x, ldx, y, ldy, blockSize, mode, primme_svds, ierr)

   Abstract interface for the callbacks |SmatrixMatvec| and |SapplyPreconditioner|.

   :param type(*) x(ldx,*) [in]: matrix with ``blockSize`` columns in column-major_ order with leading dimension ``ldx``.
   :param c_int64_t ldx: the leading dimension of the array ``x``.
   :param type(*) y(ldy,*) [out]: matrix with ``blockSize`` columns in column-major_ order with leading dimension ``ldy``.
   :param c_int64_t ldy: the leading dimension of the array ``y``.
   :param c_int blockSize [in]: number of columns in ``x`` and ``y``.
   :param c_int mode [in]: a flag.
   :param c_ptr primme_svds [in]: parameters structure created by :f:func:`primme_svds_params_create`.
   :param c_int ierr [out]: output error code; if it is set to non-zero, the current call to PRIMME will stop.

   See more details about the precision and type and dimension for `x` and `y`, and the meaning of `mode` in the documentation of the callbacks.

primme_svds_params_create
"""""""""""""""""""""""""

.. f:function:: c_ptr primme_svds_params_create()

   Allocate and initialize a parameters structure to the default values.

   After calling :f:func:`xprimme_svds` (or a variant), call :f:func:`primme_svds_params_destroy()` to release allocated resources by PRIMME.

   :return c_ptr primme_svds_params_create: pointer to a parameters structure.


primme_svds_set_method
""""""""""""""""""""""

.. f:function:: c_int primme_svds_set_method(method, methodStage1, methodStage2, primme_svds)

   Set PRIMME SVDS parameters to one of the preset configurations.

   :param integer method: (input) preset configuration to compute the singular triplets; one of

      * |PRIMME_SVDS_default|, currently set as |PRIMME_SVDS_hybrid|.
      * |PRIMME_SVDS_normalequations|, compute the eigenvectors of :math:`A^*A` or :math:`A A^*`.
      * |PRIMME_SVDS_augmented|, compute the eigenvectors of the augmented matrix, :math:`\left(\begin{array}{cc} 0 & A^* \\ A & 0 \end{array}\right)`.
      * |PRIMME_SVDS_hybrid|, start with |PRIMME_SVDS_normalequations|; use the
        resulting approximate singular vectors as initial vectors for
        |PRIMME_SVDS_augmented| if the required accuracy was not achieved.

   :param primme_preset_method methodStage1: (input) preset method to compute the eigenpairs at the first stage; see available values at :f:func:`primme_set_method`.

   :param primme_preset_method methodStage2: (input) preset method to compute the eigenpairs with
      the second stage of |PRIMME_SVDS_hybrid|; see available values at :f:func:`primme_set_method`.

   :param ptr primme_svds: (input/output) parameters structure.

   :param integer ierr: (output) if 0, successful; if negative, something went wrong.

xprimme_svds
""""""""""""

.. f:function:: c_int xprimme_svds(svals, svecs, resNorms, primme_svds)

   Solve a real or complex singular value problem.

   All arrays should be hosted on CPU. The computations are performed on CPU (see :f:func:`magma_xprimme_svds` for using GPUs).

   :param svals(*) [out]: array at least of size |SnumSvals| to store the
      computed singular values; all processes in a parallel run return this local array with the same values.
   :type svals(*): real(kind) or complex(kind)

   :param svecs(*): array at least of size (|SmLocal| + |SnLocal|) times (|SnumOrthoConst| + |SnumSvals|)
      to store column-wise the (local part for this process of the) computed left singular vectors
      and the right singular vectors.
   :type svecs(*): real(kind) or complex(kind)

   :param resNorms(*) [out]: array at least of size |SnumSvals| to store the
      residual norms of the computed triplets; all processes in parallel run return this local array with
      the same values.
   :type resNorms(*): real(kind) or complex(kind)

   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_params_create_svds`.

   :return c_int xprimme_svds: error indicator; see :ref:`error-codes-svds`.

   The arrays ``svals``, ``svecs``, and ``resNorms`` should have the same kind.

   On input, ``svecs`` should start with the content of the |SnumOrthoConst| left vectors,
   followed by the |SinitSize| left vectors, followed by the |SnumOrthoConst| right vectors and
   followed by the |SinitSize| right vectors.
 
   On return, the i-th left singular vector starts at svecs(( |SnumOrthoConst| + i - 1) \* |SmLocal| ).
   The i-th right singular vector starts at svecs(( |SnumOrthoConst| + |SinitSize| )\* |SmLocal| + ( |SnumOrthoConst| + i - 1)\* |SnLocal| ).
   The first vector has i=1.
 
   All internal operations are performed at the same precision than ``svecs`` unless the user sets |SinternalPrecision| otherwise.

   The type and precision of the callbacks depends on the type and precision of ``svecs``. See details for |SmatrixMatvec|, |SapplyPreconditioner|, |SglobalSumReal|, |SbroadcastReal|, and |SconvTestFun|.

magma_xprimme_svds
""""""""""""""""""

.. f:function:: c_int magma_xprimme_svds(svals, svecs, resNorms, primme_svds)

   Solve a real or complex singular value problem.

   Most of the computations are performed on GPU (see :c:func:`xprimme_svds` for using only the CPU).

   :param svals(*) [out]: CPU array at least of size |SnumSvals| to store the
      computed singular values; all processes in a parallel run return this local array with the same values.
   :type svals(*): real(kind) or complex(kind)

   :param svecs(*): GPU array at least of size (|SmLocal| + |SnLocal|) times (|SnumOrthoConst| + |SnumSvals|)
      to store column-wise the (local part for this process of the) computed left singular vectors
      and the right singular vectors.
   :type svecs(*): real(kind) or complex(kind)

   :param resNorms(*) [out]: CPU array at least of size |SnumSvals| to store the
      residual norms of the computed triplets; all processes in parallel run return this local array with
      the same values.
   :type resNorms(*): real(kind) or complex(kind)

   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_params_create_svds`.

   :return c_int magma_xprimme_svds: error indicator; see :ref:`error-codes-svds`.

   The arrays ``svals``, ``svecs``, and ``resNorms`` should have the same kind.

   On input, ``svecs`` should start with the content of the |SnumOrthoConst| left vectors,
   followed by the |SinitSize| left vectors, followed by the |SnumOrthoConst| right vectors and
   followed by the |SinitSize| right vectors.
 
   On return, the i-th left singular vector starts at svecs(( |SnumOrthoConst| + i - 1) \* |SmLocal| ).
   The i-th right singular vector starts at svecs(( |SnumOrthoConst| + |SinitSize| )\* |SmLocal| + ( |SnumOrthoConst| + i - 1)\* |SnLocal| ).
   The first vector has i=1.
 
   All internal operations are performed at the same precision than ``svecs`` unless the user sets |SinternalPrecision| otherwise.

   The type and precision of the callbacks depends on the type and precision of ``svecs``. See details for |SmatrixMatvec|, |SapplyPreconditioner|, |SglobalSumReal|, |SbroadcastReal|, and |SconvTestFun|.

primme_svds_params_destroy
""""""""""""""""""""""""""

.. f:function:: c_int primme_svds_params_destroy(primme_svds)

   Free memory allocated by PRIMME associated to a parameters structure created
   with :f:func:`primme_svds_params_create`.

   :param c_ptr primme_svds: parameters structure.

   :return primme_svds_params_destroy: nonzero value if the call is not successful.
    
primme_svds_set_member
""""""""""""""""""""""

.. f:function:: c_int primme_svds_set_member(primme_svds, label, value)

   Set a value in some field of the parameter structure.

   :param ptr primme_svds: (input) parameters structure.

   :param integer label: field where to set value. One of:

     | :c:member:`PRIMME_SVDS_primme                         <primme_svds_params.primme>`
     | :c:member:`PRIMME_SVDS_primmeStage2                   <primme_svds_params.primmeStage2>`
     | :c:member:`PRIMME_SVDS_m                              <primme_svds_params.m>`
     | :c:member:`PRIMME_SVDS_n                              <primme_svds_params.n>`
     | :c:member:`PRIMME_SVDS_matrixMatvec                   <primme_svds_params.matrixMatvec>`
     | :c:member:`PRIMME_SVDS_matrixMatvec_type              <primme_svds_params.matrixMatvec_type>`
     | :c:member:`PRIMME_SVDS_applyPreconditioner            <primme_svds_params.applyPreconditioner>`
     | :c:member:`PRIMME_SVDS_applyPreconditioner_type       <primme_svds_params.applyPreconditioner_type>`
     | :c:member:`PRIMME_SVDS_numProcs                       <primme_svds_params.numProcs>`
     | :c:member:`PRIMME_SVDS_procID                         <primme_svds_params.procID>`
     | :c:member:`PRIMME_SVDS_mLocal                         <primme_svds_params.mLocal>`
     | :c:member:`PRIMME_SVDS_nLocal                         <primme_svds_params.nLocal>`
     | :c:member:`PRIMME_SVDS_commInfo                       <primme_svds_params.commInfo>`
     | :c:member:`PRIMME_SVDS_globalSumReal                  <primme_svds_params.globalSumReal>`
     | :c:member:`PRIMME_SVDS_globalSumReal_type             <primme_svds_params.globalSumReal_type>`
     | :c:member:`PRIMME_SVDS_broadcastReal                  <primme_svds_params.broadcastReal>`
     | :c:member:`PRIMME_SVDS_broadcastReal_type             <primme_svds_params.broadcastReal_type>`
     | :c:member:`PRIMME_SVDS_numSvals                       <primme_svds_params.numSvals>`
     | :c:member:`PRIMME_SVDS_target                         <primme_svds_params.target>`
     | :c:member:`PRIMME_SVDS_numTargetShifts                <primme_svds_params.numTargetShifts>`
     | :c:member:`PRIMME_SVDS_targetShifts                   <primme_svds_params.targetShifts>`
     | :c:member:`PRIMME_SVDS_method                         <primme_svds_params.method>`
     | :c:member:`PRIMME_SVDS_methodStage2                   <primme_svds_params.methodStage2>`
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
     | :c:member:`PRIMME_SVDS_internalPrecision              <primme_svds_params.internalPrecision>`
     | :c:member:`PRIMME_SVDS_convTestFun                    <primme_svds_params.convTestFun>`
     | :c:member:`PRIMME_SVDS_convTestFun_type               <primme_svds_params.convTestFun_type>`
     | :c:member:`PRIMME_SVDS_convtest                       <primme_svds_params.convtest>`
     | :c:member:`PRIMME_SVDS_monitorFun                     <primme_svds_params.monitorFun>`
     | :c:member:`PRIMME_SVDS_monitorFun_type                <primme_svds_params.monitorFun_type>`
     | :c:member:`PRIMME_SVDS_monitor                        <primme_svds_params.monitor>`
     | :c:member:`PRIMME_SVDS_queue                          <primme_svds_params.queue>`
     | :c:member:`PRIMME_SVDS_stats_numOuterIterations       <primme_svds_params.stats.numOuterIterations>`
     | :c:member:`PRIMME_SVDS_stats_numRestarts              <primme_svds_params.stats.numRestarts>`
     | :c:member:`PRIMME_SVDS_stats_numMatvecs               <primme_svds_params.stats.numMatvecs>`
     | :c:member:`PRIMME_SVDS_stats_numPreconds              <primme_svds_params.stats.numPreconds>`
     | :c:member:`PRIMME_SVDS_stats_numGlobalSum             <primme_svds_params.stats.numGlobalSum>`
     | :c:member:`PRIMME_SVDS_stats_numBroadcast             <primme_svds_params.stats.numBroadcast>`
     | :c:member:`PRIMME_SVDS_stats_volumeGlobalSum          <primme_svds_params.stats.volumeGlobalSum>`
     | :c:member:`PRIMME_SVDS_stats_volumeBroadcast          <primme_svds_params.stats.volumeBroadcast>`
     | :c:member:`PRIMME_SVDS_stats_elapsedTime              <primme_svds_params.stats.elapsedTime>`
     | :c:member:`PRIMME_SVDS_stats_timeMatvec               <primme_svds_params.stats.timeMatvec>`
     | :c:member:`PRIMME_SVDS_stats_timePrecond              <primme_svds_params.stats.timePrecond>`
     | :c:member:`PRIMME_SVDS_stats_timeOrtho                <primme_svds_params.stats.timeOrtho>`
     | :c:member:`PRIMME_SVDS_stats_timeGlobalSum            <primme_svds_params.stats.timeGlobalSum>`
     | :c:member:`PRIMME_SVDS_stats_timeBroadcast            <primme_svds_params.stats.timeBroadcast>`
     | :c:member:`PRIMME_SVDS_stats_lockingIssue             <primme_svds_params.stats.lockingIssue>`

   :param value: (input) value to set.
     The allowed types are `c_int64`, `c_double`, `c_ptr`, `c_funptr` and :f:func:`procedure(primme_svds_matvec) <primme_svds_matvec>`

   :return c_int primme_svds_set_member: nonzero value if the call is not successful.

    Examples::

        type(c_ptr) :: primme_svds
        integer :: ierr
        ...
        
        integer(c_int64_t) :: m               = 100
        ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_m, m)
        ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_n, m)

        real(c_double) :: tol             = 1.0D-12
        ierr = primme_svds_set_member(primme, PRIMME_SVDS_eps, tol)

        integer(c_int64_t), parameter :: numTargetShifts = 2
        real(c_double) :: TargetShifts(numTargetShifts) = (/3.0D0, 5.1D0/)
        ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_numTargetShifts, numTargetShifts)
        ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_targetShifts, TargetShifts)

        ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_target, primme_svds_closest_abs)
        
        procedure(primme_svds_matvec) :: MV, ApplyPrecon
        ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_matrixMatvec, MV)
        
        ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_applyPreconditioner,
                                 c_funloc(ApplyPrecon))
        
        type(c_ptr) :: primme
        ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_primme, primme)
        ierr = primme_set_member(primme, PRIMME_correctionParams_precondition,
                                 1_c_int64_t)

primme_get_member
"""""""""""""""""

.. f:function:: c_int primme_svds_get_member(primme, label, value)

   Get the value in some field of the parameter structure.

   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_svds_params_create`.

   :param integer label [in]: field where to get value. One of the detailed in function :f:func:`primme_svds_set_member`.

   :param value [out]: value of the field.  The allowed types are `c_int64`, `c_double`, and `c_ptr`.

   :return c_int primme_svds_get_member: nonzero value if the call is not successful.

   Examples::

        type(c_ptr) :: primme_svds
        integer :: ierr
        ...

        integer(c_int64_t) :: m
        ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_m, m)

        real(c_double) :: aNorm
        ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_aNorm, aNorm)

        type(c_ptr) :: primme
        ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_primme, primme)
        ierr = primme_set_member(primme, PRIMME_correctionParams_precondition,
                                 1_c_int64_t)


.. include:: epilog.inc
