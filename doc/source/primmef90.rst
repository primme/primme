
.. highlight:: Fortran

.. f:currentmodule::

FORTRAN 90 Library Interface
----------------------------

.. versionadded:: 3.0

The next enumerations and functions are declared in ``primme_f90.inc``.

.. f:type:: iso_c_binding/c_ptr

   Fortran datatype for C pointers.
    
.. f:subroutine:: primme_eigs_matvec(x, ldx, y, ldy, blockSize, primme, ierr)

   Abstract interface for the callbacks |matrixMatvec|, |massMatrixMatvec|, and |applyPreconditioner|.

   :param type(*) x(ldx,*) [in]: matrix of size |nLocal| x ``blockSize`` in column-major_ order with leading dimension ``ldx``.
   :param c_int64_t ldx: the leading dimension of the array ``x``.
   :param type(*) y(ldy,*) [out]: matrix of size |nLocal| x ``blockSize`` in column-major_ order with leading dimension ``ldy``.
   :param c_int64_t ldy: the leading dimension of the array ``y``.
   :param c_int blockSize [in]: number of columns in ``x`` and ``y``.
   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_params_create`.
   :param c_int ierr [out]: output error code; if it is set to non-zero, the current call to PRIMME will stop.

   See more details about the precision and type for `x` and `y` in the documentation of the callbacks.

primme_params_create
""""""""""""""""""""

.. f:function:: c_ptr primme_params_create()

   Allocate and initialize a parameters structure to the default values.

   After calling :f:func:`xprimme` (or a variant), call :f:func:`primme_params_destroy()` to release allocated resources by PRIMME.

   :return c_ptr primme_params_create: pointer to a parameters structure.
    
primme_set_method
"""""""""""""""""

.. f:function:: c_int primme_set_method(method, primme)

   Set PRIMME parameters to one of the preset configurations.

   :param c_int method [in]: preset configuration. One of:

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

   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_params_create`.

   :return c_int primme_set_method: nonzero value if the call is not successful.
    
primme_params_destroy
"""""""""""""""""""""

.. f:function:: c_int primme_params_destroy(primme)

   Free memory allocated by PRIMME associated to a parameters structure created
   with :f:func:`primme_params_create`.

   :param c_ptr primme: parameters structure.

   :return primme_params_destroy: nonzero value if the call is not successful.
    
xprimme
"""""""

.. f:function:: c_int xprimme(evals, evecs, resNorms, primme)

   Solve a real symmetric/Hermitian standard or generalized eigenproblem.

   All arrays should be hosted on CPU. The computations are performed on CPU (see :f:func:`magma_xprimme` for using GPUs).

   :param real(kind) evals(*) [out]: array at least of size |numEvals| to store the
      computed eigenvalues; all processes in a parallel run return this local array with the same values.

   :param evecs: array at least of size |nLocal| times (|numOrthoConst| + |numEvals|) with leading dimension |ldevecs|
      to store column-wise the (local part for this process of the) computed eigenvectors.
   :rtype: real(kind) or complex(kind)

   :param real(kind) resNorms(*) [out]: array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all processes in parallel run return this local array with
      the same values.

   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_params_create`.

   :return c_int xprimme: error indicator; see :ref:`error-codes`.

   The arrays ``evals``, ``evecs``, and ``resNorms`` should have the same kind.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |ldevecs| + 1).
   The first vector has index i=1.
 
   All internal operations are performed at the same precision than ``evecs`` unless the user sets |internalPrecision| otherwise.

   The type and precision of the callbacks is also the same as ``evecs``. Although this can be changed. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

magma_xprimme
"""""""""""""

.. f:function:: c_int magma_xprimme(evals, evecs, resNorms, primme)

   Solve a real symmetric/Hermitian standard or generalized eigenproblem.

   Most of the computations are performed on GPU (see :f:func:`xprimme` for using only the CPU).

   :param real(kind) evals(*) [out]: CPU array at least of size |numEvals| to store the
      computed eigenvalues; all processes in a parallel run return this local array with the same values.

   :param evecs: GPU array at least of size |nLocal| times (|numOrthoConst| + |numEvals|) with leading dimension |ldevecs|
      to store column-wise the (local part for this process of the) computed eigenvectors.
   :rtype: real(kind) or complex(kind)

   :param real(kind) resNorms(*) [out]: CPU array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all processes in parallel run return this local array with
      the same values.

   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_params_create`.

   :return c_int magma_xprimme: error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |ldevecs| + 1).
   The first vector has index i=1.

   The arrays ``evals``, ``evecs``, and ``resNorms`` should have the same kind.

   All internal operations are performed at the same precision than ``evecs`` unless the user sets |internalPrecision| otherwise.

   The type and precision of the callbacks is also the same as ``evecs``. Although this can be changed. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.
 
xprimme_normal
""""""""""""""

.. f:function:: c_int xprimme_normal(evals, evecs, resNorms, primme)

   Solve a normal standard eigenproblem, which may not be Hermitian.

   All arrays should be hosted on CPU. The computations are performed on CPU (see :f:func:`magma_xprimme_normal` for using GPUs).

   :param complex(kind) evals(*) [out]: array at least of size |numEvals| to store the
      computed eigenvalues; all processes in a parallel run return this local array with the same values.

   :param evecs: array at least of size |nLocal| times (|numOrthoConst| + |numEvals|) with leading dimension |ldevecs|
      to store column-wise the (local part for this process of the) computed eigenvectors.
   :rtype: complex(kind)

   :param real(kind) resNorms(*) [out]: array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all processes in parallel run return this local array with
      the same values.

   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_params_create`.

   :return c_int xprimme_normal: error indicator; see :ref:`error-codes`.

   The arrays ``evals``, ``evecs``, and ``resNorms`` should have the same kind.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |ldevecs| + 1).
   The first vector has index i=1.
 
   All internal operations are performed at the same precision than ``evecs`` unless the user sets |internalPrecision| otherwise.

   The type and precision of the callbacks is also the same as ``evecs``. Although this can be changed. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

magma_xprimme_normal
""""""""""""""""""""

.. f:function:: c_int magma_xprimme_normal(evals, evecs, resNorms, primme)

   Solve a normal standard eigenproblem, which may not be Hermitian.

   Most of the arrays are stored on GPU, and also most of the computations are done on GPU (see :f:func:`xprimme` for using only the CPU).

   :param complex(kind) evals(*) [out]: CPU array at least of size |numEvals| to store the
      computed eigenvalues; all processes in a parallel run return this local array with the same values.

   :param evecs: GPU array at least of size |nLocal| times (|numOrthoConst| + |numEvals|) with leading dimension |ldevecs|
      to store column-wise the (local part for this process of the) computed eigenvectors.
   :rtype: complex(kind)

   :param real(kind) resNorms(*) [out]: CPU array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all processes in parallel run return this local array with
      the same values.

   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_params_create`.

   :return c_int magma_xprimme_normal: error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |ldevecs| + 1).
   The first vector has index i=1.

   The arrays ``evals``, ``evecs``, and ``resNorms`` should have the same kind.

   All internal operations are performed at the same precision than ``evecs`` unless the user sets |internalPrecision| otherwise.

   The type and precision of the callbacks is also the same as ``evecs``. Although this can be changed. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.
 
primme_set_member
"""""""""""""""""

.. f:function:: c_int primme_set_member(primme, label, value)

   Set a value in some field of the parameter structure.

   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_params_create`.

   :param c_int label [in]: field where to set value. One of:

      | :c:member:`PRIMME_n                                   <primme_params.n>`
      | :c:member:`PRIMME_matrixMatvec                        <primme_params.matrixMatvec>`
      | :c:member:`PRIMME_matrixMatvec_type                   <primme_params.matrixMatvec_type>`
      | :c:member:`PRIMME_applyPreconditioner                 <primme_params.applyPreconditioner>`
      | :c:member:`PRIMME_applyPreconditioner_type            <primme_params.applyPreconditioner_type>`
      | :c:member:`PRIMME_massMatrixMatvec                    <primme_params.massMatrixMatvec>`
      | :c:member:`PRIMME_massMatrixMatvec_type               <primme_params.massMatrixMatvec_type>`
      | :c:member:`PRIMME_numProcs                            <primme_params.numProcs>`
      | :c:member:`PRIMME_procID                              <primme_params.procID>`
      | :c:member:`PRIMME_commInfo                            <primme_params.commInfo>`
      | :c:member:`PRIMME_nLocal                              <primme_params.nLocal>`
      | :c:member:`PRIMME_globalSumReal                       <primme_params.globalSumReal>`
      | :c:member:`PRIMME_globalSumReal_type                  <primme_params.globalSumReal_type>`
      | :c:member:`PRIMME_broadcastReal                       <primme_params.broadcastReal>`
      | :c:member:`PRIMME_broadcastReal_type                  <primme_params.broadcastReal_type>`
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
      | :c:member:`PRIMME_iseed                               <primme_params.iseed>`
      | :c:member:`PRIMME_aNorm                               <primme_params.aNorm>`
      | :c:member:`PRIMME_BNorm                               <primme_params.BNorm>`
      | :c:member:`PRIMME_invBNorm                            <primme_params.invBNorm>`
      | :c:member:`PRIMME_eps                                 <primme_params.eps>`
      | :c:member:`PRIMME_orth                                <primme_params.orth>`
      | :c:member:`PRIMME_internalPrecision                   <primme_params.internalPrecision>`
      | :c:member:`PRIMME_printLevel                          <primme_params.printLevel>`
      | :c:member:`PRIMME_outputFile                          <primme_params.outputFile>`
      | :c:member:`PRIMME_matrix                              <primme_params.matrix>`
      | :c:member:`PRIMME_massMatrix                          <primme_params.massMatrix>`
      | :c:member:`PRIMME_preconditioner                      <primme_params.preconditioner>`
      | :c:member:`PRIMME_initBasisMode                       <primme_params.initBasisMode>`
      | :c:member:`PRIMME_projectionParams_projection         <primme_params.projectionParams.projection>`
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
      | :c:member:`PRIMME_stats_numGlobalSum                  <primme_params.stats.numGlobalSum>`
      | :c:member:`PRIMME_stats_volumeGlobalSum               <primme_params.stats.volumeGlobalSum>`
      | :c:member:`PRIMME_stats_numBroadcast                  <primme_params.stats.numBroadcast>`
      | :c:member:`PRIMME_stats_volumeBroadcast               <primme_params.stats.volumeBroadcast>`
      | :c:member:`PRIMME_stats_flopsDense                    <primme_params.stats.flopsDense>`
      | :c:member:`PRIMME_stats_numOrthoInnerProds            <primme_params.stats.numOrthoInnerProds>`
      | :c:member:`PRIMME_stats_elapsedTime                   <primme_params.stats.elapsedTime>`
      | :c:member:`PRIMME_stats_timeMatvec                    <primme_params.stats.timeMatvec>`
      | :c:member:`PRIMME_stats_timePrecond                   <primme_params.stats.timePrecond>`
      | :c:member:`PRIMME_stats_timeOrtho                     <primme_params.stats.timeOrtho>`
      | :c:member:`PRIMME_stats_timeGlobalSum                 <primme_params.stats.timeGlobalSum>`
      | :c:member:`PRIMME_stats_timeBroadcast                 <primme_params.stats.timeBroadcast>`
      | :c:member:`PRIMME_stats_timeDense                     <primme_params.stats.timeDense>`
      | :c:member:`PRIMME_stats_estimateMinEVal               <primme_params.stats.estimateMinEVal>`
      | :c:member:`PRIMME_stats_estimateMaxEVal               <primme_params.stats.estimateaxnEVal>`
      | :c:member:`PRIMME_stats_estimateLargestSVal           <primme_params.stats.estimateLargestSVal>`
      | :c:member:`PRIMME_stats_estimateBNorm                 <primme_params.stats.estimateBNorm>`
      | :c:member:`PRIMME_stats_estimateInvBNorm              <primme_params.stats.estimateInvBNorm>`
      | :c:member:`PRIMME_stats_maxConvTol                    <primme_params.stats.maxConvTol>`
      | :c:member:`PRIMME_stats_lockingIssue                  <primme_params.stats.lockingIssue>`
      | :c:member:`PRIMME_dynamicMethodSwitch                 <primme_params.dynamicMethodSwitch>`
      | :c:member:`PRIMME_convTestFun                         <primme_params.convTestFun>`
      | :c:member:`PRIMME_convTestFun_type                    <primme_params.convTestFun_type>`
      | :c:member:`PRIMME_convtest                            <primme_params.convtest>`
      | :c:member:`PRIMME_ldevecs                             <primme_params.ldevecs>`
      | :c:member:`PRIMME_ldOPs                               <primme_params.ldOPs>`
      | :c:member:`PRIMME_monitorFun                          <primme_params.monitorFun>`
      | :c:member:`PRIMME_monitorFun_type                     <primme_params.monitorFun_type>`
      | :c:member:`PRIMME_monitor                             <primme_params.monitor>`
      | :c:member:`PRIMME_queue                               <primme_params.queue>`
 
    
   :param value [in]: value to set. The allowed types are `c_int64`, `c_double`, `c_ptr`, `c_funptr` and :f:func:`procedure(primme_eigs_matvec) <primme_eigs_matvec>`

   :return c_int primme_set_member: nonzero value if the call is not successful.

    Examples::

        type(c_ptr) :: primme
        integer :: ierr
        ...
        
        integer(c_int64_t) :: n               = 100
        ierr = primme_set_member(primme, PRIMME_n, n)

        ierr = primme_set_member(primme, PRIMME_correctionParams_precondition,
                                 1_c_int64_t)

        real(c_double) :: tol             = 1.0D-12
        ierr = primme_set_member(primme, PRIMME_eps, tol)

        integer(c_int64_t), parameter :: numTargetShifts = 2
        real(c_double) :: TargetShifts(numTargetShifts) = (/3.0D0, 5.1D0/)
        ierr = primme_set_member(primme, PRIMME_numTargetShifts, numTargetShifts)
        ierr = primme_set_member(primme, PRIMME_targetShifts, TargetShifts)

        ierr = primme_set_member(primme, PRIMME_target, primme_closest_abs)
        
        procedure(primme_eigs_matvec) :: MV, ApplyPrecon
        ierr = primme_set_member(primme, PRIMME_matrixMatvec, MV)
        
        ierr = primme_set_member(primme, PRIMME_applyPreconditioner,
                                 c_funloc(ApplyPrecon))

 
primme_get_member
"""""""""""""""""

.. f:function:: c_int primme_get_member(primme, label, value)

   Get the value in some field of the parameter structure.

   :param c_ptr primme [in]: parameters structure created by :f:func:`primme_params_create`.

   :param integer label [in]: field where to get value. One of the detailed in function :f:func:`primme_set_member`.

   :param value [out]: value of the field.  The allowed types are `c_int64`, `c_double`, and `c_ptr`.

   :return c_int primme_get_member: nonzero value if the call is not successful.

   Examples::

        type(c_ptr) :: primme
        integer :: ierr
        ...

        integer(c_int64_t) :: n
        ierr = primme_get_member(primme, PRIMME_n, n)

        real(c_double) :: aNorm
        ierr = primme_get_member(primme, PRIMME_aNorm, aNorm)

        real(c_double), pointer :: shifts(:)
        type(c_ptr) :: pshifts
        ierr = primme_get_member(primme, PRIMME_ShiftsForPreconditioner, pshifts)
        call c_f_pointer(pshifts, shifts, shape=[k])


.. include:: epilog.inc
