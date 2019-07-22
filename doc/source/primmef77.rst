
.. highlight:: fortran

FORTRAN 77 Library Interface
----------------------------

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

   All arrays are stored on CPU, and also the computations are done on CPU (see :c:func:`magma_sprimme_f77` for using GPUs).
   
   :param evals(*): (output) array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): real

   :param evecs(*): (input/output) array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): real

   :param resNorms(*): (output) array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): real

   :param ptr primme: parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

cprimme_f77
"""""""""""

.. c:function:: cprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblem. The arguments have the
   same meaning as in function :c:func:`sprimme_f77`.

   All arrays are stored on CPU, and also the computations are done on CPU (see :c:func:`magma_cprimme_f77` for using GPUs).

   :param evals(*): (output) array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): real

   :param evecs(*): (input/output) array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): complex real

   :param resNorms(*): (output) array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): real

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

dprimme_f77
"""""""""""

.. c:function:: dprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a real symmetric standard eigenproblem using double precision.

   All arrays are stored on CPU, and also the computations are done on CPU (see :c:func:`magma_dprimme_f77` for using GPUs).

   :param evals(*): (output) array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): double precision

   :param evecs(*): (input/output) array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): double precision

   :param resNorms(*): (output) array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): double precision

   :param ptr primme: parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

zprimme_f77
"""""""""""

.. c:function:: zprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblem. The arguments have the
   same meaning as in function :c:func:`dprimme_f77`.

   All arrays are stored on CPU, and also the computations are done on CPU (see :c:func:`magma_zprimme_f77` for using GPUs).

   :param evals(*): (output) array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): double precision

   :param evecs(*): (input/output) array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): complex double precision

   :param resNorms(*): (output) array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): double precision

   :param ptr primme: parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

magma_sprimme_f77
"""""""""""""""""

.. c:function:: magma_sprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a real symmetric standard eigenproblem using single precision.

   Most of the arrays are stored on GPU, and also most of the computations are done on GPU (see :c:func:`sprimme_f77` for using only the CPU).
   
   :param evals(*): (output) CPU array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): real

   :param evecs(*): (input/output) GPU array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): real

   :param resNorms(*): (output) CPU array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): real

   :param ptr primme: parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

magma_cprimme_f77
"""""""""""""""""

.. c:function:: magma_cprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblem. The arguments have the
   same meaning as in function :c:func:`sprimme_f77`.

   Most of the arrays are stored on GPU, and also most of the computations are done on GPU (see :c:func:`cprimme_f77` for using only the CPU).

   :param evals(*): (output) CPU array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): real

   :param evecs(*): (input/output) GPU array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): complex real

   :param resNorms(*): (output) CPU array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): real

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

magma_dprimme_f77
"""""""""""""""""

.. c:function:: magma_dprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a real symmetric standard eigenproblem using double precision.

   Most of the arrays are stored on GPU, and also most of the computations are done on GPU (see :c:func:`dprimme_f77` for using only the CPU).

   :param evals(*): (output) CPU array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): double precision

   :param evecs(*): (input/output) GPU array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): double precision

   :param resNorms(*): (output) CPU array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): double precision

   :param ptr primme: parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

magma_zprimme_f77
"""""""""""""""""

.. c:function:: magma_zprimme_f77(evals, evecs, resNorms, primme, ierr)

   Solve a Hermitian standard eigenproblem. The arguments have the
   same meaning as in function :c:func:`dprimme_f77`.

   Most of the arrays are stored on GPU, and also most of the computations are done on GPU (see :c:func:`zprimme_f77` for using only the CPU).

   :param evals(*): (output) CPU array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): double precision

   :param evecs(*): (input/output) GPU array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): complex double precision

   :param resNorms(*): (output) CPU array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): double precision

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.


cprimme_normal_f77
""""""""""""""""""

.. c:function:: cprimme_normal_f77(evals, evecs, resNorms, primme, ierr)

   Solve a normal standard eigenproblem, which may not be Hermitian.

   All arrays are stored on CPU, and also the computations are done on CPU (see :c:func:`magma_cprimme_normal_f77` for using GPUs).

   :param evals(*): (output) array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): real

   :param evecs(*): (input/output) array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): complex real

   :param resNorms(*): (output) array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): real

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

zprimme_normal_f77
""""""""""""""""""

.. c:function:: zprimme_normal_f77(evals, evecs, resNorms, primme, ierr)

   Solve a normal standard eigenproblem, which may not be Hermitian.

   All arrays are stored on CPU, and also the computations are done on CPU (see :c:func:`magma_zprimme_normal_f77` for using GPUs).

   :param evals(*): (output) array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): double precision

   :param evecs(*): (input/output) array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): complex double precision

   :param resNorms(*): (output) array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): double precision

   :param ptr primme: parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

magma_cprimme_normal_f77
""""""""""""""""""""""""

.. c:function:: magma_cprimme_normal_f77(evals, evecs, resNorms, primme, ierr)

   Solve a normal standard eigenproblem, which may not be Hermitian.

   Most of the arrays are stored on GPU, and also most of the computations are done on GPU (see :c:func:`cprimme_normal_f77` for using only the CPU).

   :param evals(*): (output) CPU array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): real

   :param evecs(*): (input/output) GPU array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): complex real

   :param resNorms(*): (output) CPU array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): real

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

magma_zprimme_normal_f77
""""""""""""""""""""""""

.. c:function:: magma_zprimme_normal_f77(evals, evecs, resNorms, primme, ierr)

   Solve a normal standard eigenproblem, which may not be Hermitian.

   Most of the arrays are stored on GPU, and also most of the computations are done on GPU (see :c:func:`zprimme_normal_f77` for using only the CPU).

   :param evals(*): (output) CPU array at least of size |numEvals| to store the
      computed eigenvalues; all parallel calls return the same value in this array.
   :type evals(*): double precision

   :param evecs(*): (input/output) GPU array at least of size |nLocal| times (|numOrthoConst| + |numEvals|)
      to store columnwise the (local part of the) computed eigenvectors.
   :type evecs(*): complex double precision

   :param resNorms(*): (output) CPU array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all parallel calls return the same value in this array.
   :type resNorms(*): double precision

   :param ptr primme: (input) parameters structure.

   :param integer ierr: (output) error indicator; see :ref:`error-codes`.

   On input, ``evecs`` should start with the content of the |numOrthoConst| vectors,
   followed by the |initSize| vectors.
 
   On return, the i-th eigenvector starts at evecs(( |numOrthoConst| + i - 1)\* |SnLocal| + 1).
   The first vector has index i=1.
 
   The type and precision of the callbacks depends on the type and precision of `evecs`. See details for |matrixMatvec|, |massMatrixMatvec|, |applyPreconditioner|, |globalSumReal|, |broadcastReal|, and |convTestFun|.

primme_set_member_f77
"""""""""""""""""""""

.. c:function:: primme_set_member_f77(primme, label, value)

   Set a value in some field of the parameter structure.

   :param ptr primme: (input) parameters structure.

   :param integer label: field where to set value. One of:

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
      | :c:member:`PRIMME_projectionParams_projection         <primme_params.projectionParams_projection>`
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
      ``PRIMME_matrix`` and ``PRIMME_preconditioner``,
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
      ``PRIMME_matrix`` and ``PRIMME_preconditioner``,
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
