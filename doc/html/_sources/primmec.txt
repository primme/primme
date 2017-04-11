
.. role:: ccode(code) 
   :language: c

.. highlight:: c

C Library Interface
-------------------

The PRIMME interface is composed of the following functions.
To solve real symmetric and Hermitian standard eigenproblems call
respectively:

.. only:: not text

   .. parsed-literal::

      int :c:func:`sprimme <sprimme>` (float \*evals, float \*evecs, float \*resNorms,
                              primme_params \*primme)
      int :c:func:`cprimme <cprimme>` (float \*evals, PRIMME_COMPLEX_FLOAT \*evecs,
                       float \*resNorms, primme_params \*primme)
      int :c:func:`dprimme <dprimme>` (double \*evals, double \*evecs, double \*resNorms,
                              primme_params \*primme)
      int :c:func:`zprimme <zprimme>` (double \*evals, PRIMME_COMPLEX_DOUBLE \*evecs,
                       double \*resNorms, primme_params \*primme)

.. only:: text

   ::

      int sprimme(float *evals, float *evecs, float *resNorms, 
                  primme_params *primme);

      int cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, 
                  primme_params *primme);

      int dprimme(double *evals, double *evecs, double *resNorms, 
                  primme_params *primme);

      int zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, 
                  primme_params *primme);

Other useful functions:

.. only:: not text

   .. parsed-literal::

      void :c:func:`primme_initialize <primme_initialize>` (primme_params \*primme)
      int :c:func:`primme_set_method <primme_set_method>` (primme_preset_method method,
                                                           primme_params \*params)
      void :c:func:`primme_display_params <primme_display_params>` (primme_params primme)
      void :c:func:`primme_free <primme_Free>` (primme_params \*primme)

.. only:: text

   ::

      void primme_initialize(primme_params *primme);
      int primme_set_method(primme_preset_method method,
                                           primme_params *params);
      void primme_display_params(primme_params primme);
      void primme_free(primme_params primme);

PRIMME stores its data on the structure :c:type:`primme_params`.
See :ref:`guide-params` for an introduction about its fields.


Running
^^^^^^^

To use PRIMME, follow these basic steps.

#. Include::

      #include "primme.h"   /* header file is required to run primme */

#. Initialize a PRIMME parameters structure for default settings:

   .. only:: not text
   
       .. parsed-literal::

          :c:type:`primme_params` primme;
          :c:func:`primme_initialize <primme_initialize>` (&primme);

   .. only:: text
   
      ::
   
         primme_params primme;
         
         primme_initialize(&primme);
   
#. Set problem parameters (see also :ref:`guide-params`), and,
   optionally, set one of the :c:type:`preset methods <primme_preset_method>`:

   .. only:: not text

      .. parsed-literal::

         primme.\ |matrixMatvec| = LaplacianMatrixMatvec; /\* MV product \*/
         primme.\ |n| = 100;                   /\* set problem dimension \*/
         primme.\ |numEvals| = 10;       /\* Number of wanted eigenpairs \*/
         ret = :c:func:`primme_set_method <primme_set_method>` (method, &primme);
         ...

   .. only:: text

      ::

         primme.matrixMatvec = LaplacianMatrixMatvec; /* MV product */
         primme.n = 100;                   /* set problem dimension */
         primme.numEvals = 10;       /* Number of wanted eigenpairs */
         ret = primme_set_method(method, &primme);
         ...

#. Then to solve real symmetric standard eigenproblems call:

   .. only:: not text
  
      .. parsed-literal::
 
         ret = :c:func:`dprimme <dprimme>` (evals, evecs, resNorms, &primme);
   
   .. only:: text
   
      ::
   
         ret = dprimme(evals, evecs, resNorms, &primme);

   The previous is the double precision call. There is available calls for complex
   double, single and complex single; check it out :c:func:`zprimme`, :c:func:`sprimme`
   and :c:func:`cprimme`.

   The call arguments are:

   * `evals`, array to return the found eigenvalues;
   * `evecs`, array to return the found eigenvectors;
   * `resNorms`, array to return the residual norms of the found eigenpairs; and
   * `ret`, returned error code.

#. To free the work arrays in PRIMME:

   .. only:: not text
  
      .. parsed-literal::
 
         :c:func:`primme_free <primme_Free>` (&primme);
   
   .. only:: text
   
      ::
   
         primme_free(&primme);

.. _guide-params:

Parameters Guide
^^^^^^^^^^^^^^^^

PRIMME stores the data on the structure :c:type:`primme_params`, which has the next fields:
   
.. only:: not text

      | *Basic*
      | ``PRIMME_INT`` |n|,  matrix dimension.
      | ``void (*`` |matrixMatvec| ``)(...)``, matrix-vector product.
      | ``int`` |numEvals|, how many eigenpairs to find.
      | ``primme_target`` |target|, which eigenvalues to find.
      | ``int`` |numTargetShifts|, for targeting interior eigenpairs.
      | ``double *`` |targetShifts|
      | ``double`` |eps|, tolerance of the residual norm of converged eigenpairs.
      |
      | *For parallel programs*
      | ``int`` |numProcs|, number of processes
      | ``int`` |procID|,  rank of this process
      | ``PRIMME_INT`` |nLocal|,  number of rows stored in this process
      | ``void (*`` |globalSumReal| ``)(...)``, sum reduction among processes
      |
      | *Accelerate the convergence*
      | ``void (*`` |applyPreconditioner| ``)(...)``, preconditioner-vector product.
      | ``int`` |initSize|, initial vectors as approximate solutions.
      | ``int`` |maxBasisSize|
      | ``int`` |minRestartSize|
      | ``int`` |maxBlockSize|
      |
      | *User data*
      | ``void *`` |commInfo|
      | ``void *`` |matrix|
      | ``void *`` |preconditioner|
      | ``void *`` convTest
      | ``void *`` monitor
      |
      | *Advanced options*
      | ``PRIMME_INT`` |ldevecs|, leading dimension of the evecs.
      | ``int`` |numOrthoConst|, orthogonal constrains to the eigenvectors.
      | ``int`` |dynamicMethodSwitch|
      | ``int`` |locking|
      | ``PRIMME_INT`` |maxMatvecs|
      | ``PRIMME_INT`` |maxOuterIterations|
      | ``int`` |intWorkSize|
      | ``size_t`` |realWorkSize|
      | ``PRIMME_INT`` |iseed| ``[4]``
      | ``int *`` |intWork|
      | ``void *`` |realWork|
      | ``double`` |aNorm|
      | ``int`` |printLevel|
      | ``FILE *`` |outputFile|
      | ``double *`` |ShiftsForPreconditioner|
      | ``primme_init`` |initBasisMode|
      | ``struct projection_params`` :c:member:`projectionParams <primme_params.projectionParams.projection>`
      | ``struct restarting_params`` :c:member:`restartingParams <primme_params.restartingParams.scheme>`
      | ``struct correction_params`` :c:member:`correctionParams <primme_params.correctionParams.precondition>`
      | ``struct primme_stats`` :c:member:`stats <primme_params.stats.numOuterIterations>`
      | ``void (*`` |convTestFun| ``)(...)``, custom convergence criterion.
      | ``PRIMME_INT`` |ldOPS|, leading dimension to use in |matrixMatvec|.
      | ``void (*`` |monitorFun| ``)(...)``, custom convergence history.

.. only:: text

   ::

      /* Basic */
      PRIMME_INT n;                                      // matrix dimension
      void (*matrixMatvec)(...);             // matrix-vector product
      int numEvals;                    // how many eigenpairs to find
      primme_target target;              // which eigenvalues to find
      int numTargetShifts;       // for targeting interior eigenpairs
      double *targetShifts;
      double eps;            // tolerance of the converged eigenpairs
      
      /* For parallel programs */
      int numProcs;           // number of processes
      int procID;             // rank of this process 
      PRIMME_INT nLocal;      // number of rows stored in this process
      void (*globalSumReal)(...); // sum reduction among processes
      
      /* Accelerate the convergence */
      void (*applyPreconditioner)(...);     // precond-vector product
      int initSize;       // initial vectors as approximate solutions
      int maxBasisSize;
      int minRestartSize;
      int maxBlockSize;
      
      /* User data */
      void *commInfo;
      void *matrix;
      void *preconditioner;
      void *convTest;
      void *monitor;
      
      /* Advanced options */
      PRIMME_INT ldevecs; // leading dimension of the evecs
      int numOrthoConst; // orthogonal constrains to the eigenvectors
      int dynamicMethodSwitch;
      int locking;
      PRIMME_INT maxMatvecs;
      PRIMME_INT maxOuterIterations;
      int intWorkSize;
      size_t realWorkSize;
      PRIMME_INT iseed[4];
      int *intWork;
      void *realWork;
      double aNorm;
      int printLevel;
      FILE *outputFile;
      double *ShiftsForPreconditioner;
      primme_init initBasisMode;
      struct projection_params projectionParams;
      struct restarting_params restartingParams;
      struct correction_params correctionParams;
      struct primme_stats stats;
      void (*convTestFun)(...); // custom convergence criterion
      PRIMME_INT ldOPS;   // leading dimension to use in matrixMatvec
      void (*monitorFun)(...); // custom convergence history
 
PRIMME requires the user to set at least the dimension of the matrix (|n|) and
the matrix-vector product (|matrixMatvec|), as they define the problem to be solved.
For parallel programs, |nLocal|, |procID| and |globalSumReal| are also required.

In addition, most users would want to specify how many eigenpairs to find,
and provide a preconditioner (if available).

It is useful to have set all these before calling :c:func:`primme_set_method`.
Also, if users have a preference on |maxBasisSize|, |maxBlockSize|, etc, they
should also provide them into :c:type:`primme_params` prior to the
:c:func:`primme_set_method` call. This helps :c:func:`primme_set_method` make
the right choice on other parameters. It is sometimes useful to check the actual
parameters that PRIMME is going to use (before calling it) or used (on return)
by printing them with :c:func:`primme_display_params`.

Interface Description
^^^^^^^^^^^^^^^^^^^^^

The next enumerations and functions are declared in ``primme.h``.

sprimme
"""""""

.. c:function:: int sprimme (float *evals, float *evecs, float *resNorms, primme_params *primme)

   Solve a real symmetric standard eigenproblem.

   :param evals: array at least of size |numEvals| to store the
      computed eigenvalues; all processes in a parallel run return this local array with the same values.

   :param resNorms: array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all processes in parallel run return this local array with
      the same values.

   :param evecs: array at least of size |nLocal| times |numEvals|
      to store columnwise the (local part of the) computed eigenvectors.

   :param primme: parameters structure.

   :return: error indicator; see :ref:`error-codes`.

dprimme
"""""""

.. c:function:: int dprimme(double *evals, double *evecs, double *resNorms, primme_params *primme)

   Solve a real symmetric standard eigenproblem.

   :param evals: array at least of size |numEvals| to store the
      computed eigenvalues; all processes in a parallel run return this local array with the same values.

   :param resNorms: array at least of size |numEvals| to store the
      residual norms of the computed eigenpairs; all processes in parallel run return this local array with
      the same values.

   :param evecs: array at least of size |nLocal| times |numEvals|
      to store columnwise the (local part of the) computed eigenvectors.

   :param primme: parameters structure.

   :return: error indicator; see :ref:`error-codes`.

cprimme
"""""""

.. c:function:: int cprimme(float *evals, PRIMME_COMPLEX_FLOAT *evecs, float *resNorms, primme_params *primme)

   Solve a Hermitian standard eigenproblem; see function :c:func:`sprimme`.

zprimme
"""""""

.. c:function:: int zprimme(double *evals, PRIMME_COMPLEX_DOUBLE *evecs, double *resNorms, primme_params *primme)

   Solve a Hermitian standard eigenproblem; see function :c:func:`dprimme`.

primme_initialize
"""""""""""""""""

.. c:function:: void primme_initialize(primme_params *primme)

   Set PRIMME parameters structure to the default values.

   :param primme: parameters structure.

primme_set_method
"""""""""""""""""

.. c:function:: int primme_set_method(primme_preset_method method, primme_params *primme)

   Set PRIMME parameters to one of the preset configurations.

   :param method: preset configuration; one of

      | |DYNAMIC|
      | |DEFAULT_MIN_TIME|
      | |DEFAULT_MIN_MATVECS|
      | |Arnoldi|
      | |GD|
      | |GD_plusK|
      | |GD_Olsen_plusK|
      | |JD_Olsen_plusK|
      | |RQI|
      | |JDQR|
      | |JDQMR|
      | |JDQMR_ETol|
      | |STEEPEST_DESCENT|
      | |LOBPCG_OrthoBasis|
      | |LOBPCG_OrthoBasis_Window|

   :param primme: parameters structure.

   See also :ref:`methods`.

primme_display_params
"""""""""""""""""""""

.. c:function:: void primme_display_params(primme_params primme)

   Display all printable settings of ``primme`` into the file descriptor |outputFile|.

   :param primme: parameters structure.

primme_free
"""""""""""

.. c:function:: void primme_free(primme_params *primme)

   Free memory allocated by PRIMME.

   :param primme: parameters structure.

.. include:: epilog.inc
