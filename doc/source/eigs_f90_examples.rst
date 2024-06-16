
.. _f90Simple:

Simple Fortran 90 Eigs Example
------------------------------

.. code:: 

   Program primmeF90Example
   use iso_c_binding, only: c_ptr, c_int64_t, c_double, c_funloc
   implicit none

   ! Add PRIMME constants and interfaces
   include 'primme_f90.inc'
   
   ! Pointer to the PRIMME data structure used internally by PRIMME
   type(c_ptr) :: primme
   
   ! Solver Parameters
   ! NOTE: variables passed to primme_set_member should be c_int64_t, c_double,
   !       c_funptr or a function with prototype primme_eigs_matvec
   
   integer(c_int64_t) :: n               = 100
   integer(c_int64_t) :: BASISmax        = 12
   integer(c_int64_t) :: NUMEmax         = 5
   integer(c_int64_t) :: BLOCKmax        = 1
   integer(c_int64_t) :: maxMatvecs      = 300000
   real(c_double)  :: ETOL            = 1.0D-12
   integer(c_int64_t) :: printLevel      = 5
   integer(c_int64_t), parameter :: numTargetShifts = 2
   real(c_double) :: TargetShifts(numTargetShifts) = (/3.0D0, 5.1D0/)
   
   ! Eigenvalues, eigenvectors, and their residual norms
   double precision, allocatable :: evals(:)
   double precision, allocatable :: rnorms(:)
   double precision, allocatable :: evecs(:,:)
   
   ! Other vars
   integer :: ierr
   integer(c_int64_t) :: i
   real(c_double) :: epsil, aNorm
   integer(c_int64_t) :: numIts, numMatvecs
   
   ! Initialize PRIMME
   primme = primme_params_create()
   
   ! Set a few basic solver parameters
   ierr = primme_set_member(primme, PRIMME_n, n)
   ierr = primme_set_member(primme, PRIMME_numEvals, NUMEmax)
   ierr = primme_set_member(primme, PRIMME_eps, ETOL)
   ierr = primme_set_member(primme, PRIMME_target, primme_closest_abs)
   ierr = primme_set_member(primme, PRIMME_numTargetShifts, numTargetShifts)
   ierr = primme_set_member(primme, PRIMME_targetShifts, TargetShifts)
   
   ! Set matvec 
   ierr = primme_set_member_matvec_double(primme, PRIMME_matrixMatvec, MV)
   
   ! Set preconditioner  (optional)
   ierr = primme_set_member(primme, PRIMME_applyPreconditioner, c_funloc(ApplyPrecon))
   ierr = primme_set_member(primme, PRIMME_correctionParams_precondition, 1_c_int64_t)
   
   ! Set a few other solver parameters (optional)
   ierr = primme_set_member(primme, PRIMME_maxBasisSize, BASISmax)
   ierr = primme_set_member(primme, PRIMME_maxBlockSize, BLOCKmax)
   ierr = primme_set_member(primme, PRIMME_printLevel, printLevel)
   ierr = primme_set_member(primme, PRIMME_maxMatvecs, maxMatvecs)
   
   ! Set the method to be used (after n, numEvals, and precondition have
   ! been set. Also after basisSize is set, if desired.)
   ierr = primme_set_method(PRIMME_DYNAMIC, primme)
   if (ierr .lt. 0) print *, 'No preset method. Using custom settings'
   
   ! Display what parameters are used
   call primme_display_params_f77(primme)
   
   ! Allocate arrays
   allocate(evals(NUMEmax))
   allocate(rnorms(NUMEmax))
   allocate(evecs(n,NUMEmax))
   
   ! Calling the PRIMME solver
   ierr = dprimme(evals, evecs, rnorms, primme)

