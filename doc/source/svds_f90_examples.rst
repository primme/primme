.. _f90SvdsSimple:

Simple Fortran 90 SVDS Example
------------------------------

.. code::
     
    Program primmeSvdsF77Example

    use iso_c_binding, only: c_ptr, c_int64_t, c_double, c_funloc
    implicit none

    ! Add PRIMME constants and interfaces
    include 'primme_f90.inc'

    ! Pointer to the PRIMME SVDS data structure used internally by PRIMME
    type(c_ptr) primme_svds

    ! Solver Parameters
    ! NOTE: variables passed to primme_svds_set_member should be c_int64_t, c_double,
    !       c_funptr or a function with prototype primme_svds_matvec

    integer(c_int64_t) :: m               = 101
    integer(c_int64_t) :: n               = 100
    integer(c_int64_t) :: BASISmax        = 12
    integer(c_int64_t) :: NUMSmax         = 5
    integer(c_int64_t) :: BLOCKmax        = 1
    integer(c_int64_t) :: maxMatvecs      = 300000
    real(c_double)  :: STOL            = 1.0D-12
    integer(c_int64_t) :: printLevel      = 2
    integer(c_int64_t), parameter :: numTargetShifts = 1
    real(c_double) :: TargetShifts(numTargetShifts) = (/0.5D0/)
    double precision :: c

    common c

    ! Singular values, vectors and their residual norms

    double precision, allocatable :: svals(:)
    double precision, allocatable :: rnorms(:)
    double precision, allocatable :: svecs(:,:)

    ! Other vars

    integer :: ierr
    integer(c_int64_t) :: i
    real(c_double)  epsil, aNorm

    ! Initialize PRIMME
    primme_svds = primme_svds_params_create()

    c = 1.0D-4

    ! Set a few basic solver parameters
    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_m, m)
    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_n, n)
    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_numSvals, NUMSmax)
    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_eps, STOL)
    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_target, primme_svds_closest_abs)
    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_numTargetShifts, numTargetShifts)
    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_targetShifts, TargetShifts)

    ! Set matvec 
    ierr = primme_svds_set_member_matvec_double(primme_svds, PRIMME_SVDS_matrixMatvec, MV)

    ! Set preconditioner based on A^H*A (optional)
    ierr = primme_svds_set_member_matvec_double(primme_svds, PRIMME_SVDS_applyPreconditioner, ApplyPrecon)

    ! Set a few other solver parameters (optional) 

    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_maxBasisSize, BASISmax)
    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_maxBlockSize, BLOCKmax)
    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_printLevel, printLevel)
    ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_maxMatvecs, maxMatvecs)

    ! Set the method to be used (after m, n, numSvals, and precondition have
    ! been set. Also after maxBasisSize is set, if desired.)

    ierr = primme_svds_set_method(PRIMME_SVDS_default, PRIMME_DEFAULT_METHOD, PRIMME_DEFAULT_METHOD, primme_svds)
    if (ierr .lt. 0) write(*,*) 'No preset method. Using custom settings'

    ! Display what parameters are used
    call primme_svds_display_params_f77(primme_svds)

    ! Allocate arrays
    
    allocate(svals(NUMSmax))
    allocate(rnorms(NUMSmax))
    allocate(svecs(m+n,NUMSmax))

    ! Calling the PRIMME solver

    ierr = dprimme_svds(svals, svecs, rnorms, primme_svds)
    if (ierr.ne.0) then
        stop 1
    endif