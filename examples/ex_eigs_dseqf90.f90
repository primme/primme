!*******************************************************************************
!  Copyright (c) 2018, College of William & Mary                                   
!  All rights reserved.                                                            
!                                                                                  
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are met:     
!      * Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.             
!      * Redistributions in binary form must reproduce the above copyright         
!        notice, this list of conditions and the following disclaimer in the       
!        documentation and/or other materials provided with the distribution.      
!      * Neither the name of the College of William & Mary nor the
!        names of its contributors may be used to endorse or promote products      
!        derived from this software without specific prior written permission.     
!                                                                                  
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
!  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          
!  DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY       
!  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      
!  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!  
!  PRIMME: https://github.com/primme/primme
!  Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
!******************************************************************************
!
!  Example to compute the k largest eigenvalues in a 1-D Laplacian matrix.
!
!******************************************************************************

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
   
   procedure(primme_eigs_matvec) :: MV, ApplyPrecon
   
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
   ierr = primme_set_member(primme, PRIMME_matrixMatvec, MV)
   
   ! Set preconditioner  (optional)
   ierr = primme_set_member(primme, PRIMME_applyPreconditioner, c_funloc(ApplyPrecon))
   ierr = primme_set_member(primme, PRIMME_correctionParams_precondition, 1_c_int64_t)
   
   !       Set a few other solver parameters (optional) 
   
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
   if (ierr.ne.0) then
     stop 1
   endif
   
   ! Reporting results
   
   if (ierr.eq.0) then
      print *, 'PRIMME has returned successfully'
   else 
      print *, 'PRIMME returned with error: ', ierr
   endif
   
   ierr = primme_get_member(primme, PRIMME_eps, epsil)
   ierr = primme_get_member(primme, PRIMME_aNorm, aNorm)
   ierr = primme_get_member(primme, PRIMME_stats_numOuterIterations, numIts)
   ierr = primme_get_member(primme, PRIMME_stats_numMatvecs, numMatvecs)
   print '(A,E8.2,/,A,e12.5,/,A,I8,/,A,I8)',                  &
                              'Tolerance used:   ',epsil,     &
                              'Estimated norm(A):',aNorm,     &
                              'Iterations:       ',numIts,    &
                              'Matvecs:          ',numMatvecs 

   ! Reporting of evals and residuals

   do i = 1, numemax
      print '(a,i1,a,G24.16,a,E12.4)',' eval(', i, ') = ', evals(i), '    residual norm =', rnorms(i)
   end do
   ierr = primme_params_destroy(primme)
   stop
end

!-----------------------------------------------------------------------
! Supporting subroutines
!-----------------------------------------------------------------------

! ----------------------------------------------------------------
! 1-D Laplacian block matrix-vector product, Y = A * X, where
!
! - X, input dense matrix of size primme.n x blockSize;
! - Y, output dense matrix of size primme.n x blockSize;
! - A, tridiagonal square matrix of dimension primme.n with this form:
!
!       2 -1  0  0  0 ... 
!      -1  2 -1  0  0 ... 
!       0 -1  2 -1  0 ... 
!       ...
!
subroutine MV(x,ldx,y,ldy,k,primme,ierr)
   use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double
   implicit none

   include 'primme_f90.inc'
   integer(c_int64_t) ldx, ldy
   real(c_double) x(ldx,*), y(ldy,*)
   type(c_ptr), value :: primme
   integer(c_int64_t) n, i
   integer(c_int) :: k,ierr
   integer :: j

   ierr = primme_get_member(primme, PRIMME_n, n)
   do j=1,k
      do i=1,n
         y(i,j) = 0
         if (i.ge.2) then
            y(i,j) = y(i,j) - x(i-1,j)
         endif
         y(i,j) = y(i,j) + 2.*x(i,j)
         if (i.le.n-1) then
            y(i,j) = y(i,j) - x(i+1,j)
         endif
      enddo
   enddo
   ierr = 0
end

! ----------------------------------------------------------------
! This performs Y = M^{-1} * X, where
!
! - X, input dense matrix of size primme.n x blockSize;
! - Y, output dense matrix of size primme.n x blockSize;
! - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
!
subroutine ApplyPrecon(x,ldx,y,ldy,k,primme, ierr)
   use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double, c_f_pointer
   implicit none
   include 'primme_f90.inc'
   integer(c_int64_t) ldx, ldy
   real(c_double) x(ldx,*), y(ldy,*)
   type(c_ptr), value :: primme
   integer(c_int64_t) n, i
   integer(c_int) :: k,ierr
   integer :: j

   real(c_double), pointer :: shifts(:)
   type(c_ptr) :: pshifts

   ierr = primme_get_member(primme, PRIMME_n, n)
   ierr = primme_get_member(primme, PRIMME_ShiftsForPreconditioner, pshifts)
   call c_f_pointer(pshifts, shifts, shape=[k])
   do j=1,k
      do i=1,n
         y(i,j) = x(i,j)/(2.0 - shifts(j))
      enddo
   enddo
   ierr = 0
end
