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
!  Example to compute the k largest singular values in a Lauchli matrix.
!
!******************************************************************************

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
   procedure(primme_svds_matvec) MV, ApplyPrecon

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
   ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_matrixMatvec, MV)

   ! Set preconditioner based on A^H*A (optional)
   ierr = primme_svds_set_member(primme_svds, PRIMME_SVDS_applyPreconditioner, ApplyPrecon)

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

   ! Reporting results

   if (ierr.eq.0) then
      print *, 'PRIMME_SVDS has returned successfully'
   else 
      print *, 'PRIMME_SVDS returned with error: ', ierr
   endif

   ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_eps, epsil)
   ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_aNorm, aNorm)
   print '(A16,E8.2,A20,e12.5)', 'Tolerance used: ',epsil, '  Estimated norm(A):',aNorm

   ! Reporting of svals and residuals

   do i = 1, NUMSmax
      print '(a,i1,a,G24.16,a,E12.4)',' sval(', i, ') = ', svals(i), '    residual norm =', rnorms(i)
   enddo
   stop
end

!-----------------------------------------------------------------------
! Supporting subroutines
!-----------------------------------------------------------------------


! ----------------------------------------------------------------
! Lauchli-like block matrix-vector products, Y = A * X or Y = A' * X, where
! 
! - X, input dense matrix of size primme.n x blockSize or primme.m x blockSize;
! - Y, output dense matrix of size primme.m x blockSize or primme.n x blockSize;
! - A, rectangular matrix of size primme.m x primme.n with this form:
! 
!       1  1  1  1  1 ... ,  ei = 1 - (1 - c)*i/(min(m,n) - 1)
!      e0  0  0  0  0 ... 
!       0 e1  0  0  0 ... 
!       0  0 e2  0  0 ... 
!       ...
!      
subroutine MV(x,ldx,y,ldy,k,transpose,primme_svds,ierr)
   use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double
   implicit none

   intrinsic min
   include 'primme_f90.inc'
   integer(c_int64_t) :: ldx,ldy
   real(c_double) :: x(ldx,*), y(ldy,*)
   type(c_ptr), value :: primme_svds
   integer(c_int) :: k, transpose, ierr

   integer(c_int64_t) :: m, n, i
   integer j

   double precision :: c
   common c

   ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_m, m)
   ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_n, n)

   if (transpose.eq.0) then
      do j=1,k
         y(1,j) = 0
         do i=1,n
            y(1,j) = y(1,j) + x(i,j)
         enddo
         do i=2,m
            if (i-1.le.n) then
               y(i,j) = x(i-1,j)*(1.0 - (1.0-c)*(i-2)/(min(m,n)-1))
            else
               y(i,j) = 0
            endif
         enddo
      enddo
   else
      do j=1,k
         do i=1,n
            if (i+1.le.m) then
               y(i,j) = x(1,j) + x(i+1,j) * (1.0-(1.0-c)*(i-1)/(min(m,n)-1))
            else
               y(i,j) = x(1,j)
            endif
         enddo
      enddo
   endif
   ierr = 0
end

! ----------------------------------------------------------------
! This performs Y = M^{-1} * X, where
!
! - X, input dense matrix of size primme.n x blockSize;
! - Y, output dense matrix of size primme.n x blockSize;
! - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
!      
subroutine ApplyPrecon(x,ldx,y,ldy,k,mode,primme_svds,ierr)
   use iso_c_binding, only: c_ptr, c_int, c_int64_t, c_double, c_f_pointer
   implicit none

   intrinsic min
   include 'primme_f90.inc'
   integer(c_int64_t) :: ldx,ldy,m,n
   real(c_double) :: x(ldx,*), y(ldy,*)
   type(c_ptr), value :: primme_svds
   integer(c_int) :: k,mode,ierr
   integer(c_int64_t) i
   integer j

   double precision :: c, ei
   common c
   real(c_double), pointer :: shifts(:)
   type(c_ptr) :: pshifts
   integer(c_int64_t) :: numTargetShifts

   ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_m, m)
   ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_n, n)
   ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_targetShifts, pshifts)
   ierr = primme_svds_get_member(primme_svds, PRIMME_SVDS_numTargetShifts, numTargetShifts)
   call c_f_pointer(pshifts, shifts, shape=[numTargetShifts])

   if (mode.eq.PRIMME_SVDS_op_AtA) then
      do j=1,k
         do i=1,n
            if (i-1.le.m) then
               ei = 1.0 - (1.0 - c)*(i-1)/(min(m,n) - 1)
            else
               ei = 0
            endif
            y(i,j) = x(i,j)/(1.0 + ei*ei - shifts(1)*shifts(1))
         enddo
      enddo
   else if (mode.eq.PRIMME_SVDS_op_AAt) then
      do j=1,k
         y(1,j) = x(1,j)/m
         do i=2,m
            if (i-2.le.n) then
               ei = 1.0 - (1.0 - c)*(i-2)/(min(m,n) - 1)
            else
               ei = 0.0
            endif
            y(i,j) = x(i,j)/(ei*ei - shifts(1)*shifts(1))
         enddo
      enddo
   else if (mode.eq.PRIMME_SVDS_op_augmented) then
      ! If any preconditioner is available, just y = x
      y(1:m+n,1:k) = x(1:m+n,1:k)
   endif
   ierr = 0
end
