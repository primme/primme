C*******************************************************************************
C  Copyright (c) 2017, College of William & Mary                                   
C  All rights reserved.                                                            
C                                                                                  
C  Redistribution and use in source and binary forms, with or without
C  modification, are permitted provided that the following conditions are met:     
C      * Redistributions of source code must retain the above copyright
C        notice, this list of conditions and the following disclaimer.             
C      * Redistributions in binary form must reproduce the above copyright         
C        notice, this list of conditions and the following disclaimer in the       
C        documentation and/or other materials provided with the distribution.      
C      * Neither the name of the College of William & Mary nor the
C        names of its contributors may be used to endorse or promote products      
C        derived from this software without specific prior written permission.     
C                                                                                  
C  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
C  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
C  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          
C  DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY       
C  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      
C  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
C  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
C  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
C  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C  
C  PRIMME: https://github.com/primme/primme
C  Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
*******************************************************************************
*
*  Example to compute the k largest singular values in a Lauchli matrix.
*
*******************************************************************************

        Program primmeSvdsF77Example
!-----------------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Pointer to the PRIMME SVDS data structure used internally by PRIMME
!
!       Note that for 64 bit systems, pointers are 8 bytes so use:
        integer*8 primme_svds
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Problem setup
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        ! Solver Parameters
        integer*8 m,n,NUMSmax,BASISmax,BLOCKmax,maxMatvecs,
     :          printLevel, whichSvals, numTargetShifts
        real*8 TOL, c

        parameter (
     :            m               = 101,
     :            n               = 100,
     :            BASISmax        = 12,
     :            NUMSmax         = 5,
     :            BLOCKmax        = 1,
     :            maxMatvecs      = 300000,
     :            TOL             = 1.0D-12,
     :            printLevel      = 2,
     :            whichSvals      = primme_svds_closest_abs,
     :            numTargetShifts = 1
     :  )
        common c
        real*8 TargetShifts(numTargetShifts)
        data TargetShifts /0.5D0/

        external MV, ApplyPrecon

!       Singular values, vectors and their residual norms
!
        real*8 svals(NUMSmax)
        real*8 rnorms(NUMSmax)
        complex*16 svecs((m+n)*NUMSmax)

!       Other vars
!
        integer i,ierr
        real*8  epsil, aNorm

!-----------------------------------------------------------------------
!       Start executable 
!-----------------------------------------------------------------------
!
!       ----------------------------------------------------------------
!       Initialize PRIMME
!       ----------------------------------------------------------------
!
        call primme_svds_initialize_f77(primme_svds)

        c = 1.0D-4

!       Set a few basic solver parameters
        call primme_svds_set_member_f77(primme_svds, PRIMME_SVDS_m,
     :                                                          m, ierr)
        call primme_svds_set_member_f77(primme_svds, PRIMME_SVDS_n,
     :                                                          n, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                           PRIMME_SVDS_numSvals, NUMSmax, ierr)
        call primme_svds_set_member_f77(primme_svds, PRIMME_SVDS_eps,
     :                                                       TOL, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                         PRIMME_SVDS_target, whichSvals, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :            PRIMME_SVDS_numTargetShifts, numTargetShifts, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                 PRIMME_SVDS_targetShifts, TargetShifts, ierr)

!       Set matvec 
        call primme_svds_set_member_f77(primme_svds,
     :                            PRIMME_SVDS_matrixMatvec, MV, ierr)
        

!       Set preconditioner based on A^H*A (optional)
        call primme_svds_set_member_f77(primme_svds, 
     :       PRIMME_SVDS_applyPreconditioner, ApplyPrecon, ierr)
!
!       Set a few other solver parameters (optional) 
!
        call primme_svds_set_member_f77(primme_svds,
     :                      PRIMME_SVDS_maxBasisSize, BASISmax, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                      PRIMME_SVDS_maxBlockSize, BLOCKmax, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                      PRIMME_SVDS_printLevel, printLevel, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                      PRIMME_SVDS_maxMatvecs, maxMatvecs, ierr)
!
!       Set the method to be used (after m, n, numSvals, and precondition have
!       been set. Also after maxBasisSize is set, if desired.)
        call primme_svds_set_method_f77(PRIMME_SVDS_default,
     :        PRIMME_DEFAULT_METHOD, PRIMME_DEFAULT_METHOD,
     :        primme_svds, ierr)

        if (ierr .lt. 0) 
     :     write(*,*) 'No preset method. Using custom settings'

!       ----------------------------------------------------------------
!       Display what parameters are used
!       ----------------------------------------------------------------

       call primme_svds_display_params_f77(primme_svds)

!       ----------------------------------------------------------------
!       Calling the PRIMME solver
!       ----------------------------------------------------------------

        call  zprimme_svds_f77(svals, svecs, rnorms, primme_svds, ierr)
        if (ierr.ne.0) then
          stop 1
        endif

!       ----------------------------------------------------------------
!       Reporting results

        if (ierr.eq.0) then
           print *, 'PRIMME_SVDS has returned successfully'
        else 
           print *, 'PRIMME_SVDS returned with error: ', ierr
        endif

!       
!       Example of obtaining primme members from the driver:
!       NOTE: don't use primme_svds_get_member_f77, which can only be used in a callback
!
        call primme_svdstop_get_member_f77(primme_svds,
     :                      PRIMME_SVDS_eps, epsil, ierr)
        call primme_svdstop_get_member_f77(primme_svds,
     :                      PRIMME_SVDS_aNorm, aNorm, ierr)
        print '(A16,E8.2,A20,e12.5)', 'Tolerance used: ',epsil,
     :                             '  Estimated norm(A):',aNorm
!
!       Reporting of svals and residuals
!
        do i = 1, NUMSmax
           write (*, 9000) i, svals(i),rnorms(i)
        enddo
 9000   FORMAT (1x,'E(',i1,') = ',G24.16,4x,
     &         'residual norm =', E12.4)
        stop
        write(0,*) 'ERROR! No data in the file'
        stop
        end
!-----------------------------------------------------------------------
! Supporting subroutines
!-----------------------------------------------------------------------
!       ----------------------------------------------------------------


!       Lauchli-like block matrix-vector products, Y = A * X or Y = A' * X,
!       where
!      
!       - X, input dense matrix of size primme.n x blockSize or primme.m x blockSize;
!       - Y, output dense matrix of size primme.m x blockSize or primme.n x blockSize;
!       - A, rectangular matrix of size primme.m x primme.n with this form:
!      
!             1  1  1  1  1 ... ,  ei = 1 - (1 - c)*i/(min(m,n) - 1)
!            e0  0  0  0  0 ... 
!             0 e1  0  0  0 ... 
!             0  0 e2  0  0 ... 
!             ...
!      
        subroutine MV(x,ldx,y,ldy,k,transpose,primme_svds,ierr)
!       ----------------------------------------------------------------
        implicit none
        intrinsic min
        include 'primme_f77.h'
        integer*8 ldx,ldy
        complex*16 x(ldx,*), y(ldy,*)
        integer*8 primme_svds
        integer*8 m, n, i
        integer k,transpose,j, ierr
        real*8 c
        common c
        call primme_svds_get_member_f77(primme_svds, PRIMME_SVDS_m,
     :                                                          m, ierr)
        call primme_svds_get_member_f77(primme_svds, PRIMME_SVDS_n,
     :                                                          n, ierr)
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
                    y(i,j) = x(1,j)
     :                     + x(i+1,j) * (1.0-(1.0-c)*(i-1)/(min(m,n)-1))
                 else
                    y(i,j) = x(1,j)
                 endif
              enddo
           enddo
        endif
        ierr = 0
        end

!       This performs Y = M^{-1} * X, where
!      
!       - X, input dense matrix of size primme.n x blockSize;
!       - Y, output dense matrix of size primme.n x blockSize;
!       - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
!      
        subroutine ApplyPrecon(x,ldx,y,ldy,k,mode,primme_svds,ierr)
!       ----------------------------------------------------------------
        implicit none
        intrinsic min
        include 'primme_f77.h'
        integer*8 ldx,ldy,m,n
        complex*16 x(ldx,*), y(ldy,*)
        integer*8 primme_svds
        integer k,mode,i,j,ierr
        real*8 c, ei, shift
        common c
        call primme_svds_get_member_f77(primme_svds, PRIMME_SVDS_m,
     :                                                          m, ierr)
        call primme_svds_get_member_f77(primme_svds, PRIMME_SVDS_n,
     :                                                          n, ierr)
        call primme_svds_get_member_f77(primme_svds,
     :                         PRIMME_SVDS_targetShifts, shift, ierr)
        if (mode.eq.PRIMME_SVDS_op_AtA) then
           do j=1,k
              do i=1,n
                 if (i-1.le.m) then
                    ei = 1.0 - (1.0 - c)*(i-1)/(min(m,n) - 1)
                 else
                    ei = 0
                 endif
                 y(i,j) = x(i,j)/(1.0 + ei*ei - shift*shift)
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
                 y(i,j) = x(i,j)/(ei*ei - shift*shift)
              enddo
           enddo
        else if (mode.eq.PRIMME_SVDS_op_augmented) then
!          If any preconditioner is available, just y = x
           y(1:m+n,1:k) = x(1:m+n,1:k)
        endif
        ierr = 0
        end
