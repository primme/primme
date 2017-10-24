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
*  Example to compute the k largest eigenvalues in a 1-D Laplacian matrix.
*
*******************************************************************************

        Program primmeF77Example
!-----------------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Pointer to the PRIMME data structure used internally by PRIMME
!
!       Note that for 64 bit systems, pointers are 8 bytes so use:
        integer*8 primme
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Problem setup
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        ! Solver Parameters
        integer*8 n,NUMEmax,BASISmax,BLOCKmax,maxMatvecs,
     :          printLevel, method, whichEvals, numTargetShifts
        real*8 ETOL

        parameter (
     :            n               = 100,
     :            BASISmax        = 12,
     :            NUMEmax         = 5,
     :            BLOCKmax        = 1,
     :            maxMatvecs      = 300000,
     :            ETOL            = 1.0D-12,
     :            printLevel      = 5,
     :            whichEvals      = primme_smallest,
     :            numTargetShifts = 2,
     :            method          = PRIMME_DYNAMIC
     :  )
        real*8 TargetShifts(numTargetShifts)
        data TargetShifts /3.0D0, 5.1D0/

        external MV, ApplyPrecon

!       Eigenvalues, eigenvectors, and their residual norms
!
        real*8 evals(NUMEmax)
        real*8 rnorms(NUMEmax)
        complex*16 evecs(n*NUMEmax)

!       Other vars
!
        integer i,ierr
        real*8  epsil, aNorm
        integer*8 numIts, numMatvecs

!-----------------------------------------------------------------------
!       Start executable 
!-----------------------------------------------------------------------
!
!       ----------------------------------------------------------------
!       Initialize PRIMME
!       ----------------------------------------------------------------
!
        call primme_initialize_f77(primme)

!       Set a few basic solver parameters
        call primme_set_member_f77(primme, PRIMME_n, n, ierr)
        call primme_set_member_f77(primme, PRIMME_numEvals, NUMEmax,
     :                                                             ierr)
        call primme_set_member_f77(primme, PRIMME_eps, ETOL, ierr)
        call primme_set_member_f77(primme, PRIMME_target,
     :                                                 whichEvals, ierr)
        call primme_set_member_f77(primme, PRIMME_numTargetShifts, 
     :                                            numTargetShifts, ierr)
        call primme_set_member_f77(primme, PRIMME_targetShifts, 
     :                                               TargetShifts, ierr)

!       Set matvec 
        call primme_set_member_f77(primme, PRIMME_matrixMatvec,
     :                                                         MV, ierr)
        

!       Set preconditioner  (optional)
        call primme_set_member_f77(primme, 
     :       PRIMME_applyPreconditioner, ApplyPrecon, ierr)
        call primme_set_member_f77(primme, 
     :       PRIMME_correctionParams_precondition, 0, ierr)
!
!       Set a few other solver parameters (optional) 
!
        call primme_set_member_f77(primme, PRIMME_maxBasisSize, 
     :                                                   BASISmax, ierr)
        call primme_set_member_f77(primme, PRIMME_maxBlockSize,
     :                                                   BLOCKmax, ierr)
        call primme_set_member_f77(primme, PRIMME_printLevel, 
     :                                                 printLevel, ierr)
        call primme_set_member_f77(primme, PRIMME_maxMatvecs,
     :                                                 maxMatvecs, ierr)
        call primme_set_member_f77(primme, 
     :         PRIMME_restartingParams_scheme, PRIMME_thick, ierr)
!
!       Set the method to be used (after n, numEvals, and precondition have
!       been set. Also after basisSize is set, if desired.)
        call primme_set_method_f77(primme, method, ierr)

        if (ierr .lt. 0) 
     :     write(*,*) 'No preset method. Using custom settings'

!       ----------------------------------------------------------------
!       Display what parameters are used
!       ----------------------------------------------------------------

        call primme_display_params_f77(primme)

!       ----------------------------------------------------------------
!       Calling the PRIMME solver
!       ----------------------------------------------------------------

        call zprimme_f77(evals, evecs, rnorms, primme, ierr)
        if (ierr.ne.0) then
          stop 1
        endif

!       ----------------------------------------------------------------
!       Reporting results

        if (ierr.eq.0) then
           print *, 'PRIMME has returned successfully'
        else 
           print *, 'PRIMME returned with error: ', ierr
        endif

!       
!       Example of obtaining primme members from the driver:
!       NOTE: don't use primme_get_member_f77, which can only be used in a callback
!
        call primmetop_get_member_f77(primme, PRIMME_eps, epsil,
     :                                                        ierr)
        call primmetop_get_member_f77(primme, PRIMME_aNorm,
     :                                                aNorm, ierr)
        call primmetop_get_member_f77(primme,
     :           PRIMME_stats_numOuterIterations, numIts, ierr)
        call primmetop_get_member_f77(primme,
     :                PRIMME_stats_numMatvecs, numMatvecs, ierr)
        print '(A,E8.2,/,A,e12.5,/,A,I8,/,A,I8)',
     :                             'Tolerance used:   ',epsil,
     :                             'Estimated norm(A):',aNorm,
     :                             'Iterations:       ',numIts,
     :                             'Matvecs:          ',numMatvecs
!
!       Reporting of evals and residuals
!
        do i = 1, numemax
           write (*, 9000) i, evals(i),rnorms(i)
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


!       1-D Laplacian block matrix-vector product, Y = A * X, where
!      
!       - X, input dense matrix of size primme.n x blockSize;
!       - Y, output dense matrix of size primme.n x blockSize;
!       - A, tridiagonal square matrix of dimension primme.n with this form:
!      
!             2 -1  0  0  0 ... 
!            -1  2 -1  0  0 ... 
!             0 -1  2 -1  0 ... 
!             ...
!      
        subroutine MV(x,ldx,y,ldy,k,primme,ierr)
!       ----------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
        integer*8 ldx, ldy
        complex*16 x(ldx,*), y(ldy,*)
        integer*8 primme
        integer*8 n, i
        integer k,ierr,j
        call primme_get_member_f77(primme, PRIMME_n, n, ierr)
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

!       This performs Y = M^{-1} * X, where
!      
!       - X, input dense matrix of size primme.n x blockSize;
!       - Y, output dense matrix of size primme.n x blockSize;
!       - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
!      
        subroutine ApplyPrecon(x,ldx,y,ldy,k,primme, ierr)
!       ----------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
        integer*8 ldx, ldy
        complex*16 x(ldx,*), y(ldy,*)
        integer*8 primme
        integer*8 n, i
        integer k,ierr,j
        call primme_get_member_f77(primme, PRIMME_n, n, ierr)
        do j=1,k
           do i=1,n
              y(i,j) = x(i,j)/2.0
           enddo
        enddo
        ierr = 0
        end
