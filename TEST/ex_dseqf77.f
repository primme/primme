*******************************************************************************
*   PRIMME PReconditioned Iterative MultiMethod Eigensolver
*   Copyright (C) 2015 College of William & Mary,
*   James R. McCombs, Eloy Romero Alcalde, Andreas Stathopoulos, Lingfei Wu
*
*   This file is part of PRIMME.
*
*   PRIMME is free software; you can redistribute it and/or
*   modify it under the terms of the GNU Lesser General Public
*   License as published by the Free Software Foundation; either
*   version 2.1 of the License, or (at your option) any later version.
*
*   PRIMME is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*   Lesser General Public License for more details.
*
*   You should have received a copy of the GNU Lesser General Public
*   License along with this library; if not, write to the Free Software
*   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
*   02110-1301  USA
*
*******************************************************************************
*
*  Example to compute the k largest eigenvalues in a 1-D Laplacian matrix.
*
*******************************************************************************

        Program primmeF77Example
!-----------------------------------------------------------------------
        implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Pointer to the PRIMME data structure used internally by PRIMME
!
!       Note that for 64 bit systems, pointers are 8 bytes so use:
        integer*8 primme
        include 'primme_f77.h'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!       Problem setup
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        ! Solver Parameters
        integer n,NUMEmax,BASISmax,BLOCKmax,maxMatvecs,
     :          printLevel, method, whichEvals, numTargetShifts
        real*8 ETOL

        parameter (
     :            n               = 100,
     :            BASISmax        = 12,
     :            NUMEmax         = 5,
     :            BLOCKmax        = 1,
     :            maxMatvecs      = 300000,
     :            ETOL            = 1.0D-14,
     :            printLevel      = 5,
     :            whichEvals      = PRIMMEF77_smallest,
     :            numTargetShifts = 2,
     :            method          = PRIMMEF77_DYNAMIC
     :  )
        real*8 TargetShifts(numTargetShifts)
        data TargetShifts /3.0D0, 5.1D0/

        external MV, ApplyPrecon

!       Eigenvalues, eigenvectors, and their residual norms
!
        real*8   evals(NUMEmax), rnorms(NUMEmax)
        real*8   evecs(n*NUMEmax)

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
        call primme_initialize_f77(primme)

!       Set a few basic solver parameters
        call primme_set_member_f77(primme, PRIMMEF77_n, n)
        call primme_set_member_f77(primme, PRIMMEF77_numEvals, NUMEmax)
        call primme_set_member_f77(primme, PRIMMEF77_eps, ETOL)
        call primme_set_member_f77(primme, PRIMMEF77_target, whichEvals)
        call primme_set_member_f77(primme, PRIMMEF77_numTargetShifts, 
     :                                                 numTargetShifts)
        call primme_set_member_f77(primme, PRIMMEF77_targetShifts, 
     :                                                    TargetShifts)

!       Set matvec 
        call primme_set_member_f77(primme, PRIMMEF77_matrixMatvec, MV)
        

!       Set preconditioner  (optional)
        call primme_set_member_f77(primme, 
     :       PRIMMEF77_applyPreconditioner, ApplyPrecon)
        call primme_set_member_f77(primme, 
     :       PRIMMEF77_correctionParams_precondition, 0)
!
!       Set a few other solver parameters (optional) 
!
        call primme_set_member_f77(primme, PRIMMEF77_maxBasisSize, 
     :                                                        BASISmax)
        call primme_set_member_f77(primme, PRIMMEF77_maxBlockSize,
     :                                                        BLOCKmax)
        call primme_set_member_f77(primme, PRIMMEF77_printLevel, 
     :                                                      printLevel)
        call primme_set_member_f77(primme, PRIMMEF77_maxMatvecs,
     :                                                      maxMatvecs)
        call primme_set_member_f77(primme, 
     :              PRIMMEF77_restartingParams_scheme, PRIMMEF77_thick)
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

        call dprimme_f77(evals, evecs, rnorms, primme, ierr)

!       ----------------------------------------------------------------
!       Reporting results

        if (ierr.eq.0) then
           print *, 'PRIMME has returned successfully'
        else 
           print *, 'PRIMME returned with error: ', ierr
        endif

        call primme_display_stats_f77(primme)
!       
!       Example of obtaining primme members from the driver:
!       NOTE: don't use primme_get_member_f77, which can only be used in a callback
!
        call primmetop_get_member_f77(primme, PRIMMEF77_eps, epsil)
        call primmetop_get_member_f77(primme, PRIMMEF77_aNorm, aNorm)
        print '(A16,E8.2,A20,e12.5)', 'Tolerance used: ',epsil,
     :                             '  Estimated norm(A):',aNorm
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
        subroutine MV(x,y,k,primme)
!       ----------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
        real*8 x(*), y(*)
        integer*8 primme
        integer k,i,j,n
        call primme_get_member_f77(primme, PRIMMEF77_n, n)
        do j=0,k-1
           do i=1,n
              y(j*n+i) = 0
              if (i.ge.2) then
                 y(j*n+i) = y(j*n+i) - x(j*n+i-1)
              endif
              y(j*n+i) = y(j*n+i) + 2.*x(j*n+i)
              if (i.le.n-1) then
                 y(j*n+i) = y(j*n+i) - x(j*n+i+1)
              endif
           enddo
        enddo
        end

!       This performs Y = M^{-1} * X, where
!      
!       - X, input dense matrix of size primme.n x blockSize;
!       - Y, output dense matrix of size primme.n x blockSize;
!       - M, diagonal square matrix of dimension primme.n with 2 in the diagonal.
!      
        subroutine ApplyPrecon(x,y,k,primme)
!       ----------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
        real*8 x(*), y(*)
        integer*8 primme
        integer k,i,j,n
        call primme_get_member_f77(primme, PRIMMEF77_n, n)
        do j=0,k-1
           do i=1,n
              y(j*n+i) = x(j*n+i)/2.0
           enddo
        enddo
        end
