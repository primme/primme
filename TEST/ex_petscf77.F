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
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscmat.h>
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

        external generateLaplacian1D, PETScMatvec, ApplyPCPrecPETSC,
     :           par_GlobalSumDouble

!       Eigenvalues, eigenvectors, and their residual norms
!
        real*8   evals(NUMEmax), rnorms(NUMEmax)
        PetscScalar   evecs(n*NUMEmax)

!       Other vars
!
        Mat A
        PC pc
        COMMON A, pc
        PetscErrorCode ierr
        integer i,numProcs,procID,nLocal
        real*8  epsil, aNorm

!-----------------------------------------------------------------------
!       Start executable 
!-----------------------------------------------------------------------
!
        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
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
        call generateLaplacian1D(n, A, ierr)
        call primme_set_member_f77(primme, PRIMMEF77_matrixMatvec,
     :                                                     PETScMatvec)
        
!       Set parallel parameters
        call MatGetLocalSize(A, nLocal, PETSC_NULL_INTEGER, ierr)
        call primme_set_member_f77(primme, PRIMMEF77_nLocal, nLocal)
        call MPI_Comm_size(PETSC_COMM_WORLD, numProcs, ierr)
        call primme_set_member_f77(primme, PRIMMEF77_numProcs, numProcs)
        call MPI_Comm_rank(PETSC_COMM_WORLD, procID, ierr);
        call primme_set_member_f77(primme, PRIMMEF77_procID, procID)
        call primme_set_member_f77(primme, PRIMMEF77_globalSumDouble,
     :                                              par_GlobalSumDouble)

!       Set preconditioner  (optional)
        call PCCreate(PETSC_COMM_WORLD, pc, ierr)
        call PCSetType(pc, PCJACOBI, ierr)
        call PCSetOperators(pc, A, A, ierr)
        call PCSetFromOptions(pc, ierr)
        call PCSetUp(pc, ierr)
        call primme_set_member_f77(primme, 
     :       PRIMMEF77_applyPreconditioner, ApplyPCPrecPETSC)
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

        if (procID.eq.0) call primme_display_params_f77(primme)

!       ----------------------------------------------------------------
!       Calling the PRIMME solver
!       ----------------------------------------------------------------

#if defined(PETSC_USE_COMPLEX)
        call zprimme_f77(evals, evecs, rnorms, primme, ierr)
#else
        call dprimme_f77(evals, evecs, rnorms, primme, ierr)
#endif

!       ----------------------------------------------------------------
!       Reporting results

        if (procID.eq.0) then
            if (ierr.eq.0) then
              print *, 'PRIMME has returned successfully'
           else 
              print *, 'PRIMME returned with error: ', ierr
           endif

           call primme_display_stats_f77(primme)
!          
!          Example of obtaining primme members from the driver:
!          NOTE: don't use primme_get_member_f77, which can only be used in a callback
!
           call primmetop_get_member_f77(primme, PRIMMEF77_eps, epsil)
           call primmetop_get_member_f77(primme, PRIMMEF77_aNorm, aNorm)
           print '(A16,E8.2,A20,e12.5)', 'Tolerance used: ',epsil,
     :                                '  Estimated norm(A):',aNorm
!
!          Reporting of evals and residuals
!
           do i = 1, numemax
              write (*, 9000) i, evals(i),rnorms(i)
           enddo
 9000      FORMAT (1x,'E(',i1,') = ',G24.16,4x,
     &            'residual norm =', E12.4)
        endif

        call PetscFinalize(ierr)
        stop
        write(0,*) 'ERROR! No data in the file'
        stop
        end
!-----------------------------------------------------------------------
! Supporting subroutines
!-----------------------------------------------------------------------
!       ----------------------------------------------------------------


        subroutine generateLaplacian1D(n0,A,ierr)
!       ----------------------------------------------------------------
        implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
        integer n0
        PetscInt n, one, two, three
        Mat A

        PetscScalar value(3)
        PetscInt i, Istart,Iend,col(3)
        PetscErrorCode ierr

        call MatCreate(PETSC_COMM_WORLD, A, ierr)
        n = n0
        call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
        call MatSetUp(A, ierr)
        call MatGetOwnershipRange(A, Istart, Iend, ierr)
        one = 1
        two = 2
        three = 3
        do i=Istart,Iend-1
           if (i.eq.0) then
              col = (/0, 1, 0/)
              value = (/2.0, -1.0, 0.0/)
              call MatSetValues(A, one, i, two, col, value,
     :                                             INSERT_VALUES, ierr)
           else if (i.eq.n-1) then
              col = (/n-2, n-1, 0/)
              value = (/-1.0, 2.0, 0.0/)
              call MatSetValues(A, one, i, two, col, value,
     :                                             INSERT_VALUES, ierr)
           else
              col = (/i-1, i, i+1/)
              value = (/-1.0, 2.0, -1.0/)
              call MatSetValues(A, one, i, three, col, value,
     :                                             INSERT_VALUES, ierr)
           endif
           call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
           call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
        enddo
        end
        subroutine PETScMatvec(x,y,k,primme)
!       ----------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
        PetscScalar x(*), y(*)
        integer*8 primme
        integer k,j,nLocal
        Mat A
        COMMON A
        Vec xvec,yvec
        PetscErrorCode ierr

        call primme_get_member_f77(primme, PRIMMEF77_nLocal, nLocal)

        call MatCreateVecs(A, xvec, yvec, ierr)
        do j=0,k-1
           call VecPlaceArray(xvec, x(j*nLocal+1), ierr)
           call VecPlaceArray(yvec, y(j*nLocal+1), ierr)
           call MatMult(A, xvec, yvec, ierr)
           call VecResetArray(xvec, ierr)
           call VecResetArray(yvec, ierr)
        enddo
        call VecDestroy(xvec, ierr)
        call VecDestroy(yvec, ierr)
        end
        subroutine ApplyPCPrecPETSc(x,y,k,primme)
!       ----------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
        PetscScalar x(*), y(*)
        integer*8 primme
        integer k,j,nLocal
        Mat A
        PC pc
        COMMON A, pc
        Vec xvec,yvec
        PetscErrorCode ierr

        call primme_get_member_f77(primme, PRIMMEF77_nLocal, nLocal)

        call MatCreateVecs(A, xvec, yvec, ierr)
        do j=0,k-1
           call VecPlaceArray(xvec, x(j*nLocal+1), ierr)
           call VecPlaceArray(yvec, y(j*nLocal+1), ierr)
           call PCApply(pc, xvec, yvec, ierr)
           call VecResetArray(xvec, ierr)
           call VecResetArray(yvec, ierr)
        enddo
        call VecDestroy(xvec, ierr)
        call VecDestroy(yvec, ierr)
        end
        subroutine par_GlobalSumDouble(x,y,k,primme)
!       ----------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
#include <petsc/finclude/petscsys.h>
        real*8 x(*), y(*)
        integer*8 primme
        integer k, ierr
        call MPI_Allreduce(x, y, k, MPI_DOUBLE, MPI_SUM,
     :                                 PETSC_COMM_WORLD, ierr)
        end
