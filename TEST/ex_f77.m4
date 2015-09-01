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
define(`PRIMME_NUM', ifdef(`USE_PETSC', `PetscScalar', ifdef(`USE_COMPLEX', `complex*16', `real*8')))dnl

        Program primmeF77Example
!-----------------------------------------------------------------------
        implicit none
ifdef(`USE_PETSC', ``#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscmat.h>
'')dnl
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
ifdef(`USE_PETSC', `
        external generateLaplacian1D, PETScMatvec, ApplyPCPrecPETSC,
     :           par_GlobalSumDouble
', `
        external MV, ApplyPrecon
')dnl

!       Eigenvalues, eigenvectors, and their residual norms
!
        real*8   evals(NUMEmax), rnorms(NUMEmax)
        PRIMME_NUM   evecs(n*NUMEmax)

!       Other vars
!
ifdef(`USE_PETSC', `ifdef(`USE_POINTER',
`        Mat, target :: A
        PC, target :: pc
        MPI_Comm, target :: comm
', `        Mat A
        PC pc
        COMMON A, pc
')dnl
        PetscErrorCode ierr
        integer i,numProcs,procID,nLocal
', `        integer i,ierr
')dnl
        real*8  epsil, aNorm

!-----------------------------------------------------------------------
!       Start executable 
!-----------------------------------------------------------------------
!
ifdef(`USE_PETSC', `        call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
')dnl
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
ifdef(`USE_PETSC', `        call generateLaplacian1D(n, A, ierr)
ifdef(`USE_POINTER', `        call primme_set_member_f77(primme, PRIMMEF77_matrix, A)
')dnl
        call primme_set_member_f77(primme, PRIMMEF77_matrixMatvec,
     :                                                     PETScMatvec)
', `        call primme_set_member_f77(primme, PRIMMEF77_matrixMatvec, MV)
')dnl
        
ifdef(`USE_PETSC', `!       Set parallel parameters
        call MatGetLocalSize(A, nLocal, PETSC_NULL_INTEGER, ierr)
        call primme_set_member_f77(primme, PRIMMEF77_nLocal, nLocal)
        call MPI_Comm_size(PETSC_COMM_WORLD, numProcs, ierr)
        call primme_set_member_f77(primme, PRIMMEF77_numProcs, numProcs)
        call MPI_Comm_rank(PETSC_COMM_WORLD, procID, ierr);
        call primme_set_member_f77(primme, PRIMMEF77_procID, procID)
ifdef(`USE_POINTER', `        comm = PETSC_COMM_WORLD
        call primme_set_member_f77(primme, PRIMMEF77_commInfo, comm)
')dnl
        call primme_set_member_f77(primme, PRIMMEF77_globalSumDouble,
     :                                              par_GlobalSumDouble)
')dnl

!       Set preconditioner  (optional)
ifdef(`USE_PETSC', `        call PCCreate(PETSC_COMM_WORLD, pc, ierr)
        call PCSetType(pc, PCJACOBI, ierr)
        call PCSetOperators(pc, A, A, ierr)
        call PCSetFromOptions(pc, ierr)
        call PCSetUp(pc, ierr)
ifdef(`USE_POINTER', `        call primme_set_member_f77(primme, 
     :       PRIMMEF77_preconditioner, pc)
')dnl
        call primme_set_member_f77(primme, 
     :       PRIMMEF77_applyPreconditioner, ApplyPCPrecPETSC)
', `        call primme_set_member_f77(primme, 
     :       PRIMMEF77_applyPreconditioner, ApplyPrecon)
')dnl
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

       ifdef(`USE_PETSC', ` if (procID.eq.0)') call primme_display_params_f77(primme)

!       ----------------------------------------------------------------
!       Calling the PRIMME solver
!       ----------------------------------------------------------------
ifdef(`USE_PETSC', ``
#if defined(PETSC_USE_COMPLEX)
        call zprimme_f77(evals, evecs, rnorms, primme, ierr)
#else
        call dprimme_f77(evals, evecs, rnorms, primme, ierr)
#endif
'', `
        call ifdef(`USE_COMPLEX',`z', `d')primme_f77(evals, evecs, rnorms, primme, ierr)
')dnl

!       ----------------------------------------------------------------
!       Reporting results

ifdef(`USE_PETSC', ``        if (procID.eq.0) then
' define(sp, `   ')', `define(sp, `')')dnl
        sp()if (ierr.eq.0) then
        sp()   print *, 'PRIMME has returned successfully'
        sp()else 
        sp()   print *, 'PRIMME returned with error: ', ierr
        sp()endif

        sp()call primme_display_stats_f77(primme)
!       sp()
!       sp()Example of obtaining primme members from the driver:
!       sp()NOTE: don't use primme_get_member_f77, which can only be used in a callback
!
        sp()call primmetop_get_member_f77(primme, PRIMMEF77_eps, epsil)
        sp()call primmetop_get_member_f77(primme, PRIMMEF77_aNorm, aNorm)
        sp()print '(A16,E8.2,A20,e12.5)', 'Tolerance used: ',epsil,
     :  sp()                           '  Estimated norm(A):',aNorm
!
!       sp()Reporting of evals and residuals
!
        sp()do i = 1, numemax
        sp()   write (*, 9000) i, evals(i),rnorms(i)
        sp()enddo
 9000   sp()FORMAT (1x,'E(',i1,') = ',G24.16,4x,
     &  sp()       'residual norm =', E12.4)
ifdef(`USE_PETSC',`        endif

        call PetscFinalize(ierr)
')dnl
        stop
        write(0,*) 'ERROR! No data in the file'
        stop
        end
!-----------------------------------------------------------------------
! Supporting subroutines
!-----------------------------------------------------------------------
!       ----------------------------------------------------------------
changequote(`[',`]')
ifdef([USE_PETSC], [
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
ifdef([USE_POINTER], [        use iso_c_binding
])dnl
        implicit none
        include 'primme_f77.h'
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
        PRIMME_NUM x(*), y(*)
        integer*8 primme
        integer k,j,nLocal
ifdef([USE_POINTER], [        Mat, pointer :: A
        type(c_ptr) :: pA
], [        Mat A
        COMMON A
])dnl
        Vec xvec,yvec
        PetscErrorCode ierr

        call primme_get_member_f77(primme, PRIMMEF77_nLocal, nLocal)
ifdef([USE_POINTER], [        call primme_get_member_f77(primme, PRIMMEF77_matrix, pA)
        call c_f_pointer(pA, A)
])
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
ifdef([USE_POINTER], [        use iso_c_binding
])dnl
        implicit none
        include 'primme_f77.h'
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
        PRIMME_NUM x(*), y(*)
        integer*8 primme
        integer k,j,nLocal
ifdef([USE_POINTER], [        Mat, pointer :: A
        PC, pointer :: pc
        type(c_ptr) :: pA, ppc
], [        Mat A
        PC pc
        COMMON A, pc
])dnl
        Vec xvec,yvec
        PetscErrorCode ierr

        call primme_get_member_f77(primme, PRIMMEF77_nLocal, nLocal)
ifdef([USE_POINTER], [        call primme_get_member_f77(primme, PRIMMEF77_matrix, pA)
        call primme_get_member_f77(primme, PRIMMEF77_preconditioner,
     :                                                          ppc)
        call c_f_pointer(pA, A)
        call c_f_pointer(ppc, pc)
])
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
ifdef([USE_POINTER], [        use iso_c_binding
])dnl
        implicit none
        include 'primme_f77.h'
#include <petsc/finclude/petscsys.h>
        real*8 x(*), y(*)
        integer*8 primme
        integer k, ierr
ifdef([USE_POINTER], [        MPI_Comm, pointer :: comm
        type(c_ptr) :: pcomm

        call primme_get_member_f77(primme, PRIMMEF77_commInfo, pcomm)
        call c_f_pointer(pcomm, comm)
])dnl
        call MPI_Allreduce(x, y, k, MPI_DOUBLE, MPI_SUM,
     :                                 ifdef([USE_POINTER], [comm], [PETSC_COMM_WORLD]), ierr)
        end
], [
!       1-D Laplacian block matrix-vector product, Y = A * X, where
!      
!       - X, input dense matrix of size primme.n x blockSize;
!       - Y, output dense matrix of size primme.n x blockSize;
!       - A, tridiagonal square matrix of dimension primme.n with this form:
!      
!            [ 2 -1  0  0  0 ... ]
!            [-1  2 -1  0  0 ... ]
!            [ 0 -1  2 -1  0 ... ]
!             ...
!      
        subroutine MV(x,y,k,primme)
!       ----------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
        PRIMME_NUM x(*), y(*)
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
        PRIMME_NUM x(*), y(*)
        integer*8 primme
        integer k,i,j,n
        call primme_get_member_f77(primme, PRIMMEF77_n, n)
        do j=0,k-1
           do i=1,n
              y(j*n+i) = x(j*n+i)/2.0
           enddo
        enddo
        end
])dnl
