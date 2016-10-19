C*******************************************************************************
C  Copyright (c) 2016, College of William & Mary                                   
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
define(`PRIMME_NUM', ifdef(`USE_PETSC', `PetscScalar', ifdef(`USE_COMPLEX', `complex*16', `real*8')))dnl

        Program primmeSvdsF77Example
!-----------------------------------------------------------------------
        implicit none
        include 'primme_f77.h'
ifdef(`USE_PETSC', ``#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscmat.h>
'')dnl
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
ifdef(`USE_PETSC', `',`        common c')
        real*8 TargetShifts(numTargetShifts)
        data TargetShifts /0.5D0/
ifdef(`USE_PETSC', `
        external generateLauchli, PETScMatvec, ApplyPCPrecAHA,
     :           par_GlobalSum
', `
        external MV, ApplyPrecon
')dnl

!       Singular values, vectors and their residual norms
!
        real*8   svals(NUMSmax), rnorms(NUMSmax)
        PRIMME_NUM   svecs((m+n)*NUMSmax)

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
        Mat AH, AHA
        PetscErrorCode ierr
        integer i,numProcs,procID,mLocal,nLocal
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
ifdef(`USE_PETSC', `        call generateLauchli(m, n, c, A, ierr)
ifdef(`USE_POINTER', `        call primme_svds_set_member_f77(primme_svds,
     :                                   PRIMME_SVDS_matrix, A, ierr)
')dnl
        call primme_svds_set_member_f77(primme_svds,
     :                   PRIMME_SVDS_matrixMatvec, PETScMatvec, ierr)
', `        call primme_svds_set_member_f77(primme_svds,
     :                            PRIMME_SVDS_matrixMatvec, MV, ierr)
')dnl
        
ifdef(`USE_PETSC', `!       Set parallel parameters
        call MatGetLocalSize(A, mLocal, nLocal, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                              PRIMME_SVDS_mLocal, mLocal, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                              PRIMME_SVDS_nLocal, nLocal, ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD, numProcs, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                         PRIMME_SVDS_numProcs, numProcs, ierr)
        call MPI_Comm_rank(PETSC_COMM_WORLD, procID, ierr);
        call primme_svds_set_member_f77(primme_svds,
     :                              PRIMME_SVDS_procID, procID, ierr)
ifdef(`USE_POINTER', `        comm = PETSC_COMM_WORLD
        call primme_svds_set_member_f77(primme_svds,
     :                              PRIMME_SVDS_commInfo, comm, ierr)
')dnl
        call primme_svds_set_member_f77(primme_svds,
     :        PRIMME_SVDS_globalSumReal, par_GlobalSum, ierr)
')dnl

!       Set preconditioner based on A^H*A (optional)
ifdef(`USE_PETSC', `        call PCCreate(PETSC_COMM_WORLD, pc, ierr)
        call MatHermitianTranspose(A, MAT_INITIAL_MATRIX, AH, ierr)
        call MatMatMult(AH, A, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL,
     :                                                        AHA, ierr)
        call PCSetType(pc, PCJACOBI, ierr)
        call PCSetOperators(pc, AHA, AHA, ierr)
        call PCSetFromOptions(pc, ierr)
        call PCSetUp(pc, ierr)
ifdef(`USE_POINTER', `        call primme_svds_set_member_f77(primme_svds, 
     :       PRIMME_SVDS_preconditioner, pc, ierr)
')dnl
        call primme_svds_set_member_f77(primme_svds, 
     :       PRIMME_SVDS_applyPreconditioner, ApplyPCPrecAHA, ierr)
', `        call primme_svds_set_member_f77(primme_svds, 
     :       PRIMME_SVDS_applyPreconditioner, ApplyPrecon, ierr)
')dnl
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

       ifdef(`USE_PETSC', ` if (procID.eq.0)')call primme_svds_display_params_f77(primme_svds)

!       ----------------------------------------------------------------
!       Calling the PRIMME solver
!       ----------------------------------------------------------------
ifdef(`USE_PETSC', ``
#if defined(PETSC_USE_COMPLEX)
        call zprimme_svds_f77(svals, svecs, rnorms, primme_svds, ierr)
#else
        call dprimme_svds_f77(svals, svecs, rnorms, primme_svds, ierr)
#endif
'', `
        call ifdef(`USE_COMPLEX',`z', `d')primme_svds_f77(svals, svecs, rnorms, primme_svds, ierr)
')dnl

!       ----------------------------------------------------------------
!       Reporting results

ifdef(`USE_PETSC', ``        if (procID.eq.0) then
' define(sp, `   ')', `define(sp, `')')dnl
        sp()if (ierr.eq.0) then
        sp()   print *, 'PRIMME_SVDS has returned successfully'
        sp()else 
        sp()   print *, 'PRIMME_SVDS returned with error: ', ierr
        sp()endif

!       sp()
!       sp()Example of obtaining primme members from the driver:
!       sp()NOTE: don't use primme_svds_get_member_f77, which can only be used in a callback
!
        sp()call primme_svdstop_get_member_f77(primme_svds,
     :  sp()                    PRIMME_SVDS_eps, epsil, ierr)
        sp()call primme_svdstop_get_member_f77(primme_svds,
     :  sp()                    PRIMME_SVDS_aNorm, aNorm, ierr)
        sp()print '(A16,E8.2,A20,e12.5)', 'Tolerance used: ',epsil,
     :  sp()                           '  Estimated norm(A):',aNorm
!
!       sp()Reporting of svals and residuals
!
        sp()do i = 1, NUMSmax
        sp()   write (*, 9000) i, svals(i),rnorms(i)
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
        subroutine generateLauchli(m0,n0,c,A,ierr)
!       ----------------------------------------------------------------
        implicit none
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
        integer*8 m0, n0
        real*8 c
        PetscInt m, n, zeroi
        Mat A

        PetscScalar mu, oned
        PetscInt i,i_1, Istart,Iend
        PetscErrorCode ierr

        call MatCreate(PETSC_COMM_WORLD, A, ierr)
        m = m0
        n = n0
        call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, m, n, ierr)
        call MatSetUp(A, ierr)
        call MatGetOwnershipRange(A, Istart, Iend, ierr)
        zeroi = 0
        oned = 1.0
        if (Istart.eq.0) then
           do i=0,n-1
              call MatSetValue(A, zeroi, i, oned, INSERT_VALUES, ierr)
           enddo
        endif
        do i=max(1,Istart),min(Iend,n)-1
           mu = (1.0 - c)*(i-1)/(min(m,n) - 1)
           i_1 = i-1
           call MatSetValue(A, i, i_1, mu, INSERT_VALUES, ierr)
        enddo
        call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
        end
        subroutine PETScMatvec(x,ldx,y,ldy,k,transpose,primme_svds,err)
!       ----------------------------------------------------------------
ifdef([USE_POINTER], [        use iso_c_binding
])dnl
        implicit none
        include 'primme_f77.h'
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
        integer*8 ldx,ldy
        integer k,transpose
        PRIMME_NUM x(ldx,*), y(ldy,*)
        integer*8 primme_svds
        integer j,err
ifdef([USE_POINTER], [        Mat, pointer :: A
        type(c_ptr) :: pA
], [        Mat A
        PC pc
        COMMON A, pc
])dnl
        Vec xvec,yvec,aux
        PetscErrorCode ierr

ifdef([USE_POINTER], [        call primme_svds_get_member_f77(primme_svds,
     :                                  PRIMME_SVDS_matrix, pA, err)
        call c_f_pointer(pA, A)
])
#if PETSC_VERSION_LT(3,6,0)
        call MatGetVecs(A, xvec, yvec, ierr)
#else
        call MatCreateVecs(A, xvec, yvec, ierr)
#endif
        if (transpose.ne.0) then
           aux = xvec
           xvec = yvec
           yvec = aux
        endif
        do j=1,k
           call VecPlaceArray(xvec, x(1,j), ierr)
           call VecPlaceArray(yvec, y(1,j), ierr)
           if (transpose.eq.0) then
              call MatMult(A, xvec, yvec, ierr)
           else
              call MatMultHermitianTranspose(A, xvec, yvec, ierr)
           endif
           call VecResetArray(xvec, ierr)
           call VecResetArray(yvec, ierr)
        enddo
        call VecDestroy(xvec, ierr)
        call VecDestroy(yvec, ierr)
        err = 0
        end
        subroutine ApplyPCPrecAHA(x,ldx,y,ldy,k,mode,primme_svds,err)
!       ----------------------------------------------------------------
ifdef([USE_POINTER], [        use iso_c_binding
])dnl
        implicit none
        include 'primme_f77.h'
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
        integer*8 ldx,ldy,mLocal,nLocal
        integer k,mode,err
        PRIMME_NUM x(ldx,*), y(ldy,*)
        integer*8 primme_svds
        integer j
ifdef([USE_POINTER], [        Mat, pointer :: A
        PC, pointer :: pc
        type(c_ptr) :: pA, ppc
], [        Mat A
        PC pc
        COMMON A, pc
])dnl
        Vec xvec,yvec
        PetscErrorCode ierr

        call primme_svds_get_member_f77(primme_svds,
     :                              PRIMME_SVDS_mLocal, mLocal, err)
        call primme_svds_get_member_f77(primme_svds,
     :                              PRIMME_SVDS_nLocal, nLocal, err)
ifdef([USE_POINTER], [        call primme_svds_get_member_f77(primme_svds,
     :                                  PRIMME_SVDS_matrix, pA, err)
        call primme_svds_get_member_f77(primme_svds,
     :                         PRIMME_SVDS_preconditioner, ppc, err)
        call c_f_pointer(pA, A)
        call c_f_pointer(ppc, pc)
])
        if (mode.eq.PRIMME_SVDS_op_AtA) then
           call MatCreateVecs(A, xvec, PETSC_NULL_OBJECT, ierr)
           call MatCreateVecs(A, yvec, PETSC_NULL_OBJECT, ierr)
           do j=1,k
              call VecPlaceArray(xvec, x(1,j), ierr)
              call VecPlaceArray(yvec, y(1,j), ierr)
              call PCApply(pc, xvec, yvec, ierr)
              call VecResetArray(xvec, ierr)
              call VecResetArray(yvec, ierr)
           enddo
           call VecDestroy(xvec, ierr)
           call VecDestroy(yvec, ierr)
        else if (mode.eq.PRIMME_SVDS_op_AAt) then
           y(1:mLocal,1:k) = x(1:mLocal,1:k) 
        else if (mode.eq.PRIMME_SVDS_op_augmented) then
           y(1:mLocal+nLocal,1:k) = x(1:mLocal+nLocal,1:k) 
        endif
        err = 0
        end
        subroutine par_GlobalSum(x,y,k,primme_svds,ierr)
!       ----------------------------------------------------------------
        use iso_c_binding
        implicit none
        include 'primme_f77.h'
#include <petsc/finclude/petscsys.h>
        real*8, target :: x(k), y(k)
        integer*8 primme_svds
        integer k,ierr
ifdef([USE_POINTER], [        MPI_Comm, pointer :: comm
        type(c_ptr) :: pcomm

        call primme_svds_get_member_f77(primme_svds,
     :                             PRIMME_SVDS_commInfo, pcomm, ierr)
        call c_f_pointer(pcomm, comm)
])dnl
        if (c_associated(c_loc(x),c_loc(y))) then
          call MPI_Allreduce(MPI_IN_PLACE, y, k, MPIU_REAL, MPIU_SUM,
     :                                 ifdef([USE_POINTER], [comm], [PETSC_COMM_WORLD]), ierr)
        else
          call MPI_Allreduce(x, y, k, MPIU_REAL, MPIU_SUM,
     :                                 ifdef([USE_POINTER], [comm], [PETSC_COMM_WORLD]), ierr)
        endif
        end
], [
!       Lauchli-like block matrix-vector products, Y = A * X or Y = A' * X,
!       where
!      
!       - X, input dense matrix of size primme.n x blockSize or primme.m x blockSize;
!       - Y, output dense matrix of size primme.m x blockSize or primme.n x blockSize;
!       - A, rectangular matrix of size primme.m x primme.n with this form:
!      
!            [ 1  1  1  1  1 ... ],  ei = 1 - (1 - c)*i/(min(m,n) - 1)
!            [e0  0  0  0  0 ... ]
!            [ 0 e1  0  0  0 ... ]
!            [ 0  0 e2  0  0 ... ]
!             ...
!      
        subroutine MV(x,ldx,y,ldy,k,transpose,primme_svds,ierr)
!       ----------------------------------------------------------------
        implicit none
        intrinsic min
        include 'primme_f77.h'
        integer*8 ldx,ldy
        PRIMME_NUM x(ldx,*), y(ldy,*)
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
        PRIMME_NUM x(ldx,*), y(ldy,*)
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
])dnl
