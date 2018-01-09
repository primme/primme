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
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscmat.h>
        use petscvec
        use petscmat
        use petscpc
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
     :            TOL             = 1.0D-5,
     :            printLevel      = 2,
     :            whichSvals      = primme_svds_closest_abs,
     :            numTargetShifts = 1
     :  )

        real*8 TargetShifts(numTargetShifts)
        data TargetShifts /0.5D0/

        external generateLauchli, PETScMatvec, ApplyPCPrecAHA,
     :           par_GlobalSum

!       Singular values, vectors and their residual norms
!
        PetscReal svals(NUMSmax)
        PetscReal rnorms(NUMSmax)
        PetscScalar svecs((m+n)*NUMSmax)

!       Other vars
!
        Mat, target :: A
        PC, target :: pc
        MPI_Comm, target :: comm
        Mat AH, AHA
        PetscErrorCode ierr
        PetscInt mLocalpetsc,nLocalpetsc
        integer i,numProcsMpi,procIDMpi
        integer*8 numProcs,procID,mLocal,nLocal
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
        call generateLauchli(m, n, c, A, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                                   PRIMME_SVDS_matrix, A, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                   PRIMME_SVDS_matrixMatvec, PETScMatvec, ierr)
        
!       Set parallel parameters
        call MatGetLocalSize(A, mLocalpetsc, nLocalpetsc, ierr)
        mLocal = mLocalpetsc
        nLocal = nLocalpetsc
        call primme_svds_set_member_f77(primme_svds,
     :                              PRIMME_SVDS_mLocal, mLocal, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :                              PRIMME_SVDS_nLocal, nLocal, ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD, numProcsMpi, ierr)
        numProcs = numProcsMpi
        call primme_svds_set_member_f77(primme_svds,
     :                         PRIMME_SVDS_numProcs, numProcs, ierr)
        call MPI_Comm_rank(PETSC_COMM_WORLD, procIDMpi, ierr)
        procID = procIDMpi
        call primme_svds_set_member_f77(primme_svds,
     :                              PRIMME_SVDS_procID, procID, ierr)
        comm = PETSC_COMM_WORLD
        call primme_svds_set_member_f77(primme_svds,
     :                              PRIMME_SVDS_commInfo, comm, ierr)
        call primme_svds_set_member_f77(primme_svds,
     :        PRIMME_SVDS_globalSumReal, par_GlobalSum, ierr)

!       Set preconditioner based on A^H*A (optional)
        pc = tPC(0)
        call PCCreate(PETSC_COMM_WORLD, pc, ierr)
        call MatHermitianTranspose(A, MAT_INITIAL_MATRIX, AH, ierr)
        call MatMatMult(AH, A, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL,
     :                                                        AHA, ierr)
        call PCSetType(pc, PCJACOBI, ierr)
        call PCSetOperators(pc, AHA, AHA, ierr)
        call PCSetFromOptions(pc, ierr)
        call PCSetUp(pc, ierr)
        call primme_svds_set_member_f77(primme_svds, 
     :       PRIMME_SVDS_preconditioner, pc, ierr)
        call primme_svds_set_member_f77(primme_svds, 
     :       PRIMME_SVDS_applyPreconditioner, ApplyPCPrecAHA, ierr)
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

        if (procID.eq.0)call primme_svds_display_params_f77(primme_svds)

!       ----------------------------------------------------------------
!       Calling the PRIMME solver
!       ----------------------------------------------------------------

#if defined(PETSC_USE_COMPLEX) && defined(PETSC_USE_REAL_SINGLE)
        call cprimme_svds_f77(svals, svecs, rnorms, primme_svds, ierr)
#elif defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_REAL_SINGLE)
        call zprimme_svds_f77(svals, svecs, rnorms, primme_svds, ierr)
#elif !defined(PETSC_USE_COMPLEX) && defined(PETSC_USE_REAL_SINGLE)
        call sprimme_svds_f77(svals, svecs, rnorms, primme_svds, ierr)
#elif !defined(PETSC_USE_COMPLEX) && !defined(PETSC_USE_REAL_SINGLE)
        call dprimme_svds_f77(svals, svecs, rnorms, primme_svds, ierr)
#endif
        if (ierr.ne.0) then
          stop 1
        endif

!       ----------------------------------------------------------------
!       Reporting results

        if (procID.eq.0) then
            if (ierr.eq.0) then
              print *, 'PRIMME_SVDS has returned successfully'
           else 
              print *, 'PRIMME_SVDS returned with error: ', ierr
           endif

!          
!          Example of obtaining primme members from the driver:
!          NOTE: don't use primme_svds_get_member_f77, which can only be used in a callback
!
           call primme_svdstop_get_member_f77(primme_svds,
     :                         PRIMME_SVDS_eps, epsil, ierr)
           call primme_svdstop_get_member_f77(primme_svds,
     :                         PRIMME_SVDS_aNorm, aNorm, ierr)
           print '(A16,E8.2,A20,e12.5)', 'Tolerance used: ',epsil,
     :                                '  Estimated norm(A):',aNorm
!
!          Reporting of svals and residuals
!
           do i = 1, NUMSmax
              write (*, 9000) i, svals(i),rnorms(i)
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


        subroutine generateLauchli(m0,n0,c,A,ierr)
!       ----------------------------------------------------------------
        use petscmat
        implicit none
        integer*8 m0, n0
        real*8 c
        PetscInt m, n, zeroi
        Mat A

        PetscScalar mu, oned
        PetscInt i,i_1, Istart,Iend
        PetscErrorCode ierr

        A = tMat(0)
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
        use iso_c_binding
        use petscvec
        use petscmat
        use petscpc
        implicit none
        include 'primme_f77.h'
        integer*8 ldx,ldy
        integer k,transpose
        PetscScalar x(ldx,*), y(ldy,*)
        integer*8 primme_svds
        integer j,err
        Mat, pointer :: A
        type(c_ptr) :: pA
        Vec xvec,yvec,aux
        PetscErrorCode ierr

        call primme_svds_get_member_f77(primme_svds,
     :                                  PRIMME_SVDS_matrix, pA, err)
        call c_f_pointer(pA, A)

        xvec = tVec(0)
        yvec = tVec(0)
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
        use iso_c_binding
        use petscvec
        use petscmat
        use petscpc
        implicit none
        include 'primme_f77.h'
        integer*8 ldx,ldy,mLocal,nLocal
        integer k,mode,err
        PetscScalar x(ldx,*), y(ldy,*)
        integer*8 primme_svds
        integer j
        Mat, pointer :: A
        PC, pointer :: pc
        type(c_ptr) :: pA, ppc
        Vec xvec,yvec
        PetscErrorCode ierr

        call primme_svds_get_member_f77(primme_svds,
     :                              PRIMME_SVDS_mLocal, mLocal, err)
        call primme_svds_get_member_f77(primme_svds,
     :                              PRIMME_SVDS_nLocal, nLocal, err)
        call primme_svds_get_member_f77(primme_svds,
     :                                  PRIMME_SVDS_matrix, pA, err)
        call primme_svds_get_member_f77(primme_svds,
     :                         PRIMME_SVDS_preconditioner, ppc, err)
        call c_f_pointer(pA, A)
        call c_f_pointer(ppc, pc)

        if (mode.eq.PRIMME_SVDS_op_AtA) then
           xvec = tVec(0)
           yvec = tVec(0)
           call MatCreateVecs(A, xvec, PETSC_NULL_VEC, ierr)
           call MatCreateVecs(A, yvec, PETSC_NULL_VEC, ierr)
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
        use petscsys
        implicit none
        include 'primme_f77.h'
        PetscReal, target :: x(k), y(k)
        integer*8 primme_svds
        integer k,ierr
        MPI_Comm, pointer :: comm
        type(c_ptr) :: pcomm

        call primme_svds_get_member_f77(primme_svds,
     :                             PRIMME_SVDS_commInfo, pcomm, ierr)
        call c_f_pointer(pcomm, comm)
        if (c_associated(c_loc(x),c_loc(y))) then
          call MPI_Allreduce(MPI_IN_PLACE, y, k, MPIU_REAL, MPIU_SUM,
     :                                 comm, ierr)
        else
          call MPI_Allreduce(x, y, k, MPIU_REAL, MPIU_SUM,
     :                                 comm, ierr)
        endif
        end
