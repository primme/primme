/*******************************************************************************
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2017 College of William & Mary,
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
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *******************************************************************************
 * File: petscw.c
 * 
 * Purpose - PETSc wrapper.
 * 
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include "mmio.h"
#include <petscpc.h>
#include <petscmat.h>
#include <petsc/private/matimpl.h>
#include "primme_svds.h"
#include "petscw.h"

static PetscErrorCode preallocation(Mat M,PetscInt *d_nz, PetscInt *o_nz);
static PetscErrorCode loadmtx(const char* filename, Mat *M, PetscBool *pattern);
static PetscErrorCode permutematrix(Mat Ain, Mat Bin, Mat *Aout, Mat *Bout, int **permIndices);

#undef __FUNCT__
#define __FUNCT__ "readMatrixPetsc"
int readMatrixPetsc(const char* matrixFileName, PRIMME_INT *m, PRIMME_INT *n,
      PRIMME_INT *mLocal, PRIMME_INT *nLocal, int *numProcs, int *procID,
      Mat **matrix, double *fnorm_, int **perm) {

   PetscErrorCode ierr;
   PetscReal fnorm;
   PetscBool pattern;
   PetscViewer viewer;
   PetscInt m0, n0, mLocal0, nLocal0;

   PetscFunctionBegin;

   *matrix = (Mat *)primme_calloc(1, sizeof(Mat), "mat");
   if (!strcmp("mtx", &matrixFileName[strlen(matrixFileName)-3])) {  
      // coordinate format storing both lower and upper triangular parts
      ierr = loadmtx(matrixFileName, *matrix, &pattern); CHKERRQ(ierr);
   }
   else if (!strcmp("petsc", &matrixFileName[strlen(matrixFileName)-5])) {
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, matrixFileName, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
      ierr = MatCreate(PETSC_COMM_WORLD, *matrix); CHKERRQ(ierr);
      ierr = MatSetFromOptions(**matrix); CHKERRQ(ierr);
      ierr = MatLoad(**matrix, viewer); CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
   }
   else {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Could not read matrix file.");
   }
   if (fnorm_) {
      ierr = MatNorm(**matrix, NORM_FROBENIUS, &fnorm); CHKERRQ(ierr);
      *fnorm_ = fnorm;
   }

   ierr = MatGetSize(**matrix, &m0, &n0); CHKERRQ(ierr);
   *m = m0;
   *n = n0;

   if (perm && *m == *n) {
      Mat Atemp;
      ierr = permutematrix(**matrix, NULL, &Atemp, NULL, perm);CHKERRQ(ierr);
      ierr = MatDestroy(*matrix);CHKERRQ(ierr);
      **matrix = Atemp;
   }
   else if (perm) {
      *perm = NULL;
   }

   ierr = MatGetLocalSize(**matrix, &mLocal0, &nLocal0); CHKERRQ(ierr);
   *mLocal = mLocal0;
   *nLocal = nLocal0;
   MPI_Comm_size(MPI_COMM_WORLD, numProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, procID);

   PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "loadmtx"
static PetscErrorCode loadmtx(const char* filename, Mat *M, PetscBool *pattern) {
   PetscErrorCode ierr;
   FILE        *f;
   MM_typecode type;
   int         m,n,nz,i,j,k;
   PetscInt    low,high,lowj,highj,*d_nz,*o_nz;
   double      re,im;
   PetscScalar s;
   long        pos;

   PetscFunctionBegin;
   
   f = fopen(filename,"r");
   if (!f) SETERRQ2(PETSC_COMM_SELF,1,"fopen '%s': %s",filename,strerror(errno));
   
   /* first read to set matrix kind and size */
   ierr = mm_read_banner(f,&type);CHKERRQ(ierr);
   if (!mm_is_valid(type) || !mm_is_sparse(type) ||
       !(mm_is_real(type) || mm_is_complex(type) || mm_is_pattern(type) || mm_is_integer(type)))
      SETERRQ1(PETSC_COMM_SELF,1,"Matrix format '%s' not supported",mm_typecode_to_str(type)); 
#if !defined(PETSC_USE_COMPLEX)
   if (mm_is_complex(type)) SETERRQ(PETSC_COMM_SELF,1,"Complex matrix not supported in real configuration"); 
#endif
   if (pattern) *pattern = mm_is_pattern(type) ? PETSC_TRUE : PETSC_FALSE;
  
   ierr = mm_read_mtx_crd_size(f,&m,&n,&nz);CHKERRQ(ierr);
   pos = ftell(f);
   ierr = MatCreate(PETSC_COMM_WORLD,M);CHKERRQ(ierr);
   ierr = MatSetSizes(*M,PETSC_DECIDE,PETSC_DECIDE,(PetscInt)m,(PetscInt)n);CHKERRQ(ierr);
   ierr = MatSetFromOptions(*M);CHKERRQ(ierr);
   ierr = MatSetUp(*M);CHKERRQ(ierr);

   ierr = MatGetOwnershipRange(*M,&low,&high);CHKERRQ(ierr);  
   ierr = MatGetOwnershipRangeColumn(*M,&lowj,&highj);CHKERRQ(ierr);  
   ierr = PetscMalloc(sizeof(PetscInt)*(high-low),&d_nz);CHKERRQ(ierr);
   ierr = PetscMalloc(sizeof(PetscInt)*(high-low),&o_nz);CHKERRQ(ierr);
   for (i=0; i<high-low;i++) {
      d_nz[i] = (i+low>=lowj && i+low<highj) ? 1 : 0;
      o_nz[i] = (i+low>=lowj && i+low<highj) ? 0 : 1;
   }
   for (k=0;k<nz;k++) {
      ierr = mm_read_mtx_crd_entry(f,&i,&j,&re,&im,type);CHKERRQ(ierr);
      i--; j--;
      if (i!=j) {
         if (i>=low && i<high) {
            if (j>=lowj && j<highj) 
               d_nz[i-low]++;
            else
               o_nz[i-low]++;
         }
         if (j>=low && j<high && !mm_is_general(type)) {
            if (i>=low && i<high) 
               d_nz[j-low]++;
            else
               o_nz[j-low]++;        
         }
      }
   }
   ierr = preallocation(*M,d_nz,o_nz);CHKERRQ(ierr);
   ierr = PetscFree(d_nz);CHKERRQ(ierr);
   ierr = PetscFree(o_nz);CHKERRQ(ierr);
  
   /* second read to load the values */ 
   ierr = fseek(f, pos, SEEK_SET);
   if (ierr) SETERRQ1(PETSC_COMM_SELF,1,"fseek: %s",strerror(errno));
    
   re = 1.0;
   im = 0.0;
   /* Set the diagonal to zero */
   for (i=low; i<PetscMin(high,n); i++) {
      ierr = MatSetValue(*M,i,i,0.0,INSERT_VALUES);CHKERRQ(ierr);
   }
   for (k=0;k<nz;k++) {
      ierr = mm_read_mtx_crd_entry(f,&i,&j,&re,&im,type);
      i--; j--;
      if (i>=low && i<high) {
         s = re + IMAGINARY * im;
         ierr = MatSetValue(*M,i,j,s,INSERT_VALUES);CHKERRQ(ierr);
      }
      if (j>=low && j<high && i != j && !mm_is_general(type)) {
         if (mm_is_symmetric(type)) s = re + IMAGINARY * im;
         else if (mm_is_hermitian(type)) s = re - IMAGINARY * im;
         else if (mm_is_skew(type)) s = -re - IMAGINARY * im;
         else {
            SETERRQ1(PETSC_COMM_SELF,1,"Matrix format '%s' not supported",mm_typecode_to_str(type));
         }
         ierr = MatSetValue(*M,j,i,s,INSERT_VALUES);CHKERRQ(ierr);
      }
   }
   ierr = MatAssemblyBegin(*M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
   ierr = MatAssemblyEnd(*M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

   if (mm_is_symmetric(type)) { 
      ierr = MatSetOption(*M,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
   }
   if ((mm_is_symmetric(type) && mm_is_real(type)) || mm_is_hermitian(type)) { 
      ierr = MatSetOption(*M,MAT_HERMITIAN,PETSC_TRUE);CHKERRQ(ierr);
   }

   ierr = fclose(f);
   if (ierr) SETERRQ1(PETSC_COMM_SELF,1,"fclose: %s",strerror(errno));

   PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "preallocation"
static PetscErrorCode preallocation(Mat M,PetscInt *d_nz, PetscInt *o_nz) {
   PetscErrorCode ierr;
   PetscBool      isaij,ismpiaij,isseqaij;
   PetscMPIInt    size;

   PetscFunctionBegin;

   ierr = PetscObjectTypeCompare((PetscObject)M,MATAIJ,&isaij);CHKERRQ(ierr);
   ierr = PetscObjectTypeCompare((PetscObject)M,MATMPIAIJ,&ismpiaij);CHKERRQ(ierr);
   ierr = PetscObjectTypeCompare((PetscObject)M,MATSEQAIJ,&isseqaij);CHKERRQ(ierr);
   ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

   if ((isaij && size == 1) || isseqaij) {
      ierr = MatSeqAIJSetPreallocation(M,0,d_nz);CHKERRQ(ierr);
   } else if (isaij || ismpiaij) {
      ierr = MatMPIAIJSetPreallocation(M,0,d_nz,0,o_nz);CHKERRQ(ierr);
   } else {
      ierr = PetscInfo(M,"NOT using preallocation\n");CHKERRQ(ierr);
   }
  
   PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "permutematrix"
static PetscErrorCode permutematrix(Mat Ain, Mat Bin, Mat *Aout, Mat *Bout, int **permIndices)
{
   PetscErrorCode  ierr;
   MatPartitioning part;
   IS              isn, is, iscols;
   PetscInt        *nlocal,localCols,m,n;
   PetscMPIInt     size, rank;
   MPI_Comm        comm;
 
   PetscFunctionBegin;
 
   ierr = PetscObjectGetComm((PetscObject)Ain,&comm);CHKERRQ(ierr);
   ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
   ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
   ierr = MatGetSize(Ain,&m,&n);CHKERRQ(ierr);
   ierr = MatPartitioningCreate(comm,&part);CHKERRQ(ierr);
   ierr = MatPartitioningSetAdjacency(part,Ain);CHKERRQ(ierr);
   ierr = MatPartitioningSetFromOptions(part);CHKERRQ(ierr);
   /* get new processor owner number of each vertex */
   ierr = MatPartitioningApply(part,&is);CHKERRQ(ierr);
   /* get new global number of each old global number */
   ierr = ISPartitioningToNumbering(is,&isn);CHKERRQ(ierr);
   ierr = PetscMalloc(size*sizeof(int),&nlocal);CHKERRQ(ierr);
   /* get number of new vertices for each processor */
   ierr = ISPartitioningCount(is,size,nlocal);CHKERRQ(ierr);
   ierr = ISDestroy(&is);CHKERRQ(ierr);
 
   /* get old global number of each new global number */
   ierr = ISInvertPermutation(isn,nlocal[rank],&is);CHKERRQ(ierr);
   ierr = ISDestroy(&isn);CHKERRQ(ierr);
   ierr = MatPartitioningDestroy(&part);CHKERRQ(ierr);
   ierr = ISSort(is);CHKERRQ(ierr);

   /* If matrix is square, the permutation is applied to rows and columns;
      otherwise it is only applied to rows. */
   if (m == n) {
      iscols = is;
      localCols = nlocal[rank];
   } else {
      PetscInt lowj, highj;
      ierr = MatGetOwnershipRangeColumn(Ain,&lowj,&highj);CHKERRQ(ierr);  
      localCols = highj-lowj;
      ierr = ISCreateStride(comm,localCols, lowj, 1, &iscols);CHKERRQ(ierr);
   }

   /* copy permutation */
   if (permIndices) {
      const PetscInt *indices;
      PetscInt i;
      *permIndices = malloc(sizeof(int)*(nlocal[rank]+localCols));
      ierr = ISGetIndices(is, &indices);CHKERRQ(ierr);
      for (i=0; i<nlocal[rank]; i++) (*permIndices)[i] = indices[i];
      ierr = ISRestoreIndices(is, &indices);CHKERRQ(ierr);
      ierr = ISGetIndices(iscols, &indices);CHKERRQ(ierr);
      for (i=0; i<localCols; i++) (*permIndices)[i+nlocal[rank]] = indices[i];
      ierr = ISRestoreIndices(iscols, &indices);CHKERRQ(ierr);
   }
 
   ierr = PetscFree(nlocal);CHKERRQ(ierr);

   ierr = MatGetSubMatrix(Ain,is,iscols,MAT_INITIAL_MATRIX,Aout);CHKERRQ(ierr);
   if (Bin && Bout) {
      ierr = MatGetSubMatrix(Bin,is,iscols,MAT_INITIAL_MATRIX,Bout);CHKERRQ(ierr);
   }
   ierr = ISDestroy(&is);CHKERRQ(ierr);
   if (m != n) {
      ierr = ISDestroy(&iscols);CHKERRQ(ierr);
   }
 
   PetscFunctionReturn(0);
}

/******************************************************************************/
/* Matvecs column by column                                                   */
/******************************************************************************/

static void PETScMatvecGenNoBlock(void *x, PRIMME_INT ldx, void *y, PRIMME_INT ldy,
      int blockSize, int trans, Mat matrix, MPI_Comm comm) {
   int i;
   Vec xvec, yvec;
   PetscInt m, n, mLocal, nLocal;
   PetscErrorCode ierr;
   
   assert(sizeof(PetscScalar) == sizeof(SCALAR));   
   ierr = MatGetSize(matrix, &m, &n); CHKERRABORT(comm, ierr);
   ierr = MatGetLocalSize(matrix, &mLocal, &nLocal); CHKERRABORT(comm, ierr);

   #if PETSC_VERSION_LT(3,6,0)
      ierr = MatGetVecs(matrix, &xvec, &yvec); CHKERRABORT(comm, ierr);
   #else
      ierr = MatCreateVecs(matrix, &xvec, &yvec); CHKERRABORT(comm, ierr);
   #endif
   if (trans == 1) {
      Vec aux = xvec; xvec = yvec; yvec = aux;
   }
   for (i=0; i<blockSize; i++) {
      ierr = VecPlaceArray(xvec, ((SCALAR*)x) + ldx*i); CHKERRABORT(comm, ierr);
      ierr = VecPlaceArray(yvec, ((SCALAR*)y) + ldy*i); CHKERRABORT(comm, ierr);
      if (trans == 0) {
         ierr = MatMult(matrix, xvec, yvec); CHKERRABORT(comm, ierr);
      } else {
         ierr = MatMultHermitianTranspose(matrix, xvec, yvec); CHKERRABORT(comm, ierr);
      }
      ierr = VecResetArray(xvec); CHKERRABORT(comm, ierr);
      ierr = VecResetArray(yvec); CHKERRABORT(comm, ierr);
   }
   ierr = VecDestroy(&xvec); CHKERRABORT(comm, ierr);
   ierr = VecDestroy(&yvec); CHKERRABORT(comm, ierr);
}

/******************************************************************************/
/* Matvecs blocked and row-major (not useful now)                             */
/******************************************************************************/

#undef __FUNCT__
#define __FUNCT__ "MatMatMult_MPIAIJ_MPIDense0"
static PetscErrorCode MatMatMult_MPIAIJ_MPIDense0(Mat A,PetscScalar *Barray, PetscInt BN, PRIMME_INT ldb,PetscScalar *Carray, PetscInt ldc)
{
  PetscErrorCode      ierr;
  PetscInt            i,j,m=A->rmap->n,n=A->cmap->n;
  PetscScalar         *btarray,*ctarray;
  Vec                 bt,ct;
  Mat                 mA;

  PetscFunctionBegin;
  /* create MAIJ matrix mA from A -- should be done in symbolic phase */
  ierr = MatCreateMAIJ(A,BN,&mA);CHKERRQ(ierr);

  /* create vectors bt and ct to hold locally transposed arrays of B and C */
  ierr = VecCreate(PetscObjectComm((PetscObject)A),&bt);CHKERRQ(ierr);
  ierr = VecSetSizes(bt,n*BN,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetType(bt,VECSTANDARD);CHKERRQ(ierr);
  ierr = VecCreate(PetscObjectComm((PetscObject)A),&ct);CHKERRQ(ierr);
  ierr = VecSetSizes(ct,m*BN,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetType(ct,VECSTANDARD);CHKERRQ(ierr);

  /* transpose local arry of B, then copy it to vector bt */
  ierr = VecGetArray(bt,&btarray);CHKERRQ(ierr);

  for (j=0; j<BN; j++) {
    for (i=0; i<n; i++) btarray[i*BN + j] = Barray[ldb*j+i];
  }
  ierr = VecRestoreArray(bt,&btarray);CHKERRQ(ierr);

  /* compute ct = mA^T * cb */
  ierr = MatMult(mA,bt,ct);CHKERRQ(ierr);

  /* transpose local arry of ct to matrix C */
  ierr = VecGetArray(ct,&ctarray);CHKERRQ(ierr);
  for (j=0; j<BN; j++) {
    for (i=0; i<m; i++) Carray[j*ldc+i] = ctarray[i*BN + j];
  }
  ierr = VecRestoreArray(ct,&ctarray);CHKERRQ(ierr);

  ierr = VecDestroy(&bt);CHKERRQ(ierr);
  ierr = VecDestroy(&ct);CHKERRQ(ierr);
  ierr = MatDestroy(&mA);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatHermitianTransposeMatMult_MPIAIJ_MPIDense0"
static PetscErrorCode MatHermitianTransposeMatMult_MPIAIJ_MPIDense0(Mat A,PetscScalar *Barray, PetscInt BN, PRIMME_INT ldb,PetscScalar *Carray, PetscInt ldc)
{
  PetscErrorCode      ierr;
  PetscInt            i,j,m=A->rmap->n,n=A->cmap->n;
  PetscScalar         *btarray,*ctarray;
  Vec                 bt,ct;
  Mat                 mA;

  PetscFunctionBegin;
  /* create MAIJ matrix mA from A -- should be done in symbolic phase */
  ierr = MatCreateMAIJ(A,BN,&mA);CHKERRQ(ierr);

  /* create vectors bt and ct to hold locally transposed arrays of B and C */
  ierr = VecCreate(PetscObjectComm((PetscObject)A),&bt);CHKERRQ(ierr);
  ierr = VecSetSizes(bt,m*BN,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetType(bt,VECSTANDARD);CHKERRQ(ierr);
  ierr = VecCreate(PetscObjectComm((PetscObject)A),&ct);CHKERRQ(ierr);
  ierr = VecSetSizes(ct,n*BN,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetType(ct,VECSTANDARD);CHKERRQ(ierr);

  /* transpose local arry of B, then copy it to vector bt */
  ierr = VecGetArray(bt,&btarray);CHKERRQ(ierr);

  for (j=0; j<BN; j++) {
    for (i=0; i<m; i++) btarray[i*BN + j] = PetscConj(Barray[ldb*j+i]);
  }
  ierr = VecRestoreArray(bt,&btarray);CHKERRQ(ierr);

  /* compute ct = mA^T * cb */
  ierr = MatMultTranspose(mA,bt,ct);CHKERRQ(ierr);

  /* transpose local arry of ct to matrix C */
  ierr = VecGetArray(ct,&ctarray);CHKERRQ(ierr);
  for (j=0; j<BN; j++) {
    for (i=0; i<n; i++) Carray[j*ldc+i] = PetscConj(ctarray[i*BN + j]);
  }
  ierr = VecRestoreArray(ct,&ctarray);CHKERRQ(ierr);

  ierr = VecDestroy(&bt);CHKERRQ(ierr);
  ierr = VecDestroy(&ct);CHKERRQ(ierr);
  ierr = MatDestroy(&mA);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static void PETScMatvecGenRowMajor(void *x, PRIMME_INT ldx, void *y, PRIMME_INT ldy,
      int blockSize, int trans, Mat matrix, MPI_Comm comm) {
   PetscInt m, n, mLocal, nLocal;
   PetscErrorCode ierr;
   Mat X, Y;
  
   if (blockSize == 1) { 
      PETScMatvecGenNoBlock(x, ldx, y, ldy, blockSize, trans, matrix, comm);
      return;
   }

   assert(sizeof(PetscScalar) == sizeof(SCALAR));   
   ierr = MatGetSize(matrix, &m, &n); CHKERRABORT(comm, ierr);
   ierr = MatGetLocalSize(matrix, &mLocal, &nLocal); CHKERRABORT(comm, ierr);

   if (trans == 0) {
      ierr = MatMatMult_MPIAIJ_MPIDense0(matrix,(PetscScalar*)x,blockSize,ldx,(PetscScalar*)y,ldy);CHKERRABORT(comm, ierr);
   }
   else {
      ierr = MatHermitianTransposeMatMult_MPIAIJ_MPIDense0(matrix,(PetscScalar*)x,blockSize,ldx,(PetscScalar*)y,ldy);CHKERRABORT(comm, ierr);
   }
}

/******************************************************************************/
/* Matvecs blocked and column-major                                           */
/******************************************************************************/

static void PETScMatvecGenColumnMajor(void *x, PRIMME_INT ldx, void *y, PRIMME_INT ldy,
      int blockSize, int trans, Mat matrix, MPI_Comm comm) {
   PetscInt m, n, mLocal, nLocal;
   PetscErrorCode ierr;
   Mat X, Y, X0, Y0;
   int xcompact, ycompact;
  
   if (blockSize == 1) { 
      PETScMatvecGenNoBlock(x, ldx, y, ldy, blockSize, trans, matrix, comm);
      return;
   }

   assert(sizeof(PetscScalar) == sizeof(SCALAR));   
   ierr = MatGetSize(matrix, &m, &n); CHKERRABORT(comm, ierr);
   ierr = MatGetLocalSize(matrix, &mLocal, &nLocal); CHKERRABORT(comm, ierr);

   if (trans == 0) {
      ierr = MatCreateDense(comm,nLocal,PETSC_DECIDE,n,blockSize,x,&X);CHKERRABORT(comm, ierr);
      ierr = MatCreateDense(comm,mLocal,PETSC_DECIDE,m,blockSize,y,&Y);CHKERRABORT(comm, ierr);
      xcompact = nLocal == ldx;
      ycompact = mLocal == ldy;
   }
   else {
      ierr = MatCreateDense(comm,mLocal,PETSC_DECIDE,m,blockSize,x,&X);CHKERRABORT(comm, ierr);
      ierr = MatCreateDense(comm,nLocal,PETSC_DECIDE,n,blockSize,y,&Y);CHKERRABORT(comm, ierr);
      xcompact = mLocal == ldx;
      ycompact = nLocal == ldy;
   }
   ierr = MatDenseGetLocalMatrix(X, &X0);CHKERRABORT(comm, ierr);
   ierr = MatSeqDenseSetLDA(X0, (PetscInt)ldx);CHKERRABORT(comm, ierr);
   ierr = MatDenseGetLocalMatrix(Y, &Y0);CHKERRABORT(comm, ierr);
   ierr = MatSeqDenseSetLDA(Y0, (PetscInt)ldy);CHKERRABORT(comm, ierr);

   /* MatMatMult doesn't support X to be non-contiguous */
   if (xcompact) {
      X0 = X;
   }
   else {
      ierr = MatDuplicate(X, MAT_COPY_VALUES, &X0);CHKERRABORT(comm, ierr);
   }
   if (trans == 0) {
      if (ycompact) {
         ierr = MatMatMult(matrix, X0, MAT_REUSE_MATRIX, PETSC_DEFAULT, &Y); CHKERRABORT(comm, ierr);
      }
      else {
         ierr = MatMatMult(matrix, X0, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Y0); CHKERRABORT(comm, ierr);
         ierr = MatCopy(Y0, Y, SAME_NONZERO_PATTERN); CHKERRABORT(comm, ierr);
         ierr = MatDestroy(&Y0); CHKERRABORT(comm, ierr);
      }
   }
   else {
      /* A^H*X is not implemented in PETSc, do instead (A^T*X^c)^c */
      ierr = MatConjugate(X0); CHKERRABORT(comm, ierr);
      if (ycompact) {
         ierr = MatTransposeMatMult(matrix, X0, MAT_REUSE_MATRIX, PETSC_DEFAULT, &Y); CHKERRABORT(comm, ierr);
      }
      else {
         ierr = MatTransposeMatMult(matrix, X0, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Y0); CHKERRABORT(comm, ierr);
         ierr = MatCopy(Y0, Y, SAME_NONZERO_PATTERN); CHKERRABORT(comm, ierr);
         ierr = MatDestroy(&Y0); CHKERRABORT(comm, ierr);
      }
      ierr = MatConjugate(Y); CHKERRABORT(comm, ierr);
      if (xcompact) {
         ierr = MatConjugate(X0); CHKERRABORT(comm, ierr);
      }
   }
   if (!xcompact) {
      ierr = MatDestroy(&X0);CHKERRABORT(comm, ierr);
   }
   ierr = MatDestroy(&X);CHKERRABORT(comm, ierr);
   ierr = MatDestroy(&Y);CHKERRABORT(comm, ierr);
}

static void PETScMatvecGen(void *x, PRIMME_INT ldx, void *y, PRIMME_INT ldy,
      int blockSize, int trans, Mat matrix, MPI_Comm comm) {
#ifdef PETSC_MATVEC_BLOCK_ROW_MAJOR
   PETScMatvecGenRowMajor (x, ldx, y, ldy, blockSize, trans, matrix, comm);
#elif defined(PETSC_MATVEC_BLOCK_COLUMN_MAJOR)
   PETScMatvecGenColumnMajor(x, ldx, y, ldy, blockSize, trans, matrix, comm);
#else
   PETScMatvecGenNoBlock(x, ldx, y, ldy, blockSize, trans, matrix, comm);
#endif
}

void PETScMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   Mat *matrix;
   PetscInt m, n, mLocal, nLocal;
   PetscErrorCode ierr;

   matrix = (Mat *)primme->matrix;
   ierr = MatGetSize(*matrix, &m, &n);  CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   ierr = MatGetLocalSize(*matrix, &mLocal, &nLocal);  CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   assert(m == primme->n && n == primme->n && mLocal == primme->nLocal
         && nLocal == primme->nLocal);
   PETScMatvecGen(x, *ldx, y, *ldy, *blockSize, 0, *matrix, *(MPI_Comm*)primme->commInfo);
   *err = 0;
}


void PETScMatvecSVD(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, int *trans, primme_svds_params *primme_svds, int *err) {
   PetscInt m, n, mLocal, nLocal;
   Mat *matrix;
   PetscErrorCode ierr;

   matrix = (Mat *)primme_svds->matrix;
   ierr = MatGetSize(*matrix, &m, &n); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
   ierr = MatGetLocalSize(*matrix, &mLocal, &nLocal); CHKERRABORT(*(MPI_Comm*)primme_svds->commInfo, ierr);
   assert(m == primme_svds->m && n == primme_svds->n && mLocal == primme_svds->mLocal
         && nLocal == primme_svds->nLocal);
   PETScMatvecGen(x, *ldx, y, *ldy, *blockSize, *trans, *matrix, *(MPI_Comm*)primme_svds->commInfo);
   *err = 0;
}
 
static void ApplyPCPrecPETSCGen(void *x, PRIMME_INT *ldx, void *y,
      PRIMME_INT *ldy, int *blockSize, int trans, PC *pc, MPI_Comm comm) {
   int i;
   Vec xvec, yvec;
   Mat matrix;
   PetscErrorCode ierr;
   PetscInt mLocal, nLocal;
   
   ierr = PCGetOperators(pc[0],&matrix,NULL); CHKERRABORT(comm, ierr);

   assert(sizeof(PetscScalar) == sizeof(SCALAR));   
   ierr = MatGetLocalSize(matrix, &mLocal, &nLocal); CHKERRABORT(comm, ierr);
   assert(mLocal == nLocal && nLocal <= *ldx && mLocal <= *ldy);

   #if PETSC_VERSION_LT(3,6,0)
      ierr = MatGetVecs(matrix, &xvec, &yvec); CHKERRABORT(comm, ierr);
   #else
      ierr = MatCreateVecs(matrix, &xvec, &yvec); CHKERRABORT(comm, ierr);
   #endif
   for (i=0; i<*blockSize; i++) {
      ierr = VecPlaceArray(xvec, ((SCALAR*)x) + (*ldx)*i); CHKERRABORT(comm, ierr);
      ierr = VecPlaceArray(yvec, ((SCALAR*)y) + (*ldy)*i); CHKERRABORT(comm, ierr);
      if (trans == 0) {
         ierr = PCApply(*pc, xvec, yvec); CHKERRABORT(comm, ierr);
      } else if (pc[1]) {
         ierr = PCApply(pc[1], xvec, yvec); CHKERRABORT(comm, ierr);
      } else {
         ierr = PCApplyTranspose(pc[0], xvec, yvec);
      }
      ierr = VecResetArray(xvec); CHKERRABORT(comm, ierr);
      ierr = VecResetArray(yvec); CHKERRABORT(comm, ierr);
   }
   ierr = VecDestroy(&xvec); CHKERRABORT(comm, ierr);
   ierr = VecDestroy(&yvec); CHKERRABORT(comm, ierr);
}


void ApplyPCPrecPETSC(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, primme_params *primme, int *ierr) {
   ApplyPCPrecPETSCGen(x, ldx, y, ldy, blockSize, 0,
         (PC*)primme->preconditioner, *(MPI_Comm*)primme->commInfo);
   *ierr = 0;
}

void ApplyPCPrecPETSCSVD(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, int *mode, primme_svds_params *primme_svds, int *ierr) {
   int i, one=1;
   SCALAR *aux;

   if (*mode == primme_svds_op_AtA) {
      aux = (SCALAR *)primme_calloc(primme_svds->mLocal, sizeof(SCALAR), "aux");
      for(i=0; i<*blockSize; i++) {
         ApplyPCPrecPETSCGen((SCALAR*)x+(*ldx)*i, ldx, aux, ldy, &one, 1,
            primme_svds->preconditioner, *(MPI_Comm*)primme_svds->commInfo);
         ApplyPCPrecPETSCGen(aux, ldx, (SCALAR*)y+(*ldy)*i, ldy, &one, 0,
            primme_svds->preconditioner, *(MPI_Comm*)primme_svds->commInfo);
      }
      free(aux);
   }
   else if (*mode == primme_svds_op_AAt) {
      aux = (SCALAR *)primme_calloc(primme_svds->nLocal, sizeof(SCALAR), "aux");
      for(i=0; i<*blockSize; i++) {
         ApplyPCPrecPETSCGen((SCALAR*)x+(*ldx)*i, ldx, aux, ldy, &one, 0,
            primme_svds->preconditioner, *(MPI_Comm*)primme_svds->commInfo);
         ApplyPCPrecPETSCGen(aux, ldx, (SCALAR*)y+(*ldy)*i, ldy, &one, 1,
            primme_svds->preconditioner, *(MPI_Comm*)primme_svds->commInfo);
      }
      free(aux);
   }
   else if (*mode == primme_svds_op_augmented) {
      ApplyPCPrecPETSCGen((SCALAR*)x+primme_svds->nLocal, ldx, y, ldy, blockSize, 0,
         primme_svds->preconditioner, *(MPI_Comm*)primme_svds->commInfo);
      ApplyPCPrecPETSCGen(x, ldx, (SCALAR*)y+primme_svds->nLocal, ldy, blockSize, 1,
         primme_svds->preconditioner, *(MPI_Comm*)primme_svds->commInfo);
   }
   *ierr = 0;
}


void ApplyInvDavidsonDiagPrecPETSc(void *x, PRIMME_INT *ldx, void *y,
      PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   int i, j;
   double shift, d, minDenominator;
   SCALAR *xvec, *yvec;
   const int nLocal = primme->nLocal, bs = *blockSize;
   const PetscScalar *diag;
   Vec vec;
   PetscErrorCode ierr;
   
   vec = *(Vec *)primme->preconditioner;
   xvec = (SCALAR *)x;
   yvec = (SCALAR *)y;
   minDenominator = 1e-14*(primme->aNorm >= 0.0L ? primme->aNorm : 1.);

   ierr = VecGetArrayRead(vec, &diag); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   for (i=0; i<bs; i++) {
      shift = primme->ShiftsForPreconditioner[i];
      for (j=0; j<nLocal; j++) {
         d = diag[j] - shift;
         d = (fabs(d) > minDenominator) ? d : copysign(minDenominator, d);
         yvec[*ldy*i+j] = xvec[*ldx*i+j]/d;
      }
   }
   ierr = VecRestoreArrayRead(vec, &diag); CHKERRABORT(*(MPI_Comm*)primme->commInfo, ierr);
   *err = 0;
}

/******************************************************************************
 * Generates sum of square values per rows and then per columns 
 *
******************************************************************************/
#undef __FUNCT__
#define __FUNCT__ "getSumSquares"
static PetscErrorCode getSumSquares(Mat matrix, double *diag) {
   PetscErrorCode ierr;
   int i, j;
   double *sumr, *sumc;
   PetscInt n, mLocal, nLocal, low, high;
   PetscReal *aux;

   PetscFunctionBegin;

   ierr = MatGetSize(matrix, NULL, &n); CHKERRQ(ierr);
   ierr = MatGetLocalSize(matrix, &mLocal, &nLocal); CHKERRQ(ierr);
   sumr = diag; sumc = &diag[mLocal];

   ierr = PetscMalloc1(n, &aux); CHKERRQ(ierr);
   ierr = MatGetColumnNorms(matrix, NORM_2, aux); CHKERRQ(ierr);
   ierr = MatGetOwnershipRangeColumn(matrix, &low, &high);CHKERRQ(ierr);  
   for (i=low; i<high; i++) {
      sumc[i-low] = aux[i]*aux[i];
   }
   ierr = PetscFree(aux); CHKERRQ(ierr);

   ierr = MatGetOwnershipRange(matrix, &low, &high); CHKERRQ(ierr);
   for (i=low; i<high; i++) {
     PetscInt          ncols;
     const PetscInt    *cols;
     const PetscScalar *vals;

     sumr[i-low] = 0.0;
     ierr = MatGetRow(matrix, i, &ncols, &cols, &vals); CHKERRQ(ierr);
     for (j = 0; j < ncols; j++) {
       sumr[i-low] += PetscRealPart(vals[j]*PetscConj(vals[j]));
     }
     ierr = MatRestoreRow(matrix, i, &ncols, &cols, &vals); CHKERRQ(ierr);
   }

   PetscFunctionReturn(0);
}

int createInvNormalPrecPETSC(Mat matrix, double shift, double **prec) {
   int i;
   PetscInt mLocal, nLocal;
   double *diag, minDenominator=1e-14;

   MatGetLocalSize(matrix, &mLocal, &nLocal);
   diag = (double*)primme_calloc(mLocal+nLocal, sizeof(double), "diag");
   getSumSquares(matrix, diag);
   for (i=0; i<mLocal+nLocal; i++) {
      diag[i] -= shift*shift;
      if (fabs(diag[i]) < minDenominator)
         diag[i] = copysign(minDenominator, diag[i]);
   }
   *prec = diag;
   return 1;
}
