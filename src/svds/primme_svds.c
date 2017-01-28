/*******************************************************************************
 * Copyright (c) 2016, College of William & Mary
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the College of William & Mary nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * PRIMME: https://github.com/primme/primme
 * Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
 *******************************************************************************
 * File: primme_svds.c
 *
 * Purpose - front end to svd problems. 
 *
 ******************************************************************************/
 
#include <stdlib.h>   /* free, qsort */
#include <stdio.h>  
#include <string.h>  
#include <math.h>  
#include <assert.h>  
#include "numerical.h"
#include "../eigs/ortho.h"
#include "wtime.h"
#include "primme_interface.h"
#include "primme_svds_interface.h"

#define ALLOCATE_WORKSPACE_FAILURE -1
#define MALLOC_FAILURE             -3

static int primme_svds_check_input(REAL *svals, SCALAR *svecs, 
        REAL *resNorms, primme_svds_params *primme_svds);
static SCALAR* copy_last_params_from_svds(primme_svds_params *primme_svds, int stage,
      REAL *svals, SCALAR *svecs, REAL *rnorms, int *allocatedTargetShifts);
static int copy_last_params_to_svds(primme_svds_params *primme_svds, int stage,
      REAL *svals, SCALAR *svecs, REAL *rnorms, int allocatedTargetShifts);
static void applyPreconditionerSVDS(void *x, PRIMME_INT *ldx, void *y,
      PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
static void matrixMatvecSVDS(void *x_, PRIMME_INT *ldx, void *y_,
      PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
static void Num_scalInv_Smatrix(SCALAR *x, PRIMME_INT m, int n, PRIMME_INT ldx, REAL *factors,
                                       primme_svds_params *primme_svds);
static int allocate_workspace_svds(primme_svds_params *primme_svds, int allocate);
static int globalSum_Rprimme_svds(REAL *sendBuf, REAL *recvBuf, int count, 
      primme_svds_params *primme_svds);
static void convTestFunAugmented(double *eval, void *evec, double *rNorm, int *isConv,
   primme_params *primme, int *ierr);
static void convTestFunATA(double *eval, void *evec, double *rNorm, int *isConv,
   primme_params *primme, int *ierr);

/*******************************************************************************
 * Subroutine Sprimme_svds - This routine is a front end used to perform 
 *    error checking on the input parameters, perform validation, 
 *    and make the call to main_iter. 
 *
 *    Calling Sprimme_svds with all svals, svecs, resNorms set to NULL
 *    returns the int and real memory required in the following primme fields:
 *            int primme->intWorkSize : bytes of int workspace needed
 *       long int primme->realWorkSize: bytes of real workspace needed
 * 
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * evals  Contains the converged Ritz values upon return.   Should be of size 
 *        primme->numEvals.
 * 
 * evecs  The local portions of the converged Ritz vectors.  The dimension of
 *        the array is at least primme->nLocal*primme->numEvals
 *
 * resNorms  The residual norms of the converged Ritz vectors.  Should be of 
 *           size primme->numEvals
 *  
 * primme  Structure containing various solver parameters and statistics
 *         See readme.txt for INPUT/OUTPUT variables in primme
 *
 * Return Value
 * ------------
 *  0 - Success
 * -1 - Failure to allocate workspace
 * -2 - Malloc failure in allocating a permutation integer array
 * -3 - main_iter encountered a problem
 * -4 ...-32 - Invalid input (parameters or primme struct) returned 
 *             by check_input()
 * -100...-199 - PRIMME error code from first stage
 * -200...-299 - PRIMME error code from second stage
 *
 ******************************************************************************/

int Sprimme_svds(REAL *svals, SCALAR *svecs, REAL *resNorms, 
      primme_svds_params *primme_svds) {

   int ret, allocatedTargetShifts;
   SCALAR *svecs0;

   /* ------------------ */
   /* Set some defaults  */
   /* ------------------ */
   primme_svds_set_defaults(primme_svds);
   if (fabs(primme_svds->eps) == 0.0) {
      primme_svds->eps = MACHINE_EPSILON*1e4;
   }

   /* -------------------------------------------------------------- */
   /* If needed, we are ready to estimate required memory and return */
   /* -------------------------------------------------------------- */
    if (svals == NULL && svecs == NULL && resNorms == NULL)
       return allocate_workspace_svds(primme_svds, 0 /* don't allocate */);

   /* ----------------------------------------------------------- */
   /* Primme_svds_initialize must be called by users unless users */  
   /* specify all parameters in primme_svds structure. Check if   */
   /* primme_svds inputs are good for bounds, correct values etc. */
   /* ----------------------------------------------------------- */
   ret = primme_svds_check_input(svals, svecs, resNorms, primme_svds); 
   if (ret != 0) {
      return(ret);
   }

   /* ----------------------------------------------------------------------- */
   /* Compute AND allocate memory requirements for main_iter and subordinates */
   /* ----------------------------------------------------------------------- */
   ret = allocate_workspace_svds(primme_svds, 1 /* allocate */);
   if (ret != 0) {
      return ALLOCATE_WORKSPACE_FAILURE;
   }

   /* --------------- */
   /* Execute stage 1 */
   /* --------------- */

   CHKERRS((svecs0 = copy_last_params_from_svds(primme_svds, 0, NULL, svecs,
               NULL, &allocatedTargetShifts)) == NULL,
         ALLOCATE_WORKSPACE_FAILURE);

   ret = Sprimme(svals, svecs0, resNorms, &primme_svds->primme); 

   CHKERRS(copy_last_params_to_svds(primme_svds, 0, svals, svecs, resNorms,
            allocatedTargetShifts), ALLOCATE_WORKSPACE_FAILURE);

   if(ret != 0) {
      return ret - 100;
   }
   if (primme_svds->methodStage2 == primme_svds_op_none) {
      return 0;
   }

   /* --------------- */
   /* Execute stage 2 */
   /* --------------- */

   CHKERRS((svecs0 = copy_last_params_from_svds(primme_svds, 1, svals, svecs,
            resNorms, &allocatedTargetShifts)) == NULL,
         ALLOCATE_WORKSPACE_FAILURE);

   /* The value numSvals-primme->numEvals indicates how many svals */
   /* are already converged. So shift svals and resnorms that much */
   int nconv = primme_svds->numSvals - primme_svds->primmeStage2.numEvals;

   ret = Sprimme(svals+nconv, svecs0, resNorms+nconv, &primme_svds->primmeStage2);

   CHKERRS(copy_last_params_to_svds(primme_svds, 1, svals, svecs, resNorms,
         allocatedTargetShifts), ALLOCATE_WORKSPACE_FAILURE);

   if(ret != 0) {
      return ret - 200;
   }
   return 0;
}

static int comp_double(const void *a, const void *b)
{
   return *(double*)a <= *(double*)b ? -1 : 1;
}

static SCALAR* copy_last_params_from_svds(primme_svds_params *primme_svds, int stage,
      REAL *svals, SCALAR *svecs, REAL *rnorms, int *allocatedTargetShifts) {

   primme_params *primme;
   primme_svds_operator method;
   SCALAR *aux, *out_svecs = svecs;
   int n, nMax, i, cut;
   const double machEps = MACHINE_EPSILON;

   primme = stage == 0 ? &primme_svds->primme : &primme_svds->primmeStage2;
   method = stage == 0 ? primme_svds->method : primme_svds->methodStage2;

   *allocatedTargetShifts = 0;

   if (method == primme_svds_op_none) {
      primme->maxMatvecs = 1;
      return NULL;
   }

   if (!primme->matrixMatvec) {
      primme->matrixMatvec = matrixMatvecSVDS;
      primme->matrix = primme_svds;
   }
   if (!primme->applyPreconditioner) {
      primme->applyPreconditioner = applyPreconditionerSVDS;
      primme->preconditioner = primme_svds;
   }

   if (primme_svds->aNorm > 0.0) {
      switch(method) {
      case primme_svds_op_AtA:
      case primme_svds_op_AAt:
         primme->aNorm = primme_svds->aNorm*primme_svds->aNorm;
         break;
      case primme_svds_op_augmented:
         primme->aNorm = primme_svds->aNorm;
         break;
      case primme_svds_op_none:
         break;
      }
   }

   switch(method) {
   case primme_svds_op_AtA:
      primme->convTestFun = convTestFunATA;
      break;
   case primme_svds_op_AAt:
      primme->convTestFun = convTestFunATA;
      break;
   case primme_svds_op_augmented:
      primme->convTestFun = convTestFunAugmented;
      break;
   case primme_svds_op_none:
      break;
   }

   /* Set properly initial vectors. Now svecs = [Uc U0 Vc V0], where          */
   /* Uc, m x numOrthoConst, left constrain vectors;                          */
   /* U0, m x initSize, left initial vectors;                                 */
   /* Vc, n x numOrthoConst, right constrain vectors;                         */
   /* V0, n x numOrthoConst, right initial vectors.                           */

   primme->initSize = primme_svds->initSize;
   primme->numOrthoConst = primme_svds->numOrthoConst;
   n = primme_svds->initSize + primme_svds->numOrthoConst;
   nMax = max(primme_svds->initSize, primme_svds->numSvals) + primme_svds->numOrthoConst;
   switch(method) {
   case primme_svds_op_AtA:
   case primme_svds_op_AAt:
      /* Move Vc V0 to the rightmost position in svecs (aux).
         If using AA', only move Vc */
      aux = &svecs[nMax*primme_svds->mLocal];
      Num_copy_matrix_Sprimme(&svecs[primme_svds->mLocal*n], primme_svds->nLocal,
         method == primme_svds_op_AtA ? n : primme_svds->numOrthoConst,
         primme_svds->nLocal, aux, primme_svds->nLocal);
      if (method == primme_svds_op_AtA) out_svecs = aux;
      break;
   case primme_svds_op_augmented:
      /* Shuffle svecs so that svecs = [V; U] */
      assert(primme->nLocal == primme_svds->mLocal+primme_svds->nLocal);
      CHKERRS(MALLOC_PRIMME(primme->nLocal*n, &aux), NULL);
      Num_copy_Sprimme(primme->nLocal*n, svecs, 1, aux, 1);
      Num_copy_matrix_Sprimme(&aux[primme_svds->mLocal*n], primme_svds->nLocal,
         n, primme_svds->nLocal, svecs, primme->nLocal);
      Num_copy_matrix_Sprimme(aux, primme_svds->mLocal, n, primme_svds->mLocal,
         &svecs[primme_svds->nLocal], primme->nLocal);
      free(aux);

      /* Normalize the orthogonal constrains */
      Num_scal_Sprimme(primme->nLocal*primme_svds->numOrthoConst, 1./sqrt(2.),
            svecs, 1);
      break;
   case primme_svds_op_none:
      break;
   }
   primme->iseed[0] = primme_svds->iseed[0];
   primme->iseed[1] = primme_svds->iseed[1];
   primme->iseed[2] = primme_svds->iseed[2];
   primme->iseed[3] = primme_svds->iseed[3];
   primme->maxMatvecs = primme_svds->maxMatvecs;

   primme->intWork = primme_svds->intWork;
   primme->intWorkSize = primme_svds->intWorkSize;
   /* If matrixMatvecSVDS is used, it needs extra space to compute A*A' or A'*A */
   if ((primme->matrixMatvec == matrixMatvecSVDS) &&
       (method == primme_svds_op_AtA || method == primme_svds_op_AAt)) {
      cut = primme->maxBlockSize * (method == primme_svds_op_AtA ?
                     primme_svds->mLocal : primme_svds->nLocal);
   }
   else {
      cut = 0;
   }
   primme->realWork = (SCALAR*)primme_svds->realWork + cut;
   assert(primme_svds->realWorkSize >= cut*sizeof(SCALAR));
   primme->realWorkSize = primme_svds->realWorkSize - cut*sizeof(SCALAR);
 
   if ((stage == 0 && primme_svds->numTargetShifts > 0) ||
       (stage == 1 && primme->targetShifts == NULL &&
         primme_svds->target == primme_svds_closest_abs)) {
      primme->numTargetShifts = primme_svds->numTargetShifts;
      if (stage == 0 &&
            (method == primme_svds_op_AtA || method == primme_svds_op_AAt)) {
         *allocatedTargetShifts = 1;
         CHKERRS(MALLOC_PRIMME(primme_svds->numSvals, &primme->targetShifts),
            NULL);
         for (i=0; i<primme->numTargetShifts; i++) {
            primme->targetShifts[i] = 
               primme_svds->targetShifts[i]*primme_svds->targetShifts[i];
         }
      }
      else {
         primme->targetShifts = primme_svds->targetShifts;
      } 
   }
   else if (stage == 1 && primme->targetShifts == NULL &&
            primme_svds->target == primme_svds_smallest) {

      assert(method == primme_svds_op_augmented);
      *allocatedTargetShifts = 1;
      CHKERRS(MALLOC_PRIMME(primme_svds->numSvals, &primme->targetShifts),
            NULL);

      /* primme was configured to find the closest but greater values than */
      /* some shift. The eigensolver is not able to distinguish eigenvalues*/
      /* separated by less than machEps*|A|. The augmented matrix has      */
      /* |m-n| eigenpairs with value zero that don't correspond to         */
      /* singular triplets of A. To avoid to return incorrect triplets set */
      /* shifts not smaller than machEps*|A|.                              */

      for (i=0; i<primme_svds->initSize; i++) {
         primme->targetShifts[i] = max(svals[i]-rnorms[i], primme_svds->aNorm*machEps);
      }
      for ( ; i<primme_svds->numSvals; i++) {
         primme->targetShifts[i] = primme_svds->aNorm*machEps;
      }

      /* Sort the shifts in ascending order */

      qsort(primme->targetShifts, primme_svds->numSvals, sizeof(double),
            comp_double);
      primme->numTargetShifts = primme_svds->numSvals;

   }
   else if (method == primme_svds_op_augmented &&
         primme_svds->target == primme_svds_smallest &&
         primme->targetShifts == NULL) {

      CHKERRS(MALLOC_PRIMME(1, &primme->targetShifts), NULL);
      *allocatedTargetShifts = 1;
      primme->targetShifts[0] = 0.0;
      primme->numTargetShifts = 1;

   }

   /* Set an initial guess [x; A'x] or [Ax; x] if there is no initial guess   */
   /* and augmented matrix will be used                                       */

   if (method == primme_svds_op_augmented && primme->initSize <= 0) {
      int ONE = 1, NOTRANS = 0, TRANS = 1, ierr=0;
      REAL norms2_[2], norms2[2];
      if (primme_svds->m >= primme_svds->n) {
         Num_larnv_Sprimme(2, primme->iseed, primme_svds->mLocal,
               &svecs[primme_svds->nLocal]);
         CHKERRMS((primme_svds->matrixMatvec(&svecs[primme_svds->nLocal],
                     &primme_svds->mLocal, svecs, &primme_svds->nLocal, &ONE,
                     &TRANS, primme_svds, &ierr), ierr), NULL,
               "Error returned by 'matrixMatvec' %d", ierr);
      }
      else {
         Num_larnv_Sprimme(2, primme->iseed, primme_svds->nLocal, svecs);
         CHKERRMS((primme_svds->matrixMatvec(svecs, &primme_svds->nLocal,
                     &svecs[primme_svds->nLocal], &primme_svds->mLocal, &ONE,
                     &NOTRANS, primme_svds, &ierr), ierr), NULL,
               "Error returned by 'matrixMatvec' %d", ierr);
      }
      norms2_[0] = REAL_PART(Num_dot_Sprimme(primme_svds->nLocal, svecs, 1,
               svecs, 1));
      norms2_[1] = REAL_PART(Num_dot_Sprimme(primme_svds->mLocal,
               &svecs[primme_svds->nLocal], 1, &svecs[primme_svds->nLocal], 1));
      globalSum_Rprimme_svds(norms2_, norms2, 2, primme_svds);
      Num_scal_Sprimme(primme_svds->nLocal, 1.0/sqrt(norms2[0]), svecs, 1);
      Num_scal_Sprimme(primme_svds->mLocal, 1.0/sqrt(norms2[1]),
            &svecs[primme_svds->nLocal], 1);
      primme->initSize = 1;
      if (rnorms) rnorms[0] = HUGE_VAL;
      primme->initBasisMode = primme_init_user;
   }

   /* If second stage, set as numOrthoConst the first ones that pass */
   /* the convergence criterion.                                     */

   if (stage == 1) {
      for (i=0; primme->initSize > 0; i++) {
         double ev = (double)svals[i], resnorm = rnorms[i];
         int isConv=0, ierr=0;
         CHKERRMS((primme->convTestFun(&ev, NULL, &resnorm, &isConv, primme,
                     &ierr), ierr), NULL,
               "Error code returned by 'convTestFun' %d", ierr);
         if (!isConv) break;
         primme->numOrthoConst++;
         primme->initSize--;
         primme->numEvals--;
      }
   }

   /* Set locking */   

   if (primme_svds->locking >= 0) {
      primme->locking = primme_svds->locking;
   }

   return out_svecs;
}

/******************************************************************************
 * Function allocate_workspace_svds - This function computes the amount of
 *    integer and real workspace needed by the solver and possibly allocates
 *    the space 
 *
 * Input: 
 *   allocate  If 0, no allocation occurs, but the amounts of int and real 
 *                       workspaces in BYTES are returned in the primme fields 
 *                       primme.intWorkSize, and primme.realWorkSize 
 *             If  !0, and if the user-provided space is not sufficient,
 *                       allocation is also performed.
 *
 * Output
 *  primme.intWorkSize   Size of integer space allocated in bytes
 *  primme.realWorkSize  Size of real space allocated in bytes (LONG INT)
 *  *(primme.intWork)    Pointer to the integer space allocated
 *  *(primme.realWork)   Pointer to the real space allocated
 *   
 * 
 * Return value
 * ------------
 * int -  0 if (allocate == true) and the given workspaces are large enough or
 *             have been allocated successfully
 *       -1 if (allocate == true) and memory allocation has failed
 *        1 if (allocate==false) 
 *
 ******************************************************************************/

static int allocate_workspace_svds(primme_svds_params *primme_svds, int allocate) {
   primme_params primme;
   int intWorkSize=0;         /* Size of int work space */
   size_t realWorkSize=0;     /* Size of real work space */

   /* Require workspace for 1st stage */
   if (primme_svds->method != primme_svds_op_none) {
      primme = primme_svds->primme;
      Sprimme(NULL, NULL, NULL, &primme);
      intWorkSize = primme.intWorkSize;
      realWorkSize = primme.realWorkSize;
      /* If matrixMatvecSVDS is used, it needs extra space to compute A*A' or A'*A */
      if ((primme.matrixMatvec == NULL || primme.matrixMatvec == matrixMatvecSVDS) &&
          (primme_svds->method == primme_svds_op_AtA || primme_svds->method == primme_svds_op_AAt))
         realWorkSize += primme.maxBlockSize * sizeof(SCALAR) *
                           (primme_svds->method == primme_svds_op_AtA ?
                              primme_svds->mLocal : primme_svds->nLocal);
   }

   /* Require workspace for 2st stage */
   if (primme_svds->methodStage2 != primme_svds_op_none) {
      assert(primme_svds->methodStage2 != primme_svds_op_AtA &&
             primme_svds->methodStage2 != primme_svds_op_AAt);
      primme = primme_svds->primmeStage2;
      /* Check the case where all pairs from first stage are converged. */
      /* More numOrthoConst requires more memory */
      primme.numOrthoConst += primme.numEvals;
      Sprimme(NULL, NULL, NULL, &primme);
      intWorkSize = max(intWorkSize, primme.intWorkSize);
      realWorkSize = max(realWorkSize, primme.realWorkSize);
   }

   if (!allocate) {
      primme_svds->intWorkSize  = intWorkSize;
      primme_svds->realWorkSize = realWorkSize;
      return 1;
   }

   /*----------------------------------------------------------------------*/
   /* Allocate the required workspace, if the user did not provide enough  */
   /*----------------------------------------------------------------------*/
   if (primme_svds->realWork != NULL
         && primme_svds->realWorkSize < realWorkSize) {
      return -20;
   }
   else if (primme_svds->realWork == NULL) {
      primme_svds->realWorkSize = realWorkSize;
      if (primme_svds->printLevel >= 5) fprintf(primme_svds->outputFile, 
         "Allocating real workspace: %ld bytes\n", primme_svds->realWorkSize);
      CHKERRMS(MALLOC_PRIMME(realWorkSize, (char**)&primme_svds->realWork),
            MALLOC_FAILURE, "Failed to allocate %zd bytes\n", realWorkSize);
   }

   if (primme_svds->intWork != NULL && primme_svds->intWorkSize < intWorkSize) {
      return -21;
   }
   else if (primme_svds->intWork == NULL) {
      primme_svds->intWorkSize = intWorkSize;
      if (primme_svds->printLevel >= 5) fprintf(primme_svds->outputFile, 
         "Allocating integer workspace: %d bytes\n", primme_svds->intWorkSize);
      CHKERRMS(MALLOC_PRIMME(intWorkSize/sizeof(int), &primme_svds->intWork),
            MALLOC_FAILURE,
            "Failed to allocate %d bytes\n", primme_svds->intWorkSize);
   }

   return 0;
}
 
int copy_last_params_to_svds(primme_svds_params *primme_svds, int stage,
      REAL *svals, SCALAR *svecs, REAL *rnorms, int allocatedTargetShifts) {

   int trans = 1, notrans = 0;
   primme_params *primme;
   primme_svds_operator method;
   SCALAR *aux;
   REAL *norms2, *norms2_;
   int n, nMax, i, ierr;

   primme = stage == 0 ? &primme_svds->primme : &primme_svds->primmeStage2;
   method = stage == 0 ? primme_svds->method : primme_svds->methodStage2;

   if (method == primme_svds_op_none) {
      primme->maxMatvecs = 1;
      return 0;
   }

   /* Pass back the converged vectors in first stage to regular vectors */

   if (stage == 1) {
      int nconv = primme_svds->numSvals - primme->numEvals;
      primme->initSize += nconv;
      primme->numOrthoConst -= nconv;
      primme->numEvals += nconv;
   }

   /* Record performance measurements */ 
   primme_svds->stats.numOuterIterations += primme->stats.numOuterIterations;
   primme_svds->stats.numRestarts        += primme->stats.numRestarts;
   primme_svds->stats.numMatvecs         += primme->stats.numMatvecs;
   primme_svds->stats.numPreconds        += primme->stats.numPreconds;
   primme_svds->stats.elapsedTime        += primme->stats.elapsedTime;


   if (primme->aNorm > 0.0) {
      switch(method) {
      case primme_svds_op_AtA:
      case primme_svds_op_AAt:
         primme_svds->aNorm = sqrt(primme->aNorm);
         break;
      case primme_svds_op_augmented:
         primme_svds->aNorm = primme->aNorm;
         break;
      case primme_svds_op_none:
         break;
      }
   }

   if (method == primme_svds_op_AtA || method == primme_svds_op_AAt) {
      for (i=0; i<primme->initSize; i++) {
         svals[i] = sqrt(max(0.0, svals[i]));
      }
   }
 
   /* Set svecs = [Uc U Vc V] */
   nMax = max(primme_svds->initSize, primme_svds->numSvals) + primme_svds->numOrthoConst;
   primme_svds->initSize = primme->initSize;
   n = primme_svds->initSize + primme_svds->numOrthoConst;
   switch(method) {
   case primme_svds_op_AtA:
      /* Transform svecs to [Uc A*V/Sigma Vc V] */
      CHKERRMS((primme_svds->matrixMatvec(
            &svecs[primme_svds->mLocal*nMax+primme->nLocal*primme_svds->numOrthoConst],
            &primme_svds->nLocal, &svecs[primme_svds->mLocal*primme_svds->numOrthoConst],
            &primme_svds->mLocal, &primme_svds->initSize, &notrans, primme_svds,
            &ierr), ierr), -1,
         "Error returned by 'matrixMatvec' %d", ierr);
      Num_scalInv_Smatrix(&svecs[primme_svds->mLocal*primme_svds->numOrthoConst],
            primme_svds->mLocal, primme_svds->initSize, primme_svds->mLocal, svals, primme_svds);
      Num_copy_matrix_Sprimme(&svecs[primme_svds->mLocal*nMax], primme_svds->nLocal, n,
            primme_svds->nLocal, &svecs[primme_svds->mLocal*n], primme_svds->nLocal);
      break;
   case primme_svds_op_AAt:
      /* Transform svecs to [Uc U Vc A'*U/Sigma] */
      Num_copy_matrix_Sprimme(&svecs[primme_svds->mLocal*nMax], primme_svds->nLocal,
            primme_svds->numOrthoConst, primme_svds->nLocal,
            &svecs[primme_svds->mLocal*n], primme_svds->nLocal);
      CHKERRMS((primme_svds->matrixMatvec(
            &svecs[primme_svds->mLocal*primme_svds->numOrthoConst], &primme_svds->mLocal,
            &svecs[primme_svds->mLocal*n+primme->nLocal*primme_svds->numOrthoConst],
            &primme_svds->nLocal, &primme_svds->initSize, &trans, primme_svds,
            &ierr), ierr), -1,
         "Error returned by 'matrixMatvec' %d", ierr);
      Num_scalInv_Smatrix(
            &svecs[primme_svds->mLocal*n+primme->nLocal*primme_svds->numOrthoConst],
            primme_svds->nLocal, primme_svds->initSize, primme_svds->nLocal, svals, primme_svds);
      break;
   case primme_svds_op_augmented:
      assert(primme->nLocal == primme_svds->mLocal+primme_svds->nLocal);

      /* Normalize back the orthogonal constrains */
      Num_scal_Sprimme(primme->nLocal*primme_svds->numOrthoConst, sqrt(2.),
            svecs, 1);

      /* Shuffle svecs from [Vc V; Uc U] to [Uc U Vc V] */
      CHKERRS(MALLOC_PRIMME(primme->nLocal*n, &aux), -1);
      Num_copy_Sprimme(primme->nLocal*n, svecs, 1, aux, 1);
      Num_copy_matrix_Sprimme(aux, primme_svds->nLocal, n, primme->nLocal,
         &svecs[primme_svds->mLocal*n], primme_svds->nLocal);
      Num_copy_matrix_Sprimme(&aux[primme_svds->nLocal], primme_svds->mLocal, n,
         primme->nLocal, svecs, primme_svds->mLocal);
      free(aux);

      /* Normalize every column in U and V */
      CHKERRS(MALLOC_PRIMME(4*n, &norms2_), -1);
      norms2 = norms2_ + 2*n;
      for (i=0; i<n; i++) {
         norms2_[i] = REAL_PART(Num_dot_Sprimme(primme_svds->mLocal,
            &svecs[primme_svds->mLocal*i], 1, &svecs[primme_svds->mLocal*i],
            1));
      }
      for (i=0; i<n; i++) {
         norms2_[n+i] = REAL_PART(Num_dot_Sprimme(primme_svds->nLocal,
               &svecs[primme_svds->mLocal*n+primme_svds->nLocal*i], 1,
               &svecs[primme_svds->mLocal*n+primme_svds->nLocal*i], 1));
      }
      globalSum_Rprimme_svds(norms2_, norms2, 2*n, primme_svds);
      for (i=0; i<n; i++) {
         Num_scal_Sprimme(primme_svds->mLocal, 1.0/sqrt(norms2[i]),
               &svecs[primme_svds->mLocal*i], 1);
      }
      for (i=0; i<n; i++) {
         Num_scal_Sprimme(primme_svds->nLocal, 1.0/sqrt(norms2[n+i]),
               &svecs[primme_svds->mLocal*n+primme_svds->nLocal*i], 1);
      }
      free(norms2_);
      break;
   case primme_svds_op_none:
      break;
   }

   primme_svds->iseed[0] = primme->iseed[0];
   primme_svds->iseed[1] = primme->iseed[1];
   primme_svds->iseed[2] = primme->iseed[2];
   primme_svds->iseed[3] = primme->iseed[3];
   primme_svds->maxMatvecs -= primme->stats.numMatvecs;

   /* Zero references to primme workspaces to prevent to be release by primme_free */
   primme->intWork = NULL;
   primme->realWork = NULL;

   if (allocatedTargetShifts) {
      free(primme->targetShifts);
      primme->targetShifts = NULL;
   }

   /* Update residual norms when final stage */
   if (primme_svds->methodStage2 != primme_svds_op_none) {
      switch(method) {
      case primme_svds_op_AtA:
      case primme_svds_op_AAt:
         for (i=0; i<primme_svds->initSize; i++) {
            rnorms[i] = min(rnorms[i]/svals[i], primme_svds->aNorm);
         }
         break;
      case primme_svds_op_augmented:
         for (i=0; i<primme_svds->initSize; i++) {
            rnorms[i] *= sqrt(2.0);
         }
         break;
      case primme_svds_op_none:
         break;
      }
   }

   return 0;
}

/******************************************************************************
 *
 * static int primme_svds_check_input(double *svals, SCALAR *svecs, double *resNorms, 
 *                        primme_svds_params *primme_svds) 
 *
 * INPUT
 * -----
 *  svals, svecs, resNorms   Output arrays for primme
 *  primme_svds              the main structure of parameters 
 *
 * return value -   0    If input parameters in primme are appropriate
 *              -4..-19  Inappropriate input parameters were found
 *
 ******************************************************************************/
static int primme_svds_check_input(REAL *svals, SCALAR *svecs, REAL *resNorms, 
      primme_svds_params *primme_svds) {
   int ret;
   ret = 0;

   if (primme_svds == NULL)
      ret = -4;
   else if (primme_svds->n < 0 || primme_svds->m < 0 || primme_svds->nLocal < 0
            || primme_svds->mLocal < 0 || primme_svds->nLocal > primme_svds->n
            || primme_svds->mLocal > primme_svds->m) 
      ret = -5;
   else if (primme_svds->numProcs < 1)
      ret = -6;
   else if (primme_svds->matrixMatvec == NULL) 
      ret = -7;
   else if (primme_svds->applyPreconditioner == NULL && 
         primme_svds->precondition == 1) 
      ret = -8;
   else if (primme_svds->numProcs >1 && primme_svds->globalSumReal == NULL)
      ret = -9;
   else if (primme_svds->numSvals > min(primme_svds->n, primme_svds->m))
      ret = -10;
   else if (primme_svds->numSvals < 1)
      ret = -11;
   else if ( primme_svds->target != primme_svds_smallest  &&
             primme_svds->target != primme_svds_largest   &&
             primme_svds->target != primme_svds_closest_abs)
      ret = -13;
   else if ( primme_svds->method != primme_svds_op_AtA &&
             primme_svds->method != primme_svds_op_AAt &&
             primme_svds->method != primme_svds_op_augmented)
      ret = -14;
   else if ( (primme_svds->method == primme_svds_op_augmented &&
              primme_svds->methodStage2 != primme_svds_op_none) ||
             (primme_svds->method != primme_svds_op_augmented &&
              primme_svds->methodStage2 != primme_svds_op_augmented &&
              primme_svds->methodStage2 != primme_svds_op_none))
      ret = -15;
   else if (primme_svds->printLevel < 0 || primme_svds->printLevel > 5)
      ret = -16; 
   else if (svals == NULL)
      ret = -17;
   else if (svecs == NULL)
      ret = -18;
   else if (resNorms == NULL)
      ret = -19;
   /* Booked -20 and -21*/

   return ret;
   /***************************************************************************/
} /* end of check_input
   ***************************************************************************/

/**********************************************************************************
 * void MatrixATA_Matvec(void *x, void *y, int *blockSize, primme_params *primme) *
 **********************************************************************************/
static void matrixMatvecSVDS(void *x_, PRIMME_INT *ldx, void *y_,
      PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {

   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   int trans = 1, notrans = 0;
   SCALAR *x = (SCALAR*)x_, *y = (SCALAR*)y_;
   primme_svds_operator method = &primme_svds->primme == primme ?
      primme_svds->method : primme_svds->methodStage2;
   int i, bs;

   switch(method) {
   case primme_svds_op_AtA:
      for (i=0, bs=min((*blockSize-i), primme->maxBlockSize); bs>0;
               i+= bs, bs=min((*blockSize-i), primme->maxBlockSize))
      {
         primme_svds->matrixMatvec(&x[*ldx*i], ldx, primme_svds->realWork,
               &primme_svds->mLocal, &bs, &notrans, primme_svds, ierr);
         if (*ierr != 0) return;
         primme_svds->matrixMatvec(primme_svds->realWork, &primme_svds->mLocal,
            &y[*ldy*i], ldy, &bs, &trans, primme_svds, ierr);
         if (*ierr != 0) return;
      }
      break;
   case primme_svds_op_AAt:
      for (i=0, bs=min((*blockSize-i), primme->maxBlockSize); bs>0;
               i+= bs, bs=min((*blockSize-i), primme->maxBlockSize))
      {
         primme_svds->matrixMatvec(&x[*ldx*i], ldx, primme_svds->realWork,
               &primme_svds->nLocal, &bs, &trans, primme_svds, ierr);
         if (*ierr != 0) return;
         primme_svds->matrixMatvec(primme_svds->realWork, &primme_svds->nLocal,
            &y[*ldy*i], ldy, &bs, &notrans, primme_svds, ierr);
         if (*ierr != 0) return;
      }
      break;
   case primme_svds_op_augmented:
      primme_svds->matrixMatvec(&x[primme_svds->nLocal], ldx, y, ldy, blockSize,
            &trans, primme_svds, ierr);
         if (*ierr != 0) return;
      primme_svds->matrixMatvec(x, ldx, &y[primme_svds->nLocal],
         ldy, blockSize, &notrans, primme_svds, ierr);
         if (*ierr != 0) return;
      break;
   case primme_svds_op_none:
      break;
   }
}

static void applyPreconditionerSVDS(void *x, PRIMME_INT *ldx, void *y,
      PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr) {

   primme_svds_params *primme_svds = (primme_svds_params *) primme->preconditioner;
   int method = (int)(&primme_svds->primme == primme ?
                        primme_svds->method : primme_svds->methodStage2);

   primme_svds->applyPreconditioner(x, ldx, y, ldy, blockSize, &method,
         primme_svds, ierr);
}

static void Num_scalInv_Smatrix(SCALAR *x, PRIMME_INT m, int n, PRIMME_INT ldx,
      REAL *factors, primme_svds_params *primme_svds) {

   int i;
   REAL norm, norm0, factor;

   assert(ldx >= m);
   for (i=0; i<n; i++) {
      if (factors[i] > 0.0 && 1.0L/factors[i] < HUGE_VAL) {
         factor = factors[i];
      }
      else {
         norm0 = REAL_PART(Num_dot_Sprimme(m, &x[i*ldx], 1, &x[i*ldx], 1));
         globalSum_Rprimme_svds(&norm0, &norm, 1, primme_svds);
         factor = sqrt(norm);
      }
      Num_scal_Sprimme(m, 1.0/factor, &x[i*ldx], 1);
   }
}

static int globalSum_Rprimme_svds(REAL *sendBuf, REAL *recvBuf, int count, 
      primme_svds_params *primme_svds) {

   int ierr;

   if (primme_svds && primme_svds->globalSumReal) {
      CHKERRMS((primme_svds->globalSumReal(sendBuf, recvBuf, &count,
                  primme_svds, &ierr), ierr), -1,
            "Error returned by 'globalSumReal' %d", ierr);
   }
   else {
      Num_copy_Rprimme(count, sendBuf, 1, recvBuf, 1);
   }

   return 0;
}

/*******************************************************************************
 * Subroutine convTestFunATA - This routine implements primme_params.
 *    convTestFun and returns an approximate eigenpair converged when           
 *    resNorm < eps * sval * primme_svds.aNorm = eps * sqrt(eval*primme.aNorm)
 *    resNorm is close to machineEpsilon * primme.aNorm.
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * evec         The approximate eigenvector
 * eval         The approximate eigenvalue 
 * rNorm        The norm of the residual vector
 * primme       Structure containing various solver parameters
 *
 * OUTPUT PARAMETERS
 * ----------------------------------
 * isConv      if it isn't zero the approximate pair is marked as converged
 ******************************************************************************/

static void convTestFunATA(double *eval, void *evec, double *rNorm, int *isConv,
   primme_params *primme, int *ierr) {

   const double machEps = MACHINE_EPSILON;
   const double aNorm = (primme->aNorm > 0.0) ?
      primme->aNorm : primme->stats.estimateLargestSVal;
   (void)evec;  /* unused argument */
   *isConv = *rNorm < max(
               primme->eps * sqrt(fabs(*eval * aNorm)),
               machEps * 3.16 * aNorm);
   *ierr = 0;
}

/*******************************************************************************
 * Subroutine convTestFunAugmented - This routine implements primme_params.
 *    convTestFun and returns an approximate eigenpair converged when           
 *
 *       sqrt(||Av - su||^2 + ||A'u - sv||^2) < ||A|| * eps = aNorm/sqrt(2)*eps
 *
 *    However, the previous test is expensive so it is only checked after the
 *    next one passes:
 *
 *       resNorm < ||A||*eps = aNorm/sqrt(2)*eps
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * evec         The approximate eigenvector
 * eval         The approximate eigenvalue 
 * rNorm        The norm of the residual vector
 * primme       Structure containing various solver parameters
 *
 * OUTPUT PARAMETERS
 * ----------------------------------
 * isConv      if it isn't zero the approximate pair is marked as converged
 ******************************************************************************/

static void convTestFunAugmented(double *eval, void *evec_, double *rNorm,
      int *isConv, primme_params *primme, int *ierr) {

   const double machEps = MACHINE_EPSILON;
   const double aNorm = (primme->aNorm > 0.0) ?
      primme->aNorm : primme->stats.estimateLargestSVal;
   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;
   SCALAR *evec = (SCALAR*)evec_;

   /* Pre-test */

   *isConv = 
      *rNorm < max(
               primme->eps / sqrt(2.0) * aNorm,
               machEps * 3.16 * aNorm) 
      && *eval >= aNorm*machEps;

   /* Actual test */

   if (*isConv && evec) {
      int one = 1;
      SCALAR *r = (SCALAR*)malloc(sizeof(SCALAR)*primme->nLocal);
      if (r == NULL) {*ierr = 1; return;}

      /* r = [0 A';A 0] * evec = [ A'u; Av ] */

      matrixMatvecSVDS(evec, &primme->nLocal, r, &primme->nLocal, &one, primme,
            ierr);
      if (*ierr != 0) return;
      primme->stats.numMatvecs++;
      
      /* ip[0] = ||evec[0:nLocal-1]|| = ||v|| */
      /* ip[1] = ||evec[nLocal:nLocal+mLocal-1]|| = ||u|| */

      REAL ip0[2], ip[2];
      ip0[0] = REAL_PART(Num_dot_Sprimme(primme_svds->nLocal, evec, 1, evec,
               1));
      ip0[1] = REAL_PART(Num_dot_Sprimme(primme_svds->mLocal,
               &evec[primme_svds->nLocal], 1, &evec[primme_svds->nLocal], 1));
      *ierr = globalSum_Rprimme_svds(ip0, ip, 2, primme_svds);
      if (*ierr != 0) return;

      /* r[0:nLocal-1] = r[0:nLocal-1]/ip[1] - eval * evec[0:nLocal-1]/ip[0]  */
      /*               = A'u/||u|| - eval*v/||v||                             */

      Num_scal_Sprimme(primme_svds->nLocal, 1.0/ip[1], r, 1);
      Num_axpy_Sprimme(primme_svds->nLocal, -(SCALAR)*eval/ip[1], evec, 1, r,
            1);

      /* r[nLocal:end] = r[nLocal:end]/ip[1] - eval * evec[nLocal:end]/ip[0] */
      /*               = Av/||v|| - eval*u/||u||                             */

      Num_scal_Sprimme(primme_svds->mLocal, 1.0/ip[0], &r[primme_svds->nLocal],
            1);
      Num_axpy_Sprimme(primme_svds->mLocal, -(SCALAR)*eval/ip[0],
            &evec[primme_svds->nLocal], 1, &r[primme_svds->nLocal], 1);

      /* normr = sqrt(||Av - su||^2 + ||A'u - sv||^2) */

      REAL normr0, normr;
      normr0 = REAL_PART(Num_dot_Sprimme(primme->nLocal, r, 1, r, 1));
      *ierr = globalSum_Rprimme_svds(&normr0, &normr, 1, primme_svds);
      if (*ierr != 0) return;
      normr = sqrt(normr);

      /* isConv = 1 iff normr <= ||A||*eps = aNorm/sqrt(2)*eps */

      if (normr < max(aNorm/sqrt(2.0)*primme->eps, machEps * 3.16 * aNorm)) { 
         *isConv = 1;
      }
      else {
         *isConv = 0;
      }

      free(r);
   }

   *ierr = 0;
} 
