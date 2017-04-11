/*******************************************************************************
 * Copyright (c) 2017, College of William & Mary
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
#include "../eigs/const.h"
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
static void default_monitor(void *basisSvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedSvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, primme_event *event, int *stage,
      primme_svds_params *primme_svds, int *err);
static void monitor_single_stage(void *basisEvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedEvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, primme_event *event, primme_params *primme,
      int *err);
static void monitor_stage1(void *basisEvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedEvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, primme_event *event, primme_params *primme,
      int *err);
static void monitor_stage2(void *basisEvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedEvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, primme_event *event, primme_params *primme,
      int *err);

#define UPDATE_STATS(PRIMME_SVDS_STATS, OP, PRIMME_STATS) {\
   (PRIMME_SVDS_STATS).numOuterIterations OP  (PRIMME_STATS).numOuterIterations;\
   (PRIMME_SVDS_STATS).numRestarts        OP  (PRIMME_STATS).numRestarts       ;\
   /* NOTE: for the augmented and for normal equations, every matvec for the  */\
   /* eigensolver involves the direct and the transpose matrix-vector product */\
   (PRIMME_SVDS_STATS).numMatvecs         OP  (PRIMME_STATS).numMatvecs*2      ;\
   (PRIMME_SVDS_STATS).numPreconds        OP  (PRIMME_STATS).numPreconds       ;\
   (PRIMME_SVDS_STATS).numGlobalSum       OP  (PRIMME_STATS).numGlobalSum      ;\
   (PRIMME_SVDS_STATS).volumeGlobalSum    OP  (PRIMME_STATS).volumeGlobalSum   ;\
   (PRIMME_SVDS_STATS).numOrthoInnerProds OP  (PRIMME_STATS).numOrthoInnerProds;\
   (PRIMME_SVDS_STATS).elapsedTime        OP  (PRIMME_STATS).elapsedTime       ;\
   (PRIMME_SVDS_STATS).timeMatvec         OP  (PRIMME_STATS).timeMatvec        ;\
   (PRIMME_SVDS_STATS).timePrecond        OP  (PRIMME_STATS).timePrecond       ;\
   (PRIMME_SVDS_STATS).timeOrtho          OP  (PRIMME_STATS).timeOrtho         ;\
   (PRIMME_SVDS_STATS).timeGlobalSum      OP  (PRIMME_STATS).timeGlobalSum     ;\
}


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

   /* ----------------------- */
   /* Set default monitor     */
   /* ----------------------- */

   if (!primme_svds->monitorFun) {
      primme_svds->monitorFun = default_monitor;
   }

   /* ----------------------- */
   /* Reset stats             */
   /* ----------------------- */

   primme_svds->stats.numOuterIterations            = 0; 
   primme_svds->stats.numRestarts                   = 0;
   primme_svds->stats.numMatvecs                    = 0;
   primme_svds->stats.numPreconds                   = 0;
   primme_svds->stats.numGlobalSum                  = 0;
   primme_svds->stats.volumeGlobalSum               = 0;
   primme_svds->stats.numOrthoInnerProds            = 0.0;
   primme_svds->stats.elapsedTime                   = 0.0;
   primme_svds->stats.timeMatvec                    = 0.0;
   primme_svds->stats.timePrecond                   = 0.0;
   primme_svds->stats.timeOrtho                     = 0.0;
   primme_svds->stats.timeGlobalSum                 = 0.0;

   /* --------------- */
   /* Execute stage 1 */
   /* --------------- */

   if (primme_svds->eps == 0.0) {
      primme_svds->eps = MACHINE_EPSILON*1e4;
   }

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
      /* If d^2 and d'^2 are the exact and the approximate eigenvalues     */
      /* from normal equations respectively, and assuming that d' >= d,    */
      /* then d can be lower bounded as:                                   */
      /*    d'^2 - d^2 <= |r_AtA| -> sqrt(d'^2-|r_AtA|) <= d               */
      /*                             sqrt(d'^2-|r|*d')  <= d               */

      double min_val = primme_svds->aNorm*machEps; 
      for (i=0; i<primme_svds->initSize; i++) {
         primme->targetShifts[i] = max(
               sqrt(fabs(max(svals[i]-rnorms[i], 0.0)*svals[i])), min_val);
      }
      for ( ; i<primme_svds->numSvals; i++) {
         primme->targetShifts[i] = min_val;
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
      int flags[primme->initSize];
      for (i=0; primme->initSize > 0; i++) {
         /* NOTE: convTestFun at this stage expects the residual norm for the */
         /*       the augmented problem; this is why the residual norm is     */
         /*       divided by sqrt(2).                                         */
         double ev = (double)svals[i], resnorm = rnorms[i]/sqrt(2.0);
         int isConv=0, ierr=0;
         primme_svds->stats.elapsedTime = primme_wTimer(0);
         CHKERRMS((primme->convTestFun(&ev, NULL, &resnorm, &isConv, primme,
                     &ierr), ierr), NULL,
               "Error code returned by 'convTestFun' %d", ierr);
         if (!isConv) break;

         /* Report a triplet is locked */
         int ip1 = i+1;
         flags[i] = CONVERGED;
         primme_event EVENT_LOCKED = primme_event_locked;
         int ZERO = 0;
         CHKERRMS((primme_svds->monitorFun(NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, svals, &ip1, flags, rnorms, NULL, NULL,
                     &EVENT_LOCKED, &ZERO, primme_svds, &ierr), ierr), NULL,
                  "Error code returned by 'monitorFun' %d", ierr);

         primme->numOrthoConst++;
         primme->initSize--;
         primme->numEvals--;
      }
   }

   /* Set locking */   

   if (primme_svds->locking >= 0) {
      primme->locking = primme_svds->locking;
   }

   /* Set monitor */

   if (primme->monitorFun == NULL) {
      if (primme_svds->methodStage2 == primme_svds_op_none) {
         primme->monitorFun = monitor_single_stage;
      } else if (stage == 0) {
         primme->monitorFun = monitor_stage1;
      } else {
         primme->monitorFun = monitor_stage2;
      }
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
   UPDATE_STATS(primme_svds->stats, +=, primme->stats);

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
      /* ip[2:4] = u'*A*v */

      REAL ip0[4], ip[4];
      ip0[0] = REAL_PART(Num_dot_Sprimme(primme_svds->nLocal, evec, 1, evec,
               1));
      ip0[1] = REAL_PART(Num_dot_Sprimme(primme_svds->mLocal,
               &evec[primme_svds->nLocal], 1, &evec[primme_svds->nLocal], 1));
      ip0[3] = 0.0;
      *(SCALAR*)&ip0[2] = Num_dot_Sprimme(primme_svds->mLocal,
               &evec[primme_svds->nLocal], 1, &r[primme_svds->nLocal], 1);
      *ierr = globalSum_Rprimme_svds(ip0, ip, 4, primme_svds);
      if (*ierr != 0) return;
      ip[0] = sqrt(ip[0]);
      ip[1] = sqrt(ip[1]);
      SCALAR sval = *(SCALAR*)&ip[2]/ip[0]/ip[1];
      /* r[0:nLocal-1] = r[0:nLocal-1]/ip[1] - sval * evec[0:nLocal-1]/ip[0]  */
      /*               = A'u/||u|| - sval*v/||v||                             */

      Num_scal_Sprimme(primme_svds->nLocal, 1.0/ip[1], r, 1);
      Num_axpy_Sprimme(primme_svds->nLocal, -sval/ip[0], evec, 1, r, 1);

      /* r[nLocal:end] = r[nLocal:end]/ip[0] - sval * evec[nLocal:end]/ip[1] */
      /*               = Av/||v|| - sval*u/||u||                             */

      Num_scal_Sprimme(primme_svds->mLocal, 1.0/ip[0], &r[primme_svds->nLocal],
            1);
      Num_axpy_Sprimme(primme_svds->mLocal, -sval/ip[1],
            &evec[primme_svds->nLocal], 1, &r[primme_svds->nLocal], 1);

      /* normr = sqrt(||Av - su||^2 + ||A'u - sv||^2) */

      REAL normr0, normr;
      normr0 = REAL_PART(Num_dot_Sprimme(primme->nLocal, r, 1, r, 1));
      *ierr = globalSum_Rprimme_svds(&normr0, &normr, 1, primme_svds);
      if (*ierr != 0) return;
      normr = sqrt(normr);

      /* isConv = 1 iff normr <= ||A||*eps = aNorm/sqrt(2)*eps */

      if (normr < max(aNorm/sqrt(2.0)*primme->eps, machEps * 5 * aNorm)) { 
         *isConv = 1;
      }
      else {
         *isConv = 0;
      }

      free(r);
   }

   *ierr = 0;
}


/*******************************************************************************
 * Subroutine default_monitor - report iterations, #MV, residual norm,
 *    singular values, etc. at every inner/outer iteration and when some triplet
 *    converges.       
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * basisSvals   The approximate singular values of the basis
 * basisSize    The size of the basis
 * basisFlags   The state of every approximate triplet of the basis (see conv_flags)
 * iblock       Indices of the approximate triplet in the block
 * blockSize    The size of the block
 * basisNorms   The approximate residual norms of the triplet of the basis
 * numConverged The number of triplets converged in the basis and the locked triplets
 *              (this value isn't monotonic!)
 * lockedSvals  The locked singular values
 * numLocked    The number of triplets locked
 * lockedFlags  The state of each locked triplet (see conv_flags)
 * lockedNorms  The residual norms of the locked pairs
 * inner_its    The number of performed QMR iterations in the current correction equation
 * LSRes        The residual norm of the linear system at the current QMR iteration
 * event        The event reported
 * stage        0 for first stage and 1 for second stage
 * primme_svds  Structure containing various solver parameters and statistics
 *
 * OUTPUT
 * ------
 * err          Error code
 * 
 ******************************************************************************/

static void default_monitor(void *basisSvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedSvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, primme_event *event, int *stage,
      primme_svds_params *primme_svds, int *err)
{
   REAL *basisSvals = (REAL*)basisSvals_, *basisNorms = (REAL*)basisNorms_,
        *lockedSvals = (REAL*)lockedSvals_, *lockedNorms = (REAL*)lockedNorms_,
        *LSRes = (REAL*)LSRes_;
   assert(event != NULL && primme_svds != NULL && stage != NULL);

   /* Only print report if this is proc zero */
   if (primme_svds->procID == 0) {
      switch(*event) {
      case primme_event_outer_iteration:
         assert(basisSvals && basisSize && basisFlags && iblock && blockSize
                && basisNorms && numConverged);
         if (primme_svds->printLevel >= 3) {
            int i;  /* Loop variable */
            for (i=0; i < *blockSize; i++) {
               fprintf(primme_svds->outputFile, 
                     "OUT %" PRIMME_INT_P " conv %d blk %d MV %" PRIMME_INT_P " Sec %E SV %13E |r| %.3E stage %d\n",
                     primme_svds->stats.numOuterIterations, *numConverged, i,
                     primme_svds->stats.numMatvecs,
                     primme_svds->stats.elapsedTime, basisSvals[iblock[i]],
                     (double)basisNorms[iblock[i]], *stage+1);
            }
         }
         break;
      case primme_event_inner_iteration:
         assert(basisSize && iblock && basisNorms && inner_its && LSRes);
         (void)inner_its;
         if (primme_svds->printLevel >= 4) {
            fprintf(primme_svds->outputFile,
                  "INN MV %" PRIMME_INT_P " Sec %e Sval %e Lin|r| %.3e SV|r| %.3e stage %d\n",
                  primme_svds->stats.numMatvecs, primme_svds->stats.elapsedTime,
                  (double)basisSvals[iblock[0]], (double)*LSRes,
                  (double)basisNorms[iblock[0]], *stage+1);
         }
        break;
      case primme_event_restart:
         break;
      case primme_event_reset:
         break;
      case primme_event_converged:
         if ((*stage == 0 && primme_svds->printLevel >= 2)
               || (primme_svds->printLevel >= 5))
            fprintf(primme_svds->outputFile, 
                  "#Converged %d sval[ %d ]= %e norm %e Mvecs %" PRIMME_INT_P " Time %g stage %d\n",
                  *numConverged, iblock[0], basisSvals[iblock[0]],
                  basisNorms[iblock[0]], primme_svds->stats.numMatvecs,
                  primme_svds->stats.elapsedTime, *stage+1);
         break;
      case primme_event_locked:
         if (primme_svds->printLevel >= 2) { 
            fprintf(primme_svds->outputFile, 
                  "Lock striplet[ %d ]= %e norm %.4e Mvecs %" PRIMME_INT_P " Time %.4e Flag %d stage %d\n",
                  *numLocked-1, lockedSvals[*numLocked-1],
                  lockedNorms[*numLocked-1], primme_svds->stats.numMatvecs,
                  primme_svds->stats.elapsedTime, lockedFlags[*numLocked-1],
                  *stage+1);
         }
         break;
      }
      fflush(primme_svds->outputFile);
   }
   *err = 0;
}


/*******************************************************************************
 * Subroutine monitor_single_stage - report iterations, #MV, residual norm,
 *    eigenvalues, etc. at every inner/outer iteration and when some pair
 *    converges.       
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * basisEvals   The approximate eigenvalues of the basis
 * basisSize    The size of the basis
 * basisFlags   The state of every approximate pair of the basis (see conv_flags)
 * iblock       Indices of the approximate pairs in the block
 * blockSize    The size of the block
 * basisNorms   The approximate residual norms of the pairs of the basis
 * numConverged The number of pairs converged in the basis and the locked pairs
 *              (this value isn't monotonic!)
 * lockedEvals  The locked eigenvalues
 * numLocked    The number of pairs locked
 * lockedFlags  The state of each locked eigenpair (see conv_flags)
 * lockedNorms  The residual norms of the locked pairs
 * inner_its    The number of performed QMR iterations in the current correction equation
 * LSRes        The residual norm of the linear system at the current QMR iteration
 * event        The event reported
 * primme       Structure containing various solver parameters and statistics
 *
 * OUTPUT
 * ------
 * err          Error code
 * 
 ******************************************************************************/

static void monitor_single_stage(void *basisEvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedEvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, primme_event *event, primme_params *primme,
      int *err)
{
   int i;
   REAL *basisEvals = (REAL*)basisEvals_, *basisNorms = (REAL*)basisNorms_,
        *lockedEvals = (REAL*)lockedEvals_, *lockedNorms = (REAL*)lockedNorms_,
        *LSRes = (REAL*)LSRes_;
   assert(event != NULL && primme != NULL);

   REAL basisSvals[basisEvals&&basisSize?*basisSize:0],
        basisSVNorms[basisEvals&&basisSize?*basisSize:0],
        lockedSvals[lockedEvals&&numLocked?*numLocked:0],
        lockedSVNorms[lockedEvals&&numLocked?*numLocked:0];
   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;

   if (primme_svds->method == primme_svds_op_AtA
         || primme_svds->method == primme_svds_op_AAt) {
      /* sval = sqrt(abs(eval)) and SVrnorm = rnorm/sval */

      if (basisEvals && basisSize) for (i=0; i<*basisSize; i++) {
         basisSvals[i] = sqrt(fabs(basisEvals[i]));
         basisSVNorms[i] = basisNorms[i]/basisSvals[i];
      }

      if (lockedEvals && numLocked) for (i=0; i<*numLocked; i++) {
         lockedSvals[i] = sqrt(fabs(lockedEvals[i]));
         lockedSVNorms[i] = lockedNorms[i]/lockedSvals[i];
      }
   }
   else if (primme_svds->method == primme_svds_op_augmented) {
      /* SVrnorm = rnorm/sqrt(2) */

      if (basisEvals && basisSize) for (i=0; i<*basisSize; i++) {
         basisSVNorms[i] = basisNorms[i]/sqrt(2.0);
      }

      if (lockedEvals && numLocked) for (i=0; i<*numLocked; i++) {
         lockedSVNorms[i] = lockedNorms[i]/sqrt(2.0);
      }
   }

   /* When two stages, set primme_event_locked as primme_event_converged */

   primme_event event_svds = *event;
   if (primme_svds->methodStage2 != primme_svds_op_none
         && event_svds == primme_event_locked) {
      event_svds = primme_event_converged;
   }

   /* Record performance measurements */ 

   primme_svds_stats stats = primme_svds->stats;
   UPDATE_STATS(primme_svds->stats, +=, primme->stats);

   /* Call the user function report */

   int ZERO = 0;
   primme_svds->monitorFun(basisSvals, basisSize, basisFlags, iblock, blockSize,
         basisSVNorms, numConverged, lockedSvals, numLocked, lockedFlags,
         lockedSVNorms, inner_its, LSRes, &event_svds, &ZERO, primme_svds, err);
   primme_svds->stats = stats; /* restore original values */
}


/*******************************************************************************
 * Subroutine monitor_stage1 - translate monitored information from eigenvalues
 *    to singular values and call the monitor in primme_svds. Notice that
 *    because there is a second stage, the locked pairs at this stage are
 *    reported as converged. 
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * basisEvals   The approximate eigenvalues of the basis
 * basisSize    The size of the basis
 * basisFlags   The state of every approximate pair of the basis (see conv_flags)
 * iblock       Indices of the approximate pairs in the block
 * blockSize    The size of the block
 * basisNorms   The approximate residual norms of the pairs of the basis
 * numConverged The number of pairs converged in the basis and the locked pairs
 *              (this value isn't monotonic!)
 * lockedEvals  The locked eigenvalues
 * numLocked    The number of pairs locked
 * lockedFlags  The state of each locked eigenpair (see conv_flags)
 * lockedNorms  The residual norms of the locked pairs
 * inner_its    The number of performed QMR iterations in the current correction equation
 * LSRes        The residual norm of the linear system at the current QMR iteration
 * event        The event reported
 * primme       Structure containing various solver parameters and statistics
 *
 * OUTPUT
 * ------
 * err          Error code
 * 
 ******************************************************************************/

static void monitor_stage1(void *basisEvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedEvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, primme_event *event, primme_params *primme,
      int *err)
{
   REAL *basisEvals = (REAL*)basisEvals_, *basisNorms = (REAL*)basisNorms_,
        *lockedEvals = (REAL*)lockedEvals_, *lockedNorms = (REAL*)lockedNorms_,
        *LSRes = (REAL*)LSRes_;
   assert(event != NULL && primme != NULL);

   /* Ignore the converged events if locking is active and printLevel <= 4 */

   if (*event == primme_event_converged && primme->locking
         && primme->printLevel <= 4) {
      *err = 0;
      return;
   }

   /* Show locked pairs as converged pairs of the basis */

   int numLocked0 = lockedEvals&&numLocked?*numLocked:0;
   int basisSize0 = (basisEvals&&basisSize?*basisSize:0) + numLocked0;
   REAL basisSvals[basisSize0], basisSVNorms[basisSize0];
   int basisSVFlags[basisSize0], iblockSV[blockSize?*blockSize:1];
   int numConvergedSV = (numConverged?*numConverged:numLocked0);

   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;

   assert(primme_svds->method == primme_svds_op_AtA
         || primme_svds->method == primme_svds_op_AAt);

   /* sval = sqrt(abs(eval)) and SVrnorm = rnorm/sval */

   int i, j=0;
   if (lockedEvals && numLocked) for (i=0; i<*numLocked; i++, j++) {
      basisSvals[j] = sqrt(fabs(lockedEvals[i]));
      basisSVNorms[j] = lockedNorms[i]/basisSvals[i];
      basisSVFlags[j] = lockedFlags[i];
   }

   if (basisEvals && basisSize) for (i=0; i<*basisSize; i++, j++) {
      basisSvals[j] = sqrt(fabs(basisEvals[i]));
      basisSVNorms[j] = basisNorms[i]/basisSvals[i];
      basisSVFlags[j] = basisFlags ? basisFlags[i] : UNCONVERGED;
   }

   if (iblock && blockSize) for (i=0; i<*blockSize; i++) {
      iblockSV[i] = iblock[i] + numLocked0;
   }

   primme_event eventSV = *event;
   if (*event == primme_event_locked) {
      eventSV = primme_event_converged;
      iblockSV[0] = *numLocked-1;
   }
  
   /* Record performance measurements */ 

   primme_svds_stats stats = primme_svds->stats;
   UPDATE_STATS(primme_svds->stats, +=, primme->stats);

   /* Call the user function report */

   int ZERO = 0;
   primme_svds->monitorFun(basisSvals, &basisSize0, basisSVFlags, iblockSV,
         blockSize, basisSVNorms, &numConvergedSV, NULL, NULL, NULL, NULL,
         inner_its, LSRes, &eventSV, &ZERO, primme_svds, err);
   primme_svds->stats = stats; /* restore original values */
}

/*******************************************************************************
 * Subroutine monitor_stage2 - report iterations, #MV, residual norm,
 *    eigenvalues, etc. at every inner/outer iteration and when some pair
 *    converges.       
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * basisEvals   The approximate eigenvalues of the basis
 * basisSize    The size of the basis
 * basisFlags   The state of every approximate pair of the basis (see conv_flags)
 * iblock       Indices of the approximate pairs in the block
 * blockSize    The size of the block
 * basisNorms   The approximate residual norms of the pairs of the basis
 * numConverged The number of pairs converged in the basis and the locked pairs
 *              (this value isn't monotonic!)
 * lockedEvals  The locked eigenvalues
 * numLocked    The number of pairs locked
 * lockedFlags  The state of each locked eigenpair (see conv_flags)
 * lockedNorms  The residual norms of the locked pairs
 * inner_its    The number of performed QMR iterations in the current correction equation
 * LSRes        The residual norm of the linear system at the current QMR iteration
 * event        The event reported
 * primme       Structure containing various solver parameters and statistics
 *
 * OUTPUT
 * ------
 * err          Error code
 * 
 ******************************************************************************/

static void monitor_stage2(void *basisEvals_, int *basisSize, int *basisFlags,
      int *iblock, int *blockSize, void *basisNorms_, int *numConverged,
      void *lockedEvals_, int *numLocked, int *lockedFlags, void *lockedNorms_,
      int *inner_its, void *LSRes_, primme_event *event, primme_params *primme,
      int *err)
{
   REAL *basisEvals = (REAL*)basisEvals_, *basisNorms = (REAL*)basisNorms_,
        *lockedEvals = (REAL*)lockedEvals_, *lockedNorms = (REAL*)lockedNorms_,
        *LSRes = (REAL*)LSRes_;
   assert(event != NULL && primme != NULL);
   primme_svds_params *primme_svds = (primme_svds_params *) primme->matrix;

   /* Included the converged triplets after the first stage as locked */

   int numLockedExtra = lockedEvals&&numLocked ?
      primme_svds->numSvals - primme->numEvals : 0;
   int numLockedSV = (lockedEvals&&numLocked?*numLocked:0) + numLockedExtra;
   int basisSize0 = (basisEvals&&basisSize?*basisSize:0);
   REAL basisSVNorms[basisSize0],
        lockedSVNorms[numLockedSV];
   int lockedSVFlags[numLockedSV];


   /* SVrnorm = rnorm/sqrt(2) */

   int i;
   if (basisEvals && basisSize) for (i=0; i<*basisSize; i++) {
      basisSVNorms[i] = basisNorms[i]/sqrt(2.0);
   }

   lockedEvals -= numLockedExtra;
   lockedNorms -= numLockedExtra;

   for (i=0; i<numLockedExtra; i++) {
      lockedSVNorms[i] = lockedNorms[i];
      lockedSVFlags[i] = CONVERGED;
   }

   for (i=numLockedExtra; i<numLockedSV; i++) {
      lockedSVNorms[i] = lockedNorms[i]/sqrt(2.0);
      lockedSVFlags[i] = lockedFlags[i-numLockedExtra];
   }

   /* Record performance measurements */ 

   primme_svds_stats stats = primme_svds->stats;
   UPDATE_STATS(primme_svds->stats, +=, primme->stats);

   /* Call the user function report */

   int ONE = 1;
   primme_svds->monitorFun(basisEvals, basisSize, basisFlags, iblock, blockSize,
      basisSVNorms, numConverged, lockedEvals, &numLockedSV, lockedSVFlags,
      lockedSVNorms, inner_its, LSRes, event, &ONE, primme_svds, err);
   primme_svds->stats = stats; /* restore original values */
}
