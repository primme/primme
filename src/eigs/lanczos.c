/*******************************************************************************
 * Copyright (c) 2018, College of William & Mary
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
 * File: main_iter.c
 *
 * Purpose - This is the main Davidson-type iteration 
 *
 ******************************************************************************/

#ifndef THIS_FILE
#define THIS_FILE "../eigs/lanczos.c"
#endif


#include "numerical.h"
#include "template_normal.h"
#include "common_eigs.h"
/* Keep automatically generated headers under this section  */
#ifndef CHECK_TEMPLATE
#include "lanczos.h"
#include "convergence.h"
#include "correction.h"
#include "init.h"
#include "ortho.h"
#include "restart.h"
#include "solve_projection.h"
#include "update_projection.h"
#include "update_W.h"
#include "auxiliary_eigs.h"
#include "auxiliary_eigs_normal.h"
#endif

/******************************************************************************
 * Subroutine lanczos - This routine implements a more general, parallel, 
 *    block (Jacobi)-Davidson outer iteration with a variety of options.
 *
 * INPUT/OUTPUT arrays and parameters
 * ----------------------------------
 * evecs    Stores initial guesses. Upon return, it contains the converged Ritz
 *          vectors.  If locking is not engaged, then converged Ritz vectors 
 *          are copied to this array just before return.  
 *
 * primme.initSize: On output, it stores the number of converged eigenvectors. 
 *           If smaller than numEvals and locking is used, there are
 *              only primme.initSize vectors in evecs.
 *           Without locking all numEvals approximations are in evecs
 *              but only the initSize ones are converged.
 *           During the execution, access to primme.initSize gives 
 *              the number of converged pairs up to that point. The pairs
 *              are available in evals and evecs, but only when locking is used
 * ret      successful state
 * numRet   number of returned pairs in evals, evecs, and resNorms
 *
 * Return Value
 * ------------
 * error code      
 ******************************************************************************/

TEMPLATE_PLEASE
int lanczos_Sprimme(HEVAL *evals, SCALAR *evecs, PRIMME_INT ldevecs,
      HREAL *resNorms, double startTime, int *ret, int *numRet,
      primme_context ctx) {

   primme_params *primme = ctx.primme;
   PRIMME_INT nLocal = primme->nLocal;

   // Write first random vector
   Num_larnv_Sprimme(2, primme->iseed, nLocal, &evecs[0], ctx);

   // inner_prod = evecs[0]'*evecs[0]
   HSCALAR inner_prod = Num_dot_Sprimme(nLocal, &evecs[0], 1, &evecs[0], 1, ctx);

   // evecs[0] /= inner_prod
   CHKERR(Num_scal_Sprimme(nLocal, 1.0 / inner_prod, &evecs[0], 1, ctx));

   // evecs[1] = A * evecs[0]
   CHKERR(matrixMatvec_Sprimme(
         &evecs[0], nLocal, ldevecs, &evecs[1], ldevecs, 0, 1, ctx));

   // Allocate an array on the CPU
   HSCALAR *overlaps;
   int numEvals = primme->numEvals;
   CHKERR(Num_malloc_SHprimme(numEvals * numEvals, &overlaps, ctx));
 }
