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
 * File: ortho.c
 *
 * Purpose - Orthonormalizes a block of vectors, vector by vector, 
 *           against two bases and among themselves. Gram-Schmidt is used 
 *           with reorthogonalization based on Daniel's test. 
 *           For the purpose of the test, the norm of the resulting vector 
 *           is computed without synchronizations. Because of floating point,
 *           this norm is not accurate if reortho is needed, but it is more 
 *           than sufficient for the test. If at least one reortho has been
 *           performed, the actual norm must be computed. 
 *
 * Note on numerical properties
 *           The relative error in the implicit s1 is 
 *                    O( machineEps * (s0/s1)^2 ), 
 *           so when Daniel's test succeeds error(s1) = O(2*machineEps)).
 *
 *           For almost linearly dependent vectors a large error in s1
 *           might fail to detect a problem during the second reortho.
 *           Let s1 = s1hat*(1+err), with err = O(macheps*(s0/s1)^2).
 *              eg. s0/s1 = 1e8 => (s1-s1hat)/s1hat = O(1)
 *              eg. s0/s1 = 1e7 => (s1-s1hat)/s1hat = O(1e-2)
 *              eg. s0/s1 = 1e4 => (s1-s1hat)/s1hat = O(1e-8)
 *           When s0 is reset to s00 = s1, and perform again the Daniel's test
 *              s11/s00 = s11/s1hat(1+err) < 0.75 <=> s11/s1hat < .75*(1+err)
 *           if err is around 1e-2 it does not affect Daniels test a lot.
 *           Thus, if s0/s1 > 1.49e7, we must compute the actual s1.
 *
 * Note on linear dependency
 *           If the original vector norm is reduced more than machEps times,
 *           it is considered linear dependent. If R is not returned, the
 *           vector is replaced by a random vector. Otherwise the vector is
 *           zeroed.
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "numerical.h"
#include "ortho.h"
#include "const.h"
#include "globalsum.h"
#include "wtime.h"
 

/**********************************************************************
 * Function ortho - This routine orthonormalizes
 * a block of of vectors (from b1 to including b2 in basis)
 * against other vectors in the same array (from 0 to b1-1 in basis),
 * against a set of locked vectors (from 0 to numLocked-1 in locked),
 * and themselves.
 *
 * The following conditions must always be met: 
 * ldBasis > 0, nLocal > 0, b1 >= 0, b2 >= 0, b2 >= b1, numLocked >= 0, 
 * rworkSize > 0
 *
 * INPUT ARRAYS AND PARAMETERS
 * ---------------------------
 * ldBasis    Leading dimension of the basis
 * b1, b2     Range of indices of vectors to be orthonormalized
 *            (b1 can be zero, but b1 must be <= b2)
 * ldR        Leading dimension in R
 * locked     Array that holds locked vectors if they are in-core
 * ldLocked   Leading dimension of locked
 * numLocked  Number of vectors in locked
 * nLocal     Number of rows of each vector stored on this node
 * machEps    Double machine precision
 *
 * rworkSize  Length of rwork array
 * primme     Primme struct. Contains globalSumDouble and Parallelism info
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * basis   Basis vectors stored in core memory
 * R       Rotations done in the basis: input_basis = output_basis * R
 * iseed   Seeds used to generate random vectors
 * rwork   Contains buffers and other necessary work arrays
 *
 * Return Value
 * ------------
 * >0  - Insufficient workspace provided.  Positive number indicates
 *       appropriate amount.
 *  0  - success
 * -1  - some size or leading dimension was < 0
 * -2  - b1 > b2
 * -3  - A limit number of randomizations has been performed without
 *       yielding an orthogonal direction
 * 
 **********************************************************************/

TEMPLATE_PLEASE
int ortho_Sprimme(SCALAR *basis, PRIMME_INT ldBasis, SCALAR *R,
      PRIMME_INT ldR, int b1, int b2, SCALAR *locked, PRIMME_INT ldLocked,
      int numLocked, PRIMME_INT nLocal, PRIMME_INT *iseed, double machEps,
      SCALAR *rwork, size_t *rworkSize, primme_params *primme) {
              
   int i, j;                /* Loop indices */
   size_t minWorkSize;         
   int nOrth, reorth;
   int randomizations;
   int updateR;             /* update factor R */
   int messages = 1;        /* messages = 1 prints the intermediate results */
   /* TODO: replace by a dynamic criterion when to stop orthogonalizing local */
   /* vectors. Observed performance improvement when maxNumOrthos increases.  */
   int maxNumOrthos = primme?3:7; /* We let 2 reorthogonalizations before randomize */
                                  /* for nLocal length vectors, and 6 orthogonalisations */
                                  /* for the rest */
   int maxNumRandoms = 10;  /* We do not allow more than 10 randomizations */
   double tol = sqrt(2.0L)/2.0L; /* We set Daniel et al. test to .707 */
   REAL s0=0.0, s02=0.0, s1=0.0, s12=0.0;
   REAL temp;
   SCALAR *overlaps;
   double t0;

   messages = (primme && primme->procID == 0 && primme->printLevel >= 3
         && primme->outputFile);

   minWorkSize = 2*(numLocked + b2 + 1);

   /* Return memory requirement */
   if (basis == NULL) {
      *rworkSize = max(*rworkSize, minWorkSize);
      return 0;
   }

   /*----------------------------------*/
   /* input and workspace verification */
   /*----------------------------------*/
   assert(nLocal >= 0 && numLocked >= 0 && *rworkSize >= minWorkSize &&
          ldBasis >= nLocal && (numLocked == 0 || ldLocked >= nLocal) &&
          (R == NULL || ldR > b2));

   tol = sqrt(2.0L)/2.0L;

   /* Zero the columns from b1 to b2 of R */
   if (R)
      for(i=b1; i <= b2; i++)
         for (j=0; j <= i; j++)
            R[ldR*i+j] = 0.0;

   /*---------------------------------------------------*/
   /* main loop to orthogonalize new vectors one by one */
   /*---------------------------------------------------*/

   t0 = primme_wTimer(0);

   for(i=b1; i <= b2; i++) {
    
      nOrth = 0;
      reorth = 1;
      randomizations = 0;
      updateR = (R ? 1 : 0);
      s1 = 0.0;

      while (reorth) {

         if (nOrth >= maxNumOrthos) {
            /* Stop updating R when replacing one of the columns of the basis */
            /* with a random vector                                           */

            if (updateR) {
               R[ldR*i + i] = 0.0;
               updateR = 0;
            }

            if (randomizations >= maxNumRandoms) {
               return -3;
            }
            if (messages){
               fprintf(primme->outputFile, "Randomizing in ortho: %d, vector size of %" PRIMME_INT_P "\n", i, nLocal);
            }

            Num_larnv_Sprimme(2, iseed, nLocal, &basis[ldBasis*i]); 
            randomizations++;
            nOrth = 0;
         }

         nOrth++;

         if (nOrth == 1) {
            s02 = REAL_PART(Num_dot_Sprimme(nLocal, &basis[ldBasis*i], 1, 
                     &basis[ldBasis*i], 1));
            if (primme) primme->stats.numOrthoInnerProds += 1;
         }
            
         if (i > 0) {
            Num_gemv_Sprimme("C", nLocal, i, 1.0, basis, ldBasis, 
               &basis[ldBasis*i], 1, 0.0, rwork, 1);
            if (primme) primme->stats.numOrthoInnerProds += i;
         }

         if (numLocked > 0) {
            Num_gemv_Sprimme("C", nLocal, numLocked, 1.0, locked, ldLocked,
               &basis[ldBasis*i], 1, 0.0, &rwork[i], 1);
            if (primme) primme->stats.numOrthoInnerProds += numLocked;
         }

         rwork[i+numLocked] = s02;
         overlaps = &rwork[i+numLocked+1];
         CHKERR(globalSum_Sprimme(rwork, overlaps, i + numLocked + 1,
                  primme), -1);

         if (updateR) {
             Num_axpy_Sprimme(i, 1.0, overlaps, 1, &R[ldR*i], 1);
         }

         if (numLocked > 0) { /* locked array most recently accessed */
            Num_gemv_Sprimme("N", nLocal, numLocked, -1.0, locked, ldLocked, 
               &overlaps[i], 1, 1.0, &basis[ldBasis*i], 1); 
            if (primme) primme->stats.numOrthoInnerProds += numLocked;
         }

         if (i > 0) {
            Num_gemv_Sprimme("N", nLocal, i, -1.0, basis, ldBasis, 
               overlaps, 1, 1.0, &basis[ldBasis*i], 1);
            if (primme) primme->stats.numOrthoInnerProds += i;
         }
 
         if (nOrth == 1) {
            s0 = sqrt(s02 = REAL_PART(overlaps[i+numLocked]));
         }

         /* Compute the norm of the resulting vector implicitly */
         
         temp = REAL_PART(Num_dot_Sprimme(i+numLocked,overlaps,1,overlaps,1));
         s1 = sqrt(s12 = max(0.0L, s02-temp));
         
         /* s1 decreased too much. Numerical problems expected   */
         /* with its implicit computation. Compute s1 explicitly */
         
         if ( s1 < s0*sqrt(machEps) || nOrth > 1 || !primme) {  
            temp = REAL_PART(Num_dot_Sprimme(nLocal, &basis[ldBasis*i], 1, 
                                           &basis[ldBasis*i], 1));
            if (primme) primme->stats.numOrthoInnerProds += 1;
            CHKERR(globalSum_Rprimme(&temp, &s12, 1, primme), -1);
            s1 = sqrt(s12);
         }

         if (s1 <= machEps*s0) {
            if (messages) {
               fprintf(primme->outputFile, 
                 "Vector %d lost all significant digits in ortho\n", i-b1);
            }
            nOrth = maxNumOrthos;
         }
         else if (s1 <= tol*s0 || (!primme && nOrth < maxNumOrthos)) {
            /* No numerical benefit in normalizing the vector before reortho */
            s0 = s1;
            s02 = s12;
         }
         else {
            if (updateR) {
               if (!primme || nOrth == 1) {
                  temp = REAL_PART(Num_dot_Sprimme(nLocal,
                           &basis[ldBasis*i], 1, &basis[ldBasis*i], 1));
                  if (primme) primme->stats.numOrthoInnerProds += 1;
                  CHKERR(globalSum_Rprimme(&temp, &s1, 1, primme), -1);
                  s1 = sqrt(s1);
               }
               R[ldR*i + i] = s1;
            }

            Num_scal_Sprimme(nLocal, 1.0/s1, &basis[ldBasis*i], 1);
            reorth = 0;
         } 
 
      }
   }

   if (primme) primme->stats.timeOrtho += primme_wTimer(0) - t0;

   /* Check orthogonality */
   /*
   if (numLocked) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*numLocked*numLocked);
      Num_gemm_Sprimme("C", "N", numLocked, numLocked, nLocal, 1.0, locked,
            ldLocked, locked, ldLocked, 0.0, H, numLocked);
      for(i=0; i < numLocked; i++) {
         for(j=0; j < i; j++) assert(ABS(H[numLocked*i+j]) < 1e-13);
         assert(fabs(1 - ABS(H[numLocked*i+i])) < 1e-13);
      }
      free(H);
   }
   if (b2+1) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*(b2+1)*(b2+1));
      Num_gemm_Sprimme("C", "N", b2+1, b2+1, nLocal, 1.0, basis,
            ldBasis, basis, ldBasis, 0.0, H, b2+1);
      for(i=0; i < b2+1; i++) {
         for(j=0; j < i; j++) assert(ABS(H[(b2+1)*i+j]) < 1e-13);
         assert(H[(b2+1)*i+i] == 0.0 || fabs(1 - ABS(H[(b2+1)*i+i])) < 1e-13);
      }
      free(H);
   }
   if (numLocked) {
      SCALAR *H = (SCALAR*)malloc(sizeof(SCALAR)*(b2+1)*numLocked);
      Num_gemm_Sprimme("C", "N", numLocked, b2+1, nLocal, 1.0, locked,
            ldLocked, basis, ldBasis, 0.0, H, numLocked);
      for(i=0; i < b2+1; i++) {
         for(j=0; j < numLocked; j++) assert(ABS(H[numLocked*i+j]) < 1e-13);
      }
      free(H);
   }
   */

   return 0;
}

/**********************************************************************
 * Function ortho_single_iteration -- This function orthogonalizes
 *    applies ones the projector (I-QQ') on X. Optionally returns
 *    the norms ||Q'X(i)|| and ||(I-QQ')X(i)||.
 *   
 * ARRAYS AND PARAMETERS
 * ----------------
 * Q               The basis of the projector I-QQ'.
 * mQ, nQ, ldQ     Rows, columns and leading dimension of Q.
 * X               The vectors to apply the projector I-QQ'.
 * inX             Column indices to apply the projector (optional).
 * nX              Number of columns to apply the projector.
 * ldX             The leading dimension of X.
 * overlaps        The norms of Q'X(i) (optional).
 * norms           The norms of (I-QQ')X(i) (optional).
 * rwork           Auxiliary space
 * lrwork          Available rwork
 *
 ****************************************************************************/

TEMPLATE_PLEASE
int ortho_single_iteration_Sprimme(SCALAR *Q, PRIMME_INT mQ, PRIMME_INT nQ,
      PRIMME_INT ldQ, SCALAR *X, int *inX, int nX, PRIMME_INT ldX,
      REAL *overlaps, REAL *norms, SCALAR *rwork, size_t *lrwork,
      primme_params *primme) {

   int i, j, M=PRIMME_BLOCK_SIZE, m=min(M, mQ);
   SCALAR *y, *y0, *X0;
   REAL *norms0;

   /* Return memory requirement */
   if (Q == NULL) {
      *lrwork = max(*lrwork, (size_t)nQ*nX*2 + (size_t)M*nX);
      return 0;
   }

   double t0 = primme_wTimer(0);

   assert((size_t)nQ*nX*2 + (size_t)m*nX <= *lrwork);

   /* Warning: norms0 and y overlap, so don't use them at the same time */
   norms0 = (REAL*)rwork;
   y = rwork;
   y0 = y + nQ*nX;
   X0 = y0 + nQ*nX;

   /* Check if the indices of inX are contiguous */

   if (inX && nX > 0) {
      for (i=0, j=inX[0]; i<nX && j==inX[i]; i++, j++);
      if (i >= nX) {
         X = &X[inX[0]*ldX];
         inX = NULL;
      }
   }

   /* y = Q'*X */
   if (!inX) {
      Num_gemm_Sprimme("C", "N", nQ, nX, mQ, 1.0, Q, ldQ, X, ldX, 0.0, y, nQ);
   }
   else {
      Num_zero_matrix_Sprimme(y, nQ, nX, nQ);
      for (i=0, m=min(M,mQ); i < mQ; i+=m, m=min(m,mQ-i)) {
         Num_copy_matrix_columns_Sprimme(&X[i], m, inX, nX, ldX, X0, NULL, m);
         Num_gemm_Sprimme("C", "N", nQ, nX, m, 1.0, &Q[i], ldQ, X0, m, 1.0,
               y, nQ);
      }
   }
   primme->stats.numOrthoInnerProds += nQ*nX;

   /* Store the reduction of y in y0 */
   CHKERR(globalSum_Sprimme(y, y0, nQ*nX, primme), -1);
   
   /* overlaps(i) = norm(y0(:,i))^2 */
   for (i=0; i<nX; i++) {
      overlaps[i] =
         sqrt(REAL_PART(Num_dot_Sprimme(nQ, &y0[nQ*i], 1, &y0[nQ*i], 1)));
   }

   /* X = X - Q*y0; norms0(i) = norms(X(i))^2 */
   if (norms) for (i=0; i<nX; i++) norms0[i] = 0.0;
   for (i=0, m=min(M,mQ); i < mQ; i+=m, m=min(m,mQ-i)) {
      if (inX) {
         Num_copy_matrix_columns_Sprimme(&X[i], m, inX, nX, ldX, X0, NULL, m);
      }
      Num_gemm_Sprimme("N", "N", m, nX, nQ, -1.0, &Q[i], ldQ, y0, nQ, 1.0,
            inX?X0:&X[i], inX?m:ldX);
      if (inX) {
         Num_copy_matrix_columns_Sprimme(X0, m, NULL, nX, m, &X[i], inX, ldX);
      }
      if (norms) for (j=0; j<nX; j++) {
         SCALAR *v = inX ? &X0[j*m] : &X[j*ldX+i];
         norms0[j] += REAL_PART(Num_dot_Sprimme(m, v, 1, v, 1));
      }
   }

   if (norms) {
      /* Store the reduction of norms0 in norms */
      CHKERR(globalSum_Rprimme(norms0, norms, nX, primme), -1);
 
      for (i=0; i<nX; i++) norms[i] = sqrt(norms[i]);
      primme->stats.numOrthoInnerProds += nX;
   }

   primme->stats.timeOrtho += primme_wTimer(0) - t0;

   return 0;
}
