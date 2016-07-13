/*******************************************************************************
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
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
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
#include "primme.h"         
#include "numerical_z.h"
#include "ortho_z.h"
 

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

int ortho_zprimme(Complex_Z *basis, int ldBasis, Complex_Z *R, int ldR,
   int b1, int b2, Complex_Z *locked, int ldLocked, int numLocked,
   int nLocal, int *iseed, double machEps, Complex_Z *rwork, int rworkSize,
   primme_params *primme) {
              
   int i, j;                /* Loop indices */
   int count;
   int minWorkSize;         
   int nOrth, reorth;
   int randomizations;
   int messages = 0;        /* messages = 1 prints the intermediate results */
   int maxNumOrthos = primme?3:5; /* We let 2 reorthogonalizations before randomize */
                                  /* for nLocal length vectors, and 5 orthogonaliations */
                                  /* for the rest */
   int maxNumRandoms = 10;  /* We do not allow more than 10 randomizations */
   double tol = sqrt(2.0L)/2.0L; /* We set Daniel et al. test to .707 */
   double s0=0.0, s02=0.0, s1=0.0, s00=0.0;
   double temp;
   Complex_Z ztmp={+0.0e+00,+0.0e00};
   Complex_Z *overlaps;
   Complex_Z tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00}, tmone = {-1.0e+00,+0.0e00};
   FILE *outputFile;

   messages = (primme && primme->procID == 0 && primme->printLevel >= 3
         && primme->outputFile);
   outputFile = primme ? primme->outputFile : stdout;

   minWorkSize = 2*(numLocked + b2 + 1);

   /* Return memory requirement */
   if (basis == NULL) {
      return minWorkSize;
   }

   /*----------------------------------*/
   /* input and workspace verification */
   /*----------------------------------*/
   assert(nLocal >= 0 && numLocked >= 0 && rworkSize >= minWorkSize &&
          ldBasis >= nLocal && (numLocked == 0 || ldLocked >= nLocal) &&
          (R == NULL || ldR >= b2));

   tol = sqrt(2.0L)/2.0L;

   /* Zero the columns from b1 to b2 of R */
   if (R)
      for(i=b1; i <= b2; i++)
         for (j=0; j <= i; j++)
            R[ldR*i+j] = tzero;

   /*---------------------------------------------------*/
   /* main loop to orthogonalize new vectors one by one */
   /*---------------------------------------------------*/

   for(i=b1; i <= b2; i++) {
    
      nOrth = 0;
      reorth = 1;
      randomizations = 0;

      while (reorth) {

         if (nOrth >= maxNumOrthos) {
            if (randomizations >= maxNumRandoms) {
               return -3;
            }
            if (messages){
               fprintf(outputFile, "Randomizing in ortho: %d, vector size of %d\n", i, nLocal);
            }

            assert(R == NULL);
            Num_larnv_zprimme(2, iseed, nLocal, &basis[ldBasis*i]); 
            randomizations++;
            nOrth = 0;
         }

         nOrth++;

         if (nOrth == 1) {
            ztmp = Num_dot_zprimme(nLocal, &basis[ldBasis*i], 1, 
                                           &basis[ldBasis*i], 1);
         }
            
         if (i > 0) {
            Num_gemv_zprimme("C", nLocal, i, tpone, basis, ldBasis, 
               &basis[ldBasis*i], 1, tzero, rwork, 1);
         }

         if (numLocked > 0) {
            Num_gemv_zprimme("C", nLocal, numLocked, tpone, locked, ldLocked,
               &basis[ldBasis*i], 1, tzero, &rwork[i], 1);
         }

         rwork[i+numLocked] = ztmp;
         overlaps = &rwork[i+numLocked+1];
         /* In Complex, the size of the array to globalSum is twice as large */
         count = 2*(i + numLocked + 1);
         (primme ? primme->globalSumDouble : primme_seq_globalSumDouble)
            (rwork, overlaps, &count, primme);

         if (R != NULL) {
             Num_axpy_zprimme(i, tpone, overlaps, 1, &R[ldR*i], 1);
         }

         if (numLocked > 0) { /* locked array most recently accessed */
            Num_gemv_zprimme("N", nLocal, numLocked, tmone, locked, ldLocked, 
               &overlaps[i], 1, tpone, &basis[ldBasis*i], 1); 
         }

         if (i > 0) {
            Num_gemv_zprimme("N", nLocal, i, tmone, basis, ldBasis, 
               overlaps, 1, tpone, &basis[ldBasis*i], 1);
         }
 
         if (nOrth == 1) {
            s02 = overlaps[i+numLocked].r;
            s00 = s0 = sqrt(s02);
         }

         /* Compute the norm of the resulting vector implicitly */
         
         ztmp = Num_dot_zprimme(i+numLocked,overlaps,1,overlaps,1);
         temp = ztmp.r;
         s1 = sqrt(max(0.0L, s02-temp));
         
         /* s1 decreased too much. Numerical problems expected   */
         /* with its implicit computation. Compute s1 explicitly */
         
         if ( s1 < s0*sqrt(machEps) || nOrth > 1 || !primme) {  
            ztmp = Num_dot_zprimme(nLocal, &basis[ldBasis*i], 1, 
                                           &basis[ldBasis*i], 1);
            temp = ztmp.r;
            count = 1;
            (primme ? primme->globalSumDouble : primme_seq_globalSumDouble)
               (&temp, &s1, &count, primme);
            s1 = sqrt(s1);
         }

         if (R && (s1 <= machEps*s00 || (s1 <= tol*s0 && nOrth >= maxNumOrthos))) {
            if (messages) {
               fprintf(outputFile, "Zeroing column %d\n", i);
            }
            /* No randomization when computing the QR decomposition */
            Num_scal_zprimme(nLocal, tzero, &basis[ldBasis*i], 1);
            R[ldR*i + i] = tzero;
            reorth = 0;
         }
         else if (s1 <= machEps*s00) {
            if (messages) {
               fprintf(outputFile, 
                 "Vector %d lost all significant digits in ortho\n", i-b1);
            }
            nOrth = maxNumOrthos;
         }
         else if (s1 <= tol*s0 || (!primme && nOrth < maxNumOrthos)) {
            /* No numerical benefit in normalizing the vector before reortho */
            s0 = s1;
            s02 = s1*s1;
         }
         else {
            if (R != NULL) {
               if (!primme || nOrth == 1) {
                  ztmp = Num_dot_zprimme(nLocal, &basis[ldBasis*i], 1,
                        &basis[ldBasis*i], 1);   
                  temp = ztmp.r;
                  count = 1;
                  (primme ? primme->globalSumDouble : primme_seq_globalSumDouble)
                     (&temp, &s1, &count, primme);
                  s1 = sqrt(s1);
               }
               R[ldR*i + i].r = s1;
               R[ldR*i + i].i = 0.0L;
            }

            {ztmp.r = 1.0L/s1; ztmp.i = 0.0L;}
            Num_scal_zprimme(nLocal, ztmp, &basis[ldBasis*i], 1);
            reorth = 0;
         } 
 
      }
   }

   /* Check orthogonality */
   /*
   if (numLocked) {
      Complex_Z *H = (Complex_Z*)malloc(sizeof(Complex_Z)*numLocked*numLocked);
      Num_gemm_zprimme("C", "N", numLocked, numLocked, nLocal, tpone, locked,
            ldLocked, locked, ldLocked, tzero, H, numLocked);
      for(i=0; i < numLocked; i++) {
         for(j=0; j < i; j++) assert(fabs(*(double*)&H[numLocked*i+j]) < 1e-13);
         assert(fabs(1 - *(double*)&H[numLocked*i+i]) < 1e-13);
      }
      free(H);
   }
   if (b2+1) {
      Complex_Z *H = (Complex_Z*)malloc(sizeof(Complex_Z)*(b2+1)*(b2+1));
      Num_gemm_zprimme("C", "N", b2+1, b2+1, nLocal, tpone, basis,
            ldBasis, basis, ldBasis, tzero, H, b2+1);
      for(i=0; i < b2+1; i++) {
         for(j=0; j < i; j++) assert(fabs(*(double*)&H[(b2+1)*i+j]) < 1e-13);
         assert(*(double*)&H[(b2+1)*i+i] == 0.0 || fabs(1 - *(double*)&H[(b2+1)*i+i]) < 1e-13);
      }
      free(H);
   }
   if (numLocked) {
      Complex_Z *H = (Complex_Z*)malloc(sizeof(Complex_Z)*(b2+1)*numLocked);
      Num_gemm_zprimme("C", "N", numLocked, b2+1, nLocal, tpone, locked,
            ldLocked, basis, ldBasis, tzero, H, numLocked);
      for(i=0; i < b2+1; i++) {
         for(j=0; j < numLocked; j++) assert(fabs(*(double*)&H[numLocked*i+j]) < 1e-13);
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

int ortho_single_iteration_zprimme(Complex_Z *Q, int mQ, int nQ, int ldQ, Complex_Z *X,
   int *inX, int nX, int ldX, double *overlaps, double *norms, Complex_Z *rwork, int lrwork,
   primme_params *primme) {

   int i, j, M=PRIMME_BLOCK_SIZE, m=min(M, mQ), count;
   Complex_Z *y, *y0, *X0;
   double *norms0;
   Complex_Z tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00}, tmone = {-1.0e+00,+0.0e00};

   /* Return memory requirement */
   if (Q == NULL) {
      return nQ*nX*2 + M*nX;
   }

   assert(nQ*nX*2 + m*nX <= lrwork);

   /* Warning: norms0 and y overlap, so don't use them at the same time */
   norms0 = (double*)rwork;
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
      Num_gemm_zprimme("C", "N", nQ, nX, mQ, tpone, Q, ldQ, X, ldX, tzero, y, nQ);
   }
   else {
      for (i=0; i<nQ*nX; i++)
         y[i] = tzero;
      for (i=0, m=min(M,mQ); i < mQ; i+=m, m=min(m,mQ-i)) {
         Num_copy_matrix_columns_zprimme(&X[i], m, inX, nX, ldX, X0, NULL, m);
         Num_gemm_zprimme("C", "N", nQ, nX, m, tpone, &Q[i], ldQ, X0, m, tpone,
               y, nQ);
      }
   }

   /* Store the reduction of y in y0 */
   count = nQ*nX*(sizeof(Complex_Z)/sizeof(double));
   primme->globalSumDouble(y, y0, &count, primme);
   
   /* overlaps(i) = norm(y0(:,i))^2 */
   for (i=0; i<nX; i++) {
      Complex_Z ztmp = Num_dot_zprimme(nQ, &y0[nQ*i], 1, &y0[nQ*i], 1);
      overlaps[i] = sqrt(*(double*)&ztmp);
   }

   /* X = X - Q*y0; norms0(i) = norms(X(i))^2 */
   if (norms) for (i=0; i<nX; i++) norms0[i] = 0.0;
   for (i=0, m=min(M,mQ); i < mQ; i+=m, m=min(m,mQ-i)) {
      Num_gemm_zprimme("N", "N", m, nX, nQ, tmone, &Q[i], ldQ, y0, nQ, tpone,
            inX?X0:&X[i], inX?m:ldX);
      if (inX) {
         Num_copy_matrix_columns_zprimme(X0, m, NULL, nX, ldX, &X[i], inX, ldX);
      }
      if (norms) for (j=0; j<nX; j++) {
         Complex_Z *v = inX ? &X0[j*m] : &X[j*ldX+i];
         Complex_Z ztmp = Num_dot_zprimme(m, v, 1, v, 1);
         norms0[j] += *(double*)&ztmp;
      }
   }

   if (norms) {
      /* Store the reduction of norms0 in norms */
      primme->globalSumDouble(norms0, norms, &nX, primme);
 
      for (i=0; i<nX; i++) norms[i] = sqrt(norms[i]);
   }

   return 0;
}
