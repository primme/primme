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
 *           against two bases and among themselves. Gram-Scmidt is used 
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
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
 * locked     Array that holds locked vectors if they are in-core
 * ldLocked   Leading dimension of locked
 * numLocked  Number of vectors in locked
 * nLocal     Number of rows of each vector stored on this node
 * machEps    Double machine precision
 *
 * lockFile   Handle of file containing the locked vectors
 * rworkSize  Length of rwork array
 * primme     Primme struct. Contains globalSumDouble and Parallelism info
 *
 * INPUT/OUTPUT ARRAYS AND PARAMETERS
 * ----------------------------------
 * basis   Basis vectors stored in core memory
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

int ortho_zprimme(Complex_Z *basis, int ldBasis, int b1, int b2, 
   Complex_Z *locked, int ldLocked, int numLocked, int nLocal, int *iseed, 
   double machEps, Complex_Z *rwork, int rworkSize, primme_params *primme) {
              
   int i;                   /* Loop indices */
   int count;
   int returnValue;
   int minWorkSize;         
   int nOrth, reorth;
   int randomizations;
   int messages = 0;        /* messages = 1 prints the intermediate results */
   int maxNumOrthos = 2;    /* We let 2 reorthogonalizations before randomize */
   int maxNumRandoms = 10;  /* We do not allow more than 10 randomizations */
   double tol = sqrt(2.0L)/2.0L; /* We set Daniel et al. test to .707 */
   double s0, s02, s1;
   double temp;
   Complex_Z ztmp;
   Complex_Z *overlaps;
   Complex_Z tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00}, tmone = {-1.0e+00,+0.0e00};
   FILE *outputFile;

   /* messages = (primme->procID == 0 && primme->printLevel >= 5); */
   /* outputFile = primme->outputFile; */
   outputFile = stderr;

   returnValue = 0;

   /*----------------------------------*/
   /* input and workspace verification */
   /*----------------------------------*/
   if (ldBasis <= 0 || nLocal <= 0 || b1 < 0 || b2 < 0 || numLocked < 0
      || rworkSize < 0)
   {
      returnValue = -1;
   }
   else if (b1 > b2) {
      returnValue = -2;
   }

   if (returnValue != 0) {
      return(returnValue);
   }

   minWorkSize = 2*(numLocked + b2 + 1);

   if (rworkSize < minWorkSize) {
      return(minWorkSize);
   }
   
   tol = sqrt(2.0L)/2.0L;

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
               fprintf(outputFile, "Randomizing in ortho:\n");
            }

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
         (*primme->globalSumDouble)(rwork, overlaps, &count, primme);

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
            s0 = sqrt(s02);
         }

         /* Compute the norm of the resulting vector implicitly */
         
         ztmp = Num_dot_zprimme(i+numLocked,overlaps,1,overlaps,1);
         temp = ztmp.r;
         s1 = sqrt(max(0.0L, s02-temp));
         
         /* s1 decreased too much. Numerical problems expected   */
         /* with its implicit computation. Compute s1 explicitly */
         
         if ( s1 < s0*sqrt(machEps) || nOrth > 1) {  
            ztmp = Num_dot_zprimme(nLocal, &basis[ldBasis*i], 1, 
                                           &basis[ldBasis*i], 1);
            temp = ztmp.r;
            count = 1;
            (*primme->globalSumDouble)(&temp, &s1, &count, primme);
            s1 = sqrt(s1);
         }

         if (s1 <= machEps*s0) {
            if (messages) {
               fprintf(outputFile, 
                 "Vector %d lost all significant digits in ortho\n", i-b1);
            }
            nOrth = maxNumOrthos;
         }
         else if (s1 <= tol*s0) {
            if (messages) {
               fprintf(outputFile, "Reorthogonalizing: %d\n", i-b1);
            }
            /* No numerical benefit in normalizing the vector before reortho */
            s0 = s1;
            s02 = s1*s1;
         }
         else {
            {ztmp.r = 1.0L/s1; ztmp.i = 0.0L;}
            Num_scal_zprimme(nLocal, ztmp, &basis[ldBasis*i], 1);
            reorth = 0;
         } 
            
      }
   }
         
   return 0;
}


/**********************************************************************
 * Function ortho_retained_vectors -- This function orthogonalizes
 *   coefficient vectors (the eigenvectors of the projection H) that
 *   are retained from a previous iteration.  The retained coefficient
 *   vectors are orthogonalized versus the coefficient vectors of the
 *   restarted basis.  This orthogonalizes previous and current Ritz 
 *   vectors cheaply and without communication because the coefficient 
 *   vectors are stored by each process.
 *   
 * Input parameters
 * ----------------
 * currentVectors  Coefficient vectors from the current iteration to 
 *    be orthogonalized against.
 *
 * length          The dimension of each vector
 *
 * numVectors      The number of current vectors
 * 
 * 
 * Input/Output parameters
 * -----------------------
 * previousVectors The coefficient vectors obtained from the previous iteration
 *                 Note that their leading dimension is primme->maxBasisSize.
 *
 * numPrevious     The number of previous vectors to retain
 *
 * rwork           work array of size max(numVectors, numPrevious).  If it
 *                 is of size maxBasisSize, then it will always be sufficiently
 *                 large.
 *
 *
 * Return value
 * ------------
 * int   The number of previous vectors successfully orthogonalized
 *
 ****************************************************************************/

int ortho_retained_vectors_zprimme (Complex_Z *currentVectors, 
  int length, int numVectors, Complex_Z *previousVectors, int numPrevious, 
  double machEps, Complex_Z *rwork) {

   int i;       /* Loop counter                                     */
   int nOrths;  /* Number of times a vector has been orthogonalized */
   int zeroed;  /* True if the vector norm was reduced below 1e-14  */
   double norm; /* Vector norm.                                     */
   Complex_Z ztmp;/* Temp accumulation var                            */
                /* and some constants                               */
   Complex_Z tpone = {+1.0e+00,+0.0e00}, tzero = {+0.0e+00,+0.0e00}, tmone = {-1.0e+00,+0.0e00};

   /* Orthogonalize each of the numPrevious vectors against the current */
   /* vectors and amongst themselves.                                   */

   i = 0;  

   while (i < numPrevious) {
      zeroed = 0;

      /* Orthogonalize each vector twice to ensure numerical stability */

      for (nOrths = 0; nOrths < 3; nOrths++) {

         /* Orthogonalize versus the numVectors current vectors */

         Num_gemv_zprimme("C", length, numVectors, tpone, currentVectors, 
            length, &previousVectors[length*i], 1, tzero, rwork, 1);

         Num_gemv_zprimme("N", length, numVectors, tmone, currentVectors, 
            length, rwork, 1, tpone, &previousVectors[length*i], 1);

         /* Orthogonalize against the i previous vectors that have */
         /* been orthogonalized thus far.                          */

         if (i > 0) {
            Num_gemv_zprimme("C", length, i, tpone, previousVectors, 
               length, &previousVectors[length*i], 1, tzero, rwork, 1);

            Num_gemv_zprimme("N", length, i, tmone, previousVectors, 
               length, rwork, 1, tpone, &previousVectors[length*i], 1);
         }

         ztmp = Num_dot_zprimme(length, &previousVectors[length*i], 1,
                                        &previousVectors[length*i], 1);
         norm = ztmp.r;
         norm = sqrt(norm);

         if (norm < 5.0L*machEps) {
            numPrevious--;
            Num_zcopy_zprimme(length*(numPrevious-i), 
              &previousVectors[length*(i+1)], 1, &previousVectors[length*i], 1);
            zeroed = 1;
            break;
         }

            {ztmp.r = 1.0L/norm; ztmp.i = 0.0L;}
         Num_scal_zprimme(length, ztmp, &previousVectors[length*i], 1);

      } /* for 3 reorthos */

      if (!zeroed) {
         i++;
      }

   } /* main while loop */

   return numPrevious;
}
