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
 * File: inner_solve.c
 *
 * Purpose - Solves the correction equation using hermitian simplified QMR.
 *  
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "wtime.h"
#include "const.h"
#include "numerical.h"
#include "inner_solve.h"
#include "factorize.h"
#include "update_W.h"
#include "globalsum.h"
#include "auxiliary_eigs.h"

static int apply_projected_preconditioner(SCALAR *v, SCALAR *Q, PRIMME_INT ldQ,
      SCALAR *RprojectorQ, PRIMME_INT ldRprojectorQ, SCALAR *x,
      SCALAR *RprojectorX,  PRIMME_INT ldRprojectorX, int sizeRprojectorQ,
      int sizeRprojectorX, SCALAR *xKinvx, SCALAR *UDU, int *ipivot,
      SCALAR *result, SCALAR *rwork, primme_params *primme);

static int apply_skew_projector(SCALAR *Q, PRIMME_INT ldQ, SCALAR *Qhat,
      PRIMME_INT ldQhat, SCALAR *UDU, int *ipivot, int numCols, SCALAR *v,
      SCALAR *rwork, primme_params *primme);

static int apply_projected_matrix(SCALAR *v, REAL shift, SCALAR *Q, 
      PRIMME_INT ldQ, int dimQ, SCALAR *result, SCALAR *rwork,
      primme_params *primme);

static int apply_projector(SCALAR *Q, PRIMME_INT ldQ, int numCols, SCALAR *v, 
   SCALAR *rwork, primme_params *primme);

static int dist_dot(SCALAR *x, int incx,
   SCALAR *y, int incy, primme_params *primme, SCALAR *result);

static int dist_dot_real(SCALAR *x, int incx,
   SCALAR *y, int incy, primme_params *primme, REAL *result);


/*******************************************************************************
 * Function inner_solve - This subroutine solves the correction equation
 *    
 *           (I-QQ')(I-xx')(A-shift*I)(I-xx')(I-QQ')sol = -r 
 *
 *    with Q = evecs, using hermitian simplified QMR.
 *    A preconditioner may be applied to this system to accelerate convergence.
 *    The preconditioner is assumed to approximate (A-shift*I)^{-1}.  The
 *    classical JD method as described in Templates for the Solution of 
 *    Eigenvalue Problems by Bai, et. al. requires that the preconditioner is 
 *    computed orthogonally to x and evecs. This code implements all 
 *    possible variations for projectors as defined by user parameters
 *    and setup_JD_projectors(). The QMR transparently calls the resulting
 *    projected matrix and preconditioner.
 *
 *
 * Input parameters
 * ----------------
 * x           The current Ritz vector for which the correction is being solved.
 *
 * r           The residual with respect to the Ritz vector.
 *
 * evecs       The converged Ritz vectors
 *
 * evecsHat    K^{-1}*evecs where K is a hermitian preconditioner.
 *
 * UDU         The factors of the hermitian projection (evecs'*evecsHat). 
 *
 * ipivot      The pivoting for the UDU factorization
 *
 * xKinvx      The value x'*Kinv*x needed if skew-X projection
 *
 * Lprojector  Points to an array that includes all the left projector.
 *             Can be any combination of [evecs x], [evecs], [x], NULL.
 *
 * RprojectorQ Points to an array that includes the right skew projector for Q:
 *             It can be [evecsHat] or Null
 *
 * RprojectorX Points to an array that includes the right skew projector for x:
 *             It can be [Kinvx] or Null
 *
 * sizeLprojector   Number of colums of Lprojector
 *
 * sizeRprojectorQ  Number of colums of RprojectorQ
 *
 * sizeRprojectorX  Number of colums of LprojectorX
 *
 * eval        The current Ritz value 
 *
 * shift       Correction eq. shift. The closer the shift is to the target 
 *             eigenvalue, the more accurate the correction will be.
 *
 * machEps     machine precision
 *
 * rwork       Real workspace of size 
 *             4*primme->nLocal + 2*(primme->numOrthoConst+primme->numEvals)
 *
 * rworkSize   Size of the rwork array
 *
 * primme      Structure containing various solver parameters
 *
 *
 * Input/Output parameters
 * -----------------------
 * r       The residual with respect to the Ritz vector.  May be altered upon
 *         return.
 * rnorm   On input, the 2 norm of r. No need to recompute it initially.
 *         On output, the estimated 2 norm of the updated eigenvalue residual
 * touch   Parameter used in inner solve stopping criteria
 * 
 * Output parameters
 * -----------------
 * sol   The solution (correction) of the correction equation
 *
 * Return Value
 * ------------
 * Error code: 0 upon success
 *            -1 apply_projected_preconditioner failed
 *
 ******************************************************************************/

TEMPLATE_PLEASE
int inner_solve_Sprimme(SCALAR *x, SCALAR *r, REAL *rnorm, SCALAR *evecs,
      PRIMME_INT ldevecs, SCALAR *UDU, int *ipivot, SCALAR *xKinvx,
      SCALAR *Lprojector, PRIMME_INT ldLprojector, SCALAR *RprojectorQ,
      PRIMME_INT ldRprojectorQ, SCALAR *RprojectorX, PRIMME_INT ldRprojectorX,
      int sizeLprojector, int sizeRprojectorQ, int sizeRprojectorX, SCALAR *sol,
      REAL eval, REAL shift, int *touch, double machEps, SCALAR *rwork,
      size_t rworkSize, primme_params *primme) {

   int i;             /* loop variable                                       */
   int numIts;        /* Number of inner iterations                          */
   int maxIterations; /* The maximum # iterations allowed. Depends on primme */

   SCALAR *workSpace; /* Workspace needed by UDU routine */

   /* QMR parameters */

   SCALAR *g, *d, *delta, *w, *ptmp;
   REAL alpha_prev, beta, rho_prev, rho;
   REAL Theta_prev, Theta, c, sigma_prev, tau_init, tau_prev, tau; 

   /* Parameters used to dynamically update eigenpair */
   REAL Beta=0.0, Delta=0.0, Psi=0.0, Beta_prev, Delta_prev, Psi_prev, eta;
   REAL dot_sol, eval_updated, eval_prev, eres2_updated, eres_updated=0.0;
   REAL eres_prev=0.0;
   REAL Gamma_prev, Phi_prev;
   REAL Gamma=0.0, Phi=0.0;
   REAL gamma;

   REAL LTolerance, ETolerance, LTolerance_factor, ETolerance_factor;
   int isConv;
   double aNorm;

   /* -------------------------------------------*/
   /* Subdivide the workspace into needed arrays */
   /* -------------------------------------------*/

   g      = rwork;
   d      = g + primme->nLocal;
   delta  = d + primme->nLocal;
   w      = delta + primme->nLocal;
   workSpace = w + primme->nLocal; /* This needs at least 2*numOrth+NumEvals) */
   assert(rworkSize >= (size_t)primme->nLocal*4
                       + 2*(primme->numOrthoConst+primme->numEvals));

   /* -----------------------------------------*/
   /* Set up convergence criteria by Tolerance */
   /* -----------------------------------------*/

   aNorm = max(primme->stats.estimateLargestSVal, primme->aNorm);
   tau_prev = tau_init = *rnorm;       /* Assumes zero initial guess */

   /* NOTE: In any case stop when linear system residual is less than         */
   /*       max(machEps,eps)*aNorm.                                           */
   LTolerance = machEps*aNorm;
   LTolerance_factor = 1.0;
   ETolerance = 0.0;
   ETolerance_factor = 0.0;

   switch(primme->correctionParams.convTest) {
   case primme_full_LTolerance:
      /* stop when linear system residual norm is less than aNorm*eps.        */
      /* NOTE: the criterion is covered by the default values set before.     */
       break;
   case primme_decreasing_LTolerance:
      /* stop when linear system residual norm is less than relTolBase^-its   */
      LTolerance = max(LTolerance,
            pow(primme->correctionParams.relTolBase, 
               -(double)*touch));
      (*touch)++;
      break;
   case primme_adaptive:
      /* stop when estimate eigenvalue residual norm is less than aNorm*eps.  */
      /* Eigenresidual tol may not be achievable, because it iterates on      */
      /* P(A-s)P not on (A-s). But tau reflects the residual norm on P(A-s)P. */
      /* So stop when linear system residual norm or the estimate eigenvalue  */
      /* residual norm is less than aNorm*eps/1.8.                            */
      LTolerance_factor = pow(1.8, -(double)*touch);
      ETolerance_factor = pow(1.8, -(double)*touch);
      break; 
   case primme_adaptive_ETolerance:
      /* Besides the primme_adaptive criteria, stop when estimate eigenvalue  */
      /* residual norm is less than tau_init*0.1                              */
      LTolerance_factor = pow(1.8, -(double)*touch);
      ETolerance_factor = pow(1.8, -(double)*touch);
      ETolerance = tau_init*0.1;
     }
   
   /* --------------------------------------------------------*/
   /* Set up convergence criteria by max number of iterations */
   /* --------------------------------------------------------*/

   /* compute first total number of remaining matvecs */

   if (primme->maxMatvecs > 0) {
      maxIterations = primme->maxMatvecs - primme->stats.numMatvecs;
   }
   else {
      maxIterations = INT_MAX;
   }

   /* Perform primme.maxInnerIterations, but do not exceed total remaining */
   if (primme->correctionParams.maxInnerIterations > 0) {

      maxIterations = min(primme->correctionParams.maxInnerIterations, 
                          maxIterations);
   }

   /* --------------------------------------------------------*/
   /* Rest of initializations                                 */
   /* --------------------------------------------------------*/

   /* Assume zero initial guess */
   Num_copy_Sprimme(primme->nLocal, r, 1, g, 1);

   CHKERR(apply_projected_preconditioner(g, evecs, ldevecs, RprojectorQ,
           ldRprojectorQ, x, RprojectorX, ldRprojectorX, sizeRprojectorQ,
           sizeRprojectorX, xKinvx, UDU, ipivot, d, workSpace, primme), -1);

   Theta_prev = 0.0L;
   eval_prev = eval;
   CHKERR(dist_dot_real(g, 1, d, 1, primme, &rho_prev), -1);

   /* Initialize recurrences used to dynamically update the eigenpair */

   Beta_prev = Delta_prev = Psi_prev = 0.0L;
   Gamma_prev = Phi_prev = 0.0L;

   /* other initializations */
   for (i = 0; i < primme->nLocal; i++) {
      delta[i] = 0.0;
      sol[i] = 0.0;
   }

   numIts = 0;
      
   /*----------------------------------------------------------------------*/
   /*------------------------ Begin Inner Loop ----------------------------*/
   /*----------------------------------------------------------------------*/

   while (numIts < maxIterations) {

      CHKERR(apply_projected_matrix(d, shift, Lprojector, ldLprojector,
               sizeLprojector, w, workSpace, primme), -1);
      CHKERR(dist_dot_real(d, 1, w, 1, primme, &sigma_prev), -1);

      if (sigma_prev == 0.0L) {
         if (primme->printLevel >= 5 && primme->procID == 0) {
            fprintf(primme->outputFile,"Exiting because SIGMA %e\n",sigma_prev);
         }
         /* sol = r if first iteration */
         if (numIts == 0) {
            Num_copy_Sprimme(primme->nLocal, r, 1, sol, 1);
         }
         break;
      }

      alpha_prev = rho_prev/sigma_prev;
      if (fabs(alpha_prev) < machEps || fabs(alpha_prev) > 1.0L/machEps){
         if (primme->printLevel >= 5 && primme->procID == 0) {
            fprintf(primme->outputFile,"Exiting because ALPHA %e\n",alpha_prev);
         }
         /* sol = r if first iteration */
         if (numIts == 0) {
            Num_copy_Sprimme(primme->nLocal, r, 1, sol, 1);
         }
         break;
      }

      Num_axpy_Sprimme(primme->nLocal, -alpha_prev, w, 1, g, 1);

      CHKERR(dist_dot_real(g, 1, g, 1, primme, &Theta), -1);
      Theta = sqrt(Theta);
      Theta = Theta/tau_prev;
      c = 1.0L/sqrt(1+Theta*Theta);
      tau = tau_prev*Theta*c;

      gamma = c*c*Theta_prev*Theta_prev;
      eta = alpha_prev*c*c;
      for (i = 0; i < primme->nLocal; i++) {
          delta[i] = gamma*delta[i] + eta*d[i];
          sol[i] = delta[i]+sol[i];
      }
      numIts++;

      if (fabs(rho_prev) == 0.0L ) {
         if (primme->printLevel >= 5 && primme->procID == 0) {
            fprintf(primme->outputFile,"Exiting because abs(rho) %e\n",
               fabs(rho_prev));
         }
         break;
      }
      
      if (numIts > 1 && tau < LTolerance) {
         if (primme->printLevel >= 5 && primme->procID == 0) {
            fprintf(primme->outputFile, " tau < LTol %e %e\n",tau, LTolerance);
         }
         break;
      }
      if (ETolerance > 0.0 || ETolerance_factor > 0.0) {
         /* --------------------------------------------------------*/
         /* Adaptive stopping based on dynamic monitoring of eResid */
         /* --------------------------------------------------------*/

         /* Update the Ritz value and eigenresidual using the */
         /* following recurrences.                            */
      
         Delta = gamma*Delta_prev + eta*rho_prev;
         Beta = Beta_prev - Delta;
         Phi = gamma*gamma*Phi_prev + eta*eta*sigma_prev;
         Psi = gamma*Psi_prev + gamma*Phi_prev;
         Gamma = Gamma_prev + 2.0L*Psi + Phi;
        
         /* Perform the update: update the eigenvalue and the square of the  */
         /* residual norm.                                                   */
         
         CHKERR(dist_dot_real(sol, 1, sol, 1, primme, &dot_sol), -1);
         eval_updated = shift + (eval - shift + 2*Beta + Gamma)/(1 + dot_sol);
         eres2_updated = (tau*tau)/(1 + dot_sol) + 
            ((eval - shift + Beta)*(eval - shift + Beta))/(1 + dot_sol) - 
            (eval_updated - shift)*(eval_updated - shift);

         /* If numerical problems, let eres about the same as tau */
         eres_prev = eres_updated;
         if (eres2_updated < 0){
            eres_updated = sqrt( (tau*tau)/(1 + dot_sol) );
         }
         else 
            eres_updated = sqrt(eres2_updated);

         /* --------------------------------------------------------*/
         /* Stopping criteria                                       */
         /* --------------------------------------------------------*/

         if (numIts > 1 && (tau_prev <= eres_updated || eres_prev <= tau)) {
            if (primme->printLevel >= 5 && primme->procID == 0) {
               fprintf(primme->outputFile, " tau < R eres \n");
            }
            break;
         }

         if (primme->target == primme_smallest && eval_updated > eval_prev) {
            if (primme->printLevel >= 5 && primme->procID == 0) {
               fprintf(primme->outputFile, "eval_updated > eval_prev\n");
            }
            break;
         }
         else if (primme->target == primme_largest && eval_updated < eval_prev){
            if (primme->printLevel >= 5 && primme->procID == 0) {
               fprintf(primme->outputFile, "eval_updated < eval_prev\n");
            }
            break;
         }
         else if (primme->target == primme_closest_abs
               && fabs(eval-eval_updated) > tau_init+eres_updated){
            if (primme->printLevel >= 5 && primme->procID == 0) {
               fprintf(primme->outputFile, "|eval-eval_updated| > tau0+eres\n");
            }
            break;
         }
          
         if (numIts > 1 && eres_updated < ETolerance) {
            if (primme->printLevel >= 5 && primme->procID == 0) {
               fprintf(primme->outputFile, "eres < eresTol %e \n",eres_updated);
            }
            break;
         }

         /* Check if some of the next conditions is satisfied:                */
         /* a) estimate eigenvalue residual norm (eres_updated) is less       */
         /*    than eps*aNorm*Etolerance_factor                               */
         /* b) linear system residual norm is less                            */
         /*    than eps*aNorm*LTolerance_factor                               */
         /* The result is to check if eps*aNorm is less than                  */
         /* max(tau/LTolerance_factor, eres_updated/ETolerance_factor).       */

         double tol = min(tau/LTolerance_factor, eres_updated/ETolerance_factor);
         CHKERR(convTestFun_Sprimme(eval_updated, NULL, tol, &isConv, primme),
               -1);

         if (numIts > 1 && isConv) {
            if (primme->printLevel >= 5 && primme->procID == 0) {
               fprintf(primme->outputFile, " eigenvalue and residual norm "
                     "passed convergence criterion \n");
            }
            (*touch)++;
            break;
         }

         eval_prev = eval_updated;

         /* Report inner iteration */
         if (primme->monitorFun) {
            int ZERO = 0, ONE = 1;
            primme_event EVENT_INNER_ITERATION = primme_event_inner_iteration;
            int err;
            primme->stats.elapsedTime = primme_wTimer(0);
            CHKERRM((primme->monitorFun(&eval_updated, &ONE, NULL, &ZERO,
                        &ONE, &eres_updated, NULL, NULL, NULL, NULL,
                        NULL, &numIts, &tau, &EVENT_INNER_ITERATION, primme, &err),
                     err), -1, "Error returned by monitorFun: %d", err);
         }

        /* --------------------------------------------------------*/
      } /* End of if adaptive JDQMR section                        */
        /* --------------------------------------------------------*/
      else {
         /* Check if the linear system residual norm (tau) is less            */
         /* than eps*aNorm*LTolerance_factor                                  */

         CHKERR(convTestFun_Sprimme(eval, NULL, tau/LTolerance_factor, &isConv,
                  primme), -1);

         if (numIts > 1 && isConv) {
            if (primme->printLevel >= 5 && primme->procID == 0) {
               fprintf(primme->outputFile, " eigenvalue and residual norm "
                     "passed convergence criterion \n");
            }
            break;
         }

         else if (primme->monitorFun) {
            /* Report for non adaptive inner iterations */
            int ZERO = 0, ONE = 1, UNCO = UNCONVERGED;
            primme_event EVENT_INNER_ITERATION = primme_event_inner_iteration;
            int err;
            primme->stats.elapsedTime = primme_wTimer(0);
            CHKERRM((primme->monitorFun(&eval, &ONE, &UNCO, &ZERO, &ONE, rnorm,
                        NULL, NULL, NULL, NULL, NULL, &numIts, &tau,
                        &EVENT_INNER_ITERATION, primme, &err),
                     err), -1, "Error returned by monitorFun: %d", err);
         }
      }

      if (numIts < maxIterations) {

         CHKERR(apply_projected_preconditioner(g, evecs, ldevecs, RprojectorQ, 
            ldRprojectorQ, x, RprojectorX, ldRprojectorX, sizeRprojectorQ,
            sizeRprojectorX, xKinvx, UDU, ipivot, w, workSpace, primme), -1);

         CHKERR(dist_dot_real(g, 1, w, 1, primme, &rho), -1);
         beta = rho/rho_prev;
         Num_axpy_Sprimme(primme->nLocal, beta, d, 1, w, 1);
         /* Alternate between w and d buffers in successive iterations
          * This saves a memory copy. */
         ptmp = d; d = w; w = ptmp;
      
         rho_prev = rho; 
         tau_prev = tau;
         Theta_prev = Theta;

         Delta_prev = Delta;
         Beta_prev = Beta;
         Phi_prev = Phi;
         Psi_prev = Psi;
         Gamma_prev = Gamma;
      }

     /* --------------------------------------------------------*/
   } /* End of QMR main while loop                              */
     /* --------------------------------------------------------*/

   *rnorm = eres_updated;
   return 0;
}
   

/*******************************************************************************
 * Function apply_projected_preconditioner - This routine applies the
 *    projected preconditioner to a vector v by computing:
 *
 *         result = (I-Kinvx/xKinvx*x') (I - Qhat (Q'*Qhat)^{-1}Q') Kinv*v
 *
 *    First we apply the preconditioner Kinv*v, and then the two projectors 
 *    are computed one after the other.
 *    
 * Input Parameters
 * ----------------
 * v      The vector the projected preconditioner will be applied to.
 *
 * Q      The matrix evecs where evecs are the locked/converged eigenvectors
 *
 * RprojectorQ     The matrix K^{-1}Q (often called Qhat), Q, or nothing,
 *                 as determined by setup_JD_projectors.
 *
 * x               The current Ritz vector.
 *
 * RprojectorX     The matrix K^{-1}x (if needed)
 *
 * sizeRprojectorQ The number of columns in RprojectorQ
 *
 * sizeRprojectorX The number of columns in RprojectorX
 *
 * xKinvx The value x^T (Kinv*x). It is computed in the setup_JD_projectors
 *
 * UDU    The UDU decomposition of (Q'*K^{-1}*Q).  See LAPACK routine dsytrf
 *        for more details
 *
 * ipivot Permutation array indicating how the rows of the UDU decomposition
 *        have been pivoted.
 *
 * rwork  Real work array of size 2*sizeRprojectorQ=(2*orthoConst+2*numEvals)
 *
 * primme   Structure containing various solver parameters.
 *
 *
 * Output parameters
 * -----------------
 * result The result of the application.
 *
 ******************************************************************************/

static int apply_projected_preconditioner(SCALAR *v, SCALAR *Q, PRIMME_INT ldQ,
      SCALAR *RprojectorQ, PRIMME_INT ldRprojectorQ, SCALAR *x,
      SCALAR *RprojectorX,  PRIMME_INT ldRprojectorX, int sizeRprojectorQ,
      int sizeRprojectorX, SCALAR *xKinvx, SCALAR *UDU, int *ipivot,
      SCALAR *result, SCALAR *rwork, primme_params *primme) {  

   /* Place K^{-1}v in result */
   CHKERR(applyPreconditioner_Sprimme(v, primme->nLocal, primme->nLocal, result,
            primme->nLocal, 1, primme), -1);

   CHKERR(apply_skew_projector(Q, ldQ, RprojectorQ, ldRprojectorQ, UDU, ipivot,
            sizeRprojectorQ, result, rwork, primme), -1);

   CHKERR(apply_skew_projector(x, primme->nLocal, RprojectorX, ldRprojectorX,
            xKinvx, ipivot, sizeRprojectorX, result, rwork, primme), -1);

   return 0;
}

/*******************************************************************************
 * Subroutine apply_skew_projector - Apply the skew projector to a vector v:
 *
 *     v = (I-Qhat*inv(Q'Qhat)*Q') v
 *
 *   The result is placed back in v.  Q is the matrix of converged Ritz 
 *   vectors or the current Ritz vector.
 *
 * Input Parameters
 * ----------------
 * Q       The matrix of converged Ritz vectors and the current Ritz vector
 *
 * Qhat    The matrix of K^{-1}Q
 *
 * UDU     The factorization of the (Q'*Qhat) matrix
 *
 * ipivot  The pivot array for the UDU factorization
 *
 * numCols Number of columns of Q and Qhat
 *
 * rwork   Work array of size 2*numCols
 *
 * Input/Output Parameters
 * -----------------------
 * v       The vector to be skewed orthogonalized 
 * 
 ******************************************************************************/

static int apply_skew_projector(SCALAR *Q, PRIMME_INT ldQ, SCALAR *Qhat,
      PRIMME_INT ldQhat, SCALAR *UDU, int *ipivot, int numCols, SCALAR *v,
      SCALAR *rwork, primme_params *primme) {

   if (numCols > 0) {    /* there is a projector to be applied */

      SCALAR *overlaps;  /* overlaps of v with columns of Q   */
      SCALAR *workSpace; /* Used for computing local overlaps */

      overlaps = rwork;
      workSpace = overlaps + numCols;

      /* --------------------------------------------------------*/
      /* Treat the one vector case with BLAS 1 calls             */
      /* --------------------------------------------------------*/
      if (numCols == 1) {
         /* Compute workspace = Q'*v */
         CHKERR(dist_dot(Q, 1, v, 1, primme, &overlaps[0]), -1);

         /* Backsolve only if there is a skew projector */
         if (UDU != NULL) {
            CHKERRM(ABS(UDU[0]) == 0.0, -1, "Failure factorizing UDU.");
            overlaps[0] = overlaps[0]/UDU[0];
         }
         /* Compute v=v-Qhat*overlaps */
         Num_axpy_Sprimme(primme->nLocal, -overlaps[0], Qhat, 1, v, 1);
      }
      else {
         /* ------------------------------------------------------*/
         /* More than one vectors. Use BLAS 2.                    */
         /* ------------------------------------------------------*/
         /* Compute workspace = Q'*v */
         Num_gemv_Sprimme("C", primme->nLocal, numCols, 1.0, Q, ldQ, v, 1, 0.0,
               workSpace, 1);

         /* Global sum: overlaps = Q'*v */
         CHKERR(globalSum_Sprimme(workSpace, overlaps, numCols, primme), -1);

         /* --------------------------------------------*/
         /* Backsolve only if there is a skew projector */
         /* --------------------------------------------*/
         if (UDU != NULL) {
            /* Solve (Q'Qhat)^{-1}*workSpace = overlaps = Q'*v for alpha by */
            /* backsolving  with the UDU decomposition.                 */
   
            CHKERR(UDUSolve_Sprimme(UDU, ipivot, numCols, overlaps,
                     workSpace, primme), -1);

            /* Compute v=v-Qhat*workspace */
            Num_gemv_Sprimme("N", primme->nLocal, numCols, -1.0, Qhat, ldQhat,
                  workSpace, 1, 1.0, v, 1);
         }
         else  {
            /* Compute v=v-Qhat*overlaps  */
            Num_gemv_Sprimme("N", primme->nLocal, numCols, -1.0, Qhat, ldQhat,
                  overlaps, 1, 1.0, v, 1);
         } /* UDU==null */
      } /* numCols != 1 */
   } /* numCols > 0 */

   return 0;
}


/*******************************************************************************
 * Subroutine apply_projected_matrix - This subroutine applies the 
 *    projected matrix (I-Q*Q')*(A-shift*I) to a vector v by computing 
 *    (A-shift*I)v then orthogonalizing the result with Q.
 *
 * Input Parameters
 * ----------------
 * v      The vector the projected matrix will be applied to
 *
 * shift  The amount the matrix is shifted by.
 *
 * Q      The converged Ritz vectors and the current Ritz vector
 *
 * dimQ   The number of columns of Q
 * 
 * rwork  Workspace of size 2*dimQ
 *
 * primme   Structure containing various solver parameters
 *
 *
 * Output Parameters
 * -----------------
 * result The result of the application.
 *
 ******************************************************************************/

static int apply_projected_matrix(SCALAR *v, REAL shift, SCALAR *Q, 
      PRIMME_INT ldQ, int dimQ, SCALAR *result, SCALAR *rwork,
      primme_params *primme) {

   CHKERR(matrixMatvec_Sprimme(v, primme->nLocal, primme->nLocal, result,
         primme->nLocal, 0, 1, primme), -1);
   Num_axpy_Sprimme(primme->nLocal, -shift, v, 1, result, 1); 
   if (dimQ > 0)
      CHKERR(apply_projector(Q, ldQ, dimQ, result, rwork, primme), -1);

   return 0;
}
   

/*******************************************************************************
 * Subroutine apply_projector - Apply the projector (I-Q*Q') to a vector v and
 *   place the result in v.  Q is the matrix of converged Ritz vectors and the
 *   current Ritz vector.
 *
 * Input Parameters
 * ----------------
 * Q       The matrix of converged Ritz vectors and the current Ritz vector
 *
 * nLocal  The number of rows of Q and v the process has
 * 
 * numCols Number of columns of Q
 *
 * rwork   Work array of size 2*numCols
 *
 * Input/Output Parameters
 * -----------------------
 * v       The vector to be orthogonalized against Q
 * 
 ******************************************************************************/

static int apply_projector(SCALAR *Q, PRIMME_INT ldQ, int numCols, SCALAR *v, 
   SCALAR *rwork, primme_params *primme) {

   SCALAR *overlaps;  /* overlaps of v with columns of Q   */
   SCALAR *workSpace; /* Used for computing local overlaps */

   overlaps = rwork;
   workSpace = overlaps + numCols;

   Num_gemv_Sprimme("C", primme->nLocal, numCols, 1.0, Q, ldQ, v, 1, 0.0,
         workSpace, 1);
   CHKERR(globalSum_Sprimme(workSpace, overlaps, numCols, primme), -1);
   Num_gemv_Sprimme("N", primme->nLocal, numCols, -1.0, Q, ldQ, overlaps,
         1, 1.0, v, 1);

   return 0;
}


/*******************************************************************************
 * Function dist_dot - Computes dot products in parallel.
 *
 * Input Parameters
 * ----------------
 * x, y  Operands of the dot product operation
 *
 * incx  Array increment for x.  A value of 1 implies the elements are
 *       contiguous in memory.
 *
 * incy  Array increment for y.  A value of 1 implies the elements are
 *       contiguous in memory.
 *
 * primme  Structure containing various solver parameters
 *
 * result The inner product
 *
 ******************************************************************************/

static int dist_dot(SCALAR *x, int incx,
   SCALAR *y, int incy, primme_params *primme, SCALAR *result) {
                                                                                
   SCALAR temp;
                                                                                
   temp = Num_dot_Sprimme(primme->nLocal, x, incx, y, incy);
   CHKERR(globalSum_Sprimme(&temp, result, 1, primme), -1);

   return 0;
}

/*******************************************************************************
 * Function dist_dot_real - Computes dot products in parallel and return the
 *    real part.
*
 * Input Parameters
 * ----------------
 * x, y  Operands of the dot product operation
 *
 * incx  Array increment for x.  A value of 1 implies the elements are
 *       contiguous in memory.
 *
 * incy  Array increment for y.  A value of 1 implies the elements are
 *       contiguous in memory.
 *
 * primme  Structure containing various solver parameters
 *
 * Output Parameter
 * ----------------
 * result The real part of the inner product
 *
 ******************************************************************************/

static int dist_dot_real(SCALAR *x, int incx,
   SCALAR *y, int incy, primme_params *primme, REAL *result) {
                                                                                
   SCALAR temp, product;
                                                                                
   temp = Num_dot_Sprimme(primme->nLocal, x, incx, y, incy);
   CHKERR(globalSum_Sprimme(&temp, &product, 1, primme), -1);
   *result = REAL_PART(product);

   return 0;
}
