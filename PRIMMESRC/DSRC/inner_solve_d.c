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
 * File: inner_solve.c
 *
 * Purpose - Solves the correction equation using hermitian simplified QMR.
 *  
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "primme.h"
#include "wtime.h"
#include "inner_solve_d.h"
#include "inner_solve_private_d.h"
#include "factorize_d.h"
#include "numerical_d.h"


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
 * eresTol     The convergence tolerance for the eigenpair residual
 *
 * aNormEstimate Some approximate norm of A.
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

int inner_solve_dprimme(double *x, double *r, double *rnorm, 
   double *evecs, double *evecsHat, double *UDU, int *ipivot, 
   double *xKinvx, double *Lprojector, double *RprojectorQ, 
   double *RprojectorX, int sizeLprojector, int sizeRprojectorQ, 
   int sizeRprojectorX, double *sol, double eval, double shift, 
   double eresTol, double aNormEstimate, double machEps, double *rwork, 
   int rworkSize, primme_params *primme) {

   int i;             /* loop variable                                       */
   int numIts;        /* Number of inner iterations                          */
   int ret;           /* Return value used for error checking.               */
   int maxIterations; /* The maximum # iterations allowed. Depends on primme */

   double *workSpace; /* Workspace needed by UDU routine */

   /* QMR parameters */

   double *g, *d, *delta, *w, *ptmp;
   double alpha_prev, beta, rho_prev, rho;
   double Theta_prev, Theta, c, sigma_prev, tau_init, tau_prev, tau; 

   /* Parameters used to dynamically update eigenpair */
   double Beta, Delta, Psi, Beta_prev, Delta_prev, Psi_prev, eta;
   double dot_sol, eval_updated, eval_prev, eres2_updated, eres_updated, R;
   double Gamma_prev, Phi_prev;
   double Gamma, Phi;
   double gamma;

   /* The convergence criteria of the inner linear system must satisfy:       */
   /* || current residual || <= relativeTolerance * || initial residual ||    */
   /*                                               + absoluteTol             */

   double relativeTolerance; 
   double absoluteTolerance;
   double LTolerance, ETolerance;

   /* Some constants                                                          */
   double tzero = +0.0e+00;

   /* -------------------------------------------*/
   /* Subdivide the workspace into needed arrays */
   /* -------------------------------------------*/

   g      = rwork;
   d      = g + primme->nLocal;
   delta  = d + primme->nLocal;
   w      = delta + primme->nLocal;
   workSpace = w + primme->nLocal; /* This needs at least 2*numOrth+NumEvals) */
   
   /* -----------------------------------------*/
   /* Set up convergence criteria by Tolerance */
   /* -----------------------------------------*/

   if (primme->aNorm <= 0.0L) {
      absoluteTolerance = aNormEstimate*machEps;
      eresTol = eresTol*aNormEstimate;
   }
   else {
      absoluteTolerance = primme->aNorm*machEps;
   }
   tau_prev = tau_init = *rnorm;       /* Assumes zero initial guess */
   LTolerance = eresTol;

   /* Andreas: note that eigenresidual tol may not be achievable, because we */
   /* iterate on P(A-s)P not (A-s). But tau reflects linSys on P(A-s)P. */
   if (primme->correctionParams.convTest == primme_adaptive) {
      ETolerance = max(eresTol/1.8L, absoluteTolerance);
      LTolerance = ETolerance;
   }
   else if (primme->correctionParams.convTest == primme_adaptive_ETolerance) {
      LTolerance = max(eresTol/1.8L, absoluteTolerance);
      ETolerance = max(tau_init*0.1L, LTolerance);
   }
   else if (primme->correctionParams.convTest == primme_decreasing_LTolerance) {
      relativeTolerance = pow(primme->correctionParams.relTolBase, 
         (double)-primme->stats.numOuterIterations);
      LTolerance = relativeTolerance * tau_init 
                   + absoluteTolerance + eresTol;
   /*printf(" RL %e INI %e abso %e LToler %e aNormEstimate %e \n", */
   /*relativeTolerance, tau_init, absoluteTolerance,LTolerance,aNormEstimate);*/
   }
   
   /* --------------------------------------------------------*/
   /* Set up convergence criteria by max number of iterations */
   /* --------------------------------------------------------*/

   /* compute first total number of remaining matvecs */

   maxIterations = primme->maxMatvecs - primme->stats.numMatvecs;

   /* Perform primme.maxInnerIterations, but do not exceed total remaining */
   if (primme->correctionParams.maxInnerIterations > 0) {

      maxIterations = min(primme->correctionParams.maxInnerIterations, 
                          maxIterations);
   }

   /* --------------------------------------------------------*/
   /* Rest of initializations                                 */
   /* --------------------------------------------------------*/

   /* Assume zero initial guess */
   Num_dcopy_dprimme(primme->nLocal, r, 1, g, 1);

   ret = apply_projected_preconditioner(g, evecs, RprojectorQ, 
           x, RprojectorX, sizeRprojectorQ, sizeRprojectorX, 
           xKinvx, UDU, ipivot, d, workSpace, primme);

   if (ret != 0) {
      primme_PushErrorMessage(Primme_inner_solve, 
         Primme_apply_projected_preconditioner, ret, __FILE__, __LINE__, 
         primme);
      return APPLYPROJECTEDPRECONDITIONER_FAILURE;
   }
      
   Theta_prev = 0.0L;
   eval_prev = eval;
   rho_prev = dist_dot(g, 1, d, 1, primme);
      
   /* Initialize recurrences used to dynamically update the eigenpair */

   Beta_prev = Delta_prev = Psi_prev = 0.0L;
   Gamma_prev = Phi_prev = 0.0L;

   /* other initializations */
   for (i = 0; i < primme->nLocal; i++) {
      delta[i] = tzero;
      sol[i] = tzero;
   }

   numIts = 0;
      
   /*----------------------------------------------------------------------*/
   /*------------------------ Begin Inner Loop ----------------------------*/
   /*----------------------------------------------------------------------*/

   while (numIts < maxIterations) {

      apply_projected_matrix(d, shift, Lprojector, sizeLprojector, 
                             w, workSpace, primme);
      sigma_prev = dist_dot(d, 1, w, 1, primme);

      if (sigma_prev == 0.0L) {
         if (primme->printLevel >= 5 && primme->procID == 0) {
            fprintf(primme->outputFile,"Exiting because SIGMA %e\n",sigma_prev);
         }
         break;
      }

      alpha_prev = rho_prev/sigma_prev;
      if (fabs(alpha_prev) < machEps || fabs(alpha_prev) > 1.0L/machEps){
         if (primme->printLevel >= 5 && primme->procID == 0) {
            fprintf(primme->outputFile,"Exiting because ALPHA %e\n",alpha_prev);
         }
         break;
      }

      Num_axpy_dprimme(primme->nLocal, -alpha_prev, w, 1, g, 1);

      Theta = dist_dot(g, 1, g, 1, primme);
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
      
      if (tau < LTolerance) {
         if (primme->printLevel >= 5 && primme->procID == 0) {
            fprintf(primme->outputFile, " tau < LTol %e %e\n",tau, LTolerance);
         }
         break;
      }
      else if (primme->correctionParams.convTest == primme_adaptive_ETolerance
            || primme->correctionParams.convTest == primme_adaptive) {
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
         
         dot_sol = dist_dot(sol, 1, sol, 1, primme);
         eval_updated = shift + (eval - shift + 2*Beta + Gamma)/(1 + dot_sol);
         eres2_updated = (tau*tau)/(1 + dot_sol) + 
            ((eval - shift + Beta)*(eval - shift + Beta))/(1 + dot_sol) - 
            (eval_updated - shift)*(eval_updated - shift);

         /* If numerical problems, let eres about the same as tau */
         if (eres2_updated < 0){
            eres_updated = sqrt( (tau*tau)/(1 + dot_sol) );
         }
         else 
            eres_updated = sqrt(eres2_updated);

         /* --------------------------------------------------------*/
         /* Stopping criteria                                       */
         /* --------------------------------------------------------*/

         R = max(0.9878, sqrt(tau/tau_prev))*sqrt(1+dot_sol);
        
         if ( tau <= R*eres_updated || eres_updated <= tau*R ) {
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
         
         if (eres_updated < ETolerance) {    /* tau < LTol has been checked */
            if (primme->printLevel >= 5 && primme->procID == 0) {
               fprintf(primme->outputFile, "eres < eresTol %e \n",eres_updated);
            }
            break;
         }

         eval_prev = eval_updated;

         if (primme->printLevel >= 4 && primme->procID == 0) {
            fprintf(primme->outputFile,
           "INN MV %d Sec %e Eval %e Lin|r| %.3e EV|r| %.3e\n", primme->stats.
            numMatvecs, primme_wTimer(0), eval_updated, tau, eres_updated);
            fflush(primme->outputFile);
         }

        /* --------------------------------------------------------*/
      } /* End of if adaptive JDQMR section                        */
        /* --------------------------------------------------------*/
      else if (primme->printLevel >= 4 && primme->procID == 0) {
        /* Report for non adaptive inner iterations */
        fprintf(primme->outputFile,
           "INN MV %d Sec %e Lin|r| %e\n", primme->stats.numMatvecs,
           primme_wTimer(0),tau);
        fflush(primme->outputFile);
      }

      if (numIts < maxIterations) {

         ret = apply_projected_preconditioner(g, evecs, RprojectorQ, 
            x, RprojectorX, sizeRprojectorQ, sizeRprojectorX, 
            xKinvx, UDU, ipivot, w, workSpace, primme);

         if (ret != 0) {
            primme_PushErrorMessage(Primme_inner_solve, 
               Primme_apply_projected_preconditioner, ret, __FILE__, __LINE__, 
               primme);
               ret = APPLYPROJECTEDPRECONDITIONER_FAILURE;
               break;
         }
         rho = dist_dot(g, 1, w, 1, primme);
         beta = rho/rho_prev;
         Num_axpy_dprimme(primme->nLocal, beta, d, 1, w, 1);
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

static int apply_projected_preconditioner(double *v, double *Q, 
   double *RprojectorQ, double *x, double *RprojectorX, 
   int sizeRprojectorQ, int sizeRprojectorX, double *xKinvx, 
   double *UDU, int *ipivot, double *result, double *rwork, 
   primme_params *primme) {  

   int ONE = 1;
   int ret;

   if (primme->correctionParams.precondition) {
      /* Place K^{-1}v in result */
      (*primme->applyPreconditioner)(v, result, &ONE, primme);
      primme->stats.numPreconds += 1;
   }
   else {
      Num_dcopy_dprimme(primme->nLocal, v, 1, result, 1);
   }

   ret = apply_skew_projector(Q, RprojectorQ, UDU, ipivot, sizeRprojectorQ,
                           result, rwork, primme);
   if (ret != 0) {
         primme_PushErrorMessage(Primme_apply_projected_preconditioner, 
            Primme_apply_skew_projector, ret, __FILE__, __LINE__, primme);
         return APPLYSKEWPROJECTOR_FAILURE;
   }

   apply_skew_projector(x, RprojectorX, xKinvx, ipivot, sizeRprojectorX,
                           result, rwork, primme);
   if (ret != 0) {
         primme_PushErrorMessage(Primme_apply_projected_preconditioner, 
            Primme_apply_skew_projector, ret, __FILE__, __LINE__, primme);
         return APPLYSKEWPROJECTOR_FAILURE;
   }

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

static int apply_skew_projector(double *Q, double *Qhat, double *UDU, 
   int *ipivot, int numCols, double *v, double *rwork, 
   primme_params *primme) {

   int count;
   double tpone = +1.0e+00, tzero = +0.0e+00, tmone = -1.0e+00;

   if (numCols > 0) {    /* there is a projector to be applied */

      int ret;
      double *overlaps;  /* overlaps of v with columns of Q   */
      double *workSpace; /* Used for computing local overlaps */

      overlaps = rwork;
      workSpace = overlaps + numCols;

      /* --------------------------------------------------------*/
      /* Treat the one vector case with BLAS 1 calls             */
      /* --------------------------------------------------------*/
      if (numCols == 1) {
         /* Compute workspace = Q'*v */
         overlaps[0] = dist_dot(Q, 1, v, 1, primme);

         /* Backsolve only if there is a skew projector */
         if (UDU != NULL) {
            if (UDU[0] == 0.0L) {
               return UDUSOLVE_FAILURE;
            }
            overlaps[0] = overlaps[0]/UDU[0];
         }
         /* Compute v=v-Qhat*overlaps */
         Num_axpy_dprimme(primme->nLocal, -overlaps[0], Qhat, 1, v, 1);
      }
      else {
         /* ------------------------------------------------------*/
         /* More than one vectors. Use BLAS 2.                    */
         /* ------------------------------------------------------*/
         /* Compute workspace = Q'*v */
         Num_gemv_dprimme("C", primme->nLocal, numCols, tpone, Q, 
                      primme->nLocal, v, 1, tzero, workSpace, 1);

         /* Global sum: overlaps = Q'*v */
         count = numCols;
         (*primme->globalSumDouble)(workSpace, overlaps, &count, primme);   

         /* --------------------------------------------*/
         /* Backsolve only if there is a skew projector */
         /* --------------------------------------------*/
         if (UDU != NULL) {
            /* Solve (Q'Qhat)^{-1}*workSpace = overlaps = Q'*v for alpha by */
            /* backsolving  with the UDU decomposition.                 */
   
            ret = UDUSolve_dprimme(UDU, ipivot, numCols, overlaps, workSpace);
            if (ret != 0) {
               primme_PushErrorMessage(Primme_apply_skew_projector,
                  Primme_udusolve, ret, __FILE__, __LINE__, primme);
               return UDUSOLVE_FAILURE;
            }
            /* Compute v=v-Qhat*workspace */
            Num_gemv_dprimme("N", primme->nLocal, numCols, tmone, Qhat, 
                         primme->nLocal, workSpace, 1, tpone, v, 1);
         }
         else  {
            /* Compute v=v-Qhat*overlaps  */
            Num_gemv_dprimme("N", primme->nLocal, numCols, tmone, Qhat, 
                         primme->nLocal, overlaps, 1, tpone, v, 1);
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

static void apply_projected_matrix(double *v, double shift, double *Q, 
   int dimQ, double *result, double *rwork, primme_params *primme) {
   
   int ONE = 1;   /* For passing it by reference in matrixMatvec */

   (*primme->matrixMatvec)(v, result, &ONE, primme);
   Num_axpy_dprimme(primme->nLocal, -shift, v, 1, result, 1); 
   if (dimQ > 0)
      apply_projector(Q, dimQ, result, rwork, primme); 

   primme->stats.numMatvecs += 1;
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

static void apply_projector(double *Q, int numCols, double *v, 
   double *rwork, primme_params *primme) {

   int count;          /* globalSum counter                 */
   double *overlaps;  /* overlaps of v with columns of Q   */
   double *workSpace; /* Used for computing local overlaps */

   double tpone = +1.0e+00, tzero = +0.0e+00, tmone = -1.0e+00;

   overlaps = rwork;
   workSpace = overlaps + numCols;

   Num_gemv_dprimme("C", primme->nLocal, numCols, tpone, Q, primme->nLocal,
      v, 1, tzero, workSpace, 1);
   count = numCols;
   (*primme->globalSumDouble)(workSpace, overlaps, &count, primme);   
   Num_gemv_dprimme("N", primme->nLocal, numCols, tmone, Q, primme->nLocal,
      overlaps, 1, tpone, v, 1);

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
 ******************************************************************************/

static double dist_dot(double *x, int incx,
   double *y, int incy, primme_params *primme) {
                                                                                
   double temp, product;
   int count;
                                                                                
   temp = Num_dot_dprimme(primme->nLocal, x, incx, y, incy);
   count = 1;
   (*primme->globalSumDouble)(&temp, &product, &count, primme);
   return product;
                                                                                
}
