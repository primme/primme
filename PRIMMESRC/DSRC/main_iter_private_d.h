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
 * File: main_iter_private.h
 *
 * Purpose - Definitions used exclusively by main_iter.c
 *
 ******************************************************************************/

#ifndef MAIN_ITER_PRIVATE_H
#define MAIN_ITER_PRIVATE_H

/* Failure codes returned by main_iter */

#define MAX_ITERATIONS_REACHED    -1
#define INIT_FAILURE              -2
#define ORTHO_FAILURE             -3
#define SOLVE_H_FAILURE           -4
#define SOLVE_CORRECTION_FAILURE  -5
#define RESTART_FAILURE           -6
#define LOCK_VECTORS_FAILURE      -7

static void adjust_blockSize(int *iev, int *flag, int *blockSize, 
   int maxBlockSize, int *ievMax, int basisSize, int maxBasisSize, 
   int numLocked, int numConverged, int numWantedEvs, int matrixDimension);

static int retain_previous_coefficients(double *hVecs, double *previousHVecs, 
   int basisSize, int *iev, int blockSize, primme_params *primme);

void check_reset_flags_dprimme(int *flag, int *numConverged, 
   double *hVals, double *prevRitzVals, int numPrevRitzVals,
   double tol, double aNormEstimate, primme_params *primme);


static int verify_norms(double *V, double *W, double *hVecs, double *hVals, 
   int basisSize, double *resNorms, int *flag, double tol, double aNormEstimate,
   void *rwork, int *numConverged, primme_params *primme);

/*----------------------------------------------------------------------------*
 * The following are needed for the Dynamic Method Switching
 *----------------------------------------------------------------------------*/

typedef struct {
   /* Time measurements for various components of the solver */
   double MV_PR;          /* OPp operator MV+PR time                          */
   double MV;             /* MV time only                                     */
   double PR;             /* PRecond only                                     */
   double qmr_only;       /* c_q only the QMR iteration                       */
   double qmr_plus_MV_PR; /* a   QMR plus operators                           */
   double gdk_plus_MV_PR; /* b   GD plus operators                            */
   double gdk_plus_MV;    /* b_m GD plus MV (i.e., GD without correction)     */
                          /* The following two are not currently updated:     */
   double project_locked; /* c_p projection time per locked eigenvector in QMR*/
   double reortho_locked; /* c_g (usu 2*c_p) ortho per locked vector in outer */

   /* Average convergence estimates. Updated at restart/switch/convergence    */
   double gdk_conv_rate;  /* convergence rate of all (|r|/|r0|) seen for GD+k */
   double jdq_conv_rate;  /* convergence rate of all (|r|/|r0|) seen for JDQMR*/
   double JDQMR_slowdown; /* log(gdRate)/log(jdRate) restricted in [1.1, 2.5] */
   double ratio_MV_outer; /* TotalMV/outerIts (average of last two updates)   */

   /* To maintain a convergence rate of all residual reductions seen we need: */
   int    nextReset;      /* When to reset the following averaging sums       */
   double gdk_sum_logResReductions;/* Sumof all log(residual reductions) in GD*/
   double jdq_sum_logResReductions;/* Sumof all log(residual reductions) in JD*/
   double gdk_sum_MV;     /* Total num of MV resulting in these GD reductions */
   double jdq_sum_MV;     /* Total num of MV resulting in these JD reductions */
   int nevals_by_gdk;     /* Number of evals found by GD+k since last reset   */
   int nevals_by_jdq;     /* Number of evals found by JDQMR since last reset  */

   /* Variables to remember MV/its/time/resNorm/etc since last update/switch */
   int numIt_0;           /*Remembers starting outer its/MVs since switched   */
   int numMV_0;           /*   to current method, or since an epair converged */
   double timer_0;        /*Remembers starting time since switched to a method*/
                          /*   or since an epair converged with that method   */
   double time_in_inner;  /*Accumulative time spent in inner iterations       */
                          /*   since last switch or since an epair converged  */
   double resid_0;        /*First residual norm of the convergence of a method*/
                          /*   since last switch or since an epair converged */

   /* Weighted ratio of expected times, for final method recommendation */
   double accum_jdq_gdk;  /*Expected ratio of accumulative times of JDQMR/GD+k*/
   double accum_jdq;      /* Accumulates jdq_times += ratio*(gdk+MV+PR)       */
   double accum_gdk;      /* Accumulates gdk_times += gdk+MV+PR               */

} primme_CostModel;

static void initializeModel(primme_CostModel *model, primme_params *primme);
static void switch_from_JDQMR(primme_CostModel *model, primme_params *primme);
static void switch_from_GDpk (primme_CostModel *model, primme_params *primme);
static int update_statistics(primme_CostModel *model, primme_params *primme,
   double current_time, int recentConv, int calledAtRestart, int numConverged, 
   double currentResNorm, double aNormEst);
static double ratio_JDQMR_GDpk(primme_CostModel *CostModel, int numLocked,
   double estimate_slowdown, double estimate_ratio_outer_MV);
static void update_slowdown(primme_CostModel *model);

#if 0
static void displayModel(primme_CostModel *model);
#endif

#endif
