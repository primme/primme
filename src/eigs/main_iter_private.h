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
 * File: main_iter_private.h
 *
 * Purpose - Definitions used exclusively by main_iter.c
 *
 ******************************************************************************/

#ifndef MAIN_ITER_PRIVATE_H
#define MAIN_ITER_PRIVATE_H

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
static int switch_from_JDQMR(primme_CostModel *model, primme_params *primme);
static int switch_from_GDpk (primme_CostModel *model, primme_params *primme);
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
