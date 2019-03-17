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
 *   NOTE: THIS FILE IS AUTOMATICALLY GENERATED. PLEASE DON'T MODIFY
 ******************************************************************************/


#ifndef convergence_H
#define convergence_H
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_Sprimme)
#  define check_convergence_Sprimme CONCAT(check_convergence_,SCALAR_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_Rprimme)
#  define check_convergence_Rprimme CONCAT(check_convergence_,REAL_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_SHprimme)
#  define check_convergence_SHprimme CONCAT(check_convergence_,HOST_SCALAR_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_RHprimme)
#  define check_convergence_RHprimme CONCAT(check_convergence_,HOST_REAL_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_SXprimme)
#  define check_convergence_SXprimme CONCAT(check_convergence_,XSCALAR_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_RXprimme)
#  define check_convergence_RXprimme CONCAT(check_convergence_,XREAL_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_Shprimme)
#  define check_convergence_Shprimme CONCAT(check_convergence_,CONCAT(CONCAT(STEM_C,USE_ARITH(h,k)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_Rhprimme)
#  define check_convergence_Rhprimme CONCAT(check_convergence_,CONCAT(CONCAT(STEM_C,h),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_Ssprimme)
#  define check_convergence_Ssprimme CONCAT(check_convergence_,CONCAT(CONCAT(STEM_C,USE_ARITH(s,c)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_Rsprimme)
#  define check_convergence_Rsprimme CONCAT(check_convergence_,CONCAT(CONCAT(STEM_C,s),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_Sdprimme)
#  define check_convergence_Sdprimme CONCAT(check_convergence_,CONCAT(CONCAT(STEM_C,USE_ARITH(d,z)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_Rdprimme)
#  define check_convergence_Rdprimme CONCAT(check_convergence_,CONCAT(CONCAT(STEM_C,d),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_Sqprimme)
#  define check_convergence_Sqprimme CONCAT(check_convergence_,CONCAT(CONCAT(STEM_C,USE_ARITH(q,w)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_Rqprimme)
#  define check_convergence_Rqprimme CONCAT(check_convergence_,CONCAT(CONCAT(STEM_C,q),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_SXhprimme)
#  define check_convergence_SXhprimme CONCAT(check_convergence_,CONCAT(CONCAT(,USE_ARITH(h,k)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_RXhprimme)
#  define check_convergence_RXhprimme CONCAT(check_convergence_,CONCAT(CONCAT(,h),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_SXsprimme)
#  define check_convergence_SXsprimme CONCAT(check_convergence_,CONCAT(CONCAT(,USE_ARITH(s,c)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_RXsprimme)
#  define check_convergence_RXsprimme CONCAT(check_convergence_,CONCAT(CONCAT(,s),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_SXdprimme)
#  define check_convergence_SXdprimme CONCAT(check_convergence_,CONCAT(CONCAT(,USE_ARITH(d,z)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_RXdprimme)
#  define check_convergence_RXdprimme CONCAT(check_convergence_,CONCAT(CONCAT(,d),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_SXqprimme)
#  define check_convergence_SXqprimme CONCAT(check_convergence_,CONCAT(CONCAT(,USE_ARITH(q,w)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_RXqprimme)
#  define check_convergence_RXqprimme CONCAT(check_convergence_,CONCAT(CONCAT(,q),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_SHhprimme)
#  define check_convergence_SHhprimme CONCAT(check_convergence_,CONCAT(CONCAT(,USE_ARITH(s,c)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_RHhprimme)
#  define check_convergence_RHhprimme CONCAT(check_convergence_,CONCAT(CONCAT(,s),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_SHsprimme)
#  define check_convergence_SHsprimme CONCAT(check_convergence_,CONCAT(CONCAT(,USE_ARITH(s,c)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_RHsprimme)
#  define check_convergence_RHsprimme CONCAT(check_convergence_,CONCAT(CONCAT(,s),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_SHdprimme)
#  define check_convergence_SHdprimme CONCAT(check_convergence_,CONCAT(CONCAT(,USE_ARITH(d,z)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_RHdprimme)
#  define check_convergence_RHdprimme CONCAT(check_convergence_,CONCAT(CONCAT(,d),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_SHqprimme)
#  define check_convergence_SHqprimme CONCAT(check_convergence_,CONCAT(CONCAT(,USE_ARITH(q,w)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(check_convergence_RHqprimme)
#  define check_convergence_RHqprimme CONCAT(check_convergence_,CONCAT(CONCAT(,q),primme))
#endif
int check_convergence_dprimme(dummy_type_dprimme *X, PRIMME_INT ldX, int givenX, dummy_type_dprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_dprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_dprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_dprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_dprimme *blockNorms,
      dummy_type_dprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
#if !defined(CHECK_TEMPLATE) && !defined(check_practical_convergence)
#  define check_practical_convergence CONCAT(check_practical_convergence,SCALAR_SUF)
#endif
int check_practical_convergencedprimme(dummy_type_dprimme *R, PRIMME_INT ldR, dummy_type_dprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_dprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_dprimme *blockNorms,
      double tol, dummy_type_dprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_hprimme(dummy_type_hprimme *X, PRIMME_INT ldX, int givenX, dummy_type_hprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_hprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_hprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_sprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_sprimme *blockNorms,
      dummy_type_sprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencehprimme(dummy_type_hprimme *R, PRIMME_INT ldR, dummy_type_hprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_hprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_sprimme *blockNorms,
      double tol, dummy_type_sprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_kprimme(dummy_type_kprimme *X, PRIMME_INT ldX, int givenX, dummy_type_kprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_kprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_kprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_cprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_sprimme *blockNorms,
      dummy_type_sprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencekprimme(dummy_type_kprimme *R, PRIMME_INT ldR, dummy_type_kprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_kprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_sprimme *blockNorms,
      double tol, dummy_type_cprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_sprimme(dummy_type_sprimme *X, PRIMME_INT ldX, int givenX, dummy_type_sprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_sprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_sprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_sprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_sprimme *blockNorms,
      dummy_type_sprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencesprimme(dummy_type_sprimme *R, PRIMME_INT ldR, dummy_type_sprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_sprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_sprimme *blockNorms,
      double tol, dummy_type_sprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_cprimme(dummy_type_cprimme *X, PRIMME_INT ldX, int givenX, dummy_type_cprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_cprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_cprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_cprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_sprimme *blockNorms,
      dummy_type_sprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencecprimme(dummy_type_cprimme *R, PRIMME_INT ldR, dummy_type_cprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_cprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_sprimme *blockNorms,
      double tol, dummy_type_cprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_zprimme(dummy_type_zprimme *X, PRIMME_INT ldX, int givenX, dummy_type_zprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_zprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_zprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_zprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_dprimme *blockNorms,
      dummy_type_dprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencezprimme(dummy_type_zprimme *R, PRIMME_INT ldR, dummy_type_zprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_zprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_dprimme *blockNorms,
      double tol, dummy_type_zprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_magma_hprimme(dummy_type_magma_hprimme *X, PRIMME_INT ldX, int givenX, dummy_type_magma_hprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_magma_hprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_magma_hprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_sprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_sprimme *blockNorms,
      dummy_type_sprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencemagma_hprimme(dummy_type_magma_hprimme *R, PRIMME_INT ldR, dummy_type_magma_hprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_magma_hprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_sprimme *blockNorms,
      double tol, dummy_type_sprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_magma_kprimme(dummy_type_magma_kprimme *X, PRIMME_INT ldX, int givenX, dummy_type_magma_kprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_magma_kprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_magma_kprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_cprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_sprimme *blockNorms,
      dummy_type_sprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencemagma_kprimme(dummy_type_magma_kprimme *R, PRIMME_INT ldR, dummy_type_magma_kprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_magma_kprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_sprimme *blockNorms,
      double tol, dummy_type_cprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_magma_sprimme(dummy_type_magma_sprimme *X, PRIMME_INT ldX, int givenX, dummy_type_magma_sprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_magma_sprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_magma_sprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_sprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_sprimme *blockNorms,
      dummy_type_sprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencemagma_sprimme(dummy_type_magma_sprimme *R, PRIMME_INT ldR, dummy_type_magma_sprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_magma_sprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_sprimme *blockNorms,
      double tol, dummy_type_sprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_magma_cprimme(dummy_type_magma_cprimme *X, PRIMME_INT ldX, int givenX, dummy_type_magma_cprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_magma_cprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_magma_cprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_cprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_sprimme *blockNorms,
      dummy_type_sprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencemagma_cprimme(dummy_type_magma_cprimme *R, PRIMME_INT ldR, dummy_type_magma_cprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_magma_cprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_sprimme *blockNorms,
      double tol, dummy_type_cprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_magma_dprimme(dummy_type_magma_dprimme *X, PRIMME_INT ldX, int givenX, dummy_type_magma_dprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_magma_dprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_magma_dprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_dprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_dprimme *blockNorms,
      dummy_type_dprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencemagma_dprimme(dummy_type_magma_dprimme *R, PRIMME_INT ldR, dummy_type_magma_dprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_magma_dprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_dprimme *blockNorms,
      double tol, dummy_type_dprimme *VtBV, int ldVtBV, primme_context ctx);
int check_convergence_magma_zprimme(dummy_type_magma_zprimme *X, PRIMME_INT ldX, int givenX, dummy_type_magma_zprimme *R,
      PRIMME_INT ldR, int givenR, dummy_type_magma_zprimme *evecs, int numLocked,
      PRIMME_INT ldevecs, dummy_type_magma_zprimme *Bevecs, PRIMME_INT ldBevecs, dummy_type_zprimme *VtBV,
      int ldVtBV, int left, int right, int *flags, dummy_type_dprimme *blockNorms,
      dummy_type_dprimme *hVals, int *reset, int practConvCheck, primme_context ctx);
int check_practical_convergencemagma_zprimme(dummy_type_magma_zprimme *R, PRIMME_INT ldR, dummy_type_magma_zprimme *evecs,
      int evecsSize, PRIMME_INT ldevecs, dummy_type_magma_zprimme *Bevecs, PRIMME_INT ldBevecs,
      int left, int *iev, int numToProject, int *flags, dummy_type_dprimme *blockNorms,
      double tol, dummy_type_zprimme *VtBV, int ldVtBV, primme_context ctx);
#endif
