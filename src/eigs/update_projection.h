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


#ifndef update_projection_H
#define update_projection_H
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_Sprimme)
#  define update_projection_Sprimme CONCAT(update_projection_,SCALAR_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_Rprimme)
#  define update_projection_Rprimme CONCAT(update_projection_,REAL_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_SHprimme)
#  define update_projection_SHprimme CONCAT(update_projection_,HOST_SCALAR_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_RHprimme)
#  define update_projection_RHprimme CONCAT(update_projection_,HOST_REAL_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_SXprimme)
#  define update_projection_SXprimme CONCAT(update_projection_,XSCALAR_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_RXprimme)
#  define update_projection_RXprimme CONCAT(update_projection_,XREAL_SUF)
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_Shprimme)
#  define update_projection_Shprimme CONCAT(update_projection_,CONCAT(CONCAT(STEM_C,USE_ARITH(h,k)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_Rhprimme)
#  define update_projection_Rhprimme CONCAT(update_projection_,CONCAT(CONCAT(STEM_C,h),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_Ssprimme)
#  define update_projection_Ssprimme CONCAT(update_projection_,CONCAT(CONCAT(STEM_C,USE_ARITH(s,c)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_Rsprimme)
#  define update_projection_Rsprimme CONCAT(update_projection_,CONCAT(CONCAT(STEM_C,s),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_Sdprimme)
#  define update_projection_Sdprimme CONCAT(update_projection_,CONCAT(CONCAT(STEM_C,USE_ARITH(d,z)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_Rdprimme)
#  define update_projection_Rdprimme CONCAT(update_projection_,CONCAT(CONCAT(STEM_C,d),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_Sqprimme)
#  define update_projection_Sqprimme CONCAT(update_projection_,CONCAT(CONCAT(STEM_C,USE_ARITH(q,w)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_Rqprimme)
#  define update_projection_Rqprimme CONCAT(update_projection_,CONCAT(CONCAT(STEM_C,q),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_SXhprimme)
#  define update_projection_SXhprimme CONCAT(update_projection_,CONCAT(CONCAT(,USE_ARITH(h,k)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_RXhprimme)
#  define update_projection_RXhprimme CONCAT(update_projection_,CONCAT(CONCAT(,h),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_SXsprimme)
#  define update_projection_SXsprimme CONCAT(update_projection_,CONCAT(CONCAT(,USE_ARITH(s,c)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_RXsprimme)
#  define update_projection_RXsprimme CONCAT(update_projection_,CONCAT(CONCAT(,s),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_SXdprimme)
#  define update_projection_SXdprimme CONCAT(update_projection_,CONCAT(CONCAT(,USE_ARITH(d,z)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_RXdprimme)
#  define update_projection_RXdprimme CONCAT(update_projection_,CONCAT(CONCAT(,d),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_SXqprimme)
#  define update_projection_SXqprimme CONCAT(update_projection_,CONCAT(CONCAT(,USE_ARITH(q,w)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_RXqprimme)
#  define update_projection_RXqprimme CONCAT(update_projection_,CONCAT(CONCAT(,q),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_SHhprimme)
#  define update_projection_SHhprimme CONCAT(update_projection_,CONCAT(CONCAT(,USE_ARITH(s,c)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_RHhprimme)
#  define update_projection_RHhprimme CONCAT(update_projection_,CONCAT(CONCAT(,s),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_SHsprimme)
#  define update_projection_SHsprimme CONCAT(update_projection_,CONCAT(CONCAT(,USE_ARITH(s,c)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_RHsprimme)
#  define update_projection_RHsprimme CONCAT(update_projection_,CONCAT(CONCAT(,s),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_SHdprimme)
#  define update_projection_SHdprimme CONCAT(update_projection_,CONCAT(CONCAT(,USE_ARITH(d,z)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_RHdprimme)
#  define update_projection_RHdprimme CONCAT(update_projection_,CONCAT(CONCAT(,d),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_SHqprimme)
#  define update_projection_SHqprimme CONCAT(update_projection_,CONCAT(CONCAT(,USE_ARITH(q,w)),primme))
#endif
#if !defined(CHECK_TEMPLATE) && !defined(update_projection_RHqprimme)
#  define update_projection_RHqprimme CONCAT(update_projection_,CONCAT(CONCAT(,q),primme))
#endif
int update_projection_dprimme(dummy_type_dprimme *X, PRIMME_INT ldX, dummy_type_dprimme *Y,
      PRIMME_INT ldY, dummy_type_dprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_hprimme(dummy_type_hprimme *X, PRIMME_INT ldX, dummy_type_hprimme *Y,
      PRIMME_INT ldY, dummy_type_sprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_kprimme(dummy_type_kprimme *X, PRIMME_INT ldX, dummy_type_kprimme *Y,
      PRIMME_INT ldY, dummy_type_cprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_sprimme(dummy_type_sprimme *X, PRIMME_INT ldX, dummy_type_sprimme *Y,
      PRIMME_INT ldY, dummy_type_sprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_cprimme(dummy_type_cprimme *X, PRIMME_INT ldX, dummy_type_cprimme *Y,
      PRIMME_INT ldY, dummy_type_cprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_zprimme(dummy_type_zprimme *X, PRIMME_INT ldX, dummy_type_zprimme *Y,
      PRIMME_INT ldY, dummy_type_zprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_magma_hprimme(dummy_type_magma_hprimme *X, PRIMME_INT ldX, dummy_type_magma_hprimme *Y,
      PRIMME_INT ldY, dummy_type_sprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_magma_kprimme(dummy_type_magma_kprimme *X, PRIMME_INT ldX, dummy_type_magma_kprimme *Y,
      PRIMME_INT ldY, dummy_type_cprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_magma_sprimme(dummy_type_magma_sprimme *X, PRIMME_INT ldX, dummy_type_magma_sprimme *Y,
      PRIMME_INT ldY, dummy_type_sprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_magma_cprimme(dummy_type_magma_cprimme *X, PRIMME_INT ldX, dummy_type_magma_cprimme *Y,
      PRIMME_INT ldY, dummy_type_cprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_magma_dprimme(dummy_type_magma_dprimme *X, PRIMME_INT ldX, dummy_type_magma_dprimme *Y,
      PRIMME_INT ldY, dummy_type_dprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
int update_projection_magma_zprimme(dummy_type_magma_zprimme *X, PRIMME_INT ldX, dummy_type_magma_zprimme *Y,
      PRIMME_INT ldY, dummy_type_zprimme *Z, PRIMME_INT ldZ, PRIMME_INT nLocal,
      int numCols, int blockSize, int isSymmetric, primme_context ctx);
#endif
