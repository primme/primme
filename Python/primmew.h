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
 * File: primmew.h
 * 
 * Purpose - Wrap parts of PRIMME interface to make easier the SWIG description.
 * 
 ******************************************************************************/


#include <cstring>
#include <cassert>
#include <complex>

#include "../include/primme.h"

class PrimmeParams : public primme_params {
   public:

   PrimmeParams() {
      primme_initialize(static_cast<primme_params*>(this));
      correctionParams.precondition = 0;
      globalSum_set = 0;
      monitor_set = 0;
   }

   virtual ~PrimmeParams() {
      if (targetShifts) delete [] targetShifts;
      primme_free(static_cast<primme_params*>(this));
   }

   void display() {
      if (outputFile) outputFile = stdout;
      primme_display_params(*static_cast<primme_params*>(this));
      fflush(outputFile);
   }

   void set_method(primme_preset_method method) {
      primme_set_method(method, static_cast<primme_params*>(this));
   }

   void _set_targetShifts(double *targetShifts, int n) {
      if (this->targetShifts)
         delete [] this->targetShifts;
      this->targetShifts = new double[n];
      for (int i=0; i<n; i++)
         this->targetShifts[i] = targetShifts[i];
      this->numTargetShifts = n;
   }

   void _get_targetShifts(double **targetShifts, int *n) {
      *targetShifts = this->targetShifts;
      *n = this->numTargetShifts;
   }

   virtual void matvec(int len1YD, int len2YD, int ldYD, float *yd, int len1XD, int len2XD, int ldXD, float *xd)=0;
   virtual void matvec(int len1YD, int len2YD, int ldYD, std::complex<float> *yd, int len1XD, int len2XD, int ldXD, std::complex<float> *xd)=0;
   virtual void matvec(int len1YD, int len2YD, int ldYD, double *yd, int len1XD, int len2XD, int ldXD, double *xd)=0;
   virtual void matvec(int len1YD, int len2YD, int ldYD, std::complex<double> *yd, int len1XD, int len2XD, int ldXD, std::complex<double> *xd)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, float *yd, int len1XD, int len2XD, int ldXD, float *xd)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, std::complex<float> *yd, int len1XD, int len2XD, int ldXD, std::complex<float> *xd)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, double *yd, int len1XD, int len2XD, int ldXD, double *xd)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, std::complex<double> *yd, int len1XD, int len2XD, int ldXD, std::complex<double> *xd)=0;
   virtual void globalSum(int lenYD, float *yd, int lenXD, float *xd)=0;
   virtual void globalSum(int lenYD, double *yd, int lenXD, double *xd)=0;
   int globalSum_set;
   virtual void mon(int lenbasisEvals, float *basisEvals, int lenbasisFlags, int *basisFlags,
      int leniblock, int *iblock, int lenbasisNorms, float *basisNorms, int numConverged,
      int lenlockedEvals, float *lockedEvals, int lenlockedFlags, int *lockedFlags, int lenlockedNorms, float *lockedNorms,
      int inner_its, float LSRes, int event)=0;
   virtual void mon(int lenbasisEvals, double *basisEvals, int lenbasisFlags, int *basisFlags,
      int leniblock, int *iblock, int lenbasisNorms, double *basisNorms, int numConverged,
      int lenlockedEvals, double *lockedEvals, int lenlockedFlags, int *lockedFlags, int lenlockedNorms, double *lockedNorms,
      int inner_its, double LSRes, int event)=0;
   int monitor_set;
};

class PrimmeSvdsParams : public primme_svds_params {
   public:

   PrimmeSvdsParams() {
      primme_svds_initialize(static_cast<primme_svds_params*>(this));
      precondition = 0;
      globalSum_set = 0;
      monitor_set = 0;
   }

   virtual ~PrimmeSvdsParams() {
      if (targetShifts) delete [] targetShifts;
      primme_svds_free(static_cast<primme_svds_params*>(this));
   }

   void display() {
      if (outputFile) outputFile = stdout;
      primme_svds_display_params(*static_cast<primme_svds_params*>(this));
      fflush(outputFile);
   }

   void set_method(primme_svds_preset_method method,
         primme_preset_method methodStage1, primme_preset_method methodStage2) {
      primme_svds_set_method(method, methodStage1, methodStage2, static_cast<primme_svds_params*>(this));
   }

   void _set_targetShifts(double *targetShifts, int n) {
      if (this->targetShifts)
         delete [] this->targetShifts;
      this->targetShifts = new double[n];
      for (int i=0; i<n; i++)
         this->targetShifts[i] = targetShifts[i];
      this->numTargetShifts = n;
   }

   void _get_targetShifts(double **targetShifts, int *n) {
      *targetShifts = this->targetShifts;
      *n = this->numTargetShifts;
   }

   virtual void matvec(int len1YD, int len2YD, int ldYD, float *yd, int len1XD, int len2XD, int ldXD, float *xd, int transpose)=0;
   virtual void matvec(int len1YD, int len2YD, int ldYD, std::complex<float> *yd, int len1XD, int len2XD, int ldXD, std::complex<float> *xd, int transpose)=0;
   virtual void matvec(int len1YD, int len2YD, int ldYD, double *yd, int len1XD, int len2XD, int ldXD, double *xd, int transpose)=0;
   virtual void matvec(int len1YD, int len2YD, int ldYD, std::complex<double> *yd, int len1XD, int len2XD, int ldXD, std::complex<double> *xd, int transpose)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, float *yd, int len1XD, int len2XD, int ldXD, float *xd, int mode)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, std::complex<float> *yd, int len1XD, int len2XD, int ldXD, std::complex<float> *xd, int mode)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, double *yd, int len1XD, int len2XD, int ldXD, double *xd, int mode)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, std::complex<double> *yd, int len1XD, int len2XD, int ldXD, std::complex<double> *xd, int mode)=0;
   virtual void globalSum(int lenYD, float *yd, int lenXD, float *xd)=0;
   virtual void globalSum(int lenYD, double *yd, int lenXD, double *xd)=0;
   int globalSum_set;
   virtual void mon(int lenbasisSvals, float *basisSvals, int lenbasisFlags, int *basisFlags,
      int leniblock, int *iblock, int lenbasisNorms, float *basisNorms, int numConverged,
      int lenlockedSvals, float *lockedSvals, int lenlockedFlags, int *lockedFlags, int lenlockedNorms, float *lockedNorms,
      int inner_its, float LSRes, int event, int stage)=0;
   virtual void mon(int lenbasisSvals, double *basisSvals, int lenbasisFlags, int *basisFlags,
      int leniblock, int *iblock, int lenbasisNorms, double *basisNorms, int numConverged,
      int lenlockedSvals, double *lockedSvals, int lenlockedFlags, int *lockedFlags, int lenlockedNorms, double *lockedNorms,
      int inner_its, double LSRes, int event, int stage)=0;
   int monitor_set;
};
