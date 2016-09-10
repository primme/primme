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
   }

   virtual ~PrimmeParams() {
      if (targetShifts) delete [] targetShifts;
      primme_Free(static_cast<primme_params*>(this));
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
};

class PrimmeSvdsParams : public primme_svds_params {
   public:

   PrimmeSvdsParams() {
      primme_svds_initialize(static_cast<primme_svds_params*>(this));
      precondition = 0;
      globalSum_set = 0;
   }

   virtual ~PrimmeSvdsParams() {
      if (targetShifts) delete [] targetShifts;
      primme_svds_Free(static_cast<primme_svds_params*>(this));
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
};
