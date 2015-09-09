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


#include "../PRIMMESRC/COMMONSRC/primme.h"
#include <cstring>
#include <cassert>
#include <complex>

class primme_params_w : public primme_params {
   private:
   void *cx, *cy;
   int cbs;
   bool issety;

   static void mymatvec(void *x,  void *y, int *blockSize, struct primme_params *primme) {
      primme_params_w *pp = static_cast<primme_params_w*>(primme);
      pp->cx = x; pp->cy = y; pp->cbs = *blockSize;
      pp->issety = false;
      pp->matvec();
      assert(pp->issety);
   }

   static void myprevec(void *x,  void *y, int *blockSize, struct primme_params *primme) {
      primme_params_w *pp = static_cast<primme_params_w*>(primme);
      pp->cx = x; pp->cy = y; pp->cbs = *blockSize;
      pp->issety = false;
      pp->prevec();
      assert(pp->issety);
   }

   public:
   int __kind;

   primme_params_w() {
      primme_initialize(static_cast<primme_params*>(this));
      matrixMatvec = mymatvec;
      applyPreconditioner = myprevec;
   }

   virtual ~primme_params_w() {}

   void display() {
      if (outputFile) outputFile = stdout;
      primme_display_params(*static_cast<primme_params*>(this));
   }

   void set_method(primme_preset_method method) {
      primme_set_method(method, static_cast<primme_params*>(this));
   }

   virtual void matvec() {}
   virtual void prevec() {}

   void getXd(int *len1X, int *len2X, double **x) {
      *len1X = nLocal;
      *len2X = cbs;
      *x = (double*)cx;
   }
   void getXz(int *len1X, int *len2X, std::complex<double> **x) {
      *len1X = nLocal;
      *len2X = cbs;
      *x = (std::complex<double>*)cx;
   }
   void setYd(int len1Y, int len2Y, double *y) {
      assert(len1Y == nLocal && len2Y == cbs);
      memcpy(cy, y, nLocal*cbs*sizeof(double));
      issety = true;
   }
   void setYz(int len1Y, int len2Y, std::complex<double> *y) {
      assert(len1Y == nLocal && len2Y == cbs);
      memcpy(cy, y, nLocal*cbs*sizeof(std::complex<double>));
      issety = true;
   }
};
