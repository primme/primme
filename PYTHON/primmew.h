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
#include "../PRIMMESRC/SVDS/COMMONSRC/primme_svds.h"
#include <cstring>
#include <cassert>
#include <complex>

class primme_params_w : public primme_params {
   public:
   int __kind;

   primme_params_w() {
      primme_initialize(static_cast<primme_params*>(this));
      matrixMatvec = mymatvec;
      applyPreconditioner = myprevec;
      correctionParams.precondition = 0;
   }

   virtual ~primme_params_w() {}

   void display() {
      if (outputFile) outputFile = stdout;
      primme_display_params(*static_cast<primme_params*>(this));
   }

   void set_method(primme_preset_method method) {
      primme_set_method(method, static_cast<primme_params*>(this));
   }

   virtual void matvec(int len1YD, int len2YD, double *yd, int len1XD, int len2XD, double *xd)=0;
   virtual void matvec(int len1YD, int len2YD, std::complex<double> *yd, int len1XD, int len2XD, std::complex<double> *xd)=0;
   virtual void prevec(int len1YD, int len2YD, double *yd, int len1XD, int len2XD, double *xd)=0;
   virtual void prevec(int len1YD, int len2YD, std::complex<double> *yd, int len1XD, int len2XD, std::complex<double> *xd)=0;

   private:

   static void mymatvec(void *x,  void *y, int *blockSize, struct primme_params *primme) {
      primme_params_w *pp = static_cast<primme_params_w*>(primme);
      if (pp->__kind == 1)
         pp->matvec(primme->nLocal, *blockSize, (double*)x, primme->nLocal, *blockSize, (double*)y);
      else if (pp->__kind == 3)
         pp->matvec(primme->nLocal, *blockSize, (std::complex<double>*)x, primme->nLocal, *blockSize, (std::complex<double>*)y);
   }

   static void myprevec(void *x,  void *y, int *blockSize, struct primme_params *primme) {
      primme_params_w *pp = static_cast<primme_params_w*>(primme);
      if (pp->__kind == 1)
         pp->prevec(primme->nLocal, *blockSize, (double*)x, primme->nLocal, *blockSize, (double*)y);
      else if (pp->__kind == 3)
         pp->prevec(primme->nLocal, *blockSize, (std::complex<double>*)x, primme->nLocal, *blockSize, (std::complex<double>*)y);
    }
};
