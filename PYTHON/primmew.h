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

#include "../PRIMMESRC/COMMONSRC/primme.h"
#include "../PRIMMESRC/SVDS/COMMONSRC/primme_svds.h"

class PrimmeParams : public primme_params {
   public:
   int __kind;

   PrimmeParams() {
      primme_initialize(static_cast<primme_params*>(this));
      matrixMatvec = mymatvec;
      applyPreconditioner = myprevec;
      correctionParams.precondition = 0;
   }

   virtual ~PrimmeParams() {
      primme_Free(static_cast<primme_params*>(this));
   }

   void display() {
      if (outputFile) outputFile = stdout;
      primme_display_params(*static_cast<primme_params*>(this));
   }

   void set_method(primme_preset_method method) {
      primme_set_method(method, static_cast<primme_params*>(this));
   }

   virtual void matvec(int len1YD, int len2YD, int ldYD, double *yd, int len1XD, int len2XD, int ldXD, double *xd)=0;
   virtual void matvec(int len1YD, int len2YD, int ldYD, std::complex<double> *yd, int len1XD, int len2XD, int ldXD, std::complex<double> *xd)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, double *yd, int len1XD, int len2XD, int ldXD, double *xd)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, std::complex<double> *yd, int len1XD, int len2XD, int ldXD, std::complex<double> *xd)=0;

   private:

   static void mymatvec(void *x,  void *y, int *blockSize, struct primme_params *primme) {
      PrimmeParams *pp = static_cast<PrimmeParams*>(primme);
      if (pp->__kind == 1)
         pp->matvec(primme->nLocal, *blockSize, primme->nLocal, (double*)x, primme->nLocal, *blockSize, primme->nLocal, (double*)y);
      else if (pp->__kind == 3)
         pp->matvec(primme->nLocal, *blockSize, primme->nLocal, (std::complex<double>*)x, primme->nLocal, *blockSize, primme->nLocal, (std::complex<double>*)y);
   }

   static void myprevec(void *x,  void *y, int *blockSize, struct primme_params *primme) {
      PrimmeParams *pp = static_cast<PrimmeParams*>(primme);
      if (pp->__kind == 1)
         pp->prevec(primme->nLocal, *blockSize, primme->nLocal, (double*)x, primme->nLocal, *blockSize, primme->nLocal, (double*)y);
      else if (pp->__kind == 3)
         pp->prevec(primme->nLocal, *blockSize, primme->nLocal, (std::complex<double>*)x, primme->nLocal, *blockSize, primme->nLocal, (std::complex<double>*)y);
    }
};

class PrimmeSvdsParams : public primme_svds_params {
   public:
   int __kind;

   PrimmeSvdsParams() {
      primme_svds_initialize(static_cast<primme_svds_params*>(this));
      matrixMatvec = mymatvec;
      applyPreconditioner = myprevec;
      precondition = 0;
   }

   virtual ~PrimmeSvdsParams() {
      primme_svds_Free(static_cast<primme_svds_params*>(this));
   }

   void display() {
      if (outputFile) outputFile = stdout;
      primme_svds_display_params(*static_cast<primme_svds_params*>(this));
   }

   void set_method(primme_svds_preset_method method,
         primme_preset_method methodStage1, primme_preset_method methodStage2) {
      primme_svds_set_method(method, methodStage1, methodStage2, static_cast<primme_svds_params*>(this));
   }

   virtual void matvec(int len1YD, int len2YD, int ldYD, double *yd, int len1XD, int len2XD, int ldXD, double *xd, int transpose)=0;
   virtual void matvec(int len1YD, int len2YD, int ldYD, std::complex<double> *yd, int len1XD, int len2XD, int ldXD, std::complex<double> *xd, int transpose)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, double *yd, int len1XD, int len2XD, int ldXD, double *xd, int mode)=0;
   virtual void prevec(int len1YD, int len2YD, int ldYD, std::complex<double> *yd, int len1XD, int len2XD, int ldXD, std::complex<double> *xd, int mode)=0;

   private:

   static void mymatvec(void *x, int *ldx, void *y, int *ldy, int *blockSize, int *mode, struct primme_svds_params *primme_svds) {
      PrimmeSvdsParams *pp = static_cast<PrimmeSvdsParams*>(primme_svds);
      int m, n;
      if (mode == 0) {
         m = primme_svds->mLocal;
         n = primme_svds->nLocal;
      }
      else {
         n = primme_svds->mLocal;
         m = primme_svds->nLocal;
      }
      if (pp->__kind == 1)
         pp->matvec(n, *blockSize, *ldx, (double*)x, m, *blockSize, *ldy, (double*)y, *mode);
      else if (pp->__kind == 3)
         pp->matvec(n, *blockSize, *ldx, (std::complex<double>*)x, m, *blockSize, *ldy, (std::complex<double>*)y, *mode);
   }

   static void myprevec(void *x, int *ldx, void *y, int *ldy, int *blockSize, int *mode, struct primme_svds_params *primme_svds) {
      PrimmeSvdsParams *pp = static_cast<PrimmeSvdsParams*>(primme_svds);
      int m, n;
      if (mode == 0) {
         m = primme_svds->mLocal;
         n = primme_svds->nLocal;
      }
      else {
         n = primme_svds->mLocal;
         m = primme_svds->nLocal;
      }
       if (pp->__kind == 1)
         pp->prevec(n, *blockSize, *ldx, (double*)x, m, *blockSize, *ldy, (double*)y, *mode);
      else if (pp->__kind == 3)
         pp->prevec(m, *blockSize, *ldx, (std::complex<double>*)x, n, *blockSize, *ldy, (std::complex<double>*)y, *mode);
   }
};
