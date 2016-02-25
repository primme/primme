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
 * File: ioandtest.c
 *
 * Purpose - Functions to check the computed eigenpairs and to read and write
 *           eigenvectors and primme_params.
 * 
 ******************************************************************************/

#include <math.h>
#include <assert.h>
#include "primme.h"
#include "num.h"
#include "ioandtest.h"

static PRIMME_NUM primme_dot(PRIMME_NUM *x, PRIMME_NUM *y, primme_params *primme) {
   PRIMME_NUM aux, aux0;
#ifdef USE_DOUBLECOMPLEX
   int one = 2;
#else
   int one = 1;
#endif
   aux = COMPLEXV(SUF(Num_dot)(primme->nLocal, COMPLEXZ(x), 1, COMPLEXZ(y), 1));
   if (primme->globalSumDouble) {
      primme->globalSumDouble(&aux, &aux0, &one, primme);
      return aux0;
   }
   return aux;
}

static PRIMME_NUM primme_svds_dot(PRIMME_NUM *x, PRIMME_NUM *y, int trans, primme_svds_params *primme) {
   PRIMME_NUM aux, aux0;
#ifdef USE_DOUBLECOMPLEX
   int one = 2;
#else
   int one = 1;
#endif
   aux = COMPLEXV(SUF(Num_dot)(trans ? primme->nLocal : primme->mLocal, COMPLEXZ(x), 1, COMPLEXZ(y), 1));
   if (primme->globalSumDouble) {
      primme->globalSumDouble(&aux, &aux0, &one, primme);
      return aux0;
   }
   return aux;
}

#undef __FUNCT__
#define __FUNCT__ "check_solution"
int check_solution(const char *checkXFileName, primme_params *primme, double *evals,
                   PRIMME_NUM *evecs, double *rnorms, int *perm) {

   double eval0, rnorm0, prod, bound, delta;
   PRIMME_NUM *Ax, *r, *X=NULL, *h, *h0;
   int i, j, cols, retX=0, one=1;
   primme_params primme0;

   /* Read stored eigenvectors and primme_params */
   ASSERT_MSG(readBinaryEvecsAndPrimmeParams(checkXFileName, NULL, &X, primme->n, primme->n, &cols,
                                             primme->nLocal, perm, &primme0) == 0, -1, "");
   /* Check primme_params */
#  define CHECK_PRIMME_PARAM(F) \
        if (primme0. F != primme-> F ) { \
           fprintf(stderr, "Warning: discrepancy in primme." #F ", %d should be close to %d\n", primme-> F , primme0. F ); \
           retX = 1; \
        }
#  define CHECK_PRIMME_PARAM_DOUBLE(F) \
        if (fabs(primme0. F - primme-> F) > primme-> F * 1e-14) { \
           fprintf(stderr, "Warning: discrepancy in primme." #F ", %.16e should be close to %.16e\n", primme-> F , primme0. F ); \
           retX = 1; \
        }
#  define CHECK_PRIMME_PARAM_TOL(F, T) \
        if (abs(primme0. F - primme-> F ) > primme-> F * T /100+1) { \
           fprintf(stderr, "Warning: discrepancy in primme." #F ", %d should be close to %d\n", primme-> F , primme0. F ); \
           retX = 1; \
        }

   if (primme0.n) {
      CHECK_PRIMME_PARAM(n);
      CHECK_PRIMME_PARAM(numEvals);
      CHECK_PRIMME_PARAM(target);
      CHECK_PRIMME_PARAM(numTargetShifts);
      CHECK_PRIMME_PARAM(dynamicMethodSwitch);
      CHECK_PRIMME_PARAM(locking);
      CHECK_PRIMME_PARAM(numOrthoConst);
      CHECK_PRIMME_PARAM(maxBasisSize);
      CHECK_PRIMME_PARAM(minRestartSize);
      CHECK_PRIMME_PARAM(restartingParams.scheme);
      CHECK_PRIMME_PARAM(restartingParams.maxPrevRetain);
      CHECK_PRIMME_PARAM(correctionParams.precondition);
      CHECK_PRIMME_PARAM(correctionParams.robustShifts);
      CHECK_PRIMME_PARAM(correctionParams.maxInnerIterations);
      CHECK_PRIMME_PARAM(correctionParams.projectors.LeftQ);
      CHECK_PRIMME_PARAM(correctionParams.projectors.LeftX);
      CHECK_PRIMME_PARAM(correctionParams.projectors.RightQ);
      CHECK_PRIMME_PARAM(correctionParams.projectors.RightX);
      CHECK_PRIMME_PARAM(correctionParams.projectors.SkewQ);
      CHECK_PRIMME_PARAM(correctionParams.projectors.SkewX);
      CHECK_PRIMME_PARAM(correctionParams.convTest);
      CHECK_PRIMME_PARAM_DOUBLE(aNorm);
      CHECK_PRIMME_PARAM_DOUBLE(eps);
      CHECK_PRIMME_PARAM_DOUBLE(correctionParams.relTolBase);
      CHECK_PRIMME_PARAM(initSize);
      CHECK_PRIMME_PARAM_TOL(stats.numMatvecs, 40);
   }

#  undef CHECK_PRIMME_PARAM
#  undef CHECK_PRIMME_PARAM_DOUBLE
#  undef CHECK_PRIMME_PARAM_TOL

   h = (PRIMME_NUM *)primme_calloc(max(cols*2, primme->initSize), sizeof(PRIMME_NUM), "h"); h0 = &h[cols];
   Ax = (PRIMME_NUM *)primme_calloc(primme->nLocal, sizeof(PRIMME_NUM), "Ax");
   r = (PRIMME_NUM *)primme_calloc(primme->nLocal, sizeof(PRIMME_NUM), "r");

   /* Estimate the separation between eigenvalues */
   delta = primme->aNorm;
   for (i=1; i < primme->initSize; i++) {
      delta = min(delta, fabs(evals[i]-evals[i-1]));
   }

   for (i=0; i < primme->initSize; i++) {
      /* Check |V(:,0:i-1)'V(:,i)| < sqrt(machEps) */
      SUF(Num_gemv)("C", primme->nLocal, i+1, COMPLEXZV(1.0), COMPLEXZ(evecs), primme->nLocal, COMPLEXZ(&evecs[primme->nLocal*i]), 1, COMPLEXZV(0.), COMPLEXZ(h), 1);
      if (primme->globalSumDouble) {
         int cols0 = (i+1)*sizeof(PRIMME_NUM)/sizeof(double);
         primme->globalSumDouble(h, h0, &cols0, primme);
      }
      else h0 = h;
      prod = REAL_PARTZ(SUF(Num_dot)(i, COMPLEXZ(h0), 1, COMPLEXZ(h0), 1));
      prod = sqrt(prod);
      if (prod > 1e-7 && primme->procID == 0) {
         fprintf(stderr, "Warning: |EVecs[0:%d-1]'Evec[%d]| = %-3E\n", i, i, prod);
         retX = 1;
      } 
      if (fabs(sqrt(REAL_PART(h0[i]))-1) > 1e-7 && primme->procID == 0) {
         fprintf(stderr, "Warning: |Evec[%d]|-1 = %-3E\n", i, fabs(sqrt(REAL_PART(h0[i]))-1));
         retX = 1;
      } 
       
      /* Check |V(:,i)'A*V(:,i) - evals[i]| < |r|*|A| */
      primme->matrixMatvec(&evecs[primme->nLocal*i], Ax, &one, primme);
      eval0 = REAL_PART(primme_dot(&evecs[primme->nLocal*i], Ax, primme));
      if (fabs(evals[i] - eval0) > rnorms[i]*primme->aNorm && primme->procID == 0) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E should be close to %-22.1E\n", i, evals[i], eval0);
         retX = 1;
      }
      /* Check |A*V(:,i) - (V(:,i)'A*V(:,i))*V(:,i)| < |r| */
      for (j=0; j<primme->nLocal; j++) r[j] = Ax[j] - evals[i]*evecs[primme->nLocal*i+j];
      rnorm0 = sqrt(REAL_PART(primme_dot(r, r, primme)));
      if (fabs(rnorms[i]-rnorm0) > 6*max(primme->aNorm,fabs(evals[i]))*MACHINE_EPSILON && primme->procID == 0) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E, residual | %5E - %5E | <= %5E\n", i, evals[i], rnorms[i], rnorm0, 4*max(primme->aNorm,fabs(evals[i]))*MACHINE_EPSILON);
         retX = 1;
      }
      if (rnorm0 > primme->eps*primme->aNorm*sqrt((double)(i+1)) && primme->procID == 0) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E, RR residual %5E is larger than tolerance %5E\n", i, evals[i], rnorm0, primme->eps*primme->aNorm*sqrt((double)(i+1)));
         retX = 1;
      }
      /* Check angle X and V(:,i) is less than twice the max angle of the eigenvector with largest residual  */
      SUF(Num_gemv)("C", primme->nLocal, cols, COMPLEXZV(1.0), COMPLEXZ(X), primme->nLocal, COMPLEXZ(&evecs[primme->nLocal*i]), 1, COMPLEXZV(0.), COMPLEXZ(h), 1);
      if (primme->globalSumDouble) {
         int cols0 = cols*sizeof(PRIMME_NUM)/sizeof(double);
         primme->globalSumDouble(h, h0, &cols0, primme);
      }
      else h0 = h;
      prod = REAL_PARTZ(SUF(Num_dot)(cols, COMPLEXZ(h0), 1, COMPLEXZ(h0), 1));
      bound = primme->aNorm*primme->eps/delta;
      if ((sqrt(2.0)*prod+1.0)/(sqrt(2.0)*bound+1.0) < (sqrt(2.0)*prod - 1.0)/(1.0 - sqrt(2.0)*bound) && primme->procID == 0) {
         fprintf(stderr, "Warning: Eval[%d] = %-22.15E not found on X, cos angle = %5E, delta = %5E\n", i, evals[i], prod, delta);
         retX = 1;
      }
   }
   free(h);
   free(X);
   free(r);
   free(Ax);

   return retX; 
}

#undef __FUNCT__
#define __FUNCT__ "readBinaryEvecsAndPrimmeParams"
int readBinaryEvecsAndPrimmeParams(const char *fileName, PRIMME_NUM *X, PRIMME_NUM **Xout,
                                          int n, int Xcols, int *Xcolsout, int nLocal,
                                          int *perm, primme_params *primme_out) {

#  define FREAD(A, B, C, D) { ASSERT_MSG(fread(A, B, C, D) == (size_t)C, -1, "Unexpected end of file\n"); }

   FILE *f;
   PRIMME_NUM d;
   int i, j, cols;

   ASSERT_MSG((f = fopen(fileName, "rb")),
                  -1, "Could not open file %s\n", fileName);

   /* Check number size */
   /* NOTE: 2*IMAGINARY*IMAGINARY+1 is -1 in complex arith and 1 in real arith */
   FREAD(&d, sizeof(d), 1, f);
   ASSERT_MSG((int)(REAL_PART(d*(2.*IMAGINARY*IMAGINARY + 1.))) == (int)sizeof(d),
                  -1, "Mismatch arithmetic in file %s\n", fileName);
   /* Check matrix size */
   FREAD(&d, sizeof(d), 1, f);
   ASSERT_MSG(((int)REAL_PART(d)) == n,
                  -1, "Mismatch matrix size in file %s\n", fileName);

   /* Read X */
   FREAD(&d, sizeof(d), 1, f); cols = REAL_PART(d);
   if (Xcols > 0 && (X || Xout)) {
      if (!X) *Xout = X = (PRIMME_NUM*)malloc(sizeof(PRIMME_NUM)*min(cols, Xcols)*nLocal);
      if (Xcolsout) *Xcolsout = min(cols, Xcols);
      if (!perm) {
         for (i=0; i<min(cols, Xcols); i++) {
            fseek(f, (i*n + 3)*sizeof(d), SEEK_SET);
            FREAD(&X[nLocal*i], sizeof(d), nLocal, f);
         }
      }
      else {
         for (i=0; i<min(cols, Xcols); i++) {
            for (j=0; j<nLocal; j++) {
               fseek(f, (i*n + perm[j] + 3)*sizeof(d), SEEK_SET);
               FREAD(&X[nLocal*i+j], sizeof(d), 1, f);
            }
         }
      }
   }
   fseek(f, (cols*n + 3)*sizeof(d), SEEK_SET);

   /* Read primme_params */
   if (primme_out) {
      FREAD(&d, sizeof(d), 1, f);
      if ((int)REAL_PART(d) == (int)sizeof(*primme_out)) {
         FREAD(primme_out, sizeof(*primme_out), 1, f);
      }
      else
         primme_out->n = 0;
   }

   fclose(f);
   return 0;

#  undef FREAD
}

#undef __FUNCT__
#define __FUNCT__ "writeBinaryEvecsAndPrimmeParams"
int writeBinaryEvecsAndPrimmeParams(const char *fileName, PRIMME_NUM *X, int *perm,
                                           primme_params *primme) {

#  define FWRITE(A, B, C, D) { ASSERT_MSG(fwrite(A, B, C, D) == (size_t)C, -1, "Unexpected error writing on %s\n", fileName); }

   FILE *f;
   PRIMME_NUM d;
   int i, j;

   ASSERT_MSG((f = fopen(fileName, "wb")),
                  -1, "Could not open file %s\n", fileName);

   /* Write number size */
   if (primme->procID == 0) {
      /* NOTE: 2*IMAGINARY*IMAGINARY+1 is -1 in complex arith and 1 in real arith */
      d = (2.*IMAGINARY*IMAGINARY + 1.)*sizeof(d);
      FWRITE(&d, sizeof(d), 1, f);
      /* Write matrix size */
      d = primme->n;
      FWRITE(&d, sizeof(d), 1, f);
      /* Write number of columns */
      d = primme->initSize;
      FWRITE(&d, sizeof(d), 1, f);
   }

   /* Write X */
   if (!perm) {
      for (i=0; i<primme->initSize; i++) {
         fseek(f, (i*primme->n + 3)*sizeof(d), SEEK_SET);
         FWRITE(&X[primme->nLocal*i], sizeof(d), primme->nLocal, f);
      }
   }
   else {
      for (i=0; i<primme->initSize; i++) {
         for (j=0; j<primme->nLocal; j++) {
            fseek(f, (i*primme->n + perm[j] + 3)*sizeof(d), SEEK_SET);
            FWRITE(&X[primme->nLocal*i+j], sizeof(d), 1, f);
         }
      }
   }

   /* Write primme_params */
   if (primme->procID == 0) {
      fseek(f, sizeof(d)*(primme->n*primme->initSize + 3), SEEK_SET);
      d = sizeof(*primme);
      FWRITE(&d, sizeof(d), 1, f);
      FWRITE(primme, sizeof(*primme), 1, f);
   }

   fclose(f);
   return 0;

#  undef FWRITE
}

#undef __FUNCT__
#define __FUNCT__ "check_solution_svds"
int check_solution_svds(const char *checkXFileName, primme_svds_params *primme_svds, double *svals,
                        PRIMME_NUM *svecs, double *rnorms, int *perm) {

   double sval0, rnorm0, prod, delta, bound;
   PRIMME_NUM *Ax, *r, *X=NULL, *h, *h0, *U, *V;
   int i, j, cols, retX=0, one=1, notrans=0, trans=1;
   primme_svds_params primme_svds0;

   /* Read stored singular vectors and primme_svds_params */
   ASSERT_MSG(readBinaryEvecsAndPrimmeSvdsParams(checkXFileName, NULL, &X, primme_svds->m,
      primme_svds->n, primme_svds->m, &cols, primme_svds->mLocal, primme_svds->nLocal,
      perm, &primme_svds0) == 0, -1, "");

   /* Check primme_svds_params */
#  define CHECK_PRIMME_PARAM(F) \
        if (primme_svds0. F != primme_svds-> F ) { \
           fprintf(stderr, "Warning: discrepancy in primme_svds." #F ", %d should be close to %d\n", primme_svds-> F , primme_svds0. F ); \
           retX = 1; \
        }
#  define CHECK_PRIMME_PARAM_DOUBLE(F) \
        if (fabs(primme_svds0. F - primme_svds-> F) > primme_svds-> F * 1e-14) { \
           fprintf(stderr, "Warning: discrepancy in primme_svds." #F ", %.16e should be close to %.16e\n", primme_svds-> F , primme_svds0. F ); \
           retX = 1; \
        }
#  define CHECK_PRIMME_PARAM_TOL(F, T) \
        if (abs(primme_svds0. F - primme_svds-> F ) > primme_svds-> F * T /100+1) { \
           fprintf(stderr, "Warning: discrepancy in primme_svds." #F ", %d should be close to %d\n", primme_svds-> F , primme_svds0. F ); \
           retX = 1; \
        }

   if (primme_svds0.n) {
      CHECK_PRIMME_PARAM(m);
      CHECK_PRIMME_PARAM(n);
      CHECK_PRIMME_PARAM(numSvals);
      CHECK_PRIMME_PARAM(target);
      CHECK_PRIMME_PARAM(numTargetShifts);
      CHECK_PRIMME_PARAM(locking);
      CHECK_PRIMME_PARAM(numOrthoConst);
      CHECK_PRIMME_PARAM(maxBasisSize);
      CHECK_PRIMME_PARAM(maxBlockSize);
      CHECK_PRIMME_PARAM_DOUBLE(aNorm);
      CHECK_PRIMME_PARAM_DOUBLE(eps);
      CHECK_PRIMME_PARAM(initSize);
      CHECK_PRIMME_PARAM_TOL(stats.numMatvecs, 60);
      CHECK_PRIMME_PARAM(method);
      CHECK_PRIMME_PARAM(methodStage2);
   }

#  undef CHECK_PRIMME_PARAM
#  undef CHECK_PRIMME_PARAM_DOUBLE
#  undef CHECK_PRIMME_PARAM_TOL

   h = (PRIMME_NUM *)primme_calloc(cols*2, sizeof(PRIMME_NUM), "h"); h0 = &h[cols];
   Ax = (PRIMME_NUM *)primme_calloc(max(primme_svds->mLocal, primme_svds->nLocal), sizeof(PRIMME_NUM), "Ax");
   r = (PRIMME_NUM *)primme_calloc(max(primme_svds->mLocal, primme_svds->nLocal), sizeof(PRIMME_NUM), "r");

   /* Estimate the separation between eigenvalues */
   delta = primme_svds->aNorm;
   for (i=1; i < primme_svds->initSize; i++) {
      delta = min(delta, fabs(svals[i]-svals[i-1]));
   }

   U = svecs;
   V = &svecs[primme_svds->mLocal*cols];   
   for (i=0; i < primme_svds->initSize; i++) {
      /* Check normality of U(:,i) and V(:,i) */
      sval0 = REAL_PART(primme_svds_dot(&U[primme_svds->mLocal*i], &U[primme_svds->mLocal*i], 0, primme_svds));
      if (fabs(1.0 - sval0) > 1e-8 && primme_svds->procID == 0) {
         fprintf(stderr, "Warning: norm of U[%d] = %e\n", i, sval0);
         retX = 1;
      }
      sval0 = REAL_PART(primme_svds_dot(&V[primme_svds->nLocal*i], &V[primme_svds->nLocal*i], 1, primme_svds));
      if (fabs(1.0 - sval0) > 1e-8 && primme_svds->procID == 0) {
         fprintf(stderr, "Warning: norm of V[%d] = %e\n", i, sval0);
         retX = 1;
      }
      /* Check |U(:,i)'A*V(:,i) - svals[i]| < |r|*|A| */
      primme_svds->matrixMatvec(&V[primme_svds->nLocal*i], &primme_svds->nLocal, Ax, &primme_svds->mLocal, &one, &notrans, primme_svds);
      sval0 = REAL_PART(primme_svds_dot(&U[primme_svds->mLocal*i], Ax, 0, primme_svds));
      if (fabs(svals[i] - sval0) > rnorms[i]*primme_svds->aNorm && primme_svds->procID == 0) {
         fprintf(stderr, "Warning: Sval[%d] = %-22.15E should be close to %-22.1E\n", i, svals[i], sval0);
         retX = 1;
      }
      /* Check |A*V(:,i) - (U(:,i)'A*V(:,i))*U(:,i)|^2 + |A'*U(:,i) - (U(:,i)'A*V(:,i))*V(:,i)|^2 < |r|^2 */
      for (j=0; j<primme_svds->mLocal; j++) r[j] = Ax[j] - svals[i]*U[primme_svds->mLocal*i+j];
      rnorm0 = REAL_PART(primme_svds_dot(r, r, 0, primme_svds));
      primme_svds->matrixMatvec(&U[primme_svds->mLocal*i], &primme_svds->mLocal, Ax, &primme_svds->nLocal, &one, &trans, primme_svds);
      for (j=0; j<primme_svds->nLocal; j++) r[j] = Ax[j] - svals[i]*V[primme_svds->nLocal*i+j];
      rnorm0 += REAL_PART(primme_svds_dot(r, r, 1, primme_svds));
      rnorm0 = sqrt(rnorm0);
      if (rnorms[i] < rnorm0 - 4*max(primme_svds->aNorm,fabs(svals[i]))*MACHINE_EPSILON && primme_svds->procID == 0) {
         fprintf(stderr, "Warning: rnorms[%d] = %5E, but the computed residual is %5E\n", i, rnorms[i], rnorm0);
         retX = 1;
      }
      if (rnorm0 > 8*primme_svds->eps*primme_svds->aNorm*sqrt((double)(i+1)) && primme_svds->procID == 0) {
         fprintf(stderr, "Warning: Sval[%d] = %-22.15E, RR residual %5E is larger than tolerance %5E\n", i, svals[i], rnorm0, primme_svds->eps*primme_svds->aNorm*sqrt((double)(i+1)));
         retX = 1;
      }
      /* Check angle X and U(:,i) is less than twice the max angle of the eigenvector with largest residual  */
      SUF(Num_gemv)("C", primme_svds->mLocal, cols, COMPLEXZV(1.0), COMPLEXZ(X), primme_svds->mLocal, COMPLEXZ(&svecs[primme_svds->mLocal*i]), 1, COMPLEXZV(0.), COMPLEXZ(h), 1);
      if (primme_svds->globalSumDouble) {
         int cols0 = cols*sizeof(PRIMME_NUM)/sizeof(double);
         primme_svds->globalSumDouble(h, h0, &cols0, primme_svds);
      }
      else h0 = h;
      prod = REAL_PARTZ(SUF(Num_dot)(cols, COMPLEXZ(h0), 1, COMPLEXZ(h0), 1));
      bound = primme_svds->aNorm*primme_svds->eps/delta;
      if ((sqrt(2.0)*prod+1.0)/(sqrt(2.0)*bound+1.0) < (sqrt(2.0)*prod - 1.0)/(1.0 - sqrt(2.0)*bound) && primme_svds->procID == 0) {
         fprintf(stderr, "Warning: Sval[%d] = %-22.15E not found on X, cos angle = %5E, delta = %5E\n", i, svals[i], prod, delta);
         retX = 1;
      }
   }
   free(h);
   free(X);
   free(r);
   free(Ax);

   return retX; 
}

#undef __FUNCT__
#define __FUNCT__ "readBinaryEvecsAndPrimmeSvdsParams"
int readBinaryEvecsAndPrimmeSvdsParams(const char *fileName, PRIMME_NUM *X, PRIMME_NUM **Xout,
                                       int m, int n, int Xcols, int *Xcolsout, int mLocal, int nLocal,
                                       int *perm, primme_svds_params *primme_svds_out) {

#  define FREAD(A, B, C, D) { ASSERT_MSG(fread(A, B, C, D) == (size_t)C, -1, "Unexpected end of file\n"); }

   FILE *f;
   PRIMME_NUM d;
   int i, j, cols;

   ASSERT_MSG((f = fopen(fileName, "rb")),
                  -1, "Could not open file %s\n", fileName);

   /* Check number size */
   /* NOTE: 2*IMAGINARY*IMAGINARY+1 is -1 in complex arith and 1 in real arith */
   FREAD(&d, sizeof(d), 1, f);
   ASSERT_MSG((int)(REAL_PART(d*(2.*IMAGINARY*IMAGINARY + 1.))) == (int)sizeof(d),
                  -1, "Mismatch arithmetic in file %s\n", fileName);
   /* Check matrix size */
   FREAD(&d, sizeof(d), 1, f);
   ASSERT_MSG(((int)REAL_PART(d)) == m,
                  -1, "Mismatch matrix size in file %s\n", fileName);
   FREAD(&d, sizeof(d), 1, f);
   ASSERT_MSG(((int)REAL_PART(d)) == n,
                  -1, "Mismatch matrix size in file %s\n", fileName);

   /* Read X = [U V] */
   FREAD(&d, sizeof(d), 1, f); cols = REAL_PART(d);
   if (Xcols > 0 && (X || Xout)) {
      if (!X) *Xout = X = (PRIMME_NUM*)malloc(sizeof(PRIMME_NUM)*min(cols, Xcols)*(mLocal+nLocal));
      if (Xcolsout) *Xcolsout = min(cols, Xcols);
      if (!perm) {
         assert(n == nLocal && m == mLocal);
         FREAD(X, sizeof(d), mLocal*min(cols, Xcols), f);
         fseek(f, (cols*m + 4)*sizeof(d), SEEK_SET);
         FREAD(&X[mLocal*min(cols, Xcols)], sizeof(d), nLocal*min(cols, Xcols), f);
      }
      else {
         for (i=0; i<min(cols, Xcols); i++) {
            for (j=0; j<mLocal; j++) {
               fseek(f, (i*m + perm[j] + 4)*sizeof(d), SEEK_SET);
               FREAD(&X[mLocal*i+j], sizeof(d), 1, f);
            }
         }
         for (i=0; i<min(cols, Xcols); i++) {
            for (j=0; j<nLocal; j++) {
               fseek(f, (min(cols, Xcols)*m + i*n + perm[j] + 4)*sizeof(d), SEEK_SET);
               FREAD(&X[min(cols, Xcols)*mLocal + nLocal*i+j], sizeof(d), 1, f);
            }
         }
      }
   }
   fseek(f, (cols*(m+n) + 4)*sizeof(d), SEEK_SET);

   /* Read primme_params */
   if (primme_svds_out) {
      FREAD(&d, sizeof(d), 1, f);
      if ((int)REAL_PART(d) == (int)sizeof(*primme_svds_out)) {
         FREAD(primme_svds_out, sizeof(*primme_svds_out), 1, f);
      }
      else
         primme_svds_out->m = primme_svds_out-> n = 0;
   }

   fclose(f);
   return 0;

#  undef FREAD
}


#undef __FUNCT__
#define __FUNCT__ "writeBinaryEvecsAndPrimmeSvdsParams"
int writeBinaryEvecsAndPrimmeSvdsParams(const char *fileName, PRIMME_NUM *X, int *perm,
                                    primme_svds_params *primme_svds) {

#  define FWRITE(A, B, C, D) { ASSERT_MSG(fwrite(A, B, C, D) == (size_t)C, -1, "Unexpected error writing on %s\n", fileName); }

   FILE *f;
   PRIMME_NUM d;
   int i, j;

   ASSERT_MSG((f = fopen(fileName, "wb")),
                  -1, "Could not open file %s\n", fileName);

   /* Write number size */
   if (primme_svds->procID == 0) {
      /* NOTE: 2*IMAGINARY*IMAGINARY+1 is -1 in complex arith and 1 in real arith */
      d = (2.*IMAGINARY*IMAGINARY + 1.)*sizeof(d);
      FWRITE(&d, sizeof(d), 1, f);
      /* Write matrix size */
      d = primme_svds->m;
      FWRITE(&d, sizeof(d), 1, f);
      d = primme_svds->n;
      FWRITE(&d, sizeof(d), 1, f);
      /* Write number of columns */
      d = primme_svds->initSize;
      FWRITE(&d, sizeof(d), 1, f);
   }

   /* Write X = [U V] */
   if (!perm) {
      assert(primme_svds->n == primme_svds->nLocal && primme_svds->m == primme_svds->mLocal);
      FWRITE(X, sizeof(d), (primme_svds->m+primme_svds->n)*primme_svds->initSize, f);
   }
   else {
      for (i=0; i<primme_svds->initSize; i++) {
         for (j=0; j<primme_svds->mLocal; j++) {
            fseek(f, (i*primme_svds->m + perm[j] + 4)*sizeof(d), SEEK_SET);
            FWRITE(&X[primme_svds->mLocal*i+j], sizeof(d), 1, f);
         }
      }
      for (i=0; i<primme_svds->initSize; i++) {
         for (j=0; j<primme_svds->nLocal; j++) {
            fseek(f, (primme_svds->m*primme_svds->initSize + i*primme_svds->n + perm[primme_svds->mLocal+j] + 4)*sizeof(d), SEEK_SET);
            FWRITE(&X[primme_svds->mLocal*primme_svds->initSize+primme_svds->nLocal*i+j], sizeof(d), 1, f);
         }
      }
   }

   /* Write primme_svds_params */
   if (primme_svds->procID == 0) {
      fseek(f, sizeof(d)*((primme_svds->m+primme_svds->n)*primme_svds->initSize + 4), SEEK_SET);
      d = sizeof(*primme_svds);
      FWRITE(&d, sizeof(d), 1, f);
      FWRITE(primme_svds, sizeof(*primme_svds), 1, f);
   }

   fclose(f);
   return 0;

#  undef FWRITE
}
