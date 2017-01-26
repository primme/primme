/*******************************************************************************
 * Copyright (c) 2016, College of William & Mary
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
 * File: primmeR.c
 * 
 * Purpose - PRIMME R interface.
 * 
 * For details about PRIMME parameters, methods, and settings see ../readme.txt
 *
 ******************************************************************************/

#include <Rcpp.h>
#include <algorithm>
#include "primme.h"
#include "Primme_types.h"

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
//
// Auxiliary functions and macros
//

#define ASSERT(X) if(X) stop("This shouldn't happen (" #X ")");

// Template version of dprimme and zprimme

static int tprimme(double *evals, double *evecs, double *resNorms, primme_params *primme) {
   return dprimme(evals, evecs, resNorms, primme);
}
static int tprimme(double *evals, std::complex<double> *evecs, double *resNorms, primme_params *primme) {
   return zprimme(evals, evecs, resNorms, primme);
}

// Template version of dprimme_svds and zprimme_svds

static int tprimme_svds(double *svals, double *svecs, double *resNorms, primme_svds_params *primme_svds) { 
   return dprimme_svds(svals, svecs, resNorms, primme_svds);
}
static int tprimme_svds(double *svals, std::complex<double> *svecs, double *resNorms, primme_svds_params *primme_svds) {
   return zprimme_svds(svals, svecs, resNorms, primme_svds);
}


// Create a new Rcpp Matrix of type S copying the content of a Fortran matrix of
// type T.
// Arguments:
// - x: Fortran matrix from copy the values
// - m: number of rows of x and the output matrix
// - n: number of columns of x and the output matrix
// - ld: leading dimension of x
// Return: Matrix<S>

template <typename T, typename S>
S createMatrix(T *x, PRIMME_INT m, int n, PRIMME_INT ld) {
   if (ld == m) {
      return S((int)m, (int)n, x);
   }
   S mat((int)m, (int)n);
   for (int i=0; i<n; i++) {
      std::copy(&x[ld*i], &x[ld*i+m], &mat(0, i));
   }
   return mat;
}
      
// Copy the content of a Rcpp Matrix of type S into a Fortran matrix of type T.
// Arguments:
// - mat: input matrix
// - x: Fortran matrix to copy the values
// - m: number of rows of x and mat
// - n: number of columns of x and mat
// - ld: leading dimension of x
// - checkDimensions: if true, stop if mat hasn't dimensions m x n
// Return: Matrix<S>

template <typename T, typename S>
void copyMatrix(S mat, T *x, PRIMME_INT m, int n, PRIMME_INT ld,
      bool checkDimensions=true) {

   if (checkDimensions && (mat.rows() != m || mat.cols() != n))
      stop("expected matrix with different dimensions");
   if (ld == m) {
      std::copy(mat.begin(), mat.end(), x);
   }
   else {
      for (int i=0; i<n; i++) {
         std::copy(&mat(0, i), &mat(m, i), &x[ld*i]);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
//
// R wrappers around function in PRIMME
//


// NOTE: to indicate that these functions are for internal use in the package
//       the function name starts with ".".

// [[Rcpp::export(.primme_initialize)]]
PrimmeParams primme_initialize_rcpp() {
   primme_params *primme = new primme_params;
   primme_initialize(primme);
   return PrimmeParams(primme);
}

// [[Rcpp::export(.primme_free)]]
void primme_free_rcpp(PrimmeParams primme) {
   if (primme->targetShifts) delete [] primme->targetShifts;
   primme_free(primme);
}

// [[Rcpp::export(.primme_set_method)]]
void primme_set_method_rcpp(std::string methodstr, PrimmeParams primme) {
   int method;
   if (primme_constant_info(methodstr.c_str(), &method))
      stop("method isn't valid");
   primme_set_method((primme_preset_method)method, primme);
}

// [[Rcpp::export(.primme_get_member)]]
SEXP primme_get_member_rcpp(std::string labelstr, PrimmeParams primme) {
   primme_params_label label;
   const char *labelname = labelstr.c_str();
   primme_type ptype;
   int arity;
   if (primme_member_info(&label, &labelname, &ptype, &arity))
      stop("invalid label");

   switch(label) {
      // Get members with arity > 1

      case PRIMME_iseed:
      {
         IntegerVector v(4);
         std::copy(primme->iseed, &primme->iseed[4], v.begin());
         return v;
      }

      case PRIMME_targetShifts:
      {
         NumericVector v(primme->numTargetShifts);
         std::copy(primme->targetShifts, &primme->targetShifts[primme->numTargetShifts], v.begin());
         return v;
      }

      // Forbidden members
 
      case PRIMME_numProcs:
      case PRIMME_procID:
      case PRIMME_commInfo:
      case PRIMME_nLocal:
      case PRIMME_globalSumReal:
      case PRIMME_numTargetShifts:
      case PRIMME_intWorkSize:
      case PRIMME_realWorkSize:
      case PRIMME_intWork:
      case PRIMME_realWork:
      case PRIMME_outputFile:
      case PRIMME_matrix:
      case PRIMME_preconditioner:
      case PRIMME_convTestFun:
      case PRIMME_ldevecs:
      case PRIMME_ldOPs:
      case PRIMME_massMatrixMatvec:
      case PRIMME_matrixMatvec:
      case PRIMME_applyPreconditioner:
         stop("Unsupported to get this option");
         break;

      default : 
      {
         ASSERT(arity != 1);

         // Set members with type int

         if (ptype == primme_int) {
            PRIMME_INT v;
            ASSERT(primme_get_member(primme, label, &v));
            return wrap((int)v);
         }

         // Set members with type double

         else if (ptype == primme_double) {
            double v;
            ASSERT(primme_get_member(primme, label, &v));
            return wrap(v);
         }

         else {
            /* This shouldn't happen */
            ASSERT(1);
         }
      }
   }
}

// [[Rcpp::export(.primme_set_member)]]
void primme_set_member_rcpp(std::string labelstr, SEXP value, PrimmeParams primme) {
   primme_params_label label;
   const char *labelname = labelstr.c_str();
   primme_type ptype;
   int arity;
   if (primme_member_info(&label, &labelname, &ptype, &arity))
      stop("invalid label");

   switch(label) {
      // Set members with arity > 1

      case PRIMME_iseed:
      {
         IntegerVector v = as<IntegerVector>(value);
         if (v.size() != 4)
            stop("value should have four elements");
         std::copy(v.begin(), v.end(), primme->iseed);
         break;
      }

      case PRIMME_targetShifts:
      {
         NumericVector v = as<NumericVector>(value);
         if (primme->targetShifts) delete [] primme->targetShifts;
         primme->targetShifts = new double[v.size()];
         primme->numTargetShifts = v.size();
         std::copy(v.begin(), v.end(), primme->targetShifts);
         break;
      }

      // Forbidden members
 
      case PRIMME_numProcs:
      case PRIMME_procID:
      case PRIMME_commInfo:
      case PRIMME_nLocal:
      case PRIMME_globalSumReal:
      case PRIMME_numTargetShifts:
      case PRIMME_intWorkSize:
      case PRIMME_realWorkSize:
      case PRIMME_intWork:
      case PRIMME_realWork:
      case PRIMME_outputFile:
      case PRIMME_matrix:
      case PRIMME_preconditioner:
      case PRIMME_convTestFun:
      case PRIMME_ldevecs:
      case PRIMME_ldOPs:
      case PRIMME_massMatrixMatvec:
      case PRIMME_matrixMatvec:
      case PRIMME_applyPreconditioner:
         stop("Unsupported to set this option");
         break;

      default : 
      {
         ASSERT(arity != 1);

         // Set members with type int

         if (ptype == primme_int) {
            PRIMME_INT v;
            if (TYPEOF(value) == STRSXP) {
               int v0;
               if (primme_constant_info(as<std::string>(value).c_str(), &v0)) {
                  stop("Invalid value");
               }
               v = v0;
            }
            else {
               v = as<int>(value);
            }
            ASSERT(primme_set_member(primme, label, &v));
         }

         // Set members with type double

         else if (ptype == primme_double) {
            double v = as<double>(value);
            ASSERT(primme_set_member(primme, label, &v));
         }

         else {
            /* This shouldn't happen */
            ASSERT(1);
         }
      }
   }
}

// Auxiliary functions for xprimme; they return the value of
// some option in primme_params

struct getMatrixField {
   static void *get(primme_params *primme) {
      return primme->matrix;
   }
};

struct getPreconditionerField {
   static void* get(primme_params *primme) {
      return primme->preconditioner;
   }
};

struct getConvTestField {
   static void* get(primme_params *primme) {
      return primme->convtest;
   }
};

// Auxiliary function for xprimme; PRIMME wrapper around
// matrixMatvec, massMatrixMatvec and applyPreconditioner. Create a Matrix<S>
// from input vector x, call the function handler returned by F(primme) and
// copy the content of its returned Matrix<S> into the output vector y.
// Arguments:
// - T: type of PRIMME evecs
// - S: R type
// - TS: type of elements in Matrix<S>
// - F: F::get(primme) return the function pointer
// - x, ldx, y, ldy, ...: arguments of matrixMatvec, applyPreconditioner...

template <typename T, int S, typename TS, typename F>
static void matrixMatvecEigs(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, struct primme_params *primme, int *ierr)
{  
   // Create input vector
   Matrix<S,NoProtectStorage> vx =
      createMatrix<TS,Matrix<S> >((TS*)x, primme->nLocal,
            *blockSize, *ldx);

   // Call the callback
   Function *f = (Function*)F::get(primme);
   Matrix<S> vy = (*f)(vx);

   // Copy output vector
   copyMatrix(vy, (TS*)y, primme->nLocal, *blockSize, *ldy);

   *ierr = 0;
}

// Auxiliary function for xprimme; PRIMME wrapper around convTestFun.
// Create a Vector<S>
// from input vector x, call the function handler returned by F(primme) and
// copy the content of its returned Matrix<S> into the output vector y.
// Arguments:
// - T: type of PRIMME evecs
// - S: R type
// - TS: type of elements in Vector<S>
// - F: F::get(primme) return the function pointer
// - x, ldx, y, ldy, ...: arguments of convTestFun

template <typename T, int S, typename TS, typename F>
static void convTestFunEigs(double *eval, void *evec, double *rNorm, int *isconv, 
      struct primme_params *primme, int *ierr) {

   // Pass objects to R types
   Vector<S,NoProtectStorage>
      sevec(evec?primme->nLocal:0, *(TS*)evec);
   Vector<REALSXP,NoProtectStorage>
      seval(eval?1:0, *eval),
      srnorm(rNorm?1:0, *rNorm);

   // Call the callback
   Function *f = (Function*)F::get(primme);
   *isconv = as<bool>((*f)(seval, sevec, srnorm)) ? 1 : 0;

   *ierr = 0;
}

// Generic function for dprimme and zprimme
// Arguments:
// - T: type of PRIMME evecs
// - S: Matrix<S> is the type of the output eigenvectors
// - TS: type of elements in Matrix<S>
// - ortho: orthogonal constrains
// - init: initial guesses
// - A: matrix-vector product
// - B: mass matrix-vector product
// - prec: preconditioner application
// - convTest: convergence criterion

template<typename T, int S, typename TS>
static List xprimme(Matrix<S> ortho, Matrix<S> init, Function A, SEXP B,
      SEXP prec, SEXP convTest, PrimmeParams primme)
{
   if (primme->nLocal == 0) primme->nLocal = primme->n;

   // Check dimensions of ortho and init

   if (ortho.rows() != 0 && ortho.rows() != primme->nLocal)
      stop("Invalid number of rows in input matrix ortho");
   if (init.rows() != 0 && init.rows() != primme->nLocal)
      stop("Invalid number of rows in input matrix init");

   int ncols = ortho.cols() + std::max(init.cols(), primme->numEvals);

   // Allocate evals, rnorms and evecs

   NumericVector vevals(primme->numEvals);
   NumericVector vrnorms(primme->numEvals);
   Matrix<S> vevecs(primme->nLocal, ncols);
   double *evals = vevals.begin(), *rnorms = vrnorms.begin();
   T *evecs = (T*)vevecs.begin();

   // Copy initial and orthogonal constrain vectors

   primme->numOrthoConst = ortho.cols();
   copyMatrix(ortho, (TS*)evecs, primme->nLocal, ortho.cols(), primme->nLocal, false);
   primme->initSize = init.cols();
   copyMatrix(init, (TS*)&evecs[primme->nLocal*ortho.cols()], primme->nLocal, init.cols(), primme->nLocal, false);

   // Set matvec and preconditioner

   primme->matrix = &A;
   primme->matrixMatvec = matrixMatvecEigs<T, S, TS, getMatrixField>;

   Function fprec = A; /* dummy set value */
   if (prec != R_NilValue) {
      fprec = as<Function>(prec);
      primme->preconditioner = &fprec;
      primme->applyPreconditioner = matrixMatvecEigs<T, S, TS, getPreconditionerField>;
   }

   if (B != R_NilValue) {
      stop("Unsupported generalized eigenvalue problems, for now");
   }

   Function fconvTest = A; /* dummy set value */
   if (convTest != R_NilValue) {
      fconvTest = as<Function>(convTest);
      primme->convtest = &fconvTest;
      primme->convTestFun = convTestFunEigs<T, S, TS, getConvTestField>;
   }

   // Call xprimme
   int ret = tprimme(evals, evecs, rnorms, primme);

   // Return only the eigenvectors
   SubMatrix<S> revecs(vevecs, Range(0, primme->nLocal-1), Range(ortho.cols(), ortho.cols()+primme->initSize-1));

   // Return
   return List::create(
      Named("ret") = IntegerVector::create(ret),
      Named("evals") = vevals,
      Named("evecs") = revecs,
      Named("rnorms") = vrnorms);
}

// [[Rcpp::export(.dprimme)]]
List dprimme_rcpp(NumericMatrix ortho, NumericMatrix init, Function A, SEXP B, SEXP prec, SEXP convTest, PrimmeParams primme) {
   return xprimme<double, REALSXP, double>(ortho, init, A, B, prec, convTest, primme);
}

// [[Rcpp::export(.zprimme)]]
List zprimme_rcpp(ComplexMatrix ortho, ComplexMatrix init, Function A, SEXP B, SEXP prec, SEXP convTest, PrimmeParams primme) {
   return xprimme<PRIMME_COMPLEX_DOUBLE, CPLXSXP, Rcomplex>(ortho, init, A, B, prec, convTest, primme);
}

// [[Rcpp::export(.primme_svds_initialize)]]
PrimmeSvdsParams primme_svds_initialize_rcpp() {
   primme_svds_params *primme_svds = new primme_svds_params;
   primme_svds_initialize(primme_svds);
   return PrimmeSvdsParams(primme_svds);
}

// [[Rcpp::export(.primme_svds_free)]]
void primme_svds_free_rcpp(PrimmeSvdsParams primme_svds) {
   if (primme_svds->targetShifts) delete [] primme_svds->targetShifts;
   primme_svds_free(primme_svds);
}

// [[Rcpp::export(.primme_svds_set_method)]]
void primme_svds_set_method_rcpp(std::string methodstr,
      std::string methodStage1str, std::string methodStage2str,
      PrimmeSvdsParams primme_svds) {

   int method, methodStage1, methodStage2;
   if (primme_svds_constant_info(methodstr.c_str(), &method))
      stop("method isn't valid");
   if (primme_constant_info(methodStage1str.c_str(), &methodStage1))
      stop("methodStage1 isn't valid");
   if (primme_constant_info(methodStage2str.c_str(), &methodStage2))
      stop("methodStage2 isn't valid");
   primme_svds_set_method((primme_svds_preset_method)method,
      (primme_preset_method)methodStage1, (primme_preset_method)methodStage2,
      primme_svds);
}

// [[Rcpp::export(.primme_svds_get_member)]]
SEXP primme_svds_get_member_rcpp(std::string labelstr,
      PrimmeSvdsParams primme_svds) {

   primme_svds_params_label label;
   const char *labelname = labelstr.c_str();
   primme_type ptype;
   int arity;
   if (primme_svds_member_info(&label, &labelname, &ptype, &arity))
      stop("invalid label");

   switch(label) {
      // Get members with arity > 1

      case PRIMME_SVDS_iseed:
      {
         IntegerVector v(4);
         std::copy(primme_svds->iseed, &primme_svds->iseed[4], v.begin());
         return v;
      }

      case PRIMME_SVDS_targetShifts:
      {
         NumericVector v(primme_svds->numTargetShifts);
         std::copy(primme_svds->targetShifts,
               &primme_svds->targetShifts[primme_svds->numTargetShifts],
               v.begin());
         return v;
      }

      // Get primme_params references

      case PRIMME_SVDS_primme:
      { 
         return PrimmeParams(&primme_svds->primme,
               false /* don't delete when collected by garbage collector */);
      }
      case PRIMME_SVDS_primmeStage2:
      { 
         return PrimmeParams(&primme_svds->primmeStage2,
               false /* don't delete when collected by garbage collector */);
      }

      // Forbidden members
 
      case PRIMME_SVDS_numProcs: 
      case PRIMME_SVDS_procID: 
      case PRIMME_SVDS_mLocal: 
      case PRIMME_SVDS_nLocal: 
      case PRIMME_SVDS_commInfo:
      case PRIMME_SVDS_globalSumReal:
      case PRIMME_SVDS_numTargetShifts:
      case PRIMME_SVDS_intWorkSize:
      case PRIMME_SVDS_realWorkSize:
      case PRIMME_SVDS_intWork:
      case PRIMME_SVDS_realWork:
      case PRIMME_SVDS_matrix:
      case PRIMME_SVDS_preconditioner:
      case PRIMME_SVDS_outputFile:
      case PRIMME_SVDS_matrixMatvec:
      case PRIMME_SVDS_applyPreconditioner:
         stop("Unsupported to get this option");
         break;

      default : 
      {
         ASSERT(arity != 1);

         // Set members with type int

         if (ptype == primme_int) {
            PRIMME_INT v;
            ASSERT(primme_svds_get_member(primme_svds, label, &v));
            return wrap((int)v);
         }

         // Set members with type double

         else if (ptype == primme_double) {
            double v;
            ASSERT(primme_svds_get_member(primme_svds, label, &v));
            return wrap(v);
         }

         else {
            /* This shouldn't happen */
            ASSERT(1);
         }
      }
   }
}

// [[Rcpp::export(.primme_svds_set_member)]]
void primme_svds_set_member_rcpp(std::string labelstr, SEXP value,
      PrimmeSvdsParams primme_svds) {

   primme_svds_params_label label;
   const char *labelname = labelstr.c_str();
   primme_type ptype;
   int arity;
   if (primme_svds_member_info(&label, &labelname, &ptype, &arity))
      stop("invalid label");

   switch(label) {
      // Set members with arity > 1

      case PRIMME_SVDS_iseed:
      {
         IntegerVector v = as<IntegerVector>(value);
         if (v.size() != 4)
            stop("value should have four elements");
         std::copy(v.begin(), v.end(), primme_svds->iseed);
         break;
      }

      case PRIMME_SVDS_targetShifts:
      {
         NumericVector v = as<NumericVector>(value);
         if (primme_svds->targetShifts) delete [] primme_svds->targetShifts;
         primme_svds->targetShifts = new double[v.size()];
         primme_svds->numTargetShifts = v.size();
         std::copy(v.begin(), v.end(), primme_svds->targetShifts);
         break;
      }

      // Forbidden members
 
      case PRIMME_SVDS_primme: 
      case PRIMME_SVDS_primmeStage2:
      case PRIMME_SVDS_numProcs: 
      case PRIMME_SVDS_procID: 
      case PRIMME_SVDS_mLocal: 
      case PRIMME_SVDS_nLocal: 
      case PRIMME_SVDS_commInfo:
      case PRIMME_SVDS_globalSumReal:
      case PRIMME_SVDS_numTargetShifts:
      case PRIMME_SVDS_intWorkSize:
      case PRIMME_SVDS_realWorkSize:
      case PRIMME_SVDS_intWork:
      case PRIMME_SVDS_realWork:
      case PRIMME_SVDS_matrix:
      case PRIMME_SVDS_preconditioner:
      case PRIMME_SVDS_outputFile:
      case PRIMME_SVDS_matrixMatvec:
      case PRIMME_SVDS_applyPreconditioner:
         stop("Unsupported to set this option");
         break;

      default : 
      {
         ASSERT(arity != 1);

         // Set members with type int

         if (ptype == primme_int) {
            PRIMME_INT v;
            if (TYPEOF(value) == STRSXP) {
               int v0;
               if (primme_svds_constant_info(as<std::string>(value).c_str(), &v0)) {
                  stop("Invalid value");
               }
               v = v0;
            }
            else {
               v = as<int>(value);
            }
            ASSERT(primme_svds_set_member(primme_svds, label, &v));
         }

         // Set members with type double

         else if (ptype == primme_double) {
            double v = as<double>(value);
            ASSERT(primme_svds_set_member(primme_svds, label, &v));
         }

         else {
            /* This shouldn't happen */
            ASSERT(1);
         }
      }
   }
}

// Auxiliary functions for mexFunction_xprimme_svds; they return the next
// information based on the mode and primme_svds_params:
// Arguments:
// - mode: the corresponding input value in primme_svds.matrixMatvec and
//         primme_svds.applyPreconditioner.
// - primme_svds: primme_svds_params
// - mx: return the number of rows of input vectors x
// - my: return the number of rows of output vectors y
// - str: return the corresponding string for mode (notransp/transp or
//        AHA/AAH/aug).

struct getSvdsForMatrix {
   static void get(int transpose, primme_svds_params *primme_svds,
         PRIMME_INT *mx, PRIMME_INT *my, void **AFUN, const char **str) {
      *AFUN = primme_svds->matrix;
      if (transpose == 0) { /* Doing y <- A * x */
         *mx = primme_svds->nLocal;
         *my = primme_svds->mLocal;
         *str = "n";
      }
      else { /* Doing y <- A' * x */
         *mx = primme_svds->mLocal;
         *my = primme_svds->nLocal;
         *str = "c";
      }
   }
};

struct getSvdsForPreconditioner {
   static void get(int mode, primme_svds_params *primme_svds,
         PRIMME_INT *mx, PRIMME_INT *my, void **AFUN, const char **str) {
      *AFUN = primme_svds->preconditioner;
      if (mode == primme_svds_op_AtA) {
         /* Preconditioner for A^t*A */
         *mx = *my = primme_svds->nLocal;
         *str = "AHA";
      }
      else if (mode == primme_svds_op_AAt) {
         /* Preconditioner for A*A^t */
         *mx = *my = primme_svds->mLocal;
         *str = "AAH";
      }
      else if (mode == primme_svds_op_augmented) {
         /* Preconditioner for [0 A^t; A 0] */
         *mx = *my = primme_svds->mLocal + primme_svds->nLocal;
         *str = "aug";
      }
      else {
         stop("Unsupported preconditioner type");
      }
   }
};


// Auxiliary function for xprimme_svds; PRIMME SVDS wrapper around
// matrixMatvec and applyPreconditioner. Create a Matrix<S>
// from input vector x, call the function handler returned by F(primme_svds) and
// copy the content of its returned Matrix<S> into the output vector y.
// Arguments:
// - T: type of PRIMME evecs
// - S: R type
// - TS: type of elements in Matrix<S>
// - F: F::get(primme_svds) return the function pointer
// - x, ldx, y, ldy, ...: arguments of matrixMatvec, applyPreconditioner...

template <typename T, int S, typename TS, typename F>
static void matrixMatvecSvds(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy,
      int *blockSize, int *mode, struct primme_svds_params *primme_svds,
      int *ierr)
{  
   // Get numbers of rows of x and y
   PRIMME_INT mx, my;
   const char *str;
   void *fp;
   F::get(*mode, primme_svds, &mx, &my, &fp, &str);

   // Create input vector
   Matrix<S> vx =
      createMatrix<TS,Matrix<S> >((TS*)x, mx, *blockSize,
            *ldx);

   // Call the callback
   Function *f = (Function*)fp;
   Matrix<S> vy = (*f)(vx, wrap(str));

   // Copy output vector
   copyMatrix(vy, (TS*)y, my, *blockSize, *ldy);

   *ierr = 0;
}

// Generic function for dprimme_svds and zprimme_svds
// Arguments:
// - T: type of PRIMME evecs
// - S: Matrix<S> is the type of the output eigenvectors
// - TS: type of elements in Matrix<S>
// - ortho: orthogonal constrains
// - init: initial guesses
// - A: matrix-vector product
// - B: mass matrix-vector product
// - prec: preconditioner application
// - convTest: convergence criterion

template<typename T, int S, typename TS>
static List xprimme_svds(Matrix<S> orthol, Matrix<S> orthor, Matrix<S> initl,
      Matrix<S> initr, Function A, SEXP prec, PrimmeSvdsParams primme_svds)
{
   if (primme_svds->mLocal == 0) primme_svds->mLocal = primme_svds->m;
   if (primme_svds->nLocal == 0) primme_svds->nLocal = primme_svds->n;

   // Check dimensions of ortho and init

   if (orthol.rows() != 0 && orthol.rows() != primme_svds->mLocal)
      stop("Invalid number of rows in input matrix orthol");
   if (orthor.rows() != 0 && orthor.rows() != primme_svds->nLocal)
      stop("Invalid number of rows in input matrix orthor");
   if (orthol.cols() != orthor.cols())
      stop("orthol and orthor should have the same number of columns");
   if (initl.rows() != 0 && initl.rows() != primme_svds->mLocal)
      stop("Invalid number of rows in input matrix initl");
   if (initr.rows() != 0 && initr.rows() != primme_svds->nLocal)
      stop("Invalid number of rows in input matrix initr");
   if (initl.cols() != initr.cols())
      stop("initl and initr should have the same number of columns");

   // Allocate svals, rnorms and svecs

   NumericVector vsvals(primme_svds->numSvals);
   NumericVector vrnorms(primme_svds->numSvals);
   double *svals = vsvals.begin(), *rnorms = vrnorms.begin();
   T *svecs =
      new T[(primme_svds->mLocal+primme_svds->nLocal)
      * (orthol.cols()+std::max(initl.cols(), (int)primme_svds->numSvals))];

   // Copy initial and orthogonal constrain vectors

   primme_svds->numOrthoConst = orthol.cols();
   primme_svds->initSize = initl.cols();
   copyMatrix(orthol, (TS*)svecs, primme_svds->mLocal, orthol.cols(), primme_svds->mLocal, false);
   copyMatrix(initl, (TS*)&svecs[primme_svds->mLocal*orthol.cols()], primme_svds->mLocal, initl.cols(), primme_svds->mLocal, false);
   int ncols = orthol.cols() + initl.cols();
   copyMatrix(orthor, (TS*)&svecs[primme_svds->mLocal*ncols], primme_svds->nLocal, orthor.cols(), primme_svds->nLocal, false);
   copyMatrix(initr, (TS*)&svecs[primme_svds->mLocal*ncols+primme_svds->nLocal*orthor.cols()], primme_svds->nLocal, initr.cols(), primme_svds->nLocal, false);

   // Set matvec and preconditioner

   primme_svds->matrix = &A;
   primme_svds->matrixMatvec = matrixMatvecSvds<T, S, TS, getSvdsForMatrix>;

   Function fprec = A; /* dummy set value */
   if (prec != R_NilValue) {
      fprec = as<Function>(prec);
      primme_svds->preconditioner = &fprec;
      primme_svds->applyPreconditioner = matrixMatvecSvds<T, S, TS, getSvdsForPreconditioner>;
   }

   // Call xprimme_svds
   int ret = tprimme_svds(svals, svecs, rnorms, primme_svds);

   // Return only the singular vectors (not the ortho)
   Matrix<S> svecsl(primme_svds->mLocal, primme_svds->initSize, (TS*)&svecs[primme_svds->mLocal*orthol.cols()]);
   Matrix<S> svecsr(primme_svds->nLocal, primme_svds->initSize, (TS*)&svecs[primme_svds->mLocal*(orthol.cols()+primme_svds->initSize)+primme_svds->nLocal*orthor.cols()]);

   // Return
   return List::create(
      Named("ret") = IntegerVector::create(ret),
      Named("svals") = vsvals,
      Named("svecsu") = svecsl,
      Named("svecsv") = svecsr,
      Named("rnorms") = vrnorms);
}

// [[Rcpp::export(.dprimme_svds)]]
List dprimme_svds_rcpp(NumericMatrix orthol, NumericMatrix orthor, NumericMatrix initl, NumericMatrix initr, Function A, SEXP prec, PrimmeSvdsParams primme_svds) {
   return xprimme_svds<double, REALSXP, double>(orthol, orthor, initl, initr, A, prec, primme_svds);
}

// [[Rcpp::export(.zprimme_svds)]]
List zprimme_svds_rcpp(ComplexMatrix orthol, ComplexMatrix orthor, ComplexMatrix initl, ComplexMatrix initr, Function A, SEXP prec, PrimmeSvdsParams primme_svds) {
   return xprimme_svds<PRIMME_COMPLEX_DOUBLE, CPLXSXP, Rcomplex>(orthol, orthor, initl, initr, A, prec, primme_svds);
}
