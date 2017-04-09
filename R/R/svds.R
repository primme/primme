#******************************************************************************
# Copyright (c) 2016, College of William & Mary
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the College of William & Mary nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# PRIMME: https://github.com/primme/primme
# Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u
#******************************************************************************
# File: svds.R
# 
# Purpose - Driver to compute singular values and vectors.
# 
#*****************************************************************************

#' Find a few singular values and vectors on large, sparse matrix
#'
#' Compute a few singular triplets from a specified region (the largest, the
#' smallest, the closest to a point) on a matrix using PRIMME [1].
#' Only the matrix-vector product of the matrix is required. The used method is
#' usually faster than a direct method (such as \code{\link{svd}}) if
#' seeking few singular values and the matrix-vector product is cheap. For
#' accelerating the convergence consider to use preconditioning  and/or
#' educated initial guesses.
#'
#' @param A matrix or a function with signature f(x, trans) that returns
#'        \code{A \%*\% x} when \code{trans == "n"} and
#'        \code{t(Conj(A)) \%*\% x} when \code{trans == "c"}.
#' @param NSvals number of singular triplets to seek.
#' @param which which singular values to find:
#'    \describe{
#'       \item{\code{"L"}}{the largest values;}
#'       \item{\code{"S"}}{the smallest values;}
#'       \item{vector of numbers}{the closest values to these points.}
#'    }
#' @param tol a triplet \eqn{(\sigma,u,v)} is marked as converged when
#'    \eqn{\sqrt{\|Av - \sigma u\|^2+\|A^*u - \sigma v\|^2} \le tol\|A\|}{sqrt(||A*v - sigma*u||^2 + ||A'*u - \sigma*v||^2}
#'    is smaller than \eqn{tol*||A||}, or close to the minimum tolerance that
#'    the selected method can achieve.
#' @param u0 matrix whose columns are educated guesses of the left singular
#'        vectors to find.
#' @param v0 matrix whose columns are educated guesses of the right singular
#'        vectors to find.
#' @param orthou find left singular vectors orthogonal to the space spanned by
#'        the columns of this matrix; useful to avoid finding some triplets or
#'        to find new solutions.
#' @param orthov find right singular vectors orthogonal to the space spanned by
#'        the columns of this matrix.
#' @param prec preconditioner used to accelerated the convergence; it is a named
#'        list of matrices or functions such as \code{solve(prec[[mode]],x)} or
#'        \code{prec[[mode]](x)} return an approximation of \eqn{OP^{-1} x},
#'        where
#'        \tabular{cc}{
#'          \code{mode}  \tab \eqn{OP} \cr
#'          \code{"AHA"} \tab \eqn{A^*A} \cr
#'          \code{"AAH"} \tab \eqn{A A^*} \cr
#'          \code{"aug"} \tab \eqn{[0 A; A^* 0]}
#'        }
#'        The three values haven't to be set. It is recommended to set
#'        \code{"AHA"} for matrices with nrow > ncol; \code{"AAH"} for
#'        matrices with nrow < ncol; and additionally \code{"aug"} for
#'        \code{tol} < 1e-8.
#' @param isreal whether A \%*\% x always returns real number and not complex.
#' @param ... other PRIMME options (see details).
#' @return list with the next elements
#'    \describe{
#'       \item{\code{d}}{the singular values \eqn{\sigma_i}}
#'       \item{\code{u}}{the left singular vectors \eqn{u_i}}
#'       \item{\code{v}}{the right singular vectors \eqn{v_i}}
#'       \item{\code{rnorms}}{the residual vector norms
#'          \eqn{\sqrt{\|Av - \sigma u\|^2+\|A^*u - \sigma v\|^2}}{sqrt(||A*v - sigma*u||^2 + ||A'*u - \sigma*v||^2}}
#'       \item{\code{stats$numMatvecs}}{matrix-vector products performed}
#'       \item{\code{stats$numPreconds}}{number of preconditioner applications performed}
#'       \item{\code{stats$elapsedTime}}{time expended by the eigensolver}
#'       \item{\code{stats$timeMatvec}}{time expended in the matrix-vector products}
#'       \item{\code{stats$timePrecond}}{time expended in applying the preconditioner}
#'       \item{\code{stats$estimateANorm}}{estimation of the norm of A}
#'    }
#'
#' @details
#' Optional arguments to pass to PRIMME eignesolver (see further details at
#' [2]):
#'
#' \describe{
#'    \item{\code{aNorm}}{estimation of norm-2 of A, used in convergence test
#'       (if not provided, it is estimated as the largest eigenvalue in 
#'       magnitude seen)}
#'    \item{\code{maxBlockSize}}{maximum block size (like in subspace iteration
#'       or LOBPCG)}
#'    \item{\code{printLevel}}{message level reporting, from 0 (no output) to 5
#'       (show all)} 
#'    \item{\code{locking}}{1, hard locking; 0, soft locking}
#'    \item{\code{maxBasisSize}}{maximum size of the search subspace}
#'    \item{\code{minRestartSize}}{ minimum Ritz vectors to keep in restarting}
#'    \item{\code{maxMatvecs}}{ maximum number of matrix vector multiplications}
#'    \item{\code{iseed}}{ an array of four numbers used as a random seed}
#'    \item{\code{method}}{which equivalent eigenproblem to solve
#'       \describe{
#'          \item{\code{"primme_svds_normalequation"}}{\eqn{A^*A} or \eqn{AA^*}}
#'          \item{\code{"primme_svds_augmented"}}{ \eqn{[0 A^*;A 0]}}
#'          \item{\code{"primme_svds_hybrid"}}{ first normal equations and
#'                      then augmented (default)}
#'       }                   
#'    }
#'    \item{\code{locking}}{1, hard locking; 0, soft locking}
#'    \item{\code{primmeStage1, primmeStage2}}{list with options for the first
#'       and the second stage solver; see \code{\link{eigs_sym}}}
#' }
#'
#' If \code{method} is \code{"primme_svds_normalequation"}, the minimum
#' tolerance that can be achieved is \eqn{\|A\|\epsilon/\sigma}, where \eqn{\epsilon}
#' is the machine precision. If \code{method} is \code{"primme_svds_augmented"}
#' or \code{"primme_svds_hybrid"}, the minimum tolerance is \eqn{\|A\|\epsilon}.
#' However it may not return triplets with singular values smaller than
#' \eqn{\|A\|\epsilon}.
#'
#' @references
#' [1]  L. Wu, E. Romero and A. Stathopoulos, \emph{PRIMME_SVDS: A High-Performance
#'      Preconditioned SVD Solver for Accurate Large-Scale Computations},
#'      arXiv:1607.01404
#'
#' [2] \url{http://www.cs.wm.edu/~andreas/software/doc/svdsc.html#parameters-guide}
#'
#' @seealso
#' \code{\link{svd}} for computing all singular triplets;
#' \code{\link{eigs_sym}} for computing a few eigenvalues and vectors
#'    from a symmetric/Hermitian matrix.
#'
#' @examples
#' A <- diag(1:5,10,5)  # the singular values of this matrix are 1:10 and the
#'                         # left and right singular vectors are the columns of
#'                         # diag(1,100,10) and diag(10), respectively
#' r <- svds(A, 3);
#' r$d # the three largest singular values on A
#' r$u # the corresponding approximate left singular vectors
#' r$v # the corresponding approximate right singular vectors
#' r$rnorms # the corresponding residual norms
#' r$stats$numMatvecs # total matrix-vector products spend
#'
#' r <- svds(A, 3, "S") # compute the three smallest values
#'
#' r <- svds(A, 3, 2.5) # compute the three closest values to 2.5
#'
#' A <- diag(1:500,500,100)   # we use a larger matrix to amplify the difference
#' r <- svds(A, 3, 2.5, tol=1e-3); # compute the values with 
#' r$rnorms                               # residual norm <= 1e-3*||A||
#'
#' # Build the diagonal squared preconditioner
#' # and see how reduce the number matrix-vector products
#' P <- diag(colSums(A^2))
#' svds(A, 3, "S", tol=1e-3)$stats$numMatvecs
#' svds(A, 3, "S", tol=1e-3, prec=list(AHA=P))$stats$numMatvecs
#' 
#' # Passing A and the preconditioner as functions
#' Af <- function(x,mode) if (mode == "n") A%*%x else crossprod(A,x);
#' P = colSums(A^2);
#' PAHAf <- function(x) x / P;
#' r <- svds(Af, 3, "S", tol=1e-3, prec=list(AHA=PAHAf), m=500, n=100)
#'
#' # Passing initial guesses
#' v0 <- diag(1,100,4) + matrix(rnorm(400), 100, 4)/100;
#' svds(A, 4, "S", tol=1e-3)$stats$numMatvecs
#' svds(A, 4, "S", tol=1e-3, v0=v0)$stats$numMatvecs
#' 
#' # Passing orthogonal constrain, in this case, already compute singular vectors
#' r <- svds(A, 4, "S", tol=1e-3); r$d
#' svds(A, 4, "S", tol=1e-3, orthov=r$v)$d
#'
#' @useDynLib PRIMME
#' @importFrom Rcpp evalCpp
#' @export

svds <- function(A, NSvals, which="L", tol=1e-6, u0=NULL, v0=NULL,
      orthou=NULL, orthov=NULL, prec=NULL, isreal=NULL, ...) {

   # Extra arguments are considered PRIMME options
   opts <- list(...);

   # If A is a function, check that n is defined
   if (is.function(A)) {
      if (!.is.wholenumber(opts$n) || !.is.wholenumber(opts$m))
         stop("matrix dimension not set (set 'm' and 'n')");
      Af <- A;
      isreal_suggestion <- FALSE;
   }
   else if (length(dim(A)) != 2) {
      stop("A should be a matrix or a function")
   }
   else {
      opts$m <- nrow(A);
      opts$n <- ncol(A);
      if (is.matrix(A)) {
         Af <- A;
         isreal_suggestion <- is.numeric(A);
      }
      else if (any(c("dmatrix", "dgeMatrix", "dgCMatrix", "dsCMatrix") %in% class(A))) {
        Af <- A;
        isreal_suggestion <- TRUE;
      }
      else if (any(c("zmatrix", "zgeMatrix", "zgCMatrix", "zsCMatrix") %in% class(A))) {
        Af <- A;
        isreal_suggestion <- FALSE;
      }
      else {
         Af <- function(x,trans)
            if (trans == "n") A %*% x else Conj(t(crossprod(Conj(x),A)));
         isreal_suggestion <- FALSE;
      }
   }

   # Check nsvals and set the option
   if (!.is.wholenumber(NSvals) || NSvals > min(opts$m, opts$n))
      stop("NSvals should be an integer not greater than the smallest dimension of the matrix");
   opts$numSvals <- NSvals

   # Check target at set the option
   targets = list(L="primme_svds_largest",
                  S="primme_svds_smallest");
   if (is.character(which) && which %in% names(targets)) {
      opts$target <- targets[[which]];
   }
   else if (is.numeric(which)) {
      opts$targetShifts <- which;
      opts$target <- "primme_svds_closest_abs";
   }
   else {
      stop("target should be numeric or L or S");
   }
 
   # Check tol and set the option
   if (!is.numeric(tol) || tol < 0 || tol >= 1) {
      stop("tol should be a positive number smaller than 1 or a function")
   }
   else {
      opts$eps <- tol;
   }

   # Check u0,v0 and orthou,orthov is a matrix of proper dimensions
   check_uv <- function(u, v, u_name, v_name) {
      if (!is.null(u) && (!is.matrix(u) || nrow(u) != opts$m)) {
         stop(paste(u_name, "should be NULL or a matrix with the same number of rows as A"))
      }
      else if (!is.null(u) && is.null(v)) {
         v <- Af(u, "c");
      }

      if (!is.null(v) && (!is.matrix(v) || nrow(v) != opts$n)) {
         stop(paste(v_name, "should be NULL or a matrix with the same number of rows as columns A has"))
      }
      else if (!is.null(v) && is.null(u)) {
         u <- Af(v, "n");
      }

      if (is.null(u) && is.null(v))
         list(u=matrix(nrow=0, ncol=0), v=matrix(nrow=0, ncol=0))
      else 
         list(u=u, v=v);
   }

   ortho <- check_uv(orthou, orthov, "orthou", "orthov");
   init <- check_uv(u0, v0, "u0", "v0");

   # Check that prec is a function or a square matrix with proper dimensions
   if (is.null(prec) || is.function(prec)) {
      precf <- prec
   }
   else if (!is.list(prec)) {
      stop("prec should be a function or a list(AHA=...,AAH=...,aug=...)")
   }
   else if (!is.null(prec$AHA) && !is.function(prec$AHA)
         && (!is.matrix(prec$AHA) || length(dim(prec$AHA)) != 2
            || ncol(prec$AHA) != opts$n || nrow(prec$AHA) != opts$n)) {
      stop("prec$AHA should be a function or a square matrix of dimension ncol(A)")
   }
   else if (!is.null(prec$AAH) && !is.function(prec$AAH)
         && (!is.matrix(prec$AAH) || length(dim(prec$AAH)) != 2
            || ncol(prec$AAH) != opts$m || nrow(prec$AAH) != opts$m)) {
      stop("prec$AAH should be a function or a square matrix of dimension nrow(A)")
   }
   else if (!is.null(prec$aug) && !is.function(prec$aug)
         && (!is.matrix(prec$aug) || length(dim(prec$aug)) != 2
            || ncol(prec$aug) != opts$m+opts$n || nrow(prec$aug) != opts$m+opts$n)) {
      stop("prec$aug should be a function or a square matrix of dimension nrow(A)+ncol(A)")
   }
   else {
      tofunc <- function(x)
         if (is.function(x)) x
         else if (is.null(x)) identity
         else function(v) solve(x, v);
     precf <- function(x, mode) tofunc(prec[[mode]])(x);
   }

   # Extract method* from opts
   methodStage1 <- opts[["methodStage1"]];
   opts$methodStage1 <- NULL;
   methodStage2 <- opts[["methodStage2"]];
   opts$methodStage2 <- NULL;
   method <- opts[["method"]];
   opts$method <- NULL;

   # Process isreal
   if (!is.null(isreal) && !is.logical(isreal)) {
      stop("isreal should be logical");
   }
   else if (is.null(isreal)) {
      isreal <- isreal_suggestion;
   }

   # Initialize PRIMME SVDS
   primme_svds <- .primme_svds_initialize();

   # Set options
   for (x in names(opts)) {
      if (is.list(opts[[x]])) {
         primme <- .primme_svds_get_member(x, primme_svds);
         for (primmex in names(opts[[x]]))
            .primme_set_member(primmex, opts[[x]][[primmex]], primme);
      }
      else {
         .primme_svds_set_member(x, opts[[x]], primme_svds);
      }
   }

   # Set method
   if (!is.null(method) || !is.null(methodStage1) || !is.null(methodStage2)) {
      if (is.null(method)) method <- "primme_svds_default";
      if (is.null(methodStage1)) methodStage1 <- "PRIMME_DEFAULT_METHOD";
      if (is.null(methodStage2)) methodStage2 <- "PRIMME_DEFAULT_METHOD";
      .primme_svds_set_method(method, methodStage1, methodStage2, primme_svds);
   }

   # Call PRIMME SVDS
   r <- if (!isreal)
      .zprimme_svds(ortho$u, ortho$v, init$u, init$v, Af, precf, primme_svds)
   else
      .dprimme_svds(ortho$u, ortho$v, init$u, init$v, Af, precf, primme_svds)

   # Get stats
   r$stats$numMatvecs <- .primme_svds_get_member("stats_numMatvecs", primme_svds)
   r$stats$numPreconds <- .primme_svds_get_member("stats_numPreconds", primme_svds)
   r$stats$elapsedTime <- .primme_svds_get_member("stats_elapsedTime", primme_svds)
   r$stats$estimateANorm <- .primme_svds_get_member("aNorm", primme_svds)
   r$stats$timeMatvec <- .primme_svds_get_member("stats_timeMatvec", primme_svds)
   r$stats$timePrecond <- .primme_svds_get_member("stats_timePrecond", primme_svds)
   
   # Free PRIMME SVDS structure
   .primme_svds_free(primme_svds);

   # Return values, vectors, residuals norms and stats if no error;
   # stop otherwise
   if (r$ret != 0)
      stop(.getSvdsErrorMsg(r$ret))
   else
      r$ret <- NULL;
   r
}


.getSvdsErrorMsg <- function(n) {
   l <- list(
      "0" = "success",
      "1" = "reported only amount of required memory",
      "-1" = "failed in allocating int or real workspace",
      "-2" = "malloc failed in allocating a permutation integer array",
      "-3" = "main_iter() encountered problem; the calling stack of the functions where the error occurred was printed in 'stderr'",
      "-4" = "primme_svds is NULL",
      "-5" = "Wrong value for m or n or mLocal or nLocal",
      "-6" = "Wrong value for numProcs",
      "-7" = "matrixMatvec is not set",
      "-8" = "applyPreconditioner is not set but precondition == 1 ",
      "-9" = "numProcs >1 but globalSumDouble is not set",
      "-10" = "Wrong value for numSvals, it's larger than min(m, n)",
      "-11" = "Wrong value for numSvals, it's smaller than 1",
      "-13" = "Wrong value for target",
      "-14" = "Wrong value for method",
      "-15" = "Not supported combination of method and methodStage2",
      "-16" = "Wrong value for printLevel",
      "-17" = "svals is not set",
      "-18" = "svecs is not set",
      "-19" = "resNorms is not set",
      "-20" = "not enough memory for realWork",
      "-21" = "not enough memory for intWork");
   if (n >= -100)
      l[[as.character(n)]]
   else if (n >= -200)
      paste("Error from PRIMME first stage:", .getEigsErrorMsg(n+100))
   else if (n >= -300)
      paste("Error from PRIMME second stage:", .getEigsErrorMsg(n+200))
   else
      "Unknown error code";
}
