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
# File: eigs.R
# 
# Purpose - Driver to compute eigenvalues and eigenvectors.
# 
#*****************************************************************************

#' Find a few eigenvalues and vectors on large, sparse Hermitian matrix
#'
#' Compute a few eigenpairs from a specified region (the largest, the smallest,
#' the closest to a point) on a symmetric/Hermitian matrix using PRIMME [1].
#' Only the matrix-vector product of the matrix is required. The used method is
#' usually faster than a direct method (such as \code{\link{eigen}}) if
#' seeking a few eigenpairs and the matrix-vector product is cheap. For
#' accelerating the convergence consider to use preconditioning and/or educated
#' initial guesses.
#'
#' @param A symmetric/Hermitian matrix or a function with signature f(x) that
#'        returns \code{A \%*\% x}.
#' @param NEig number of eigenvalues and vectors to seek.
#' @param which which eigenvalues to find:
#'    \describe{
#'       \item{\code{"LA"}}{the largest (rightmost) values;}
#'       \item{\code{"SA"}}{the smallest (leftmost) values;}
#'       \item{\code{"LM"}}{the farthest from \code{targetShifts};}
#'       \item{\code{"SM"}}{the closest to \code{targetShifts};}
#'       \item{\code{"CLT"}}{the closest but left to \code{targetShifts};}
#'       \item{\code{"CGT"}}{the closest but greater than \code{targetShifts};}
#'       \item{vector of numbers}{the closest values to these points.}
#'    }
#' @param tol the convergence tolerance:
#'    \eqn{\|A x - x\lambda\| \le tol\|A\|}{||A*x - x*lambda|| <= tol*||A||}.
#' @param targetShifts return the closest eigenvalues to these points as
#'        indicated by \code{target}.
#' @param x0 matrix whose columns are educated guesses of the eigenvectors to
#'        to find.
#' @param ortho find eigenvectors orthogonal to the space spanned by the
#'        columns of this matrix; useful to avoid finding some eigenvalues or
#'        to find new solutions.
#' @param prec preconditioner used to accelerated the convergence; usually it
#'        is an approximation of the inverse of \eqn{A - \sigma I} if finding
#'        the closest eigenvalues to \eqn{\sigma}. If it is a matrix
#'        it is used as prec \%*\% x; otherwise it is used as prec(x).
#' @param isreal whether A \%*\% x always returns real number and not complex.
#' @param ... other PRIMME options (see details).
#' @return list with the next elements
#'    \describe{
#'       \item{\code{values}}{the eigenvalues \eqn{\lambda_i}}
#'       \item{\code{vectors}}{the eigenvectors \eqn{x_i}}
#'       \item{\code{rnorms}}{the residual vector norms
#'          \eqn{\|A x_i - \lambda_i x_i\|}{||A*x_i - lambda_i*x_i||}.}
#'       \item{\code{stats$numMatvecs}}{number of matrix-vector products performed}
#'       \item{\code{stats$numPreconds}}{number of preconditioner applications performed}
#'       \item{\code{stats$elapsedTime}}{time expended by the eigensolver}
#'       \item{\code{stats$timeMatvec}}{time expended in the matrix-vector products}
#'       \item{\code{stats$timePrecond}}{time expended in applying the preconditioner}
#'       \item{\code{stats$estimateMinEval}}{estimation of the smallest eigenvalue of A}
#'       \item{\code{stats$estimateMaxEval}}{estimation of the largest eigenvalue of A}
#'       \item{\code{stats$estimateANorm}}{estimation of the norm of A}
#'    }
#'
#' @details
#' Optional arguments to pass to PRIMME eigensolver (see further details at [2]):
#'
#' \describe{
#' \item{\code{method}}{ used by the solver, one of:
#'    \describe{
#'    \item{\code{"DYNAMIC"}}{                  switches dynamically between DEFAULT_MIN_TIME and DEFAULT_MIN_MATVECS}
#'    \item{\code{"DEFAULT_MIN_TIME"}}{         best method for light matrix-vector product}
#'    \item{\code{"DEFAULT_MIN_MATVECS"}}{      best method for heavy matrix-vector product or preconditioner}
#'    \item{\code{"Arnoldi"}}{                  an Arnoldi not implemented efficiently}
#'    \item{\code{"GD"}}{                       classical block Generalized Davidson }
#'    \item{\code{"GD_plusK"}}{                 GD+k block GD with recurrence restarting}
#'    \item{\code{"GD_Olsen_plusK"}}{           GD+k with approximate Olsen preconditioning}
#'    \item{\code{"JD_Olsen_plusK"}}{           GD+k, exact Olsen (two preconditioner applications per step)}
#'    \item{\code{"RQI"}}{                      Rayleigh Quotient Iteration, also Inverse Iteration
#'                                              if \code{targetShifts} is provided}
#'    \item{\code{"JDQR"}}{                     original block, Jacobi Davidson}
#'    \item{\code{"JDQMR"}}{                    our block JDQMR method (similar to JDCG)}
#'    \item{\code{"JDQMR_ETol"}}{               slight, but efficient JDQMR modification}
#'    \item{\code{"STEEPEST_DESCENT"}}{         equivalent to GD(\code{maxBlockSize},2*\code{maxBlockSize})}
#'    \item{\code{"LOBPCG_OrthoBasis"}}{        equivalent to GD(\code{neig},3*\code{neig})+\code{neig}}
#'    \item{\code{"LOBPCG_OrthoBasis_Window"}}{ equivalent to GD(\code{maxBlockSize},3*\code{maxBlockSize})+\code{maxBlockSize} when neig>\code{maxBlockSize}}
#'    }}
#'    \item{\code{aNorm}}{estimation of norm-2 of A, used in convergence test (if not
#'        provided, it is estimated as the largest eigenvalue in magnitude
#'        seen).}
#'    \item{\code{maxBlockSize}}{maximum block size (like in subspace iteration or
#'        LOBPCG).}
#'    \item{\code{printLevel}}{message level reporting, from 0 (no output) to 5 (show all).} 
#'    \item{\code{locking}}{1, hard locking; 0, soft locking.}
#'    \item{\code{maxBasisSize}}{maximum size of the search subspace.}
#'    \item{\code{minRestartSize}}{ minimum Ritz vectors to keep in restarting.}
#'    \item{\code{maxMatvecs}}{ maximum number of matrix vector multiplications.}
#'    \item{\code{maxit}}{ maximum number of outer iterations.}
#'    \item{\code{scheme}}{ the restart scheme (thick restart by default).}
#'    \item{\code{maxPrevRetain}}{ number of approximate eigenvectors retained from
#'          previous iteration, that are kept after restart.}
#'    \item{\code{robustShifts}}{ set to true to avoid stagnation.}
#'    \item{\code{maxInnerIterations}}{ maximum number of inner QMR iterations.}
#'    \item{\code{LeftQ}}{ use the locked vectors in the left projector.}
#'    \item{\code{LeftX}}{ use the approx. eigenvector in the left projector.}
#'    \item{\code{RightQ}}{ use the locked vectors in the right projector.}
#'    \item{\code{RightX}}{ use the approx. eigenvector in the right projector.}
#'    \item{\code{SkewQ}}{ use the preconditioned locked vectors in the right projector.}
#'    \item{\code{SkewX}}{ use the preconditioned approximate eigenvector in the right
#'                projector.}
#'    \item{\code{relTolBase}}{ a legacy from classical JDQR (recommend not use).}
#'    \item{\code{iseed}}{ an array of four numbers used as a random seed.}
#' }
#'
#' @references
#' [1] A. Stathopoulos and J. R. McCombs \emph{PRIMME: PReconditioned Iterative
#'     MultiMethod Eigensolver: Methods and software description}, ACM
#'     Transaction on Mathematical Software Vol. 37, No. 2, (2010)
#'     21:1-21:30.
#'
#' [2] \url{http://www.cs.wm.edu/~andreas/software/doc/primmec.html#parameters-guide}
#'
#' @seealso
#' \code{\link{eigen}} for computing all values;
#' \code{\link{svds}} for computing a few singular values
#'
#' @examples
#' A <- diag(1:10)  # the eigenvalues of this matrix are 1:10 and the
#'                  # eigenvectors are the columns of diag(10)
#' r <- eigs_sym(A, 3);
#' r$values  # the three largest eigenvalues on diag(1:10)
#' r$vectors # the corresponding approximate eigenvectors
#' r$rnorms # the corresponding residual norms
#' r$stats$numMatvecs # total matrix-vector products spend
#'
#' r <- eigs_sym(A, 3, 'SA') # compute the three smallest values
#'
#' r <- eigs_sym(A, 3, 2.5) # compute the three closest values to 2.5
#'
#' r <- eigs_sym(A, 3, 2.5, tol=1e-3); # compute the values with
#' r$rnorms                                    # residual norm <= 1e-3*||A||
#'
#' # Build a Jacobi preconditioner (too convenient for a diagonal matrix!)
#' # and see how reduce the number matrix-vector products
#' A <- diag(1:1000)   # we use a larger matrix to amplify the difference
#' P <- diag(diag(A) - 2.5)
#' eigs_sym(A, 3, 2.5, tol=1e-3)$stats$numMatvecs
#' eigs_sym(A, 3, 2.5, tol=1e-3, prec=P)$stats$numMatvecs
#' 
#' # Passing A and the preconditioner as functions
#' Af <- function(x) (1:100) * x; # = diag(1:100) %*% x
#' Pf <- function(x) x / (1:100 - 2.5); # = solve(diag(1:100 - 2.5), x)
#' r <- eigs_sym(Af, 3, 2.5, tol=1e-3, prec=Pf, n=100)
#'
#' # Passing initial guesses
#' A <- diag(1:1000)   # we use a larger matrix to amplify the difference
#' x0 <- diag(1,1000,4) + matrix(rnorm(4000), 1000, 4)/100;
#' eigs_sym(A, 4, "SA", tol=1e-3)$stats$numMatvecs
#' eigs_sym(A, 4, "SA", tol=1e-3, x0=x0)$stats$numMatvecs
#' 
#' # Passing orthogonal constrain, in this case, already compute eigenvectors
#' r <- eigs_sym(A, 4, "SA", tol=1e-3); r$values
#' eigs_sym(A, 4, "SA", tol=1e-3, ortho=r$vectors)$values
#' 
#' @useDynLib PRIMME
#' @importFrom Rcpp evalCpp
#' @export

eigs_sym <- function(A, NEig=1, which="LA", targetShifts=NULL, tol=1e-6,
      x0=NULL, ortho=NULL, prec=NULL, isreal=NULL, ...) {

   # Extra arguments are considered PRIMME options
   opts <- list(...);

   # If A is a function, check that n is defined. Otherwise check that
   # A is a square matrix and get its dimension
   if (is.function(A)) {
      if (!.is.wholenumber(opts$n))
         stop("matrix dimension not set (set 'n')");
      Af <- A;
      isreal_suggestion <- FALSE;
   }
   else if (length(dim(A)) != 2 || ncol(A) != nrow(A)) {
      stop("A should be a square matrix or a function")
   }
   else {
      opts$n <- nrow(A);
      isreal_suggestion <-
         if (is.matrix(A)) is.numeric(A)
         else (inherits(A, "Matrix") && substr(class(A), 0, 1) == "d");
      if ((is.null(isreal) || isreal == isreal_suggestion) && (
               is.matrix(A) ||
               any(c("dmatrix", "dgeMatrix", "dgCMatrix", "dsCMatrix") %in% class(A)) ||
               any(c("zmatrix", "zgeMatrix", "zgCMatrix", "zsCMatrix") %in% class(A)) )) {
         Af <- A;
      }
      else {
         Af <- function(x) A %*% x;
      }
   }

   # Check NEig and set the option
   if (!.is.wholenumber(NEig) || NEig > opts$n) {
      stop("NEig should be an integer not greater than the matrix dimension");
   } else {
      opts$numEvals <- NEig
   }

   # Check target at set the option
   targets = list(LA="primme_largest",
         LM="primme_largest_abs",
         SA="primme_smallest",
         CGT="primme_closest_geq",
         CLT="primme_closest_leq",
         SM="primme_closest_abs");
   if (is.character(which) && which %in% names(targets)) {
      if (which %in% list("CGT", "CLT", "SM") && is.null(targetShifts))
         targetShifts <- 0.0;
      opts$target <- targets[[which]];
      opts$targetShifts <- targetShifts;
   }
   else if (is.numeric(which)) {
      opts$targetShifts <- which;
      opts$target <- "primme_closest_abs";
   }
   else {
      stop("target should be numeric or LA, LM, SA, CGT, CLT or SM");
   }
 
   # Check tol and set the option
   if (is.function(tol)) {
      convTest <- tol
   }
   else if (!is.numeric(tol) || tol < 0 || tol >= 1) {
      stop("tol should be a positive number smaller than 1 or a function")
   }
   else {
      opts$eps <- tol;
      convTest <- NULL;
   }

   # Check x0 is a matrix of proper dimensions
   if (is.null(x0))
      x0 <- matrix(nrow=0, ncol=0)
   else if (!is.matrix(x0) || nrow(x0) != opts$n)
      stop("x0 should be NULL or a matrix with the same number of rows as A");

   # Check ortho is a matrix of proper dimensions
   if (is.null(ortho))
      ortho <- matrix(nrow=0, ncol=0)
   else if (!is.matrix(ortho) || nrow(ortho) != opts$n)
      stop("ortho should be NULL or a matrix with the same number of rows as A");

   # Check that prec is a function or a square matrix with proper dimensions
   if (is.null(prec) || is.function(prec)) {
      precf <- prec
   }
   else if (!is.matrix(prec) || length(dim(prec)) != 2 || ncol(prec) != opts$n || nrow(prec) != opts$n) {
      stop("prec should be a square matrix with the same dimension as A or a function")
   }
   else {
      precf <- function(x) solve(prec, x)
   }

   # Extract method from opts
   method <- opts$method;
   opts$method <- NULL;

   # Process isreal
   if (!is.null(isreal) && !is.logical(isreal)) {
      stop("isreal should be logical");
   }
   else if (is.null(isreal)) {
      isreal <- isreal_suggestion;
   }

   # Initialize PRIMME
   primme <- .primme_initialize();

   # Set options
   for (x in names(opts))
      .primme_set_member(x, opts[[x]], primme);

   # Set method
   if (!is.null(method))
      .primme_set_method(method, primme);

   # Call PRIMME
   r <- if (!isreal)
      .zprimme(ortho, x0, Af, NULL, precf, convTest, primme)
   else
      .dprimme(ortho, x0, Af, NULL, precf, convTest, primme);

   # Get stats
   r$stats$numMatvecs <- .primme_get_member("stats_numMatvecs", primme)
   r$stats$numPreconds <- .primme_get_member("stats_numPreconds", primme)
   r$stats$elapsedTime <- .primme_get_member("stats_elapsedTime", primme)
   r$stats$estimateMinEval <- .primme_get_member("stats_estimateMinEVal", primme)
   r$stats$estimateMaxEval <- .primme_get_member("stats_estimateMaxEVal", primme)
   r$stats$estimateANorm <- .primme_get_member("stats_estimateLargestSVal", primme)
   r$stats$timeMatvec <- .primme_get_member("stats_timeMatvec", primme)
   r$stats$timePrecond <- .primme_get_member("stats_timePrecond", primme)
   
   # Free PRIMME structure
   .primme_free(primme);

   # Return values, vectors, residuals norms and stats if no error;
   # stop otherwise
   if (r$ret != 0)
      stop(.getEigsErrorMsg(r$ret))
   else
      r$ret <- NULL;
   r
}

.getEigsErrorMsg <- function(n) {
   l <- list(
      "0"= "success",
      "1"= "reported only amount of required memory",
      "-1"= "failed in allocating int or real workspace",
      "-2"= "malloc failed in allocating a permutation integer array",
      "-3"= "main_iter() encountered problem; the calling stack of the functions where the error occurred was printed in 'stderr'",
      "-4"= "argument 'primme' is NULL",
      "-5"= "'n' < 0 or 'nLocal' < 0 or 'nLocal' > 'n'",
      "-6"= "'numProcs' < 1",
      "-7"= "'matrixMatvec' is NULL",
      "-8"= "'applyPreconditioner' is NULL and 'precondition' is not NULL",
      "-9"= "'not used",
      "-10"= "'numEvals' > 'n'",
      "-11"= "'numEvals' < 0",
      "-12"= "'eps' > 0 and 'eps' < machine precision",
      "-13"= "'target' is not properly defined",
      "-14"= "'target' is one of 'primme_largest_abs', 'primme_closest_geq', 'primme_closest_leq' or 'primme_closest_abs' but 'numTargetShifts' <= 0 (no shifts)",
      "-15"= "'target' is one of 'primme_largest_abs', 'primme_closest_geq', 'primme_closest_leq' or 'primme_closest_abs' but 'targetShifts' is NULL  (no shifts array)",
      "-16"= "'numOrthoConst' < 0 or 'numOrthoConst' > 'n'. (no free dimensions left)",
      "-17"= "'maxBasisSize' < 2",
      "-18"= "'minRestartSize' < 0 or 'minRestartSize' shouldn't be zero",
      "-19"= "'maxBlockSize' < 0 or 'maxBlockSize' shouldn't be zero",
      "-20"= "'maxPrevRetain' < 0",
      "-21"= "'scheme' is not one of *primme_thick* or *primme_dtr*",
      "-22"= "'initSize' < 0",
      "-23"= "'locking' == 0 and 'initSize' > 'maxBasisSize'",
      "-24"= "'locking' and 'initSize' > 'numEvals'",
      "-25"= "'maxPrevRetain' + 'minRestartSize' >= 'maxBasisSize'",
      "-26"= "'minRestartSize' >= 'n'",
      "-27"= "'printLevel' < 0 or 'printLevel' > 5",
      "-28"= "'convTest' is not one of 'primme_full_LTolerance', 'primme_decreasing_LTolerance', 'primme_adaptive_ETolerance' or 'primme_adaptive'",
      "-29"= "'convTest' == 'primme_decreasing_LTolerance' and 'relTolBase' <= 1",
      "-30"= "'evals' is NULL, but not 'evecs' and 'resNorms'",
      "-31"= "'evecs' is NULL, but not 'evals' and 'resNorms'",
      "-32"= "'resNorms' is NULL, but not 'evecs' and 'evals'",
      "-33"= "'locking' == 0 and 'minRestartSize' < 'numEvals'",
      "-34"= "'ldevecs' is less than 'nLocal'",
      "-35"= "'ldOPs' is non-zero and less than 'nLocal'",
      "-36"= "not enough memory for realWork",
      "-37"= "not enough memory for intWork",
      "-38"= "'locking' == 0 and 'target' is 'primme_closest_leq' or 'primme_closet_geq'");
   l[[as.character(n)]];
}

#' @keywords internal
.is.wholenumber <-
         function(x, tol = .Machine$double.eps^0.5)
            is.integer(x) || (is.numeric(x) && abs(x - round(x)) < tol)
