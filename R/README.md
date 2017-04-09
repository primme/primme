<!-- README.md is generated from README.Rmd. Please edit that file -->
PRIMME
======

This package is an R interface to PRIMME, a C library for computing a few eigenvalues and their corresponding eigenvectors of a real symmetric or complex Hermitian matrix. It can also compute singular values and vectors of a square or rectangular matrix. It can find largest, smallest, or interior singular/eigenvalues and can use preconditioning to accelerate convergence. It is especially optimized for large, difficult problems, and can be a useful tool for both non-experts and experts.

Use the following two references to cite this package:

-   A. Stathopoulos and J. R. McCombs *PRIMME: PReconditioned Iterative MultiMethod Eigensolver: Methods and software description*, ACM Transaction on Mathematical Software Vol. 37, No. 2, (2010), 21:1-21:30.

-   L. Wu, E. Romero and A. Stathopoulos, *PRIMME\_SVDS: A High-Performance Preconditioned SVD Solver for Accurate Large-Scale Computations*, arXiv:1607.01404

Installation Instructions
=========================

We are currently working to put the PRIMME package on CRAN. Meanwhile, to install the latest version:

``` r
library(devtools)
install_github("primme/primme", subdir="R")
```

Usage
=====

Load the package as usual:

``` r
library(PRIMME)
```

Eigenvalue problems
-------------------

The next example computes the three largest eigenvalues of the matrix `A`, which in this case is a dense diagonal matrix. It shows all the eigenvalues `values`, the eigenvectors `vectors`, the residual norms `rnorms` and some stats, such as the time `stats$elapsedTime` and the number of matrix vector multiplications performed `stats$numMatvecs`:

``` r
A <- diag(1:10) 
r <- eigs_sym(A, 3);
r
#> $values
#> [1] 10  9  8
#> 
#> $vectors
#>                [,1]          [,2]          [,3]
#>  [1,] -1.371868e-16  2.381771e-16 -2.252380e-16
#>  [2,]  6.980780e-17  2.866614e-17  1.292433e-16
#>  [3,] -2.299606e-16  1.269985e-16 -5.609795e-17
#>  [4,] -1.960802e-16  2.701925e-17 -2.503660e-17
#>  [5,] -4.857749e-17  2.462784e-16 -1.267024e-16
#>  [6,] -2.577747e-16 -2.079437e-16 -1.676989e-16
#>  [7,] -8.486805e-17 -7.529381e-16 -2.007853e-15
#>  [8,] -1.120694e-15 -2.065498e-15 -1.000000e+00
#>  [9,] -6.414024e-16 -1.000000e+00  1.560476e-15
#> [10,]  1.000000e+00 -2.987345e-16 -1.484140e-15
#> 
#> $rnorms
#> [1] 5.155287e-15 4.440892e-15 4.740049e-15
#> 
#> $stats
#> $stats$numMatvecs
#> [1] 10
#> 
#> $stats$numPreconds
#> [1] 0
#> 
#> $stats$elapsedTime
#> [1] 0.0006670952
#> 
#> $stats$estimateMinEval
#> [1] 1
#> 
#> $stats$estimateMaxEval
#> [1] 10
#> 
#> $stats$estimateANorm
#> [1] 10
#> 
#> $stats$timeMatvec
#> [1] 0.0004501343
#> 
#> $stats$timePrecond
#> [1] 0
```

The next examples show how to compute eigenvalues in other parts of the spectrum:

``` r
A <- diag(1:10)

r <- eigs_sym(A, 3, 'SA'); # compute the three smallest values
r$values
#> [1] 1 2 3

r <- eigs_sym(A, 3, 5.1); # compute the three closest values to 5.1
r$values
#> [1] 5 6 4
```

In some cases, a larger convergence tolerance may suffice:

``` r
A <- diag(1:5000)

r <- eigs_sym(A, 10, 'SA');
r$stats$numMatvecs
#> [1] 1146

r <- eigs_sym(A, 10, 'SA', tol=1e-3); 
r$stats$numMatvecs
#> [1] 409
```

Preconditioners, if available can reduce the time/matrix-vector multiplications significantly (see TODO):

``` r
# A is a tridiagonal
A <- diag(1:5000)
for(i in 1:4999) {A[i,i+1]<-1; A[i+1,i]<-1}

r <- eigs_sym(A, 10, 'SA');
r$stats$numMatvecs
#> [1] 1179
r$stats$elapsedTime
#> [1] 5.297332

# Jacobi preconditioner
P = diag(A);
r <- eigs_sym(A, 10, 'SA', prec=function(x)x/P);
r$stats$numMatvecs
#> [1] 51
r$stats$elapsedTime
#> [1] 0.2373209
```

Dense matrices, sparse matrices, and functions that return the matrix-vector product can be passed as the matrix problem `A`:

``` r
r <- eigs_sym(diag(1:10), 1); # dense matrix
library(Matrix)
r <- eigs_sym(Matrix(diag(1:10), sparse=TRUE), 1); # sparse matrix
Afun = function(x) matrix(1:10)*x;  # function that does diag(1:10) %*% x
r <- eigs_sym(Afun, 1, n=10); # n is the matrix dimension corresponding to Afun
```

The next benchmark function extends `rbenchmark` to return besides the time, the number of matrix-vector multiplications and the maximum residual norm among all returned eigenpairs.

``` r
library(knitr)

bench_eigs <- function(..., A, environment=parent.frame()) {
   arguments = match.call()[-1]
   if (!is.null(names(arguments)))
      arguments = arguments[!names(arguments) %in% c("A", "environment")]
   testRes <- function(s,v)
      sapply(1:ncol(v), function(i)
         base::norm(A%*%v[,i]-v[,i]*s[i],"2"));
   labels <- (if (!is.null(names(arguments))) names else as.character)(arguments) 
   data.frame(row.names=NULL, test=labels, t(mapply(function(test) {
      r_t <- system.time(r <- eval(test, environment));
      if (!"values" %in% names(r)) r$values <- r$d;
      if (!"vectors" %in% names(r)) r$vectors <- r$u;
      resNorm <- max(testRes(r$values, r$vectors))
      matvecs <- if ("mprod" %in% names(r)) r$mprod
                 else if ("nops" %in% names(r)) r$nops
                 else if ("stats" %in% names(r)) r$stats$numMatvecs
                 else "--";
      list(time=r_t[3], matvecs=matvecs, rnorm=resNorm)
   }, arguments)))
}
```

PRIMME eigs\_sym is based on Davidson-type methods and they may be faster than Lanczos/Arnoldi based method (e.g., svd, RSpectra and irlba) in difficult problems that eigenpairs take many iterations to convergence or an efficient preconditioner is available.

``` r
library(RSpectra, warn.conflicts=FALSE, pos=5)
library(irlba, pos=5)
library(svd, pos=5)

Ad <- diag(1:12000);
for(i in 1:11999) {Ad[i,i+1]<-1; Ad[i+1,i]<-1}
set.seed(1)
r <- bench_eigs(
   PRIMME=PRIMME::eigs_sym(Ad,2,tol=1e-5),
   irlba=partial_eigen(Ad,2,tol=1e-5),
   RSpectra=RSpectra::eigs_sym(Ad,2,tol=1e-5),
   trlan=svd::trlan.eigen(Ad,2,opts=list(tol=1e-5)),
   A=Ad
)
kable(r, digits=2, caption="2 largest eigenvalues on dense matrix")
```

| test     | time    | matvecs | rnorm        |
|:---------|:--------|:--------|:-------------|
| PRIMME   | 13.312  | 500     | 0.1129397    |
| irlba    | 86.767  | --      | 0.04308973   |
| RSpectra | 57.78   | 2192    | 9.512001e-07 |
| trlan    | 324.302 | --      | 0.1197901    |

``` r
Ad <- diag(1:6000);
for(i in 1:5999) {Ad[i,i+1]<-1; Ad[i+1,i]<-1}
P <- diag(Ad);
set.seed(1)
r <- bench_eigs(
   PRIMME=PRIMME::eigs_sym(Ad,5,'SM',tol=1e-7),
   "PRIMME Prec"=PRIMME::eigs_sym(Ad,5,'SM',tol=1e-7,prec=function(x)x/P),
   RSpectra=RSpectra::eigs_sym(Ad,5,'SM',tol=1e-7),
   A=Ad
)
kable(r, digits=2, caption="5 eigenvalues closest to zero on dense matrix")
```

| test        | time  | matvecs | rnorm        |
|:------------|:------|:--------|:-------------|
| PRIMME      | 3.852 | 555     | 0.0005940415 |
| PRIMME Prec | 0.284 | 42      | 0.0004805318 |
| RSpectra    | 9.067 | 1433    | 4.884529e-08 |

By default PRIMME tries to guess the best configuration, but a little hint can help sometimes. The next example sets the preset method `'PRIMME_DEFAULT_MIN_TIME'` that takes advantage of very light matrix-vector products.

``` r
As <- as(sparseMatrix(i=1:50000,j=1:50000,x=1:50000),"dgCMatrix");
for(i in 1:49999) {As[i,i+1]<-1; As[i+1,i]<-1}
P = 1:50000; # Jacobi preconditioner of As
set.seed(1)
r <- bench_eigs(
   "PRIMME defaults"=PRIMME::eigs_sym(As,40,'SM',tol=1e-10),
   "PRIMME min time"=PRIMME::eigs_sym(As,40,'SM',tol=1e-10,method='PRIMME_DEFAULT_MIN_TIME'),
   "PRIMME Prec"=PRIMME::eigs_sym(As,40,'SM',tol=1e-10,prec=function(x)x/P),
   RSpectra=RSpectra::eigs_sym(As,40,'SM',tol=1e-10,opts=list(maxitr=9999)),
   A=As
)
kable(r, digits=2, caption="40 eigenvalues closest to zero on dense matrix")
```

| test            | time   | matvecs | rnorm        |
|:----------------|:-------|:--------|:-------------|
| PRIMME defaults | 73.319 | 13436   | 4.991904e-06 |
| PRIMME min time | 8.332  | 18945   | 4.935499e-06 |
| PRIMME Prec     | 2.153  | 311     | 4.411372e-06 |
| RSpectra        | 14.508 | 4343    | 4.224989e-09 |

Singular value problems
-----------------------

For SVD problems, the package provides a similar interface:

``` r
A <- diag(1:10, 20,10) # rectangular matrix of dimension 20x10
r <- svds(A, 3); # compute the three largest singular values
r
#> $d
#> [1] 10  9  8
#> 
#> $u
#>                [,1]          [,2]          [,3]
#>  [1,] -1.005532e-17 -2.363039e-17 -2.054460e-18
#>  [2,] -1.258396e-17 -8.675434e-18  3.439066e-17
#>  [3,] -7.359263e-18 -3.292225e-17 -8.656467e-18
#>  [4,]  3.835071e-17  4.091713e-17  6.938004e-18
#>  [5,] -1.440351e-17 -2.958001e-17  4.547859e-19
#>  [6,]  7.167005e-17  1.429539e-16  8.137773e-20
#>  [7,] -5.629758e-17  1.196127e-17  3.508478e-16
#>  [8,] -2.642821e-17 -1.260746e-16 -1.000000e+00
#>  [9,]  5.819616e-16 -1.000000e+00  1.740527e-16
#> [10,]  1.000000e+00  1.131234e-15 -4.132377e-17
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00
#> [12,]  0.000000e+00  0.000000e+00  0.000000e+00
#> [13,]  0.000000e+00  0.000000e+00  0.000000e+00
#> [14,]  0.000000e+00  0.000000e+00  0.000000e+00
#> [15,]  0.000000e+00  0.000000e+00  0.000000e+00
#> [16,]  0.000000e+00  0.000000e+00  0.000000e+00
#> [17,]  0.000000e+00  0.000000e+00  0.000000e+00
#> [18,]  0.000000e+00  0.000000e+00  0.000000e+00
#> [19,]  0.000000e+00  0.000000e+00  0.000000e+00
#> [20,]  0.000000e+00  0.000000e+00  0.000000e+00
#> 
#> $v
#>                [,1]          [,2]          [,3]
#>  [1,] -1.005532e-16 -2.126735e-16 -1.643568e-17
#>  [2,] -6.291981e-17 -3.903946e-17  1.375627e-16
#>  [3,] -2.453088e-17 -9.876675e-17 -2.308391e-17
#>  [4,]  9.587679e-17  9.206354e-17  1.387601e-17
#>  [5,] -2.880701e-17 -5.324401e-17  7.276574e-19
#>  [6,]  1.194501e-16  2.144308e-16  1.085036e-19
#>  [7,] -8.042512e-17  1.537878e-17  4.009689e-16
#>  [8,] -3.303526e-17 -1.418340e-16 -1.000000e+00
#>  [9,]  6.466240e-16 -1.000000e+00  1.547135e-16
#> [10,]  1.000000e+00  1.018111e-15 -3.305902e-17
#> 
#> $rnorms
#> [1] 6.280370e-15 6.978189e-15 7.850462e-15
#> 
#> $stats
#> $stats$numMatvecs
#> [1] 20
#> 
#> $stats$numPreconds
#> [1] 0
#> 
#> $stats$elapsedTime
#> [1] 0.0001549721
#> 
#> $stats$estimateANorm
#> [1] 10
#> 
#> $stats$timeMatvec
#> [1] 2.861023e-06
#> 
#> $stats$timePrecond
#> [1] 0
```

The next examples show how to compute the smallest singular values and how to specify some tolerance:

``` r
A <- diag(1:100, 500,100)

r <- svds(A, 3, 'S'); # compute the three smallest values
r$d
#> [1] 1 2 3

r <- svds(A, 3, 'S', tol=1e-5);
r$rnorms # this is should be smaller than ||A||*tol
#> [1] 0.0009489247 0.0007254660 0.0010843363
```

The next example shows the use of a diagonal preconditioner based on \(A^*A\) (see TODO):

``` r
A <- rbind(rep(1,n=100), diag(1:100, 500,100))
r <- svds(A, 3, 'S');
r$stats$numMatvecs
#> [1] 702

P <- colSums(A^2);  # Jacobi preconditioner of Conj(t(A))%*%A
r <- svds(A, 3, 'S', prec=list(AHA=function(x)x/P));
r$stats$numMatvecs
#> [1] 44
```

The next benchmark function extends `rbenchmark` to return besides the time, the number of matrix-vector multiplications and the maximum residual norm of the returned triplets.

``` r
bench_svds <- function(..., A, environment=parent.frame()) {
   arguments = match.call()[-1]
   if (!is.null(names(arguments)))
      arguments = arguments[!names(arguments) %in% c("A", "environment")]
   testRes <- function(u,s,v)
      sapply(1:ncol(u), function(i)
         base::norm(rbind(A%*%v[,i]-u[,i]*s[i], Conj(t(as.matrix(Conj(t(u[,i]))%*%A)))-v[,i]*s[i]),"2"));
   labels <- (if (!is.null(names(arguments))) names else as.character)(arguments) 
   data.frame(row.names=NULL, test=labels, t(mapply(function(test) {
      r_t <- system.time(r <- eval(test, environment));
      if (is.null(r$v)) r$v <- sapply(1:ncol(r$u), function(i) crossprod(A,r$u[,i])/base::norm(crossprod(A,r$u[,i]),"2"));
      resNorm <- max(testRes(r$u, r$d, r$v))
      matvecs <- if ("mprod" %in% names(r)) r$mprod
                 else if ("nops" %in% names(r)) r$nops
                 else if ("stats" %in% names(r)) r$stats$numMatvecs
                 else "--";
      list(time=r_t[3], matvecs=matvecs, rnorm=resNorm)
   }, arguments)))
}
```

PRIMME svds may perform as good as similar methods in the packages svd, RSpectra and irlba in solving few singular values.

``` r
Ad <- matrix(rnorm(6000*6000),6000)
set.seed(1)
r <- bench_svds(
   PRIMME=PRIMME::svds(Ad,2,tol=1e-5),
   irlba=irlba(Ad,2,tol=1e-5),
   RSpectra=RSpectra::svds(Ad,2,tol=1e-5),
   trlan=trlan.svd(Ad,2,opts=list(tol=1e-5)),
   propack=propack.svd(Ad,2,opts=list(tol=1e-5,maxiter=99999)),
   A=Ad
)
kable(r, digits=2, caption="2 largest singular values on dense matrix")
```

| test     | time   | matvecs | rnorm        |
|:---------|:-------|:--------|:-------------|
| PRIMME   | 3.02   | 232     | 0.001487547  |
| irlba    | 4.364  | 342     | 0.001719602  |
| RSpectra | 10.816 | 636     | 2.995926e-09 |
| trlan    | 6.058  | --      | 0.001331501  |
| propack  | 3.072  | --      | 0.001757105  |

PRIMME can take advantage of a light matrix-vector product:

``` r
As <- as(sparseMatrix(i=1:50000,j=1:50000,x=1:50000),"dgCMatrix");
r <- bench_svds(
   PRIMME=PRIMME::svds(As,40,tol=1e-5),
   irlba=irlba(As,40,tol=1e-5,maxit=5000,work=100),
   RSpectra=RSpectra::svds(As,40,tol=1e-5),
   A=As
)
kable(r, digits=2, caption="40 largest singular values on sparse matrix")
```

| test     | time   | matvecs | rnorm        |
|:---------|:-------|:--------|:-------------|
| PRIMME   | 3.7    | 12216   | 0.4924661    |
| irlba    | 12.657 | 4244    | 1.708491     |
| RSpectra | 14.198 | 4236    | 5.444241e-06 |

And for now it is the only package that supports computing the smallest singular values:

``` r
# Get LargeReFile from UF matrix collection
tf <- tempfile();
download.file('http://www.cise.ufl.edu/research/sparse/MM/Stevenson/LargeRegFile.tar.gz',tf);
td <- tempdir();
untar(tf, exdir=td);
As <- as(readMM(paste(td,'LargeRegFile/LargeRegFile.mtx',sep='/')), "dgCMatrix");
unlink(tf)
unlink(td, recursive=TRUE)

P <- colSums(As^2);  # Jacobi preconditioner of Conj(t(A))%*%A
r <- bench_svds(
   PRIMME=PRIMME::svds(As,5,'S',tol=1e-10),
   "PRIMME Prec"=PRIMME::svds(As,5,'S',tol=1e-10,prec=list(AHA=function(x)x/P)),
   A=As
)
kable(r, digits=2, caption="5 smallest singular values on sparse matrix")
```

| test        | time    | matvecs | rnorm        |
|:------------|:--------|:--------|:-------------|
| PRIMME      | 438.589 | 25528   | 2.715613e-07 |
| PRIMME Prec | 19.636  | 1086    | 2.906932e-07 |

TODO
====

-   Optimize the application of preconditioner when it is passed as a dense or sparse matrix. When solving small problems the overhead of calling the R function that applies the preconditioner can dominate over the reduction of iterations:

``` r
# A is a tridiagonal
A <- diag(1:1000)
for(i in 1:999) {A[i,i+1]<-1; A[i+1,i]<-1}

r <- eigs_sym(A, 10, 'SA');
r$stats$numMatvecs
#> [1] 698
r$stats$elapsedTime
#> [1] 0.06746888

# Jacobi preconditioner
P = diag(diag(A));
r <- eigs_sym(A, 10, 'SA', prec=P);
r$stats$numMatvecs
#> [1] 58
r$stats$elapsedTime
#> [1] 1.027588
```

-   Add support for matrices distributed among processes.
