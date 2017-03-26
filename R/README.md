<!-- README.md is generated from README.Rmd. Please edit that file -->
PRIMME
======

This package is an interface to PRIMME, a C library for computing a few eigenvalues and their corresponding eigenvectors of a real symmetric or complex Hermitian matrix. It can also compute singular values and vectors. It can find largest, smallest, or interior eigenvalues or singular values and can use preconditioning to accelerate convergence. It is a useful tool for both non-experts and experts.

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

The next example computes the three largest eigenvalues of the matrix `A`, which in this case is a dense diagonal matrix. It shows all the eigenvalues `evals`, the eigenvectors `evecs`, the residual norms `rnorms` and some stats, such as the time `stats$elapsedTime` and the number of matrix vector multiplications performed `stats$numMatvecs`:

``` r
A <- diag(1:10) 
r <- primme.eigs_symm(A, 3);
r
#> $evals
#> [1] 10  9  8
#> 
#> $evecs
#>                [,1]          [,2]          [,3]
#>  [1,] -7.907476e-17  1.440058e-16  1.683224e-17
#>  [2,] -2.158240e-17 -4.239230e-17 -2.406929e-17
#>  [3,] -3.016793e-17 -1.602451e-16 -1.266348e-16
#>  [4,] -1.908196e-17  2.211772e-17 -2.532696e-16
#>  [5,] -1.526557e-16 -1.196959e-16  2.359224e-16
#>  [6,]  6.005938e-16 -2.159731e-16 -1.040834e-16
#>  [7,] -9.470506e-17  3.174544e-16 -1.960238e-15
#>  [8,] -5.657909e-16  3.992032e-16 -1.000000e+00
#>  [9,]  2.491104e-15 -1.000000e+00 -5.095750e-16
#> [10,]  1.000000e+00  2.169617e-15 -6.078579e-16
#> 
#> $rnorms
#> [1] 4.440892e-15 4.440892e-15 4.440892e-15
#> 
#> $stats
#> $stats$numMatvecs
#> [1] 10
#> 
#> $stats$elapsedTime
#> [1] 0.0003149509
#> 
#> $stats$estimateMinEval
#> [1] 1
#> 
#> $stats$estimateMaxEval
#> [1] 10
#> 
#> $stats$estimateANorm
#> [1] 10
```

The next examples show how to compute eigenvalues in other parts of the spectrum:

``` r
A <- diag(1:10)

r <- primme.eigs_symm(A, 3, 'SA'); # compute the three smallest values
r$evals
#> [1] 1 2 3

r <- primme.eigs_symm(A, 3, 5.1); # compute the three closest values to 5.1
r$evals
#> [1] 5 6 4
```

In some cases, a larger convergence tolerance may suffice:

``` r
A <- diag(1:1000)

r <- primme.eigs_symm(A, 10, 'SA');
r$stats$numMatvecs
#> [1] 643
r$stats$elapsedTime
#> [1] 1.645482

r <- primme.eigs_symm(A, 10, 'SA', tol=1e-4); 
r$stats$numMatvecs
#> [1] 359
r$stats$elapsedTime
#> [1] 0.9441411
```

Preconditioners, if available can reduce the time/matrix-vector multiplications significantly (see TODO):

``` r
# A is a tridiagonal
A <- diag(1:1000)
for(i in 1:999) {A[i,i+1]<-1; A[i+1,i]<-1}

r <- primme.eigs_symm(A, 10, 'SA');
r$stats$numMatvecs
#> [1] 565

# Jacobi preconditioner
P = diag(diag(A));
r <- primme.eigs_symm(A, 10, 'SA', prec=P);
r$stats$numMatvecs
#> [1] 58
```

Dense matrices, sparse matrices, and functions that return the matrix-vector product can be passed as the matrix problem `A`:

``` r
r <- primme.eigs_symm(diag(1:10), 1); # dense matrix
require(Matrix)
r <- primme.eigs_symm(Matrix(diag(1:10), sparse=TRUE), 1); # sparse matrix
Afun = function(x) matrix(1:10)*x;  # function that does diag(1:10) %*% x
r <- primme.eigs_symm(Afun, 1, n=10); # n is the matrix dimension corresponding to Afun
```

Singular value problems
-----------------------

For SVD problems, the package provides a similar interface:

``` r
A <- diag(1:10, 20,10) # rectangular matrix of dimension 20x10
r <- primme.svds(A, 3); # compute the three largest singular values
r
#> $svals
#> [1] 10  9  8
#> 
#> $svecsu
#>                [,1]          [,2]          [,3]
#>  [1,] -2.342300e-17 -6.664832e-18 -4.472334e-19
#>  [2,] -3.604701e-17  2.200328e-17  4.239230e-17
#>  [3,]  2.289022e-18 -4.674718e-17 -4.878910e-17
#>  [4,]  3.568651e-17  1.611847e-17  1.040834e-17
#>  [5,] -2.662394e-17 -2.523782e-17  2.981556e-17
#>  [6,]  2.483909e-16 -7.976114e-17  6.960578e-17
#>  [7,] -1.623183e-16  6.806230e-17 -2.533442e-16
#>  [8,]  9.121977e-17 -3.640058e-16 -1.000000e+00
#>  [9,]  1.465345e-16 -1.000000e+00  1.934088e-16
#> [10,]  1.000000e+00  3.636782e-16  2.432348e-17
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
#> $svecsv
#>                [,1]          [,2]          [,3]
#>  [1,] -2.342300e-16 -5.998349e-17 -3.577867e-18
#>  [2,] -1.802351e-16  9.901476e-17  1.695692e-16
#>  [3,]  7.630073e-18 -1.402416e-16 -1.301043e-16
#>  [4,]  8.921629e-17  3.626656e-17  2.081668e-17
#>  [5,] -5.324788e-17 -4.542807e-17  4.770490e-17
#>  [6,]  4.139848e-16 -1.196417e-16  9.280771e-17
#>  [7,] -2.318833e-16  8.750867e-17 -2.895362e-16
#>  [8,]  1.140247e-16 -4.095065e-16 -1.000000e+00
#>  [9,]  1.628161e-16 -1.000000e+00  1.719189e-16
#> [10,]  1.000000e+00  3.273104e-16  1.945878e-17
#> 
#> $rnorms
#> [1] 6.332677e-15 1.191533e-14 7.850462e-15
#> 
#> $stats
#> $stats$numMatvecs
#> [1] 20
#> 
#> $stats$elapsedTime
#> [1] 0.0071311
#> 
#> $stats$estimateANorm
#> [1] 10
```

The next examples show how to compute the smallest singular values and how to specify some tolerance:

``` r
A <- diag(1:100, 500,100)

r <- primme.svds(A, 3, 'S'); # compute the three smallest values
r$svals
#> [1] 1 2 3

r <- primme.svds(A, 3, 'S', tol=1e-5);
r$rnorms # this is should be smaller than ||A||*tol
#> [1] 0.0010734258 0.0009443993 0.0011803952
```

The next example shows the use of a diagonal preconditioner based on \(A^*A\) (see TODO):

``` r
A <- rbind(rep(1,n=100), diag(1:100, 500,100))
r <- primme.svds(A, 3, 'S');
r$stats$numMatvecs
#> [1] 764

P = diag(diag(crossprod(A)));
r <- primme.svds(A, 3, 'S', prec=list(AHA=P));
r$stats$numMatvecs
#> [1] 44
```

TODO
====

-   Optimize the application of preconditioner when it is passed as a dense or sparse matrix. When solving small problems the overhead of calling the R function that applies the preconditioner can dominate over the reduction of iterations:

``` r
# A is a tridiagonal
A <- diag(1:1000)
for(i in 1:999) {A[i,i+1]<-1; A[i+1,i]<-1}

r <- primme.eigs_symm(A, 10, 'SA');
r$stats$numMatvecs
#> [1] 571
r$stats$elapsedTime
#> [1] 1.565808

# Jacobi preconditioner
P = diag(diag(A));
r <- primme.eigs_symm(A, 10, 'SA', prec=P);
r$stats$numMatvecs
#> [1] 64
r$stats$elapsedTime
#> [1] 4.666244
```

-   Add support for matrices distributed among processes.
