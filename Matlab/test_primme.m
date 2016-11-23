%  Copyright (c) 2016, College of William & Mary
%  All rights reserved.
%  
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%      * Redistributions of source code must retain the above copyright
%        notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in the
%        documentation and/or other materials provided with the distribution.
%      * Neither the name of College of William & Mary nor the
%        names of its contributors may be used to endorse or promote products
%        derived from this software without specific prior written permission.
%  
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%  DISCLAIMED. IN NO EVENT SHALL COLLEGE OF WILLIAM & MARY BE LIABLE FOR ANY
%  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%  
%  PRIMME: https://github.com/primme/primme
%  Contact: Andreas Stathopoulos, a n d r e a s _at_ c s . w m . e d u


% Compute the 6 largest eigenvalues of a matrix with tolerance 1e-6

A = diag(1:50);
ops = struct();
ops.tol = 1e-6; % residual norm tolerance 
k = 6;          % number of eigenvalues
evals = primme_eigs(A, k, 'LA', ops);

assert(norm(evals - (50:-1:50-k+1)') < 1e-6*norm(A))

% Compute in single arithmetic the 6 eigenvalues closest to 25.2
% using DEFAULT_MIN_TIME

[evecs, evals, rnorms, stats] = primme_eigs(single(A), k, 25.2, ops, ...
                                            'DEFAULT_MIN_TIME');

[a,p] = sort(abs((1:50) - 25.2));
assert(all(abs(diag(evals) - p(1:k)') < 1e-6*stats.estimateAnorm))
assert(all(rnorms < 1e-6*stats.estimateAnorm));

% Compute the 6 smallest eigenvalues and vectors of a matrix defined by
% the matrix-vector product

fun = @(x)A*x;
matrix_dim = 50;
[evecs, evals] = primme_eigs(fun, matrix_dim, k, 'SA', ops);

assert(norm(diag(evals) - (1:k)') < 1e-6*norm(A))
for i=1:k
  assert(norm(A*evecs(:,i) - evecs(:,i)*evals(i,i)) < 1e-6*norm(A))
end

% Compute the 6 closest eigenvalues to 30.5 using the Jacobi preconditioner
% (too much convenient for a diagonal matrix)

Adiag = diag(A);
Pfun = @(x)((Adiag - 30.5).\x);
evals = primme_eigs(A, k, 30.5, [], [], Pfun);
 
% Compute the 6 largest singular values of a matrix with tolerance 1e-6

A = diag(1:50); A(200,1) = 0; % rectangular matrix of size 200x50
ops = struct();
ops.eps = 1e-6; % residual norm tolerance 
k = 6;          % number of singular values
svals = primme_svds(A, k, 'L', ops);

assert(norm(svals - (50:-1:50-k+1)') < 1e-6*norm(A))

% Compute the 6 smallest singular values and vectors of a matrix defined by
% the matrix-vector product

% Mathwork's MATLAB doesn't allow to define a function here... The next trick
% is an easy bypass to do different things without using 'if'.

structA = struct('notransp', A, 'transp', A');
fun = @(x,mode)structA.(mode)*x;
matrix_dim_m = 200;
matrix_dim_n = 50;

[svecsl, svals, svecsr] = primme_svds(fun, matrix_dim_m, matrix_dim_n, k, 'S', ops);

assert(norm(diag(svals) - (1:k)') < 1e-6*norm(A))
for i=1:k
  assert(norm(A*svecsr(:,i) - svecsl(:,i)*svals(i,i)) < 1e-6*norm(A))
end

% Compute the 10 singular values closest to 30.5 using the normal equation
% approach, with a tolerance of 1e-4. Also set the method (DEFAULT_MIN_TIME) and
% the block size (2) to the solver that computes the eigenvalues on A'*A.

ops = struct();
ops.tol = 1e-4; % set tolerance
ops.method = 'primme_svds_normalequations'; % set solver method
ops.primme.method = 'DEFAULT_MIN_TIME'; % set first stage solver method
ops.primme.maxBlockSize = 2; % set block size for first stage
[svecsl, svals, svecsr] = primme_svds(A, 10, 'S', ops);

% Compute the 5 smallest singular values, using a preconditioner only
% for the first stage

Pstruct = struct('AHA', diag(A'*A), ...
                 'AAH', ones(200, 1), 'aug', ones(250, 1));
Pfun = @(x,mode)Pstruct.(mode).\x;
svals = primme_svds(A, 5, 'S', [], Pfun);

disp('Success');
