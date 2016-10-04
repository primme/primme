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
ops.eps = 1e-6; % residual norm tolerance 
k = 6;          % number of eigenvalues
evals = primme_eigs(A, k, 'LA', ops);

assert(norm(evals - (50:-1:50-k+1)') < 1e-6*norm(A))

% Compute the 6 largest eigenvalues and vectors of a matrix defined by
% the matrix-vector product

fun = @(x)A*x;
matrix_dim = 50;
[evecs, evals] = primme_eigs(fun, matrix_dim, 6, 'LA', ops);

assert(norm(diag(evals) - (50:-1:50-k+1)') < 1e-6*norm(A))
for i=1:k
  assert(norm(A*evecs(:,i) - evecs(:,i)*evals(i,i)) < 1e-6*norm(A))
end

disp('Success');
