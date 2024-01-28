% Matlab PRIMME comparisons run on a core i7-9700k with 32GB of RAM available
A = sparse(diag(1:12000) + diag(ones(11999,1),1) + diag(ones(11999,1), -1));
ops = struct();
ops.tol =1e-5;

fprintf("----------------\n");
fprintf("Eigenvalue Tests\n");
fprintf("----------------\n");

fprintf("Strongly diagonally dominant tridiagional matrix, finding the 2 largest eigenvalues\n\n");
fprintf("test\ttime\tmatvecs\trnorms\n");
% Primme Performance test
tic;
[evecs, evals, rnorms, stats] = primme_eigs(A, 2, 'LA', ops);
t = toc();
fprintf("Primme:\t%.3f\t%d\t%.3f\n", t, stats.numMatvecs, max(rnorms));

% eigs Performance test
tic;
[evecs, evals] = eigs(A, 2, 'largestreal', ops);
t = toc();
rnorms = vecnorm((A*evecs) - (evecs * evals));
fprintf("eigs:\t%.3f\t---\t%.3f\n\n", t, max(rnorms));

fprintf("---------------------------\n");
fprintf("Primme Preconditioned Tests\n");
fprintf("---------------------------\n");

fprintf("Strongly diagonally dominant tridiagional matrix, finding the 2 smallest eigenvalues\n\n");
fprintf("test\t\ttime\tmatvecs\trnorms\n");
% Primme Performance test
tic;
[evecs, evals, rnorms, stats] = primme_eigs(A, 2, 'SM', ops);
t = toc();
fprintf("Primme:\t\t%.3f\t%d\t%.3f\n", t, stats.numMatvecs, max(rnorms));

diagInv = diag(((1:12000)), 0);
% Primme Preconditioned Performance test
tic;
[evecs, evals, rnorms, stats] = primme_eigs(A, 2, 'SM', ops, [], diagInv);
t = toc();
fprintf("Precondition:\t%.3f\t%d\t%.3f\n\n", t, stats.numMatvecs, max(rnorms));

fprintf("---------------------\n");
fprintf("Singular Value Tests\n");
fprintf("---------------------\n");

fprintf("A = randn(6000,6000);\t finding the largest 2 singular values\n\n");
A = randn(6000, 6000);
fprintf("test\ttime\tmatvecs\trnorms\n");

%Primme Performance test
tic;
[u, svals, v, rnorms, stats] = primme_svds(A, 2, 'L', ops);
t = toc();
fprintf("Primme:\t%.3f\t%d\t%.3f\n", t, stats.numMatvecs, max(rnorms));

%svds Performance test
tic;
[u, svals, v] = svds(A, 2, 'largest', ops);
t = toc();
rnorms = sqrt(vecnorm(A*v - u * svals).^2 + vecnorm(A'*u - v * svals).^2); 
fprintf("svds:\t%.3f\t---\t%.3f\n\n", t, max(rnorms));

fprintf("--------------------------\n");
fprintf("Large Singular Value Tests\n");
fprintf("--------------------------\n");

fprintf("A = randn(6000,6000);\t finding the largest 100 singular values\n\n");
fprintf("test\ttime\tmatvecs\trnorms\n");
%Primme Performance test
ops.maxBlockSize = 100;
tic;
[u, svals, v, rnorms, stats] = primme_svds(A, 100, 'L', ops);
t = toc();
fprintf("Primme:\t%.3f\t%d\t%.3f\n", t, stats.numMatvecs, max(rnorms));

ops = struct();
ops.tol = 1e-5;
%svds Performance test
tic;
[u, svals, v] = svds(A, 100, 'largest', ops);
t = toc();
rnorms = sqrt(vecnorm(A*v - u * svals).^2 + vecnorm(A'*u - v * svals).^2); 
fprintf("svds:\t%.3f\t---\t%.3f\n", t, max(rnorms));
