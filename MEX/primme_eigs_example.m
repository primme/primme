
function primme_eigs_example

% Global variables used in test_eigs
global Matvec_counter A;

% Matrix problem
A = delsq(numgrid('C', 200));
% Number of eigenvalues
numEvals = 5;   
% Target: largest eigenvalues
target = 'LA';
primme_start = tic;
[primme_V,primme_D,norms,primmeout] = primme_eigs(A, numEvals, target);
primme_time_elapsed = toc(primme_start)
primme_evalues = diag(primme_D)'
primme_numMatvec = primmeout(3)


Matvec_counter = 0;
n = size(A,1);
opts.tol = 1e-12;
opts.issym = 1;
opts.isreal = 1;
eigs_start = tic;

[V ,D] = eigs(@test_eigs, n , numEvals, 'LA', opts);
eigs_time_elapsed = toc(eigs_start)
eigs_evalues = diag(D)'
eigs_numMatvec = Matvec_counter
 
end

function y = test_eigs(x)
   global Matvec_counter A;
   Matvec_counter = Matvec_counter + 1;
   y = A*x;
end

