% Show primme_svds using Afun and Pfun for smallest singular triplets

clear,clc

% A = delsq(numgrid('C',100));
load deter4.mat;
A = Problem.A; 
[m,n] = size(A);
global primmeA;
primmeA = A;    

% Compute the 5 smallest singular values, using a preconditioner only
% for the first stage

Pstruct = struct('AHA', diag(A'*A), ...
                 'AAH', ones(size(A,1), 1), 'aug', ones(size(A,1)+size(A,2), 1));
Pfun = @(x,mode)Pstruct.(mode).\x;
svals = primme_svds(A, 5, 'S', [], Pfun);

% Use ILU to generate a preconditioner for A. 
% ILU_time = tic;
% [L,U] = ilu(A,struct('type','ilutp','droptol',1e-3,'thresh', 1.0));
% precond_runtime = toc(ILU_time);

% Use RIF to generate a preconditioner directly for ATA or AAT. 
RIF_begin = tic;
rifthresh = 1E-3;
zrifthresh = 1E-9;
rifnnz = 0;
if m >= n
    C = A;
else
    C = A';
end
global LL;
LL = RIF(C,0,rifthresh,zrifthresh,rifnnz);
precond_runtime = toc(RIF_begin)

K = 5;
sigma = 'S';
opts = struct('tol',1E-6,'method','primme_svds_normalequations');
primme_svds_time = tic;
[primme_U, primme_S, primme_V, Resnorm, primmeout]= ...
    primme_svds(@(x,tflag)Afun(x,tflag), m, n, K, sigma, opts);
primme_svds_runtime = toc(primme_svds_time)   
matvec = primmeout.numMatvecs
primme_S = diag(primme_S)
sval = primme_S';
for j = 1:length(primme_S)
    res_U(j) = norm(A'*primme_U(:,j) - sval(j)*primme_V(:,j));
    res_V(j) = norm(A*primme_V(:,j) - sval(j)*primme_U(:,j));
    res_OAAO(j) = sqrt((res_U(j)^2 + res_V(j)^2)/...
        (norm(primme_V(:,1))^2 + norm(primme_U(:,1))^2));
    res_ATA(j) = norm(A'*(A*primme_V(:,j)) - (sval(j)^2)*primme_V(:,j));          
end
primme_svds_Norm.res_U = res_U;
primme_svds_Norm.res_V = res_V;       
primme_svds_Norm.res_OAAO = res_OAAO;
primme_svds_Norm.res_ATA = res_ATA;
primme_svds_Norm

primme_svds_time = tic;
[primme_U, primme_S, primme_V, Resnorm, primmeout]= ...
    primme_svds(@(x,tflag)Afun(x,tflag), m, n, K, sigma, opts, @(x,model)Pfun(x,model));
primme_svds_runtime = toc(primme_svds_time) + precond_runtime   
matvec = primmeout.numMatvecs
primme_S = diag(primme_S)
sval = primme_S';
for j = 1:length(primme_S)
    res_U(j) = norm(A'*primme_U(:,j) - sval(j)*primme_V(:,j));
    res_V(j) = norm(A*primme_V(:,j) - sval(j)*primme_U(:,j));
    res_OAAO(j) = sqrt((res_U(j)^2 + res_V(j)^2)/...
        (norm(primme_V(:,1))^2 + norm(primme_U(:,1))^2));
    res_ATA(j) = norm(A'*(A*primme_V(:,j)) - (sval(j)^2)*primme_V(:,j));          
end
primme_svds_Norm.res_U = res_U;
primme_svds_Norm.res_V = res_V;       
primme_svds_Norm.res_OAAO = res_OAAO;
primme_svds_Norm.res_ATA = res_ATA;
primme_svds_Norm