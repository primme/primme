% Compare primme_svds against svds for largest singular triplets

clear,clc
nodes = 15:15:150;
matrixSize = zeros(length(nodes),1);
primme_svds_runtime = zeros(length(nodes),1);
svds_runtime = zeros(length(nodes),1);
matvec = zeros(length(nodes),1);
for i = 1:10
    nd = nodes(i);
    A = delsq(numgrid('C',nd));
    [m,n] = size(A);
    matrixSize(i) = m;
    K = 100;
    sigma = 'L';
    opts = struct('tol',1E-10,'v0',{ones(min(m,n),1)./sqrt(m)},'p',15);
    primme_svds_time = tic;
    [primme_U, primme_S, primme_V, Resnorm, primmeout]= ...
        primme_svds(A, K, sigma, opts);
    primme_svds_runtime(i) = toc(primme_svds_time);   
    matvec(i) = primmeout.numMatvecs;
    primme_S = diag(primme_S);
    
    opts.maxit = 1000;
    opts.p = K+2; %3*K;
    sigma = 'largest';
    svds_time = tic;
    [U,S,V,FLAG] = svds(A, K, sigma, opts);
    svds_runtime(i) = toc(svds_time);   
    S = diag(S);
end
primme_svds_runtime
svds_runtime

set(0,'defaultaxesfontsize',18);
h1 = figure(1);
plot(matrixSize,primme_svds_runtime,'r*-','LineWidth',2,'MarkerSize',10); 
hold on
plot(matrixSize,svds_runtime,'bo-','LineWidth',2,'MarkerSize',10); 
hold off
title(['Find ',num2str(K),' largest singular values'])
legend('primme\_svds', 'svds','Location','northwest')
xlabel('Size of Matrix')
ylabel('Runtime')
filename = ['svds-',num2str(opts.tol),'-num',num2str(K),'-',sigma];
saveas(h1,filename,'fig');
saveas(h1,filename,'epsc');