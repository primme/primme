% Compare primme_svds against svds for smallest singular triplets

% clear,clc
nodes = 5:5:25;
% nodes = 100;
matrixSize = zeros(length(nodes),1);
primme_svds_runtime = zeros(length(nodes),1);
svds_runtime = zeros(length(nodes),1);
matvec = zeros(length(nodes),1);
for i = 1:length(nodes)
    nd = nodes(i);
%     A = delsq(numgrid('C',nd)); % 2D laplacian
    [~,~,A]=laplacian([nd, nd, nd]); % 3D laplacian
%     A = gallery('neumann',nd); % miss zero singular value
    ILU_time = tic;
    [L,U] = ilu(A,struct('type','ilutp','droptol',1e-3,'thresh', 1.0));
    precond_runtime = toc(ILU_time);

    [m,n] = size(A);
    matrixSize(i) = m;
    K = 6;
    sigma = 'S';
    rng('default');
    opts = struct('tol',1E-14,'v0',{ones(min(m,n),1)./sqrt(m)},'p',35);
    primme_svds_time = tic;
    [primme_U, primme_S, primme_V, Resnorm, primmeout]= ...
        primme_svds(A, K, sigma, opts);
    primme_svds_runtime(i) = toc(primme_svds_time) + precond_runtime;   
    matvec(i) = primmeout.numMatvecs;
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
    
    opts.maxit = 1000;
    opts.v0 = ones(size(A,2),1)./sqrt(m);
    sigma = 'smallest';
    svds_time = tic;
    [U,S,V,FLAG] = svds(A, K, sigma, opts);
    svds_runtime(i) = toc(svds_time);   
    S = diag(S)
    sval = S';
    for j = 1:length(S)
        res_U(j) = norm(A'*U(:,j) - sval(j)*V(:,j));
        res_V(j) = norm(A*V(:,j) - sval(j)*U(:,j));
        res_OAAO(j) = sqrt((res_U(j)^2 + res_V(j)^2)/...
            (norm(V(:,1))^2 + norm(U(:,1))^2));
        res_ATA(j) = norm(A'*(A*V(:,j)) - (sval(j)^2)*V(:,j));          
    end
    svds_Norm.res_U = res_U;
    svds_Norm.res_V = res_V;       
    svds_Norm.res_OAAO = res_OAAO;
    svds_Norm.res_ATA = res_ATA;
    svds_Norm
end
primme_svds_runtime
svds_runtime

set(0,'defaultaxesfontsize',18);
h1 = figure(1);
plot(matrixSize,primme_svds_runtime,'r*-','LineWidth',2,'MarkerSize',10); 
hold on
plot(matrixSize,svds_runtime,'bo-','LineWidth',2,'MarkerSize',10); 
hold off
title(['Find ',num2str(K),' smallest singular values'])
legend('primme\_svds(ILU)', 'svds(QR)','Location','northwest')
xlabel('Size of Matrix')
ylabel('Runtime')
filename = ['svds-',num2str(opts.tol),'-num',num2str(K),'-',sigma];
saveas(h1,filename,'fig');
saveas(h1,filename,'epsc');