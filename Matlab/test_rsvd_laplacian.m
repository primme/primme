% Compare primme_svds against rsvd for largest singular triplets

clear,clc
addpath('~/Documents/Work/GITHUB/randSVDPack/matlab_code')

nodes = 15:15:150;
for i = 1:10
    nd = nodes(i);
    A = delsq(numgrid('C',nd));
    [m,n] = size(A);
    matrixSize(i) = m;
    K = 100;
    sigma = 'L';
    opts = struct('tol',1E-1,'v0',{ones(min(m,n),1)./sqrt(m)},'p',15);
    primme_svds_time = tic;
    [primme_U, primme_S, primme_V, Resnorm, primmeout]= ...
        primme_svds(A, K, sigma, opts);
    runtime_primmesvds(i) = toc(primme_svds_time);   
    matvec(i) = primmeout.numMatvecs;
    sigma_primmesvds = diag(primme_S);
    perror_primmesvds(i) = normest(A - primme_U*primme_S*primme_V')/normest(A);
    res_U = normest(A'*primme_U - primme_V*primme_S);
    res_V = norm(A*primme_V - primme_U*primme_S);
    reserror_primmesvds(i) = sqrt((res_U^2 + res_V^2)/(norm(primme_U)^2 + norm(primme_V)^2));
    
    % set params for variants of random svds
    k = 100;
    p = 50;
    q = 2;
    s = 1;

    fprintf('rsvd 1..\n');
    if length(A) - k < p
        p = length(A) - k;
    end
    start_rsvd1 = tic;
    [U,S1,V] = rsvd_version1(A,k,p,q,s);
    whos M U Sigma V
    sigma_rsvd1 = diag(S1);
    runtime_rsvd1(i) = toc(start_rsvd1);
    perror_rsvd1(i) = normest(A - U*S1*V')/normest(A); 
    reserror_rsvd1(i) = normest(A*V - U*S1);

    fprintf('rsvd 2..\n');
    start_rsvd2 = tic;
    [U,S2,V] = rsvd_version2(A,k,p,q,s);
    whos M U Sigma V
    sigma_rsvd2 = diag(S2);
    runtime_rsvd2(i) = toc(start_rsvd2);
    perror_rsvd2(i) = normest(A - U*S2*V')/normest(A);
    reserror_rsvd2(i) = normest(A*V - U*S2);
    
    fprintf('rsvd 3..\n');
    kstep = 20;
    q = 2;
    s = 1;
    start_rsvd3 = tic;
    [U,S3,V] = rsvd_version3(A,K,kstep,q,s);
    sigma_rsvd3 = diag(S3);
    runtime_rsvd3(i) = toc(start_rsvd3);
    perror_rsvd3(i) = normest(A - U*S3*V')/normest(A);
    reserror_rsvd3(i) = normest(A*V - U*S3);
end
runtime_primmesvds
runtime_rsvd1
runtime_rsvd2
runtime_rsvd3

set(0,'defaultaxesfontsize',18);
h1 = figure(1);
plot(matrixSize,runtime_primmesvds,'r*-','LineWidth',2,'MarkerSize',10); 
hold on
plot(matrixSize,runtime_rsvd1,'bo-.','LineWidth',2,'MarkerSize',10);
plot(matrixSize,runtime_rsvd2,'bs--','LineWidth',2,'MarkerSize',10); 
plot(matrixSize,runtime_rsvd3,'bx-','LineWidth',2,'MarkerSize',10); 
hold off
title(['Laplacian2D: Find ',num2str(K),' largest singular values'])
legend('primme\_svds', 'rsvd1', 'rsvd2', 'rsvd3', 'Location','northwest')
xlabel('Size of Matrix')
ylabel('Runtime (Seconds)')
filename = ['Lap_runtime_svds-1E-1','-num',num2str(K),'-',sigma];
saveas(h1,filename,'fig');
saveas(h1,filename,'png');

set(0,'defaultaxesfontsize',18);
h2 = figure(2);
plot(matrixSize,reserror_primmesvds,'r*-','LineWidth',2,'MarkerSize',10); 
hold on
plot(matrixSize,reserror_rsvd1,'bo-.','LineWidth',2,'MarkerSize',10);
plot(matrixSize,reserror_rsvd2,'bs--','LineWidth',2,'MarkerSize',10); 
plot(matrixSize,reserror_rsvd3,'bx-','LineWidth',2,'MarkerSize',10); 
hold off
title(['Laplacian2D: Find ',num2str(K),' largest singular values'])
legend('primme\_svds', 'rsvd1', 'rsvd2', 'rsvd3', 'Location','southeast')
xlabel('Size of Matrix')
ylabel('Residual Norm of Rank-K Approx')
filename = ['Lap_reserror_svds-1E-1','-num',num2str(K),'-',sigma];
saveas(h2,filename,'fig');
saveas(h2,filename,'png');