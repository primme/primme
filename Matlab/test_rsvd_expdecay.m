% test the different rsvd algorithms

clear,clc
addpath('~/Documents/Work/GITHUB/randSVDPack/matlab_code')
make_matrix = 1;
save_matrix = 0;

for i = 1:10
    if make_matrix == 1
        fprintf('making matrix..\n');
        m = 2000*i; n = 1000*i; p = min(m,n);
        S = logspace(0,-2,p);
        S = diag(S);
        mat_filename = ['M_m', num2str(m/1000), 'k_n', num2str(n/1000), 'k.mat']
        A = make_matrix1(m,n,S,mat_filename,save_matrix);
    else 
        load('M_m10k_n5k.mat');
        m = size(A,1);
        n = size(A,2);
    end
    matrixSize(i) = m;

    % set params for variants of random svds
    k = 10;
    p = 20;
    q = 2;
    s = 1;

    fprintf('rsvd 1..\n');
    start_rsvd1 = tic;
    [U,S1,V] = rsvd_version1(A,k,p,q,s);
    whos M U Sigma V
    sigma_rsvd1 = diag(S1);
    runtime_rsvd1(i) = toc(start_rsvd1);
    perror_rsvd1(i) = norm(A - U*S1*V')/norm(A); 
    reserror_rsvd1(i) = norm(A*V - U*S1);

    fprintf('rsvd 2..\n');
    start_rsvd2 = tic;
    [U,S2,V] = rsvd_version2(A,k,p,q,s);
    whos M U Sigma V
    sigma_rsvd2 = diag(S2);
    runtime_rsvd2(i) = toc(start_rsvd2);
    perror_rsvd2(i) = norm(A - U*S2*V')/norm(A);
    reserror_rsvd2(i) = norm(A*V - U*S2);

    fprintf('rsvd 3..\n');
    kstep = 20;
    start_rsvd3 = tic;
    [U,S3,V] = rsvd_version3(A,k,kstep,q,s);
    sigma_rsvd3 = diag(S3);
    runtime_rsvd3(i) = toc(start_rsvd3);
    perror_rsvd3(i) = norm(A - U*S3*V')/norm(A);
    reserror_rsvd3(i) = norm(A*V - U*S3);

    % set params for primme_svds
    fprintf('primme_svds ...\n');
    opts.tol = 1e-1;
    opts.p = 15;
    opts.primme.maxBlockSize = 1;
    opts.disp = 0;
    start_primmesvds = tic;
    [U,PS,V] = primme_svds(A,k,'L',opts);
    sigma_primmesvds = diag(PS);
    runtime_primmesvds(i) = toc(start_primmesvds);
    perror_primmesvds(i) = norm(A - U*PS*V')/norm(A);
    res_U = norm(A'*U - V*PS);
    res_V = norm(A*V - U*PS);
    reserror_primmesvds(i) = sqrt((res_U^2 + res_V^2)/(norm(U)^2 + norm(V)^2));
end
runtime_rsvd1
runtime_rsvd2
runtime_rsvd3
runtime_primmesvds

set(0,'defaultaxesfontsize',18);
h3 = figure(3);
plot(matrixSize,runtime_primmesvds,'r*-','LineWidth',2,'MarkerSize',10); 
hold on
plot(matrixSize,runtime_rsvd1,'bo-.','LineWidth',2,'MarkerSize',10);
plot(matrixSize,runtime_rsvd2,'bs--','LineWidth',2,'MarkerSize',10); 
plot(matrixSize,runtime_rsvd3,'bx-','LineWidth',2,'MarkerSize',10); 
hold off
title(['ExpDecay: Find ',num2str(k),' largest singular values'])
legend('primme\_svds', 'rsvd1', 'rsvd2', 'rsvd3', 'Location','northwest')
xlabel('Size of Matrix')
ylabel('Runtime (Seconds)')
filename = ['expdecay_runtime_svds-1E-1','-num',num2str(k),'-L'];
saveas(h3,filename,'fig');
saveas(h3,filename,'png');

set(0,'defaultaxesfontsize',18);
h4 = figure(4);
plot(matrixSize,reserror_primmesvds,'r*-','LineWidth',2,'MarkerSize',10); 
hold on
plot(matrixSize,reserror_rsvd1,'bo-.','LineWidth',2,'MarkerSize',10);
plot(matrixSize,reserror_rsvd2,'bs--','LineWidth',2,'MarkerSize',10); 
plot(matrixSize,reserror_rsvd3,'bx-','LineWidth',2,'MarkerSize',10); 
hold off
title(['ExpDecay: Find ',num2str(k),' largest singular values'])
legend('primme\_svds', 'rsvd1', 'rsvd2', 'rsvd3', 'Location','northwest')
xlabel('Size of Matrix')
ylabel('Residual Norm of Rank-K Approx')
filename = ['expdecay_reserror_svds-1E-1','-num',num2str(k),'-L'];
saveas(h4,filename,'fig');
saveas(h4,filename,'png');
