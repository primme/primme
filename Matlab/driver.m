% TEST DRIVER
% 1. Generate test matrix A = SDU + N/xi
%       where
%       S - n x m, S(i,j) approx N(0,1)
%       D - m x m, diagonal gives linearly diminishing signal singular values
%       U - m x d, signal row space U*U' = I
%       N - n x m, N(i,j) approx N(0,1)
%       xi- scalar, for 1 <= xi <= sqrt(d/m) signal recoverable, noise does not dominate
%
% 2. Compute SVD of A with svds
%
% 3. Compute B (sketch of A) with fastFrequentDirections
%       requires fastFD.m which also contains fastRankReduce.m
%
% 4. Calculate metrics for comparison
%       a) sines of left & right sing. vecs
%       b) projection error of left & right sing. vecs
%          for k = 1 & k = 1:10:
%               || Uk - Us*Us'*Uk ||,   Uk the left sing. vecs of A,
%                                       Us the left sing. vecs of svd(A*Vb) ([~,~,Vb]=svd(B))
%               || Vk - Vb*Vb'*Vk ||,   Vk the right sing. vecs of A,
%                                       Vb the right sing. vecs of B
%       c) timing data
%
% 5. Output in format '<primme/ffd>_d###_m###_s###_#.mat' where the last values
%       before .mat refers to the tolerance passed to primme
%
%% TEST ONLY VALS
n = 1e3;
dd = [300 400 500];
%dd = 300;
mm = [10 20 50];
%mm = 10;
ells = [10:10:200];

%% MODEL RUN
% n = 1e5; 1e6?
% dd = [1e4 5e4 1e5];      % test columns
% mm = [10 20 50];         % rank of signal


%% 
opts = struct();
opts.disp = 1;
l = zeros(10,1);
r = zeros(10,1);

singvecs = struct();
singvecs.rows = n;
singvecs.cols = dd;
singvecs.ranks = mm;
singvecs.ells = ells;

projerr = struct();
projerr.rows = n;
projerr.cols = dd;
projerr.ranks = mm;
projerr.ells = ells;

time = struct();
time.rows = n;
time.cols = dd;
time.ranks = mm;
time.ells = ells;

for w = 1:numel(dd)          % iter over col sizes d
d = dd(w);
    
    for x = 1:numel(mm)      % iter over rank sizes m
	m = mm(x);

        % MATRIX GEN A = SDU + N/zeta
        sprintf('Setting up the problem: N = %d, d = %d, m = %d.', n, d, m);
        % S
        A = normrnd(0,1,n,m);
        % D
        D = zeros(m,1);
        for i = 1:m
            D(i) = 1-(i-1)/m;
        end
        D = diag(D);
        % U
        Q = normrnd(0,1,d,m);
        [Q,~] = qr(Q,0);

        A = A*(D*Q');   % Assuming m << d
        clear D Q;
        save('tmpA.mat','A','-v7.3'); % noiseless A
        
        % Add noise
        zetas = [0.5 (1+sqrt(d/m))/2 sqrt(d/m) 10*sqrt(d/m) ];
        
        for y = 1:numel(zetas)
        zeta = zetas(y);
            
            if y > 1
                load('tmpA.mat'); % fetch A after 1st iter
            end
            for i = 1:d
                A(:,i) = A(:,i) + normrnd(0,1,n,1)/zeta;
            end

            % PRIMME or SVDS
            opts.tol = 1e-4;
            opts.maxit = d;	% Do not iterate more than the full rank of A (d)
            tic
            [u,s,v] = svds(A,m+2,'L',opts);

            % Restrict everything to the rank m we try to recover
            u = u(:,1:m); s = s(1:m,1:m); v = v(:,1:m);
            time(w,x,y).svd = toc;

            % Compute Em = ||A-UU'A||_F the exact projection error. U are the "exact" singular vectors
            k = m; 	       % Note u has m columns. 
            P = u(:,1:k)'*A;   % this is k x d. We don't want to build uk*P because is of size(A)
            fnorm = 0;
            for i=1:n	   % compute the Frobenius norm row by row
               fnorm = fnorm + norm(A(i,:)-u(i,1:k)*P,'fro')^2;
            end
            fnorm = sqrt(fnorm);

            sprintf('Run method for various ells and evaluate its efficacy');
            	% ell = 5*m;    % FFD's bound m/(ell-m) => want ell about 5 times as the intrinsic model rank m
            for z = 1:length(ells)
                ell = ells(z);

                % FFD
                tic;
                B = fastFD(A,ell); % ALT: Randomized Method Here
                time(w,x,y,z).approx = toc;
                [u1,s1,v1]=svd(B,'econ');
                [u2,s2,v2]=svd(A*v1,'econ');
                % the Sketch and its SVD are then Sketch = A*v1*v1' = u2 * s2 * (v1*v2)^T

                % Evaluations

                % 1. Principal angles between the m sing.vectors of A and the ell-dim space of the sketch
                
                singvecs(w,x,y,z).left = sin(acos(min(1-eps,svd(u2'*u))));
                singvecs(w,x,y,z).right = sin(acos(min(1-eps,svd(v1'*v))));

               

                % 2. LEk = ||A-us*us'A||_F and REk = ||A-Avs*vs'|| left and right projection errors for us,vs from FFD 
                % IMPORTANT: 
                %     us,vs have ell > m columns. We consider only the m largest of them for our projection error
                %     It is possible that for k > m, the sketch projection error is less than the exact projection error for m.

                k = min(m,size(u2,2)); 
                P = u2(:,1:k)'*A;   % this is k x d. We don't want to build uk*P because is of size(A)
                Q = A*v1(:,1:k);    % this is n x k. We don't want to build Q*vk' because is of size(A)
                lpe = 0;
                rpe = 0;
                for i=1:n	   % compute the Frobenius norm row by row
                  lpe = lpe + norm(A(i,:)-u2(i,1:k)*P,'fro')^2;
                  rpe = rpe + norm(A(i,:)-Q(i,:)*v1(:,1:k)','fro')^2;
                end
                projerr(w,x,y,z).left = sqrt(lpe)/fnorm;
                projerr(w,x,y,z).right = sqrt(rpe)/fnorm;
                clear Q;
            end
        end
    end
end
save('testOut.mat','n','dd','mm','zetas','singvecs','projerr','-v7.3');
save('timeOut.mat','time','-v7.3');
clear time;
