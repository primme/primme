% TEST DRIVER for Model I from the RSVDPack paper
% our goal: compute m desired number of singular vectors
% 1. Generate test matrix A = UDV'
%       where
%       S - n x m, orthogonal matrix S*S' = I
%       D - m x m, diagonal gives linearly diminishing signal singular values
%       U - m x d, orthogonal matrix U*U' = I
%
% 2. Run four methods
%	SVDS, PSVDS, RSVD, FD (calling fastFD.m)
%
% 4. Calculate metrics for comparison
%	 a) for all methods compute residual norms
%    b) compute sines of left & right sing. vecs (as computed by SVDS) with the spaces returned by RSVD and FD
%    c) projection error of left & right sing. vecs ||A-us*us'A||_F/||A-uu'A|| where u are the m exact sing.vecs
%    d) timing data
%
% 5. Save the data for each method in a separate structure/file
%    We can parallelize the execution for various methods, by having multiple drivers, each set with a different
%    DO_method = 1, and run them in parallel.
%
% 6. Similarly, copy and paste this driverIV and generate the rest of the models.

clear,clc
% Decide which methods to run
DO_PSVDS = 1;
DO_RSVD  = 1;
DO_FD    = 1;
DO_SVDS  = 1;  % SVDS info will always be run or read from disk. if DO_SVDS==0, then the info struct is not written to disk.

% TEST parameters for model IV
nn = [25e3 50e3]; % this size is a reasonable starting matrix size
% nn = [2e3]; % testing DELETE
dd = nn/5;     % We fix dd = 1/5 nn.
% For each (n,d) we explore three ranks, mm(rank,case)=[10, sqrt(d), d/20]; Limit max rank (and evals we compute) to 200
mm = floor([10*ones(size(dd)); sqrt(dd); min(dd/20,200)]); 

% PRIMME & SVDS pars
psvdsopts = struct();
psvdsopts.tol = 1e-2;
% psvdsopts.disp = 3;

% SVDS pars
svdsopts = struct();
svdsopts.tol = 1e-4;

% Four super structures for four different methods
% Frequent Directions------------------
if (DO_FD)
  FD = struct();
  FD(1).rows = nn;
  FD(1).cols = dd; 	% dd is not needed since = nn/5, but for generality
  FD(1).ranks = mm;
  % Other fields:
  % FD(w,x,y,z).sines.left
  % FD(w,x,y,z).sines.right
  % FD(w,x,y,z).ProjErr.left
  % FD(w,x,y,z).ProjErr.right
  % FD(w,x,y,z).residuals
  % FD(w,x,y,z).time
  % FD.zetas = zetas;   % 1 element, a matrix
  % FD.ells = ells;     % 1 element, a vector
end

% Randomized SVD------------------
if (DO_RSVD)
  RSVD = struct();
  RSVD(1).rows = nn;
  RSVD(1).cols = dd;
  RSVD(1).ranks = mm;
  % Other fields
  % RSVD(w,x,y).sines.left
  % RSVD(w,x,y).sines.right
  % RSVD(w,x,y).ProjErr.left
  % RSVD(w,x,y).ProjErr.right
  % RSVD(w,x,y).residuals
  % RSVD(w,x,y).time
  % RSVD.zetas = zetas;   % 1 element, a matrix
  % RSVD.ells = ells;     % 1 element, a vector
end

% MATLAB's SVDS------------------
if (DO_SVDS)
  SVDS = struct();
  SVDS(1).rows = nn;
  SVDS(1).cols = dd;
  SVDS(1).ranks = mm;
  % Other fields
  % SVDS(w,x,y).residuals
  % SVDS(w,x,y).time
  % SVDS.zetas = zetas;   % 1 element, a matrix
  % SVDS.ells = ells;     % 1 element, a vector
end

% PRIMME_SVDS------------------
if (DO_PSVDS)
  PSVDS = struct();
  PSVDS(1).rows = nn;
  PSVDS(1).cols = dd;
  PSVDS(1).ranks = mm;
  % Other fields  
  % PSVDS(w,x,y).sines.left
  % PSVDS(w,x,y).sines.right
  % PSVDS(w,x,y).ProjErr.left
  % PSVDS(w,x,y).ProjErr.right
  % PSVDS(w,x,y).residuals
  % PSVDS(w,x,y).time
  % PSVDS.zetas = zetas;   % 1 element, a matrix
  % PSVDS.ells = ells;     % 1 element, a vector
end
count = 1;
for w = 1:numel(nn) % iterate over different matrix dimensions A(n,d)
   n = nn(w);
   d = n/5;

   psvdsopts.maxit = 2*d; % Do not iterate PRIMME more than d (= max rank(A))

   for x = 1:size(mm,1) % iter over rank sizes of signal: m
        m = mm(x,w);

        fprintf('\n Setting up the problem: N = %d, d = %d, m = %d \n', n, d, m)

        % For this (n,d) and m, set the appropriate ells for FD _only_:
        ells = floor([m : max(m/2,10) : min(5*m,500)]);

        % Building the exp-decay model matrix for various zetas
        zetas = [-0.5 -2 -3.5];
        Zetas(count,:) = zetas;
        count = count + 1;

        for y = 1:numel(zetas)
            zeta = zetas(y);
            p = min(n,d);
            S = logspace(0,zeta,p);
            S = diag(S);
            A = makeRSVDmodel(n,d,S,'',0); % no reading from disk. Generate A with the same seeds.

            %--------------- SVDS --------------%
            % If SVDS has been run before load the U,S,V space for m vectors
            % otherwise, call SVDS(A,m,tol=1e-4)
            disp('Run SVDS method')
            svdinfo = ['svdInfo_',num2str(n),'_',num2str(m),'_',num2str(d),'_',...
                        'z',num2str(y),'.mat'];
            if exist(svdinfo,'file')
                load(svdinfo,'u','s','v');
            else
                tic
                [u,s,v] = svds(A,m,'L',svdsopts);
                SVDS(w,x,y).time = toc;
%                 save(svdinfo,'u','s','v');
                res = zeros(m,1);
                for j=1:m
                    res(j) = norm(A'*u(:,j)-s(j,j)*v(:,j));
                end
                SVDS(w,x,y).residuals = res;
            end
            % Restrict everything to the rank m we try to recover

            % Compute Em = ||A-UU'A||_F the exact projection error. U are the "exact" singular vectors
            k = m; 	       % Note u has m columns.
            P = u(:,1:k)'*A;   % this is k x d. We don't want to build uk*P because is of size(A)
            fnorm = 0;
            for i=1:n % compute the Frobenius norm row by row
               fnorm = fnorm + norm(A(i,:)-u(i,1:k)*P,'fro')^2;
            end
            fnorm = sqrt(fnorm);

            %--------------- PRIMME_SVDS --------------%
            % Call primme_svds(A,m,tol=1e-1) and store its residual norms
            if (DO_PSVDS)
               disp('Run PRIMME_SVDS method')
               tic;
               [up1,sp1,vp1,res] = primme_svds(A,m,'L',psvdsopts);
	       	   PSVDS(w,x,y).time = toc;
               PSVDS(w,x,y).residuals = res'/sqrt(2);
               clear res;
               
               % Evaluations
               % 1. Principal angles between the m sing.vectors of A and the ell-dim space of the sketch
               PSVDS(w,x,y).sine.left = sin(acos(min(1-eps,svd(up1'*u))));
               PSVDS(w,x,y).sine.right = sin(acos(min(1-eps,svd(vp1'*v))));

               % 2. LEk = ||A-us*us'A||_F and REk = ||A-Avs*vs'|| left and right projection errors for us,vs from FFD
               % IMPORTANT:
               %     us,vs have ell > m columns. We consider only the m largest of them for our projection error
               %     It is possible that for k > m, the sketch projection error is less than the exact projection error for m.

               k = min(m,size(up1,2));
               P = up1(:,1:k)'*A;   % this is k x d. We don't want to build uk*P because is of size(A)
               Q = A*vp1(:,1:k);    % this is n x k. We don't want to build Q*vk' because is of size(A)
               lpe = 0;
               rpe = 0;
               for i=1:n % compute the Frobenius norm row by row
                 lpe = lpe + norm(A(i,:)-up1(i,1:k)*P,'fro')^2;
                 rpe = rpe + norm(A(i,:)-Q(i,:)*vp1(:,1:k)','fro')^2;
               end
               PSVDS(w,x,y).ProjErr.left = sqrt(lpe)/fnorm;
               PSVDS(w,x,y).ProjErr.right = sqrt(rpe)/fnorm;
            end
            
            %-------------- RSVD --------------%
            if (DO_RSVD)
               disp('Run RSVD method')
               tic;
               pp = 20; % oversampling to ell
               qq = 2; % the iteration of the power method 
               ss = 1; % orthogonalizations after each multiplication with A
               % Choose one of three variants of RSVD
               [ur1,sr1,vr1] = rsvd_version1(A,m,pp,qq,ss);				   
               RSVD(w,x,y).time = toc;
               len = min(m,size(vr1,2));
               for i=len:-1:1 % implicitly report from large svals to small ones
                   nrm(len-i+1) = norm(A*vr1(:,i) - ur1(:,i)*sr1(i,i));
               end
               RSVD(w,x,y).residuals = nrm;
               clear nrm;

               % Evaluations
               % 1. Principal angles between the m sing.vectors of A and the ell-dim space of the sketch
               RSVD(w,x,y).sine.left = sin(acos(min(1-eps,svd(ur1'*u))));
               RSVD(w,x,y).sine.right = sin(acos(min(1-eps,svd(vr1'*v))));

               % 2. LEk = ||A-us*us'A||_F and REk = ||A-Avs*vs'|| left and right projection errors for us,vs from FFD
               % IMPORTANT:
               %     us,vs have ell > m columns. We consider only the m largest of them for our projection error
               %     It is possible that for k > m, the sketch projection error is less than the exact projection error for m.

               k = min(m,size(ur1,2));
               P = ur1(:,1:k)'*A;   % this is k x d. We don't want to build uk*P because is of size(A)
               Q = A*vr1(:,1:k);    % this is n x k. We don't want to build Q*vk' because is of size(A)
               lpe = 0;
               rpe = 0;
               for i=1:n % compute the Frobenius norm row by row
                 lpe = lpe + norm(A(i,:)-ur1(i,1:k)*P,'fro')^2;
                 rpe = rpe + norm(A(i,:)-Q(i,:)*vr1(:,1:k)','fro')^2;
               end
               RSVD(w,x,y).ProjErr.left = sqrt(lpe)/fnorm;
               RSVD(w,x,y).ProjErr.right = sqrt(rpe)/fnorm;
            end
                
            if (DO_FD)
               fprintf('Run FD method for n = %d, d = %d, m = %d, zeta = %d with various ells and evaluate its efficacy\n',n,d,m,zeta)
               %----------------------------------
               % For various space sizes ell
               for z = 1:length(ells)
                   ell = ells(z);

                   %--------------- FD -----------------%
                   tic;
                   B = fastFD(A,ell);
                   FD(w,x,y,z).time = toc;

                    % Obtain sketch (ASSUMES that B approximates the RIGHT singular space)
                   [u1,s1,v1]=svd(B,'econ');
                   [u2,s2,v2]=svd(A*v1,'econ');
                   vtmp = v1*v2;
                   for i=1:min(m,size(vtmp,2))
                       nrm(i) = norm(A'*u2(:,i) - vtmp(:,i)*s2(i,i));
                   end
                   FD(w,x,y,z).residuals = nrm;
                   clear nrm vtmp;
                   % the Sketch and its SVD are then Sketch = A*v1*v1' = u2 * s2 * (v1*v2)^T

                   % Evaluations
                   FD(w,x,y,z).ells = ells;

                   % 1. Principal angles between the m sing.vectors of A and the ell-dim space of the sketch
                   FD(w,x,y,z).sine.left = sin(acos(min(1-eps,svd(u2'*u))));
                   FD(w,x,y,z).sine.right = sin(acos(min(1-eps,svd(v1'*v))));

                   % 2. LEk = ||A-us*us'A||_F and REk = ||A-Avs*vs'|| left and right projection errors for us,vs from FFD
                   % IMPORTANT:
                   %     us,vs have ell > m columns. We consider only the m largest of them for our projection error
                   %     It is possible that for k > m, the sketch projection error is less than the exact projection error for m.

                   k = min(m,size(u2,2));
                   P = u2(:,1:k)'*A;   % this is k x d. We don't want to build uk*P because is of size(A)
                   Q = A*v1(:,1:k);    % this is n x k. We don't want to build Q*vk' because is of size(A)
                   lpe = 0;
                   rpe = 0;
                   for i=1:n% compute the Frobenius norm row by row
                     lpe = lpe + norm(A(i,:)-u2(i,1:k)*P,'fro')^2;
                     rpe = rpe + norm(A(i,:)-Q(i,:)*v1(:,1:k)','fro')^2;
                   end
                   FD(w,x,y,z).ProjErr.left = sqrt(lpe)/fnorm;
                   FD(w,x,y,z).ProjErr.right = sqrt(rpe)/fnorm;
                end % for z
            end % if DO_FD
            
            clear P Q;
            %------------- Other methods? --------------%                
        end % for y
    end % for x
end % for w

if (DO_SVDS),  SVDS(1).zetas = Zetas;   end;
if (DO_FD),    FD(1).zetas = Zetas;   end;
if (DO_RSVD),  RSVD(1).zetas = Zetas;  end;
if (DO_PSVDS), PSVDS(1).zetas = Zetas;  end;
if (DO_SVDS), save('out_SVDS_I.mat','SVDS','-v7.3'); end;
if (DO_PSVDS), save('out_PSVDS_I.mat','PSVDS','-v7.3');end;
if (DO_RSVD), save('out_RSVD_I.mat','RSVD','-v7.3');end;
if (DO_FD), save('out_FD_I.mat','FD','-v7.3');end;
