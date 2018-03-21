function A = makeFDmodel(varargin)
%
%   A = makeFDmodel(n,m,d,zeta,...)
%
% Generate test matrix A = SDU + N/zeta, where
%       S - n x m, S(i,j) \sim N(0,1)
%       D - m x m, diagonal with D_ii =  1-(i-1)/m
%       U - m x d, signal row space U*U' = I
%       N - n x m, N(i,j) \sim N(0,1)
%       zeta- scalar, zeta <= 1, the full rank noise dominates 
%                 1 < zeta <= sqrt(d/m), signal noisy but recoverable
%                     zeta > sqrt(d/m), A is a strongly low rank matrix
% Input:
% n     column size
% m     the rank of the signal S
% d     the number of columns
% zeta  the scaling of the noise
%   The following arguments are optional
% readWrite: determines how B=S*D*U is generated:
%       readWrite == 0 B is generated (default)
%       readWrite == 1 B is read from a file tmpA.mat 
%       otherwise      B is generated and saved in tmpA.mat
% signalSeed 
%       the seed used to generate the signal S
% noiseSeed 
%       the seed used to generate the noise N

% Check and read input
narginchk(4,7);
n = varargin{1};
m = varargin{2};
d = varargin{3};
zeta = varargin{4};
if nargin >4 
   readWrite = varargin{5};
else
   readWrite = 0;
end

if readWrite == 1
   load('tmpA.mat'); % fetch noiseless signal A
   if (n~=size(A,1) || d~=size(A,2))
      fprintf('Matrix in tmpA.mat has dimensions (%d,%d) not (%d,%d)\n',size(A,1),size(A,2),n,d); 
      return;
   end
else
   % MATRIX GEN A = SDU + N/zeta
   % S
   if nargin >5 
      signalSeed = varargin{6};
      oldSeed = rng;
      rng(signalSeed);
   end
   A = normrnd(0,1,n,m);
   % D
   D = zeros(m,1);
   for i = 1:m
      D(i) = 1-(i-1)/m;
   end
   D = diag(D);
   % U
   U = normrnd(0,1,d,m);
   [U,~] = qr(U,0);
   
   A = A*(D*U');   % Assuming m << d
   clear D U;
   if nargin >5 
      rng(oldSeed);
   end

   if (readWrite ~=0)
      save('tmpA.mat','A','-v7.3'); % save noiseless A
   end
end

% Add noise at level 1/zeta
if nargin > 6 
   noiseSeed = varargin{7};
   oldSeed = rng;
   rng(noiseSeed); % Set the seed so different zetas can be applied to the same noise
end
for i = 1:d
   A(:,i) = A(:,i) + normrnd(0,1,n,1)/zeta;
end
if nargin > 6 
   rng(oldSeed); % Revert to original seed
end
