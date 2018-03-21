% randomized low rank svd using eigendecomp of BBt version 
function [U,Sigma,V] = rsvd_version1(A,k,p,q,s)
    m = size(A,1);
    n = size(A,2);
    l = k + p;

    R = randn(n,l);
    Y = A*R; % m \times n * n \times l = m \times l

    for j=1:q
        if mod(2*j-2,s) == 0
            [Y,~] = qr(Y,0);
        end
        Z = A'*Y;

        if mod(2*j-1,s) == 0
            [Z,~] = qr(Z,0);
        end
        Y = A*Z;
    end    
    [Q,~] = qr(Y,0);

    B = Q'*A;

    BBt = B*B'; 
    BBt=0.5*(BBt+BBt'); % make sure it's symmetric

    [Uhat,D] = eig(BBt);
    Sigma = sqrt(D);

    U = Q*Uhat;

    V = zeros(n,l);
    for j=1:l
        v_j = 1/Sigma(j,j) * (B' * Uhat(:,j));
        V(:,j) = v_j;
    end

    % take last k vectors and values, due to eig routine value ordering
    U = U(:,(end-k+1):end);
    Sigma = Sigma((end-k+1):end,(end-k+1):end);
    V = V(:,(end-k+1):end);
end

