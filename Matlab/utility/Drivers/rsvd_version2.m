% randomized low rank SVD using QR of B^T version 
function [U,Sigma,V] = rsvd_version2(A,k,p,q,s)
    m = size(A,1);
    n = size(A,2);
    l = k + p;

    R = randn(n,l);
    Y = A*R; % m \times n * n \times k = m \times k

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


    %B = Q'*A; % l \times m * m \times n = l \times n
    %Bt = B'; % n \times l
    Bt = A'*Q;

    [Qhat,Rhat] = qr(Bt,0);

    % Rhat is l \times l
    whos Qhat Rhat

    [Uhat,Sigmahat,Vhat] = svd(Rhat);

    U = Q*Vhat;
    Sigma = Sigmahat;
    V = Qhat*Uhat;
    
    % take first k components
    U = U(:,1:k);
    Sigma = Sigma(1:k,1:k);
    V = V(:,1:k);
end

