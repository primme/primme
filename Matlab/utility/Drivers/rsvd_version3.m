% randomized low rank SVD via QB algorithm 
function [U,Sigma,V] = rsvd_version3(A,k,kstep,q,s)
    m = size(A,1);
    n = size(A,2);
    nstep = round(k/kstep) + 1;
    kval = kstep*nstep;

    while kval > min(m,n)
        kstep = kstep - 1;
        kval = kstep*nstep;
    end

    [Q,B] = randpbQB(A,q,s,kstep,nstep);

    BBt = B*B';
    BBt=0.5*(BBt+BBt'); % make sure it's symmetric
    [Uhat,D] = eig(BBt);

    Sigma = sqrt(D);
    U = Q*Uhat;

    V = zeros(n,kval);
    for j=1:kval
        v_j = 1/Sigma(j,j) * (B' * Uhat(:,j));
        V(:,j) = v_j;
    end

    % take last k vectors and values, due to eig routine value ordering
    U = U(:,(end-k+1):end);
    Sigma = Sigma((end-k+1):end,(end-k+1):end);
    V = V(:,(end-k+1):end);
end

