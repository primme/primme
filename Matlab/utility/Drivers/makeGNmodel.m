function [ A ] = model(m, n, r, nsignal, rkins, savefilename,save_matrix)
% A simple model for human genotype data matrices in Matlab
%
% Inputs: 
% m: number of rows (gene markers)
% n: number of columns (patients)
% r: number of subblocks to model population admixing
% nsignal: number of entries to represent signal
% rkins: fraction of columns to duplicate

% Outputs: A: a dense matrix of size m x n with matrix elements 0, 1 or 2

    A = zeros(m, n);
    % Model population admixing by randomly setting a subblock to the same 
    % value, k
    for i = 1:r
        k = randi(3)-1;
        r1_lower = randi(m);
        r1_upper = randi(m);
        if r1_lower > r1_upper
            temp = r1_lower;
            r1_lower = r1_upper;
            r1_upper = temp;
        end
        r2_lower = randi(n);
        r2_upper = randi(n);
        if r2_lower > r2_upper
            temp = r2_lower;
            r2_lower = r2_upper;
            r2_upper = temp;
        end
        A(r1_lower:r1_upper, r2_lower:r2_upper) = k;
    end
        
    % Model signal
    for i = 1:nsignal
        A(randi(m),randi(n)) = randi(3)-1;
    end
        
    % Model kinship by duplicating rows
    nkins = floor(rkins*m);
    for i = 1:nkins
       A(randi(m),:) = A(randi(m),:);
    end
    
    if save_matrix == 1
        save(savefilename,'A');
    end
end
