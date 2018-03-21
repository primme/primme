function M = makeRSVDmodel(m,n,S,mat_filename,save_matrix)
% makes a matrix of given dimensions and with singular values 
% specified by the diagonal matrix S, writes result to filename
    fprintf('making matrix..\n');
    p = min(m,n);
    if m >= n
       [U, temp] = qr(randn(m,n),0);
       [V, temp] = qr(randn(n));
    else
       [U, temp] = qr(randn(m));
       [V, temp] = qr(randn(n,m),0);
    end

%     whos U S V
    M = U*S*V';

    if save_matrix
        fprintf('writing to mat file..\n');
        save(mat_filename,'M','S');
    end
end
